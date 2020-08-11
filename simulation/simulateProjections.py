"""
Marten Chaillet's cryoET simulator - updated version of the simulator used in SHREC2019 and SHREC2020
Contributions from Ilja Gubins and Gijs van der Schot.
Original simulator (used in SHREC2019) written by Gijs van der Schot, which was loosely based on the simulator in the
TOM toolbox for matlab.
"""

# IO related modules
from pytom.gui.mrcOperations import *
from pytom.gui.guiFunctions import loadstar, datatype
import pytom.tompy.io
import configparser
import logging
import os
import datetime
import sys
from tqdm import tqdm
import time

# Plotting
import matplotlib
# use Qt5Agg to prevent conflict with tkinter in pylab.
matplotlib.use('Qt5Agg')
from pylab import *

# math
from pytom.basic.files import *
import numpy as xp
import random
import constant_dictionaries as phys


class ConfigLogger(object):
    """
    Facilitates writing the conf file to a .log file in the output_folder for reference of settings.
    """
    def __init__(self, log):
        self.__log = log

    def __call__(self, config):
        self.__log.info("Config:")
        config.write(self)

    def write(self, data):
        # stripping the data makes the output nicer and avoids empty lines
        line = data.strip()
        self.__log.info(line)


def create_gold_marker(diameter, voxel_size, solvent_potential):
    """
    From Rahman 2018 (International Journal of Biosensors and Bioelectronics).
    Volume of unit cell gold is 0.0679 nm^3 with 4 atoms per unit cell.
    Volume of gold bead is 4/3 pi r^3.
    """
    from pytom.tompy.tools import create_sphere
    from pytom.tompy.filter import gaussian3d

    # constants
    unit_cell_volume = 0.0679 # nm^3
    atoms_per_unit_cell = 4
    C = 2 * xp.pi * phys.constants['h_bar']**2 / (phys.constants['el'] * phys.constants['me']) * 1E20  # A^2

    # dimension of gold box, always add 5 nm to the sides
    dimension = int(xp.ceil(diameter / voxel_size * 1E-9)) * 3
    sphere = create_sphere((dimension,)*3, radius=(diameter*0.5)*1E-9/voxel_size)
    rounded_sphere = gaussian3d(sphere, sigma=0.3)

    # values transformed to occupied volume per voxel from 1 nm**3 per voxel to actual voxel size
    solvent_correction = rounded_sphere * solvent_potential
    gold_atoms = (rounded_sphere / unit_cell_volume) * atoms_per_unit_cell

    # interaction potential
    gold_scattering_factors = xp.array(phys.scattering_factors['AU']['g'])
    gold_potential = gold_atoms * gold_scattering_factors[0:5].sum() * C / 1000 # 1000 A^3 = 1 nm^3

    return gold_potential - solvent_correction


def generate_model(particleFolder, output_folder, model_ID, listpdbs, pixelSize = 1, size=1024, thickness=200,
                   solvent_potential=4.5301, sigma_structural=0.2, numberOfParticles=1000, placement_size=512, retries=5000,
                   add_gold_markers=True, number_of_markers=20):
    # IMPORTANT: We assume the particle models are in the desired voxel spacing for the pixel size of the simulation!

    from voltools import transform

    # outputs
    X, Y, Z = size, size, thickness
    cell = xp.zeros((X, Y, Z))

    occupancy_bbox_mask = xp.zeros_like(cell)
    occupancy_accurate_mask = xp.zeros_like(cell)
    class_bbox_mask = xp.zeros_like(cell)
    class_accurate_mask = xp.zeros_like(cell)
    ground_truth_txt_file = ''

    # load pdb volumes and pad them
    volumes = []
    for pdb in listpdbs:
        try:
            # First find the voxel size of the map...
            # rotate
            # then bin
            vol = pytom.tompy.io.read_mrc(f'{particleFolder}/{pdb}_{pixelSize*1E10:.2f}A_solvent-{solvent_potential:.3f}V.mrc')
            dx, dy, dz = vol.shape
            vol2 = xp.zeros((dx*2, dy*2, dz*2), dtype=xp.float32)
            vol2[dx//2:-dx//2, dy//2:-dy//2, dz//2:-dz//2] = vol
            volumes.append(vol2)
        except Exception as ee:
            print(ee)
            raise Exception('Could not open pdb ', pdb)

    # attributes
    number_of_classes = len(listpdbs)
    save_path = f'{output_folder}/model_{model_ID}'
    dims = [v.shape for v in volumes] # TODO solve this!!!!!!!!!
    particles_by_class = [0, ] * number_of_classes
    particle_nr = 1
    default_tries_left = retries
    skipped_particles = 0

    difference = size - placement_size
    if not difference:
        loc_x_start = loc_y_start = 0
        loc_x_end = loc_y_end = size
    else:
        loc_x_start = loc_y_start = difference // 2
        loc_x_end = loc_y_end = int(size - xp.ceil(difference / 2))

    # GOLD MARKERS WILL ALSO BE COUNTER TOWARDS TOTAL PARTICLE NUMBER
    # There are also added as an additional class
    if add_gold_markers:
        number_of_classes += 1
        particles_by_class += [0]

        # select a gold particle of either 5, 10, or 15 nm diameter
        size_list = [5,10,15] # in nm

        # class id of gold markers is always the same
        cls_id = number_of_classes - 1

        for _ in tqdm(range(number_of_markers), desc='Placing gold markers'):
            # select a random size for the gold marker
            marker_size = size_list[xp.random.randint(0,len(size_list))]
            # create the gold marker in other function
            gold_marker = create_gold_marker(marker_size, pixelSize, solvent_potential)
            dimensions = gold_marker.shape

            threshold = 0.01
            gold_marker[gold_marker < threshold] = 0

            u = xp.random.uniform(0.0, 1.0, (2,))
            theta = xp.arccos(2 * u[0] - 1)
            phi = 2 * xp.pi * u[1]
            psi = xp.random.uniform(0.0, 2 * xp.pi)
            p_angles = xp.rad2deg([theta, phi, psi])

            # randomly mirror and rotate particle
            try:
                gold_marker = transform(gold_marker, rotation=p_angles,
                                        rotation_order='szxz', interpolation='filt_bspline', device='cpu')
            except Exception as e:
                print(e)
                print('Something went wrong while rotating?')
                continue

            threshold = 0.001
            gold_marker[gold_marker < threshold] = 0

            # thresholded particle
            accurate_particle_occupancy = gold_marker > 0

            # find random location for the particle
            xx, yy, zz = gold_marker.shape
            tries_left = default_tries_left
            while tries_left > 0:
                loc_x = xp.random.randint(loc_x_start + xx // 2 + 1, loc_x_end - xx // 2 - 1)
                loc_y = xp.random.randint(loc_y_start + yy // 2 + 1, loc_y_end - yy // 2 - 1)
                loc_z = xp.random.randint(zz // 2 + 1, Z - zz // 2 - 1)

                tries_left -= 1

                # calculate coordinates of bbox for the newly rotated particle
                bbox_x = [loc_x - dimensions[0] // 2, loc_x + dimensions[0] // 2 + dimensions[0] % 2]
                bbox_y = [loc_y - dimensions[1] // 2, loc_y + dimensions[1] // 2 + dimensions[1] % 2]
                bbox_z = [loc_z - dimensions[2] // 2, loc_z + dimensions[2] // 2 + dimensions[2] % 2]

                # create masked occupancy mask
                masked_occupancy_mask = occupancy_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1],
                                        bbox_z[0]:bbox_z[1]]
                masked_occupancy_mask = masked_occupancy_mask * accurate_particle_occupancy

                # if the location fits (masked occupancy pixel-wise mask is empty), break the loop and use this location
                if masked_occupancy_mask.sum() == 0:
                    break

            # however if still can't fit, ignore this particle (also adds variance in how many particles are
            # actually put)
            if tries_left < 1:
                skipped_particles += 1
                continue

            # populate occupancy volumes
            occupancy_bbox_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] = particle_nr
            occupancy_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
                accurate_particle_occupancy * particle_nr

            # populate class masks
            class_bbox_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] = (cls_id + 1)
            class_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
                accurate_particle_occupancy * (cls_id + 1)

            # populate density volume
            cell[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += gold_marker

            # update stats
            particle_nr += 1
            particles_by_class[cls_id] += 1

            # update text
            ground_truth_txt_file += f'fiducial {int(loc_x - loc_x_start)} {int(loc_y - loc_y_start)} {int(loc_z)}' \
                                     f'NaN NaN NaN\n'

    for _ in tqdm(range(numberOfParticles), desc='Placing particles'):

        # select random class but correct for artifact class if adding gold particles
        cls_id = xp.random.randint(0, number_of_classes - 1) if add_gold_markers else xp.random.randint(0, number_of_classes)

        # generate random rotation
        # to be properly random, use uniform sphere sampling
        # https://math.stackexchange.com/a/442423/72032
        # https://en.wikipedia.org/wiki/Rotation_matrix#Uniform_random_rotation_matrices
        # http://corysimon.github.io/articles/uniformdistn-on-sphere/
        u = xp.random.uniform(0.0, 1.0, (2,))
        theta = xp.arccos(2 * u[0] - 1)
        phi = 2 * xp.pi * u[1]
        psi = xp.random.uniform(0.0, 2 * xp.pi)
        p_angles = xp.rad2deg([theta, phi, psi])

        # randomly mirror and rotate particle
        try:
            particle = volumes[cls_id]
            if xp.random.randint(2): # Generate true/false randomly
                # Mirror the particle to cover both left and right handedness of the proteins
                particle = xp.flip(particle, axis=xp.random.randint(3))
            rotated_particle = transform(particle, rotation=p_angles,
                                         rotation_order='szxz', interpolation='filt_bspline', device='cpu')
        except Exception as e:
            print(e)
            print('Something went wrong while rotating?')
            continue

        # remove particle rotation artifacts
        # threshold = min(volumes[cls_id][volumes[cls_id] > 0]) / 10
        # the approach above doesn't work well with PDB particles, there are often voxels with values of ^-11
        threshold = 0.001
        rotated_particle[rotated_particle < threshold] = 0

        # thresholded rotated particle
        accurate_particle_occupancy = rotated_particle > 0

        # find random location for the particle
        xx, yy, zz = rotated_particle.shape
        tries_left = default_tries_left
        while tries_left > 0:
            loc_x = xp.random.randint(loc_x_start + xx // 2 + 1, loc_x_end - xx // 2 - 1)
            loc_y = xp.random.randint(loc_y_start + yy // 2 + 1, loc_y_end - yy // 2 - 1)
            loc_z = xp.random.randint(zz // 2 + 1, Z - zz // 2 - 1)

            tries_left -= 1

            # calculate coordinates of bbox for the newly rotated particle
            bbox_x = [loc_x - dims[cls_id][0] // 2, loc_x + dims[cls_id][0] // 2 + dims[cls_id][0] % 2]
            bbox_y = [loc_y - dims[cls_id][1] // 2, loc_y + dims[cls_id][1] // 2 + dims[cls_id][1] % 2]
            bbox_z = [loc_z - dims[cls_id][2] // 2, loc_z + dims[cls_id][2] // 2 + dims[cls_id][2] % 2]

            # create masked occupancy mask
            masked_occupancy_mask = occupancy_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]]
            masked_occupancy_mask = masked_occupancy_mask * accurate_particle_occupancy

            # if the location fits (masked occupancy pixel-wise mask is empty), break the loop and use this location
            if masked_occupancy_mask.sum() == 0:
                break

        # however if still can't fit, ignore this particle (also adds variance in how many particles are actually put)
        if tries_left < 1:
            skipped_particles += 1
            continue

        # populate occupancy volumes
        occupancy_bbox_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] = particle_nr
        occupancy_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
            accurate_particle_occupancy * particle_nr

        # populate class masks
        class_bbox_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] = (cls_id + 1)
        class_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
            accurate_particle_occupancy * (cls_id + 1)

        # populate density volume
        cell[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += rotated_particle

        # update stats
        particle_nr += 1
        particles_by_class[cls_id] += 1

        # update text
        ground_truth_txt_file += f'{listpdbs[cls_id]} {int(loc_x - loc_x_start)} {int(loc_y - loc_y_start)} {int(loc_z)} ' \
                                 f'{p_angles[0]:.4f} {p_angles[1]:.4f} {p_angles[2]:.4f}\n'

    # add solvent background potential
    print('Adding background solvent potential')
    cell += solvent_potential

    # Add structural nois
    print('Adding structural noise to grand model cell')
    noisy_cell = cell + xp.random.normal(0, sigma_structural, cell.shape)

    # save grandmodels
    print('Saving grandmodels')
    pytom.tompy.io.write(f'{save_path}/grandmodel_noisefree_original.mrc', cell)
    pytom.tompy.io.write(f'{save_path}/grandmodel_original.mrc', noisy_cell)

    # save class masks
    print('Saving class volumes')
    pytom.tompy.io.write(f'{save_path}/class_bbox.mrc', class_bbox_mask)
    pytom.tompy.io.write(f'{save_path}/class_mask.mrc', class_accurate_mask)

    # save occupancy masks
    print('Saving occupancy volumes')
    pytom.tompy.io.write(f'{save_path}/occupancy_bbox.mrc', occupancy_bbox_mask)
    pytom.tompy.io.write(f'{save_path}/occupancy_mask.mrc', occupancy_accurate_mask)

    # save particle text file
    with open(f'{save_path}/particle_locations.txt', 'w') as f:
        f.write(ground_truth_txt_file)

    # reporting
    print(f'Total number of particles in the tomogram: {particle_nr - 1}\n'
          f'Skipped {skipped_particles} particles ({default_tries_left} random location retries)\n'
          f'Particles by class: {particles_by_class}')

    # # debug: inspect all generated volumes
    # import napari
    # with napari.gui_qt():
    #     viewer = napari.Viewer()
    #     viewer.add_image(occupancy_bbox_mask, name='occupancy bbox mask', interpolation='bicubic')
    #     viewer.add_image(occupancy_accurate_mask, name='occupancy accurate mask', interpolation='bicubic')
    #     viewer.add_image(class_bbox_mask, name='class bbox mask', interpolation='bicubic')
    #     viewer.add_image(class_accurate_mask, name='class accurate mask', interpolation='bicubic')
    #     viewer.add_image(cell, name='cell', interpolation='bicubic')

    # trying to reduce intermediate memory usage
    del cell, noisy_cell, class_accurate_mask, class_bbox_mask, occupancy_accurate_mask, occupancy_bbox_mask

    return


def parallel_rotate_model(volume, outname, angle):
    print(f'Starting rotation process for angle {angle}')
    sys.stdout.flush()
    from voltools import transform
    # volume = pytom.tompy.io.read_mrc(filename)
    rotated_volume = transform(volume, rotation=(0, angle, 0), rotation_order='sxyz', interpolation='filt_bspline', device='cpu')
    pytom.tompy.io.write(outname, rotated_volume)
    print(f'Process for angle {angle} is finished ({outname})')
    sys.stdout.flush()
    return True


def create_rotation_model(output_folder, model_ID):

    # Load grandmodel
    save_path = f'{output_folder}/model_{model_ID}'
    if os.path.exists(f'{save_path}/grandmodel_original.mrc'):
        grandcell = pytom.tompy.io.read_mrc(f'{save_path}/grandmodel_original.mrc')
    else:
        raise Exception(f'create_rotation_model expects grandmodel be created before ({save_path}/grandmodel_original.mrc)')

    print(f'Successfully loaded grandcell (shape: {grandcell.shape})')

    size = grandcell.shape[0]
    height = grandcell.shape[2]
    if grandcell.shape[2] >= heightBox:
        raise Exception('Your model is larger than than the box that we will rotate (heightBox parameter)')

    volume = xp.zeros((size, size, heightBox), dtype=xp.float32)
    offset = (heightBox - height) // 2

    if (heightBox - height) % 2:
        volume[:, :, offset + 1:-offset] = grandcell[:, :, :]
    else:
        volume[:, :, offset:-offset] = grandcell[:, :, :]

    dir = f'{save_path}/rotations'
    filename = f'{dir}/rotated_volume_0.mrc'
    pytom.tompy.io.write(filename, volume)

    print(f'Saved initial rotation volume at {filename}')
    return volume


def wavelength_eV2m(V):
    # OLD FUNCTION
    # h / sqrt( 2 * me * el * ev * 1)
    # return 6.625 * 10 ** (-34) / ((2 * 9.1 * 10 ** (-31)) * 1.6 * (10 ** (-19)) * ev * 1) ** 0.5

    # NEW FUNCTION
    # function calculates the electron wavelength given a certain accelaration voltage V
    h = phys.constants["h"]
    e = phys.constants["el"]
    m = phys.constants["me"]
    c = phys.constants["c"]

    # matlab original: lambda = h/sqrt(e*V*m*(e/m*V/c^2 + 2 ));
    Lambda = h/xp.sqrt(e*V*m*(e/m*V/c**2 + 2 ))

    return Lambda


def create_complex_CTF(image_shape, pix_size, Dz, voltage=300E3, Cs=2.7E-3, sigma_decay_CTF=0.4, amplitude_contrast=0.07):
    """
    Adapated from Vulovic et al., 2013. Returns a complex contrast transfer function. Dimensions of input image or
    volume should be equal. Only phase part of CTF.

    @param vol:
    @type vol:
    @param pix_size:
    @type pix_size:
    @param Dz:
    @type Dz:
    @param voltage:
    @type voltage:
    @param Cs:
    @type Cs:
    @param sigma_decay_ctf:
    @type sigma_decay_ctf:
    @param amplitude_contrast:
    @type amplitude_contrast:

    @return:
    @rtype:

    @author: Marten Chaillet
    """

    Lambda = wavelength_eV2m(voltage)

    Ny = 1 / (2 * pix_size)

    if len(image_shape) > 2:
        R, Y, Z = xp.meshgrid(xp.arange(-Ny, Ny, 2 * Ny / image_shape[0]), xp.arange(-Ny, Ny, 2 * Ny / image_shape[1]),
                              xp.arange(-Ny, Ny, 2 * Ny / image_shape[2]))
        r = xp.sqrt(R ** 2 + Y ** 2 + Z ** 2)
    else:
        R, Y = xp.meshgrid(xp.arange(-Ny, Ny, 2. * Ny / (image_shape[0])), xp.arange(-Ny, Ny, 2. * Ny / (image_shape[1])))
        r = xp.sqrt(R ** 2 + Y ** 2)

    # print(r)
    complex_CTF = xp.exp( -1j * xp.pi / 2 * (Cs * (Lambda ** 3) * (r ** 4) - 2 * Dz * Lambda * (r ** 2)) )

    # TODO Should amplitude contrast be in??? Discuss with Friedrich and Gijs
    # ctf.real = ctf.real * amplitude_contrast
    # ctf.imag = ctf.imag * (1-amplitude_contrast)

    if sigma_decay_CTF:
        decay = xp.exp(-(r / (sigma_decay_CTF * Ny)) ** 2)
        complex_CTF = complex_CTF * decay

    return complex_CTF


def transmission_function(sliced_potential, voltage, dz):
    # wavelength
    Lambda = wavelength_eV2m(voltage)
    # relative mass
    relmass = phys.constants["me"] + phys.constants["el"] * voltage / (phys.constants["c"] ** 2)
    # sigma_transfer
    sig_transfer = 2 * xp.pi * relmass * phys.constants["el"] * Lambda / (phys.constants["h"] ** 2)

    return xp.exp(1j * sig_transfer * sliced_potential * dz)


def fresnel_propagator(imageSize, pixelSize, voltage, dz):
    Lambda = wavelength_eV2m(voltage)

    xwm = pixelSize * (imageSize)  # pixelsize for multislice * size sample
    q_true_pix_m = 1 / xwm  # spacing in Fourier space
    x_line = xp.arange(-(imageSize - 1) / 2, (imageSize - 1) / 2 + 1, 1)
    [xv, yv] = xp.meshgrid(x_line, x_line)
    q_m = xp.sqrt(xv**2 + yv**2) * q_true_pix_m  # frequencies in Fourier domain

    return xp.exp(-1j * xp.pi * Lambda * (q_m ** 2) * dz)


def gradient_image(size, factor, angle=0, center_shift=0):
    """
    Creates an image with a gradient of values rotated along angle. Factor determines the strength of the gradient.
    @param size:
    @param factor:
    @param angle:
    @param shift:
    @return:
    """
    from scipy.ndimage import rotate
    max_rotation_radius = (size/2) / xp.cos(45 * xp.pi / 180)
    extension = int(xp.ceil(max_rotation_radius - size/2))
    left = 1 - factor
    right = 1 + factor
    step = (right-left) / size
    values = xp.arange(left - extension * step + center_shift * step,
                       right + extension * step + center_shift * step, step)
    image = xp.repeat(values[xp.newaxis, :], size + 2*extension, axis=0)
    return rotate(image, angle, reshape=False)[extension:size+extension, extension:size+extension]


def microscope(noisefree_projections, angles, ice_thickness=200, dose=80, pixel_size=1E-9,
               binning=1, voltage=300E3, camera_type='K2SUMMIT', camera_folder=''):
    """
    TODO add parameters for camera type and folder with detector data
    Inspired by InSilicoTEM (Vulovic et al., 2013)
    @author: Marten Chaillet
    """
    import detector
    # noisefree_projections is a stack of images
    # Assuming images are same size in x and y
    size = noisefree_projections.shape[0]
    n_images = noisefree_projections.shape[2]

    print('number of projections is ', n_images)

    projections = xp.zeros((size, size, n_images), dtype=xp.float32)

    dqe = detector.create_detector_response(camera_type, 'DQE', xp.zeros((size, size)), voltage=voltage,
                                            folder=camera_folder)
    mtf = detector.create_detector_response(camera_type, 'MTF', xp.zeros((size, size)), voltage=voltage,
                                            folder=camera_folder)
    ntf = xp.sqrt(mtf ** 2 / dqe) # square root because dqe = mtf^2 / ntf^2
    mtf_shift = xp.fft.ifftshift(mtf)
    ntf_shift = xp.fft.ifftshift(ntf)

    # NUMBER OF ELECTRONS PER PIXEL
    dose_per_tilt = dose / n_images  # electrons per square A per tilt
    dose_per_pixel = dose_per_tilt * (pixel_size*1E10)**2 / binning**2 # from square A to square nm (10A pixels)
    print(f'number of electrons per pixel (before binning and absorption): {dose_per_pixel}')

    for n, angle in enumerate(angles):
        projection = noisefree_projections[:, :, n]
        # Fourier transform and multiply with sqrt(dqe) = mtf/ntf
        projection_fourier = xp.fft.fftn(xp.fft.ifftshift(projection))
        projection_fourier = projection_fourier * mtf_shift / ntf_shift
        # Convert back to real space
        projection = xp.real(xp.fft.fftshift(xp.fft.ifftn(projection_fourier)))
        projection[projection<0] = 0
        # Draw from poisson distribution and scale by camera's conversion factor
        # conversion_factor = 100  # in ADU/e- , this is an arbitrary unit. Value taken from Vulovic et al., 2010
        absorption_factor = xp.exp( - (ice_thickness / phys.mean_free_path[voltage]) / xp.cos(angle * xp.pi / 180) )
        print(f'number of electrons per pixel for tilt angle {angle:.2f} degrees (before binning): {dose_per_pixel * absorption_factor:.2f}')
        poisson_mean = projection * ( dose_per_pixel * absorption_factor )
        projection_poisson = xp.zeros(projection.shape)
        for _ in range(binning**2):
            projection_poisson += ( xp.random.poisson(lam=poisson_mean) / binning**2 )
            # imshow(projection_poisson)
            # show()

        # Image values are now in ADU
        # Apply the camera's noise transfer function to the noisy image
        projection_fourier = xp.fft.fftn(xp.fft.ifftshift(projection_poisson)) * ntf_shift

        # Add additional noise from digitization process, less relevant for modern day cameras.
        # readout noise standard deviation can be 7 ADUs, from Vulovic et al., 2010
        # sigma_readout = 7
        # readsim = xp.random.normal(0, sigma_readout, projection.shape) # readout noise has a gaussian distribution
        # darksim = 0     # dark current noise has a poisson distribution, usually an order of magnitude smaller than
        #                 # readout noise and can hence be neglected

        # Add readout noise and dark noise in real space
        projection = xp.real(xp.fft.fftshift(xp.fft.ifftn(projection_fourier))) # + readsim + darksim

        projections[:,:,n] = projection

    return projections


def generate_projections(angles, output_folder, model_ID, image_size= 512, pixel_size=1E-9, binning=1, ice_thickness=200E-9,
                         dose=80, voltage=300E3, spherical_aberration=2.7E-3, multislice=False, msdz=5E-9,
                         amplitudeContrast=0.07, defocus=2E-6, sigma_decay_CTF=0.4, camera_type='K2SUMMIT',
                         camera_folder='', random_gradient=False, random_rotation=False):
    """

    @param angles:
    @type angles:
    @param output_folder:
    @type output_folder:
    @param model_ID:
    @type model_ID
    @param image_size:
    @type image_size:
    @param pixel_size:
    @type pixel_size:
    @param binning:
    @type binning:
    @param ice_thickness:
    @type ice_thickness:
    @param dose:
    @type dose:
    @param voltage:
    @type voltage:
    @param sphericalAberration:
    @type sphericalAberration:
    @param multislice:
    @type multislice
    @param msdz:
    @type msdz:
    @param amplitudeContrast:
    @type amplitudeContrast
    @param defocus:
    @type defocus:
    @param sigmaDecayCTF:
    @type sigmaDecayCTF:
    @param camera_type:
    @type camera_type:
    @param camera_folder:
    @type camera_folder:
    @param random_gradient:
    @type random_gradient:
    @param random_rotation:
    @type random_rotation:

    @return:
    @rtype:

    @author: Marten Chaillet
    """
    from pytom.basic.datatypes import DATATYPE_ALIGNMENT_RESULTS as dar
    from pytom.basic.datatypes import fmtAlignmentResults, HEADER_ALIGNMENT_RESULTS
    from pytom.gui.guiFunctions import savestar
    from scipy.ndimage import rotate
    # NOTE; Parameter defocus specifies the defocus at the bottom of the model!

    # grab model
    save_path = f'{output_folder}/model_{model_ID}'
    grandcell = pytom.tompy.io.read_mrc(f'{save_path}/rotations/rotated_volume_{0}.mrc')

    imageSize = grandcell.shape[0]//2
    box_size = grandcell.shape[0]
    box_height = grandcell.shape[2]

    zheight = box_height * pixel_size  # thickness of rotation volume in nm
    defocus -= (zheight / 2) # center defocus value at tilt angle

    ileft = (box_size - image_size)//2
    iright = -int(xp.ceil((box_size - image_size)/2))

    # create a 3d array for the exit wave of each image
    noisefree_projections = xp.zeros((image_size, image_size, len(angles)))
    # get the contrast transfer function
    complex_CTF = create_complex_CTF((image_size, image_size), pixel_size, defocus, voltage=voltage, Cs=spherical_aberration,
                  sigma_decay_CTF=sigma_decay_CTF, amplitude_contrast=amplitudeContrast)

    # Check if msdz is viable, else correct it
    if msdz % pixel_size != 0:
        # Round the msdz to the closest integer mutltiple of the pixel size
        round_up = xp.round(msdz % pixel_size / pixel_size)
        if round_up:
            msdz += (pixel_size - msdz % pixel_size)
        else:
            msdz -= (msdz % pixel_size)
        print(f'We adjusted msdz to {msdz*1E9:.3f} nm to make it an integer multiple of pixel size.')

    for n, angle in enumerate(angles):

        filename = f'{save_path}/rotations/rotated_volume_{int(angle)}.mrc'
        rotated_volume = pytom.tompy.io.read_mrc(filename)  #* 1E1

        if not multislice:
            print('Simulating projection (without ms) from tilt angle ', angle)
            # remove .get() if fully cupy
            projected_tilt_image = rotated_volume[ileft:iright, ileft:iright, :].sum(axis=2)#.get()
            projected_tilt_image = xp.abs(xp.fft.ifftn(xp.fft.fftn(projected_tilt_image)*complex_CTF))**2
            #TODO is abs()**2 the way to go for exit wave field of simple projection?

        else:
            print('Simulating projection (with ms) from tilt angle ', angle)

            # Determine the number of slices
            if msdz > zheight:
                n_slices = 1
                msdz = zheight
                print('The multislice step can not be larger than volume thickness. Use only one slice.\n')
            elif msdz < pixel_size:
                n_slices = box_height
                msdz = pixel_size
                print('The multislice steps can not be smaller than the pixel size. Use the pixel size instead for the step size.\n')
            else:
                n_slices = int(xp.ceil(xp.around(zheight/msdz,3)))

            print('Number of slices: ', n_slices)
            # Allocate space for multislice projection
            projected_potent_ms = xp.zeros((image_size, image_size, n_slices), dtype=complex)

            px_per_slice = int(msdz/pixel_size)
            num_px_last_slice = box_height % px_per_slice  # consider last slice might have a different number of pixels

            # Project potential for each slice (phase grating)
            for ii in range(n_slices):
                # remove .get() if fully cupy
                projected_potent_ms[:, :, ii] = rotated_volume[ileft:iright, ileft:iright,
                                                ii * px_per_slice: (ii + 1) * px_per_slice].mean(axis=2)  # .get()

            # calculate the transmission function for each slice
            psi_t = transmission_function(projected_potent_ms, voltage, msdz)

            # calculate the fresnel propagator (identical for same dz)
            propagator = fresnel_propagator(image_size, pixel_size, voltage, msdz)

            # Wave propagation with MULTISLICE method
            psi_multislice = xp.zeros((image_size,image_size), dtype=complex) + 1 # should be complex datatype

            # imshow(psi_t[:, :, 0].real)
            # show()
            # imshow(psi_t[:, :, 0].imag)
            # show()
            # imshow(xp.abs(psi_multislice) ** 2)
            # show()

            for ii in range(n_slices-min(1,num_px_last_slice)):
                waveField = xp.fft.fftn( xp.fft.ifftshift(psi_multislice) * xp.fft.ifftshift(psi_t[:, :, ii]) )
                psi_multislice = xp.fft.fftshift( xp.fft.ifftn((waveField * xp.fft.ifftshift(propagator) )) )

                if ii%50 == 0:
                    fig, ax = plt.subplots(1, 3)
                    ax[0].imshow(psi_t[:,:,ii].real, vmin=0, vmax=2)
                    ax[1].imshow(psi_t[:, :, ii].imag, vmin=0, vmax=2)
                    ax[2].imshow(xp.abs(psi_multislice)**2, vmin=0, vmax=2)
                    show()

            # Calculate propagation through last slice in case the last slice contains a different number of pixels
            if num_px_last_slice:
                msdz_end = num_px_last_slice * pixel_size
                psi_t[:, :, -1] = transmission_function(projected_potent_ms[:, :, -1], voltage, msdz_end)
                propagator_end = fresnel_propagator(image_size, pixel_size, voltage, msdz_end)
                waveField = xp.fft.fftn( xp.fft.ifftshift(psi_multislice) * xp.fft.ifftshift(psi_t[:, :, -1]) )
                psi_multislice = xp.fft.fftshift( xp.fft.ifftn( waveField * xp.fft.ifftshift(propagator_end) ) )

            imshow(xp.abs(psi_multislice)**2)
            show()

            # Multiple by CTF for microscope effects on electron wave
            wave_CTF = xp.fft.ifftshift(complex_CTF) * xp.fft.fftn(xp.fft.ifftshift(psi_multislice) )
            # Intensity in image plane is obtained by taking the absolute square of the wave function
            projected_tilt_image = xp.abs(xp.fft.fftshift(xp.fft.ifftn(wave_CTF))) ** 2

            imshow(projected_tilt_image)
            show()

        noisefree_projections[:,:,n] = projected_tilt_image

    # Add random gradient (artifact) to all projected intensities
    if random_gradient:
        random_factor = xp.random.uniform(.0, .2)
        random_angle = xp.random.uniform(-180, 180)
        random_shift = xp.random.uniform(-image_size/2, image_size/2)
        intensity_gradient = gradient_image(image_size, random_factor, angle=random_angle, center_shift=random_shift)
        # Apply the gradient to each projections
        for i in range(len(angles)):
            noisefree_projections[:,:,i] *= intensity_gradient

    # Add random rotations if desired
    if random_rotation:
        random_angles = xp.random.uniform(-2, 2, size=len(angles))
        for i in range(len(angles)):
            noisefree_projections[:, :, i] = rotate(noisefree_projections[:, :, i],
                                                    -random_angles[i], reshape=False)
            # remove rotation artifacts smaller than 0
            noisefree_projections[:, :, i][noisefree_projections[:, :, i] < 0] = 0

    # First save the noisefree projections
    pytom.tompy.io.write(f'{save_path}/projections_noisefree.mrc',
                         noisefree_projections)  # noisefree_projections.astype(xp.float32)?

    # Now add microscope effects to the projections
    projections = microscope(noisefree_projections, angles, ice_thickness=ice_thickness, dose=dose, pixel_size=pixel_size,
                             binning=binning, voltage=voltage, camera_type=camera_type, camera_folder=camera_folder)
    pytom.tompy.io.write(f'{save_path}/projections.mrc', projections)

    # Create a folder for individual noisy projections (needed for reconstruction algorithm)
    if not os.path.exists(f'{save_path}/projections/'):
        os.mkdir(f'{save_path}/projections/')
    for i in range(len(angles)):
        # These files are needed for the reconstruction
        pytom.tompy.io.write(f'{save_path}/projections/synthetic_{i+1}.mrc',
                             projections[:, :, n].squeeze())

    # len(angles) is the number of files that we have
    alignment = xp.zeros(len(angles), dtype=dar)
    alignment['TiltAngle'] = angles
    alignment['Magnification'] = xp.repeat(1.0, len(angles))
    if random_rotation:
        alignment['InPlaneRotation'] = random_angles

    for i in range(len(angles)):
        alignment['FileName'][i] = f'{save_path}/projections/synthetic_{i+1}.mrc'

    # Then the alignment file, because the noisyProjections folder needs to have been created
    alignment_file = f'{save_path}/projections/alignment_simulated.txt'
    savestar(alignment_file, alignment, fmt=fmtAlignmentResults, header=HEADER_ALIGNMENT_RESULTS)
    return


def reconstruct_tomogram(start_idx, end_idx, size_reconstruction, angles, output_folder, model_ID, weighting=-1):
    """
    Reconstruction of simulated tilt series into a tomogram. First creates an alignemnt file, then uses weighted back
    projection to make a reconstruction.

    @param prefix: file name prefix
    @type prefix: string
    @param suffix: file extension
    @type suffix: basestring
    @param start_idx: starting index number of first tilt image in file names
    @type start_idx: int
    @param end_idx: ending index of last tilt image in file names
    @type end_idx: int
    @param vol_size: size of the reconstruction volume
    @type vol_size: [int, int, int]
    @param angles: list of angles of the tilt projections
    @type angles: [float, float, ...]
    @param output_folder: output directory path
    @type output_folder: string
    @param model_ID: number of the model
    @type model_ID: int
    @param weighting: weighting factor, either -1, 0, or 1
    @type weighting: int

    @return: -
    @rtype: None

    @author: Gijs van der Schot, Marten Chaillet
    """
    from pytom.reconstruction.reconstructionStructures import Projection, ProjectionList
    save_path = f'{output_folder}/model_{model_ID}'

    alignment_file = f'{save_path}/projections/alignment_simulated.txt'
    # prefix = f'{save_path}/projections/simu_'
    # suffix = '.mrc'
    vol_size = [size_reconstruction, ] * 3

    projections = ProjectionList()
    # IMPORTANT: angles *-1 to get the right reconstrunction relative to the orignal model!
    for i in range(start_idx, end_idx+1):
        # p = Projection(prefix+str(i)+suffix, tiltAngle= -1 * angles[i-1])
        p = Projection(f'{save_path}/projections/synthetic_{i+1}.mrc', tiltAngle=-1 * angles[i - 1])
        projections.append(p)

    outputname = f'{save_path}/reconstruction.em'

    # IF EM alignment file provided, filters applied and reconstruction will be identical.
    vol = projections.reconstructVolume(dims=vol_size, reconstructionPosition=[0,0,0], binning=1,
                                        applyWeighting=weighting, alignResultFile=alignment_file)
    vol.write(outputname)
    os.system(f'em2mrc.py -f {outputname} -t {os.path.dirname(outputname)}')
    os.system(f'rm {outputname}')

    # Crop the models to the ouput reconstruction for training data annotation
    cell = pytom.tompy.io.read_mrc(f'{save_path}/grandmodel_original.mrc')
    reconstruction = pytom.tompy.io.read_mrc(f'{save_path}/reconstruction.mrc')

    difference = cell.shape[0] - reconstruction.shape[0]
    del reconstruction

    if difference:
        ileft = difference // 2
        iright = - int(xp.ceil(difference/2))
        print('-- Cropping models')
        pytom.tompy.io.write(f'{save_path}/grandmodel.mrc', cell[ileft:iright, ileft:iright, :])
        cell = pytom.tompy.io.read_mrc(f'{save_path}/grandmodel_noisefree_original.mrc')
        pytom.tompy.io.write(f'{save_path}/grandmodel_noisefree.mrc', cell[ileft:iright, ileft:iright, :])
        cell = pytom.tompy.io.read_mrc(f'{save_path}/class_bbox.mrc')
        pytom.tompy.io.write(f'{save_path}/class_bbox.mrc', cell[ileft:iright, ileft:iright, :])
        cell = pytom.tompy.io.read_mrc(f'{save_path}/class_mask.mrc')
        pytom.tompy.io.write(f'{save_path}/class_mask.mrc', cell[ileft:iright, ileft:iright, :])
        cell = pytom.tompy.io.read_mrc(f'{save_path}/occupancy_bbox.mrc')
        pytom.tompy.io.write(f'{save_path}/occupancy_bbox.mrc', cell[ileft:iright, ileft:iright, :])
        cell = pytom.tompy.io.read_mrc(f'{save_path}/occupancy_mask.mrc')
        pytom.tompy.io.write(f'{save_path}/occupancy_mask.mrc', cell[ileft:iright, ileft:iright, :])
    else:
        pytom.tompy.io.write(f'{save_path}/grandmodel.mrc', cell[:, :, :])
        cell = pytom.tompy.io.read_mrc(f'{save_path}/grandmodel_noisefree_original.mrc')
        pytom.tompy.io.write(f'{save_path}/grandmodel_noisefree.mrc', cell[:, :, :])

    del cell
    return


if __name__ == '__main__':

    # Read config
    config = configparser.ConfigParser()
    try:
        if len(sys.argv) > 1:
            config_given = sys.argv[1]
            if config_given and os.path.exists(config_given):
                print(f'\nLoading a given configuration file: {config_given}')
                config.read_file(open(config_given))
        else:
            print(f'\nLoading default configuration file: pytom/simulation/simulation.conf')
            config.read_file(open('simulation.conf'))
    except Exception as e:
        print(e)
        raise Exception('Could not open config file.')

    print('Configuration sections:', config.sections())

    # Set simulation parameters
    try:
        output_folder = config['General']['OutputFolder']
        model_ID = int(config['General']['ModelID'])
        seed = int(config['General']['Seed'])
        pixel_size = float(config['General']['PixelSize']) * 1E-9
        metadata = loadstar(config['General']['MetaFile'], dtype=datatype)

        print(f'Generating model {model_ID} in folder {output_folder}')
    except Exception as e:
        print(e)
        raise Exception('Missing general parameters in config file.')

    if 'GenerateModel' in config.sections():
        try:
            # We assume the particle models are in the desired voxel spacing for the pixel size of the simulation!
            particleFolder = config['GenerateModel']['ParticleFolder']
            listpdbs = eval(config['GenerateModel']['Models'])
            placement_size = int(config['GenerateModel']['PlacementSize'])
            size = int(config['GenerateModel']['Size'])
            sigma_structural = float(config['GenerateModel']['SigmaStructuralNoise'])
            add_gold_markers = config['GenerateModel'].getboolean('GoldMarkers')
            if add_gold_markers:
                m_range = number_of_markers = config['GenerateModel']['NumberOfMarkers'].split('-')
                if len(m_range) > 1:
                    xp.random.seed(seed)
                    random.seed(seed)
                    number_of_markers = xp.random.randint(int(m_range[0]), int(m_range[1]))
                else:
                    number_of_markers = int(m_range[0])
            else:
                number_of_markers = 0

            # Randomly vary solvent potential
            solvent_potential = float(config['GenerateModel']['SolventPotential'])
            vary_solvent = config['GenerateModel'].getboolean('VarySolvent')
            if vary_solvent:
                xp.random.seed(seed)
                random.seed(seed)
                factor = xp.random.randint(-1, 4) * 0.1
                solvent_potential = solvent_potential + factor * solvent_potential

            # parse range of ice thickness, provided in nm
            thicc_range = config['GenerateModel']['Thickness'].split('-')
            if len(thicc_range) > 1:
                xp.random.seed(seed)
                random.seed(seed)

                thickness = xp.random.uniform(float(thicc_range[0]), float(thicc_range[1])) * 1E-9
            else:
                thickness = float(thicc_range[0]) * 1E-9
            thickness = int(thickness / pixel_size)
            if thickness % 2 == 1:
                thickness += 1

            # parse range of number of particles
            p_range = config['GenerateModel']['NumberOfParticles'].split('-')
            if len(p_range) > 1:
                xp.random.seed(seed)
                random.seed(seed)

                numberOfParticles = xp.random.randint(int(p_range[0]), int(p_range[1]))
            else:
                numberOfParticles = int(p_range[0])
        except Exception as e:
            print(e)
            raise Exception('Missing generate model parameters in config file.')

    if 'Rotation' in config.sections():
        try:
            heightBox = int(config['Rotation']['HeightBox'])
            nodes = int(config['Rotation']['Nodes'])
        except Exception as e:
            print(e)
            raise Exception('Missing rotation parameters in config file.')

    if 'GenerateProjections' in config.sections():
        try:
            image_size = int(config['GenerateProjections']['ImageSize'])
            multislice = config['GenerateProjections'].getboolean('MultiSlice')
            if multislice:
                msdz = float(config['GenerateProjections']['MultiSliceSize']) * 1E-9 # multislice step size in nm
            else:
                msdz = 0
            sigma_decay_CTF = float(config['GenerateProjections']['SigmaDecayCTF'])
            camera_type = config['GenerateProjections']['Camera']
            camera_folder = config['GenerateProjections']['CameraFolder']
            binning = int(config['GenerateProjections']['Binning'])
            random_gradient = config['GenerateProjections'].getboolean('RandomGradient')
            random_rotation = config['GenerateProjections'].getboolean('RandomRotation')

            # parse range of dose
            dose_range = config['GenerateProjections']['ElectronDose'].split('-') # dose over full tilt series (not per tilt angle)
            if len(dose_range) > 1:
                xp.random.seed(seed)
                random.seed(seed)

                dose = xp.random.randint(int(dose_range[0]), int(dose_range[1]))
            else:
                dose = int(dose_range[0])


            # get all relevant parameters for the simulation from the meta file
            angles = metadata['TiltAngle']  # specified in degrees
            voltage = metadata['Voltage'][0] * 1E3  # voltage in keV
            spherical_aberration = metadata['SphericalAberration'][0] * 1E-3  # spherical aberration in mm
            amplitude_contrast = metadata['AmplitudeContrast'][0]  # fraction of amplitude contrast

            # defocus = 3E-6 # We want to get these from the metafile eventually
            defocus_range = config['GenerateProjections']['Defocus'].split('-') # defocus in um
            if len(defocus_range) > 1:
                xp.random.seed(seed)
                random.seed(seed)

                defocus = xp.random.uniform(float(defocus_range[0]),float(defocus_range[1])) * 1E-6
            else:
                defocus = float(defocus_range[0]) * 1E-6

            # pixelSize = 1E-9 # We want to get this from the metafile eventually
            # pixelSize = metadata['PixelSpacing'][0] * 1E-9 # pixel size in nm
        except Exception as e:
            print(e)
            raise Exception('Missing generate projection parameters.')

    if 'ReconstructTomogram' in config.sections():
        try:
            start = int(config['ReconstructTomogram']['StartIdx'])
            end = int(config['ReconstructTomogram']['EndIdx'])
            weighting = int(config['ReconstructTomogram']['Weighting'])
            sizeRecon = int(config['ReconstructTomogram']['SizeRecon'])
            angles = metadata['TiltAngle']  # specified in degrees
        except Exception as e:
            print(e)
            raise Exception('Missing tomogram reconstruction parameters in config file.')

    # Create directories and logger
    if not os.path.exists(os.path.join(output_folder, f'model_{model_ID}')):
        os.mkdir(os.path.join(output_folder, f'model_{model_ID}'))

    if os.path.exists(f'{output_folder}/model_{model_ID}/simulator.log'):
        os.remove(f'{output_folder}/model_{model_ID}/simulator.log')

    logging.basicConfig(filename='{}/model_{}/simulator-{date:%Y-%m-%d_%H:%M:%S}.log'.format(output_folder, model_ID,
                                                                                             date=datetime.datetime.now()), level=logging.INFO)
    config_logger = ConfigLogger(logging)
    config_logger(config)

    logging.info('Values of parameters that we randomly vary per simulation (only for generate model and generate projections):')
    if 'GenerateModel' in config.sections():
        logging.info(f'ice thickness = {thickness*pixel_size*1E9:.2f}nm')
        logging.info(f'solvent potential = {solvent_potential}V')
        logging.info(f'# of particles = {numberOfParticles}')
    if 'GenerateProjections' in config.sections():
        logging.info(f'defocus = {defocus*1E6:.2f}um')
        logging.info(f'dose = {dose} e-/A^2')

    # Generate a grand model
    if 'GenerateModel' in config.sections():
        # set seed for random number generation
        xp.random.seed(seed)
        random.seed(seed)

        print('\n- Generating grand model')
        generate_model(particleFolder, output_folder, model_ID, listpdbs, pixelSize=pixel_size, size=size,
                       thickness=thickness, placement_size=placement_size, solvent_potential=solvent_potential,
                       numberOfParticles=numberOfParticles, sigma_structural=sigma_structural,
                       add_gold_markers=add_gold_markers, number_of_markers=number_of_markers)

    # Generated rotated grand model versions
    if 'Rotation' in config.sections():
        print(f'\n- Rotating model')

        # If needed, creating rotations directory, removing leftover initial rotation volume
        dir = f'{output_folder}/model_{model_ID}/rotations'
        if not (os.path.exists(dir)):
            os.mkdir(dir)
        filename = f'{dir}/rotated_volume_0.mrc'
        if os.path.exists(filename):
            print(f'Found leftover rotated volume 0, removing it (path: {filename}')
            os.remove(filename)

        #  Create initial rotation volume (0 degree rotation)
        volume = create_rotation_model(output_folder, model_ID)

        # parallel rotate
        sys.stdout.flush()
        time.sleep(5)
        # sleep helps with memory usage
        # my theory (ilja): it gives time for garbage collector to clean up

        from joblib import Parallel, delayed
        verbosity = 55 # set to 55 for debugging, 11 to see progress, 0 to turn off output

        # joblib automatically memory maps a numpy array to child processes
        print(f'Rotating the model with {nodes} processes')
        results = Parallel(n_jobs=nodes, verbose=verbosity, prefer="threads")\
            (delayed(parallel_rotate_model)(volume, f'{dir}/rotated_volume_{int(ang)}.mrc', ang)
             for ang in angles if ang != 0.)

        sys.stdout.flush()
        if all(results):
            print('All rotation processes finished successfully')
        else:
            print(f'{results.count(None)} rotation processes did not finish successfully')

    # Generate noise-free projections
    if 'GenerateProjections' in config.sections():
        # set seed for random number generation
        xp.random.seed(seed)
        random.seed(seed)
        # Grab the ice thickness from the initial model
        ice_thickness = (
                    pytom.tompy.io.read_mrc(f'{output_folder}/model_{model_ID}/grandmodel_original.mrc').shape[2] *
                    pixel_size)

        print('\n- Generating projections')
        generate_projections(angles, output_folder, model_ID, image_size=image_size, pixel_size=pixel_size, binning=binning,
                             ice_thickness=ice_thickness, dose=dose, voltage=voltage, spherical_aberration=spherical_aberration,
                             multislice=multislice, msdz=msdz, amplitudeContrast=amplitude_contrast, defocus=defocus,
                             sigma_decay_CTF=sigma_decay_CTF, camera_type=camera_type, camera_folder=camera_folder,
                             random_gradient=random_gradient, random_rotation=random_rotation)

    # Reconstruct tomogram
    if 'ReconstructTomogram' in config.sections():
        print('\n- Reconstructing tomogram')
        reconstruct_tomogram(start, end, sizeRecon, angles, output_folder, model_ID, weighting=weighting)

