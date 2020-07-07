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
from pytom.reconstruction.reconstructionStructures import *
from pytom.basic.files import *
import numpy as xp
import random
import constant_dictionaries as phys

class ConfigLogger(object):
    """
    Facilitates writing the conf file to a .log file in the outputFolder for reference of settings.
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


def generate_model(particleFolder, outputFolder, modelID, listpdbs, pixelSize = 1, size=1024, thickness=200,
                   solvent_potential=4.5301, sigma_structural=0.2, numberOfParticles=1000, placementSize=512, retries=5000):
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
    save_path = f'{outputFolder}/model_{modelID}'
    dims = [v.shape for v in volumes]
    particles_by_class = [0, ] * number_of_classes
    particle_nr = 1
    default_tries_left = retries
    skipped_particles = 0
    loc_x_start = loc_y_start = (size - placementSize) - (placementSize // 2)
    loc_x_end = loc_y_end = (size - placementSize) + (placementSize // 2)

    for _ in tqdm(range(numberOfParticles), desc='Placing particles'):

        # select random class
        cls_id = xp.random.randint(0, number_of_classes)

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
        ground_truth_txt_file += f'{listpdbs[cls_id]} {loc_x:.4f} {loc_y - 256:.4f} {loc_z - 256:.4f} ' \
                                 f'{p_angles[0]:.4f} {p_angles[1]:.4f} {p_angles[2]:.4f}\n'

    # add solvent background potential
    print('Adding background solvent potential')
    cell += solvent_potential

    # Add structural nois
    print('Adding structural noise to grand model cell')
    noisy_cell = cell + xp.random.normal(0, sigma_structural, cell.shape)

    # save grandmodels
    print('Saving grandmodels')
    pytom.tompy.io.write(f'{save_path}/grandmodel_noisefree.mrc', cell)
    pytom.tompy.io.write(f'{save_path}/grandmodel.mrc', noisy_cell)
    pytom.tompy.io.write(f'{save_path}/grandmodel_noisefree_cropped.mrc', cell[placementSize // 2:-placementSize // 2,
                                                                  placementSize // 2:-placementSize // 2, :])
    pytom.tompy.io.write(f'{save_path}/grandmodel_cropped.mrc', noisy_cell[placementSize // 2:-placementSize // 2,
                                                        placementSize // 2:-placementSize // 2, :])

    # save class masks
    print('Saving class volumes')
    pytom.tompy.io.write(f'{save_path}/class_bbox.mrc', class_bbox_mask[placementSize//2:-placementSize//2,
                                                                  placementSize//2:-placementSize//2, :])
    pytom.tompy.io.write(f'{save_path}/class_mask.mrc', class_accurate_mask[placementSize//2:-placementSize//2,
                                                                  placementSize//2:-placementSize//2, :])

    # save occupancy masks
    print('Saving occupancy volumes')
    pytom.tompy.io.write(f'{save_path}/occupancy_bbox.mrc', occupancy_bbox_mask[placementSize//2:-placementSize//2,
                                                                  placementSize//2:-placementSize//2, :])
    pytom.tompy.io.write(f'{save_path}/occupancy_mask.mrc', occupancy_accurate_mask[placementSize//2:-placementSize//2,
                                                                  placementSize//2:-placementSize//2, :])

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


def create_rotation_model(outputFolder, modelID):

    # Load grandmodel
    save_path = f'{outputFolder}/model_{modelID}'
    if os.path.exists(f'{save_path}/grandmodel.mrc'):
        grandcell = pytom.tompy.io.read_mrc(f'{save_path}/grandmodel.mrc')
    else:
        raise Exception(f'create_rotation_model expects grandmodel be created before ({save_path}/grandmodel.mrc)')

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


def calcCTF(vol, pix_size, Dz, voltage=300E3, Cs=2.7E-3, sigma_decay_ctf=0.4, amplitude_contrast=0.07):
    '''
    TODO UPDATE DESCRIPTION
    %TOM_CREATE_CTF calculates 2D or 3D CTF (pure phase contrast)
    %
    %   Note: only tested for even dimensions!
    %
    %   [ctf_out amplitude] = tom_create_ctf(Dz, vol, pix_size, voltage, Cs, sigma)
    %
    %PARAMETERS
    %  INPUT
    @param Dz: Defocus (>0 underfocus, <0 overfocus) (in m)
    @type Dz: float
    %   vol      : Volume (or image)
    %   pix_size : pixel size (in nm) (default: 0.72 nm)
    %   voltage  : accelerating Voltage (in eV) (default: 300 kV)
    %   Cs       : sperical aberration
    %   Sigma    : envelope of ctf (optional). If chosen decay
    %              ctf ~exp(-(freq/sigma)^2) is assumed. SIGMA is in units of Nyqvist.
    %              => ctf(sigma) = 1/e
    %
    %  OUTPUT
    %   ctf_out  : output containing the centrosymmetric ctf. It can be used as
    %              a fourier filter.
    %
    %EXAMPLE
    %   im = tom_emread('proteasome.em');
    %   ctf = tom_create_ctf(-4.4,im.Value,1.0240, 120,2);
    %   tom_imagesc(ctf)
    %
    %REFERENCES
    %
    %SEE ALSO
    %   TOM_CTF TOM_FOURIER TOM_IFOURIER
    %
    %
    %   created by FF 09/14/04
    %
    %   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
    %   Journal of Structural Biology, 149 (2005), 227-234.
    %
    %   Copyright (c) 2004-2007
    %   TOM toolbox for Electron Tomography
    %   Max-Planck-Institute of Biochemistry
    %   Dept. Molecular Structural Biology
    %   82152 Martinsried, Germany
    %   http://www.biochem.mpg.de/tom
    '''

    Lambda = wavelength_eV2m(voltage)

    Ny = 1 / (2 * pix_size)

    if len(vol.shape) > 2:
        R, Y, Z = xp.meshgrid(xp.arange(-Ny, Ny, 2 * Ny / vol.shape[0]), xp.arange(-Ny, Ny, 2 * Ny / vol.shape[1]),
                              xp.arange(-Ny, Ny, 2 * Ny / vol.shape[2]))
        r = xp.sqrt(R ** 2 + Y ** 2 + Z ** 2)
    else:
        R, Y = xp.meshgrid(xp.arange(-Ny, Ny, 2. * Ny / (vol.shape[0])), xp.arange(-Ny, Ny, 2. * Ny / (vol.shape[1])))
        r = xp.sqrt(R ** 2 + Y ** 2)

    # print(r)
    ctf = xp.exp( -1j * xp.pi / 2 * (Cs * (Lambda ** 3) * (r ** 4) - 2 * Dz * Lambda * (r ** 2)) )

    ctf.real = ctf.real * amplitude_contrast
    ctf.imag = ctf.imag * (1-amplitude_contrast)

    if sigma_decay_ctf:
        decay = xp.exp(-(r / (sigma_decay_ctf * Ny)) ** 2)
        ctf = ctf * decay

    return ctf


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


def microscope(noisefree_projections, outputFolder, modelID, dose=80, pixelsize=1E-9, voltage=300E3,
               camera_type='K2SUMMIT', camera_folder=''):
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

    if not os.path.exists(f'{outputFolder}/model_{modelID}/noisyProjections/'):
        os.mkdir(f'{outputFolder}/model_{modelID}/noisyProjections/')

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
    dose_per_pixel = dose_per_tilt * (pixelsize*1E9*10)**2  # from square A to square nm (10A pixels)
    print(f'number of electrons per pixel: {dose_per_pixel}')

    for n in range(n_images):
        projection = noisefree_projections[:, :, n]
        # Fourier transform and multiply with sqrt(dqe) = mtf/ntf
        projection_fourier = xp.fft.fftn(xp.fft.ifftshift(projection))
        projection_fourier = projection_fourier * mtf_shift / ntf_shift
        # Convert back to real space
        projection = xp.real(xp.fft.fftshift(xp.fft.ifftn(projection_fourier)))
        # Draw from poisson distribution and scale by camera's conversion factor
        conversion_factor = 100  # in ADU/e- , this is an arbitrary unit. Value taken from Vulovic et al., 2010
        projection_poisson = xp.random.poisson(lam=(projection * dose_per_pixel)) # + conversion_factor
        # Image values are now in ADU
        # Apply the camera's noise transfer function to the noisy image
        projection_fourier = xp.fft.fftn(xp.fft.ifftshift(projection_poisson)) * ntf_shift

        # Add additional noise from digitization process, less relevant for modern day cameras.
        # readout noise standard deviation can be 7 ADUs, from Vulovic et al., 2010
        # sigma_readout = 7
        # readsim = xp.random.normal(0, sigma_readout, projection.shape) # readout noise has a gaussian distribution
        # darksim = 0     # dark current noise has a poisson distribution, usually an order of magnitude smaller than readout
        #                 # noise and can hence be neglected

        # Add readout noise and dark noise in real space
        projection = xp.real(xp.fft.fftshift(xp.fft.ifftn(projection_fourier))) # + readsim + darksim

        projections[:,:,n] = projection

    return projections


def generate_projections(angles, outputFolder, modelID, dose=80, pixelSize=1E-9, voltage=300E3, sphericalAberration=2.7E-3,
                         multislice=False, msdz=5E-9, amplitudeContrast=0.07, defocus=2E-6, sigmaDecayCTF=0.4,
                         camera_type='K2SUMMIT', camera_folder=''):
    # NOTE; Parameter defocus specifies the defocus at the bottom of the model!

    grandcell = pytom.tompy.io.read_mrc(f'{outputFolder}/model_{modelID}/rotations/rotated_volume_{0}.mrc')

    imageSize = grandcell.shape[0]//2
    heightBox = grandcell.shape[2]

    zheight = heightBox * pixelSize  # thickness of volume in nm?
    defocus -= (zheight / 2) # center defocus value at tilt angle

    noisefree_projections = xp.zeros((imageSize, imageSize, len(angles)))

    #TODO allow for different defocus per tilt image. Now single defocus for all images
    ctf = calcCTF(xp.zeros((imageSize,imageSize)), pixelSize, defocus, voltage=voltage, Cs=sphericalAberration,
                  sigma_decay_ctf=sigmaDecayCTF, amplitude_contrast=amplitudeContrast)

    # Check if msdz is viable, else correct it
    if msdz % pixelSize != 0:
        # Round the msdz to the closest integer mutltiple of the pixel size
        round_up = xp.round(msdz % pixelSize / pixelSize)
        if round_up:
            msdz += (pixelSize - msdz % pixelSize)
        else:
            msdz -= (msdz % pixelSize)
        print(f'We adjusted msdz to {msdz*1E9:.3f} nm to make it an integer multiple of pixel size.')

    for n, angle in enumerate(angles):

        filename = f'{outputFolder}/model_{modelID}/rotations/rotated_volume_{int(angle)}.mrc'
        rotated_volume = pytom.tompy.io.read_mrc(filename)  #* 1E1

        if not multislice:
            print('Simulating projection (without ms) from tilt angle ', angle)
            projected_tilt_image = rotated_volume[imageSize//2:-imageSize//2,imageSize//2:-imageSize//2,:].sum(axis=2)#.get() # remove .get() if fully cupy
            projected_tilt_image = xp.abs(xp.fft.ifftn(xp.fft.fftn(projected_tilt_image)*ctf))**2
            #TODO is abs()**2 the way to go for exit wave field of simple projection?

        else:
            print('Simulating projection (with ms) from tilt angle ', angle)

            # Determine the number of slices
            if msdz > zheight:
                n_slices = 1
                msdz = zheight
                print('The multislice step can not be larger than volume thickness. Use only one slice.\n')
            elif msdz < pixelSize:
                n_slices = heightBox
                msdz = pixelSize
                print('The multislice steps can not be smaller than the pixel size. Use the pixel size instead for the step size.\n')
            else:
                n_slices = int(xp.ceil(xp.around(zheight/msdz,3)))

            print('Number of slices: ', n_slices)
            # Allocate space for multislice projection
            projected_potent_ms = xp.zeros((imageSize,imageSize, n_slices), dtype=complex)

            px_per_slice = int(msdz/pixelSize)
            num_px_last_slice = heightBox % px_per_slice  # consider last slice might have a different number of pixels

            # Project potential for each slice (phase grating)
            for ii in range(n_slices):
                projected_potent_ms[:, :, ii] = rotated_volume[imageSize//2:-imageSize//2,imageSize//2:-imageSize//2,
                                                ii*px_per_slice : (ii+1)*px_per_slice].mean(axis=2) #.get() # remove .get() if fully cupy

            # calculate the transmission function for each slice
            psi_t = transmission_function(projected_potent_ms, voltage, msdz)

            # calculate the fresnel propagator (identical for same dz)
            P = fresnel_propagator(imageSize, pixelSize, voltage, msdz)

            # Wave propagation with MULTISLICE method
            psi_multislice = xp.zeros((imageSize,imageSize), dtype=complex) + 1 # should be complex datatype
            for ii in range(n_slices-min(1,num_px_last_slice)):
                waveField = xp.fft.fftn( xp.fft.ifftshift(psi_multislice) * xp.fft.ifftshift(psi_t[:, :, ii]) )
                psi_multislice = xp.fft.fftshift( xp.fft.ifftn((waveField * xp.fft.ifftshift(P) )) )

            # Calculate propagation through last slice in case the last slice contains a different number of pixels
            if num_px_last_slice:
                msdz_end = num_px_last_slice * pixelSize
                psi_t[:, :, -1] = transmission_function(projected_potent_ms[:, :, -1], voltage, msdz_end)
                P_end = fresnel_propagator(imageSize, pixelSize, voltage, msdz_end)
                waveField = xp.fft.fftn( xp.fft.ifftshift(psi_multislice) * xp.fft.ifftshift(psi_t[:, :, -1]) )
                psi_multislice = xp.fft.fftshift( xp.fft.ifftn( waveField * xp.fft.ifftshift(P_end) ) )

            # Multiple by CTF for microscope effects on electron wave
            waveCTF = xp.fft.ifftshift(ctf) * xp.fft.fftn(xp.fft.ifftshift(psi_multislice) )
            # Intensity in image plane is obtained by taking the absolute square of the wave function
            projected_tilt_image = xp.abs(xp.fft.fftshift(xp.fft.ifftn(waveCTF))) ** 2

        noisefree_projections[:,:,n] = projected_tilt_image

    pytom.tompy.io.write(f'{outputFolder}/model_{modelID}/projections_noisefree.mrc',
                         noisefree_projections) # noisefree_projections.astype(xp.float32)?

    # Now add the microscope effects to the projections
    projections = microscope(noisefree_projections, outputFolder, modelID, dose=dose, pixelsize=pixelSize,
                             voltage=voltage, camera_type=camera_type, camera_folder=camera_folder)

    pytom.tompy.io.write(f'{outputFolder}/model_{modelID}/projections.mrc', projections)

    for n in range(projections.shape[2]):
        # These files are needed for the reconstruction
        pytom.tompy.io.write(f'{outputFolder}/model_{modelID}/noisyProjections/simulated_proj_{n+1}.mrc',
                             projections[:, :, n].squeeze())
    return


def reconstruct_tomogram(prefix, suffix, start_idx, end_idx, vol_size, angles, outputFolder, modelID, weighting=-1):
    from pytom.reconstruction.reconstructionStructures import Projection, ProjectionList
    projections = ProjectionList()

    for i in range(start_idx, end_idx+1):
        # IMPORTANT: angles *-1 to get the right reconstrunction relative to the orignal model!
        p = Projection(prefix+str(i)+suffix, tiltAngle= -1 * angles[i-1])
        projections.append(p)

    outputname = os.path.join(outputFolder, f'model_{modelID}/reconstruction.em')

    # IF EM alignment file provided, filters applied and reconstruction will be identical.
    vol = projections.reconstructVolume(dims=vol_size, reconstructionPosition=[0,0,0], binning=1, applyWeighting=weighting)
    vol.write(outputname)
    os.system(f'em2mrc.py -f {outputname} -t {os.path.dirname(outputname)}')
    os.system(f'rm {outputname}')
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
        outputFolder = config['General']['OutputFolder']
        modelID = int(config['General']['ModelID'])
        seed = int(config['General']['Seed'])
        pixelSize = float(config['General']['PixelSize']) * 1E-9
        metadata = loadstar(config['General']['MetaFile'], dtype=datatype)

        print(f'Generating model {modelID} in folder {outputFolder}')
    except Exception as e:
        print(e)
        raise Exception('Missing general parameters in config file.')

    if 'GenerateModel' in config.sections():
        try:
            # We assume the particle models are in the desired voxel spacing for the pixel size of the simulation!
            particleFolder = config['GenerateModel']['ParticleFolder']
            listpdbs = eval(config['GenerateModel']['Models'])
            placementSize = int(config['GenerateModel']['PlacementSize'])
            size = int(config['GenerateModel']['Size'])
            sigma_structural = float(config['GenerateModel']['SigmaStructuralNoise'])

            # Randomly vary solvent potential
            solvent_potential = float(config['GenerateModel']['SolventPotential'])
            xp.random.seed(seed)
            random.seed(seed)
            factor = xp.random.randint(-1, 4) * 0.1
            solvent_potential = solvent_potential + factor * solvent_potential

            # parse range of ice thickness
            thicc_range = config['GenerateModel']['Thickness'].split('-')
            if len(thicc_range) > 1:
                xp.random.seed(seed)
                random.seed(seed)

                thickness = xp.random.randint(int(thicc_range[0]), int(thicc_range[1])) # Thickness is provided in number of voxels
            else:
                thickness = int(thicc_range[0])
            if thickness % 2 == 1:
                thickness += 1 # make it even to prevent mismatch of reconstruction and model

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
            multislice = config['GenerateProjections'].getboolean('MultiSlice')
            if multislice:
                msdz = float(config['GenerateProjections']['MultiSliceSize']) * 1E-9 # multislice step size in nm
            sigmaDecayCTF = float(config['GenerateProjections']['SigmaDecayCTF'])
            camera_type = config['GenerateProjections']['Camera']
            camera_folder = config['GenerateProjections']['CameraFolder']

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
            sphericalAberration = metadata['SphericalAberration'][0] * 1E-3  # spherical aberration in mm
            amplitudeContrast = metadata['AmplitudeContrast'][0]  # fraction of amplitude contrast

            # defocus = 3E-6 # We want to get these from the metafile eventually
            defocus_range = config['GenerateProjections']['Defocus'].split('-') # defocus in um
            if len(defocus_range) > 1:
                xp.random.seed(seed)
                random.seed(seed)

                defocus = xp.random.uniform(float(defocus_range[0]),float(defocus_range[1])) * 1E-6
            else:
                defocus = float(defocus_range) * 1E-6

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
    if not os.path.exists(os.path.join(outputFolder, f'model_{modelID}')):
        os.mkdir(os.path.join(outputFolder, f'model_{modelID}'))

    if os.path.exists(f'{outputFolder}/model_{modelID}/simulator.log'):
        os.remove(f'{outputFolder}/model_{modelID}/simulator.log')

    logging.basicConfig(filename='{}/model_{}/simulator-{date:%Y-%m-%d_%H:%M:%S}.log'.format(outputFolder, modelID,
                                                                                             date=datetime.datetime.now()), level=logging.INFO)
    config_logger = ConfigLogger(logging)
    config_logger(config)

    logging.info('Values of parameters that we randomly vary per simulation (only for generate model and generate projections):')
    if 'GenerateModel' in config.sections():
        logging.info(f'ice thickness = {thickness*pixelSize*1E9:.2f}nm')
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
        generate_model(particleFolder, outputFolder, modelID, listpdbs, pixelSize=pixelSize, size=size,
                       thickness=thickness, placementSize=placementSize, solvent_potential=solvent_potential,
                       numberOfParticles=numberOfParticles, sigma_structural=sigma_structural)

    # Generated rotated grand model versions
    if 'Rotation' in config.sections():
        print(f'\n- Rotating model')

        # If needed, creating rotations directory, removing leftover initial rotation volume
        dir = f'{outputFolder}/model_{modelID}/rotations'
        if not (os.path.exists(dir)):
            os.mkdir(dir)
        filename = f'{dir}/rotated_volume_0.mrc'
        if os.path.exists(filename):
            print(f'Found leftover rotated volume 0, removing it (path: {filename}')
            os.remove(filename)

        #  Create initial rotation volume (0 degree rotation)
        volume = create_rotation_model(outputFolder, modelID)

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

        print('\n- Generating projections')
        generate_projections(angles, outputFolder, modelID, dose=dose, pixelSize=pixelSize,
                             voltage=voltage, sphericalAberration=sphericalAberration,multislice=multislice, msdz=msdz,
                             amplitudeContrast=amplitudeContrast, defocus=defocus, sigmaDecayCTF=sigmaDecayCTF,
                             camera_type=camera_type, camera_folder=camera_folder)

    # Reconstruct tomogram
    if 'ReconstructTomogram' in config.sections():
        print('\n- Reconstructing tomogram')
        prefix = os.path.join(outputFolder, f'model_{modelID}/noisyProjections/simulated_proj_')
        suffix = '.mrc'
        vol_size = [sizeRecon, sizeRecon, sizeRecon]
        reconstruct_tomogram(prefix, suffix, start, end, vol_size, angles, outputFolder, modelID, weighting=weighting)

