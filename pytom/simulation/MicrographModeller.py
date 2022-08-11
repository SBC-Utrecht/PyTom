"""
Marten Chaillet's cryoET simulator - updated version of the simulator used in SHREC2019 and SHREC2020
Contributions from Ilja Gubins and Gijs van der Schot.
Original simulator (used in SHREC2019) written by Gijs van der Schot, which was loosely based on the simulator in the
TOM toolbox for matlab.
"""

# IO related modules
# from pytom.gui.mrcOperations import *
import configparser
import tracemalloc
import logging
import os
import datetime
import sys
from tqdm import tqdm

# math
# from pytom.basic.files import *
import numpy as xp
import random
import pytom.simulation.physics as physics

# Plotting, use Qt5Agg to prevent conflict with tkinter in pylab on cluster
# import matplotlib
# matplotlib.use('Qt5Agg')
# import matplotlib.pylab as plt


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


def downscale_class_mask(volume, binning, order=0):
    """
    Downscale class mask for binned reconstructions.

    @param volume: the mask to scale down, 3d array of ints
    @type  volume: L{np.ndarray}
    @param binning: number of times to bin (demagnify)
    @type  binning: L{int}
    @param order: interpolation order for scipy
    @type  order: L{int}

    @return: binned volume, 3d array of ints
    @rtype: L{np.ndarray}

    @author: Marten Chaillet
    """
    from scipy.ndimage import zoom

    if binning == 1:
        return volume

    if volume.dtype != 'int':
        print('Dtype of class mask/occupancy mask was not int, forcing to int.')
        volume.astype(int)

    return zoom(volume, 1/binning, order=order)


def draw_range(range, datatype, name):
    """
    Input parsing from config file. This parses possible ranges of values and randomly samples in the range. In case
    the input is just a single value that will be used instead and no random selection will be done. This allows users
    to dynamically specify a range or single value depending on the needs.

    @param range: list of two values to select in between
    @type  range: L{list} - [L{float},] * 2
    @param datatype: desired type for the parameter, either int or float
    @type  datatype: L{str}
    @param name: name of the simulation parameter
    @type  name: L{str}

    @return: a single value that was selected from the range, or single parsed value
    @rtype:  L{int} or L{float}
    """
    if type(range) == list and len(range) == 2:
        xp.random.seed(seed)
        random.seed(seed)
        if datatype == int:
            return xp.random.randint(range[0], range[1])
        elif datatype == float:
            return xp.random.uniform(range[0], range[1])
    elif type(range) == list and len(range) == 1:
        if datatype == int:
            return int(range[0])
        elif datatype == float:
            return float(range[0])
    elif type(range) == float or type(range) == int:
        if datatype == int:
            return int(range)
        elif datatype == float:
            return float(range)
    else:
        print(f'invalid data range or input type for parameter {name}')
        sys.exit(0)


def motion_blur(model, spacing, sigma):

    from pytom.simulation.microscope import fourier_grids, ctf_grids

    blurring_filter = xp.flip(xp.fft.ifftshift(xp.exp(- 2 * xp.pi ** 2 * sigma ** 2 *
                                       ctf_grids(fourier_grids(model.shape, 1 / (2 * spacing), reduced=True))[1]),
                                               axes=(0,1)), axis=2)

    return xp.fft.irfftn(xp.fft.rfftn(model) * blurring_filter).real


def generate_model(particle_folder, save_path, listpdbs, listmembranes, pixel_size=1.0,
                   size=1024, thickness=200,
                   solvent_potential=physics.V_WATER, solvent_factor=1.0, number_of_particles=1000,
                   placement_size=512, retries=5000, number_of_markers=0,
                   absorption_contrast=False, voltage=300E3, number_of_membranes=0, sigma_motion_blur=.0,
                   particle_flipping=None):
    """
    Generate a grand model of a cryo-EM sample. Particles, membranes and gold markers will be randomly rotated before
    being randomly placed in the volume. The program attempts to place the specified numbers, but only takes a max of
    tries before moving on to the next particle. The model will be saved to the save_path. In case the absorption
    contrast flag is set to true, an imaginary part will also be generated and saved as a separate file. Corresponding
    ground truth information will be saved to the output folder as well. These are a locations list, class mask, and
    occupancy mask.

    The particle that are loaded from the particle folder can be macromolecules or membrane pieces. The files in this
    folder are assumed as .mrc format and should have been previously generated with pytom/simulation/potential.py.
    These files can also store real and imaginary parts of the potential. Both have a strict file name type.
    real: {pdb}_{pixel_size:.2f}A_solvent-{solvent_potential*solvent_factor:.3f}_real.mrc
    imaginary: {pdb}_{pixel_size:.2f}A_solvent-{solvent_potential*solvent_factor:.3f}_imag_{voltage*1E-3:.0f}V.mrc'
    potential.py will automatically generate the files according to this format.

    @param particle_folder: folder where particles are stored that are to be placed in the sample
    @type  particle_folder: L{str}
    @param save_path: simulation project folder
    @type  save_path: L{str}
    @param listpdbs: a list of pdb identifiers to be placed in sample, structures need to present in particle folder in .mrc format (generated by potential.py)
    @type  listpdbs: L{list} - [L{str},] * n
    @param listmembranes: list of membranes model to place in simulation, should also be present in particle folder
    @type  listmembranes: L{list} - [L{str},] * n
    @param pixel_size: size of voxels in model in A, needed for getting particles at correct size, and generating gold markers
    @type  pixel_size: L{float}
    @param size: x, y size of box and consequently projections
    @type  size: L{int}
    @param thickness: z height of box, corresponds to ice layer thickness
    @type  thickness: L{int}
    @param solvent_potential: background solvent electrostatic potential, this is also a parameter for loading electrostatic potential particle folder
    @type  solvent_potential: L{float}
    @param solvent_factor: factor to relatively increase the solvent background, needs to correspond to generated electrostatic potential in particle folder
    @type  solvent_factor: L{float}
    @param number_of_particles: number of particles that we will try to place in the sample
    @type  number_of_particles: L{int}
    @param placement_size: specify placement region within the full volume
    @type  placement_size: L{int}
    @param retries: number of times to try and randomly place each particle
    @type  retries: L{int}
    @param number_of_markers: number of gold markers that we will try to place
    @type  number_of_markers: L{int}
    @param absorption_contrast: whether to include an imaginary part of the sample for absorption contrast
    @type  absorption_contrast: L{float}
    @param voltage: voltage of electron beam in eV, absorption contrast depends on this, and is therefore a variable for the file name of the particle
    @type  voltage: L{float}
    @param number_of_membranes: number of membranes that we will try to place
    @type  number_of_membranes: L{int}
    @param beam_damage_snr: signal to noise ratio of beam damage, will be applied as a normal distribution over the sample
    @type  beam_damage_snr: L{float}

    @return: - (model is written to save_path)
    @rtype:  None

    @author: Gijs van der Schot, Ilja Gubins, Marten Chaillet
    """
    # IMPORTANT: We assume the particle models are in the desired voxel spacing for the pixel size of the simulation!
    from pytom.simulation.potential import create_gold_marker
    from pytom.voltools import transform
    from pytom.agnostic.io import read_mrc, write

    # outputs
    X, Y, Z = size, size, thickness
    cell_real = xp.zeros((X, Y, Z))
    if absorption_contrast: cell_imag = xp.zeros_like(cell_real)

    # occupancy_bbox_mask = xp.zeros_like(cell_real)
    occupancy_accurate_mask = xp.zeros_like(cell_real)
    # class_bbox_mask = xp.zeros_like(cell_real)
    class_accurate_mask = xp.zeros_like(cell_real)
    ground_truth_txt_file = ''

    # load pdb volumes and pad them
    volumes_real = []
    if absorption_contrast: volumes_imag = []
    for pdb in listpdbs:
        try:
            # file paths for particles
            vol_base_name   = f'{pdb}_{pixel_size:.2f}A_solvent-{solvent_potential*solvent_factor:.3f}V'
            filename_real       = os.path.join(particle_folder, f'{vol_base_name}_real.mrc')
            filename_imag       = os.path.join(particle_folder, f'{vol_base_name}_imag_{voltage*1E-3:.0f}V.mrc')
            # load the particle
            vol_real = read_mrc(filename_real)
            if absorption_contrast:
                vol_imag = read_mrc(filename_imag)
                # make sure real and imaginary part are the same size
                if vol_real.shape != vol_imag.shape:
                    print(f'real and imaginary interaction potential not the same shape for {pdb}, skipping model')
                    continue
                volumes_imag.append(vol_imag)
            volumes_real.append(vol_real)
        except Exception as ee:
            print(ee)
            raise Exception(f'Could not open pdb {pdb}, skipping the model')

    # attributes
    number_of_classes = len(listpdbs)
    names_of_classes = listpdbs
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

    # Add large cell structures, such as membranes first!
    if number_of_membranes:
        number_of_classes += 1
        names_of_classes.append('vesicles')
        particles_by_class += [0]

        # class id of cell structures is always the same
        cls_id = number_of_classes - 1

        for _ in tqdm(range(number_of_membranes), desc='Placing membranes and micelles'):

            # Give membranes a numbered names to randomly index them
            membrane_model = listmembranes[xp.random.randint(0, len(listmembranes))]

            # file names membranes
            vol_base_name   = f'{membrane_model}_{pixel_size:.2f}A_solvent-{solvent_potential*solvent_factor:.3f}V'
            filename_real       = os.path.join(particle_folder, f'{vol_base_name}_real.mrc')
            filename_imag       = os.path.join(particle_folder, f'{vol_base_name}_imag_{voltage*1E-3:.0f}V.mrc')
            # load the vesicle
            membrane_vol_real = read_mrc(filename_real)
            if absorption_contrast:
                membrane_vol_imag = read_mrc(filename_imag)
                if membrane_vol_real.shape != membrane_vol_imag.shape:
                    print(f'skipped mebrane model {membrane_model} due to real and imaginary shape not matching')
                    continue

            u = xp.random.uniform(0.0, 1.0, (2,))
            theta = xp.arccos(2 * u[0] - 1)
            phi = 2 * xp.pi * u[1]
            psi = xp.random.uniform(0.0, 2 * xp.pi)
            p_angles = xp.rad2deg([theta, phi, psi])

            # randomly mirror and rotate particle
            try:
                membrane_vol_real = transform(membrane_vol_real, rotation=p_angles,
                                        rotation_order='szxz', interpolation='linear', device='cpu')
                if absorption_contrast: membrane_vol_imag = transform(membrane_vol_imag, rotation=p_angles,
                                        rotation_order='szxz', interpolation='linear', device='cpu')
                # do something to remove edge artifacts of rotation? linear instead of filt_bspline
            except Exception as e:
                print(e)
                print('Something went wrong while rotating?')
                continue

            threshold = 0.001
            membrane_vol_real[membrane_vol_real < threshold] = 0
            if absorption_contrast: membrane_vol_imag[membrane_vol_imag < threshold] = 0

            # thresholded particle
            accurate_particle_occupancy = membrane_vol_real > 0

            # find random location for the particle
            # allow membranes to be placed half outside of the grandcell
            xx, yy, zz = membrane_vol_real.shape
            tries_left = default_tries_left
            # x_cut_left, x_cut_right = 0,0
            x_cut_left, x_cut_right, y_cut_left, y_cut_right, z_cut_low, z_cut_high = 0,0,0,0,0,0
            while tries_left > 0:
                # loc_x = xp.random.randint(loc_x_start + xx // 2 + 1, loc_x_end - xx // 2 - 1)
                loc_x = xp.random.randint(loc_x_start, loc_x_end)
                loc_y = xp.random.randint(loc_y_start, loc_y_end)
                loc_z = xp.random.randint(0, Z)

                tries_left -= 1

                x_cut_left = abs(loc_x - xx // 2) if (loc_x - xx // 2) < 0 else 0
                x_cut_right = abs(X - (loc_x + xx // 2 + xx % 2)) if (loc_x + xx // 2 + xx % 2) > Y else 0
                # adjust for bbox indexing for membrane sticking out of box
                # for x limit between x_start and x_end
                y_cut_left = abs(loc_y - yy//2) if (loc_y - yy//2) < 0 else 0
                y_cut_right = abs(Y - (loc_y + yy//2 + yy%2)) if (loc_y + yy//2 + yy%2) > Y else 0

                z_cut_low = abs(loc_z - zz//2) if (loc_z - zz//2) < 0 else 0
                z_cut_high = abs(Z - (loc_z + zz//2 + zz%2)) if (loc_z + zz//2 + zz%2) > Z else 0

                # crop the occupancy by removing overhang
                accurate_particle_occupancy_crop = accurate_particle_occupancy[x_cut_left:xx-x_cut_right,
                                                            y_cut_left:yy-y_cut_right, z_cut_low:zz-z_cut_high]

                # calculate coordinates of bbox for the newly rotated particle
                # bbox_x = [loc_x - xx // 2, loc_x + xx // 2 + xx % 2]
                bbox_x = [loc_x - xx // 2 + x_cut_left, loc_x + xx // 2 + xx % 2 - x_cut_right]
                bbox_y = [loc_y - yy // 2 + y_cut_left, loc_y + yy // 2 + yy % 2 - y_cut_right]
                bbox_z = [loc_z - zz // 2 + z_cut_low, loc_z + zz // 2 + zz % 2 - z_cut_high]

                # print(bbox_x, bbox_y, bbox_z)
                # print(xx,yy, zz)
                # print(y_cut_left, y_cut_right, z_cut_low, z_cut_high)

                # create masked occupancy mask
                masked_occupancy_mask = occupancy_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1],
                                        bbox_z[0]:bbox_z[1]]
                masked_occupancy_mask = masked_occupancy_mask * accurate_particle_occupancy_crop

                # if the location fits (masked occupancy pixel-wise mask is empty), break the loop and use this location
                if masked_occupancy_mask.sum() == 0:
                    break

            # however if still can't fit, ignore this particle (also adds variance in how many particles are
            # actually put)
            if tries_left < 1:
                skipped_particles += 1
                continue

            # populate occupancy volumes
            # occupancy_bbox_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] = particle_nr
            occupancy_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
                accurate_particle_occupancy_crop * particle_nr

            # populate class masks
            # class_bbox_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] = (cls_id + 1)
            class_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
                accurate_particle_occupancy_crop * (cls_id + 1)

            # populate density volume
            cell_real[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
                membrane_vol_real[x_cut_left:xx-x_cut_right, y_cut_left:yy-y_cut_right, z_cut_low:zz-z_cut_high]
            if absorption_contrast: cell_imag[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
                membrane_vol_imag[x_cut_left:xx-x_cut_right, y_cut_left:yy-y_cut_right, z_cut_low:zz-z_cut_high]

            # update stats
            particle_nr += 1
            particles_by_class[cls_id] += 1

            # update text
            ground_truth_txt_file += f'vesicle {int(loc_x - loc_x_start)} {int(loc_y - loc_y_start)} {int(loc_z)} ' \
                                     f'NaN NaN NaN\n'

    # GOLD MARKERS WILL ALSO BE COUNTED TOWARDS TOTAL PARTICLE NUMBER
    # There are also added as an additional class
    if number_of_markers:
        number_of_classes += 1
        names_of_classes.append('fiducials')
        particles_by_class += [0]

        # class id of gold markers is always the same
        cls_id = number_of_classes - 1

        for _ in tqdm(range(number_of_markers), desc='Placing gold markers'):

            # create the gold marker in other function
            if absorption_contrast:
                gold_real, gold_imag = create_gold_marker(pixel_size, solvent_potential, oversampling=2,
                                                            solvent_factor=solvent_factor,
                                                            imaginary=True, voltage=voltage)
            else:
                gold_real = create_gold_marker(pixel_size, solvent_potential, oversampling=2,
                                                 solvent_factor=solvent_factor)

            u = xp.random.uniform(0.0, 1.0, (2,))
            theta = xp.arccos(2 * u[0] - 1)
            phi = 2 * xp.pi * u[1]
            psi = xp.random.uniform(0.0, 2 * xp.pi)
            p_angles = xp.rad2deg([theta, phi, psi])

            # randomly mirror and rotate particle
            try:
                gold_real = transform(gold_real, rotation=p_angles,
                                        rotation_order='szxz', interpolation='linear', device='cpu')
                if absorption_contrast: gold_imag = transform(gold_imag, rotation=p_angles,
                                        rotation_order='szxz', interpolation='linear', device='cpu')
                # do something to remove edge artifacts of rotation? linear instead of filt_bspline
            except Exception as e:
                print(e)
                print('Something went wrong while rotating?')
                continue

            threshold = 0.001
            gold_real[gold_real < threshold] = 0
            if absorption_contrast: gold_imag[gold_imag < threshold] = 0

            # thresholded particle
            accurate_particle_occupancy = gold_real > 0

            # find random location for the particle
            xx, yy, zz = gold_real.shape
            tries_left = default_tries_left
            while tries_left > 0:
                loc_x = xp.random.randint(loc_x_start + xx // 2 + 1, loc_x_end - xx // 2 - 1)
                loc_y = xp.random.randint(loc_y_start + yy // 2 + 1, loc_y_end - yy // 2 - 1)
                loc_z = xp.random.randint(zz // 2 + 1, Z - zz // 2 - 1)

                tries_left -= 1

                # calculate coordinates of bbox for the newly rotated particle
                bbox_x = [loc_x - xx // 2, loc_x + xx // 2 + xx % 2]
                bbox_y = [loc_y - yy // 2, loc_y + yy // 2 + yy % 2]
                bbox_z = [loc_z - zz // 2, loc_z + zz // 2 + zz % 2]

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
            # occupancy_bbox_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] = particle_nr
            occupancy_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
                accurate_particle_occupancy * particle_nr

            # populate class masks
            # class_bbox_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] = (cls_id + 1)
            class_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
                accurate_particle_occupancy * (cls_id + 1)

            # populate density volume
            cell_real[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += gold_real
            if absorption_contrast: cell_imag[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += gold_imag

            # update stats
            particle_nr += 1
            particles_by_class[cls_id] += 1

            # update text
            ground_truth_txt_file += f'fiducial {int(loc_x - loc_x_start)} {int(loc_y - loc_y_start)} {int(loc_z)} ' \
                                     f'NaN NaN NaN\n'

    for _ in tqdm(range(number_of_particles), desc='Placing particles'):

        # select random class but correct for artifact class if adding gold particles

        cls_id = xp.random.randint(0, number_of_classes - bool(number_of_markers) - bool(number_of_membranes))

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
            particle_real = volumes_real[cls_id]
            if absorption_contrast: particle_imag = volumes_imag[cls_id]

            # Code for randomly flipping particles to add both handedness to the model, useful for generalizing
            # training datasets
            if particle_flipping == 'Random':
                if xp.random.randint(2): # Generate true/false randomly
                    # Mirror the particle to cover both left and right handedness of the proteins
                    # TODO if particle is flipped it should be written down in the particle locations txt file
                    ax = xp.random.randint(3)
                    particle_real = xp.flip(particle_real, axis=ax)
                    if absorption_contrast: particle_imag = xp.flip(particle_imag, axis=ax)
            elif particle_flipping == 'Yes':
                particle_real = xp.flip(particle_real, axis=0)  # the axis to flip does not matter for mirroring
                if absorption_contrast: particle_imag = xp.flip(particle_imag, axis=0)

            # Rotate particle to the specified orientation
            rotated_particle_real = transform(particle_real, rotation=p_angles,
                                         rotation_order='szxz', interpolation='filt_bspline', device='cpu')
            if absorption_contrast: rotated_particle_imag = transform(particle_imag, rotation=p_angles,
                                         rotation_order='szxz', interpolation='filt_bspline', device='cpu')
        except Exception as e:
            print(e)
            print('Something went wrong while rotating?')
            continue

        # remove particle rotation artifacts
        # threshold = min(volumes[cls_id][volumes[cls_id] > 0]) / 10
        # the approach above doesn't work well with PDB particles, there are often voxels with values of ^-11
        threshold = 0.01
        rotated_particle_real[rotated_particle_real < threshold] = 0
        if absorption_contrast: rotated_particle_imag[rotated_particle_imag < threshold]

        # thresholded rotated particle
        accurate_particle_occupancy = rotated_particle_real > 0

        # find random location for the particle
        xx, yy, zz = rotated_particle_real.shape
        tries_left = default_tries_left
        while tries_left > 0:
            loc_x = xp.random.randint(loc_x_start + xx // 2 + 1, loc_x_end - xx // 2 - 1)
            loc_y = xp.random.randint(loc_y_start + yy // 2 + 1, loc_y_end - yy // 2 - 1)
            loc_z = xp.random.randint(zz // 2 + 1, Z - zz // 2 - 1)

            tries_left -= 1

            # calculate coordinates of bbox for the newly rotated particle
            # bbox_x = [loc_x - dims[cls_id][0] // 2, loc_x + dims[cls_id][0] // 2 + dims[cls_id][0] % 2]
            # bbox_y = [loc_y - dims[cls_id][1] // 2, loc_y + dims[cls_id][1] // 2 + dims[cls_id][1] % 2]
            # bbox_z = [loc_z - dims[cls_id][2] // 2, loc_z + dims[cls_id][2] // 2 + dims[cls_id][2] % 2]
            bbox_x = [loc_x - xx // 2, loc_x + xx // 2 + xx % 2]
            bbox_y = [loc_y - yy // 2, loc_y + yy // 2 + yy % 2]
            bbox_z = [loc_z - zz // 2, loc_z + zz // 2 + zz % 2]

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
        # occupancy_bbox_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] = particle_nr
        occupancy_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
            accurate_particle_occupancy * particle_nr

        # populate class masks
        # class_bbox_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] = (cls_id + 1)
        class_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
            accurate_particle_occupancy * (cls_id + 1)

        # populate density volume
        cell_real[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += rotated_particle_real
        if absorption_contrast: cell_imag[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
                                                    rotated_particle_imag

        # update stats
        particle_nr += 1
        particles_by_class[cls_id] += 1

        # update text
        ground_truth_txt_file += f'{listpdbs[cls_id]} {int(loc_x - loc_x_start)} {int(loc_y - loc_y_start)} {int(loc_z)} ' \
                                 f'{p_angles[0]:.4f} {p_angles[1]:.4f} {p_angles[2]:.4f}\n'

    # add motion blur to a certain resolution (like a b_factor)
    cell_real = motion_blur(cell_real, pixel_size, sigma_motion_blur)
    cell_imag = motion_blur(cell_imag, pixel_size, sigma_motion_blur)

    # grandmodel names
    filename_gm_real    = os.path.join(save_path, 'grandmodel.mrc')
    filename_gm_imag    = os.path.join(save_path, 'grandmodel_imag.mrc')
    filename_cm         = os.path.join(save_path, 'class_mask.mrc')
    filename_om         = os.path.join(save_path, 'occupancy_mask.mrc')
    filename_loc        = os.path.join(save_path, 'particle_locations.txt')
    filename_con        = os.path.join(save_path, 'class_conversion_table.txt')

    # save grandmodels
    print('Saving grandmodels')
    write(filename_gm_real, cell_real)
    if absorption_contrast: write(filename_gm_imag, cell_imag)
    # save class masks
    print('Saving class volumes')
    # pytom.agnostic.io.write(f'{save_path}/class_bbox.mrc', class_bbox_mask)
    write(filename_cm, class_accurate_mask)
    # save occupancy masks
    print('Saving occupancy volumes')
    # pytom.agnostic.io.write(f'{save_path}/occupancy_bbox.mrc', occupancy_bbox_mask)
    write(filename_om, occupancy_accurate_mask)
    # save particle text file
    with open(filename_loc, 'w') as f:
        f.write(ground_truth_txt_file)
    # save conversion table from pdb to class
    with open(filename_con, 'w') as f:
        conversion_table = 'background 0\n'
        for i in range(number_of_classes):
            conversion_table += f'{names_of_classes[i]} {i+1}\n'
        f.write(conversion_table)

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
    #     viewer.add_image(cell_real, name='cell_real', interpolation='bicubic')

    # trying to reduce intermediate memory usage
    del cell_real, class_accurate_mask, occupancy_accurate_mask  # , class_bbox_mask, occupancy_bbox_mask
    if absorption_contrast: del cell_imag
    return


def create_ice_layer(shape, angle, width, value=1.0, sigma=0.0):
    """
    Efficiently create an ice layer at specified angle within the volume shape.
    Assumes the angle rotates perpendicular to the y-axis of the volume.

    @param shape: shape of volume to add the ice to, (x,y,z)
    @type  shape: L{tuple}, (L{int},) * 3
    @param angle: rotation angle of ice layer
    @type  angle: L{float}
    @param width: width of the ice layer in number of pixels
    @type  width: L{int}
    @param value: value to fill ice layer with
    @type  value: L{float}
    @param sigma: std of gaussian filter for smoothing edges in number of pixels (can be float)
    @type  sigma: L{float}

    @return: ice layer box, 3d array of floats
    @rtype:  L{np.ndarray}

    @author: Marten Chaillet
    """
    from scipy.ndimage import gaussian_filter
    assert xp.abs(angle) <= 90, print('rotation angle of ice layer cannot be larger than +- 90 degrees.')

    # get size
    xsize = shape[0]
    ysize = shape[1]
    zsize = shape[2]

    # create coordinates for x and z
    x = xp.arange(-xsize / 2, xsize / 2, 1, dtype=xp.float32)
    z = xp.arange(-zsize / 2, zsize / 2, 1, dtype=xp.float32)
    zm = xp.tile(z[xp.newaxis, :], [xsize, 1])

    # draw x values dependent on angle
    xline = x * xp.tan(-angle * xp.pi / 180)
    nwidth = width / xp.cos(angle * xp.pi / 180)
    xmin = xline - nwidth / 2
    xmax = xline + nwidth / 2

    # gradient for min and max value of x
    square_min = xp.tile(xmin[:, xp.newaxis], [1, zsize])
    square_max = xp.tile(xmax[:, xp.newaxis], [1, zsize])

    # Create smooth edge for side 1 of the layer
    grad1 = zm - square_min
    range1 = xp.abs(grad1.min() - grad1.max())
    c1 = (range1 / grad1.shape[0]) / 0.5
    grad1[grad1 > c1] = c1
    grad1[grad1 < -c1] = -c1
    grad1 = (grad1 - grad1.min()) / (grad1.max() - grad1.min())
    # Create smooth edge for side 2 of the layer
    grad2 = square_max - zm
    range2 = xp.abs(grad2.min() - grad2.max())
    c2 = (range2 / grad2.shape[0]) / 0.5
    grad2[grad2 > c2] = c2
    grad2[grad2 < -c2] = -c2
    # prevent division by zero if ice layer fills up the whole volume
    if (grad2.max() - grad2.min()) == 0:
        grad2 = grad2 / grad2.max()
    else:
        grad2 = (grad2 - grad2.min()) / (grad2.max() - grad2.min())

    # if a sigma is provided apply gaussian filter
    if not (sigma == 0.0):
        layer = gaussian_filter(grad1 * grad2 * value, sigma=sigma)
    else:
        layer = grad1 * grad2 * value

    # before returning tile the 2d layer to a volume
    return xp.tile(layer[:, xp.newaxis, :], [1, ysize, 1])


def microscope_single_projection(noisefree_projection, dqe, mtf, dose, pixel_size, oversampling=1):
    """
    Apply detector functions to the electron probability wave in the image plane and sample from Poisson distribution
    to obtain electron counts per pixel.

    Projections is first multiplied with mtf/ntf in fourier space. In real space it is then multiplied with dose. Then
    draw from Poisson distribution oversampling**2 times. Finally:
    coarse_projection = (1/oversampling**2) * sum_{i}^{N} F-1( F(count_i) * ntf )

    @param noisefree_projection: input projection giving probability of detecting electron in pixel, 2d array of floats
    @type  noisefree_projection: L{np.ndarray}
    @param dqe: Fourier space detector DQE, 2d array of floats
    @type  dqe: L{np.ndarray}
    @param mtf: Fourier space detector MTF, 2d array of floats
    @type  mtf: L{np.ndarray}
    @param dose: electron dose for this projection in e-/A^2
    @type  dose: L{float}
    @param pixel_size: pixel size in m
    @type  pixel_size: L{float}
    @param oversampling: oversampling times of pixel size was binned, influences the sampling from poisson distribution
    @type  oversampling: L{int}

    @return: projection after detection process, 2d array of floats
    @rtype:  L{np.ndarray}

    @author: Marten Chaillet
    """
    # from scipy.ndimage import shift

    ntf = xp.sqrt(mtf ** 2 / dqe)  # square root because dqe = mtf^2 / ntf^2
    ntf = xp.maximum(ntf, 1E-7)  # ensure we do not divide by zero
    mtf_shift = xp.fft.ifftshift(mtf)
    ntf_shift = xp.fft.ifftshift(ntf)

    # NUMBER OF ELECTRONS PER PIXEL
    dose_per_pixel = dose * (pixel_size*1E10)**2 / oversampling**2 # from square A to square nm (10A pixels)
    print(f'Number of electrons per pixel (before oversampling and sample absorption): {dose_per_pixel:.2f}')

    # Fourier transform and multiply with sqrt(dqe) = mtf/ntf
    projection_fourier = xp.fft.fftn(noisefree_projection) * mtf_shift / ntf_shift
    # projection_fourier = projection_fourier
    # Convert back to real space
    projection = xp.real(xp.fft.ifftn(projection_fourier))
    projection[projection<0] = 0

    # Apply poissonian noise
    projection_poisson = xp.zeros(projection.shape)
    for _ in range(oversampling**2):
        # Sample from the poisson distribution for oversampling^2 to account for the coarse graining of the simulation.
        # Normally in oversampling multiple pixels are averaged, obtaining a mean value of multiple poisson distributed
        # variables.
        poisson_intermediate = xp.random.poisson(lam=projection * dose_per_pixel)
        projection_poisson += xp.real(xp.fft.ifftn(xp.fft.fftn(poisson_intermediate) * ntf_shift)) / oversampling**2

    # Add additional noise from digitization process, less relevant for modern day cameras.
    # readout noise standard deviation can be 7 ADUs, from Vulovic et al., 2010
    # sigma_readout = 7
    # readsim = xp.random.normal(0, sigma_readout, projection.shape) # readout noise has a gaussian distribution
    # darksim = 0     # dark current noise has a poisson distribution, usually an order of magnitude smaller than
    #                 # readout noise and can hence be neglected

    # Add readout noise and dark noise in real space
    # return projection_poisson + readsim + darksim

    return projection_poisson


def parallel_project(grandcell, frame, image_size, pixel_size, msdz, n_slices, ctf, dose, dqe, mtf, voltage,
                     oversampling=1, translation=(.0,.0,.0), rotation=(.0,.0,.0), scale=(1., 1., 1.),
                     solvent_potential=physics.V_WATER, solvent_absorption=.0, ice_thickness_voxels=None, beam_damage_snr=0):
    """
    Project grandcell to create a frame/tilt. The grandcell will first be transformed according to the rotation and
    translation parameters, before projection. Microscope functions such as CTF, DQE, MTF need to be supplied. This
    function will create both a noisefree projection, and a projection with camera noise. The noisefree projection has
    been projected and convoluted with CTF, it is the electron wave at the image plane. The projection with noise is
    the noisefree image fed through the detection process, where DQE/MTF are applied and the probability of electron
    dose is sampled from poisonnian noise.

    @param grandcell: grandmodel that needs to be projected, 3d array of np.float32 or np.complex64
    @type  grandcell: L{np.ndarray}
    @param frame: index of the frame in the series
    @type  frame: L{int}
    @param image_size: size of the projection image
    @type  image_size: L{int}
    @param pixel_size: pixel size of image in m, 1 nm would be 1e-9
    @type  pixel_size: L{float}
    @param msdz: slice (step) size for multislice method in m, 5e-9 is often a good value
    @type  msdz: L{float}
    @param n_slices: total number of multislice steps along the z dimension of volume
    @type  n_slices: L{int}
    @param ctf: Fourier space CTF function, 3d array of complex values
    @type  ctf: {np.ndarray}
    @param dose: specific dose for this single tilt/frame, in e-/A^2
    @type  dose: L{float}
    @param dqe: DQE function in fourier space, 3d array of floats
    @type  dqe: L{np.ndarray}
    @param mtf: MTF function in fourier space, 3d array of floats
    @type  mtf: L{np.ndarray}
    @param voltage: voltage of electron beam in eV
    @type  voltage: L{float}
    @param oversampling: number of times pixel size was binned for dose correction
    @type  oversampling: L{int}
    @param translation: translation in x,y,z
    @type  translation: L{tuple} - (L{float},) * 3
    @param rotation: rotation along x,y,z angles, for tilt-series it suffices to specify y, i.e. (.0, tilt angle, .0)
    @type  rotation: L{tuple} - (L{float},) * 3
    @param solvent_potential: solvent potential used to fill in background ice
    @type  solvent_potential: L{float}
    @param solvent_absorption: absorption contrast to fill in background ice
    @type  solvent_absorption: L{float}
    @param ice_thickness_voxels: thickness of ice layers in voxels
    @type  ice_thickness_voxels: L{int}

    @return: (noisefree projection, projection)
    @rtype:  L{tuple} - (L{np.ndarray}, L{np.ndarray})

    @author: Marten Chaillet
    """
    from pytom.voltools import transform
    from pytom.simulation.microscope import transmission_function, fresnel_propagator

    print('Transforming sample for tilt/frame ', frame)

    sample = grandcell.copy()

    # think of a rotating cube in a larger box
    # sample box height = image_size * sin(tilt_angle) + ice_height * sin(90 - tilt_angle)
    max_tilt_radians = abs(rotation[1]) * xp.pi / 180
    max_tilt_radians_opp = (90 - abs(rotation[1])) * xp.pi / 180
    rotation_height = int(xp.ceil(xp.sin(max_tilt_radians) * image_size +
                                           xp.sin(max_tilt_radians_opp) * ice_thickness_voxels))
    print(f'Reduced rotation height for relevant specimens: {rotation_height}')
    if rotation_height % 2: rotation_height += 1
    diff = sample.shape[2]-rotation_height
    i = diff // 2

    # model parameter is a complex volume or real volume depending on the addition of absorption contrast
    # first transform the volume
    if sample.dtype == 'complex64':
        transform(sample[...,i:-i].real, translation=translation, rotation=rotation, rotation_order='sxyz',
                  scale=scale, interpolation='filt_bspline', device='cpu', output=sample[...,i:-i].real)
        transform(sample[...,i:-i].imag, translation=translation, rotation=rotation, rotation_order='sxyz', scale=scale,
                  interpolation='filt_bspline', device='cpu', output=sample[...,i:-i].imag)
        # remove rotation artifacts
        threshold = 0.001
        sample.real[sample.real < threshold] = 0
        sample.imag[sample.imag < threshold] = 0
    elif sample.dtype == 'float32':
        transform(sample[...,i:-i], translation=translation, rotation=rotation, rotation_order='sxyz', scale=scale,
                                interpolation='filt_bspline', device='cpu', output=sample[...,i:-i])
        # remove rotation artifacts
        threshold = 0.001
        sample[sample < threshold] = 0
    else:
        print('Invalid dtype of sample.')
        sys.exit(0)

    print('Simulating projection with multislice method for frame/tilt ', frame)
    # Start with projection, first prepare parameters
    box_size = sample.shape[0]
    box_height = sample.shape[2]

    # add the ice layer to the sample
    if sample.dtype == 'complex64' and ice_thickness_voxels is not None:
        sample.imag += create_ice_layer(sample.shape, rotation[1], ice_thickness_voxels, value=solvent_absorption,
                                        sigma=0.0)

    if n_slices==box_height and image_size==box_size:
        projected_potent_ms = sample
        px_per_slice = 1
        num_px_last_slice = box_height % px_per_slice
    else:
        ileft = (box_size - image_size) // 2
        iright = -int(xp.ceil((box_size - image_size) / 2))

        # Allocate space for multislice projection
        projected_potent_ms = xp.zeros((image_size, image_size, n_slices), dtype=xp.complex64)

        px_per_slice = int(msdz/pixel_size)
        num_px_last_slice = box_height % px_per_slice  # consider last slice might have a different number of pixels

        # Project potential for each slice (phase grating)
        for ii in range(n_slices):
            if ileft==0 and iright==0:
                projected_potent_ms[:, :, ii] = sample[:,:, ii * px_per_slice: (ii + 1) * px_per_slice].mean(axis=2) #.get()
            else:
                projected_potent_ms[:, :, ii] = sample[ileft:iright, ileft:iright,
                                                ii * px_per_slice: (ii + 1) * px_per_slice].mean(axis=2) #.get()
    # at this point we no longer need the sample, because all information in now contained in the projected slices
    # free the memory to accomodate space
    del sample

    # calculate the transmission function for each slice
    psi_t = transmission_function(projected_potent_ms, voltage, msdz)

    # calculate the fresnel propagator (identical for same dz)
    propagator = fresnel_propagator(image_size, pixel_size, voltage, msdz)

    # Wave propagation with MULTISLICE method, psi_multislice is complex
    psi_multislice = xp.zeros((image_size,image_size), dtype=xp.complex64) + 1 # +1 for initial probability

    # Loop over all the slices, except the last one if it has a different slice size
    for ii in range(n_slices-min(1,num_px_last_slice)):
        wave_field = xp.fft.fftn( xp.fft.ifftshift(psi_multislice) * xp.fft.ifftshift(psi_t[:, :, ii]) )
        psi_multislice = xp.fft.fftshift( xp.fft.ifftn((wave_field * xp.fft.ifftshift(propagator) )) )

    # Calculate propagation through last slice in case the last slice contains a different number of pixels
    if num_px_last_slice:
        msdz_end = num_px_last_slice * pixel_size
        psi_t[:, :, -1] = transmission_function(projected_potent_ms[:, :, -1], voltage, msdz_end)
        propagator_end = fresnel_propagator(image_size, pixel_size, voltage, msdz_end)
        wave_field = xp.fft.fftn( xp.fft.ifftshift(psi_multislice) * xp.fft.ifftshift(psi_t[:, :, -1]) )
        psi_multislice = xp.fft.fftshift( xp.fft.ifftn( wave_field * xp.fft.ifftshift(propagator_end) ) )

    # Multiple by CTF for microscope effects on electron wave
    wave_ctf = xp.fft.ifftshift(ctf) * xp.fft.fftn(xp.fft.ifftshift(psi_multislice) )
    # Intensity in image plane is obtained by taking the absolute square of the wave function
    noisefree_projection = xp.abs(xp.fft.fftshift(xp.fft.ifftn(wave_ctf))) ** 2

    # TODO Beam damage might be better modelled as a Gaussian drop off that increases in the order of projecting
    if beam_damage_snr > 0.0:
        sigma_signal = noisefree_projection[noisefree_projection>0.1].std()
        # snr = sigma_signal**2 / sigma_noise**2
        sigma_damage = xp.sqrt(sigma_signal**2 / beam_damage_snr)
        print(f'Standard deviation in beam damage SNR set to {sigma_damage} for SNR of {beam_damage_snr}')
        beam_noise = xp.random.normal(0, scale=sigma_damage, size=noisefree_projection.shape)
    else:
        beam_noise = 0

    # DEBUGGING
    # test = xp.log(xp.abs(xp.fft.fftshift(xp.fft.fftn(noisefree_projection))))
    # r1, m1 = radial_average(test)
    # fig, ax = plt.subplots()
    # ax.plot(r1, m1)
    # ax.legend()
    # plt.savefig(f'{folder}/radial.png')

    # Apply the microscope function
    projection = microscope_single_projection(noisefree_projection + beam_noise, dqe, mtf, dose, pixel_size,
                                              oversampling=oversampling)
    # Write the projection to the projection folder
    # pytom.agnostic.io.write(f'{folder}/synthetic_{frame+1}.mrc', projection)
    # Return noisefree and projection as tuple for writing as mrc stack in higher function
    return (noisefree_projection, projection)


def generate_tilt_series_cpu(save_path,
                             angles,
                             nodes=1,
                             image_size=None,
                             rotation_box_height=None,
                             pixel_size=1E-9,
                             oversampling=1,
                             dose=80,
                             voltage=300E3,
                             spherical_aberration=2.7E-3,
                             chromatic_aberration=2.7E-3,
                             energy_spread=0.7,
                             illumination_aperture=0.030E-3,
                             objective_diameter=100E-6,
                             focus_length=4.7E-3,
                             astigmatism=0.0E-9,
                             astigmatism_angle=0,
                             msdz=1E-9,
                             defocus=2E-6,
                             sigma_shift=0.,
                             sigma_angle_in_plane_rotation=0.,
                             sigma_magnification=0.,
                             sigma_tilt_angle=0.,
                             camera_type='K2SUMMIT',
                             camera_folder='',
                             solvent_potential=physics.V_WATER,
                             absorption_contrast=False,
                             beam_damage_snr=0,
                             grandcell=None):
    """
    Creating a tilt series for the initial grand model by rotating the sample along a set of tilt angles. For each angle
    the projection process of the microscope will be simulated. Calculation of each projection will be done on CPU
    nodes, as specified by nodes parameter. Computational cost can quickly increase for small pixel sizes with no
    oversampling factor. Sufficient RAM memory needs to be available to store each instance of the grandmodel for projection.

    @param save_path: simulation project folder
    @type  save_path: L{str}
    @param angles: tilt angles for the tilt series
    @type  angles: L{list} - [L{float},]
    @param nodes: number of CPU nodes to parallelize simulation over
    @type  nodes: L{int}
    @param image_size: number of pixels on x and y dimension
    @type  image_size: L{int}
    @param rotation_box_height: specify box height to rotate model in, if none given box height will be adjusted to encompass grandmodel at highest tilt angle
    @type  rotation_box_height: L{int}
    @param pixel_size: pixel size in m, default 1e-9 (i.e. 1 nm)
    @type  pixel_size: L{float}
    @param oversampling: number of times pixel size is binned for dose, mtf, and dqe correction
    @type  oversampling: L{int}
    @param dose: electron dose over the full series in e-/A^2, will be equally divided over number of frames, default 80
    @type  dose: L{float}
    @param voltage: voltage of electron beam in eV, default 300E3
    @type  voltage: L{float}
    @param spherical_aberration: spherical aberration of optical system in m, default 2.7e-3
    @type  spherical_aberration: L{float}
    @param chromatic_aberration: chromatic aberration of optical system in m, default 2.7e-3
    @type  chromatic_aberration: L{float}
    @param energy_spread: energy spread of electron beam as a fraction, default 0.7
    @type  energy_spread: L{float}
    @param illumination_aperture: aperture of electron beam in m, default 0.030e-3
    @type  illumination_aperture: L{float}
    @param objective_diameter: diameter of objective lens in m, default 100e-6
    @type  objective_diameter: L{float}
    @param focus_length: distance to focal point in m, default 4.7e-3
    @type  focus_length: L{float}
    @param astigmatism: strength of astigmatism in m, default 0.0e-9 (i.e. no astigmatism)
    @type  astigmatism: L{float}
    @param astigmatism_angle: direction of the astigmatism as an angle, default 0
    @type  astigmatism_angle: L{float}
    @param msdz: slice (step) size for the multislice method in m, default 5e-9
    @type  msdz: L{float}
    @param defocus: defocus of the sample in m, default 2e-6 (2 um), negative value would indicate overfocus
    @type  defocus: L{float}
    @param sigma_shift: standard deviation of random image shift in x and y, in A
    @type  sigma_shift: L{float}
    @param camera_type: type of detector device, e.g. K2SUMMIT
    @type  camera_type: L{str}
    @param camera_folder: folder where camera data is stored
    @type  camera_folder: L{str}
    @param solvent_potential: average background solvent electrostatic potential, default 4.5301
    @type  solvent_potential: L{float}
    todo solvent potential should be a boolean similar to absorption contrast??
    @param absorption_contrast: flag for including absorption contrast, grandmodel should have been generated with absorption contrast
    @type  absorption_contrast: L{bool}
    @param grandcell: optional parameter for passing grandcell directly to function, 3d array of floats or complex
    values
    @type  grandcell: L{np.ndarray}

    @return: - (projection are stored in save_path)
    @rtype:  None

    @author: Marten Chaillet
    """
    from pytom.basic.datatypes import DATATYPE_ALIGNMENT_RESULTS as dar
    from pytom.basic.datatypes import DATATYPE_METAFILE as dmf
    from pytom.basic.datatypes import fmtAlignmentResults, HEADER_ALIGNMENT_RESULTS, FMT_METAFILE, HEADER_METAFILE
    from pytom.gui.guiFunctions import savestar
    from pytom.simulation.microscope import create_detector_response, create_complex_ctf
    from pytom.simulation.microscope import convert_defocus_astigmatism_to_defocusU_defocusV
    from pytom.agnostic.io import read_mrc, write
    from joblib import Parallel, delayed
    # NOTE; Parameter defocus specifies the defocus at the bottom of the model!

    if grandcell is None:
        # grab model
        filename_gm_real = os.path.join(save_path, 'grandmodel.mrc')
        filename_gm_imag = os.path.join(save_path, 'grandmodel_imag.mrc')
        grandcell = read_mrc(filename_gm_real)
        if absorption_contrast:
            # set dtype to be complex64 to save memory
            xp_type = xp.complex64
            grandcell = grandcell.astype(xp_type)
            grandcell.imag = read_mrc(filename_gm_imag)
            # calculate the absorption for amorphous ice at the specified voltage
            solvent_amplitude = physics.potential_amplitude(physics.AMORPHOUS_ICE_DENSITY, physics.WATER_MW, voltage)
            print(f'solvent absorption = {solvent_amplitude:.3f}')
        else:
            # set dtype as float32 to save memory
            xp_type = xp.float32
            grandcell = grandcell.astype(xp_type)
            solvent_amplitude = 0.0
    else:
        if grandcell.dtype == complex:
            solvent_amplitude = physics.potential_amplitude(physics.AMORPHOUS_ICE_DENSITY, physics.WATER_MW, voltage)
            # set dtype to be complex64 to save memory
            xp_type = xp.complex64
            grandcell = grandcell.astype(xp_type)
        else:
            # set dtype as float32 to save memory
            xp_type = xp.float32
            grandcell = grandcell.astype(xp_type)
            solvent_amplitude = 0.0

    # extract variables of grandcell shape
    box_size = grandcell.shape[0]
    box_height = grandcell.shape[2]

    # set the image size, can be smaller than grandcell size
    if image_size is None:
        image_size = box_size
    else:
        # if provided, assert it is valid
        assert image_size <= box_size, 'Specified projection image size is invalid as it is larger than the model ' \
                                       'dimension.'

    assert box_size % 2 == 0 and image_size % 2 == 0, print(' - Simulations work only with even box sizes for easier '
                                                            'cropping and scaling.')

    # For sample set the arrays specifically to np.complex64 datatype to save memory space
    if rotation_box_height is not None:
        # Use the specified rotation box height
        rotation_volume = xp.zeros((box_size, box_size, rotation_box_height), dtype=xp_type)  # complex or real
        offset = (rotation_box_height - box_height) // 2

        if (rotation_box_height - box_height) % 2:
            rotation_volume[:, :, offset + 1:rotation_box_height-offset] = grandcell[:, :, :]
        else:
            rotation_volume[:, :, offset:rotation_box_height-offset] = grandcell[:, :, :]
    else:
        # Calculate maximal box height to contain all the information for the rotation
        # think of rotating parallel lines in a box
        # rotation_height = ice_height / cos(tilt_angle) + image_size * tan(tilt_angle)
        max_tilt = max([abs(a) for a in angles])
        max_tilt_radians = max_tilt * xp.pi / 180
        rotation_box_height = int(xp.ceil(xp.tan(max_tilt_radians) * image_size +
                                          box_height / xp.cos(max_tilt_radians)))
        if rotation_box_height % 2: rotation_box_height += 1
        # Place grandcell in rotation volume
        rotation_volume = xp.zeros((box_size, box_size, rotation_box_height), dtype=xp_type)
        offset = (rotation_box_height - box_height) // 2

        if (rotation_box_height - box_height) % 2:
            rotation_volume[:, :, offset + 1:rotation_box_height-offset] = grandcell[:, :, :]
        else:
            rotation_volume[:, :, offset:rotation_box_height-offset] = grandcell[:, :, :]
    del grandcell

    # adjust defocus because the value specifies defocus at the bottom of the box. the input expects defocus at the
    # center of the sample, therefore subtract half the box size.
    zheight = rotation_box_height * pixel_size  # thickness of rotation volume in nm
    defocus -= (zheight / 2)  # center defocus value at tilt angle

    # Check if msdz is viable, else correct it
    if msdz % pixel_size != 0:
        # Round the msdz to the closest integer mutltiple of the pixel size
        round_up = xp.round(msdz % pixel_size / pixel_size)
        if round_up:
            msdz += (pixel_size - msdz % pixel_size)
        else:
            msdz -= (msdz % pixel_size)
        print(f'We adjusted msdz to {msdz*1E9:.3f} nm to make it an integer multiple of pixel size.')

    # Determine the number of slices
    if msdz > zheight:
        n_slices = 1
        msdz = zheight
        print('The multislice step can not be larger than volume thickness. Use only one slice.\n')
    elif msdz < pixel_size:
        n_slices = box_height
        msdz = pixel_size
        print(
            'The multislice steps can not be smaller than the pixel size. Use the pixel size instead for the step size.\n')
    else:
        n_slices = int(xp.ceil(xp.around(zheight / msdz, 3)))
    print('Number of slices for multislice: ', n_slices)

    # determine dose per frame
    dose_per_tilt = dose / len(angles)

    # create motion blur by introducing random  x,y translation for each tilt
    translations = []
    magnifications = []
    in_plane_rotations = []
    for i in range(len(angles)):
        # generate random translation
        x = xp.random.normal(0, sigma_shift) / (pixel_size*1E10)
        y = xp.random.normal(0, sigma_shift) / (pixel_size*1E10)
        translations.append((x, y, 0.0))

        # generate random in plane rotation
        in_plane_rotations.append(xp.random.normal(0, sigma_angle_in_plane_rotation))

        # add magnifications
        if sigma_magnification != 0:
            gamma_a = 1. / (sigma_magnification ** 2)  # a = 1 / sigma**2
            gamma_b = 1. / gamma_a                     # b = mu / a
            # generate random magnification
            mag = xp.random.gamma(gamma_a, gamma_b)
            magnifications.append((mag, mag, 1.))
        else:
            magnifications.append((1., 1., 1.))

        # random dz for tilt angle
        angles[i] += xp.random.normal(0, sigma_tilt_angle)

    # defocus_series = [xp.random.normal(defocus, 0.2E-6) for a in angles]
    ctf_series = []
    dz_series, ast_series, ast_angle_series = [], [], []
    for x in angles:
        # todo currently input astigmatism is overriden by these options
        # todo add these options for frame series
        dz = xp.random.normal(defocus, 0.2e-6)
        ast = xp.random.normal(astigmatism, 0.1e-6)  # introduce astigmastism with 100 nm variation
        ast_angle = xp.random.normal(astigmatism_angle, 5)  # vary angle randomly around a 40 degree angle
        ctf = create_complex_ctf((image_size, image_size), pixel_size, dz, voltage=voltage,
                                     Cs=spherical_aberration, Cc=chromatic_aberration, energy_spread=energy_spread,
                                     illumination_aperture=illumination_aperture, objective_diameter=objective_diameter,
                                     focus_length=focus_length, astigmatism=ast,
                                     astigmatism_angle=ast_angle, display=False)
        ctf_series.append(ctf)
        dz_series.append(dz)
        ast_series.append(ast)
        ast_angle_series.append(ast_angle)

    dqe = create_detector_response(camera_type, 'DQE', image_size, voltage=voltage,
                                            folder=camera_folder, oversampling=oversampling)
    mtf = create_detector_response(camera_type, 'MTF', image_size, voltage=voltage,
                                            folder=camera_folder, oversampling=oversampling)

    # joblib automatically memory maps a numpy array to child processes
    print(f'Projecting the model with {nodes} processes')

    verbosity = 11  # set to 55 for debugging, 11 to see progress, 0 to turn off output
    results = Parallel(n_jobs=nodes, verbose=verbosity, prefer="threads") \
        (delayed(parallel_project)(rotation_volume, i, image_size, pixel_size, msdz, n_slices, ctf,
                                   dose_per_tilt, dqe, mtf, voltage, oversampling=oversampling, translation=translation,
                                   rotation=(.0, angle, in_plane_rotation), scale=magnification,
                                   solvent_potential=solvent_potential,
                                   solvent_absorption=solvent_amplitude, ice_thickness_voxels=box_height,
                                   beam_damage_snr=beam_damage_snr)
         for i, (angle, in_plane_rotation,
                 translation, magnification, ctf) in enumerate(zip(angles, in_plane_rotations,
                                                                   translations, magnifications, ctf_series)))

    sys.stdout.flush()

    if results.count(None) == 0:
        print('All projection processes finished successfully')
    else:
        print(f'{results.count(None)} rotation processes did not finish successfully')

    # write (noisefree) projections as mrc stacks
    filename_nf = os.path.join(save_path, 'noisefree_projections.mrc')
    filename_pr = os.path.join(save_path, 'projections.mrc')
    write(filename_nf, xp.stack([n for (n, p) in results], axis=2))
    write(filename_pr, xp.stack([p for (n, p) in results], axis=2))

    # todo create alignment and misalignment file, in reconstruction choice for alignment and misalignment
    # store alignment information
    # len(angles) is the number of files that we have
    alignment                       = xp.zeros(len(angles), dtype=dar)
    # IMPORTANT: get the inverse of each parameter for correct reconstruction
    alignment['TiltAngle']          = -1. * xp.array(angles)
    alignment['Magnification']      = 1. / xp.array([x for (x, y, z) in magnifications])
    alignment['AlignmentTransX']    = -1 * xp.array([x for (x, y, z) in translations])
    alignment['AlignmentTransY']    = -1 * xp.array([y for (x, y, z) in translations])
    alignment['InPlaneRotation']    = -1 * xp.array(in_plane_rotations)
    for i in range(len(angles)):
        alignment['FileName'][i]    = os.path.join(save_path, 'projections', f'synthetic_{i+1}.mrc')

    # write the alignment file
    filename_align                      = os.path.join(save_path, 'alignment_simulated.txt')
    savestar(filename_align, alignment, fmt=fmtAlignmentResults, header=HEADER_ALIGNMENT_RESULTS)

    # write meta file containing exactly all varied parameters in the simulation
    # get defocusU defocusV type defocus parameters
    defocusU, defocusV = convert_defocus_astigmatism_to_defocusU_defocusV(xp.array(dz_series), xp.array(ast_series))
    metafile                        = xp.zeros(len(angles), dtype=dmf)
    metafile['DefocusU']            = defocusU * 1e6
    metafile['DefocusV']            = defocusV * 1e6
    metafile['DefocusAngle']        = xp.array(ast_angle_series)
    metafile['Voltage']             = xp.array([voltage * 1e-3, ] * len(angles))
    metafile['SphericalAberration'] = xp.array([spherical_aberration * 1e3, ] * len(angles))
    metafile['PixelSpacing']        = xp.array([pixel_size * 1e10, ] * len(angles))
    metafile['TiltAngle']           = xp.array(angles)
    metafile['InPlaneRotation']     = xp.array(in_plane_rotations)
    metafile['TranslationX']        = xp.array([x for (x, y, z) in translations])
    metafile['TranslationY']        = xp.array([y for (x, y, z) in translations])
    metafile['Magnification']       = xp.array([x for (x, y, z) in magnifications])
    for i in range(len(angles)):
        alignment['FileName'][i]    = os.path.join(save_path, 'projections', f'synthetic_{i+1}.mrc')

    savestar(os.path.join(save_path, 'simulation.meta'), metafile, fmt=FMT_METAFILE, header=HEADER_METAFILE)

    # translations, dz_series, ast_series, ast_angle_series
    # header_varied_params = ''
    # dt_varied_params = [('TranslationX', 'f4'),
    #                               ('TranslationY', 'f4'),
    #                               ('Defocus', 'f4'),
    #                               ('Astigmatism', 'f4'),
    #                               ('AstigmatismAngle', 'f4')]
    # units_varied_params = ['px', 'px', 'um', 'um', 'degrees']
    #
    # for n, h in enumerate(dt_varied_params):
    #     header_varied_params += '{} {}\n'.format(h[0], '({})'.format(units_varied_params[n]) * (
    #             units_varied_params[n] != ''))
    #
    # fmt_varied_params = '%15.10f %15.10f %15.10f %15.10f %15.10f'
    #
    # varied_params = xp.zeros(len(angles), dtype=dt_varied_params)
    # varied_params['TranslationX'] = xp.array([x for (x, y, z) in translations])
    # varied_params['TranslationY'] = xp.array([y for (x, y, z) in translations])
    # varied_params['Defocus'] = (xp.array(dz_series) + zheight/2) * 1e6
    # varied_params['Astigmatism'] = xp.array(ast_series) * 1e6
    # varied_params['AstigmatismAngle'] = xp.array(ast_angle_series)
    #
    # savestar(os.path.join(save_path, 'varied_parameters.txt'), varied_params, fmt=fmt_varied_params,
    #          header=header_varied_params)

    return


def generate_frame_series_cpu(save_path, n_frames=20, nodes=1, image_size=None, pixel_size=1E-9,
                              oversampling=1, dose=80, voltage=300E3, spherical_aberration=2.7E-3,
                              chromatic_aberration=2.7E-3, energy_spread=0.7, illumination_aperture=0.030E-3,
                              objective_diameter=100E-6, focus_length=4.7E-3, astigmatism=0.0E-9, astigmatism_angle=0,
                              msdz=5E-9, defocus=2E-6, mean_shift=0.0, camera_type='K2SUMMIT', camera_folder='',
                              solvent_potential=physics.V_WATER, absorption_contrast=False, beam_damage_snr=0,
                              sigma_angle_in_plane_rotation=0.,
                              sigma_magnification=0.,
                              grandcell=None):
    """
    Creating a frame series for the initial grand model by applying stage drift (translation) for each frame, and
    calculating the sample projection in the microscope. Calculation of each projection will be done on CPU nodes, as
    specified by nodes parameter.

    todo add an option for incremental degradation of frame-series to simulate increased exposure

    @param save_path: simulation project folder
    @type  save_path: L{str}
    @param n_frames: number of frames to simulate
    @type  n_frames: L{int}
    @param nodes: number of CPU nodes to parallelize simulation over
    @type  nodes: L{int}
    @param image_size: number of pixels on x and y dimension
    @type  image_size: L{int}
    @param pixel_size: pixel size in m, default 1e-9 (i.e. 1 nm)
    @type  pixel_size: L{float}
    @param oversampling: number of times pixel size is binned for dose, mtf, and dqe correction
    @type  oversampling: L{int}
    @param dose: electron dose over the full series in e-/A^2, will be equally divided over number of frames, default 80
    @type  dose: L{float}
    @param voltage: voltage of electron beam in eV, default 300E3
    @type  voltage: L{float}
    @param spherical_aberration: spherical aberration of optical system in m, default 2.7e-3
    @type  spherical_aberration: L{float}
    @param chromatic_aberration: chromatic aberration of optical system in m, default 2.7e-3
    @type  chromatic_aberration: L{float}
    @param energy_spread: energy spread of electron beam as a fraction, default 0.7
    @type  energy_spread: L{float}
    @param illumination_aperture: aperture of electron beam in m, default 0.030e-3
    @type  illumination_aperture: L{float}
    @param objective_diameter: diameter of objective lens in m, default 100e-6
    @type  objective_diameter: L{float}
    @param focus_length: distance to focal point in m, default 4.7e-3
    @type  focus_length: L{float}
    @param astigmatism: strength of astigmatism in m, default 0.0e-9 (i.e. no astigmatism)
    @type  astigmatism: L{float}
    @param astigmatism_angle: direction of the astigmatism as an angle, default 0
    @type  astigmatism_angle: L{float}
    @param msdz: slice (step) size for the multislice method in m, default 5e-9
    @type  msdz: L{float}
    @param defocus: defocus of the sample in m, default 2e-6 (2 um), negative value would indicate overfocus
    @type  defocus: L{float}
    @param mean_shift: mean global shift of image in A, this will be spread over frames and in a random direction, default 10
    @type  mean_shift: L{float}
    @param camera_type: type of detector device, e.g. K2SUMMIT
    @type  camera_type: L{str}
    @param camera_folder: folder where camera data is stored
    @type  camera_folder: L{str}
    @param solvent_potential: average background solvent electrostatic potential, default 4.5301
    @type  solvent_potential: L{float}
    todo solvent potential should be a boolean similar to absorption contrast??
    @param absorption_contrast: flag for including absorption contrast, grandmodel should have been generated with absorption contrast
    @type  absorption_contrast: L{bool}
    @param grandcell: optional parameter for passing grandcell directly to function, 3d array of floats or complex
    values
    @type  grandcell: L{np.ndarray}

    @return: - (projection are stored in save_path)
    @rtype:  None

    @author: Marten Chaillet
    """
    from pytom.basic.datatypes import DATATYPE_ALIGNMENT_RESULTS as dar
    from pytom.basic.datatypes import DATATYPE_METAFILE as dmf
    from pytom.basic.datatypes import fmtAlignmentResults, HEADER_ALIGNMENT_RESULTS, FMT_METAFILE, HEADER_METAFILE
    from pytom.gui.guiFunctions import savestar
    from pytom.simulation.microscope import create_detector_response, create_complex_ctf, \
        convert_defocus_astigmatism_to_defocusU_defocusV
    from pytom.agnostic.io import read_mrc, write
    from joblib import Parallel, delayed

    # print('this function is not working because it does not create a metafile and does not correctly pass '
    #       'magnification to parallel project')
    # sys.exit(0)

    # NOTE; Parameter defocus specifies the defocus at the bottom of the model!
    if grandcell is None:
        # grab model
        filename_gm_real = os.path.join(save_path, 'grandmodel.mrc')
        filename_gm_imag = os.path.join(save_path, 'grandmodel_imag.mrc')
        grandcell = read_mrc(filename_gm_real)
        if absorption_contrast:
            # set dtype to be complex64 to save memory
            grandcell = grandcell.astype(xp.complex64)
            grandcell.imag = read_mrc(filename_gm_imag)
            # calculate the absorption for amorphous ice at the specified voltage
            solvent_amplitude = physics.potential_amplitude(physics.AMORPHOUS_ICE_DENSITY, physics.WATER_MW, voltage)
            print(f'solvent absorption = {solvent_amplitude:.3f}')
        else:
            # set dtype as float32 to save memory
            grandcell = grandcell.astype(xp.float32)
            solvent_amplitude = 0.0
    else:
        if grandcell.dtype == complex:
            solvent_amplitude = physics.potential_amplitude(physics.AMORPHOUS_ICE_DENSITY, physics.WATER_MW, voltage)
            # set dtype to be complex64 to save memory
            grandcell = grandcell.astype(xp.complex64)
        else:
            # set dtype as float32 to save memory
            grandcell = grandcell.astype(xp.float32)

    # extract size
    box_size = grandcell.shape[0]
    box_height = grandcell.shape[2]

    if image_size is None:
        image_size = box_size

    # confirm image_size is valid
    assert image_size <= box_size, 'Specified projection image size is invalid as it is larger than the model dimension.'

    # adjust defocus
    zheight = box_height * pixel_size  # thickness of rotation volume in nm
    defocus -= (zheight / 2)  # center defocus value at tilt angle

    # Check if msdz is viable, else correct it
    if msdz % pixel_size != 0:
        # Round the msdz to the closest integer mutltiple of the pixel size
        round_up = xp.round(msdz % pixel_size / pixel_size)
        if round_up:
            msdz += (pixel_size - msdz % pixel_size)
        else:
            msdz -= (msdz % pixel_size)
        print(f'We adjusted msdz to {msdz*1E9:.3f} nm to make it an integer multiple of pixel size.')

    # Determine the number of slices
    if msdz > zheight:
        n_slices = 1
        msdz = zheight
        print('The multislice step can not be larger than volume thickness. Use only one slice.\n')
    elif msdz < pixel_size:
        n_slices = box_height
        msdz = pixel_size
        print(
            'The multislice steps can not be smaller than the pixel size. Use the pixel size instead for the step size.\n')
    else:
        n_slices = int(xp.ceil(xp.around(zheight / msdz, 3)))
    print('Number of slices for multislice: ', n_slices)

    # determine dose per frame
    dose_per_frame = dose / n_frames

    # First generate stage drift and in-plane rotation, stage drift is a set of correlated translations across the
    # number of frames. In MotionCorr2 paper accumulated motion for 20S proteasome dataset is 11A across the whole
    # frame series.
    # First generate global motion and global direction of motion.
    global_motion = xp.random.normal(mean_shift, 3)  # normal around mean 10 A and std 3A
    average_motion_per_frame = global_motion / n_frames
    global_angle = xp.random.uniform(0,360)  # random angle from uniform
    translations, cumulative_translations, translations_voxel = [], [], []
    magnifications = []
    in_plane_rotations = []
    x, y = 0, 0
    for i in range(n_frames):
        # randomly vary the motion per frame and angle
        motion_i = xp.random.normal(average_motion_per_frame, average_motion_per_frame/2)
        angle_i = xp.random.normal(global_angle, 20)
        # decompose motion into x and y translation
        y_i = xp.sin(angle_i * xp.pi / 180) * motion_i
        x_i = xp.cos(angle_i * xp.pi / 180) * motion_i
        # only seem to need cumulative_translations
        translations.append((x_i, y_i, 0)) # append the translation for the frame as a tuple
        y += y_i
        x += x_i
        cumulative_translations.append((x,y, 0)) # translation for z coordinate as we are referring to volumes
        translations_voxel.append((x*1E-10 / pixel_size, y*1E-10 / pixel_size, 0))

        # generate random in plane rotation
        in_plane_rotations.append(xp.random.normal(0, sigma_angle_in_plane_rotation))

        # add magnifications
        if sigma_magnification != 0:
            gamma_a = 1. / (sigma_magnification ** 2)  # a = 1 / sigma**2
            gamma_b = 1. / gamma_a  # b = mu / a
            # generate random magnification
            mag = xp.random.gamma(gamma_a, gamma_b)
            magnifications.append((mag, mag, 1.))
        else:
            magnifications.append((1., 1., 1.))

    # write motion trajectory to a png file for debugging
    # fig, ax = plt.subplots(2)
    # ax[0].plot([x for (x,y,z) in cumulative_translations], [y for (x,y,z) in cumulative_translations], label='trajectory')
    # ax[0].set_xlabel('x (A)')
    # ax[0].set_ylabel('y (A)')
    # ax[0].legend()
    # ax[1].plot([x for (x,y,z) in translations_voxel], [y for (x,y,z) in translations_voxel], label='trajectory')
    # ax[1].set_xlabel('x (voxels)')
    # ax[1].set_ylabel('y (voxels)')
    # ax[1].legend()
    # plt.savefig(f'{save_path}/global_motion.png')
    # plt.close()

    # defocus_series = [xp.random.normal(defocus, 0.2E-6) for a in angles]
    ctf_series = []
    dz_series, ast_series, ast_angle_series = [], [], []
    for x in range(n_frames):
        # todo currently input astigmatism is overriden by these options
        # todo add these options for frame series
        dz = xp.random.normal(defocus, 0.2e-6)
        ast = xp.random.normal(astigmatism, 0.1e-6)  # introduce astigmastism with 100 nm variation
        ast_angle = xp.random.normal(astigmatism_angle, 5)  # vary angle randomly around a 40 degree angle
        ctf = create_complex_ctf((image_size, image_size), pixel_size, dz, voltage=voltage,
                                 Cs=spherical_aberration, Cc=chromatic_aberration, energy_spread=energy_spread,
                                 illumination_aperture=illumination_aperture, objective_diameter=objective_diameter,
                                 focus_length=focus_length, astigmatism=ast,
                                 astigmatism_angle=ast_angle, display=False)
        ctf_series.append(ctf)
        dz_series.append(dz)
        ast_series.append(ast)
        ast_angle_series.append(ast_angle)

    dqe = create_detector_response(camera_type, 'DQE', image_size, voltage=voltage,
                                            folder=camera_folder, oversampling=oversampling)
    mtf = create_detector_response(camera_type, 'MTF', image_size, voltage=voltage,
                                            folder=camera_folder, oversampling=oversampling)

    # joblib automatically memory maps a numpy array to child processes
    print(f'Projecting the model with {nodes} processes')

    verbosity = 55  # set to 55 for debugging, 11 to see progress, 0 to turn off output
    results = Parallel(n_jobs=nodes, verbose=verbosity, prefer="threads") \
        (delayed(parallel_project)(grandcell, frame, image_size, pixel_size, msdz, n_slices, ctf,
                                   dose_per_frame, dqe, mtf, voltage, oversampling=oversampling, translation=shift,
                                   rotation=(.0, .0, in_plane_rotation),
                                   scale=magnification, solvent_potential=solvent_potential,
                                   solvent_absorption=solvent_amplitude, ice_thickness_voxels=box_height,
                                   beam_damage_snr=beam_damage_snr)
        for frame, (in_plane_rotation,
                shift, magnification, ctf) in enumerate(zip(in_plane_rotations,
                                                                  translations_voxel, magnifications, ctf_series)))

    sys.stdout.flush()
    if results.count(None) == 0:
        print('All projection processes finished successfully')
    else:
        print(f'{results.count(None)} rotation processes did not finish successfully')

    # write (noisefree) projections as mrc stacks
    filename_nf = os.path.join(save_path, 'noisefree_projections.mrc')
    filename_pr = os.path.join(save_path, 'projections.mrc')
    write(filename_nf, xp.stack([n for (n,p) in results], axis=2))
    write(filename_pr, xp.stack([p for (n, p) in results], axis=2))

    # todo create alignment and misalignment file, in reconstruction choice for alignment and misalignment
    # store alignment information
    # len(angles) is the number of files that we have
    alignment = xp.zeros(n_frames, dtype=dar)
    # IMPORTANT: get the inverse of each parameter for correct reconstruction
    alignment['TiltAngle'] = xp.array([.0] * n_frames)
    alignment['Magnification'] = 1. / xp.array([x for (x, y, z) in magnifications])
    alignment['AlignmentTransX'] = -1 * xp.array([x for (x, y, z) in translations_voxel])
    alignment['AlignmentTransY'] = -1 * xp.array([y for (x, y, z) in translations_voxel])
    alignment['InPlaneRotation'] = -1 * xp.array(in_plane_rotations)
    for i in range(n_frames):
        alignment['FileName'][i] = os.path.join(save_path, 'projections', f'synthetic_{i+1}.mrc')

    # write the alignment file
    filename_align = os.path.join(save_path, 'alignment_simulated.txt')
    savestar(filename_align, alignment, fmt=fmtAlignmentResults, header=HEADER_ALIGNMENT_RESULTS)

    # write meta file containing exactly all varied parameters in the simulation
    # get defocusU defocusV type defocus parameters
    defocusU, defocusV = convert_defocus_astigmatism_to_defocusU_defocusV(xp.array(dz_series), xp.array(ast_series))
    metafile = xp.zeros(n_frames, dtype=dmf)
    metafile['DefocusU'] = defocusU * 1e6
    metafile['DefocusV'] = defocusV * 1e6
    metafile['DefocusAngle'] = xp.array(ast_angle_series)
    metafile['Voltage'] = xp.array([voltage * 1e-3, ] * n_frames)
    metafile['SphericalAberration'] = xp.array([spherical_aberration * 1e3, ] * n_frames)
    metafile['PixelSpacing'] = xp.array([pixel_size * 1e10, ] * n_frames)
    metafile['TiltAngle'] = xp.array([.0] * n_frames)
    metafile['InPlaneRotation'] = xp.array(in_plane_rotations)
    metafile['TranslationX'] = xp.array([x for (x, y, z) in translations_voxel])
    metafile['TranslationY'] = xp.array([y for (x, y, z) in translations_voxel])
    metafile['Magnification'] = xp.array([x for (x, y, z) in magnifications])
    for i in range(n_frames):
        alignment['FileName'][i] = os.path.join(save_path, 'projections', f'synthetic_{i+1}.mrc')

    savestar(os.path.join(save_path, 'simulation.meta'), metafile, fmt=FMT_METAFILE, header=HEADER_METAFILE)

    # # Store translations as reference for model
    # # len(angles) is the number of files that we have
    # alignment                       = xp.zeros(n_frames, dtype=dar)
    # alignment['AlignmentTransX']    = xp.array([x for (x,y,z) in cumulative_translations])
    # alignment['AlignmentTransY']    = xp.array([y for (x,y,z) in cumulative_translations])
    # alignment['Magnification']      = xp.repeat(1.0, n_frames)
    # for i in range(n_frames):
    #     alignment['FileName'][i]    = os.path.join(save_path, 'projections', f'synthetic_{i+1}.mrc')
    #
    # # Write the alignment file as a text file
    # filename_align                      = os.path.join(save_path, 'alignment_simulated.txt')
    # savestar(filename_align, alignment, fmt=fmtAlignmentResults, header=HEADER_ALIGNMENT_RESULTS)

    return


def FSS(fimage1, fimage2, numberBands, verbose=False):
    """
    algorithm FSS = Fourier Shell Scaling
    Scale the values of fimage1 to the values in fimage2 per band in fourier space. The mean of each band is calculated
    for both images. The values of the bands in image1 are divided by dividing by mean1 and multiplying by mean2.
    M{band = band * (m2 / m1)}

    @param fimage1: the fourier amplitudes to be scaled, 2d array
    @type  fimage1: L{np.ndarray}
    @param fimage2: the example fourier amplitudes, 2d array
    @type  fimage2: L{np.ndarray}
    @param numberBands: number of rings to scale, determines precision
    @type  numberBands: L{int}
    @param verbose: be talkative
    @type  verbose: L{bool}

    @return: the scaled fourier image, 2D array
    @type:   L{np.ndarray}

    @author: Marten Chaillet
    """
    from pytom.agnostic.correlation import meanUnderMask
    from pytom.simulation.support import bandpass_mask

    assert fimage1.shape == fimage2.shape, "volumes not of same size"
    assert len(set(fimage1.shape)) == 1, "volumes are not perfect cubes"

    if verbose:
        print(f'shape of images is: {fimage1.shape}')

    increment = int(fimage1.shape[0] / 2 * 1 / numberBands)
    band = [-1, -1]

    output = xp.zeros(fimage1.shape)

    for i in range(0, fimage1.shape[0] // 2, increment):

        band[0] = i
        band[1] = i + increment

        if verbose:
            print('Band : ', band)

        if band[1] >= fimage1.shape[0] // 2:
            bandpass = bandpass_mask(fimage1.shape, 0, high=band[0])
            bandpass = (bandpass == 0) * 1
        else:
            bandpass = bandpass_mask(fimage1.shape, band[0], band[1])

        if i == 0:
            # remove center point from mask
            c = bandpass.shape[0] // 2
            bandpass[c, c] = 0
            # scale center separately as center adjustment has large influence on the image
            output[c, c] = fimage1[c, c] / (fimage1[c, c] / fimage2[c, c])

        n = bandpass.sum()
        # get mean amplitude of each band
        m1 = meanUnderMask(fimage1, bandpass, p=n)
        m2 = meanUnderMask(fimage2, bandpass, p=n)

        # scale the values inside the band of volume1
        outband = fimage1 * bandpass / (m1 / m2)

        output += outband

    del outband, bandpass
    return output


def scale_image(image1, image2, numberBands):
    """
    Scale the amplitudes of image1 to those of image2 in fourier space. This function handles Fourier transforms and
    passes amplitude (absolute) to the Fourier shell scaling algorithm. Upon return the scaled amplitudes are recombined
    with the phase information of image1. The output is returned in real space.

    @param image1: input image to be scaled, 2d array
    @type  image1: L{np.ndarray}
    @param image2: example image for scaling, 2d array
    @type  image2: L{np.ndarray}
    @param numberBands: number of rings to use for scaling, determines the sharpness of scaling
    @type  numberBands: L{int}

    @return: scaled image1, 2d array
    @rtype:  L{np.ndarray}

    @author: Marten Chaillet
    """
    assert image1.shape == image2.shape, "volumes not of same size"
    assert len(set(image1.shape)) == 1, "volumes are not perfect cubes"

    # scale the amplitudes of 1 to those of 2 using fourier shells
    scaled_amp = FSS(xp.abs(xp.fft.fftshift(xp.fft.fftn(image1))), xp.abs(xp.fft.fftshift(xp.fft.fftn(image2))),
                     numberBands, verbose=False)
    # construct the output volume with the scaled amplitudes and phase infomation of volume 1
    fout = xp.fft.ifftn(xp.fft.ifftshift(scaled_amp) * xp.exp(1j * xp.angle(xp.fft.fftn(image1))))

    return fout.real


def parallel_scale(number, projection, example, pixel_size, example_pixel_size, oversampling, make_even_factor):
    """
    Function that prepares projection and example for fourier shell scaling by ensuring same shape and pixel size.
    It is used in a parallel CPU call.

    @param number: number of the projection image in the series, purely passed for output reasons
    @type  number: L{int}
    @param projection: simulated projection, 2d array
    @type  projection: L{np.ndarray}
    @param example: experimental projection example for scaling, 2d array
    @type  example: L{np.ndarray}
    @param pixel_size: pixel size of simulation in A
    @type  pixel_size: L{float}
    @param example_pixel_size: pixel size experimental in A
    @type  example_pixel_size: L{float}
    @param oversampling: oversampling factor
    @type  oversampling: L{int}
    @param make_even_factor: force the output projection size to be divisible by 2*make_even_factor, i.e. a 1 as input forces the size to be divisble by two
    @type  make_even_factor: L{int}

    @return: (number, scaled image), scaled image is 2d array
    @rtype:  L{tuple} - (L{int}, L{np.ndarray})

    @author: Marten Chaillet
    """
    from pytom.agnostic.transform import resize

    # print(projection.shape, example.shape)

    print(f' -- scaling projection {number+1}')
    if pixel_size != (example_pixel_size * oversampling):
        # magnify or de-magnify if the pixel size does not match yet
        print('(de)magnifying pixel size')
        example = resize(example, (example_pixel_size * oversampling) / pixel_size, interpolation='Spline').squeeze()
        #todo squeeze() can be removed after updating pytom

    # prevent issues later on with oversampling in case experimental and simulated pixel size do not match
    if example.shape[0] % (2*make_even_factor):
        example = example[:-(example.shape[0] % (2*make_even_factor)), :-(example.shape[0] % (2*make_even_factor))]

    if projection.shape != example.shape:
        # crop the largest
        if projection.shape[0] > example.shape[0]:
            cut = (projection.shape[0] - example.shape[0]) // 2
            projection = projection[cut:-cut, cut:-cut]
        else:
            cut = (example.shape[0] - projection.shape[0]) // 2
            example = example[cut:-cut, cut:-cut]

    print('using FSS to scale amplitudes')
    return number, scale_image(projection, example, projection.shape[0] // 4)


def scale_projections(save_path, pixel_size, example_folder, example_pixel_size, oversampling, nodes,
                      make_even_factor):
    """
    Scale the amplitudes of a simulated tilt/frame-series by the amplitudes from experimental images in Fourier space.
    Experimental tilt/frame-series should have been captured under same settings to avoid a mismatch of microscope
    functions, especially CTF rings. This tool is mainly intended to adjust for the increased low spatial frequency
    signal from energy filters.

    If the experimental image has a different pixel size than the simulation, the images will be resized. That is why
    the experimental pixel size is needed as argument here. The oversampling factor is also needed to downsample the
    the data to the binned simulation pixel size.

    @param save_path: simulation project folder
    @type  save_path: L{str}
    @param pixel_size: pixel size of simulated projections in A
    @type  pixel_size: L{float}
    @param example_folder: folder where experimental frame/tilt-series are stored as mrc stacks, a stack will be selected at random
    @type  example_folder: L{str}
    @param example_pixel_size: pixel size of experimental image in A
    @type  example_pixel_size: L{float}
    @param oversampling: the oversampling factor of the simulation
    @type  oversampling: L{int}
    @param nodes: number of CPU nodes to use to parallelize scaling per tilt/frame
    @type  nodes: L{int}
    @param make_even_factor: force the output projection size to be divisible by 2*make_even_factor, i.e. a 1 as input forces the size to be divisble by two
    @type  make_even_factor: L{int}

    @return: - (scaled projections are saved to save_path)
    @rtype:  None

    @author: Marten Chaillet
    """
    from pytom.agnostic.transform import resize
    from pytom.simulation.support import reduce_resolution_fourier
    from pytom.agnostic.io import read_mrc, write
    from joblib import Parallel, delayed

    # generate list of all the possible example projections in the provided parameter example_folder
    files = [f for f in os.listdir(example_folder) if os.path.isfile(os.path.join(example_folder, f))]
    random_file = xp.random.randint(0, len(files))
    print(f'Selected {files[random_file]} for experimental amplitude scaling.')

    filename_pr = os.path.join(save_path, 'projections.mrc')
    filename_ex = os.path.join(example_folder, files[random_file])
    projections         = read_mrc(filename_pr)
    example_projections = read_mrc(filename_ex)

    # assert projections.shape[2] == example_projections.shape[2], 'not enough or too many example projections'
    sim_size = projections.shape
    exp_size = example_projections.shape
    if sim_size[2] > exp_size[2]:  # make sure number of projection angles match
        diff = sim_size[2] - exp_size[2]
        new = xp.zeros((exp_size[0], exp_size[1], sim_size[2]))
        new[..., diff//2: - (diff//2+diff%2)] = example_projections
        new[..., :diff//2] = example_projections[..., :diff//2]
        new[..., -(diff//2+diff%2):] = example_projections[..., -(diff//2+diff%2):]
        example_projections = new
    elif sim_size[2] < exp_size[2]:
        diff = exp_size[2] - sim_size[2]
        example_projections = example_projections[..., diff//2: -(diff//2 + diff%2)]

    # joblib automatically memory maps a numpy array to child processes
    print(f'Scaling projections with {nodes} processes')

    verbosity = 55  # set to 55 for debugging, 11 to see progress, 0 to turn off output
    results = Parallel(n_jobs=nodes, verbose=verbosity, prefer="threads") \
        (delayed(parallel_scale)(i, projections[:, :, i].squeeze(),
                                 resize(reduce_resolution_fourier(example_projections[:, :, i], 1, 2*oversampling),
                                        1/oversampling, interpolation='Spline').squeeze(),
                                 pixel_size, example_pixel_size, oversampling, make_even_factor)
         for i in range(projections.shape[2]))
    # todo squeeze() can be removed after updating pytom

    sys.stdout.flush()

    if results.count(None) == 0:
        print('All scaling processes finished successfully')
    else:
        print(f'{results.count(None)} scaling processes did not finish successfully')

    new_projections = xp.dstack(tuple([r for (i,r) in results]))

    filename_scaled = os.path.join(save_path, 'projections_scaled.mrc')
    write(filename_scaled, new_projections)

    return


def reconstruct_tomogram(save_path, weighting=-1, reconstruction_bin=1,
                         filter_projections=False, use_scaled_projections=False):
    """
    Reconstruction of simulated tilt series into a tomogram. Uses weighted backprojection to make a reconstruction with
    -1, 0, or 1 as a weighting option. -1 is ramp weighting, 0 is no weighting, and 1 is exact weighting.

    Flags can be set for applying a slight low-pass filter (0.9 * image_width) and for using the fourier amplitude
    experimental scaled projections.

    In case binning option is set, or in case scaled projections are reduced in size compared to original projections,
    then the reconstruction will have a different size than the ground truth annotation. In this case we created copies
    of the annotation at the reduced size/scaling.

    @param save_path: simulation project folder
    @type  save_path: L{str}
    @param weighting: type of weighting for WBP (-1, 0, or 1)
    @type  weighting: L{int}
    @param reconstruction_bin: number of times reconstruction should be binned
    @type  reconstruction_bin: L{int}
    @param filter_projections: flag for low-pass filtering images lightly (0.9 width in fourier space)
    @type  filter_projections: L{bool}
    @param use_scaled_projections: flag to reconstruct projections that have amplitude scaled to experimental images
    @type  use_scaled_projections: L{bool}

    @return: -
    @rtype:  None

    @author: Gijs van der Schot, Marten Chaillet
    """
    from pytom.reconstruction.reconstructionStructures import Projection, ProjectionList
    from pytom.agnostic.transform import resize
    from pytom.simulation.support import reduce_resolution_real
    from pytom.agnostic.io import read_mrc, write

    # create folder for individual projections, reconstruction algorithm reads single projections from a folder
    projection_folder = os.path.join(save_path, 'projections')
    if not os.path.exists(projection_folder):
        os.mkdir(projection_folder)

    # select which projections to use, scaled or original
    if use_scaled_projections:
        filename_sc = os.path.join(save_path, 'projections_scaled.mrc')
        projections = read_mrc(filename_sc)
    else:
        filename_pr = os.path.join(save_path, 'projections.mrc')
        projections = read_mrc(filename_pr)

    # possibly apply low pass filter, before storing as single projections to reconstruction folder
    if filter_projections:
        for i in range(projections.shape[2]):
            # 2.3 corresponds to circular filter with width 0.9 of half of the image
            projection_scaled = reduce_resolution_real(projections[:, :, i].squeeze(), 1.0, 2.0 * reconstruction_bin)
            filename_single = os.path.join(save_path, 'projections', f'synthetic_{i+1}.mrc')
            write(filename_single, projection_scaled)
    else:
        for i in range(projections.shape[2]):
            filename_single = os.path.join(save_path, 'projections', f'synthetic_{i+1}.mrc')
            write(filename_single, projections[:, :, i])

    # size and volume shape of the reconstruction
    size_reconstruction = projections.shape[0] // reconstruction_bin
    vol_size = [size_reconstruction, ] * 3
    # alignment file name for reconstruction
    filename_align = os.path.join(save_path, 'alignment_simulated.txt')

    if reconstruction_bin == 1:
        filename_output = os.path.join(save_path, 'reconstruction.em')
    else:
        filename_output = os.path.join(save_path, f'reconstruction_bin{reconstruction_bin}.em')
    # IF EM alignment file provided, filters applied and reconstruction will be identical.
    projections = ProjectionList()
    projections.load_alignment(filename_align)
    vol = projections.reconstructVolume(dims=vol_size, reconstructionPosition=[0, 0, 0], binning=reconstruction_bin,
                                        weighting=weighting)
    vol.write(filename_output)
    os.system(f'em2mrc.py -f {filename_output} -t {os.path.dirname(filename_output)}')
    os.system(f'rm {filename_output}')

    # only if binning during reconstruction is used or the scaled projections are used, ground truth annotation needs
    # to be updated
    if reconstruction_bin > 1 or use_scaled_projections:

        # Adjust ground truth data to match cuts after scaling or after binning for reconstruction
        filename_gm = os.path.join(save_path, 'grandmodel.mrc')
        filename_cm = os.path.join(save_path, 'class_mask.mrc')
        filename_om = os.path.join(save_path, 'occupancy_mask.mrc')
        # todo| New names of grandmodel and class masks should not contain the word bin if they are only fourier
        #  shell scaled. In that case use different identifier (e.g. _scaled)
        filename_gm_bin = os.path.join(save_path, f'grandmodel_bin{reconstruction_bin}.mrc')
        filename_cm_bin = os.path.join(save_path, f'class_mask_bin{reconstruction_bin}.mrc')
        filename_om_bin = os.path.join(save_path, f'occupancy_mask_bin{reconstruction_bin}.mrc')

        cell = read_mrc(filename_gm)
        # find cropping indices
        lind = (cell.shape[0] - reconstruction_bin * size_reconstruction)//2
        rind = cell.shape[0] - lind
        cell = resize(reduce_resolution_real(cell[lind:rind, lind:rind, :], 1, 2*reconstruction_bin), 1/reconstruction_bin,
                      interpolation='Spline')
        write(filename_gm_bin, cell)
        # bin class mask and bbox needed for training
        cell = read_mrc(filename_cm)
        write(filename_cm_bin, downscale_class_mask(cell[lind:rind,lind:rind,:], reconstruction_bin))
        # bin occupancy mask as well
        cell = read_mrc(filename_om)
        write(filename_om_bin, downscale_class_mask(cell[lind:rind,lind:rind,:], reconstruction_bin))

        # create particle locations binned file
        adjusted_ground_truth = ''
        filename_loc = os.path.join(save_path, 'particle_locations.txt')
        filename_loc_bin = os.path.join(save_path, f'particle_locations_bin{reconstruction_bin}.txt')
        with open(filename_loc, 'r') as fin:
            line = fin.readline()
            while line:
                data = line.split()
                data[1] = int(data[1]) // reconstruction_bin - lind
                data[2] = int(data[2]) // reconstruction_bin - lind
                data[3] = int(data[3]) // reconstruction_bin - lind
                # only add the location back if the center of particle is still inside the box after binning
                if 0 <= data[1] < rind and 0 <= data[2] < rind and 0 <= data[3] < rind:
                    data[1] = str(data[1])
                    data[2] = str(data[2])
                    data[3] = str(data[3])
                    adjusted_ground_truth += ' '.join(data) + '\n'
                line = fin.readline()
        with open(filename_loc_bin, 'w') as fout:
            fout.write(adjusted_ground_truth)

    return


if __name__ == '__main__':
    # ------------------------------Import functions used in main-------------------------------------------------------
    # loadstar is for reading .meta files containing data collection parameters (tilt angles, etc.).
    # literal_eval is used for passing arguments from the config file.
    from pytom.gui.guiFunctions import loadstar, datatype
    from ast import literal_eval
    # Use tracemalloc to record the peak memory usage of the program
    tracemalloc.start()

    # --------------------------------------Read config-----------------------------------------------------------------
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

    # ----------------------------------------Set simulation parameters-------------------------------------------------
    try:
        output_folder           = config['General']['OutputFolder']
        simulator_mode          = config['General']['Mode']
        device                  = config['General']['Device']
        nodes                   = config['General'].getint('Nodes')
        model_ID                = config['General'].getint('ModelID')
        seed                    = config['General'].getint('Seed')
        pixel_size              = config['General'].getfloat('PixelSize') * 1E-10 # pixel_size in nm
        oversampling            = config['General'].getint('Oversampling')  # oversampling is used for correcting
        # poisson statistics and camera DQE and MTF functions
        solvent_potential       = config['General'].getfloat('SolventConstant')
        absorption_contrast     = config['General'].getboolean('AbsorptionContrast')
        voltage                 = config['General'].getfloat('Voltage') * 1E3  # voltage in keV
        # voltage and pixelsize are needed for model generation and projection, thus general parameters

        # ensure simulator mode and device are valid options
        if (simulator_mode in ['TiltSeries', 'FrameSeries']) or (device in ['CPU', 'GPU']):
            print(f'Generating model {model_ID} on {device} in folder {output_folder}')
        else:
            print('Invalid entry for simulator mode or device in config.')
            sys.exit(0)
    except Exception as e:
        print(e)
        raise Exception('Missing general parameters in config file.')

    if 'GenerateModel' in config.sections():
        try:
            # We assume the particle models are in the desired voxel spacing for the pixel size of the simulation!
            particle_folder     = config['GenerateModel']['ParticleFolder']
            listpdbs            = literal_eval(config['GenerateModel']['Models'])
            listmembranes       = literal_eval(config['GenerateModel']['MembraneModels'])
            size                = config['GenerateModel'].getint('Size')
            placement_size      = config['GenerateModel'].getint('PlacementSize')
            # parse range of ice thickness, provided in nm
            thickness           = draw_range(literal_eval(config['GenerateModel']['Thickness']), float, 'Thickness') * 1E-9
            thickness_voxels    = int(thickness / pixel_size) # calculate thickness in number of voxels!
            # make even number to solve tomogram reconstruction mismatch bug
            thickness_voxels -= (thickness_voxels % 4)
            # gold markers
            number_of_markers   = draw_range(literal_eval(config['GenerateModel']['NumberOfMarkers']), int,
                                           'NumberOfMarkers')
            # parse range of number of particles
            number_of_particles = draw_range(literal_eval(config['GenerateModel']['NumberOfParticles']), int,
                                             'NumberOfParticles')
            number_of_membranes = draw_range(literal_eval(config['GenerateModel']['NumberOfMembranes']), int,
                                             'NumberOfMembranes')
            sigma_motion_blur   = config['GenerateModel'].getfloat('SigmaMotionBlur')  # in A units
            particle_flipping   = config['GenerateModel']['Mirror']

            # TODO add parameter for meta mode of random variation or stick exactly to input values
        except Exception as e:
            print(e)
            raise Exception('Missing generate model parameters in config file.')

    if 'Microscope' in config.sections():
        try:
            camera                  = config['Microscope']['Camera']
            try:
                camera_folder       = config['Microscope']['CameraFolder']
            except Exception as e:
                camera_folder       = os.path.join(os.path.split(os.path.realpath(__file__))[0], 'detectors')
            # beam damage SNR
            beam_damage_snr         = draw_range(literal_eval(config['Microscope']['BeamDamageSNR']), float,
                                                 'BeamDamageSNR')
            defocus                 = draw_range(literal_eval(config['Microscope']['Defocus']), float, 'Defocus') * 1E-6
            electron_dose           = draw_range(literal_eval(config['Microscope']['ElectronDose']), float, 'ElectronDose')
            spherical_aberration    = config['Microscope'].getfloat('SphericalAberration') * 1E-3
            chromatic_aberration    = config['Microscope'].getfloat('ChromaticAberration') * 1E-3
            energy_spread           = config['Microscope'].getfloat('EnergySpread')
            illumination_aperture   = config['Microscope'].getfloat('IlluminationAperture') * 1E-3
            objective_diameter      = config['Microscope'].getfloat('ObjectiveDiameter') * 1E-6
            focus_length            = config['Microscope'].getfloat('FocalDistance') * 1E-3
            astigmatism             = config['Microscope'].getfloat('Astigmatism') * 1E-9
            astigmatism_angle       = draw_range(literal_eval(config['Microscope']['AstigmatismAngle']), float,
                                                 'AstigmatismAngle')
        except Exception as e:
            print(e)
            raise Exception('Missing microscope parameters in config file.')

    if simulator_mode in config.sections():
        try:
            # first read common parameters between tilt and frame series
            image_size              = config[simulator_mode].getint('ImageSize')
            msdz                    = config[simulator_mode].getfloat('MultisliceStep') * 1E-9
            # random translations between frames/tilts in A
            translation_shift       = draw_range(literal_eval(config[simulator_mode]['TranslationalShift']), float,
                                                 'TranslationalShift')
            # mode specific parameters
            if simulator_mode == 'TiltSeries':
                metadata            = loadstar(config['TiltSeries']['MetaFile'], dtype=datatype)
                angles              = metadata['TiltAngle'] # in degrees
            elif simulator_mode == 'FrameSeries':
                number_of_frames    = config['FrameSeries'].getint('NumberOfFrames')
        except Exception as e:
            print(e)
            raise Exception(f'Missing {simulator_mode} parameters in config file.')

    if 'ScaleProjections' in config.sections():
        try:
            example_folder      = config['ScaleProjections']['ExampleFolder']
            example_pixel_size  = config['ScaleProjections'].getfloat('ExamplePixelSize')
            # If experimental and simulated projections have different size, we need to crop. This should be done with
            # care if the option oversampling is set for reconstructions, because in that case the ground truth data needs
            # to be binned and cropped as well. Uneven size of the volume means the ground truth data will be shifted
            # by half a pixel compared to the reconstruction. This options makes sure that this not happen.
            make_even_factor    = config['ScaleProjections'].getint('EvenSizeFactor')
        except Exception as e:
            print(e)
            raise Exception('Missing experimental projection scaling parameters.')

    if 'TomogramReconstruction' in config.sections():
        try:
            weighting               = config['TomogramReconstruction'].getint('Weighting')
            reconstruction_bin      = config['TomogramReconstruction'].getint('Binning')
            filter_projections      = config['TomogramReconstruction'].getboolean('FilterProjections')
            use_scaled_projections  = config['TomogramReconstruction'].getboolean('UseScaledProjections')
        except Exception as e:
            print(e)
            raise Exception('Missing tomogram reconstruction parameters in config file.')

    # --------------------------------------Create directories and logger-----------------------------------------------
    save_path = os.path.join(output_folder, f'model_{model_ID}')
    if not os.path.exists(save_path):
        os.mkdir(save_path)

    logging.basicConfig(filename='{}/simulator-{date:%Y-%m-%d_%H:%M:%S}.log'.format(save_path,
                                                                date=datetime.datetime.now()), level=logging.INFO)
    config_logger = ConfigLogger(logging)
    config_logger(config)

    logging.info('Values of parameters that we randomly vary per simulation (only for generate model and generate projections):')
    if 'GenerateModel' in config.sections():
        logging.info(f'model thickness = {thickness_voxels*pixel_size*1E9:.2f}nm (adjusted to be an even number of voxels)')
        logging.info(f'# of particles = {number_of_particles}')
        logging.info(f'# of markers = {number_of_markers}')
        logging.info(f'# of membranes = {number_of_membranes}')
    if 'Microscope' in config.sections():
        logging.info(f'defocus = {defocus*1E6:.2f}um')
        logging.info(f'electron dose = {electron_dose} e-/A^2')

    # ----------------------------------------Execute simulation--------------------------------------------------------
    if 'GenerateModel' in config.sections():
        # set seed for random number generation
        xp.random.seed(seed)
        random.seed(seed)

        print('\n- Generating grand model')
        generate_model(particle_folder, save_path, listpdbs, listmembranes,
                       pixel_size           =pixel_size * 1E10,
                       size                 =size,
                       thickness            =thickness_voxels,
                       placement_size       =placement_size,
                       solvent_potential    =solvent_potential,
                       number_of_particles  =number_of_particles,
                       number_of_markers    =number_of_markers,
                       absorption_contrast  =absorption_contrast,
                       voltage              =voltage,
                       number_of_membranes  =number_of_membranes,
                       sigma_motion_blur    =sigma_motion_blur,
                       particle_flipping    =particle_flipping)

    if simulator_mode in config.sections() and simulator_mode == 'TiltSeries':
        # set seed for random number generation
        xp.random.seed(seed)
        random.seed(seed)
        # Grab the ice thickness from the initial model in case program is only executed for projections
        print('\n- Generating projections')
        if device == 'CPU':
            generate_tilt_series_cpu(save_path, angles,
                                      nodes                 =nodes,
                                      image_size            =image_size,
                                      rotation_box_height   =None, # will automatically calculate fitting size if None
                                      pixel_size            =pixel_size,
                                      oversampling          =oversampling,
                                      dose                  =electron_dose,
                                      voltage               =voltage,
                                      spherical_aberration  =spherical_aberration,
                                      chromatic_aberration  =chromatic_aberration,
                                      energy_spread         =energy_spread,
                                      illumination_aperture =illumination_aperture,
                                      objective_diameter    =objective_diameter,
                                      focus_length          =focus_length,
                                      astigmatism           =astigmatism,
                                      astigmatism_angle     =astigmatism_angle,
                                      msdz                  =msdz,
                                      defocus               =defocus,
                                      sigma_shift           =translation_shift,
                                      camera_type           =camera,
                                      camera_folder         =camera_folder,
                                      solvent_potential     =solvent_potential,
                                      absorption_contrast   =absorption_contrast,
                                      beam_damage_snr       =beam_damage_snr)
        elif device == 'GPU':
            print('This option needs to be implemented.')
            sys.exit(0)
        else:
            print('Invalid device type.')
            sys.exit(0)
    elif simulator_mode in config.sections() and simulator_mode == 'FrameSeries':
        xp.random.seed(seed)
        random.seed(seed)
        print('\n- Generate frame series projections')
        if device == 'CPU':
            generate_frame_series_cpu(save_path,
                                      n_frames              =number_of_frames,
                                      nodes                 =nodes,
                                      image_size            =image_size,
                                      pixel_size            =pixel_size,
                                      oversampling          =oversampling,
                                      dose                  =electron_dose,
                                      voltage               =voltage,
                                      spherical_aberration  =spherical_aberration,
                                      chromatic_aberration  =chromatic_aberration,
                                      energy_spread         =energy_spread,
                                      illumination_aperture =illumination_aperture,
                                      objective_diameter    =objective_diameter,
                                      focus_length          =focus_length,
                                      astigmatism           =astigmatism,
                                      astigmatism_angle     =astigmatism_angle,
                                      msdz                  =msdz,
                                      defocus               =defocus,
                                      mean_shift            =translation_shift,
                                      camera_type           =camera,
                                      camera_folder         =camera_folder,
                                      solvent_potential     =solvent_potential,
                                      absorption_contrast   =absorption_contrast,
                                      beam_damage_snr       =beam_damage_snr)
        elif device == 'GPU':
            print('This option needs to be implemented.')
            sys.exit(0)
        else:
            print('Invalid device type.')
            sys.exit(0)

    if 'ScaleProjections' in config.sections():
        # set seed for random number generation
        xp.random.seed(seed)
        random.seed(seed)
        print('\n- Scaling projections with experimental data')
        scale_projections(save_path, pixel_size * 1E10, example_folder,
                                            example_pixel_size, oversampling, nodes, make_even_factor)

    if 'TomogramReconstruction' in config.sections():
        print('\n- Reconstructing tomogram')
        reconstruct_tomogram(save_path,
                             weighting=weighting,
                             reconstruction_bin=reconstruction_bin,
                             filter_projections=filter_projections,
                             use_scaled_projections=use_scaled_projections)

    current, peak = tracemalloc.get_traced_memory()
    print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
    tracemalloc.stop()

