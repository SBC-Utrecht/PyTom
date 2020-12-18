"""
Marten Chaillet's cryoET simulator - updated version of the simulator used in SHREC2019 and SHREC2020
Contributions from Ilja Gubins and Gijs van der Schot.
Original simulator (used in SHREC2019) written by Gijs van der Schot, which was loosely based on the simulator in the
TOM toolbox for matlab.
"""

# IO related modules
from pytom.gui.mrcOperations import *
import pytom.tompy.io
import configparser
import tracemalloc
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
# use Agg for plotting without display, while only writing to files
# matplotlib.use('Agg')
import matplotlib.pylab as plt

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


def draw_range(range, datatype, name):
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


def create_ellipse(size, mj, mn1, mn2, smooth=0, cutoff_SD=3):
    """
    Generate an ellipse defined by 3 radii along x,y,z - parameters mj, mn1, mn2.
    """
    X,Y,Z = xp.meshgrid(xp.arange(size/1), xp.arange(size/1), xp.arange(size/1))

    X -= size/2-0.5
    Y -= size/2-0.5
    Z -= size/2-0.5

    R = xp.sqrt( (X/mj)**2 + (Y/mn1)**2 + (Z/mn2)**2)

    # print(R.max(), R.min())

    out = xp.zeros((size,size,size),dtype=xp.float32)
    out[ R <= 1] = 1

    if smooth:
        R2 = R.copy()
        R2[R <= 1] = 1
        sphere = xp.exp(-1 * ((R2-1)/smooth)**2)
        sphere[sphere <= xp.exp(-cutoff_SD**2/2.)] = 0
        out = sphere

    return out


def addNoise(noise_size, dim):
    """
    Correlated noise to create density deformations.
    """
    from numpy.fft import ifftn, fftn, fftshift
    from numpy.random import random

    if noise_size == 0:
        # 0 value will not work
        noise_size = 1

    noise_no_norm = abs(ifftn(fftshift(fftn(random((noise_size, noise_size, noise_size)))), [dim] * 3))
    noise = 0.2 * noise_no_norm / abs(noise_no_norm).max()

    return 1 + (noise - noise.mean())


def create_gold_marker(voxel_size, solvent_potential, binning=1, solvent_factor=1.0, imaginary=False, voltage=300E3):
    """
    From Rahman 2018 (International Journal of Biosensors and Bioelectronics).
    Volume of unit cell gold is 0.0679 nm^3 with 4 atoms per unit cell.
    Volume of gold bead is 4/3 pi r^3.
    """
    from pytom.tompy.tools import create_sphere
    from potential import reduce_resolution, bin
    # from pytom.tompy.filter import gaussian3d

    assert (type(binning) is int) and (binning >= 1), print('Stop gold marker creation binning factor is not a positive'
                                                            ' integer.')

    # select a random size for the gold marker in nm
    diameter = xp.random.uniform(low=4.0, high=10.0)

    # constants
    unit_cell_volume = 0.0679 # nm^3
    atoms_per_unit_cell = 4
    C = 2 * xp.pi * phys.constants['h_bar']**2 / (phys.constants['el'] * phys.constants['me']) * 1E20  # nm^2
    voxel_size_nm = voxel_size*1E9 / binning
    voxel_volume = voxel_size_nm**3

    # dimension of gold box, always add 5 nm to the sides
    dimension = int(xp.ceil(diameter / voxel_size_nm)) * 3
    # sigma half of radius?
    r = 0.8 * ((diameter * 0.5) / voxel_size_nm) # fraction of radius due to extension with exponential smoothing
    ellipse = True
    if ellipse:
        r2 = r * xp.random.uniform(0.8, 1.2)
        r3 = r * xp.random.uniform(0.8, 1.2)
        bead = create_ellipse(dimension, r, r2, r3, smooth=2)
    else:
        bead = create_sphere((dimension,)*3, radius=r)

    bead *= addNoise(int(r*0.75), dimension) * addNoise(int(r*0.25), dimension)
    # SIGMA DEPENDENT ON VOXEL SIZE
    # rounded_sphere = gaussian3d(sphere, sigma=(1 * 0.25 / voxel_size_nm))
    bead[bead < 0.9] = 0 # remove too small values
    # add random noise to gold particle to prevent perfect CTF ringing around the particle.
    # random noise also dependent on voxel size maybe?
    # rounded_sphere = (rounded_sphere > 0) * (rounded_sphere * xp.random.normal(1, 0.3, rounded_sphere.shape))
    # rounded_sphere[rounded_sphere < 0] = 0

    if imaginary:
        solvent_amplitude = potential_amplitude(0.93, 18, voltage) * solvent_factor
        gold_amplitude = potential_amplitude(19.3, 197, voltage)
        gold_imaginary = bead * (gold_amplitude - solvent_amplitude)
        # filter and bin
        gold_imaginary = bin(gold_imaginary, binning)
        # gold_imaginary = bin(reduce_resolution(gold_imaginary, voxel_size_nm, voxel_size_nm*2), binning)

    # values transformed to occupied volume per voxel from 1 nm**3 per voxel to actual voxel size
    solvent_correction = bead * (solvent_potential * solvent_factor)
    unit_cells_per_voxel = (bead * voxel_volume / unit_cell_volume)
    gold_atoms = unit_cells_per_voxel * atoms_per_unit_cell

    # interaction potential
    gold_scattering_factors = xp.array(phys.scattering_factors['AU']['g'])
    # gold_scattering_factors[0:5].sum() == 10.57
    # C and scattering factor are in A units thus divided by 1000 A^3 = 1 nm^3 to convert
    gold_potential = gold_atoms * gold_scattering_factors[0:5].sum() * C / voxel_volume / 1000
    gold_real = gold_potential - solvent_correction
    # filter and bin
    # gold_real = bin(reduce_resolution(gold_real, voxel_size_nm, voxel_size_nm * 2), binning)
    gold_real = bin(gold_real, binning)

    if imaginary:
        return gold_real, gold_imaginary
    else:
        return gold_real


def generate_model(particleFolder, output_folder, model_ID, listpdbs, listmembranes, pixel_size = 1E-10, binning=1,
                   size=1024, thickness=200,
                   solvent_potential=4.5301, solvent_factor=1.0, numberOfParticles=1000,
                   placement_size=512, retries=5000, number_of_markers=0,
                   absorption_contrast=False, voltage=300E3, number_of_membranes=0):
    # IMPORTANT: We assume the particle models are in the desired voxel spacing for the pixel size of the simulation!
    from scipy.ndimage import gaussian_filter
    from voltools import transform

    # outputs
    X, Y, Z = size, size, thickness
    cell = xp.zeros((X, Y, Z))
    if absorption_contrast: cell_imag = xp.zeros_like(cell)

    occupancy_bbox_mask = xp.zeros_like(cell)
    occupancy_accurate_mask = xp.zeros_like(cell)
    class_bbox_mask = xp.zeros_like(cell)
    class_accurate_mask = xp.zeros_like(cell)
    ground_truth_txt_file = ''

    # load pdb volumes and pad them
    volumes = []
    if absorption_contrast: volumes_imag = []
    for pdb in listpdbs:
        try:
            # First find the voxel size of the map...
            # rotate
            # then bin
            if absorption_contrast:
                vol_real = pytom.tompy.io.read_mrc(f'{particleFolder}/{pdb}_{pixel_size*1E10:.2f}A_solvent-'
                                              f'{solvent_potential*solvent_factor:.3f}V_real.mrc')
                vol_imag = pytom.tompy.io.read_mrc(f'{particleFolder}/{pdb}_{pixel_size*1E10:.2f}A_solvent-'
                                                   f'{solvent_potential*solvent_factor:.3f}V_imag.mrc')
                assert vol_real.shape == vol_imag.shape, print(f'real and imaginary interaction potential not the same '
                                                               f'shape for {pdb}')
                # dx, dy, dz = vol_real.shape
                # vol2_real = xp.zeros((dx*2, dy*2, dz*2), dtype=xp.float32) # do not double
                # vol2_imag = xp.zeros_like(vol2_real)
                # vol2_real[dx // 2:-dx // 2, dy // 2:-dy // 2, dz // 2:-dz // 2] = vol_real
                # vol2_imag[dx // 2:-dx // 2, dy // 2:-dy // 2, dz // 2:-dz // 2] = vol_imag
                volumes.append(vol_real)
                volumes_imag.append(vol_imag)
            else:
                vol = pytom.tompy.io.read_mrc(f'{particleFolder}/{pdb}_{pixel_size*1E10:.2f}A_solvent-'
                                              f'{solvent_potential*solvent_factor:.3f}V_real.mrc')
                # dx, dy, dz = vol.shape
                # vol2 = xp.zeros((dx*2, dy*2, dz*2), dtype=xp.float32)
                # vol2[dx//2:-dx//2, dy//2:-dy//2, dz//2:-dz//2] = vol
                volumes.append(vol)
        except Exception as ee:
            print(ee)
            raise Exception('Could not open pdb ', pdb)

    # attributes
    number_of_classes = len(listpdbs)
    save_path = f'{output_folder}/model_{model_ID}'
    dims = [v.shape for v in volumes] # TODO this parameter is not needed, shapes can be extracted later
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
        particles_by_class += [0]

        # class id of cell structures is always the same
        cls_id = number_of_classes - 1

        for _ in tqdm(range(number_of_membranes), desc='Placing membranes and micelles'):

            # Give membranes a numbered names to randomly index them
            membrane_model = listmembranes[xp.random.randint(0, len(listmembranes))]

            # create the gold marker in other function
            if absorption_contrast:
                membrane_vol = pytom.tompy.io.read_mrc(f'{particleFolder}/{membrane_model}_{pixel_size*1E10:.2f}A_'
                                                   f'solvent-{solvent_potential*solvent_factor:.3f}V_real.mrc')
                membrane_vol_imag = pytom.tompy.io.read_mrc(f'{particleFolder}/{membrane_model}_{pixel_size*1E10:.2f}A_'
                                        f'solvent-{solvent_potential*solvent_factor:.3f}V_imag_{voltage*1E-3:.0f}V.mrc')
            else:
                membrane_vol = pytom.tompy.io.read_mrc(f'{particleFolder}/{membrane_model}_{pixel_size*1E10:.2f}A_'
                                                   f'solvent-{solvent_potential*solvent_factor:.3f}V_real.mrc')

            u = xp.random.uniform(0.0, 1.0, (2,))
            theta = xp.arccos(2 * u[0] - 1)
            phi = 2 * xp.pi * u[1]
            psi = xp.random.uniform(0.0, 2 * xp.pi)
            p_angles = xp.rad2deg([theta, phi, psi])

            # randomly mirror and rotate particle
            try:
                membrane_vol = transform(membrane_vol, rotation=p_angles,
                                        rotation_order='szxz', interpolation='linear', device='cpu')
                if absorption_contrast: membrane_vol_imag = transform(membrane_vol_imag, rotation=p_angles,
                                        rotation_order='szxz', interpolation='linear', device='cpu')
                # do something to remove edge artifacts of rotation? linear instead of filt_bspline
            except Exception as e:
                print(e)
                print('Something went wrong while rotating?')
                continue

            threshold = 0.001
            membrane_vol[membrane_vol < threshold] = 0
            if absorption_contrast: membrane_vol_imag[membrane_vol_imag < threshold] = 0

            # thresholded particle
            accurate_particle_occupancy = membrane_vol > 0

            # find random location for the particle
            # allow membranes to be placed half outside of the grandcell
            dimensions = membrane_vol.shape
            xx, yy, zz = dimensions
            tries_left = default_tries_left
            # x_cut_left, x_cut_right = 0,0
            y_cut_left, y_cut_right, z_cut_low, z_cut_high = 0,0,0,0
            while tries_left > 0:
                loc_x = xp.random.randint(loc_x_start + xx // 2 + 1, loc_x_end - xx // 2 - 1)
                # loc_y = xp.random.randint(loc_y_start + yy // 2 + 1, loc_y_end - yy // 2 - 1)
                # loc_z = xp.random.randint(zz // 2 + 1, Z - zz // 2 - 1)
                # loc_x = xp.random.randint(loc_x_start, loc_x_end)
                loc_y = xp.random.randint(loc_y_start, loc_y_end)
                loc_z = xp.random.randint(0, Z)

                tries_left -= 1

                # adjust for bbox indexing for membrane sticking out of box
                # for x limit between x_start and x_end
                # if loc_x - xx//2 < 0: x_cut_left = abs(loc_x - xx//2)
                if loc_y - yy//2 < 0: y_cut_left = abs(loc_y - yy//2)
                if loc_z - zz//2 < 0: z_cut_low = abs(loc_z - zz//2)
                # if loc_x + xx//2 + xx%2 > loc_x_end: x_cut_right = loc_x - (loc_x + xx//2 + xx%2)
                if loc_y + yy//2 + yy%2 > Y: y_cut_right = abs(Y - (loc_y + yy//2 + yy%2))
                if loc_z + zz//2 + zz%2 > Z: z_cut_high = abs(Z - (loc_z + zz//2 + zz%2))

                # crop the occupancy by removing overhang
                accurate_particle_occupancy_crop = accurate_particle_occupancy[:,
                                                            y_cut_left:yy-y_cut_right, z_cut_low:zz-z_cut_high]

                # calculate coordinates of bbox for the newly rotated particle
                bbox_x = [loc_x - xx // 2, loc_x + xx // 2 + xx % 2]
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
            occupancy_bbox_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] = particle_nr
            occupancy_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
                accurate_particle_occupancy_crop * particle_nr

            # populate class masks
            class_bbox_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] = (cls_id + 1)
            class_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
                accurate_particle_occupancy_crop * (cls_id + 1)

            # populate density volume
            cell[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
                membrane_vol[:, y_cut_left:yy-y_cut_right, z_cut_low:zz-z_cut_high]
            cell_imag[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
                membrane_vol_imag[:, y_cut_left:yy-y_cut_right, z_cut_low:zz-z_cut_high]

            # update stats
            particle_nr += 1
            particles_by_class[cls_id] += 1

            # update text
            ground_truth_txt_file += f'fiducial {int(loc_x - loc_x_start)} {int(loc_y - loc_y_start)} {int(loc_z)} ' \
                                     f'NaN NaN NaN\n'

    # GOLD MARKERS WILL ALSO BE COUNTER TOWARDS TOTAL PARTICLE NUMBER
    # There are also added as an additional class
    if number_of_markers:
        number_of_classes += 1
        particles_by_class += [0]

        # class id of gold markers is always the same
        cls_id = number_of_classes - 1

        for _ in tqdm(range(number_of_markers), desc='Placing gold markers'):

            # create the gold marker in other function
            if absorption_contrast:
                gold_marker, gold_imag = create_gold_marker(pixel_size, solvent_potential, binning=binning,
                                                            solvent_factor=solvent_factor,
                                                            imaginary=True, voltage=voltage)
            else:
                gold_marker = create_gold_marker(pixel_size, solvent_potential, binning=binning,
                                                 solvent_factor=solvent_factor)

            u = xp.random.uniform(0.0, 1.0, (2,))
            theta = xp.arccos(2 * u[0] - 1)
            phi = 2 * xp.pi * u[1]
            psi = xp.random.uniform(0.0, 2 * xp.pi)
            p_angles = xp.rad2deg([theta, phi, psi])

            # randomly mirror and rotate particle
            try:
                gold_marker = transform(gold_marker, rotation=p_angles,
                                        rotation_order='szxz', interpolation='linear', device='cpu')
                if absorption_contrast: gold_imag = transform(gold_imag, rotation=p_angles,
                                        rotation_order='szxz', interpolation='linear', device='cpu')
                # do something to remove edge artifacts of rotation? linear instead of filt_bspline
            except Exception as e:
                print(e)
                print('Something went wrong while rotating?')
                continue

            threshold = 0.001
            gold_marker[gold_marker < threshold] = 0
            if absorption_contrast: gold_imag[gold_imag < threshold] = 0

            # thresholded particle
            accurate_particle_occupancy = gold_marker > 0

            # find random location for the particle
            dimensions = gold_marker.shape
            xx, yy, zz = dimensions
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
            if absorption_contrast: cell_imag[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += gold_imag

            # update stats
            particle_nr += 1
            particles_by_class[cls_id] += 1

            # update text
            ground_truth_txt_file += f'fiducial {int(loc_x - loc_x_start)} {int(loc_y - loc_y_start)} {int(loc_z)} ' \
                                     f'NaN NaN NaN\n'

    for _ in tqdm(range(numberOfParticles), desc='Placing particles'):

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
            particle = volumes[cls_id]
            if absorption_contrast: particle_imag = volumes_imag[cls_id]
            if xp.random.randint(2): # Generate true/false randomly
                # Mirror the particle to cover both left and right handedness of the proteins
                ax = xp.random.randint(3)
                particle = xp.flip(particle, axis=ax)
                if absorption_contrast: particle_imag = xp.flip(particle_imag, axis=ax)
            rotated_particle = transform(particle, rotation=p_angles,
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
        threshold = 0.001
        rotated_particle[rotated_particle < threshold] = 0
        if absorption_contrast: rotated_particle_imag[rotated_particle_imag < threshold]

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

            # TODO is the dims variable needed here? Could just use xx,yy,zz because the shape did not change
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
        if absorption_contrast: cell_imag[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
                                                    rotated_particle_imag

        # update stats
        particle_nr += 1
        particles_by_class[cls_id] += 1

        # update text
        ground_truth_txt_file += f'{listpdbs[cls_id]} {int(loc_x - loc_x_start)} {int(loc_y - loc_y_start)} {int(loc_z)} ' \
                                 f'{p_angles[0]:.4f} {p_angles[1]:.4f} {p_angles[2]:.4f}\n'

    # # add solvent background potential
    # print('Adding background solvent potential')
    # cell += solvent_potential
    #
    # # Add structural nois
    # print('Adding structural noise to grand model cell')
    # noisy_cell = cell #+ xp.random.normal(0, sigma_structural, cell.shape)
    # pytom.tompy.io.write(f'{save_path}/grandmodel_noisefree_original.mrc', cell)

    # motion blur is simply a gaussian in fourier space to the required frequency component
    # sigma motion is roughly 8A
    sigma_motion = 8
    gaussian_filter(cell, sigma = (sigma_motion/(2 * (pixel_size * 1E10)))/2.35, output=cell)
    gaussian_filter(cell_imag, sigma=(40 / (2 * (pixel_size * 1E10))) / 2.35, output=cell_imag)

    # save grandmodels
    print('Saving grandmodels')
    pytom.tompy.io.write(f'{save_path}/grandmodel_original.mrc', cell)
    if absorption_contrast: pytom.tompy.io.write(f'{save_path}/grandmodel_original_imag.mrc', cell_imag)

    # save class masks
    print('Saving class volumes')
    pytom.tompy.io.write(f'{save_path}/class_bbox_original.mrc', class_bbox_mask)
    pytom.tompy.io.write(f'{save_path}/class_mask_original.mrc', class_accurate_mask)

    # save occupancy masks
    print('Saving occupancy volumes')
    pytom.tompy.io.write(f'{save_path}/occupancy_bbox_original.mrc', occupancy_bbox_mask)
    pytom.tompy.io.write(f'{save_path}/occupancy_mask_original.mrc', occupancy_accurate_mask)

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
    del cell, class_accurate_mask, class_bbox_mask, occupancy_accurate_mask, occupancy_bbox_mask
    if absorption_contrast: del cell_imag
    return


def parallel_rotate_model(volume, outname, angle):
    print(f'Starting rotation process for angle {angle}')
    sys.stdout.flush()
    from voltools import transform
    # volume = pytom.tompy.io.read_mrc(filename)
    rotated_volume = transform(volume, rotation=(0, angle, 0), rotation_order='sxyz', interpolation='filt_bspline', device='cpu')
    threshold = 0.001
    rotated_volume[rotated_volume < threshold] = 0
    pytom.tompy.io.write(outname, rotated_volume)
    print(f'Process for angle {angle} is finished ({outname})')
    sys.stdout.flush()
    return True


def create_rotation_model(output_folder, model_ID, heightBox, imaginary=False):

    # Load grandmodel
    save_path = f'{output_folder}/model_{model_ID}'
    if imaginary:
        file = f'{save_path}/grandmodel_original_imag.mrc'
    else:
        file = f'{save_path}/grandmodel_original.mrc'

    if os.path.exists(file):
        grandcell = pytom.tompy.io.read_mrc(file)
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
    if imaginary:
        filename = f'{dir}/rotated_volume_0_imag.mrc'
    else:
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


def potential_amplitude(rho, molecular_weight, voltage):
    """
    Calculate the inelastic Mean Free Path for different environments.
    @param rho:
    @param molecular_weight:
    @param voltage:
    @return:
    """

    # if Mw == 18 % water
    #     ZO = 8;
    #     sigma_inH = 8.8e-6 * beta2_100. * log(beta2. * (U0 + E0) / (E_loss / 2)) / (
    #                 beta2. * log(beta2_100. * (100e3 + E0) / (E_loss / 2)));
    #     sigma_inO = 1.5 * 1e-6 * ZO ^ 0.5. / beta2. * log(beta2. * (U0 + E0) / (E_loss / 2));
    #     sigma_in = 2 * sigma_inH + sigma_inO;
    #
    # elseif
    # Mw == 12 % carbon
    # ZC = 6;
    # sigma_in = 1.5 * 1e-6 * ZC ^ 0.5. / beta2. * log(beta2. * (U0 + E0) / (E_loss / 2));
    #
    # elseif
    # Mw == 7.2 % protein
    # sigma_in = 0.82 * 1e-4 * beta2_100. * log(beta2. * (U0 + E0) / (E_loss / 2)) / (
    #             beta2. * log(beta2_100. * (100e3 + E0) / (E_loss / 2)));
    # end
    #
    # Nconc = rho * 1000 * nc.Na / (Mw / 1000);
    # Lambda_in = 1. / (Nconc * sigma_in * 1e-18);

    Lambda = wavelength_eV2m(voltage)
    relative_mass = phys.constants["me"] + phys.constants["el"] * voltage / (phys.constants["c"] ** 2)
    sigma_transfer = 2 * xp.pi * relative_mass * phys.constants["el"] * Lambda / (phys.constants["h"] ** 2)
    E0 = phys.constants['me'] * phys.constants['c'] ** 2 / phys.constants['el'] # rest mass energy
    E_loss = 20 # eV mean plasmon loss

    beta2 = 1 - (E0 / (voltage + E0)) ** 2
    beta2_100 = 1 - (E0 / (100E3 + E0)) ** 2

    if molecular_weight == 18: # water molecule
        # rho for amorphous ice is 0.93 g/cm^3
        ZO = 8
        sigma_inelastic_H = ( 8.8E-6 * beta2_100 * xp.log(beta2 * (voltage + E0) / (E_loss / 2)) /
                              (beta2 * xp.log(beta2_100 * (100E3 + E0) / (E_loss / 2))) )
        sigma_inelastic_O = 1.5 * 1E-6 * ZO ** 0.5 / beta2 * xp.log(beta2 * (voltage+E0) / (E_loss/2))
        sigma_inelastic = 2 * sigma_inelastic_H + sigma_inelastic_O
    elif molecular_weight == 12: # carbon film
        ZC = 6
        sigma_inelastic = 1.5 * 1E-6 * ZC ** 0.5 / beta2 * xp.log(beta2 * (voltage + E0) / (E_loss / 2))
    elif molecular_weight == 7.2 or molecular_weight == 734.1: # protein or lipid
        # rho for proteins is assumed to be 1.35 g/cm^3
        # fractional composition is 0.492, 0.313, 0.094, and 0.101 for elements H, C, N, and O, respectively
        sigma_inelastic = (0.82 * 1E-4 * beta2_100 * xp.log(beta2 * (voltage + E0) / (E_loss / 2)) /
                           (beta2 * xp.log(beta2_100 * (100E3 + E0) / (E_loss / 2))) )
    elif molecular_weight == 197: # gold particles
        ZAu = 79
        sigma_inelastic = 1.5 * 1E-6 * ZAu ** 0.5 / beta2 * xp.log(beta2 * (voltage + E0) / (E_loss / 2))
    # elif molecular_weight == 734.1:
    #     # DPPC lipid composition C40H80NO8P, molar mass 734.053 g/mol
    #     # rho is 0.92 g/cm^3 for most vegetable oils, find a better reference!
    #     ZC = 6
    #     ZN = 7
    #     ZO = 8
    #     ZP = 15
    #     sigma_inelastic_H = (8.8E-6 * beta2_100 * xp.log(beta2 * (voltage + E0) / (E_loss / 2)) /
    #                          (beta2 * xp.log(beta2_100 * (100E3 + E0) / (E_loss / 2))))
    #     sigma_inelastic_C = 1.5 * 1E-6 * ZC ** 0.5 / beta2 * xp.log(beta2 * (voltage + E0) / (E_loss / 2))
    #     sigma_inelastic_N = 1.5 * 1E-6 * ZN ** 0.5 / beta2 * xp.log(beta2 * (voltage + E0) / (E_loss / 2))
    #     sigma_inelastic_O = 1.5 * 1E-6 * ZO ** 0.5 / beta2 * xp.log(beta2 * (voltage + E0) / (E_loss / 2))
    #     sigma_inelastic_P = 1.5 * 1E-6 * ZP ** 0.5 / beta2 * xp.log(beta2 * (voltage + E0) / (E_loss / 2))
    #     sigma_inelastic = 40 * sigma_inelastic_C + 80 * sigma_inelastic_H + sigma_inelastic_N + 8 * sigma_inelastic_O + \
    #         sigma_inelastic_P

    concentration = rho * 1000 * phys.constants['na'] / (molecular_weight / 1000)
    lamda_inelastic = ( 1. / (concentration * sigma_inelastic * 1E-18) )

    return 1.0 / (2 * sigma_transfer * lamda_inelastic)


def create_ice_layer(shape, angle, width, value=1.0, sigma=0.0):
    """
    Efficiently create an ice layer at specified angle within the volume shape.
    Assumes the angle rotates perpendicular to the y-axis of the volume.

    @param volume: volume to add the ice to
    @param angle: rotation angle of ice layer
    @param width: width of the ice layer in number of pixels
    @param value: value of ice layer
    @param sigma: std of gaussian filter for smoothing edges

    @author: Marten Chaillet
    """
    from scipy.ndimage import gaussian_filter
    assert xp.abs(angle) <= 90, print('rotation angle of ice layer cannot be larger than +- 90 degrees.')

    xsize = shape[0]
    ysize = shape[1]
    zsize = shape[2]

    x = xp.arange(-xsize / 2, xsize / 2, 1, dtype=xp.float32)
    z = xp.arange(-zsize / 2, zsize / 2, 1, dtype=xp.float32)
    zm = xp.tile(z[xp.newaxis, :], [xsize, 1])

    xline = x * xp.tan(-angle * xp.pi / 180)
    nwidth = width / xp.cos(angle * xp.pi / 180)
    xmin = xline - nwidth / 2
    xmax = xline + nwidth / 2

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
    grad2 = (grad2 - grad2.min()) / (grad2.max() - grad2.min())

    # if a sigma is provided apply gaussian filter
    if not (sigma == 0.0):
        layer = gaussian_filter(grad1 * grad2 * value, sigma=sigma)
    else:
        layer = grad1 * grad2 * value

    # before returning tile the 2d layer to a volume
    return xp.tile(layer[:, xp.newaxis, :], [1, ysize, 1])


def radial_average(image):
    assert len(set(image.shape)) == 1, 'differently size dimension, cannot perform radial averaging'
    size = image.shape[0]
    center = (size - 1) / 2
    x, y = xp.meshgrid(xp.arange(size), xp.arange(size))
    R = xp.sqrt((x - center) ** 2 + (y - center) ** 2)

    f = lambda r: image[(R >= r - .5) & (R < r + .5)].mean()
    r = xp.linspace(1, size // 2, num=size // 2)
    mean = xp.vectorize(f)(r)

    # plot it
    # fig, ax = plt.subplots()
    # ax.plot(r, mean)
    # plt.show()
    return r, mean


def create_complex_CTF(image_shape, pix_size, Dz, voltage=300E3, Cs=2.7E-3, sigma_decay_CTF=0.4, display_CTF=False):
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

    # Decay function
    if sigma_decay_CTF:
        decay = xp.exp(-(r / (sigma_decay_CTF * Ny)) ** 2)
        CTF = complex_CTF * decay
    else:
        CTF = complex_CTF

    # visual inspection of CTF
    if display_CTF:
        r1, m1 = radial_average(CTF.real)
        r2, m2 = radial_average(CTF.imag)
        fig, ax = plt.subplots()

        ax.plot(r1, m1, label='real')
        ax.plot(r2, m2, label='imaginary')
        ax.legend()
        show()

    return CTF


def create_complex_CTF_ext(image_shape, pixel_size, defocus, voltage=300E3, Cs=2.7E-3, Cc=2.7E-3,
                           energy_spread=0.7, illumination_aperture=0.030E-3, objective_diameter=100E-6,
                           focus_length=4.7E-3, astigmatism=0.0, astigmatism_angle=0.0, display_CTF=False):
    """
        # parameters for extended CTF function (InSilicoTEM)
    chromatic_abberation    = 2.7E-3 # C_c
    energy_spread           = 0.7 # deltaE
    illumination_aperture   = 0.030E-3 # a_i
    objective_diameter      = 100E-6 #
    focus_length            = 4.7E-3 # focal distance
    astigmatism in 0.0E-9
    astigmatism angle in degrees

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
    from scipy.ndimage import gaussian_filter

    assert len(set(image_shape)) == 1, print('invalid input image/volume for create CTF, dimensions need to be equal.')

    image_size = image_shape[0]

    Lambda = wavelength_eV2m(voltage)

    q_true = 1 / (image_size * pixel_size)

    dfsh = defocus + astigmatism
    dfln = defocus - astigmatism
    inratioqqs = xp.sqrt(dfsh/defocus)
    inratioqlq = xp.sqrt(defocus/dfln)

    x = xp.arange(-(image_size//2), abs(-image_size//2), 1)
    xx, yy = xp.meshgrid(x, x)
    xdot = xx * xp.cos(astigmatism_angle*xp.pi/180) - yy * xp.sin(astigmatism_angle*xp.pi/180)
    ydot = xx * xp.sin(astigmatism_angle*xp.pi/180) - yy * xp.cos(astigmatism_angle*xp.pi/180)
    q = xp.sqrt((xdot / inratioqlq) ** 2 + (ydot * inratioqqs) ** 2) * q_true
    qsym = xp.sqrt(xdot ** 2 + ydot ** 2) * q_true

    # print(r)
    c = 0.5 * xp.pi * (Cs * (Lambda ** 3) * (qsym ** 4) - 2 * defocus * Lambda * (q ** 2))
    complex_CTF = xp.cos(c) - 1j * xp.sin(c)

    # chromatic envelope
    H = Cc * energy_spread / voltage
    nom_chrom = xp.pi * Lambda * q ** 2 * H
    denom_chrom = 4 * xp.sqrt(xp.log(2))
    Kc = xp.exp(- (nom_chrom / denom_chrom) ** 2)
    # spatial envelope
    nums = (xp.pi * Cs * Lambda ** 2 * q ** 3 - xp.pi * defocus * q) ** 2 * illumination_aperture ** 2
    Ks = xp.exp(- nums / xp.log(2))
    # full envelope
    K = Kc * Ks

    # aperture function
    A = xp.ones(image_shape)
    qmax = 2 * xp.pi * objective_diameter / (Lambda * focus_length)
    A[q > qmax] = 0
    A = gaussian_filter(A, sigma=3)

    CTF = complex_CTF * K * A

    # print(complex_CTF.shape, K.shape, A.shape)
    # visual inspection of CTF
    if display_CTF:
        r1, m1 = radial_average(CTF.real)
        r2, m2 = radial_average(CTF.imag)
        fig, ax = plt.subplots()

        ax.plot(r1, m1, label='real')
        ax.plot(r2, m2, label='imaginary')
        ax.legend()
        plt.savefig('ctf.png')

    return CTF


def transmission_function(sliced_potential, voltage, dz):
    # wavelength
    Lambda = wavelength_eV2m(voltage)
    # relative mass
    relative_mass = phys.constants["me"] + phys.constants["el"] * voltage / (phys.constants["c"] ** 2)
    # sigma_transfer
    sigma_transfer = 2 * xp.pi * relative_mass * phys.constants["el"] * Lambda / (phys.constants["h"] ** 2)

    return xp.exp(1j * sigma_transfer * sliced_potential * dz)


def fresnel_propagator(imageSize, pixel_size, voltage, dz):
    Lambda = wavelength_eV2m(voltage)

    xwm = pixel_size * (imageSize)  # pixelsize for multislice * size sample
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
               binning=1, voltage=300E3, camera_type='K2SUMMIT', camera_folder='', correction=False):
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
        # Absorption factor by sample
        if correction:
            absorption_factor = xp.exp( - (ice_thickness / phys.mean_free_path[voltage]) / xp.cos(angle * xp.pi / 180) )
        else:
            absorption_factor = 1
        print(f'number of electrons per pixel for tilt angle {angle:.2f} degrees (before binning): {dose_per_pixel * absorption_factor:.2f}')

        # Apply poissonian noise
        poisson_mean = projection * ( dose_per_pixel * absorption_factor )
        projection_poisson = xp.zeros(projection.shape)
        for _ in range(binning**2):
            projection_poisson += ( xp.random.poisson(lam=poisson_mean) / binning**2 )
            # imshow(projection_poisson)
            # show()
        # projection_poisson = xp.random.poisson(lam=poisson_mean)

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
                         defocus=2E-6, sigma_damage=0.0, camera_type='K2SUMMIT',
                         camera_folder='', random_gradient=False, random_rotation=False, solvent_potential=4.5301,
                         solvent_factor=1.0, absorption_contrast=False):
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
    @param sigma_damage:
    @type sigma_damage:
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

    # imageSize = grandcell.shape[0]//2
    box_size = grandcell.shape[0]
    box_height = grandcell.shape[2]

    del grandcell # delete volume here?

    zheight = box_height * pixel_size  # thickness of rotation volume in nm
    defocus -= (zheight / 2) # center defocus value at tilt angle

    ileft = (box_size - image_size)//2
    iright = -int(xp.ceil((box_size - image_size)/2))

    if absorption_contrast:
        solvent_amplitude = potential_amplitude(0.93, 18, voltage) * solvent_factor
        print(f'solvent absorption = {solvent_amplitude:.3f}')

    # create a 3d array for the exit wave of each image
    noisefree_projections = xp.zeros((image_size, image_size, len(angles)))
    # get the contrast transfer function
    complex_CTF, alt_ctf = create_complex_CTF_ext((image_size, image_size), pixel_size, defocus, voltage=voltage,
                                         Cs=spherical_aberration, display_CTF=False)

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

        filename = f'{save_path}/rotations/rotated_volume_{int(angle)}'
        rotated_volume = pytom.tompy.io.read_mrc(f'{filename}.mrc')  #* 1E1
        if absorption_contrast: rotated_volume_imag = pytom.tompy.io.read_mrc(f'{filename}_imag.mrc')

        if not multislice:
            print('Simulating projection (without ms) from tilt angle ', angle)
            # remove .get() if fully cupy
            projected_tilt_image = rotated_volume[ileft:iright, ileft:iright, :].sum(axis=2)#.get()
            projected_tilt_image = xp.abs(xp.fft.ifftn(xp.fft.fftn(projected_tilt_image)*complex_CTF))**2
            #TODO is abs()**2 the way to go for exit wave field of simple projection?

        else:
            print('Simulating projection (with ms) from tilt angle ', angle)

            # Rotation artifacts are removed in rotate_model() !

            add_ice_layer=True
            if add_ice_layer:
                # add ice layer with some smoothing to prevent sharp edges
                ice = create_ice_layer(rotated_volume.shape, angle, ice_thickness//pixel_size,
                                       value=1.0, sigma=0.0)
                rotated_volume += (ice * (solvent_potential * solvent_factor))
                if absorption_contrast:
                    rotated_volume_imag += (ice * solvent_amplitude)
                del ice # remove to free up memory
            else: # might be good to also have an option for not adding the ice layer, purely simulation of rotation and projection
                if absorption_contrast:
                    rotated_volume_imag += solvent_amplitude

            if absorption_contrast:
                rotated_volume = rotated_volume.astype(complex)
                rotated_volume = rotated_volume.real + 1j * rotated_volume_imag

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
                if ileft==0 and iright==0:
                    projected_potent_ms[:, :, ii] = rotated_volume[:,:,
                                                    ii * px_per_slice: (ii + 1) * px_per_slice].mean(axis=2)  # .get()
                else:
                    projected_potent_ms[:, :, ii] = rotated_volume[ileft:iright, ileft:iright,
                                                    ii * px_per_slice: (ii + 1) * px_per_slice].mean(axis=2)  # .get()

            # calculate the transmission function for each slice
            psi_t = transmission_function(projected_potent_ms, voltage, msdz)

            # calculate the fresnel propagator (identical for same dz)
            propagator = fresnel_propagator(image_size, pixel_size, voltage, msdz)

            # Wave propagation with MULTISLICE method
            psi_multislice = xp.zeros((image_size,image_size), dtype=complex) + 1 # should be complex datatype

            for ii in range(n_slices-min(1,num_px_last_slice)):
                waveField = xp.fft.fftn( xp.fft.ifftshift(psi_multislice) * xp.fft.ifftshift(psi_t[:, :, ii]) )
                psi_multislice = xp.fft.fftshift( xp.fft.ifftn((waveField * xp.fft.ifftshift(propagator) )) )

            # Calculate propagation through last slice in case the last slice contains a different number of pixels
            if num_px_last_slice:
                msdz_end = num_px_last_slice * pixel_size
                psi_t[:, :, -1] = transmission_function(projected_potent_ms[:, :, -1], voltage, msdz_end)
                propagator_end = fresnel_propagator(image_size, pixel_size, voltage, msdz_end)
                waveField = xp.fft.fftn( xp.fft.ifftshift(psi_multislice) * xp.fft.ifftshift(psi_t[:, :, -1]) )
                psi_multislice = xp.fft.fftshift( xp.fft.ifftn( waveField * xp.fft.ifftshift(propagator_end) ) )

            # Multiple by CTF for microscope effects on electron wave
            wave_CTF = xp.fft.ifftshift(complex_CTF) * xp.fft.fftn(xp.fft.ifftshift(psi_multislice) )
            # Intensity in image plane is obtained by taking the absolute square of the wave function
            projected_tilt_image = xp.abs(xp.fft.fftshift(xp.fft.ifftn(wave_CTF))) ** 2

        # imshow(projected_tilt_image)
        # show()
        #
        # test = xp.log(xp.abs(xp.fft.fftshift(xp.fft.fftn(projected_tilt_image))))
        #
        # imshow(test)
        # show()
        #
        # _ = radial_average(test)

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
    if absorption_contrast:
        projections = microscope(noisefree_projections, angles, ice_thickness=ice_thickness, dose=dose,
                                 pixel_size=pixel_size, binning=binning, voltage=voltage, camera_type=camera_type,
                                 camera_folder=camera_folder, correction=False)
    else:
        projections = microscope(noisefree_projections, angles, ice_thickness=ice_thickness, dose=dose,
                                 pixel_size=pixel_size, binning=binning, voltage=voltage, camera_type=camera_type,
                                 camera_folder=camera_folder, correction=True)
    pytom.tompy.io.write(f'{save_path}/projections.mrc', projections)

    # Create a folder for individual noisy projections (needed for reconstruction algorithm)
    if not os.path.exists(f'{save_path}/projections/'):
        os.mkdir(f'{save_path}/projections/')
    for i in range(len(angles)):
        # These files are needed for the reconstruction
        pytom.tompy.io.write(f'{save_path}/projections/synthetic_{i+1}.mrc',
                             projections[:, :, i].squeeze())

    # len(angles) is the number of files that we have
    alignment = xp.zeros(len(angles), dtype=dar)
    # correct angle by -1 to get the right reconstruction
    alignment['TiltAngle'] = -1.0 * xp.array(angles)
    alignment['Magnification'] = xp.repeat(1.0, len(angles))
    if random_rotation:
        alignment['InPlaneRotation'] = random_angles

    for i in range(len(angles)):
        alignment['FileName'][i] = f'{save_path}/projections/synthetic_{i+1}.mrc'

    # Then the alignment file, because the noisyProjections folder needs to have been created
    alignment_file = f'{save_path}/projections/alignment_simulated.txt'
    savestar(alignment_file, alignment, fmt=fmtAlignmentResults, header=HEADER_ALIGNMENT_RESULTS)
    return


def microscope_single_projection(noisefree_projection, dqe, mtf, dose, pixel_size, binning=1):
    """
    Inspired by InSilicoTEM (Vulovic et al., 2013)
    @author: Marten Chaillet
    """
    from scipy.ndimage import shift

    ntf = xp.sqrt(mtf ** 2 / dqe) # square root because dqe = mtf^2 / ntf^2
    ntf = xp.maximum(ntf, 1E-7) # ensure we do not divide by zero
    mtf_shift = xp.fft.ifftshift(mtf)
    ntf_shift = xp.fft.ifftshift(ntf)

    # NUMBER OF ELECTRONS PER PIXEL
    dose_per_pixel = dose * (pixel_size*1E10)**2 / binning**2 # from square A to square nm (10A pixels)
    print(f'Number of electrons per pixel (before binning and sample absorption): {dose_per_pixel}')

    # Fourier transform and multiply with sqrt(dqe) = mtf/ntf
    projection_fourier = xp.fft.fftn(xp.fft.ifftshift(noisefree_projection)) * mtf_shift / ntf_shift
    # projection_fourier = projection_fourier
    # Convert back to real space
    projection = xp.real(xp.fft.fftshift(xp.fft.ifftn(projection_fourier)))
    projection[projection<0] = 0
    # Draw from poisson distribution and scale by camera's conversion factor
    # conversion_factor = 100  # in ADU/e- , this is an arbitrary unit. Value taken from Vulovic et al., 2010

    # Apply poissonian noise
    # poisson_mean = projection * dose_per_pixel
    sigma_motion = 1 / (pixel_size*1E9)
    projection_poisson = xp.zeros(projection.shape)
    for _ in range(binning**2):
        translation = (xp.random.normal(0, sigma_motion), xp.random.normal(0, sigma_motion))
        # projection_poisson += (
        #         xp.random.poisson(lam=shift(projection, translation, mode='nearest') * dose_per_pixel) / binning ** 2)
        poisson_intermediate = xp.random.poisson(lam=shift(projection, translation, mode='nearest') * dose_per_pixel)
        projection_poisson += xp.real(xp.fft.fftshift(xp.fft.ifftn(xp.fft.fftn(xp.fft.ifftshift(poisson_intermediate))
                                                                   * ntf_shift))) / binning**2
        # imshow(projection_poisson)
        # show()

    # Image values are now in ADU
    # Apply the camera's noise transfer function to the noisy image
    # projection_fourier = xp.fft.fftn(xp.fft.ifftshift(projection_poisson)) * ntf_shift

    # Add additional noise from digitization process, less relevant for modern day cameras.
    # readout noise standard deviation can be 7 ADUs, from Vulovic et al., 2010
    # sigma_readout = 7
    # readsim = xp.random.normal(0, sigma_readout, projection.shape) # readout noise has a gaussian distribution
    # darksim = 0     # dark current noise has a poisson distribution, usually an order of magnitude smaller than
    #                 # readout noise and can hence be neglected

    # Add readout noise and dark noise in real space
    # return xp.real(xp.fft.fftshift(xp.fft.ifftn(projection_fourier))) # + readsim + darksim
    return projection_poisson


def parallel_project(grandcell, folder, frame, image_size, pixel_size, msdz, n_slices, ctf, dose, dqe, mtf, voltage,
                     binning=1, translation=(.0,.0,.0), rotation=(.0,.0,.0), solvent_potential=.0,
                     solvent_absorption=.0, sigma_damage=.0, ice_voxels=None):
    """
    Only use multislice for this.
    @param model:
    @param frame:
    @param translation:
    @param rotation:
    @return:
    """
    from voltools import transform

    print('Transforming sample for tilt/frame ', frame)

    sample = grandcell.copy()

    max_tilt_radians = abs(rotation[1]) * xp.pi / 180
    max_tilt_radians_opp = (90 - abs(rotation[1])) * xp.pi / 180
    rotation_height = int(xp.ceil(xp.sin(max_tilt_radians) * image_size +
                                           xp.sin(max_tilt_radians_opp) * ice_voxels))
    print(f'Reduced rotation height for relevant specimens: {rotation_height}')
    if rotation_height % 2: rotation_height += 1
    diff = sample.shape[2]-rotation_height
    i = diff // 2

    # model parameter is a complex volume or real volume depending on the addition of absorption contrast
    # first transform the volume
    if sample.dtype == 'complex64':
        transform(sample[...,i:-i].real, translation=translation, rotation=rotation, rotation_order='sxyz',
                                interpolation='filt_bspline', device='cpu', output=sample[...,i:-i].real)
        transform(sample[...,i:-i].imag, translation=translation, rotation=rotation, rotation_order='sxyz',
                                interpolation='filt_bspline', device='cpu', output=sample[...,i:-i].imag)
        # remove rotation artifacts
        threshold = 0.001
        sample.real[sample.real < threshold] = 0
        sample.imag[sample.imag < threshold] = 0
    elif sample.dtype == 'float32':
        transform(sample[...,i:-i], translation=translation, rotation=rotation, rotation_order='sxyz',
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
    if ice_voxels is None: # this does not seem like the best condition 'or not sample.dtype=='complex64' '
        ice_layer = 1
    else:
        ice_layer = create_ice_layer(sample.shape, rotation[1], ice_voxels, value=1.0, sigma=0.0)
        # ice layer datatype at this point should be np.float32
        # print(f'data type of ice: {ice_layer.dtype}')
    # apply structural deterioration due to beam damage via random noise with increment based on frame number
    # incremental damage based on frame number
    # if not (sigma_damage == 0.0):
    #     # this does not add noise to imaginary part, but maybe not needed
    #     sample += ( ice_layer * xp.random.normal(0, sigma_damage, sample.shape) )

    if sample.dtype == 'complex64':
        sample.real += ( ice_layer * solvent_potential ) # Maybe we can remove this? Are we only interested in absorption of ice layer?
        # ice_layer *= solvent_absorption
        sample.imag += (ice_layer * solvent_absorption )
    else:
        sample += ( ice_layer * solvent_potential )

    del ice_layer # free memory

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

    # Blur the image at this point for some motion blur?
    # from scipy.ndimage import gaussian_filter
    # gaussian_filter(psi_multislice.real, sigma=(20/(2 * (pixel_size * 1E10)))/2.35, output=psi_multislice.real)

    # Multiple by CTF for microscope effects on electron wave
    wave_ctf = xp.fft.ifftshift(ctf) * xp.fft.fftn(xp.fft.ifftshift(psi_multislice) )
    # Intensity in image plane is obtained by taking the absolute square of the wave function
    noisefree_projection = xp.abs(xp.fft.fftshift(xp.fft.ifftn(wave_ctf))) ** 2

    # wave_ctf = xp.fft.ifftshift(ctf.imag) * xp.fft.fftn(xp.fft.ifftshift(xp.abs(psi_multislice) ** 2))
    # noisefree_projection = xp.abs(xp.fft.fftshift(xp.fft.ifftn(wave_ctf)))

    # noisefree_projection = xp.abs(psi_multislice)**2

    # # Blur the image at this point for some motion blur?
    # from scipy.ndimage import gaussian_filter
    # gaussian_filter(noisefree_projection, sigma=2, output=noisefree_projection)

    # DEBUGGING
    # test = xp.log(xp.abs(xp.fft.fftshift(xp.fft.fftn(noisefree_projection))))
    # r1, m1 = radial_average(test)
    # fig, ax = plt.subplots()
    # ax.plot(r1, m1)
    # ax.legend()
    # plt.savefig(f'{folder}/radial.png')

    # Apply the microscope function
    projection = microscope_single_projection(noisefree_projection, dqe, mtf, dose, pixel_size, binning=binning)
    # Write the projection to the projection folder
    pytom.tompy.io.write(f'{folder}/synthetic_{frame+1}.mrc', projection)
    # Return noisefree and projection as tuple for writing as mrc stack in higher function
    return (noisefree_projection, projection)


def generate_tilt_series_cpu(output_folder, model_ID, angles, nodes=1, image_size=None, rotation_box_height=None,
                             pixel_size=1E-9, binning=1, dose=80, voltage=300E3, spherical_aberration=2.7E-3,
                             chromatic_aberration=2.7E-3, energy_spread=0.7, illumination_aperture=0.030E-3,
                             objective_diameter=100E-6, focus_length=4.7E-3, astigmatism=0.0E-9, astigmatism_angle=0,
                             msdz=1E-9, defocus=2E-6, sigma_damage=0.0, camera_type='K2SUMMIT', camera_folder='',
                             solvent_potential=4.5301, absorption_contrast=False):
    from pytom.basic.datatypes import DATATYPE_ALIGNMENT_RESULTS as dar
    from pytom.basic.datatypes import fmtAlignmentResults, HEADER_ALIGNMENT_RESULTS
    from pytom.gui.guiFunctions import savestar
    import detector
    from joblib import Parallel, delayed
    # NOTE; Parameter defocus specifies the defocus at the bottom of the model!
    save_path = f'{output_folder}/model_{model_ID}'
    # create folder for individual projections
    projection_folder = f'{save_path}/projections'
    if not os.path.exists(projection_folder):
        os.mkdir(projection_folder)

    # grab model
    grandcell = pytom.tompy.io.read_mrc(f'{save_path}/grandmodel_original.mrc')
    if absorption_contrast:
        xp_type = xp.complex64
        grandcell = grandcell.astype(xp_type)
        grandcell.imag = pytom.tompy.io.read_mrc(f'{save_path}/grandmodel_original_imag.mrc')
        # calculate the absorption for amorphous ice at the specified voltage
        solvent_amplitude = potential_amplitude(0.93, 18, voltage)
        print(f'solvent absorption = {solvent_amplitude:.3f}')
    else:
        xp_type = xp.float32
        grandcell = grandcell.astype(xp_type)
        solvent_amplitude = 0.0

    box_size = grandcell.shape[0]
    box_height = grandcell.shape[2]

    # For sample set the arrays specifically to np.complex64 datatype to save memory space
    if rotation_box_height is not None:
        # Use the specified rotation box height
        rotation_volume = xp.zeros((box_size, box_size, rotation_box_height), dtype=xp_type) # TODO datatype should not be complex if working with only real potential
        offset = (rotation_box_height - box_height) // 2

        if (rotation_box_height - box_height) % 2:
            rotation_volume[:, :, offset + 1:-offset] = grandcell[:, :, :]
        else:
            rotation_volume[:, :, offset:-offset] = grandcell[:, :, :]
    else:
        # Calculate maximal box height to contain all the information for the rotation
        max_tilt = max([abs(a) for a in angles])
        max_tilt_radians = max_tilt * xp.pi / 180
        rotation_box_height = int(xp.ceil(xp.tan(max_tilt_radians) * image_size +
                                          thickness_voxels / xp.cos(max_tilt_radians)))
        if rotation_box_height % 2: rotation_box_height += 1
        # Place grandcell in rotation volume
        rotation_volume = xp.zeros((box_size, box_size, rotation_box_height), dtype=xp_type)
        offset = (rotation_box_height - box_height) // 2

        if (rotation_box_height - box_height) % 2:
            rotation_volume[:, :, offset + 1:-offset] = grandcell[:, :, :]
        else:
            rotation_volume[:, :, offset:-offset] = grandcell[:, :, :]
    del grandcell

    if image_size is None:
        image_size = box_size

    # confirm image_size is valid
    assert image_size <= box_size, 'Specified projection image size is invalid as it is larger than the model dimension.'

    # adjust defocus
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
    sigma_motion = 0.0 # nm
    translations = []
    for i in angles:
        x = xp.random.normal(0, sigma_motion) / (pixel_size*1E9)
        y = xp.random.normal(0, sigma_motion) / (pixel_size*1E9)
        # print((x,y,0.0))
        translations.append((x,y,0.0))

    # defocus_series = [xp.random.normal(defocus, 0.2E-6) for a in angles]
    ctf_series = []
    for i in angles:
        dz = xp.random.normal(defocus, 0.2E-6)
        ctf = create_complex_CTF_ext((image_size, image_size), pixel_size, dz, voltage=voltage,
                                     Cs=spherical_aberration, Cc=chromatic_aberration, energy_spread=energy_spread,
                                     illumination_aperture=illumination_aperture, objective_diameter=objective_diameter,
                                     focus_length=focus_length, astigmatism=astigmatism,
                                     astigmatism_angle=astigmatism_angle, display_CTF=False)
        ctf_series.append(ctf)

    dqe = detector.create_detector_response(camera_type, 'DQE', xp.zeros((image_size, image_size)), voltage=voltage,
                                            folder=camera_folder)
    mtf = detector.create_detector_response(camera_type, 'MTF', xp.zeros((image_size, image_size)), voltage=voltage,
                                            folder=camera_folder)

    # joblib automatically memory maps a numpy array to child processes
    print(f'Projecting the model with {nodes} processes')

    verbosity = 55  # set to 55 for debugging, 11 to see progress, 0 to turn off output
    results = Parallel(n_jobs=nodes, verbose=verbosity, prefer="threads") \
        (delayed(parallel_project)(rotation_volume, projection_folder, i, image_size, pixel_size, msdz, n_slices, ctf,
                                   dose_per_tilt, dqe, mtf, voltage, binning=binning, translation=translation,
                                   rotation=(.0, angle, .0), solvent_potential=solvent_potential,
                                   solvent_absorption=solvent_amplitude, sigma_damage=sigma_damage,
                                   ice_voxels=box_height)
         for i, (angle, translation, ctf) in enumerate(zip(angles, translations, ctf_series)))

    sys.stdout.flush()
    if results.count(None) == 0:
        print('All projection processes finished successfully')
    else:
        print(f'{results.count(None)} rotation processes did not finish successfully')

    # write (noisefree) projections as mrc stacks
    pytom.tompy.io.write(f'{save_path}/noisefree_projections.mrc', xp.stack([n for (n, p) in results], axis=2))
    pytom.tompy.io.write(f'{save_path}/projections.mrc', xp.stack([p for (n, p) in results], axis=2))

    # store alignment information
    # len(angles) is the number of files that we have
    alignment = xp.zeros(len(angles), dtype=dar)
    # correct angle by -1 to get the right reconstruction
    alignment['TiltAngle'] = -1.0 * xp.array(angles)
    alignment['Magnification'] = xp.repeat(1.0, len(angles))
    for i in range(len(angles)):
        alignment['FileName'][i] = f'{save_path}/projections/synthetic_{i+1}.mrc'
    # write the alignment file
    alignment_file = f'{save_path}/projections/alignment_simulated.txt'
    savestar(alignment_file, alignment, fmt=fmtAlignmentResults, header=HEADER_ALIGNMENT_RESULTS)
    return


def generate_frame_series_cpu(output_folder, model_ID, n_frames=20, nodes=1, image_size=None, pixel_size=1E-9,
                              binning=1, dose=80, voltage=300E3, spherical_aberration=2.7E-3,
                              chromatic_aberration=2.7E-3, energy_spread=0.7, illumination_aperture=0.030E-3,
                              objective_diameter=100E-6, focus_length=4.7E-3, astigmatism=0.0E-9, astigmatism_angle=0,
                              msdz=1E-9, defocus=2E-6, sigma_damage=0.0, camera_type='K2SUMMIT', camera_folder='',
                              solvent_potential=4.5301, absorption_contrast=False):
    """
    Creating a frame series for the initial grand model by applying stage drift (translation) and some rotations for
    each frame, and the calculating the sample projection in the microscope. Additionally apply increasing random noise
    to the images to imitate particle degradation through beam interaction. I cannot do this incremental as I want
    to parallelize the rotation+projection, but the overall effect on signal should be identical.

    @param output_folder:
    @param model_ID:
    @param n_frames:
    @param image_size:
    @param pixel_size:
    @param binning:
    @param dose:
    @param voltage:
    @param spherical_aberration:
    @param multislice:
    @param msdz:
    @param defocus:
    @param sigma_decay_CTF:
    @param camera_type:
    @param camera_folder:
    @param random_gradient:
    @param solvent_potential:
    @param solvent_factor:
    @param absorption_contrast:
    @return:
    """
    from pytom.basic.datatypes import DATATYPE_ALIGNMENT_RESULTS as dar
    from pytom.basic.datatypes import fmtAlignmentResults, HEADER_ALIGNMENT_RESULTS
    from pytom.gui.guiFunctions import savestar
    import detector
    from joblib import Parallel, delayed

    save_path = f'{output_folder}/model_{model_ID}'

    # NOTE; Parameter defocus specifies the defocus at the bottom of the model!
    # create folder for individual projections
    projection_folder = f'{save_path}/projections'
    if not os.path.exists(projection_folder):
        os.mkdir(projection_folder)

    # grab model
    grandcell = pytom.tompy.io.read_mrc(f'{save_path}/grandmodel_original.mrc')
    if absorption_contrast:
        grandcell = grandcell.astype(xp.complex64)
        grandcell.imag = pytom.tompy.io.read_mrc(f'{save_path}/grandmodel_original_imag.mrc')
        # calculate the absorption for amorphous ice at the specified voltage
        solvent_amplitude = potential_amplitude(0.93, 18, voltage)
        print(f'solvent absorption = {solvent_amplitude:.3f}')
    else:
        grandcell = grandcell.astype(xp.float32)
        solvent_amplitude = 0.0

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
    global_motion = xp.random.normal(10, 3) # normal around mean 10 A and std 3A
    average_motion_per_frame = global_motion / n_frames
    global_angle = xp.random.uniform(0,360) # random angle from uniform
    translations, cumulative_translations, translations_voxel = [], [], []
    x, y = 0, 0
    for i in range(n_frames):
        # randomly vary the motion per frame and angle
        motion_i = xp.random.normal(average_motion_per_frame, average_motion_per_frame/2)
        angle_i = xp.random.normal(global_angle, 20)
        # decompose motion into x and y translation
        y_i = xp.sin(angle_i * xp.pi / 180) * motion_i
        x_i = xp.cos(angle_i * xp.pi / 180) * motion_i
        # TODO decide if translations can be removed, only seem to need cumulative_translations
        translations.append((x_i, y_i, 0)) # append the translation for the frame as a tuple
        y += y_i
        x += x_i
        cumulative_translations.append((x,y, 0)) # translation for z coordinate as we are referring to volumes
        translations_voxel.append((x*1E-10 / pixel_size, y*1E-10 / pixel_size, 0))

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

    # get the contrast transfer function
    ctf = create_complex_CTF_ext((image_size, image_size), pixel_size, defocus, voltage=voltage,
                                   Cs=spherical_aberration, Cc=chromatic_aberration, energy_spread=energy_spread,
                                   illumination_aperture=illumination_aperture, objective_diameter=objective_diameter,
                                   focus_length=focus_length, astigmatism=astigmatism,
                                   astigmatism_angle=astigmatism_angle, display_CTF=False)

    dqe = detector.create_detector_response(camera_type, 'DQE', xp.zeros((image_size, image_size)), voltage=voltage,
                                            folder=camera_folder)
    mtf = detector.create_detector_response(camera_type, 'MTF', xp.zeros((image_size, image_size)), voltage=voltage,
                                            folder=camera_folder)

    # joblib automatically memory maps a numpy array to child processes
    print(f'Projecting the model with {nodes} processes')

    verbosity = 55  # set to 55 for debugging, 11 to see progress, 0 to turn off output
    results = Parallel(n_jobs=nodes, verbose=verbosity, prefer="threads") \
        (delayed(parallel_project)(grandcell, projection_folder, frame, image_size, pixel_size, msdz, n_slices, ctf,
                                   dose_per_frame, dqe, mtf, voltage, binning=binning, translation=shift, rotation=(.0,.0,.0),
                                   solvent_potential=solvent_potential, solvent_absorption=solvent_amplitude,
                                   sigma_damage=sigma_damage, ice_voxels=None)
         for frame, shift in enumerate(translations_voxel))

    sys.stdout.flush()
    if results.count(None) == 0:
        print('All projection processes finished successfully')
    else:
        print(f'{results.count(None)} rotation processes did not finish successfully')

    # write (noisefree) projections as mrc stacks
    pytom.tompy.io.write(f'{save_path}/noisefree_projections.mrc', xp.stack([n for (n,p) in results], axis=2))
    pytom.tompy.io.write(f'{save_path}/projections.mrc', xp.stack([p for (n, p) in results], axis=2))

    # Store translations as reference for model
    # len(angles) is the number of files that we have
    alignment = xp.zeros(n_frames, dtype=dar)
    alignment['AlignmentTransX'] = xp.array([x for (x,y,z) in cumulative_translations])
    alignment['AlignmentTransY'] = xp.array([y for (x,y,z) in cumulative_translations])
    alignment['Magnification'] = xp.repeat(1.0, n_frames)
    for i in range(n_frames):
        alignment['FileName'][i] = f'{projection_folder}/synthetic_{i+1}.mrc'

    # Write the alignment file as a text file
    alignment_file = f'{projection_folder}/alignment_simulated.txt'
    savestar(alignment_file, alignment, fmt=fmtAlignmentResults, header=HEADER_ALIGNMENT_RESULTS)
    return


def reconstruct_tomogram(size_reconstruction, output_folder, model_ID, weighting=-1, crop=False, binning=1):
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
    from pytom.basic.datatypes import DATATYPE_ALIGNMENT_RESULTS as dar
    from pytom.gui.guiFunctions import loadstar

    save_path = f'{output_folder}/model_{model_ID}'
    alignment_file = f'{save_path}/projections/alignment_simulated.txt'
    alignment_data = loadstar(alignment_file, dtype=dar)
    vol_size = [size_reconstruction, ] * 3

    projections = ProjectionList()
    # IMPORTANT: angles *-1 to get the right reconstrunction relative to the orignal model!
    # for d in alignment_data:
    #     # In the alignment file the angles were already multiplied by -1
    #     p = Projection(d[5], tiltAngle= d[2])
    #     # Old method: here we multiply by -1 to get right angles
    #     # p = Projection(f'{save_path}/projections/synthetic_{i+1}.mrc', tiltAngle=-1 * angles[i - 1])
    #     projections.append(p)

    outputname = f'{save_path}/reconstruction.em'
    # IF EM alignment file provided, filters applied and reconstruction will be identical.
    vol = projections.reconstructVolume(dims=vol_size, reconstructionPosition=[0,0,0], binning=binning,
                                        applyWeighting=weighting, alignResultFile=alignment_file)
    vol.write(outputname)
    os.system(f'em2mrc.py -f {outputname} -t {os.path.dirname(outputname)}')
    os.system(f'rm {outputname}')

    # Flag for checking cropping?!
    if crop:
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
            # cell = pytom.tompy.io.read_mrc(f'{save_path}/grandmodel_noisefree_original.mrc')
            # pytom.tompy.io.write(f'{save_path}/grandmodel_noisefree.mrc', cell[ileft:iright, ileft:iright, :])
            cell = pytom.tompy.io.read_mrc(f'{save_path}/class_bbox_original.mrc')
            pytom.tompy.io.write(f'{save_path}/class_bbox.mrc', cell[ileft:iright, ileft:iright, :])
            cell = pytom.tompy.io.read_mrc(f'{save_path}/class_mask_original.mrc')
            pytom.tompy.io.write(f'{save_path}/class_mask.mrc', cell[ileft:iright, ileft:iright, :])
            cell = pytom.tompy.io.read_mrc(f'{save_path}/occupancy_bbox_original.mrc')
            pytom.tompy.io.write(f'{save_path}/occupancy_bbox.mrc', cell[ileft:iright, ileft:iright, :])
            cell = pytom.tompy.io.read_mrc(f'{save_path}/occupancy_mask_original.mrc')
            pytom.tompy.io.write(f'{save_path}/occupancy_mask.mrc', cell[ileft:iright, ileft:iright, :])
        else:
            pytom.tompy.io.write(f'{save_path}/grandmodel.mrc', cell[:, :, :])
            # os.system(f'cp {save_path}/grandmodel_noisefree_original.mrc {save_path}/grandmodel_noisefree.mrc')
            os.system(f'cp {save_path}/class_bbox_original.mrc {save_path}/class_bbox.mrc')
            os.system(f'cp {save_path}/class_mask_original.mrc {save_path}/class_mask.mrc')
            os.system(f'cp {save_path}/occupancy_bbox_original.mrc {save_path}/occupancy_bbox.mrc')
            os.system(f'cp {save_path}/occupancy_mask_original.mrc {save_path}/occupancy_mask.mrc')

        del cell
    return


if __name__ == '__main__':
    from pytom.gui.guiFunctions import loadstar, datatype
    from ast import literal_eval
    # Use tracemalloc to record the peak memory usage of the program
    tracemalloc.start()

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
        output_folder           = config['General']['OutputFolder']
        simulator_mode          = config['General']['Mode']
        device                  = config['General']['Device']
        nodes                   = config['General'].getint('Nodes')
        model_ID                = config['General'].getint('ModelID')
        seed                    = config['General'].getint('Seed')
        pixel_size              = config['General'].getfloat('PixelSize') * 1E-9 # pixel_size in nm
        binning                 = config['General'].getint('Binning')
        solvent_potential       = config['General'].getfloat('SolventConstant')
        absorption_contrast     = config['General'].getboolean('AbsorptionContrast')
        voltage                 = config['General'].getfloat('Voltage') * 1E3  # voltage in keV
        # voltage and pixelsize are needed for model generation and projection, thus general parameters
        # ensure simulator mode and device are valid options
        if (simulator_mode in ['TiltSeries', 'FrameSeries']) or (device in ['CPU', 'GPU']):
            print(f'Generating model {model_ID} in folder {output_folder}')
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
            if thickness_voxels % 2 == 1: thickness_voxels += 1
            # gold markers
            number_of_markers   = draw_range(literal_eval(config['GenerateModel']['NumberOfMarkers']), int,
                                           'NumberOfMarkers')
            # parse range of number of particles
            number_of_particles = draw_range(literal_eval(config['GenerateModel']['NumberOfParticles']), int,
                                             'NumberOfParticles')
            number_of_membranes = draw_range(literal_eval(config['GenerateModel']['NumberOfMembranes']), int,
                                             'NumberOfMembranes')
        except Exception as e:
            print(e)
            raise Exception('Missing generate model parameters in config file.')

    if 'Microscope' in config.sections():
        try:
            camera                  = config['Microscope']['Camera']
            camera_folder           = config['Microscope']['CameraFolder']
            defocus                 = draw_range(literal_eval(config['Microscope']['Defocus']), float, 'Defocus') * 1E-6
            electron_dose           = draw_range(literal_eval(config['Microscope']['ElectronDose']), float, 'ElectronDose')
            spherical_aberration    = config['Microscope'].getfloat('SphericalAberration') * 1E-3
            chromatic_aberration    = config['Microscope'].getfloat('ChromaticAberration') * 1E-3
            energy_spread           = config['Microscope'].getfloat('EnergySpread')
            illumination_aperture   = config['Microscope'].getfloat('IlluminationAperture') * 1E-3
            objective_diameter      = config['Microscope'].getfloat('ObjectiveDiameter') * 1E-6
            focus_length            = config['Microscope'].getfloat('FocalDistance') * 1E-3
            astigmatism             = config['Microscope'].getfloat('Astigmatism') * 1E-9
            astigmatism_angle       = config['Microscope'].getfloat('AstigmatismAngle')
        except Exception as e:
            print(e)
            raise Exception('Missing microscope parameters in config file.')

    if simulator_mode in config.sections():
        try:
            # first read common parameters between tilt and frame series
            image_size      = config[simulator_mode].getint('ImageSize')
            msdz            = config[simulator_mode].getfloat('MultisliceStep') * 1E-9
            beam_damage     = config[simulator_mode].getfloat('BeamDamage')
            # mode specific parameters
            if simulator_mode == 'TiltSeries':
                metadata            = loadstar(config['TiltSeries']['MetaFile'], dtype=datatype)
                angles              = metadata['TiltAngle'] # in degrees
            elif simulator_mode == 'FrameSeries':
                number_of_frames    = config['FrameSeries'].getint('NumberOfFrames')
        except Exception as e:
            print(e)
            raise Exception(f'Missing {simulator_mode} parameters in config file.')

    if 'TomogramReconstruction' in config.sections():
        try:
            weighting           = config['TomogramReconstruction'].getint('Weighting')
            reconstruction_size = config['TomogramReconstruction'].getint('ReconstructionSize')
            reconstruction_bin  = config['TomogramReconstruction'].getint('Binning')
        except Exception as e:
            print(e)
            raise Exception('Missing tomogram reconstruction parameters in config file.')

    # Create directories and logger
    if not os.path.exists(os.path.join(output_folder, f'model_{model_ID}')):
        os.mkdir(os.path.join(output_folder, f'model_{model_ID}'))

    logging.basicConfig(filename='{}/model_{}/simulator-{date:%Y-%m-%d_%H:%M:%S}.log'.format(output_folder, model_ID,
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

    # Generate a grand model
    if 'GenerateModel' in config.sections():
        # set seed for random number generation
        xp.random.seed(seed)
        random.seed(seed)

        print('\n- Generating grand model')
        generate_model(particle_folder, output_folder, model_ID, listpdbs, listmembranes,
                       pixel_size           =pixel_size,
                       binning              =binning,
                       size                 =size,
                       thickness            =thickness_voxels,
                       placement_size       =placement_size,
                       solvent_potential    =solvent_potential,
                       numberOfParticles    =number_of_particles,
                       number_of_markers    =number_of_markers,
                       absorption_contrast  =absorption_contrast,
                       voltage              =voltage,
                       number_of_membranes  =number_of_membranes)

    if simulator_mode == 'TiltSeries':
        # set seed for random number generation
        xp.random.seed(seed)
        random.seed(seed)
        # Grab the ice thickness from the initial model in case program is only executed for projections
        print('\n- Generating projections')
        if device == 'CPU':
            generate_tilt_series_cpu(output_folder, model_ID, angles,
                                      nodes                 =nodes,
                                      image_size            =image_size,
                                      rotation_box_height   =None, # will automatically calculate fitting size if None
                                      pixel_size            =pixel_size,
                                      binning               =binning,
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
                                      sigma_damage          =beam_damage,
                                      camera_type           =camera,
                                      camera_folder         =camera_folder,
                                      solvent_potential     =solvent_potential,
                                      absorption_contrast   =absorption_contrast)
        elif device == 'GPU':
            print('This option needs to be implemented.')
            sys.exit(0)
        else:
            print('Invalid device type.')
            sys.exit(0)
    elif simulator_mode == 'FrameSeries':
        xp.random.seed(seed)
        random.seed(seed)
        print('\n- Generate frame series projections')
        if device == 'CPU':
            generate_frame_series_cpu(output_folder, model_ID,
                                      n_frames              =number_of_frames,
                                      nodes                 =nodes,
                                      image_size            =image_size,
                                      pixel_size            =pixel_size,
                                      binning               =binning,
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
                                      sigma_damage          =beam_damage,
                                      camera_type           =camera,
                                      camera_folder         =camera_folder,
                                      solvent_potential     =solvent_potential,
                                      absorption_contrast   =absorption_contrast)
        elif device == 'GPU':
            print('This option needs to be implemented.')
            sys.exit(0)
        else:
            print('Invalid device type.')
            sys.exit(0)
    else:
        print('Invalid simulation mode specified in config file.')
        sys.exit(0)

    if 'TomogramReconstruction' in config.sections():
        print('\n- Reconstructing tomogram')
        reconstruct_tomogram(reconstruction_size, output_folder, model_ID, weighting=weighting, binning=reconstruction_bin)

    current, peak = tracemalloc.get_traced_memory()
    print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
    tracemalloc.stop()

