#!/usr/bin/env pytom
import sys, os, numpy as np
from scipy.spatial.distance import cdist
from pytom.agnostic.io import read
from pytom.angles.angleFnc import matToZYZ
from pytom.basic.structures import ParticleList
from pytom.tools.script_helper import ScriptHelper2, ScriptOption2
from pytom.tools.parse_script_options import parse_script_options2
from pytom.simulation.membrane import Vector
from pytom.voltools.utils import rotation_matrix

import matplotlib

try:  # try to load this backend for the cluster
    matplotlib.use('Qt5Agg')
except:
    pass

import matplotlib.pyplot as plt


def neighbour_position_3d(tomograms, coordinates, rotations, neighbourhood=4, center_class=None,
                          plane_norm=[0, 0, 1], reference=[0, 0, 1], estimate_membrane_plane=False):
    # get set of the unique tomogram names to loop over
    unique_tomograms = np.unique(tomograms)

    # plane adjustment
    reference = Vector(reference)
    plane_norm = Vector(plane_norm)
    rm_adjust = reference.get_rotation(plane_norm)

    # list of output coordinates
    relative_coords = []
    n_center_particles = 0

    if center_class is None:
        center_class = np.array([True, ] * coordinates.shape[0])

    for tomo in unique_tomograms:
        # select tomogram from dataframe
        tomo_particles = (tomograms == tomo)
        n_part_in_tomo = tomo_particles.sum()
        tomo_coordinates = coordinates[tomo_particles]
        tomo_rotations = rotations[tomo_particles]
        tomo_centers = center_class[tomo_particles]

        if n_part_in_tomo < 2:
            # skip this tomogram because the particles has no neighbours
            continue
        else:
            n_center_particles += tomo_centers.sum()

        # calculate distance matrix
        dist_matrix = cdist(tomo_coordinates, tomo_coordinates)
        max_dist = np.max(dist_matrix)
        dist_matrix[dist_matrix == 0] = max_dist + 1
        # print(dist_matrix.min())

        for i in range(n_part_in_tomo):

            if not tomo_centers[i]:
                continue

            loop = neighbourhood if (n_part_in_tomo - 1) >= neighbourhood else n_part_in_tomo - 1

            for _ in range(loop):

                j = np.argmin(dist_matrix[i])
                distance = dist_matrix[i][j]  # ; could write out distance? or check distance?
                dist_matrix[i][j] = max_dist + 1

                if distance < 1000:  # assume distance larger than 1000 A will not be relevant
                    # get center and neighbour coordinate
                    coord_p = tomo_coordinates[i]
                    coord_n = tomo_coordinates[j]

                    # vector from center particle to neighbour
                    v = Vector(coord_n - coord_p)

                    # voltools rotation matrix
                    rm = rotation_matrix(rotation=-tomo_rotations[i], rotation_order='rzyz').T
                    v.rotate(rm[:3, :3])

                    # append to results
                    relative_coords.append(list(v.get()))

    relative_coords = np.array(relative_coords).T
    print('---> you have ', n_center_particles, ' particle centers')
    print('---> they have ', relative_coords[1].shape, ' neighbours')

    if estimate_membrane_plane:
        svd = np.linalg.svd(relative_coords[:, 0:25000] -
                            np.mean(relative_coords[:, 0:25000], axis=1, keepdims=True))
        left = svd[0]
        plane_norm = Vector(left[:, -1])
        rm_adjust = reference.get_rotation(plane_norm)

        print('---> estimated membrane plane vector ', left[:, -1])

    return np.dot(relative_coords.T, rm_adjust).T, plane_norm.get()


def density_plot(data, hist_limits, hist_voxel_size, fig_size=(5, 5), vrange=None, probability=True, plane='xy'):
    """
    @param data: (3, N) shape array, x, y, z coordinates with A units
    @param hist_limits: tuple of two elements giving min and max value for all axis
    @param hist_voxel_size: voxel size of histogram in A units
    """
    assert data.shape[0] == 3, 'data does not have x, y, z as the second axis'

    plane_to_axis = {'xy': 2,
                     'xz': 1,
                     'yz': 0}

    # calculate bins for the dimensions
    n_bins = int(abs(hist_limits[1] - hist_limits[0]) / hist_voxel_size)

    hist_3d, hist_3d_edges = np.histogramdd(data.T, bins=n_bins,
                                            range=[(hist_limits[0], hist_limits[1]), ] * 3, density=False)
    if probability:
        hist_3d /= hist_3d.sum()

    fig, ax = plt.subplots(figsize=fig_size)

    if vrange is not None:
        h = ax.imshow(np.rot90(hist_3d.sum(axis=plane_to_axis[plane])),
                      cmap='plasma', vmin=vrange[0], vmax=vrange[1])
    else:
        h = ax.imshow(np.rot90(hist_3d.sum(axis=plane_to_axis[plane])),
                      cmap='plasma')

    # add colorbar and axis
    plt.colorbar(h, ax=ax)
    ax.set_xlabel('Relative x-coordinate $(\AA)$')
    ax.set_ylabel('Relative y-coordinate $(\AA)$')

    ticks = [0, hist_3d.shape[0] // 2, hist_3d.shape[0] - 1]
    xlabels = [hist_limits[0], 0, hist_limits[1]]
    ylabels = [hist_limits[1], 0, hist_limits[0]]

    ax.set_xticks(ticks)
    ax.set_xticklabels(xlabels)

    ax.set_yticks(ticks)
    ax.set_yticklabels(ylabels)

    return fig, ax, hist_3d, hist_3d_edges


def scatter_plot(data, axis_limits, fig_size=(5, 5), plane='xy', msize=0.01):
    """
    @param data: (3, N) shape array, x, y, z coordinates with A units
    @param axis_limits: tuple of two elements giving min and max value for all axis
    """
    assert data.shape[0] == 3, 'data does not have x, y, z as the second axis'

    fig, ax = plt.subplots(figsize=fig_size)

    if plane == 'xy':
        ax.scatter(data[0], data[1], s=msize)
    elif plane == 'xz':
        ax.scatter(data[0], data[2], s=msize)
    elif plane == 'yz':
        ax.scatter(data[1], data[2], s=msize)
    else:
        print('invalid plane')
        sys.exit(0)

    # add colorbar and axis
    ax.set_xlabel(f'Relative {plane[0]}-coordinate $(\AA)$')
    ax.set_ylabel(f'Relative {plane[1]}-coordinate $(\AA)$')

    ax.set_xlim(*axis_limits)
    ax.set_ylim(*axis_limits)

    return fig, ax


def plot_3d(data, plane_norm, axis_limits, msize=0.01):
    # ======= plot
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111, projection='3d')

    if data.shape[1] > 20000:
        display = data[:, np.random.choice(data.shape[1], size=20000)]
    else:
        display = data

    ax.scatter(display[0], display[1], display[2], s=msize*100)
    ax.quiver(0, 0, 0, *plane_norm, length=100, color='red', label='membrane plane')
    ax.quiver(0, 0, 0, 0, 0, 1, length=100, color='blue', label='plane correction')

    ax.set_xlabel('x $(\AA)$')
    ax.set_xlim(*axis_limits)
    ax.set_ylabel('y $(\AA)$')
    ax.set_ylim(*axis_limits)
    ax.set_zlabel('z $(\AA)$')
    ax.set_zlim(*axis_limits)
    return fig, ax


def get_classifier(subtomos, center_subtomos):

    classifier = []

    for s in subtomos:
        if s in center_subtomos:
            classifier.append(1)
        else:
            classifier.append(0)

    return classifier


if __name__ == '__main__':
    helper = ScriptHelper2(
        sys.argv[0].split('/')[-1],  # script name
        description='Inspect neighbour distribution around particles in 3d. Can be used to inspect polysome '
                    'assocation of ribosomes, for example.',
        authors='Marten Chaillet',
        options=[
            ScriptOption2(['-f', '--file'], '.star (relion) or .xml (pytom) particle list ', 'file', 'required'),
            ScriptOption2(['--center-class'], 'Particle list with subset of main file that will be used as center '
                                              'class. Currently only supported for star file',
                          'file', 'optional'),
            ScriptOption2(['-d', '--destination'], 'Folder where output graphs are stored.', 'directory',
                          'optional', '.'),
            ScriptOption2(['-o', '--output-name'], 'Base name of png files to write as output (do not write '
                                                   'extension).', 'string', 'optional'),
            ScriptOption2(['-s', '--pixel-size'], 'The pixel spacing of coordinates in the particle list (in A). If '
                                                  'a .star file is provided, script will attempt to read it '
                                                  'automatically.', 'float', 'optional'),
            ScriptOption2(['--particle-diameter'], 'Diameter of the particle for restricting plot sizes (in A).',
                          'float', 'optional'),
            # ScriptOption2(['-v', '--views'], 'List of output views to produce, default is xy and xz. Options are xy, '
            #                                  'xz and yz.', 'string,string', 'optional', ['xy', 'xz']),
            ScriptOption2(['--estimate-plane'], 'Fit plane to neighbour positions and rotate it into the xy '
                                                'plane. Can be used to rotate points on membrane into xy.',
                          'no arguments', 'optional'),
            ScriptOption2(['--plane-norm'], 'Provide a normal for the xy plane rotation, instead of estimating the '
                                            'plane. Script echos estimated plane norm in previous calculations, '
                                            'so you can put it back here to get a consistent rotation over multiple '
                                            'particle lists. Default is no rotation.', 'float,float,float', 'optional',
                          [0, 0, 1]),
            ScriptOption2(['--neighbourhood'], 'Number of neighbours to plot around each particle.', 'int',
                          'optional', 4),
            ScriptOption2(['--density'], 'Plot side views as density map instead of point cloud.', 'no arguments',
                          'optional'),
            ScriptOption2(['--view3d'], 'Opens interactive 3d plot point cloud that can be rotated and zoomed.',
                          'no arguments', 'optional'),
            ScriptOption2(['--marker-size'], 'Size of markers for 2d and 3d scatter plots, increase for small lists.',
                          'float', 'optional', 0.001)
        ])

    options = parse_script_options2(sys.argv[1:], helper)

    input_file, center_file, destination, output_name, pixel_size, particle_diameter, estimate_plane, plane_norm, \
        neighbourhood, plot_density, plot_3d_interactive, marker_size = options
    views_list = ['xy', 'xz', 'yz']

    # set output path if not specified with a name
    if output_name is None:
        output_name = os.path.splitext(os.path.split(input_file)[1])[0]  # base input name will be used
    output_path = os.path.join(destination, output_name)

    subtomogram_names = None
    # read the particle list
    if os.path.splitext(input_file)[1] == '.xml':
        if pixel_size is None:
            print('WARNING: pixel size not provided with pytom particle lists, axis of plots should not be '
                  'interpreted with distances. Setting pixel size to 1A.')
            pixel_size = 1

        # read the .xml particle list
        particle_list = ParticleList()
        particle_list.fromXMLFile(input_file)

        # create list for storing the data
        tomograms, coordinates, rotations = [], [], []

        # loop over the particles
        for particle in particle_list:
            pick_position, shift, rotation = particle.getPickPosition(), particle.getShift(), particle.getRotation()

            # shift position
            pick_position + shift.toVector()

            # get the tomogram
            tomogram = pick_position.getOriginFilename()
            # get coordinate
            coordinate = pick_position.toVector()
            # get the angles in zyz notation
            # put them in RELION notation by taking tranpose of matrix and negative (other option is to do this with
            # relion angles of course)
            zyz_angles = matToZYZ(rotation.toMatrix())  # returns a list

            # add the information to the lists
            tomograms.append(tomogram)
            coordinates.append(coordinate)
            rotations.append(zyz_angles)

        # convert arrays to numpy
        coordinates = np.array(coordinates) * pixel_size
        rotations = np.array(rotations)
        tomograms = np.array(tomograms)

    else:

        data = read(input_file)

        assert np.all(data['PixelSize'] == data['PixelSize'][0]), 'pixel size not identical over dataset'
        if pixel_size is None:
            pixel_size = data['PixelSize'][0]
            print(f'No pixel size provided, reading from star file ({pixel_size}A).')

        # get subtomogram names
        subtomogram_names = data['ImageName']

        # tomograms are in micrograph name folder
        tomograms = np.array([os.path.split(name)[1].split('.')[0] for name in data['MicrographName']])
        # OriginXAngst, etc, should already be in A according to RELION docs
        # see https://relion.readthedocs.io/en/release-3.1/Reference/Conventions.html
        shifts = np.array([data['OriginXAngst'], data['OriginYAngst'], data['OriginZAngst']]).T
        coordinates = np.array([data['CoordinateX'],
                                data['CoordinateY'],
                                data['CoordinateZ']]).T * pixel_size + shifts
        rotations = np.array([data['AngleRot'], data['AngleTilt'], data['AnglePsi']]).T

    center_class_index = None

    if center_file is not None and subtomogram_names is not None:
        data_center = read(center_file)

        # get center subtomo names
        subtomogram_names_center = data_center['ImageName']

        center_class_index = np.array(get_classifier(subtomogram_names, subtomogram_names_center))

    # run the neighbour rotation
    relative_coordinates, plane_norm = neighbour_position_3d(tomograms, coordinates, rotations,
                                                             neighbourhood=neighbourhood,
                                                             center_class=center_class_index,
                                                             plane_norm=plane_norm,
                                                             estimate_membrane_plane=estimate_plane)

    # determine plot limits
    if particle_diameter is not None:
        axis_size = particle_diameter * 3
        lim = int(axis_size - axis_size / 2)
        limits = (-lim, lim)
    else:
        limits = relative_coordinates.min(), relative_coordinates.max()

    # make interactive 3d plot
    if plot_3d_interactive:
        plot_3d(relative_coordinates, plane_norm, limits, msize=marker_size)
        plt.show()

    # make 2d plots
    if plot_density:
        for view in views_list:
            density_plot(relative_coordinates, limits, 15, fig_size=(5, 5), probability=True, plane=view)
            plt.tight_layout()
            plt.savefig(output_path + '_density_' + view + '.png', dpi=300)
    else:
        for view in views_list:
            scatter_plot(relative_coordinates, limits, fig_size=(5, 5), plane=view, msize=marker_size)
            plt.tight_layout()
        plt.savefig(output_path + '_scatter_' + view + '.png', dpi=300)
