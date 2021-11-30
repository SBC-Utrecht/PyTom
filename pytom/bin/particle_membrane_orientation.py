#!/usr/bin/env pytom

"""
Determine orientation of membrane-bound particles in relation to the segmentation map of a membrane.

Author: Marten Chaillet
"""
import numpy as np
import sys


def convert_to_mesh(volume, cutoff=0.2, mesh_detail=2, display=False):
    from skimage import measure

    verts, faces, normals, values = measure.marching_cubes(volume, level=cutoff, step_size=mesh_detail)

    if display:
        try:
            import matplotlib.pyplot as plt
        except Exception as e:
            import matplotlib
            matplotlib.use('Qt5Agg')
            import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        # Display resulting triangular mesh using Matplotlib. This can also be done
        # with mayavi (see skimage.measure.marching_cubes_lewiner docstring).
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')

        # Fancy indexing: `verts[faces]` to generate a collection of triangles
        mesh = Poly3DCollection(verts[faces])
        mesh.set_edgecolor('k')
        ax.add_collection3d(mesh)

        ax.set_xlabel(f"x-axis: {volume.shape[0]}")
        ax.set_ylabel(f"y-axis: {volume.shape[1]}")
        ax.set_zlabel(f"z-axis: {volume.shape[2]}")

        ax.set_xlim(0, volume.shape[0])  # a = 6 (times two for 2nd ellipsoid)
        ax.set_ylim(0, volume.shape[1])  # b = 10
        ax.set_zlim(0, volume.shape[2])  # c = 16

        plt.tight_layout()
        plt.show()

    # invert normals because they point inwards
    return verts, faces, normals, values


def find_outliers(distances, nstd=3):
    norm_distances = (distances - distances.mean()) / distances.std()
    return np.logical_and(norm_distances >= -nstd, norm_distances <= nstd)


def find_orientations(plist, segmentation, cutoff, mesh_detail, reference_normal):
    from pytom.simulation.membrane import Vector
    from pytom_numpy import vol2npy
    from pytom.voltools.utils import rotation_matrix

    # get the triangular mesh
    verts, faces, normals, values = convert_to_mesh(segmentation, cutoff, mesh_detail)

    # load reference normal and make sure it is actually a normal
    unit_vector = Vector(reference_normal)
    unit_vector.normalize()

    distances, orientations = [], []
    particle_arrows, membrane_arrows = [], []
    for p in plist:
        # get rotation matrix and convert to axis-angle
        rotation = p.getRotation().toVector()  # z1, z2, x
        matrix = rotation_matrix(rotation=(rotation[0], rotation[2], rotation[1]), rotation_order='rzxz')
        # rotate unit vector by rotation matrix...
        particle_normal = Vector(unit_vector.get())  # copy reference normal
        particle_normal.rotate(matrix[:3, :3])

        # angle, vector = matToAxisAngle(rotation_matrix)
        # particle_normal = Vector(*vector)

        # get the coordinates of the particle
        coordinates = np.array(p.getPickPosition().toVector())

        # distance to each vertex in the triangle mesh
        distance = np.sqrt(np.sum(np.subtract(verts, coordinates) ** 2, axis=1))
        # distance = np.sqrt((verts[:, 0] - coordinates[0]) ** 2 +
        #                    (verts[:, 1] - coordinates[1]) ** 2 + (verts[:, 2] - coordinates[2]) ** 2)
        min_distance = np.min(distance)
        distances.append(min_distance)
        min_distance_idx = np.argmin(distance)

        # find all triangles that have closest_point index in them
        vertex_normal = Vector(normals[min_distance_idx])

        # get the difference angle
        difference_angle = particle_normal.angle(vertex_normal, degrees=True)
        orientations.append(difference_angle)

        particle_arrows.append((coordinates, particle_normal.get()))
        membrane_arrows.append((verts[min_distance_idx], vertex_normal.get()))

    return np.array(distances), np.array(orientations), particle_arrows, membrane_arrows


def write_arrow_bild(p_arrows, m_arrows, outlier_filter, filename):
    with open(filename, 'w') as stream:
        for (arrow, normal), inlier in zip(p_arrows, outlier_filter):
            if inlier:
                stream.write('.color red\n')
                stream.write(f'.arrow {arrow[0]:.2f} {arrow[1]:.2f} {arrow[2]:.2f} '
                             f'{arrow[0] + 15 * normal[0]:.2f} '
                             f'{arrow[1] + 15 * normal[1]:.2f} '
                             f'{arrow[2] + 15 * normal[2]:.2f}\n')
        for (arrow, normal), inlier in zip(m_arrows, outlier_filter):
            if inlier:
                stream.write('.color blue\n')
                stream.write(f'.arrow {arrow[0]:.2f} {arrow[1]:.2f} {arrow[2]:.2f} '
                             f'{arrow[0] + 15 * normal[0]:.2f} '
                             f'{arrow[1] + 15 * normal[1]:.2f} '
                             f'{arrow[2] + 15 * normal[2]:.2f}\n')
    # except :
    #     print('Unable to write error file.')
    #     return


if __name__ == '__main__':
    from pytom.tools.script_helper import ScriptHelper2, ScriptOption2
    from pytom.tools.parse_script_options import parse_script_options2
    from pytom.agnostic.io import read
    from pytom.basic.structures import ParticleList

    helper = ScriptHelper2(
        sys.argv[0].split('/')[-1],  # script name
        description='Find the orientations of particles with respect to a membrane. The particles should be in the '
                    'pytom particle list (xml) format. Membrane volume should be a segmented map of the membrane in '
                    'the tomogram, take care sizes correspond.',
        authors='Marten Chaillet',
        options=[
            ScriptOption2(['-p', '--particle_list'], 'XML type particle list.', 'file', 'required'),
            ScriptOption2(['-s', '--segmentation_file'], 'Membrane map in mrc/em/rec format.', 'file',
                          'required'),
            ScriptOption2(['-o', '--output_name'], 'Histogram plot output name. If not provided plot to screen.',
                          'string', 'optional'),
            ScriptOption2(['-c', '--cutoff'], 'Cutoff value for converting membrane model to a triangular mesh.',
                          'float', 'optional', 0.2),
            ScriptOption2(['-m', '--mesh_detail'], 'Detail of the mesh, i.e. how fine it should be sampled from the '
                                                   'volume.', 'int', 'optional', 2),
            ScriptOption2(['-n', '--template_normal'], 'Direction of normal vector of template that the orientations '
                                                       'are relative to. It will not affect the distribution of '
                                                       'populations that you find, only their angular difference with '
                                                       'the template.', 'float,float,float', 'optional', [.0, .0, 1.]),
            ScriptOption2(['-f', '--filter_nstd'], 'Number of standard deviations to filter outliers based on '
                                                   'distance to membrane mesh.', 'int', 'optional', 3)])

    options = parse_script_options2(sys.argv[1:], helper)

    particle_list_file, segmentation_file, output_file, cutoff, mesh_detail, template_normal, nstd = options

    segmentation = read(segmentation_file)
    particle_list = ParticleList()
    particle_list.fromXMLFile(particle_list_file)

    distances, orientations, p_arrows, m_arrows = find_orientations(particle_list, segmentation, cutoff, mesh_detail,
                                                                    template_normal)

    outlier_filter = find_outliers(distances, nstd)
    distances = distances[outlier_filter]
    orientations = orientations[outlier_filter]

    # ======== create plot
    try:
        import matplotlib.pyplot as plt
    except Exception as e:
        import matplotlib
        matplotlib.use('Qt5Agg')
        import matplotlib.pyplot as plt

    # do some plotting
    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    ax[0].hist(distances)
    # ax[0].set_label('distance')
    ax[0].set_xlabel('distance (voxels)')
    ax[0].set_ylabel('number of particles')
    ax[1].hist(orientations)
    # ax[1].set_label('orientation')
    ax[1].set_xlabel('angle (degrees)')

    plt.tight_layout()

    if output_file is not None:
        plt.savefig(output_file + '.png', dpi=200, format='png')
        print(f"wrote {output_file + '.png'}")
        write_arrow_bild(p_arrows, m_arrows, outlier_filter, output_file + '_vectors.bild')
        print(f"wrote {output_file + '_vectors.bild'}")
    else:
        plt.show()

