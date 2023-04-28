#!/usr/bin/env pytom
"""
this is code from a bachelor project. at this point, the code still has fundamental flaws.
when unit tests have been supplied, this disclaimer can be removed.
25.12.2019 - FF
"""
"""
Created on Oct 01, 2019

@author: dschulte
"""

"""
create a function that gets the projection alignment transformations (alignmentresults.txt) and uses those to find the 
2d position on the original data given the tilt angle and the 3d position of the point of interest, also it should
provide a transformation matrix to be used to transform a template of the particle to the same point of view ad the 2d
position, for this it also needs the 3d rotation of the particle.
"""

import numpy as np

def test(index=0):
    from pytom.agnostic.io import read
    from pytom.agnostic.transform import rotate3d

    path_raw_projection = "/data2/dschulte/BachelorThesis/Data/VPP2/03_Tomographic_Reconstruction/tomogram_000/sorted/sorted_29.em"
    path_aligned_projection = "/data2/dschulte/BachelorThesis/Data/VPP2/03_Tomographic_Reconstruction/tomogram_000/alignment/marker_0001_-60.0,60.0/sorted_aligned_30.em"
    path_template = "/data2/dschulte/BachelorThesis/Data/VPP2/05_Subtomogram_Analysis/combo_reduced.em"
    tilt_angle = -1.9989999533
    raw_projection = read(path_raw_projection)
    aligned_projection = read(path_aligned_projection)
    template = read(path_template)
    dim = aligned_projection.shape[0]

    from pytom.basic.structures import ParticleList

    particlelist1 = ParticleList()
    particlelist1.fromXMLFile("/data2/dschulte/BachelorThesis/Data/VPP2/04_Particle_Picking/Picked_Particles/combined_reduced_extracted/particleList_combined_reduced_tomogram_010_WBP.xml")

    align_results = (-15.3044300079, 3.7495634556, -184.3835906982, 1.0000053644) # -2
    #align_results = (-1.2395387888, 4.9647006989, -184.3754882812, 1.0000000000) # 0
    #pick_position = (278.12318382089245 * 8, 222.6395540890773 * 8, 268.97256780848085 * 8) #0
    #pick_position = (381.21906883806 * 8, 153.61353397521387 * 8, 246.8315433927568 * 8) #74
    #particle_rotation = (26.442828473505173, 44.44149544840194, 58.160958298848676) #0
    #particle_rotation = (85.2456894956599, 30.815061362336394, 9.543300915975514) #74

    pick_position = particlelist1[index].getPickPosition().toVector()
    particle_rotation = (particlelist1[index].getRotation().getZ1(), particlelist1[index].getRotation().getX(), particlelist1[index].getRotation().getZ2())

    print(pick_position, particle_rotation)

    align_transformation, raw_position, aligned_position, template_transformation = combine_trans_projection(align_results, pick_position, particle_rotation, tilt_angle, dim, 8, template.shape[0])

    d = 100
    raw_patch = raw_projection[int(raw_position[0]-d):int(raw_position[0]+d), int(raw_position[1]-d):int(raw_position[1]+d)].squeeze()
    raw_patch = raw_patch / np.mean(raw_patch)

    aligned_patch = aligned_projection[int(aligned_position[0]-d):int(aligned_position[0]+d), int(aligned_position[1]-d):int(aligned_position[1]+d)].squeeze()
    aligned_patch = aligned_patch / (np.mean(aligned_patch))

    aligned_raw = matrix_apply_to_2d(aligned_projection.squeeze(), align_transformation)
    aligned_raw_patch = aligned_raw[int(raw_position[0]-d):int(raw_position[0]+d), int(raw_position[1]-d):int(raw_position[1]+d)].squeeze()
    aligned_raw_patch = aligned_raw_patch / (np.mean(aligned_raw_patch))

    transformed_template = matrix_apply_to_3d_3x3(template, template_transformation)
    print(np.mean(transformed_template), np.mean(template))
    template_2d = transformed_template.sum(axis=2)
    template_2d = template_2d / np.mean(template_2d)
    print(template_2d.shape)

    template_vol = rotate3d(template, phi=particle_rotation[0], the=particle_rotation[1], psi=particle_rotation[2])
    template_vol = rotate3d(template_vol, the=tilt_angle)
    template_vol_2d = template_vol.sum(axis=2)

    from pytom.reconstruction.reconstruct_local_alignment import normalised_cross_correlation_numpy, find_sub_pixel_max_value_2d

    raw_cc = normalised_cross_correlation_numpy(raw_patch, template_2d)
    rx, ry, _ = find_sub_pixel_max_value_2d(raw_cc)

    aligned_cc = normalised_cross_correlation_numpy(aligned_patch, template_vol_2d)
    tx, ty, _ = find_sub_pixel_max_value_2d(aligned_cc)

    import pylab as pl
    import scipy

    f, ax = pl.subplots(2, 3,figsize=(15,10))

    for i in range(2):
        for j in range(3):
            ax[i][j].axis('off')

    ax[0][0].set_title('Raw Data Particle')
    ax[0][0].imshow(scipy.ndimage.gaussian_filter(raw_patch, 3))
    ax[0][1].set_title('Template Transformed to Raw Data')
    ax[0][1].imshow(template_2d)
    ax[0][2].set_title('Cross Correlation')
    ax[0][2].imshow(raw_cc)
    ax[0][2].text(0.05, 0.05, 'Peak\nx: {0:.2f}\ny: {0:.2f}'.format(rx, ry), transform=ax[0][2].transAxes, color='white')
    ax[1][0].set_title('Aligned Data Particle')
    ax[1][0].imshow(scipy.ndimage.gaussian_filter(aligned_patch, 3))
    ax[1][1].set_title('Template Aligned to Aligned Data')
    ax[1][1].imshow(template_vol_2d)
    ax[1][2].set_title('Cross Correlation')
    ax[1][2].imshow(aligned_cc)
    ax[1][2].text(0.05, 0.05, 'Peak\nx: {0:.2f}\ny: {0:.2f}'.format(tx, ty), transform=ax[1][2].transAxes, color='white')
    f.tight_layout()
    pl.show()


def combine_trans_projection(tx, ty, rot, mag, x, y, z, phi, the, psi, tiltangle, dim, binning, particle_dim=200):
    """
    Combines all transfomrations of the raw o aligned and the position and rotations of a particle into one single matrix
    to effeciently and with only one round of interpolation be able to compare a template with raw data.

    @param tiltangle: The tiltangle
    @param dim: The dimensions of the raw data (assumes a square)
    @param binning: The binning factor of the coordinates of the particle position
    @param: particle_dim: The dimension of the template (assumes a cube)
    @return: (align transformations, particle position in raw data, particle position in aligned data, matrix for template)

    @author: Douwe Schulte and Gijs van der Schot
    """
    from numpy import cos, sin, pi

    # Calculates the inverse transformation matrix of the projection alignment transformations
    alpha = -rot * pi/180
    c = cos(alpha)
    s = sin(alpha)

    rotate    = np.matrix([[c, s, 0],  [-s, c, 0],   [0, 0, 1]])
    magnify   = np.matrix([[mag, 0, 0], [0, mag, 0], [0, 0, 1]])
    translate = np.matrix([[1, 0, tx],  [0, 1, ty],  [0, 0, 1]])

    align_transformations = np.linalg.inv(rotate * magnify * translate)

    # Map the 3D position to a 2D position on the projection of the tiltangle
    x = x * binning
    y = y * binning
    z = z * binning

    aligned_y = y  # assume the rotation axis is around y
    aligned_x = (cos(tiltangle * pi / 180) * (x - dim / 2) - sin(tiltangle * pi / 180) * (z - dim / 2)) + dim / 2

    # Use the projection alignment transformations to map this 2D position to a 2D position on the raw projections
    aligned_pos = np.matrix([[aligned_x - dim/2], [aligned_y - dim/2], [1]])
    raw_pos = align_transformations * aligned_pos

    # Calculate the rotation matrix for the template, a combination of the particle rotation and the tilt angle
    template_3d_rotation = generate_rotation_matrix(0, tiltangle, 0) * generate_rotation_matrix(phi, the, psi) * matrix_rotate_3d_z(rot) * matrix_magnify_3d(mag)

    # Merge this matrix with the projection transformations
    merged_matrix = template_3d_rotation

    return (align_transformations, (raw_pos.item(0, 0) + dim/2, raw_pos.item(1, 0) + dim/2), (aligned_x, aligned_y), merged_matrix)


def matrix_rotate_3d_x(deg):
    """Creates a 3d 3x3 rotation matrix for a deg turn on the x axis"""
    from numpy import cos, sin, pi
    rad_x = -deg * pi/180
    c_x = cos(rad_x)
    s_x = sin(rad_x)
    return np.matrix([[1, 0, 0], [0, c_x, -s_x], [0, s_x, c_x]])


def matrix_rotate_3d_y(deg):
    """Creates a 3d 3x3 rotation matrix for a deg turn on the y axis"""
    from numpy import cos, sin, pi
    rad_y = -deg * pi/180
    c_y = cos(rad_y)
    s_y = sin(rad_y)
    return np.matrix([[c_y, 0, s_y], [0, 1, 0], [-s_y, 0, c_y]])


def matrix_rotate_3d_z(deg):
    """Creates a 3d 3x3 rotation matrix for a deg turn on the z axis"""
    from numpy import cos, sin, pi
    rad_z = -deg * pi/180
    c_z = cos(rad_z)
    s_z = sin(rad_z)
    return np.matrix([[c_z, -s_z, 0], [s_z, c_z, 0], [0, 0, 1]])


def matrix_translate_3d(tx, ty, tz):
    """Creates a 3d 4x4 affine transformation matrix for a 3d translation"""
    return np.matrix([[1, 0, 0, tx], [0, 1, 0, ty], [0, 0, 1, tz], [0, 0, 0, 1]])


def matrix_magnify_3d(f):
    """Creates a 3d 3x3 rotation matrix for a magnification in every axis"""
    return np.matrix([[f, 0, 0], [0, f, 0], [0, 0, f]])


def matrix_2d_to_3d(matrix):
    """Calculates the 3d affine transformation matrix from the given 2d affine transformation matrix"""
    return np.matrix([
        [matrix.item(0, 0), matrix.item(0, 1), 0, matrix.item(0, 2)],
        [matrix.item(1, 0), matrix.item(1, 1), 0, matrix.item(1, 2)],
        [0,                 0,                 1, 0                ],
        [matrix.item(2, 0), matrix.item(2, 1), 0, matrix.item(2, 2)]])


def matrix_apply_to_3d_4x4(vol, matrix):
    from scipy import mgrid

    # Calculate the new coordinates of every point
    grid = mgrid[0.:vol.shape[0], 0.:vol.shape[1], 0.:vol.shape[2]]
    temp = grid.reshape((3, grid.size / 3))
    # Add the fourth dimension (just 1s but needed for the computations)
    newrow = np.ones(grid.size / 3)
    temp = np.vstack([temp, newrow])
    # Use the matrix to calculate the new positions of every point
    temp = np.dot(matrix, temp)
    # Delete the fourth dimension
    temp = np.delete(temp, 3, axis=0)
    temp = np.array(temp)
    grid = np.reshape(temp, (3, vol.shape[0], vol.shape[1], vol.shape[2]))

    from scipy.ndimage.interpolation import map_coordinates
    d = map_coordinates(vol, grid, order=3)

    return d


def matrix_apply_to_3d_3x3(vol, matrix):
    """Applies a given 3d 3x3 rotation matrix to the given volume, rotating around the center"""
    from scipy import mgrid

    cx = vol.shape[0]/2
    cy = vol.shape[1]/2
    cz = vol.shape[2]/2

    # Calculate the new coordinates of every point
    grid = mgrid[-cx:vol.shape[0]-cx, -cy:vol.shape[1]-cy, -cz:vol.shape[2]-cz]
    temp = grid.reshape((3, grid.size / 3))
    # Add the fourth dimension (just 1s but needed for the computations)
    # Use the matrix to calculate the new positions of every point
    temp = np.dot(matrix, temp)
    # Delete the fourth dimension
    temp = np.array(temp)
    grid = np.reshape(temp, (3, vol.shape[0], vol.shape[1], vol.shape[2]))

    grid[0] += cx
    grid[1] += cy
    grid[2] += cz

    from scipy.ndimage.interpolation import map_coordinates
    d = map_coordinates(vol, grid, order=3)

    return d


def matrix_apply_to_2d(data, matrix):
    """Applies a given 2d 2x2 rotation matrix to the given array, rotating around the center"""
    from scipy import mgrid

    cx = data.shape[0] / 2
    cy = data.shape[1] / 2

    # Calculate the new coordinates of every point
    grid = mgrid[-cx:data.shape[0]-cx, -cy:data.shape[1]-cy]
    temp = grid.reshape((2, grid.size / 2))
    # Add the fourth dimension (just 1s but needed for the computations)
    newrow = np.ones(grid.size / 2)
    temp = np.vstack([temp, newrow])
    # Use the matrix to calculate the new positions of every point
    temp = np.dot(matrix, temp)
    # Delete the fourth dimension
    temp = np.delete(temp, 2, axis=0)
    temp = np.array(temp)
    grid = np.reshape(temp, (2, data.shape[0], data.shape[1]))

    grid[0] += cx
    grid[1] += cy

    from scipy.ndimage.interpolation import map_coordinates
    d = map_coordinates(data, grid, order=3)

    return d


def generate_rotation_matrix(phi, the, psi):
    """Creates a 3d 3x3 rotation matrix with the given rotations in ZXZ notation."""
    # Transfer the angle to Euclidean
    phi = -float(phi) * np.pi / 180.0
    the = -float(the) * np.pi / 180.0
    psi = -float(psi) * np.pi / 180.0
    sin_alpha = np.sin(phi)
    cos_alpha = np.cos(phi)
    sin_beta = np.sin(the)
    cos_beta = np.cos(the)
    sin_gamma = np.sin(psi)
    cos_gamma = np.cos(psi)

    # Calculate inverse rotation matrix
    Inv_R = np.zeros((3, 3), dtype='float32')

    Inv_R[0, 0] = cos_alpha * cos_gamma - cos_beta * sin_alpha \
                  * sin_gamma
    Inv_R[0, 1] = -cos_alpha * sin_gamma - cos_beta * sin_alpha \
                  * cos_gamma
    Inv_R[0, 2] = sin_beta * sin_alpha

    Inv_R[1, 0] = sin_alpha * cos_gamma + cos_beta * cos_alpha \
                  * sin_gamma
    Inv_R[1, 1] = -sin_alpha * sin_gamma + cos_beta * cos_alpha \
                  * cos_gamma
    Inv_R[1, 2] = -sin_beta * cos_alpha

    Inv_R[2, 0] = sin_beta * sin_gamma
    Inv_R[2, 1] = sin_beta * cos_gamma
    Inv_R[2, 2] = cos_beta
    #Inv_R[3, 3] = 1

    return np.matrix(Inv_R)
