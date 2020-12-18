"""
Membrane geometrical structures

Depedencies: pyvista

Author: Marten Chaillet
"""
# essential
import numpy as xp
import pyvista
from pytom.tompy.io import write

# plotting
import matplotlib
matplotlib.use('Qt5Agg')
from pylab import *
from mpl_toolkits.mplot3d import Axes3D # <--- This is important for 3d plotting


class Vector:
    # Class can be used as both a 3d coordinate, and a vector
    def __init__(self, vx, vy, vz):
        self.axis = xp.array([vx, vy, vz])

    def get_x(self):
        return self.axis[0]

    def set_x(self, x):
        self.axis[0] = x

    def get_y(self):
        return self.axis[1]

    def set_y(self, y):
        self.axis[1] = y

    def get_z(self):
        return self.axis[2]

    def set_z(self, z):
        self.axis[2] = z

    def show(self):
        print(self.axis)

    def cross(self, other):
        return Vector(self.axis[1] * other.axis[2] - self.axis[2] * other.axis[1],
                      self.axis[2] * other.axis[0] - self.axis[0] * other.axis[2],
                      self.axis[0] * other.axis[1] - self.axis[1] * other.axis[0])

    def dot(self, other):
        # return the dot product of vectors v1 and v2, of form (x,y,z)
        # dot product of two vectors is zero if they are perpendicular
        return self.axis[0] * other.axis[0] + self.axis[1] * other.axis[1] + self.axis[2] * other.axis[2]

    def magnitude(self):
        # calculate the magnitude (length) of vector p
        return xp.sqrt(xp.sum(self.axis ** 2))

    def normalize(self):
        self.axis = self.axis / self.magnitude()

    def angle(self, other, degrees=False):
        # returns angle in radians
        angle = xp.arccos(self.dot(other) / (self.magnitude() * other.magnitude()))
        if degrees:
            return angle * 180 / xp.pi
        else:
            return angle


def matrix_from_axis_angle(ax, angle):
    # ax input is of type Vector
    # angle input should be in radians
    # calculate the rotation matrix from and axis-angle (euler) representation
    ax.normalize()

    x, y, z = ax.axis
    c = xp.cos(angle)
    s = xp.sin(angle)
    t = 1.0 - c

    m00 = c + x * x * t
    m11 = c + y * y * t
    m22 = c + z * z * t

    tmp1 = x * y * t
    tmp2 = z * s
    m10 = tmp1 + tmp2
    m01 = tmp1 - tmp2
    tmp1 = x * z * t
    tmp2 = y * s
    m20 = tmp1 - tmp2
    m02 = tmp1 + tmp2
    tmp1 = y * z * t
    tmp2 = x * s
    m21 = tmp1 + tmp2
    m12 = tmp1 - tmp2
    return xp.array([[m00, m01, m02], [m10, m11, m12], [m20, m21, m22]])


def z_axis_rotation_matrix(angle):
    m00 = xp.cos(angle*xp.pi/180)
    m01 = - xp.sin(angle*xp.pi/180)
    m10 = xp.sin(angle*xp.pi/180)
    m11 = xp.cos(angle*xp.pi/180)
    return xp.array([[m00,m01,0], [m10,m11,0], [0,0,1]])


def get_rotation(vector, reference):
    # get the rotations to rotate vector on top of reference
    # both inputs are of type Vector
    vector.normalize()
    reference.normalize()

    axis = vector.cross(reference)
    angle = vector.angle(reference)

    return matrix_from_axis_angle(axis, angle)


def random_point_ellipsoid(a, b, c):
    # a,b, and c are paremeters of the ellipsoid
    # generating random (x,y,z) points on ellipsoid
    u = xp.random.rand()
    v = xp.random.rand()
    theta = u * 2.0 * xp.pi
    phi = xp.arccos(2.0 * v - 1.0)
    sinTheta = xp.sin(theta)
    cosTheta = xp.cos(theta)
    sinPhi = xp.sin(phi)
    cosPhi = xp.cos(phi)
    rx = a * sinPhi * cosTheta
    ry = b * sinPhi * sinTheta
    rz = c * cosPhi
    return xp.array([rx, ry, rz])


def sample_points_ellipsoid(number, a=2, b=3, c=4):
    points = random_point_ellipsoid(a, b, c)
    for i in range(number - 1):
        point = random_point_ellipsoid(a, b, c)
        points = xp.vstack((points, point))
    return points


def triangulate(points, alpha):
    # points is a 3D numpy array (n_points, 3) coordinates of a sphere
    cloud = pyvista.PolyData(points)
    # cloud.plot()

    volume = cloud.delaunay_3d(alpha=alpha)
    shell = volume.extract_geometry()
    # shell.plot()
    return shell


def centroid(triangle):
    # cetroid of the three 3D points a, b, and  c
    # a,b, and c are numpy array of length 3
    return (1 / 3) * (triangle[0] + triangle[1] + triangle[2])


def rotate_point(point, rotation_matrix):
    new_point = xp.matmul(rotation_matrix, point.reshape(-1, 1))
    return new_point.reshape(1, -1).squeeze()


def rotate_triangle(triangle, matrix):
    return xp.vstack((rotate_point(triangle[i], matrix) for i in range(3)))


def shift_point(point, shift):
    return point + shift


def shift_triangle(triangle, shift):
    rtriangle = triangle
    rtriangle[0] = shift_point(rtriangle[0], shift)
    rtriangle[1] = shift_point(rtriangle[1], shift)
    rtriangle[2] = shift_point(rtriangle[2], shift)
    return rtriangle


def display_points(points, zlim=0):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(points[:, 0], points[:, 1], points[:, 2])
    #     fig.show()
    if zlim:
        ax.set_zlim3d(-zlim, zlim)
    return


def display_vectors(v1, v2):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot([0, v1[0]], [0, v1[1]], [0, v1[2]], color='blue')
    ax.plot([0, v2[0]], [0, v2[1]], [0, v2[2]], color='red')
    ax.set_xlim3d(-1, 1)
    ax.set_ylim3d(-1, 1)
    ax.set_zlim3d(-1, 1)
    #     fig.show()
    return


def display_triangle_normal(triangle, normal):
    # first center triangle around origin
    center = centroid(triangle)
    centered_triangle = shift_triangle(triangle, -center)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(centered_triangle[:, 0], centered_triangle[:, 1], centered_triangle[:, 2])
    ax.plot([0, normal[0]], [0, normal[1]], [0, normal[2]], color='blue')
    ax.set_xlim3d(-.5, .5)
    ax.set_ylim3d(-.5, .5)
    ax.set_zlim3d(-.5, .5)
    return


def boundary_box_from_pdb(filename):
    try:
        with open(filename, 'r') as pdb:
            line = pdb.readline()
            while line.split()[0] != 'CRYST1':
                # note: if we do not encounter CRYST1 in the file, we go to the except statement.
                line = pdb.readline()
        return float(line.split()[1]), float(line.split()[2]), float(line.split()[3])
    except Exception as e:
        print(e)
        raise Exception('Could not read pdb file.')


def sign(p1, p2, p3):
    """
    Determine on which side of the line formed by p2 and p3, p1 lays.
    """
    return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])


def point_in_triangle(pt, triangle):
    """
    All points are numpy arrays of length 2, so points in 2D.
    @param pt:
    @param v1:
    @param v2:
    @param v3:
    @return:
    """
    v1 = triangle[0]
    v2 = triangle[1]
    v3 = triangle[2]
    d1 = sign(pt, v1, v2)
    d2 = sign(pt, v2, v3)
    d3 = sign(pt, v3, v1)

    has_neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
    has_pos = (d1 > 0) or (d2 > 0) or (d3 > 0)

    return not(has_neg and has_pos)


def membrane_potential(surface_mesh, voxel_size, membrane_pdb, solvent, voltage):
    """

    @param ellipsoid_mesh: this is a pyvista surface
    @param voxel_size:
    @param membrane_pdb:
    @param solvent_potential:
    @return:
    """
    from potential import read_structure, iasa_integration

    reference = Vector(.0, .0, 1.0)

    membrane_x, membrane_y, membrane_z, membrane_e, membrane_b, membrane_o = [], [], [], [], [], []

    atom_count = 0

    for icell in range(surface_mesh.n_cells):
        print(f'Triangle {icell+1} out of {surface_mesh.n_cells}.')

        triangle = surface_mesh.extract_cells(icell).points
        # print(triangle)
        normal = Vector(*surface_mesh.cell_normals[icell])

        center = centroid(triangle)
        triangle1 = shift_triangle(triangle, -center)

        matrix1 = get_rotation(normal, reference)
        triangle2 = rotate_triangle(triangle1, matrix1)

        # add a random rotation of the triangle in the x-y plane to rotate the membrane's coordinates to other locations
        angle = xp.random.uniform(0,360)
        matrix2 = z_axis_rotation_matrix(angle)
        triangle3 = rotate_triangle(triangle2, matrix2)

        # apply a shift to the points so that the coordinates are all above the origin
        shift = xp.array([xp.min(triangle3[:, 0]), xp.min(triangle3[:, 1]), 0])
        triangle_sample = shift_triangle(triangle3, -shift)
        # print(xp.round(triangle_sample, decimals=2))

        # READ THE STRUCTURE AND EXTEND IT ONLY ONCE
        # membrane pdb should have solvent deleted at this point
        x_coordinates, y_coordinates, z_coordinates, elements, b_factors, occupancies = read_structure(membrane_pdb)
        z_coordinates = list(xp.array(z_coordinates) - sum(z_coordinates)/len(z_coordinates))
        # get the periodice boundary box from the pdb file
        x_bound, y_bound, z_bound = boundary_box_from_pdb(membrane_pdb)

        xmax = xp.max(triangle_sample[:,0])
        ymax = xp.max(triangle_sample[:,1])

        xext = int(xp.ceil(xmax / x_bound))
        yext = int(xp.ceil(ymax / y_bound))
        newx, newy, newz, newe, newb, newo = [], [], [], [], [], []
        for i in range(xext):
            for j in range(yext):
                newx += list(xp.array(x_coordinates) + i*x_bound)
                newy += list(xp.array(y_coordinates) + j*y_bound)
                newz += z_coordinates
                newe += elements
                newb += b_factors
                newo += occupancies
        # at this points our triangle fits onto the structure
        # now delete all the atoms that fall outside of the triangle
        i = 0
        while i < len(newx):
            coordinate = xp.array([newx[i], newy[i]])
            if not point_in_triangle(coordinate, triangle_sample[:,:2]):
                newx.pop(i)
                newy.pop(i)
                newz.pop(i)
                newe.pop(i)
                newb.pop(i)
                newo.pop(i)
            else:
                i += 1
        # now we have all the atoms needed to build the membrane model for the current triangle
        # lets first rotate the atoms to the original orientation of the triangle

        # applying the rotation to each atom individually might me a bit slow...
        n_atoms = len(newe)
        matrix1_t = matrix1.T
        matrix2_t = matrix2.T
        for i in range(n_atoms):
            atom = xp.array([newx[i], newy[i], newz[i]])
            catom = rotate_point(shift_point(atom, shift), matrix2_t)
            atom_surface = shift_point(rotate_point(catom, matrix1_t), center)
            newx[i] = atom_surface[0]
            newy[i] = atom_surface[1]
            newz[i] = atom_surface[2]

        # add positions to larger array that describes atom coordinates on the full surface
        membrane_x += newx
        membrane_y += newy
        membrane_z += newz
        membrane_e += newe
        membrane_b += newb
        membrane_o += newo

        atom_count += len(newe)
        print(f'Number of atoms added {len(newe)} to a total count of {atom_count}.')

    structure = (membrane_x, membrane_y, membrane_z, membrane_e, membrane_b, membrane_o)
    # pass directly to iasa_integration
    potential = iasa_integration('', voxel_size, solvent_exclusion=True, V_sol=solvent,
                                 absorption_contrast=True, voltage=voltage, density=0.92, molecular_weight=734.1,
                                 structure_tuple=structure)
    return potential


if __name__ == '__main__':

    folder = '/data2/mchaillet/structures'

    size_factor = 3

    # automatically scale these points
    N = 300 * size_factor^3# number of points
    a, b, c = (x*size_factor for x in (xp.random.randint(180, 280), xp.random.randint(180, 280),
                                       xp.random.randint(180, 280)))
    alpha = 2000 * size_factor
    voxel = 2.62 * 4 # A
    pdb = f'{folder}/pdb/lipid/dppc128_dehydrated.pdb'
    solvent = 4.5301
    voltage = 300E3

    # generate an ellipsoid and triangulate it
    points = sample_points_ellipsoid(N, a=a, b=b, c=c)
    surface = triangulate(points, alpha)

    # fill the triangles with lipid molecules and calculate potential for it
    volume = membrane_potential(surface, voxel, pdb, solvent, voltage)
    real = volume[0]
    imag = volume[1]

    from potential import reduce_resolution, bin

    name = 'bilayer'
    size = f'{a*2/10:.0f}x{b*2/10:.0f}x{c*2/10:.0f}nm' # double the values of the ellipsoid radii for actual size

    real_fil = reduce_resolution(real, voxel, 2 * voxel)
    imag_fil = reduce_resolution(imag, voxel, 2 * voxel)

    write(f'{folder}/potential/membrane/{name}_{voxel:.2f}A_{size}_4.53V_real.mrc', real)
    write(f'{folder}/potential/membrane/{name}_{voxel:.2f}A_{size}_4.53V_imag_300V.mrc', imag)

    binning = 2

    real_bin = bin(reduce_resolution(real, voxel, binning * voxel * 2), binning) # *2 still?
    imag_bin = bin(reduce_resolution(imag, voxel, binning * voxel * 2), binning)

    write(f'{folder}/potential/membrane/{name}_{voxel*binning:.2f}A_{size}_4.53V_real.mrc', real_bin)
    write(f'{folder}/potential/membrane/{name}_{voxel*binning:.2f}A_{size}_4.53V_imag_300V.mrc', imag_bin)

