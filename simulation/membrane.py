"""
Membrane geometrical structures

Depedencies: pyvista

Author: Marten Chaillet
"""

# essential
import numpy as xp
import pytom.simulation.physics as physics

# plotting
# TODO remove plotting import
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


def distance_nonsqrt(p1, p2):
    return sum( (p1 - p2) ** 2 )


def distance(p1, p2):
    dd = (p1 - p2)**2
    return np.sqrt(sum(dd))


def get_point_ellipsoid(a, b, c, theta, phi):
    sinTheta = xp.sin(theta)
    cosTheta = xp.cos(theta)
    sinPhi = xp.sin(phi)
    cosPhi = xp.cos(phi)
    # rx = a * sinPhi * cosTheta
    # ry = b * sinPhi * sinTheta
    # rz = c * cosPhi
    rx = a * cosPhi * cosTheta
    ry = b * cosPhi * sinTheta
    rz = c * sinPhi
    return xp.array([rx, ry, rz])


def random_point_ellipsoid(a, b, c):
    # a,b, and c are paremeters of the ellipsoid
    # generating random (x,y,z) points on ellipsoid
    u = xp.random.rand()
    v = xp.random.rand()
    theta = u * 2.0 * xp.pi
    phi = xp.arccos(2.0 * v - 1.0) - xp.pi / 2
    # phi = u * 2.0 * xp.pi
    # theta = xp.arccos(2.0 * v - 1.0)
    return get_point_ellipsoid(a, b, c, theta, phi)


def ellipsoid_intersection(thetaphi, a, b, c, x, y, z):
    theta, phi = thetaphi
    # phi, theta = thetaphi
    costheta = xp.cos(theta)
    sintheta = xp.sin(theta)
    cosphi = xp.cos(phi)
    sinphi = xp.sin(phi)
    eq1 = ((a**2 - b**2) * costheta * sintheta * cosphi
           - x * a * sintheta
           + y * b * costheta)
    eq2 = ((a**2 * costheta**2 + b**2 * sintheta**2 - c**2) * sinphi * cosphi
           - x * a * sinphi * costheta
           - y * b * sinphi * sintheta
           + z * c * cosphi)
    return eq1, eq2


def place_back_to_ellipsoid(point, a, b, c, newton_iter=3):

    # TODO Missing r parameter in function?
    # TODO Function does not optimize for 'strong' ellipsoids, possibly because of bad initial guess for theta phi

    from scipy.optimize import newton

    # def fprime(thetaphi, a, b, c, x, y, z):
    #     theta, phi = thetaphi
    #     costheta = xp.cos(theta)
    #     sintheta = xp.sin(theta)
    #     cosphi = xp.cos(phi)
    #     sinphi = xp.sin(phi)
    #     a2, b2 = a**2, b**2
    #     a11 = (a2 - b2) * (costheta**2 - sintheta**2) * cosphi - x * a * costheta - y * b * sintheta
    #     a12 = - (a2 - b2) * costheta * sintheta * sinphi
    #     a21 = - 2 * (a2 - b2) * costheta * sintheta * sinphi * cosphi + x * a * sinphi * sintheta - y * b * sinphi * \
    #           costheta
    #     a22 = (a2 * costheta**2 + b2 * sintheta**2 - c**2) * (cosphi**2 - sinphi**2) - x * a * cosphi * costheta - y \
    #           * b * cosphi * sintheta - z * c * sinphi
    #     return [[a11, a12], [a21, a22]]

    x, y, z = point
    thetaphi0 = [xp.arctan2(a*y, b*x), xp.arctan2(z, c * xp.sqrt((x/a)**2 + (y/b)**2))]
    angles_ellipsoid = newton(ellipsoid_intersection, thetaphi0, args=(a, b, c, x, y, z), maxiter=newton_iter)
    return get_point_ellipsoid(a, b, c, angles_ellipsoid[0], angles_ellipsoid[1])


def test_place_back_to_ellipsoid(size=100, a=1, b=1, c=1, iter=3):
    # xyz points between -1 and 1
    distances = xp.zeros((size, size, size))
    for i,x in zip(range(size), xp.arange(-a, a, (2.*a)/size)):
        for j,y in zip(range(size), xp.arange(-b, b, (2.*b)/size)):
            for k,z in zip(range(size), xp.arange(-c, c, (2.*c)/size)):
                point = xp.array([x,y,z])
                ellipsoid_point = place_back_to_ellipsoid(point, a, b, c, newton_iter=iter)
                distances[i,j,k] = distance(point, ellipsoid_point)
    # visualize
    slice_x = distances[size//2,:,:]
    slice_y = distances[:,size//2,:]
    slice_z = distances[:,:,size//2]
    fig, ax = plt.subplots(1, 3)
    ax[0].imshow(slice_x)
    ax[1].imshow(slice_y)
    ax[2].imshow(slice_z)
    plt.show()
    return


class DistanceMatrix:
    def __init__(self, points):
        from scipy.spatial.distance import cdist
        self.matrix = cdist(points, points, metric='sqeuclidean')  # squared euclidean because we compare distance
        # remove the points correlating with themselves
        self.upper = np.max(self.matrix)
        self.matrix[self.matrix == 0] = self.upper

    def update(self, points, new_point_index):
        dist_update = np.sum((points - points[new_point_index]) ** 2, axis=1)  # squared euclidean (see above)
        # remove point correlating with itself
        dist_update[dist_update == 0] = self.upper
        self.matrix[new_point_index, :] = dist_update
        self.matrix[:, new_point_index] = dist_update

    def shortest_distance(self):
        return xp.unravel_index(self.matrix.argmin(), self.matrix.shape)


def equilibrate_ellipsoid(points, a=2, b=3, c=4, maxiter=10000, factor=0.01, display=False):

    dmatrix = DistanceMatrix(points)

    for x in range(maxiter):
        # get the indices of the points that form the closest pair
        minp1, minp2 = dmatrix.shortest_distance()

        # move closest pair away from each other
        p1 = points[minp1].copy()
        p2 = points[minp2].copy()
        p1_new = p1 - factor * (p2 - p1)
        p2_new = p1 + (1+factor) * (p2 - p1)

        # use newton optimization to place the points back on the ellipsoid
        points[minp1] = place_back_to_ellipsoid(p1_new, a, b, c)
        points[minp2] = place_back_to_ellipsoid(p2_new, a, b, c)

        # update distance matrix with the new points
        dmatrix.update(points, minp1)
        dmatrix.update(points, minp2)

        # placed_back = xp.vstack((placed_back, points[minp1]))
        # placed_back = xp.vstack((placed_back, points[minp2]))

        if display:
            print(p1, p2)
            print(p1_new, p2_new)
            print(points[minp1], points[minp2])
            plt.close()
            display_points_3d(points)

    # display_points_3d(placed_back)

    return points


def sample_points_ellipsoid(number, a=2, b=3, c=4, evenly=True, maxiter=10000, factor=0.01, display=False):
    points = random_point_ellipsoid(a, b, c)
    for i in range(number-1):
        points = xp.vstack((points, random_point_ellipsoid(a, b, c)))
    if evenly:
        return equilibrate_ellipsoid(points, a, b, c, maxiter=maxiter, factor=factor, display=display)
    else:
        return points


def get_point_ellipse(a, b, theta):
    rx = a * xp.cos(theta)
    ry = b * xp.sin(theta)
    return xp.array([rx, ry])


def random_point_ellipse(a, b):
    u = xp.random.rand()
    theta = u * 2.0 * xp.pi
    return get_point_ellipse(a, b, theta)


def place_back_to_ellipse(point, a, b):

    from scipy.optimize import newton

    def f(theta, a, b, x, y):
        costheta = xp.cos(theta)
        sintheta = xp.sin(theta)
        return (a**2 - b**2) * costheta * sintheta - x * a * sintheta + y * b * costheta

    def fprime(theta, a, b, x, y):
        costheta = xp.cos(theta)
        sintheta = xp.sin(theta)
        return (a ** 2 - b ** 2) * (costheta**2 - sintheta**2) - x * a * costheta - y * b * sintheta

    x, y = point
    theta0 = xp.arctan2(a*y, b*x)
    angle_ellipse = newton(f, theta0, fprime=fprime, args=(a, b, x, y), maxiter=5)
    return get_point_ellipse(a, b, angle_ellipse)


def equilibrate_ellipse(points, a=2, b=3, maxiter=10000, factor=0.01, display=False):
    for x in range(maxiter):
        minp1, minp2 = 0, 1
        mind = distance_nonsqrt(points[minp1], points[minp2])
        maxd = mind
        # find closest two points
        for i in range(points.shape[0]-1):
            for j in range(i+1, points.shape[0]):
                d = distance_nonsqrt(points[i], points[j])
                if d < mind:
                    minp1, minp2 = i, j
                    mind = d
                if d > maxd:
                    maxd = d

        p1 = points[minp1].copy()
        p2 = points[minp2].copy()
        p1_new = p1 - factor * (p2 - p1)
        p2_new = p1 + (1+factor) * (p2 - p1)

        # use newton optimization to place the points back on the ellipsoid
        points[minp1] = place_back_to_ellipse(p1_new, a, b)
        points[minp2] = place_back_to_ellipse(p2_new, a, b)

        if display:
            plt.close()
            display_points_2d(points)

    return points


def sample_points_ellipse(number, a=2, b=3, evenly=True, maxiter=1000, factor=0.01, display=False):
    points = random_point_ellipse(a, b)
    for i in range(number-1):
        points = xp.vstack((points, random_point_ellipse(a, b)))
    if evenly:
        return equilibrate_ellipse(points, a, b, maxiter=maxiter, factor=factor, display=display)
    else:
        return points


def triangulate(points, alpha):
    import pyvista as pv
    # pyvista can also directly generate an ellipsoid
    # ellipsoid = pv.ParametricEllipsoid(10, 5, 5)
    # this returns a surface as pyvista.PolyData
    # delaunay 3d should work directly on this

    # points is a 3D numpy array (n_points, 3) coordinates of a sphere
    cloud = pv.PolyData(points)
    # cloud.plot()

    # reconstructs the surface from a set of points on an assumed solid surface
    # DataSetFilters.reconstruct_surface()
    # delaunay_3d should also work on a surface


    # for noise search for "pyvista perlin noise 3d"

    volume = cloud.delaunay_3d(alpha=alpha, progress_bar=True)
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


def display_points_2d(points):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(points[:, 0], points[:, 1])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.show()
    return


def display_points_3d(points, zlim=0):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(points[:, 0], points[:, 1], points[:, 2])
    #     fig.show()
    if zlim:
        ax.set_zlim3d(-zlim, zlim)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()
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


def point_array_sign(point_array, p2, p3):
    return (point_array[:, 0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (point_array[:, 1] - p3[1])


def point_array_in_triangle(point_array, triangle):

    v1 = triangle[0]
    v2 = triangle[1]
    v3 = triangle[2]

    d1 = sign(point_array, v1, v2)
    d2 = sign(point_array, v2, v3)
    d3 = sign(point_array, v3, v1)

    has_neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
    has_pos = (d1 > 0) or (d2 > 0) or (d3 > 0)

    return not(has_neg and has_pos)


def rotation_matrix_to_affine_matrix(rotation_matrix):
    m = np.identity(4)
    m[0:3, 0:3] = rotation_matrix
    return m


def membrane_potential_old(surface_mesh, voxel_size, membrane_pdb, solvent, voltage):
    """

    @param ellipsoid_mesh: this is a pyvista surface
    @param voxel_size:
    @param membrane_pdb:
    @param solvent_potential:
    @return:
    """
    #todo function is very slow

    from potential import read_structure, iasa_integration

    # READ THE STRUCTURE AND EXTEND IT ONLY ONCE
    # membrane pdb should have solvent deleted at this point
    x_coordinates, y_coordinates, z_coordinates, elements, b_factors, occupancies = read_structure(membrane_pdb)
    z_coordinates = list(xp.array(z_coordinates) - sum(z_coordinates) / len(z_coordinates))
    # get the periodice boundary box from the pdb file
    x_bound, y_bound, z_bound = boundary_box_from_pdb(membrane_pdb)

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
            if not point_in_triangle(coordinate, triangle_sample[:, :2]):
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
        # n_atoms = len(newe)
        # matrix1_t = matrix1.T
        # matrix2_t = matrix2.T

        atoms = xp.array([newx, newy, newz])
        atoms[0, :] += shift[0]
        atoms[1, :] += shift[1]
        atoms[2, :] += shift[2]
        # inverse of rotation matrix for the rotating the points back to original triangle position
        # inverse of atoms for correct order of applying rotation matrices
        ratoms = (matrix1.T @ (matrix2.T @ atoms))
        ratoms[0, :] += center[0]
        ratoms[1, :] += center[1]
        ratoms[2, :] += center[2]
        # newx = list(ratoms[0, :])
        # newy = list(ratoms[1, :])
        # newz = list(ratoms[2, :])

        # for i in range(n_atoms):
        #     atom = xp.array([newx[i], newy[i], newz[i]])
        #     catom = rotate_point(shift_point(atom, shift), matrix2_t)
        #     atom_surface = shift_point(rotate_point(catom, matrix1_t), center)
        #     newx[i] = atom_surface[0]
        #     newy[i] = atom_surface[1]
        #     newz[i] = atom_surface[2]

        # add positions to larger array that describes atom coordinates on the full surface
        membrane_x += list(ratoms[0, :])
        membrane_y += list(ratoms[1, :])
        membrane_z += list(ratoms[2, :])
        membrane_e += newe
        membrane_b += newb
        membrane_o += newo

        atom_count += len(newe)
        print(f'Number of atoms added {len(newe)} to a total count of {atom_count}.')

    structure = (membrane_x, membrane_y, membrane_z, membrane_e, membrane_b, membrane_o)
    # pass directly to iasa_integration
    potential = iasa_integration('', voxel_size, solvent_exclusion=True, V_sol=solvent,
                                 absorption_contrast=True, voltage=voltage, density=physics.PROTEIN_DENSITY,
                                 molecular_weight=physics.PROTEIN_MW, structure_tuple=structure)
    return potential


def membrane_potential(surface_mesh, voxel_size, membrane_pdb, solvent, voltage):
    """

    @param ellipsoid_mesh: this is a pyvista surface
    @param voxel_size:
    @param membrane_pdb:
    @param solvent_potential:
    @return:
    """
    #todo function is very slow

    from potential import read_structure, iasa_integration
    from pytom.voltools.utils import translation_matrix

    # READ THE STRUCTURE AND EXTEND IT ONLY ONCE
    # membrane pdb should have solvent deleted at this point
    x_coordinates, y_coordinates, z_coordinates, elements, b_factors, occupancies = read_structure(membrane_pdb)
    z_coordinates = list(xp.array(z_coordinates) - sum(z_coordinates) / len(z_coordinates))
    # get the periodice boundary box from the pdb file
    x_bound, y_bound, z_bound = boundary_box_from_pdb(membrane_pdb)

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

        xmax = xp.max(triangle_sample[:, 0])
        ymax = xp.max(triangle_sample[:, 1])

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

        atoms = xp.array([newx, newy, newz])
        newe = xp.array(elements)
        newb = xp.array(b_factors)
        newo = xp.array(occupancies)

        locs = point_array_in_triangle(atoms, triangle_sample[:, :2])
        atoms = atoms[locs]
        newe = newe[locs]
        newb = newb[locs]
        newo = newo[locs]

        # at this points our triangle fits onto the structure
        # now delete all the atoms that fall outside of the triangle
        # i = 0
        # while i < len(newx):
        #     coordinate = xp.array([newx[i], newy[i]])
        #     if not point_in_triangle(coordinate, triangle_sample[:, :2]):
        #         newx.pop(i)
        #         newy.pop(i)
        #         newz.pop(i)
        #         newe.pop(i)
        #         newb.pop(i)
        #         newo.pop(i)
        #     else:
        #         i += 1
        # now we have all the atoms needed to build the membrane model for the current triangle
        # lets first rotate the atoms to the original orientation of the triangle

        # applying the rotation to each atom individually might me a bit slow...
        # n_atoms = len(newe)
        # matrix1_t = matrix1.T
        # matrix2_t = matrix2.T
        t2 = translation_matrix(translation=-shift)  # voltools implements matrices as inverse operations, here we want
        t1 = translation_matrix(translation=-center)  # forward transformation
        affine_matrix = xp.matmul(xp.matmul(t1, matrix1.T), xp.matmul(matrix2.T, t2))

        ratoms = xp.matmul(affine_matrix, atoms)

        # atoms[0, :] += shift[0]
        # atoms[1, :] += shift[1]
        # atoms[2, :] += shift[2]
        # # inverse of rotation matrix for the rotating the points back to original triangle position
        # inverse of atoms for correct order of applying rotation matrices
        # np.matmul(matrix1.T, matrix2.T)
        # ratoms = (matrix1.T @ (matrix2.T @ atoms))
        # ratoms[0, :] += center[0]
        # ratoms[1, :] += center[1]
        # ratoms[2, :] += center[2]
        # newx = list(ratoms[0, :])
        # newy = list(ratoms[1, :])
        # newz = list(ratoms[2, :])

        # for i in range(n_atoms):
        #     atom = xp.array([newx[i], newy[i], newz[i]])
        #     catom = rotate_point(shift_point(atom, shift), matrix2_t)
        #     atom_surface = shift_point(rotate_point(catom, matrix1_t), center)
        #     newx[i] = atom_surface[0]
        #     newy[i] = atom_surface[1]
        #     newz[i] = atom_surface[2]

        # add positions to larger array that describes atom coordinates on the full surface
        membrane_x += list(ratoms[0, :])
        membrane_y += list(ratoms[1, :])
        membrane_z += list(ratoms[2, :])
        membrane_e += list(newe)
        membrane_b += list(newb)
        membrane_o += list(newo)

        atom_count += len(newe)
        print(f'Number of atoms added {len(newe)} to a total count of {atom_count}.')

    structure = (membrane_x, membrane_y, membrane_z, membrane_e, membrane_b, membrane_o)
    # pass directly to iasa_integration
    potential = iasa_integration('', voxel_size, solvent_exclusion=True, V_sol=solvent,
                                 absorption_contrast=True, voltage=voltage, density=physics.PROTEIN_DENSITY,
                                 molecular_weight=physics.PROTEIN_MW, structure_tuple=structure)
    return potential


if __name__ == '__main__':

    import os
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.tompy.io import write
    from pytom.simulation.support import reduce_resolution
    from pytom.tompy.transform import resize

    # syntax is ScriptOption([short, long], description, requires argument, is optional)
    options = [ScriptOption(['-s', '--size_factor'], '1 corresponds to a vesicle which has an average diameter of'
                                                     '45 nm across.', True, False),
               ScriptOption(['-d', '--destination'], 'Folder where output should be stored.', True, False),
               ScriptOption(['-m', '--membrane_pdb'], 'Membrane file, default '
                                                      '/data2/mchaillet/structures/pdb/lipid/dppc128_dehydrated.pdb'
                            , True, True),
               ScriptOption(['-h', '--help'], 'Help.', False, True)]

    helper = ScriptHelper(sys.argv[0].split('/')[-1], description='Create a template from the specified structure file '
                                                                  'with options for ctf correction and filtering. \n'
                                                                  'Script has dependencies on pytom and chimera.',
                          authors='Marten Chaillet', options=options)
    if len(sys.argv) == 2:
        print(helper)
        sys.exit()
    try:
        size_factor, folder, input_membrane, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()

    if help:
        print(helper)
        sys.exit()

    if size_factor: size_factor = float(size_factor)
    if input_membrane: pdb = input_membrane
    else: pdb = '/data2/mchaillet/structures/pdb/lipid/dppc128_dehydrated.pdb'

    # automatically scale these points
    N = int(100 * size_factor**2.2)  # number of points
    a, b, c = (x*size_factor for x in (xp.random.randint(180, 280), xp.random.randint(180, 280),
                                       xp.random.randint(180, 280)))
    alpha = 2000 * size_factor
    voxel = 5 # A

    solvent = 4.5301
    voltage = 300E3

    # generate an ellipsoid and triangulate it
    print('ellipsoid parameters: ' , a,b,c)
    points = sample_points_ellipsoid(N, a=a, b=b, c=c, evenly=True, maxiter=50000, factor=0.1)
    surface = triangulate(points, alpha)

    # fill the triangles with lipid molecules and calculate potential for it
    volume = membrane_potential(surface, voxel, pdb, solvent, voltage)
    real = volume[0]
    imag = volume[1]

    name = 'bilayer'
    size = f'{a*2/10:.0f}x{b*2/10:.0f}x{c*2/10:.0f}nm' # double the values of the ellipsoid radii for actual size

    real_fil = reduce_resolution(real, voxel, 2 * voxel)
    imag_fil = reduce_resolution(imag, voxel, 2 * voxel)

    write(os.path.join(folder, f'{name}_{voxel:.2f}A_{size}_solvent-4.530V_real.mrc'), real)
    write(os.path.join(folder, f'{name}_{voxel:.2f}A_{size}_solvent-4.530V_imag_300V.mrc'), imag)

    binning = 2

    real_bin = resize(reduce_resolution(real, voxel, binning * voxel * 2), 1/binning, interpolation='Spline')
    imag_bin = resize(reduce_resolution(imag, voxel, binning * voxel * 2), 1/binning, interpolation='Spline')

    write(os.path.join(folder, f'{name}_{voxel*binning:.2f}A_{size}_solvent-4.530V_real.mrc'), real_bin)
    write(os.path.join(folder, f'{name}_{voxel*binning:.2f}A_{size}_solvent-4.530V_imag_300V.mrc'), imag_bin)

