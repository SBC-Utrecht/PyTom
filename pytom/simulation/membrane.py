"""
Membrane geometrical structures

Depedencies: pyvista

Author: Marten Chaillet
"""

# essential
import numpy as np
import pytom.simulation.physics as physics
from numba import jit


class Vector:
    # Class can be used as both a 3d coordinate, and a vector
    def __init__(self, coordinates):
        assert len(coordinates) == 3, 'Invalid axis list for a 3d vector, input does not contain 3 coordinates.'
        self._axis = np.array(coordinates)
        self._zero_vector = np.all(self._axis==0)

    def get(self):
        return self._axis

    def show(self):
        print(self._axis)

    def copy(self):
        return Vector(self.get())

    def inverse(self):
        self._axis *= -1

    def cross(self, other):
        return Vector([self._axis[1] * other._axis[2] - self._axis[2] * other._axis[1],
                       self._axis[2] * other._axis[0] - self._axis[0] * other._axis[2],
                       self._axis[0] * other._axis[1] - self._axis[1] * other._axis[0]])

    def dot(self, other):
        # return the dot product of vectors v1 and v2, of form (x,y,z)
        # dot product of two vectors is zero if they are perpendicular
        return self._axis[0] * other._axis[0] + self._axis[1] * other._axis[1] + self._axis[2] * other._axis[2]

    def magnitude(self):
        # calculate the magnitude (length) of vector p
        return np.sqrt(np.sum(self._axis ** 2))

    def normalize(self):
        if not self._zero_vector:
            self._axis = self._axis / self.magnitude()

    def angle(self, other, degrees=False):
        # returns angle in radians
        if self._zero_vector or other._zero_vector:
            angle = 0
        else:
            angle = np.arccos(self.dot(other) / (self.magnitude() * other.magnitude()))
        if degrees:
            return angle * 180 / np.pi
        else:
            return angle

    def rotate(self, rotation_matrix):
        self._axis = np.dot(self._axis, rotation_matrix)

    def _get_orthogonal_unit_vector(self):
        # A vector orthogonal to (a, b, c) is (-b, a, 0), or (-c, 0, a) or (0, -c, b).
        if self._zero_vector:
            return Vector([1, 0, 0])  # zero vector is orthogonal to any vector
        else:
            if self._axis[2] != 0:
                x, y = 1, 1
                z = (- 1 / self._axis[2]) * (x * self._axis[0] + y * self._axis[1])
            elif self._axis[1] != 0:
                x, z = 1, 1
                y = (- 1 / self._axis[1]) * (x * self._axis[0] + z * self._axis[2])
            else:
                y, z = 1, 1
                x = (- 1 / self._axis[0]) * (y * self._axis[1] + z * self._axis[2])
            orth = Vector([x, y, z])
            orth.normalize()
            np.testing.assert_allclose(self.dot(orth), 0, atol=1e-7, err_msg='something went wrong in finding ' \
                                                                             'perpendicular vector')
            return orth

    def get_rotation(self, other):
        """
        Get rotation to rotate other onto self. Take the transpose to rotate self onto other.
        :param other:
        :type other:
        :return:
        :rtype:
        """
        if self._zero_vector or other._zero_vector:
            return np.identity(3)

        nself, nother = self.copy(), other.copy()
        nself.normalize()
        nother.normalize()

        if nself.dot(nother) > 0.99999:  # if the vectors are parallel
            return np.identity(3)  # return identity
        elif nself.dot(nother) < -0.99999:  # if the vectors are opposite
            axis = nself._get_orthogonal_unit_vector()  # return 180 along whatever axis
            angle = np.pi  # 180 degrees rotation around perpendicular vector
        else:
            axis = nself.cross(nother)
            axis.normalize()
            angle = nself.angle(nother)

        x, y, z = axis.get()
        c = np.cos(angle)
        s = np.sin(angle)
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
        return np.array([[m00, m01, m02], [m10, m11, m12], [m20, m21, m22]])


def z_axis_rotation_matrix(angle):
    m00 = np.cos(angle*np.pi/180)
    m01 = - np.sin(angle*np.pi/180)
    m10 = np.sin(angle*np.pi/180)
    m11 = np.cos(angle*np.pi/180)
    return np.array([[m00,m01,0], [m10,m11,0], [0,0,1]])


def distance_nonsqrt(p1, p2):
    return sum((p1 - p2) ** 2)


def distance(p1, p2):
    dd = (p1 - p2)**2
    return np.sqrt(sum(dd))


def get_point_ellipsoid(a, b, c, theta, phi):
    sinTheta = np.sin(theta)
    cosTheta = np.cos(theta)
    sinPhi = np.sin(phi)
    cosPhi = np.cos(phi)
    # rx = a * sinPhi * cosTheta
    # ry = b * sinPhi * sinTheta
    # rz = c * cosPhi
    rx = a * cosPhi * cosTheta
    ry = b * cosPhi * sinTheta
    rz = c * sinPhi
    return np.array([rx, ry, rz])


def random_point_ellipsoid(a, b, c):
    # a,b, and c are paremeters of the ellipsoid
    # generating random (x,y,z) points on ellipsoid
    u = np.random.rand()
    v = np.random.rand()
    theta = u * 2.0 * np.pi
    phi = np.arccos(2.0 * v - 1.0) - np.pi / 2
    # phi = u * 2.0 * np.pi
    # theta = np.arccos(2.0 * v - 1.0)
    return get_point_ellipsoid(a, b, c, theta, phi)


@jit(nopython=True)
def get_root_ellipse(r0, z0, z1, g, maxiter=20):
    # use bisection method to find root ellipse
    n0 = r0 * z0
    s0, s1 = z1 - 1, 0 if g < 0 else np.sqrt(n0 ** 2 + z1 ** 2) - 1
    s = 0
    for _ in range(maxiter):
        s = (s0 + s1) / 2
        if s == s0 or s == s1:
            break
        ratio0, ratio1 = n0 / (s + r0), z1 / (s + 1)
        g = ratio0 ** 2 + ratio1 ** 2 - 1
        if g > 0:
            s0 = s
        elif g < 0:
            s1 = s
        else:
            break
    return s


@jit(nopython=True)
def get_root_ellipsoid(r0, r1, z0, z1, z2, g, maxiter=20):
    # use bisection method to find root ellipsoid
    n0, n1 = r0 * z0, r1 * z1
    s0, s1 = z2 - 1, 0 if g < 0 else np.sqrt(n0 ** 2 + n1 ** 2 + z2 ** 2) - 1
    s = 0
    for _ in range(maxiter):
        s = (s0 + s1) / 2
        if s == s0 or s == s1:
            break
        ratio0, ratio1, ratio2 = n0 / (s + r0), n1 / (s + r1), z2 / (s + 1)
        g = ratio0 ** 2 + ratio1 ** 2 + ratio2 ** 2 - 1
        if g > 0:
            s0 = s
        elif g < 0:
            s1 = s
        else:
            break
    return s


@jit(nopython=True)
def distance_point_ellipse_quadrant(e0, e1, y0, y1, maxiter=20, epsilon=0):
    """
    from https://www.geometrictools.com/Documentation/DistancePointEllipseEllipsoid.pdf
    e0 >= e1 > 0, y0 >= 0, y1 >= 0  (Y is in the first quadrant)

    :param e0:
    :type e0:
    :param e1:
    :type e1:
    :param y0:
    :type y0:
    :param y1:
    :type y1:
    :return:
    :rtype:
    """
    if y1 > epsilon:
        if y0 > epsilon:
            z0, z1 = y0 / e0, y1 / e1
            g = z0 ** 2 + z1 ** 2 - 1
            if g != 0:
                r0 = (e0 / e1) ** 2
                sbar = get_root_ellipse(r0, z0, z1, g, maxiter=maxiter)
                x0, x1 = r0 * y0 / (sbar + r0), y1 / (sbar + 1)
                distance = np.sqrt((x0 - y0) ** 2 + (x1 - y1) ** 2)
            else:
                x0, x1, distance = y0, y1, 0
        else:  # y0 == 0
            x0, x1, distance = 0, e1, abs(y1 - e1)
    else:  # y1 == 0
        numer0, denom0 = e0 * y0, e0 ** 2 - e1 ** 2
        if numer0 < denom0:
            xde0 = numer0 / denom0
            x0, x1 = e0 * xde0, e1 * np.sqrt(1 - xde0 ** 2)
            distance = np.sqrt((x0 - y0) ** 2 + x1 ** 2)
        else:
            x0, x1, distance = e0, 0, abs(y0 - e0)
    return x0, x1, distance  # return point, distance


@jit(nopython=True)
def distance_point_ellipsoid_octant(e0, e1, e2, y0, y1, y2, maxiter=20, epsilon=0):
    """
    from https://www.geometrictools.com/Documentation/DistancePointEllipseEllipsoid.pdf
    e0 >= e1 >= e2 > 0, y0 >= 0, y1 >= 0, y2 >= 0  (Y is in the first octant)

    :param e0:
    :type e0:
    :param e1:
    :type e1:
    :param e2:
    :type e2:
    :param y0:
    :type y0:
    :param y1:
    :type y1:
    :param y2:
    :type y2:
    :return:
    :rtype:
    """
    if y2 > epsilon:
        if y1 > epsilon:
            if y0 > epsilon:
                z0, z1, z2 = y0 / e0, y1 / e1, y2 / e2
                g = z0 ** 2 + z1 ** 2 + z2 ** 2 - 1
                if g != 0:
                    r0, r1 = (e0 / e2) ** 2, (e1 / e2) ** 2
                    sbar = get_root_ellipsoid(r0, r1, z0, z1, z2, g, maxiter=maxiter)
                    x0, x1, x2 = r0 * y0 / (sbar + r0), r1 * y1 / (sbar + r1), y2 / (sbar + 1)
                    distance = np.sqrt((x0 - y0) ** 2 + (x1 - y1) ** 2 + (x2 - y2) ** 2)
                else:
                    x0, x1, x2, distance = y0, y1, y2, 0
            else:  # y0 == 0
                x0 = 0
                x1, x2, distance = distance_point_ellipse_quadrant(e1, e2, y1, y2)
        else:  # y1 == 0
            if y0 > epsilon:
                x1 = 0
                x0, x2, distance = distance_point_ellipse_quadrant(e0, e2, y0, y2)
            else:  # y0 == 0
                x0, x1, x2, distance = 0, 0, e2, abs(y2 - e2)
    else:  # y2 == 0
        denom0, denom1 = e0 ** 2 - e2 ** 2, e1 ** 2 - e2 ** 2
        numer0, numer1 = e0 * y0, e1 * y1
        computed = False
        if numer0 < denom0 and numer1 < denom1:
            xde0, xde1 = numer0 / denom0 , numer1 / denom1
            xde0sqr, xde1sqr = xde0 ** 2, xde1 ** 2
            discr = 1 - xde0sqr - xde1sqr
            if discr > 0:
                x0, x1, x2 = e0 * xde0, e1 * xde1, e2 * np.sqrt(discr)
                distance = np.sqrt((x0 - y0) ** 2 + (x1 - y1) ** 2 + x2 ** 2)
                computed = True
        if not computed:
            x2 = 0
            x0, x1, distance = distance_point_ellipse_quadrant(e0, e1, y0, y1)
    return x0, x1, x2, distance


def place_back_to_ellipsoid(point, a, b, c, maxiter=20):
    # This is not a problem as we can anyway rotate the ellipsoid afterwards to place them in simulations,
    # so the directions of a, b and c do not matter.
    assert a >= b >= c > 0, "finding point currently only possible on ellipsoid with a >= b >= c > 0"

    # find rotation to move point to the first quadrant
    v1, v2 = Vector(point), Vector(abs(point))
    octant_rotation_matrix = v1.get_rotation(v2)  # this is the rotation of v2 onto v1

    # find rotation of ellipsoid parameters so that a >= b >= c > 0 holds
    x, y, z, _ = distance_point_ellipsoid_octant(a, b, c, *v2.get(), maxiter=maxiter)
    intersection = Vector([x, y, z])
    intersection.rotate(octant_rotation_matrix)
    return intersection.get()


def test_place_back_to_ellipsoid(size=11, a=10, b=3, c=1, iterations=20):
    import matplotlib
    matplotlib.use('Qt5Agg')
    import matplotlib.pyplot as plt
    # xyz points between -a/b/c and a/b/c
    xx, yy, zz = map(int, [size*a, size*b, size*c])
    distances = np.zeros((xx, yy, zz))
    for i, x in enumerate(np.linspace(-a, a, xx)):
        for j, y in enumerate(np.linspace(-b, b, yy)):
            for k, z in enumerate(np.linspace(-c, c, zz)):
                point = np.array([x, y, z])
                ellipsoid_point = place_back_to_ellipsoid(point, a, b, c, maxiter=iterations)
                distances[i,j,k] = distance(point, ellipsoid_point)
    # visualize
    slice_x = distances[xx//2,:,:]
    slice_y = distances[:,yy//2,:]
    slice_z = distances[:,:,zz//2]
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
        return np.unravel_index(self.matrix.argmin(), self.matrix.shape)


def equilibrate_ellipsoid(points, a=2, b=3, c=4, maxiter=10000, factor=0.01, display=False):

    dmatrix = DistanceMatrix(points)

    for x in range(maxiter):
        if x % 1000 == 0: print(f'equilibrator iteration {x}')
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

        if display:
            import matplotlib
            matplotlib.use('Qt5Agg')
            import matplotlib.pyplot as plt
            print(p1, p2)
            print(p1_new, p2_new)
            print(points[minp1], points[minp2])
            plt.close()
            display_points_3d(points)

    return points


def sample_points_ellipsoid(number, a=2, b=3, c=4, evenly=True, maxiter=10000, factor=0.01, display=False):
    points = random_point_ellipsoid(a, b, c)
    for i in range(number-1):
        points = np.vstack((points, random_point_ellipsoid(a, b, c)))
    if evenly:
        return equilibrate_ellipsoid(points, a, b, c, maxiter=maxiter, factor=factor, display=display)
    else:
        return points


def get_point_ellipse(a, b, theta):
    rx = a * np.cos(theta)
    ry = b * np.sin(theta)
    return np.array([rx, ry])


def random_point_ellipse(a, b):
    u = np.random.rand()
    theta = u * 2.0 * np.pi
    return get_point_ellipse(a, b, theta)


def place_back_to_ellipse(point, a, b):

    from scipy.optimize import newton

    def f(theta, a, b, x, y):
        costheta = np.cos(theta)
        sintheta = np.sin(theta)
        return (a**2 - b**2) * costheta * sintheta - x * a * sintheta + y * b * costheta

    def fprime(theta, a, b, x, y):
        costheta = np.cos(theta)
        sintheta = np.sin(theta)
        return (a ** 2 - b ** 2) * (costheta**2 - sintheta**2) - x * a * costheta - y * b * sintheta

    x, y = point
    theta0 = np.arctan2(a*y, b*x)
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
            import matplotlib
            matplotlib.use('Qt5Agg')
            import matplotlib.pyplot as plt
            plt.close()
            display_points_2d(points)

    return points


def sample_points_ellipse(number, a=2, b=3, evenly=True, maxiter=1000, factor=0.01, display=False):
    points = random_point_ellipse(a, b)
    for i in range(number-1):
        points = np.vstack((points, random_point_ellipse(a, b)))
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
    new_point = np.matmul(rotation_matrix, point.reshape(-1, 1))
    return new_point.reshape(1, -1).squeeze()


def rotate_triangle(triangle, matrix):
    return np.vstack((rotate_point(triangle[i], matrix) for i in range(3)))


def shift_point(point, shift):
    return point + shift


def shift_triangle(triangle, shift):
    rtriangle = triangle
    rtriangle[0] = shift_point(rtriangle[0], shift)
    rtriangle[1] = shift_point(rtriangle[1], shift)
    rtriangle[2] = shift_point(rtriangle[2], shift)
    return rtriangle


def display_points_2d(points):
    import matplotlib
    matplotlib.use('Qt5Agg')
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(points[:, 0], points[:, 1])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.show()
    return


def display_points_3d(points, zlim=0):
    import matplotlib
    matplotlib.use('Qt5Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # <--- This is important for 3d plotting
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
    import matplotlib
    matplotlib.use('Qt5Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # <--- This is important for 3d plotting
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
    import matplotlib
    matplotlib.use('Qt5Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # <--- This is important for 3d plotting
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

    d1 = point_array_sign(point_array, v1, v2)
    d2 = point_array_sign(point_array, v2, v3)
    d3 = point_array_sign(point_array, v3, v1)

    # for numpy arrays
    has_neg = np.any([d1 < 0, d2 < 0, d3 < 0], axis=0)
    has_pos = np.any([d1 > 0, d2 > 0, d3 > 0], axis=0)

    return np.invert(np.all([has_neg, has_pos], axis=0))


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
    z_coordinates = list(np.array(z_coordinates) - sum(z_coordinates) / len(z_coordinates))
    # get the periodice boundary box from the pdb file
    x_bound, y_bound, z_bound = boundary_box_from_pdb(membrane_pdb)

    reference = Vector([.0, .0, 1.0])

    membrane_x, membrane_y, membrane_z, membrane_e, membrane_b, membrane_o = [], [], [], [], [], []

    atom_count = 0

    for icell in range(surface_mesh.n_cells):
        print(f'Triangle {icell+1} out of {surface_mesh.n_cells}.')

        triangle = surface_mesh.extract_cells(icell).points
        # print(triangle)
        normal = Vector(surface_mesh.cell_normals[icell])

        center = centroid(triangle)
        triangle1 = shift_triangle(triangle, -center)

        matrix1 = normal.get_rotation(reference)
        triangle2 = rotate_triangle(triangle1, matrix1)

        # add a random rotation of the triangle in the x-y plane to rotate the membrane's coordinates to other locations
        angle = np.random.uniform(0,360)
        matrix2 = z_axis_rotation_matrix(angle)
        triangle3 = rotate_triangle(triangle2, matrix2)

        # apply a shift to the points so that the coordinates are all above the origin
        shift = np.array([np.min(triangle3[:, 0]), np.min(triangle3[:, 1]), 0])
        triangle_sample = shift_triangle(triangle3, -shift)
        # print(np.round(triangle_sample, decimals=2))

        xmax = np.max(triangle_sample[:,0])
        ymax = np.max(triangle_sample[:,1])

        xext = int(np.ceil(xmax / x_bound))
        yext = int(np.ceil(ymax / y_bound))
        newx, newy, newz, newe, newb, newo = [], [], [], [], [], []
        for i in range(xext):
            for j in range(yext):
                newx += list(np.array(x_coordinates) + i*x_bound)
                newy += list(np.array(y_coordinates) + j*y_bound)
                newz += z_coordinates
                newe += elements
                newb += b_factors
                newo += occupancies
        # at this points our triangle fits onto the structure
        # now delete all the atoms that fall outside of the triangle
        i = 0
        while i < len(newx):
            coordinate = np.array([newx[i], newy[i]])
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

        atoms = np.array([newx, newy, newz])
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
        #     atom = np.array([newx[i], newy[i], newz[i]])
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
    potential = iasa_integration('', voxel_size, solvent_exclusion='masking', V_sol=solvent,
                                 absorption_contrast=True, voltage=voltage, density=physics.PROTEIN_DENSITY,
                                 molecular_weight=physics.PROTEIN_MW, structure_tuple=structure)
    return potential


def membrane_potential(surface_mesh, voxel_size, membrane_pdb, solvent, voltage, cores=1, gpu_id=None):
    """

    @param ellipsoid_mesh: this is a pyvista surface
    @param voxel_size:
    @param membrane_pdb:
    @param solvent_potential:
    @return:
    """

    from pytom.simulation.potential import read_structure, iasa_integration_parallel, iasa_integration_gpu
    from pytom.voltools.utils import translation_matrix
    from threadpoolctl import threadpool_info, threadpool_limits

    # READ THE STRUCTURE AND EXTEND IT ONLY ONCE
    # membrane pdb should have solvent deleted at this point
    x_coordinates, y_coordinates, z_coordinates, \
        elements, b_factors, occupancies = map(np.array, read_structure(membrane_pdb))
    z_coordinates -= z_coordinates.mean()
    n_atoms_box = len(elements)
    # get the periodice boundary box from the pdb file
    x_bound, y_bound, z_bound = boundary_box_from_pdb(membrane_pdb)

    reference = Vector([.0, .0, 1.0])

    membrane_x, membrane_y, membrane_z, membrane_e, membrane_b, membrane_o = [], [], [], [], [], []

    atom_count = 0

    user_apis = []
    for lib in threadpool_info():
        user_apis.append(lib['user_api'])
    if 'blas' in user_apis:
        print('BLAS in multithreading user apis, will try to limit number of threads for numpy.dot')

    for icell in range(surface_mesh.n_cells):


        triangle = surface_mesh.extract_cells(icell).points
        # print(triangle)
        normal = Vector(surface_mesh.cell_normals[icell])

        center = centroid(triangle)
        triangle1 = shift_triangle(triangle, -center)

        matrix1 = normal.get_rotation(reference)
        triangle2 = rotate_triangle(triangle1, matrix1)

        # add a random rotation of the triangle in the x-y plane to rotate the membrane's coordinates to other locations
        angle = np.random.uniform(0,360)
        matrix2 = z_axis_rotation_matrix(angle)
        triangle3 = rotate_triangle(triangle2, matrix2)

        # apply a shift to the points so that the coordinates are all above the origin
        shift = np.array([np.min(triangle3[:, 0]), np.min(triangle3[:, 1]), 0])
        triangle_sample = shift_triangle(triangle3, -shift)
        # print(np.round(triangle_sample, decimals=2))

        xmax = np.max(triangle_sample[:, 0])
        ymax = np.max(triangle_sample[:, 1])

        xext = int(np.ceil(xmax / x_bound))
        yext = int(np.ceil(ymax / y_bound))

        atoms = np.empty((xext * yext * n_atoms_box, 4))  # this should be (xext*yext, 3) for other multiplication
        # definition
        newe = np.zeros(xext * yext * n_atoms_box, dtype='<U1')
        newb, newo = (np.zeros(xext * yext * n_atoms_box),) * 2

        for i in range(xext):
            for j in range(yext):
                index = i * yext + j
                # add ones to 4th coordinate position for affine transformation
                atoms[index * n_atoms_box: (index + 1) * n_atoms_box, :] = np.array([x_coordinates + i * x_bound,
                                                                                     y_coordinates + j * y_bound,
                                                                                     z_coordinates,
                                                                                     np.ones(n_atoms_box)]).T
                newe[index * n_atoms_box: (index + 1) * n_atoms_box] = elements.copy()
                newb[index * n_atoms_box: (index + 1) * n_atoms_box] = b_factors.copy()
                newo[index * n_atoms_box: (index + 1) * n_atoms_box] = occupancies.copy()
                # newx += list(np.array(x_coordinates) + i*x_bound)
                # newy += list(np.array(y_coordinates) + j*y_bound)
                # newz += z_coordinates
                # newe += elements
                # newb += b_factors
                # newo += occupancies

        locs = point_array_in_triangle(atoms[:, :2], triangle_sample[:, :2])
        atoms = atoms[locs]
        newe = newe[locs]
        newb = newb[locs]
        newo = newo[locs]

        # at this points our triangle fits onto the structure
        # now delete all the atoms that fall outside of the triangle
        # i = 0
        # while i < len(newx):
        #     coordinate = np.array([newx[i], newy[i]])
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
        matrix1 = rotation_matrix_to_affine_matrix(matrix1)
        matrix2 = rotation_matrix_to_affine_matrix(matrix2)
        t2 = translation_matrix(translation=-shift)  # voltools implements matrices as inverse operations, here we want
        t1 = translation_matrix(translation=-center)  # forward transformation
        # TODO np.dot does the same as matmul in this situation but allows multithreading
        affine_matrix = np.matmul(np.matmul(t1, matrix1.T), np.matmul(matrix2.T, t2))  # transpose for the inverse
        # rotations

        if 'blas' in user_apis:
            with threadpool_limits(limits=cores, user_api='blas'):
                ratoms = np.dot(atoms, affine_matrix.T)
        else:
            ratoms = np.dot(atoms, affine_matrix.T)

        # ratoms = np.matmul(affine_matrix, atoms)

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
        #     atom = np.array([newx[i], newy[i], newz[i]])
        #     catom = rotate_point(shift_point(atom, shift), matrix2_t)
        #     atom_surface = shift_point(rotate_point(catom, matrix1_t), center)
        #     newx[i] = atom_surface[0]
        #     newy[i] = atom_surface[1]
        #     newz[i] = atom_surface[2]

        # TODO It might be much faster to write these coordinates to the end of a file. Because append will force new
        # TODO allocation of memory each time it is called.
        # add positions to larger array that describes atom coordinates on the full surface
        membrane_x += list(ratoms[:, 0])
        membrane_y += list(ratoms[:, 1])
        membrane_z += list(ratoms[:, 2])
        membrane_e += list(newe)
        membrane_b += list(newb)
        membrane_o += list(newo)

        atom_count += len(newe)

        if icell % 1000 == 0: print(f'At triangle {icell+1} out of {surface_mesh.n_cells}. Current atom count is '
                                    f'{atom_count}.')

    structure = (membrane_x, membrane_y, membrane_z, membrane_e, membrane_b, membrane_o)
    # pass directly to iasa_integration
    if gpu_id is not None:
        potential = iasa_integration_gpu('', voxel_size, solvent_exclusion='masking', V_sol=solvent,
                                              absorption_contrast=True, voltage=voltage,
                                              density=physics.PROTEIN_DENSITY,
                                              molecular_weight=physics.PROTEIN_MW, structure_tuple=structure,
                                              gpu_id=gpu_id)
    else:
        potential = iasa_integration_parallel('', voxel_size, solvent_exclusion='masking', V_sol=solvent,
                                     absorption_contrast=True, voltage=voltage, density=physics.PROTEIN_DENSITY,
                                     molecular_weight=physics.PROTEIN_MW, structure_tuple=structure, cores=cores)

    return potential


if __name__ == '__main__':
    import time
    import os
    import sys
    from pytom.tools.script_helper import ScriptHelper2, ScriptOption2
    from pytom.tools.parse_script_options import parse_script_options2
    from pytom.agnostic.io import write
    from pytom.simulation.support import reduce_resolution_fourier
    from pytom.agnostic.transform import resize

    start = time.time()

    # syntax is ScriptOption([short, long], description, requires argument, is optional)
    helper = ScriptHelper2(
        sys.argv[0].split('/')[-1],  # script name
        description='Generate a membrane structure for simulation',
        authors='Marten Chaillet',
        options = [ScriptOption2(['-r', '--radius_factor'], '1 corresponds to a vesicle which has an average diameter '
                                                          'of 45 nm across.', 'float', 'optional', 1.),
                   ScriptOption2(['-s', '--spacing'], 'Voxel spacing.', 'float', 'optional', 5),
                   ScriptOption2(['-d', '--destination'], 'Folder where output should be stored.', 'directory',
                                 'required'),
                   ScriptOption2(['-m', '--membrane_pdb'], 'Membrane file, default '
                                                           '/path/to/pdb/lipid/dppc128_dehydrated.pdb',
                                 'file', 'required'),
                   ScriptOption2(['-x', '--solvent'], 'Solvent background potential', 'float', 'optional',
                                 physics.V_WATER),
                   ScriptOption2(['-v', '--voltage'], 'Voltage for absorption contrast ', 'float', 'optional', 300),
                   ScriptOption2(['-c', '--cores'], 'Number of cores to use for numpy dot operations and later iasa '
                                                    'integration.', 'int', 'optional', 1),
                   ScriptOption2(['-g', '--gpuID'], 'GPU id to execute on', 'int', 'optional')])

    options = parse_script_options2(sys.argv[1:], helper)
    size_factor, voxel, folder, input_membrane, solvent, voltage, cores, gpuID = options
    voltage *= 1E3

    # automatically scale these points
    N = int(100 * size_factor**2.2)  # number of points
    a, b, c = sorted((x*size_factor for x in (np.random.randint(180, 280), np.random.randint(180, 280),
                                       np.random.randint(180, 280))), reverse=True)
    alpha = 2000 * size_factor
    voxel = 5

    # generate an ellipsoid and triangulate it
    print('Ellipsoid parameters: ' , a, b, c)
    points = sample_points_ellipsoid(N, a=a, b=b, c=c, evenly=True, maxiter=50000, factor=0.1)
    surface = triangulate(points, alpha)

    # fill the triangles with lipid molecules and calculate potential for it
    volume = membrane_potential(surface, voxel, input_membrane, solvent, voltage, cores=cores, gpu_id=gpuID)

    name = 'bilayer'
    size = f'{a*2/10:.0f}x{b*2/10:.0f}x{c*2/10:.0f}nm'  # double the values of the ellipsoid radii for actual size

    real_fil = reduce_resolution_fourier(volume.real, voxel, 2 * voxel)
    imag_fil = reduce_resolution_fourier(volume.imag, voxel, 2 * voxel)

    write(os.path.join(folder, f'{name}_{size}_{voxel:.2f}A_solvent-4.530V_real.mrc'), real_fil)
    write(os.path.join(folder, f'{name}_{size}_{voxel:.2f}A_solvent-4.530V_imag_300V.mrc'), imag_fil)

    # binning = 2
    #
    # real_bin = resize(reduce_resolution_fourier(volume.real, voxel, binning * voxel * 2), 1/binning,
    #                   interpolation='Spline')
    # imag_bin = resize(reduce_resolution_fourier(volume.imag, voxel, binning * voxel * 2), 1/binning,
    #                   interpolation='Spline')
    #
    # write(os.path.join(folder, f'{name}_{voxel*binning:.2f}A_{size}_solvent-4.530V_real.mrc'), real_bin)
    # write(os.path.join(folder, f'{name}_{voxel*binning:.2f}A_{size}_solvent-4.530V_imag_300V.mrc'), imag_bin)

    end = time.time()

    print('\n Time elapsed: ', end-start, '\n')

