#! /usr/bin/env python
# -*- coding: utf-8 -*-

from pytom.gpu.initialize import xp

from pytom.tompy.filter import gaussian3d
from numpy.random import standard_normal


def create_sphere(size, radius=-1, sigma=0, num_sigma=2, center=None, gpu=False):
    """Create a 3D sphere volume.

    @param size: size of the resulting volume.
    @param radius: radius of the sphere inside the volume.
    @param sigma: sigma of the Gaussian.
    @param center: center of the sphere.

    @return: sphere inside a volume.
    """


    if len(size) == 1:
        size = (size, size, size)
    assert len(size) == 3

    if center is None:
        center = [size[0]//2, size[1]//2, size[2]//2]
    if radius == -1:
        radius = xp.min(size)//2

    sphere = xp.zeros(size)
    [x,y,z] = xp.mgrid[0:size[0], 0:size[1], 0:size[2]]
    r = xp.sqrt((x-center[0])**2+(y-center[1])**2+(z-center[2])**2)
    sphere[r<=radius] = 1

    if sigma > 0:
        ind = xp.logical_and(r>radius, r<radius+num_sigma*sigma)
        sphere[ind] = xp.exp(-((r[ind] - radius)/sigma)**2/2)

    return sphere

def create_circle(size, radius=-1, sigma=0, num_sigma=3, center=None):
    """Create a 3D sphere volume.

    @param size: size of the resulting volume.
    @param radius: radius of the sphere inside the volume.
    @param sigma: sigma of the Gaussian.
    @param center: center of the sphere.

    @return: sphere inside a volume.
    """


    if len(size) == 1:
        size = (size, size)
    assert len(size) == 2

    if center is None:
        center = [size[0]//2, size[1]//2]
    if radius == -1:
        radius = xp.min(size)//2

    circle = xp.zeros(size)
    [x,y] = xp.mgrid[0:size[0], 0:size[1]]
    r = xp.sqrt((x-center[0])**2+(y-center[1])**2)
    circle[r<=radius] = 1

    if sigma > 0:
        ind = xp.logical_and(r>radius, r<=radius+num_sigma*sigma)
        circle[ind] = xp.exp(-((r[ind] - radius)/sigma)**2/2)

    return circle

def prepare_mask(v, threshold, smooth):
    """Prepare a mask according to the given volume.
    Everything above the given threshold will be set to 1.

    @param v: input volume.
    @param threshold: threshold.
    @param smooth: sigma of Gaussian.

    @return: mask.
    """
    from pytom.tompy.filter import gaussian3d
    ind = xp.where(v>threshold)
    mask = xp.zeros(v.shape)
    mask[ind] = 1
    
    return gaussian3d(mask, smooth)

def add_noise(data, snr=0.1, m=0):
    """Add gaussian noise to the given volume.
    @param data: input volume.
    @param snr: SNR of the added noise.
    @param m: mean of the added noise.

    @return The image with gaussian noise	
    """
    vs = xp.var(data)
    vn = vs/snr
    sd = xp.sqrt(vn)
    s = xp.ndarray(data.shape)
    noise = sd*standard_normal(s.shape)+m
    t = data + noise
    return t

def paste_in_center(volume, volume2, gpu=False):
    if 0:
        pass#raise Exception('pasteCenter not defined for gpu yet.')
    else:
        l,l2 = len(volume.shape), len(volume.shape)
        assert l == l2
        for i in range(l):
            assert volume.shape[i] <= volume2.shape[i]

        if len(volume.shape) == 3:
            sx,sy,sz = volume.shape
            SX, SY, SZ = volume2.shape
            volume2[SX//2-sx//2:SX//2+sx//2+sx%2,SY//2-sy//2:SY//2+sy//2+sy%2,SZ//2-sz//2:SZ//2+sz//2+sz%2 ] = volume
            return volume2

        if len(volume.shape) == 2:
            sx,sy = volume.shape
            SX, SY = volume2.shape
            volume2[SX//2-sx//2:SX//2+sx//2+sx%2,SY//2-sy//2:SY//2+sy//2+sy%2] = volume
            return volume2

def rotation_matrix_x(angle):
    """Return the 3x3 rotation matrix around x axis.

    @param angle: rotation angle around x axis (in degree).

    @return: rotation matrix.
    """
    angle = xp.deg2rad(angle)
    mtx = xp.matrix(xp.zeros((3,3)))
    mtx[1,1] = xp.cos(angle)
    mtx[2,1] = xp.sin(angle)
    mtx[2,2] = xp.cos(angle)
    mtx[1,2] = -xp.sin(angle)
    mtx[0,0] = 1

    return mtx

def rotation_matrix_y(angle):
    """Return the 3x3 rotation matrix around y axis.

    @param angle: rotation angle around y axis (in degree).

    @return: rotation matrix.
    """
    angle = xp.deg2rad(angle)
    mtx = xp.matrix(xp.zeros((3,3)))
    mtx[0,0] = xp.cos(angle)
    mtx[2,0] = -xp.sin(angle)
    mtx[2,2] = xp.cos(angle)
    mtx[0,2] = xp.sin(angle)
    mtx[1,1] = 1

    return mtx

def rotation_matrix_z(angle):
    """Return the 3x3 rotation matrix around z axis.

    @param angle: rotation angle around z axis (in degree).

    @return: rotation matrix.
    """
    angle = xp.deg2rad(angle)
    mtx = xp.matrix(xp.zeros((3,3)))
    mtx[0,0] = xp.cos(angle)
    mtx[1,0] = xp.sin(angle)
    mtx[1,1] = xp.cos(angle)
    mtx[0,1] = -xp.sin(angle)
    mtx[2,2] = 1

    return mtx

def rotation_matrix_zxz(angle):
    """Return the 3x3 rotation matrix of an Euler angle in ZXZ convention.
    Note the order of the specified angle should be [Phi, Psi, Theta], or [Z1, Z2, X] in Pytom format.

    @param angle: list of [Phi, Psi, Theta] in degree.

    @return: rotation matrix.
    """
    assert len(angle) == 3

    z1 = angle[0]
    z2 = angle[1]
    x = angle[2]

    zm1 = rotation_matrix_z(z1)
    xm = rotation_matrix_x(x)
    zm2= rotation_matrix_z(z2)
    
    res = zm2 * (xm * zm1)

    return res

def rotation_distance(ang1, ang2):
    """Given two angles (lists), calculate the angular distance (degree).

    @param ang1: angle 1. [Phi, Psi, Theta], or [Z1, Z2, X] in Pytom format.
    @param ang2: angle 2. [Phi, Psi, Theta], or [Z1, Z2, X] in Pytom format.
    
    @return: rotation distance in degree.
    """
    mtx1 = rotation_matrix_zxz(ang1)
    mtx2 = rotation_matrix_zxz(ang2)
    res = xp.multiply(mtx1, mtx2) # elementwise multiplication
    trace = xp.sum(res)
    
    from math import pi, acos
    temp=0.5*(trace-1.0)
    if temp >= 1.0:
        return 0.0
    if temp <= -1.0:
        return 180
    return acos(temp)*180/pi

def euclidian_distance(pos1, pos2):
    return xp.linalg.norm(xp.array(pos1)-xp.array(pos2))

def volumesSameSize(v0, v1):
    if len(v0.shape) != len(v1.shape):
        return False
    return xp.array([abs(v0.shape[pos] - v1.shape[pos]) == 0 for pos in range(len(v0.shape))]).sum() == len(v0.shape)

def taper_edges(image, width):
    """
    taper edges of image (or volume) with cos function

    @param image: input image (or volume)
    @type image: ndarray
    @param width: width of edge
    @type width: int

    @return: image with smoothened edges, taper_mask
    @rtype: array-like

    @author: GvdS
    """

    dims = list(image.shape) + [0]
    val = xp.cos(xp.arange(1, width + 1) * xp.pi / (2. * (width)))
    taperX = xp.ones((dims[0]), dtype=xp.float32)
    taperY = xp.ones((dims[1]))
    taperX[:width] = val[::-1]
    taperX[-width:] = val
    taperY[:width] = val[::-1]
    taperY[-width:] = val
    if dims[2] > 1:
        taperZ = xp.ones((dims[2]))
        taperZ[:width] = val[::-1]
        taperZ[-width:] = val
        Z, X, Y = xp.meshgrid(taperX, taperY, taperZ)
        taper_mask = X * (X < Y) * (X < Z) + Y * (Y <= X) * (Y < Z) + Z * (Z <= Y) * (Z <= X)
    else:
        X, Y = xp.meshgrid(taperX, taperY)
        taper_mask = X * (X < Y) + Y * (Y <= X)

    return image * taper_mask, taper_mask

def resize(*args, **kwargs):
    pass

def determineRotationCenter(particle, binning):
    """
    determineRotationCenter:
    @param particle: The particle
    @type particle: Either L{pytom_volume.vol} or string specifying the particle file name
    @param binning: Binning factor
    @return: [centerX,centerY,centerZ]

    @author: GvdS
    """
    if particle.__class__ == str:
        from pytom.tompy.io import read
        particle = read(particle)

    centerX = particle.shape[0] / 2.0 * (1.0 / float(binning))
    centerY = particle.shape[1] / 2.0 * (1.0 / float(binning))
    centerZ = particle.shape[2] / 2.0 * (1.0 / float(binning))
    #
    # if binning > 1:
    #   centerX = centerX - 0.25*(binning-1)
    #    centerY = centerY - 0.25*(binning-1)
    #    centerZ = centerZ - 0.25*(binning-1)

    return [centerX, centerY, centerZ]

def invert_WedgeSum( invol, r_max=None, lowlimit=0., lowval=0.):
    """
    invert wedge sum - avoid division by zero and boost of high frequencies

    @param invol: input volume
    @type invol: L{pytom_volume.vol} or L{pytom_volume.vol_comp}
    @param r_max: radius
    @type r_max: L{int}
    @param lowlimit: lower limit - all values below this value that lie in the specified radius will be replaced \
                by lowval
    @type lowlimit: L{float}
    @param lowval: replacement value
    @type lowval: L{float}

    @author: FF
    """
    from math import sqrt
    if not r_max:
        r_max=invol.shape[1]//2-1

    dx,dy,dz= invol.shape
    if dz != dx:
        X, Y, Z = meshgrid(arange(-dx // 2, dx // 2 + dx % 2), arange(-dy // 2, dy // 2 + dy % 2), arange(0,dz//2+1))
        invol = xp.fft.fftshift(invol,axes=(0,1))
    else:
        X, Y, Z = meshgrid(arange(-dx // 2, dx // 2 + dx % 2), arange(-dy // 2, dy // 2 + dy % 2), arange(-dz // 2, dz // 2 + dz % 2))
    R = sqrt(X ** 2 + Y ** 2 + Z**2).astype(int)

    invol_out= invol.copy()
    invol_out[invol < lowlimit] =  lowval
    invol_out = 1/invol_out
    invol_out[R >= rmax] = 0

    if dx != dz:
        invol_out = xp.fft.fftshift(invol_out, axes=(0,1))

    return invol_out