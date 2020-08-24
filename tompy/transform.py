#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import scipy

def rotate_axis(data, angle, axis='z'):
    """Rotate the volume around certain axis.

    @param data: input volume.
    @param angle: angle to rotate in degrees (counter-clockwise).
    @param axis: axis to rotate around.

    @return: rotated volume.
    """
    if axis == 'z':
        a = (0, 1)
    elif axis == 'y':
        a = (0, 2)
        angle = -angle
    elif axis == 'x':
        a = (1, 2)
    else:
        raise ValueError("Invalid input argument! Choose between 'x', 'y' or 'z'.")

    from scipy.ndimage.interpolation import rotate
    res = rotate(data, angle, a, reshape=False, mode='constant')
    return res


def rotate3d(data, phi=0, psi=0, the=0, center=None, order=2):
    """Rotate a 3D data using ZXZ convention (phi: z1, the: x, psi: z2).

    @param data: data to be rotated.
    @param phi: 1st rotate around Z axis, in degree.
    @param psi: 3rd rotate around Z axis, in degree.
    @param the: 2nd rotate around X axis, in degree.
    @param center: rotation center.

    @return: the data after rotation.
    """
    # Figure out the rotation center
    if center is None:
        cx = data.shape[0] // 2
        cy = data.shape[1] // 2
        cz = data.shape[2] // 2
    else:
        assert len(center) == 3
        (cx, cy, cz) = center

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



    from scipy import mgrid
    grid = mgrid[-cx:data.shape[0]-cx, -cy:data.shape[1]-cy, -cz:data.shape[2]-cz]
    temp = grid.reshape((3, grid.size // 3))
    temp = np.dot(Inv_R, temp)
    grid = np.reshape(temp, grid.shape)
    grid[0] += cx
    grid[1] += cy
    grid[2] += cz

    # Interpolation
    from scipy.ndimage import map_coordinates
    d = map_coordinates(data, grid, order=order)

    return d


def translate3d(data, dx=0, dy=0, dz=0, order=2):
    """Translate the data in real space.

    @param data: data to be shifted.
    @param dx: translation along x-axis.
    @param dy: translation along y-axis.
    @param dz: translation along z-axis.
    
    @return: the data after translation.
    """
    if dx == 0 and dy == 0 and dz == 0:
        return data

    # from scipy.ndimage.interpolation import shift
    # res = shift(data, [dx, dy, dz])
    # return res
    from scipy import mgrid
    grid = mgrid[0.:data.shape[0], 0.:data.shape[1], 0.:data.shape[2]]
    grid[0] -= dx
    grid[1] -= dy
    grid[2] -= dz
    from scipy.ndimage import map_coordinates
    d = map_coordinates(data, grid, order=order)

    return d


def translate3d_f(data, dx=0, dy=0, dz=0):
    """Translate the data using Fourier shift theorem.
    """
    if dx == 0 and dy == 0 and dz == 0:
        return data

    sx = data.shape[0]
    sy = data.shape[1]
    sz = data.shape[2]

    xx, yy, zz = np.indices((sx, sy, sz/2+1))

    xx[np.where(xx >= sx/2)] -= sx
    yy[np.where(yy >= sy/2)] -= sy

    # Fourier shift theorem
    shift = np.exp(-2j*np.pi/sx*xx*dx) * np.exp(-2j*np.pi/sy*yy*dy) * np.exp(-2j*np.pi/sz*zz*dz)

    fdata = rfft(data)

    res = irfft(fdata * shift, data.shape)

    return res


def transform3d(data, m, order=2):
    """Transform 3D data using 3x3 transformation matrix.

    @param data: data.
    @param m: 3x3 transformation matrix.

    @return: data after transformation.
    """
    from scipy import mgrid
    grid = mgrid[0.:data.shape[0], 0.:data.shape[1], 0.:data.shape[2]]
    temp = grid.reshape((3, grid.size // 3))
    temp = np.dot(m, temp)
    grid = np.reshape(temp, grid.shape)

    from scipy.ndimage import map_coordinates
    d = map_coordinates(data, grid, order=order)
    return d


def resize(data, x, y, z):
    """Resize the data in real space.

    @param data: input data.
    @param x: resized dimension x.
    @param y: resized dimension y.
    @param z: resized dimension z.

    @return: resized data.
    """
    s = data.shape
    from scipy import mgrid, array
    from scipy.ndimage import map_coordinates
    grid = mgrid[0:s[0]-1:x * 1j, 0:s[1]-1:y * 1j, 0:s[2]-1:z * 1j]
    d = map_coordinates(data, grid, order=2)
    return d


def cut_from_projection(proj, center, size):
    """Cut out a subregion out from a 2D projection.

    @param proj: 2D projection.
    @param center: cutting center.
    @param size: cutting size.

    @return: subprojection.
    """
    from scipy import mgrid
    from scipy.ndimage import map_coordinates

    if len(proj.shape) > 2 and proj.shape[2] > 1:
        raise Exception('We assume that projections come from a 3D object, thus your projection cannot be a 3D object itself')

    # adjusted to python3
    grid = mgrid[center[0]-size[0]//2:center[0]+size[0]-size[0]//2-1:size[0]*1j, 
                 center[1]-size[1]//2:center[1]+size[1]-size[1]//2-1:size[1]*1j]
    v = map_coordinates(proj.squeeze(), grid, order=2)
    return v

def projTiltY(vol, tiltangle, projSize=None, center=None):
    """
    project volume along z, tilted around y
    @param vol: volume
    @type vol: numpy array
    @param tiltangle
    @type tiltangle: C{float}
    @param center: (x,z) rotation center (default Nx//2, Nz//2, first index=0)
    @type center: L{numpy.array}
    @return: projection
    @rtype: numpy array
    """
    dims = vol.shape
    assert dims[2] == 1, "input volume must be 3-dim"
    if projSize is None:
        projSize = 2*[0]
        projSize[0] = dims[0]
        projSize[1] = dims[1]
    else:
        assert len(projSize) == 2, "projection Size must be 2D vector"
        assert projSize[1] <= dims[1], "projection y-Size must not exceed dimension of volume"
        assert projSize[0] <= dims[0], "projection x-Size must not exceed dimension of volume"

    if center is None:
        cx = dims[0] // 2
        cz = dims[2] // 2
    else:
        assert len(center) == 2, "center olume must be 2D vector"
        (cx, cz) = center

    # define rotation around y-axis - inverse because you start from the projection
    print('tiltangle: ',tiltangle)
    the = -float(tiltangle) * np.pi / 180.0
    sin_beta = np.sin(the)
    cos_beta = np.cos(the)

    # mapping
    rotmat = np.zeros((2, 2), dtype='float32')
    projmat = np.zeros((2, 2), dtype='float32')
    projrotmat = np.zeros((2, 2), dtype='float32')
    rotmat[0, 0] = cos_beta
    rotmat[1, 1] = cos_beta
    rotmat[0, 1] = -sin_beta
    rotmat[1, 0] = sin_beta
    projmat[0,0] = 1.
    projrotmat[0,0] = projmat[0,0]*rotmat[0, 0]
    projrotmat[0,1] = projmat[0,0]*rotmat[0, 1]

    from scipy import mgrid
    grid = mgrid[-cx:dims[0]-cx, -cz:dims[2]-cz]
    temp = grid.reshape((2, grid.size // 2))
    temp = np.dot(projrotmat, temp)
    grid = np.reshape(temp, grid.shape)
    grid[0] += cx
    projCoordsInXZSlice = grid[0]

    # split into 2D problems
    projgrid = mgrid[-cx:dims[0]-cx]
    # linear interpolation weights
    #i1 = projCoordsInXZSlice.round()
    i1 = projCoordsInXZSlice.astype(int)
    w1 = 1. - np.abs(projCoordsInXZSlice - i1)
    i2 = (1+projCoordsInXZSlice).astype(int)
    #w2 = np.abs(projCoordsInXZSlice - i1)
    w2 = 1-w1
    # make lists - can be re-used for each slice
    indexlist = []
    weightlist = []
    for ix in range(0,projSize[0]):
        itmp1 = np.where(i1 == ix)
        wtmp1 = w1[itmp1]
        itmp2 = np.where(i2 == ix)
        wtmp2 = w2[itmp2]
        itmpx = (np.concatenate((itmp1[0],itmp2[0])), np.concatenate((itmp1[1],itmp2[1])))
        indexlist.append(itmpx)
        weightlist.append( np.concatenate((w1[itmp1], w2[itmp2])))

    projection = np.zeros((projSize[0], projSize[1]), dtype='float32')
    projline = np.zeros( projSize[0], dtype='float32')
    # loop over y
    for iy in range(0,projSize[1]):
        projline[:] = 0.
        volslice = vol[:,iy,:]
        for ix in range(0,projSize[0]):
            itmp = indexlist[ix]
            projline[ix] = np.sum( volslice[itmp]*weightlist[ix] )
        projection[:,iy] = projline

    return projection



################################
#    Fourier Relevant Stuff    #
################################

def rfft(data):
    """Do 3D Fourier transformaion of data of real numbers.

    @param input data.

    @return: the data after transformation.
    """
    return np.fft.rfftn(data)


def irfft(data, s=None):
    """Do 3D Inverse Fourier transformaion to get back data of real numbers.

    @param data: input data.
    
    @return: the data after ifft (without the need to scale).
    """
    return np.fft.irfftn(data, s)


def fft(data):
    """Do 3D Fourier transformaion of data (real or complex numbers).

    @param input data.

    @return: the data after transformation.
    """
    return np.fft.fftn(data)


def ifft(data):
    """Do 3D Inverse Fourier transformaion to get back data (real or complex numbers).

    @param data: input data.
    
    @return: the data after ifft (without the need to scale).
    """
    return np.fft.ifftn(data)


def fftshift(data):
    return np.fft.fftshift(data)


def ifftshift(data):
    return np.fft.ifftshift(data)


def conv3d(data, kernel):
    """Do 3D convolution.

    @param: input data.
    @param kernel: the kernel to convolve.

    @return: data after convolution.
    """
    from scipy.ndimage.filters import convolve
    d = convolve(data, kernel)
    return d


def fourier_reduced2full(data, isodd=False):
    """Return an Hermitian symmetried data.
    """
    # get the original shape
    sx = data.shape[0]
    sy = data.shape[1]
    if isodd:
        sz = (data.shape[2]-1)*2+1
    else:
        sz = (data.shape[2]-1)*2

    res = np.zeros((sx, sy, sz), dtype=data.dtype)
    res[:, :, 0:data.shape[2]] = data

    # calculate the coodinates accordingly
    szz = sz - data.shape[2]
    x, y, z = np.indices((sx, sy, szz))
    ind = [np.mod(sx-x, sx), np.mod(sy-y, sy), szz-z]

    # do the complex conjugate of the second part
    res[:, :, data.shape[2]:] = np.ma.conjugate(data[tuple(ind)])

    return res


def fourier_full2reduced(data):
    return data[:,:,0:data.shape[2]//2+1]



def fourier_filter(data, fltr, human=True):
    if human:
        fltr = ifftshift(fltr)
        fltr = fourier_full2reduced(fltr)

    fd = rfft(data)
    res = irfft(fd * fltr, data.shape)

    return res


def resiz2e(volume, factor, interpolation='Fourier'):
    """
    resize volume in real or Fourier space
    @param volume: input volume
    @type volume: L{pytom_volume.vol}
    @param factor: a factor > 0. Factors < 1 will de-magnify the volume, factors > 1 will magnify.
    @type factor: L{float}
    @param interpolation: Can be 'Fourier' (default), 'Spline', 'Cubic' or 'Linear'
    @type interpolation: L{str}

    @return: The re-sized volume
    @rtype: L{pytom_volume.vol}
    @author: FF
    """

    if (interpolation == 'Spline') or (interpolation == 'Cubic') or (interpolation == 'Linear'):
        return scale(volume, factor, interpolation='Spline')
    else:
        fvol = xp.fft.rfftn(volume)
        newfvol = resizeFourier(fvol=fvol, factor=factor, isodd=volume.shape[2]%2)
        newvol = (xp.fft.irfftn(newfvol, s=[newfvol.shape[0],]*len(newfvol.shape)))

        return newvol, newfvol
