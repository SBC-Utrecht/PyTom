#! /usr/bin/env python
# -*- coding: utf-8 -*-
from pytom.gpu.initialize import xp, device
# Typing imports
from pytom.gpu.initialize import xpt

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
    import numpy as xp

    if center is None:
        cx = data.shape[0] // 2
        cy = data.shape[1] // 2
        cz = data.shape[2] // 2
    else:
        assert len(center) == 3
        (cx, cy, cz) = center

    if phi == 0 and psi==0 and the==0:
        return data

    # Transfer the angle to Euclidean
    phi = -float(phi) * xp.pi / 180.0
    the = -float(the) * xp.pi / 180.0
    psi = -float(psi) * xp.pi / 180.0
    sin_alpha = xp.sin(phi)
    cos_alpha = xp.cos(phi)
    sin_beta = xp.sin(the)
    cos_beta = xp.cos(the)
    sin_gamma = xp.sin(psi)
    cos_gamma = xp.cos(psi)

    # Calculate inverse rotation matrix
    Inv_R = xp.zeros((3, 3), dtype='float32')

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


    grid = xp.mgrid[-cx:data.shape[0]-cx, -cy:data.shape[1]-cy, -cz:data.shape[2]-cz]
    temp = grid.reshape((3, grid.size // 3))
    temp = xp.dot(Inv_R, temp)
    grid = xp.reshape(temp, grid.shape)
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

    xx, yy, zz = xp.indices((sx, sy, sz/2+1))

    xx[xp.where(xx >= sx/2)] -= sx
    yy[xp.where(yy >= sy/2)] -= sy

    # Fourier shift theorem
    shift = xp.exp(-2j*xp.pi/sx*xx*dx) * xp.exp(-2j*xp.pi/sy*yy*dy) * xp.exp(-2j*xp.pi/sz*zz*dz)

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
    temp = xp.dot(m, temp)
    grid = xp.reshape(temp, grid.shape)

    from scipy.ndimage import map_coordinates
    d = map_coordinates(data, grid, order=order)
    return d

def resize_scipy(data, x, y, z):
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

def cut_from_projection(proj, center, size, device=2):
    """Cut out a subregion out from a 2D projection.

    @param proj: 2D projection.
    @param center: cutting center.
    @param size: cutting size.

    @return: subprojection.
    """
    import pytom.voltools as vt
    if len(proj.shape) > 2 and proj.shape[2] > 1:
        raise Exception('We assume that projections come from a 3D object, thus your projection cannot be a 3D object itself')

    import numpy as np
    from scipy.ndimage import map_coordinates
    try:
        proj = proj.squeeze().get()
    except:
        proj = proj.squeeze()

    #v =
    # adjusted to python3
    grid = np.mgrid[center[0]-size[0]//2:center[0]+size[0]-size[0]//2-1:size[0]*1j,
                    center[1]-size[1]//2:center[1]+size[1]-size[1]//2-1:size[1]*1j]


    v = map_coordinates(proj, grid, order=2)

    return xp.array(v)

def scale(volume, factor, interpolation='Spline'):
    """
    scale: Scale a volume by a factor in REAL SPACE - see also resize function for more accurate operation in Fourier \
    space
    @param volume: input volume
    @type volume: L{pytom_volume.vol}
    @param factor: a factor > 0. Factors < 1 will de-magnify the volume, factors > 1 will magnify.
    @type factor: L{float}
    @param interpolation: Can be Spline (default), Cubic or Linear
    @type interpolation: L{str}

    @return: The scaled volume
    @author: Thomas Hrabe
    """

    from pytom.voltools import transform
    from pytom.agnostic.tools import paste_in_center
    if factor <= 0:
        raise RuntimeError('Scaling factor must be > 0!')

    interpolation_dict = {'Spline': 'filt_bspline',
                          'Linear': 'linear',
                          'Cubic': 'filt_bspline',
                          'filt_bspline': 'filt_bspline',
                          'linear':'linear'}
    interpolation = interpolation_dict[interpolation]

    volume = volume.squeeze()

    size_x = volume.shape[0]
    size_y = volume.shape[1]
    size_z = 1
    newSizeX = int(xp.floor(size_x * float(factor) + 0.5))
    newSizeY = int(xp.floor(size_y * float(factor) + 0.5))
    newSizeZ = 1

    if len(volume.shape) == 3:
        size_z = volume.shape[2]
        newSizeZ = int(xp.floor(size_z * factor + 0.5))
        scaleF = [1/factor, 1/factor, 1/factor]
    else:
        scaleF = [1/factor, 1/factor, 1]
        volume = xp.expand_dims(volume,2)

    if factor ==1:
        rescaledVolume = volume
    elif factor > 1:
        newVolume = xp.zeros((newSizeX, newSizeY, newSizeZ),dtype=volume.dtype)
        newVolume = paste_in_center(volume, newVolume)
        rescaledVolume = xp.zeros_like(newVolume)
        transform(newVolume, scale=scaleF, output=rescaledVolume, device=device, interpolation=interpolation)
    else:
        rescaledVolumeFull = xp.zeros_like(volume)
        transform(volume, scale=scaleF, output=rescaledVolumeFull, device=device, interpolation=interpolation)
        rescaledVolume = xp.zeros((newSizeX, newSizeY, newSizeZ), dtype=volume.dtype)
        rescaledVolume = paste_in_center(rescaledVolumeFull, rescaledVolume)

    return rescaledVolume

def resize(volume, factor, interpolation='Fourier'):
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
    import numpy
    ss = len(volume.shape)
    org_shape = volume.shape
    org_size = volume.size
    volume = volume.squeeze()

    if (interpolation == 'Spline') or (interpolation == 'Cubic') or (interpolation == 'Linear'):
        return scale(volume, factor, interpolation=interpolation)
    else:
        fvol = xp.fft.rfftn(volume)
        outsize= tuple((numpy.around(numpy.array(volume.shape)*factor,0)).astype(int))
        newfvol = resizeFourier(fvol=fvol, factor=factor, isodd=volume.shape[-1]%2)

        newvol = xp.fft.irfftn(newfvol, s=outsize)

        if ss != len(newvol.shape):
            newvol = xp.expand_dims(newvol,2)

        newvol *=  newvol.size/org_size

        return newvol

def resizeFourier(fvol, factor, isodd=False):
    """
    resize Fourier transformed by factor

    @param fvol: Fourier transformed of a volume  - reduced complex
    @type fvol: L{pytom_volume.vol_comp}
    @param factor:  a factor > 0. Factors < 1 will de-magnify the volume, factors > 1 will magnify.
    @type factor: float
    @return: resized Fourier volume (deruced complex)
    @rtype: L{pytom_volume.vol_comp}
    @author: FF
    """

    fvol = fvol.squeeze()
    if factor == 1: return fvol
    oldFNx = fvol.shape[0]
    oldFNy = fvol.shape[1]


    # new dims in real and Fourier space
    newNx = newFNx = int(xp.floor(oldFNx * factor + 0.5))
    newNy = newFNy = int(xp.floor(oldFNy * factor + 0.5 ))

    # check 3D images
    if len(fvol.shape) == 3:
        oldNz = int((fvol.shape[2] - 1)*2 + 1*isodd)
        newNz = int(xp.floor(oldNz * factor + 0.5))
        newFNz = newNz // 2 + 1
        oldFNz = fvol.shape[2]

        scf = 1#oldNz **3 / (newNx * newNy * newNz)
        newxIsEven = newFNx%2
        newyIsEven = newFNy%2

        fvol_center_scaled = xp.fft.fftshift(fvol,axes=(0,1)) * scf
        newfvol = xp.zeros((newFNx, newFNy, newFNz), dtype=fvol.dtype)
        if factor >= 1:
            # Effectively zero-padding
            newfvol[newFNx // 2 - oldFNx // 2 + newxIsEven :newFNx // 2 + oldFNx // 2 + isodd + newxIsEven,
                    newFNy // 2 - oldFNy // 2 + newyIsEven :newFNy // 2 + oldFNy // 2 + isodd + newyIsEven,
                    :oldFNz] = fvol_center_scaled
                    # newFNz // 2 - oldFNz // 2 :newFNz // 2 + oldFNz // 2 + 1] = fvol_center_scaled
        else:
            # Effectively cropping
            newfvol = fvol_center_scaled[oldFNx // 2 - newFNx // 2 - newxIsEven: oldFNx // 2 + newFNx // 2 + isodd ,
                                         oldFNy // 2 - newFNy // 2 - newyIsEven: oldFNy // 2 + newFNy // 2 + isodd ,
                                         :newFNz]
        newfvol = xp.fft.fftshift(newfvol, axes=(0,1))

    else:
        newNz = 1
        scf = 1. # / (newNx * newNy * newNz)
        newFNz = 1
        oldNy = int((fvol.shape[1] - 1)*2 + 1*isodd)
        newNy = int(xp.floor(oldNy * factor + 0.5))
        newFNy = newNy // 2 + 1

        fvol_center_scaled = xp.fft.fftshift(fvol,axes=(0))
        newfvol = xp.zeros((newFNx, newFNy), dtype=fvol.dtype) *scf

        newxIsEven = newFNx % 2
        newyIsEven = newFNy % 2


        if factor >= 1:
            newfvol[newFNx // 2 - oldFNx // 2 + newxIsEven :newFNx // 2 + oldFNx // 2 + isodd + newxIsEven, :oldFNy] = fvol_center_scaled
        else:
            newfvol = fvol_center_scaled[oldFNx // 2 - newFNx // 2 - newxIsEven: oldFNx // 2 + newFNx // 2 + isodd, :newFNy  ]

        newfvol = xp.fft.fftshift(newfvol, axes=(0))
        #newfvol = xp.expand_dims(newfvol,2)

    return newfvol

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

def shiftFourier(volume, shift=None, grid=None):
    PI = xp.pi
    sizex, sizey, sizez = volume.shape
    shift_x, shift_y, shift_z = shift

    if grid is None:
        Y,X,Z = xp.meshgrid(xp.arange(volume.shape[0]*1.),xp.arange(volume.shape[1]*1.),xp.arange(volume.shape[2]*1.))
        X -= volume.shape[0] // 2 - 0.
        Y -= volume.shape[1] // 2 - 0.
        Z -= volume.shape[2] // 2 - 0.
        print(X.min(), X.max())
    else:
        X,Y,Z,R = grid


    import matplotlib
    try:
        matplotlib.use('Qt5Agg')
    except:
        pass
    from pylab import imshow, show



    fshift  = (xp.cos(2 * PI / sizex * X * shift_x) - xp.sin(2 * PI / sizex * X * shift_x)*1j).astype(xp.complex64)
    # fshift *= (xp.cos(2 * PI / sizey * Y * shift_y) - xp.sin(2 * PI / sizey * Y * shift_y)*1j)
    # fshift *= (xp.cos(2 * PI / sizez * Z * shift_z) - xp.sin(2 * PI / sizez * Z * shift_z)*1j)

    # f2shift = (xp.cos(2 * PI / sizex * X * -shift_x) - xp.sin(2 * PI / sizex * X * -shift_x) * 1j)
    # f2shift *= (xp.cos(2 * PI / sizey * Y * -shift_y) - xp.sin(2 * PI / sizey * Y * -shift_y) * 1j)
    # f2shift *= (xp.cos(2 * PI / sizez * Z * -shift_z) - xp.sin(2 * PI / sizez * Z * -shift_z) * 1j)

    if 0:
        ff = xp.abs(fshift).real[:, 2, :]
        imshow(ff)
        show()

    return volume*fshift

def shift(volume, shiftX, shiftY, shiftZ, imethod='fourier', twice=False, grid=None):
    """
    shift: Performs a shift on a volume
    @param volume: the volume
    @param shiftX: shift in x direction
    @param shiftY: shift in y direction
    @param shiftZ: shift in z direction
    @param imethod: Select interpolation method. Real space : linear, cubic, spline . Fourier space: fourier
    @param twice: Zero pad volume into a twice sized volume and perform calculation there.
    @return: The shifted volume.
    @author: Yuxiang Chen and Thomas Hrabe
    """
    if imethod == 'fourier':
        fvolume = xp.fft.fftshift(xp.fft.fftn(volume))
        destFourier = shiftFourier(fvolume, shift=[shiftX, shiftY, shiftZ], grid=grid)

        v2 = (xp.fft.ifftn(xp.fft.fftshift(destFourier)))
        print(v2.imag.max())
        return (v2.real)

    else:
        from pytom.voltools import transform

        # now results should be consistent with python2
        centerX = int(volume.shape[0] // 2)
        centerY = int(volume.shape[1] // 2)
        centerZ = int(volume.shape[2] // 2)

        res = xp.zeros_like(volume, dtype=xp.float32)
        transform(volume, output=res, center=[centerX, centerY, centerZ], translation=[shiftX, shiftY, shiftZ],
                  interpolation='filt_bspline', device=device)

        return res

################################
#    Fourier Relevant Stuff    #
################################

def rfft(data):
    """Do 3D Fourier transformaion of data of real numbers.

    @param input data.

    @return: the data after transformation.
    """
    return xp.fft.rfftn(data)

def irfft(data, s=None):
    """Do 3D Inverse Fourier transformaion to get back data of real numbers.

    @param data: input data.
    
    @return: the data after ifft (without the need to scale).
    """
    return xp.fft.irfftn(data, s)

def fft(data):
    """Do 3D Fourier transformaion of data (real or complex numbers).

    @param input data.

    @return: the data after transformation.
    """
    return xp.fft.fftn(data)

def ifft(data):
    """Do 3D Inverse Fourier transformaion to get back data (real or complex numbers).

    @param data: input data.
    
    @return: the data after ifft (without the need to scale).
    """
    return xp.fft.ifftn(data)

def fftshift(data):
    return xp.fft.fftshift(data)

def ifftshift(data):
    return xp.fft.ifftshift(data)

def conv3d(data, kernel):
    """Do 3D convolution.

    @param: input data.
    @param kernel: the kernel to convolve.

    @return: data after convolution.
    """
    from scipy.ndimage.filters import convolve
    d = convolve(data, kernel)
    return d

def fourier_reduced2full(data, isodd=False, reduced_axis=2) -> xpt.NDArray:
    """Return an Hermitian symmetried data.
    Only defined for volumes
    """
    # get the output shape
    new_shape = [0, ] * 3
    for i, s in enumerate(data.shape):
        if i == reduced_axis:
            if isodd:
                new_shape[i] = (s - 1) * 2 + 1
            else:
                new_shape[i] = (s - 1) * 2
        else:
            new_shape[i] = s

    sx, sy, sz = new_shape
    res = xp.zeros(tuple(new_shape), dtype=data.dtype)

    if reduced_axis == 2:
        res[:, :, 0:data.shape[2]] = data

        # calculate the coodinates accordingly
        szz = sz - data.shape[2]
        x, y, z = xp.indices((sx, sy, szz))
        ind = [xp.mod(sx - x, sx), xp.mod(sy - y, sy), szz - z]

        # do the complex conjugate of the second part
        res[:, :, data.shape[2]:] = xp.conj(data[tuple(ind)])
    elif reduced_axis == 1:
        res[:, 0:data.shape[1], :] = data

        # calculate the coodinates accordingly
        syy = sy - data.shape[1]
        x, y, z = xp.indices((sx, syy, sz))
        ind = [xp.mod(sx - x, sx), syy - y, xp.mod(sz - z, sz)]

        # do the complex conjugate of the second part
        res[:, data.shape[1]:, :] = xp.conj(data[tuple(ind)])
    else:
        res[0:data.shape[0], :, :] = data

        # calculate the coodinates accordingly
        sxx = sx - data.shape[0]
        x, y, z = xp.indices((sxx, sy, sz))
        ind = [sxx - x, xp.mod(sy - y, sy), xp.mod(sz - z, sz)]

        # do the complex conjugate of the second part
        res[data.shape[0]:, :, :] = xp.conj(data[tuple(ind)])

    return res


def fourier_full2reduced(data, reduced_axis=-1):
    assert len(data.shape) in [2, 3], 'functionality only defined for arrays with at least two dimensions'
    if reduced_axis == 0:
        return data[0:data.shape[0] // 2 + 1, :, :] if len(data.shape) == 3 else data[0:data.shape[0] // 2 + 1, :]
    elif reduced_axis == 1 and len(data.shape) == 3:
        return data[:, 0:data.shape[1] // 2 + 1, :]
    else:
        return data[..., 0:data.shape[-1]//2+1]


def fourier_filter(data, fltr, human=True):
    if human:
        fltr = ifftshift(fltr)
        fltr = fourier_full2reduced(fltr)

    return irfft(rfft(data) * fltr, data.shape).real.astype(xp.float32)

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

def npyTOsf(data, r, b, m_x, m_y, m_z):
    from pytom.agnostic.interpolation import splineInterpolation
    import time

    t = time.time()

    res = xp.zeros([4*b*b])


    for i in range(4*b*b):
        jj = i//(4*b)
        kk = i % (4*b)
    #for jj in range(0, 2*b):
    #    for kk in range(0,2*b):
        the = xp.pi * (2 * jj+1) / (4 * b)
        phi = xp.pi * kk / b
        x = r * xp.cos(phi) * xp.sin(the)
        y = r * xp.sin(phi) * xp.sin(the)
        z = r * xp.cos(the)
        res[jj * 2 * b + kk] = splineInterpolation(data, x+m_x, y+m_y, z+m_z)

    print(time.time()-t)
    return res
