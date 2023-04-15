#! /usr/bin/env python
# -*- coding: utf-8 -*-

from pytom.gpu.initialize import xp, device
import numpy as np
from pytom.agnostic.transform import resize as RESIZE

epsilon = 1E-8


def create_sphere(size, radius=-1, sigma=0, num_sigma=2, center=None, gpu=False):
    """Create a 3D sphere volume.
    @param size: size of the resulting volume.
    @param radius: radius of the sphere inside the volume.
    @param sigma: sigma of the Gaussian.
    @param center: center of the sphere.
    @return: sphere inside a volume.
    """


    if size.__class__ == float or len(size) == 1:
        size = (size, size, size)
    assert len(size) == 3
    if center is None:
        center = [size[0]//2, size[1]//2, size[2]//2]
    if radius == -1:
        radius = xp.min(size)//2

    sphere = xp.zeros(size, dtype=xp.float32)
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


    if size.__class__ == float or len(size) == 1:
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
        ind = xp.logical_and(r>radius, r<radius+num_sigma*sigma)
        circle[ind] = xp.exp(-((r[ind] - radius)/sigma)**2/2)

    return circle


def create_grid(shape, center=None):
    cx,cy,cz = [shape[0]//2, shape[1]//2, shape[2]//2] if center is None else center

    Y, X, Z = xp.meshgrid(xp.arange(shape[0]), xp.arange(shape[1]), xp.arange(shape[2]))
    X -= cx
    Y -= cy
    Z -= cz

    return (X,Y,Z)

def prepare_mask(v, threshold, smooth):
    """Prepare a mask according to the given volume.
    Everything above the given threshold will be set to 1.

    @param v: input volume.
    @param threshold: threshold.
    @param smooth: sigma of Gaussian.

    @return: mask.
    """
    from pytom.agnostic.filter import gaussian3d
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
    noise = sd*np.standard_normal(s.shape)+m
    t = data + noise
    return t

def paste_in_center(volume, volume2, gpu=False):
    if 0:
        pass#raise Exception('pasteCenter not defined for gpu yet.')
    else:
        l,l2 = len(volume.shape), len(volume2.shape)
        if l != l2:
            raise ValueError(f"Can't overlap a volume with dimension {l} with a volume of dimension {l2}")

        if not (all(i <= j for i,j in zip(volume.shape, volume2.shape)) or
                all(i >= j for i,j in zip(volume.shape, volume2.shape))):
            raise ValueError(f"Do not know how to overlap shape {volume.shape} with shape {volume2.shape}")

        if len(volume.shape) == 3:
            sx,sy,sz = volume.shape
            SX, SY, SZ = volume2.shape
            if SX <= sx:
                volume2[:,:,:] = volume[sx//2-SX//2:sx//2+SX//2+SX%2,
                                        sy//2-SY//2:sy//2+SY//2+SY%2,
                                        sz//2-SZ//2:sz//2+SZ//2+SZ%2 ]
            else:
                volume2[SX//2-sx//2:SX//2+sx//2+sx%2,
                        SY//2-sy//2:SY//2+sy//2+sy%2,
                        SZ//2-sz//2:SZ//2+sz//2+sz%2 ] = volume
            return volume2

        if len(volume.shape) == 2:
            sx,sy = volume.shape
            SX, SY = volume2.shape
            if SX <= sx:
                volume2[:,:,:] = volume[sx//2-SX//2:sx//2+SX//2+SX%2,
                                        sy//2-SY//2:sy//2+SY//2+SY%2]
            else:
                volume2[SX//2-sx//2:SX//2+sx//2+sx%2,
                        SY//2-sy//2:SY//2+sy//2+sy%2] = volume
            return volume2


def euclidian_distance(pos1, pos2):
    return xp.linalg.norm(xp.array(pos1)-xp.array(pos2))


def volumesSameSize(v0, v1):
    if len(v0.shape) != len(v1.shape):
        return False
    return xp.array([abs(v0.shape[pos] - v1.shape[pos]) == 0 for pos in range(len(v0.shape))]).sum() == len(v0.shape)


def taper_edges(image, width, taper_mask=None):
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

    if taper_mask is None:
        width = int(round(width))
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
            X, Y = xp.meshgrid(taperY, taperX)
            taper_mask = X * (X < Y) + Y * (Y <= X)

    return image * taper_mask, taper_mask

def resize(*args, **kwargs):
    return RESIZE(*args, **kwargs)


def determineRotationCenter(particle, binning):
    """
    determineRotationCenter:
    @param particle: The particle
    @type particle: Either L{pytom.lib.pytom_volume.vol} or string specifying the particle file name
    @param binning: Binning factor
    @return: [centerX,centerY,centerZ]

    @author: GvdS
    """
    if particle.__class__ == str:
        from pytom.agnostic.io import read
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
    @type invol: L{pytom.lib.pytom_volume.vol} or L{pytom.lib.pytom_volume.vol_comp}
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
        X, Y, Z = xp.meshgrid(xp.arange(-dx // 2, dx // 2 + dx % 2), xp.arange(-dy // 2, dy // 2 + dy % 2), xp.arange(0,dz))
        invol = xp.fft.fftshift(invol,axes=(0,1))
    else:
        X, Y, Z = xp.meshgrid(xp.arange(-dx // 2, dx // 2 + dx % 2), xp.arange(-dy // 2, dy // 2 + dy % 2), xp.arange(-dz // 2, dz // 2 + dz % 2))

    R = xp.sqrt(X ** 2 + Y ** 2 + Z**2).astype(xp.int32)

    invol_out = invol.copy().astype(xp.float32)
    invol_out[invol < lowlimit] = lowval
    invol_out = 1. / invol_out
    invol_out[R >= r_max] = 0

    if dx != dz:
        invol_out = xp.fft.fftshift(invol_out, axes=(0,1))

    return invol_out

def alignVolumesAndFilterByFSC(vol1, vol2, mask=None, nband=None, iniRot=None, iniTrans=None, interpolation='linear',
                               fsc_criterion=0.143, verbose=0):
    """
    align two volumes, compute their FSC, and filter by FSC
    @param vol1: volume 1
    @param vol2: volume 2
    @mask: mask volume
    @type mask: L{pytom.lib.pytom_volume.vol}
    @param nband: Number of bands
    @type nband: L{int}
    @param iniRot: initial guess for rotation
    @param iniTrans: initial guess for translation
    @param interpolation: interpolation type - 'linear' (default) or 'spline'
    @param fsc_criterion: filter -> 0 according to resolution criterion
    @type fsc_criterion: float
    @param verbose: verbose level (0=mute, 1 some output, 2=talkative)
    @type verbose: int
    @type interpolation: str
    @return: (filvol1, filvol2, fsc, fsc_fil, optiRot, optiTrans) i.e., filtered volumes, their FSC, the corresponding\
        filter that was applied to the volumes, and the optimal rotation and translation of vol2 with respect to vol1\
        note: filvol2 is NOT rotated and translated!
    @author: FF
    """
    from pytom.agnostic.correlation import FSC
    from pytom.agnostic.filter import filter_volume_by_profile
    from pytom.agnostic.structures import Alignment
    from pytom.agnostic.correlation import nxcc
    from pytom.voltools import transform

    # filter volumes prior to alignment according to SNR
    fsc = FSC(volume1=vol1, volume2=vol2, numberBands=nband)
    fil = design_fsc_filter(fsc=fsc, fildim=int(vol2.shape[2]//2))
    #filter only one volume so that resulting CCC is weighted by SNR only once
    filvol2 = filter_volume_by_profile(volume=vol2, profile=fil)
    # align vol2 to vol1
    if verbose == 2:
        alignment = Alignment(vol1=vol1, vol2=filvol2, score=nxcc, mask=mask,
                              iniRot=iniRot, iniTrans=iniTrans, opti='fmin_powell', interpolation=interpolation,
                              verbose=verbose)
    else:
        alignment = Alignment(vol1=vol1, vol2=filvol2, score=nxcc, mask=mask,
                              iniRot=iniRot, iniTrans=iniTrans, opti='fmin_powell', interpolation=interpolation,
                              verbose=False)
    optiScore, optiRot, optiTrans = alignment.localOpti( iniRot=iniRot, iniTrans=iniTrans)
    if verbose:
        from pytom.angles.angleFnc import differenceAngleOfTwoRotations
        from pytom.basic.structures import Rotation
        diffAng = differenceAngleOfTwoRotations(rotation1=Rotation(0,0,0), rotation2=optiRot)
        print("Alignment densities: Rotations: %2.3f, %2.3f, %2.3f; Translations: %2.3f, %2.3f, %2.3f " % (optiRot[0],
                                    optiRot[1], optiRot[2], optiTrans[0], optiTrans[1], optiTrans[2]))
        print("Orientation difference: %2.3f deg" % diffAng)
    vol2_alig = xp.zeros(vol2.shape,dtype=xp.float32)
    transform(vol2, output=vol2_alig, rotation=(optiRot[0], optiRot[2], optiRot[1]), rotation_order='rzxz',
              center=(int(vol2.shape[0]//2),int(vol2.shape[1]//2),int(vol2.shape[2]//2)), interpolation='filt_bspline',
              translation=(optiTrans[0], optiTrans[1], optiTrans[2]), device=device)
    # finally compute FSC and filter of both volumes
    if not nband:
        nband = int(vol2.sizeX()/2)
    fsc = FSC(volume1=vol1, volume2=vol2_alig, numberBands=nband)
    fil = design_fsc_filter(fsc=fsc, fildim=int(vol2.shape[0]//2), fsc_criterion=fsc_criterion)
    filvol1 = filter_volume_by_profile(volume=vol1, profile=fil)
    #filvol2 = filter_volume_by_profile( volume=vol2_alig, profile=fil)
    filvol2 = filter_volume_by_profile(volume=vol2, profile=fil)

    return (filvol1, filvol2, fsc, fil, optiRot, optiTrans)


def design_fsc_filter(fsc, fildim=None, fsc_criterion=0.143):
    """
    design spectral filter to weight by SNR of frequency
    @param fsc: input fsc
    @type fsc: 1-d list
    @param fildim: filter dimension
    @type fildim: int
    @return: filter
    @rtype: list
    @author: FF
    """
    from math import sqrt
    from pytom.basic.resolution import getResolutionBandFromFSC
    if not fildim:
        fildim = len(fsc)
    nband = len(fsc)
    if fsc_criterion != 0.0:
        resolutionBand = getResolutionBandFromFSC(fsc, criterion=fsc_criterion)
        smooth = max(resolutionBand/5,2)
    else:
        resolutionBand = len(fsc)
        smooth = 1
    #print "filter: fsc_criterion %2.3f, resolutionBand %d" % (fsc_criterion, resolutionBand)
    # filter by sqrt(FSC)
    fsc_fil = len(fsc)*[0.]
    for ii in range(0,len(fsc)):
        fscval = fsc[ii]
        if fscval > 0.:
            #fsc_fil[ii] = sqrt(fsc[ii])
            if ii <= resolutionBand:
                fsc_fil[ii] = sqrt(fsc[ii])
            elif (ii > resolutionBand) and (ii <= resolutionBand + smooth):
                fsc_fil[ii] = sqrt(fsc[ii]) * (((resolutionBand + smooth) - ii)/smooth)**2  # squared filter
            else:
                fsc_fil[ii] = 0.
        else:
            fsc_fil[ii] = 0.
    #design filter
    fil = xp.array(fildim*[0.])
    if nband != len(fil):
        shrinkfac = 1./(len(fil)/nband)
    for ii in range(len(fil)-1):
        # linear resample fsc if nband ~= size(fil)
        if nband != len(fil):
            ilow = int(xp.floor(shrinkfac*ii))
            dlow = shrinkfac*ii - ilow
            ihi  = ilow+1
            dhi  = ihi - shrinkfac*ii
            fil[ii] = fsc_fil[ilow]*dhi + fsc_fil[ihi]*dlow
        else:
            fil[ii] = fsc_fil[ii]
    return fil



def subvolume(volume, sub_startX, sub_startY, sub_startZ, stepSizeX, stepSizeY, stepSizeZ):
    from pytom.lib.pytom_volume import vol, subvolume

    if volume.__class__ != vol:
        return volume[sub_startX:sub_startX+stepSizeX, sub_startY:sub_startY+stepSizeY, sub_startZ:sub_startZ+stepSizeZ]
    else:
        return subvolume(volume, sub_startX, sub_startY, sub_startZ, stepSizeX, stepSizeY, stepSizeZ)

def putSubVolume(subvolume, volume, startX, startY, startZ):
    from pytom.lib.pytom_volume import vol, putSubVolume
    from pytom.lib.pytom_numpy import vol2npy

    if volume.__class__ == vol and subvolume.__class__ == vol:
        putSubVolume(subvolume, volume, startX, startY, startZ)
    elif volume.__class__ == vol:
        volume = vol2npy(volume).copy()
        sx,sy,sz = subvolume.shape
        volume[startX:startX+sx, startY:startY+sy, startZ:startZ+sz] = subvolume[:,:,:]
    elif subvolume.__class__ == vol:
        subvolume = vol2npy(subvolume).copy()
        sx,sy,sz = subvolume.shape
        volume[startX:startX+sx, startY:startY+sy, startZ:startZ+sz] = subvolume[:,:,:]
    else:
        sx, sy, sz = subvolume.shape
        volume[startX:startX + sx, startY:startY + sy, startZ:startZ + sz] = subvolume[:, :, :]

# ================================ ROTATION MATRICES ==============================


def rotation_matrix_x(angle):
    """Return the 3x3 rotation matrix around x axis.

    @param angle: rotation angle around x axis (in degree).

    @return: rotation matrix.
    """
    angle = np.deg2rad(angle)
    mtx = np.zeros((3, 3))
    mtx[1, 1] = np.cos(angle)
    mtx[2, 1] = np.sin(angle)
    mtx[2, 2] = np.cos(angle)
    mtx[1, 2] = -np.sin(angle)
    mtx[0, 0] = 1

    return mtx


def rotation_matrix_y(angle):
    """Return the 3x3 rotation matrix around y axis.

    @param angle: rotation angle around y axis (in degree).

    @return: rotation matrix.
    """
    angle = np.deg2rad(angle)
    mtx = np.zeros((3, 3))
    mtx[0, 0] = np.cos(angle)
    mtx[2, 0] = -np.sin(angle)
    mtx[2, 2] = np.cos(angle)
    mtx[0, 2] = np.sin(angle)
    mtx[1, 1] = 1

    return mtx


def rotation_matrix_z(angle):
    """Return the 3x3 rotation matrix around z axis.

    @param angle: rotation angle around z axis (in degree).

    @return: rotation matrix.
    """
    angle = np.deg2rad(angle)
    mtx = np.zeros((3, 3))
    mtx[0, 0] = np.cos(angle)
    mtx[1, 0] = np.sin(angle)
    mtx[1, 1] = np.cos(angle)
    mtx[0, 1] = -np.sin(angle)
    mtx[2, 2] = 1

    return mtx


def rotation_matrix_zxz(zxz):
    """Return the 3x3 rotation matrix of an Euler angle in ZXZ convention.
    Note the order of the specified angle should be [Phi, Theta, Psi], or [Z1, X, Z2].
    Rotation matrix multiplied in order mat(Z2) * mat(X) * mat(Z1).

    @param angle: list of [Phi, Theta, Psi] in degree.

    @return: rotation matrix.
    """
    assert len(zxz) == 3

    z1, x, z2 = zxz

    zm1 = rotation_matrix_z(z1)
    xm = rotation_matrix_x(x)
    zm2 = rotation_matrix_z(z2)

    res = np.dot(zm2, np.dot(xm, zm1))

    return res


def rotation_matrix_zyz(zyz):
    """Return the 3x3 rotation matrix of an Euler angle in ZYZ convention.
    Note the order of the specified angle should be [Phi, Theta, Psi], or [Z1, X, Z2].
    Rotation matrix multiplied in order mat(Z2) * mat(X) * mat(Z1)

    @param angle: list of [Phi, Theta, Psi] in degree.

    @return: rotation matrix.
    """
    assert len(zyz) == 3

    z1, y, z2 = zyz

    zm1 = rotation_matrix_z(z1)
    xm = rotation_matrix_y(y)
    zm2 = rotation_matrix_z(z2)

    res = np.dot(zm2, np.dot(xm, zm1))

    return res


def rotation_distance(ang1, ang2, rotation_order='zxz'):
    """Given two angles (lists), calculate the angular distance (degree).

    @param ang1: angle 1. [Phi, Psi, Theta], or [Z1, Z2, X] in Pytom format.
    @param ang2: angle 2. [Phi, Psi, Theta], or [Z1, Z2, X] in Pytom format.

    @return: rotation distance in degree.
    """
    mtx1 = rotation_matrix_zxz(ang1)
    mtx2 = rotation_matrix_zxz(ang2)
    res = np.dot(mtx1, np.linalg.inv(mtx2))  # elementwise multiplication
    trace = np.sum(res)

    temp = 0.5 * (trace - 1.0)
    if temp >= 1.0:
        return 0.0
    if temp <= -1.0:
        return 180
    return np.arccos(temp) * 180 / np.pi


def rotation_matrix(angles, rotation_order='zxz', multiplication='post'):

    assert len(angles) == 3, "should provide 3 angles"
    assert multiplication in ['pre', 'post'], "multiplication can only be pre or post"

    mat_dict = {'x': rotation_matrix_x, 'y': rotation_matrix_y, 'z': rotation_matrix_z}

    mtxs = []
    for angle, rot in zip(angles, rotation_order):
        mtxs.append(mat_dict[rot](-angle) if multiplication == 'post' else mat_dict[rot](angle))

    if multiplication == 'post':
        return np.dot(np.dot(mtxs[0], mtxs[1]), mtxs[2])
    else:
        return np.dot(mtxs[2], np.dot(mtxs[1], mtxs[0]))


def convert_angles(angles, rotation_order='zxz', return_order='zyz', multiplication='post'):
    # get the rotation matrix with the input order
    m = rotation_matrix(angles, rotation_order=rotation_order, multiplication=multiplication)
    # get the angles with the specified output order
    return mat2ord(m, return_order=return_order, multiplication=multiplication)

# ================================= mat2... =========================================
# these are all defined for post-multiplication with the matrix: dot(mat, Rc)


def mat2ord(rotation_matrix, return_order='zyz', multiplication='post'):
    assert multiplication in ['pre', 'post'], "multiplication can only be pre or post"
    assert len(rotation_matrix.shape) == 2 and all([s == 3 for s in rotation_matrix.shape]), \
        "invalid rotation matrix shape"

    return_funcs = {'xyz': mat2xyz, 'xzy': mat2xzy, 'yxz': mat2yxz, 'yzx': mat2yzx, 'zxy': mat2zxy, 'zyx': mat2zyx,
                    'xyx': mat2xyx, 'xzx': mat2xzx, 'yxy': mat2yxy, 'yzy': mat2yzy, 'zxz': mat2zxz, 'zyz': mat2zyz}

    # if 'pre' multiplication, invert the matrix
    res = return_funcs[return_order](np.linalg.inv(rotation_matrix)) if \
        multiplication == 'pre' else return_funcs[return_order](rotation_matrix)

    # always take negative of angles
    return tuple([-r for r in res])  # if multiplication == 'post' else res


def mat2xyz(r):
    y = np.rad2deg(np.arcsin(r[0, 2]))

    if np.abs(90 - y) < epsilon or np.abs(90 + y) < epsilon:
        z = np.rad2deg(np.arctan2(r[1,0], r[1,1]))
        x = 0
    else:
        x = np.rad2deg(np.arctan2(-r[1, 2], r[2, 2]))
        z = np.rad2deg(np.arctan2(-r[0, 1], r[0, 0]))

    return (x, y, z)


def mat2xzy(r):
    z = np.rad2deg(np.arcsin(-r[0, 1]))

    if np.abs(90 - z) < epsilon or np.abs(90 + z) < epsilon:
        y = np.rad2deg(np.arctan2(-r[2, 0], r[2, 2]))
        x = 0
    else:
        x = np.rad2deg(np.arctan2(r[2, 1], r[1, 1]))
        y = np.rad2deg(np.arctan2(r[0, 2], r[0, 0]))

    return (x, z, y)


def mat2yxz(r):
    x = np.rad2deg(np.arcsin(-r[1, 2]))

    if np.abs(90 - x) < epsilon or np.abs(90 + x) < epsilon:
        z = np.rad2deg(np.arctan2(-r[0, 1], r[0,0]))
        y = 0
    else:
        y = np.rad2deg(np.arctan2(r[0, 2], r[2, 2]))
        z = np.rad2deg(np.arctan2(r[1, 0], r[1, 1]))

    return (y, x, z)


def mat2yzx(r):
    z = np.rad2deg(np.arcsin(r[1, 0]))

    if np.abs(90- z) < epsilon or np.abs(90 + z) < epsilon:
        x = np.rad2deg(np.arctan2(r[2, 1], r[2, 2]))
        y = 0
    else:
        y = np.rad2deg(np.arctan2(-r[2, 0], r[0, 0]))
        x = np.rad2deg(np.arctan2(-r[1, 2], r[1, 1]))

    return (y, z, x)


def mat2zxy(r):
    x = np.rad2deg(np.arcsin(r[2, 1]))

    if np.abs(90 - x) < epsilon or np.abs(90 + x) < epsilon:
        y = np.rad2deg(np.arctan2(r[0, 2], r[0, 0]))
        z = 0
    else:
        z = np.rad2deg(np.arctan2(-r[0, 1], r[1, 1]))
        y = np.rad2deg(np.arctan2(-r[2, 0], r[2, 2]))

    return (z, x, y)


def mat2zyx(r):
    y = np.rad2deg(np.arcsin(-r[2, 0]))

    if np.abs(90 - y) < epsilon or np.abs(90 + y) < epsilon:
        x = np.rad2deg(np.arctan2(-r[1, 2], r[1, 1]))
        z = 0
    else:
        z = np.rad2deg(np.arctan2(r[1, 0], r[0, 0]))
        x = np.rad2deg(np.arctan2(r[2, 1], r[2, 2]))

    return (z, y, x)


def mat2xyx(r):
    y = np.rad2deg(np.arccos(r[0, 0]))

    if np.abs(180 - y) < epsilon or np.abs(y) < epsilon:
        x1 = np.rad2deg(np.arctan2(-r[1, 2], r[1, 1]))
        x0 = 0
    else:
        x0 = np.rad2deg(np.arctan2(r[1, 0], -r[2, 0]))
        x1 = np.rad2deg(np.arctan2(r[0, 1], r[0, 2]))

    return (x0, y, x1)


def mat2xzx(r):
    z = np.rad2deg(np.arccos(r[0, 0]))

    if np.abs(180 - z) < epsilon or np.abs(z) < epsilon:
        x1 = np.rad2deg(np.arctan2(r[2, 1], r[2, 2]))
        x0 = 0
    else:
        x0 = np.rad2deg(np.arctan2(r[2, 0], r[1, 0]))
        x1 = np.rad2deg(np.arctan2(r[0, 2], -r[0, 1]))

    return (x0, z, x1)


def mat2yxy(r):
    x = np.rad2deg(np.arccos(r[1, 1]))

    if np.abs(180 - x) < epsilon or np.abs(x) < epsilon:
        y1 = np.rad2deg(np.arctan2(r[0, 2], r[0, 0]))
        y0 = 0
    else:
        y0 = np.rad2deg(np.arctan2(r[0, 1], r[2, 1]))
        y1 = np.rad2deg(np.arctan2(r[1, 0], -r[1, 2]))

    return (y0, x, y1)


def mat2yzy(r):
    z = np.rad2deg(np.arccos(r[1, 1]))

    if np.abs(180 - z) < epsilon or np.abs(z) < epsilon:
        y1 = np.rad2deg(np.arctan2(-r[2, 0], r[2, 2]))
        y0 = 0
    else:
        y0 = np.rad2deg(np.arctan2(r[2, 1], -r[0, 1]))
        y1 = np.rad2deg(np.arctan2(r[1, 2], r[1, 0]))

    return (y0, z, y1)


def mat2zxz(r):
    x = np.rad2deg(np.arccos(r[2,2]))

    if np.abs(180 - x) < epsilon or np.abs(x) < epsilon:
        z1 = np.rad2deg(np.arctan2(-r[0,1], r[0,0]))
        z0 = 0
    else:
        z0 = np.rad2deg(np.arctan2(r[0,2], -r[1,2]))
        z1 = np.rad2deg(np.arctan2(r[2,0], r[2,1]))

    return (z0, x, z1)


def mat2zyz(r):
    y = np.rad2deg(np.arccos(r[2,2]))

    if np.abs(180 - y) < epsilon or np.abs(y) < epsilon:
        z1 = np.rad2deg(np.arctan2(r[1,0], r[1,1]))
        z0 = 0
    else:
        z0 = np.rad2deg(np.arctan2(r[1,2], r[0,2]))
        z1 = np.rad2deg(np.arctan2(r[2,1], -r[2,0]))

    return (z0, y, z1)


# ========================================================================================

def zxz2zyz(z0, x, z1):
    from numpy import cos, sin, deg2rad, arctan2 as atan2, rad2deg, pi, abs

    cz0 = cos(deg2rad(z0))
    cz1 = cos(deg2rad(z1))
    cx  = cos(deg2rad(x))
    sx  = sin(deg2rad(x))
    sz0 = sin(deg2rad(z0))
    sz1 = sin(deg2rad(z1))

    r02 = sx * sz0
    r10 = cz1* sz0 + cx*cz0*sz1
    r11 = cx*cz0*cz1-sz0*sz1
    r12 = -sx*cz0
    r20 = sx * sz1
    r21 = sx* cz1

    y = x

    if abs(pi - x) < epsilon or abs(x) < epsilon:
        z1 = rad2deg(atan2(r10, r11))
        z0 = 0
    else:
        z0 = rad2deg(atan2(r12, r02))
        z1 = rad2deg(atan2(r21,-r20))

    return (z0,y,z1)

def zyz2zxz(z0, y, z1):
    from numpy import cos, sin, deg2rad, zeros, arccos as acos, arctan2 as atan2, rad2deg, arcsin as asin, abs, pi

    cz0 = cos(deg2rad(z0))
    cz1 = cos(deg2rad(z1))
    cy = cos(deg2rad(y))
    sy = sin(deg2rad(y))
    sz0 = sin(deg2rad(z0))
    sz1 = sin(deg2rad(z1))

    r00 = cy * cz0 * cz1 - sz0 * sz1
    r01 = -cz1 * sz0 - cy * cz0 * sz1
    r02 = sy * cz0
    r12 = sy * sz0
    r20 = -sy * cz1
    r21 = sy * sz1

    x = y

    if abs(pi-x) < epsilon or abs(x) < epsilon:
        z1 = rad2deg(atan2(-r01, r00))
        z0 = 0
    else:
        z0 = rad2deg(atan2(r02, -r12))
        z1 = rad2deg(atan2(r20, r21))

    return (z0, x, z1)

def is1D(volume):
    '''This function checks if a volume is 1D
    @param volume: ndarray
    @type volume: ndarray numpy or cupy
    @return: returns boolean indicating if volume is 1D
    @rtype: bool'''

    return (volume.ndim == 1)

def is2D(volume):
    '''This function checks if a volume is 2D
        @param volume: ndarray
        @type volume: ndarray numpy or cupy
        @return: returns boolean indicating if volume is 2D
        @rtype: bool'''

    return (volume.ndim == 2)

def is3D(volume):
    '''This function checks if a volume is 3D
        @param volume: ndarray
        @type volume: ndarray numpy or cupy
        @return: returns boolean indicating if volume is 3D
        @rtype: bool'''

    return (volume.ndim == 3)

def convert_operation_order_str2list(operation_string):
    '''This function convert a string of operations to a list of numbers as used in general_transform_crop
        @param operation_string: a string of S R T in any order. Each letter should appear exactly once.
        @type volume: str
        @return: list of 3 ints
        @rtype: list'''

    a = {'R': 0, 'T': 1, 'S': 2}
    order = [2, 1, 0]
    for n, operation in enumerate(operation_string):
        order[a[operation]] = n
    return order

def convert_operation_order_list2str(operation_list):
    '''This function convert a string of operations to a list of numbers as used in general_transform_crop
        @param operation_list: a list of 2,1,0 in any order. Each number should appear exactly once.
        @type volume: list
        @return: str of RTS in any order
        @rtype: str'''

    t=['R', 'T', 'S']
    order = ['', '', '']
    for n, operation in enumerate(operation_list):
        order[operation] = t[n]

    return ''.join(order)
