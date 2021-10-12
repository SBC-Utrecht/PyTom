#! /usr/bin/env python
# -*- coding: utf-8 -*-

from pytom.gpu.initialize import xp, device

from pytom.agnostic.filter import gaussian3d
from numpy.random import standard_normal

from pytom.agnostic.transform import resize as RESIZE

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

def rotation_matrix_zyz(angle):
    """Return the 3x3 rotation matrix of an Euler angle in ZXZ convention.
    Note the order of the specified angle should be [Phi, Psi, Theta], or [Z1, Z2, X] in Pytom format.

    @param angle: list of [Phi, Psi, Theta] in degree.

    @return: rotation matrix.
    """
    assert len(angle) == 3

    z1 = angle[0]
    z2 = angle[1]
    y = angle[2]

    zm1 = rotation_matrix_z(z1)
    xm = rotation_matrix_y(y)
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
    @type particle: Either L{pytom_volume.vol} or string specifying the particle file name
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
    @type mask: L{pytom_volume.vol}
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
    from pytom_volume import transformSpline, vol
    from pytom.agnostic.correlation import FSC
    from pytom.agnostic.filter import filter_volume_by_profile
    from pytom.agnostic.structures import Alignment
    from pytom.agnostic.correlation import nxcc
    from pytom.voltools import transform

    assert isinstance(object=vol1, class_or_type_or_tuple=vol), "alignVolumesAndFilterByFSC: vol1 must be of type vol"
    assert isinstance(object=vol2, class_or_type_or_tuple=vol), "alignVolumesAndFilterByFSC: vol2 must be of type vol"
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
    transform(vol2, output=vol2_alig, rotation=[optiRot[0], optiRot[2], optiRot[1]], rotation_order='rzxz',
              center=[int(vol2.shape[0]//2),int(vol2.shape[1]//2),int(vol2.shape[2]//2)], interpolation='filt_bspline',
              translation=[optiTrans[0], optiTrans[1], optiTrans[2]], device=device)
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
    from pytom_volume import vol, subvolume

    if volume.__class__ != vol:
        return volume[sub_startX:sub_startX+stepSizeX, sub_startY:sub_startY+stepSizeY, sub_startZ:sub_startZ+stepSizeZ]
    else:
        return subvolume(volume, sub_startX, sub_startY, sub_startZ, stepSizeX, stepSizeY, stepSizeZ)

def putSubVolume(subvolume, volume, startX, startY, startZ):
    from pytom_volume import vol, putSubVolume
    from pytom_numpy import vol2npy

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


def convert_angles(angles, rotation_order='rzxz', return_order='rzyz'):
    from pytom.voltools.utils.matrices import rotation_matrix
    import numpy as np

    mat_dict = {'x':rotation_matrix_x, 'y':rotation_matrix_y, 'z': rotation_matrix_z}

    # shift direction of rotation to match CCW definition used in mat2??? calculations.
    # angles = -1* np.array(angles)

    m = rotation_matrix(angles,rotation_order=rotation_order)

    return_funcs = {'xyz':mat2xyz, 'xzy':mat2xzy, 'yxz':mat2yxz, 'yzx':mat2yzx, 'zxy':mat2zxy, 'zyx':mat2zyx,
                    'xyx':mat2xyx, 'xzx':mat2xzx, 'yxy':mat2yxy, 'yzy':mat2yzy, 'zxz':mat2zxz, 'zyz':mat2zyz}

    angs = return_funcs[return_order[-3:]](m)

    if return_order[0] == 's':
        angs = np.array(angs)[::-1]

    return angs

def mat2xyz(rotation_matrix):
    #y = asin(r02), x = atan2(−r12, r22).z = atan2(−r01, r00)
    from numpy import pi, abs, rad2deg, arctan2 as atan2, arccos as acos, arcsin as asin
    r = rotation_matrix
    y = rad2deg(asin(r[0, 2]))

    epsilon = 1E-4

    if abs(90 - y) < epsilon or abs(90 + y) < epsilon:
        z = rad2deg(atan2(r[1,0], r[1,1]))
        x = 0
    else:
        x = rad2deg(atan2(-r[1, 2], r[2, 2]))
        z = rad2deg(atan2(-r[0, 1], r[0, 0]))

    return (x, y, z)

def mat2xzy(rotation_matrix):
    from numpy import pi, abs, rad2deg, arctan2 as atan2, arccos as acos, arcsin as asin
    r = rotation_matrix
    z = rad2deg(asin(-r[0, 1]))

    if abs(90 - z) < 1E-4 or abs(90 + z) < 1E-4:
        y = rad2deg(atan2(-r[2, 0], r[2, 2]))
        x = 0
    else:
        x = rad2deg(atan2(r[2, 1], r[1, 1]))
        y = rad2deg(atan2(r[0, 2], r[0, 0]))

    return (x, z, y)

def mat2yxz(rotation_matrix):
    from numpy import pi, abs, rad2deg, arctan2 as atan2, arccos as acos, arcsin as asin
    r = rotation_matrix
    x = rad2deg(asin(-r[1, 2]))

    if abs(90 - x) < 1E-4 or abs(90 + x) < 1E-4:
        z = rad2deg(atan2(-r[0, 1], r[0,0]))
        y = 0
    else:
        y = rad2deg(atan2(r[0, 2], r[2, 2]))
        z = rad2deg(atan2(r[1, 0], r[1, 1]))

    return (y, x, z)

def mat2yzx(rotation_matrix):
    from numpy import pi, abs, rad2deg, arctan2 as atan2, arccos as acos, arcsin as asin
    r = rotation_matrix

    z = rad2deg(asin(r[1, 0]))

    if abs(90- z) < 1E-8 or abs(90 + z) < 1E-8:
        x = rad2deg(atan2(r[2, 1], r[2, 2]))
        y = 0
    else:
        y = rad2deg(atan2(-r[2, 0], r[0, 0]))
        x = rad2deg(atan2(-r[1, 2], r[1, 1]))

    return (y, z, x)

def mat2zxy(rotation_matrix):
    from numpy import pi, abs, rad2deg, arctan2 as atan2, arccos as acos, arcsin as asin
    r = rotation_matrix
    x = rad2deg(asin(r[2, 1]))

    if abs(90 - x) < 1E-8 or abs(90 + x) < 1E-8:
        y = rad2deg(atan2(r[0, 2], r[0, 0]))
        z = 0
    else:
        z = rad2deg(atan2(-r[0, 1], r[1, 1]))
        y = rad2deg(atan2(-r[2, 0], r[2, 2]))

    return (z, x, y)

def mat2zyx(rotation_matrix):
    from numpy import pi, abs, rad2deg, arctan2 as atan2, arccos as acos, arcsin as asin
    r = rotation_matrix
    y = rad2deg(asin(-r[2, 0]))

    if abs(90 - y) < 1E-8 or abs(90 + y) < 1E-8:
        x = rad2deg(atan2(-r[1, 2], r[1, 1]))
        z = 0
    else:
        z = rad2deg(atan2(r[1, 0], r[0, 0]))
        x = rad2deg(atan2(r[2, 1], r[2, 2]))

    return (z, y, x)

def mat2xyx(rotation_matrix):
    from numpy import pi, abs, rad2deg, arctan2 as atan2, arccos as acos, arcsin as asin
    r = rotation_matrix
    y = rad2deg(acos(r[0, 0]))

    if abs(180 - y) < 1E-8 or abs(y) < 1E-8:
        x1 = rad2deg(atan2(-r[1, 2], r[1, 1]))
        x0 = 0
    else:
        x0 = rad2deg(atan2(r[1, 0], -r[2, 0]))
        x1 = rad2deg(atan2(r[0, 1], r[0, 2]))

    return (x0, y, x1)

def mat2xzx(rotation_matrix):
    from numpy import pi, abs, rad2deg, arctan2 as atan2, arccos as acos, arcsin as asin
    r = rotation_matrix
    z = rad2deg(acos(r[0, 0]))

    if abs(180 - z) < 1E-8 or abs(z) < 1E-8:
        x1 = rad2deg(atan2(r[2, 1], r[2, 2]))
        x0 = 0
    else:
        x0 = rad2deg(atan2(r[2, 0], r[1, 0]))
        x1 = rad2deg(atan2(r[0, 2], -r[0, 1]))

    return (x0, z, x1)

def mat2yxy(rotation_matrix):
    from numpy import pi, abs, rad2deg, arctan2 as atan2, arccos as acos, arcsin as asin
    r = rotation_matrix
    x = rad2deg(acos(r[1, 1]))

    if abs(180 - x) < 1E-8 or abs(x) < 1E-8:
        y1 = rad2deg(atan2(r[0, 2], r[0, 0]))
        y0 = 0
    else:
        y0 = rad2deg(atan2(r[0, 1], r[2, 1]))
        y1 = rad2deg(atan2(r[1, 0], -r[1, 2]))

    return (y0, x, y1)

def mat2yzy(rotation_matrix):
    from numpy import pi, abs, rad2deg, arctan2 as atan2, arccos as acos, arcsin as asin
    r = rotation_matrix
    z = rad2deg(acos(r[1, 1]))

    if abs(180 - z) < 1E-8 or abs(z) < 1E-8:
        y1 = rad2deg(atan2(-r[2, 0], r[2, 2]))
        y0 = 0
    else:
        y0 = rad2deg(atan2(r[2, 1], -r[0, 1]))
        y1 = rad2deg(atan2(r[1, 2], r[1, 0]))

    return (y0, z, y1)

def mat2zxz(rotation_matrix):
    from numpy import pi, abs, rad2deg, arctan2 as atan2, arccos as acos
    r = rotation_matrix
    x = rad2deg(acos(r[2,2]))

    if abs(180 - x) < 1E-8 or abs(x) < 1E-8:
        z1 = rad2deg(atan2(-r[0,1], r[0,0]))
        z0 = 0
    else:
        z0 = rad2deg(atan2(r[0,2], -r[1,2]))
        z1 = rad2deg(atan2(r[2,0], r[2,1]))

    return (z0, x, z1)

def mat2zyz(rotation_matrix):
    from numpy import pi, abs, rad2deg, arctan2 as atan2, arccos as acos
    r = rotation_matrix

    y = rad2deg(acos(r[2,2]))

    if abs(180 - y) < 1E-8 or abs(y) < 1E-8:
        z1 = rad2deg(atan2(r[1,0], r[1,1]))
        z0 = 0
    else:
        z0 = rad2deg(atan2(r[1,2], r[0,2]))
        z1 = rad2deg(atan2(r[2,1], -r[2,0]))

    return (z0, y, z1)

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

    if abs(pi - x) < 1E-8 or abs(x) < 1E-8:
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

    if abs(pi-x) < 1E-8 or abs(x) < 1E-8:
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