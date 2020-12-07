#! /usr/bin/env python
# -*- coding: utf-8 -*-

from pytom.gpu.initialize import xp, device

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
        X, Y = xp.meshgrid(taperY, taperX)
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
    from pytom.tompy.correlation import FSC
    from pytom.tompy.filter import filter_volume_by_profile
    from pytom.tompy.structures import Alignment
    from pytom.tompy.correlation import nxcc
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