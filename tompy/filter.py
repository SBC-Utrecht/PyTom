#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
basic filters operating on numpy arrays
"""
from pytom.gpu.initialize import xp
import scipy
import numpy as np

def normalize(v):
    """Normalize the data according to standard deviation

    @param v: input volume.

    @return: Normalized volume.
    """
    m = np.mean(v)
    v = v-m
    s = np.std(v)
    v = v/s
    return v

def bandpass(v, low=0, high=-1, sigma=0):
    """Do a bandpass filter on a given volume.

    @param v: input volume.
    @param low: low frequency in k-space.
    @param high: high frequency in k-space.
    @param sigma: smoothness factor.

    @return: bandpass filtered volume.
    """
    assert low >= 0, "lower limit must be >= 0"

    from pytom.tompy.tools import create_sphere

    if high == -1:
        high = np.min(v.shape)/2
    assert low < high, "upper bandpass must be > than lower limit"

    if low == 0:
        mask = create_sphere(v.shape, high, sigma)
    else:
        # BUG! TODO
        # the sigma
        mask = create_sphere(v.shape, high, sigma) - create_sphere(v.shape, low, sigma)

    from pytom.tompy.transform import fourier_filter

    res = fourier_filter(v, mask, True)

    return res


def median3d(data, size=3):
    """Median filter.

    @param data: data to be filtered.
    @param size: size of the median filter.
    @type size: C{int}

    @return: filtered image.
    """
    from scipy.ndimage.filters import median_filter
    assert type(size) == int, "median3d: size must be integer"
    d = median_filter(data, size)
    return d


def gaussian3d(data, sigma=0):
    """Gaussian filter.

    @param data: data to be filtered.
    @param sigma: sigma of Gaussian.

    @return: filtered data.
    """
    from scipy.ndimage.filters import gaussian_filter
    d = gaussian_filter(data, sigma)
    return d


class Wedge(object):
    """Defines the missing wedge filter in Fourier space."""
    def __init__(self):
        super(Wedge, self).__init__()

    def apply(self, data):
        raise NotImplementedError('Abstract method! Please overwrite it.')

    def toSphericalFunc(self, bw, radius):
        raise NotImplementedError('Abstract method! Please overwrite it.')


class GeneralWedge(Wedge):
    """General wedge."""
    def __init__(self, wedge_vol, half=True, isodd=False):
        """Initialize a general wedge with given Fourier volume.

        @param half: if true, the given volume is in half and zero frequency at the corner.
                     Otherwise, the given volume is full and zero frequency in the center.
        @param isodd: when half is true, this field tells the z dimension of the full Fourier volume is odd-/even-sized.
        """
        self.set_wedge_volume(wedge_vol, half)

    def set_wedge_volume(self, wedge_vol, half=True, isodd=False):
        if half:
            self._volume = wedge_vol

            # human understandable version with 0-freq in the center
            from transform import fourier_reduced2full, fftshift
            self._whole_volume = fftshift(fourier_reduced2full(self._volume, isodd))
        else:
            self._whole_volume = wedge_vol

            from transform import fourier_full2reduced, ifftshift
            self._volume = fourier_full2reduced(ifftshift(self._whole_volume))

    def apply(self, data, rotation=None):
        if rotation is not None: # rotate the wedge first
            assert len(rotation) == 3
            from transform import rotate3d, fourier_full2reduced, ifftshift
            filter_vol = rotate3d(self._whole_volume, rotation[0], rotation[1], rotation[2], order=1)
            filter_vol = fourier_full2reduced(ifftshift(filter_vol))
        else:
            filter_vol = self._volume

        from transform import fourier_filter
        res = fourier_filter(data, filter_vol, False)

        return res

    def toSphericalFunc(self, bw, radius):
        assert(bw<=128)

        # start sampling
        from vol2sf import vol2sf
        sf = vol2sf(self._whole_volume, radius, bw)
        
        return sf


class SingleTiltWedge(Wedge):
    """Missing wedge of single tilt geometry. Assume Y axis is the rotation axis."""
    def __init__(self, start_ang=30, end_ang=30, smooth=0):
        super(SingleTiltWedge, self).__init__()
        self.start_ang = start_ang
        self.end_ang = end_ang
        self._volume_shape = None # store the wedge volume shape (whole!)
        self._volume = None # store the wedge volume in k-space (only half!)

        self._bw = None # store the bandwidth of the spherical function
        self._sf = None # store the spherical function

        self.smooth=smooth

    def _create_wedge_volume(self, size, cutOffRadius=None):
        if cutOffRadius is None: cutOffRadius = size[0]//2

        self._volume = create_wedge(self.start_ang, self.end_ang, cutOffRadius, size[0], size[1], size[2], self.smooth)
        self._volume_shape = size

    def apply(self, data, rotation=None):
        """
        @param rotation: apply rotation to the wedge first
        """
        # if no missing wedge
        if self.start_ang == -90 and self.end_ang == 90:
            return data

        if self._volume is not None and np.array_equal(self._volume_shape, data.shape):
            pass
        else:
            self._create_wedge_volume(data.shape)

        if rotation is not None: # rotate the wedge first
            assert len(rotation) == 3
            from pytom.tompy.transform import rotate3d, fourier_reduced2full, fourier_full2reduced, fftshift, ifftshift
            isodd = self._volume_shape[2] % 2
            filter_vol = fftshift(fourier_reduced2full(self._volume, isodd))
            filter_vol = rotate3d(filter_vol, rotation[0], rotation[1], rotation[2], order=1) # linear interp!
            filter_vol = fourier_full2reduced(ifftshift(filter_vol))
        else:
            filter_vol = self._volume

        from pytom.tompy.transform import fourier_filter
        res = fourier_filter(data, filter_vol, False)

        return res

    def returnWedgeVolume(self, size, rotation=None):
        """Return the wedge volume in full size and zero in the center
        @param size: size of wedge
        @type size: C{list}
        @param rotation: rotation (3-dim vector if Euler angles)
        @type rotation: C{list}
        @return: wedge volume
        @rtype: Numpy array
        """
        assert len(size) == 3, "returnWedgeVolume: size must be 3-dim list"
        
        if self._volume is not None and np.array_equal(self._volume_shape, size):
            pass
        else:
            self._create_wedge_volume(size)

        from pytom.tompy.transform import rotate3d, fourier_reduced2full, fftshift
        isodd = self._volume_shape[2] % 2
        wedge_vol = fftshift(fourier_reduced2full(self._volume, isodd))
        if rotation is not None: # rotate the wedge first
            assert len(rotation) == 3
            wedge_vol = rotate3d(wedge_vol, rotation[0], rotation[1], rotation[2], order=1)

        return wedge_vol

    def toSphericalFunc(self, bw, radius=None, threshold=0.5):
        """Convert the wedge from k-space to a spherical function. \
        currently some hard-coded parameters in - bw <=128, r=45 for max bw, default vol 100
        
        @param bw: bandwidth of the spherical function (must be <=128).
        @param radius: radius in k-space. For general Wedge, not used for SingleTiltWedge.
        @param threshold: threshold, above which the value there would be set to 1.

        @return: a spherical function in numpy.array - default 100x100x100 if no self.vol defined
        """
        assert(bw<=128), "toSphericalFunc: bw currently limited to <= 128"

        # if no missing wedge
        if self.start_ang == -90 and self.end_ang == 90:
            self._sf = np.ones((4*bw**2,))
            return self._sf

        r = 45 # this radius and the volume size should be sufficient for sampling b <= 128
        if self._volume is None or np.min(self._volume.shape) < 100:
            self._create_wedge_volume((100,100,100))
        
        if self._bw == bw and self._sf is not None:
            return self._sf
        else:
            self._bw = bw
        
        from pytom.tompy.transform import fourier_reduced2full, fftshift
        isodd = self._volume_shape[2] % 2
        filter_vol = fftshift(fourier_reduced2full(self._volume, isodd))

        # start sampling
        from math import pi, sin, cos
        res = []
        
        for j in range(2*bw):
            for k in range(2*bw):
                the = pi*(2*j+1)/(4*bw) # (0,pi)
                phi = pi*k/bw # [0,2*pi)
                
                # this part actually needs interpolation
                x = int(cos(phi)*sin(the)*r+50)
                y = int(sin(phi)*sin(the)*r+50)
                z = int(cos(the)*r+50)
                
                # if the value is bigger than the threshold, we include it
                if filter_vol[x,y,z] > threshold:
                    res.append(1.0)
                else:
                    res.append(0.0)
        
        # store it so that we don't have to recompute it next time
        self._sf = np.array(res)
        
        return self._sf


def create_wedge(wedgeAngle1, wedgeAngle2, cutOffRadius, sizeX, sizeY, sizeZ, smooth=0):
    '''This function returns a wedge object. For speed reasons it decides whether to generate a symmetric or assymetric wedge.
    @param wedgeAngle1: angle of wedge1 in degrees
    @type wedgeAngle1: int
    @param wedgeAngle2: angle of wedge2 in degrees
    @type wedgeAngle2: int
    @param cutOffRadius: radius from center beyond which the wedge is set to zero.
    @type cutOffRadius: int
    @param sizeX: the size of the box in x-direction.
    @type sizeX: int
    @param sizeY: the size of the box in y-direction.
    @type sizeY: int
    @param sizeZ: the size of the box in z-direction.
    @type sizeZ: int
    @param smooth: smoothing parameter that defines the amount of smoothing  at the edge of the wedge.
    @type smooth: float
    @return: 3D array determining the wedge object.
    @rtype: ndarray of np.float64'''
    print(wedgeAngle1, wedgeAngle2)
    if wedgeAngle1 == wedgeAngle2:
        return create_symmetric_wedge(wedgeAngle1, wedgeAngle2, cutOffRadius, sizeX, sizeY, sizeZ, smooth)
    else:
        return create_asymmetric_wedge(wedgeAngle1, wedgeAngle2, cutOffRadius, sizeX, sizeY, sizeZ, smooth)

def create_symmetric_wedge(angle1, angle2, cutoffRadius, sizeX, sizeY, sizeZ, smooth):
    '''This function returns a symmetric wedge object.
    @param angle1: angle of wedge1 in degrees
    @type angle1: int
    @param angle2: angle of wedge2 in degrees
    @type angle2: int
    @param cutOffRadius: radius from center beyond which the wedge is set to zero.
    @type cutOffRadius: int
    @param sizeX: the size of the box in x-direction.
    @type sizeX: int
    @param sizeY: the size of the box in y-direction.
    @type sizeY: int
    @param sizeZ: the size of the box in z-direction.
    @type sizeZ: int
    @param smooth: smoothing parameter that defines the amount of smoothing  at the edge of the wedge.
    @type smooth: float
    @return: 3D array determining the wedge object.
    @rtype: ndarray of np.float64'''

    range_angle1Smooth = smooth / np.sin(angle1 * np.pi / 180.)
    range_angle2Smooth = smooth / np.sin(angle2 * np.pi / 180.)
    wedge = np.zeros((sizeX, sizeY, sizeZ // 2 + 1), dtype=np.float64)

    z, y, x = np.meshgrid(np.abs(np.arange(-sizeX // 2 + sizeX % 2, sizeX // 2 + sizeX % 2, 1.)),
                          np.abs(np.arange(-sizeY // 2 + sizeY % 2, sizeY // 2 + sizeY % 2, 1.)),
                          np.arange(0, sizeZ // 2 + 1, 1.))
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    with np.errstate(all='ignore'):
        wedge[np.tan(angle1 * np.pi / 180) < y / x] = 1

    wedge[sizeX // 2, :, 0] = 1

    if smooth:
        area = np.abs(x - (y / np.tan(angle1 * np.pi / 180))) <= range_angle1Smooth
        strip = 1 - (np.abs(x - (y / np.tan(angle1 * np.pi / 180.))) * np.sin(angle1 * np.pi / 180.) / smooth)
        wedge += (strip * area * (1 - wedge))

    wedge[r > cutoffRadius] = 0
    return np.fft.fftshift(wedge, axes=(0, 1))

def create_asymmetric_wedge(angle1, angle2, cutoffRadius, sizeX, sizeY, sizeZ, smooth):
    '''This function returns an asymmetric wedge object.
    @param angle1: angle of wedge1 in degrees
    @type angle1: int
    @param angle2: angle of wedge2 in degrees
    @type angle2: int
    @param cutOffRadius: radius from center beyond which the wedge is set to zero.
    @type cutOffRadius: int
    @param sizeX: the size of the box in x-direction.
    @type sizeX: int
    @param sizeY: the size of the box in y-direction.
    @type sizeY: int
    @param sizeZ: the size of the box in z-direction.
    @type sizeZ: int
    @param smooth: smoothing parameter that defines the amount of smoothing  at the edge of the wedge.
    @type smooth: float
    @return: 3D array determining the wedge object.
    @rtype: ndarray of np.float64'''

    range_angle1Smooth = smooth / np.sin(angle1 * np.pi / 180.)
    range_angle2Smooth = smooth / np.sin(angle2 * np.pi / 180.)
    wedge = np.zeros((sizeX, sizeY, sizeZ // 2 + 1))

    z, y, x = np.meshgrid(np.arange(-sizeX // 2 + sizeX % 2, sizeX // 2 + sizeX % 2),
                          np.arange(-sizeY // 2 + sizeY % 2, sizeY // 2 + sizeY % 2), np.arange(0, sizeZ // 2 + 1))
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)

    wedge[np.tan(angle1 * np.pi / 180) < y / x] = 1
    wedge[np.tan(-angle2 * np.pi / 180) > y / x] = 1
    wedge[sizeX // 2, :, 0] = 1

    if smooth:
        area = np.abs(x - (y / np.tan(angle1 * np.pi / 180))) <= range_angle1Smooth
        strip = 1 - (np.abs(x - (y / np.tan(angle1 * np.pi / 180.))) * np.sin(angle1 * np.pi / 180.) / smooth)
        wedge += (strip * area * (1 - wedge) * (y > 0))

        area2 = np.abs(x + (y / np.tan(angle2 * np.pi / 180))) <= range_angle2Smooth
        strip2 = 1 - (np.abs(x + (y / np.tan(angle2 * np.pi / 180.))) * np.sin(angle2 * np.pi / 180.) / smooth)
        wedge += (strip2 * area2 * (1 - wedge) * (y <= 0))

    wedge[r > cutoffRadius] = 0
    return wedge
