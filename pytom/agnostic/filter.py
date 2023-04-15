#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
basic filters operating on numpy arrays
"""
from pytom.gpu.initialize import xp, device
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

def bandpass_circle(image, low=0, high=-1, sigma=0, ff=1):
    """Do a bandpass filter on a given volume.

    @param volume: input volume.
    @param low: low frequency in k-space.
    @param high: high frequency in k-space.
    @param sigma: smoothness factor.

    @return: bandpass filtered volume.
    """
    from pytom.agnostic.transform import fourier_filter
    assert low >= 0, "lower limit must be >= 0"

    from pytom.agnostic.tools import create_sphere, create_circle

    if high == -1:
        high = np.min(image.shape)/2
    assert low < high, "upper bandpass must be > than lower limit"

    if low == 0:
        mask = create_circle(image.shape, high, sigma, num_sigma=5)
    else:
        # BUG! TODO
        # the sigma
        mask = create_circle(image.shape, high, sigma, num_sigma=5) - create_circle(image.shape, max(0,low-sigma*2),
                                                                               sigma, num_sigma=5)

    if ff is None:
        res = applyFourierFilterFull(image, xp.fft.fftshift(mask))
    else:
        res = applyFourierFilterFull(image, xp.fft.fftshift(mask) * ff)

    return res

def bandpass(volume, low=0, high=-1, sigma=0, returnMask=False, mask=None, fourierOnly=False):
    """Do a bandpass filter on a given volume.

    @param volume: input volume.
    @param low: low frequency in k-space.
    @param high: high frequency in k-space.
    @param sigma: smoothness factor.

    @return: bandpass filtered volume.
    """
    assert low >= 0, "lower limit must be >= 0"

    from pytom.agnostic.tools import create_sphere

    if high == -1:
        high = np.min(volume.shape)/2
    assert low < high, "upper bandpass must be > than lower limit"

    if mask is None:
        if low == 0:
            mask = create_sphere(volume.shape, high, sigma)
        else:
            # BUG! TODO
            # the sigma
            mask = create_sphere(volume.shape, high, sigma) - create_sphere(volume.shape, max(0,low-sigma*2), sigma)

    from pytom.agnostic.transform import fourier_filter
    if fourierOnly:
        fvolume = xp.fft.fftshift(volume)
        sx,sy,sz = mask.shape

        fcrop = fvolume[max(0, sx // 2 - high - 2):min(sx, sx // 2 + high + 2),
                        max(0, sy // 2 - high - 2):min(sy, sy // 2 + high + 2),
                        max(0, sz // 2 - high - 2):min(sz, sz // 2 + high + 2) ]
        mcrop =    mask[max(0, sx // 2 - high - 2):min(sx, sx // 2 + high + 2),
                        max(0, sy // 2 - high - 2):min(sy, sy // 2 + high + 2),
                        max(0, sz // 2 - high - 2):min(sz, sz // 2 + high + 2)]

        res = fourierMult(xp.fft.fftshift(fcrop), mcrop, True)

    else:
        res = fourier_filter(volume, mask, True)
    if returnMask:
        return res, mask
    else:
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


def wiener_filtered_ctf(ctf, ssnr, highpass=None, phaseflipped=False):
    """
    Highpass and ssnr might be generated as follows:

    # create sampling points for highpass and snr
    d_hp = [xp.abs(xp.arange(-1, 1, 2 / size)) for size in shape]
    grids = xp.meshgrid(*d_hp)
    points_hp = xp.sqrt(sum([g ** 2 for g in grids]))
    points_ssnr = - points_hp

    # calculate highpass
    highpass = points_hp / highpassnyquist
    highpass[highpass > 1] = 1
    highpass = 1 - xp.cos(highpass * xp.pi)

    # calculate spectral SNR
    ssnr = xp.exp(points_snr * snrfalloff * 100 / spacing_angstrom) * 10 ** (3 * deconvstrength)


    @param ctf: ctf function
    @type  ctf: L{np.ndarray}
    @param ssnr: array with spectral signal to noise ratio, same shape as ctf
    @type  ssnr: L{np.ndarray}
    @param highpass: Pass optional high pass filter, should filter to ~ (0.02 nm)^-1, same shape a ctf
    @type  highpass: L{np.ndarray}
    @param phaseflipped: Flip phases of CTF
    @type  phaseflipped: L{bool}

    @return: Wiener filter
    @rtype: L{np.ndarray}
    """
    # apply highpass to ssnr
    if highpass is not None:
        ssnr *= highpass

    # flip phases if needed
    if phaseflipped:
        ctf = xp.abs(ctf)

    # create the wiener_like filter
    with np.errstate(all='ignore'):  # division by 0 is not a problem in 1/snr because the inf will be used for the
        wiener = ctf / (ctf * ctf + 1 / ssnr)  # next division
    wiener[xp.isnan(wiener)] = 0
    return wiener


def SSNR(shape, spacing_angstrom, snrfalloff, deconvstrength, r=None):
    """

    @param shape: shape tuple
    @type  shape: L{tuple}
    @param spacing_angstrom: spacing in angstrom
    @type  spacing_angstrom: L{float}
    @param snrfalloff: fall off of gaussian ssnr function
    @type  snrfalloff: L{float}
    @param deconvstrength: strength at 0 spatial frequency, 10 ** (3 * deconvstrength)
    @type  deconvstrength: L{float}
    @param r: optional radial grid
    @type  r: L{np.ndarray}

    @return: ssnr function
    @rtype:  L{np.ndarray}
    """
    # todo what about astigmatic ssnr?
    if r is None:
        d = [xp.abs(xp.arange(size) / (size // 2) - 1) for size in shape]
        grids = xp.meshgrid(*d, indexing='ij')
        r = - xp.sqrt(sum([g ** 2 for g in grids]))

    return xp.exp(r * snrfalloff * 100 / spacing_angstrom) * 10 ** (3 * deconvstrength)


def highpass_ramp(shape, highpassnyquist, r=None):
    """

    @param shape: shape tuple
    @type  shape: L{tuple}
    @param highpassnyquist: high pass frequency, good value is 0.02
    @type  highpassnyquist: L{float}
    @param r: optional radial grid
    @type  r: L{np.ndarray}

    @return: high pass ramp filter
    @rtype:  L{np.ndarray}
    """
    if r is None:
        d = [xp.abs(xp.arange(size) / (size // 2) - 1) for size in shape]
        grids = xp.meshgrid(*d, indexing='ij')
        r = xp.sqrt(sum([g ** 2 for g in grids]))
    highpass = r / highpassnyquist
    highpass[highpass > 1] = 1
    return 1 - xp.cos(highpass * xp.pi)


def wiener_like_filter_1d(shape, spacing_angstrom, defocus, snrfalloff, deconvstrength, highpassnyquist, voltage=300e3,
                  spherical_aberration=2.7e-3, amplitude_contrast=0.07, phaseflipped=False, phase_shift=0):
    """
    # todo should defocus input be in m or in um?

    @param shape:
    @type  shape:
    @param spacing_angstrom:
    @type  spacing_angstrom:
    @param defocus:
    @type  defocus:
    @param voltage:
    @type  voltage:
    @param spherical_aberration:
    @type  spherical_aberration:
    @param amplitude_contrast:
    @type  amplitude_contrast:
    @param snrfalloff:
    @type  snrfalloff:
    @param deconvstrength:
    @type  deconvstrength:
    @param highpassnyquist:
    @type  highpassnyquist:
    @param phaseflipped:
    @type  phaseflipped:
    @param phaseshift: phasehift input in degrees
    @type  phaseshift:

    @return:
    @rtype:

    @author: Marten Chaillet
    """
    from pytom.simulation.microscope import create_ctf_1d
    size = max(shape) // 2

    points_hp = xp.linspace(.0, 1., num=size)
    # points_hp = xp.arange(0, 1 + 1 / (max(shape) - 1), 1 / (max(shape) - 1))
    highpass = points_hp / highpassnyquist
    highpass[highpass > 1] = 1
    highpass = 1 - xp.cos(highpass * xp.pi)
    points_snr = xp.linspace(.0, -1., num=size)
    # points_snr = xp.arange(0, -(1 + 1 / (max(shape) - 1)), -1 / (max(shape) - 1))
    ssnr = xp.exp(points_snr * snrfalloff * 100 / spacing_angstrom) * 10 ** (3 * deconvstrength)

    ctf = - create_ctf_1d(size, spacing_angstrom * 1e-10, defocus, amplitude_contrast=amplitude_contrast,
                        voltage=voltage, Cs=spherical_aberration, phase_shift_deg=phase_shift, bfactor=.0)

    wiener = wiener_filtered_ctf(ctf, ssnr, highpass=highpass, phaseflipped=phaseflipped)

    # now expand the wiener filter to the input shape
    # x = xp.arange(0, shape[0]) - shape[0] // 2
    # y = xp.arange(0, shape[1]) - shape[1] // 2
    # z = xp.arange(0, shape[2]) - shape[2] // 2
    # xx, yy, zz = xp.meshgrid(x, y, z)
    # xx /= (shape[0] // 2)
    # yy /= (shape[1] // 2)
    # zz /= (shape[2] // 2)
    # r = xp.sqrt(xx ** 2, yy ** 2, zz ** 2)
    # r[r > 1] = 1
    # r = xp.fft.ifftshift(r)

    return profile2FourierVol(wiener, dim=shape)


def wiener_like_filter(shape, spacing_angstrom, defocus, snrfalloff, deconvstrength, highpassnyquist, voltage=300e3,
                  spherical_aberration=2.7e-3, amplitude_contrast=0.07, phaseflipped=False, phase_shift=0):
    """
    # todo should defocus input be in m or in um?
    # todo function should have an option for reduced

    @param shape:
    @type  shape:
    @param spacing_angstrom:
    @type  spacing_angstrom:
    @param defocus:
    @type  defocus:
    @param voltage:
    @type  voltage:
    @param spherical_aberration:
    @type  spherical_aberration:
    @param amplitude_contrast:
    @type  amplitude_contrast:
    @param snrfalloff:
    @type  snrfalloff:
    @param deconvstrength:
    @type  deconvstrength:
    @param highpassnyquist:
    @type  highpassnyquist:
    @param phaseflipped:
    @type  phaseflipped:
    @param phaseshift: phasehift input in degrees
    @type  phaseshift:

    @return:
    @rtype:

    @author: Marten Chaillet
    """
    from pytom.simulation.microscope import create_ctf, create_ctf_1d

    # calculate highpass
    highpass = highpass_ramp(shape, highpassnyquist)

    # calculate spectral SNR
    ssnr = SSNR(shape, spacing_angstrom, snrfalloff, deconvstrength)

    # calculate ctf
    ctf = - create_ctf(shape, spacing_angstrom * 1e-10, defocus, amplitude_contrast, voltage, spherical_aberration,
                       phase_shift_deg=phase_shift)
    print(f'ctf shape is {ctf.shape}')
    # todo add astigmatism option

    wiener = wiener_filtered_ctf(ctf, ssnr, highpass=highpass, phaseflipped=phaseflipped)

    # if reduced:
    #     z = wiener.shape[-1]
    #     wiener = wiener[...,:z//2+1]

    return wiener


class Wedge(object):
    # TODO class not used => remove
    """Defines the missing wedge filter in Fourier space."""
    def __init__(self):
        super(Wedge, self).__init__()

    def apply(self, data):
        raise NotImplementedError('Abstract method! Please overwrite it.')

    def toSphericalFunc(self, bw, radius):
        raise NotImplementedError('Abstract method! Please overwrite it.')


class GeneralWedge(Wedge):
    # TODO class not used => remove
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

        from pytom.agnostic.transform import fourier_filter
        res = fourier_filter(data, filter_vol, False)

        return res

    def toSphericalFunc(self, bw, radius):
        assert(bw<=128)

        # start sampling
        from vol2sf import vol2sf
        sf = vol2sf(self._whole_volume, radius, bw)
        
        return sf


class SingleTiltWedge(Wedge):
    # TODO class not used => remove
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
            from pytom.agnostic.transform import rotate3d, fourier_reduced2full, fourier_full2reduced, fftshift, ifftshift
            isodd = self._volume_shape[2] % 2
            filter_vol = fftshift(fourier_reduced2full(self._volume, isodd))
            filter_vol = rotate3d(filter_vol, rotation[0], rotation[1], rotation[2], order=1) # linear interp!
            filter_vol = fourier_full2reduced(ifftshift(filter_vol))
        else:
            filter_vol = self._volume

        from pytom.agnostic.transform import fourier_filter
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

        from pytom.agnostic.transform import rotate3d, fourier_reduced2full, fftshift
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
        
        from pytom.agnostic.transform import fourier_reduced2full, fftshift
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


def create_wedge(wedgeAngle1, wedgeAngle2, cutOffRadius, sizeX, sizeY, sizeZ, smooth=0, rotation=None):
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
    # TODO update so that uneven values for sizeZ still match with pytom.lib.pytom.lib.pytom_volume

    import numpy as np

    if cutOffRadius < 1:
        cutOffRadius = sizeX // 2

    if wedgeAngle1 == wedgeAngle2:
        return create_symmetric_wedge(wedgeAngle1, wedgeAngle2, cutOffRadius, sizeX, sizeY, sizeZ, smooth, rotation).astype(np.float32)
    else:
        return create_asymmetric_wedge(wedgeAngle1, wedgeAngle2, cutOffRadius, sizeX, sizeY, sizeZ, smooth, rotation).astype(np.float32)

def create_symmetric_wedge(angle1, angle2, cutoffRadius, sizeX, sizeY, sizeZ, smooth, rotation=None):
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
    wedge = xp.zeros((sizeX, sizeY, sizeZ // 2 + 1), dtype=xp.float64)
    if rotation is None:
        # numpy meshgrid by default returns indexing with cartesian coordinates (xy)
        # shape N, M, P returns meshgrid with M, N, P (see numpy meshgrid documentation)
        # the naming here is therefore weird
        z, y, x = xp.meshgrid(xp.abs(xp.arange(-sizeY // 2 + sizeY % 2, sizeY // 2 + sizeY % 2, 1.)),
                              xp.abs(xp.arange(-sizeX // 2 + sizeX % 2, sizeX // 2 + sizeX % 2, 1.)),
                              xp.arange(0, sizeZ // 2 + 1, 1.))

    else:
        # here its different again, but result is correct.
        cx,cy,cz = [s//2 for s in (sizeX,sizeY,sizeZ)]
        grid = xp.mgrid[-cx:sizeX - cx, -cy:sizeY - cy, :sizeZ // 2 + 1]

        phi, the, psi = rotation

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

        temp = grid.reshape((3, grid.size // 3))
        temp = xp.dot(Inv_R, temp)
        grid = xp.reshape(temp, grid.shape)

        y = abs(grid[0, :, :, :])
        z = abs(grid[1, :, :, :])
        x = abs(grid[2, :, :, :])

    r = xp.sqrt(x ** 2 + y ** 2 + z ** 2)
    if angle1 > 1E-3:
        range_angle1Smooth = smooth / xp.sin(angle1 * xp.pi / 180.)

        with np.errstate(all='ignore'):
            wedge[xp.tan(xp.float32(angle1) * xp.pi / xp.float32(180.)) <= y / x] = 1

        if rotation is None:
            wedge[sizeX // 2, :, 0] = 1
        else:
            phi,the,psi = rotation
            if phi < 1E-6 and psi < 1E-6 and the<1E-6:
                wedge[sizeX // 2, :, 0] = 1

        if smooth:
            area = xp.abs(x - (y / xp.tan(angle1 * xp.pi / 180))) < range_angle1Smooth
            strip = 1 - ((xp.abs((x) - ((y) / xp.tan(angle1 * xp.pi / 180.)))) * xp.sin(angle1 * xp.pi / 180.) / smooth)
            wedge += (strip * area * (1 - wedge))

    else:
        wedge += 1
    wedge[r > cutoffRadius] = 0
    return xp.fft.ifftshift(wedge, axes=(0, 1))  # TODO should be ifftshift, because centered is shifted to corner

def create_asymmetric_wedge(angle1, angle2, cutoffRadius, sizeX, sizeY, sizeZ, smooth, rotation=None):
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
    @rtype: ndarray of xp.float64'''
    range_angle1Smooth = smooth / xp.sin(angle1 * xp.pi / 180.)
    range_angle2Smooth = smooth / xp.sin(angle2 * xp.pi / 180.)
    wedge = xp.zeros((sizeX, sizeY, sizeZ // 2 + 1))

    if rotation is None:
        # see comment above with symmetric wedge function about meshgrid
        z, y, x = xp.meshgrid(xp.arange(-sizeY // 2 + sizeY % 2, sizeY // 2 + sizeY % 2),
                              xp.arange(-sizeX // 2 + sizeX % 2, sizeX // 2 + sizeX % 2),
                              xp.arange(0, sizeZ // 2 + 1))

    else:
        cx, cy, cz = [s // 2 for s in (sizeX, sizeY, sizeZ)]
        grid = xp.mgrid[-cx:sizeX - cx, -cy:sizeY - cy, :sizeZ // 2 + 1]

        phi, the, psi = rotation

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

        temp = grid.reshape((3, grid.size // 3))
        temp = xp.dot(Inv_R, temp)
        grid = xp.reshape(temp, grid.shape)

        y = grid[0, :, :, :]
        z = grid[1, :, :, :]
        x = grid[2, :, :, :]

    r = xp.sqrt(x ** 2 + y ** 2 + z ** 2)


    with np.errstate(all='ignore'):
        wedge[xp.tan(angle1 * xp.pi / 180) < y / x] = 1
        wedge[xp.tan(-angle2 * xp.pi / 180) > y / x] = 1
    wedge[sizeX // 2, :, 0] = 1

    if smooth:
        area = xp.abs(x - (y / xp.tan(angle1 * xp.pi / 180))) <= range_angle1Smooth
        strip = 1 - (xp.abs(x - (y / xp.tan(angle1 * xp.pi / 180.))) * xp.sin(angle1 * xp.pi / 180.) / smooth)
        wedge += (strip * area * (1 - wedge) * (y > 0))

        area2 = xp.abs(x + (y / xp.tan(angle2 * xp.pi / 180))) <= range_angle2Smooth
        strip2 = 1 - (xp.abs(x + (y / xp.tan(angle2 * xp.pi / 180.))) * xp.sin(angle2 * xp.pi / 180.) / smooth)
        wedge += (strip2 * area2 * (1 - wedge) * (y <= 0))

    wedge[r > cutoffRadius] = 0

    return xp.fft.ifftshift(wedge, axes=(0, 1))  # TODO should be ifftshift, because centered is shifted to corner

def circle_filter(sizeX, sizeY, radiusCutoff):
    """
    circleFilter: NEEDS Documentation
    @param sizeX: NEEDS Documentation
    @param sizeY: NEEDS Documentation
    @param radiusCutoff: NEEDS Documentation
    """
    X, Y = xp.meshgrid(xp.arange(-sizeX//2 + sizeX%2, sizeX//2+sizeX%2), xp.arange(-sizeY//2+sizeY%2, sizeY//2+sizeY%2))
    R = xp.sqrt(X**2 + Y**2)

    filter = xp.zeros((sizeX, sizeY), dtype=xp.float32)
    filter[R <= radiusCutoff] = 1

    return filter

def ellipse_filter(sizeX, sizeY, radiusCutoffX, radiusCutoffY):
    """
    circleFilter: NEEDS Documentation
    @param sizeX: NEEDS Documentation
    @param sizeY: NEEDS Documentation
    @param radiusCutoff: NEEDS Documentation
    """
    X, Y = xp.meshgrid(xp.arange(-sizeY//2+sizeY%2, sizeY//2+sizeY%2), xp.arange(-sizeX//2 + sizeX%2, sizeX//2+sizeX%2))
    R = xp.sqrt((X/radiusCutoffX)**2 + (Y/radiusCutoffY)**2)

    filter = xp.zeros((sizeX, sizeY), dtype=xp.float32)
    #print(filter.shape, R.shape)
    filter[R <= 1] = 1

    return filter

def ramp_filter(sizeX, sizeY, crowtherFreq=None, N=None):
    """
    rampFilter: Generates the weighting function required for weighted backprojection - y-axis is tilt axis

    @param sizeX: size of weighted image in X
    @param sizeY: size of weighted image in Y
    @param crowtherFreq: size of weighted image in Y
    @return: filter volume

    """
    if crowtherFreq is None: crowtherFreq = sizeX//2
    N = 0 if N is None else 1/N

    rampLine = xp.abs(xp.arange(-sizeX//2, sizeX//2)) / crowtherFreq + N
    # should be: rampLine = xp.abs(xp.arange(-sizeX // 2, sizeX // 2)) / crowtherFreq + N
    rampLine[rampLine > 1] = 1

    rampfilter = xp.column_stack(([(rampLine), ] * (sizeY)))

    return rampfilter

def exact_filter(tilt_angles, tiltAngle, sX, sY, sliceWidth, arr=[]):
    """
    exactFilter: Generates the exact weighting function required for weighted backprojection - y-axis is tilt axis
    Reference : Optik, Exact filters for general geometry three dimensional reconstuction, vol.73,146,1986.
    @param tilt_angles: list of all the tilt angles in one tilt series
    @param titlAngle: tilt angle for which the exact weighting function is calculated
    @param sizeX: size of weighted image in X
    @param sizeY: size of weighted image in Y

    @return: filter volume

    """
    # Calculate the relative angles in radians.
    sampling = xp.sin(xp.abs((xp.array(tilt_angles) - tiltAngle) * xp.pi / 180.))
    smallest_sampling = xp.min(sampling[sampling > 0.001])

    if sliceWidth / smallest_sampling > sX // 2:  # crowther crit can be to nyquist freq (i.e. sX // 2)
        sliceWidth = smallest_sampling * (sX // 2)  # adjust sliceWidth if too large

    crowtherFreq = min(sX // 2, int(xp.ceil(sliceWidth / smallest_sampling)))
    arrCrowther = xp.abs(xp.arange(-crowtherFreq, min(sX // 2, crowtherFreq + 1)))

    # as in the paper: 1 - frequency / overlap_frequency
    # where frequency = arrCrowther, and 1 / overlap_frequency = sampling/sliceWidth
    wfuncCrowther = 1. / xp.clip(1 - ((sampling / sliceWidth)[:, xp.newaxis] * arrCrowther) ** 2, 0, 2).sum(axis=0)

    # Create full width weightFunc
    wfunc = xp.ones((sX, sY), dtype=xp.float32)

    # row_stack is not implemented in cupy
    weightingFunc = xp.tile(wfuncCrowther, (sY, 1)).T

    wfunc[sX // 2 - crowtherFreq:sX // 2 + min(sX // 2, crowtherFreq + 1), :] = weightingFunc

    return wfunc

def rotateWeighting(weighting, rotation, mask=None, binarize=False):
    """
    rotateWeighting: Rotates a frequency weighting volume around the center. If the volume provided is reduced complex, it will be rescaled to full size, ftshifted, rotated, iftshifted and scaled back to reduced size.
    @param weighting: A weighting volume in reduced complex convention
    @type weighting: cupy or numpy array
    @param rotation: rotation angles in zxz order
    @type rotation: list
    @param mask:=None is there a rotation mask? A mask with all = 1 will be generated otherwise. Such mask should be \
        provided anyway.
    @type mask: cupy or numpy ndarray
    @return: weight as reduced complex volume
    @rtype: L{pytom.lib.pytom_volume.vol_comp}
    """
    from pytom.lib.pytom_volume import vol, limit, vol_comp
    from pytom.lib.pytom_volume import rotate
    from pytom.voltools import transform
    assert type(weighting) == vol or type(weighting) == vol_comp, "rotateWeighting: input neither vol nor vol_comp"
    from pytom.agnostic.transform import fourier_reduced2full, fourier_full2reduced

    weighting = fourier_reduced2full(weighting, isodd=weighting.shape[0]%2 == 1)
    weighting = xp.fft.fftshift(weighting)

    weightingRotated = xp.zeros_like(weighting)

    transform(weighting, output=weightingRotated, rotation=rotation, rotation_order='rzxz', device=device, interpolation='filt_bspline')

    if not mask is None:
        weightingRotated *= mask


    weightingRotated = xp.fft.fftshift(weightingRotated)
    returnVolume = fourier_full2reduced(weightingRotated)

    if binarize:
        returnVolume[returnVolume < 0.5] = 0
        returnVolume[returnVolume >= 0.5] = 1

    return returnVolume

def profile2FourierVol(profile, dim=None, reduced=False):
    """
    create Volume from 1d radial profile, e.g., to modulate signal with \
    specific function such as CTF or FSC. Simple linear interpolation is used\
    for sampling.

    @param profile: profile
    @type profile: 1-d L{pytom.lib.pytom_volume.vol} or 1-d python array
    @param dim: dimension of (cubic) output
    @type dim: L{int}
    @param reduced: If true reduced Fourier representation (N/2+1, N, N) is generated.
    @type reduced: L{bool}

    @return: 3-dim complex volume with spherically symmetrical profile
    @rtype: L{pytom.lib.pytom_volume.vol}
    @author: FF
    """

    if dim is None:
        try:
            dim = [2 * profile.shape[0],]*3
        except:
            dim = [2 * len(profile),]*3

    is3D = (len(dim) ==  3)

    nx,ny = dim[:2]
    if reduced:
        if is3D:
            nz = int(dim[2] // 2) + 1
        else:
            ny = int(ny // 2) + 1
    else:
        if is3D:
            nz = dim[2]



    try:
        r_max = profile.shape[0] - 1
    except:
        r_max = len(profile) - 1

    if len(dim) ==3:
        if reduced:
            X, Y, Z = xp.meshgrid(xp.arange(-nx // 2, nx // 2 + nx % 2), xp.arange(-ny // 2, ny // 2 + ny % 2),
                               xp.arange(0, nz))
        else:
            X, Y, Z = xp.meshgrid(xp.arange(-nx // 2, nx // 2 + nx % 2), xp.arange(-ny // 2, ny // 2 + ny % 2),
                               xp.arange(-nz // 2, nz // 2 + nz % 2))
        R = xp.sqrt(X ** 2 + Y ** 2 + Z ** 2)

    else:
        if reduced:
            X, Y = xp.meshgrid(xp.arange(-nx // 2, ny // 2 + ny % 2), xp.arange(0, ny))
        else:
            X, Y = xp.meshgrid(xp.arange(-nx // 2, nx // 2 + nx % 2), xp.arange(-ny // 2, ny // 2 + ny % 2))
        R = xp.sqrt(X ** 2 + Y ** 2 )

    IR = xp.floor(R).astype(xp.int64)
    valIR_l1 = IR.copy()
    valIR_l2 = valIR_l1 + 1
    val_l1, val_l2 = xp.zeros_like(X, dtype=xp.float64), xp.zeros_like(X, dtype=xp.float64)

    l1 = R - IR.astype(xp.float32)
    l2 = 1 - l1

    try:
        profile = xp.array(profile)
    except:
        import numpy
        profile = xp.array(numpy.array(profile))

    for n in xp.arange(r_max):
        val_l1[valIR_l1 == n] = profile[n]
        val_l2[valIR_l2 == n + 1] = profile[n + 1]

    val_l1[IR == r_max] = profile[n + 1]
    val_l2[IR == r_max] = profile[n + 1]

    val_l1[R > r_max] = 0
    val_l2[R > r_max] = 0

    fkernel = l2 * val_l1 + l1 * val_l2

    if reduced:
        fkernel = xp.fft.fftshift(fkernel, axes=(0, 1))
    else:
        fkernel = xp.fft.fftshift(fkernel)

    return fkernel

def filter_volume_by_profile(volume, profile):
    """
    filter volume by 1-d profile
    @param volume: volume
    @type volume: L{pytom.lib.pytom_volume.vol}
    @param profile: 1-d profile
    @type profile: L{pytom.lib.pytom_volume.vol}
    @return: outvol
    @rtype: L{pytom.lib.pytom_volume.vol}
    @author: FF
    """
    from pytom.agnostic.filter import applyFourierFilter, applyFourierFilterFull

    if volume.shape[0] != volume.shape[-1]:
        reduced = True
        convolute = applyFourierFilter
    else:
        reduced = False
        convolute = applyFourierFilterFull

    kernel = profile2FourierVol(profile=profile, dim=volume.shape, reduced=reduced)
    outvol = convolute(volume, kernel)
    return outvol

def applyFourierFilter(particle, filter):
    return xp.fft.irfftn(xp.fft.rfftn(particle, particle.shape) * filter).real.astype(xp.float32)

def applyFourierFilterFull(particle, filter):
    return xp.fft.ifftn(xp.fft.fftn(particle) * filter).real.astype(xp.float32)

def fourierMult(fvolume, filter, human=False):
    from pytom.agnostic.transform import fourier_full2reduced

    if human:
        filter = xp.fft.fftshift(filter)
        filter = fourier_full2reduced(filter)
        fvolume = fourier_full2reduced(fvolume)

    fvolume *= filter

    return fvolume
