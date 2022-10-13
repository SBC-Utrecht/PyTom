from pytom.gpu.initialize import xp, device


def fwhm_to_sigma(fwhm):
    return fwhm / (2 * xp.sqrt(2 * xp.log(2)))


def sigma_to_fwhm(sigma):
    return sigma * (2 * xp.sqrt(2 * xp.log(2)))


def hwhm_to_sigma(hwhm):
    return hwhm / (xp.sqrt(2 * xp.log(2)))


def sigma_to_hwhm(sigma):
    return sigma * (xp.sqrt(2 * xp.log(2)))


def create_gaussian_low_pass(shape, spacing, resolution, reduced=False):
    """
    NOTE: MOVE FUNCTION TO TOMPY.FILTER
    Create a 2D or 3D Gaussian low-pass filter with cutoff (or HWHM). This value will be converted to the proper
    sigma for the Gaussian function to acquire the desired cutoff.

    @param shape: shape tuple with x,y or x,y,z dimension
    @type  shape: L{tuple} -> (L{int},) * 3 or L{tuple} -> (L{int},) * 2
    @param hwhm: half width at half maximum for gaussian low-pass filter, the cutoff value for the filter
    @type  hwhm: L{int}
    @param center: center of the Gaussian, will be calculated if not provided
    @type  center: L{tuple} -> (L{int},) * 3 or L{tuple} -> (L{int},) * 2

    @return: sphere/circle in square volume/image, 3d or 2d array dependent on input shape
    @rtype:  L{np.ndarray}

    @author: Marten Chaillet
    """
    from pytom.simulation.microscope import normalised_grid
    assert len(shape) == 2 or len(shape) == 3, "filter can only be created in 2d or 3d"

    # 2 * spacing / resolution is cutoff in fourier space
    # then convert cutoff (hwhm) to sigma for gaussian function
    sigma_fourier = hwhm_to_sigma(2 * spacing / resolution)

    return xp.exp(-normalised_grid(shape, reduced=reduced) ** 2 / (2 * sigma_fourier ** 2))


def reduce_resolution_fourier(input, spacing, resolution):
    """
    NOTE: MOVE FUNCTION TO agnostic.FILTER
    Apply scipy gaussian filter in fourier space.

    @param input: input to be filtered, either 2d or 3d array (however scipy will be able to handle higher
    dimensionality as well
    @type  input: L{np.ndarray}
    @param spacing: spacing of each pixel/voxel in relative units
    @type  spacing: L{float}
    @param resolution: desired resolution after filtering. maximal resolution is 2 * spacing, and thus resolution
    value of 2 * spacing will not filter the image
    @type  resolution: L{float}

    @return: filtered input, 2d or 3d array
    @rtype:  L{np.ndarray}

    @author: Marten Chaillet
    """
    gaussian_filter = create_gaussian_low_pass(input.shape, spacing, resolution)
    if 'gpu' in device:
        from pytom.agnostic.filter import applyFourierFilterFull
        return applyFourierFilterFull(input, xp.fft.ifftshift(gaussian_filter))
    else:
        from pytom.agnostic.transform import fourier_filter
        return fourier_filter(input, gaussian_filter, human=True)


def reduce_resolution_real(input, spacing, resolution):
    """
    NOTE: MOVE FUNCTION TO agnostic.FILTER
    Apply scipy gaussian filter in real space.

    @param input: input to be filtered, either 2d or 3d array (however scipy will be able to handle higher
    dimensionality as well
    @type  input: L{np.ndarray}
    @param spacing: spacing of each pixel/voxel in relative units
    @type  spacing: L{float}
    @param resolution: desired resolution after filtering. maximal resolution is 2 * spacing, and thus resolution
    value of 2 * spacing will not filter the image
    @type  resolution: L{float}

    @return: filtered input, 2d or 3d array
    @rtype:  L{np.ndarray}

    @author: Marten Chaillet
    """
    assert 'cpu' in device, 'reducing resolution with gaussian kernel is done through scipy which for now only works ' \
                            'on numpy arrays, and not cupy arrays'
    #TODO might be done in newer cupy versions via cupyx.scipy

    from scipy.ndimage import gaussian_filter

    # here it needs to be fwhm to sigma, while in fourier space hwhm to sigma
    return gaussian_filter(input, sigma=fwhm_to_sigma(resolution / (2 * spacing)))


def gradient_image(size, factor, angle=0, center_shift=0):
    """
    Creates an image with a gradient of values rotated along angle. Factor determines the strength of the gradient.

    @param size: size of the image, x and y size are equal
    @type  size: L{int}
    @param factor: strength of gradient, value between 0 and 1, where 0 is no gradient and 1 a gradient from 0 to 2
    @type  factor L{float}
    @param angle: angle to rotate the gradient by
    @type  angle: L{float}
    @param center_shift: whether to shift the gradient from perfect center, if no shift applied center value will
    always be equal to 1
    @type  center_shift: L{float}

    @return: image, a 2d array of floats
    @rtype:  L{np.ndarray}

    @author: Marten Chaillet
    """
    from scipy.ndimage import rotate
    max_rotation_radius = (size/2) / xp.cos(45 * xp.pi / 180)
    extension = int(xp.ceil(max_rotation_radius - size/2))
    left = 1 - factor
    right = 1 + factor
    step = (right-left) / size
    values = xp.arange(left - extension * step + center_shift * step,
                       right + extension * step + center_shift * step, step)
    image = xp.repeat(values[xp.newaxis, :], size + 2*extension, axis=0)
    return rotate(image, angle, reshape=False)[extension:size+extension, extension:size+extension]


def create_circle(shape, radius=-1, sigma=0, center=None):
    """
    Create a circle with radius in an image of shape.

    @param shape: shape of the image
    @type  shape: L{tuple} -> (L{int},) * 2
    @param radius: radius of the circle
    @type  radius: L{float}
    @param sigma: smooth gaussian edge of circle
    @type  sigma: L{float}
    @param center: center of the circle in the image
    @type  center: L{tuple} -> (L{int},) * 2

    @return: image with circle, 2d array of floats
    @rtype:  L{np.ndarray}

    @author: Marten Chaillet
    """
    assert len(shape) == 2

    if center is None:
        center = [shape[0] / 2, shape[1] / 2]
    if radius == -1:
        radius = xp.min(shape) / 2

    sphere = xp.zeros(shape)
    [x, y] = xp.mgrid[0:shape[0], 0:shape[1]]
    r = xp.sqrt((x - center[0]) ** 2 + (y - center[1]) ** 2)
    sphere[r <= radius] = 1

    if sigma > 0:
        ind = xp.logical_and(r > radius, r < radius + 2 * sigma)
        sphere[ind] = xp.exp(-((r[ind] - radius) / sigma) ** 2 / 2)

    return sphere


def bandpass_mask(shape, low=0, high=-1):
    """
    Return 2d bandpass mask in shape. Mask is created by subtracting a circle with radius low from a circle with
    radius high.

    @param shape: shape of image
    @type  shape: L{tuple} -> (L{int},) * 2
    @param low: inner radius of band
    @type  low: L{float}
    @param high: outer radius of band
    @type  high: L{float}

    @return: an image with bandpass, 2d array of ints
    @rtype:  L{np.ndarray}

    @author: Marten Chaillet
    """
    assert low >= 0, "lower limit must be >= 0"

    if high == -1:
        high = xp.min(shape) / 2
    assert low < high, "upper bandpass must be > than lower limit"

    if low == 0:
        mask = create_circle(shape, high, sigma=0)
    else:
        mask = create_circle(shape, high, sigma=0) - create_circle(shape, low, sigma=0)

    return mask


def bin_volume(potential, factor):
    """
    Bin the input volume (potential) factor times.

    @param potential: input volume, 3d array
    @type  potential: L{np.ndarray}
    @param factor: integer multiple of 1, number of times to bin (or downsample)
    @type  factor: L{int}

    @return: downsampled input, 3d array
    @rtype:  L{np.ndarray}

    @author: Marten Chaillet
    """
    assert type(factor) is int and factor >= 1, print('non-valid binning factor, should be integer above 1')

    if factor == 1:
        return potential

    size = potential.shape
    s = [(x % factor) // 2 for x in size]
    d = [(x % factor)  % 2 for x in size]
    # print(size, s, d)
    # s = (potential.shape[0] % factor) // 2
    # d = (potential.shape[0] % factor) % 2

    potential = potential[s[0]:size[0] - s[0] - d[0], s[1]:size[1] - s[1] - d[1], s[2]:size[2] - s[2] - d[2]]
    # potential = potential[s:potential.shape[0] - s - d, s:potential.shape[0] - s - d, s:potential.shape[0] - s - d]

    size = potential.shape if potential.shape != size else size
    # ds = int(potential.shape[0]//factor)
    ds = [int(x // factor) for x in size]
    # image_size = potential.shape[0]

    # binned = potential.reshape(ds, image_size // ds,
    #                            ds, image_size // ds, ds, image_size // ds).mean(-1).mean(1).mean(-2)
    binned = potential.reshape(ds[0], size[0] // ds[0], ds[1], size[1] // ds[1], ds[2], size[2] // ds[2]).mean(-1).mean(1).mean(-2)

    return binned


def create_ellipse(size, mj, mn1, mn2, smooth=0, cutoff_SD=3):
    """
    Generate an ellipse defined by 3 radii along x,y,z - parameters mj, mn1, mn2. Ellipse is generated in a square
    volume with each dimension has same size.

    @param size: length of dimensions
    @type  size: L{int}
    @param mj: major radius
    @type  mj: L{float}
    @param mn1: minor radius 1
    @type  mn1: L{float}
    @param mn2: minor radius 2
    @type  mn2: L{float}
    @param smooth: gaussian smoothing of ellips edges, where smooth is the sigma of gaussian function
    @type  smooth: L{float}
    @param cutoff_SD: number of standard deviations to determine gaussian smoothing to
    @type  cutoff_SD: L{int}

    @return: square volume with ellipse, 3d array of floats
    @rtype:  L{np.ndarray}

    @author: Marten Chaillet
    """
    X,Y,Z = xp.meshgrid(xp.arange(size/1), xp.arange(size/1), xp.arange(size/1))

    X -= size/2-0.5
    Y -= size/2-0.5
    Z -= size/2-0.5

    R = xp.sqrt( (X/mj)**2 + (Y/mn1)**2 + (Z/mn2)**2)

    # print(R.max(), R.min())

    out = xp.zeros((size,size,size),dtype=xp.float32)
    out[ R <= 1] = 1

    if smooth:
        R2 = R.copy()
        R2[R <= 1] = 1
        sphere = xp.exp(-1 * ((R2-1)/smooth)**2)
        sphere[sphere <= xp.exp(-cutoff_SD**2/2.)] = 0
        out = sphere

    return out


def add_correlated_noise(noise_size, dim):
    """
    Add correlated noise to create density deformations.

    @param noise_size: strength of the correlation in noise, in number of pixels
    @type  noise_size: L{int}
    @param dim: dimension of the volume to create the correlated noise in
    @type  dim: L{int}

    @return: volume with noise, 3d array of floats
    @rtype:  L{np.ndarray}

    @author: Marten Chaillet
    """
    from numpy.fft import ifftn, fftn, fftshift
    from numpy.random import random

    if noise_size == 0:
        # 0 value will not work
        noise_size = 1

    noise_no_norm = abs(ifftn(fftshift(fftn(random((noise_size, noise_size, noise_size)))), [dim] * 3))
    noise = 0.2 * noise_no_norm / abs(noise_no_norm).max()

    return 1 + (noise - noise.mean())


def add_white_noise(volume, SNR=1):
    """
    Adds white (normal distributed) noise to a volume.

    @param volume: the volume to be noise, 3d array of floats
    @type  volume: L{np.ndarray}
    @param SNR: signal to noise ratio of output, assuming input to contain no noise
    @type  SNR: L{float}

    @return: noisy input volume with specified SNR, 3d array of floats
    @rtype:  L{np.ndarray}

    @author: Thomas Hrabe
    """

    if (SNR < 0):
        return volume

    from math import sqrt
    from pytom_volume import vol, mean, variance, gaussianNoise

    m = mean(volume)
    s = sqrt(variance(volume, False) / SNR)  # SNR = Var(signal) / Var(noise)

    #    s = sqrt(variance(volume,False)/SNR)-variance(volume,False)
    #
    #    if s<0:
    #        return volume
    #    elif s==0:
    #        s =1

    noise = vol(volume.sizeX(), volume.sizeY(), volume.sizeZ())

    gaussianNoise(noise, m, s)  # s is actually the std

    result = volume + noise

    return result
