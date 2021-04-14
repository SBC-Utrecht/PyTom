import numpy as xp


def create_gaussian_low_pass(shape, cutoff, center=None):
    """
    NOTE: MOVE FUNCTION TO TOMPY.FILTER
    Create a 3D mask with gaussian edges.

    @param shape: shape tuple.
    @param hwhm: cutoff value for the filter.
    @param center: center of the Gaussian. If not provided will be calculated.

    @return: sphere inside a volume.
    """
    assert len(shape) == 2 or len(shape) == 3, "filter can only be created in 2d or 3d"
    assert len(set(shape)) == 1, "shape needs to have equal sizes"

    # full width at half maximum is two times the half width at half maximum (or cutoff)
    c = 2  # or can be set to np.sqrt(2) for butterworth filter
    sigma_cutoff = cutoff / xp.sqrt(2 * xp.log(c))

    if center is None:
        # center = [(shape[0]-1)/2, (shape[1]-1)/2, (shape[2]-1)/2]
        center = [s//2 for s in shape]

    grid = xp.mgrid[0:shape[0], 0:shape[1], 0:shape[2]] if len(shape) == 3 else xp.mgrid[0:shape[0], 0:shape[1]]
    # r = xp.sqrt(sum((x-center[0])**2+(y-center[1])**2+(z-center[2])**2)
    r = xp.sqrt(sum([(g-c)**2 for (g,c) in zip(grid, center)]))

    filter = xp.exp(-r ** 2 / (2 * sigma_cutoff ** 2))

    return filter


def reduce_resolution_fourier(input, spacing, resolution):
    """
    NOTE: MOVE FUNCTION TO TOMPY.FILTER
    Apply scipy gaussian filter in fourier space.

    @param input:
    @param spacing:
    @param resolution:

    @return: filtered volume
    """
    from scipy.ndimage import fourier_gaussian

    # 2.35 is the factor for the full width half maximum
    result = fourier_gaussian(xp.fft.fftn(input), sigma=(resolution / (2 * spacing)) / 2.35)
    return xp.fft.ifftn(result).real


def reduce_resolution(input, spacing, resolution):
    """
    NOTE: MOVE FUNCTION TO TOMPY.FILTER
    Apply scipy gaussian filter in real space.

    @param input:
    @param spacing:
    @param resolution:

    @return: filtered volume
    """
    from scipy.ndimage import gaussian_filter

    # 2.35 is the factor for the full width half maximum
    result = gaussian_filter(input, sigma=(resolution / (2 * spacing)) / 2.35)
    return result


def gradient_image(size, factor, angle=0, center_shift=0):
    """
    Creates an image with a gradient of values rotated along angle. Factor determines the strength of the gradient.
    @param size:
    @param factor:
    @param angle:
    @param shift:
    @return:
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


def create_circle(size, radius=-1, sigma=0, center=None):
    """
    Create a sphere in image of size with radius.
    """
    if len(size) == 1:
        size = (size, size)
    assert len(size) == 2

    if center is None:
        center = [size[0] / 2, size[1] / 2]
    if radius == -1:
        radius = xp.min(size) / 2

    sphere = xp.zeros(size)
    [x, y] = xp.mgrid[0:size[0], 0:size[1]]
    r = xp.sqrt((x - center[0]) ** 2 + (y - center[1]) ** 2)
    sphere[r <= radius] = 1

    if sigma > 0:
        ind = xp.logical_and(r > radius, r < radius + 2 * sigma)
        sphere[ind] = xp.exp(-((r[ind] - radius) / sigma) ** 2 / 2)

    return sphere


def bandpass_mask(shape, low=0, high=-1):
    """
    Return 2d bandpass mask.
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
    TODO change name to bin_volume instead of bin
    @param potential:
    @param factor: integer value
    @return:
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
    Generate an ellipse defined by 3 radii along x,y,z - parameters mj, mn1, mn2.
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
    Correlated noise to create density deformations.
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
    add Adds white noise to volume
    @param volume: A volume
    @param SNR: Signal to Noise ratio of result
    @type SNR: int or float > 0
    @return: Volume containing noise with SNR == SNR
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
