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

    input = xp.fft.fftn(input)
    # 2.35 is the factor for the full width half maximum
    result = fourier_gaussian(input, sigma=(resolution / (2 * spacing)) / 2.35)
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

