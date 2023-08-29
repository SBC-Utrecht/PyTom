"""This module defines agnostic correlation functions

Created: unknown date
@author: FF

restructured by @sroet in May 2023
"""
from pytom.gpu.initialize import xp
from pytom.agnostic.filter import bandpass
from pytom.agnostic.io import read
from pytom.agnostic.macros import volumesSameSize
from pytom.agnostic.normalise import (
    meanUnderMask,
    meanVolUnderMask,
    stdUnderMask,
    stdVolUnderMask,
    normaliseUnderMask,
    mean0std1,
)
from pytom.agnostic.tools import create_circle, create_sphere
from pytom.basic.structures import Mask

# Typing imports
from typing import Optional, Tuple, Any, List, Union, cast
from pytom.gpu.initialize import xpt


def flcf(
    volume: xpt.NDArray[float],
    template: xpt.NDArray[float],
    mask: Optional[xpt.NDArray[float]] = None,
    std_v: Optional[float] = None,
) -> xpt.NDArray[float]:
    """Fast local correlation function

    @param volume: target volume
    @param template: template to be searched.
                     It can have a smaller size then target volume.
    @param mask: template mask. If not given, a default sphere mask will be used.
    @param std_v: standard deviation of the target volume under mask, which do not need
                  to be calculated again when the mask is identical.

    @return: the local correlation function
    """

    if mask is None:
        radius = volume.shape[0] // 2 - 3
        if len(volume.shape) == 2:
            mask = create_circle(volume.shape, radius, 3)
        elif len(volume.shape) == 3:
            mask = create_sphere(volume.shape, radius, 3)

    p = mask.sum()
    mean_t = meanUnderMask(template, mask, p=p)
    temp = ((template - mean_t) / stdUnderMask(template, mask, mean_t, p=p)) * mask

    if std_v is None:
        mean_v = meanVolUnderMask(volume, mask)
        std_v = stdVolUnderMask(volume, mask, mean_v)

    res = (
        xp.fft.fftshift(
            xp.fft.ifftn(xp.conj(xp.fft.fftn(temp)) * xp.fft.fftn(volume))
        ).real
        / std_v
        / p
    )
    return res


def xcc(
    volume: xpt.NDArray[float],
    template: xpt.NDArray[float],
    mask: Optional[xpt.NDArray[float]] = None,
    volume_is_normalized: bool = False,
) -> float:
    """
    xcc: Calculates the cross correlation coefficient in real space
    @param volume: A volume
    @type volume:  L{xpt.NDArray}
    @param template: A template that is searched in volume.
                     Must have the same size as volume.
    @type template:  L{xpt.NDArray}
    @param mask: mask to constrain correlation
    @type mask: L{xpt.NDArray}
    @param volume_is_normalized: only used for compatibility with nxcc - not used
    @type volume_is_normalized: L{bool}
    @return: An unscaled value
    @raise exception: Raises a runtime error if volume and template have a different
                      size.
    @author: Thomas Hrabe
    """

    if not volumesSameSize(volume, template):
        raise RuntimeError("Volume and template must have same size!")

    if mask:  # mask is given
        result = mask * volume * template
    else:
        result = volume * template

    cc: float = result.sum()

    cc = cc / float(volume.size)

    return cc


def nxcc(
    volume: xpt.NDArray[float],
    template: xpt.NDArray[float],
    mask: Optional[xpt.NDArray[float]] = None,
    volume_is_normalized: bool = False,
) -> float:
    """
    nxcc: Calculates the normalized cross correlation coefficient in real space
    @param volume: A volume
    @type volume:  L{xpt.NDArray}
    @param template: A template that is searched in volume.
                     Must have the same size as volume.
    @type template:  L{xpt.NDArray}
    @param mask: mask to constrain correlation
    @type mask: L{xpt.NDArray}
    @param volume_is_normalized: speed up if volume is already normalized
    @type volume_is_normalized: L{bool}
    @return: A value between -1 and 1
    @raise exception: Raises a runtime error if volume and template have a different
                      size.
    @author: Thomas Hrabe
    @change: flag for pre-normalized volume, FF
    """

    if not volumesSameSize(volume, template):
        raise RuntimeError("Volume and template must have same size!")

    if mask is None:
        if not volume_is_normalized:
            v = mean0std1(volume, True)
        else:
            v = volume
        t = mean0std1(template, True)
        p = volume.size
        result = v * t
    else:
        from pytom.agnostic.normalise import normaliseUnderMask

        if not volume_is_normalized:
            (v, p) = normaliseUnderMask(volume, mask)
            (t, p) = normaliseUnderMask(template, mask, p)
            t = t * mask  # multiply with the mask
            result = v * t
        else:
            (t, p) = normaliseUnderMask(template, mask)
            t = t * mask  # multiply with the mask
            result = volume * t

    ncc = result.sum()
    ncc = ncc / float(p)

    return float(ncc)


def xcf(
    volume: Union[xpt.NDArray[float], xpt.NDArray[complex]],
    template: Union[xpt.NDArray[float], xpt.NDArray[complex]],
    mask: Any = None,
    std_v: Any = None,
) -> xpt.NDArray[float]:
    """
    XCF: returns the non-normalised cross correlation function. The xcf
    result is scaled only by the square of the number of elements.

    @param volume : The search volume
    @type volume: L{xpt.NDArray}
    @param template : The template searched (this one will be used for conjugate complex
                                             multiplication)
    @type template: L{xpt.NDArray}
    @param mask: Will be unused, only for compatibility reasons with FLCF
    @param std_v: Will be unused, only for compatibility reasons with FLCF
    @return: XCF volume
    @rtype: L{xpt.NDArray}
    @author: Thomas Hrabe
    """

    # if mask:
    #    volume = volume * mask
    # 	 template = template * mask

    # determine fourier transforms of volumes
    if template.dtype in (xp.float64, xp.float32):
        fvolume = xp.fft.fftn(volume)
    else:
        fvolume = volume

    if template.dtype in (xp.float32, xp.float64):
        ftemplate = xp.fft.fftn(template)
    else:
        ftemplate = template

    # perform element wise - conjugate multiplication
    ftemplate = xp.conj(ftemplate)
    fresult = fvolume * ftemplate

    # transform back to real space
    result = abs(xp.fft.ifftn(fresult))

    # xp.fft.iftshift(result)

    return result


def xcf_mult(
    volume: xpt.NDArray[float],
    template: xpt.NDArray[float],
    mask: Any = None,
    std_v: Any = None,
) -> xpt.NDArray[float]:
    """
    XCF: returns the non-normalised cross correlation function. The xcf
    result is scaled only by the square of the number of elements.

    @param volume : The search volume
    @type volume: L{pytom.lib.pytom_volume.vol}
    @param template : The template searched (this one will be used for conjugate complex
                                             multiplication)
    @type template: L{pytom.lib.pytom_volume.vol}
    @param mask: Will be unused, only for compatibility reasons with FLCF
    @param std_v: Will be unused, only for compatibility reasons with FLCF
    @return: XCF volume
    @rtype: L{pytom.lib.pytom_volume.vol}
    @author: Thomas Hrabe
    """

    # if mask:
    #    volume = volume * mask
    # 	 template = template * mask

    # determine fourier transforms of volumes
    fvolume = xp.fft.fftn(volume)
    ftemplate = xp.fft.fftn(template)

    # perform element wise - conjugate multiplication
    ftemplate = xp.conj(ftemplate)
    fresult = fvolume * ftemplate

    # transform back to real space
    result = abs(xp.fft.ifftn(fresult))

    # xp.fft.iftshift(result)

    return result


def norm_xcf(
    volume: xpt.NDArray[float],
    template: xpt.NDArray[float],
    mask: Optional[xpt.NDArray[float]] = None,
    std_v: Optional[float] = None,
    gpu: bool = False,
) -> xpt.NDArray[float]:
    """
    nXCF: returns the normalised cross correlation function. Autocorrelation
    of two equal objects would yield a max nxcf peak of 1.

    @param volume: The search volume
    @param template: The template searched (this one will be used for conjugate complex
                                            multiplication)
    @type template: L{xpt.NDArray}
    @param mask: template mask. If not given, a default sphere mask will be generated
                 which has the same size with the given template.
    @type mask: L{xpt.NDArray}
    @param std_v: Will be unused, only for compatibility reasons with FLCF
    @return: the calculated norm_xcf volume
    @rtype: L{xpt.NDArray}
    @author: Thomas Hrabe
    @change: masking of template implemented
    """

    if mask is not None:
        result = xcf_mult(
            normaliseUnderMask(volume=volume, mask=mask, p=None),
            normaliseUnderMask(volume=template, mask=mask, p=None),
            mask,
        )
    else:
        result = xcf_mult(
            mean0std1(volume, True), mean0std1(template, True), mask, True
        )

    return result


# Band correlation


def band_cc(
    volume: xpt.NDArray[float],
    reference: xpt.NDArray[float],
    band: Tuple[float, float],
    verbose: bool = False,
    shared: None = None,
    index: None = None,
) -> float:
    """
    band_cc: Determines the normalised correlation coefficient within a band
    @param volume: The volume
    @type volume: L{xpt.NDArray}
    @param reference: The reference
    @type reference: L{xpt.NDArray}
    @param band: [a,b] - specify the lower and upper end of band.
    @return: The correlation coefficient of the two volumes in the specified band.
    @rtype: L{float}
    @author: GS
    """

    if index is not None:
        print(index)

    # from pytom.agnostic.filter import vol_comp

    if verbose:
        print("lowest freq : ", band[0], " highest freq", band[1])

    vf, m = bandpass(volume, band[0], band[1], returnMask=True, fourierOnly=True)
    rf: xpt.NDArray = bandpass(
        reference, band[0], band[1], mask=m, fourierOnly=True
    )  # ,vf[1])

    vf = vf.astype(xp.complex128)
    cc_volume = rf.astype(vf.dtype)

    cc_volume = cc_volume * xp.conj(vf)

    cc: float = cc_volume.sum().real

    v = vf
    r = rf

    abs_v = xp.abs(v)
    abs_r = xp.abs(r)

    sum_v = xp.sum(abs_v**2)
    sum_r = xp.sum(abs_r**2)

    sum_v = xp.abs(sum_v)
    sum_r = xp.abs(sum_r)

    if sum_v == 0:
        sum_v = 1

    if sum_r == 0:
        sum_r = 1

    cc = cc / (xp.sqrt(sum_v * sum_r))

    # numerical errors will be punished with nan
    nan_treshold = 1.1
    if abs(cc) > nan_treshold:
        cc = float("nan")

    return float(cc)


# Fourier Shell Correlation and helper functions


def fsc(
    volume1: xpt.NDArray[float],
    volume2: xpt.NDArray[float],
    number_bands: Optional[int] = None,
    mask: Optional[Union[xpt.NDArray[float], Mask, str]] = None,
    verbose: bool = False,
    filename: Optional[str] = None,
) -> List[float]:
    """
    FSC - Calculates the Fourier Shell Correlation for two volumes
    @param volume1: volume one
    @type volume1: L{xpt.NDArray}
    @param volume2: volume two
    @type volume2: L{xpt.NDArray}
    @param number_bands: number of shells for FSC
    @type number_bands: int
    @param mask: mask
    @type mask: L{xpt.NDArray}
    @param verbose: flag to activate printing of info
    @type verbose: L{bool}
    @param filename: write FSC to ascii file if specified
    @type filename: string

    @return: Returns a list of cc values
    @author:GS
    @rtype: list[floats]
    """
    import time

    time.time()

    if not volumesSameSize(volume1, volume2):
        raise RuntimeError("Volumes must have the same size!")

    number_bands = volume1.shape[0] // 2 if number_bands is None else number_bands

    if mask is not None:
        if mask.__class__ == xp.array([0]).__class__:
            volume1 = volume1 * mask
            volume2 = volume2 * mask

        elif isinstance(mask, Mask):
            mask = mask.getVolume()
            volume1 = volume1 * mask
            volume2 = volume2 * mask

        elif isinstance(mask, str):
            mask = read(mask)
            volume1 = volume1 * mask
            volume2 = volume2 * mask

        else:
            raise RuntimeError(
                "FSC: Mask must be a volume OR a Mask object OR a string path to a mask"
            )

    fsc_result = []

    increment = int(volume1.shape[0] / 2 * 1 / number_bands)
    import time

    fvolume1 = xp.fft.fftn(volume1)
    fvolume2 = xp.fft.fftn(volume2)

    for n, i in enumerate(range(0, volume1.shape[0] // 2, increment)):
        band_l = i
        band_r = i + increment
        band = (band_l, band_r)
        if verbose:
            print("Band : ", band)

        res = band_cc(fvolume1, fvolume2, band, verbose)

        if i == 0 and increment == 1:
            # force a 1 for correlation of the zero frequency
            res = 1

        if verbose:
            print("Correlation ", res)

        fsc_result.append(res)

    if filename:
        f = open(filename, "w")
        for item in fsc_result:
            f.write("%s\n" % item)
        f.close()

    return fsc_result


def determine_resolution(
    fsc: List[float],
    resolution_criterion: float,
    verbose: bool = False,
    randomized_fsc: Optional[List[float]] = None,
) -> Tuple[float, float, int]:
    """
    determine_resolution: Determines frequency and band where correlation drops below
                          the resolution_criterion. Uses linear interpolation between
                          two positions
    @param fsc: The fsc list determined by L{pytom.basic.correlation.FSC}
    @param resolution_criterion: A value between 0 and 1
    @param verbose: Bool that activate writing of info, default=False
    @param randomized_fsc: A value that sets the start of the calculation of randomized
                           FSC. (0-1).
    @return: (resolution,interpolated_band,number_bands)
    @author: Thomas Hrabe
    @todo: Add test!
    """

    fsc_arr: xpt.NDArray = xp.array(fsc)
    number_bands = len(fsc)

    band = number_bands

    randomized_fsc_arr: xpt.NDArray
    if randomized_fsc is None:
        randomized_fsc_arr = xp.ones_like(fsc_arr) * (fsc_arr.min() - 0.1)
    else:
        randomized_fsc_arr = xp.array(randomized_fsc)
    if len(randomized_fsc_arr) < number_bands:
        raise ValueError("Input randomized_fsc should be at least as long as the fsc")

    for i in range(number_bands):
        if fsc_arr[i] < resolution_criterion and fsc_arr[i] > randomized_fsc_arr[i]:
            band = i - 1  # select the band that is still larger than criterion
            break

    if verbose:
        print("Band detected at ", band)

    if band == -1:
        raise RuntimeError("Please check your resolution criterion or you FSC!")

    elif band < number_bands:
        fsc1 = fsc_arr[band]
        fsc2 = fsc_arr[band + 1]

        rfsc1 = randomized_fsc_arr[band]
        rfsc2 = randomized_fsc_arr[band + 1]

        try:
            if fsc2 < rfsc2:
                interpolated_band = (fsc1 - rfsc1) / (rfsc2 - rfsc1 + fsc1 - fsc2)
                pass
            else:
                interpolated_band = (resolution_criterion - fsc1) / (fsc2 - fsc1) + band

        except ZeroDivisionError:
            interpolated_band = band

    else:
        interpolated_band = band

    if verbose:
        print("Band interpolated to ", interpolated_band)

    resolution = (interpolated_band + 1) / float(number_bands)

    if resolution < 0:
        resolution = 1
        interpolated_band = number_bands
        print(
            "Warning: PyTom determined a resolution < 0 for your data. "
            'Please check "mass" in data is positive or negative for all cubes.'
        )
        print(f"Warning: Setting resolution to 1 and {interpolated_band}")
        print("")

    return resolution, interpolated_band, number_bands


def calc_fsc_true(
    fsc_t: xpt.NDArray[float], fsc_n: xpt.NDArray[float], ring_thickness: int = 1
) -> xpt.NDArray:
    """Calculates the true FSC as defined in Henderson
    @param fsc_t: array with FSC values without randomized phases.
    @type fsc_t: ndarray
    @param fsc_n: array with FSC values without randomized phases.
    @type fsc_n: ndarray

    @return: return the fsc_true array
    @rtype: ndarray
    """

    from numpy import zeros_like

    if len(fsc_t) != len(fsc_n):
        raise Exception("FSC arrays are not of equal length.")
    fsc_true = zeros_like(fsc_t)
    steps = 0
    for i in range(len(fsc_t)):
        if abs(fsc_t[i] - fsc_n[i]) < 1e-1:
            fsc_true[i] = fsc_t[i]
        elif steps < ring_thickness:
            fsc_true[i] = fsc_t[i]
            steps += 1
        else:
            fsc_true[i] = (fsc_t[i] - max(0, fsc_n[i])) / (1 - max(0, fsc_n[i]))

    return fsc_true


def generate_random_phases_3d(
    shape: Union[Tuple[int, int], Tuple[int, int, int]], reduced_complex: bool = True
) -> xpt.NDArray[float]:
    """This function returns a set of random phases (between -pi and pi),
       optionally centrosymmetric
    @shape: shape of array
    @type: tuple
    @reduced_complex: is the shape in reduced complex format
    @type: bool, default=True

    @return: returns an array with values between -pi and pi
    @rtype: ndarray
    """

    # Using 'cast' to get arround mypy errors for now
    if len(shape) == 3:
        dx, dy, dz = cast(Tuple[int, int, int], shape)
    else:
        dx, dy = cast(Tuple[int, int], shape)
        dz = max(dx, dy)

    rnda = (xp.random.ranf(shape) * xp.pi * 2) - xp.pi

    cc = dx // 2
    ccc = (dx - 1) // 2

    loc = (dz // 2) * (reduced_complex is False)
    centralslice = rnda[:, :, loc]

    centralslice[cc, cc - ccc : cc] = centralslice[cc, -ccc:][::-1] * -1
    centralslice[cc - ccc : cc, cc - ccc :] = (
        xp.rot90(centralslice[-ccc:, cc - ccc :], 2) * -1
    )

    rnda[:, :, loc] = xp.fft.fftshift(centralslice) if reduced_complex else centralslice

    return rnda


def randomize_phase_beyond_freq(volume: xpt.NDArray, frequency: int) -> xpt.NDArray:
    """This function randomizes the phases beyond a given frequency,
    while preserving the Friedel symmetry.
    @param volume: target volume
    @type volume: L{xpt.NDArray}
    @param frequency: frequency in pixel beyond which phases are randomized.
    @type frequency: int

    @return: 3d volume.
    @rtype: L{xpt.NDArray}
    @author: GvdS"""

    dims = len(volume.shape)
    if dims not in {2, 3}:
        raise Exception("Invalid volume dimensions: either supply 3D or 2D ndarray")

    ft = xp.fft.rfftn(volume)
    phase = xp.angle(ft)

    amplitude = xp.abs(ft)
    rnda = generate_random_phases_3d(
        amplitude.shape, reduced_complex=True
    )  # (ranf((dx,dy,dz//2+1)) * pi * 2) - pi
    # TODO: why does 2D have the extra terms "+ dx %2" and "+ dy %2" and casting to int?
    if dims == 2:
        dx, dy = volume.shape
        x, y = xp.meshgrid(
            xp.arange(-dx // 2, dx // 2 + dx % 2), xp.arange(-dy // 2, dy // 2 + dy % 2)
        )
        rf = xp.sqrt(x**2 + y**2).astype(int)
        r = xp.fft.fftshift(rf)[:, : dy // 2 + 1]
    else:
        dx, dy, dz = volume.shape
        x, y, z = xp.meshgrid(
            xp.arange(-dx // 2, dx // 2),
            xp.arange(-dy // 2, dy // 2),
            xp.arange(-dz // 2, dz // 2),
        )
        rf = xp.sqrt(x**2 + y**2 + z**2)  # .astype(int)
        r = xp.fft.fftshift(rf)[:, :, : dz // 2 + 1]
        # centralslice = fftshift(rnda[:,:,0])
        # cc = dx//2
        # ccc= (dx-1)//2
        # centralslice[cc, cc-ccc:cc] = centralslice[cc,-ccc:][::-1]*-1
        # centralslice[cc-ccc:cc,cc-ccc:] = rot90(centralslice[-ccc:,cc-ccc:],2)*-1
        # rnda[:,:,0] = fftshift(centralslice)

    rnda[r <= frequency] = 0
    phase[r > frequency] = 0
    phase += rnda
    image = xp.fft.irfftn((amplitude * xp.exp(1j * phase)), s=volume.shape)

    if xp.abs(image.imag).sum() > 1e-8:
        raise Exception(
            "Imaginary part is non-zero. Failed to centro-summetrize the phases."
        )

    return image.real


def soc(
    volume: xpt.NDArray[float],
    reference: xpt.NDArray[float],
    mask: Optional[xpt.NDArray[float]] = None,
    std_v: Optional[float] = None,
) -> xpt.NDArray[float]:
    """
    soc : Second Order Correlation. Correlation of correlation peaks.
    @param volume: The volume
    @type volume:  L{xpt.NDArray}
    @param reference: The reference / template
    @type reference:  L{xpt.NDArray}
    @author: Thomas Hrabe
    """

    reference_peak = flcf(reference, reference, mask)
    peaks = flcf(volume, reference, mask)

    return flcf(peaks, reference_peak, mask)
