from pytom.gpu.initialize import xp, device

def mean0std1(volume, copyFlag=False):
    """
    mean0std1: normalises input volume to mean 0 and std 1. Procedure is performed inplace if copyFlag is unspecified!!!
    @param volume: Data containing either an image or a volume
    @param copyFlag: If True a copy of volume will be returned. False unless specified otherwise.
    @return: If copyFlag == True, then return a normalised copy.
    @author: Thomas Hrabe
    """
    from pytom.tools.maths import epsilon

    if volume.dtype == xp.complex64:
        volume = xp.fft.ifftshift(xp.fft.ifftn(volume))

    if not copyFlag:
        volume -= volume.mean()
        volumeStd = volume.std()

        # if volumeStd < epsilon:
        #     raise ValueError(
        #         'pytom_normalise.mean0std1 : The standard deviation is too low for division! ' + str(volumeStd))

        volume /= volumeStd if volumeStd > epsilon else 1
    else:
        volumeCopy = volume.copy()
        volumeStd = volume.std()

        # if volumeStd < epsilon:
        #     raise ValueError(
        #         'pytom_normalise.mean0std1 : The standard deviation is too low for division! ' + str(volumeStd))

        # volumeCopy.shiftscale(-1.*volumeMean,1)

        volumeStd = 1 if volumeStd < epsilon else volumeStd

        return (volumeCopy - volume.mean()) / volumeStd

def normaliseUnderMask(volume, mask, p=None):
    """
    normalize volume within a mask - take care: only normalization, but NOT multiplication with mask!

    @param volume: volume for normalization
    @type volume: pytom volume
    @param mask: mask
    @type mask: pytom volume
    @param p: sum of gray values in mask (if pre-computed)
    @type p: C{int} or C{float}
    @return: volume normalized to mean=0 and std=1 within mask, p
    @rtype: C{list}
    @author: FF
    """
    from pytom.agnostic.correlation import meanUnderMask, stdUnderMask
    # from math import sqrt
    if not p:
        p = xp.sum(mask)
    # meanT = sum(volume) / p
    ## subtract mean and mask
    # res = mask * (volume - meanT)
    # stdT = sum(res*res) / p
    # stdT = sqrt(stdT)
    # res = res / stdT

    meanT = meanUnderMask(volume, mask, p)

    stdT = stdUnderMask(volume, mask, meanT, p)
    res = (volume - meanT) / stdT
    return (res, p)

def subtractMeanUnderMask(volume, mask):
    """
    subtract mean from volume/image under mask

    @param volume: volume/image
    @type volume: L{pytom_volume.vol}
    @param mask: mask
    @type mask: L{pytom_volume.vol}
    @return: volume/image
    @rtype: L{pytom_volume.vol}
    """
    # npix = volume.sizeX() * volume.sizeY() * volume.sizeZ()
    normvol = volume * mask
    normvol = xp.sum(normvol) / xp.sum(mask)
    return normvol

def meanUnderMask(volume, mask=None, p=1, gpu=False):
    """
    meanValueUnderMask: Determines the mean value under a mask
    @param volume: The volume
    @type volume:  L{pytom_volume.vol}
    @param mask:  The mask
    @type mask:  L{pytom_volume.vol}
    @param p: precomputed number of voxels in mask
    @type p: float
    @return: A value (scalar)
    @rtype: single
    @change: support None as mask, FF 08.07.2014
    """

    return (volume*mask).sum() / p

def stdUnderMask(volume, mask, meanValue, p=1, gpu=False):
    """
    stdValueUnderMask: Determines the std value under a mask

    @param volume: input volume
    @type volume:  L{pytom_volume.vol}
    @param mask: mask
    @type mask:  L{pytom_volume.vol}
    @param p: non zero value numbers in the mask
    @type p: L{float} or L{int}
    @return: A value
    @rtype: L{float}
    @change: support None as mask, FF 08.07.2014
    """
    return (meanUnderMask(volume**2, mask, p) - meanValue**2)**0.5

def meanVolUnderMask(volume, mask):
    """
    meanUnderMask: calculate the mean volume under the given mask (Both should have the same size)
    @param volume: input volume
    @type volume:  L{numpy.ndarray} L{cupy.ndarray}
    @param mask: mask
    @type mask:  L{numpy.ndarray} L{cupy.ndarray}
    @param p: non zero value numbers in the mask
    @type p: L{int} or L{float}
    @param gpu: Boolean that indicates if the input volumes stored on a gpu.
    @type gpu: L{bool}
    @return: the calculated mean volume under mask
    @rtype:  L{numpy.ndarray} L{cupy.ndarray}
    @author: Gijs van der Schot
    """
    res = xp.fft.fftshift(xp.fft.ifftn(xp.fft.fftn(volume) * xp.conj(xp.fft.fftn(mask)))) / mask.sum()
    return res.real

def stdVolUnderMask(volume, mask, meanV):
    """
    stdUnderMask: calculate the std volume under the given mask
    @param volume: input volume
    @type volume:  L{numpy.ndarray} L{cupy.ndarray}
    @param mask: mask
    @type mask:  L{numpy.ndarray} L{cupy.ndarray}
    @param p: non zero value numbers in the mask
    @type p: L{int}
    @param meanV: mean volume under mask, which should already been caculated
    @type meanV:  L{numpy.ndarray} L{cupy.ndarray}
    @return: the calculated std volume under mask
    @rtype:  L{numpy.ndarray} L{cupy.ndarray}
    @author: GvdS
    """

    meanV2 = meanV * meanV
    vol2 = volume * volume
    var = meanVolUnderMask(vol2, mask) - meanV2
    var[var<1E-09] = 1

    return var**0.5

def meanVolUnderMaskPlanned(volume, mask, ifftnP, fftnP, plan):
    """
    meanUnderMask: calculate the mean volume under the given mask (Both should have the same size)
    @param volume: input volume
    @type volume:  L{numpy.ndarray} L{cupy.ndarray}
    @param mask: mask
    @type mask:  L{numpy.ndarray} L{cupy.ndarray}
    @param p: non zero value numbers in the mask
    @type p: L{int} or L{float}
    @param gpu: Boolean that indicates if the input volumes stored on a gpu.
    @type gpu: L{bool}
    @return: the calculated mean volume under mask
    @rtype:  L{numpy.ndarray} L{cupy.ndarray}
    @author: Gijs van der Schot
    """
    volume_fft = fftnP(volume.astype(xp.complex64), plan=plan)
    mask_fft = fftnP(mask.astype(xp.complex64), plan=plan)
    res = xp.fft.fftshift(ifftnP(volume_fft * xp.conj(mask_fft), plan=plan)) / mask.sum()
    return res.real

def stdVolUnderMaskPlanned(volume, mask, meanV, ifftnP, fftnP, plan):
    """
    stdUnderMask: calculate the std volume under the given mask
    @param volume: input volume
    @type volume:  L{numpy.ndarray} L{cupy.ndarray}
    @param mask: mask
    @type mask:  L{numpy.ndarray} L{cupy.ndarray}
    @param p: non zero value numbers in the mask
    @type p: L{int}
    @param meanV: mean volume under mask, which should already been caculated
    @type meanV:  L{numpy.ndarray} L{cupy.ndarray}
    @return: the calculated std volume under mask
    @rtype:  L{numpy.ndarray} L{cupy.ndarray}
    @author: GvdS
    """

    meanV2 = meanV * meanV
    vol2 = volume * volume
    var = meanVolUnderMaskPlanned(vol2, mask, ifftnP, fftnP, plan) - meanV2
    var[var<1E-09] = 1

    return var**0.5