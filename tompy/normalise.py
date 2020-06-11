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
        volume = xp.fft.ifftshift(ifftn(volume))

    if not copyFlag:
        volume -= volume.mean()
        volumeStd = volume.std()

        if volumeStd < epsilon:
            raise ValueError(
                'pytom_normalise.mean0std1 : The standard deviation is too low for division! ' + str(volumeStd))

        volume /= volumeStd
    else:
        volumeCopy = volume.copy()
        volumeStd = volume.std()

        if volumeStd < epsilon:
            raise ValueError(
                'pytom_normalise.mean0std1 : The standard deviation is too low for division! ' + str(volumeStd))

        # volumeCopy.shiftscale(-1.*volumeMean,1)
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
    from pytom.tompy.correlation import meanUnderMask, stdUnderMask
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