"""
procedures to normalize according to mean and standard deviation
"""
def mean0std1(volume,copyFlag=False):
    """
    mean0std1: normalises input volume to mean 0 and std 1. Procedure is performed inplace if copyFlag is unspecified!!!
    @param volume: Data containing either an image or a volume
    @param copyFlag: If True a copy of volume will be returned. False unless specified otherwise. 
    @return: If copyFlag == True, then return a normalised copy.
    @author: Thomas Hrabe
    """
    import pytom_volume
    from math import sqrt
    from pytom.tools.maths import epsilon
    
    if volume.__class__ == pytom_volume.vol_comp:
        from pytom.basic.fourier import ifft,iftshift
        volume = iftshift(ifft(volume))

    if not copyFlag:
        volumeMean = pytom_volume.mean(volume)
        volume.shiftscale(-volumeMean,1)
        volumeStd = sqrt(pytom_volume.variance(volume,False))
    
        if volumeStd < epsilon:
            raise ValueError('pytom_normalise.mean0std1 : The standard deviation is too low for division! ' + str(volumeStd))
    
        volume.shiftscale(0,1./float(volumeStd))
    else:
        volumeCopy = pytom_volume.vol(volume.sizeX(),volume.sizeY(),volume.sizeZ())
        volumeCopy.copyVolume(volume)
        volumeMean = pytom_volume.mean(volumeCopy)
        volumeCopy.shiftscale(-1.*volumeMean,1)
        volumeStd = sqrt(pytom_volume.variance(volumeCopy,False))
    
        if volumeStd < epsilon:
            raise ValueError('pytom_normalise.mean0std1 : The standard deviation is too low for division! ' + str(volumeStd))
    
        #volumeCopy.shiftscale(-1.*volumeMean,1)
        volumeCopy.shiftscale(0,1/float(volumeStd))
         
        return volumeCopy


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
    from pytom.basic.correlation import meanValueUnderMask, stdValueUnderMask
    #from math import sqrt
    if not p:
        from pytom_volume import sum
        p = sum(mask)
    #meanT = sum(volume) / p
    ## subtract mean and mask
    #res = mask * (volume - meanT)
    #stdT = sum(res*res) / p
    #stdT = sqrt(stdT)
    #res = res / stdT


    meanT = meanValueUnderMask(volume, mask, p)
    
    stdT = stdValueUnderMask(volume, mask, meanT, p)
    res = (volume - meanT)/stdT
    return (res,p)


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
    from pytom_volume import sum as sumvol
    #npix = volume.sizeX() * volume.sizeY() * volume.sizeZ()
    normvol = volume*mask
    normvol = sumvol(normvol) / sumvol(mask)
    return normvol

