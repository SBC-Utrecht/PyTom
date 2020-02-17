
from pytom.tompy.tools import paste_in_center, create_sphere
from pytom.gpu.initialize import xp

def meanUnderMask(volume, mask=None, p=None, gpu=False):
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

    return (volume*mask).sum() / mask.sum()

def stdUnderMask(volume, mask, meanValue, p=None, gpu=False):
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
    return (meanUnderMask(volume**2, mask, mask.sum(), gpu=gpu) - meanValue**2)**0.5

def meanVolUnderMask(volume, mask, gpu=False, xp=None):
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
    res = xp.fft.fftshift(xp.fft.irfftn(xp.fft.rfftn(volume) * xp.conj(xp.fft.rfftn(mask)))) / mask.sum()
    return res.real

def stdVolUnderMask(volume, mask, meanV, gpu=False):
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
    var = meanUnderMask(vol2, mask, gpu=gpu) - meanV2
    var[var<1E-09] = 1

    return var**0.5

def FLCF(volume, template, mask=None, stdV=None, gpu=False, mempool=None, pinned_mempool=None):
    '''Fast local correlation function

    @param volume: target volume
    @param template: template to be searched. It can have smaller size then target volume.
    @param mask: template mask. If not given, a default sphere mask will be used.
    @param stdV: standard deviation of the target volume under mask, which do not need to be calculated again when the mask is identical.

    @return: the local correlation function
    '''

    meanT = meanUnderMask(template, mask, gpu=gpu)
    temp = ((template - meanT ) / stdUnderMask(template, mask, meanT,gpu=gpu)) * mask

    res =  xp.fft.fftshift(xp.fft.ifftn(xp.conj(xp.fft.fftn(temp)) * xp.fft.fftn(volume))).real / stdV / mask.sum()
    meanT = None
    temp=None
    mempool.free_all_blocks()
    pinned_mempool.free_all_blocks()

def xcc(volume, template, mask=None, volumeIsNormalized=False, gpu=False):
    """
    xcc: Calculates the cross correlation coefficient in real space
    @param volume: A volume
    @type volume:  L{pytom_volume.vol}
    @param template: A template that is searched in volume. Must be of same size as volume.
    @type template:  L{pytom_volume.vol}
    @param mask: mask to constrain correlation
    @type mask: L{pytom_volume.vol}
    @param volumeIsNormalized: only used for compatibility with nxcc - not used
    @type volumeIsNormalized: L{bool}
    @return: A unscaled value
    @raise exception: Raises a runtime error if volume and template have a different size.
    @author: Thomas Hrabe
    """


    from pytom.tompy.macros import volumesSameSize

    if not volumesSameSize(volume, template):
        raise RuntimeError('Volume and template must have same size!')

    if mask:  # mask is given
        result = mask * volume * template
    else:
        result = volume * template

    cc = result.sum()

    cc = cc / float(volume.size)

    return cc

def nxcc(volume, template, mask=None, volumeIsNormalized=False):
    """
    nxcc: Calculates the normalized cross correlation coefficient in real space
    @param volume: A volume
    @type volume:  L{pytom_volume.vol}
    @param template: A template that is searched in volume. Must be of same size as volume.
    @type template:  L{pytom_volume.vol}
    @param mask: mask to constrain correlation
    @type mask: L{pytom_volume.vol}
    @param volumeIsNormalized: speed up if volume is already normalized
    @type volumeIsNormalized: L{bool}
    @return: A value between -1 and 1
    @raise exception: Raises a runtime error if volume and template have a different size.
    @author: Thomas Hrabe
    @change: flag for pre-normalized volume, FF
    """

    from pytom.tompy.macros import volumesSameSize

    if not volumesSameSize(volume, template):
        raise RuntimeError('Volume and template must have same size!')

    if not mask:
        from pytom.tompy.normalise import mean0std1
        if not volumeIsNormalized:
            v = mean0std1(volume, True)
        t = mean0std1(template, True)
        p = volume.size
        result = v * t
    else:
        from pytom.tompy.normalise import normaliseUnderMask
        if not volumeIsNormalized:
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

    return ncc


def mean0std1(volume, copyFlag=False):
    """
    mean0std1: normalises input volume to mean 0 and std 1. Procedure is performed inplace if copyFlag is unspecified!!!
    @param volume: Data containing either an image or a volume
    @param copyFlag: If True a copy of volume will be returned. False unless specified otherwise.
    @return: If copyFlag == True, then return a normalised copy.
    @author: Thomas Hrabe
    """
    from pytom.tools.maths import epsilon


    volume = xp.fft.fftshift(xp.fft.ifftn(volume))

    volumeMean = volume.mean()
    volume -= volumeMean
    volumeStd = volume.std()

    volume /= volumeStd if volumeStd > epsilon else 1

    return volume

def xcf(volume, template, mask=None, stdV=None, ):
    """
    XCF: returns the non-normalised cross correlation function. The xcf
    result is scaled only by the square of the number of elements.

    @param volume : The search volume
    @type volume: L{pytom_volume.vol}
    @param template : The template searched (this one will be used for conjugate complex multiplication)
    @type template: L{pytom_volume.vol}
    @param mask: Will be unused, only for compatibility reasons with FLCF
    @param stdV: Will be unused, only for compatibility reasons with FLCF
    @return: XCF volume
    @rtype: L{pytom_volume.vol}
    @author: Thomas Hrabe
    """

    # if mask:
    #    volume = volume * mask
    #	 template = template * mask

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

    #xp.fft.iftshift(result)

    return result

def xcf_mult(volume, template, mask, stdV=None, ):
    """
    XCF: returns the non-normalised cross correlation function. The xcf
    result is scaled only by the square of the number of elements.

    @param volume : The search volume
    @type volume: L{pytom_volume.vol}
    @param template : The template searched (this one will be used for conjugate complex multiplication)
    @type template: L{pytom_volume.vol}
    @param mask: Will be unused, only for compatibility reasons with FLCF
    @param stdV: Will be unused, only for compatibility reasons with FLCF
    @return: XCF volume
    @rtype: L{pytom_volume.vol}
    @author: Thomas Hrabe
    """

    # if mask:
    #    volume = volume * mask
    #	 template = template * mask

    # determine fourier transforms of volumes
    fvolume = xp.fft.fftn(volume)
    ftemplate = xp.fft.fftn(template)

    # perform element wise - conjugate multiplication
    ftemplate = xp.conj(ftemplate)
    fresult = fvolume * ftemplate

    # transform back to real space
    result = abs(xp.fft.ifftn(fresult))

    #xp.fft.iftshift(result)

    return result

def nXcf(volume, template, mask=None, stdV=None, gpu=False):
    """
    nXCF: returns the normalised cross correlation function. Autocorrelation
    of two equal objects would yield a max nxcf peak of 1.

    @param volume: The search volume
    @param template: The template searched (this one will be used for conjugate complex multiplication)
    @type template: L{pytom_volume.vol}
    @param mask: template mask. If not given, a default sphere mask will be generated which has the same size with the given template.
    @type mask: L{pytom_volume.vol}
    @param stdV: Will be unused, only for compatibility reasons with FLCF
    @return: the calculated nXcf volume
    @rtype: L{pytom_volume.vol}
    @author: Thomas Hrabe
    @change: masking of template implemented
    """

    if not (mask is None):
        result = xcf_mult(normaliseUnderMask(volume=volume, mask=mask, p=None),
                     normaliseUnderMask(volume=template, mask=mask, p=None), mask)
    else:
        result = xcf_mult(mean0std1(volume), mean0std1(template), mask)


    n = result.size

    result /= float(n * n)

    return result


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

    #from math import sqrt
    if not p:
        p = mask.sum()
    #meanT = sum(volume) / p
    ## subtract mean and mask
    #res = mask * (volume - meanT)
    #stdT = sum(res*res) / p
    #stdT = sqrt(stdT)
    #res = res / stdT

    meanT = meanUnderMask(volume, mask, p)

    stdT = stdUnderMask(volume, mask, meanT, p)
    res = (volume - meanT)/stdT
    return res

def bandCC(volume,reference,band,verbose = False, gpu=False):
    """
    bandCC: Determines the normalised correlation coefficient within a band
    @param volume: The volume
    @type volume: L{pytom_volume.vol}
    @param reference: The reference
    @type reference: L{pytom_volume.vol}
    @param band: [a,b] - specify the lower and upper end of band.
    @return: First parameter - The correlation of the two volumes in the specified band. 
             Second parameter - The bandpass filter used.
    @rtype: List - [float,L{pytom_freqweight.weight}]
    @author: Thomas Hrabe    
    """



    from pytom.tompy.filter import bandpass
    from pytom.tompy.correlation import xcf
    #from pytom.tompy.filter import vol_comp
    
    if verbose:
        print('lowest freq : ', band[0],' highest freq' , band[1])
        
    vf = bandpass(volume,band[0],band[1])
    rf = bandpass(reference,band[0],band[1])#,vf[1])
    
    #ccVolume = vol_comp(rf[0].shape[0],rf[0].shape[1],rf[0].shape[2])
    #ccVolume.copyVolume(rf[0])

    ccVolume = rf.astype('complex128')

    ccVolume = ccVolume * xp.conj(vf[0])
    #pytom_volume.conj_mult(ccVolume,vf[0])
    
    cc = ccVolume.sum()

    cc = cc.real
    
    v = vf[0]
    r = rf[0]
    
    absV = xp.abs(v)
    absR = xp.abs(r)

    sumV = xp.sum(absV**2)
    sumR = xp.sum(absR**2)
    
    sumV = xp.abs(sumV)
    sumR = xp.abs(sumR)
    
    if sumV == 0:
        sumV =1
        
    if sumR == 0:
        sumR =1
        
    cc = cc / (xp.sqrt(sumV*sumR))
    
    #numerical errors will be punished with nan
    if abs(cc) > 1.1 :
        cc = float('nan')
    
    return [cc,vf[1]]
    
def weightedXCC(volume,reference,numberOfBands,wedgeAngle=-1, gpu=False):
    """
    weightedXCC: Determines the band weighted correlation coefficient for a volume and reference. Notation according Steward/Grigorieff paper
    @param volume: A volume
    @type volume: L{pytom_volume.vol}
    @param reference: A reference of same size as volume
    @type reference: L{pytom_volume.vol}
    @param numberOfBands: Number of bands
    @param wedgeAngle: A optional wedge angle
    @return: The weighted correlation coefficient
    @rtype: float
    @author: Thomas Hrabe
    """

    if gpu:
        import cupy as xp
    else:
        import numpy as xp

    from pytom.tompy.transform import fft, fourier_reduced2full
    from pytom.basic.structures import WedgeInfo

    result = 0
    numberVoxels = 0

    #volume.write('vol.em');
    #reference.write('ref.em');


    wedge = WedgeInfo(wedgeAngle)
    wedgeVolume = wedge.returnWedgeVolume(volume.shape[0],volume.shape[1],volume.shape[2])

    increment = int(volume.shape[0]/2 * 1/numberOfBands)
    band = [0,100]
    for i in range(0,volume.shape[0]/2, increment):

        band[0] = i
        band[1] = i + increment

        r = bandCC(volume, reference, band, gpu=gpu)
        cc = r[0]

        #print cc;
        filter = r[1]

        #get bandVolume
        bandVolume = filter.getWeightVolume(True)

        filterVolumeReduced = bandVolume * wedgeVolume
        filterVolume = fourier_reduced2full(filterVolumeReduced)

        #determine number of voxels != 0
        N = (filterVolume != 0).sum()

        w = xp.sqrt(1/float(N))

        #add to number of total voxels
        numberVoxels=numberVoxels + N
        #print 'w',w;
        #print 'cc',cc;
        #print 'N',N;

        cc2 = cc*cc
        #print 'cc2',cc2;
        if cc <= 0.0:
            cc = cc2
        else:
            cc = cc2/(cc+w)

        #print 'cc',cc;
        cc = cc *cc *cc; #no abs
        #print 'cc',cc;

        #add up result
        result = result + cc*N

    return result*(1/float(numberVoxels))

def FSCSum(volume,reference,numberOfBands,wedgeAngle=-1, gpu=False):
    """
    FSCSum: Determines the sum of the Fourier Shell Correlation coefficient for a volume and reference. 
    @param volume: A volume
    @type volume: L{pytom_volume.vol}
    @param reference: A reference of same size as volume
    @type reference: L{pytom_volume.vol}
    @param numberOfBands: Number of bands
    @param wedgeAngle: A optional wedge angle
    @return: The sum FSC coefficient
    @rtype: float  
    @author: Thomas Hrabe   
    """    

    if gpu:
        import cupy as xp
    else:
        import numpy as xp

    from pytom.tompy.correlation import bandCC


    result = 0
    numberVoxels = 0
    
    #volume.write('vol.em');
    #reference.write('ref.em');
    fvolume = xp.fft.fftn(volume)
    freference = xp.fft.fftn(reference)
    numelem = volume.size
    
    fvolume = shiftscale(fvolume, 0, 1/float(numelem))
    freference = shiftscale(fvolume, 0, 1/float(numelem))

    #print '-----'
    for i in range(numberOfBands):
        #process bandCorrelation
        band = []
        band[0] = i*volume.shape[0]/numberOfBands
        band[1] = (i+1)*volume.shape[0]/numberOfBands
        
        r = bandCC(fvolume, freference, band, gpu=gpu)
        cc = r[0]
        #print cc
        result = result + cc
    #print '-----'
    
    return result*(1/float(numberOfBands))

def bandCF(volume,reference,band=[0,100]):
    """
    bandCF:
    @param volume: The volume
    @param reference: The reference
    @param band: [a,b] - specify the lower and upper end of band. [0,1] if not set.
    @return: First parameter - The correlation of the two volumes in the specified ring. 
             Second parameter - The bandpass filter used.
    @rtype: List - [L{pytom_volume.vol},L{pytom_freqweight.weight}]
    @author: Thomas Hrabe   
    @todo: does not work yet -> test is disabled
    """

    if gpu:
        import cupy as xp
    else:
        import numpy as xp

    import pytom_volume
    from math import sqrt
    from pytom.basic import fourier
    from pytom.basic.filter import bandpassFilter
    from pytom.basic.correlation import nXcf
    
    vf = bandpassFilter(volume,band[0],band[1],fourierOnly=True)
    rf = bandpassFilter(reference,band[0],band[1],vf[1],fourierOnly=True)
    
    v = pytom_volume.reducedToFull(vf[0])
    r = pytom_volume.reducedToFull(rf[0])
    
    absV = pytom_volume.abs(v)
    absR = pytom_volume.abs(r)
    
    pytom_volume.power(absV,2)
    pytom_volume.power(absR,2)
    
    sumV = abs(pytom_volume.sum(absV))
    sumR = abs(pytom_volume.sum(absR))
    
    if sumV == 0:
        sumV = 1
        
    if sumR == 0:
        sumR = 1
    
    pytom_volume.conjugate(rf[0])
    
    fresult = vf[0] * rf[0]

    #transform back to real space
    result = fourier.ifft(fresult)
    
    fourier.iftshift(result)
    
    result.shiftscale(0,1/float(sqrt(sumV*sumR)))
    
    return [result,vf[1]]

def weightedXCF(volume,reference,numberOfBands,wedgeAngle=-1, gpu=False):
    """
    weightedXCF: Determines the weighted correlation function for volume and reference
    @param volume: A volume 
    @param reference: A reference 
    @param numberOfBands:Number of bands
    @param wedgeAngle: A optional wedge angle
    @return: The weighted correlation function 
    @rtype: L{pytom_volume.vol} 
    @author: Thomas Hrabe 
    @todo: does not work yet -> test is disabled
    """

    if gpu:
        import cupy as xp
    else:
        import numpy as xp

    from pytom.tompy.correlation import bandCF
    from pytom.tompy.transforms import fourier_reduced2full
    from math import sqrt
    import pytom_freqweight
    
    result = xp.zeros_like(volume)

    q = 0
    
    if wedgeAngle >=0:
        wedgeFilter = pytom_freqweight.weight(wedgeAngle, 0, volume.sizeX(),volume.sizeY(),volume.sizeZ())
        wedgeVolume = wedgeFilter.getWeightVolume(True)
    else:
        wedgeVolume = xp.ones_like(volume)
        
    w = xp.sqrt(1/float(volume.size))
    
    numberVoxels = 0
    
    for i in range(numberOfBands):
        """
        notation according Steward/Grigorieff paper
        """
        band = [0,0]
        band[0] = i*volume.shape[0]/numberOfBands
        band[1] = (i+1)*volume.shape[0]/numberOfBands
        
        r = bandCF(volume,reference,band, gpu=gpu)
        
        cc = r[0]
                
        filter = r[1]
        #get bandVolume
        bandVolume = filter.getWeightVolume(True)
            
        filterVolumeReduced = bandVolume * wedgeVolume
        filterVolume = fourier_reduced2full(filterVolumeReduced)
        #determine number of voxels != 0    
        N = filterVolume[abs(filterVolume) < 1].sum()
            
        #add to number of total voxels
        numberVoxels = numberVoxels + N
                 
        cc2 = r[0].copy()
                
        cc2 = cc2**2
        
        cc = shiftscale(cc, w,1)
        ccdiv = cc2/(cc)
                
        ccdiv = ccdiv**3
        
        #abs(ccdiv); as suggested by grigorief
        ccdiv = shiftscale(ccdiv, 0, N)
        
        result = result + ccdiv
    
    result = shiftscale(result, 0,1/float(numberVoxels))
    
    return result

def FSC(volume1, volume2, numberBands, mask=None, verbose=False, filename=None, gpu=False):
    """
    FSC - Calculates the Fourier Shell Correlation for two volumes
    @param volume1: volume one
    @type volume1: L{pytom_volume.vol}
    @param volume2: volume two
    @type volume2: L{pytom_volume.vol}
    @param numberBands: number of shells for FSC
    @type numberBands: int
    @param filename: write FSC to ascii file if specified
    @type filename: string

    @return: Returns a list of cc values 
    @author: Thomas Hrabe  
    @rtype: list[floats]
    """

    if gpu:
        import cupy as xp
    else:
        import numpy as xp

    from pytom.tompy.correlation import bandCC
    from pytom.basic.structures import Mask
    from pytom.tompy.io import read


    if not volume1.shape == volume2.shape:
        raise RuntimeError('Volumes must have the same size!')
    print(mask.__class__, xp.array([]).__class__)
    if not mask is None:
        if mask.__class__ == xp.array([]).__class__:
            volume1 = volume1 * mask
            volume2 = volume2 * mask
          
        elif mask.__class__ == Mask:
            mask = mask.getVolume()
            volume1 = volume1 * mask
            volume2 = volume2 * mask

        elif mask.__class__ == str:
            mask = read(mask)
            volume1 = volume1 * mask
            volume2 = volume2 * mask 

        else:
            raise RuntimeError('FSC: Mask must be a volume OR a Mask object OR a string path to a mask!')  
        
    fscResult = []
    band = [-1,-1]
    
    increment = int(volume1.shape[0]/2 * 1/numberBands)
    
    for i in range(0,volume1.shape[0]//2, increment):
        
        band[0] = i
        band[1] = i + increment
        
        if verbose:
            print('Band : ' ,band)
            
        res = bandCC(volume1,volume2,band,verbose, gpu=gpu)
        
        if i == 0 and increment == 1:
            #force a 1 for correlation of the zero frequency 
            res[0] = 1
  
        if verbose:
            print('Correlation ' ,res[0])

        fscResult.append(res[0])

    if filename:
        f = open(filename,'w')
        for item in fscResult:
            f.write("%s\n" % item)
        f.close()

    return fscResult

def determineResolution(fsc,resolutionCriterion, verbose=False, randomizedFSC=None, gpu=False):
    """
    determineResolution: Determines frequency and band where correlation drops below the resolutionCriterion. Uses linear interpolation between two positions
    @param fsc: The fsc list determined by L{pytom.basic.correlation.FSC}
    @param resolutionCriterion: A value between 0 and 1
    @return: [resolution,interpolatedBand,numberBands] 
    @author: Thomas Hrabe 
    @todo: Add test! 
    """

    if gpu:
        import cupy as xp
    else:
        import numpy as xp

    numberBands = len(fsc)
    
    band = numberBands

    if randomizedFSC is None:
        randomizedFSC = xp.ones_like(fsc)*(fsc.min()-0.1)

    for i in range(numberBands):
        if fsc[i] < resolutionCriterion and fsc[i] > randomizedFSC[i]:
            band = i-1  #select the band that is still larger than criterion
            break
    
    if verbose:
        print('Band detected at ', band)
    
    if band == -1:
        raise RuntimeError("Please check your resolution criterion or you FSC!")
            
    elif band < numberBands:
        fsc1 = fsc[band]
        fsc2 = fsc[band+1]

        rfsc1 = randomizedFSC[band]
        rfsc2 = randomizedFSC[band+1]

        try:
            if fsc2 < rfsc2:
                interpolatedBand = (fsc1-rfsc1)/(rfsc2-rfsc1+fsc1-fsc2)
                pass
            else:
                interpolatedBand = (resolutionCriterion-fsc1)/(fsc2-fsc1)+band
        
        except ZeroDivisionError:
            interpolatedBand = band
        
    else:
        interpolatedBand = band
        
    if verbose:
        print('Band interpolated to ', interpolatedBand)
        
    resolution = (interpolatedBand+1) / float(numberBands)
    
    if resolution < 0 :
        resolution = 1
        interpolatedBand = numberBands
        print('Warning: PyTom determined a resolution < 0 for your data. Please check "mass" in data is positive or negative for all cubes.')
        print('Warning: Setting resolution to 1 and ',interpolatedBand)
        print('')
        
    return [resolution, interpolatedBand, numberBands]

def soc(volume,reference,mask=None, stdV=None, gpu=False):
    """
    soc : Second Order Correlation. Correlation of correlation peaks.
    @param volume: The volume
    @type volume:  L{pytom_volume.vol}
    @param reference: The reference / template
    @type reference:  L{pytom_volume.vol}
    @author: Thomas Hrabe   
    """

    if gpu:
        import cupy as xp
    else:
        import numpy as xp

    referencePeak = FLCF(reference, reference, mask, gpu=gpu)
    peaks = FLCF(volume, reference, mask, gpu=gpu)
    
    return FLCF(peaks, referencePeak, mask, gpu=gpu)

def isBorderVoxel(volume, coordinates):

    shape = volume.shape
    for n,pos in enumerate(coordinates):
        if pos < 1 or pos >= shape[n]-1:
            return True
    return False

def subPixelPeakParabolic(scoreVolume, coordinates, verbose=False, gpu=False):
    """
    quadratic interpolation of three adjacent samples
    @param scoreVolume: The score volume
    @param coordinates: [x,y,z] coordinates where the sub pixel peak will be determined
    @param verbose: be talkative
    @type verbose: bool
    @return: Returns [peakValue,peakCoordinates] with sub pixel accuracy
    """
    if gpu:
        import cupy as xp
    else:
        import numpy as xp

    if isBorderVoxel(scoreVolume, coordinates):
        if verbose:
            print("subPixelPeakParabolic: peak near borders - no interpolation done")
        return [scoreVolume[coordinates], coordinates]

    peakCoordinates = coordinates
    l = len(coordinates)
    (x, p1, a1) = qint(ym1=scoreVolume[coordinates-[1,0,0][:l]],
                     y0=scoreVolume[coordinates],
                     yp1=scoreVolume[coordinates+[1,0,0][:l]])
    (y, p2, a2) = qint(ym1=scoreVolume[coordinates - [0,1,0][:l]],
                     y0=scoreVolume[coordinates],
                     yp1=scoreVolume[coordinates+ [0,1,0][:l] ])
    if l > 2:
        (z, p3, a3) = qint(ym1=scoreVolume.getV(coordinates-[0,0,1]),
                         y0=scoreVolume.getV(coordinates),
                         yp1=scoreVolume.getV(coordinates+[0,0,1]))
        peakCoordinates[0] += x
        peakCoordinates[1] += y
        peakCoordinates[2] += z
        peakValue = p1 + a2*y**2 + a3*z**2
    else:
        peakCoordinates[0] += x
        peakCoordinates[1] += y
        peakValue = p1 + a2*y**2

    return [peakValue,peakCoordinates]

def qint(ym1, y0, yp1):
    """
    quadratic interpolation of three adjacent samples
     
    [p,y,a] = qint(ym1,y0,yp1)
     
    returns the extremum location p, height y, and half-curvature a
    of a parabolic fit through three points.
    Parabola is given by y(x) = a*(x-p)^2+b,
    where y(-1)=ym1, y(0)=y0, y(1)=yp1.
    @param ym1: y(-1)
    @type ym1: float
    @param y0: y(0)
    @type y0: float
    @param yp1: y(+1)
    @type yp1: float
    @return: peak-ordinate, peak-value
    """
    p = (yp1 - ym1)/(2*(2*y0 - yp1 - ym1))
    y = y0 - 0.25*(ym1-yp1)*p
    a = 0.5*(ym1 - 2*y0 + yp1)
    #b = y0 - a*(p**2)
    return p, y, a
            
def subPixelPeak(scoreVolume, coordinates, cubeLength=8, interpolation='Spline', verbose=False):
    """
    subPixelPeak: Will determine the sub pixel area of peak. Utilizes spline, fourier or parabolic interpolation.

    @param verbose: be talkative
    @type verbose: L{str}
    @param scoreVolume: The score volume
    @param coordinates: [x,y,z] coordinates where the sub pixel peak will be determined
    @param cubeLength: length of cube - only used for Spline and Fourier interpolation
    @type cubeLength: int (even)
    @param interpolation: interpolation type: 'Spline', 'Quadratic', or 'Fourier'
    @type interpolation: str
    @return: Returns [peakValue,peakCoordinates] with sub pixel accuracy

    last change: 02/07/2013 FF: 2D functionality added
    """
    assert type(interpolation) == str, 'subPixelPeak: interpolation must be str'
    if (interpolation.lower() == 'quadratic') or (interpolation.lower() == 'parabolic'):
        (peakValue,peakCoordinates) = subPixelPeakParabolic(scoreVolume=scoreVolume, coordinates=coordinates, verbose=verbose)
        return [peakValue,peakCoordinates]

    if gpu:
        import cupy as xp
    else:
        import numpy as xp

    from pytom_volume import vol,subvolume,rescaleSpline,peak
    from pytom.basic.transformations import resize
  
    #extend function for 2D
    twoD = (scoreVolume.shape) == 2

    cubeStart = cubeLength//2
    sizeX = scoreVolume.sizeX()
    sizeY = scoreVolume.sizeY()
    sizeZ = scoreVolume.sizeZ()
    
    if twoD:
        if (coordinates[0]-cubeStart < 1 or coordinates[1]-cubeStart < 1) or\
            (coordinates[0]-cubeStart + cubeLength >= sizeX or coordinates[1]-cubeStart + cubeLength >= sizeY):
            if verbose:
                print ("SubPixelPeak: position too close to border for sub-pixel")
            return [scoreVolume(coordinates[0],coordinates[1],coordinates[2]),coordinates]

        subVolume = subvolume(scoreVolume,coordinates[0]-cubeStart,coordinates[1]-cubeStart,0,cubeLength,cubeLength,1)
    else:
        if (coordinates[0]-cubeStart < 1 or coordinates[1]-cubeStart < 1 or coordinates[2]-cubeStart < 1) or \
                (coordinates[0]-cubeStart + cubeLength >= sizeX or coordinates[1]-cubeStart + cubeLength >= sizeY or \
                 coordinates[2]-cubeStart + cubeLength >= sizeZ):
            if verbose:
                print ("SubPixelPeak: position too close to border for sub-pixel")
            return [scoreVolume(coordinates[0],coordinates[1],coordinates[2]),coordinates]

        subVolume = subvolume(scoreVolume,coordinates[0]-cubeStart,coordinates[1]-cubeStart,coordinates[2]-cubeStart,
                              cubeLength,cubeLength,cubeLength)
    
    #size of interpolated volume
    scaleSize = 10*cubeLength

    #ratio between interpolation area and large volume
    scaleRatio = 1.0 * cubeLength / scaleSize
    
    #resize into bigger volume
    if interpolation=='Spline':
        if twoD:
            subVolumeScaled = vol(scaleSize,scaleSize,1)
        else:
            subVolumeScaled = vol(scaleSize,scaleSize,scaleSize)
        rescaleSpline(subVolume,subVolumeScaled)
    else:
        subVolumeScaled = resize(volume=subVolume, factor=10)[0]
    
    peakCoordinates = peak(subVolumeScaled)
    
    peakValue = subVolumeScaled(peakCoordinates[0],peakCoordinates[1],peakCoordinates[2])
    
    #calculate sub pixel coordinates of interpolated peak
    peakCoordinates[0] = peakCoordinates[0]*scaleRatio - cubeStart + coordinates[0]
    peakCoordinates[1] = peakCoordinates[1]*scaleRatio - cubeStart + coordinates[1]
    if twoD:
        peakCoordinates[2] = 0
    else:
        peakCoordinates[2] = peakCoordinates[2]*scaleRatio - cubeStart + coordinates[2]
    if ( peakCoordinates[0] > scoreVolume.sizeX() or peakCoordinates[1] > scoreVolume.sizeY() or
            peakCoordinates[2] > scoreVolume.sizeZ() ):
        if verbose:
            print ("SubPixelPeak: peak position too large :( return input value")
        #something went awfully wrong here. return regular value 
        return [scoreVolume(coordinates[0],coordinates[1],coordinates[2]),coordinates]
    
    return [peakValue,peakCoordinates]

def dev(volume, template, mask=None, volumeIsNormalized=False):
    """
    dev: Calculates the squared deviation of volume and template in real space
    @param volume: A volume
    @type volume:  L{pytom_volume.vol}
    @param template: A template that is searched in volume. Must be of same size as volume.
    @type template:  L{pytom_volume.vol}
    @param mask: mask to constrain correlation
    @type mask: L{pytom_volume.vol}
    @param volumeIsNormalized: speed up if volume is already normalized
    @type volumeIsNormalized: L{bool}
    @return: deviation
    @raise exception: Raises a runtime error if volume and template have a different size.
    @author: FF
    """

    if gpu:
        import cupy as xp
    else:
        import numpy as xp

    from pytom_volume import vol,sum
    from pytom.tools.macros import volumesSameSize

    assert type(volume) == vol, "dev: volume has to be of type vol!"
    assert type(template) == vol, "dev: template has to be of type vol!"
    if not volumesSameSize(volume,template):
        raise RuntimeError('Volume and template must have same size!')

    if not mask:
        p = volume.numelem()
        result = volume-template
    else:
        assert type(mask) == vol, "dev: mask has to be of type vol!"
        p = sum(mask)
        result = (volume-template)*mask

    deviat = sum(result**2)
    deviat = deviat / float(p)

    return deviat

def generate_random_phases_3d(shape, reduced_complex=True):
    '''this function generates a set of random phases (between -pi and pi), optionally centrosymmetric
    @shape: shape of array
    @type: tuple
    @reduced_complex: is the shape in reduced complex format
    @type: bool, default=True

    @return: returns an array with values between -pi and pi
    @rtype: ndarray
    '''
    from numpy.random import ranf
    from numpy import pi, rot90
    from numpy.fft import fftshift

    if len(shape) == 3:
        dx, dy, dz = shape
    else:
        dx, dy = shape
        dz = max(dx, dy)

    rnda = (ranf(shape) * pi * 2) - pi

    cc = dx // 2
    ccc = (dx - 1) // 2

    loc = (dz // 2) * (reduced_complex == False)
    centralslice = rnda[:, :, loc]

    centralslice[cc, cc - ccc:cc] = centralslice[cc, -ccc:][::-1] * -1
    centralslice[cc - ccc:cc, cc - ccc:] = rot90(centralslice[-ccc:, cc - ccc:], 2) * -1

    rnda[:, :, loc] = fftshift(centralslice) if reduced_complex else centralslice

    return rnda

def randomizePhaseBeyondFreq(volume, frequency):
    '''This function randomizes the phases beyond a given frequency, while preserving the Friedel symmetry.
    @param volume: target volume
    @type volume: 3d ndarray
    @param frequency: frequency in pixel beyond which phases are randomized.
    @type frequency: int

    @return: 3d volume.
    @rtype: 3d ndarray
    @author: GvdS'''

    from numpy.fft import ifftn, fftshift, fftn, rfftn, irfftn
    from numpy import angle, abs, exp, meshgrid, arange, sqrt, pi, rot90



    threeD = len(volume.shape) == 3
    twoD = len(volume.shape) == 2

    if not twoD and not threeD:
        raise Exception('Invalid volume dimensions: either supply 3D or 2D ndarray')

    if twoD:
        dx, dy = volume.shape
    else:
        dx, dy, dz = volume.shape

    ft = rfftn(volume)
    phase = angle(ft)

    amplitude = abs(ft)
    rnda = generate_random_phases_3d(amplitude.shape, reduced_complex=True)  # (ranf((dx,dy,dz//2+1)) * pi * 2) - pi

    if twoD:
        X, Y = meshgrid(arange(-dx // 2, dx // 2 + dx % 2), arange(-dy // 2, dy // 2 + dy % 2))
        RF = sqrt(X ** 2 + Y ** 2).astype(int)
        R = fftshift(RF)[:, :dy // 2 + 1]
    else:
        X, Y, Z = meshgrid(arange(-dx // 2, dx // 2), arange(-dy // 2, dy // 2),
                           arange(-dz // 2, dz // 2))
        RF = sqrt(X ** 2 + Y ** 2 + Z ** 2)  # .astype(int)
        R = fftshift(RF)[:, :, :dz // 2 + 1]
        # centralslice = fftshift(rnda[:,:,0])
        # cc = dx//2
        # ccc= (dx-1)//2
        # centralslice[cc, cc-ccc:cc] = centralslice[cc,-ccc:][::-1]*-1
        # centralslice[cc-ccc:cc,cc-ccc:] = rot90(centralslice[-ccc:,cc-ccc:],2)*-1
        # rnda[:,:,0] = fftshift(centralslice)

    rnda[R <= frequency] = 0
    phase[R > frequency] = 0
    phase += rnda
    image = irfftn((amplitude * exp(1j * phase)), s=volume.shape)

    if abs(image.imag).sum() > 1E-8:
        raise Exception('Imaginary part is non-zero. Failed to centro-summetrize the phases.')

    return image.real

def calc_FSC_true(FSC_t, FSC_n):
    '''Calculates the true FSC as defined in Henderson
    @param FSC_t: array with FSC values without randomized phases.
    @type FSC_t: ndarray
    @param FSC_n: array with FSC values without randomized phases.
    @type FSC_n: ndarray

    @return: return the FSC_true array
    @rtype: ndarray
    '''

    from numpy import zeros_like

    if len(FSC_t) != len(FSC_n):
        raise Exception('FSC arrays are not of equal length.')
    FSC_true = zeros_like(FSC_t)
    for i in range(len(FSC_t)):

        if abs(FSC_t[i] - FSC_n[i]) < 1e-1:
            FSC_true[i] = FSC_t[i]
        else:
            FSC_true[i] = (FSC_t[i] - FSC_n[i]) / (1 - FSC_n[i])

    return FSC_true
