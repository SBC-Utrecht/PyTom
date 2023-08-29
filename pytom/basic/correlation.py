def xcc(volume,template,mask=None, volume_is_normalized=False):
    """
    xcc: Calculates the cross correlation coefficient in real space
    @param volume: A volume
    @type volume:  L{pytom.lib.pytom_volume.vol}
    @param template: A template that is searched in volume. Must be of same size as volume.
    @type template:  L{pytom.lib.pytom_volume.vol}
    @param mask: mask to constrain correlation
    @type mask: L{pytom.lib.pytom_volume.vol}
    @param volume_is_normalized: only used for compatibility with nxcc - not used
    @type volume_is_normalized: L{bool}
    @return: A unscaled value
    @raise exception: Raises a runtime error if volume and template have a different size.  
    @author: Thomas Hrabe 
    """
    from pytom.lib.pytom_volume import sum
    from pytom.tools.macros import volumesSameSize
    
    if not volumesSameSize(volume,template):
        raise RuntimeError('Volume and template must have same size!')
   
    if mask: # mask is given
        result = mask * volume * template
    else:
        result = volume * template
    
    cc = sum(result)
    
    cc = cc / float(volume.numelem())
    
    return cc 
    
def nxcc(volume, template, mask=None, volume_is_normalized=False):
    """
    nxcc: Calculates the normalized cross correlation coefficient in real space
    @param volume: A volume
    @type volume:  L{pytom.lib.pytom_volume.vol}
    @param template: A template that is searched in volume. Must be of same size as volume.
    @type template:  L{pytom.lib.pytom_volume.vol}
    @param mask: mask to constrain correlation
    @type mask: L{pytom.lib.pytom_volume.vol}
    @param volume_is_normalized: speed up if volume is already normalized
    @type volume_is_normalized: L{bool}
    @return: A value between -1 and 1
    @raise exception: Raises a runtime error if volume and template have a different size.
    @author: Thomas Hrabe 
    @change: flag for pre-normalized volume, FF
    """

    from pytom.lib.pytom_volume import vol,sum,limit
    from pytom.tools.macros import volumesSameSize
    
    if not volumesSameSize(volume,template):
        raise RuntimeError('Volume and template must have same size!')
    
    if not mask:
        from pytom.basic.normalise import mean0std1
        if not volume_is_normalized:
           v = mean0std1(volume, True)
        t = mean0std1(template, True)
        p = volume.numelem()
        result = v*t
    else:
        from pytom.lib.pytom_numpy import vol2npy
        from pytom.basic.normalise import normaliseUnderMask
        if not volume_is_normalized:
            (v,p) = normaliseUnderMask(volume, mask)
            (t,p) = normaliseUnderMask(template, mask, p)
            t = t * mask # multiply with the mask
            result = v * t
        else:
            (t,p) = normaliseUnderMask(template,mask)
            t = t * mask # multiply with the mask
            result = volume * t
    
    ncc = sum(result)
    ncc = ncc / float(p)

    return ncc 

    
def xcf(volume, template, mask=None, std_v=None):
    """
    XCF: returns the non-normalised cross correlation function. The xcf 
    result is scaled only by the square of the number of elements.

    @param volume : The search volume
    @type volume: L{pytom.lib.pytom_volume.vol}
    @param template : The template searched (this one will be used for conjugate complex multiplication)
    @type template: L{pytom.lib.pytom_volume.vol}
    @param mask: changed: will be used if specified
    @type mask: L{pytom.lib.pytom_volume.vol}
    @param std_v: Will be unused, only for compatibility reasons with flcf 
    @return: XCF volume
    @rtype: L{pytom.lib.pytom_volume.vol}
    @author: Thomas Hrabe
    """
    import pytom.lib.pytom_volume as pytom_volume
    from pytom.basic import fourier

    if mask != None:
        volume = volume * mask
        template = template * mask

    #determine fourier transforms of volumes
    if volume.__class__ == pytom_volume.vol:
        from pytom.basic.fourier import fft
        fvolume = fft(volume)
    else:
        fvolume = volume
    
    if template.__class__ == pytom_volume.vol:
        from pytom.basic.fourier import fft
        ftemplate = fft(template)
    else:
        ftemplate = template
    
    #perform element wise - conjugate multiplication    
    pytom_volume.conjugate(ftemplate)
    fresult = fvolume * ftemplate

    #transform back to real space
    result = fourier.ifft(fresult)
    
    fourier.iftshift(result)

    n = result.numelem()
    if mask:
        n1 = pytom_volume.sum(mask)
        # real -> FFT requires n1, FFT -> real n
        result.shiftscale(0,1/float(n1*n))
    else:
        result.shiftscale(0,1/float(n*n))
    
    return result

def norm_xcf(volume,template,mask=None, std_v=None):
    """
    nXCF: returns the normalised cross correlation function. Autocorrelation 
    of two equal objects would yield a max nxcf peak of 1.

    @param volume: The search volume
    @param template: The template searched (this one will be used for conjugate complex multiplication)
    @type template: L{pytom.lib.pytom_volume.vol}
    @param mask: template mask. If not given, a default sphere mask will be generated which has the same size with the given template.
    @type mask: L{pytom.lib.pytom_volume.vol}
    @param std_v: Will be unused, only for compatibility reasons with flcf
    @return: the calculated norm_xcf volume
    @rtype: L{pytom.lib.pytom_volume.vol}
    @author: Thomas Hrabe
    @change: masking of template implemented
    """
    from pytom.basic.normalise import mean0std1

    if mask == None:
        result = xcf(mean0std1(volume,True),mean0std1(template,True), mask=None, std_v=None)
    else:
        from pytom.basic.normalise import normaliseUnderMask
        result = xcf(normaliseUnderMask(volume=volume, mask=mask, p=None)[0],
                     normaliseUnderMask(volume=template, mask=mask, p=None)[0],
                     mask=mask, std_v=None)
    #n = result.numelem()
    #result.shiftscale(0,1/float(n*n))

    return result


def meanValueUnderMask(volume, mask=None, p=None):
    """
    meanValueUnderMask: Determines the mean value under a mask
    @param volume: The volume
    @type volume:  L{pytom.lib.pytom_volume.vol}
    @param mask:  The mask
    @type mask:  L{pytom.lib.pytom_volume.vol}
    @param p: precomputed number of voxels in mask
    @type p: float
    @return: A value (scalar)
    @rtype: single
    @change: support None as mask, FF 08.07.2014
    """
    from pytom.lib.pytom_volume import vol, sum
  
    if not volume.__class__ == vol:
        raise TypeError("meanValueUnderMask: Volume parameter must be of type vol!")
    
    if mask:
        if not p:
            p = sum(mask)
        resV = volume*mask
        res = sum(resV)/p
    else:
        res = sum(volume)/(volume.size_x()*volume.size_y()*volume.size_z())
    
    return res

def stdValueUnderMask(volume, mask, meanValue, p=None):
    """
    stdValueUnderMask: Determines the std value under a mask

    @param volume: input volume
    @type volume:  L{pytom.lib.pytom_volume.vol}
    @param mask: mask
    @type mask:  L{pytom.lib.pytom_volume.vol}
    @param p: non zero value numbers in the mask
    @type p: L{float} or L{int}
    @return: A value
    @rtype: L{float}
    @change: support None as mask, FF 08.07.2014
    """
    from pytom.lib.pytom_volume import sum
    from pytom.lib.pytom_volume import vol, power, variance
    
    assert volume.__class__ == vol
    if mask:
        if not p:
            p = sum(mask)
    else:
        p = volume.size_x()*volume.size_y()*volume.size_z()
    
    squareM = meanValue**2
    squareV = vol(volume.size_x(), volume.size_y(), volume.size_z())
    squareV.copyVolume(volume)
    power(squareV, 2)
    res = meanValueUnderMask(squareV, mask, p)
    res = res - squareM
    if res >= 0:
        res = res**0.5
    else:
        print("Res = %.6f < 0 => standard deviation determination fails :(")
        print("   something went terribly wrong and program has to stop")
        raise ValueError('Program stopped in stdValueUnderMask')

    return res

def meanUnderMask(volume, mask, p):
    """
    meanUnderMask: calculate the mean volume under the given mask (Both should have the same size)
    @param volume: input volume
    @type volume:  L{pytom.lib.pytom_volume.vol}
    @param mask: mask
    @type mask:  L{pytom.lib.pytom_volume.vol}
    @param p: non zero value numbers in the mask
    @type p: L{int} or L{float}
    @return: the calculated mean volume under mask
    @rtype:  L{pytom.lib.pytom_volume.vol}
    @author: Yuxiang Chen
    """
    size = volume.numelem()
    
    from pytom.basic.fourier import fft, ifft, iftshift
    from pytom.lib.pytom_volume import conjugate
    # for some reason, this has to be conjugated. (Otherwise the asym mask won't work)
    fMask = fft(mask)
    conjugate(fMask)
    
    result = iftshift(ifft(fMask*fft(volume)))/(size*p)

    return result

def stdUnderMask(volume, mask, p, meanV):
    """
    stdUnderMask: calculate the std volume under the given mask
    @param volume: input volume
    @type volume:  L{pytom.lib.pytom_volume.vol}
    @param mask: mask
    @type mask:  L{pytom.lib.pytom_volume.vol}
    @param p: non zero value numbers in the mask
    @type p: L{int}
    @param meanV: mean volume under mask, which should already been caculated
    @type meanV:  L{pytom.lib.pytom_volume.vol}
    @return: the calculated std volume under mask
    @rtype:  L{pytom.lib.pytom_volume.vol}
    @author: Yuxiang Chen
    """
    from pytom.basic.fourier import fft, ifft, iftshift
    from pytom.lib.pytom_volume import vol, power, limit
    
    copyV = vol(volume.size_x(), volume.size_y(), volume.size_z())
    copyV.copyVolume(volume)
    power(copyV, 2)  #calculate the square of the volume
    
    copyMean = vol(meanV.size_x(), meanV.size_y(), meanV.size_z())
    copyMean.copyVolume(meanV)
    power(copyMean, 2)

    result = meanUnderMask(copyV, mask, p) - copyMean

    # from pytom.lib.pytom_volume import abs
    # abs(result)
    limit(result, 1e-09, 1, 0, 0, True, False)  # this step is needed to set all those value (close to 0) to 1
    power(result, 0.5)

    return result

def flcf(volume, template, mask=None, std_v=None, wedge=1):
    '''
    Created on Apr 13, 2010
    FLCF: compute the fast local correlation function
    This functions works only when the mask is in the middle.
    
    @param volume: target volume
    @type volume: L{pytom.lib.pytom_volume.vol}
    @param template: template to be searched. It can have smaller size then target volume.
    @type template: L{pytom.lib.pytom_volume.vol}
    @param mask: template mask. If not given, a default sphere mask will be generated which has the same size with the given template.
    @type mask: L{pytom.lib.pytom_volume.vol}
    @param std_v: standard deviation of the target volume under mask, which do not need to be calculated again when the mask is identical.
    @type std_v: L{pytom.lib.pytom_volume.vol}
    @return: the local correlation function
    @rtype: L{pytom.lib.pytom_volume.vol}
    
    @author: Yuxiang Chen
    '''
    from pytom.lib.pytom_volume import vol, pasteCenter, conjugate, sum
    from pytom.basic.fourier import fft, ifft, iftshift

    if volume.__class__ != vol and template.__class__ != vol:
        raise RuntimeError('Wrong input type!')
    
    if volume.size_x()<template.size_x() or volume.size_y()<template.size_y() or volume.size_z()<template.size_z():
        raise RuntimeError('Template size is bigger than the target volume')

    # generate the mask 
    if mask.__class__ != vol:
        from pytom.lib.pytom_volume import initSphere
        mask = vol(template.size_x(), template.size_y(), template.size_z())
        mask.setAll(0)
        initSphere(mask, template.size_x()/2,0,0,template.size_x()/2, template.size_y()/2, template.size_z()/2)
    else:
        if template.size_x()!=mask.size_x() and template.size_y()!=mask.size_y() and template.size_z()!=mask.size_z():
            raise RuntimeError('Template and mask size are not consistent!')
    
    # calculate the non-zeros
    p = sum(mask)

    # normalize the template under mask
    meanT = meanValueUnderMask(template, mask, p)
    stdT = stdValueUnderMask(template, mask, meanT, p)

    temp = (template - meanT)/stdT
    temp = temp * mask

    # construct both the template and the mask which has the same size as target volume
    tempV = temp
    if volume.size_x() != temp.size_x() or volume.size_y() != temp.size_y() or volume.size_z() != temp.size_z():
        tempV = vol(volume.size_x(), volume.size_y(), volume.size_z())
        tempV.setAll(0)
        pasteCenter(temp, tempV)
    
    maskV = mask
    if volume.size_x() != mask.size_x() or volume.size_y() != mask.size_y() or volume.size_z() != mask.size_z():
        maskV = vol(volume.size_x(), volume.size_y(), volume.size_z())
        maskV.setAll(0)
        pasteCenter(mask, maskV)
    
    # calculate the mean and std of volume
    if std_v.__class__ != vol:
        meanV = meanUnderMask(volume, maskV, p)
        std_v = stdUnderMask(volume, maskV, p, meanV)

    size = volume.numelem()
    fT = fft(tempV)
    conjugate(fT)
    result = iftshift(ifft(fT*fft(volume)))/std_v

    result.shiftscale(0, 1/(size*p))
    
    return result


def band_cc(volume,reference,band,verbose = False):
    """
    band_cc: Determines the normalised correlation coefficient within a band
    @param volume: The volume
    @type volume: L{pytom.lib.pytom_volume.vol}
    @param reference: The reference
    @type reference: L{pytom.lib.pytom_volume.vol}
    @param band: [a,b] - specify the lower and upper end of band.
    @return: First parameter - The correlation of the two volumes in the specified band. 
             Second parameter - The bandpass filter used.
    @rtype: List - [float,L{pytom.lib.pytom_freqweight.weight}]
    @author: Thomas Hrabe    
    """
    import pytom.lib.pytom_volume as pytom_volume
    from pytom.basic.filter import bandpassFilter
    from math import sqrt
    
    if verbose:
        print('lowest freq : ', band[0],' highest freq' , band[1])
        
    vf = bandpassFilter(volume,band[0],band[1],fourierOnly=True)
    rf = bandpassFilter(reference,band[0],band[1],vf[1],fourierOnly=True)
    
    ccVolume = pytom_volume.vol_comp(rf[0].size_x(),rf[0].size_y(),rf[0].size_z())
    ccVolume.copyVolume(rf[0])
    
    pytom_volume.conj_mult(ccVolume,vf[0])
    
    cc = pytom_volume.sum(ccVolume)

    cc = cc.real
    
    v = vf[0]
    r = rf[0]
    
    absV = pytom_volume.abs(v)
    absR = pytom_volume.abs(r)
    
    pytom_volume.power(absV,2)
    pytom_volume.power(absR,2)
    
    sumV = pytom_volume.sum(absV)
    sumR = pytom_volume.sum(absR)
    
    sumV = abs(sumV)
    sumR = abs(sumR)
    
    if sumV == 0:
        sumV =1
        
    if sumR == 0:
        sumR =1
        
    cc = cc / (sqrt(sumV*sumR)) 
    
    #numerical errors will be punished with nan
    if abs(cc) > 1.1 :
        cc = float('nan')
    
    return [cc,vf[1]];
    
def weighted_xcc(volume,reference,number_of_bands,wedge_angle=-1):
        """
        weighted_xcc: Determines the band weighted correlation coefficient for a volume and reference. Notation according Steward/Grigorieff paper
        @param volume: A volume
        @type volume: L{pytom.lib.pytom_volume.vol}
        @param reference: A reference of same size as volume
        @type reference: L{pytom.lib.pytom_volume.vol}
        @param number_of_bands: Number of bands
        @param wedge_angle: A optional wedge angle
        @return: The weighted correlation coefficient
        @rtype: float  
        @author: Thomas Hrabe   
        """    
        import pytom.lib.pytom_volume as pytom_volume
        from pytom.basic.fourier import fft
        from math import sqrt
        import pytom.lib.pytom_freqweight as pytom_freqweight
        result = 0
        numberVoxels = 0
        
        #volume.write('vol.em');
        #reference.write('ref.em');
        fvolume = fft(volume)
        freference = fft(reference)
        numelem = volume.numelem()
        
        fvolume.shiftscale(0,1/float(numelem))
        freference.shiftscale(0,1/float(numelem))
        from pytom.basic.structures import WedgeInfo
        wedge = WedgeInfo(wedge_angle)
        wedgeVolume = wedge.returnWedgeVolume(volume.size_x(),volume.size_y(),volume.size_z())
        
        increment = int(volume.size_x()/2 * 1/number_of_bands)
        band = [0,100]
        for i in range(0,volume.size_x()/2, increment):
        
            band[0] = i
            band[1] = i + increment
    
            r = band_cc(volume,reference,band)
            cc = r[0]
            
            #print cc;
            filter = r[1]
            #get bandVolume
            bandVolume = filter.getWeightVolume(True)
            
            filterVolumeReduced = bandVolume * wedgeVolume
            filterVolume = pytom_volume.reducedToFull(filterVolumeReduced)
            
            #determine number of voxels != 0    
            N = pytom_volume.numberSetVoxels(filterVolume)
            
            w = sqrt(1/float(N))
            
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
    
    
    
def fsc_sum(volume,reference,number_of_bands,wedge_angle=-1):
    """
    fsc_sum: Determines the sum of the Fourier Shell Correlation coefficient for a volume and reference. 
    @param volume: A volume
    @type volume: L{pytom.lib.pytom_volume.vol}
    @param reference: A reference of same size as volume
    @type reference: L{pytom.lib.pytom_volume.vol}
    @param number_of_bands: Number of bands
    @param wedge_angle: A optional wedge angle
    @return: The sum FSC coefficient
    @rtype: float  
    @author: Thomas Hrabe   
    """    
    
    from pytom.basic.correlation import band_cc
    import pytom.lib.pytom_volume as pytom_volume
    from pytom.basic.fourier import fft
    from math import sqrt
    import pytom.lib.pytom_freqweight as pytom_freqweight
    result = 0
    numberVoxels = 0
    
    #volume.write('vol.em');
    #reference.write('ref.em');
    fvolume = fft(volume)
    freference = fft(reference)
    numelem = volume.numelem()
    
    fvolume.shiftscale(0,1/float(numelem))
    freference.shiftscale(0,1/float(numelem))

    #print '-----'
    for i in range(number_of_bands):
        #process bandCorrelation
        band = []
        band[0] = i*volume.size_x()/number_of_bands 
        band[1] = (i+1)*volume.size_x()/number_of_bands
        
        r = band_cc(fvolume,freference,band)
        cc = r[0]
        #print cc
        result = result + cc
    #print '-----'
    
    return result*(1/float(number_of_bands))

def band_cf(volume,reference,band=[0,100]):
    """
    band_cf:
    @param volume: The volume
    @param reference: The reference
    @param band: [a,b] - specify the lower and upper end of band. [0,1] if not set.
    @return: First parameter - The correlation of the two volumes in the specified ring. 
             Second parameter - The bandpass filter used.
    @rtype: List - [L{pytom.lib.pytom_volume.vol},L{pytom.lib.pytom_freqweight.weight}]
    @author: Thomas Hrabe   
    @todo: does not work yet -> test is disabled
    """
    import pytom.lib.pytom_volume as pytom_volume
    from math import sqrt
    from pytom.basic import fourier
    from pytom.basic.filter import bandpassFilter
    from pytom.basic.correlation import norm_xcf
    
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

def weighted_xcf(volume,reference,number_of_bands,wedge_angle=-1):
    """
    weighted_xcf: Determines the weighted correlation function for volume and reference
    @param volume: A volume 
    @param reference: A reference 
    @param number_of_bands:Number of bands
    @param wedge_angle: A optional wedge angle
    @return: The weighted correlation function 
    @rtype: L{pytom.lib.pytom_volume.vol} 
    @author: Thomas Hrabe 
    @todo: does not work yet -> test is disabled
    """
    from pytom.basic.correlation import band_cf
    import pytom.lib.pytom_volume as pytom_volume
    from math import sqrt
    import pytom.lib.pytom_freqweight as pytom_freqweight

    result = pytom_volume.vol(volume.size_x(),volume.size_y(),volume.size_z())
    result.setAll(0)
    cc2 = pytom_volume.vol(volume.size_x(),volume.size_y(),volume.size_z())
    cc2.setAll(0)
    q = 0
    
    if wedge_angle >=0:
        wedgeFilter = pytom_freqweight.weight(wedge_angle,0,volume.size_x(),volume.size_y(),volume.size_z())
        wedgeVolume = wedgeFilter.getWeightVolume(True)
    else:
        wedgeVolume = pytom_volume.vol(volume.size_x(), volume.size_y(), int(volume.size_z()/2+1))
        wedgeVolume.setAll(1.0)
        
    w = sqrt(1/float(volume.size_x()*volume.size_y()*volume.size_z()))
    
    numberVoxels = 0
    
    for i in range(number_of_bands):
        """
        notation according Steward/Grigorieff paper
        """
        band = [0,0]
        band[0] = i*volume.size_x()/number_of_bands 
        band[1] = (i+1)*volume.size_x()/number_of_bands
        
        r = band_cf(volume,reference,band)
        
        cc = r[0]
                
        filter = r[1]
        #get bandVolume
        bandVolume = filter.getWeightVolume(True)
            
        filterVolumeReduced = bandVolume * wedgeVolume
        filterVolume = pytom_volume.reducedToFull(filterVolumeReduced)
        #determine number of voxels != 0    
        N = pytom_volume.numberSetVoxels(filterVolume)
            
        #add to number of total voxels
        numberVoxels = numberVoxels + N
                 
        cc2.copyVolume(r[0])
                
        pytom_volume.power(cc2,2)
        
        cc.shiftscale(w,1)
        ccdiv = cc2/(cc)
                
        pytom_volume.power(ccdiv,3)
        
        #abs(ccdiv); as suggested by grigorief
        ccdiv.shiftscale(0,N)
        
        result = result + ccdiv
    
    result.shiftscale(0,1/float(numberVoxels))
    
    return result

def fsc(volume1,volume2,number_bands,mask=None,verbose=False, filename=None):
    """
    FSC - Calculates the Fourier Shell Correlation for two volumes
    @param volume1: volume one
    @type volume1: L{pytom.lib.pytom_volume.vol}
    @param volume2: volume two
    @type volume2: L{pytom.lib.pytom_volume.vol}
    @param number_bands: number of shells for FSC
    @type number_bands: int
    @param filename: write FSC to ascii file if specified
    @type filename: string

    @return: Returns a list of cc values 
    @author: Thomas Hrabe  
    @rtype: list[floats]
    """
    
    from pytom.basic.correlation import band_cc
    from pytom.tools.macros import volumesSameSize
    from pytom.basic.structures import Mask
    import pytom.lib.pytom_volume as pytom_volume
    
    if not volumesSameSize(volume1, volume2):
        raise RuntimeError('Volumes must have the same size!')
    
    if mask:
        if mask.__class__ == pytom_volume.vol:
            volume1 = volume1 * mask
            volume2 = volume2 * mask
          
        elif mask.__class__ == Mask:
            mask = mask.getVolume()
            volume1 = volume1 * mask
            volume2 = volume2 * mask
        elif mask.__class__ == str:
            mask = pytom_volume.read(mask)
            volume1 = volume1 * mask
            volume2 = volume2 * mask 
        else:
            raise RuntimeError('FSC: Mask must be a volume OR a Mask object OR a string path to a mask!')  
        
    fscResult = []
    band = [-1,-1]
    if type(number_bands) == tuple:
        number_bands = tuple[0]
    
    increment = int(volume1.size_x()/2 * 1/number_bands)
    
    for i in range(0,volume1.size_x()//2, increment):
        
        band[0] = i
        band[1] = i + increment
        
        if verbose:
            print('Band : ' ,band)
            
        res = band_cc(volume1,volume2,band,verbose)
        
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

def determine_resolution(fsc,resolution_criterion,verbose=False):
    """
    determine_resolution: Determines frequency and band where correlation drops below the resolution_criterion. Uses linear interpolation between two positions
    @param fsc: The fsc list determined by L{pytom.basic.correlation.FSC}
    @param resolution_criterion: A value between 0 and 1
    @return: [resolution,interpolatedBand,number_bands] 
    @author: Thomas Hrabe 
    @todo: Add test! 
    """
    number_bands = len(fsc)
    
    band = number_bands

    for i in range(number_bands):
        if fsc[i] < resolution_criterion:     
            band = i-1  #select the band that is still larger than criterion
            break
    
    if verbose:
        print('Band detected at ', band)
    
    if band == -1:
        raise RuntimeError("Please check your resolution criterion or you FSC!")
            
    elif band < number_bands:
        fsc1 = fsc[band]
        fsc2 = fsc[band+1]
        
        try:
            interpolatedBand = (resolution_criterion-fsc1)/(fsc2-fsc1)+band            
        
        except ZeroDivisionError:
            interpolatedBand = band
        
    else:
        interpolatedBand = band
        
    if verbose:
        print('Band interpolated to ', interpolatedBand)
        
    resolution = (interpolatedBand+1) / float(number_bands)
    
    if resolution < 0 :
        resolution = 1
        interpolatedBand = number_bands
        print('Warning: PyTom determined a resolution < 0 for your data. Please check "mass" in data is positive or negative for all cubes.')
        print('Warning: Setting resolution to 1 and ',interpolatedBand)
        print('')
        
    return [resolution,interpolatedBand,number_bands]


def soc(volume,reference,mask=None, std_v=None):
    """
    soc : Second Order Correlation. Correlation of correlation peaks.
    @param volume: The volume
    @type volume:  L{pytom.lib.pytom_volume.vol}
    @param reference: The reference / template
    @type reference:  L{pytom.lib.pytom_volume.vol}
    @author: Thomas Hrabe   
    """
    referencePeak = flcf(reference,reference,mask)
    peaks = flcf(volume,reference,mask)
    
    return flcf(peaks,referencePeak,mask)


def sub_pixel_peak_parabolic(score_volume, coordinates, verbose=False):
    """
    quadratic interpolation of three adjacent samples
    @param score_volume: The score volume
    @param coordinates: [x,y,z] coordinates where the sub pixel peak will be determined
    @param verbose: be talkative
    @type verbose: bool
    @return: Returns [peakValue,peakCoordinates] with sub pixel accuracy
    """
    from pytom.lib.pytom_volume import vol
    assert type(score_volume) == vol, "score_volume must be vol!"
    if (coordinates[0] == 0 or coordinates[0] == score_volume.size_x()-1 or
                coordinates[1] == 0 or coordinates[1] == score_volume.size_y()-1 or
                coordinates[2] == 0 or coordinates[2] == score_volume.size_z()-1 ):
        if verbose:
            print("sub_pixel_peak_parabolic: peak near borders - no interpolation done")
        return [score_volume.getV(coordinates[0], coordinates[1], coordinates[2]), coordinates]
    peakCoordinates = coordinates
    (x, p1, a1) = qint(ym1=score_volume.getV(coordinates[0]-1, coordinates[1], coordinates[2]), 
                     y0=score_volume.getV(coordinates[0], coordinates[1], coordinates[2]), 
                     yp1=score_volume.getV(coordinates[0]+1, coordinates[1], coordinates[2]))
    (y, p2, a2) = qint(ym1=score_volume.getV(coordinates[0], coordinates[1]-1, coordinates[2]), 
                     y0=score_volume.getV(coordinates[0], coordinates[1], coordinates[2]), 
                     yp1=score_volume.getV(coordinates[0], coordinates[1]+1, coordinates[2]))
    if score_volume.size_z() != 1:
        (z, p3, a3) = qint(ym1=score_volume.getV(coordinates[0], coordinates[1], coordinates[2]-1), 
                         y0=score_volume.getV(coordinates[0], coordinates[1], coordinates[2]), 
                         yp1=score_volume.getV(coordinates[0], coordinates[1], coordinates[2]+1))
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
            
def sub_pixel_peak(score_volume, coordinates, cube_length=8, interpolation='Spline', verbose=False):
    """
    sub_pixel_peak: Will determine the sub pixel area of peak. Utilizes spline, fourier or parabolic interpolation.

    @param verbose: be talkative
    @type verbose: L{str}
    @param score_volume: The score volume
    @param coordinates: [x,y,z] coordinates where the sub pixel peak will be determined
    @param cube_length: length of cube - only used for Spline and Fourier interpolation
    @type cube_length: int (even)
    @param interpolation: interpolation type: 'Spline', 'Quadratic', or 'Fourier'
    @type interpolation: str
    @return: Returns [peakValue,peakCoordinates] with sub pixel accuracy

    last change: 02/07/2013 FF: 2D functionality added
    """
    assert type(interpolation) == str, 'sub_pixel_peak: interpolation must be str'
    if (interpolation.lower() == 'quadratic') or (interpolation.lower() == 'parabolic'):
        (peakValue,peakCoordinates) = sub_pixel_peak_parabolic(score_volume=score_volume, coordinates=coordinates, verbose=verbose)
        return [peakValue,peakCoordinates]
    from pytom.lib.pytom_volume import vol,subvolume,rescaleSpline,peak
    from pytom.basic.transformations import resize
  
    #extend function for 2D
    twoD = False
    if score_volume.size_z() == 1:
        twoD = True

    cubeStart = cube_length//2
    size_x = score_volume.size_x()
    size_y = score_volume.size_y()
    size_z = score_volume.size_z()
    
    if twoD:
        if (coordinates[0]-cubeStart < 1 or coordinates[1]-cubeStart < 1) or\
            (coordinates[0]-cubeStart + cube_length >= size_x or coordinates[1]-cubeStart + cube_length >= size_y):
            if verbose:
                print ("SubPixelPeak: position too close to border for sub-pixel")
            return [score_volume(coordinates[0],coordinates[1],coordinates[2]),coordinates]

        subVolume = subvolume(score_volume,coordinates[0]-cubeStart,coordinates[1]-cubeStart,0,cube_length,cube_length,1)
    else:
        if (coordinates[0]-cubeStart < 1 or coordinates[1]-cubeStart < 1 or coordinates[2]-cubeStart < 1) or \
                (coordinates[0]-cubeStart + cube_length >= size_x or coordinates[1]-cubeStart + cube_length >= size_y or \
                 coordinates[2]-cubeStart + cube_length >= size_z):
            if verbose:
                print ("SubPixelPeak: position too close to border for sub-pixel")
            return [score_volume(coordinates[0],coordinates[1],coordinates[2]),coordinates]

        subVolume = subvolume(score_volume,coordinates[0]-cubeStart,coordinates[1]-cubeStart,coordinates[2]-cubeStart,
                              cube_length,cube_length,cube_length)
    
    #size of interpolated volume
    scaleSize = 10*cube_length

    #ratio between interpolation area and large volume
    scaleRatio = 1.0 * cube_length / scaleSize
    
    #resize into bigger volume
    if interpolation=='Spline':
        if twoD:
            subVolumeScaled = vol(scaleSize,scaleSize,1)
        else:
            subVolumeScaled = vol(scaleSize,scaleSize,scaleSize)
        rescaleSpline(subVolume,subVolumeScaled)
    else:
        subVolumeScaled = resize(volume=subVolume, factor=10)[0]

    # this returns the index of the peak
    peakCoordinates = peak(subVolumeScaled)
    
    peakValue = subVolumeScaled(peakCoordinates[0],peakCoordinates[1],peakCoordinates[2])
    
    #calculate sub pixel coordinates of interpolated peak
    peakCoordinates[0] = peakCoordinates[0]*scaleRatio - cubeStart + coordinates[0]
    peakCoordinates[1] = peakCoordinates[1]*scaleRatio - cubeStart + coordinates[1]
    if twoD:
        peakCoordinates[2] = 0
    else:
        peakCoordinates[2] = peakCoordinates[2]*scaleRatio - cubeStart + coordinates[2]
    if ( peakCoordinates[0] > score_volume.size_x() or peakCoordinates[1] > score_volume.size_y() or
            peakCoordinates[2] > score_volume.size_z() ):
        if verbose:
            print ("SubPixelPeak: peak position too large :( return input value")
        #something went awfully wrong here. return regular value 
        return [score_volume(coordinates[0],coordinates[1],coordinates[2]),coordinates]
    
    return [peakValue,peakCoordinates]

