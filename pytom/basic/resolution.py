'''
Created on Feb 14, 2011

@author: hrabe
'''


def saveForFSC(newReferenceName,particleList,verbose=False):
    """
    saveForFSC: Split particles into even / odd and calculate the Fourier Shell Correlation out of these two sets.
    @param newReferenceName: Name of current alignment result
    @type newReferenceName: str
    @param particleList:
    @type particleList: L{pytom.basic.structures.ParticleList} 
    """
    if len(particleList) < 2:
        raise RuntimeError('ParticleList must have at least 2 elements to run saveForFSC!')
    
    from pytom.basic.structures import ParticleList
    from pytom.alignment.alignmentFunctions import average
    
    even = ParticleList('/')
    odd = ParticleList('/')
    
    for particleCounter in range(len(particleList)):
        
        if particleCounter % 2 == 0:
            even.append(particleList[particleCounter])
        else:
            odd.append(particleList[particleCounter])
    
    if verbose:
        print('averaging even:')
    average(even,newReferenceName + 'even.em',verbose)
    if verbose:
        print('averaging odd:')
        
    average(odd,newReferenceName + 'odd.em',verbose)

            

def bandToAngstrom(band,pixelSize,numberOfBands,upscale=1):
    """
    bandToAngstrom: Calculate resolution in Angstrom
    @param band: The current resolution in pixels (python convention, i.e., start from 0)
    @param pixelSize: Pixelsize in original sized subtomogram 
    @param numberOfBands: Number of bands used
    @param upscale: Was band determined for downscaled tomograms? Scale band up to original size.
    @author: Thomas Hrabe 
    """
    
    if band < 0:
        raise RuntimeError('pytom.basic.resolution.bandToAngstrom: Band parameter < 0. Check your parameters!')
    
    frequency = band * upscale
    
    resolution = 2 * pixelSize * numberOfBands / (frequency+1)
    
    return resolution
    
def angstromToBand(resolution, pixelSize,numberOfBands,scale = 1):
    """
    angstromToBand: 
    @param resolution:
    @param pixelSize:
    @param numberOfBands:   
    @param scale: 
    """
    
    if resolution <= 0:
        raise RuntimeError('pytom.basic.resolution.angstromToBand: Resolution <= 0. Check your parameters!')
    
    band = 2 * pixelSize * numberOfBands / (resolution * scale) -1
    
    return band
    
def angleFromResolution(resolution,particleDiameter=-1):
    """
    angleFromResolution: determines the smallest possible rotation angle for a current resolution (see parameters) 
    @param resolution: Resolution (in Angstrom)
    @param particleDiameter: Diameter of object. (in Angstrom)
    @author: Thomas Hrabe    
    """
    from pytom.angles.angle import rad2deg
    
    if resolution <= 0:
        raise RuntimeError('Resolution parameter is invalid')
    
    if particleDiameter <= 0 :
        raise RuntimeError('ParticleDiameter seems to be wrong!')
    
    angle = rad2deg * resolution / particleDiameter
    
    return angle

def getResolutionBandFromFSC(fsc, criterion=0.143):
    """
    @param fsc: fourier shell correlation (1-d)
    @type fsc: L{list}
    @param criterion: fsc value for resolution determination
    @type criterion: L{float}
    @return: band at crossing of criterion
    @rtype: L{float}
    @author: FF
    """
    crossed = False
    resband = len(fsc)-1
    for (ii,fscval) in enumerate(fsc[:-4]):
        if fscval < criterion and crossed==False:
            resband = ii
            crossed = True
        if crossed and fscval > criterion:
            crossed = False
    # interpolate crossing
    fsc1 = fsc[resband]
    fsc2 = fsc[resband-1]
    return resband - 1. + (fsc2-criterion)/(fsc2-fsc1)

def write_fsc2Ascii(fsc, filename, fsc_rand=None, fsc_corr=None):
    """
    write fsc array to ascii file
    @param fsc: Fourier shell correlation
    @type fsc: 1-d array
    @param filename: Name of generated ascii file
    @type filename: L{str}
    """
    fh = open( filename, "w")

    fsc_rand = ['',]*len(fsc) if fsc_rand is None else fsc_rand
    fsc_corr = ['',]*len(fsc) if fsc_corr is None else fsc_corr

    for fscval, fscrandval, fsccorrval in zip(fsc, fsc_rand, fsc_corr):
        fh.write(f'{fscval}\t{fscrandval}\t{fsccorrval}\n')
    fh.close()

def read_fscFromAscii(filename):
    """
    read fsc array from ascii file
    @param filename: Name of generated ascii file
    @type filename: L{str}
    @return: fsc
    @rtype: 1-d array
    """
    fsc = []
    fh = open( filename, "r")
    tmp = str(1.0)
    while len(tmp) > 0:
        tmp = fh.readline()
        if len(tmp) > 0:
            fsc.append(float(tmp))
    return fsc




