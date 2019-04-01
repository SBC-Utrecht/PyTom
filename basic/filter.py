def profile2FourierVol( profile, dim=None, reduced=False):
    """
    create Volume from 1d radial profile, e.g., to modulate signal with \
    specific function such as CTF or FSC. Simple linear interpolation is used\
    for sampling.

    @param profile: profile
    @type profile: 1-d L{pytom_volume.vol} or 1-d python array
    @param dim: dimension of (cubic) output
    @type dim: L{int}
    @param reduced: If true reduced Fourier representation (N/2+1, N, N) is generated.
    @type reduced: L{bool}

    @return: 3-dim complex volume with spherically symmetrical profile
    @rtype: L{pytom_volume.vol}
    @author: FF
    """
    from pytom_volume import vol
    from math import sqrt

    if not dim:
        if type(profile) == vol:
            dim = 2*profile.sizeX()
        else:
            dim = 2*len(profile)

    if reduced:
        nx = int(dim/2)+1
    else:
        nx = dim

    fkernel = vol( nx, dim, dim)
    if type(profile) == vol:
        r_max = profile.sizeX()-1
    else:
        r_max = len(profile)-1

    if reduced:
        centx = 0
    else:
        centx = int(dim/2)
        centy = centx
        centz = centx
    for ix in range(0,nx):
        for iy in range(0,dim):
            if reduced:
                if iy < nx-1:
                    centy=0
                else:
                    centy=dim-1
            for iz in range(0,dim):
                if reduced:
                    if iz < nx-1:
                        centz=0
                    else:
                        centz=dim-1
                r = (ix-centx)**2 + (iy-centy)**2 + (iz-centz)**2
                r = sqrt(r)
                if r > r_max:
                    fkernel.setV( 0., ix, iy, iz)
                else:
                    ir = int(r)
                    l1 = r-ir
                    if l1 > 0.:
                        l2 = 1.-l1
                        if type(profile) == vol:
                            val = l2*profile.getV(ir,0,0)+l1*profile.getV(ir+1,0,0)
                        else:
                            val = l2*profile[ir]+l1*profile[ir+1]
                    else:
                        if type(profile) == vol:
                            val = profile.getV(ir,0,0)
                        else:
                            val = profile[ir]
                    fkernel.setV( val, ix, iy, iz)

    return fkernel

def filter_volume_by_profile( volume, profile):
    """
    filter volume by 1-d profile
    @param volume: volume
    @type volume: L{pytom_volume.vol}
    @param profile: 1-d profile
    @type profile: L{pytom_volume.vol}
    @return: outvol
    @rtype: L{pytom_volume.vol}
    @author: FF
    """
    from pytom.basic.filter import profile2FourierVol
    from pytom.basic.fourier import convolute, powerspectrum

    kernel = profile2FourierVol( profile=profile, dim=volume.sizeX(), reduced=False)
    outvol = convolute(v=volume, k=kernel, kernel_in_fourier=True)
    return outvol


def gridCTF(x_array, y_array, z_array):
    """
    function similar to ndgrid in matlab

    @type x_array: 1-dim array
    @type y_array: 1-dim array
    @type z_array: 1-dim array
    @return: 3-dim volumes with x, y, z values in x,y,z dimension, respectively
    """
    from pytom_volume import vol
    
    cubeSizeX = len(x_array)  
    cubeSizeY = len(y_array) 
    cubeSizeZ = len(z_array) 
    
    r = vol(cubeSizeX,cubeSizeY,cubeSizeZ)
    r.setAll(0.0)
    
    for i in range( cubeSizeZ ):        
        for j in range( cubeSizeY ):
            for k in range( cubeSizeX ):
                
                r.setV(x_array[k], k, j, i)    
    
    y = vol(cubeSizeX,cubeSizeY,cubeSizeZ)
    y.setAll(0.0)
    
    for i in range( cubeSizeZ ):        
        for j in range( cubeSizeY ):
            for k in range( cubeSizeX ):                
                y.setV(y_array[k], j, k, i) 
    
    
    z = vol(cubeSizeX,cubeSizeY,cubeSizeZ)
    z.setAll(0.0)
    
    for i in range( cubeSizeZ ):        
        for j in range( cubeSizeY ):
            for k in range( cubeSizeX ):                
                z.setV(z_array[i], k, j, i)
    
    return r, y, z


def volCTF(defocus, x_dim, y_dim, z_dim, pixel_size=None, voltage=None, Cs=None, sigma=None):
    """
    @param defocus: defocus in mu m
    @type defocus: L{float}
    @param x_dim: dimension of volume in x
    @param y_dim: dimension of volume in y
    @param z_dim: dimension of volume in z
    """
    from pytom_volume import vol, power
    from pytom.tools.macros import frange
    import math
    
    
    if Cs == None:
        Cs = 2*(10**(-3))
    else:
        Cs = Cs*(10**(-3))
    
    if voltage == None:
        voltage = 300000
    else:
        voltage = voltage * 1000
    
    if pixel_size == None:
        pixel_size = 0.72*(10**(-9))
    else:
        pixel_size = pixel_size*(10**(-9))
        
        
    Dz = defocus * (10**(-6))    
    Ny = 1/(2*pixel_size)
    
    voltagest = voltage*(1+voltage/1022000)
    lam = ((150.4/voltagest)**0.5)*(10**(-10));
    
#     x_array = frange(-Ny, Ny-(Ny/x_dim), 2*Ny/x_dim)
#     y_array = frange(-Ny, Ny-(Ny/y_dim), 2*Ny/y_dim)
#     z_array = frange(-Ny, Ny-(Ny/z_dim), 2*Ny/z_dim)
    
    x_array = frange(-Ny, Ny, 2*Ny/x_dim)
    y_array = frange(-Ny, Ny, 2*Ny/y_dim)
    z_array = frange(-Ny, Ny, 2*Ny/z_dim)
    
    r, y, z = gridCTF(x_array, y_array, z_array)
    
    power(r, 2)
    power(y, 2)
    power(z, 2)
    r  = r + y + z
    power(r, 0.5)
    
    r4 = vol(r)
    power(r4, 4)
    r2 = vol(r)
    power(r2, 2)
    
    vol = ( r4*Cs*(lam**3) - r2*2*Dz*lam )*math.pi/2
        
    for i in range(vol.sizeX()):
        for j in range(vol.sizeY()):
            for k in range(vol.sizeZ()):
                vol(math.sin(vol(i, j, k)) , i, j, k)    
    
    return vol

def convolutionCTF(volume, defocus, pixelSize=None, voltage=None, Cs=None, sigma=None):
    """
    convolutionCTF:
    @param defocus: Negative value = underfocus (in microns)
    @param pixelSize: Size of pixels / voxels (in Anstroms)
    @param voltage: 
    @param Cs:
    @param sigma:    
    """
    from pytom_volume import subvolume, complexRealMult
    from pytom.basic.fourier import ftshift, fft, ifft
    from pytom.basic.filter import volCTF
    
    dimX = volume.sizeX()
    dimY = volume.sizeY()
    dimZ = volume.sizeZ()

    ctf = volCTF(defocus, dimX, dimY, dimZ, pixelSize, voltage, Cs, sigma)
    filterCTF = subvolume(ftshift(ctf,inplace=False), 0, 0, 0, dimX, dimY, (dimZ/2)+1)
    
    filteredVolume = ifft( complexRealMult(fft(volume), filterCTF) )
    
    return filteredVolume


def fourierFilterShift(filter):
    """
    fourierFilterShift: NEEDS Documentation
    @param filter: NEEDS Documentation 
    """
    from pytom_volume import vol
    
    widthX = filter.sizeX()
    centerX = filter.sizeX()//2
    boxX = filter.sizeX()//2
    
    widthY = filter.sizeY()
    centerY = filter.sizeY()//2
    boxY = filter.sizeY()//2
    
    shifted_filter = vol(widthX, widthY, 1)
    shifted_filter.setAll(0.0)
    
    for i in range(widthX):
        rx = (boxX-i)%widthX
        
        for j in range(widthY):
            ry = (boxY-j)%widthY
                                                            
            shifted_filter.setV(filter.getV(i, j, 0), rx, ry, 0)
                
    return shifted_filter

def fourierFilterShift_ReducedComplex(filter):
    """
    fourierFilterShift: NEEDS Documentation
    @param filter: NEEDS Documentation
    """
    from pytom_volume import vol

    widthX = filter.sizeX()
    centerX = filter.sizeX()//2
    boxX = filter.sizeX()//2

    widthY = filter.sizeY()
    centerY = filter.sizeY()//2
    boxY = filter.sizeY()//2

    shifted_filter = vol(widthX, widthY, 1)
    shifted_filter.setAll(0.0)

    for i in range(widthX):
        rx = (boxX-i)%widthX

        for j in range(widthY):
            ry = (boxY-j)%widthY

            shifted_filter.setV(filter.getV(i, j, 0), rx, widthY-j-1, 0)

    return shifted_filter

def circleFilter(sizeX,sizeY, radiusCutoff):
    """
    circleFilter: NEEDS Documentation
    @param sizeX: NEEDS Documentation 
    @param sizeY: NEEDS Documentation
    @param radiusCutoff: NEEDS Documentation
    """
    from pytom_volume import vol
    
    centerX = sizeX//2
        
    centerY = sizeY//2
    sizeY = (sizeY//2) +1
        
    filter_vol = vol(sizeX, sizeY, 1)
    filter_vol.setAll(0.0)
    
    for i in range(sizeX):        
        for j in range(sizeY):
                        
            radius = ((i-centerX)*(i-centerX)+(j-centerY)*(j-centerY))**0.5
                                                
            if radius <= radiusCutoff:
                filter_vol.setV(1.0, i, j, 0)
                
    return filter_vol


def rampFilter( sizeX, sizeY):
    """
    rampFilter: Generates the weighting function required for weighted backprojection - y-axis is tilt axis

    @param sizeX: size of weighted image in X
    @param sizeY: size of weighted image in Y

    @return: filter volume
    
    """
    from pytom_volume import vol
    from math import exp
    import time

    s = time.time()
    centerX = sizeX//2
    
    centerY = sizeY//2
    sizeY = (sizeY//2) +1
    
    Ny = sizeX//2
        
    filter_vol = vol(sizeX, sizeY, 1)
    filter_vol.setAll(0.0)
    
    for i in range(sizeX):        
        distX = abs(float(i-centerX))
        ratio = distX/Ny
        for j in range(sizeY):            
            filter_vol.setV(ratio, i, j, 0)
    print('ramp filter takes: ', time.time()-s, sizeX, sizeY)
    return filter_vol


def exactFilter(tilt_angles, tiltAngle, sX, sY, sliceWidth, arr=[]):
    """
    exactFilter: Generates the exact weighting function required for weighted backprojection - y-axis is tilt axis
    Reference : Optik, Exact filters for general geometry three dimensional reconstuction, vol.73,146,1986.
    @param tilt_angles: list of all the tilt angles in one tilt series
    @param titlAngle: tilt angle for which the exact weighting function is calculated
    @param sizeX: size of weighted image in X
    @param sizeY: size of weighted image in Y

    @return: filter volume

    """

    from pytom_numpy import npy2vol, vol2npy
    import copy
    from numpy import array, matrix, sin, pi, arange, float32, column_stack, asfortranarray, clip, argmin, ones, ceil
    from pylab import   imshow,show
    from time import time

    s = time()
    sY = sX//2+1

    #arr = matrix(abs(arange(-sX // 2, sX // 2)))
    #w = 1. / (clip(1 - array(matrix(abs(sin((tilt_angles - tiltAngle) * pi / 180.))).T * arr), 0, 2)).sum(axis=0)

    diffAngles = (tilt_angles-tiltAngle)*pi/180.
    spacing = min(abs(diffAngles)[abs(diffAngles)>0.01]) #closest angle to tiltAngle
    crowtherFreq = int(ceil(1 / sin( spacing )))
    arrC = matrix(abs(arange(-crowtherFreq, crowtherFreq)))
    wfuncC = 1. / (clip(1 - array(matrix(abs(sin(diffAngles))).T * arrC), 0, 2)).sum(axis=0)
    wfunc = ones((sX,sY,1),dtype=float32)
    wfunc[sX//2-crowtherFreq:sX//2+crowtherFreq,:, 0] = column_stack(([(wfuncC), ] * (sY))).astype(float32)

    # Create pytom vol object of wfunc.
    fortran_wfunc = asfortranarray(wfunc)
    weightFunc = npy2vol(fortran_wfunc,3)
    rr =copy.deepcopy(vol2npy(weightFunc))

    if abs(rr).sum() < 0.1:
        print('npy2vol failed: no weighting for image.')

    print('Computational time exact weighting: ', time()-s)

    return weightFunc


def rotateWeighting(weighting, z1, z2, x, mask=None, isReducedComplex=None, returnReducedComplex=False, binarize=False):
    """
    rotateWeighting: Rotates a frequency weighting volume around the center. If the volume provided is reduced complex, it will be rescaled to full size, ftshifted, rotated, iftshifted and scaled back to reduced size.
    @param weighting: A weighting volume
    @type weighting: L{pytom_volume.vol}
    @param z1: Z1 rotation angle
    @type z1: float
    @param z2: Z2 rotation angle
    @type z2: float
    @param x: X rotation angle
    @type x: float
    @param mask:=None is there a rotation mask? A mask with all = 1 will be generated otherwise. Such mask should be \
        provided anyway.
    @type mask: L{pytom_volume.vol}
    @param isReducedComplex: Either set to True or False. Will be determined otherwise
    @type isReducedComplex: bool
    @param returnReducedComplex: Return as reduced complex? (Default is False)
    @type returnReducedComplex: bool
    @param binarize: binarize weighting
    @type binarize: bool
    @return: weight as reduced complex volume
    @rtype: L{pytom_volume.vol_comp}
    """
    from pytom_volume import vol, limit, vol_comp
    from pytom_volume import rotate
    assert type(weighting) == vol or  type(weighting) == vol_comp, "rotateWeighting: input neither vol nor vol_comp"
    
    isReducedComplex = isReducedComplex or int(weighting.sizeX()/2)+1 == weighting.sizeZ();

    if isReducedComplex:
        #scale weighting to full size
        from pytom_fftplan import fftShift
        from pytom_volume import reducedToFull
        weighting = reducedToFull(weighting)
        fftShift(weighting, True)

    if not mask:
        mask = vol(weighting.sizeX(),weighting.sizeY(),weighting.sizeZ())
        mask.setAll(1)

    weightingRotated = vol(weighting.sizeX(),weighting.sizeY(),weighting.sizeZ())

    rotate(weighting,weightingRotated,z1,z2,x)
    weightingRotated = weightingRotated * mask
    
    if returnReducedComplex:
        from pytom_fftplan import fftShift
        from pytom_volume import fullToReduced
        fftShift(weightingRotated,True)
        returnVolume = fullToReduced(weightingRotated)
    else:
        returnVolume = weightingRotated

    if binarize:
        limit(returnVolume,0.5,0,0.5,1,True,True)
    
    return returnVolume


    
    
def wedgeFilter(volume,angle,radius=0,angleIsHalf=True,fourierOnly=False):
    """
    wedgeFilter: performs a single wedge filtering of volume.
    @param volume: The volume filtered.
    @param angle: Opening angle of the wedge.
    @param radius: cutoff radius (means lowpass filter radius - 0 by default)
    @param angleIsHalf: True if angle represents the whole wedge. (True by default)
    @return: The wedge filtered volume,the filter and the filtered volume in fourier space  
    @author: Thomas Hrabe
    """
    import pytom_freqweight
    
    if angle <= 0:
        return volume
    
    if not angleIsHalf:
        angle = abs(angle/2) 
    
    wf = pytom_freqweight.weight(angle,radius,volume.sizeX(),volume.sizeY(),volume.sizeZ())
    
    
    return filter(volume,wf,fourierOnly)
    
    
def bandpassFilter(volume, lowestFrequency, highestFrequency, bpf=None,
                   smooth=0,fourierOnly=False):
    """
    bandpassFilter: Performs bandpass filtering of volume.
    @param volume: The volume to be filtered
    @type volume: L{pytom_volume.vol}
    @param lowestFrequency: The lowest frequency in the filter as absolute pixel value
    @type lowestFrequency: L{float}
    @param highestFrequency: The highest frequency in the filter as absolute pixel value
    @type highestFrequency: L{float} 
    @param bpf: Bandpassfilter C++ object if it is already existing 
    @param smooth: Adjusts the size of gaussian falloff around each band end.
    @return: The bandpass filtered volume, the filter and the filtered volume in fourier space
    @author: Thomas Hrabe
    """

    import pytom_freqweight
    import pytom_volume
    
    if volume.__class__ == pytom_volume.vol:
        from pytom.basic.fourier import fft
        fvolume = fft(volume)
    else:
        fvolume = volume
 
    if not bpf:
        bpf = pytom_freqweight.weight(lowestFrequency,highestFrequency,
                                      fvolume.sizeX(),fvolume.sizeY(),fvolume.sizeZ(),smooth)    

    return filter(fvolume,bpf,fourierOnly)      


def lowpassFilter(volume,band,smooth=0,fourierOnly=False):
    """
    lowpassFilter: 
    @param volume: The volume filtered
    @param band: Upper end of band(in pixels)
    @param smooth: Adjusts the size of gaussian falloff around each band end.
    @return: The lowpass filtered volume,the filter and the filtered volume in fourier space
    @author: Thomas Hrabe 
    """
        
    import pytom_freqweight
    import pytom_volume
    
    if volume.__class__ == pytom_volume.vol:
        from pytom.basic.fourier import fft
        fvolume = fft(volume)
    else:
        fvolume = volume
        
    lpf = pytom_freqweight.weight(0.0,band,fvolume.sizeX(),fvolume.sizeY(),fvolume.sizeZ(),smooth)
    
    return filter(fvolume,lpf,fourierOnly)
    
def highpassFilter(volume,band,smooth=0,fourierOnly=False):   
    """
    highpassFilter:
    @param volume:  The volume filtered
    @param band: Lower end of band (in pixels)
    @param smooth: Adjusts the size of gaussian falloff around each band end.
    @return: The highpass filtered volume,the filter and the filtered volume in fourier space
    @author: Thomas Hrabe 
    """
    import pytom_freqweight
    import pytom_volume
    
    if volume.__class__ == pytom_volume.vol:
        from pytom.basic.fourier import fft
        fvolume = fft(volume)
    else:
        fvolume = volume
    
    hpf = pytom_freqweight.weight(band,fvolume.sizeX(),fvolume.sizeX(),fvolume.sizeY(),fvolume.sizeZ(),smooth)
       
    return filter(fvolume,hpf,fourierOnly)

def filter(volume,filterObject,fourierOnly=False):
    """
    filter: A generic filter method.
    @param volume: The volume to be filtered
    @param filterObject: A filter object (either wedgeFilter, bandpassFilter, ...)
    @return: The filtered volume,the filter and the filtered volume in fourier space
    @author: Thomas Hrabe
    """
        
    from pytom.basic.fourier import fft, ifft
    import pytom_volume
    
    if volume.__class__ == pytom_volume.vol:
        fvolume = fft(volume)
        noPixels = volume.sizeX() * volume.sizeY() * volume.sizeZ()
    else:
        fvolume = pytom_volume.vol_comp(volume.sizeX(),volume.sizeY(),volume.sizeZ())
        fvolume.copyVolume(volume)
        if (volume.sizeZ() == 1) and (volume.sizeY() == 1):
            noPixels = 2 *(volume.sizeX()-1)
        elif volume.sizeZ() == 1:
            noPixels = volume.sizeX() *2 *(volume.sizeY()-1)
        else:
            noPixels = volume.sizeX() * volume.sizeY() * 2*(volume.sizeZ()-1)
    filterObject.apply(fvolume)
    
    if not fourierOnly:
        result = ifft(fvolume)/noPixels
        return [result,filterObject,fvolume]   
    else:
        return [fvolume,filterObject]


def gaussian_filter(vol, sigma):
#    # construct the Gaussian kernel
#    import numpy as np
#    from scipy import mgrid, exp
#    sizex = vol.sizeX()
#    sizey = vol.sizeY()
#    sizez = vol.sizeZ()
#    
#    [x, y, z] = mgrid[-(sizex/2):(sizex+1)/2, -(sizey/2):(sizey+1)/2, -(sizez/2):(sizez+1)/2]
#    g = exp(-(x**2+y**2+z**2)/(2*float(sigma)**2))
#    g /= g.sum()
#    
#    # transfer to vol obj and do the filtering
#    from pytom_numpy import npy2vol
#    kernel = npy2vol(np.array(g, dtype='float32', order='F'), 3)
#    
#    from pytom.basic.fourier import convolute
#    res = convolute(vol, kernel, False)

    import numpy as np
    from pytom_numpy import vol2npy, npy2vol
    from scipy.ndimage.filters import gaussian_filter
    v = vol2npy(vol)
    r = gaussian_filter(v, sigma)
    res = npy2vol(np.array(r, dtype='float32', order='F'), 3)
    
    return res



