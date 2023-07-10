import numpy as np
from pytom.basic.fourier import radialaverage, ftshift, iftshift
from pytom.lib.pytom_volume import reducedToFull, fullToReduced, vol, max
from math import sqrt


def profile2FourierVol( profile, dim=None, reduced=False):
    """
    create Volume from 1d radial profile, e.g., to modulate signal with \
    specific function such as CTF or FSC. Simple linear interpolation is used\
    for sampling.

    @param profile: profile
    @type profile: 1-d L{pytom.lib.pytom_volume.vol} or 1-d python array
    @param dim: dimension of (cubic) output
    @type dim: L{int}
    @param reduced: If true reduced Fourier representation (N/2+1, N, N) is generated.
    @type reduced: L{bool}

    @return: 3-dim complex volume with spherically symmetrical profile
    @rtype: L{pytom.lib.pytom_volume.vol}
    @author: FF
    """
    from pytom.lib.pytom_volume import vol
    from math import sqrt

    if not dim:
        if type(profile) == vol:
            dim = 2*profile.size_x()
        else:
            dim = 2*len(profile)

    if reduced:
        nx = int(dim/2)+1
    else:
        nx = dim

    fkernel = vol( nx, dim, dim)
    if type(profile) == vol:
        r_max = profile.size_x()-1
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
    @type volume: L{pytom.lib.pytom_volume.vol}
    @param profile: 1-d profile
    @type profile: L{pytom.lib.pytom_volume.vol}
    @return: outvol
    @rtype: L{pytom.lib.pytom_volume.vol}
    @author: FF
    """
    from pytom.basic.filter import profile2FourierVol
    from pytom.basic.fourier import convolute, powerspectrum

    kernel = profile2FourierVol( profile=profile, dim=volume.size_x(), reduced=False)
    outvol = convolute(v=volume, k=kernel, kernel_in_fourier=True)
    return outvol

def gridCTF(x_array, y_array, z_array):
    """
    function similar to ndgrid in matlab

    @type x_array: 1-dim array
    @type y_array: 1-dim array
    @type z_array: 1-dim array
    @return: 3-dim volumes with x, y, z values in x,y,z dimension, respectively
    @rtype: L{pytom.lib.pytom_volume.vol}
    """
    from pytom.lib.pytom_volume import vol
    
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
    @return: 3-dim volumes with x, y, z values in x,y,z dimension, respectively
    @rtype: L{pytom.lib.pytom_volume.vol}
    """
    from pytom.lib.pytom_volume import vol, power
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
        
    for i in range(vol.size_x()):
        for j in range(vol.size_y()):
            for k in range(vol.size_z()):
                vol(math.sin(vol(i, j, k)) , i, j, k)    
    
    return vol

def convolutionCTF(volume, defocus, pixelSize=None, voltage=None, Cs=None, sigma=None):
    """
    convolutionCTF:
    @param volume: input volume to be convolved with CTF
    @type volume: L{pytom.lib.pytom_volume.vol}
    @param defocus: Negative value = underfocus (in microns)
    @param pixelSize: Size of pixels / voxels (in Anstroms)
    @param voltage: 
    @param Cs:
    @param sigma: 
    @return: CTF filtered volume
    @rtype L{pytom.lib.pytom_volume.vol}
    @author: FF
    """
    from pytom.lib.pytom_volume import subvolume, complexRealMult
    from pytom.basic.fourier import ftshift, fft, ifft
    from pytom.basic.filter import volCTF
    
    dimX = volume.size_x()
    dimY = volume.size_y()
    dimZ = volume.size_z()

    ctf = volCTF(defocus, dimX, dimY, dimZ, pixelSize, voltage, Cs, sigma)
    filterCTF = subvolume(ftshift(ctf,inplace=False), 0, 0, 0, dimX, dimY, (dimZ/2)+1)
    
    filteredVolume = ifft( complexRealMult(fft(volume), filterCTF) )
    
    return filteredVolume

def fourierFilterShift(filter):
    """
    fourierFilterShift: NEEDS Documentation
    @param filter: NEEDS Documentation 
    """
    from pytom.lib.pytom_volume import vol
    
    widthX = filter.size_x()
    centerX = filter.size_x()//2
    boxX = filter.size_x()//2
    
    widthY = filter.size_y()
    centerY = filter.size_y()//2
    boxY = filter.size_y()//2
    
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
    from pytom.lib.pytom_volume import vol

    widthX = filter.size_x()
    centerX = filter.size_x()//2
    boxX = filter.size_x()//2

    widthY = filter.size_y()
    centerY = filter.size_y()//2
    boxY = filter.size_y()//2

    shifted_filter = vol(widthX, widthY, 1)
    shifted_filter.setAll(0.0)

    for i in range(widthX):
        rx = (boxX-i)%widthX

        for j in range(widthY):
            ry = (boxY-j)%widthY

            shifted_filter.setV(filter.getV(i, j, 0), rx, widthY-j-1, 0)

    return shifted_filter

def circleFilter(size_x,size_y, radiusCutoff):
    """
    circleFilter: NEEDS Documentation
    @param size_x: NEEDS Documentation 
    @param size_y: NEEDS Documentation
    @param radiusCutoff: NEEDS Documentation
    """
    from pytom.lib.pytom_volume import vol
    
    centerX = size_x//2
        
    centerY = size_y//2
    size_y = (size_y//2) +1
        
    filter_vol = vol(size_x, size_y, 1)
    filter_vol.setAll(0.0)
    
    for i in range(size_x):        
        for j in range(size_y):
                        
            radius = ((i-centerX)*(i-centerX)+(j-centerY)*(j-centerY))**0.5
                                                
            if radius <= radiusCutoff:
                filter_vol.setV(1.0, i, j, 0)
                
    return filter_vol

def rampFilter( size_x, size_y, crowtherFreq=None, N=None):
    """
    rampFilter: Generates the weighting function required for weighted backprojection - y-axis is tilt axis

    @param size_x: size of weighted image in X
    @param size_y: size of weighted image in Y

    @return: filter volume
    
    """
    from pytom.lib.pytom_volume import vol
    from math import exp

    centerX = size_x//2
    
    centerY = size_y//2
    size_y = (size_y//2) +1

    N = 0 if N is None else 1 / N

    if crowtherFreq is None:
        Ny = size_x//2
    else:
        Ny = crowtherFreq

    filter_vol = vol(size_x, size_y, 1)
    filter_vol.setAll(0.0)

    # TODO this is very slow...
    for i in range(size_x):        
        distX = abs(float(i-centerX))
        ratio = min(1, (distX/Ny)+N)
        for j in range(size_y):            
            filter_vol.setV(ratio, i, j, 0)

    return filter_vol

def exactFilter(tilt_angles, tiltAngle, sX, sY, sliceWidth, arr=[]):
    """
    exactFilter: Generates the exact weighting function required for weighted backprojection - y-axis is tilt axis
    Reference : Optik, Exact filters for general geometry three dimensional reconstuction, vol.73,146,1986.
    @param tilt_angles: list of all the tilt angles in one tilt series
    @param titlAngle: tilt angle for which the exact weighting function is calculated
    @param size_x: size of weighted image in X
    @param size_y: size of weighted image in Y

    @return: filter volume

    """
    # Using Friedel Symmetry in Fourier space.
    sY = sY//2+1

    # Calculate the relative angles in radians.
    sampling = np.sin(np.abs((np.array(tilt_angles) - tiltAngle) * np.pi / 180.))
    smallest_sampling = np.min(sampling[sampling > 0.001])

    if sliceWidth / smallest_sampling > sX // 2:  # crowther crit can be to nyquist freq (i.e. sX // 2)
        sliceWidth = smallest_sampling * (sX // 2)  # adjust sliceWidth if too large

    crowtherFreq = min(sX // 2, int(np.ceil(sliceWidth / smallest_sampling)))
    arrCrowther = np.abs(np.arange(-crowtherFreq, min(sX // 2, crowtherFreq + 1)))

    # as in the paper: 1 - frequency / overlap_frequency
    # where frequency = arrCrowther, and 1 / overlap_frequency = sampling/sliceWidth
    wfuncCrowther = 1. / np.clip(1 - ((sampling / sliceWidth)[:, np.newaxis] * arrCrowther) ** 2, 0, 2).sum(axis=0)

    # Create full with weightFunc
    wfunc = np.ones((sX, sY, 1), dtype=float)
    wfunc[sX//2-crowtherFreq:sX//2+min(sX//2,crowtherFreq+1),:, 0] = np.tile(wfuncCrowther, (sY, 1)).T

    weightFunc = vol(sX, sY, 1)
    weightFunc.setAll(0.0)

    # use npy2vol instead??
    # weightFunc = npy2vol(array(wfunc, dtype='float32', order='F'), 3)
    for ix in range(0, sX):
        for iy in range(0, sY):
            weightFunc.setV(float(wfunc[ix][iy]), ix, iy, 0)
        
    return weightFunc

def rotateFilter(tilt_angles, tiltAngle, sX, sY, sliceWidth, arr=[]):
    """
    rotate filter function
    @param tilt_angles: ...
    @return: filter volume
    """
    from numpy import zeros_like, ones, column_stack, sin, abs, zeros, pi, ceil, floor, float32, array
    from scipy.ndimage import rotate
    from pytom.lib.pytom_volume import vol
    from pytom.basic.files import read
    from pytom.lib.pytom_numpy import npy2vol

    tilt_angles = array(tilt_angles)

    smallest = abs(tilt_angles - tiltAngle)
    smallest[smallest.argmin()] = 1000
    smallest = smallest.min()

    wFunc = ones((sX))
    size = min(sX, int(2 / sin(smallest * pi / 180.)))

    a = zeros((size, size))

    sY = sY // 2 + 1

    dame = zeros_like(a)
    dame[size // 2, :] = 1.

    out = zeros_like(a)

    for angle in (tilt_angles - tiltAngle):

        if abs(angle) > 0.0001:
            dame2 = rotate(dame, angle, axes=(0, 1), reshape=False)
        else:
            dame2 = dame.copy()
        out += dame2

    wFunc[sX // 2 - int(floor(size // 2)):sX // 2 + int(ceil(size / 2))] = 1 / out[size // 2, :]

    wfunc = column_stack( ([(wFunc), ] * (sY)) ).astype(float32)

    weightFunc = vol(sX, sY, 1)
    weightFunc.setAll(0.0)

    for ix in range(0, sX):
        for iy in range(0, sY):
            # print(ix,iy)
            weightFunc.setV(float(wfunc[ix][iy]), ix, iy, 0)

    return weightFunc

def rotateWeighting(weighting, z1, z2, x, mask=None, isReducedComplex=None, returnReducedComplex=False, binarize=False):
    """
    rotateWeighting: Rotates a frequency weighting volume around the center. If the volume provided is reduced complex, it will be rescaled to full size, ftshifted, rotated, iftshifted and scaled back to reduced size.
    @param weighting: A weighting volume
    @type weighting: L{pytom.lib.pytom_volume.vol}
    @param z1: Z1 rotation angle
    @type z1: float
    @param z2: Z2 rotation angle
    @type z2: float
    @param x: X rotation angle
    @type x: float
    @param mask:=None is there a rotation mask? A mask with all = 1 will be generated otherwise. Such mask should be \
        provided anyway.
    @type mask: L{pytom.lib.pytom_volume.vol}
    @param isReducedComplex: Either set to True or False. Will be determined otherwise
    @type isReducedComplex: bool
    @param returnReducedComplex: Return as reduced complex? (Default is False)
    @type returnReducedComplex: bool
    @param binarize: binarize weighting
    @type binarize: bool
    @return: weight as reduced complex volume
    @rtype: L{pytom.lib.pytom_volume.vol_comp}
    """
    from pytom.lib.pytom_volume import vol, limit, vol_comp, rotate
    assert type(weighting) == vol or  type(weighting) == vol_comp, "rotateWeighting: input neither vol nor vol_comp"
    
    isReducedComplex = isReducedComplex or int(weighting.size_x()/2)+1 == weighting.size_z();

    if isReducedComplex:
        #scale weighting to full size
        from pytom.lib.pytom_fftplan import fftShift
        from pytom.lib.pytom_volume import reducedToFull
        weighting = reducedToFull(weighting)
        fftShift(weighting, True)

    if not mask:
        mask = vol(weighting.size_x(),weighting.size_y(),weighting.size_z())
        mask.setAll(1)

    weightingRotated = vol(weighting.size_x(),weighting.size_y(),weighting.size_z())

    rotate(weighting,weightingRotated,z1,z2,x)
    weightingRotated = weightingRotated * mask
    
    if returnReducedComplex:
        from pytom.lib.pytom_fftplan import fftShift
        from pytom.lib.pytom_volume import fullToReduced
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
    @type volume: L{pytom.lib.pytom_volume.vol}
    @param angle: Opening angle of the wedge.
    @type angle: C{float}
    @param radius: cutoff radius (means lowpass filter radius - 0 by default)
    @type radius: C{Int}
    @param angleIsHalf: True if angle represents the whole wedge. (True by default)
    @return: The wedge filtered volume,the filter and the filtered volume in fourier space  
    @author: Thomas Hrabe
    """
    import pytom.lib.pytom_freqweight as pytom_freqweight
    
    if angle <= 0:
        return volume
    
    if not angleIsHalf:
        angle = abs(angle/2) 
    
    wf = pytom_freqweight.weight(angle,radius,volume.size_x(),volume.size_y(),volume.size_z())
    
    
    return filter(volume,wf,fourierOnly)

def bandpassFilter(volume, lowestFrequency, highestFrequency, bpf=None, smooth=0,fourierOnly=False):
    """
    bandpassFilter: Performs bandpass filtering of volume.
    @param volume: The volume to be filtered
    @type volume: L{pytom.lib.pytom_volume.vol}
    @param lowestFrequency: The lowest frequency in the filter as absolute pixel value
    @type lowestFrequency: L{float}
    @param highestFrequency: The highest frequency in the filter as absolute pixel value
    @type highestFrequency: L{float} 
    @param bpf: Bandpassfilter C++ object if it is already existing 
    @param smooth: Adjusts the size of gaussian falloff around each band end.
    @return: The bandpass filtered volume, the filter and the filtered volume in fourier space
    @author: Thomas Hrabe
    """

    import pytom.lib.pytom_freqweight as pytom_freqweight
    import pytom.lib.pytom_volume as pytom_volume
    
    if volume.__class__ == pytom_volume.vol:
        from pytom.basic.fourier import fft
        fvolume = fft(volume)
    else:
        fvolume = volume
 
    if not bpf:
        bpf = pytom_freqweight.weight(lowestFrequency,highestFrequency,
                                      fvolume.size_x(),fvolume.size_y(),fvolume.size_z(),smooth)    

    return filter(fvolume,bpf,fourierOnly)      

def lowpassFilter(volume,band,smooth=0,fourierOnly=False):
    """
    lowpassFilter: 
    @param volume: The volume filtered
    @type volume: L{pytom.lib.pytom_volume.vol}
    @param band: Upper end of band(in pixels)
    @param smooth: Adjusts the size of gaussian falloff around each band end.
    @return: The lowpass filtered volume,the filter and the filtered volume in fourier space
    @author: Thomas Hrabe 
    """
        
    import pytom.lib.pytom_freqweight as pytom_freqweight
    import pytom.lib.pytom_volume as pytom_volume
    
    if volume.__class__ == pytom_volume.vol:
        from pytom.basic.fourier import fft
        fvolume = fft(volume)
    else:
        fvolume = volume
        
    lpf = pytom_freqweight.weight(0.0,band,fvolume.size_x(),fvolume.size_y(),fvolume.size_z(),smooth)
    
    return filter(fvolume,lpf,fourierOnly)
    
def highpassFilter(volume,band,smooth=0,fourierOnly=False):   
    """
    highpassFilter:
    @param volume:  The volume filtered
    @type volume: L{pytom.lib.pytom_volume.vol}
    @param band: Lower end of band (in pixels)
    @param smooth: Adjusts the size of gaussian falloff around each band end.
    @return: The highpass filtered volume,the filter and the filtered volume in fourier space
    @author: Thomas Hrabe 
    """
    import pytom.lib.pytom_freqweight as pytom_freqweight
    import pytom.lib.pytom_volume as pytom_volume
    
    if volume.__class__ == pytom_volume.vol:
        from pytom.basic.fourier import fft
        fvolume = fft(volume)
    else:
        fvolume = volume
    
    hpf = pytom_freqweight.weight(band,fvolume.size_x(),fvolume.size_x(),fvolume.size_y(),fvolume.size_z(),smooth)
       
    return filter(fvolume,hpf,fourierOnly)

def filter(volume,filterObject,fourierOnly=False):
    """
    filter: A generic filter method.
    @param volume: The volume to be filtered
    @type volume: L{pytom.lib.pytom_volume.vol}
    @type volume: L{pytom.lib.pytom_volume.vol}
    @param filterObject: A filter object (either wedgeFilter, bandpassFilter, ...)
    @return: The filtered volume,the filter and the filtered volume in fourier space
    @author: Thomas Hrabe
    """
    from pytom.basic.fourier import fft, ifft
    import pytom.lib.pytom_volume as pytom_volume
    
    if volume.__class__ == pytom_volume.vol:
        fvolume = fft(volume)
        noPixels = volume.size_x() * volume.size_y() * volume.size_z()
    else:
        fvolume = pytom_volume.vol_comp(volume.size_x(),volume.size_y(),volume.size_z())
        fvolume.copyVolume(volume)
        if (volume.size_z() == 1) and (volume.size_y() == 1):
            noPixels = 2 *(volume.size_x()-1)
        elif volume.size_z() == 1:
            noPixels = volume.size_x() *2 *(volume.size_y()-1)
        else:
            noPixels = volume.size_x() * volume.size_y() * 2*(volume.size_z()-1)
    filterObject.apply(fvolume)
    
    if not fourierOnly:
        result = ifft(fvolume)/noPixels
        return [result,filterObject,fvolume]   
    else:
        return [fvolume,filterObject]

def gaussian_filter(vol, sigma):
    """
    @param vol: input volume
    @type vol: L{pytom.lib.pytom_volume.vol}
    @param sigma: xxx
    @type sigma: xxx
    @return: resulting volume
    @rtype: L{pytom.lib.pytom_volume.vol}
    """
    import numpy as np
    from pytom.lib.pytom_numpy import vol2npy, npy2vol
    from scipy.ndimage.filters import gaussian_filter
    v = vol2npy(vol)
    r = gaussian_filter(v, sigma)
    res = npy2vol(np.array(r, dtype='float32', order='F'), 3)
    
    return res

def wiener_filter(fsc, even_ctf_sum, odd_ctf_sum, fudge=1):
    """
    Generate relion type wiener filters for the even and odd averages.

    @param fsc: fsc curve of even and odd averages
    @type fsc: L{list}
    @param even_ctf_sum: even ctf sum volume, normalized over the number of particles in the half set, reduced complex
    @type even_ctf_sum: L{pytom.lib.pytom_volume.vol}
    @param odd_ctf_sum: odd ctf sum volume, normalized over the number of particles in the half set, reduced complex
    @type odd_ctf_sum: L{pytom.lib.pytom_volume.vol}
    @param fudge: relion fudge factor (T)
    @type fudge: L[float}
    @return: returns two filters for the odd and even half map filter1, filter2
    @rtype: L{[pytom.lib.pytom_volume.vol, pytom.lib.pytom_volume.vol]}
    """
    assert even_ctf_sum.size_x() == even_ctf_sum.size_y() and even_ctf_sum.size_x() // 2 + 1 == even_ctf_sum.size_z(), \
        'even volume does not have correct dimensions'
    assert even_ctf_sum.size_x() == odd_ctf_sum.size_x() and even_ctf_sum.size_y() == odd_ctf_sum.size_y() and \
           even_ctf_sum.size_z() == odd_ctf_sum.size_z(), 'even and odd ctf volumes not matching'

    # get a radial profile of the summed wedges
    wedge_ctf_sum = even_ctf_sum + odd_ctf_sum
    radial = radialaverage(wedge_ctf_sum, isreduced=True)

    # get full complex version of the ctf volumes
    even_ctf_full = ftshift(reducedToFull(even_ctf_sum), inplace=False)
    odd_ctf_full = ftshift(reducedToFull(odd_ctf_sum), inplace=False)
    # sum_ctf_full = ftshift(reducedToFull(wedge_ctf_sum), inplace=False)
    sx, sy, sz = even_ctf_full.size_x(), even_ctf_full.size_y(), even_ctf_full.size_z()
    centerx, centery, centerz = sx // 2 + 1, sy // 2 + 1, sz // 2 + 1

    # find the first zero crossing of the ctf
    fsc = np.array(fsc)
    zero_crossings = np.where(np.diff(fsc > 0))[0] + 1
    if len(zero_crossings) != 0:
        rmax = zero_crossings[0]
    else:
        rmax = sx // 2 - 1
    fsc[rmax:] = 0  # set to 0 after first crossing to prevent artifacts

    # turn radial average profile and fsc profile into volumes
    fsc_vol = profile2FourierVol(fsc, dim=sx)
    radial_vol = profile2FourierVol(radial, dim=sx)

    even_filter = vol(sx, sy, sz)
    even_filter.setAll(0.)
    odd_filter = vol(sx, sy, sz)
    odd_filter.setAll(0.)
    # final_average_filter = vol(sx, sy, sz)
    # final_average_filter.setAll(0.)

    for x in range(sx):
        for y in range(sy):
            for z in range(sz):
                r = sqrt((x - centerx) ** 2 + (y - centery) ** 2 + (z - centerz) ** 2)
                if r > (rmax - 0.5):  # after this needs to be zero
                    continue
                elif r == 0:  # center should always be 1
                    even_filter.setV(0, x, y, z)
                    odd_filter.setV(0, x, y, z)
                else:
                    fsc_val = fsc_vol.getV(x, y, z)
                    snr = fsc_val / (1 - fsc_val)
                    tau2 = snr / radial_vol.getV(x, y, z)
                    even_filter.setV(1 / (even_ctf_full.getV(x, y, z) + 1 / (fudge * tau2)), x, y, z)
                    odd_filter.setV(1 / (odd_ctf_full.getV(x, y, z) + 1 / (fudge * tau2)), x, y, z)
                    # fsc_prime = sqrt(2 * fsc_val / (1 + fsc_val))
                    # snr = fsc_prime / (1 - fsc_prime)

    even_filter.setV(max(even_filter), centerx, centery, centerz)
    odd_filter.setV(max(odd_filter), centerx, centery, centerz)

    return fullToReduced(iftshift(even_filter, inplace=False)), fullToReduced(iftshift(odd_filter, inplace=False))





