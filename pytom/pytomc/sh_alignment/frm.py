'''
Created on Sep 16, 2011

@author: yuxiangchen
'''

import numpy as np
import swig_frm # the path of this swig module should be set correctly in $PYTHONPATH
from .vol2sf import vol2sf


def enlarge2(corr):
    """Enlarge a volume to make each dimension twice as large using 3rd order spline interpolation.

    Parameters
    ----------
    corr: 3D volume to enlarge.
          numpy.array

    Returns
    -------
    A new 3D volume.
    numpy.array
    """
    nx, ny, nz = corr.shape
    corr = corr.reshape((nx*ny*nz,))
    res = np.zeros((nx*2*ny*2*nz*2,), dtype='double')

    swig_frm.enlarge2(corr, nx, ny, nz, res)

    return res.reshape(nx*2, ny*2, nz*2)

def create_wedge_sf(start, end, b, valIn=1.0, valOut=0.0):
    """Create a single tilt wedge-shaped spherical function.

    Parameters
    ----------
    start: Starting angle (in degree) of acquistion, e.g. -60.
           Integer or float.

    end: Ending angle (in degree) of acquistion, e.g. 60.
         Integer or float.

    b: Bandwidth of the spherical function.
       Integer.

    valIn: Value inside the sampled region.
           Integer or float. Default is 1.

    valOut: Value outside the sampled region.
            Integer or float. Default is 0.

    Returns
    -------
    A spherical function.
    numpy.array
    """
    from math import pi, sin, cos, atan
    
    res = []
    angle_range = [-end, -start]
    
    for j in range(2*b):
        for k in range(2*b):
            the = pi*(2*j+1)/(4*b) # (0,pi)
            phi = pi*k/b # [0,2*pi)
            
            x = cos(phi)*sin(the)
            z = cos(the)
            
            if x == 0.0:
                if z > 0:
                    angle = 90
                else:
                    angle = -90
            else:
                angle = atan(z/x)*180/pi
            
            if angle >= angle_range[0] and angle <= angle_range[1]:
                res.append(valIn)
            else:
                res.append(valOut)
    
    return np.array(res)

def frm_corr_peaks(f, g, npeaks=10, norm=False):
    """Do the correlation between two spherical functions and return the top N correlation peaks and the corresponding angles.

    Parameters
    ----------
    f: Spherical function f.
       numpy.array

    g: Spherical function g.
       numpy.array

    npeaks: Number of peaks returned.
            Integer, default 10.

    norm: Normalize both functions or not.
          Boolean, default False.

    Returns
    -------
    [[score, [phi, psi, theta]], ...]
    """
    f = np.array(f, dtype='double')
    g = np.array(g, dtype='double')
    
    # normalize the two signals
    if norm is True and np.std(f) != 0 and np.std(g) !=0:
        f = (f - np.mean(f))/np.std(f)
        g = (g - np.mean(g))/np.std(g)
    
    peaks = np.zeros(4*npeaks)
    
    if swig_frm.frm(f, g, peaks) != 0:
        raise RuntimeError('Something is wrong during FRM!')
    
    res = []
    for i in range(npeaks):
        if peaks[i*4] != 0: # reorder the result angles in pytom convention (phi, psi, the)
            psi = peaks[i*4+1]
            the = peaks[i*4+2]
            phi = peaks[i*4+3]
            
#            # make the angle inside meaningful range
#            while phi < 0 or phi > 360:
#                phi -= 360*sign(phi)
#            while psi < 0 or psi > 360:
#                psi -= 360*sign(psi)
#            while the < 0 or the > 360:
#                the -= 360*sign(the)
            
            res.append([peaks[i*4], [phi, psi, the]])
    
    return res

def frm_corr(f, g):
    """Do the correlation between two spherical functions and return the correlation as 3d array.

    Parameters
    ----------
    f: Spherical function f.
       numpy.array

    g: Spherical function g.
       numpy.array

    Returns
    -------
    numpy.array
    """
    f = np.array(f, dtype='double')
    g = np.array(g, dtype='double')
    
    b = len(f)**0.5/2
    b = int(np.floor(b))
    if b <= 2:
        raise RuntimeError('Bandwidth too small: %d' % b)
    c = np.zeros(16*b**3, dtype='double')
    
    if swig_frm.frm_corr(f, g, c) != 0:
        raise RuntimeError('Something is wrong during FRM!')
    
    c = c[::2] # retrieve the real part only
    c = np.reshape(c, (2*b, 2*b, 2*b), order='C')
    
    return c

def frm_ncorr(f, g):
    """Obsolete.
    """
    f = f - np.mean(f) # the offset is corrected here!
    g = g - np.mean(g)
    
    c = frm_corr(f, g)/(frm_corr(f, f).max()*frm_corr(g, g).max())**0.5 # the scaling is corrected here!

#    m = numpy.ones(len(f))
#    c = frm_corr(f, g)/(frm_corr(f**2, m) * frm_corr(g**2, m))**0.5
    
    return c
    

def frm_fourier_corr(fr, fi, gr, gi, return_real=False):
    """The Fourier space version of frm_corr.

    Parameters
    ----------
    fr: The real part of spherical function f.
        numpy.array

    fi: The imaginary part of spherical function f.
        numpy.array

    gr: The real part of spherical function g.
        numpy.array

    gi: The imaginary part of spherical function g.
        numpy.array

    return_real: Return the real part of the result of not.
                 Boolean, default False.

    Returns
    -------
    The correlation function.
    numpy.array
    """
    fr = np.array(fr, dtype='double')
    fi = np.array(fi, dtype='double')
    gr = np.array(gr, dtype='double')
    gi = np.array(gi, dtype='double')
    
    b = len(fr)**0.5/2
    b = int(np.floor(b))
    if b <= 2:
        raise RuntimeError('Bandwidth too small: %d' % b)
    c = np.zeros(16*b**3, dtype='double')
    
    if swig_frm.frm_fourier_corr(fr, fi, gr, gi, c) != 0:
        raise RuntimeError('Something is wrong during FRM!')
    
    if return_real: # return the real part only
        c = c[::2]
        c = np.reshape(c, (2*b, 2*b, 2*b), order='C')
    else:
        c = np.ndarray(shape=(8*b**3,), buffer=c, dtype='complex')
        c = np.reshape(c, (2*b, 2*b, 2*b), order='C')
        
#        re = c[::2]
#        re = np.reshape(re, (2*b, 2*b, 2*b), order='C')
#        im = c[1::2]
#        im = np.reshape(im, (2*b, 2*b, 2*b), order='C')
#        c = re+im*1j
    
    return c

def frm_idx2angle(bw, i, j, k):
    """Transfer the index from correlation volume to the actual Euler angles in ZXZ convention (degree).
    Note the order of the returned Euler angle is [Phi, Psi, Theta], or [Z1, Z2, X] in Pytom format.

    Parameters
    ----------
    bw: Bandwidth of the spherical harmonics.
        Integer

    i: First index of the correlation volume.
       Integer

    j: Second index of the correlation volume.
       Integer

    k: Second index of the correlation volume.
       Integer

    Returns
    -------
    Euler angle in degrees.
    """
    from math import pi, floor
    
    phe=i*pi/bw
    phc=j*pi/bw
    ohm=k*pi/bw
    
    phi = pi-ohm; phi -= 2.0*pi*floor(phi/2.0/pi) # make it inside [0, 2*PI]
    the = pi-phc
    psi = -phe; psi -= 2.0*pi*floor(psi/2.0/pi)
    
    # return in pytom way, the symbols here are different, ignore the naming!
    return [psi*180.0/pi, phi*180.0/pi, the*180.0/pi]


def frm_angle2idx(bw, phi, psi, the):
    """Transfer the Euler angle (ZXZ, degree) to the index in the correlation volume.
    """
    from math import floor
    phi -= floor(phi/360.)*360.
    psi -= floor(psi/360.)*360.
    the -= floor(the/180.)*180.

    from math import pi

    phi = phi*pi/180
    psi = psi*pi/180
    the = the*pi/180

    i = int(round(-(phi*bw/pi)))
    j = int(round((pi-the)*bw/pi))
    k = int(round((pi-psi)*bw/pi))

    if i < 0:
        i += 2*bw
    if k < 0:
        k += 2*bw

    assert i > -1 and i < 2*bw and j > -1 and j < 2*bw and k > -1 and k < 2*bw

    return i, j, k


def frm_vol(v1, v2, b, radius=None):
    """Calculate the correlation function in terms of Euler angle of two 3D volumes.

    Parameters
    ----------
    v1: Volume Nr. 1
        pytom_volume.vol

    v2: Volume Nr. 2 / Reference
        pytom_volume.vol

    b: Bandwidth of spherical harmonics.
       Integer

    radius: Radius of the volume involved during the calculation.
            Integer, default is half of the volume size.

    Returns
    -------
    Correlation volume.
    numpy.array
    """
    if not radius:
        radius = v1.sizeX()/2
    
    res = np.zeros((2*b, 2*b, 2*b))
    for r in range(1, radius+1):
        corr = frm_corr(vol2sf(v1, r, b), vol2sf(v2, r, b))
        res += corr*(r**2) # should multiply by r**2
    
    return res

def frm_constrained_vol(v1, m1, v2, m2, b, radius=None):
    """Calculate the real space constrained correlation of two volumes and return the correlation volume.

    Parameters
    ----------
    v1: Volume Nr. 1
        pytom_volume.vol

    m1: Spherical mask of volume Nr. 1
        numpy.array

    v2: Volume Nr. 2 / Reference
        pytom_volume.vol

    m2: Spherical mask of v2
        numpy.array

    b: Bandwidth of spherical harmonics.
       Integer

    radius: Radius of the volume involved during the calculation.
            Integer, default is half of the volume size.

    Returns
    -------
    Correlation volume.
    numpy.array
    """
    if not radius:
        radius = v1.sizeX()/2
    
    res = np.zeros((2*b, 2*b, 2*b))
    for r in range(1, radius+1):
        corr = frm_constrained_corr(vol2sf(v1, r, b), m1, vol2sf(v2, r, b), m2)
        res += corr*(r**2) # should multiply by r**2
    
    return res


def frm_find_best_angle(corr, b=None):
    """Find the Euler angle corresponding to the peak of a correlation volume.

    Parameters
    ----------
    corr: Correlation volume.
          numpy.array

    b: Bandwidth of spherical harmonics.
       Integer. Automatically determined, if not specified.

    Returns
    -------
    Euler angle [Phi, Psi, Theta] in degrees.
    """
    if not b:
        b = corr.shape[0]/2
    
    idx = np.where(corr == corr.max())
    ang = frm_idx2angle(b, idx[0][0], idx[1][0], idx[2][0])
    return ang


def find_subpixel_peak_position(corr, pos):
    """Find the peak position in a correlation volume in sub-pixel accuracy around the given position using interpolation.

    Parameters
    ----------
    corr: Correlation volume.
          numpy.array

    pos: Position around which to search.
         List, [x, y, z].

    Returns
    -------
    The sub-pixel peak position [x', y', z'] and the interpolated score.
    """
    try:
        x = pos[0]
        y = pos[1]
        z = pos[2]
        
        dx = (corr[x+1][y][z] - corr[x-1][y][z])/2
        dy = (corr[x][y+1][z] - corr[x][y-1][z])/2
        dz = (corr[x][y][z+1] - corr[x][y][z-1])/2
        
        dxx = corr[x+1][y][z] + corr[x-1][y][z] - 2*corr[x][y][z]
        dyy = corr[x][y+1][z] + corr[x][y-1][z] - 2*corr[x][y][z]
        dzz = corr[x][y][z+1] + corr[x][y][z-1] - 2*corr[x][y][z]
        
        dxy = (corr[x+1][y+1][z] + corr[x-1][y-1][z] - corr[x+1][y-1][z] - corr[x-1][y+1][z])/4
        dxz = (corr[x+1][y][z+1] + corr[x-1][y][z-1] - corr[x+1][y][z-1] - corr[x-1][y][z+1])/4
        dyz = (corr[x][y+1][z+1] + corr[x][y-1][z-1] - corr[x][y-1][z+1] - corr[x][y+1][z-1])/4
        
        aa = np.array([[dxx, dxy, dxz], [dxy, dyy, dyz], [dxz, dyz, dzz]])
        bb = np.array([-dx, -dy, -dz])
        detx, dety, detz = np.linalg.solve(aa, bb)
        
        peak = corr[x][y][z] + dx*detx + dy*dety + dz*detz + dxx*detx**2/2 + dyy*dety**2/2 + dzz*detz**2/2 + detx*dety*dxy + detx*detz*dxz + dety*detz*dyz
    except:
        detx, dety, detz = 0, 0, 0
        peak = corr[x][y][z]

    if np.abs(detx) >= 1 or np.abs(dety) >= 1 or np.abs(detz) >= 1: # if enter this condition, something is not right
        detx, dety, detz = 0, 0, 0
        peak = corr[x][y][z]

    return [x+detx, y+dety, z+detz], peak

def frm_find_best_angle_interp(corr, b=None):
    """Interpolated version of frm_find_best_angle.

    Parameters
    ----------
    corr: Correlation volume.
          numpy.array

    b: Bandwidth of spherical harmonics.
       Integer. Automatically determined, if not specified.

    Returns
    -------
    Euler angle [Phi, Psi, Theta] in degrees.
    """
    if not b:
        b = corr.shape[0]/2
    
    idx = np.where(corr == corr.max())

    # interpolation for subpixel accuracy
    pos, peak = find_subpixel_peak_position(corr, [idx[0][0], idx[1][0], idx[2][0]])
    ang = frm_idx2angle(b, pos[0], pos[1], pos[2])
    
    return ang, peak


def frm_find_best_angluar_match(v1, v2, b, radius=None):
    """Calculate the best angular match of two volumes, assuming no translation and no missing wedge.

    Parameters
    ----------
    v1: Volume Nr. 1
        pytom_volume.vol

    v2: Volume Nr. 2 / Reference
        pytom_volume.vol

    b: Bandwidth of spherical harmonics.
       Integer

    radius: Radius of the volume involved during the calculation.
            Integer, default is half of the volume size.

    Returns
    -------
    Best Euler angle [Phi, Psi, Theta] to rotate v2 to match v1.
    """
    res = frm_vol(v1, v2, b, radius)
    
    return frm_find_best_angle(res, b)


def frm_constrained_corr(f, mf, g, mg, norm=False, return_score=True):
    """Calculate the correlation as a function of Euler angle between two spherical functions.

    Parameters
    ----------
    f: Spherical function Nr. 1
       numpy.array

    mf: Spherical mask of f
        numpy.array

    g: Spherical function Nr. 2
       numpy.array

    mg: Spherical mask of g
        numpy.array

    norm: Normalize the calculation or not.
          Boolean, default is False.

    return_score: Return the correlation score or return the intermediate result (numerator, denominator1, denominator2).
                  Boolean, default is True.

    Returns
    -------
    If return_score is set to True, return the correlation function; otherwise return the intermediate result.
    """
    f = np.array(f, dtype='double')
    mf = np.array(mf, dtype='double')
    g = np.array(g, dtype='double')
    mg = np.array(mg, dtype='double')
    
    if norm:
        dumb1 = frm_corr(f*mf, g*mg)
        dumb2 = frm_corr(mf, g*mg)
        dumb3 = frm_corr(f*mf, mg)
        dumb4 = frm_corr(mf, mg)
        dumb5 = frm_corr(f**2*mf, mg)
        dumb6 = frm_corr(mf, g**2*mg)
        
        fmean = dumb3 / dumb4
        gmean = dumb2 / dumb4
        
        numerator = dumb1 - fmean*dumb2 - gmean*dumb3 + fmean*gmean*dumb4
        denominator1 = dumb5 + fmean**2*dumb4 - 2*fmean*dumb3
        denominator2 = dumb6 + gmean**2*dumb4 - 2*gmean*dumb2
    else:
        numerator = frm_corr(f*mf, g*mg)
        denominator1 = frm_corr(f**2*mf, mg)
        denominator2 = frm_corr(mf, g**2*mg)
    
    if return_score:
        denominator = (np.abs(denominator1) * np.abs(denominator2))**0.5    
        res = numerator / denominator
    else:
        res = (numerator, np.abs(denominator1), np.abs(denominator2))
    
    return res

def frm_fourier_constrained_corr(fr, fi, mf, gr, gi, mg, return_real=False, norm=False, return_score=True):
    """The Fourier space version of frm_constrained_corr.

    Parameters
    ----------
    fr: Real part of spherical function Nr. 1
        numpy.array

    fi: Imaginary part of spherical function Nr. 1
        numpy.array

    mf: Spherical mask of f
        numpy.array

    gr: Real part of spherical function Nr. 2
        numpy.array

    gi: Imaginary part of spherical function Nr. 2
        numpy.array

    mg: Spherical mask of g
        numpy.array

    norm: Normalize the calculation or not.
          Boolean, default is False.

    return_score: Return the correlation score or return the intermediate result (numerator, denominator1, denominator2).
                  Boolean, default is True.

    Returns
    -------
    If return_score is set to True, return the correlation function; otherwise return the intermediate result.
    """
    if fr.__class__ != np.array:
        fr = np.array(fr, dtype='double')
    if fi.__class__ != np.array:
        fi = np.array(fi, dtype='double')
    if mf.__class__ != np.array:
        mf = np.array(mf, dtype='double')
    if gr.__class__ != np.array:
        gr = np.array(gr, dtype='double')
    if gi.__class__ != np.array:
        gi = np.array(gi, dtype='double')
    if mg.__class__ != np.array:
        mg = np.array(mg, dtype='double')
    
    mi = np.zeros(mf.shape) # the imag part of the masks, all zero
    
    if norm:
        dumb1 = frm_fourier_corr(fr*mf, fi*mf, gr*mg, gi*mg)
        dumb2 = frm_fourier_corr(mf, mi, gr*mg, gi*mg)
        dumb3 = frm_fourier_corr(fr*mf, fi*mf, mg, mi)
        dumb4 = frm_fourier_corr(mf, mi, mg, mi)
        dumb5 = frm_fourier_corr((fr**2+fi**2)*mf, mi, mg, mi) # be careful with this part! abs^2!
        dumb6 = frm_fourier_corr(mf, mi, (gr**2+gi**2)*mg, mi)
        
        fmean = dumb3 / dumb4
        gmean = dumb2 / dumb4
        
        numerator = dumb1 - fmean*dumb2 - gmean*dumb3 + fmean*gmean*dumb4
        denominator1 = dumb5 + fmean**2*dumb4 - 2*fmean*dumb3
        denominator2 = dumb6 + gmean**2*dumb4 - 2*gmean*dumb2
    else:
        numerator = frm_fourier_corr(fr*mf, fi*mf, gr*mg, gi*mg)
        denominator1 = frm_fourier_corr((fr**2+fi**2)*mf, mi, mg, mi)
        denominator2 = frm_fourier_corr(mf, mi, (gr**2+gi**2)*mg, mi)

    if return_score:
        denominator = (np.abs(denominator1) * np.abs(denominator2))**0.5
        # denominator = (denominator1 * denominator2)**0.5
        res = numerator / denominator
        
        if return_real:
            res = np.real(res)
    else:
        if return_real:
            res = (np.real(numerator), np.abs(denominator1), np.abs(denominator2))
        else:
            res = (numerator, np.abs(denominator1), np.abs(denominator2))
    
    return res


def np_transfer_idx(idx, shape):
    """Transfer the flattened numpy index into 3d index.
    """
    x = idx / (shape[1]*shape[2])
    y = (idx - x*shape[1]*shape[2]) / shape[2]
    z = idx - x*shape[1]*shape[2] - y*shape[2]
    return [x, y, z]


def frm_find_topn_angles_interp(corr, n=5, dist=3.0):
    """Find the Euler angles corresponding to the top N peaks in the correlation function using interpolation, given the minimal distance between them.
    
    Parameters
    ----------
    corr: Correlation function.
          numpy.array

    n: Number of peaks.
       Integer, default 5.

    dist: Euler angle distance threshold in bandwidth unit, e.g. if bandwidth is b, then the angle distance is dist/b*180
          Integer or float, default 3.

    Returns
    -------
    List: [(phi, psi, theta, peak_value), ...]
    """
    from pytom.tools.maths import rotation_distance

    b = corr.shape[0]/2
    
    res = []
    sort = corr.argsort(axis=None)
    
    for i in range(len(sort)-1, -1, -1):
        x,y,z = np_transfer_idx(sort[i], corr.shape)
        # ang = frm_idx2angle(b, x, y, z)
        pos, peak = find_subpixel_peak_position(corr, [x, y, z])
        ang = frm_idx2angle(b, pos[0], pos[1], pos[2])
        
        # eliminate the nearby solution
        for a in res:
            if rotation_distance(ang, a[0]) < float(dist)/b*180.0:
                break
        else:
            # res.append((ang, corr[x,y,z]))
            res.append((ang, peak))
        
        if len(res) >= n:
            break
    
    return res


def frm_find_topn_angles_interp2(corr, n=5, dist_cut=3.0):
    """2nd version of frm_find_topn_angles_interp, much faster!

    Parameters
    ----------
    corr: Correlation function.
          numpy.array

    n: Number of peaks.
       Integer, default 5.

    dist: Euler angle distance threshold in bandwidth unit, e.g. if bandwidth is b, then the angle distance is dist/b*180
          Integer or float, default 3.

    Returns
    -------
    List: [(phi, psi, theta, peak_value), ...]
    """
    if dist_cut < 1.0:
        dist_cut = 1.0

    nx, ny, nz = corr.shape
    corr = corr.reshape((nx*ny*nz,))
    b = nx//2

    peaks = np.zeros((n*4,), dtype='double')
    
    if swig_frm.find_topn_angles(corr, b, peaks, dist_cut) != 0:
        raise RuntimeError('Error happens in finding peaks!')
    
    peaks = peaks.reshape((n, 4))
    res = []
    for a in peaks:
        res.append(([a[1], a[3], a[2]], a[0])) # reorder it for return

    return res


def rt2tr(translation, rotation):
    """Change the transformation order from 1T2R to 1R2T.
    """
    from pytom.tools.maths import TransformationMatrix
    from pytom.angles.angleFnc import matToZXZ
    m = TransformationMatrix(rotation, [0,0,0])*TransformationMatrix([0,0,0], translation)
    return (matToZXZ(m.getRotationMatrix()).toList(), m.getTranslation())


def sph_correlate_ps(f, mf, g, mg, to_calculate=0):
    """Calculate the correlation of two (power spectrum) spherical functions.

    Parameters
    ----------
    f: Spherical function Nr. 1
       numpy.array

    mf: Spherical mask of f
        numpy.array

    g: Spherical function Nr. 2
       numpy.array

    mg: Spherical mask of g
        numpy.array

    to_calculate: 1 - calculate the numerator
                  2 - calculate the numerator and the denominator1
                  else - calculate all

    Returns
    -------
    Tuple: (numerator, denominator1, denominator2)
    """
    f = np.array(f, dtype='double')
    mf = np.array(mf, dtype='double')
    g = np.array(g, dtype='double')
    mg = np.array(mg, dtype='double')
    
    numerator = frm_corr(f*mf, g*mg)

    if to_calculate == 1: # calculate only the numerator
        return (numerator, None, None)
    elif to_calculate == 2: # calculate the numerator and denominator1
        denominator1 = frm_corr(f**2*mf, mg)
        return (numerator, np.abs(denominator1), None)
    else: # calculate all
        denominator1 = frm_corr(f**2*mf, mg)
        denominator2 = frm_corr(mf, g**2*mg)
        return (numerator, np.abs(denominator1), np.abs(denominator2))


def sph_correlate_fourier(fr, fi, mf, gr, gi, mg, to_calculate=0):
    """Calculate the correlation of two (Fourier space) spherical functions.

    Parameters
    ----------
    fr: Real part of spherical function Nr. 1
        numpy.array

    fi: Imaginary part of spherical function Nr. 1
        numpy.array

    mf: Spherical mask of f
        numpy.array

    gr: Real part of spherical function Nr. 2
        numpy.array

    gi: Imaginary part of spherical function Nr. 2
        numpy.array

    mg: Spherical mask of g
        numpy.array

    to_calculate: 1 - calculate the numerator
                  2 - calculate the numerator and the denominator1
                  else - calculate all

    Returns
    -------
    Tuple: (numerator, denominator1, denominator2)
    """
    if fr.__class__ != np.array:
        fr = np.array(fr, dtype='double')
    if fi.__class__ != np.array:
        fi = np.array(fi, dtype='double')
    if mf.__class__ != np.array:
        mf = np.array(mf, dtype='double')
    if gr.__class__ != np.array:
        gr = np.array(gr, dtype='double')
    if gi.__class__ != np.array:
        gi = np.array(gi, dtype='double')
    if mg.__class__ != np.array:
        mg = np.array(mg, dtype='double')
    
    mi = np.zeros(mf.shape) # the imag part of the masks, all zero
    
    numerator = frm_fourier_corr(fr*mf, fi*mf, gr*mg, gi*mg)

    if to_calculate == 1: # calculate only the numerator
        return (np.real(numerator), None, None)
    elif to_calculate == 2: # calculate the numerator and denominator1
        denominator1 = frm_fourier_corr((fr**2+fi**2)*mf, mi, mg, mi)
        return (np.real(numerator), np.abs(denominator1), None)
    else: # calculate all
        denominator1 = frm_fourier_corr((fr**2+fi**2)*mf, mi, mg, mi)
        denominator2 = frm_fourier_corr(mf, mi, (gr**2+gi**2)*mg, mi)
        return (np.real(numerator), np.abs(denominator1), np.abs(denominator2))


def get_adaptive_bw(k, b=None):
    """Get the bandwidth adaptively according to the radius/frequency.

    Parameters
    ----------
    k: Radius/frequency.
       Integer

    b: Bandwidth range of spherical harmonics.
       None -> [4, 64]
       List -> [b_min, b_max]
       Integer -> [b, b]

    Returns
    -------
    Integer.
    """
    from math import log, ceil, pow

    if not b: # set the bandwidth adaptively
        b_min = 4
        b_max = 64
    elif b.__class__ == tuple or b.__class__ == list:
        b_min = b[0]
        b_max = b[1]
    elif isinstance(b, int): # fixed bandwidth
        b_min = b
        b_max = b
    else:
        raise RuntimeError("Argument b is not valid!")

    bw = int(pow(2, int(ceil(log(2*k, 2)))))
    if bw < b_min:
        return b_min
    elif bw > b_max:
        return b_max
    else:
        return bw

def frm_correlate(vf, wf, vg, wg, b, max_freq, weights=None, ps=False, denominator1=None, denominator2=None, return_score=True):
    """Calculate the correlation of two volumes as a function of Euler angle.

    Parameters
    ----------
    vf: Volume Nr. 1
        pytom_volume.vol

    wf: Mask of vf in Fourier space.
        pytom.basic.structures.Wedge

    vg: Volume Nr. 2 / Reference
        pytom_volume.vol

    wg: Mask of vg in Fourier space.
        pytom.basic.structures.Wedge

    b: Bandwidth range of spherical harmonics.
       None -> [4, 64]
       List -> [b_min, b_max]
       Integer -> [b, b]

    max_freq: Maximal frequency involved in calculation.
              Integer.

    weights: Obsolete.

    ps: Calculation based on only the power spectrum of two volumes or not.
        Boolean. Default is False.

    denominator1: If the denominator1 is provided or not. If yes, do not have to re-calculate it again.
                  This field is used out of computation effeciency consideration.
                  Default is None, not provided.

    denominator2: If the denominator2 is provided or not. If yes, do not have to re-calculate it again.
                  This field is used out of computation effeciency consideration.
                  Default is None, not provided.

    return_score: Return the correlation score or return the intermediate result (numerator, denominator1, denominator2).
                  Boolean, default is True.

    Returns
    -------
    If return_score is set to True, return the correlation function; otherwise return the intermediate result.
    """
    if not weights: # weights, not used yet
        weights = [1 for i in range(max_freq)]

    from pytom.basic.fourier import fft, ifft, ftshift, iftshift
    from pytom_volume import vol, reducedToFull, abs, real, imag, rescale
    from .vol2sf import vol2sf
    from math import log, ceil, pow

    # IMPORTANT!!! Should firstly do the IFFTSHIFT on the volume data (NOT FFTSHIFT since for odd-sized data it matters!),
    # and then followed by the FFT.
    vf = ftshift(reducedToFull(fft(iftshift(vf, inplace=False))), inplace=False)
    vg = ftshift(reducedToFull(fft(iftshift(vg, inplace=False))), inplace=False)

    if ps: # power spectrum only
        ff = abs(vf)
        ff = real(ff)
        gg = abs(vg)
        gg = real(gg)
    else: # use spline intepolation on the real/imaginary parts. Can be done better, but now it suffices.
        vfr = real(vf)
        vfi = imag(vf)
        vgr = real(vg)
        vgi = imag(vg)

    numerator = None
    if denominator1 is not None and denominator2 is not None:
        to_calculate = 1
    elif denominator1 is None and denominator2 is not None:
        to_calculate = 2
    else:
        to_calculate = 0

    _last_bw = 0
    # might be a better idea to start from 2 due to the bad interpolation around 0 frequency!
    # this can be better solved by NFFT!
    for r in range(1, max_freq+1):
        # calculate the appropriate bw
        bw = get_adaptive_bw(r, b)

        # construct the wedge masks accordingly
        # now this part has been shifted to Pytom
        # if _last_bw != bw:
        #     # mf = create_wedge_sf(wf[0], wf[1], bw)
        #     # mg = create_wedge_sf(wg[0], wg[1], bw)
        #     mf = wf.toSphericalFunc(bw)
        #     mg = wg.toSphericalFunc(bw)
        mf = wf.toSphericalFunc(bw, r)
        mg = wg.toSphericalFunc(bw, r)

        if ps:
            corr1, corr2, corr3 = sph_correlate_ps(vol2sf(ff, r, bw), mf, vol2sf(gg, r, bw), mg, to_calculate)
        else:
            corr1, corr2, corr3 = sph_correlate_fourier(vol2sf(vfr, r, bw), vol2sf(vfi, r, bw), mf, vol2sf(vgr, r, bw), vol2sf(vgi, r, bw), mg, to_calculate)
        
        if _last_bw != bw: # size is different, have to do enlarge
            if numerator is None:
                numerator = np.zeros((2*bw, 2*bw, 2*bw), dtype='double')
                if to_calculate == 1:
                    pass
                elif to_calculate == 2:
                    denominator1 = np.zeros((2*bw, 2*bw, 2*bw), dtype='double')
                else:
                    denominator1 = np.zeros((2*bw, 2*bw, 2*bw), dtype='double')
                    denominator2 = np.zeros((2*bw, 2*bw, 2*bw), dtype='double')
            else:
                numerator = enlarge2(numerator)
                if to_calculate == 1:
                    pass
                elif to_calculate == 2:
                    denominator1 = enlarge2(denominator1)
                else:
                    denominator1 = enlarge2(denominator1)
                    denominator2 = enlarge2(denominator2)

        numerator += corr1*(r**2)*weights[r-1]
        if to_calculate == 1:
            pass
        elif to_calculate == 2:
            denominator1 += corr2*(r**2)*weights[r-1]
        else:
            denominator1 += corr2*(r**2)*weights[r-1]
            denominator2 += corr3*(r**2)*weights[r-1]

        _last_bw = bw

    if return_score:
        res = numerator/(denominator1 * denominator2)**0.5
        return res
    else:
        return (numerator, denominator1, denominator2)


def frm_align(vf, wf, vg, wg, b, max_freq, peak_offset=None, mask=None, weights=None, position=None, num_seeds=5):
    """Find the best alignment (translation & rotation) of volume vg (reference) to match vf.
    For details, please check the paper.

    Parameters
    ----------
    vf: Volume Nr. 1
        pytom_volume.vol

    wf: Mask of vf in Fourier space.
        pytom.basic.structures.Wedge. If none, no missing wedge.

    vg: Volume Nr. 2 / Reference
        pytom_volume.vol

    wg: Mask of vg in Fourier space.
        pytom.basic.structures.Wedge. If none, no missing wedge.

    b: Bandwidth range of spherical harmonics.
       None -> [4, 64]
       List -> [b_min, b_max]
       Integer -> [b, b]

    max_freq: Maximal frequency involved in calculation.
              Integer.

    peak_offset: The maximal offset which allows the peak of the score to be.
                 Or simply speaking, the maximal distance allowed to shift vg to match vf.
                 This parameter is needed to prevent shifting the reference volume out of the frame.
                 pytom_volume.vol / Integer. By default is half of the volume radius.

    mask: Mask volume for vg in real space.
          pytom_volume.vol

    weights: Obsolete.

    position: If the translation is already known or not. If provided, no translational search will be conducted.
              List: [x,y,z], default None.

    num_seeds: Number of threads for the expectation maximization procedure. The more the better, yet slower.
               Integer, default is 5.

    Returns
    -------
    (The best translation and rotation (Euler angle, ZXZ convention [Phi, Psi, Theta]) to transform vg to match vf.
    (best_translation, best_rotation, correlation_score)
    """
    from pytom_volume import vol, rotateSpline, peak
    from pytom.basic.transformations import shift
    from pytom.basic.correlation import FLCF
    from pytom.basic.filter import lowpassFilter
    from pytom.basic.structures import Mask, SingleTiltWedge, Rotation
    from pytom_volume import initSphere
    from pytom_numpy import vol2npy

    if vf.sizeX()!=vg.sizeX() or vf.sizeY()!=vg.sizeY() or vf.sizeZ()!=vg.sizeZ():
        raise RuntimeError('Two volumes must have the same size!')

    if wf is None:
        wf = SingleTiltWedge(0)
    else: # apply wedge respectively
        vf = wf.apply(vf)
    if wg is None:
        wg = SingleTiltWedge(0)
    else: # apply wedge respectively
        vg = wg.apply(vg)

    if peak_offset is None:
        peak_offset = vol(vf.sizeX(), vf.sizeY(), vf.sizeZ())
        initSphere(peak_offset, vf.sizeX()/4, 0,0, vf.sizeX()/2,vf.sizeY()/2,vf.sizeZ()/2)
    elif isinstance(peak_offset, int):
        peak_radius = peak_offset
        peak_offset = vol(vf.sizeX(), vf.sizeY(), vf.sizeZ())
        initSphere(peak_offset, peak_radius, 0,0, vf.sizeX()/2,vf.sizeY()/2,vf.sizeZ()/2)
    elif peak_offset.__class__ == vol:
        pass
    else:
        raise RuntimeError('Peak offset is given wrong!')

    # cut out the outer part which normally contains nonsense
    m = vol(vf.sizeX(), vf.sizeY(), vf.sizeZ())
    initSphere(m, vf.sizeX()/2, 0,0, vf.sizeX()/2,vf.sizeY()/2,vf.sizeZ()/2)
    vf = vf*m
    vg = vg*m
    if mask is None:
        mask = m
    else:
        vg = vg*mask

    if position is None: # if position is not given, we have to find it ourself
        # first roughtly determine the orientation (only according to the energy info)
        # get multiple candidate orientations
        numerator, denominator1, denominator2 = frm_correlate(vf, wf, vg, wg, b, max_freq, weights, True, None, None, False)
        score = numerator/(denominator1 * denominator2)**0.5
        res = frm_find_topn_angles_interp2(score, num_seeds, get_adaptive_bw(max_freq, b)/16.)
    else:
        # the position is given by the user
        vf2 = shift(vf, -position[0]+vf.sizeX()/2, -position[1]+vf.sizeY()/2, -position[2]+vf.sizeZ()/2, 'fourier')
        score = frm_correlate(vf2, wf, vg, wg, b, max_freq, weights, ps=False)
        orientation, max_value = frm_find_best_angle_interp(score)

        return position, orientation, max_value

    # iteratively refine the position & orientation
    from pytom.tools.maths import euclidianDistance
    max_iter = 10 # maximal number of iterations
    mask2 = vol(mask.sizeX(), mask.sizeY(), mask.sizeZ()) # store the rotated mask
    vg2 = vol(vg.sizeX(), vg.sizeY(), vg.sizeZ())
    lowpass_vf = lowpassFilter(vf, max_freq, max_freq/10.)[0]
    # apply vg's wedge to vf
    lowpass_vf = wg.apply(lowpass_vf)

    max_position = None
    max_orientation = None
    max_value = -1.0
    for i in range(num_seeds):
        old_pos = [-1, -1, -1]
        lm_pos = [-1, -1, -1]
        lm_ang = None
        lm_value = -1.0
        orientation = res[i][0] # initial orientation
        for j in range(max_iter):
            rotateSpline(vg, vg2, orientation[0], orientation[1], orientation[2]) # first rotate
            rotateSpline(mask, mask2, orientation[0], orientation[1], orientation[2]) # rotate the mask as well
            vg2 = wf.apply(vg2) # then apply the wedge
            vg2 = lowpassFilter(vg2, max_freq, max_freq/10.)[0]
            vf2 = wg.apply(lowpass_vf, Rotation(orientation)) # apply vg's wedge to vf with rotation!
            score = FLCF(vf2, vg2, mask2) # find the position
            pos = peak(score, peak_offset)
            pos, val = find_subpixel_peak_position(vol2npy(score), pos)
            if val > lm_value:
                lm_pos = pos
                lm_ang = orientation
                lm_value = val
        
            if euclidianDistance(lm_pos, old_pos) <= 1.0:
                # terminate this thread
                if lm_value > max_value:
                    max_position = lm_pos
                    max_orientation = lm_ang
                    max_value = lm_value

                break
            else:
                old_pos = lm_pos

            # here we shift the target, not the reference
            # if you dont want the shift to change the energy landscape, use fourier shift
            vf2 = shift(vf, -lm_pos[0]+vf.sizeX()/2, -lm_pos[1]+vf.sizeY()/2, -lm_pos[2]+vf.sizeZ()/2, 'fourier')
            score = frm_correlate(vf2, wf, vg, wg, b, max_freq, weights, False, denominator1, denominator2, True)
            # orientation, val = frm_find_topn_angles_interp(score, 1, 3)[0]
            orientation, val = frm_find_best_angle_interp(score)
            
        else: # no converge after the specified iteration, still get the best result as we can
            if lm_value > max_value:
                max_position = lm_pos
                max_orientation = lm_ang
                max_value = lm_value

        # print max_value # for show the convergence of the algorithm

    return max_position, max_orientation, max_value

