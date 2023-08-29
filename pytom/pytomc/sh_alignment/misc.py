from .frm import *

def frm_determine_orientation_rscore(vf, wf, vg, wg, b, radius=None, weights=None):
    """Obsolete.
    """
    if not radius: # set the radius
        radius = vf.size_x()/2
    if not weights: # set the weights
        weights = [1 for i in range(radius)]
    
    if not b: # set the bandwidth adaptively
        b_min = 4
        b_max = 128
    elif b.__class__ == tuple or b.__class__ == list:
        b_min = b[0]
        b_max = b[1]
    elif isinstance(b, int): # fixed bandwidth
        b_min = b
        b_max = b
    else:
        raise RuntimeError("Argument b is not valid!")
    
    from pytom.basic.fourier import fft, ifft, ftshift, iftshift
    from pytom.lib.pytom_volume import vol, reducedToFull, rescale, abs, real
    from .vol2sf import vol2sf
    from pytom.lib.pytom_numpy import vol2npy
    from math import log, ceil, pow

    # IMPORTANT!!! Should firstly do the IFFTSHIFT on the volume data (NOT FFTSHIFT since for odd-sized data it matters!),
    # and then followed by the FFT.
    vf = ftshift(reducedToFull(fft(iftshift(vf, inplace=False))), inplace=False)
    vg = ftshift(reducedToFull(fft(iftshift(vg, inplace=False))), inplace=False)

    ff = abs(vf)
    ff = real(ff)
    gg = abs(vg)
    gg = real(gg)
    
    get_bw = lambda x: int(pow(2, int(ceil(log(2*x, 2)))))
    
    res = None
    _last_bw = 0
    for r in range(1, radius+1):
        # calculate the appropriate bw
        bw = get_bw(r)
        if bw < b_min:
            bw = b_min
        if bw > b_max:
            bw = b_max
            
        # construct the wedge masks accordingly
        if _last_bw != bw:
            mf = create_wedge_sf(wf[0], wf[1], bw)
            mg = create_wedge_sf(wg[0], wg[1], bw)
        
        corr = frm_constrained_corr(vol2sf(ff, r, bw), mf, vol2sf(gg, r, bw), mg, norm=False, return_score=True)
        
        if _last_bw != bw:
            if res is None:
                res = np.zeros((2*bw, 2*bw, 2*bw), dtype='double')
            else:
                res = enlarge2(res)
        
        res += corr*(r**2)*weights[r-1]
        
        _last_bw = bw

    return frm_find_topn_angles_interp(res)

def frm_fourier_adaptive_wedge_vol_rscore(vf, wf, vg, wg, b, radius=None, weights=None):
    """Obsolete.
    """
    if not radius: # set the radius
        radius = vf.size_x()/2
    if not weights: # set the weights
        weights = [1 for i in range(radius)]

    if not b: # set the bandwidth adaptively
        b_min = 4
        b_max = 128
    elif b.__class__ == tuple or b.__class__ == list:
        b_min = b[0]
        b_max = b[1]
    elif isinstance(b, int): # fixed bandwidth
        b_min = b
        b_max = b
    else:
        raise RuntimeError("Argument b is not valid!")

    from pytom.basic.fourier import fft, ifft, ftshift, iftshift
    from pytom.lib.pytom_volume import vol, reducedToFull, real, imag, rescale
    from .vol2sf import vol2sf
    from pytom.lib.pytom_numpy import vol2npy
    from math import log, ceil, pow

    # IMPORTANT!!! Should firstly do the IFFTSHIFT on the volume data (NOT FFTSHIFT since for odd-sized data it matters!),
    # and then followed by the FFT.
    vf = ftshift(reducedToFull(fft(iftshift(vf, inplace=False))), inplace=False)
    vg = ftshift(reducedToFull(fft(iftshift(vg, inplace=False))), inplace=False)

    vfr = real(vf)
    vfi = imag(vf)
    vgr = real(vg)
    vgi = imag(vg)

    get_bw = lambda x: int(pow(2, int(ceil(log(2*x, 2)))))

    res = None
    _last_bw = 0
    for r in range(1, radius+1):
        # calculate the appropriate bw
        bw = get_bw(r)
        if bw < b_min:
            bw = b_min
        if bw > b_max:
            bw = b_max

        # construct the wedge masks accordingly
        if _last_bw != bw:
            mf = create_wedge_sf(wf[0], wf[1], bw)
            mg = create_wedge_sf(wg[0], wg[1], bw)

        corr = frm_fourier_constrained_corr(vol2sf(vfr, r, bw), vol2sf(vfi, r, bw), mf, vol2sf(vgr, r, bw), vol2sf(vgi, r, bw), mg, True, False, True)
        
        if _last_bw != bw:
            if res is None:
                res = np.zeros((2*bw, 2*bw, 2*bw), dtype='double')
            else:
                res = enlarge2(res)

        res += corr*(r**2)*weights[r-1]

        _last_bw = bw

    return res

def frm_align_vol_rscore(vf, wf, vg, wg, b, radius=None, mask=None, peak_offset=None, weights=None, position=None):
    """Obsolete.
    """
    from pytom.lib.pytom_volume import vol, rotateSpline, peak, initSphere
    from pytom.basic.transformations import shift
    from pytom.basic.correlation import xcf
    from pytom.basic.filter import lowpassFilter
    from pytom.basic.structures import Mask
    from pytom.lib.pytom_numpy import vol2npy

    if vf.size_x()!=vg.size_x() or vf.size_y()!=vg.size_y() or vf.size_z()!=vg.size_z():
        raise RuntimeError('Two volumes must have the same size!')

    if mask is None:
        mask = vol(vf.size_x(), vf.size_y(), vf.size_z())
        initSphere(mask, vf.size_x()/2, 0,0, vf.size_x()/2,vf.size_y()/2,vf.size_z()/2)
    elif mask.__class__ == vol:
        pass
    elif mask.__class__ == Mask:
        mask = mask.getVolume()
    elif isinstance(mask, int):
        mask_radius = mask
        mask = vol(vf.size_x(), vf.size_y(), vf.size_z())
        initSphere(mask, mask_radius, 0,0, vf.size_x()/2,vf.size_y()/2,vf.size_z()/2)
    else:
        raise RuntimeError('Given mask has wrong type!')

    if peak_offset is None:
        peak_offset = vol(vf.size_x(), vf.size_y(), vf.size_z())
        initSphere(peak_offset, vf.size_x()/2, 0,0, vf.size_x()/2,vf.size_y()/2,vf.size_z()/2)
    elif isinstance(peak_offset, int):
        peak_radius = peak_offset
        peak_offset = vol(vf.size_x(), vf.size_y(), vf.size_z())
        initSphere(peak_offset, peak_radius, 0,0, vf.size_x()/2,vf.size_y()/2,vf.size_z()/2)
    elif peak_offset.__class__ == vol:
        pass
    else:
        raise RuntimeError('Peak offset is given wrong!')

    # cut out the outer part which normally contains nonsense
    vf = vf*mask

    # # normalize them first
    # from pytom.basic.normalise import mean0std1
    # mean0std1(vf)
    # mean0std1(vg)

    if position is None: # if position is not given, we have to find it ourself
        # first roughtly determine the orientation (only according to the energy info)
        # get multiple candidate orientations
        orientations = frm_determine_orientation_rscore(vf, wf, vg, wg, b, radius, weights)
    else:
        # the position is given by the user
        vf2 = shift(vf, -position[0]+vf.size_x()/2, -position[1]+vf.size_y()/2, -position[2]+vf.size_z()/2, 'spline')
        res = frm_fourier_adaptive_wedge_vol_rscore(vf2, wf, vg, wg, b, radius, weights)
        orientation, max_value = frm_find_best_angle_interp(res)

        return position, orientation, max_value

    # iteratively refine the position & orientation
    from pytom.basic.structures import WedgeInfo
    from pytom.tools.maths import euclidianDistance
    max_iter = 10 # maximal number of iterations
    wedge = WedgeInfo([90+wf[0], 90-wf[1]])
    old_pos = [-1, -1, -1]
    vg2 = vol(vg.size_x(), vg.size_y(), vg.size_z())
    lowpass_vf = lowpassFilter(vf, radius, 0)[0]
    for i in range(max_iter):
        peak_value = 0.0
        position = None
        for orientation in orientations:
            orientation = orientation[0]

            rotateSpline(vg, vg2, orientation[0], orientation[1], orientation[2]) # first rotate
            vg2 = wedge.apply(vg2) # then apply the wedge
            vg2 = lowpassFilter(vg2, radius, 0)[0]
            res = xcf(lowpass_vf, vg2) # find the position
            pos = peak(res, peak_offset)
            # val = res(pos[0], pos[1], pos[2])
            pos, val = find_subpixel_peak_position(vol2npy(res), pos)
            if val > peak_value:
                position = pos
                peak_value = val

        if euclidianDistance(position, old_pos) <= 1.0:
            break
        else:
            old_pos = position

        # here we shift the target, not the reference
        # if you dont want the shift to change the energy landscape, use fourier shift
        vf2 = shift(vf, -position[0]+vf.size_x()/2, -position[1]+vf.size_y()/2, -position[2]+vf.size_z()/2, 'fourier')
        res = frm_fourier_adaptive_wedge_vol_rscore(vf2, wf, vg, wg, b, radius, weights)
        orientations = frm_find_topn_angles_interp(res)

    return position, orientations[0][0], orientations[0][1]


# ------------------------------------------------------------------------
def bart_align_vol(vf, wf, vg, wg, b, radius=None, peak_offset=None):
    """Implementation of Bartesaghi's approach for alignment. For detail, please check the paper.

    Parameters
    ----------
    vf: The volume you want to match.
        pytom.lib.pytom_volume.vol

    wf: The single tilt wedge information of volume vf.
        [missing_wedge_angle1, missing_wedge_angle2]. Note this is defined different with frm_align im frm.py!

    vg: The reference volume.
        pytom.lib.pytom_volume.vol

    wg: The single tilt wedge information of volume vg.
        [missing_wedge_angle1, missing_wedge_angle2]. Note this is defined different with frm_align im frm.py!

    b: The bandwidth of spherical harmonics.
       Integer in the range [4, 64]

    radius: The maximal radius in the Fourier space, which is equal to say the maximal frequency involved in calculation.
            Integer. By default is half of the volume size.

    peak_offset: The maximal offset which allows the peak of the score to be.
                 Or simply speaking, the maximal distance allowed to shift vg to match vf.
                 This parameter is needed to prevent shifting the reference volume out of the frame.
                 Integer. By default is half of the volume size.

    Returns
    -------
    The best translation and rotation (Euler angle, ZXZ convention [Phi, Psi, Theta]) to transform vg to match vf.
    (best_translation, best_rotation, correlation_score)
    """
    from pytom.lib.pytom_volume import vol, rotateSpline, max, peak, initSphere
    from pytom.basic.correlation import norm_xcf
    from pytom.basic.filter import lowpassFilter
    from pytom.basic.structures import WedgeInfo

    if not radius: # set the radius
        radius = vf.size_x()/2

    if peak_offset is None:
        peak_offset = vol(vf.size_x(), vf.size_y(), vf.size_z())
        initSphere(peak_offset, vf.size_x()/2, 0,0, vf.size_x()/2,vf.size_y()/2,vf.size_z()/2)
    elif isinstance(peak_offset, int):
        peak_radius = peak_offset
        peak_offset = vol(vf.size_x(), vf.size_y(), vf.size_z())
        initSphere(peak_offset, peak_radius, 0,0, vf.size_x()/2,vf.size_y()/2,vf.size_z()/2)
    elif peak_offset.__class__ == vol:
        pass
    else:
        raise RuntimeError('Peak offset is given wrong!')
    
    from pytom.basic.fourier import fft, ifft, ftshift, iftshift
    from pytom.lib.pytom_volume import vol, reducedToFull, rescale, abs, real
    from .vol2sf import vol2sf
    from pytom.lib.pytom_numpy import vol2npy
    from math import log, ceil, pow

    # IMPORTANT!!! Should firstly do the IFFTSHIFT on the volume data (NOT FFTSHIFT since for odd-sized data it matters!),
    # and then followed by the FFT.
    ff = abs(ftshift(reducedToFull(fft(iftshift(vf, inplace=False))), inplace=False))
    ff = real(ff)
    gg = abs(ftshift(reducedToFull(fft(iftshift(vg, inplace=False))), inplace=False))
    gg = real(gg)
    
    sf = None
    sg = None
    mf = create_wedge_sf(wf[0], wf[1], b)
    mg = create_wedge_sf(wg[0], wg[1], b)

    for r in range(3, radius+1): # Should start from 3 since the interpolation in the first 2 bands is not accurate.
        if sf is None:
            sf = vol2sf(ff, r, b)
            sg = vol2sf(gg, r, b)
        else:
            sf += vol2sf(ff, r, b)
            sg += vol2sf(gg, r, b)
    
    corr = frm_constrained_corr(sf, mf, sg, mg)
    ang, val = frm_find_best_angle_interp(corr)

    tmp = vol(vg.size_x(),vg.size_y(),vg.size_z())
    rotateSpline(vg, tmp, ang[0], ang[1], ang[2])
    wedge_f = WedgeInfo(90+wf[0], 90-wf[1])
    wedge_g = WedgeInfo(90+wg[0], 90-wg[1])
    cc = norm_xcf(lowpassFilter(wedge_g.apply(vf), radius, 0)[0], lowpassFilter(wedge_f.apply(tmp), radius, 0)[0])
    pos = peak(cc, peak_offset)
    pos, score = find_subpixel_peak_position(vol2npy(cc), pos)

    return (pos, ang, score)

# ------------------------------------------------------------------------
def frm_fourier_adaptive_wedge_vol(vf, wf, vg, wg, b, radius=None, weights=None, r_score=False, norm=False):
    """Auxiliary function for xu_align_vol. Assuming no translation, return the correlation function of Euler angle.

    Parameters
    ----------
    vf: The volume you want to match.
        pytom.lib.pytom_volume.vol

    wf: The single tilt wedge information of volume vf.
        [missing_wedge_angle1, missing_wedge_angle2]. Note this is defined different with frm_align im frm.py!

    vg: The reference volume.
        pytom.lib.pytom_volume.vol

    wg: The single tilt wedge information of volume vg.
        [missing_wedge_angle1, missing_wedge_angle2]. Note this is defined different with frm_align im frm.py!

    b: The adaptive bandwidth of spherical harmonics.
       List [min_bandwidth, max_bandwidth], min_bandwidth, max_bandwidth in the range [4, 64].
       Or integer, which would then mean to use fixed bandwidth: min_bandwidth = max_bandwidth = integer.

    radius: The maximal radius in the Fourier space, which is equal to say the maximal frequency involved in calculation.
            Integer. By default is half of the volume size.

    weights: Obsolete.

    r_score: Obsolete.

    norm: Obsolete.

    Returns
    -------
    The correlation function of Euler angle.
    """
    if not radius: # set the radius
        radius = vf.size_x()/2
    if not weights: # set the weights
        weights = [1 for i in range(radius)]

    if not b: # set the bandwidth adaptively
        b_min = 4
        b_max = 128
    elif b.__class__ == tuple or b.__class__ == list:
        b_min = b[0]
        b_max = b[1]
    elif isinstance(b, int): # fixed bandwidth
        b_min = b
        b_max = b
    else:
        raise RuntimeError("Argument b is not valid!")

    from pytom.basic.fourier import fft, ifft, ftshift, iftshift
    from pytom.lib.pytom_volume import vol, reducedToFull, real, imag, rescale
    from .vol2sf import vol2sf
    from pytom.lib.pytom_numpy import vol2npy
    from math import log, ceil, pow

    # IMPORTANT!!! Should firstly do the IFFTSHIFT on the volume data (NOT FFTSHIFT since for odd-sized data it matters!),
    # and then followed by the FFT.
    vf = ftshift(reducedToFull(fft(iftshift(vf, inplace=False))), inplace=False)
    vg = ftshift(reducedToFull(fft(iftshift(vg, inplace=False))), inplace=False)

    vfr = real(vf)
    vfi = imag(vf)
    vgr = real(vg)
    vgi = imag(vg)

    get_bw = lambda x: int(pow(2, int(ceil(log(2*x, 2)))))

    numerator = None
    denominator1 = None
    denominator2 = None
    _last_bw = 0
    for r in range(1, radius+1):
        # calculate the appropriate bw
        bw = get_bw(r)
        if bw < b_min:
            bw = b_min
        if bw > b_max:
            bw = b_max

        # construct the wedge masks accordingly
        if _last_bw != bw:
            mf = create_wedge_sf(wf[0], wf[1], bw)
            mg = create_wedge_sf(wg[0], wg[1], bw)

        corr1, corr2, corr3 = frm_fourier_constrained_corr(vol2sf(vfr, r, bw), vol2sf(vfi, r, bw), mf, vol2sf(vgr, r, bw), vol2sf(vgi, r, bw), mg, True, norm, False)
        
        if _last_bw != bw:
            if numerator is None:
                numerator = np.zeros((2*bw, 2*bw, 2*bw), dtype='double')
                denominator1 = np.zeros((2*bw, 2*bw, 2*bw), dtype='double')
                denominator2 = np.zeros((2*bw, 2*bw, 2*bw), dtype='double')
            else:
                numerator = enlarge2(numerator)
                denominator1 = enlarge2(denominator1)
                denominator2 = enlarge2(denominator2)

        numerator += corr1*(r**2)*weights[r-1]
        denominator1 += corr2*(r**2)*weights[r-1]
        denominator2 += corr3*(r**2)*weights[r-1]

        _last_bw = bw

    res = numerator/(denominator1 * denominator2)**0.5

    return res

def frm_determine_orientation(vf, wf, vg, wg, b, radius=None, weights=None, r_score=False, norm=False):
    """Auxiliary function for xu_align_vol. Find the angle to rotate vg to match vf, using only their power spectrums.

    Parameters
    ----------
    vf: The volume you want to match.
        pytom.lib.pytom_volume.vol

    wf: The single tilt wedge information of volume vf.
        [missing_wedge_angle1, missing_wedge_angle2]. Note this is defined different with frm_align im frm.py!

    vg: The reference volume.
        pytom.lib.pytom_volume.vol

    wg: The single tilt wedge information of volume vg.
        [missing_wedge_angle1, missing_wedge_angle2]. Note this is defined different with frm_align im frm.py!

    b: The adaptive bandwidth of spherical harmonics.
       List [min_bandwidth, max_bandwidth], min_bandwidth, max_bandwidth in the range [4, 64].
       Or integer, which would then mean to use fixed bandwidth: min_bandwidth = max_bandwidth = integer.

    radius: The maximal radius in the Fourier space, which is equal to say the maximal frequency involved in calculation.
            Integer. By default is half of the volume size.

    weights: Obsolete.

    r_score: Obsolete.

    norm: Obsolete.

    Returns
    -------
    The angle (Euler angle, ZXZ convention [Phi, Psi, Theta]) to rotate vg to match vf.
    """
    if not radius: # set the radius
        radius = vf.size_x()/2
    if not weights: # set the weights
        weights = [1 for i in range(radius)]
    
    if not b: # set the bandwidth adaptively
        b_min = 4
        b_max = 128
    elif b.__class__ == tuple or b.__class__ == list:
        b_min = b[0]
        b_max = b[1]
    elif isinstance(b, int): # fixed bandwidth
        b_min = b
        b_max = b
    else:
        raise RuntimeError("Argument b is not valid!")
    
    from pytom.basic.fourier import fft, ifft, ftshift, iftshift
    from pytom.lib.pytom_volume import vol, reducedToFull, rescale, abs, real
    from .vol2sf import vol2sf
    from pytom.lib.pytom_numpy import vol2npy
    from math import log, ceil, pow

    # IMPORTANT!!! Should firstly do the IFFTSHIFT on the volume data (NOT FFTSHIFT since for odd-sized data it matters!),
    # and then followed by the FFT.
    vf = ftshift(reducedToFull(fft(iftshift(vf, inplace=False))), inplace=False)
    vg = ftshift(reducedToFull(fft(iftshift(vg, inplace=False))), inplace=False)

    ff = abs(vf)
    ff = real(ff)
    gg = abs(vg)
    gg = real(gg)
    
    get_bw = lambda x: int(pow(2, int(ceil(log(2*x, 2)))))
    
    numerator = None
    denominator1 = None
    denominator2 = None
    _last_bw = 0
    for r in range(1, radius+1):
        # calculate the appropriate bw
        bw = get_bw(r)
        if bw < b_min:
            bw = b_min
        if bw > b_max:
            bw = b_max
            
        # construct the wedge masks accordingly
        if _last_bw != bw:
            mf = create_wedge_sf(wf[0], wf[1], bw)
            mg = create_wedge_sf(wg[0], wg[1], bw)
        
        corr1, corr2, corr3 = frm_constrained_corr(vol2sf(ff, r, bw), mf, vol2sf(gg, r, bw), mg, norm, return_score=False)
        
        if _last_bw != bw:
            if numerator is None:
                numerator = np.zeros((2*bw, 2*bw, 2*bw), dtype='double')
                denominator1 = np.zeros((2*bw, 2*bw, 2*bw), dtype='double')
                denominator2 = np.zeros((2*bw, 2*bw, 2*bw), dtype='double')
            else:
                numerator = enlarge2(numerator)
                denominator1 = enlarge2(denominator1)
                denominator2 = enlarge2(denominator2)
        
        numerator += corr1*(r**2)*weights[r-1]
        denominator1 += corr2*(r**2)*weights[r-1]
        denominator2 += corr3*(r**2)*weights[r-1]
        
        _last_bw = bw
    
    res = numerator/(denominator1 * denominator2)**0.5

    return frm_find_topn_angles_interp2(res)


def xu_align_vol(vf, wf, vg, wg, b, radius=None, mask=None, peak_offset=None):
    """Implementation of Xu's approach for alignment. For detail, please check the paper.

    Parameters
    ----------
    vf: The volume you want to match.
        pytom.lib.pytom_volume.vol

    wf: The single tilt wedge information of volume vf.
        [missing_wedge_angle1, missing_wedge_angle2]. Note this is defined different with frm_align im frm.py!

    vg: The reference volume.
        pytom.lib.pytom_volume.vol

    wg: The single tilt wedge information of volume vg.
        [missing_wedge_angle1, missing_wedge_angle2]. Note this is defined different with frm_align im frm.py!

    b: The adaptive bandwidth of spherical harmonics.
       List [min_bandwidth, max_bandwidth], min_bandwidth, max_bandwidth in the range [4, 64].
       Or integer, which would then mean to use fixed bandwidth: min_bandwidth = max_bandwidth = integer.

    radius: The maximal radius in the Fourier space, which is equal to say the maximal frequency involved in calculation.
            Integer. By default is half of the volume size.

    peak_offset: The maximal offset which allows the peak of the score to be.
                 Or simply speaking, the maximal distance allowed to shift vg to match vf.
                 This parameter is needed to prevent shifting the reference volume out of the frame.
                 Integer. By default is half of the volume size.

    Returns
    -------
    The best translation and rotation (Euler angle, ZXZ convention [Phi, Psi, Theta]) to transform vg to match vf.
    (best_translation, best_rotation, correlation_score)
    """
    from pytom.lib.pytom_volume import vol, rotateSpline, peak, initSphere
    from pytom.basic.transformations import shift
    from pytom.basic.correlation import norm_xcf
    from pytom.basic.filter import lowpassFilter
    from pytom.basic.structures import Mask
    from pytom.lib.pytom_numpy import vol2npy

    if vf.size_x()!=vg.size_x() or vf.size_y()!=vg.size_y() or vf.size_z()!=vg.size_z():
        raise RuntimeError('Two volumes must have the same size!')

    if mask is None:
        mask = vol(vf.size_x(), vf.size_y(), vf.size_z())
        initSphere(mask, vf.size_x()/2, 0,0, vf.size_x()/2,vf.size_y()/2,vf.size_z()/2)
    elif mask.__class__ == vol:
        pass
    elif mask.__class__ == Mask:
        mask = mask.getVolume()
    elif isinstance(mask, int):
        mask_radius = mask
        mask = vol(vf.size_x(), vf.size_y(), vf.size_z())
        initSphere(mask, mask_radius, 0,0, vf.size_x()/2,vf.size_y()/2,vf.size_z()/2)
    else:
        raise RuntimeError('Given mask has wrong type!')

    if peak_offset is None:
        peak_offset = vol(vf.size_x(), vf.size_y(), vf.size_z())
        initSphere(peak_offset, vf.size_x()/2, 0,0, vf.size_x()/2,vf.size_y()/2,vf.size_z()/2)
    elif isinstance(peak_offset, int):
        peak_radius = peak_offset
        peak_offset = vol(vf.size_x(), vf.size_y(), vf.size_z())
        initSphere(peak_offset, peak_radius, 0,0, vf.size_x()/2,vf.size_y()/2,vf.size_z()/2)
    elif peak_offset.__class__ == vol:
        pass
    else:
        raise RuntimeError('Peak offset is given wrong!')

    # cut out the outer part which normally contains nonsense
    vf = vf*mask

    position = None
    if position is None: # if position is not given, we have to find it ourself
        # first roughtly determine the orientation (only according to the energy info)
        # get multiple candidate orientations
        orientations = frm_determine_orientation(vf, wf, vg, wg, b, radius, None, None, False)
    else:
        # the position is given by the user
        vf2 = shift(vf, -position[0]+vf.size_x()/2, -position[1]+vf.size_y()/2, -position[2]+vf.size_z()/2, 'spline')
        res = frm_fourier_adaptive_wedge_vol(vf2, wf, vg, wg, b, radius, None, None, False)
        orientation, max_value = frm_find_best_angle_interp(res)

        return position, orientation, max_value

    
    from pytom.basic.structures import WedgeInfo
    from pytom.tools.maths import euclidianDistance
    max_iter = 1 # maximal number of iterations
    wedge = WedgeInfo([90+wf[0], 90-wf[1]])
    old_pos = [-1, -1, -1]
    vg2 = vol(vg.size_x(), vg.size_y(), vg.size_z())
    lowpass_vf = lowpassFilter(vf, radius, 0)[0]
    
    peak_value = 0.0
    position = None
    ang = None
    for orientation in orientations:
        orientation = orientation[0]

        rotateSpline(vg, vg2, orientation[0], orientation[1], orientation[2]) # first rotate
        vg2 = wedge.apply(vg2) # then apply the wedge
        vg2 = lowpassFilter(vg2, radius, 0)[0]
        res = norm_xcf(lowpass_vf, vg2) # find the position
        pos = peak(res, peak_offset)
        # val = res(pos[0], pos[1], pos[2])
        pos, val = find_subpixel_peak_position(vol2npy(res), pos)
        if val > peak_value:
            position = pos
            ang = orientation
            peak_value = val

    return position, ang, peak_value

