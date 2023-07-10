

def extractPeaks(volume, reference, rotations, scoreFnc=None, mask=None, maskIsSphere=False, wedgeInfo=None, debug=False,
                 **kwargs):
    '''
    Created on May 17, 2010
    @param volume: target volume
    @type volume: L{pytom.lib.pytom_volume.vol}
    @param reference: reference
    @type reference: L{pytom.lib.pytom_volume.vol}
    @param rotations: rotation angle list
    @type rotations: L{pytom.angles.globalSampling.GlobalSampling}
    @param scoreFnc: score function that is used
    @type scoreFnc: L{pytom.basic.correlation}
    @param mask: mask volume
    @type mask: L{pytom.lib.pytom_volume.vol}
    @param maskIsSphere: flag to indicate whether the mask is sphere or not
    @type maskIsSphere: boolean
    @param wedgeInfo: wedge information
    @type wedgeInfo: L{pytom.basic.structures.WedgeInfo}
    @return: both the score volume and the corresponding rotation index volume
    @rtype: L{pytom.lib.pytom_volume.vol}
    @author: chen
    '''
    nodeName = kwargs.get('nodeName', '')
    verbose = kwargs.get('verboseMode', True)
    if verbose not in [True, False]:
        verbose = True
    moreInfo = kwargs.get('moreInfo', False)
    if moreInfo not in [True, False]:
        moreInfo = False
    
    from pytom.basic.correlation import flcf
    from pytom.basic.structures import WedgeInfo, Wedge
    from pytom.lib.pytom_volume import vol, pasteCenter
    from pytom.lib.pytom_volume import rotateSpline as rotate  # for more accuracy
    from pytom.lib.pytom_volume import updateResFromIdx

    if scoreFnc == None:
        scoreFnc = flcf
    
    # only FLCF needs mask
    if scoreFnc == flcf:
        if mask.__class__ != vol: # construct a sphere mask by default
            from pytom.lib.pytom_volume import initSphere
            mask = vol(reference.size_x(), reference.size_y(), reference.size_z())
            mask.setAll(0)
            initSphere(mask, reference.size_x()/2,0,0,reference.size_x()/2, reference.size_x()/2,reference.size_x()/2)
            maskIsSphere = True
    
    # result volume which stores the score
    result = vol(volume.size_x(), volume.size_y(), volume.size_z())
    result.setAll(-1)
    
    # result orientation of the peak value (index)
    orientation = vol(volume.size_x(), volume.size_y(), volume.size_z())
    orientation.setAll(0)
    
    currentRotation = rotations.nextRotation()
    index = 0
    
    if verbose == True:
        from pytom.tools.ProgressBar import FixedProgBar
        max = rotations.numberRotations()-1
        prog = FixedProgBar(0, max, nodeName)
    if moreInfo:
        sumV = vol(volume.size_x(), volume.size_y(), volume.size_z())
        sumV.setAll(0)
        sqrV = vol(volume.size_x(), volume.size_y(), volume.size_z())
        sqrV.setAll(0)
    else:
        sumV = None
        sqrV = None

    if wedgeInfo.__class__ == WedgeInfo or wedgeInfo.__class__ == Wedge:
        print('Applied wedge to volume')
        volume = wedgeInfo.apply(volume)

    while currentRotation != [None,None,None]:
        if verbose == True:
            prog.update(index)
        from pytom.basic.files import read
        # rotate the reference
        ref = vol(reference.size_x(),reference.size_y(),reference.size_z())
        rotate(reference, ref, currentRotation[0], currentRotation[1], currentRotation[2])

        if debug: ref.write('rot_cpu.em')

        # apply wedge
        if wedgeInfo.__class__ == WedgeInfo or wedgeInfo.__class__ == Wedge:
            ref = wedgeInfo.apply(ref)

        if debug: ref.write('wedge_rot_cpu.em')

        # rotate the mask if it is asymmetric
        if scoreFnc == flcf:
            if maskIsSphere == False: # if mask is not a sphere, then rotate it
                m = vol(mask.size_x(),mask.size_y(),mask.size_z())
                rotate(mask, m, currentRotation[0], currentRotation[1], currentRotation[2])
            else:
                m = mask
        
        # compute the score
        # if mask is sphere and it is the first run,
        # compute the standard deviation of the volume under mask for late use
        if scoreFnc == flcf and index == 0 and maskIsSphere == True:
            # compute standard deviation of the volume under mask
            maskV = m
            if volume.size_x() != m.size_x() or volume.size_y() != m.size_y() or volume.size_z() != m.size_z():
                maskV = vol(volume.size_x(), volume.size_y(), volume.size_z())
                maskV.setAll(0)
                pasteCenter(m, maskV)
            from pytom.lib.pytom_volume import sum
            p = sum(m);
            from pytom.basic.correlation import meanUnderMask, stdUnderMask
            meanV = meanUnderMask(volume, maskV, p)
            std_v = stdUnderMask(volume, maskV, p, meanV)
            if debug:
                volume.write('volume_cpu.mrc')
                meanV.write('meanV_cpu.mrc')

        if scoreFnc == flcf:
            if maskIsSphere == True:
                score = scoreFnc(volume, ref, m, std_v, wedge=1)
            else:
                score = scoreFnc(volume, ref, m)
        else: # not FLCF, so doesn't need mask as parameter and perhaps the reference should have the same size
            _ref = vol(volume.size_x(), volume.size_y(), volume.size_z())
            _ref.setAll(0)
            pasteCenter(ref, _ref)
            
            score = scoreFnc(volume, _ref)

        if debug: score.write('cmap_cpu2.mrc')

        # update the result volume and the orientation volume
        updateResFromIdx(result, score, orientation, index)
        
        if moreInfo:
            sumV = sumV + score
            sqrV = sqrV + score*score
        
        currentRotation = rotations.nextRotation()
        index = index+1

    return [result, orientation, sumV, sqrV]


def create_structured_wedge(tilt_angles, angle2=None, cutoffRadius=10, size_x=10, size_y=10, size_z=10, smooth=0, rotation=None, c=1):
    from pytom.gpu.initialize import xp
    print(tilt_angles)
    z, y, x = xp.meshgrid(xp.arange(-size_y // 2 + size_y % 2, size_y // 2 + size_y % 2),
                          xp.arange(-size_x // 2 + size_x % 2, size_x // 2 + size_x % 2),
                          xp.arange(0, size_z // 2 + 1))

    r = xp.sqrt((x*size_x/size_z) ** 2 + (y) ** 2 + (z*size_x/size_y) ** 2)

    tot = xp.zeros_like(z)

    for i in tilt_angles:
        alpha = (i+90)*xp.pi/180
        if abs(i) < 0.001:
            tot += 1*(abs(x) < c)
        else:
            tot += (abs(xp.sin(alpha)*(x-y/xp.tan(alpha))) <= c)*1

    tot[tot>1] = 1
    tot[r > cutoffRadius] = 0
    return xp.fft.fftshift(tot, axes=(0, 1))


def templateMatchingGPU(volume, reference, rotations, scoreFnc=None, mask=None, maskIsSphere=False, wedgeInfo=None, padding=True, jobid=0, **kwargs):
    '''
    Created on May 17, 2020
    @param volume: target volume
    @type volume: numpy ndarray
    @param reference: reference
    @type reference: numpy ndarray
    @param rotations: rotation angle list
    @type rotations: L{pytom.angles.globalSampling.GlobalSampling}
    @param scoreFnc: score function that is used (currently not used)
    @type scoreFnc: L{pytom.basic.correlation}
    @param mask: mask volume
    @type mask: numpy ndarray
    @param maskIsSphere: flag to indicate whether the mask is sphere or not. (currently assumed to be true)
    @type maskIsSphere: boolean
    @param wedgeInfo: wedge information
    @type wedgeInfo: L{pytom.basic.structures.WedgeInfo}
    @param padding: pad the volume to a size that is fast for fft's. Example: 127 should be padded to 128.
    @type padding: boolean
    @return: both the score volume and the corresponding rotation index volume
    @rtype: two numpy ndarray's
    @author: GvdS
    '''

    from pytom.agnostic.filter import create_wedge
    from pytom.gpu.gpuStructures import TemplateMatchingGPU
    from pytom.tools.calcFactors import calc_fast_gpu_dimensions
    import numpy as np
    from pytom.gpu.initialize import xp, device

    if not kwargs['gpuID'] is None:
        import cupy as xp

    xp.cuda.Device(kwargs['gpuID']).use()

    angles = rotations[:]
    SX,SY,SZ = volume.shape
    sx,sy,sz = reference.shape
    angle = wedgeInfo.getWedgeAngle()

    # reference -= reference.max()
    print('\n\nWEDGEANGLE', angle)
    if angle.__class__ != list:
        w1 = angle
        w2 = angle
    else:
        w1, w2 = angle

    if w1 > 1E-3 or w2 > 1E-3:
        # cutoff was previously sx // 2 - 1
        # replace with Wedge.convert2numpy() and the returnWedgeVolume
        cutoff = wedgeInfo._wedgeObject._cutoffRadius
        smooth = wedgeInfo._wedgeObject._smooth
        wedge = create_wedge(w1, w2, cutoff, sx, sy, sz, smooth).get().astype(np.float32)
        wedge_tomogram = create_wedge(w1, w2, cutoff, SX, SY, SZ, smooth).get().astype(np.float32)

        # convolve the search volume with the wedge
        volume = np.real(np.fft.irfftn(np.fft.rfftn(volume) * wedge_tomogram, s=volume.shape))

        del wedge_tomogram
        print('Wedge filter applied to volume')
    else:
        wedge = np.ones((sx,sy,sz//2+1),dtype='float32')
    
    print('dim of tomogram ', volume.shape)
    dimx, dimy, dimz = volume.shape
    if padding:
        cx = max(reference.shape[0], calc_fast_gpu_dimensions(SX - 2, 4000)[0])
        cy = max(reference.shape[1], calc_fast_gpu_dimensions(SY - 2, 4000)[0])
        cz = max(reference.shape[2], calc_fast_gpu_dimensions(SZ - 2, 4000)[0])
        voluNDAs = np.zeros([cx, cy, cz], dtype=np.float32)
        voluNDAs[:min(cx, SX), :min(cy, SY), :min(cz, SZ)] = volume[:min(cx, SX), :min(cy, SY),:min(cz, SZ)]
        volume = voluNDAs

    print(f'dimensions of template and mask: {reference.shape} {mask.shape} ')

    input = (volume, reference, mask, wedge, angles, maskIsSphere)

    tm_process = TemplateMatchingGPU(jobid, kwargs['gpuID'], input=input)
    tm_process.start()

    import time
    sleep_time, max_sleep_time = 0, 3600 * 12  # set max runtime for gpu to 12 hours
    while tm_process.is_alive() and sleep_time < max_sleep_time:
        time.sleep(1)
        sleep_time += 1

    if tm_process.completed:
        print(f'Templated matching completed successfully on {kwargs["gpuID"]}')
        if padding:
            # change shape to original size
            scrs, angs = np.zeros((SX, SY, SZ), dtype=np.float32), np.zeros((SX, SY, SZ), dtype=np.float32)
            scrs[:min(cx, SX), :min(cy, SY), :min(cz, SZ)] = tm_process.plan.scores.get()[:min(cx, SX),
                                                                   :min(cy, SY), :min(cz, SZ)]
            angs[:min(cx, SX), :min(cy, SY), :min(cz, SZ)] = tm_process.plan.angles.get()[:min(cx, SX),
                                                                   :min(cy, SY), :min(cz, SZ)]
            return [scrs, angs, None, None]
        else:
            return [tm_process.plan.scores.get(), tm_process.plan.angles.get(), None, None]
    else:
        if sleep_time >= max_sleep_time:
            print(f'Job terminated on {kwargs["gpuID"]} due to time limit (exceeded {max_sleep_time} sec).')

        raise Exception('failed template matching')
