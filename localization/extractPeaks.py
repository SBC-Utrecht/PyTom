
#def updateResult(resVol, newVol, orientVol, index):
#    '''
#    Created on May 17, 2010
#    @param resVol: result score volume that is about to update
#    @type resVol: L{pytom_volume.vol}
#    @param newVol: the new score volume that is about to compare to the result volume
#    @type newVol: L{pytom_volume.vol}
#    @param orientvol: orientation volume that stores the index of the rotation information
#    @type orientVol: L{pytom_volume.vol}
#    @param index: index of the loop
#    @type index: int
#    @author: chen
#    '''
#    
#    for i in xrange(resVol.sizeX()):
#        for j in xrange(resVol.sizeY()):
#            for k in xrange(resVol.sizeZ()):
#                if resVol.getV(i,j,k) < newVol.getV(i,j,k):
#                    resVol.setV(newVol.getV(i,j,k), i,j,k)
#                    orientVol.setV(index, i,j,k)
    
def extractPeaks(volume, reference, rotations, scoreFnc=None, mask=None, maskIsSphere=False, wedgeInfo=None, **kwargs):
    '''
    Created on May 17, 2010
    @param volume: target volume
    @type volume: L{pytom_volume.vol}
    @param reference: reference
    @type reference: L{pytom_volume.vol}
    @param rotations: rotation angle list
    @type rotations: L{pytom.angles.globalSampling.GlobalSampling}
    @param scoreFnc: score function that is used
    @type scoreFnc: L{pytom.basic.correlation}
    @param mask: mask volume
    @type mask: L{pytom_volume.vol}
    @param maskIsSphere: flag to indicate whether the mask is sphere or not
    @type maskIsSphere: boolean
    @param wedgeInfo: wedge information
    @type wedgeInfo: L{pytom.basic.structures.WedgeInfo}
    @return: both the score volume and the corresponding rotation index volume
    @rtype: L{pytom_volume.vol}
    @author: chen
    '''
#    from pytom.tools.timing import timing
#    t = timing(); t.start()
    
    # parse the parameters
    nodeName = kwargs.get('nodeName', '')
    verbose = kwargs.get('verboseMode', True)
    if verbose not in [True, False]:
        verbose = True
    moreInfo = kwargs.get('moreInfo', False)
    if moreInfo not in [True, False]:
        moreInfo = False
    
    from pytom.basic.correlation import FLCF
    from pytom.basic.structures import WedgeInfo, Wedge
    from pytom_volume import vol, pasteCenter
    from pytom_volume import rotateSpline as rotate # for more accuracy
    from pytom_volume import updateResFromIdx
    from pytom.basic.files import write_em

    if scoreFnc == None:
        scoreFnc = FLCF
    
    # only FLCF needs mask
    if scoreFnc == FLCF:
        if mask.__class__ != vol: # construct a sphere mask by default
            from pytom_volume import initSphere;
            mask = vol(reference.sizeX(), reference.sizeY(), reference.sizeZ());
            mask.setAll(0);
            initSphere(mask, reference.sizeX()/2,0,0,reference.sizeX()/2,
	        reference.sizeX()/2,reference.sizeX()/2);
            maskIsSphere = True
    
    # result volume which stores the score
    result = vol(volume.sizeX(), volume.sizeY(), volume.sizeZ())
    result.setAll(-1)
    
    # result orientation of the peak value (index)
    orientation = vol(volume.sizeX(), volume.sizeY(), volume.sizeZ())
    orientation.setAll(0)
    
    currentRotation = rotations.nextRotation()
    index = 0
    
    if verbose == True:
        from pytom.tools.ProgressBar import FixedProgBar
        max = rotations.numberRotations()-1
        prog = FixedProgBar(0, max, nodeName)
    if moreInfo:
        sumV = vol(volume.sizeX(), volume.sizeY(), volume.sizeZ())
        sumV.setAll(0)
        sqrV = vol(volume.sizeX(), volume.sizeY(), volume.sizeZ())
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
        
        # rotate the reference
        ref = vol(reference.sizeX(),reference.sizeY(),reference.sizeZ())
        rotate(reference, ref, currentRotation[0], currentRotation[1], currentRotation[2])

        # apply wedge
        if wedgeInfo.__class__ == WedgeInfo or wedgeInfo.__class__ == Wedge:
            ref = wedgeInfo.apply(ref)

        # rotate the mask if it is asymmetric
        if scoreFnc == FLCF:
            if maskIsSphere == False: # if mask is not a sphere, then rotate it
                m = vol(mask.sizeX(),mask.sizeY(),mask.sizeZ())
                rotate(mask, m, currentRotation[0], currentRotation[1], currentRotation[2])
            else:
                m = mask
        
        # compute the score
        # if mask is sphere and it is the first run, compute the standard deviation of the volume under mask for late use
        if scoreFnc == FLCF and index == 0 and maskIsSphere == True:
            # compute standard deviation of the volume under mask
            maskV = m
            if volume.sizeX() != m.sizeX() or volume.sizeY() != m.sizeY() or volume.sizeZ() != m.sizeZ():
                maskV = vol(volume.sizeX(), volume.sizeY(), volume.sizeZ())
                maskV.setAll(0)
                pasteCenter(m, maskV)
            from pytom_volume import sum
            p = sum(m);
            from pytom.basic.correlation import meanUnderMask, stdUnderMask
            meanV = meanUnderMask(volume, maskV, p);
            stdV = stdUnderMask(volume, maskV, p, meanV);

        # ref.write('template_cpu.em')

        if scoreFnc == FLCF:
            if maskIsSphere == True:
                score = scoreFnc(volume, ref, m, stdV, wedge=1)
            else:
                score = scoreFnc(volume, ref, m)
        else: # not FLCF, so doesn't need mask as parameter and perhaps the reference should have the same size
            _ref = vol(volume.sizeX(), volume.sizeY(), volume.sizeZ())
            _ref.setAll(0)
            pasteCenter(ref, _ref)
            
            score = scoreFnc(volume, _ref)

        # update the result volume and the orientation volume
        updateResFromIdx(result, score, orientation, index)
        
        if moreInfo:
            sumV = sumV + score
            sqrV = sqrV + score*score
        
        currentRotation = rotations.nextRotation()
        index = index+1
        
#    if moreInfo:
#        sumV = sumV/rotations.numberRotations()
#        sqrV = sqrV/rotations.numberRotations()

#    time = t.end(); print 'The overall execution time: %f' % time

    return [result, orientation, sumV, sqrV]

def extractPeaksGPU(volume, reference, rotations, scoreFnc=None, mask=None, maskIsSphere=False, wedgeInfo=None, padding=True, **kwargs):
    from pytom_numpy import vol2npy
    from pytom.tompy.filter import create_wedge, applyFourierFilter
    from pytom.tompy.io import write
    from pytom.gpu.gpuStructures import TemplateMatchingGPU
    from pytom.tools.calcFactors import calc_fast_gpu_dimensions
    import time
    import numpy as np
    from pytom_freqweight import weight




    angles = rotations[:]
    #volume = wedgeInfo.apply(volume)

    volume, ref, mask = [vol2npy(vol) for vol in (volume, reference, mask)]
    SX,SY,SZ = volume.shape
    sx,sy,sz = ref.shape
    angle = wedgeInfo.getWedgeAngle()

    if angle.__class__ != list:
        w1 = w2 = angle
    else:
        w1,w1 = angle

    if w1 > 1E-3 or w2 > 1E-3:
        print('Wedge applied to volume')
        cutoff = wedgeInfo._wedgeObject._cutoffRadius if wedgeInfo._wedgeObject._cutoffRadius > 1E-3 else sx//2
        smooth = wedgeInfo._wedgeObject._smooth
        wedge = create_wedge(w1, w2, cutoff, sx, sy, sz, smooth).astype(np.complex64)
        wedgeVolume = create_wedge(w1, w2, SX//2-2, SX, SY, SZ, smooth).astype(np.float32)
        volume = np.real(np.fft.irfftn(np.fft.rfftn(volume) * wedgeVolume.get()))

    else:
        wedge = np.ones((sx,sy,sz//2+1),dtype='float32')
        #wedgeVolume = np.ones_like(volume,dtype='float32')

    if padding:
        dimx, dimy, dimz = volume.shape

        cx = max(ref.shape[0], calc_fast_gpu_dimensions(dimx - 2, 4000)[0])
        cy = max(ref.shape[1], calc_fast_gpu_dimensions(dimy - 2, 4000)[0])
        cz = max(ref.shape[2], calc_fast_gpu_dimensions(dimz - 2, 4000)[0])
        voluNDAs = np.zeros([cx, cy, cz], dtype=np.float32)
        voluNDAs[:min(cx, dimx), :min(cy, dimy), :min(cz, dimz)] = volume[:min(cx, dimx), :min(cy, dimy),
                                                                   :min(cz, dimz)]
        volume = voluNDAs
        print(f'dimensions of volume: {ref.shape} {mask.shape} ')

    input = (volume, ref, mask, wedge, angles, volume.shape)

    tm_process = TemplateMatchingGPU(0, kwargs['gpuID'][0], input=input)
    tm_process.start()
    while tm_process.is_alive():
        time.sleep(1)

    if tm_process.completed:
        print('Templated matching completed successfully')
        return [tm_process.plan.scores.get(), tm_process.plan.angles.get(), None, None]
    else:
        raise Exception('failed template matching')
