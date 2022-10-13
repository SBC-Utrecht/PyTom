def simpleSimulation(volume,rotation,shiftV,wedgeInfo=None,SNR=0.1,mask=None):
    """
    simpleSimulation: Simulates an ET by applying rotation,shift,wedge and noise to an volume
    
    @param volume: the volume used for simulations
    @param rotation: the rotation applied to volume
    @param shiftV: shift vector applied to volume
    @param wedgeInfo: wedge applied to volume
    @param SNR: noise level applied to volume
    @param mask: Apodisation mask 
    @return: a simple cryo em simulation of volume 
    """
    from pytom_volume import vol,rotate,shift,initSphere
    from pytom.simulation.support import add_white_noise
    
    if not rotation == [0,0,0]:
        #print '---ROTATE---'
        #print 'EMSimulation simpleSimulation: in rotation 1 ' + str(rotation)
        rotatedCopy = vol(volume.sizeX(),volume.sizeY(),volume.sizeZ())
        rotate(volume,rotatedCopy,rotation[0],rotation[1],rotation[2])
    else:
        rotatedCopy = vol(volume.sizeX(),volume.sizeY(),volume.sizeZ())
        rotatedCopy.copyVolume(volume)
    
    #print 'EMSimulation simpleSimulation: after rotation ' 
    
    if not mask:
        #print 'EMSimulation simpleSimulation: in mask 1' 
        mask = vol(volume.sizeX(),volume.sizeY(),volume.sizeZ())
        initSphere(mask,volume.sizeX()//2-1,0,0, volume.sizeX()//2,
	    volume.sizeX()//2, volume.sizeX()//2)
        maskedCopy = rotatedCopy * mask
    if not mask.__class__ == vol:
        #print 'EMSimulation simpleSimulation: in mask 2'
        
        mask = mask.getVolume(rotation)
        maskedCopy = rotatedCopy * mask
        
    else:
        #print 'EMSimulation simpleSimulation: in mask 3'
        maskedCopy = rotatedCopy * mask
    
    #print 'EMSimulation simpleSimulation:  after mask'
    
    if not shiftV == [0,0,0]:
        #print '--SHIFT---'
        shiftedCopy = vol(volume.sizeX(),volume.sizeY(),volume.sizeZ())
        shift(maskedCopy,shiftedCopy,shiftV[0],shiftV[1],shiftV[2])
    else:
        shiftedCopy = vol(volume.sizeX(),volume.sizeY(),volume.sizeZ())
        shiftedCopy.copyVolume(maskedCopy)
        
    if (shiftV == [0,0,0]) and (rotation==[0,0,0]):
        #no shift and no rotation -> simply take the original volume
        c = vol(maskedCopy.sizeX(),volume.sizeY(),volume.sizeZ())
        c.copyVolume(maskedCopy)
        noisyCopy = add_white_noise(c,SNR)
    else:
        noisyCopy = add_white_noise(shiftedCopy,SNR)
    
    if wedgeInfo:
        #print '---WEDGE---'
        result = wedgeInfo.apply(noisyCopy)
    else:
        result = noisyCopy
    
    #print 'EMSimulation simpleSimulation: end function'
        
    if result.__class__ == list :
        return result[0]
    else:
        return result
