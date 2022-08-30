'''
Created on Jan 27, 2010

@author: Thomas Hrabe, FF
'''
analytWedge=False


def invert_WedgeSum( invol, r_max=None, lowlimit=0., lowval=0.):
    """
    invert wedge sum - avoid division by zero and boost of high frequencies

    @param invol: input volume
    @type invol: L{pytom_volume.vol} or L{pytom_volume.vol_comp}
    @param r_max: radius
    @type r_max: L{int}
    @param lowlimit: lower limit - all values below this value that lie in the specified radius will be replaced \
                by lowval
    @type lowlimit: L{float}
    @param lowval: replacement value
    @type lowval: L{float}

    @author: FF
    """
    from math import sqrt
    if not r_max:
        r_max=invol.sizeY()/2-1

    # full representation with origin in center
    if invol.sizeZ() == invol.sizeX():
        centX1 = int(invol.sizeX()/2)
        centY1 = int(invol.sizeY()/2)
        centZ  = int(invol.sizeZ()/2)
        for ix in range(0,invol.sizeX()):
            for iy in range(0,invol.sizeY()):
                for iz in range(0,invol.sizeZ()):
                    dx = (ix-centX1)**2
                    dy = (iy-centY1)**2
                    dz = (iz-centZ)**2
                    r = sqrt(dx+dy+dz)
                    if r < r_max:
                        v = invol.getV( ix, iy, iz)
                        if v < lowlimit:
                            v = 1./ lowval
                        else:
                            v = 1./ v
                    else:
                        v = 0.
                    invol.setV( v, ix, iy, iz)
    else:
        centX1 = 0
        centX2 = invol.sizeX()-1
        centY1 = 0
        centY2 = invol.sizeY()-1
        centZ  = 0
        for ix in range(0,invol.sizeX()):
            for iy in range(0,invol.sizeY()):
                for iz in range(0,invol.sizeZ()):
                    d1 = (ix-centX1)**2
                    d2 = (ix-centX2)**2
                    dx = min(d1,d2)
                    d1 = (iy-centY1)**2
                    d2 = (iy-centY2)**2
                    dy = min(d1,d2)
                    dz = (iz-centZ)**2
                    r = sqrt(dx+dy+dz)
                    if r < r_max:
                        v = invol.getV( ix, iy, iz)
                        if v < lowlimit:
                            v = 1./lowval
                        else:
                            v = 1./v
                    else:
                        v = 0.
                    invol.setV( v, ix, iy, iz)


def cloneMaximisationResult(maxResult,symmetry):
    """
    cloneMaximisationResult: Clones a Maximization Result according to particles nfold symmetry -> We get n-1 results that will contribute to the later average.
    @param maxResult:
    @param symmetry:   
    """
    
    resultClones = []
    
    if symmetry.isOneFold():
        return resultClones
    
    originalRotation = maxResult.getRotation()
    
    symmetry.setPsi(originalRotation.getPsi())
    symmetry.setTheta(originalRotation.getTheta())
    
    phi = originalRotation.getPhi()
    angleList = symmetry.getAngleList(phi)
    
    rotation = angleList.nextRotation()
    
    while not rotation == [None,None,None]:
        
        result = maxResult.copy()
        result.setRotation(rotation)
        resultClones.append(result)
        
        rotation = angleList.nextRotation()
        
    return resultClones 


def applySymmetryOnParticleList(particleList,symmetry):
    """
    applySymmetryOnParticleList
    @deprecated: use L{pytom.basic.structures.Symmetry.apply} instead!
    """
    
    if symmetry.isOneFold():
        return particleList

    from pytom.basic.structures import ParticleList,Rotation
    newList = ParticleList(particleList.getDirectory())

    for i in range(len(particleList)):
        particle = particleList[i] 
        
        originalRotation = particle.getRotation()
        
        symmetry.setPsi(originalRotation.getPsi())
        symmetry.setTheta(originalRotation.getTheta())
        
        phi = originalRotation.getPhi()
        angleList = symmetry.getAngleList(phi)
    
        rotation = angleList.nextRotation()
        rotation = angleList.nextRotation()
    
        while not rotation == [None,None,None]:
        
            p2 = particle.copy()
        
            p2.setRotation(Rotation(rotation))
            
            newList.append(p2)
        
            rotation = angleList.nextRotation()
        
    return newList + particleList


def applySymmetryToVolume(volume,symmetryObject,wedgeInfo):
    """
    applySymmetryToVolume
    @deprecated: use L{pytom.basic.structures.Symmetry.apply} instead!
    """
    from pytom_volume import read,rotate,shift,vol,initSphere,complexDiv
    from pytom_freqweight import weight
    from pytom.basic.fourier import fft,ifft,ftshift
    from pytom.basic.filter import filter
    from pytom.alignment.structures import ExpectationResult
    from pytom.basic.structures import Reference,Symmetry
    from pytom.tools.maths import epsilon
    
    
    if not volume.__class__ == vol:
        raise Exception('You must provide a volume as first parameter to applySymmetryToObject')
    
    if not symmetryObject.__class__ == Symmetry:
        raise Exception('You must provide a Symmetry object as second parameter to applySymmetryToObject')
    
    result = vol(volume.sizeX(),volume.sizeY(),volume.sizeZ())
    result.copyVolume(volume)
    
    sizeX = volume.sizeX()
    sizeY = volume.sizeY()
    sizeZ = volume.sizeZ()
    
    rot = vol(sizeX,sizeY,sizeZ)
    
    wedgeSum = vol(sizeX,sizeY,sizeZ/2+1)
    wedgeSum.setAll(0.0)
    
    angleList = symmetryObject.getAngleList()
    
    rotation = angleList.nextRotation() 
    
    while not rotation == [None,None,None]:
        
          rotate(volume,rot,rotation[0],rotation[1],rotation[2])          
          result = result + rot
          
          rotation = angleList.nextRotation()
    

    result.shiftscale(0.0,1/float(angleList.numberRotations()))
    return result


def _disrtibuteAverageMPI(particleList,averageName,showProgressBar = False,verbose=False,
                          createInfoVolumes = False,setParticleNodesRatio = 3,sendEndMessage = False):
    """
    _distributeAverageMPI : Distributes averaging to multiple MPI nodes.
    @param particleList: The particles
    @param averageName: Filename of new average 
    @param verbose: Prints particle information. Disabled by default. 
    @param createInfoVolumes: Create info data (wedge sum, inverted density) too? False by default. 
    @return: A new Reference object
    @rtype: L{pytom.basic.structures.Reference}
    @author: Thomas Hrabe
    """
    
    import pytom_mpi
    from pytom.alignment.structures import ExpectationJob
    from pytom.parallel.parallelWorker import ParallelWorker
    from pytom.parallel.alignmentMessages import ExpectationJobMsg
    from pytom_volume import read,complexDiv,complexRealMult
    from pytom.basic.fourier import fft,ifft
    from pytom.basic.filter import lowpassFilter
    from pytom.basic.structures import Reference
    
    import os 
    import sys
    numberOfNodes = pytom_mpi.size()
    
    particleNodesRatio = float(len(particleList)) / float(numberOfNodes)

    splitFactor = numberOfNodes
    
    if particleNodesRatio < setParticleNodesRatio:
        #make sure each node gets at least 20 particles. 
        splitFactor = len(particleList) / setParticleNodesRatio
    
    splitLists = particleList.splitNSublists(splitFactor)

    msgList = []
    avgNameList = []
    preList = []
    wedgeList = []
    
    for i in range(len(splitLists)):
        plist = splitLists[i]
        avgName = averageName + '_dist' +str(i) + '.em'
        avgNameList.append(avgName)
        
        preList.append(averageName + '_dist' +str(i) + '-PreWedge.em')
        wedgeList.append(averageName + '_dist' +str(i) + '-WedgeSumUnscaled.em')
        
        job = ExpectationJob(plist,avgName)
        message = ExpectationJobMsg(0,i)
        message.setJob(job)
        msgList.append(message)
        
    #distribute averaging
    worker = ParallelWorker()
    worker.fillJobList(msgList)
    
    worker.parallelWork(True,sendEndMessage)
    
    #collect results
    result = read(preList[0])
    wedgeSum = read(wedgeList[0])
    
    for i in range(1,len(preList)):
        r = read(preList[i])
        result += r
    
        w = read(wedgeList[i])
        wedgeSum += w 
        
    result.write(averageName[:len(averageName)-3]+'-PreWedge.em')
    wedgeSum.write(averageName[:len(averageName)-3] + '-WedgeSumUnscaled.em')

    invert_WedgeSum( invol=wedgeSum, r_max=result.sizeX()/2-2., lowlimit=.05*len(particleList),
                     lowval=.05*len(particleList))
    fResult = fft(result)
    r = complexRealMult(fResult,wedgeSum)

    result = ifft(r)
    result.shiftscale(0.0,1/float(result.sizeX()*result.sizeY()*result.sizeZ()))

    # do a low pass filter
    result = lowpassFilter(result, result.sizeX()/2-2, (result.sizeX()/2-1)/10.)[0]

    result.write(averageName)
    
    # clean results
    for i in range(0,len(preList)):
        os.system('rm ' + avgNameList[i])
        os.system('rm ' + preList[i])
        os.system('rm ' + wedgeList[i])
        
    return Reference(averageName,particleList)    


def distributeAverage(particleList,averageName,showProgressBar = False,verbose=False,createInfoVolumes = False,sendEndMessage = False):
    """
    distributeAverage : Distributes averaging to multiple nodes
    @param particleList: The particles
    @param averageName: Filename of new average 
    @param verbose: Prints particle information. Disabled by default. 
    @param createInfoVolumes: Create info data (wedge sum, inverted density) too? False by default. 
    @return: A new Reference object
    @rtype: L{pytom.basic.structures.Reference}
    @author: Thomas Hrabe
    """

    import pytom_mpi
    
    mpiInitialized = pytom_mpi.isInitialised()
    mpiAvailable = False
    if not mpiInitialized:
        try:
            pytom_mpi.init()
            
            if pytom_mpi.size() > 1:
                mpiAvailable = True
            
        except:
            print('Could not initialize MPI properly! Running in sequential mode!')
        
    if mpiAvailable:
        if pytom_mpi.rank() == 0:
            return _disrtibuteAverageMPI(particleList,averageName,showProgressBar,verbose,createInfoVolumes,sendEndMessage)
        else:
            from pytom.alignment.ExMaxAlignment import ExMaxWorker
            worker = ExMaxWorker()
            worker.parallelRun(False)
            
    else:
        print('MPI not available')
        return average(particleList,averageName,showProgressBar,verbose,createInfoVolumes)


def average( particleList, averageName, showProgressBar=False, verbose=False,
        createInfoVolumes=False, weighting=False, norm=False, gpuID=None):
    """
    average : Creates new average from a particleList
    @param particleList: The particles
    @param averageName: Filename of new average 
    @param verbose: Prints particle information. Disabled by default. 
    @param createInfoVolumes: Create info data (wedge sum, inverted density) too? False by default.
    @param weighting: apply weighting to each average according to its correlation score
    @param norm: apply normalization for each particle
    @return: A new Reference object
    @rtype: L{pytom.basic.structures.Reference}
    @author: Thomas Hrabe
    @change: limit for wedgeSum set to 1% or particles to avoid division by small numbers - FF
    """
    from pytom_volume import read,vol,reducedToFull,limit, complexRealMult
    from pytom.basic.filter import lowpassFilter, rotateWeighting
    from pytom_volume import transformSpline as transform
    from pytom.basic.fourier import convolute
    from pytom.basic.structures import Reference
    from pytom.basic.normalise import mean0std1
    from pytom.tools.ProgressBar import FixedProgBar
    from math import exp
    import os

    if len(particleList) == 0:
        raise RuntimeError('The particle list is empty. Aborting!')
    
    if showProgressBar:
        progressBar = FixedProgBar(0,len(particleList),'Particles averaged ')
        progressBar.update(0)
        numberAlignedParticles = 0
    
    result = []
    wedgeSum = []
    
    newParticle = None
    # pre-check that scores != 0
    if weighting:
        wsum = 0.
        for particleObject in particleList:
            wsum += particleObject.getScore().getValue()
        if wsum < 0.00001:
            weighting = False
            print("Warning: all scores have been zero - weighting not applied")

    for particleObject in particleList:
        
        if verbose:
            print(particleObject)

        if not os.path.exists(particleObject.getFilename()): continue
        particle = read(particleObject.getFilename())
        if norm: # normalize the particle
            mean0std1(particle) # happen inplace
        
        wedgeInfo = particleObject.getWedge()
        # apply its wedge to itself
        particle = wedgeInfo.apply(particle)
        
        if result == []:
            sizeX = particle.sizeX() 
            sizeY = particle.sizeY()
            sizeZ = particle.sizeZ()
            
            newParticle = vol(sizeX,sizeY,sizeZ)
            
            centerX = sizeX/2 
            centerY = sizeY/2 
            centerZ = sizeZ/2 
            
            result = vol(sizeX,sizeY,sizeZ)
            result.setAll(0.0)
            if analytWedge:
                wedgeSum = wedgeInfo.returnWedgeVolume(wedgeSizeX=sizeX, wedgeSizeY=sizeY, wedgeSizeZ=sizeZ)
            else:
                # > FF bugfix
                wedgeSum = wedgeInfo.returnWedgeVolume(sizeX,sizeY,sizeZ)
                # < FF
                # > TH bugfix
                #wedgeSum = vol(sizeX,sizeY,sizeZ)
                # < TH
                #wedgeSum.setAll(0)
            assert wedgeSum.sizeX() == sizeX and wedgeSum.sizeY() == sizeY and wedgeSum.sizeZ() == sizeZ/2+1, \
                    "wedge initialization result in wrong dims :("
            wedgeSum.setAll(0)

        ### create spectral wedge weighting
        rotation = particleObject.getRotation()
        rotinvert = rotation.invert()
        if analytWedge:
            # > analytical buggy version
            wedge = wedgeInfo.returnWedgeVolume(sizeX,sizeY,sizeZ,False, rotinvert)
        else:
            # > FF: interpol bugfix
            wedge = rotateWeighting( weighting=wedgeInfo.returnWedgeVolume(sizeX,sizeY,sizeZ,False),
                                     z1=rotinvert[0], z2=rotinvert[1], x=rotinvert[2], mask=None,
                                     isReducedComplex=True, returnReducedComplex=True)
            # < FF
            # > TH bugfix
            #wedgeVolume = wedgeInfo.returnWedgeVolume(wedgeSizeX=sizeX, wedgeSizeY=sizeY, wedgeSizeZ=sizeZ,
            #                                    humanUnderstandable=True, rotation=rotinvert)
            #wedge = rotate(volume=wedgeVolume, rotation=rotinvert, imethod='linear')
            # < TH

        ### shift and rotate particle
        shiftV = particleObject.getShift()
        newParticle.setAll(0)
            
        transform(particle,newParticle,-rotation[1],-rotation[0],-rotation[2],
                  centerX,centerY,centerZ,-shiftV[0],-shiftV[1],-shiftV[2],0,0,0)
        
        if weighting:
            weight = 1.-particleObject.getScore().getValue()
            #weight = weight**2
            weight = exp(-1.*weight)
            result = result + newParticle * weight
            wedgeSum = wedgeSum + wedge * weight
        else:
            result = result + newParticle
            wedgeSum = wedgeSum + wedge
        
        if showProgressBar:
            numberAlignedParticles = numberAlignedParticles + 1
            progressBar.update(numberAlignedParticles)

    ###apply spectral weighting to sum
    result = lowpassFilter(result, sizeX//2-1, 0.)[0]
    #if createInfoVolumes:
    result.write(averageName[:len(averageName)-3]+'-PreWedge.em')
    wedgeSum.write(averageName[:len(averageName)-3] + '-WedgeSumUnscaled.em')
        
    invert_WedgeSum( invol=wedgeSum, r_max=sizeX/2-2., lowlimit=.05*len(particleList), lowval=.05*len(particleList))
    
    if createInfoVolumes:
        wedgeSum.write(averageName[:len(averageName)-3] + '-WedgeSumInverted.em')
        
    result = convolute(v=result, k=wedgeSum, kernel_in_fourier=True)

    # do a low pass filter
    #result = lowpassFilter(result, sizeX/2-2, (sizeX/2-1)/10.)[0]
    result.write(averageName)
    
    if createInfoVolumes:
        resultINV = result * -1
        #write sign inverted result to disk (good for chimera viewing ... )
        resultINV.write(averageName[:len(averageName)-3]+'-INV.em')
    newReference = Reference(averageName,particleList)
    
    return newReference

def average2(particleList, weighting=False, norm=False, determine_resolution=False,
             mask=None, binning=1, verbose=False):
    """
    2nd version of average function. Will not write the averages to the disk. Also support internal \
    resolution determination.
    """
    from pytom_volume import read, vol, complexDiv, complexRealMult
    from pytom_volume import transformSpline as transform
    from pytom.basic.fourier import fft, ifft, convolute
    from pytom.basic.normalise import mean0std1
    from pytom.tools.ProgressBar import FixedProgBar
    from pytom.basic.filter import lowpassFilter, rotateWeighting
    from math import exp
    
    if len(particleList) == 0:
        raise RuntimeError('The particlelist provided is empty. Aborting!')
    
    if verbose:
        progressBar = FixedProgBar(0,len(particleList),'Particles averaged ')
        progressBar.update(0)
        numberAlignedParticles = 0
    
    even = None
    odd = None
    wedgeSum_even = None
    wedgeSum_odd = None
    newParticle = None
    
    is_odd = True
    for particleObject in particleList:
        particle = read(particleObject.getFilename(), 0,0,0,0,0,0,0,0,0, binning,binning,binning)
        if norm:
            mean0std1(particle)
        wedgeInfo = particleObject.getWedge()
        
        # apply its wedge to itself
        particle = wedgeInfo.apply(particle)
        
        if odd is None: # initialization
            sizeX = particle.sizeX() 
            sizeY = particle.sizeY()
            sizeZ = particle.sizeZ()
            
            newParticle = vol(sizeX,sizeY,sizeZ)
            
            centerX = sizeX/2 
            centerY = sizeY/2 
            centerZ = sizeZ/2 
            
            odd = vol(sizeX,sizeY,sizeZ)
            odd.setAll(0.0)
            even = vol(sizeX,sizeY,sizeZ)
            even.setAll(0.0)
            
            wedgeSum_odd = wedgeInfo.returnWedgeVolume(sizeX,sizeY,sizeZ)
            wedgeSum_odd.setAll(0)
            wedgeSum_even = wedgeInfo.returnWedgeVolume(sizeX,sizeY,sizeZ)
            wedgeSum_even.setAll(0)
        

        # create spectral wedge weighting
        rotation = particleObject.getRotation()
        rotinvert =  rotation.invert()
        if analytWedge:
            # > original buggy version
            wedge = wedgeInfo.returnWedgeVolume(sizeX,sizeY,sizeZ,False, rotinvert)
            # < original buggy version
        else:
            # > FF: interpol bugfix
            wedge = rotateWeighting( weighting=wedgeInfo.returnWedgeVolume(sizeX,sizeY,sizeZ,False),
                                     z1=rotinvert[0], z2=rotinvert[1], x=rotinvert[2], mask=None,
                                     isReducedComplex=True, returnReducedComplex=True)
            # < FF
            # > TH bugfix
            #wedgeVolume = wedgeInfo.returnWedgeVolume(wedgeSizeX=sizeX, wedgeSizeY=sizeY, wedgeSizeZ=sizeZ,
            #                                          humanUnderstandable=True, rotation=rotinvert)
            #wedge = rotate(volume=wedgeVolume, rotation=rotinvert, imethod='linear')
            # < TH
        if is_odd:
            wedgeSum_odd = wedgeSum_odd + wedge
        else:
            wedgeSum_even = wedgeSum_even + wedge
        
        # shift and rotate particle
        shiftV = particleObject.getShift()
        newParticle.setAll(0)
        transform(particle,newParticle,-rotation[1],-rotation[0],-rotation[2],
                  centerX,centerY,centerZ,-shiftV[0]/binning,
                  -shiftV[1]/binning,-shiftV[2]/binning,0,0,0)

        if is_odd:
            if weighting:
                weight = 1. - particleObject.getScore().getValue()
                #weight = weight**2
                weight = exp(-1.*weight)
                odd = odd + newParticle * weight
            else:
                odd = odd + newParticle
        else:
            if weighting:
                weight = 1. - particleObject.getScore().getValue()
                #weight = weight**2
                weight = exp(-1.*weight)
                even = even + newParticle * weight
            else:
                even = even + newParticle
        
        is_odd = not is_odd
        
        if verbose:
            numberAlignedParticles = numberAlignedParticles + 1
            progressBar.update(numberAlignedParticles)

    # determine resolution if needed
    fsc = None
    if determine_resolution:
        # apply spectral weighting to sum
        f_even = fft(even)
        w_even = complexDiv(f_even, wedgeSum_even)
        w_even = ifft(w_even)        
        w_even.shiftscale(0.0,1/float(sizeX*sizeY*sizeZ))
        
        f_odd = fft(odd)
        w_odd = complexDiv(f_odd, wedgeSum_odd)
        w_odd = ifft(w_odd)        
        w_odd.shiftscale(0.0,1/float(sizeX*sizeY*sizeZ))
        
        from pytom.basic.correlation import FSC
        fsc = FSC(w_even, w_odd, sizeX/2, mask, verbose=False)
    
    # add together
    result = even+odd
    wedgeSum = wedgeSum_even+wedgeSum_odd

    invert_WedgeSum( invol=wedgeSum, r_max=sizeX/2-2., lowlimit=.05*len(particleList), lowval=.05*len(particleList))
    #wedgeSum.write(averageName[:len(averageName)-3] + '-WedgeSumInverted.em')
    result = convolute(v=result, k=wedgeSum, kernel_in_fourier=True)
    # do a low pass filter
    #result = lowpassFilter(result, sizeX/2-2, (sizeX/2-1)/10.)[0]
    
    return (result, fsc)


def alignTwoVolumes(particle,reference,angleObject,mask,score,preprocessing,progressBar=False):
    """
    alignTwoVolumes: align two volumes with respect to each other
    @param particle: volume to be aligned
    @param reference: reference volume (not translated and rotated)
    @param angleObject: angular sampling
    @param mask: mask volume
    @param score: score type
    @param preprocessing: preprocessing parameters
    @param progressBar: show progress bar
    @type progressBar: bool
    """
    
    from pytom.alignment.structures import GrowingAverageInterimResult
    from pytom.basic.structures import Rotation,Shift
    
    partVolume = particle.getVolume()
    
    refVolume = reference.getVolume() 

    refWeight = reference.getWeighting()

    peak = bestAlignment(partVolume,refVolume,refWeight,particle.getWedge(),
            angleObject,score,mask,preprocessing)
    
    score.setValue(peak.getScoreValue())
    
    return GrowingAverageInterimResult(particle,reference,Rotation(peak.getRotation()),Shift(peak.getShift()),score)


def _rotateWedgeReference(reference,rotation,wedgeInfo,mask,rotationCenter):
    """
    _rotateShiftWedgeParticle: Wrapper for Rotation, Shift and WedgeWeighting of Reference
    @param reference: The reference
    @param rotation: The rotation 
    @type rotation: L{pytom.basic.structures.Rotation}
    @param wedgeInfo: Wedge info object
    @type wedgeInfo: L{pytom.basic.structures.Wedge}
    @param mask: The mask object (a volume) or None
    @type mask: L{pytom_volume.vol}
    @return:
    @change: support mask == None, FF
    """
    from pytom_volume import vol, transform as transform #developers: your can also import transformSpline for more accurate rotation!
    
    rotatedVolume = vol(reference.sizeX(),reference.sizeY(),reference.sizeZ())
    transform(reference,rotatedVolume,rotation[0],rotation[1],rotation[2],rotationCenter[0],rotationCenter[1],rotationCenter[2],0,0,0,0,0,0)
   
    if mask:
        return wedgeInfo.apply(rotatedVolume) * mask
    else:
        return wedgeInfo.apply(rotatedVolume)


def bestAlignment(particle, reference, referenceWeighting, wedgeInfo, rotations,
         scoreObject=0, mask=None, preprocessing=None, progressBar=False, binning=1,
         bestPeak=None, verbose=False):
    """
    bestAlignment: Determines best alignment of particle relative to the reference
    @param particle: A particle
    @type particle: L{pytom_volume.vol}
    @param reference: A reference
    @type reference: L{pytom_volume.vol}
    @param referenceWeighting: Fourier weighting of the reference (sum of wedges for instance)
    @type referenceWeighting: L{pytom.basic.structures.vol}
    @param wedgeInfo: What does the wedge look alike?
    @type wedgeInfo: L{pytom.basic.structures.Wedge}
    @param rotations: All rotations to be scanned
    @type rotations: L{pytom.angles.AngleObject}
    @param scoreObject: 
    @type scoreObject: L{pytom.score.score.Score}
    @param mask: real-space mask for correlation function
    @type mask: L{pytom.basic.structures.Particle}
    @param preprocessing: Class storing preprocessing of particle and reference such as bandpass
    @type preprocessing: L{pytom.alignment.preprocessing.Preprocessing}
    @param progressBar: Display progress bar of alignment. False by default.
    @param binning: binning factor - e.g. binning=2 reduces size by FACTOR of 2
    @type binning: int or float
    @param bestPeak: Initialise best peak with old values.   
    @param verbose: Print out infos. Writes CC volume to disk!!! Default is False  
    @return: Returns the best rotation for particle and the corresponding scoring result.
    @author: Thomas Hrabe
    """
    from pytom.basic.correlation import subPixelPeak, subPixelPeakParabolic
    from pytom.alignment.structures import Peak
    from pytom_volume import peak, vol, vol_comp
    from pytom.basic.filter import filter,rotateWeighting
    from pytom.basic.structures import Rotation, Shift, Particle, Mask
    from pytom.angles.angle import AngleObject
    from pytom.alignment.preprocessing import Preprocessing
    from pytom.basic.transformations import resize, resizeFourier
    binningType = 'Fourier' # or 'Fourier'

    assert isinstance(rotations, AngleObject), "bestAlignment: rotations must be " \
                                                                             "AngleObject!"
    currentRotation = rotations.nextRotation()
    if currentRotation == [None,None,None]:
        raise Exception('bestAlignment: No rotations are sampled! Something is wrong with input rotations')

    assert particle.__class__ == vol, "particle not of type vol"
    assert reference.__class__ == vol, "reference not of type vol"
    assert (referenceWeighting.__class__ == vol or referenceWeighting.__class__ == str), \
        "referenceWeighting not volume or str"
    if mask:
        assert mask.__class__ == Mask, "Mask not of type Mask"
        m = mask.getVolume()

    if scoreObject == 0 or not scoreObject:
        from pytom.basic.score import xcfScore
        scoreObject = xcfScore()
    # fix binning
    if binning == 0:
        binning = 1
    if binning != 1:
        particleUnbinned = vol(particle.sizeX(), particle.sizeY(), particle.sizeZ())
        particleUnbinned.copyVolume(particle)
        particle = resize(volume=particle, factor=1./binning, interpolation=binningType)
        if type(particle) == tuple:
            particle = particle[0]
        referenceUnbinned = vol(reference.sizeX(), reference.sizeY(), reference.sizeZ())
        referenceUnbinned.copyVolume(reference)
        reference = resize(volume=reference, factor=1./binning, interpolation=binningType)
        if type(reference) == tuple:
            reference = reference[0]
        if mask:
            m = resize(volume=m, factor=1./binning, interpolation='Spline')
        if not referenceWeighting.__class__ == str:
            referenceWeightingUnbinned = vol_comp(referenceWeighting.sizeX(), referenceWeighting.sizeY(),
                                                  referenceWeighting.sizeZ())
            referenceWeightingUnbinned.copyVolume(referenceWeighting)
            if binning != 1:
                referenceWeighting = resizeFourier(fvol=referenceWeighting, factor=1./binning)
    centerX, centerY, centerZ = int(particle.sizeX()/2), int(particle.sizeY()/2), int(particle.sizeZ()/2)

    # create buffer volume for transformed particle 
    particleCopy = vol(particle.sizeX(),particle.sizeY(),particle.sizeZ())
    particle = wedgeInfo.apply(particle) #apply wedge to itself
    if preprocessing is None:
        preprocessing = Preprocessing()
    preprocessing.setTaper( taper=particle.sizeX()/10.)
    particle = preprocessing.apply(volume=particle, bypassFlag=True)  # filter particle to some resolution
    particleCopy.copyVolume(particle)

    if mask:
        from pytom_volume import sum
        from pytom.basic.correlation import meanUnderMask, stdUnderMask
        p = sum(m)
        meanV = meanUnderMask(particle, m, p)
        stdV = stdUnderMask(particle, m, p, meanV)
    else:
        meanV = None
        stdV = None
    cntr = 0
    while currentRotation != [None,None,None]:
        if mask:
            m = mask.getVolume(currentRotation)
            if binning != 1:
                m = resize(volume=m, factor=1./binning, interpolation='Spline')
            #update stdV if mask is not a sphere
            # compute standard deviation volume really only if needed
            # print('Mask is sphere: ', mask.isSphere(), scoreObject._type)

            if (not mask.isSphere()) and (scoreObject._type=='FLCFScore'):
                if 1: print('recalc meanV en stdV')
                meanV   = meanUnderMask(particle, m, p)
                stdV    = stdUnderMask(particle, m, p, meanV)
        else:
            m = None
        
        simulatedVol = _rotateWedgeReference(reference, currentRotation, wedgeInfo, m, [centerX, centerY, centerZ])
        simulatedVol = preprocessing.apply(volume=simulatedVol, bypassFlag=True)
        
        #weight particle
        if not referenceWeighting.__class__ == str:
            from pytom_freqweight import weight
            weightingRotated = rotateWeighting(weighting=referenceWeighting, z1=currentRotation[0],
                                               z2=currentRotation[1], x=currentRotation[2], isReducedComplex=True,
                                               returnReducedComplex=True, binarize=False)
            particleCopy.copyVolume(particle)
            r = list(filter(particleCopy, weight(weightingRotated)))
            particleCopy = r[0]


        scoringResult = scoreObject.score(particleCopy, simulatedVol, m, stdV)

        pk = peak(scoringResult)

        # with subPixelPeak
        [peakValue,peakPosition] = subPixelPeak(scoreVolume=scoringResult, coordinates=pk,
                                               interpolation='Quadratic', verbose=False)
        #[peakValue,peakPosition] = subPixelPeakParabolic(scoreVolume=scoringResult, coordinates=pk, verbose=False)

        # determine shift relative to center
        shiftX = (peakPosition[0] - centerX) * binning
        shiftY = (peakPosition[1] - centerY) * binning
        shiftZ = (peakPosition[2] - centerZ) * binning

        #NANs would fail this test.
        assert peakValue == peakValue, "peakValue seems to be NaN"
        
        newPeak = Peak(peakValue, Rotation(currentRotation), Shift(shiftX, shiftY, shiftZ))
        
        if verbose:
            print('Rotation: z1=%3.1f, z2=%3.1f, x=%3.1f; Dx=%2.2f, Dy=%2.2f, Dz=%2.2f, CC=%2.3f' % \
                  (currentRotation[0], currentRotation[1], currentRotation[2], shiftX, shiftY, shiftZ, peakValue))

        if bestPeak is None:
            bestPeak = newPeak

            if verbose:
                scoringResult.write('BestScore.em')
        if bestPeak < newPeak:
            bestPeak = newPeak

            if verbose:
                scoringResult.write('BestScore.em')

        currentRotation = rotations.nextRotation()

    # repeat ccf for binned sampling to get translation accurately
    if binning != 1:
        m = mask.getVolume(bestPeak.getRotation())
        centerX, centerY, centerZ = int(particleUnbinned.sizeX()/2), int(particleUnbinned.sizeY()/2), \
                                    int(particleUnbinned.sizeZ()/2)
        simulatedVol = _rotateWedgeReference(referenceUnbinned, bestPeak.getRotation(), wedgeInfo, m,
                                             [centerX, centerY, centerZ])
        simulatedVol = preprocessing.apply(volume=simulatedVol, bypassFlag=True)
        if mask and scoreObject._type=='FLCFScore':
            p = sum(m)
            meanV = meanUnderMask(volume=particleUnbinned, mask=m, p=p)
            stdV  = stdUnderMask(volume=particleUnbinned, mask=m, p=p, meanV=meanV)
        scoreObject._peakPrior.reset_weight()
        scoringResult = scoreObject.score(particle=particleUnbinned, reference=simulatedVol, mask=m, stdV=stdV)
        pk = peak(scoringResult)
        [peakValue,peakPosition] = subPixelPeak(scoreVolume=scoringResult, coordinates=pk, interpolation='Quadratic',
                                                verbose=False)
        shiftX = (peakPosition[0] - centerX)
        shiftY = (peakPosition[1] - centerY)
        shiftZ = (peakPosition[2] - centerZ)
        bestPeak = Peak(peakValue, bestPeak.getRotation(), Shift(shiftX, shiftY, shiftZ))
    if verbose:
            print('BestAlignment: z1=%3.1f, z2=%3.1f, x=%3.1f; Dx=%2.2f, Dy=%2.2f, Dz=%2.2f, CC=%2.3f' % \
                  (bestPeak.getRotation()[0], bestPeak.getRotation()[1], bestPeak.getRotation()[2],
                   bestPeak.getShift()[0], bestPeak.getShift()[1], bestPeak.getShift()[2], bestPeak.getScoreValue()))
    rotations.reset()
    scoreObject._peakPrior.reset_weight()
    return bestPeak


def bestAlignmentGPU(particle, rotations, plan, preprocessing=None, wedgeInfo=None, isSphere=True,
                     rotation_order='rzxz', max_shift=40, profile=False, interpolation_factor=0.1):
    """
    bestAlignment: Determines best alignment of particle relative to the reference
    @param particle: A particle
    @type particle: L{pytom_volume.vol}
    @param reference: A reference
    @type reference: L{pytom_volume.vol}
    @param referenceWeighting: Fourier weighting of the reference (sum of wedges for instance)
    @type referenceWeighting: L{pytom.basic.structures.vol}
    @param wedgeInfo: What does the wedge look alike?
    @type wedgeInfo: L{pytom.basic.structures.Wedge}
    @param rotations: All rotations to be scanned
    @type rotations: L{pytom.angles.AngleObject}
    @param scoreObject:
    @type scoreObject: L{pytom.score.score.Score}
    @param mask: real-space mask for correlation function
    @type mask: L{pytom.basic.structures.Particle}
    @param preprocessing: Class storing preprocessing of particle and reference such as bandpass
    @type preprocessing: L{pytom.alignment.preprocessing.Preprocessing}
    @param progressBar: Display progress bar of alignment. False by default.
    @param binning: binning factor - e.g. binning=2 reduces size by FACTOR of 2
    @type binning: int or float
    @param bestPeak: Initialise best peak with old values.
    @param verbose: Print out infos. Writes CC volume to disk!!! Default is False
    @return: Returns the best rotation for particle and the corresponding scoring result.
    @author: Thomas Hrabe
    """
    from pytom.basic.structures import Rotation, Shift
    from pytom.alignment.structures import Peak

    centerCoordinates = tuple([size//2 for size in plan.volume.shape])

    # create buffer volume for transformed particle
    plan.volume = plan.cp.array(particle, dtype=plan.cp.float32) * plan.taperMask
    plan.updateWedge(wedgeInfo)
    plan.wedgeParticle()
    plan.calc_stdV()
    currentRotation = rotations.nextRotation()
    if currentRotation == [None, None, None]:
        raise Exception('bestAlignment: No rotations are sampled! Something is wrong with input rotations')

    bestPeak = Peak(float(-100000.), Rotation(currentRotation), Shift([0,0,0]))

    while currentRotation != [None, None, None]:
        plan.rotatedRef *= 0
        # If not spherical mask, recalculate respective arrays
        if not isSphere:
            plan.mask.transform(output=plan.rotatedMask,
                                rotation=(currentRotation[0], currentRotation[2], currentRotation[1]),
                                center=centerCoordinates, rotation_order=rotation_order)
            plan.p = plan.rotatedMask.sum()
            plan.mask_fft = plan.fftnP(plan.rotatedMask.astype(plan.cp.complex64),plan=plan.fftplan)
            plan.calc_stdV()

        # Rotate reference
        plan.referenceTex.transform(rotation=(float(currentRotation[0]), float(currentRotation[2]),
                                              float(currentRotation[1])),
                                    center=centerCoordinates, rotation_order=rotation_order, output=plan.rotatedRef)

        # Apply wedge filter to rotated reference
        plan.wedgeRotatedRef()

        # Normalize rotated wedged reference
        plan.normalizeVolume()

        # Calculate normalized crosscorrelation
        plan.cross_correlation()

        # Find subPixelPeak
        # [peakValue, peakShifts] = plan.subPixelMax3D(ignore_border=border, k=0.1, profile=profile)
        peakValue, peakShifts = plan.subPixelMaxSpline()
        newPeak = Peak(float(peakValue), Rotation(currentRotation), Shift(*peakShifts))

        # Save peak if correlation is better
        if newPeak > bestPeak:
            bestPeak = newPeak

        # Update current rotation
        currentRotation = rotations.nextRotation()

    # Update respective averages, saves time as one does not have to allocate particles to gpu again.
    # plan.addParticleAndWedgeToSum(particle, bestPeak, centerCoordinates)

    # TODO GPU should also repeat ccf for unbinned sample at this point
    return bestPeak


def compareTwoVolumes(particle,reference,referenceWeighting,wedgeInfo,rotations,shift,
         scoreObject=0,mask=None,preprocessing=[],binning=1,verbose=False):
    """
    compare: Compares two volumes
    @param particle: A particle
    @param reference: A reference
    @param referenceWeighting: Fourier weighting of the reference (sum of wedges for instance)
    @param wedgeInfo: What does the wedge look alike?
    @type wedgeInfo: L{pytom.basic.structures.WedgeInfo}
    @param rotations: All rotations to be scanned
    @type rotations: Must be a L{pytom.angles.angleList.OneAngleList}
    @param scoreObject: 
    @type scoreObject: L{pytom.score.score.Score}
    @param mask: 
    @param preprocessing: Class storing preprocessing of particle and reference such as bandpass
    @type preprocessing: L{pytom.alignment.preprocessing.Preprocessing}
    @param binning: Is binning applied (Disabled when == 1)
    @return: Returns the correlation coefficient of two volumes.
    @author: Thomas Hrabe
    """

    from pytom_volume import vol,transformSpline
    from pytom.basic.filter import filter,rotateWeighting
    from pytom.angles.angleList import OneAngleList
    
    assert particle.__class__ == vol, "particle not of type vol"
    assert reference.__class__ == vol, "reference not of type vol"
    assert referenceWeighting.__class__ == vol or referenceWeighting.__class__ == str
    assert rotations.__class__ == OneAngleList or len(rotations) == 1
    
    if scoreObject == 0:
        from pytom.basic.score import nxcfScore
        scoreObject = nxcfScore()
    
    particleCopy = vol(particle.sizeX(),particle.sizeY(),particle.sizeZ())
    
    #process particle   
    particle = wedgeInfo.apply(particle)
    particle = preprocessing.apply(particle,True)
    particleCopy.copyVolume(particle)
    
    #manipulate reference
    currentRotation = rotations.nextRotation()
    
    sizeX = particle.sizeX() 
    sizeY = particle.sizeY()
    sizeZ = particle.sizeZ()
    
    centerX = sizeX/2.0 
    centerY = sizeY/2.0 
    centerZ = sizeZ/2.0 
    
    
    simulatedVol= vol(sizeX,sizeY,sizeZ)

    transformSpline(reference,simulatedVol,currentRotation[0],currentRotation[1],currentRotation[2],
                    centerX,centerY,centerZ,0,0,0,shift[0]/binning,shift[1]/binning,shift[2]/binning)
    
    simulatedVol = wedgeInfo.apply(simulatedVol)
    m = mask.getVolume(currentRotation)
    simulatedVol = simulatedVol * m
    
    simulatedVol = preprocessing.apply(simulatedVol,True)

    if not referenceWeighting.__class__ == str:
        from pytom_freqweight import weight
        
        weightingRotated = rotateWeighting(weighting=referenceWeighting, z1=currentRotation[0], z2=currentRotation[1],
                                           x=currentRotation[2], isReducedComplex=None, returnReducedComplex=True,
                                           binarize=False)
        if verbose:
            particleCopy.info('pc')
            weightingRotated.info('wr')
        
        r = list(filter(particleCopy,weight(weightingRotated)))            
        particleCopy = r[0]
        
    scoringResult = scoreObject.scoringCoefficient(particleCopy,simulatedVol,m)
    
    if scoringResult.__class__ == list:
        scoringResult = scoringResult[0]
    
    
    assert scoringResult == scoringResult, "scoringResult possibly NaN"
    
    return scoringResult 


def alignmentProgressToHTML(filename,iterationNumber,score,resolution,angularIncrement,angleDistance,shiftDistance):
    """
    alignmentProgressToHTML
    """
    from pytom.tools.files import getPytomPath,readStringFile  

    if angleDistance == 0:
        angleDistance= [0,0]
        shiftDistance= [0,0]
    
    html = '<html>\n'
    html = html + '<title>Angular Progress for iteration ' + str(iterationNumber) + '</title>\n'
    html = html + '<body><center><table>\n'
    html = html + '<tr><td>Sum of scores</td><td>' + str(score) + ' </td></tr>\n'
    html = html + '<tr><td>Resolution in angstrom</td><td>' + str(resolution) + ' </td></tr>\n' 
    html = html + '<tr><td>AngularIncrement for next iteration </td><td>' + str(angularIncrement)  + '\n</td></tr>'
    html = html + '<tr><td>Mean angular distance after alignment</td><td>Mean: ' + str(angleDistance[0])  + '</td><td>Std: ' + str(angleDistance[1]) + '\n</td></tr>'
    html = html + '<tr><td>Mean shift distance after alignment </td><td>Mean: ' + str(shiftDistance[0])  + '</td><td>Std: ' + str(shiftDistance[1]) + '\n</td></tr>'
    html = html + '</table></center></body></html>'
    
    
    file = open(filename, "w")
    file.write(html)
    file.close()


def FRMAlignmentWrapper(particle,wedgeParticle, reference, wedgeReference,bandwidth, highestFrequency, mask = None , peakPrior = None):
    """
    FRMAlignmentWrapper: Wrapper for frm_align to handle PyTom objects.
    @param particle: The particle 
    @type particle: L{pytom.basic.structures.Particle}
    @param wedgeParticle: Wedge object of particle
    @type wedgeParticle: L{pytom.basic.structures.Wedge} 
    @param reference: Reference used for alignment
    @type reference: L{pytom.basic.structures.Reference}
    @param wedgeReference: Information about reference wedge 
    @type wedgeReference: L{pytom.basic.structures.Wedge}
    @param bandwidth: The bandwidth of the spherical harmonics - lowestBand used, highestBand used 
    @type bandwidth: [lowestBand,highestBand]
    @param highestFrequency: Highest frequency for lowpass filter in fourierspace
    @type highestFrequency: int
    @param mask: Mask that is applied to the particle
    @type mask: L{pytom.basic.structures.Mask}
    @param peakPrior: Maximum distance of peak from origin
    @type peakPrior: L{pytom.score.score.PeakPrior} or an integer
    @return: Returns a list of [L{pytom.basic.structures.Shift}, L{pytom.basic.structures.Rotation}, scoreValue]
    """
    from pytom.basic.structures import Particle,Reference,Mask,Wedge,Shift,Rotation
    from pytom.basic.score import PeakPrior
    from sh_alignment.frm import frm_align
    
    if particle.__class__ == Particle:
        particle = particle.getVolume()
    
    if reference.__class__ == Reference:
        reference = reference.getVolume()
    
    if mask.__class__ == Mask:
        mask = mask.getVolume()
    else:
        mask = None
        
    if bandwidth.__class__ == list and len(bandwidth) != 2:
        raise RuntimeError('Bandwidth parameter must be a list of two integers!')

    if peakPrior and peakPrior.__class__ == PeakPrior:
        peakPrior = int(peakPrior.getRadius())
    
    if not peakPrior or peakPrior < 0.0001:
        peakPrior = None
        
    pos, angle, score = frm_align(particle, wedgeParticle.getWedgeObject(), reference*mask, wedgeReference.getWedgeObject(), bandwidth, int(highestFrequency), peakPrior)
            
    return [Shift([pos[0]-particle.sizeX()/2, pos[1]-particle.sizeY()/2, pos[2]-particle.sizeZ()/2]), Rotation(angle), score]

