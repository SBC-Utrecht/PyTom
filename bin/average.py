#!/usr/bin/env pytom
from pytom.basic.structures import PyTomClass
from pytom.basic.structures import ParticleList
from pytom.angles.localSampling import LocalSampling
from tompy.mpi import MPI
import os
import sys
analytWedge=False


def splitParticleList(particleList, setParticleNodesRatio=3, numberOfNodes=10):
    """
    @param particleList: The particle list
    @param setParticleNodesRatio: minimum number of particles per node
    @type setParticleNodesRatio: L{int}
    @return: list of particle lists, splitFactor (number of processors or smaller for few particles)
    @rtype: list, L{int}
    @author: FF
    """

    particleNodesRatio = float(len(particleList)) / float(numberOfNodes)
    splitFactor = numberOfNodes
    #make sure each node gets at least setParticleNodesRatio particles.
    if particleNodesRatio < setParticleNodesRatio:
        splitFactor = len(particleList) / int(setParticleNodesRatio)
    splitLists = particleList.splitNSublists(splitFactor)  # somehow ...
    return splitLists

def average( particleList, averageName, showProgressBar=False, verbose=False,
        createInfoVolumes=False, weighting=False, norm=False):
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
        print('go', result, not result)
        if 0 and verbose:
            print(particleObject)

    
        if not os.path.exists(particleObject.getFilename()):
            continue
        particle = read(particleObject.getFilename())
        if norm: # normalize the particle
            mean0std1(particle) # happen inplace
        
        wedgeInfo = particleObject.getWedge()
        # apply its wedge to itself
        particle = wedgeInfo.apply(particle)

        if not result:
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
        rotinvert =  rotation.invert()
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
        print(sizeX)
    ###apply spectral weighting to sum

    result = lowpassFilter(result, sizeX/2-1, 0.)[0]
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
                            v = 1./lowval
                        else:
                            v = 1./v
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

def averageParallel(particleList,averageName, showProgressBar=False, verbose=False,
                    createInfoVolumes=False, weighting=None, norm=False,
                    setParticleNodesRatio=3,cores=6):
    """
    compute average using parfor
    @param particleList: The particles
    @param averageName: Filename of new average
    @param verbose: Prints particle information. Disabled by default.
    @param createInfoVolumes: Create info data (wedge sum, inverted density) too? False by default.
    @param weighting: weight particles by exp CC in average
    @type weighting: bool
    @param setParticleNodesRatio: minimum number of particles per node
    @type setParticleNodesRatio: L{int}
    @return: A new Reference object
    @rtype: L{pytom.basic.structures.Reference}
    @author: FF

    """
    from pytom_volume import read, complexRealMult
    from pytom.basic.fourier import fft,ifft
    from pytom.basic.filter import lowpassFilter
    from pytom.basic.structures import Reference
    from pytom.alignment.alignmentFunctions import invert_WedgeSum

    import os

    splitLists = splitParticleList(particleList, setParticleNodesRatio=setParticleNodesRatio, numberOfNodes=cores)
    splitFactor = len(splitLists)

    avgNameList = []
    preList = []
    wedgeList = []
    for ii in range(splitFactor):
        avgName = averageName + '_dist' + str(ii) + '.em'
        avgNameList.append(avgName)
        preList.append(averageName + '_dist' + str(ii) + '-PreWedge.em')
        wedgeList.append(averageName + '_dist' + str(ii) + '-WedgeSumUnscaled.em')

    #reference = average(particleList=plist, averageName=xxx, showProgressBar=True, verbose=False,
    # createInfoVolumes=False, weighting=weighting, norm=False)
    from multiprocessing import Process

    procs = []
    for i in range(splitFactor):
        proc = Process(target=average,args=(splitLists[i],avgNameList[i], showProgressBar,verbose,createInfoVolumes, weighting, norm) )
        procs.append(proc)
        proc.start()

    import time
    while procs:
        procs = [proc for proc in procs if proc.is_alive()]
        time.sleep(.1)

    #averageList = mpi.parfor( average, list(zip(splitLists, avgNameList, [showProgressBar]*splitFactor,
    #                                       [verbose]*splitFactor, [createInfoVolumes]*splitFactor,
    #                                            [weighting]*splitFactor, [norm]*splitFactor)), verbose=True)
    
    #collect results from files
    unweiAv = read(preList[0])
    wedgeSum = read(wedgeList[0])
    os.system('rm ' + wedgeList[0])
    os.system('rm ' + avgNameList[0])
    os.system('rm ' + preList[0])
    for ii in range(1,splitFactor):
        av = read(preList[ii])
        unweiAv += av
        os.system('rm ' + preList[ii])
        w = read(wedgeList[ii])
        wedgeSum += w
        os.system('rm ' + wedgeList[ii])
        os.system('rm ' + avgNameList[ii])

    
    if createInfoVolumes:
        unweiAv.write(averageName[:len(averageName)-3]+'-PreWedge.em')
        wedgeSum.write(averageName[:len(averageName)-3] + '-WedgeSumUnscaled.em')

    # convolute unweighted average with inverse of wedge sum
    invert_WedgeSum( invol=wedgeSum, r_max=unweiAv.sizeX()/2-2., lowlimit=.05*len(particleList),
                     lowval=.05*len(particleList))
    fResult = fft(unweiAv)
    r = complexRealMult(fResult,wedgeSum)
    unweiAv = ifft(r)
    unweiAv.shiftscale(0.0,1/float(unweiAv.sizeX()*unweiAv.sizeY()*unweiAv.sizeZ()))
    # low pass filter to remove artifacts at fringes
    unweiAv = lowpassFilter(volume=unweiAv, band=unweiAv.sizeX()/2-2, smooth=(unweiAv.sizeX()/2-1)/10.)[0]

    unweiAv.write(averageName)

    return Reference(averageName,particleList)

def run(fname, outname, cores=6):


    even = ParticleList()

    even.fromXMLFile(fname)

    aa = averageParallel(particleList=even,
                         averageName=outname,
                         showProgressBar=True, verbose=False, createInfoVolumes=False,
                         weighting=False, norm=False,
                         setParticleNodesRatio=3, cores=cores)

    #evenAverage.getVolume().write('Even.em')
    



if __name__=='__main__':
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options

    options = [ScriptOption(['-p','--particleList'], 'Particle List', True, False),
               ScriptOption(['-a','--averageName'], 'Filename of output average.', True, False),
               ScriptOption(['-c','--numberOfCores'],'Number of Cores used for average calculation', True, False),
               ScriptOption(['-w','--weighting'],'Weight particles by exp CC in average. False by default.', False, True),
               ScriptOption(['-v','--verbose'],'Print particle information. False by default.', False, True),
               ScriptOption(['-s','--showProgressBar'],'Show progress bar. False by default.', False, True),
               ScriptOption(['-i','--createInfoVolumes'],'Create Info data (wedge sum, inverted density) too? False by default.', False, True),
               ScriptOption(['-n','--normalize'],'Normalize average. False by default.', False, True),
               ScriptOption(['-h', '--help'], 'Help.', False, True)]


    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Convert coordinate list to particle list.',
                          authors='Friedrich Foerster',
                          options=options)
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        plName, outname, cores, weighting, verbose, showProgressBar, createInfoVol, norm, help = dd = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()

    if help is True:
        print(helper)
        sys.exit()

    try: cores = int(cores)
    except: cores = 1

    if not os.path.exists(plName):
        print('Please provide an existing particle list')
        sys.exit()

    even = ParticleList()
    even.fromXMLFile(plName)

    averageParallel(particleList=even,
                    averageName=outname,
                    showProgressBar=showProgressBar, verbose=verbose, createInfoVolumes=createInfoVol,
                    weighting=weighting, norm=norm,
                    setParticleNodesRatio=3, cores=cores)
    
