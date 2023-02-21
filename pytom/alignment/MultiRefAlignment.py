'''
Created on Jan 2, 2012

@author: hrabe
'''
from pytom.cluster.mcoEXMX import classifyParticleList,distributeExpectation


def multiRef_EXMXAlign(multiRefJob,doFinalize=True,verbose=False):
    """
    multiRef_EXMXAlign: Performs multi reference alignment on a particle list
    @param multiRefJob: The multi reference alignment job 
    @param doFinalize: Send finalize msgs to workers or not. Default is true  
    @param verbose: Default is false
    """
    import pytom.lib.pytom_mpi as pytom_mpi
    
    if doFinalize:
        pytom_mpi.init()
     
    if pytom_mpi.rank() == 0:
        
        from pytom.alignment.ExMaxAlignment import ExMaxManager
        from pytom.tools.files import checkDirExists
        from os import mkdir
        from pytom.basic.resolution import bandToAngstrom,angstromToBand,angleFromResolution
        
        particleList = multiRefJob.getParticleList()
        initialParticleList = particleList
        previousParticleList = initialParticleList
        
        
        destinationDirectory = multiRefJob.getDestinationDirectory()
        numberIterations = multiRefJob.getNumberIterations()
        numberClasses = multiRefJob.getNumberClasses()
         
        exMaxJob = multiRefJob.getExMaxJob()
        p = particleList[0]
        pVol = p.getVolume()
        cubeSize = pVol.sizeX()
        
        preprocessing = exMaxJob.getPreprocessing()
        sampleInfo = exMaxJob.getSampleInformation()
        
        if verbose:
            print(multiRefJob)
            
        if not checkDirExists(destinationDirectory):
            raise IOError('Destination directory ' + destinationDirectory + ' not found!')
        
        try:
            particleLists = particleList.splitByClass()
            if len(particleLists) <= 1:
                raise Exception()
            
        except Exception:
            from pytom.cluster.clusterFunctions import randomiseParticleListClasses
            
            if numberClasses:
                if verbose:
                    print('Randomizing particle list')
                    
                particleList = randomiseParticleListClasses(particleList,numberClasses)
                particleList.toXMLFile(destinationDirectory + '/RandomisedParticleList.xml')
                particleLists = particleList.splitByClass()
                
            else:
                raise RuntimeError('The particle list provided is not pre-classified and you did not set numberClasses for a random seed!')
        
        iteration = 0
        converged = False
        
        while iteration < numberIterations and (not converged):
            
            if verbose:
                print('Running iteration ' + str(iteration) + ' of ' +str(numberIterations))
                
            iterationDirectory = destinationDirectory + '/' + str(iteration) + '/'
            
            if not checkDirExists(iterationDirectory):
                mkdir(iterationDirectory)
            
            
            #determine resolution of all classes 
            maxRes = 0
            minRes = 1000000
            
            if not checkDirExists(iterationDirectory+'resolution/'):
                mkdir(iterationDirectory+'resolution/')
            
            for classIterator in range(len(particleLists)):
                currentParticleList = particleLists[classIterator]
                if len(currentParticleList) > 1:
                    [resNyquist,resolutionBand,numberBands] = currentParticleList.determineResolution( criterion= exMaxJob.getFSCCriterion(), numberBands = cubeSize / 2, mask= exMaxJob.getMask(), keepHalfsetAverages = False, halfsetPrefix=iterationDirectory +'resolution/' + 'class'+str(classIterator)+'_fsc-', verbose=verbose )
                else:
                    continue
                resolutionAngstrom = bandToAngstrom(resolutionBand,sampleInfo.getPixelSize(),numberBands,1 )
                #resolutionAngstrom = bandToAngstrom(resolutionBand,sampleInfo.getPixelSize(),numberBands,exMaxJob.getBinning() )
                
                if resolutionBand > maxRes:
                    maxRes = resolutionBand
                if resolutionBand < minRes:
                    minRes = resolutionBand
                    
                if verbose:
                    print('Class ',classIterator,' - current resolution :' + str(resolutionAngstrom) + ' Angstrom')
            
            #set highest frequency according to user specification
            band = maxRes
            if not multiRefJob.getUseMaxResolution():
                band = minRes
            
            if band == numberBands:   
                #determineResolution returns numberBands for filter if fsc result is invalid. in that case, use nyquist /2 as filter setting
                print('Warning MultiRefAlignment.py: LL 114')
                print('Warning: Resolution determined for all classes was invalid. Will use Nyquist/2 for current iteration') 
                band = numberBands / 2
                
            preprocessing.setHighestFrequency(band)
            exMaxJob.setPreprocessing(preprocessing)
            
            
            alignmentLists = [None] * len(particleLists)
            #generate cluster centers
            referenceList = distributeExpectation(particleLists,iterationDirectory,'clusterCenter'+ str(iteration),verbose,exMaxJob.getSymmetry())
            
            for classIterator in range(len(particleLists)):
                
                classDirectory = iterationDirectory + 'class' + str(classIterator) + '/'

                #determine distance for all particles 
                refinementDirectory = classDirectory + 'refinement/'
                
                if verbose:
                    print(refinementDirectory)
                    
                if not checkDirExists(refinementDirectory):
                    mkdir(refinementDirectory)
                 
                exMaxJob.setParticleList(particleList)
                exMaxJob.setReference(referenceList[classIterator])
                exMaxJob.setDestination(refinementDirectory)
                
                #run refinement
                manager = ExMaxManager(exMaxJob)
            
                manager.distributeAlignment(verbose)

                alignmentLists[classIterator] = manager.getAlignmentList()    
                alignmentLists[classIterator].toXMLFile(iterationDirectory + 'AlignmentList'+ str(classIterator) +'.xml')


            #perform classification here    
            if verbose:
                print('Classifying after iteration ' + str(iteration))
                
            particleList = classifyParticleList(initialParticleList,alignmentLists,verbose)
            particleList.toXMLFile(iterationDirectory + 'classifiedParticles.xml')
            particleLists = particleList.splitByClass()
            
            
            difference = previousParticleList.classDifference(particleList)
            converged = multiRefJob.getEndThreshold() >= difference[3] 
            
            #set up for next round!
            previousParticleList = particleList
            iteration = iteration + 1
        
        if doFinalize:    
            manager.parallelEnd()
            pytom_mpi.finalise()
               
        return [particleList,alignmentLists]
        
    else:
        from pytom.alignment.ExMaxAlignment import ExMaxWorker
        worker = ExMaxWorker()
        worker.parallelRun()
        pytom_mpi.finalise()


