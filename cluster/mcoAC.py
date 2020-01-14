'''
Created on Feb 8, 2011

@author: hrabe
'''


def classifyParticleList(particleList,alignmentLists,criterion,temperature,verbose=False):
    """
    classifyParticleList:
    @param particleList:
    @param alignmentLists:
    @param criterion: 
    @param temperature:
    @param verbose: Print classification process. (Default is False)   
    """

    from pytom.basic.structures import ParticleList
    from pytom.cluster.mcoEXMXStructures import SwapList
    
    for alignmentList in alignmentLists:
        alignmentList.sortByParticleList(particleList)
    
    temperature.initializeTemperature(alignmentLists)
    
    priorClassNumber = len(particleList.splitByClass())
    
    iteration = 0
    clusterSizeSame = False
    
    #prevent that annealing deletes clusters. 
    #leave at least one particle in class! repeat for 10 iterations if criterion.getAllowClassCollapse() == False!
    
    while iteration < 10 and not clusterSizeSame:
        swapList = SwapList()
        classifiedParticleList = ParticleList(particleList.getDirectory())
        
        for particleIndex in range(len(particleList)):
            
            particle = particleList[particleIndex]

            results = []
            
            for alignmentIterator in range(len(alignmentLists)):
                
                alignmentList = alignmentLists[alignmentIterator]
                results.append(alignmentList[particleIndex])
                
            #results = sorted(results, key=lambda MaximisationResult: MaximisationResult.getScore().getValue(),reverse=True)    
            
            if verbose:
                print('')
                print('')
                print('')
                print('----------------------')
                print(results)
            
            [bestResult, bestCluster,swapList] = criterion.apply(results,temperature,swapList,False)
            
            if bestResult == []:
                raise RuntimeError('Could not determine best cluster for particle ' + particle.getFilename())
            
            classifiedParticle = bestResult.toParticle()
            classifiedParticle.setClass(bestCluster)
            classifiedParticle.setRotation(particle.getRotation())
            classifiedParticle.setShift(particle.getShift())
            
            if verbose:
                print(classifiedParticle)
                
            classifiedParticleList.append(classifiedParticle)
    
        iteration += 1

        clusterSizeSame = priorClassNumber == len(classifiedParticleList.splitByClass(False))
        clusterSizeSame = clusterSizeSame or criterion.getAllowClassCollapse()
    
    return [classifiedParticleList,swapList]
    
    
    
def classifyParticleListThreads(particleList,alignmentLists,criterion,temperature,verbose=False):
    """
    classifyParticleListThreads: Same as above, but distributes classifyParticleList to threads
    @param particleList:
    @param alignmentLists:
    @param criterion: 
    @param temperature:
    @param verbose: Print classification process. (Default is False)   
    """

    from pytom.basic.structures import ParticleList
    from pytom.cluster.mcoEXMXStructures import SwapList
    
    for alignmentList in alignmentLists:
        alignmentList.sortByParticleList(particleList)
    
    temperature.initializeTemperature(alignmentLists)
    
    priorClassNumber = len(particleList.splitByClass())
    
    iteration = 0
    clusterSizeSame = False
    
    #prevent that annealing deletes clusters. 
    #leave at least one particle in class! repeat for 10 iterations if criterion.getAllowClassCollapse() == False!
    
    while iteration < 10 and not clusterSizeSame:
        swapList = SwapList()
        classifiedParticleList = ParticleList(particleList.getDirectory())
        
        for particleIndex in range(len(particleList)):
            
            particle = particleList[particleIndex]

            results = []
            
            for alignmentIterator in range(len(alignmentLists)):
                
                alignmentList = alignmentLists[alignmentIterator]
                results.append(alignmentList[particleIndex])
                
            #results = sorted(results, key=lambda MaximisationResult: MaximisationResult.getScore().getValue(),reverse=True)    
            
            if verbose:
                print('')
                print('')
                print('')
                print('----------------------')
                print(results)
            
            [bestResult, bestCluster,swapList] = criterion.apply(results,temperature,swapList,False)
            
            if bestResult == []:
                raise RuntimeError('Could not determine best cluster for particle ' + particle.getFilename())
            
            classifiedParticle = bestResult.toParticle()
            classifiedParticle.setClass(bestCluster)
            classifiedParticle.setRotation(particle.getRotation())
            classifiedParticle.setShift(particle.getShift())
            
            if verbose:
                print(classifiedParticle)
                
            classifiedParticleList.append(classifiedParticle)
    
        iteration += 1

        clusterSizeSame = priorClassNumber == len(classifiedParticleList.splitByClass(False))
        clusterSizeSame = clusterSizeSame or criterion.getAllowClassCollapse()
    
    return [classifiedParticleList,swapList]
    

def mcoAC(annealingJob,doFinalize=True,verbose=False):
    """
    mcoAC: Performs mcoAC clustering on particleList 
    @param annealingJob:  
    @param doFinalize: Send finalize messages to workers or not. Default is true. Should be false when this process is integrated into another parallel process.   
    @param verbose: Default is false
    """
    import pytom_mpi
    
    if doFinalize:
        pytom_mpi.init()
     
    if pytom_mpi.rank() == 0:
        
        from pytom.cluster.mcoEXMX import mcoEXMX
        
        from pytom.tools.files import checkDirExists
        from os import mkdir,system
        
        particleList = annealingJob.getParticleList()
        
        if len(particleList) == 0:
            raise RuntimeError('Particle list is empty! Abort!')
        
        initialParticleList = particleList
        previousParticleList = initialParticleList
        
        destinationDirectory = annealingJob.getDestinationDirectory()
        numberIterations = annealingJob.getNumberIterations()
        numberClasses = annealingJob.getNumberClasses() 

        if verbose:
            print(annealingJob)
            
        if not checkDirExists(destinationDirectory):
            raise IOError('Destination directory ' + destinationDirectory + ' not found!')
        
        try:
            particleLists = particleList.splitByClass()
            if len(particleLists) <= 1 or (len(particleLists) == 1 and len(particleLists[0]) == len(particleList)):
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
        
        bestScoreSum = None
        bestParticleList = None
        
        while (not annealingJob.cooledDown()) and (not converged):
            
            if verbose:
                print('Running iteration ' + str(iteration) + ' of ' +str(numberIterations))
        
            iterationDirectory = destinationDirectory + '/' + str(iteration) + '/'
            mcoEXMXDirectory = iterationDirectory + 'mcoEXMX/'
            
            
            if not checkDirExists(iterationDirectory):
                mkdir(iterationDirectory)
            
            annealingJob.setDestinationDirectory(mcoEXMXDirectory)
            annealingJob.setParticleList(previousParticleList)
            
            annealingJob.setNumberIterations(annealingJob.getLocalIncrement())
            
            annealingJob.toXMLFile(iterationDirectory + 'annealingJob.xml')
            
            #run local refinement
            [pl, alignmentLists] = mcoEXMX(annealingJob, doFinalize = False, verbose = verbose)
            annealingJob.setNumberIterations(numberIterations)
            
            #store currently best solution
            if iteration == 0 or (bestScoreSum < pl.sumOfScores() and len(pl.splitByClass()) > 1):
                bestScoreSum = pl.sumOfScores() 
                bestParticleList = pl
            
            bestParticleList.toXMLFile(iterationDirectory + 'currentBestParticleList.xml')
            
            #perform classification here    
            [particleList , swapList] = classifyParticleList(initialParticleList,alignmentLists,annealingJob.getCriterion(),annealingJob.getTemperature().copy(),verbose)
            
            #save iteration results to disk 
            particleList.toXMLFile(iterationDirectory + 'classifiedParticles.xml')
            swapList.toXMLFile(iterationDirectory + 'swapList.xml')
            
            #if not verbose mode, delete mcoEXMX files
            if not verbose:
                system('rm -rf ' + mcoEXMXDirectory)
            
            #print number class swaps
            difference = previousParticleList.classDifference(particleList)
            converged = annealingJob.getEndThreshold() >= difference[3]  
            
            previousParticleList = particleList
                
            #set up for new round
            annealingJob.decreaseTemperature()
            particleLists = particleList.splitByClass()
            
            iteration = iteration + 1
            
            print('Annealing iteration ' + str(iteration) + ' finished!')
            
        if doFinalize: 
            from pytom.alignment.ExMaxAlignment import ExMaxManager
               
            manager = ExMaxManager(annealingJob.getExMaxJob())
            manager.parallelEnd()
            pytom_mpi.finalise()
            
        return bestParticleList
        
    else:
        from pytom.cluster.mcoEXMXStructures import MCOEXMXWorker
        worker = MCOEXMXWorker()
        worker.setDoAlignment(annealingJob.getDoAlignment())
        worker.setFRMBandwidth(annealingJob.getFRMBandwidth())
        worker.parallelRun(verbose = verbose)
        pytom_mpi.finalise()

