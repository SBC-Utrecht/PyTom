'''
Created on Dec 22, 2010

@author: hrabe
'''

def classifyParticleList(particleList,alignmentLists,verbose=False):   
    """
    classifyParticleList: Classifies input particle list according to scores determined in alignment lists. 
    @param particleList:
    @param alignmentLists:
    @param verbose:
    @return: Will return particle list in same order as input but with assigned class memberships.
    """
    
    from pytom.basic.structures import ParticleList
    classifiedParticleList = ParticleList(particleList.getDirectory())
    
    for alignmentList in alignmentLists:
        alignmentList.sortByParticleList(particleList)
    
    
    for particleIndex in range(len(particleList)):
        
        particle = particleList[particleIndex]
        score = particle.getScore()
        
        if verbose:
            print(particle.getFilename())
        
        if not score:
            bestScore= -999999999
        else:
            bestScore = particle.getScore().getWorstValue()
            
        bestResult = []
        bestCluster = []
        
        for alignmentIterator in range(len(alignmentLists)):
            alignmentList = alignmentLists[alignmentIterator]
            
            result = alignmentList[particleIndex]
            
            score = result.getScore()
                
            value = float(score.getValue());

            if verbose:
                print(value)
                
            if value >= bestScore:
                bestResult = result
                bestCluster = alignmentIterator
                bestScore = value
        
        if bestResult == []:
            raise Exception('Could not determine best cluster for particle ' + particle.getFilename())
        
        classifiedParticle = bestResult.toParticle()
        classifiedParticle.setClass(bestCluster)
        
        if verbose:
            print(classifiedParticle)
            
        classifiedParticleList.append(classifiedParticle)
    
    return classifiedParticleList

def distributeExpectation(particleLists,iterationDirectory,averagePrefix,verbose=False,symmetry = None):
    """
    distributeExpectation: Distributes particle expectation (averaging) to multiple workers. Required by many algorithms such as MCOEXMX  
    @param particleLists: list of particleLists 
    @param iterationDirectory:
    @param averagePrefix:  
    @param verbose:
    @param symmetry:  
    """
    import pytom_mpi
    from pytom.tools.files import checkDirExists
    from pytom.parallel.alignmentMessages import ExpectationJobMsg,ExpectationResultMsg
    from pytom.alignment.structures import ExpectationJob
    from pytom.basic.structures import Reference,ReferenceList
    from os import mkdir
    
    if not pytom_mpi.isInitialised():
            pytom_mpi.init()
            
    mpi_myid = pytom_mpi.rank()
        
    if not mpi_myid == 0:
        raise RuntimeError('This function (distributeExpectation) can only be processed by mpi_id = 0! ID == ' + str(mpi_myid) + ' Aborting!')

    if not checkDirExists(iterationDirectory):
            raise IOError('The iteration directory does not exist. ' + iterationDirectory)
        
    mpi_numberNodes = pytom_mpi.size()
    
    if mpi_numberNodes <= 1:
        raise RuntimeError('You must run clustering with openMPI on multiple CPUs')
    
    listIterator = 0
        
    referenceList = ReferenceList()
    
    #distribute jobs to all nodes
    for i in range(1,mpi_numberNodes):
        
        if verbose:
                print('Starting first job distribute step')
                
        if listIterator < len(particleLists):
            
            if not checkDirExists(iterationDirectory + 'class' + str(listIterator) + '/'):
                mkdir(iterationDirectory + 'class' + str(listIterator) + '/')
            
            averageName = iterationDirectory + 'class' + str(listIterator) + '/' + averagePrefix + '-' +str(listIterator) + '.em'
            
            if not symmetry.isOneFold():
                newPl = symmetry.apply(particleLists[listIterator])
                job = ExpectationJob(newPl,averageName)
            else:
                job = ExpectationJob(particleLists[listIterator],averageName)
                
            newReference = Reference(averageName,particleLists[listIterator])
            
            referenceList.append(newReference)
            
            jobMsg = ExpectationJobMsg(0,str(i))
            jobMsg.setJob(job)

            pytom_mpi.send(str(jobMsg),i)
            if verbose:
                print(jobMsg)
            
            listIterator = listIterator + 1
            
            
    finished = False

    #there are more jobs than nodes. continue distributing and collect results
    receivedMsgCounter = 0
    
    while not finished:

        #listen and collect
        mpi_msgString = pytom_mpi.receive() 
        
        if verbose:
            print(mpi_msgString)

        jobResultMsg = ExpectationResultMsg('','')   
        jobResultMsg.fromStr(mpi_msgString)
        
        receivedMsgCounter = receivedMsgCounter + 1
        
        #send new job to free node
        if listIterator < len(particleLists):
            
            if not checkDirExists(iterationDirectory + 'class' + str(listIterator) + '/'):
                mkdir(iterationDirectory + 'class' + str(listIterator) + '/')
                
            averageName = iterationDirectory + 'class' + str(listIterator) + '/' + averagePrefix + '-' +str(listIterator) + '.em'
            job = ExpectationJob(particleLists[listIterator],averageName)
            newReference = Reference(averageName,particleLists[listIterator])
            referenceList.append(newReference)
            
            
            jobMsg = ExpectationJobMsg(0,str(jobResultMsg.getSender()))
            jobMsg.setJob(job)
            
            pytom_mpi.send(str(jobMsg),i)
            if verbose:
                print(jobMsg)
            
            listIterator = listIterator + 1
            
        finished = listIterator >= len(particleLists) and receivedMsgCounter == len(particleLists) 
    
    return referenceList
    
def mcoEXMX(mcoEMJob,doFinalize=True,verbose=False):
    """
    mcoEXMX: Perfomrs kmeans clustering on particleList
    @param mcoEMJob: The clustering job 
    @param doFinalize: Send finalize msgs to workers or not. Default is true  
    @param verbose: Default is false
    """
    import pytom_mpi
    
    if doFinalize:
        pytom_mpi.init()
     
    if pytom_mpi.rank() == 0:
        
        from pytom.alignment.ExMaxAlignment import ExMaxManager
        from pytom.tools.files import checkDirExists
        from os import mkdir
        from pytom.basic.plot import plotClassSizes
        from builtins import min as minList
        
        particleList = mcoEMJob.getParticleList()

        if len(particleList) == 0:
            raise RuntimeError('Particle list is empty! Abort!')

        destinationDirectory = mcoEMJob.getDestinationDirectory()
        numberIterations = mcoEMJob.getNumberIterations()
        numberClasses = mcoEMJob.getNumberClasses() 
        exMaxJob = mcoEMJob.getExMaxJob()
        
        if verbose:
            print(mcoEMJob)
            
        if not checkDirExists(destinationDirectory):
            raise IOError('Destination directory ' + destinationDirectory + ' not found!')
        
        try:
            particleLists = particleList.splitByClass()
            if len(particleLists) < 1 or (len(particleLists) == 1 and len(particleLists[0]) == len(particleList)):
                raise Exception()
            
        except Exception:
            from pytom.cluster.clusterFunctions import randomiseParticleListClasses
            
            if numberClasses:
                if verbose:
                    print('Randomising particle list')
                    
                particleList = randomiseParticleListClasses(particleList,numberClasses)
                particleList.toXMLFile(destinationDirectory + '/RandomisedParticleList.xml')
                particleLists = particleList.splitByClass()
                
            else:
                raise RuntimeError('The particle list provided is not pre-classified and you did not set numberClasses for a random seed!')
        
        initialParticleList = particleList
        previousParticleList = initialParticleList
        
        iteration = 0
        converged = False
        doAdaptiveResolution = mcoEMJob.getAdaptiveResolution()
        
        if doAdaptiveResolution:
            preProcessing = exMaxJob.getPreprocessing()
            highestFrequency = preProcessing.getHighestFrequency() 
            resolutionList = [highestFrequency] * len(particleLists) 
        
        while iteration < numberIterations and (not converged):
            
            if verbose:
                print('Running iteration ' + str(iteration) + ' of ' +str(numberIterations))
        
            iterationDirectory = destinationDirectory + '/' + str(iteration) + '/'
            
            if not checkDirExists(iterationDirectory):
                mkdir(iterationDirectory)
            
            #referenceList = ReferenceList()
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
                
                #use adaptive resolution -> update lowpass filter

                if doAdaptiveResolution and len(resolutionList) > 0 and resolutionList[classIterator] > 0:
#                    preProcessing = exMaxJob.getPreprocessing()
#                    resolution = resolutionList[classIterator]
#                    preProcessing.setHighestFrequency(resolution)
#                    exMaxJob.setPreprocessing(preProcessing)

                    preProcessing = exMaxJob.getPreprocessing()
                    resolution = minList(resolutionList) *1.1
                    preProcessing.setHighestFrequency(resolution)
                    exMaxJob.setPreprocessing(preProcessing)
                    
                    
                
                #run refinement
                manager = ExMaxManager(exMaxJob)
                exMaxJob.toXMLFile(classDirectory + 'Job.xml')
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
            converged = mcoEMJob.getEndThreshold() >= difference[3] 
            
            #determine resolution in each class
            if doAdaptiveResolution:
                resolutionList = [-1] * len(particleLists)
                
                for classIterator in range(len(particleLists)):
                    classList = particleLists[classIterator]
                    
                    if len(classList) == 1:
                        #if there is only one particle in that class, override resolution
                        print('Class ', classIterator , ' has only 1 particle! Will be assigned the lowest resolution determined.')
                        continue
                    
                    className = classList[0].getClass()
                    
                    v = classList[0].getVolume()
                    cubeSize = v.sizeX()
                    
                    resolution = classList.determine_resolution(criterion=0.5,number_bands = cubeSize / 2,mask=exMaxJob.getMask(),verbose=False,plot='',keepHalfsetAverages = False,halfsetPrefix='class' + str(className), parallel=True)
                    #resolution = [Resolution in Nyquist , resolution in band, number_bands]
                    resolutionList[classIterator] = resolution[1]
                    print('Resolution for class ', classIterator , ' determined to ', resolution[1], ' pixels. Class size is ', len(classList) , ' particles')
            
                #get lowest resolution determined for classes with more than 1 particle
                min = 999999999999999
                for classIterator in range(len(particleLists)):
                    if min >= resolutionList[classIterator] and resolutionList[classIterator] >= 0:
                        min = resolutionList[classIterator]
                        
                #set resolution for all classes with only 1 particle to lowest resolution    
                for classIterator in range(len(particleLists)):
                    if resolutionList[classIterator] < 0:
                        resolutionList[classIterator] = min
                    
            #set up for next round!
            previousParticleList = particleList
            iteration = iteration + 1
        
        if doFinalize:    
            manager.parallelEnd()
            pytom_mpi.finalise()
               
        return [particleList,alignmentLists]
        
    else:
        from pytom.cluster.mcoEXMXStructures import MCOEXMXWorker
        worker = MCOEXMXWorker()
        worker.setDoAlignment(mcoEMJob.getDoAlignment())
        worker.setFRMBandwidth(mcoEMJob.getFRMBandwidth())
        worker.parallelRun()
        pytom_mpi.finalise()
        
#--------------------------------------------------------------------------------

def clusterKMeans(data, k, scale=True, plot=False):
    """
    clusterKMeans: k-means clustering function
    @param data: data set to cluster
    @type data: list of list or {pytom.cluster.kmeansStructures.KMDataSet}
    @param k: the number of clusters
    @type k: integer
    @param scale: unify the data along each dimension or not
    @type scale: boolean
    @param plot: plot the clustering procedure or not (if the data's dimension is more than 2, you can specified the plotting dimension as e.g. [0,1])
    @type plot: boolean or list
    @return: the clustered data set
    @rtype: {pytom.cluster.kmeans.KMDataSet}
    """
    assert (k>1)
    from pytom.cluster.kmeansStructures import KMDataSet
    
    if data.__class__ == KMDataSet:
        kmdata = data
    else:
        kmdata = KMDataSet(data)
    
    if scale:
        from pytom.cluster.kmeans import scaleKMDataSet
        kmdata = scaleKMDataSet(kmdata)
    
    if len(kmdata.clusters) == 0:
        kmdata.initializeClusters(k)
    
    # initialize the plotting and plot the original data set
    colors = ['b', 'g', 'r', 'c', 'y', 'm']
    if k > len(colors):
        plot = False
    if plot:
        # determine the plot dimension
        if plot.__class__ == list:
            dim_x = plot[0]; dim_y = plot[1]
        else:
            dim_x = 0; dim_y = 1
        
        import pylab
        pylab.ion() # interactive mode on
        plt_centroids=[]
        plt_points=[]
        for i in range(k):
            c = kmdata.clusters[i].centroid
            tmp, = pylab.plot(c[dim_x],c[dim_y], color=colors[i],marker='*',markersize=12)
            plt_centroids.append(tmp)
            
        for i in range(kmdata.num):
            d = kmdata.getData(i)
            tmp, = pylab.plot(d[dim_x],d[dim_y], color='k',marker='+',markersize=12)
            plt_points.append(tmp)
            
        pylab.draw()
    
    # use l2 distance here
    from pytom.tools.maths import euclidianDistance
    stop = 1
    while(stop > 0.01):
        changedNum = 0.
        for i in range(kmdata.num):
            if kmdata.assignCluster(i, euclidianDistance):
                changedNum = changedNum+1
        stop = changedNum/kmdata.num
        
        # update the centroid of each cluster
        kmdata.updateClusters()
        
        # handle empty cluster
        for i in kmdata.findEmptyClusters():
#            print 'Empty cluster'
            # remove the empty cluster
            kmdata.removeCluster(i)
            # find the cluster with the largest SSE and split
            kmdata.sortClusters(euclidianDistance)
            stop = stop + len(kmdata.clusters[-1].members)/kmdata.num
            kmdata.splitCluster(-1)
        
#        print stop
        
        # plot
        if plot:
            for i in range(k):
                c = kmdata.clusters[i].centroid
                plt_centroids[i].set_xdata(c[dim_x])
                plt_centroids[i].set_ydata(c[dim_y])
            
            for i in range(kmdata.num):
                c = 0
                for j in range(k):
                    if(i in kmdata.clusters[j].members):
                        c = j
                        break
                plt_points[i].set_color(colors[c])
            
            pylab.draw()
            
#            import sys
#            sys.stdin.readline()
    
    return kmdata

def readKMFile(datafile, clusterfile=None):
    """
    readKMFile: read the file used for storing the data set for k-means clustering
    """
    
    # read the data file firstly
    f = open(datafile)
    
    from pytom.cluster.kmeansStructures import KMDataSet, Cluster
    kmdata = KMDataSet()
    
    for line in f:
        tmp = line.split(' ')
        kmdata.data.append([float(n) for n in tmp])
        
    f.close()
    kmdata.refresh()
    
    # read the label file
    if clusterfile:
        f = open(clusterfile)
        
        class_labels = {}
        i = 0
        for line in f:
            class_id = int(line)
            if class_id == -1:
                pass
            else:
                if class_id in list(class_labels.keys()):
                    class_labels[class_id].members.append(i)
                else:
                    c = Cluster(kmdata, None, name=str(class_id))
                    c.members.append(i)
                    kmdata.clusters.append(c)
                    class_labels[class_id] = c
            i = i + 1
        
        f.close()
    
    return kmdata

def readClassFile(classfile):
    class_labels=[]
    f=open(classfile, 'r')
    for line in f:
        class_labels.append(int(line))
    f.close()
    
    return class_labels

def writeKMResult(result, filename):
    """
    writeKMResult: write the clustering result to disk
    @param result: k-means result
    @type result: {pytom.cluster.kmeansStructures.KMDataSet}
    @param filename: file name
    @type filename: string
    """
    f = open(filename, 'w')
    for i in range(result.num):
        c = None
        for j in range(len(result.clusters)):
            if i in result.clusters[j].members:
                c = result.clusters[j]
                c.name = str(j) # set the name of the cluster according to the index
                break
        
        if c:
            line = c.name+'\n'
        else:
            line = '-1\n'
        f.write(line)
    f.close()

def scaleKMDataSet(kmdata):
    """
    scaleKMDataSet: Scale the KMData in each dimension.
    """
#    from scipy import std
#    scale = []
#    for i in xrange(kmdata.dim):
#        scale.append(std([d[i] for d in kmdata.data]))
#    
#    ndata = []
#    for d in kmdata.data:
#        nd = []
#        for i in xrange(kmdata.dim):
#            nd.append(float(d[i])/scale[i])
#        ndata.append(nd)
    from pytom.tools.maths import normalize2D
    ndata = normalize2D(kmdata.data, byRow=False)
    
    from pytom.cluster.kmeansStructures import KMDataSet
    res = KMDataSet(ndata)
    return res

def clusterBisectKMeans(data, k, scale=True):
    assert (k>1)
    bisect_kmdata = clusterKMeans(data, 2, scale)
    
    from pytom.tools.maths import euclidianDistance
    while len(bisect_kmdata.clusters)<k:
        bisect_kmdata.sortClusters(euclidianDistance)
        bisect_kmdata.splitCluster(-1)
    
    from pytom.cluster.kmeansStructures import KMDataSet
    centroids = [c.centroid for c in bisect_kmdata.clusters]
    kmdata = KMDataSet(data)
    kmdata.initializeClusters(k, centroids)

    kmdata = clusterKMeans(kmdata, k, scale)
    
    return kmdata
