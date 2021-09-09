'''
@todo: ADD UNIT TESTS!!!
Created on Mar 23, 2010
A correlation matrix is one basic mean for the classification / clustering of a particle list. 
@author: hrabe
'''

        
def calculateCorrelationVector(particle,particleList,mask,particleIndex,applyWedge=True,binningFactor=0,lowestFrequency=0,highestFrequency=1):
    """
    calculateCorrelationVector: 
    @param particle: the current particle
    @param particleList: all the other particles
    @param mask: corellation under a mask
    @param applyWedge: Apply constrained correlation? True by default. Set to false for rotation classification for instance.
    @param binningFactor: Binning of data when read. Default is 0 (larger values according to libtomc notation)   
    @param lowestFrequency: Lowest frequency for bandpass
    @param highestFrequency: Highest frequency for bandpass
    """
    from pytom.basic.structures import Particle,ParticleList,Mask 
    from pytom.cluster.correlationMatrixStructures import CorrelationVector
    from pytom_volume import vol,rotate,shift
    from pytom.basic.filter import bandpassFilter,filter
    
    assert particle.__class__ == Particle
    assert particleList.__class__ == ParticleList
    assert mask.__class__ == Mask or mask.__class__ == str
    
    from pytom.score.score import FLCFScore
    from pytom.tools.memory import read as readBuffer
    from pytom_volume import read
    
    if mask.__class__ == str:
        mask = read(mask,0,0,0,0,0,0,0,0,0,int(binningFactor),int(binningFactor),int(binningFactor))
    else:
        mask = read(mask.getFilename(),0,0,0,0,0,0,0,0,0,int(binningFactor),int(binningFactor),int(binningFactor))
    
    #sc = nxcfScore()
    sc = FLCFScore()
    
    correlationVector = CorrelationVector(particle,particleList,particleIndex)
    particleRotation = particle.getRotation()
    particleShift = particle.getShift()
    
    particleShifted = None
    otherParticleShifted = None
    
    particleRotated = None
    otherParticleRotated = None
    
    bandpassFilterObject = None
    
    for particleIndex in range(len(particleList)):
        
        otherParticle = particleList[particleIndex]
        
        #read from disk
        
        particleVolume = read(particle.getFilename(),0,0,0,0,0,0,0,0,0,int(binningFactor),int(binningFactor),int(binningFactor))
        
        #otherParticleVolume = read(otherParticle.getFilename(),binningX=binningFactor,binningY=binningFactor,binningZ=binningFactor)
        otherParticleVolume = read(otherParticle.getFilename(),0,0,0,0,0,0,0,0,0,int(binningFactor),int(binningFactor),int(binningFactor))
        
        #initialise memory for buffer volumes
        if not particleRotated:
            particleRotated = vol(particleVolume.sizeX(),particleVolume.sizeY(),particleVolume.sizeZ())
            otherParticleRotated = vol(particleVolume.sizeX(),particleVolume.sizeY(),particleVolume.sizeZ())
            particleShifted = vol(particleVolume.sizeX(),particleVolume.sizeY(),particleVolume.sizeZ())
            otherParticleShifted = vol(particleVolume.sizeX(),particleVolume.sizeY(),particleVolume.sizeZ())
            
        #get settings
        rotationOtherParticle = otherParticle.getRotation()
        otherParticleShift = otherParticle.getShift()
        
        #apply shift determined
        if particleShift.getX() != 0 or particleShift.getY() != 0 or particleShift.getZ() != 0: 
            shift(particleVolume,particleShifted,-particleShift[0],-particleShift[1],-particleShift[2]);
        else:
            particleShifted = particleVolume
            
        if otherParticleShift.getX() != 0 or otherParticleShift.getY() != 0 or otherParticleShift.getZ() != 0:
            shift(otherParticleVolume,otherParticleShifted,-otherParticleShift[0],-otherParticleShift[1],-otherParticleShift[2]);
        else:
            otherParticleShifted = otherParticleVolume
        
        
        #apply rotation determined    
        if particleRotation.getZ1() != 0 or particleRotation.getX() != 0 or particleRotation.getZ2() != 0:
            rotate(otherParticleShifted,otherParticleRotated,-particleRotation.getZ2(),-particleRotation.getZ1(),-particleRotation.getX())
            otherParticleRotated = otherParticleRotated * mask
        else:
            otherParticleRotated = otherParticleVolume
            
        if rotationOtherParticle.getZ1() != 0 or rotationOtherParticle.getX() != 0 or rotationOtherParticle.getZ2() != 0:  
            rotate(particleShifted,particleRotated,-rotationOtherParticle.getZ2(),-rotationOtherParticle.getZ1(),-rotationOtherParticle.getX())
            particleRotated = particleRotated * mask
        else:
            particleRotated = particleVolume

        applyWedge = particleRotation.getZ1() != 0 or particleRotation.getX() != 0 or particleRotation.getZ2() != 0
        applyWedge = applyWedge or rotationOtherParticle.getZ1() != 0 or rotationOtherParticle.getX() != 0 or rotationOtherParticle.getZ2() != 0
        
        #apply cross wedge
        if applyWedge:
            
            particleWedge = particle.getWedge()
            otherParticleWedge = otherParticle.getWedge()
            
            particleWedge.setRotation(particleRotation)
            otherParticleWedge.setRotation(rotationOtherParticle)
            
            particleVolume = otherParticleWedge.apply(particleRotated)
            otherParticleVolume = particleWedge.apply(otherParticleRotated)
    
        #apply bandpass
        if 0 <= lowestFrequency < 1 and 0 < highestFrequency <= 1:
            if not bandpassFilterObject:
                bandpassString = str(lowestFrequency) + ':' + str(highestFrequency)+';'
                
                r = bandpassFilter(particleVolume,bandpassString)
                
                particleVolume = r[0]
                bandpassFilterObject = r[1]
                #print 'calculateCorrelationVector : Init bandpass'
            else:
                r = list(filter(particleVolume,bandpassFilterObject))
                particleVolume = r[0]
                
                #print 'calculateCorrelationVector : existing bandpass'
                
            r = list(filter(otherParticleVolume,bandpassFilterObject))
            otherParticleVolume = r[0]
        
        #particleVolume.write('p.em')
        #otherParticleVolume.write('op.em')
        #assert False
        
        #apply scoring here
        value = sc.scoringCoefficient(otherParticleVolume,particleVolume,mask)
        
        if value != value:
            print('Error during calculation of correlationMatrix! Check files written to disk (error_part.em,error_otherPart.em). ABORTING!')
            print(particle.getFilename(),otherParticle.getFilename())
            particleVolume.write('error_part.em')
            otherParticleVolume.write('error_otherPart.em')
            
            assert False
            
        correlationVector.append(value)
    
    return correlationVector
    
class CMManager():
    """
    CMManager: Supervises correlation matrix calculation on multiple processes. 
    """
    
    
    def __init__(self,job=None):
        from pytom.cluster.correlationMatrixStructures import CorrelationMatrixJob
        from pytom_volume import vol
        
        if not job or job.__class__ != CorrelationMatrixJob:
            raise ParameterError('Please create instances of CMManager only with CorrelationMatrixJob objects')
        elif job.__class__ == CorrelationMatrixJob:
            self._particleList = None 
            self.fromJob(job)
            
        matrixSize = len(self._particleList)    
        
        self._correlationMatrix = vol(matrixSize,matrixSize,1) 
        self._correlationMatrix.setAll(1)
        
        
    def fromJob(self,job):
        """
        fromJob: Initialise this object with a CorrelationMatrixJob
        @param job: 
        @type job: L{pytom.cluster.correlationMatrixStructures.CorrelationMatrixJob} 
        """
        from pytom.cluster.correlationMatrixStructures import CorrelationMatrixJob
        
        assert job.__class__ == CorrelationMatrixJob
        
        self._particleList = job.getParticleList()
        self._mask = job.getMask()
        self._resultMatrixName = job.getResultMatrixName()
        self._applyWedge = job.getApplyWedge()
        self._binningFactor = job.getBinning()
        self._lowestFrequency = job.getLowestFrequency()
        self._highestFrequency = job.getHighestFrequency()
        
    def _setMatrixValuesFromVector(self,particleIndex,vector):
        """
        _setMatrixValuesFromVector: Fills matrix entries with values in vector. Starts at [particleIndex,particleIndex] and copies values in vector to the row and col in matrix. 
        @param particleIndex: 
        @param vector:  
        """
        from pytom.cluster.correlationMatrixStructures import CorrelationVector
    
        if len(vector) == 0:
            return
    
        assert vector.__class__ == CorrelationVector
        assert particleIndex < self._correlationMatrix.sizeX()
        assert len(vector) <= self._correlationMatrix.sizeX()
        
        self._correlationMatrix.setV(1,particleIndex,particleIndex,0)
        
        for vectorIndex in range(len(vector)):
            
            value = vector[vectorIndex]
            
            self._correlationMatrix.setV(value,particleIndex+vectorIndex+1,particleIndex,0)
            self._correlationMatrix.setV(value,particleIndex,particleIndex+vectorIndex+1,0)
        
    def distributeCalculation(self,mpi_myid,verbose=False):
        """
        distributeCalculation: Distribute calculation of matrix to multiple nodes.
        """
        import pytom_mpi
        from pytom.cluster.correlationMatrixStructures import CorrelationVectorJob
        from pytom.parallel.clusterMessages import CorrelationVectorJobMessage,CorrelationVectorMessage
        from pytom.tools.ProgressBar import FixedProgBar
        
        if not mpi_myid == 0:
            raise Exception('This function (distributeCalculation) can only be processed by mpi_id = 0! ID == ' + mpi_myid.__str__() +' Aborting!')
        
        mpi_myname = 'node_'+mpi_myid.__str__()
        mpi_numberNodes = pytom_mpi.size()
        
        particleIndex = 0
        
        progressBar = FixedProgBar(0,len(self._particleList),'Particles correlated ')
        progressBar.update(0)
        
        
        #distribute on all nodes
        for nodeIndex in range(1,mpi_numberNodes):
            
            if particleIndex < len(self._particleList):
                
                particle = self._particleList[particleIndex]
                
                reducedParticleList = self._particleList[particleIndex+1:]
                
                job = CorrelationVectorJob(particle,reducedParticleList,self._mask,particleIndex,self._applyWedge,self._binningFactor,self._lowestFrequency,self._highestFrequency)
            
                jobMsg = CorrelationVectorJobMessage(str(mpi_myid),str(nodeIndex))
                jobMsg.setJob(job)
                
                if verbose:
                    print(jobMsg)
                
                pytom_mpi.send(str(jobMsg),nodeIndex)
                    
                particleIndex = particleIndex + 1
                    
        numberVectorsReceived = 0
        
        finished = numberVectorsReceived > len(self._particleList)
        
        while not finished:
            
            #listen until numberVectorsReceived > len(self._particleList) and continue distributing
            mpi_msgString = pytom_mpi.receive()
            
            if verbose:
                print(mpi_msgString)
            
            correlationVectorMsg = CorrelationVectorMessage()
            correlationVectorMsg.fromStr(mpi_msgString)
            
            assert correlationVectorMsg.__str__() == mpi_msgString
            
            
            vector = correlationVectorMsg.getVector()
            self._setMatrixValuesFromVector(vector.getParticleIndex(), vector)
            
            self._savePreliminaryResult()
            
            #print 'Result received from ' + correlationVectorMsg.getSender().__str__() + ' and matrix saved to disk.'
            
            numberVectorsReceived = numberVectorsReceived + 1
            
            if particleIndex < len(self._particleList):
                #print 'Send particle number :' , particleIndex
                particle = self._particleList[particleIndex]
                
                reducedParticleList = self._particleList[particleIndex+1:]
                
                job = CorrelationVectorJob(particle,reducedParticleList,self._mask,particleIndex,self._applyWedge,self._binningFactor,self._lowestFrequency,self._highestFrequency)
            
                jobMsg = CorrelationVectorJobMessage(mpi_myid.__str__(),correlationVectorMsg.getSender().__str__())
                jobMsg.setJob(job)
                
                pytom_mpi.send(jobMsg.__str__(),int(correlationVectorMsg.getSender()))
                    
                particleIndex = particleIndex + 1
            
            #update progress bar
            progressBar.update(numberVectorsReceived)
            
            finished = numberVectorsReceived >= len(self._particleList)
        
    def calculateMatrix(self):
        """
        calculateMatrix: Perform calculation of matrix only on this node.
        """
        from pytom.cluster.correlationMatrixStructures import CorrelationVectorJob
        
        for particleIndex in range(len(self._particleList)-1):

            particle = self._particleList[particleIndex]
            reducedParticleList = self._particleList[particleIndex+1:]
            
            job = CorrelationVectorJob(particle,reducedParticleList,self._mask,particleIndex,self._applyWedge,self._binningFactor,self._lowestFrequency,self._highestFrequency)
            
            worker = CMWorker(job)
            
            correlationVector = worker.run()
            
            self._setMatrixValuesFromVector(particleIndex, correlationVector)
    
    def _savePreliminaryResult(self):
        self.saveMatrix()
        
    def saveMatrix(self):
        """
        saveMatrix: Save correlation matrix to disk. 
        """
        
        self._correlationMatrix.write(self._resultMatrixName)
    
    def parallelEnd(self):
        """
        parallelEnd : Sends status message = end to all workers. All workers will terminate upon receiving this message.
        @author: Thomas Hrabe
        """
        import pytom_mpi
        from pytom.parallel.messages import StatusMessage
        
        if not pytom_mpi.isInitialised():
            pytom_mpi.init()
        
        mpi_myid = pytom_mpi.rank()
        
        mpi_numberNodes = pytom_mpi.size()
        for i in range(1,mpi_numberNodes):
            msg = StatusMessage(mpi_myid.__str__(),i.__str__())
            msg.setStatus("End")
            pytom_mpi.send(msg.__str__(),i)
     
        
class CMWorker():
    """
    CMWorker: Performs all the calculation for CMManager
    """
    
    def __init__(self,job=None):
        from pytom.cluster.correlationMatrixStructures import CorrelationVectorJob
        
        if not job or job.__class__ != CorrelationVectorJob:
            raise ParameterError('Please create instances of CMWorker only with CorrelationVectorJob objects')
        elif job.__class__ == CorrelationVectorJob:
            self._particleList = None
            self._particle = None
            self._mask = None 
            self.fromJob(job)
            
    def fromJob(self,job):
        from pytom.cluster.correlationMatrixStructures import CorrelationVectorJob
        
        assert job.__class__ == CorrelationVectorJob
        
        self._particle = job.getParticle()
        self._particleList = job.getParticleList()
        self._mask = job.getMask()
        self._particleIndex = job.getParticleIndex()
        self._applyWedge = job.getApplyWedge()
        self._binningFactor = job.getBinning()
        self._lowestFrequency = job.getLowestFrequency()
        self._highestFrequency = job.getHighestFrequency()
        
    def dumpMsg2Log(self,logfile,msg):
        """
        dumpMsg2Log:
        @param logfile:
        @param msg:
        @author: Thomas Hrabe 
        """
        f = open(logfile,'a')
        f.write(msg.__str__())
        f.close()
       
    def run(self):
        """
        run: Calculates correlation vector
        """
        correlationVector = calculateCorrelationVector(self._particle,self._particleList,self._mask,self._particleIndex,self._applyWedge,self._binningFactor,self._lowestFrequency,self._highestFrequency)
        return correlationVector
    
    
def distributedCorrelationMatrix(job,verbose=False):
    """
    distributedCorrelationMatrix: Performs calculation of correlation matrix either on multiple processes or sequentially.
    """
    import pytom_mpi
    
    pytom_mpi.init()
    
    if pytom_mpi.size() >1:
        
        mpi_myid = pytom_mpi.rank()
        
        if mpi_myid == 0:
            manager = CMManager(job)
            
            manager.distributeCalculation(mpi_myid,verbose)
            
            manager.parallelEnd()
            
            manager.saveMatrix()
            
        else:
            
            from pytom.parallel.clusterMessages import CorrelationVectorMessage,CorrelationVectorJobMessage
            from pytom.parallel.messages import StatusMessage,MessageError
            
            end = False
            while not end:
                mpi_msg = pytom_mpi.receive()
                
                if verbose:
                    print(mpi_msg)
                    
                try:
                    msg = CorrelationVectorJobMessage()
                    msg.fromStr(mpi_msg)
                    
                    worker = CMWorker(msg.getJob())
                    #worker.dumpMsg2Log('node'+str(mpi_myid)+'.log', msg.__str__())
                    resultVector = worker.run()
                    resultMessage = CorrelationVectorMessage(mpi_myid,0)
                    resultMessage.setVector(resultVector)
                    
                    #worker.dumpMsg2Log('node'+mpi_myid.__str__()+'.log', resultMessage.__str__())
                    if verbose and False:
                        print(resultMessage)
                        
                    pytom_mpi.send(resultMessage.__str__(),0)
                    
                except(MessageError,RuntimeError,IndexError):
                    msg = StatusMessage('', '')
                    msg.fromStr(mpi_msg)
                    if msg.getStatus() == 'End':
                        end = True
            print('Node ' + mpi_myid.__str__() + ' finished')
    else:
        print('Sequential Processing! Running on one machine only!')
        manager = CMManager(job)
        
        manager.calculateMatrix()
        
        manager.saveMatrix()
    
    pytom_mpi.finalise()
        
def decomposeMatrix(matrixFile):
    """
    decomposeMatrix: Will perform analysis on calculated matrix using scipy.linalg tools
    @param matrixFile: Correlation matrix in em format
    @type matrixFile: Either the matrix as volume or a string   
    @return: [eigenValues,eigenVectors] of the matrix
    """
    
    from scipy import linalg
    from numpy import matrix
    
    if matrixFile.__class__ == str:
        from pytom_volume import read
        corrMatrix = read(matrixFile)
    else:
        corrMatrix = matrixFile 
        
    matStr = '' #copy all matrix values into string
    try:
        #try to use shortest path by converting from tom volume to numpy matrix
        from pytom_numpy import vol2npy
        
        mat = matrix(vol2npy(corrMatrix))
        
    except ImportError:
        #numpy did not work. copy each entry into string and proceed, takes long time
        
        for x in range(corrMatrix.sizeX()):
            for y in range(corrMatrix.sizeY()):
                matStr = matStr + ' ' + str(corrMatrix.getV(x,y,0))
            if x < (corrMatrix.sizeX() -1):      
                matStr = matStr + ';'
    
            mat = matrix(matStr) #generate numpy.matrix from string

    return linalg.eigh(mat) #determine values
    

def scaleEigenVectors(eigenValues,eigenVectors):
    """
    scaleEigenVectors: Scales eigenvectors according sqrt(eigenvalue * eigenvector) = eigenvector
    @param eigenValues:
    @param eigenVectors:
    @return:  scaled eigenVectors of the matrix 
    """
    from numpy import sqrt
    
    for i in range(len(eigenValues)):
        eigenVectors[i] = sqrt(abs(eigenValues[i])) * eigenVectors[i]
    
    return eigenVectors
    
    
    
    
def clusterMatrix(matrix,clusterAlgorithm='Hierarchical',numberClusters = 5,numberEigenVectors=None,verbose=False):
    """
    clusterMatrix: Will perform cluster analysis on correlation matrix
    @param matrix: The correlation matrix
    @param clusterAlgorithm: The clustering algorithm. 'Hierarchical' or 'K-Means'. Hierarchical is default
    @param numberClusters: Into how many clusters do you want to split the data? 
    @param numberEigenVectors: Number < len(EigenVectors) . len(EigenVectors) is used otherwise 
    @param verbose: Show plots (if matplotlib is installed). False is default  
    @return: Cluster affiliation for each original observation. 
    @rtype: List. Length = number original observations 
    """
    
    [eigenValues,eigenVectors] = decomposeMatrix(matrix)
    
    numberEigenVectors = numberEigenVectors or len(eigenValues) 
    
    eigenVectors = scaleEigenVectors(eigenValues,eigenVectors)
    
    if clusterAlgorithm == 'Hierarchical':
        from scipy.cluster.hierarchy import  average,dendrogram,fcluster
        from scipy.spatial.distance import pdist
        
        distances = pdist(eigenVectors)
        if verbose:
            print(distances)
    
        linkage = average(distances)
        if verbose:
            print(linkage)
        
        if verbose:
            dendrogram(linkage)
            
        return fcluster(linkage,50,criterion='maxclust',depth=2)
         
    elif clusterAlgorithm == 'K-Means':
        
        from scipy.cluster.vq import kmeans2
        
        #reduce array size according to numberEigenVectors    
        print(eigenVectors[:,0:numberEigenVectors])
        #perform kmeans
        [centroids,clusterLabels] = kmeans2(eigenVectors[:,0:numberEigenVectors],numberClusters,10,'random')
        
        return clusterLabels
    