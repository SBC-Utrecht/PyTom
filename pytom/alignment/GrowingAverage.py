'''
Created on May 3, 2010

@author: hrabe
'''
from pytom.basic.structures import PyTomClass


class GAWorker(PyTomClass):
    """
    GAWorker: Processes growing average job.
    """
    def __init__(self,mpiid):
        self._mpiId = mpiid

    def dumpMsg2Log(self,logfile,msg):
        """
        dumpMsg2Log:
        @param logfile:
        @param msg:
        @author: Thomas Hrabe 
        """
        from pytom.tools.files import dump2TextFile
        
        dump2TextFile(logfile,msg)

    def run(self):
        
        import pytom.lib.pytom_mpi as pytom_mpi
        from pytom.parallel.messages import StatusMessage,MessageError
        from pytom.basic.exceptions import ParameterError
        from pytom.basic.structures import PyTomClassError
        
        if not pytom_mpi.isInitialised():
            pytom_mpi.init()

        end = False

        while not end:
            try:
                mpi_msgString = pytom_mpi.receive()
                msg = GrowingAverageJobMessage('','')
                msg.fromStr(mpi_msgString)
                
                
                self.fromJob(msg.getJob())
                self._run()
                
                returnMessage = GrowingAverageResultMessage(self._mpiId,0)
                
                pytom_mpi.send(returnMessage.__str__(),0)
                
            except (MessageError,PyTomClassError,ParameterError,IndexError):
                        #as StatusMessage and finish
                        msg = StatusMessage('','')
                        msg.fromStr(mpi_msgString)
                        if msg.getStatus() == 'End':
                            end = True


    def _run(self):
        """
        run: Starts growing average alignment
        """
        
        from pytom.basic.structures import Reference
        from pytom.lib.pytom_volume import read
        from pytom_fftplan import fftShift
        
        #create reference object - as self.reference and weighting on disk
        reference = self._particleList[self._startParticleNumber]
        referenceFile = reference.getFilename()

        refVolume = read(referenceFile);
        
        #parse filename
        pos1 = referenceFile.rfind('/')
        pos2 = referenceFile.rfind('.em')
        
        r = referenceFile[pos1+1:pos2]
        
        #self._reference = Reference(referenceFile,r + '-StartWeight.em')
        self._reference = Reference(referenceFile)
        
        wedgeInfo = reference.getWedgeInfo()
        wedgeVolume = wedgeInfo.returnWedgeVolume(refVolume.sizeX(),refVolume.sizeY(),refVolume.sizeZ(),False)
        
        #wedgeVolume.setAll(1);
        wedgeVolume.write(r + '-StartWeight.em')
        
        self._startGrowingAverageLoop()
        
    def _startGrowingAverageLoop(self):
        """
        _startGrowingAverageLoop: Loops over all particles in self.particleList and computes growing average
        """
        
        from pytom.alignment.structures import GrowingAverageInterimResult
        from pytom.basic.structures import Rotation,Shift
        from pytom.alignment.alignmentFunctions import alignTwoVolumes
        numberParticles = len(self._particleList)
        
        self._resultList = []
        
        result = GrowingAverageInterimResult(self._particleList[self._startParticleNumber],self._reference,Rotation(0,0,0),Shift(0,0,0),self._score)
    
        #init result list with first element 
        self._resultList.append(result)
        self._createNewAverage(self._startParticleNumber)
        
        for particleIterator in range(numberParticles):
            
            if particleIterator == self._startParticleNumber:
                continue
                    
            print('Running particle no ' + str(particleIterator)) 
                    
            particle = self._particleList[particleIterator]
            
            #determine alignment of current reference with particle
            result = alignTwoVolumes(particle,self._reference,self._angleObject,self._maskFile,self._score,self._preprocessing)
            
            self._angleObject.reset()
            
            self._resultList.append(result)
            
            self._createNewAverage(particleIterator)
            
    def _createNewAverage(self,iteration):
        """
        _createNewAverage:
        @param iteration: 
        """
        from pytom.alignment.alignmentFunctions import average
        filename = self._destinationDirectory + str(iteration) + '.em'
        self._reference = average(self._resultList,filename)
        
        self._resultList[len(self._resultList)-1].toXMLFile(filename+'.xml')
            
            
    def fromJob(self,job):
        """
        fromJob: Initializes this object.
        @param job: 
        @type job: L{pytom.alignment.structures.GrowingAverageJob}   
        """
        
        from pytom.alignment.structures import GrowingAverageJob
        
        if not job.__class__ == GrowingAverageJob:
            
            from pytom.basic.exceptions import ParameterError
            raise ParameterError('You must provide a GrowingAverageJob object to initialize GrowingAverage.')
        
        self._particleList = job.particleList
        self._angleObject = job.angleObject
        self._startParticleNumber = job.startParticleNumber
        self._maskFile = job.maskFile
        self._score = job.score
        self._destinationDirectory = job.destinationDirectory
        self._preprocessing = job.preprocessing
        
class GAManager(PyTomClass):
    
    def __init__(self,particleClassLists,score,angleObject,mask,destinationDirectory,preprocessing):
        """
        @param particleClassLists: List ( [] ) of particleLists
        @param score: Scoring object
        @param angleObject: Stores all scanning angles
        @param mask: Scoring mask 
        @param destinationDirectory: Directory to save results to     
        """
        self.particleClassLists = particleClassLists
        self.score = score
        self.angleObject = angleObject
        self.mask = mask
        self.destinationDirectory = destinationDirectory
        self.preprocessing = preprocessing
        
    
    def sequentialGA(self,verbose=False):
        """
        sequentialGA: Sequential growing average procedure
        @todo: not tested yet
        """
        import os
        from pytom.tools.ProgressBar import FixedProgBar
        
        jobCounter = 0
        
        classIterator = 0
        numberClasses = len(self.particleClassLists)
        
        progressBar = FixedProgBar(0,len(self.particleClassLists),'Jobs finished ')
        progressBar.update(0)
        numberJobsDone = 1

        
        
        for classIterator in range(numberClasses):
            #distribute Jobs, make sure all available nodes are working 
            

            print('Sending class ' +classIterator.__str__())
            
            os.system('mkdir ' + self.destinationDirectory + '/class' + classIterator.__str__())
            
            job = GrowingAverageJob(self.particleClassLists[classIterator],self.angleObject,self.mask,self.score,0,self.destinationDirectory+ '/class' + str(classIterator),self.preprocessing)
                
            if verbose:
                print(job)
            
            
            worker = GAWorker(-1)
            worker.fromJob(job)
            worker._run()
            
            #update progress bar
            progressBar.update(numberJobsDone)
            numberJobsDone = numberJobsDone +1
            
                
def growingAverage(particleClassLists,score,angleObject,mask,destinationDirectory,preprocessing,verbose=False):
    
    import pytom.lib.pytom_mpi as pytom_mpi
    
    
    if not pytom_mpi.isInitialised():
        pytom_mpi.init()
    
    
    
    if pytom_mpi.size() > 1:
        
        mpi_myid = pytom_mpi.rank()
    
        if mpi_myid == 0:
            
            
            manager = GAManager(particleClassLists,score,angleObject,mask,destinationDirectory,preprocessing)
            
            manager.parallelGA(verbose)
            manager.parallelEnd()
        
        else:
            w = GAWorker(mpi_myid)
            w.run()
    
    else:
        print('Processing in sequential mode')
        
        manager = GAManager(particleClassLists,score,angleObject,mask,destinationDirectory,preprocessing)
        
        manager.sequentialGA(verbose)
        
        
    pytom_mpi.finalise()
    
    
def growingAverageNew(particleList=None,angleObject=None,maskFile=None,scoreObject=None,startClassNumber=0,destinationDirectory='.',preprocessing = None,binning=1,verbose=False):
    """
    
    """
    
    from pytom.alignment.alignmentFunctions import bestAlignment
    from pytom.basic.structures import Reference,Particle,Rotation,ParticleList
    from pytom.alignment.preprocessing import Preprocessing
    
    if not preprocessing:
        preprocessing = Preprocessing()
    
    
    numberOfClasses = len(particleList.splitByClass())
    if verbose:
        print('Processing ' + str(numberOfClasses) + ' classes.')
        print('Generating start average')
        
    startAverageList = particleList.particlesFromClass(float(startClassNumber))
    
    startAverageList.average(destinationDirectory + '/GA_it0.em',progressBar=verbose)
    
    currentReference = Reference(destinationDirectory + '/GA_it0.em')
    
    growingAverageParticleList = ParticleList(particleList.getDirectory())
    
    for p in startAverageList:
        p.setRotation(Rotation(0,0,0))
        
        growingAverageParticleList.append(p)
    
    for i in range(2,numberOfClasses):
        
        currentParticleList = particleList.particlesFromClass(float(i))
        
        if verbose:
            print('Generating ' + str(i) + '. class average')
            
        currentParticleList.average(destinationDirectory + '/CA_it'+str(i)+'.em',progressBar=verbose)
        
        currentParticle = Particle(destinationDirectory + '/CA_it'+str(i)+'.em',wedgeInfo=currentParticleList[0].getWedgeInfo())
        
        if verbose:
            print('Running alignment iteration ' + str(i))
            print(currentParticle)
            print(currentReference)    
        
        currentPeak = bestAlignment(currentParticle.getVolume(),currentReference.getVolume(),currentReference.getWeighting(),currentParticle.getWedgeInfo(),angleObject,scoreObject,maskFile,preprocessing=preprocessing,binning=binning)
        
        if verbose:
            print('Parameters determined:')
            print(currentPeak)
            
        for p in currentParticleList:
            p.setRotation(currentPeak.getRotation())
            p.setShift(currentPeak.getShift())
            
            growingAverageParticleList.append(p)
        
        if verbose:
            print('Generating growing average ' + str(i))
            
        growingAverageParticleList.average(destinationDirectory + '/GA_it'+ str(i) +'.em',progressBar=verbose)
        
        currentReference = Reference(destinationDirectory + '/GA_it'+ str(i) +'.em')
        angleObject.reset()
        
class GrowingAverageJob(PyTomClass):
    """
    GrowingAverageJob:
    @ivar particleList: List of particles to be aligned
    @ivar angleObject: Angle object L{pytom.angles.AngleObject}
    @ivar startParticleNumber: Number of start particle (default 0)
    @ivar maskFile: Mask used for appedizing 
    """
    
    def __init__(self,particleList=None,angleObject=None,maskFile=None,scoreObject=None,startClassNumber=0,destinationDirectory='.',preprocessing = None):
        
        from pytom.tools.files import checkDirExists
        from pytom.angles.angleList import AngleList
        from pytom.basic.structures import ParticleList
        
        self._particleList = particleList or ParticleList('/')
        self._angleObject = angleObject or AngleList()
        self._startClassNumber = startClassNumber
        self._maskFile = maskFile or None
        
        if preprocessing:
            self._preprocessing = preprocessing
        else:
            from pytom.alignment.preprocessing import Preprocessing
            self._preprocessing = Preprocessing()
        
        
        if self._maskFile.__class__ == str:
            from pytom.basic.structures import Mask 
            self._maskFile = Mask(self.maskFile)
            
        self._score = scoreObject
        
        if not checkDirExists(destinationDirectory):
            raise Exception('Destination directory ' + destinationDirectory + ' does not exist.')
        
        if not destinationDirectory[len(destinationDirectory)-1] == '/':
            destinationDirectory = destinationDirectory + '/'
        
        self._destinationDirectory = destinationDirectory
        
    def toXML(self):
        """
        toXML : Compiles a XML file from result object
        rtype : L{lxml.etree._Element}
        @author: Thomas Hrabe
        """ 
        from lxml import etree
        
        job_element = etree.Element('GrowingAverageJob',StartClassNumber = str(self._startClassNumber), DestinationDirectory = str(self._destinationDirectory))
        
        job_element.append(self._maskFile.toXML())
        
        job_element.append(self._particleList.toXML())
        
        job_element.append(self._angleObject.toXML())
        
        job_element.append(self._score.toXML())
        
        job_element.append(self._preprocessing.toXML())
        
        return job_element
        
    def fromXML(self,xmlObj):
        """
        fromXML:
        @param xmlObj: A xml object  
        @type xmlObj: L{lxml.etree._Element}
        @author: Thomas Hrabe 
        """
        
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element :
            from pytom.basic.exceptions import ParameterError
            raise ParameterError('Is not a lxml.etree._Element! You must provide a valid XMLobject.')
        
        if xmlObj.tag == 'GrowingAverageJob':
            job_element = xmlObj
        else:
            from pytom.basic.exceptions import ParameterError
            raise ParameterError('Is not a GrowingAverageJobXML! You must provide a valid GrowingAverageJobXML object.')
        
        from pytom.angles.angle import AngleObject
        from pytom.basic.score import fromXML as scoreFromXML
        from pytom.basic.structures import ParticleList
        from pytom.alignment.preprocessing import Preprocessing
        
        self._startClassNumber = int(job_element.get('StartClassNumber'))
        
        mask = job_element.xpath('Mask')[0]
        from pytom.basic.structures import Mask
        self._maskFile = Mask('')
        self._maskFile.fromXML(mask)
         
        self._destinationDirectory = job_element.get('DestinationDirectory')
        
        particleXML = job_element.xpath('ParticleList')[0]
        
        self._particleList = ParticleList('/',[])
        self._particleList.fromXML(particleXML)
        
        angleXML = job_element.xpath('Angles')[0]
        
        ang = AngleObject()
        self._angleObject = ang.fromXML(angleXML)
        
        scoreXML = job_element.xpath('Score')[0]
        
        self._score = scoreFromXML(scoreXML)
        
        self._preprocessing = Preprocessing()
        preprocessingXML = job_element.xpath('Preprocessing')[0]
        self._preprocessing.fromXML(preprocessingXML)
        
class GrowingAverageInterimResult(PyTomClass):
    """
    GrowingAverageInterimResult:
    """
    
    
    def __init__(self,particle,reference,rotation,shift,score):
        self.particle = particle
        self.reference = reference
        self.rotation = rotation
        self.shift = shift
        self.score = score
    
    def getFilename(self):
        return self.particle.getFilename()
      
    def getWedgeInfo(self):
        return self.particle.getWedgeInfo()
    
    def getRotation(self):
        return self.rotation
    
    def getShift(self):
        return self.shift
    
    def toXML(self):
        """
        toXML : Compiles a XML file from result object
        rtype : L{lxml.etree._Element}
        @author: Thomas Hrabe
        """ 
        
        from lxml import etree
        
        result_element = etree.Element('GrowingAverageInterimResult')
        
        result_element.append(self.particle.toXML())
        result_element.append(self.reference.toXML())
        
        result_element.append(self.rotation.toXML())
        result_element.append(self.shift.toXML())
        
        result_element.append(self.score.toXML())
        return result_element
    
    def fromXML(self,xmlObj):
        """
        fromXML:
        @param xmlObj: A xml object  
        @type xmlObj: L{lxml.etree._Element}
        @author: Thomas Hrabe 
        """ 
        
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element :
            raise Exception('Is not a lxml.etree._Element! You must provide a valid XMLobject.')
        
        if xmlObj.tag == 'GrowingAverageInterimResult':
            result_element = xmlObj
        else:
            Exception('XML object is not a GrowingAverageInterimResult! You must provide a valid GrowingAverageInterimResultXML object.')
        
        from pytom.basic.structures import Particle,Reference
        
        particleXML = result_element.xpath('/GrowingAverageInterimResult/Particle')[0]
        
        self.particle = Particle('')
        self.particle.fromXML(particleXML)
        
        referenceXML = result_element.xpath('/GrowingAverageInterimResult/Result')[0]
        self.reference = Reference('')
        self.reference.fromXML(referenceXML)
        

from pytom.parallel.messages import Message
       
class GrowingAverageJobMessage(Message):

    def setJob(self,job):
        self.job = job;
        
    def getJob(self):
        return self.job;
    
    def fromXML(self,xmlObj):
        """
        fromXML : Assigns values to job attributes from XML object
        @param xmlObj: A xml object  
        @author: Thomas Hrabe 
        """
        from lxml.etree import _Element,tostring;
        
        if xmlObj.__class__ != _Element :
            
            from pytom.basic.exceptions import ParameterError
            raise ParameterError('Is not a lxml.etree._Element! You must provide a valid XML-MaximisationJobMsg object.')
            
        main = xmlObj.xpath('/GrowingAverageJobMsg')
    
        message = main[0].xpath('Message')
        message = message[0]
        
        self._sender = int(message.get('Sender'))
        self._recipient = int(message.get('Recipient'))
        self._timestamp = message.get('Timestamp')
        
        jobDescription = main[0].xpath('GrowingAverageJob')
        
        j = GrowingAverageJob()
        j.fromXML(jobDescription[0])
        self.setJob(j)
    
    def toXML(self):
        """
        toXML : Compiles a XML file from object
        @author: Thomas Hrabe
        """    
        from lxml import etree
        
        jobMsg = etree.Element("GrowingAverageJobMsg")
        
        messageElement = etree.Element("Message",Sender = self._sender.__str__(), Recipient = self._recipient.__str__(), Timestamp = self._timestamp.__str__())
        
        jobMsg.append(messageElement)
        
        jobElement = self.job.toXML()
        jobMsg.append(jobElement)
        
        return jobMsg

class GrowingAverageResultMessage(Message):
    
    def fromXML(self,xmlObj):
        """
        fromXML : Assigns values to job attributes from XML object
        @param xmlObj: A xml object  
        @author: Thomas Hrabe 
        """
        from lxml.etree import _Element,tostring
        
        if xmlObj.__class__ != _Element :
            
            from pytom.basic.exceptions import ParameterError
            raise ParameterError('Is not a lxml.etree._Element! You must provide a valid XML-MaximisationJobMsg object.')
        
        main = xmlObj.xpath('/GrowingAverageResultMsg')
    
        message = main[0].xpath('Message')
        message = message[0]
    
        self._sender = int(message.get('Sender'))
        self._recipient = int(message.get('Recipient'))
        self._timestamp = message.get('Timestamp')

    
    def toXML(self):
        """
        toXML : Compiles a XML file from object
        @author: Thomas Hrabe
        """    
        from lxml import etree
        
        resultMsg = etree.Element("GrowingAverageResultMsg")
        
        messageElement = etree.Element("Message",Sender = self._sender.__str__(), Recipient = self._recipient.__str__(), Timestamp = self._timestamp.__str__())
        
        resultMsg.append(messageElement)
        
        return resultMsg