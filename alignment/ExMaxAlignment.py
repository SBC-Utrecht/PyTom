'''
Created on Apr 19, 2010

@author: hrabe
'''

from pytom.basic.structures import PyTomClass

class ExMaxJob(PyTomClass):
    """
    ExMaxJob : stores all settings needed for Expectation Maximisation based alignment
    """    
    def __init__(self,particleList=None,destination=None,reference=None,score=None,
                 rotations=None,mask=None,symmetry=None,
                 numberRefinementRounds=None,numberIterations=None,preprocessing=None,
                 excludeThreshold=-1,binning=1,sampleInformation=None,fscCriterion=0.5,
                 adaptiveResolution=True,adaptiveOffset=0.1,angleFactor=0.5):
        """
        __init__ :
        @param particleList:
        @type particleList: L{pytom.basic.structures.ParticleList}
        @param destination:
        @type destination: str
        @param reference: 
        @type reference: L{pytom.basic.structures.Reference}
        @param score:
        @type score: L{pytom.score.score.Score}
        @param rotations:     
        @type rotations: L{pytom.angles.angle}
        @param mask:
        @type mask: str
        @param symmetry:
        @type symmetry:  Child of L{pytom.basic.structures.Symmetry}
        @param numberRefinementRounds: not used (why here anyway ....?)
        @param numberIterations: 
        @param preprocessing: 
        @type preprocessing: L{pytom.alignment.preprocessing.Preprocessing}
        @param excludeThreshold: Include only the N% best particles for further processing 
        @type excludeThreshold:  int
        @param binning: Integer determining the binning (libtomc notation) 
        @param sampleInformation:
        @type sampleInformation: L{pytom.basic.SampleInformation}
        @param fscCriterion: Fourier shell criterion for determining the resolution. Define at which value resolution will be determined.0.5 is default.
        @type fscCriterion: float
        @param adaptiveResolution: Apply bandpass filter to current resolution and adjust refinement angle by current resolution? Default is True 
        @param adaptiveOffset: Offset for bandpass to include more signal in processing. Default is 0.1, can be set to anything. 
        @param angleFactor: Scaling of adaptive angle for angular sampling. Default is 0.5. new_angle = angleFactor * resolution_angle   
        @todo: Add default parameter support here!
        """
        from pytom.basic.structures import ParticleList,SampleInformation,PointSymmetry
        from os import sep as pathSeperator
        
        if particleList.__class__ == str:
            self._particleList = ParticleList(particleList,[])
        else:
            self._particleList = particleList
            
        self._destination = destination
        
        if self._destination.__class__ == str and not self._destination[len(self._destination)-1] == pathSeperator:
            self._destination =  self._destination + pathSeperator
        
        self._reference = reference
        self._score = score
        self._rotations = rotations
        
        if mask.__class__ == str:
            from pytom.basic.structures import Mask
            self._mask = Mask(mask)
        else:
            self._mask = mask

        self._symmetry = symmetry if symmetry is not None else PointSymmetry(1)
        self._numberIterations = numberIterations
        self._preprocessing = preprocessing
        self._excludeThreshold = excludeThreshold
        self._binning = binning
        if not sampleInformation:
            self._sampleInformation = SampleInformation(1,1)
        else:
            self._sampleInformation = sampleInformation
    
        self._fscCriterion = fscCriterion
        self._adaptiveResolution = adaptiveResolution
        self._adaptiveOffset = adaptiveOffset
        self._angleFactor = angleFactor 
    
    def getBinning(self):
        return self._binning
    
    def getClassificationParameter(self):
        '''
        @deprecated: Use getExcludeThreshold instead
        '''
        return self._excludeThreshold
    
    def getExcludeThreshold(self):
        return self._excludeThreshold
    
    def getParticleList(self):
        return self._particleList
    
    def getParticlePath(self):
        return self._particlePath
    
    def getDestination(self):
        return self._destination

    def setDestination(self,destination):
        self._destination = destination
        
    def getReference(self):
        return self._reference
    
    def setReference(self,reference):
        self._reference = reference
        
    def getRotations(self):
        return self._rotations
    
    def getMask(self):
        return self._mask
    
    def getScore(self):
        return self._score

    def getNumberIterations(self):
        return self._numberIterations
    
    def getPreprocessing(self):
        return self._preprocessing
    
    def getFSCCriterion(self):
        return self._fscCriterion
    
    def getAdaptiveResolution(self):
        return float(self._adaptiveResolution)
    
    def getAdaptiveOffset(self):
        return float(self._adaptiveOffset)
    
    def getAngleFactor(self):
        return float(self._angleFactor)
    
    def getSampleInformation(self):
        return self._sampleInformation
    
    def setAngleFactor(self,factor):
        
        if not factor.__class__ == float:
            raise RuntimeError('Parameter for angle factor must be a float!')
        
        if not (0 <= factor <= 1):
            raise RuntimeError('Angle factor must be 0 <= factor <= 1')
        
        self._angleFactor = factor
        
    def setAdaptiveOffset(self,offset):
        
        if not offset.__class__ == float:
            raise RuntimeError('Adaptive offset must be a float!')
        
        if offset < 0 or offset > 1:
            raise RuntimeError('Adaptive offset must be in [0,1]')
        
        self._adaptiveOffset = offset
        
    def setFSCCriterion(self,fscCriterion):
        self._fscCriterion = fscCriterion
        
    def setAdaptiveResolution(self,adaptiveResolution):
        self._adaptiveResolution = adaptiveResolution
        
    def setPreprocessing(self,preprocessing = None):
        from pytom.alignment.preprocessing import Preprocessing
        if not preprocessing or preprocessing.__class__ != Preprocessing:
            return 
        self._preprocessing = preprocessing
        
    def getSymmetry(self):
        return self._symmetry
    
    def setRotations(self,rotations):
        self._rotations = rotations
    
    def setParticleList(self,particleList):
        self._particleList = particleList

    def setFSCOffset(self,value):
        self._fscOffset = float(value)
        
    def fromXML(self,xmlObj=-1):
        """
        fromXML : Assigns values to job attributes from XML object
        @param xmlObj: A xml object  
        @type xmlObj: L{lxml.etree._Element}
        @author: Thomas Hrabe 
        """
        from lxml.etree import _Element,tostring
        
        if xmlObj.__class__ != _Element:
            raise TypeError('You must provide a valid XML-ExpectationMaximisationJob object.')
        
        from pytom.angles.angle import AngleObject
        from pytom.score.score import fromXML as fromXMLScore
        from pytom.alignment.preprocessing import Preprocessing
        from pytom.basic.structures import ParticleList,Symmetry,Reference,SampleInformation
        
        if xmlObj.tag == 'ExpectationMaximisationJob':
            job = xmlObj
        else:
            job = xmlObj.xpath('ExpectationMaximisationJob')
            
            if len(job) == 0:
                raise TypeError("This XML is not an ExpectationMaximisationJob.")
            
            job = job[0]
         
        jobDescription = job.xpath('Description')    
        jobDescription = jobDescription[0]

        particleList_element = jobDescription.xpath('ParticleList')[0]
        pl = ParticleList('/',[])
        pl.fromXML(particleList_element)
        self._particleList = pl
        
        self._destination = jobDescription.get('Destination')
        if jobDescription.get('ExcludeThreshold') is None:
            self._excludeThreshold = int(float(jobDescription.get('ResultClassification')) * 100)
        else:
            self._excludeThreshold = int(float(jobDescription.get('ExcludeThreshold')))
        
        
        r = jobDescription.xpath('Reference')    
        self._reference = Reference('')
        self._reference.fromXML(r[0])
    
        mask_element = jobDescription.xpath('Mask')[0]
        from pytom.basic.structures import Mask
        self._mask = Mask()
        self._mask.fromXML(mask_element)
        
        self._numberIterations = int(jobDescription.get('NumberIterations'))
        self._binning = int(jobDescription.get('Binning'))
        
        if jobDescription.get('FSCCriterion') == None:
            self._fscCriterion = 0.5
        else:
            self._fscCriterion = float(jobDescription.get('FSCCriterion'))
        
        if jobDescription.get('AdaptiveResolution') == None:
            self._adaptiveResolution = True
        else:
            self._adaptiveResolution = jobDescription.get('AdaptiveResolution') == 'True'
            
        if jobDescription.get('AdaptiveOffset') == None:
            self._adaptiveOffset = 0.1
        else:
            self._adaptiveOffset = float(jobDescription.get('AdaptiveOffset'))
        
        if jobDescription.get('AngleFactor') == None:
            self._angleFactor = 0.5
        else:
            self._angleFactor = float(jobDescription.get('AngleFactor'))
        
        score = jobDescription.xpath('Score')
        self._score = fromXMLScore(score[0])
        
        angles = jobDescription.xpath('Angles')
        ang = AngleObject()
        self._rotations = ang.fromXML(angles[0])
        
        symmetryXML = jobDescription.xpath('Symmetry')
        symmetryFactory = Symmetry()
        self._symmetry  = symmetryFactory.fromXML(symmetryXML[0])
           
        p = Preprocessing()
        preObj = jobDescription.xpath('Preprocessing')
        p.fromXML(preObj[0])
        p.setTaper( taper=0) # later set hard-coded to dim/10
        self._preprocessing = p
        
        siObj = jobDescription.xpath('SampleInformation')
        if len(siObj) > 0:
            si = SampleInformation()
            si.fromXML(siObj[0])
            self._sampleInformation = si
        else:
            self._sampleInformation = SampleInformation()
        
    def toXML(self):
        """
        toXML : Compiles a XML object from job object
        @rtype: L{lxml.etree._Element}
        @return: XML Object
        @author: Thomas Hrabe
        """
        from lxml import etree

        jobElement = etree.Element("ExpectationMaximisationJob")
        descriptionElement = etree.Element("Description")
        
        descriptionElement.append(self._particleList.toXML())
        descriptionElement.set("Destination",self._destination)
        
        descriptionElement.append(self._reference.toXML())

        descriptionElement.append(self._mask.toXML())
        
        descriptionElement.set("NumberIterations",str(self._numberIterations))
        descriptionElement.set("ExcludeThreshold",str(self._excludeThreshold))
        descriptionElement.set("Binning", str(self._binning))
        descriptionElement.set("FSCCriterion",str(self._fscCriterion))
        descriptionElement.set("AdaptiveResolution",str(self._adaptiveResolution))
        descriptionElement.set("AdaptiveOffset",str(self._adaptiveOffset))
        descriptionElement.set("AngleFactor",str(self._angleFactor))
        
        descriptionElement.append(self._score.toXML())
        
        descriptionElement.append(self._rotations.toXML())
        
        descriptionElement.append(self._preprocessing.toXML())
        descriptionElement.append(self._symmetry.toXML())
        
        descriptionElement.append(self._sampleInformation.toXML())
        
        jobElement.append(descriptionElement)
        
        assert jobElement.__class__ == etree._Element
        return jobElement
        
    
    def toHTMLFile(self,filename):
        """
        toHTMLFile: Overrides parent method and stores ExpectationMaximisationJob to HMTL
        @param filename: HTML filename
        """
        from pytom.tools.files import getPytomPath

        super(self.__class__,self).toHTMLFile(filename,xsltFile=getPytomPath() + '/xslt/ExpMaxJob.xsl')


    def check(self):
        """
        check: Performs logical check on self to make sure that job settings are correct. Will display error message if not.
        """
        
        from pytom.tools.files import checkFileExists,checkDirExists
        from pytom_volume import read
        
        self._particleList.check()
        
        self._mask.check()
        
        self._reference.check()
        
        if not checkDirExists(self._destination):
            raise ValueError('Result destination path not found! ' + self._destination)
        

class ExMaxWorker(object):        
    """
    ExMaxWorker: Will perform all operations neccessary to align two particles. 
    I.e, determines the orientation and shift of a particle according to a reference for a certain number of rotations.
    It is driven mainly by job objects defining the next operation. Accepted job structures are L{pytom.alignment.structures.MaximisationJob} or L{pytom.alignment.structures.ExpectationJob} 
    @author: Thomas Hrabe
    """
    
    def run(self):
        """
        run
        """
        from pytom.alignment.structures import MaximisationJob,ExpectationJob
        #trying to resolve old naming twist: change expectation and maximization naming
        if self._jobType == MaximisationJob:
            #print 'exMax worker run: expectation'
            return self._expectation()
        elif self._jobType == ExpectationJob:
            #print 'exMax worker run: maximization'
            return self._maximisation()
    
    def parallelRun(self,verbose=False):
        """
        parallelRun: Run the worker in parallel mode
        @param verbose: Print debug messages (False by default) 
        """
    
        import pytom_mpi
        from pytom.parallel.alignmentMessages import MaximisationJobMsg,MaximisationResultMsg,ExpectationJobMsg,ExpectationResultMsg
        from pytom.parallel.messages import StatusMessage,MessageError
        from pytom.basic.exceptions import ParameterError
        from pytom.basic.structures import PyTomClassError
        from pytom.tools.files import dumpMsg2Log,checkFileExists
        import os
        
        end = False
        
        if not pytom_mpi.isInitialised():
            pytom_mpi.init()
        
        mpi_id = pytom_mpi.rank()
        
        mpi_myname = 'node_' + str(mpi_id)
        
        logFileName = mpi_myname+'.log'
        
        while not end:
            
            #listen for messages
            mpi_msgString = pytom_mpi.receive()
            
            if verbose:
                print('Worker received:')
                print(mpi_msgString)
                dumpMsg2Log(logFileName,mpi_msgString)
            
            #parse message
            try:
                #as job and start processing
                
                msg = MaximisationJobMsg()
                msg.fromStr(mpi_msgString)
                
                self.fromJob(msg.getJob())
                
                result = self.run()
                resultMsg = MaximisationResultMsg(str(mpi_id),'0')
                resultMsg.setResult(result)
                
                pytom_mpi.send(str(resultMsg),0)
                
                if verbose:
                    dumpMsg2Log(logFileName,str(resultMsg))
                    
            except (MessageError,PyTomClassError,ParameterError,TypeError):
                try:
                    msg = ExpectationJobMsg()
                    msg.fromStr(mpi_msgString)

                    self.fromJob(msg.getJob())
                    
                    result = self.run()
                    
                    resultMsg = ExpectationResultMsg(str(mpi_id),'0')                    
                    resultMsg.setResult(result)                    
                    pytom_mpi.send(str(resultMsg),0)
                    
                    if verbose:
                        dumpMsg2Log(logFileName,str(resultMsg))
                    
                except (MessageError,PyTomClassError,ParameterError,TypeError):
                    try:
                        #as StatusMessage and finish
                        msg = StatusMessage('','')
                        msg.fromStr(mpi_msgString)
                        if msg.getStatus() == 'End':
                            end = True
                        
                    except :
                        print(mpi_msgString)
                        raise RuntimeError('Error parsing message. Message either unknown or invalid.')
        if verbose:              
            os.system('rm ' + logFileName)
        
        
    def _maximisation(self):
        """
        _maximisation : Performs the maximization step for particle alignment (i.e. creating a new reference)
        @rtype: L{pytom.alignment.structures.ExpectationResult}
        @return: An ExpectationResult object
        @author: Thomas Hrabe
        """ 
        from pytom.alignment.structures import ExpectationResult
        
        self._particleList.average(self._newAverageName)
        
        return ExpectationResult(self._newAverageName)
        
    def _expectation(self):
        """
        _expectation : Performs the expectation step alignment for a particle and a reference defined in the previous constructor
        @return: [endShift,endRotation] contains the final shift and final rotation        
        @author: Thomas Hrabe
        """
        from pytom.alignment.structures import MaximisationResult,Peak
        from pytom.basic.structures import Particle,Reference
        from pytom.alignment.alignmentFunctions import bestAlignment
        
        
        if self._particle.__class__ == Particle:
            from pytom_volume import read
            particleFile    = self._particle.getFilename()
            # changed FF: binning now solely done in bestAlignment
            particle        = read(particleFile,0,0,0,0,0,0,0,0,0,
                                   1, 1, 1)
                                   #self._binning, self._binning, self._binning)
        else:
            raise TypeError('Provided particle must be of type pytom.basic.structures.Particle')
         
        if self._scoreObject.__class__ == str:
            from pytom.score.score import Score
            string              = self._scoreObject
            self._scoreObject   = Score() 
            self._scoreObject.fromStr(string)
        # changed FF: mask does not need to know about binning
        #self._mask.setBinning(self._binning)
        
        # determine if the current reference is an average of the particle under scrutiny
        #if reference was generated by the current particle 
        #and
        #an unweighted particle average exists
        #and
        #a sum of all particle wedges exists
        #then subtract the particle from the current reference and proceed with the unbiased reference. 
        #read the particle as is from disk otherwise  
        if (self._reference.wasGeneratedBy(self._particle) and self._reference.hasPreWedge() and
                self._reference.hasWeighting()):
            [reference, self._referenceWeighting] = self._reference.subtractParticle(particle=self._particle, binning=1)
                #self._particle,self._binning)
            referenceObject = self._reference
        else:
            if self._reference.__class__ == Reference:
                from pytom_volume import read
                referenceObject     = self._reference
                referenceFile       = self._reference.getReferenceFilename()
                # changed FF: binning only in bestAlignment function
                reference           = read(referenceFile,0,0,0,0,0,0,0,0,0,
                                           1, 1, 1)
                                           #self._binning,self._binning,self._binning)
           
            if self._referenceWeighting == str and len(self._referenceWeighting) > 0:
                # changed FF: binning only in bestAlignment function
                from pytom_volume import read
                self._referenceWeighting        = read(self._referenceWeighting)
                #if self._binning == 1:
                #    from pytom_volume import read
                #    self._referenceWeighting        = read(self._referenceWeighting)
                #else:
                #    from pytom.basic.files import readSubvolumeFromFourierspaceFile
                #    self._referenceWeighting        = readSubvolumeFromFourierspaceFile(
                #                                    self._referenceWeightingFile,reference.sizeX(),
                #                                    reference.sizeY(),reference.sizeZ())
        pScore = self._particle.getScore()
        
        if not pScore:
            pScore = self._scoreObject
            
        # if Particle has a old and HIGH score value from localization, alignment 
        # will always take that one -> set peak = None, the angle will be scanned anyway!  
        peak = None

        self._rotationList.reset()
        peak = bestAlignment(particle=particle, reference=reference, referenceWeighting='',
                             wedgeInfo=self._particle.getWedge(), rotations=self._rotationList,
                             scoreObject=self._scoreObject, mask=self._mask, preprocessing=self._preprocessing,
                             progressBar=False, binning=self._binning, bestPeak=peak, verbose=False)
        
        if peak.__class__ != Peak:
            print('Alignment result for file ' , particleFile , ' was invalid. Using old alignment result for the next round.')
            s = pScore
            
            if not s:
                s = self._scoreObject
            peak = Peak(s.getValue(),self._particle.getRotation(),self._particle.getShift())
            
            
        self._scoreObject.setValue(peak.getScoreValue())
       
        #revert binning
        shift = peak.getShift()
        # changed FF: binning only in bestAlignment function
        #shift.scale(self._binning)

        result = MaximisationResult(self._particle, referenceObject, self._scoreObject, shift,
                                    peak.getRotation(), self._rotationList)
        
        #release object manually
        self._scoreObject._peakPrior.reset()
        return result
    
    def fromJob(self,job):
        """
        fromJob : Determines whether this object is a expectation job or maximisation job and sets the appropriate attributes. 
        For maximisation, the rotated wedge of the particle will be used if it is rotated. The global wedge will be used otherwise. 
        @param job: The job structure.
        @type job: either L{pytom.alignment.strutures.MaximisationJob} or L{pytom.alignment.strutures.ExpectationJob}  
        """
        from pytom.alignment.structures import MaximisationJob,ExpectationJob
        self._jobType = job.__class__
        
        if job.__class__ == MaximisationJob:
            self._reference = job.reference
            self._referenceWeighting = job.reference.getWeightingFilename()
            self._particle = job.particle
            self._scoreObject = job.score

            self._rotationList = job.rotations
            self._mask = job.mask
            self._numberRefinementRounds = job.numberRefinementRounds
            self._preprocessing = job.preprocessing
            self._binning = job.binning
            
        elif job.__class__ == ExpectationJob:
            self._particleList = job.getParticleList()
            self._newAverageName = job.getNewAverageName()


class ExMaxManager(PyTomClass):
    """
    ExMaxManager: Will distribute jobs among a worker pool and take care of other organisation tasks. Driven by ExpectationMaximisationJob.
    """

    def __init__(self,alignmentJob):
        """
        __init__
        """
        from pytom.alignment.ExMaxAlignment import ExMaxJob
        
        assert alignmentJob.__class__ == ExMaxJob
        
        self._fromExpMaxJob(alignmentJob)
        
        
    def _fromExpMaxJob(self,alignmentJob):
        """
        _fromExpMaxJob:
        """
        
        from pytom.alignment.structures import AlignmentList
        
        self._particleList = alignmentJob.getParticleList()
        self._reference = alignmentJob.getReference()
        self._mask = alignmentJob.getMask()
        self._numberIterations = alignmentJob.getNumberIterations()
        self._score = alignmentJob.getScore() 
        self._preprocessing = alignmentJob.getPreprocessing()
        self._symmetry = alignmentJob.getSymmetry()
        self._angleList = alignmentJob.getRotations()
        self._destination = alignmentJob.getDestination()
        self._binning = alignmentJob.getBinning()
        self._alignmentList = AlignmentList()
        self._alignmentList.clear()
        
    def getAlignmentList(self):
        return self._alignmentList
    
    def setAlignmentList(self,alignmentList): 
        from pytom.alignment.structures import AlignmentList
        
        assert alignmentList.__class__ == AlignmentList
        
        self._alignmentList = alignmentList
        
    def sequentialAlignment(self,verbose=False):
        """
        sequentialAlignment: Perform alignment on one node if no parallel environment available
        """
        
        from pytom.alignment.structures import MaximisationJob
        from pytom.tools.ProgressBar import FixedProgBar
        
        progressBar = FixedProgBar(0,len(self._particleList),'Particles aligned ')
        progressBar.update(0)
        
        for particle in self._particleList:
            
            startRotation = particle.getRotation()
            
            customisedAngles = self._angleList.setStartRotation(startRotation)
            job = MaximisationJob(particle=particle, reference=self._reference, score=self._score.copy(),
                                  rotations=customisedAngles,mask=self._mask, numberRefinementRounds=0,
                                  preprocessing=self._preprocessing, binning=self._binning)
            if verbose:
                print(job)
                
            worker = ExMaxWorker()
            worker.fromJob(job)
            result = worker.run()
            
            if verbose:
                print(result)
          
            self._alignmentList.append(result)
            progressBar.update(len(self._alignmentList))
        
    def distributeAlignment(self,verbose=False):
        """
        distributeAlignment: Distribute job in worker pool
        """
        
        import pytom_mpi
        from pytom.alignment.structures import MaximisationJob,ExpectationJob
        from pytom.tools.files import checkFileExists,checkDirExists
        from pytom.parallel.alignmentMessages import MaximisationJobMsg,MaximisationResultMsg
        from pytom.tools.ProgressBar import FixedProgBar
        
        if not pytom_mpi.isInitialised():
            pytom_mpi.init()
            
        mpi_myid = pytom_mpi.rank()
        
        if not mpi_myid == 0:
            raise RuntimeError('This function (distributeAlignment) can only be processed by mpi_id = 0! ID == ' + str(mpi_myid) +' Aborting!')
        
        if not checkDirExists(self._destination):
            raise IOError('The destination directory is not valid. ' + self._destination)
        
        if not checkFileExists(self._reference.getReferenceFilename()):
            raise IOError('The reference can not be found in the path : ' + self._reference.getReferenceFilename())
        
        mpi_numberNodes = pytom_mpi.size()
        
        particleIterator =0
        
        progressBar = FixedProgBar(0,len(self._particleList),'Particles aligned ')
        
        #distribute jobs to all nodes
        if verbose:
            print('Starting first job distribute step')
            
        for i in range(1,mpi_numberNodes):
                
            if particleIterator < len(self._particleList):
                particle = self._particleList[particleIterator]

                #determine rotation for particle
                startRotation = particle.getRotation()
                customisedAngles = self._angleList.setStartRotation(startRotation)
                
                #create job message for particle
                job = MaximisationJob(particle=particle, reference=self._reference, score=self._score,
                                      rotations=customisedAngles, mask=self._mask, numberRefinementRounds=0,
                                      preprocessing=self._preprocessing, binning=self._binning)
                job.check()
                jobMsg = MaximisationJobMsg(str(mpi_myid),str(i))
                jobMsg.setJob(job)
                
                #send job to worker
                pytom_mpi.send(str(jobMsg),i)
                if verbose:
                    print(jobMsg)
                
                particleIterator += 1
                
        finished = False
        
        if verbose:
            print('Starting job distribute and collect result step')

        progressBar.update(0)
        
        numberJobsDone = 1
        #there are more jobs than nodes. continue distributing and collect results
        while not finished:
            
            #listen and collect
            mpi_msgString = pytom_mpi.receive() 
            
            if verbose:
                print(mpi_msgString)
            
            jobResultMsg = MaximisationResultMsg('','')   
            jobResultMsg.fromStr(mpi_msgString)
            
            self._alignmentList.append(jobResultMsg.result)
            
            #send new job to free node
            if particleIterator < len(self._particleList):
                jobMsg = MaximisationJobMsg(str(mpi_myid),jobResultMsg.getSender())
                
                particle = self._particleList[particleIterator]
    
                #determine rotation for particle
                startRotation = particle.getRotation()
                customisedAngles = self._angleList.setStartRotation(startRotation)
                
                job = MaximisationJob(particle=particle, reference=self._reference, score=self._score,
                                      rotations=customisedAngles, mask=self._mask, numberRefinementRounds=0,
                                      preprocessing=self._preprocessing, binning=self._binning)
                job.check()
                
                jobMsg.setJob(job)
                
                pytom_mpi.send(jobMsg.__str__(),int(jobResultMsg.getSender()))
                
                if verbose:
                    print(jobMsg)
                    
                particleIterator = particleIterator + 1
                
            #all particles have been distributed
            finished = particleIterator >= len(self._particleList)
            #all alignment results have arrived
            finished = (len(self._particleList) == len(self._alignmentList) and finished)
    
            #update progress bar
            
            progressBar.update(numberJobsDone)
                
            numberJobsDone += 1
        
    
    def _saveForFSC(self,newReferenceName):
        """
        saveForFSC: Split particles into even / odd and calculate the Fourier Shell Correlation out of these two sets.
        @param newReferenceName: Name of current alignment result
        @type newReferenceName: str
        """
        from pytom.alignment.structures import ExpectationJob
        from pytom.alignment.aligner import ExMax
        from math import ceil
        
        exJobOdd = ExpectationJob(None, newReferenceName+'odd.em')
        exJobEve = ExpectationJob(None, newReferenceName+'even.em')
        
        resultCounter = 0
        
        results = self._alignmentList.getList()
        
        for result in results:
            
            if resultCounter % 2 == 0:
                exJobEve.appendMaximisationResult(result)
            else:
                exJobOdd.appendMaximisationResult(result)
            
            resultCounter = resultCounter + 1
        
        worker = ExMaxWorker()
        worker.fromJob(exJobEve)
        worker.run()

        ll =  exJobOdd.getParticleList()

        if len(ll)== 0:
            exJobOdd.setParticleList(exJobEve.getParticleList())
        
        worker.fromJob(exJobOdd)
        worker.run()
    
    def parallelEnd(self):
        """
        parallelEnd : Sends status message = end to all workers. All workers will terminate upon receiving this message.
        @author: Thomas Hrabe
        """
        import pytom_mpi
        from pytom.parallel.messages import StatusMessage
        
        mpi_numberNodes = pytom_mpi.size()
        mpi_myid = pytom_mpi.rank()
        
        for i in range(1,mpi_numberNodes):
            msg = StatusMessage(str(mpi_myid),str(i))
            msg.setStatus("End")
            pytom_mpi.send(str(msg),i)
    
    def determineResolution(self,filename,numberShells = None , resolutionCriterion=0.125):
        """
        determineResolution: determines resolution of current alignment. Sets the L{pytom.alignment.preprocessing} bandpass to current resolution + x.
        @param filename: Filenames for FSC 
        @param numberShells: How many shells do we use
        @param resolutionCriterion: 
        @author: Thomas Hrabe
        """
        from pytom.basic.correlation import FSC,determineResolution
        from pytom_volume import read
        from pytom.score.score import RScore
        self._saveForFSC(self._destination + filename)
        
        if self._score.__class__ == RScore:
            #do not set preprocessing bandpass for RScore. 
            return
        
        odd  = read(self._destination + filename+'odd.em')
        even = read(self._destination + filename+'even.em')
        
        if not numberShells:
            numberShells = int(odd.sizeX())
            """@ivar numberShells: max. number of shells"""
        
        fsc = FSC(odd,even,numberShells)
        
        resolution = determineResolution(fsc,resolutionCriterion)
        
        return resolution
     
                    
def parallelStart(exMaxJob,verbose,sendFinishMessage = True):
    """
    parallelStart:
    @param exMaxJob: Object defining alignment job
    @type exMaxJob: L{pytom.alignment.ExMaxAlignment.ExMaxJob}
    @param verbose: Will print debug information if True. (False by default)
    @type verbose: bool
    @param sendFinishMessage: Send finish message to worker nodes. True by default, False for multi ref alignment for instance 
    @type sendFinishMessage: bool
    """
    from pytom.basic.structures import Reference
    mpiWorks = True
    try: 
        import pytom_mpi
    
        if not pytom_mpi.isInitialised():
            if verbose:
                print('Initializing mpi')
            pytom_mpi.init()
            
    except Exception:
        mpiWorks = False
        
    if not mpiWorks or pytom_mpi.size() == 1:
        sequentialStart(exMaxJob,verbose,sendFinishMessage)
        return
    
    mpi_myid = pytom_mpi.rank()
    if pytom_mpi.size() >1:
        
        if mpi_myid == 0:
            from pytom_volume import read
            from pytom.basic.filter import lowpassFilter
            from pytom.alignment.alignmentFunctions import _disrtibuteAverageMPI,alignmentProgressToHTML
            
            from pytom.alignment.structures import AlignmentList
            from pytom.basic.resolution import bandToAngstrom,angstromToBand,angleFromResolution
            from math import ceil
            
            iteration = 1
            resultsEqual = False
        
            r = exMaxJob.getReference()
            rVol = r.getVolume()
            
            cubeSize = rVol.sizeX() 
            
            symmetry = exMaxJob.getSymmetry()
            
            while iteration <= exMaxJob.getNumberIterations() and not resultsEqual:
                #save current job
                exMaxJob.toXMLFile(exMaxJob.getDestination() + 'ExMaxJob-' + 
                                   str(iteration) + '.xml')
                
                exMaxJob.toHTMLFile(exMaxJob.getDestination() + 
                                    'ExMaxJob-' + str(iteration) + '.html')
                
                manager = ExMaxManager(exMaxJob)
                print('Distribute jobs iteration ' + str(iteration))
                manager.distributeAlignment(verbose)
                
                alignmentList = manager.getAlignmentList()
                
                #dump alignmentList and exMaxJob to file
                print('Save results iteration ' + str(iteration))                
                #alignmentList.toXMLFile(exMaxJob.getDestination() + 'AlignmentList-' + str(iteration) + '.xml')
                
                newAllParticleList = alignmentList.toParticleList()
                newAllParticleList.sortByScore('descending')
                newAllParticleList.toXMLFile(exMaxJob.getDestination() + 'ParticleList-' + str(iteration) + '.xml')
                
                #save current average
                manager.setAlignmentList(alignmentList)
                
                
                if exMaxJob.getExcludeThreshold() > 0:
                    newParticleList = newAllParticleList.listFromBestScorePercentage(exMaxJob.getExcludeThreshold())
                    #reduce particleList to N% for averaging, determine resolution. 
                    #the next round however will all particles again!
                else:
                    newParticleList = newAllParticleList
                
                print('')
                print('Write average iteration ', iteration)
                
                if symmetry.getNFold() == 1:
                    newReference = _disrtibuteAverageMPI(newParticleList,
                                                         exMaxJob.getDestination() + str(iteration) + '.em',False,verbose)
                else:
                    if verbose:
                        print('Applying ', symmetry.getNFold(), 'fold symmetry to particles')
                    
                    symmetrizedParticleList = symmetry.apply(newParticleList)
                    newReference = _disrtibuteAverageMPI(symmetrizedParticleList,
                                                         exMaxJob.getDestination() + str(iteration) + '.em',False,
                                                         verbose)
                #update reference
                
                exMaxJob.setReference(newReference)
                print('')
                print('Determine new parameters iteration ' + str(iteration))
                #set new refinement angles and update rotation strategy 
                
                sampleInfo = exMaxJob.getSampleInformation()
                if exMaxJob.getAdaptiveResolution():
                    
                    if verbose:
                        print('Determine resolution.')
                    if symmetry.getNFold() == 1:
                        #here we get the current resolution in nyquist / the band  and returns the number of bands used
                        [resNyquist,resolutionBand,numberBands]  = newParticleList.determineResolution( 
                                                                   criterion=exMaxJob.getFSCCriterion(), numberBands = cubeSize / 2, 
                                                                   mask=exMaxJob.getMask(), keepHalfsetAverages = True, 
                                                                   halfsetPrefix=exMaxJob.getDestination() + 'fsc-'+str(iteration)+'-', 
                                                                   verbose=verbose )
                    else:
                        if verbose:
                            print('Applying ', symmetry.getNFold(), 'fold symmetry to resolution')
                        [resNyquist,resolutionBand,numberBands]  = symmetrizedParticleList.determineResolution(
                                                                   criterion=exMaxJob.getFSCCriterion(), numberBands = cubeSize / 2, 
                                                                   mask=exMaxJob.getMask(), keepHalfsetAverages = True, 
                                                                   halfsetPrefix=exMaxJob.getDestination() + 'fsc-'+str(iteration)+'-', 
                                                                   verbose=verbose )
                    
                    #here we convert the current resolution to angstroms.    
                    resolutionAngstrom = bandToAngstrom( resolutionBand,
                                                         sampleInfo.getPixelSize(),numberBands,1 )
		        #sampleInfo.getPixelSize(),numberBands,exMaxJob.getBinning() )

                    if verbose:
                        print('Current resolution :' + str(resolutionAngstrom) + ' Angstrom')
                        
                    #set new angular values
                    try:
                        refinementAngle = angleFromResolution(resolutionAngstrom,sampleInfo.getParticleDiameter())
                        refinementAngle = refinementAngle * exMaxJob.getAngleFactor()
                        oldRotations = exMaxJob.getRotations() 
                        newRotations = oldRotations.focusRotation(refinementAngle=refinementAngle)
                        exMaxJob.setRotations(newRotations)
                        #update preprocessing bandpass (adaptive bandpass filter if enabled). note the  exMaxJob.getAdaptiveOffset()!!!
                        preprocessing = exMaxJob.getPreprocessing()

                        #determine better resolution (in ANGSTROM)        
                        newResolution = resolutionAngstrom * (1 - exMaxJob.getAdaptiveOffset()) 
                        #convert upper value to BAND
                        newResolutionBand = int( ceil( angstromToBand( 
                                                newResolution,sampleInfo.getPixelSize(),
                                                numberBands,1) ))
                        #update job
                        preprocessing.setHighestFrequency(newResolutionBand)
                        exMaxJob.setPreprocessing(preprocessing)
                    except RuntimeError:
                        print('Warning for adaptive angular sampling and adaptive bandpass')
                        print('Angle and filter will be unchanged due to error. Occurs when determined resolution is invalid due to too few particles.')
                        print('Current resolution :' + str(resolutionAngstrom) + ' Angstrom')
                        rotations = exMaxJob.getRotations()
                        refinementAngle = rotations.getIncrement() 
                        
                else:
                    rotations = exMaxJob.getRotations()
                    refinementAngle = rotations.getIncrement() 
                    
                    oldRotations = rotations 
                    newRotations = oldRotations.focusRotation(refinementAngle=refinementAngle)
                    exMaxJob.setRotations(newRotations)
                    
                    resolutionAngstrom = bandToAngstrom(cubeSize/2,sampleInfo.getPixelSize(),cubeSize/2,1)
                
                if exMaxJob.getAdaptiveResolution():
                    #filter new reference according to current resolution and write to file
                    newReferenceVolume = read(exMaxJob.getDestination() + str(iteration) + '.em')
                    filtered = lowpassFilter(newReferenceVolume,ceil(resolutionBand))
                    filteredReference = filtered[0]
                    filteredReference.write(exMaxJob.getDestination() + str(iteration) + '_ANG' + str(resolutionAngstrom) + '.em')
                    filteredReferenceInv = filteredReference * -1
                    filteredReferenceInv.write(exMaxJob.getDestination() + str(iteration) + '_ANG' + str(resolutionAngstrom) + '_INV.em')
 
                #update particle list. determined angles will be set to each particle
                
                if exMaxJob.getClassificationParameter() > 0:
                    newParticleList = newAllParticleList
                
                exMaxJob.setParticleList(newParticleList)
                
                #save to XML and HTML
                exMaxJob.toXMLFile(exMaxJob.getDestination() + 'ExMaxJob-' + str(iteration) + '.xml')
                exMaxJob.toHTMLFile(exMaxJob.getDestination() + 'ExMaxJob-' + str(iteration) + '.html')
                
                print('Iteration ' + str(iteration) + ' finished')
                iteration = iteration +1
                
                print('')
                print('')
                
                if (not ((iteration <= exMaxJob.getNumberIterations()) and (not resultsEqual))) and sendFinishMessage:
                    print('Send finish message to other nodes')
                    manager.parallelEnd()
                    
        else:            
            worker = ExMaxWorker()

            worker.parallelRun()
            
        if sendFinishMessage:    
            pytom_mpi.finalise()
   
def sequentialStart(exMaxJob,verbose,sendFinishMessage = True):     
    """
    sequentialStart: same as parallelStart but processing will take place on one node only
    @param exMaxJob: Object defining alignment job
    @type exMaxJob: L{pytom.alignment.ExMaxAlignment.ExMaxJob}
    @param verbose: Will print debug information if True. (False by default)
    @param sendFinishMessage: Send finish message to worker nodes. True by default, False for multi ref alignment for instance 
    """
        
    from pytom_volume import read
    from pytom.basic.filter import lowpassFilter
    from pytom.alignment.alignmentFunctions import average,alignmentProgressToHTML
    from pytom.alignment.structures import AlignmentList
    from pytom.basic.resolution import bandToAngstrom,angstromToBand,angleFromResolution
    from math import ceil
    
    iteration = 1
    resultsEqual = False

    r = exMaxJob.getReference()
    rVol = r.getVolume()
    
    cubeSize = rVol.sizeX() 
    
    symmetry = exMaxJob.getSymmetry()
    
    while iteration <= exMaxJob.getNumberIterations() and not resultsEqual:
        #save current job
        exMaxJob.toXMLFile(exMaxJob.getDestination() + 'ExMaxJob-' + str(iteration) + '.xml')
        exMaxJob.toHTMLFile(exMaxJob.getDestination() + 'ExMaxJob-' + str(iteration) + '.html')
        
        manager = ExMaxManager(exMaxJob)
        
        manager.sequentialAlignment(verbose)
        
        alignmentList = manager.getAlignmentList()
        
        #dump alignmentList and exMaxJob to file
        print('Save results iteration ' + str(iteration))
        newAllParticleList = alignmentList.toParticleList()
        newAllParticleList.sortByScore('descending')
        newAllParticleList.toXMLFile(exMaxJob.getDestination() + 'ParticleList-' + str(iteration) + '.xml')                
        
        
        #save current average
        manager.setAlignmentList(alignmentList)
        
        
        if exMaxJob.getExcludeThreshold() > 0:
            newParticleList = newAllParticleList.listFromBestScorePercentage(exMaxJob.getExcludeThreshold())
            #reduce particleList to N% for averaging, determine resolution. 
            #the next round however will all particles again!
        else:
            newParticleList = newAllParticleList
        
        print('')
        print('Write average iteration ', iteration)
        
        if symmetry.getNFold() == 1:            
            newParticleList,exMaxJob.getDestination() + str(iteration) + '.em',False,verbose
            newReference = average(newParticleList,exMaxJob.getDestination() + str(iteration) + '.em',False,verbose)
        else:
            if verbose:
                print('Applying ', symmetry.getNFold(), 'fold symmetry to particles')
            symmetrizedParticleList = symmetry.apply(newParticleList)
            newReference = average(symmetrizedParticleList,exMaxJob.getDestination() + str(iteration) + '.em',False,verbose)
            
        #update reference
        exMaxJob.setReference(newReference)
        print('')
        print('Determine new parameters iteration ' + str(iteration))
        #set new refinement angles and update rotation strategy 
        
        sampleInfo = exMaxJob.getSampleInformation()
        if exMaxJob.getAdaptiveResolution():
            
            if verbose:
                print('Determine resolution.')
            if symmetry.getNFold() == 1:
                [resNyquist,resolutionBand,numberBands]  = newParticleList.determineResolution(
                    criterion=exMaxJob.getFSCCriterion(), numberBands=cubeSize/2, mask=exMaxJob.getMask(),
                    keepHalfsetAverages = True, halfsetPrefix=exMaxJob.getDestination() + 'fsc-'+str(iteration)+'-',
                    verbose=verbose )
            else:
                if verbose:
                    print('Applying ', symmetry.getNFold(), 'fold symmetry to resolution')
                [resNyquist,resolutionBand,numberBands]  = symmetrizedParticleList.determineResolution(
                    criterion=exMaxJob.getFSCCriterion(), numberBands =cubeSize/2, mask=exMaxJob.getMask(),
                    keepHalfsetAverages=True, halfsetPrefix=exMaxJob.getDestination() + 'fsc-'+str(iteration)+'-',
                    verbose=verbose )
                
            resolutionAngstrom = bandToAngstrom( resolutionBand,sampleInfo.getPixelSize(),numberBands,1 )
            #resolutionAngstrom = bandToAngstrom( resolutionBand,sampleInfo.getPixelSize(),numberBands,exMaxJob.getBinning() )

            if verbose:
                print('Current resolution :' + str(resolutionAngstrom) + ' Angstrom')
                
            #set new angular values
            try:
                refinementAngle = angleFromResolution(resolutionAngstrom,sampleInfo.getParticleDiameter())
                refinementAngle = refinementAngle * exMaxJob.getAngleFactor()
                oldRotations = exMaxJob.getRotations() 
                newRotations = oldRotations.focusRotation(refinementAngle=refinementAngle)
                exMaxJob.setRotations(newRotations)
                #update preprocessing bandpass (adaptive bandpass filter if enabled). note the  exMaxJob.getAdaptiveOffset()!!!
                preprocessing = exMaxJob.getPreprocessing()
                newResolution = resolutionAngstrom * (1 - exMaxJob.getAdaptiveOffset()) #determine better resolution (in ANGSTROM)        
                newResolutionBand = int( ceil( angstromToBand(resolution=newResolution,
                                                              pixelSize=sampleInfo.getPixelSize(),
                                                              numberOfBands=numberBands, scale=1))) # FF: fixed bug!
                                #numberOfBands=numberBands, scale=exMaxJob.getBinning()))) #convert upper value to BANDS
                #update job
                preprocessing.setHighestFrequency(newResolutionBand)
                exMaxJob.setPreprocessing(preprocessing)
            except RuntimeError:
                print('Warning for adaptive angular sampling and adaptive bandpass')
                print('Angle and filter will be unchanged due to error. Occurs when determined resolution is invalid due to too few particles.')
                print('Current resolution :' + str(resolutionAngstrom) + ' Angstrom')
                rotations = exMaxJob.getRotations()
                refinementAngle = rotations.getIncrement() 
                
        else:
            rotations = exMaxJob.getRotations()
            refinementAngle = rotations.getIncrement()
            oldRotations = rotations 
            newRotations = oldRotations.focusRotation(refinementAngle=refinementAngle)
            exMaxJob.setRotations(newRotations)
            
            resolutionAngstrom = bandToAngstrom(cubeSize/2,sampleInfo.getPixelSize(),cubeSize/2,1)
            #resolutionAngstrom = bandToAngstrom(cubeSize/2,sampleInfo.getPixelSize(),cubeSize/2,exMaxJob.getBinning())

        if exMaxJob.getAdaptiveResolution():
            #filter new reference according to current resolution and write to file
            newReferenceVolume = read(exMaxJob.getDestination() + str(iteration) + '.em')
            filtered = lowpassFilter(newReferenceVolume,ceil(resolutionBand))
            filteredReference = filtered[0]
            filteredReference.write(exMaxJob.getDestination()+str(iteration)+'_ANG'+str(resolutionAngstrom)+'.em')
            filteredReferenceInv = filteredReference * (-1)
            filteredReferenceInv.write(exMaxJob.getDestination()+str(iteration)+'_ANG'+str(resolutionAngstrom)+'_INV.em')


        #update particle list. determined angles will be set to each particle
        
        if exMaxJob.getClassificationParameter() > 0:
            newParticleList = newAllParticleList
        
        exMaxJob.setParticleList(newParticleList)
        
        #save to XML and HTML
        exMaxJob.toXMLFile(exMaxJob.getDestination() + 'ExMaxJob-' + str(iteration) + '.xml')
        exMaxJob.toHTMLFile(exMaxJob.getDestination() + 'ExMaxJob-' + str(iteration) + '.html')
        
        print('Iteration ' + str(iteration) + ' finished')
        iteration = iteration +1
        
        print('')
        print('')

        
        
def resume(job,verbose=False):
    """
    resume: Will resume a ExMaxJob
    @param job: Either a string pointing to Job XML file or a L{pytom.alignment.ExMaxAlignment.ExMaxJob} object
    @type job: str or L{pytom.alignment.ExMaxAlignment.ExMaxJob} 
    @param verbose: Will print debug information if True. (False by default)
    """
    from pytom.alignment.ExMaxAlignment import ExMaxJob
    
    if job.__class__ == str:
        em = ExMaxJob(0,0,0,0,0,0,0,0,0,0,0,0)
        em.fromXMLFile(job)
        job = em
        
    if job.__class__ == ExMaxJob:
        parallelStart(job,verbose)
        
    else:
        print('Could not resume job. Aborting!')
        
