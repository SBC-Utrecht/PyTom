# -*- coding: utf-8 -*-
"""
G{importgraph}
@author: Thomas Hrabe
"""

from pytom.basic.structures import PyTomClass
from pytom.alignment.ExMaxAlignment import ExMaxJob
def fromUnicode(code):

    return code.__str__()

class ExpectationMaximisationJob(ExMaxJob):
    """
    ExpectationMaximisationJob: Class specification refactored to ExMaxAlingment. This class is required for backward compatibility...
    """


class MaximisationJob(PyTomClass):
    """
    MaximisationJob : Stores all infos needed for a maximisation job
    """    
    def __init__(self,particle='',reference='',score='',rotations='',
            mask='',numberRefinementRounds=1,preprocessing='',binning=1):
        """
        @param particle: particle to be aligned
        @type particle: L{pytom.basic.structures.Particle} or str
        @param reference: reference for particle alignment
        @param score: type of score used for alignment
        @param rotations: rotations scanned in search
        @param mask: mask on reference
        @param numberRefinementRounds: iteration number
        @param preprocessing: preprocessing parameters
        @param binning: Binning Factor
        @type binning: int
        """
       
        from lxml.etree import _Element
        
        if particle.__class__ == _Element:
            self.fromXML(particle)
        else:
            if particle.__class__ == str:
                from pytom.basic.structures import Particle
                self.particle = Particle(particle)
            else:
                self.particle = particle
                
            self.reference = reference
            self.score = score
            self.rotations = rotations
            self.mask = mask   
            self.numberRefinementRounds = numberRefinementRounds
            self.preprocessing = preprocessing
            self.binning = binning
            
    def fromXML(self,xmlObj):
        """
        fromXML : Assigns values to job attributes from XML object
        @param xmlObj: A xml object  
        @author: Thomas Hrabe 
        """
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element:
            from pytom.basic.exceptions import ParameterError
            raise ParameterError('You must provide a valid XML-MaximisationJob object.')
        
        if xmlObj.tag == "JobDescription":
            jobDescription = xmlObj
        else:  
            jobDescription = xmlObj.xpath('JobDescription')
            
            if len(jobDescription) == 0:
                from pytom.basic.structures import PyTomClassError
                raise PyTomClassError("This XML is not an JobDescription.")
            
            jobDescription = jobDescription[0]

        from pytom.angles.angle import AngleObject
        from pytom.basic.score import fromXML as fromXMLScore
        from pytom.alignment.preprocessing import Preprocessing
        from pytom.basic.structures import Mask,Particle,Reference,ReferenceList
        
        self.binning = int(jobDescription.get('Binning'))
        
        particle_element = jobDescription.xpath('Particle')[0]
        p = Particle('')
        p.fromXML(particle_element)
        self.particle = p
        
        r = jobDescription.xpath('Reference')
        
        if len(r) > 0:
            ref = Reference('')
            ref.fromXML(r[0])
            self.reference = ref
        else:
            r = jobDescription.xpath('ReferenceList')
            ref = ReferenceList()
            ref.fromXML(r[0])
            self.reference = ref
                    
        mask = jobDescription.xpath('Mask')[0]
        self.mask = Mask('')
        self.mask.fromXML(mask)
        
        self.numberRefinementRounds = jobDescription.get('NumberRefinementRounds')
        self.numberRefinementRounds = int(self.numberRefinementRounds)
        
        score = jobDescription.xpath('Score')
        self.score = fromXMLScore(score[0])
        
        angles = jobDescription.xpath('Angles')
        ang = AngleObject()
        self.rotations = ang.fromXML(angles[0])
        
        preObj = xmlObj.xpath('Preprocessing')
        if len(preObj) == 0:
            self.preprocessing = Preprocessing()
        else:
            p = Preprocessing()
            p.fromXML(preObj[0])
            self.preprocessing = p
        
        
    def toXML(self):        
        """
        toXML : Compiles a XML file from job object
        @author: Thomas Hrabe
        """    
        
        from lxml import etree
             
        jobElement = etree.Element("JobDescription")
        
        jobElement.set("NumberRefinementRounds",str(self.numberRefinementRounds))
        
        jobElement.set("Binning",str(self.binning))
        
        jobElement.append(self.particle.toXML())

        jobElement.append(self.reference.toXML())

        jobElement.append(self.mask.toXML())
        jobElement.append(self.rotations.toXML())
        
        if self.score.__class__ == str:
            jobElement.append(self.score)
        else:
            jobElement.append(self.score.toXML())
        
        preObj = self.preprocessing.toXML()
        jobElement.append(preObj)
        
        return jobElement


    def check(self): 
        """
        check: Performs check on self whether all settings were sane. Paths and Files exists 
        """
        
        from pytom.tools.files import checkFileExists

        returnValue = checkFileExists(self.particle.getFilename())
        if not returnValue:
            raise IOError(str(self.particle.getFilename()) + ' not found!')
            
        returnValue = returnValue and checkFileExists(self.reference.getFilename())
        
        if not returnValue:
            raise IOError(str(self.reference) + ' not found!')
        
        return returnValue

  
class MaximisationResult(PyTomClass):
    """
    MaximisationResult : Stores results of one maximisation process
    """

    def __init__(self,particle='',reference=-1.0,score=-1.0,shift=-1.0,rotation=-1.0,angleObject=-1):
        
        from pytom.basic.structures import Particle,Reference,Shift,Rotation
        from numpy import long

        if particle.__class__ == str:
            self._particle = Particle(particle)
        elif particle.__class__ == Particle:
            self._particle = particle
        else:
            self._particle = Particle()

        if reference.__class__ == str:
            self._reference = Reference(reference)
        elif reference.__class__ == Reference:
            self._reference = reference
        else:
            self._reference = Reference()            
        
        if shift.__class__ == list:
            self._shift = Shift(shift)
        elif shift.__class__ == float:
            self._shift = Shift()
        else:
            self._shift = shift
        
        if rotation.__class__ == list:
            self._rotation = Rotation(rotation)
        elif rotation.__class__ == Rotation:
            self._rotation = rotation
        else:
            self._rotation = Rotation()
            
        if score.__class__ == float:
            from pytom.basic.score import xcfScore
            self._score = xcfScore()
        else:
            self._score = score
            
        if angleObject.__class__ == float or isinstance(angleObject, (int, long)):
            from pytom.angles.angleList import AngleList
            self._angleObject = AngleList()
        else:
            self._angleObject = angleObject
            
    def toParticle(self):
        """
        toParticle: Converts this object to a Particle object.
        @return:
        @rtype: L{pytom.basic.structures.Particle}
        """
        particle = self._particle
        particle.setRotation(self._rotation)
        particle.setShift(self._shift)
        particle.setScore(self._score)
        
        return particle
    
    def getParticle(self):
        return self._particle
    
    def getShift(self):
        from pytom.basic.structures import Shift
        
        if self._shift.__class__ == list:
            return Shift(self._shift)
        else:
            return self._shift
    
    def setShift(self,shift):
        from pytom.basic.structures import Shift
        
        assert shift.__class__ == Shift
        
        self._shift = shift
        
    def getAngleObject(self):    
        return self._angleObject
    
    def getRotation(self):
        
        from pytom.basic.structures import Rotation
        
        if self._rotation.__class__ == list:
            return Rotation(self._rotation[0],self._rotation[1],self._rotation[2])
        else:
            return self._rotation.copy()
    
    def getScore(self):
        """
        getScore: Returns score object
        """
        return self._score
    
    def setRotation(self,rotation):
        """
        setRotation:
        @param rotation: 
        """
        from pytom.basic.structures import Rotation
        
        if rotation.__class__ == list:
            rotation = Rotation(rotation)
            
        self._rotation = rotation
        
    def fromXML(self,xmlObj):
        """
        fromXML : Assigns values to result attributes from XML object
        @param xmlObj: A xml object  
        @author: Thomas Hrabe 
        """
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element:
            from pytom.basic.exceptions import ParameterError
            raise ParameterError('Is not a lxml.etree._Element! You must provide a valid XML object.') 
        
        from pytom.basic.score import fromXML as fromXMLScore 
        from pytom.angles.angle import AngleObject
        
        if xmlObj.tag == "Result":
            result = xmlObj
        else:  
            result = xmlObj.xpath('Result')
            
            if len(result) == 0:
                raise PyTomClassError("This XML is not an MaximisationResult. No Result provided.")
            
            result = result[0]
        from pytom.basic.structures import Particle,Reference,Rotation,Shift
        particle_element = result.xpath('Particle')[0]
        p = Particle('')
        p.fromXML(particle_element)
        self._particle = p

        r = result.xpath('Reference')
        ref = Reference('')
        ref.fromXML(r[0])
        self._reference = ref
        
        scoreXML = result.xpath('Score')[0]
        self._score = fromXMLScore(scoreXML)
        
        shiftXML = result.xpath('Shift')[0]    
        self._shift = Shift()
        self._shift.fromXML(shiftXML)
        
        rotationXML = result.xpath('Rotation')[0]
        self._rotation = Rotation()    
        self._rotation.fromXML(rotationXML)
        
        angleElement = result.xpath('Angles')
        ang = AngleObject()
        self._angleObject = ang.fromXML(angleElement[0])
    
        
    def toXML(self):
        """
        toXML : Compiles a XML from result object
        @author: Thomas Hrabe
        """    
        from lxml import etree
        
        resultElement = etree.Element("Result")
        
        resultElement.append(self._particle.toXML())
        
        if self._reference.hasGeneratedByInfo():
            from pytom.basic.structures import Reference
            newRef = Reference(self._reference.getReferenceFilename())
            resultElement.append(newRef.toXML())
        else:
            resultElement.append(self._reference.toXML())
        
        resultElement.append(self._shift.toXML())
        
        resultElement.append(self._rotation.toXML())
                
        resultElement.append(self._score.toXML())
        
        resultElement.append(self._angleObject.toXML())
        
        return resultElement
    
    def copy(self):
        return MaximisationResult(self._particle,self._reference,self._score,self._shift,self._rotation,self._angleObject)
                                  
    
class ExpectationJob(PyTomClass):
    
    def __init__(self, particleList = None , newAverageName = ''):
        """
        @param particleList: particle list
        @type particleList: L{pytom.basic.structures.ParticleList}
        @param newAverageName: name of output average
        @type newAverageName: L{str}
        """
        from pytom.basic.structures import ParticleList
        
        self._particleList = particleList or ParticleList('/',[])
        self._newAverageName = newAverageName
        
    def appendMaximisationResult(self,result):
        """
        @param result: maximization result
        @type result: L{pytom.alignment.structures.MaximisationResult}
        """
        from pytom.alignment.structures import MaximisationResult
        
        if not result.__class__ == MaximisationResult:
            raise RuntimeError('The object you are appending must be a MaximisationResult!')
        
        self._particleList.append(result.toParticle())
    
    def appendParticle(self,particle):
        """
        @param particle: particle
        @type particle: L{pytom.basic.structures.Particle}
        """
        self._particleList.append(particle)

    def setParticleList(self, particleList):
        """
        @param particleList: particle list
        @type particleList: L{pytom.basic.structures.ParticleList}
        """
        self._particleList = particleList
        
    def getParticleList(self):
        return self._particleList
    
    def getNewAverageName(self):
        return self._newAverageName
        
    def fromXML(self,xmlObj):
        """
        fromXML : Assigns values to job attributes from XML object
        @param xmlObj: A xml object
        @type xmlObj: L{lxml.etree._Element}

        @author: Thomas Hrabe 
        """
        from lxml.etree import _Element
        from pytom.basic.structures import ParticleList

        if xmlObj.__class__ != _Element:
            from pytom.basic.exceptions import ParameterError
            raise ParameterError('Is not a lxml.etree._Element! You must provide a valid XML object.')

        if not xmlObj.tag == 'ExpectationJob':
            jobElement = xmlObj.xpath('ExpectationJob')
        
            if len(jobElement) == 0:
                from pytom.basic.exceptions import ParameterError
                raise ParameterError('You must provide a valid XML-ExpectationJob object.')
            else:
                jobElement = jobElement[0]
        else:
            jobElement = xmlObj
             
        self._newAverageName = jobElement.get('AverageName').__str__()
        
        particleListXML = xmlObj.xpath('ParticleList')
        self._particleList = ParticleList('/')
        if len(particleListXML) > 0:
            self._particleList.fromXML(particleListXML[0])
        
    def toXML(self):
        """
        toXML : Compiles a XML file from job object
        @author: Thomas Hrabe
        """ 
        from lxml import etree
        
        expJob = etree.Element("ExpectationJob", AverageName = self._newAverageName)
        expJob.append(self._particleList.toXML())
        return expJob
        
            
class ExpectationResult(PyTomClass):
    def __init__(self,resultName):
        self.resultFileName = resultName
        self.wedgeSumName = resultName[:len(resultName)-3] + '-WedgeSum.em'
        
    def getResult(self):
        return self.resultFileName
    
    def getWedgeSumName(self):
        return self.wedgeSumName
    
    def fromXML(self,xmlObj):
        """
        fromXML : Assigns values to result attributes from XML object
        @param xmlObj: A xml object  
        @author: Thomas Hrabe 
        """
        from lxml.etree import _Element

        if xmlObj.__class__ != _Element :
            from pytom.basic.exceptions import ParameterError
            raise ParameterError('Is not a lxml.etree._Element! You must provide a valid XML object.')
        
        self.resultFileName = xmlObj.get('ResultFileName')
        self.wedgeSumName = xmlObj.get('WedgeSumName')
        
    def toXML(self):
        """
        toXML : Compiles a XML file from result object
        @author: Thomas Hrabe
        """ 
        from lxml import etree
        
        result = etree.Element("ExpectationResult", ResultFileName = self.resultFileName, WedgeSumName = self.wedgeSumName)
          
        return result

    
class AlignmentList(PyTomClass):
    """
    AlignmentList : Stores alignment of particles. Similar to the av3 - motif list, rotations and shifts are stored for each individual. 
    The list of L{pytom.alignment.structures.MaximisationResult} stores the information.
    """
    def __init__(self):
        """
        __init__: Initialise
        """
        self._alignmentList = []
        
    def fromXML(self,xmlObj):
        """
        @param xmlObj: A xml object  
        @type xmlObj: L{lxml.etree._Element}
        @author: Thomas Hrabe 
        """
        
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element :
            from pytom.basic.exceptions import ParameterError
            raise ParameterError('Is not a lxml.etree._Element! You must provide a valid XML-WedgeInfo object.')
   
        from pytom.alignment.structures import MaximisationResult
        
        if xmlObj.tag == 'AlignmentList':
            list_element = xmlObj
        else:
            list_element = xmlObj.xpath('AlignmentList')
            list_element = list_element[0]
        
        self._alignmentList = []
        
        for child in list_element.iter('Result'):
            result = MaximisationResult()
            result.fromXML(child)
            self._alignmentList.append(result)
            
    def toXML(self):
        """
        toXML : Compiles a XML file from result object
        rtype : L{lxml.etree._Element}
        @author: Thomas Hrabe
        """ 
        from lxml import etree
        
        list_element = etree.Element('AlignmentList')
        
        for result in self._alignmentList:
            list_element.append(result.toXML())
        
        return list_element
    
    def append(self,object):
        """
        append: Appends object to the alignment list. 
        @param object: 
        @type object: L{pytom.alignment.structures.MaximisationResult}
        @author: Thomas Hrabe  
        """
        from pytom.alignment.structures import MaximisationResult
        assert object.__class__ == MaximisationResult
        
        self._alignmentList.append(object.copy())
    
    def getList(self):
        """
        getList: getter for list
        @return: List of alignments
        @rtype: L{list}
        """
        return self._alignmentList
    
    def setList(self,list):
        from pytom.alignment.structures import MaximisationResult
        assert list[0].__class__ == MaximisationResult
        
        self._alignmentList = list
        
    def getResultByParticle(self,particle,alXML = None):
        """
        getResultByParticle: Returns alignment result of specified particle
        @param particle: The particle
        @type particle: Either str of L{pytom.basic.structures.Particle}  
        @rtype: L{pytom.alignment.structures.MaximizationResult}
        """
        from pytom.basic.structures import Particle
        
        if alXML == None:
            alXML = self.toXML()
        
        if particle.__class__ == Particle:
            particleName = particle.getFilename()
        else:
            particleName = particle
            
        #r = alXML.xpath('/AlignmentList/Result/Particle[@Filename="' + particleName + '"]/..') 
        
        #if len(r) == 0:
        #    raise Exception('Particle ' + particle.getFilename() + ' was not found in this alignment list!')
        
        result = None
        found = False
        for res in self._alignmentList:
            if found:
                break
            
            p = res.getParticle()
            
            if p.getFilename() == particleName: 
                result = MaximisationResult()
                result.fromStr(str(res))
                found = True
               
        if not found:
            raise RuntimeError('Particle ' + particle.getFilename() + ' was not found in this alignment list!')
        
        return result

    def sortByParticleList(self,particleList):
        """
        sortByParticleList: Sorts particles according to particleList order. Results might be in different order after parallel processing
        @param particleList:
        @type particleList: L{pytom.basic.structures.ParticleList} 
        """
        
        newList = []
        selfXML = self.toXML()
        
        for particle in particleList:
            result = self.getResultByParticle(particle, selfXML)
            newList.append(result)
            
        self._alignmentList = newList
        
    def toMotif(self,filename):
        """
        toMotif:
        @param filename: 
        @type filename: L{str}
        @todo: add unit test
        """
        from pytom_volume import vol
        l = self.len()
        motif = vol(20,l,1)
        motif.setAll(0)
        for i in range(l):
            
            #correlation coefficient
            motif.setV(1,0,i,0)
            
            motif.setV(i+1,3,i,0)
            
            res = self._alignmentList[i]
            shift = res.getShift()
            motif.setV(shift[0],10,i,0)
            motif.setV(shift[1],11,i,0)
            motif.setV(shift[2],12,i,0)
            
            rotation = res.getRotation()
            motif.setV(rotation[0],16,i,0)
            motif.setV(rotation[1],17,i,0)
            motif.setV(rotation[2],18,i,0)
            
        motif.write(filename)


    def classifyResultScores(self,peakPercentage = 10):
        """
        classifyList: Classifies results according to their score value. Any result will remain unless it's score value is greater than the mean score value / percentage
        @param peakPercentage: value of score. The formula is mean * peakPercentage / 100. Default value will be 10. Negative values are also valid. 
        @return: alignmentList of 'good' results. (if no result matches the criterion, all results will be returned)  
        """
        from pytom.alignment.structures import MaximisationResult
        
        #get all score values from results through xpath
        values = self.xpath('/AlignmentList/Result/Score/@Value')
        
        mean = 0
        
        for value in values:
            mean = mean + float(value)
            
        #determine mean score value
        mean = mean / len(values)
    
        classificationValue = mean * peakPercentage / 100
        
        #select all results above threshold through xpath
        goodResults = self.xpath('/AlignmentList/Result[Score/@Value >= ' + str(classificationValue) + ']')
        
        newResults = []
        
        for resultXML in goodResults:
            result = MaximisationResult()
            result.fromXML(resultXML)
            
            newResults.append(result)
        
        newAl = AlignmentList()
        
        if len(newResults) > 0:
            newAl.setList(newResults)
        else:
            newAl = self
            
        return newAl
        

    def sumOfScores(self):
        """
        sumOfScores: Returns average score value for current alignment
        return:
        """
        sum = 0
        
        for i in xrange(0,len(self)):
            
            mxRes = self[i] 
            sc = mxRes.getScore()
            sum = sum + float(sc.getValue())
    
        return sum / len(self)
    
    def distanceOfRotations(self,oldAlignmentList):
        """
        distanceOfRotations:Determines distance of current rotations from previous rotations 
        @param oldAlignmentList: 
        """
        from pytom.tools.maths import listMean,listStd
        
        distance = []
        numberResults = len(self) 
        
        oldListXML = oldAlignmentList.toXML()
        
        for i in xrange(0,numberResults):
            newResult = self[i] 
            particle = newResult.getParticle()
            
            try:
                #oldResult = oldAlignmentList.getResultByParticle(particle)

                oldResultXML = oldListXML.xpath('/AlignmentList/Result/Particle[@Filename="' + particle.getFilename() + '"]/..')
                oldResult = MaximisationResult()
                oldResult.fromXML(oldResultXML[0])
                
            except:
                numberResults = numberResults - 1
                continue
            
            newRotation = newResult.getRotation()
            oldRotation = oldResult.getRotation()
            
            newQuaternion = newRotation.toQuaternion()
            oldQuaternion = oldRotation.toQuaternion()
            
            d = newQuaternion.distance(oldQuaternion)
            
            distance.append(d)
            
        return [listMean(distance),listStd(distance)]
    
    def distanceOfShifts(self,oldAlignmentList):
        """
        distanceOfShifts: Determines distance of current shifts from previous shifts 
        @param oldAlignmentList: 
        """
        from pytom.tools.maths import euclidianDistance,listMean,listStd
        
        distance = []
        numberResults = len(self) 
        
        oldListXML = oldAlignmentList.toXML()
        
        for i in xrange(0,numberResults):
            
            newResult = self[i] 
            particle = newResult.getParticle()
            
            try:
                oldResultXML = oldListXML.xpath('/AlignmentList/Result/Particle[@Filename="' + particle.getFilename() + '"]/..')
                oldResult = MaximisationResult()
                oldResult.fromXML(oldResultXML[0])
            except:
                numberResults = numberResults - 1
                continue
            
            newShift = newResult.getShift()
            oldShift = oldResult.getShift()
            
            d = euclidianDistance(newShift.toVector(),oldShift.toVector())
            
            distance.append(d)
        
        return [listMean(distance),listStd(distance)]
       
    def toParticleList(self):
        """
        toParticleList: Converts this object to a Particle List
        @return: A particle list that can be used for further processing
        @rtype: L{pytom.basic.structures.ParticleList} 
        """
        from pytom.basic.structures import ParticleList
        pl = ParticleList('/')
        
        for result in self._alignmentList:
            
            p = result.toParticle()
            pl.append(p)
            
        return pl
            
    def __len__(self):
        return len(self._alignmentList)
    
    def __getitem__(self,key):
        return self._alignmentList[key]
    
    def __setitem__(self,key,value):
        if isinstance(key, (int, long)):
            if key < len(self._alignmentList):
                    self._alignmentList[key] = value
            else:
                IndexError('Index is too large!')
                
    def clear(self):
        self._alignmentList = []


class GrowingAverageJob(PyTomClass):
    """
    GrowingAverageJob:
    @ivar particleList: List of particles to be aligned
    @ivar angleObject: Angle object L{pytom.angles.AngleObject}
    @ivar startParticleNumber: Number of start particle (default 0)
    @ivar maskFile: Mask used for appedizing 
    """
    
    def __init__(self,particleList=None,angleObject=None,maskFile=None,scoreObject=None,startParticleNumber=0,destinationDirectory='.',preprocessing = None):
        
        from pytom.tools.files import checkDirExists
        from pytom.angles.angleList import AngleList
        from pytom.basic.structures import ParticleList
        
        self.particleList = particleList or ParticleList('/')
        self.angleObject = angleObject or AngleList()
        self.startParticleNumber = startParticleNumber
        self.maskFile = maskFile or None
        
        if preprocessing:
            self.preprocessing = preprocessing
        else:
            from pytom.alignment.preprocessing import Preprocessing
            self.preprocessing = Preprocessing()
        
        
        if self.maskFile.__class__ == str:
            from pytom.basic.structures import Mask 
            self.maskFile = Mask(self.maskFile)
            
        self.score = scoreObject
        
        if not checkDirExists(destinationDirectory):
            raise Exception('Destination directory ' + destinationDirectory + ' does not exist.')
        
        if not destinationDirectory[len(destinationDirectory)-1] == '/':
            destinationDirectory = destinationDirectory + '/'
        
        self.destinationDirectory = destinationDirectory
        
    def toXML(self):
        """
        toXML : Compiles a XML file from result object
        rtype : L{lxml.etree._Element}
        @author: Thomas Hrabe
        """ 
        from lxml import etree
        
        job_element = etree.Element('GrowingAverageJob',StartParticleNumber = str(self.startParticleNumber), DestinationDirectory = str(self.destinationDirectory))
        
        job_element.append(self.maskFile.toXML())
        
        job_element.append(self.particleList.toXML())
        
        job_element.append(self.angleObject.toXML())
        
        job_element.append(self.score.toXML())
        
        job_element.append(self.preprocessing.toXML())
        
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
        
        self.startParticleNumber = int(job_element.get('StartParticleNumber'))
        
        mask = job_element.xpath('Mask')[0]
        from pytom.basic.structures import Mask
        self.maskFile = Mask('')
        self.maskFile.fromXML(mask)
         
        self.destinationDirectory = job_element.get('DestinationDirectory')
        
        particleXML = job_element.xpath('ParticleList')[0]
        
        self.particleList = ParticleList('/',[])
        self.particleList.fromXML(particleXML)
        
        angleXML = job_element.xpath('Angles')[0]
        ang = AngleObject()
        self.angleObject = ang.fromXML(angleXML)
        
        scoreXML = job_element.xpath('Score')[0]
        
        self.score = scoreFromXML(scoreXML)
        
        self.preprocessing = Preprocessing()
        preprocessingXML = job_element.xpath('Preprocessing')[0]
        self.preprocessing.fromXML(preprocessingXML)


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
        
        from pytom.basic.structures import Particle
        
        particleXML = result_element.xpath('/GrowingAverageInterimResult/Particle')[0]
        
        self.particle = Particle('')
        self.particle.fromXML(particleXML)
        
        from pytom.basic.structures import Reference
        referenceXML = result_element.xpath('/GrowingAverageInterimResult/Result')[0]
        self.reference = Reference('')
        self.reference.fromXML(referenceXML)


class Peak(PyTomClass):
    def __init__(self,scoreValue,rotation,shift):
        """
        Peak of a correlation search
        @param scoreValue: value of score
        @type scoreValue: L{float}
        @param rotation: Rotation (Orientation) corresponding to max score value
        @type rotation: L{pytom.basic.structures.Rotation}
        @param shift: Translation corresponding to max score value
        @type shift: L{pytom.basic.structures.Shift}
        """
        
        from pytom.basic.structures import Rotation,Shift
        assert rotation.__class__ == Rotation
        assert shift.__class__ == Shift
        
        self._scoreValue = scoreValue
        self._rotation = rotation
        self._shift = shift
        
    def getShiftLength(self):
        return self._shift.getLength()
    
    def getScoreValue(self):
        return self._scoreValue
    
    def getShift(self):
        return self._shift
    
    def getRotation(self):
        return self._rotation
    
    def __lt__(self,otherPeak):
        return self._scoreValue < otherPeak.getScoreValue()
    
    def __le__(self,otherPeak):
        return self._scoreValue <= otherPeak.getScoreValue()
    
    def __eq__(self,otherPeak):
        return self._scoreValue == otherPeak.getScoreValue()
    
    def __ge__(self,otherPeak):
        return self._scoreValue >= otherPeak.getScoreValue()
    
    def __gt__(self,otherPeak):
        return self._scoreValue > otherPeak.getScoreValue()
    
    def __ne__(self,otherPeak):
        return not self == otherPeak
       
    def toXML(self):
        """
        toXML : Compiles a XML file from result object
        rtype : L{lxml.etree._Element}
        @author: Thomas Hrabe
        """ 
        
        from lxml import etree
        
        peak_element = etree.Element('Peak',Value = str(self._scoreValue))
        
        peak_element.append(self._rotation.toXML())
        peak_element.append(self._shift.toXML())
        
        return peak_element
