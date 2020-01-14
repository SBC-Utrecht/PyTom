'''
Created on Feb 8, 2011

@author: hrabe
'''


from pytom.basic.structures import PyTomClass
from pytom.cluster.mcoEXMXStructures import MCOEXMXJob,SwapList


class RandomJump(PyTomClass):
    """
    RandomJump
    """
    
    def __init__(self,result,temperature,oldScore,newScore,randomNumber,criterion):
        self._result = result
        self._temparature = temperature
        self._oldScore = oldScore
        self._newScore = newScore
        self._randomNumber = randomNumber
        self._criterion = criterion
        
    def toXML(self):
        from lxml import etree
        
        jumpElement = etree.Element('Jump',OldScore=str(self._oldScore),NewScore=str(self._newScore),RandomNumber = str(self._randomNumber),Criterion = str(self._criterion))
    
        jumpElement.append(self._result.toXML())
        jumpElement.append(self._temparature.toXML())
    
        return jumpElement
    
def criterionFromXML(xmlObj):
    
    criterionName = str(xmlObj.get('Name'))
    
    returnCriterion = None
    
    if criterionName == 'MetropolisCriterion':
        returnCriterion = MetropolisCriterion(allowClassCollapse=xmlObj.get('AllowClassCollapse') == 'True') 
    elif criterionName == 'ThresholdAcceptance':
        returnCriterion = ThresholdAcceptance(allowClassCollapse=bool(xmlObj.get('AllowClassCollapse')))   
    
    return returnCriterion
    
class AnnealingCriterion(PyTomClass):
    """
    AnnealingCriterion: Abstract class for annealing criteria
    """
    def __init__(self,name,allowClassCollapse):
        self._name = name
        
        self._allowClassCollapse = allowClassCollapse
        
    def setAllowClassCollapse(self,allowClassCollapse):
        
        if not allowClassCollapse.__class__ == bool:
            raise RuntimeError('Parameter must be of type bool!')
        
        self._allowClassCollapse = allowClassCollapse
        
    def getAllowClassCollapse(self):
        return self._allowClassCollapse
    
    def getName(self):
        return self._name
    
    def toXML(self):
        from lxml import etree
        
        return etree.Element('AnnealingCriterion',Name=str(self._name),AllowClassCollapse=str(self._allowClassCollapse))
    
    def fromXML(self,xmlObj):
        criterion = criterionFromXML(xmlObj)
        
        self._name = criterion.getName()
        self._allowClassCollapse = criterion.getAllowClassCollapse()
    def apply(self,results,temperature,swapList = SwapList(),verbose=False):
        """
        apply: abstract function 
        @param results:
        @param temperature:
        @param swapList: 
        @return: [bestResult,bestCluster]   
        """
        raise RuntimeError('This is a abstract class method! pytom.cluster.annealingStructures.AnnealingCriterion')
        

class ThresholdAcceptance(AnnealingCriterion):
    """
    ThresholdAcceptance: Threshold Accepting: A General Purpose Optimization Algorithm Appearing Superior to Simulated Annealing - GUNTER DUECK AND TOBIAS SCHEUER
    """
    
    
    def __init__(self,allowClassCollapse=False):
        super(self.__class__,self).__init__('ThresholdAcceptance',allowClassCollapse)
    
    def apply(self,results,temperature,swapList = SwapList(),verbose=False):    
        """
        apply: Apply randomized  
        @param results:
        @param temperature:
        @param swapList: 
        @return: [bestResult,bestCluster]   
        """
        
        bestResult = results[0]         
        clusterIterator = 0
        bestCluster = 0
        
        #loop over all results
        for result in results:
            if verbose:
                print('')
                print(bestResult)
                print(result)
            
            difference = bestResult.getScore().getValue() - result.getScore().getValue()
            
            if difference <= 0:
                #accept new score if it is better than the old one...
                bestResult = result
                bestCluster = clusterIterator
                continue
            
            if verbose:
                print('Difference : ' + str(difference))
                print('Temperature : ' + str(temperature))
                
                
            if difference > temperature.getTemperature():
                if verbose:
                    print('Swap!')
                                  
                jump = RandomJump(result,temperature,bestResult.getScore().getValue(),result.getScore().getValue(),-1,-1)
                
                swapList.append(jump)
                
                bestResult = result
                bestCluster = clusterIterator
                
            clusterIterator = clusterIterator + 1
            
        return [bestResult,bestCluster,swapList]
        
        
class MetropolisCriterion(AnnealingCriterion):
    """
    MetropolisCriterion: As defined in http://mathworld.wolfram.com/SimulatedAnnealing.html
    """
    
    def __init__(self,allowClassCollapse=False):
        super(self.__class__,self).__init__('MetropolisCriterion',allowClassCollapse)
        
        
    def apply(self,results,temperature,swapList = SwapList(),verbose = False):
        """
        apply: Apply randomized  
        @param results:
        @param temperature:
        @param swapList: 
        @return: [bestResult,bestCluster]   
        """
        
        from math import e
        import random
    
        random.seed()
        
        #sort results in descending order
        
        bestResult = results[0]         
        bestCluster = 0
        
        #allows only one jump per particle!
        didJumpBefore = False
        
        #loop over all results
        for resultIndex in range(len(results)):
            result = results[resultIndex]
            
            if verbose:
                print('')
                print(bestResult)
                print(result)
            
            difference = bestResult.getScore().getValue() - result.getScore().getValue()
            
            if difference <= 0.0:
                #accept new score if it is better than the old one...
                bestResult = result
                bestCluster = resultIndex
                continue
            
            if temperature.getTemperature() <= 0:
                continue
            
            if verbose:
                print('Difference : ' + str(difference))
                print('Temperature : ' + str(temperature))
                
            criterion = e**(- difference / temperature.getTemperature())
            
            if verbose:
                print('Criterion : ' + str(criterion))
            
            rand = random.random()
            
            if verbose:
                print('Random : '+ str(rand))
                
            if criterion > rand and not didJumpBefore:
                if verbose:
                    print('Swap!')
                                  
                jump = RandomJump(result,temperature,bestResult.getScore().getValue(),result.getScore().getValue(),rand,criterion)
                
                swapList.append(jump)
                
                bestResult = result
                bestCluster = resultIndex
                
                didJumpBefore = True
            
        return [bestResult,bestCluster,swapList]
            
def temperatureFromXML(xmlObj):
    
    from lxml.etree import _Element
    
    if xmlObj.__class__ != _Element :
        raise Exception('Is not a lxml.etree._Element! You must provide a valid XMLobject.')
    
    temperatureName = str(xmlObj.get('Name'))
    
    temperature = float(xmlObj.get('Temperature'))
    temperatureStep = float(xmlObj.get('TemperatureStep'))
    
    temperatureObject = None
    
    if temperatureName == 'AnnealingTemperature':
        temperatureObject = AnnealingTemperature(temperature,temperatureStep)
    elif temperatureName == 'SigmaTemperature':
        temperatureObject = SigmaTemperature(temperature,temperatureStep)
        
    return temperatureObject
     
class AnnealingTemperature(PyTomClass):
    """
    AnnealingTemperature: Defines the annealing temperature including decrement step
    """
    
    def __init__(self,temperature,temperatureStep=1):
    
        
        self._temperature = float(temperature)
        self._temperatureStep = float(temperatureStep)
        
    def getTemperature(self): 
        if self._temperature == 0:
            raise RuntimeError('Annealing temperature is 0, do something!')
        
        return self._temperature
    
    def decreaseTemperature(self):
        self._temperature = self._temperature - self._temperatureStep
    
    def numberIterations(self):
        """
        numberIterations: Returns the number of iterations expected for current cooling rate
        """
        from math import ceil
        
        return ceil(float(self._temperature) / self._temperatureStep)
        
        
    def cooledDown(self):
        """
        cooledDown: Returns True if temperature <= 0
        """
        return self._temperature <= 0
    
    def initializeTemperature(self,alignmentLists = None):
        """
        initializeTemperature:
        """
        pass
    
    def toXML(self):
        from lxml import etree
        
        temperatureElement = etree.Element('AnnealingTemperature',Name = 'AnnealingTemperature',Temperature=str(self._temperature),TemperatureStep=str(self._temperatureStep))
    
        return temperatureElement
    
        
class SigmaTemperature(AnnealingTemperature): 
    """
    SigmaTemperature: Annealing temperature will be set according to sigma scalled by current temperature
    """
    def toXML(self):
        from lxml import etree
        
        temperatureElement = etree.Element('AnnealingTemperature',Name = 'SigmaTemperature',Temperature=str(self._temperature),TemperatureStep=str(self._temperatureStep))
    
        return temperatureElement
    
    def initializeTemperature(self,alignmentLists = None):
        """
        initializeTemperature:
        @param alignmentLists: list of alignmentLists 
        """
        
        if alignmentLists == None or len(alignmentLists) == 0:
            raise RuntimeError('AlignmentLists parameter was empty!')
    
        from pytom.alignment.structures import AlignmentList
        from pytom.tools.maths import listStd
        
        if not alignmentLists[0].__class__ == AlignmentList:
            raise RuntimeError('Please provide a list of AlignmentLists here!')
    
        scores = []
        
        for alignmentList in alignmentLists:
            for resultIterator in range(len(alignmentList)):
                
                result = alignmentList[resultIterator]
        
                scores.append(result.getScore().getValue())
        
        scoresSTD  = listStd(scores) 
        
        self._temperature = self._temperature * scoresSTD
        
        
class MCOACJob(MCOEXMXJob):
    """
    MCOACJob: Specifies an annealed clustering job, extends MCOEXMXJOB
    """
    
    def __init__(self,particleList=None,destinationDirectory=None,mask=None,score=None,preprocessing=None,wedgeInfo=None,binning=None,sampleInformation=None,numberClasses=None,temperature=None,criterion=None,endThreshold=None,localSearchIncrement=None,symmetry=None):
        """
        __init__: Call parent constructor 
        @param temperature: Annealing temperature
        @type temperature: L{pytom.cluster.annealingStructures.AnnealingTemperature}
        @param criterion: 
        @type criterion: Child of L{pytom.cluster.annealingStructures.AnnealingCriterion} 
        """
        if not temperature:
            numberIterations = 0
        else:
            numberIterations = temperature.numberIterations()
        super(self.__class__,self).__init__(particleList,numberIterations,destinationDirectory,mask,score,preprocessing,wedgeInfo,binning,sampleInformation,numberClasses,endThreshold,symmetry)
        
        self._temperature = temperature
        self._criterion = criterion
        
        if not localSearchIncrement and temperature:
            localSearchIncrement = float(temperature.numberIterations())
        else:
            localSearchIncrement = numberIterations
            
        self._localSearchIncrement = float(localSearchIncrement)
        
    def toXML(self):
        from lxml import etree
        
        xml = super(self.__class__,self).toXML()
        xml.append(self._temperature.toXML())
        xml.append(self._criterion.toXML())
        
        xml.tag = 'MCOACJob'
        xml.set('LocalSearchIncrement',str(float(self._localSearchIncrement)))
        
        return xml
    
    def fromXML(self,xmlObj):
        
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element :
            raise Exception('Is not a lxml.etree._Element! You must provide a valid XMLobject.')
        
        temperatureXML = xmlObj.xpath('AnnealingTemperature')
        self._temperature = temperatureFromXML(temperatureXML[0])
        
        self._criterion = criterionFromXML(xmlObj.xpath('AnnealingCriterion')[0])
        self._localSearchIncrement = float(xmlObj.get('LocalSearchIncrement'))
        
        xmlObj.tag = 'MCOEXMXJob'
        super(self.__class__, self).fromXML(xmlObj)
            
    def getTemperature(self):
        return self._temperature
    
    def getCriterion(self):
        return self._criterion
    
    def getAllowClassCollapse(self):
        return self._allowClassCollapse
    
    def getLocalIncrement(self):
        return self._localSearchIncrement
    
    def setAllowClassCollapse(self,value):
        self._allowClassCollapse = value
    
    def decreaseTemperature(self):
        """
        decreaseTemperature: Decreases temperature of annealing temperature object
        """
        self._temperature.decreaseTemperature()
        
    def cooledDown(self):
        """
        cooledDown: Returns True if temperature <= 0
        """
        return self._temperature.cooledDown()
    
    def copy(self):
        """
        copy: Returns a copy of this object
        """
        new = MCOACJob()
        new.fromStr(str(self))
        return new