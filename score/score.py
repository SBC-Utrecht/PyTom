from pytom.basic.structures import PyTomClass

def Vol_G_Val(volume,value):
    """
    Vol_GE_Val: returns True when peak in volume greater than value
    @param volume: A volume
    @type volume: L{pytom_volume.vol}
    @param value: A value
    @type value: L{pytom_volume.vol}
    @return: True if peak in volume > value
    @rtype: boolean
    @author: Thomas Hrabe
    """
    import pytom_volume
    
    if value.__class__ == pytom_volume.vol:
        p = pytom_volume.peak(value)
        value = value.getV(p[0],p[1],p[2])
    
    if volume.__class__ == pytom_volume.vol:    
        p = pytom_volume.peak(volume)
        volume = volume.getV(p[0],p[1],p[2])
        
    return volume > value


def weightedCoefficient(self,volume,reference,mask=None,stdV=None):
    """
    weightedCoefficient: Determines the peak coefficient of the scoring function. 
    The distance from the center contributes to the peak value. Must be activated by hand. 
    @param volume: A volume.
    @type volume: L{pytom_volume.vol}
    @param reference: A reference.
    @type reference: L{pytom_volume.vol}
    @param mask: A mask.
    @type mask: L{pytom_volume.vol}
    @param stdV: Deviation volume of volume  
    @type stdV: L{pytom_volume.vol}
    @return: The highest coefficient determined.
    @author: Thomas Hrabe
    """

    resFunction = self.scoringFunction(volume,reference,mask,stdV)
    resFunction = self._peakPrior.apply(resFunction)
    
    return resFunction

def peakCoef(self,volume,reference,mask=None):
    """
    peakCoef: Determines the coefficient of the scoring function.
    @param volume: A volume.
    @type volume: L{pytom_volume.vol}
    @param reference: A reference.
    @type reference: L{pytom_volume.vol}
    @return: The highest coefficient determined.
    @author: Thomas Hrabe, FF
    """
    from pytom_volume import peak
    from pytom.tools.maths import euclidianDistance
    from pytom.basic.correlation import subPixelPeak
    
    if mask is None:
        resFunction = self.scoringFunction(volume,reference)
    else:
        resFunction = self.scoringFunction(volume,reference,mask)

    # change FF: 07.01.2020
    #centerX = resFunction.sizeX()//2 -1
    #centerY = resFunction.sizeY()//2 -1
    #centerZ = resFunction.sizeZ()//2 -1
    
    pcoarse = peak(resFunction)
    
    p = subPixelPeak(scoreVolume=resFunction, coordinates=pcoarse)
    #p = subPixelPeak(scoreVolume=resFunction, coordinates=[centerX,centerY,centerZ])
    
    #if euclidianDistance([centerX,centerY,centerZ],p[1]) <= 1.4142135623730951:
    #    c = p[0]

    return p[0]

def fromXML(xmlObj):
    """
    fromXML: Returns a score object depending on the string constant
    @param xmlObj: The DOM object provided
    @return: The score object 
    @author: Thomas Hrabe
    """
    
    type = xmlObj.get('Type')
    value = xmlObj.get('Value')
    
    if value == 'NAN':
        value = []
        
    if type == 'xcfScore':
        score = xcfScore()
    elif type == 'nxcfScore':
        score = nxcfScore()
    elif type == 'FLCFScore':
        score = FLCFScore()
    # make FLCF default score
    elif type == 'undefined':
        score = FLCFScore()
    elif type == 'SOCScore':
        score = SOCScore()
    elif type == 'RScore':    
        score = RScore()
        nb = xmlObj.get('NumberBands')
        wa = xmlObj.get('WedgeAngle')
        nb = nb
        wa = wa
        
        if nb.__class__ == str:
            nb = int(nb)
        if wa.__class__ == str:
            wa = int(wa)
                
        score.initAttributes(nb, wa)
    elif type == 'FSCScore':    
        score = FSCScore()
        nb = xmlObj.get('NumberBands')
        wa = xmlObj.get('WedgeAngle')
        nb = nb
        wa = wa
    
        if nb.__class__ == str:
            nb = int(nb)
        if wa.__class__ == str:
            wa = int(wa)
                
        score.initAttributes(nb, wa)
    elif type == 'FRMScore':
        from pytom.frm.FRMAlignment import FRMScore
        score = FRMScore()
    else:
        from pytom.basic.exceptions import ParameterError
        raise ParameterError('Type ' +type+' not available in pytom.alignment.score!')
    
    score.setValue(value)
    
    prElement = xmlObj.xpath('PeakPrior')

    if prElement == None:
        prElement = xmlObj.xpath('DistanceFunction')
        
    if len(prElement) >0:
        prObject = PeakPrior()
        prObject.fromXML(prElement[0])
        score._peakPrior = prObject
    else:  
        score._peakPrior = PeakPrior()

    removeAutocorrelation = xmlObj.get('RemoveAutocorr')
    if removeAutocorrelation == 'True':
        removeAutocorrelation = True
    else:
        removeAutocorrelation = False
    score.setRemoveAutocorrelation(flag=removeAutocorrelation)
        
    return score


def fromXMLFile(filename):
    """
    get score object from File
    """
    from lxml import etree
    from pytom.tools.files import readStringFile

    lines = readStringFile(filename)
    xmlObj = etree.fromstring(lines)
    score = fromXML(xmlObj)
    return score


class Score:
    """
    Score: Template class used for scoring alignments.  
    @author: Thomas Hrabe
    G{classtree}
    """
    def ctor(self,scoringFunction = 0,scoringCoefficient = 0,scoringCriterion = 0,scoreValue=[], removeAutocorr=False):
        """
        ctor : The constructor of this class
        @param scoringFunction: Will be assigned to self.scoringFunction - the true scoring function (XCF,nXCF,FLCF...)
        @param scoringCoefficient: Will be assigned to self.scoringCoefficient - function for finding the best scoring value.  
        @param scoringCriterion: Will be assigned to self.scoringCriterion  - function determining the best scoring value.
        @param scoreValue: Best score value determined by this scoring object.
        @param removeAutocorr: remove autocorrelation of average and particle?
        @type removeAutocorr: L{bool}
        @author: Thomas Hrabe, FF
        @change: removeAutocorr - FF
        """
        self.scoringFunction = scoringFunction
        self.scoringCoefficient = scoringCoefficient
        self.scoringCriterion = scoringCriterion
        self.setPeakPrior()
        self.scoreValue = scoreValue
        self.weightedCoefficient = weightedCoefficient
        self.setRemoveAutocorrelation( flag=removeAutocorr)
        self._type = 'undefined'
        
    def score(self,particle,reference,mask=None,stdV=None):
        """
        returns weighted Coefficient
        @param volume: A volume.
        @type volume: L{pytom_volume.vol}
        @param reference: A reference.
        @type reference: L{pytom_volume.vol}
        @param mask: A mask.
        @type mask: L{pytom_volume.vol}
        @param stdV: Deviation volume of volume  
        @type stdV: L{pytom_volume.vol}
        @return: The highest coefficient determined.
        """
        return self.weightedCoefficient(self,particle,reference,mask,stdV)
        
    def getScoreFunc(self):
        """
        getScoreFunc: return the used score function
        @author: chen
        """
        return self.scoringFunction
    
    def getValue(self):
        """
        getValue: Returns current score value
        @rtype: score coefficient
        """
        return self.scoreValue
        
    def setValue(self,value):
        """
        setValue: Sets score value determined by this scoring object.
        @param value:  The value.
        """
        if value.__class__ ==  list or value == 'NAN' or value == None:
            self.scoreValue = 'NAN'
        else:
            self.scoreValue = float(value)
    
    def getDistanceFunction(self):
        """
        @deprecated: Use getPeakPrior instead!
        """
        print('Deprecated: Use getPeakPrior instead!')
        return self._peakPrior
    
    def setDistanceFunction(self,filename='',mean=0,deviation=-1):
        """
        @deprecated: Use setPeakPrior instead!
        """
        print('Deprecated: Use setPeakPrior instead!')
        self.setPeakPrior(filename,mean,deviation)   
         
    def getPeakPrior(self):
        return self._peakPrior
    
    def setPeakPrior(self,filename='',radius=0,smooth=-1):
        """
        setPeakPrior : Defines the peak distance function for the scoring object. Shape 
        @param filename:
        @type filename: Either L{pytom.alignment.score.DistanceFunction} or a filename to a binary mask
        @param radius: Radius of ones around center
        @param smooth: Width of smoothed (1 -> 0) area past radius 
        @author: Thomas Hrabe 
        """
        if filename.__class__ == PeakPrior:
            self._peakPrior = filename
        else:
            self._peakPrior = PeakPrior(filename,radius,smooth)

    def setRemoveAutocorrelation(self, flag=False):
        """
        remove autocorrelation of average and particles, i.e., subtract particle from average prior to cc?
        @param flag: flag True/False
        @type flag: L{bool}
        """
        self._removeAutocorrelation = flag

    def getRemoveAutocorrelation(self):
        """
        @rtype: L{bool}
        """
        return self._removeAutocorrelation
    
    def toXML(self):
        """
        toXML : Compiles a XML file from this object 
        @author: Thomas Hrabe
        """
        from lxml import etree
        
        if self.scoreValue.__class__ == list:
            self.scoreValue = 'NAN'
            
        score_element = etree.Element("Score", Type = self._type, Value = str(self.scoreValue),
                                      RemoveAutocorr=str(self.getRemoveAutocorrelation()))
        
        prElement = self._peakPrior.toXML()
        
        score_element.append(prElement)
        
        return score_element

    def toXMLFile(self, filename):
        """
        write score to XML file
        """
        import pytom

        self._filename = filename
        versionString = '<!-- PyTom Version: ' + pytom.__version__ + ' -->\n'


        file = open(filename, "w")
        file.write(versionString + str(self))
        file.close()

    def fromXMLFile(self, filename):
        """
        fromXMLFile: Configures object according to file. Do NOT overload this function in any child!
        @param filename: Absolute / relative path to file
        @type filename: L{str}
        @author: Thomas Hrabe
        """
        from pytom.tools.files import readStringFile

        self._filename = filename
        lines = readStringFile(filename)

        self.fromStr(lines)
    
    def __str__(self):
        """
        __str__ : Prints object as XML String
        @author: Thomas Hrabe
        """    
        from lxml.etree import tostring
        
        doc = self.toXML()

        xmlstring =  tostring(doc, pretty_print=True).decode('utf-8')

        return xmlstring

    def getWorstValue(self):
        raise RuntimeError('This function must be overridden by child!')
        
    def fromStr(self,string):
        """
        fromStr : Creates object from XML String
        @param string: The XML String 
        @type string: str
        @author: Thomas Hrabe
        """    
        from lxml import etree
        
        root = etree.fromstring(string)
        
        self = fromXML(root)
    
    def getType(self):
        return self._type
    
    def __eq__(self,other):
        if super(other.__class__,other) == Score:
            return str(self) == str(other)
        else:
            return False
        
    def __ne__(self,other):
        return not self == other
        
    def __cmp__(self,otherScore):
        
        if not self._type == otherScore.getType():
            raise TypeError('Score types must be identical for comparison')
        
        if self.scoreValue < otherScore.getValue():
            return -1
        elif self.scoreValue == otherScore.getValue():
            return 0
        elif self.scoreValue > otherScore.getValue():
            return 1
    
    def __lt__(self, otherScore):
        self.__cmp__(otherScore)

    def __gt__(self, otherScore):
        if self.scoreValue > otherScore.getValue():
            return -1
        elif self.scoreValue == otherScore.getValue():
            return 0
        elif self.scoreValue < otherScore.getValue():
            return 1

    
    def copy(self):    
        selfAsXML = self.toXML()
        return fromXML(selfAsXML)
    
class xcfScore(Score):    
    """
    xcfScore: Uses the non-normalised correlation function for scoring
    @author: Thomas Hrabe
    """
    coefFnc = peakCoef

    def __init__(self,value=None):
        """
        __init__ : Assigns the xcf as scoringFunction, peakCoef as scoringCoefficient and Vol_G_Val as scoringCriterion
        @param value: Current value of score 
        """
        from pytom.basic.correlation import xcf,xcc
        self.ctor(xcf,xcc,Vol_G_Val)
        self._type = 'xcfScore'
        #if value and (isinstance(value, (int, long)) or value.__class__ == float):
        if value and (value.__class__ == int or value.__class__ == int or value.__class__ == float):
            self.setValue(value)
        else:
            self.setValue(self.getWorstValue())
            
    def getWorstValue(self):
        return -10000000000
    
class nxcfScore(Score):    
    """
    nxcfScore: Uses the normalised correlation function for scoring
    @author: Thomas Hrabe
    """
    
    def __init__(self,value=None):
        """
        __init__ : Assigns the normalized xcf (nxcf) as scoringFunction, peakCoef as scoringCoefficient and Vol_G_Val as scoringCriterion
        @param value: Current value of score
        """
        from pytom.basic.correlation import nXcf,nxcc
        self.ctor(nXcf,nxcc,Vol_G_Val)   
        self._type = 'nxcfScore'
        #if value and (isinstance(value, (int, long)) or value.__class__ == float):
        if value and (value.__class__ == int or value.__class__ == int or value.__class__ == float):
            self.setValue(value)
        else:
            self.setValue(self.getWorstValue())
            
    def getWorstValue(self):
        return -10000000000

class FLCFScore(Score):    
    """
    FLCFScore: Uses the FLCF correlation function for scoring
    @author: Thomas Hrabe
    """
    # coefFnc = peakCoef
    def __init__(self,value=None):
        """
        __init__ : Assigns the fast local correlation as scoringFunction, peakCoef as scoringCoefficient and Vol_G_Val as scoringCriterion
        @param value: Current value of score
        """
        from pytom.basic.correlation import FLCF
        from pytom.score.score import peakCoef

        self.ctor(FLCF, peakCoef, Vol_G_Val)
        self._type = 'FLCFScore'
        
        #if value and (isinstance(value, (int, long)) or value.__class__ == float):
        if value and (value.__class__ == int or value.__class__ == int or value.__class__ == float):
            self.setValue(value)
        else:
            self.setValue(self.getWorstValue())
            
    def getWorstValue(self):
        return -10000000000
    
def wXCCWrapper(self,volume,reference,mask=None):
    from pytom.basic.correlation import weightedXCC
        
    if self.getNumberOfBands() == 0:
        raise RuntimeError('RScore: Number of bands is Zero! Abort.')
    
    return weightedXCC(volume,reference,self.getNumberOfBands(),self.getWedgeAngle())

def wXCFWrapper(self,volume,reference,mask=None):
    from pytom.basic.correlation import weightedXCF
        
    if self.getNumberOfBands() == 0:
        raise RuntimeError('RScore: Number of bands is Zero! Abort.')
    
    return weightedXCF(volume,reference,self.bands,self.wedgeAngle)

class RScore(Score):
    """
    RScore: Uses the weighted correlation function for scoring. See 
    Stewart, A. 2004 Ultramicroscopy - Noise bias in the refinement of structures derived from single particles 
    for more info. Implementation of this class is a little more complicated. wXCFWrapper and wXCCWrapper 
    @author: Thomas Hrabe
    """
    
    def __init__(self,value=None):
        from pytom.basic.correlation import nXcf,weightedXCC
        self.ctor(nXcf,weightedXCC,Vol_G_Val)
        self._type = 'RScore'
        self._numberOfBands=0
        self._wedgeAngle =-1

        if value:
            self.setValue(value)
        else:
            self.setValue(self.getWorstValue())
            
    def getWorstValue(self):
        return -10000000000
    
    def getNumberOfBands(self):
        return self._numberOfBands
    
    def getWedgeAngle(self):
        return self._wedgeAngle
    
    
    
    def initAttributes(self,numberOfBands=10,wedgeAngle=-1):
        """
        initBands: Must be called prior to the scoring. Initializes the bands used for wxcf.
        @param numberOfBands:   
        @param wedgeAngle: 
        @author: Thomas Hrabe
        """
        self._numberOfBands = numberOfBands
        self._wedgeAngle = wedgeAngle
        
    def toXML(self,value=-10000000):
        """
        toXML : Compiles a XML file from this object
        @author: Thomas Hrabe
        """
        from lxml import etree
            
        score_element = etree.Element("Score",Type=self._type,Value = str(self.scoreValue))
        
        score_element.set('NumberBands',str(self._numberOfBands))
        score_element.set('WedgeAngle',str(self._wedgeAngle))
        
        return score_element

def FSCWrapper(self,volume,reference):
    from pytom.basic.correlation import FSCSum
        
    if self.numberOfBands == 0:
        from pytom.basic.exceptions import ParameterError
        raise ParameterError('Bands attribute is empty. Abort.')
        
    return FSCSum(volume,reference,self.bands,self.wedgeAngle)

def FSFWrapper(self,volume,reference):
    from pytom.basic.correlation import weightedXCF
        
    if self.numberOfBands == 0:
        from pytom.basic.exceptions import ParameterError
        raise ParameterError('Bands attribute is empty. Abort.')
        
    return weightedXCF(volume,reference,self.bands,self.wedgeAngle)   
  
class FSCScore(Score):
    """
    FSCScore: Uses the Sum of the Fourier Shell Correlation function for scoring. 
    @author: Thomas Hrabe
    """
    FSC = FSCWrapper
    FSF = FSFWrapper
    
    def __init__(self):
        from pytom.basic.correlation import nXcf
        self.ctor(nXcf,self.FSC,Vol_G_Val)
        self._type = 'FSCScore'
        self._numberOfBands = 0
        self._wedgeAngle = 0
        
    def initAttributes(self,numberOfBands=10,wedgeAngle=-1):
        """
        initBands: Must be called prior to the scoring. Initialises the bands used for wxcf.
        @param numberOfBands:   
        @param wedgeAngle: 
        @author: Thomas Hrabe
        """
        self._numberOfBands = numberOfBands
        self._bands = []
        
        for i in range(self._numberOfBands):
            self.bands.append([float(i)/self._numberOfBands*1/2,float(i+1)/self._numberOfBands*1/2])
            
        self._wedgeAngle = wedgeAngle
        
    def toXML(self,value=-10000000):
        """
        toXML : Compiles a XML file from this object
        @author: Thomas Hrabe
        """
        from lxml import etree
            
        score_element = etree.Element("Score",Type=self._type,Value = str(value))
        
        score_element.set('NumberBands',str(self._numberOfBands))
        score_element.set('WedgeAngle',str(self._wedgeAngle))
        
        return score_element

    def getWorstValue(self):
        return -10000000000
        
class PeakPrior(PyTomClass):
    """
    PeakPrior: Weights each correlation coefficient according it's peaks distance from the center
    """
    
    def __init__(self,filename='',radius=0.0,smooth=-1):
        """
        __init__ : Initialises attributes. Set either filename or mean > 0 and std > 0
        @param filename: Initialize either via existing file or  
        @param radius: Radius of ones around center
        @param smooth: Width of smoothed (1 -> 0) area past radius 
        @author: Thomas Hrabe  
        """
        self._weight = None
        self._filename = filename   
        self._radius = radius
        self._smooth = smooth
    
    def getFileName(self):
        return self._filename
    
    def getRadius(self):
        return self._radius
    
    def getSmooth(self):
        return self._smooth()
    
    def isInitialized(self):
        from pytom_volume import vol
        return self._weight.__class__ == vol
             
    def apply(self,volume):
        """
        apply: Applies weighting defined in this object to value. The return value can be modified if needed. 
        @param volume: A volume
        @return: self.weight * volume
        @author: Thomas Hrabe 
        """
        from pytom.tools.macros import volumesSameSize 
        from pytom.tools.files import checkFileExists
        
        if not self.isInitialized() and (not checkFileExists(self._filename)):
            self.initVolume(volume.sizeX(),volume.sizeY(),volume.sizeZ())
        elif not self.isInitialized():
            self.fromFile()
        
        assert volumesSameSize(self._weight,volume)#make sure both have same size

        from pytom_numpy import vol2npy

        return (self._weight * volume) - (self._weight <= 0.000001)*2
    
    def fromFile(self,filename=None):
        from pytom_volume import read
        from pytom.tools.files import checkFileExists
        
        filename = filename or self._filename
        
        if checkFileExists(filename):
            self._weight = read(filename)
        else:
            raise RuntimeError('PeakPrior: File ' + filename + ' not found')


    def reset_weight(self):
        self._weight = None


    def reset(self):
        self._mean = 0.0
        self._deviation = 0.0
        self._weight = None
        del(self._weight)

        
    def initVolume(self,sizeX,sizeY,sizeZ):
        """
        initVolume:
        @param sizeX:
        @param sizeY:
        @param sizeZ:
        @return: 
        @author: Thomas Hrabe 
        """
        from pytom_volume import vol,initSphere
        
        self._weight = vol(sizeX,sizeY,sizeZ)
        if self._radius > 0 or self._smooth > 0:
            initSphere(self._weight,self._radius,self._smooth,0.0, 
	        sizeX/2, sizeY/2, sizeZ/2)
        else:
            self._weight.setAll(1)
        
    def toXML(self):
        """
        toXML : Compiles a XML object from distanceFunction object
        @author: Thomas Hrabe
        """
        from lxml import etree
        
        prElement = etree.Element("PeakPrior", Filename=self._filename, Radius = str(float(self._radius)),
                                  Smooth = str(float(self._smooth)))
        
        return prElement
    
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

        if xmlObj.tag == 'PeakPrior':
            self._filename = xmlObj.get('Filename')
            self._radius = float(xmlObj.get('Radius'))
            self._smooth = float(xmlObj.get('Smooth'))
            
        elif xmlObj.tag == 'DistanceFunction':
            self._filename = xmlObj.get('Filename')
            self._radius = float(xmlObj.get('Mean'))
            self._smooth = float(xmlObj.get('Deviation'))

        else:
            Exception('Is not a PeakPrior XML object!')
        
        
        if self._filename != '':
            self.fromFile(self._filename)
        
    def __str__(self):
        """
        __str__ : Prints object as XML String
        @author: Thomas Hrabe
        """    
        from lxml import etree
        
        tree = self.toXML()
        
        return etree.tostring(tree,pretty_print=True).decode("utf-8")[:-1]
    
class SOCScore(Score):    
    """
    SOCScore: Uses the Second Order Correlation function for scoring
    @author: Thomas Hrabe
    """
    coefFnc = peakCoef
    
    def __init__(self):
        from pytom.basic.correlation import soc
        self.ctor(soc,self.coefFnc,Vol_G_Val)   
        self._type = 'SOCScore'

    def getWorstValue(self):
        return -10000000000.0
    
class MFCScore(Score):
    """
    MFCScore : Uses the Mutual Correlation Function for scoring
    @todo: Implementation
    @author: Thomas Hrabe
    """
    coefFnc = peakCoef
    
    def getWorstValue(self):
        return -10000000000
    
class POFScore(Score):
    """
    POFScore : Uses the Phase Only Correlation Function for scoring
    @todo: Implementation
    @author: Thomas Hrabe
    """
    coefFnc = peakCoef
    
    
    def getWorstValue(self):
        return -10000000000
    
