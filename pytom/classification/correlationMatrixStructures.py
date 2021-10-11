'''
Created on Mar 24, 2010

@author: hrabe
'''

from pytom.basic.structures import PyTomClass 

class CorrelationVector(PyTomClass):
    """
    CorrelationVector: Stores a vector of correlation values of one particle with many others from particle list
    @todo: Add unittest
    """
    
    def __init__(self,particle=None,particleList=None,particleIndex=None):
        self._correlations = []
        self._particle = particle
        self._particleList = particleList
        self._particleIndex = particleIndex
        
    def append(self,value):
        self._correlations.append(value)
    
    def getParticleIndex(self):
        return self._particleIndex
    
    def __len__(self):
        return len(self._correlations)
    
    def __getitem__(self,key):
        """
        """
        if isinstance(key, int):
            if key < len(self):
                return self._correlations[key]
            else:
                raise IndexError('Index out of range.')
        else:
            assert False
            
    def toXML(self): 
        """
        toXML : Compiles a XML file from result object
        rtype : L{lxml.etree._Element}
        @author: Thomas Hrabe
        """        
        from lxml import etree

        vectorElement = etree.Element("CorrelationVector",ParticleIndex = self._particleIndex.__str__())
        
        vectorElement.append(self._particle.toXML())
        vectorElement.append(self._particleList.toXML())
        
        for i in range(len(self._correlations)):
            
            valueElement = etree.Element("Correlation",Index=i.__str__(),Value=self._correlations[i].__str__())
            
            vectorElement.append(valueElement)

        return vectorElement

    def fromXML(self,xmlObj):
        """
        fromXML : Assigns values to job attributes from XML object
        @param xmlObj: A xml object  
        @type xmlObj: L{lxml.etree._Element}
        @author: Thomas Hrabe 
        """
        from lxml.etree import _Element,tostring
        
        if xmlObj.__class__ != _Element:
            raise ParameterError('You must provide a valid XML-CorrelationVector object.')
        
        from pytom.basic.structures import Particle,ParticleList
        
        self._particleIndex = int(xmlObj.get('ParticleIndex'))
        
        particleObject = xmlObj.xpath('Particle')
        self._particle = Particle('.')
        self._particle.fromXML(particleObject[0])
        
        particleListObject = xmlObj.xpath('ParticleList')
        self._particleList = ParticleList('/')
        self._particleList.fromXML(particleListObject[0])
        
        values = xmlObj.xpath('Correlation')
        
        self._correlations = [0  for _ in range(len(values))]
        
        for v in values:
            
            index = int(v.get('Index'))
            self._correlations[index] = float(v.get('Value'))
             
class CorrelationVectorJob(PyTomClass):
    """
    CorrelationVectorJob: All settings needed for a correlation vector job. Explore class for more info.
    """
    
    def __init__(self,particle=None,particleList=None,mask=None,particleIndex = None,applyWedge = True,binningFactor=0,lowestFrequency=-1,highestFrequency=-1):
        """
        __init__:
        @param particle: Particle
        @type particle: pytom.alignment.structures.Particle
        @param particleList: ParticleList of all particles will be correlated with self._particle
        @type particleList: pytom.alignment.structures.ParticleList
        @param mask: Mask used for correlation
        @type mask: str
        @param applyWedge: Apply wedge during correlation. True by default, disable for rotation classification.
        @type applyWedge: Bool  
        @param binningFactor: Binning factor accroding to libtomc definition. 0 by default.
        @type binningFactor: unsigned int
        @param lowestFrequency: Lowest frequency for bandpass in nyquist
        @type lowestFrequency: float
        @param highestFrequency: Highest frequency for bandpass in nyquist
        @type highestFrequency: float  
        """
        from pytom.basic.structures import Particle,ParticleList,Mask
        
        if particle and particle.__class__ != Particle:
            raise ParameterError('You must provide a Particle object.')
        
        if particleList and particleList.__class__ != ParticleList:
            raise ParameterError('You must provide a ParticleList object.')
        
        self._particle = particle
        self._particleList = particleList
        
        if not mask:
            mask = Mask()
        elif mask.__class__ == str:
            mask = Mask(mask)
        elif mask.__class__ != Mask:
            from pytom.basic.exceptions import ParameterError
            raise ParameterError('Mask must be a string or Mask object!')
        
        self._mask = mask
        self._particleIndex = particleIndex
        self._applyWedge = applyWedge
        self._binningFactor = binningFactor
        self._lowestFrequency = lowestFrequency
        self._highestFrequency = highestFrequency
        
    def getMask(self):
        return self._mask
    
    def getParticle(self):
        return self._particle
    
    def getParticleList(self):
        return self._particleList
    
    def getParticleIndex(self):
        return self._particleIndex
    
    def getApplyWedge(self):
        return self._applyWedge
    
    def getBinning(self):
        return self._binningFactor
    
    def getLowestFrequency(self):
        return self._lowestFrequency
    
    def getHighestFrequency(self):
        return self._highestFrequency
    
    def toXML(self):
        """
        toXML : Compiles a XML object from result object
        rtype : L{lxml.etree._Element}
        @author: Thomas Hrabe
        """
        
        from lxml import etree
        
        jobElement = etree.Element('CorrelationVectorJob',ParticleIndex = self._particleIndex.__str__(),ApplyWedge=self._applyWedge.__str__(),Binning = self._binningFactor.__str__(),LowestFrequency = str(self._lowestFrequency),HighestFrequency = str(self._highestFrequency))
        
        jobElement.append(self._particle.toXML())
        jobElement.append(self._particleList.toXML())
        jobElement.append(self._mask.toXML())
        
        return jobElement
        
    def fromXML(self,xmlObj):
        """
        fromXML : Assigns values to job attributes from XML object
        @param xmlObj: A xml object  
        @type xmlObj: L{lxml.etree._Element}
        @author: Thomas Hrabe 
        """
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element:
            raise ParameterError('You must provide a valid XML-CorrelationVectorJob object.')
        
        
        from pytom.basic.structures import Particle,ParticleList,Mask
        
        particleObject = xmlObj.xpath('Particle')
        self._particle = Particle('.')
        self._particle.fromXML(particleObject[0])
        
        particleListObject = xmlObj.xpath('ParticleList')
        self._particleList = ParticleList('/')
        self._particleList.fromXML(particleListObject[0])
    
        maskObject = xmlObj.xpath('Mask')[0]
        self._mask = Mask()
        self._mask.fromXML(maskObject)
        
        self._particleIndex = xmlObj.get('ParticleIndex')
        self._applyWedge = xmlObj.get('ApplyWedge') == 'True'
        self._binningFactor = float(xmlObj.get('Binning'))
        
        self._lowestFrequency = float(xmlObj.get('LowestFrequency'))
        self._highestFrequency = float(xmlObj.get('HighestFrequency'))
        
class CorrelationMatrixJob(PyTomClass):
    """
    CorrelationVectorJob: Represents all settings for a correlation matrix job. Explore the class for further information
    """
    
    def __init__(self,particleList=None,mask='',resultMatrixName='',applyWedge = True,binningFactor=0,lowestFrequency=-1,highestFrequency=-1):
        """
        __init__:
        @param particleList: ParticleList of all particles that will be correlated 
        @type particleList: L{pytom.basic.structures.ParticleList}
        @param mask: Mask used for correlation
        @type mask: str or L{pytom.basic.structures.Mask}
        @param resultMatrixName: Result filename
        @type resultMatrixName: str  
        @param applyWedge: Apply wedge weighting if available?
        @type applyWedge: boolean, False by default  
        @param binningFactor: Binning factor accroding to libtomc definition. 0 by default.
        @type binningFactor: unsigned int
        @param lowestFrequency: Lowest frequency for bandpass in nyquist
        @type lowestFrequency: float
        @param highestFrequency: Highest frequency for bandpass in nyquist
        @type highestFrequency: float
        """
        from pytom.basic.structures import ParticleList,Mask
        from pytom.tools.files import checkFileExists
        
        if not particleList:
            particleList = ParticleList('/')
        elif particleList.__class__ != ParticleList:
            raise ParameterError('You must provide a ParticleList object!')
        
        self._particleList = particleList
        
        if mask.__class__ == str:
            mask = Mask(mask)
        elif mask.__class__ != Mask:
            raise ParameterError('Mask must be a string or Mask object!')
        
        self._mask = mask
        self._resultMatrixName = resultMatrixName
        self._applyWedge = applyWedge
        self._binningFactor = binningFactor
        self._lowestFrequency = lowestFrequency
        self._highestFrequency = highestFrequency
        
    def getParticleList(self):    
        return self._particleList
    
    def getMask(self):
        return self._mask
    
    def getResultMatrixName(self):
        return self._resultMatrixName
    
    def getApplyWedge(self):
        return self._applyWedge
    
    def getBinning(self):
        return self._binningFactor
    
    def getLowestFrequency(self):
        return self._lowestFrequency
    
    def getHighestFrequency(self):
        return self._highestFrequency
    
    def toXML(self):
        """
        toXML : Compiles a XML object from result object
        rtype : L{lxml.etree._Element}
        @author: Thomas Hrabe
        """
        from lxml import etree
        
        jobElement = etree.Element('CorrelationMatrixJob',ResultMatrixName = self._resultMatrixName, ApplyWedge=self._applyWedge.__str__(),Binning=self._binningFactor.__str__(),LowestFrequency = str(self._lowestFrequency),HighestFrequency = str(self._highestFrequency))
        
        jobElement.append(self._particleList.toXML())
        
        jobElement.append(self._mask.toXML())
        
        return jobElement
        
    def fromXML(self,xmlObj):
        """
        fromXML : Assigns values to job attributes from XML object
        @param xmlObj: A xml object  
        @type xmlObj: L{lxml.etree._Element}
        @author: Thomas Hrabe 
        """
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element:
            raise ParameterError('You must provide a valid XML-CorrelationMatrixJob object.')
        
        
        from pytom.basic.structures import ParticleList,Mask
        
        particleListObject = xmlObj.xpath('ParticleList')
        self._particleList = ParticleList('/')
        self._particleList.fromXML(particleListObject[0]) 
        
        maskObject = xmlObj.xpath('Mask')[0]
        self._mask = Mask()
        self._mask.fromXML(maskObject)
        
        self._resultMatrixName = xmlObj.get('ResultMatrixName')
        self._applyWedge = xmlObj.get('ApplyWedge') == 'True'
        self._binningFactor = float(xmlObj.get('Binning'))
        
        self._lowestFrequency = float(xmlObj.get('LowestFrequency'))
        self._highestFrequency = float(xmlObj.get('HighestFrequency'))
        
    def toHTMLFile(self,filename):
        """
        toHTMLFile: Overrides parent method and stores CorrelationMatrixJob to HMTL
        @param filename: HTML filename
        """
        from lxml import etree
        import io
        from pytom.tools.files import getPytomPath,readStringFile
        
        pytomPath = getPytomPath()
        
        xsltString = readStringFile(pytomPath + '/xslt/CorrelationMatrix.xsl')
        
        xsltStringIO = io.StringIO(xsltString)
        
        selfHTML = self.xsltTransform(xsltStringIO)
        
        htmlString = etree.tostring(selfHTML,pretty_print=True)

        file = open(filename, "w")

        file.write(htmlString)
    
        file.close()
        