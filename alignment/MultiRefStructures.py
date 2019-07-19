'''
Created on Jan 2, 2012

@author: hrabe
'''
from pytom.cluster.mcoEXMXStructures import MCOEXMXJob

class MCOEXMXAlignJob(MCOEXMXJob):
    """
    MCOEXMXAlign:
    Multi reference alignment job based on the MCOEXMX clustering. Difference to L{pytom.cluster.mcoEXMXStructures.MCOEXMXJob} is that 
    anykind of angle list L{pytom.angles.AngleObject} will be looped to determine the best alignment of each particle and reference.    
    """
    
    def __init__(self,particleList,numberIterations,destinationDirectory,mask,score,preprocessing,wedgeInfo,binning,sampleInformation,numberClasses,endThreshold,symmetry = None,rotationList=None,adaptiveResolution=True,fscCriterion = 0.5, adaptiveOffset=0.1,angleFactor=0.5,useMaxResolution=True):
        
        if rotationList == None:
           self = MCOEXMXJob(particleList,numberIterations,destinationDirectory,mask,score,preprocessing,wedgeInfo,binning,sampleInformation,numberClasses,endThreshold,symmetry) 
        else:
            super(self.__class__,self).__init__(particleList,numberIterations,destinationDirectory,mask,score,preprocessing,wedgeInfo,binning,sampleInformation,numberClasses,endThreshold,symmetry)
            
            self._exMaxJob.setAdaptiveResolution(adaptiveResolution)
            self._exMaxJob.setFSCCriterion(fscCriterion)
            self._exMaxJob.setAdaptiveOffset(adaptiveOffset)
            self._exMaxJob.setAngleFactor(angleFactor)
            self._exMaxJob.setRotations(rotationList)
                        
            self._useMaxResolution = useMaxResolution
    
    def getUseMaxResolution(self):
        return self._useMaxResolution
    
    def setUseMaxResolution(self,useMaxResolution):
        self._useMaxResolution = useMaxResolution
        
    def toXML(self):
        from lxml import etree
        
        job_element = super(self.__class__,self).toXML()
        job_element.tag = 'MCOEXMXAlignJob'
        
        job_element.set('UseMaxResolution',str(self._useMaxResolution))
            
        return job_element
    
    def fromXML(self,xmlObj):
        
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element :
            raise TypeError('Is not a lxml.etree._Element! You must provide a valid XMLobject.')
       
        if not xmlObj.tag == 'MCOEXMXAlignJob':
            raise TypeError('You must provide a MCOEXMXAlignJob XML object!')
        
        self.useMaxResolution = xmlObj.get('UseMaxResolution') == 'True'
        
        xmlObj.tag = 'MCOEXMXJob'
        super(self.__class__, self).fromXML(xmlObj)
        
        
    
    
    