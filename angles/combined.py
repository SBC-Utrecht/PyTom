'''
Created on Nov 21, 2011

@author: hrabe
'''


from pytom.angles.angle import AngleObject


class GlobalLocalCombined(AngleObject):
    """
    GlobalLocalCombined: Switches between global and local sampling lists. Starts with global sampling, refines with local sampling, repeats global sampling again
    """
    
    def __init__(self,globalSampling=None,localSampling=None,startGlobal = True):
        """
        """
        from pytom.angles.globalSampling import GlobalSampling
        from pytom.angles.localSampling import LocalSampling
        if not globalSampling.__class__ == GlobalSampling and globalSampling:
            raise TypeError('GlobalLocalCombined : globalSampling must be of type GlobalSampling or None.')
        else:
            self._globalSampling = globalSampling
        
        if not localSampling.__class__ == LocalSampling and localSampling:
            raise TypeError('GlobalLocalCombined : localSampling must be of type LocalSampling or None.')
        else:
            self._localSampling = localSampling
        
        self._currentSamplingIsGlobal = startGlobal
        
    def toXML(self):
        from lxml import etree
               
        angles_element = etree.Element("Angles",Type = 'Combined',currentIsGlobal = str(self._currentSamplingIsGlobal))
        angles_element.append(self._globalSampling.toXML())
        angles_element.append(self._localSampling.toXML())
        
        return angles_element
    
    def fromXML(self,xmlObj):
        from lxml.etree import _Element
        from pytom.angles.globalSampling import GlobalSampling
        from pytom.angles.localSampling import LocalSampling
        
        if xmlObj.__class__ != _Element:
            raise TypeError('You must provide a valid XML-AngleEquidistant object.')
        
        self._currentSamplingIsGlobal = bool(xmlObj.get('currentIsGlobal')) 
        
        globalXML = xmlObj.xpath("Angles[@Type='GlobalSampling']")
        if len(globalXML) == 0:
            globalXML = xmlObj.xpath("Angles[@Type='FromEMFile']")
        globalSampling = GlobalSampling()
        globalSampling.fromXML(globalXML[0])
        self._globalSampling = globalSampling
        
        localXML = xmlObj.xpath("Angles[@Type='LocalSampling']")
        if len(localXML) == 0:
            localXML = xmlObj.xpath("Angles[@Type='EquidistantList']")
        localSampling = LocalSampling()
        localSampling.fromXML(localXML[0])
        self._localSampling = localSampling
    
    def nextRotation(self):
        if self._currentSamplingIsGlobal:
            return self._globalSampling.nextRotation()
        else:
            return self._localSampling.nextRotation()
        
    def focusRotation(self,rotation=None,refinementAngle=None):
        """
        focusRotation: create list of rotations centered around chosen rotation \
	(takes self._numberShells in addition to specified parameters

        @param rotation: rotation that defines center of local sampling
	@type rotation: list
        @param refinementAngle: angular increment for search
	@type refinementAngle: float (or int)
	@return: LocalSampling
        """
        if self._currentSamplingIsGlobal:
            from pytom.angles.localSampling import LocalSampling
            if not rotation:
                rotation = [0,0,0]
        
            if refinementAngle == None:
                refinementAngle = 10
                
            newLocalSampling = LocalSampling(self._localSampling.getNumberShells(),refinementAngle,rotation[0],rotation[1],rotation[2])
            
            return GlobalLocalCombined(self._globalSampling,newLocalSampling, startGlobal = False)    
        else:
            return GlobalLocalCombined(self._globalSampling,self._localSampling, startGlobal = True)
        
    def reset(self):
        self._globalSampling.reset()
        self._localSampling.reset()
        
    def setStartRotation(self,startRotation):    
        from pytom.basic.structures import Rotation
        if not startRotation.__class__ == Rotation:
            raise TypeError('Angles.setStartRotation requires startRotation to be of type Rotation')
    
        if self._currentSamplingIsGlobal:
            from pytom.angles.localSampling import LocalSampling
            newLocalSampling = LocalSampling(self._localSampling.getNumberShells(), self._localSampling.getIncrement(), startRotation[0], startRotation[1] , startRotation[2])
            return GlobalLocalCombined(self._globalSampling,newLocalSampling, False)
        else:
            return GlobalLocalCombined(self._globalSampling,self._localSampling, True)
            
    def getIncrement(self):
        """
        getIncrement:
        """
        if self._currentSamplingIsGlobal:
            return self._globalSampling.getIncrement()
        else:
            return self._localSampling.getIncrement()
        
