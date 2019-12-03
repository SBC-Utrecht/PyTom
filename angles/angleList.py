from pytom.angles.angle import AngleObject
from numpy import arange

class AngleList(AngleObject):
    """
    AngleList : Derives from AngleObject. Can be used to perform calculations on arbitrary rotation lists.
    @ivar _rotationList: The list stored
    @ivar _currentIndex: Current index used in nextRotation
    @author: Thomas Hrabe
    """
    
    _rotationList = []
    _currentIndex = 0
    
    def __init__(self,rotationList = None):
        """
        @param rotationList: Either [] or a list of rotations - [ [phi,psi,theta],  [phi,psi,theta], [phi,psi,theta] ,...]
        """
        self._rotationList = rotationList or []
        
    def append(self,phi,psi=[],theta=[]):
        """
        append: Appends a combination of Euler Angles to the current list
        @param phi: Euler Angle or Rotation object or list object
        @param psi: Euler Angle
        @param theta: Euler Angle
        """
        from pytom.basic.structures import Rotation
        
        if phi.__class__ == Rotation:
            self._rotationList.append(phi)
        elif phi.__class__ == list:
            self._rotationList.append(Rotation(phi[0],phi[1],phi[2]))
        else:
            self._rotationList.append(Rotation(phi,psi,theta))
        
    def nextRotation(self):
        """
        nextRotation: Returns the next rotation for the current object

        @return: rotation
        @rtype: L{pytom.basic.structures.Rotation}

        @author: Thomas Hrabe 
        """
        if self._currentIndex >= self.numberRotations():
            return [None,None,None]
        
        rotation = self._rotationList[self._currentIndex]
        self._currentIndex = self._currentIndex + 1
        
        if rotation.__class__ == list:
            return rotation
        else:
            return rotation.toList()

    def toXML(self):
        from lxml import etree
        
        job_element = etree.Element('Angles',Type='AngleList')
        
        for i in range(self.numberRotations()):
            
            rotation = self._rotationList[i]
            if rotation.__class__ == list:
                rotationElement = etree.Element('Rotation',Index = str(i),Z1 = str(rotation[0]),Z2 = str(rotation[1]),X = str(rotation[2]))
            else:
                rotationElement = rotation.toXML()
                rotationElement.set('Index',i.__str__())
                
            job_element.append(rotationElement)
        
        return job_element
    
    def fromXML(self,xmlObj):
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element :
            raise RuntimeError('Is not a lxml.etree._Element! You must provide a valid XMLobject.')
        
        if xmlObj.tag == 'Angles':
            angleList_element = xmlObj
        else:
            RuntimeError('Is not an AngleList! You must provide a valid AngleList object.')
    
        from pytom.basic.structures import Rotation
    
        rotations = angleList_element.xpath('Rotation') # modified by chen
        
        for rotation in rotations:
            rot = Rotation(0,0,0)
            rot.fromXML(rotation)
            self._rotationList.append(rot)
            
    def numberRotations(self):
        return len(self._rotationList)

    def __getitem__(self,key):
        if isinstance(key, int):
            
            if key < self.numberRotations():
                return self._rotationList[key]
            else:
                raise IndexError('Index out of range.')
            
        elif key.__class__ == slice:
            
            start = key.start
            stop = key.stop
            step = key.step
            if not start:
                start = 0

            if not stop:
                stop = self.numberRotations()

            if not step:
                step = 1
            if stop >= 922337203685477580:
                stop = self.numberRotations()
                
            rotations = []
            
            for r in self._rotationList[start:stop:step]:
                rotations.append(r)
                
            return rotations
        else:
            assert False
    
    def reset(self):
        """
        reset: Resets self._currentIndex
        """
        self._currentIndex = 0
        
    def setStartRotation(self,rotation):
        """
        """
        from pytom.basic.structures import Rotation
        
        if rotation.__class__ == list:
            rotation = Rotation(rotation)
        
        if not rotation.__class__ == Rotation:
            raise TypeError('You must provide a Rotation to this function!')
        
        if self._rotationList[0] == rotation:
            return self
         
        self.append(rotation)
        
        if len(self) > 1:
            self._rotationList[len(self)-1] = self[0]
            self._rotationList[0] = rotation
    
        return self

    
class OneAngleList(AngleList):
    
    def __init__(self,rotationList = None):
        """
        @param rotationList: [phi,psi,theta] or one Rotation 
        """
        from pytom.basic.structures import Rotation
        
        self._rotationList = rotationList or []
         
        if self._rotationList.__class__ == list and len(self._rotationList) == 3:
            #make sure its only one object
            self._rotationList = [Rotation(rotationList[0],rotationList[1],rotationList[2])]
        elif self._rotationList.__class__ == list and len(self._rotationList) > 1:
            self._rotationList = [Rotation(rotationList[0])]
        elif self._rotationList.__class__ == Rotation:
            self._rotationList = [self._rotationList]
            
    def append(self,phi,psi=[],theta=[]):
        """
        append: Appends a combination of Euler Angles to the current list
        @param phi: Euler Angle or Rotation object or list object
        @param psi: Euler Angle
        @param theta: Euler Angle
        """
        from pytom.basic.structures import Rotation
        
        if phi.__class__ == Rotation:
            self._rotationList = [phi]
        elif phi.__class__ == list:
            self._rotationList = [Rotation(phi[0],phi[1],phi[2])]
        else:
            self._rotationList = [Rotation(phi,psi,theta)]
    
    def setStartRotation(self,rotation):
        """
        """
        from pytom.basic.structures import Rotation
        
        if rotation.__class__ == list:
            rotation = Rotation(rotation)
        
        if not rotation.__class__ == Rotation:
            raise TypeError('You must provide a Rotation to this function!')
        
        self._rotationList = [rotation]
        
        return self

    def toXML(self):
        from lxml import etree
        
        job_element = etree.Element('Angles',Type='OneAngleList')
        
        
        rotation = self._rotationList[0]
        
        if rotation.__class__ == list:
            rotationElement = etree.Element('Rotation',Index = str(i),Z1 = str(rotation[0]),Z2 = str(rotation[1]),X = str(rotation[2]))
        else:
            rotationElement = rotation.toXML()
            rotationElement.set('Index',str(0))
            
        job_element.append(rotationElement)
    
        return job_element
    
    def fromXML(self,xmlObj):
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element :
            raise RuntimeError('Is not a lxml.etree._Element! You must provide a valid XMLobject.')
        
        if xmlObj.tag == 'Angles':
            angleList_element = xmlObj
        else:
            RuntimeError('Is not an AngleList! You must provide a valid AngleList object.')
    
        from pytom.basic.structures import Rotation
    
        rotations = angleList_element.xpath('Rotation') # modified by chen
        
        for rotation in rotations:
            rot = Rotation(0,0,0)
            rot.fromXML(rotation)
            self._rotationList.append(rot)

def angleRange(start,end,inc):

    dif = (end-start)
    if dif <0:
        r = list(arange(start,360,inc))
        r1 = list(arange(0,end,inc))
        r.extend(r1)
    else:
        r = list(arange(start,start+dif,inc))
        
        
    return r
    
class EulerAngleList(AngleObject):
    """
    EulerAngleList : Represents a list of eulerian angles.  
    """
    
    def __init__(self,phiStart=0,psiStart=0,thetaStart=0,phiInc=10,psiInc=10,
            thetaInc=10,phiEnd=360,psiEnd=360,thetaEnd=360):
        """
        __init__: Creates a simple list stored in the class attributes for each eulerian angle.
        @param phiStart: first phi value (in deg)
        @param psiStart: first psi value
        @param thetaStart: first theta value
	@param phiInc: phi increment
	@param psiInc: psi increment
	@param thetaInc: theta increment
        @param phiEnd: last phi value (in deg)
        @param psiEnd: last psi value
        @param thetaEnd: last theta value

        @author: Thomas Hrabe
        """
        
        from pytom.angles.angleList import angleRange
        
        self._eulerStart=[phiStart,psiStart,thetaStart]
        self._eulerInc=[phiInc,psiInc,thetaInc]
        self._eulerEnd=[phiEnd,psiEnd,thetaEnd]
        
        self._thetaRange = angleRange(thetaStart,thetaEnd,thetaInc)
        self._psiRange = angleRange(psiStart,psiEnd,psiInc)
        self._phiRange = angleRange(phiStart,phiEnd,phiInc)
        
        for i in range(0,len(self._phiRange)):
            self._phiRange[i] = self._phiRange[i] % 360

        for i in range(0,len(self._thetaRange)):
            self._thetaRange[i] = self._thetaRange[i] % 180     
        
        for i in range(0,len(self._psiRange)):
            self._psiRange[i] = self._psiRange[i] % 360
            
    
        self._currentPhiIndex   =0
        self._currentPsiIndex   =0
        self._currentThetaIndex =0
        
        self._finished = False

        
    def nextRotation(self):
        """
        nextRotation : 
        @return: Either [None,None,None] if all rotations were scanned or [phi psi theta] for the next rotation
        @author: Thomas Hrabe
        """
        
        if self._finished :
            return [None,None,None]
        
        
        phi = self._phiRange[self._currentPhiIndex]
        psi = self._psiRange[self._currentPsiIndex]
        theta = self._thetaRange[self._currentThetaIndex]
        
        self._currentThetaIndex = self._currentThetaIndex + 1
        
        if self._currentThetaIndex >= len(self._thetaRange):
            self._currentThetaIndex = 0
            self._currentPsiIndex = self._currentPsiIndex + 1
            if self._currentPsiIndex >= len(self._psiRange):
                self._currentPsiIndex = 0
                self._currentPhiIndex = self._currentPhiIndex + 1
                if self._currentPhiIndex >= len(self._phiRange):
                    self._finished = True
                    
                    return [None,None,None]                  
            
        return [phi, psi, theta]
    
    def _findClosestAngle(self,angle,angleRange,angleName):
        """
        __findClosestAngle:
        @param angle:
        @param angleRange:
        @param angleName:
        @return:
        @author: Thomas Hrabe   
        """
        dist = 1000000
        
        rangeLength = len(angleRange)
        
        if rangeLength == 1:
            return 0
        
        for iterator in range(rangeLength):
            if abs(angle - angleRange[iterator]) < dist:
                dist = abs(angle - angleRange[iterator])
                pos = iterator 
        increment = abs(angleRange[0] - angleRange[1])
        
        if dist > increment:
            print(angleName , ' value : ' , angle , ' ' , angleRange)
            raise AngleError(' The ' + angleName + ' value was not found in the range of this object!')
        
        return pos
        
    def focusRotation(self,rotation): 
        """
        focusRotation: create list of rotations centered around chosen rotation \
	(takes self._numberShells in addition to specified parameters

        @param rotation: rotation that defines center of local sampling
	@type rotation: list
	@return: EulerAngleList
        """    
        phi = self._findClosestAngle(rotation[0],self._phiRange,'phi')
        
        psi = self._findClosestAngle(rotation[1],self._psiRange,'psi')
        
        theta = self._findClosestAngle(rotation[2],self._thetaRange,'theta')
        
        if phi>0:
            phiStart=self._phiRange[phi-1]
        else:
            phiStart=(self._phiRange[phi]-self._eulerInc[0])%360
                    
        if phi < len(self._phiRange)-2:
            phiEnd=self._phiRange[phi+1]
        else:
            phiEnd=(self._phiRange[phi]+self._eulerInc[0])%360

        if psi>0:
            psiStart=self._psiRange[psi-1]
        else:
            psiStart=(self._psiRange[psi]-self._eulerInc[1])%360
                    
        if psi < len(self._psiRange)-2:
            psiEnd=self._psiRange[psi+1]
        else:
            psiEnd=(self._psiRange[psi]+self._eulerInc[1])%360
                    
        if theta>0:
            thetaStart=self._thetaRange[theta-1]
        else:
            thetaStart=(self._thetaRange[theta]-self._eulerInc[2])%180
                    
        if theta < len(self._thetaRange)-2:
            thetaEnd=self._thetaRange[theta+1]
        else:
            thetaEnd=(self._thetaRange[theta]+self._eulerInc[2])%180
        
        phiInc = self._eulerInc[0]/2
        #if phiInc <1:
        #    phiInc =1
        psiInc = self._eulerInc[1]/2
        #if psiInc <1:
        #    psiInc =1
        thetaInc = self._eulerInc[2]/2
        #if thetaInc <1:
        #    thetaInc =1
        
        if len(self._phiRange) == 1:
            phiStart = 0
            phiInc   = 2
            phiEnd   = 1
        
        if len(self._psiRange) == 1:
            psiStart = 0
            psiInc   = 2
            psiEnd   = 1
            
        if len(self._thetaRange) == 1:
            thetaStart = 0
            thetaInc   = 2
            thetaEnd   = 1

        if phiInc == 0:
            phiInc = 1
            
        if psiInc == 0:
            psiInc = 1
            
        if thetaInc == 0:
            thetaInc = 1
        
        
        ang = EulerAngleList(phiStart,psiStart,thetaStart,phiInc,psiInc,thetaInc,phiEnd,psiEnd,thetaEnd)
        return ang
    
    def toXML(self):
        """
        toXML : Compiles a XML attribute from the Object
        @author: Thomas Hrabe
        """
        from lxml import etree
               
        angles_element = etree.Element("Angles",Type = 'Eulerian')
        
        start_element = etree.Element("Start",phi = self._eulerStart[0].__str__(),psi = self._eulerStart[1].__str__(), theta = self._eulerStart[2].__str__())
        
        increment_element = etree.Element( "Increment",phi = self._eulerInc[0].__str__(),psi = self._eulerInc[1].__str__(),theta = self._eulerInc[2].__str__())
        
        end_element = etree.Element("End",phi = self._eulerEnd[0].__str__(), psi = self._eulerEnd[1].__str__(),theta = self._eulerEnd[2].__str__())
        
        angles_element.append(start_element)
        angles_element.append(increment_element)
        angles_element.append(end_element)
        
        
        return angles_element
    
    def fromXML(self,xmlObj):
        """
        fromXML: Sets attributes of self to those defined by xmlObj
        @param xmlObj: A XML attribute as generated by toXML 
        @author: Thomas Hrabe
        """
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element:
            raise ParameterError('You must provide a valid XML-EulerAngleList object.')
        
        startNode = xmlObj.xpath('Start')
        startNode = startNode[0]
        self._eulerStart[0] = int(startNode.get("phi"))
        self._eulerStart[1] = int(startNode.get("psi"))
        self._eulerStart[2] = int(startNode.get("theta"))
        
        incrementNode = xmlObj.xpath('Increment')
        incrementNode = incrementNode[0]
        self._eulerInc[0] = int(incrementNode.get("phi"))
        self._eulerInc[1] = int(incrementNode.get("psi"))
        self._eulerInc[2] = int(incrementNode.get("theta"))
        
        endNode = xmlObj.xpath('End')
        endNode = endNode[0]
        self._eulerEnd[0] = int(endNode.get("phi"))
        self._eulerEnd[1] = int(endNode.get("psi"))
        self._eulerEnd[2] = int(endNode.get("theta"))
        
        self._currentPhiIndex   =0
        self._currentPsiIndex   =0
        self._currentThetaIndex =0
        

        self._phiRange = angleRange(self._eulerStart[0],self._eulerEnd[0],self._eulerInc[0])
        self._psiRange = angleRange(self._eulerStart[1],self._eulerEnd[1],self._eulerInc[1])
        self._thetaRange = angleRange(self._eulerStart[2],self._eulerEnd[2],self._eulerInc[2])
        
        
        self.finished = False

        
    def numberRotations(self):
        """
        numberRotations: Returns the total number of rotations for the current object
        @author: Thomas Hrabe 
        """
        otherObj = EulerAngleList(self._eulerStart[0],self._eulerStart[1],self._eulerStart[2],self._eulerInc[0],self._eulerInc[1],self._eulerInc[2],self._eulerEnd[0],self._eulerEnd[1],self._eulerEnd[2])
        ang = [0,0,0]
        counter =0
        
        while not ang == [None,None,None]:
            ang = otherObj.nextRotation()
            if not ang == [None,None,None]:
                counter = counter +1
            
        return counter
    
class AngleError(Exception):    
    def __init__(self,value):
        self.value = value
    def __str__(self):
        print(self.value)
