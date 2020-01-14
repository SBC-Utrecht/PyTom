'''
Created on May 2, 2013

@author: thrabe
'''
from pytom.angles.angleList import AngleList
from pytom.angles.angle import AngleObject

class LocalSampling(AngleObject):
    """
    LocalSampling: Angular sampling identical to the AV3 package. \
    See http://pytom.org/doc/pytom/alignment.html for more information.
    """
        
    def __init__(self,shells=3,increment=3,z1Start=0.0,z2Start=0.0,xStart=0.0):
        """
        @param shells: Number of shells to be scanned
        @param increment: Angular increment used
        @param z1Start: Start angle for Z1 rotation 
        @param z2Start: Start angle for Z2 rotation 
        @param xStart: Start angle for X rotation 
        @author: Thomas Hrabe     
        """
        from pytom.basic.structures import Rotation
        self._shells = float(shells)

        if increment == 0.0:
            raise ValueError('LocalSampling : Increment is 0!')
        else:
            self._increment = float(increment)
            
        if shells == 0:
            raise ValueError('LocalSampling : Shells is 0!')
        
        self.setStartRotation(Rotation(z1=z1Start, z2=z2Start, x=xStart))
        # initialize final rotation around z-axis of REFERENCE
        self.reset()

        
    def setStartRotation(self,startRotation):
        """
        setStartRotation: Sets start rotation.
        @param startRotation: The start rotation.
        @type startRotation: L{pytom.basic.structures.Rotation}
        @return: Reference to current object
        @author: FF
        """
        from pytom.basic.structures import Rotation

        self._startZ1 = float(startRotation[0])
        self._startZ2 = float(startRotation[1])
        self._startX  = float(startRotation[2])
        
        self._startRotation = Rotation( z1=self._startZ1, z2=self._startZ2, x=self._startX, paradigm='ZXZ')

        self._startMatrix = self._startRotation.toMatrix()
        self.reset()
    
        return self
    
    def getNumberShells(self):
        """
        get number of shells for alignment

        @return: number of shells
        @rtype: L{int}
        """
        return self._shells

    def setNumberShells(self, shells):
        """
        set number of shells
        @param shells: number of shells
        @type shells: L{int}
        """
        self._shells = shells
    
    def getIncrement(self):
        """
        get angular increment
        @return: angular increment
        @rtype: L{float}
        """
        return self._increment

    def setIncrement(self, increment):
        """
        set angular increment
        @param increment: new increment in deg
        @type increment: L{float}
        @author: FF
        """
        self._increment = increment
    
    def nextRotation(self):
        """
        nextRotation : 
        @return: [z1 z2 x] for the next rotation or [None,None,None] after all rotations were sampled
        @author: Friedrich Foerster
        @change: Local Rotation had a bug causing too large rotations in Phi
        @date: 07/07/2014
        """
        if self._finished:
            return [None,None,None]
        
        from math import sin,ceil,pi,sqrt,atan2#,modf
        from pytom.basic.structures import Rotation
        from pytom.angles.angleFnc import matToZXZ
        
        phi = self._currentZ1
        if self._currentX == 0:
            npsi=1
            dpsi=360.
        else:
            dpsi = self._increment /sin(float(self._currentX*self._increment)/180.*pi)
            npsi = ceil(360./dpsi)
            #make dpsi equidistant again
            dpsi = 360./npsi

        localRotation = Rotation( z1=phi-self._currentZ2*dpsi, z2=self._currentZ2*dpsi,
                                  x=self._currentX*self._increment, paradigm='ZXZ')

        globalMatrix = localRotation.toMatrix() * self._startMatrix
        
        [phi,psi,theta] = matToZXZ(globalMatrix)

        if self._currentZ2 >= npsi-1:
            self._currentZ2 =0
            if self._currentX >= ceil(self._shells/2):
                self._currentX = 0
                if self._currentZ1 >= self._shells*self._increment:
                    self._finished = True
                    return [self._startZ1,self._startZ2,self._startX]
                else:
                    self._currentZ1 = self._currentZ1 + self._increment
            else:
                self._currentX = self._currentX +1
        else:
            self._currentZ2 = self._currentZ2 + 1

        return [phi%360,psi%360,theta%360]
    
    
    def focusRotation(self,rotation=None,refinementAngle=10):
        """
        focusRotation: Will focus this object on a given rotation by tightening the scan area. 
        @param rotation: The rotation to focus on
        @type rotation: [phi,psi,theta] list 
        @return: New LocaLSampling object
        @rtype: L{pytom.angles.localSampling.LocalSampling}
        @author: Thomas Hrabe 
        """
        if not rotation:
            rotation = [0,0,0]
        
        if not refinementAngle:
            refinementAngle = 10
            
        return LocalSampling(self._shells,refinementAngle,rotation[0],rotation[1],rotation[2])
    
    def distanceFunction(self,rotation):
        """
        distanceFunction: determines the distance of of given rotation to old rotation
        @param rotation: rotation
        @type rotation: L{pytom.basic.structures.Rotation}
        @return: proposed increment for next iteration
        @author: FF
        """
        from pytom.tools.maths import rotation_distance
        increment =  rotation_distance(ang1=self._startRotation, ang2=rotation)
        
        #from math import sqrt,ceil
        #sqrt2 = 1.4142135623730951
        ##must modulo 360 each value for avoiding negative values
        #deltaTheta  = (self._startX % 360) - (rotation[2] % 360)
        #deltaPhi    = (self._startZ1 % 360)   - (rotation[0] % 360)
        #
        #increment   = sqrt(pow(deltaTheta,2)+pow(deltaPhi,2)) / sqrt2
        
        return increment
        
    def toXML(self):
        from lxml import etree

        try:
            angles_element = etree.Element("Angles",Type = 'AV3Sampling')
            angles_element.set("Increment", str(self._increment))
            angles_element.set("Shells", str(self._shells))
            angles_element.set("Phi_old", str(self._startZ1))
            angles_element.set("Psi_old", str(self._startZ2))
            angles_element.set("Theta_old", str(self._startX))
            angles_element.set("ShellsParameter", str(self._shellsParameter))
            angles_element.set("IncrementParameter", str(self._incrementParameter))
        except:
            angles_element = etree.Element("Angles", Type='LocalSampling')
            angles_element.set("Increment", str(self._increment))
            angles_element.set("Shells", str(self._shells))
            angles_element.set("StartZ1", str(self._startZ1))
            angles_element.set("StartZ2", str(self._startZ2))
            angles_element.set("StartX", str(self._startX))

        return angles_element
        
    def fromXML(self,xmlObj=-1):
        """
        fromXML : Assigns values to job attributes from XML object
        @param xmlObj: A xml object  
        @author: Thomas Hrabe 
        """
        from lxml.etree import _Element
        from pytom.basic.structures import Rotation
        
        if xmlObj.__class__ != _Element or xmlObj.get('Type') not in  ['LocalSampling','Equidistant','AV3Sampling']:
            raise TypeError('You must provide a valid XML-LocalSampling object.')
        
        if xmlObj.get('Type') == 'Equidistant':
            xmlObj = xmlObj.xpath('Parameters')[0]

        try:    
            self._increment = float(xmlObj.get('Increment'))
        except TypeError:
            self._increment = float(xmlObj.get('AngleIncrement'))
            
        try:
            self._shells  = float(xmlObj.get('Shells'))
        except TypeError:
            self._shells  = float(xmlObj.get('NumberShells'))
            
        try:
            self._startZ1    = float(xmlObj.get('StartZ1'))
            self.av3 = False
        except:
            self._startZ1    = float(xmlObj.get('Phi_old'))
            self.av3 = True
        try:
            self._startZ2    = float(xmlObj.get('StartZ2'))
        except:
            self._startZ2    = float(xmlObj.get('Psi_old'))
        try:
            self._startX    = float(xmlObj.get('StartX'))
        except:
            self._startX  = float(xmlObj.get('Theta_old'))

        if self.av3:
            try:
                self._shellsParameter = int(xmlObj.get('ShellsParameter'))
            except TypeError:
                raise Exception('No ShellsParameter Defined')
            try:
                self._incrementParameter = int(xmlObj.get('IncrementParameter'))
            except TypeError:
                raise Exception('No IncrementParameter Defined')

        self.reset()

        self.setStartRotation(startRotation=Rotation(z1=self._startZ1, z2=self._startZ2, x=self._startX))
    
    def numberRotations(self):
        """
        numberRotations: Returns the total number of rotations for the current object
        @author: Thomas Hrabe 
        """

        ang = [0,0,0]
        counter =0
        
        while not ang == [None,None,None]:
            ang = self.nextRotation()
            counter = counter +1
            
        self.reset()
        
        return counter


    def reset(self):
        """
        reset: Resets the object that nextRotation would return the same sequence of results again
        @author: FF
        """
        self._currentZ1 = -1.*self._shells*self._increment
        self._currentZ2 = 0
        self._currentX = 0
        self._finished = False
        #self._startRotation = None

    def copy(self,rotation):
        """
        copy: Copies the current angle object to a new one with same settings but different starting rotation
        @param rotation: The other starting rotation   
        @author: Thomas Hrabe 
        """
        if rotation.__class__ == list:
            return LocalSampling(self._shells,self._increment,rotation[0],rotation[1],rotation[2])
        else:
            return LocalSampling(self._shells,self._increment,rotation.getPhi(),rotation.getPsi(),rotation.getTheta())

class AV3Sampling(LocalSampling):
    """
    AV3Sampling: Deprecated and for backward compatibility
    @deprecated: 
    """ 
    
    def __init__(self,shells=3,increment=3,z1Start=0.0,z2Start=0.0,xStart=0.0):
        """
        call LocalSampling instead ...
        """
        super(self.__class__,self).__init__(shells,increment,z1Start,z2Start,xStart)


class ExtendedInplaneSampling(AngleList):    
    """
    ExtendedInplaneSampling: Samples in the proximity of a specific rotation but samples all inplane rotations. \
    See http://pytom.org/doc/pytom/alignment.html for more information.
    @author: Thomas Hrabe
    """

    def __init__(self,numberShells = 0, angleIncrement = 1, startZ1 = 0, startZ2 =0, 
            startX=0,shellIncrement = None):
        """
        @param numberShells: Number of shells around start rotation
        @param angleIncrement: Angular increment for inplane rotation
        @param startZ1: Z1 position of start rotation 
        @param startZ2: Z2 position of start rotation 
        @param startX: X position of start rotation
        @param shellIncrement: Distance of shell from start rotation or previous \
        shell. Will be set to angleIncrement if not specified. 
        """
        self._rotationList = []

        self._numberShells = int(numberShells)
        self._angleIncrement = float(angleIncrement)
        self._startZ1 = float(startZ1)
        self._startZ2 = float(startZ2)
        self._startX = float(startX)
        
        if shellIncrement:
            self._shellIncrement = shellIncrement
        else:
            self._shellIncrement = angleIncrement
        
    def _doAllInplaneRotations(self,z2,x):
        
        from pytom.basic.structures import Rotation
        from pytom.tools.macros import frange
        rotationList = []
        
        for z1 in frange(0,360,self._angleIncrement):
            rotationList.append(Rotation((self._startZ1 + z1)%360,z2,x))
        
        return rotationList
    
    def _doShellRotations(self,currentShell):
        from math import sin,cos
        from pytom.tools.maths import epsilon
        from pytom.angles.angle import deg2rad
        from pytom.tools.macros import frange
        rotationList = []
        
        
        for angle in frange(0,360,self._shellIncrement):
            psi = currentShell * sin(deg2rad * angle) * self._shellIncrement
            if abs(psi) <=epsilon:
                psi = 0.0
                
            theta = currentShell * cos(deg2rad * angle) * self._shellIncrement
            if abs(theta) <=epsilon:
                theta = 0.0
                
            rotationList = rotationList + self._doAllInplaneRotations((psi+ self._startZ2 ), (theta + self._startX) )
            
        return rotationList
    
    def _initRotationList(self):
        
        for currentShell in range(self._numberShells+1):
            
            if currentShell == 0.0:
                self._rotationList = self._rotationList + self._doAllInplaneRotations(self._startZ2, self._startX)
            else: 
                self._rotationList = self._rotationList + self._doShellRotations(currentShell)
            
    def getIncrement(self):
        return self._angleIncrement
    
    def getNumberShells(self):
        return self._numberShells
    
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
        
        if rotation is None:
            rotation = [0,0,0]
        
        if refinementAngle is None:
            refinementAngle = 10
        
        return LocalSampling(self._numberShells,refinementAngle,
            rotation[0],rotation[1],rotation[2])
        
    def setStartRotation(self,startRotation):
        """
        setStartRotation: Sets start rotation.
        @param startRotation: The start rotation.  
        @return: Reference to current object
        """
        
        return LocalSampling(self._numberShells, self._angleIncrement, startRotation[0], startRotation[1] , startRotation[2])
        
    def copy(self,rotation):
        from pytom.basic.structures import Rotation
        
        if rotation.__class__ == Rotation:
            return LocalSampling(self._numberShells,self._angleIncrement,rotation.getZ1(),rotation.getZ2(),rotation.getX())
        elif rotation.__class__ == list:
            return LocalSampling(self._numberShells,self._angleIncrement,rotation[0],rotation[1],rotation[2])
        else:
            raise RuntimeError('You must provide either a Rotation or a list [z1,x,z2]/[phi,psi,theta] to ExtendedInplaneSampling.focusRotation') 
    
    def toXML(self):
        from lxml import etree
        
        angles_element = etree.Element('Angles',Type='ExtendedInplaneSampling',NumberShells = str(self._numberShells),AngleIncrement = str(self._angleIncrement),StartZ1 = str(self._startZ1), StartZ2 = str(self._startZ2), StartX = str(self._startX),ShellIncrement = str(self._shellIncrement))
        
        return angles_element
    
    def fromXML(self,xmlObj):
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element :
            raise RuntimeError('Is not a lxml.etree._Element! You must provide a valid XMLobject.')
        
        if xmlObj.tag == 'Angles':
            angleListElement = xmlObj
        else:
            RuntimeError('You must provide a valid ExtendedInplaneSampling object.')
        
        self._rotationList = []

        self._numberShells = None
        self._angleIncrement = None
        self._startZ1 = None
        self._startX = None
        self._startZ2 = None
        
        self._numberShells = int(angleListElement.get('NumberShells'))
        self._angleIncrement = float(angleListElement.get('AngleIncrement'))
        
        
        if angleList_element.get('ShellIncrement') is None:
            self._shellIncrement = self._angleIncrement
        else:
            self._shellIncrement = float(angleListElement.get('ShellIncrement'))
        
        self._startZ1 = float(angleListElement.get('StartZ1'))
        self._startX = float(angleListElement.get('StartX'))
        self._startZ2 = float(angleListElement.get('StartZ2'))
        
        
    def nextRotation(self):
        if len(self._rotationList) == 0:
            self._initRotationList()
            
        return super(LocalSampling,self).nextRotation()    
    
    def numberRotations(self):
        if len(self._rotationList) == 0:
            self._initRotationList()
        return len(self._rotationList)
    
    
