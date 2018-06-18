'''
Created on Dec 11, 2009

@author: hrabe
'''
from pytom.basic.structures import PyTomClass

class Quaternion(PyTomClass):
    """
    Quaternion: Quaternion represent rotations in SE(3). Used for rotation distance metrics and so forth. Notation according
    'Effective Sampling and Distance Metrics for 3D Rigid Body Path Planning , Kuffner J. 2004'
    """
    
    def __init__(self,x,y,z,theta=0,isEulerRotation=True):
        """
        __init__ : Initialize - see isEulerRotation for how2
        @param x: x coordinate of rotation vector OR z1 rotation in Euler ZXZ paradigm in degrees
        @param y: y coordinate of rotation vector OR z2 rotation in Euler ZXZ paradigm in degrees
        @param z: z coordinate of rotation vector OR x rotation in Euler ZXZ paradigm in degrees
        @param theta: rotation around the vector (default 0)
        @param isEulerRotation: Default is True. Assumes that Quaternion object will be initialized from Euler rotation where x=z1,y=z2,z=x. Furthermore, all angles MUST be in degrees!      
        """
        from math import sin,cos
        
        if isEulerRotation:
            #convert to unit sphere rotation first
            from pytom.angles.angleFnc import zxzToAxisAngle
            [theta,[x,y,z]] = zxzToAxisAngle(x,y,z)
            
        from pytom.angles.angle import deg2rad
        
        self._w = cos(theta*deg2rad/2)
        self._x = x*sin(theta*deg2rad/2)
        self._y = y*sin(theta*deg2rad/2)
        self._z = z*sin(theta*deg2rad/2)
        
    def getW(self):
        return self._w
    
    def getX(self):
        return self._x
    
    def getY(self):
        return self._y
    
    def getZ(self):
        return self._z
    
    def setW(self,value):
        self._w = value
        
    def setX(self,value):    
        self._x = value
        
    def setY(self,value):    
        self._y = value
        
    def setZ(self,value):    
        self._z = value
        
    def distance2(self,otherQuaternion):
        """
        distance2: Calculates distance between this and another Quaternion according to algorithm 5 in source paper
        @param otherQuaternion: The other quaternion  
        """
        from pytom.angles.angle import rad2deg
        from math import acos
        
        lamb = self._w * otherQuaternion.getW()
        lamb = lamb + self._x*otherQuaternion.getX()
        lamb = lamb + self._y*otherQuaternion.getY()
        lamb = lamb + self._z*otherQuaternion.getZ()
    
        if lamb < 0:
            lamb = self._w * (-1 * otherQuaternion.getW())
            lamb = lamb + self._x*(-1 * otherQuaternion.getX())
            lamb = lamb + self._y*(-1 * otherQuaternion.getY())
            lamb = lamb + self._z*(-1 * otherQuaternion.getZ())
            
        return abs(1-lamb)
    
    
    def flip(self):
        """
        flip: Inverts the quaternion values. Returns self as a inversed quaternion.
        """
        q = Quaternion(0,0,0,0)
        q.setW(-self.getW())
        q.setX(-self.getX())
        q.setY(-self.getY())
        q.setZ(-self.getZ())
        
        return q
        
    def distance(self,otherQuaternion):
        """
        distance: Determines distance according to dotProduct and arcus cosinus
        @param otherQuaternion: The other quaternion
        @rtype: Angle in degrees
        @author: Luis Kuhn
        """
        from math import acos
        from pytom.angles.angle import rad2deg
        from pytom.tools.maths import epsilon
        
        if self == otherQuaternion:
            return 0.0
        
        d1 = self.dotProduct(otherQuaternion)
        
        if d1 >= 1.0+epsilon:
            d1 = 1.0
        elif d1 <= -1.0-epsilon:
            d1 = -1.0
            
        try: 
            d1 = acos(d1) * rad2deg * 2
        except:
            print ('Warning: Quarternion distance - acos can not be determined for ' + str(d1) + ' returning 0 instead!')
            print 'Quaternion values: '
            print self
            print otherQuaternion
            return 0.0
        
        d2 = self.dotProduct(otherQuaternion.flip())
        if d2 >= 1.0+epsilon:
            d2 = 1.0
        elif d2 <= -1.0-epsilon:
            d2 = -1.0
            
        try:
            d2 = acos(d2) * rad2deg * 2
        except:
            print 'Warning: Quarternion distance - acos can not be determined for ' + str(d2) + ' returning 0 instead!'
            print 'Quaternion values: '
            print self
            print otherQuaternion
            return 0.0
        
        dist = d1
        if d2 < d1:
            dist = d2
            
        return dist
        
    
    def toXML(self):
        """
        toXML: Converts this object to XML
        @return: XML object
        @author: Thomas Hrabe
        """
        from lxml import etree
        
        quaternion_element = etree.Element('Quaternion',x = str(self._x),y = str(self._y),z = str(self._z),w = str(self._w))
        
        return quaternion_element
        
    def fromXML(self,xmlObj=None):
        """
        fromXML: Initialises this object through the provided XML object
        @param xmlObj: The XML object. An exception is thrown if xmlObj is invalid.
        @author: Thomas Hrabe 
        """
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element :
            raise Exception('Is not a lxml.etree._Element! You must provide a valid XML-Quaternion object.')
        
        quaternion_element = xmlObj
        
        self._x = float(quaternion_element.get('x'))
        self._y = float(quaternion_element.get('y'))
        self._z = float(quaternion_element.get('z'))
        self._w = float(quaternion_element.get('w'))
    
    def toRotation(self):
        """
        toRotation:
        """
        from pytom.angles.angleFnc import psiThetaFromCenterVector
        
        
    def dotProduct(self,otherQuaternion):
        """
        dotProduct: performs the dot product (inner product) of two quaternions. Everything greater than 1 and lower than -1 will be set to 1 or -1.
        @param otherQuaternion: An other quaternion 
        @return: Value between -1 and 1
        """
        assert otherQuaternion.__class__ == Quaternion

        from pytom.tools.maths import epsilon
    
        sc  = self._w * otherQuaternion.getW()
        sc += self._x * otherQuaternion.getX()
        sc += self._y * otherQuaternion.getY()
        sc += self._z * otherQuaternion.getZ()
        
        return sc
        
        
        
        
    