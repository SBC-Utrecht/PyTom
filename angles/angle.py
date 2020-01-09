import math
rad2deg = 180.0/math.pi
deg2rad = math.pi/180.0

from pytom.basic.structures import PyTomClass

def fromStr(string):
        """
        fromStr: Creates any AngleObject child depending on the XML string specified
        @param string:
        @author: Thomas Hrabe 
        """
        from lxml import etree

        root = etree.fromstring(string)

        return fromXML(root)
    
def fromXML(xmlObj):
    """
    fromXML
    @deprecated: Use AngleObject.fromXML instead!!!
    """
    print('This function is deprecated. Use the method of AngleObject.fromXML instead!')
    from pytom.angles.angle import AngleObject
    ang = AngleObject()
    return ang.fromXML(xmlObj)
    

class AngleObject(PyTomClass):
    """
    AngleObject: Template class for all angle objects. Methods defined here must be overloaded by others. The angles stored are in degrees ( pytom sublayer converts them to radiants autmatically). 
    G{classtree}
    """    
    
    def __init__(self):
        pass

    def toXML(self):
        """
        overload toXML function, just for the sake of it
        """
        raise Exception("AngleObject: Function toXML should not be called, instead overloaded")
    
    def fromXML(self,xmlObj, verbose=False):
        """
        fromXML: Creates any AngleObject child depending xmlObj
        @param xmlObj: a XML DOM object
        @param verbose: verbose mode?
        @type verbose: L{bool}
        @return: an Angle List as defined by xmlObj
        @author: Thomas Hrabe
        G{callgraph}
        """
        angType = xmlObj.get('Type')

        if angType == 'Eulerian':
            from pytom.angles.angleList import EulerAngleList
            ang = EulerAngleList()
            ang.fromXML(xmlObj)
        elif angType == 'Equidistant':
            from pytom.angles.localSampling import LocalSampling
            ang = LocalSampling()
            ang.fromXML(xmlObj)
        elif angType == 'AV3Sampling':
            from pytom.angles.localSampling import AV3Sampling
            ang = AV3Sampling()
            ang.fromXML(xmlObj)
        elif angType == 'ExtendedInplaneSampling':
            from pytom.angles.localSampling import ExtendedInplaneSampling
            ang = ExtendedInplaneSampling()
            ang.fromXML(xmlObj)
        elif angType == 'GlobalSampling':
            from pytom.angles.globalSampling import GlobalSampling
            ang = GlobalSampling()
            ang.fromXML(xmlObj)
        elif angType == 'FromEMFile':
            from pytom.angles.globalSampling import GlobalSampling
            ang = GlobalSampling()
            ang.fromXML(xmlObj)
        elif angType == 'AngleList':
            from pytom.angles.angleList import AngleList
            ang = AngleList()
            ang.fromXML(xmlObj)
        elif angType == 'OneAngleList':
            from pytom.angles.angleList import OneAngleList
            ang = OneAngleList()
            ang.fromXML(xmlObj)
        elif angType == 'LocalSampling':
            from pytom.angles.localSampling import LocalSampling
            ang = LocalSampling()
            ang.fromXML(xmlObj)
        #elif angType == 'RestrictedInplaneLocalSampling':
        #    from pytom.angles.localSampling import  InplaneLocalSampling
        #    ang = RestrictedInplaneEuidistantList()
        #    ang.fromXML(xmlObj)
        elif angType == 'EquidistantList':
            from pytom.angles.localSampling import LocalSampling
            ang = LocalSampling()
            ang.fromXML(xmlObj)
        #elif angType == 'RestrictedInplaneEuidistantList':
        #    from pytom.angles.localSampling import RestrictedInplaneLocalSampling
        #    ang = RestrictedInplaneEuidistantList()
        #    ang.fromXML(xmlObj)
        elif angType == 'Combined':
            from pytom.angles.combined import GlobalLocalCombined
            ang = GlobalLocalCombined()
            ang.fromXML(xmlObj)
        else:
            raise TypeError('Type ' +angType+' not available in Angles.')

        if verbose:
            print("AngleObject.fromXML: Returned AngleObject: "+str(ang))

        self.ang = ang

        return ang


    def fromXMLFile(self, filename):
        """
        fromXMLFile: Configures object according to file. Do NOT overload this function in any child!
        @param filename: Absolute / relative path to file
        @type filename: L{str}
        @author: Thomas Hrabe
        """
        from pytom.angles.globalSampling import GlobalSampling
        from pytom.tools.files import readStringFile
        from lxml import etree

        self._filename = filename
        lines = readStringFile(filename)
        root = etree.fromstring(lines)
        angObject = self.fromXML(root)
        
        if not type(angObject) == GlobalSampling:
            angObject._filename = filename

        return angObject

    
    def setRotationCenter(self,x,y,z):
        """
        setRotationCenter: Set rotation center in pixels
        @param x: Rotation center in x direction
        @param y: Rotation center in y direction
        @param z: Rotation center in z direction
        """
        self._rotationCenterX = x
        self._rotationCenterY = y
        self._rotationCenterZ = z
        
    def getRotationCenter(self):
        """
        getRotationCenter: Get rotation center in pixels
        @return: Returns rotation center as [x,y,z] if values are set or [0,0,0]
        """
        if hasattr(self,'_rotationCenterX'):
            return [self._rotationCenterX,self._rotationCenterY,self._rotationCenterZ]
        else:
            return [0,0,0]
        
        
    def setStartRotation(self,startRotation=None):
        """
        setStartRotation: Sets start rotation.
        @param startRotation: The start rotation. 
        @return: Reference to current object
        @rtype: L{pytom.angles.angle.AngleObject}
        @author: Thomas Hrabe 
        """
        print('This is an abstract method, You must override this method in your subclass')
        
        assert False
        
    def oldRotationsToXML(self):
        """
        oldRotationsToXML : compiles XML object from oldRotationList so that they will be scanned in the next round,too.
        """
        from lxml import etree
        oldElement = etree.Element("OldRotations")
        
        for i in range(len(self._oldRotations)):
            rotation = self._oldRotations[i]
            rotElement = etree.Element('Rotation',phi=str(rotation[0]),psi=str(rotation[1]),theta=str(rotation[2])) 
            oldElement.append(rotElement)
            
        return oldElement
    
    def oldRotationsFromXML(self,xmlObj):
        """
        oldRotationsFromXML: Creates self.oldRotations list from xmlObj
        @param xmlObj: A XML attribute as generated by toXML 
        @author: Thomas Hrabe
        """
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element:
            raise TypeError('You must provide a valid XML-EulerAngleList object.')
        
        oldElement = xmlObj.xpath('OldRotations')
        
        if len(oldElement) == 0:
            oldElement = xmlObj.xpath('/*/OldRotations')
        
        oldElement = oldElement[0]

        for rot in oldElement.iter('Rotation'):
            phi = float(rot.get('phi'))
            psi = float(rot.get('psi'))
            the = float(rot.get('theta'))
            self._oldRotations.append([phi,psi,the])
    
    def getOldRotations(self):
        return self._oldRotations
    
    def nextRotation(self):
        """
        nextRotation: Returns the next rotation for the current object
        @author: Thomas Hrabe 
        """
        print('This is a virtual class and thus virtual member function')
        print(self.__class__)
        assert False
    
    def getIncrement(self):
        """
        getIncrement:
        @return: Current angular increment
        """
        print('This is a virtual class and thus virtual member function')
        print(self.__class__)
        assert False
        
    def focusRotation(self,rotation,refinementAngle=None): 
        """
        focusRotation: Returns new object which focuses on the rotation provided
        @param rotation: The rotation to focus on
        @return: New angle object
        @author: Thomas Hrabe 
        """
        print('This is a virtual class and thus virtual member function')
        print(self.__class__)
        assert False
        
    def numberRotations(self):
        """
        numberRotations: Returns the total number of rotations for the current object.
        @author: Thomas Hrabe 
        @deprecated: Use __len__ instead
        """
        print('This is a virtual class and thus virtual member function')
        assert False

    def __len__(self):
        return self.numberRotations()
    
    def printRotations(self):
        
        rotations = self.getAllRotations()
        
        for i in range(len(rotations)):
            print(rotations[i])
        
    def getAllRotations(self):
        """
        getAllRotations: Returns a list of all rotations this object will return
        @return: list of all rotations 
        @rtype: L{list}
        """
        
        rotations = []
        ang = [0,0,0]
        
        while not ang == [None,None,None]:
            ang = self.nextRotation()
            rotations.append(ang)
            
        self.reset()
        
        return rotations
        
    def reset(self):
        """
        reset: Resets the object that nextRotation would return the same same sequence of results again   
        @author: Thomas Hrabe 
        """
        print('This is a virtual class and thus virtual member function')
        print(self.__class__)
        assert False
    
    def __getitem__(self,key):
        """
        __getitem__: Allow array style and slice access to rotations
        @return: A rotation
        @rtype: L{pytom.basic.structures.Rotation}
        """
                
        if isinstance(key, int): 
            
            index = 0
            rotation = self.nextRotation()
            
            while not index == key and not rotation == [None,None,None]:
                index = index + 1
                rotation = self.nextRotation()
            
            self.reset()
            
            return rotation
        else:
            assert False
        
        
    def copy(self,rotation=None):
        """
        copy: Copies the current angle object to a new one with same settings but different starting rotation
        @param rotation: The other starting rotation   
        @author: Thomas Hrabe 
        """
        print(self.__class__)
        print('This is a virtual class and thus virtual member function')
        assert False
        
    """
    def distanceFunction(self,rotation):
        
        distanceFunction: determines the distance of of given rotation to old rotation
        @param rotation:    
        @author: Thomas Hrabe 
        
        print 'This is a virtual class and thus virtual member function'
        assert False
    """
    
    def toEMFile(self,filename, inputIsRadians=False):
        """
        toEMFile: Will save all rotations to EM file. The EM file can be read out by L{pytom.angles.fromFile.AngleListFromEM}
        @param filename: Name of file to save to 
        @param inputIsRadians: 
        """
        from pytom_volume import vol
        
        no = self.numberRotations()
         
        anglesVolume = vol(3,no,1) 
        
        rotations = self.getAllRotations()
        
        i = 0
        for rotation in rotations:
            if rotation == [None, None, None]:
                break
            if inputIsRadians:
                anglesVolume.setV(rotation[0],0,i,0)
                anglesVolume.setV(rotation[1],1,i,0)
                anglesVolume.setV(rotation[2],2,i,0)
            else:
                from pytom.angles.angle import deg2rad
                anglesVolume.setV(rotation[0]*deg2rad,0,i,0)
                anglesVolume.setV(rotation[1]*deg2rad,1,i,0)
                anglesVolume.setV(rotation[2]*deg2rad,2,i,0)               
            
            i = i+1
            
        anglesVolume.write(filename)
    
    def showRotations(self,volume,directory,mask=None):
        """
        showRotations: Will store x,y,z slice through volume after each rotation was applied 
        """
        from pytom.tools.toImage import volumeToPNG
        from pytom.tools.files import dump2TextFile
        from pytom_volume import vol,rotateSpline
        from pytom.tools.files import getPytomPath,readStringFile
        from pytom.frontend.htmlDefines import headResponse200,copyright
    
        pytomPath = getPytomPath()
        cssDefines = readStringFile(pytomPath + '/frontend/htmlDefines/CSS.css')
    
        html = '<html>' + cssDefines
        html = html + '\n<title>PyTom Alignment</title>\n'
        html = html + '<body><center><h1 type="Main">List of Rotations</h1></center>\n'
        html = html + '</br></br></br></br></br></br></br></br><center>\n'
        html = html + '<table>\n'
        html = html + '<tr><td>Z1 (Phi) Z2 (Psi) X(Theta)</td><td>x-Slice</td><td>y-Slice</td><td>z-Slice</td></tr>'
        
        rotated = vol(volume.sizeX(),volume.sizeY(),volume.sizeZ())
        
        if not mask:
            mask = vol(volume.sizeX(),volume.sizeY(),volume.sizeZ())
            mask.setAll(1)
            
        volume = volume * mask    
        
        self.reset()
        
        rotation = self.nextRotation()
        rotationCounter = 0
        
        while not rotation == [None,None,None]:
        
            rotateSpline(volume,rotated,rotation[0],rotation[1],rotation[2])
            html = html + '<tr><td>' + str(rotation[0]) + ' ' + str(rotation[1]) + ' ' + str(rotation[2]) +'</td>/n'
            volumeToPNG(rotated,directory + str(rotationCounter) + '-x.png',rotated.sizeX()/2,'x')
            html = html + '<td><img src="' + directory + str(rotationCounter) + '-x.png"/></td>'
            volumeToPNG(rotated,directory + str(rotationCounter) + '-y.png',rotated.sizeY()/2,'y')
            html = html + '<td><img src="' + directory + str(rotationCounter) + '-y.png"/></td>'
            volumeToPNG(rotated,directory + str(rotationCounter) + '-z.png',rotated.sizeZ()/2,'z')
            html = html + '<td><img src="' + directory + str(rotationCounter) + '-z.png"/></td>'
            
            rotation = self.nextRotation()
            rotationCounter = rotationCounter + 1 
        
        html = html + '</body></html>'
        
        dump2TextFile(directory + 'rotations.html',html)
        
    def __eq__(self,other):    
        return (str(self) == str(other))    
    
    
    def rotationDistanceMatrix(self):
        """
        rotationDistanceMatrix: Generate a matrix of rotation distances. The 
        distance is determined by L{pytom.angles.quaternions.Quaternion.distance}.
        @return: L{pytom_volume.vol} of rotation distances
        @author: Thomas Hrabe
        """
        from pytom_volume import vol
        from pytom.tools.maths import rotation_distance
        numberOfRotations = len(self)

        distanceMatrix = vol(numberOfRotations,numberOfRotations,1)
        distanceMatrix.setAll(0)

        for i in range(len(self)):
            if i < len(self):
                for j in range(i+1,len(self)):
                    ri = self[i]
                    rj = self[j]

                    d = rotation_distance(ri,rj)

                    distanceMatrix(d,i,j,0)
                    distanceMatrix(d,j,i,0)

            distanceMatrix(0,i,i,0)

        return distanceMatrix
    
    
