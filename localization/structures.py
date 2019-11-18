'''
Created on May 21, 2010

@author: chen
'''

from pytom.basic.structures import PyTomClass


class Volume(PyTomClass):
    """
    Volume: Class for target & result volume
    """
    def __init__(self, filename='', subregion=[0,0,0,0,0,0], sampling=[0,0,0], binning=[0,0,0]):
        self._filename = filename
        self.subregion = subregion
        self.sampling = sampling
        self.binning = binning
        
    def getFilename(self):
        """
        getFilename: Get filename of the volume
        """
        return self._filename
    
    def setFilename(self,filename):
        """
        setFilename: Set the filename of the volume
        """
        self._filename = filename
    
    def getVolume(self, subregion=[0,0,0,0,0,0], sampling=[0,0,0], binning=[0,0,0]):
        """
        getVolume: Get the volume according to the given parameters.\
        Note the subregion, sampling and binning parameters are not the same as \
	those attributes in the class.

        @param subregion: subregion [startX, startY, startZ, sizeX, sizeY, sizeZ]
        @type subregion: list
        @param sampling: sampling [factorX, factorY, factorZ]
        @type sampling: list
        @param binning: binning [factorX, factorY, factorZ]
        @type binning: list
        
        @rtype: L{pytom_volume.vol}
        """
        from pytom_volume import read
        from pytom.tools.files import checkFileExists
        
        if not checkFileExists(self._filename):
            raise Exception('File ' + self._filename + ' does not exist!')
        
        return read(self._filename, subregion[0], subregion[1], subregion[2],
                    subregion[3], subregion[4], subregion[5],
                    sampling[0], sampling[1], sampling[2],
                    binning[0], binning[1], binning[2])
     
    def toXML(self): 
        """
        toXML : Compiles a XML file from object
        rtype : L{lxml.etree._Element}
        @author: chen
        """ 
        
        from lxml import etree
        e = etree.Element('Volume',Filename=self._filename, Subregion=str(self.subregion),
                          Sampling=str(self.sampling), Binning=str(self.binning))
        
        return e
    
    def fromXML(self,xmlObj):
        """
        @param xmlObj: A xml object  
        @type xmlObj: L{lxml.etree._Element}
        @author: chen
        """
        
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element :
            from pytom.basic.exceptions import ParameterError
            raise ParameterError('Is not a lxml.etree._Element! You must provide a valid XML-Volume object.')
    
        if xmlObj.tag == 'Volume':
            e = xmlObj
        else:
            e = xmlObj.xpath('Volume')
            e = e[0]
    
        self._filename = e.get('Filename')
        
        s = e.get('Subregion')
        if s != None and s!= 'None':
            s = s[1:-1]
            print(s)
            self.subregion = [int(float(i)) for i in s.split(',')]

        s = e.get('Sampling')
        if s != None and s!= 'None':
            s = s[1:-1]        
            self.sampling = [int(i) for i in s.split(',')]

        s = e.get('Binning')
        if s != None and s!= 'None':
            s = s[1:-1]
            self.binning = [int(i) for i in s.split(',')]
    
    def copy(self):
        """
        copy: Copies self to new Particle object
        """
        p = Volume()
        pXML = self.toXML()
        
        p.fromXML(pXML)
        return p


class Orientation(Volume):
    """
    Orientation: class for the result orientation
    """
    
    def toXML(self): 
        """
        toXML : Compiles a XML file from object
        rtype : L{lxml.etree._Element}
        @author: chen
        """ 
        
        from lxml import etree
        e = etree.Element('Orientation',Filename = self._filename)
        return e
    
    def fromXML(self,xmlObj):
        """
        @param xmlObj: A xml object  
        @type xmlObj: L{lxml.etree._Element}
        @author: chen
        """
        
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element :
            raise TypeError('Is not a lxml.etree._Element! You must provide a valid XML-Volume object.')
    
        if xmlObj.tag == 'Orientation':
            e = xmlObj
        else:
            e = xmlObj.xpath('Orientation')
            e = e[0]
    
        self._filename = e.get('Filename')
        

class FoundParticle(PyTomClass):
    """
    FoundParticle: Class for getting the particle according to the result score volume
    """
    def __init__(self, pos=None, orient=None, score=None, filename=''):
        """
        @param pos: pick position of the particle
        @type pos: L{pytom.basic.structures.PickPosition}
        @param orient: orientation of the pick particle
        @type orient: L{pytom.basic.structures.Rotation}
        @param score: score vaule of the picked particle
        @type score: L{pytom.score.score.Score}
        @param filename: filename of the found particle, if any
        @type filename: string
        """
        self.pos = pos
        self.orient = orient
        self.score = score
        self.filename = filename
        
    def getVolume(self):
        """
        getVolume: Get volume of the found particle
        
        @rtype: L{pytom_volume.vol}
        """
        from pytom_volume import read
        from pytom.tools.files import checkFileExists
        
        if not checkFileExists(self._filename):
            raise Exception('File ' + self._filename + ' does not exist!')
        
        return read(self.filename)
    
    def getPos(self):
        return [ int(i) for i in self.pos.toVector()]
    
    def toXML(self):
        from lxml import etree
        
        particle_element = etree.Element('FoundParticle', Filename = self.filename)
        
        particle_element.append(self.pos.toXML())
            
        particle_element.append(self.orient.toXML())
        
        particle_element.append(self.score.toXML())
        
#        class_element = etree.Element('Filename', Name = str(self.filename))
#        particle_element.append(class_element)
        
        return particle_element
    
    def fromXML(self, xmlObj):
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element:
            from pytom.basic.exceptions import ParameterError
            raise ParameterError('You must provide a valid XML object.')
        
        if xmlObj.tag == "FoundParticle":
            main = xmlObj
        else:  
            main = xmlObj.xpath('FoundParticle')
            if len(main) == 0:
                from pytom.basic.exceptions import PyTomClassError
                raise PyTomClassError("This XML is not a FoundParticle.")
            main = main[0]
            
        from pytom.basic.structures import PickPosition, Rotation
        from pytom.score.score import fromXML as fromXMLScore
        
        self.filename = main.get('Filename')
        
        p = main.xpath('PickPosition')[0]
        self.pos = PickPosition()
        self.pos.fromXML(p)
        
        o = main.xpath('Rotation')[0]
        self.orient = Rotation()
        self.orient.fromXML(o)
        
        s = main.xpath('Score')[0]
        self.score = fromXMLScore(s)
        
        
    def toParticle(self):
        """
        toParticle: Converts this object to a Particle object.
        @rtype: L{pytom.basic.structures.Particle}
        """
        from pytom.basic.structures import Particle
        return Particle(filename=self.filename, rotation=self.orient, pickPosition=self.pos, score=self.score)
    
    def fromParticle(self, p):
        self.filename = p.getFilename()
        self.orient = p.getRotation()
        self.pos = p.getPickPosition()
        self.score = p.getScore()


class ClassifiedParticle(PyTomClass):
    """
    ClassifiedParticle: Class for storing the info about classified result particle
    """
    def __init__(self, particle=None, classified=None, correct=None):
        """
        @param particle: classified particle
        @type particle: L{pytom.localization.structures.FoundParticle or simulation.SimulatedParticle or simulation.IdentifiedParticle}
        @param classified: classified result
        @type classified: string
        @param correct: correct classified or not
        @type correct: boolean
        """
        self.particle=particle
        self.classified=classified
        self.correct=correct
    
    def getPos(self):
        if not self.particle:
            return None
        
        from pytom.localization.simulation import SimulatedParticle, IdentifiedParticle
        if self.particle.__class__ == SimulatedParticle:
            return self.particle.pos
        if self.particle.__class__ == IdentifiedParticle:
            return self.particle.f_particle.pos.toVector()
        if self.particle.__class__ == FoundParticle:
            return self.particle.pos.toVector()
    
    def getScore(self):
        if not self.particle:
            return None
        
        from pytom.localization.simulation import IdentifiedParticle
        if self.particle.__class__ == IdentifiedParticle:
            return float(self.particle.f_particle.score.scoreValue)
        if self.particle.__class__ == FoundParticle:
            return float(self.particle.score.scoreValue)
    
    def getRot(self):
        pass
    
    def toXML(self):
        from lxml import etree
        
        particle_element = etree.Element('ClassifiedParticle', Classified=str(self.classified), Correct=str(self.correct))
        
        particle_element.append(self.particle.toXML())
        
        return particle_element
    
    def fromXML(self, xmlObj):
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element:
            from pytom.basic.exceptions import ParameterError
            raise ParameterError('You must provide a valid XML object.')
        
        if xmlObj.tag == "ClassifiedParticle":
            main = xmlObj
        else:  
            main = xmlObj.xpath('ClassifiedParticle')
            if len(main) == 0:
                from pytom.basic.exceptions import PyTomClassError
                raise PyTomClassError("This XML is not a IdentifiedParticle.")
            main = main[0]
            
        self.classified = main.get('Classified')
        
        correct = main.get('Correct')
        if correct == 'True':
            self.correct = True
        elif correct == 'False':
            self.correct = False
        else:
            self.correct = None
        
        f_particle = main.xpath('FoundParticle')
        if len(f_particle) > 0:
            self.particle = FoundParticle()
            self.particle.fromXML(f_particle[0])
        
        s_particle = main.xpath('SimulatedParticle')
        if len(s_particle) > 0:
            from pytom.localization.simulation import SimulatedParticle
            self.particle = SimulatedParticle()
            self.particle.fromXML(s_particle[0])
        
        i_particle = main.xpath('IdentifiedParticle')
        if len(i_particle) > 0:
            from pytom.localization.simulation import IdentifiedParticle
            self.particle = IdentifiedParticle()
            self.particle.fromXML(i_particle[0])


class ParticleList(PyTomClass):
    """
    ParticleList: Class for storing all kinds of particles
    """
    def __init__(self):
        self.pl = []
    
    def append(self, p):
        self.pl.append(p)
    
    def toXML(self):
        from lxml import etree
        root = etree.Element('ParticleList')
        for p in self.pl:
            root.append(p.toXML())
            
        return root
    
    def fromXML(self, xmlObj):
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element:
            from pytom.basic.exceptions import ParameterError
            raise ParameterError('You must provide a valid XML object.')
        
        if xmlObj.tag == "ParticleList":
            main = xmlObj
        else:  
            main = xmlObj.xpath('ParticleList')
            if len(main) == 0:
                from pytom.basic.exceptions import PyTomClassError
                raise PyTomClassError("This XML is not a ParticleList.")
            main = main[0]
            
        for p in main.xpath('FoundParticle'):
            particle = FoundParticle()
            particle.fromXML(p)
            self.pl.append(particle)
        
        for p in main.xpath('ClassifiedParticle'):
            particle = ClassifiedParticle()
            particle.fromXML(p)
            self.pl.append(particle)
            
        for p in main.xpath('SimulatedParticle'):
            from pytom.localization.simulation import SimulatedParticle
            particle = SimulatedParticle()
            particle.fromXML(p)
            self.pl.append(particle)
            
        for p in main.xpath('IdentifiedParticle'):
            from pytom.localization.simulation import IdentifiedParticle
            particle = IdentifiedParticle()
            particle.fromXML(p)
            self.pl.append(particle)
        
        for p in main.xpath('Particle'): # now also support reading class Particle, but will convert to FoundParticle
            from pytom.basic.structures import Particle
            pp = Particle()
            pp.fromXML(p)
            
            # convert to FoundParticle
            particle = FoundParticle()
            particle.fromParticle(pp)
            self.pl.append(particle)


def readParticleFile(filename):
    """
    readParticleFile: Read the particle file from the disk
    @param filename: Filename of particles
    @type filename: string
    
    @return: List of particles
    @rtype: list
    """
    pl = ParticleList()
    pl.fromXMLFile(filename)
    
    return pl.pl

