from pytom.basic.structures import PyTomClass

class Preprocessing(PyTomClass):
    """
    Preprocessing: Defines procedures used for preprocessing before each alignment 
    step such as - bandpass filtering, - prerotation of volume, - rotation of weighting.
    The highest frequency in the bandpass filter will never be lower than 0.1, even 
    when specified lower. 
    Check code in apply function.
    """
    
    
    def __init__(self,lowestFrequency = None, highestFrequency = None, smooth=0,
                 prerotation=-1, weightingFile='', substractParticle=False, taper=0):
        """
        @param lowestFrequency: lowest frequency of bandpass filter
        @type lowestFrequency: int
        @param highestFrequency: highest frequency of bandpass filter
        @type highestFrequency: int
        @param smooth: smoothing width for bandpass
        @type smooth: float
        @param prerotation: ???
        @type prerotation: ???
        @param weightingFile: ???
        @type weightingFile: ???
        @param substractParticle: subtract particle from reference
        @type substractParticle: bool
        @param taper: taper width for edge tapering
        @type taper: float
        """
        self.setBandpass(lowestFrequency, highestFrequency, smooth)
        self.setRotation(prerotation)
        self.setWeighting(weightingFile)
        self.setSubstractParticle(substractParticle)
        self.setTaper(taper=taper)
        self._taperMask = None

    def setTaper(self, taper=0):
        """
        set width for apodizing
        @param taper: width of tapered edge
        @type taper: float
        """
        
        self._taper = taper
        
    def getTaper(self):
        """
        get width for apodizing
        
        @return: taper
        @rtype: float
        """
        
        return self._taper

    def getHighestFrequency(self): 
        return self._highestFrequency

    def getLowestFrequency(self): 
        return self._lowestFrequency

    def getBandpassSmooth(self):
        return self._bandpassSmooth
    
    def getSubstractParticle(self):
        return self._substractParticle
    
    def setHighestFrequency(self,highestFrequency):
        self._highestFrequency = highestFrequency
        
    def setLowestFrequency(self,lowestFrequency):
        self._lowestFrequency  = lowestFrequency
    
    def setSubstractParticle(self,on):
        self._substractParticle = on
        
    def apply(self, volume, bypassFlag=False, downscale=1, particle=None):
        """
        apply: Performs preprocessing of volume and reference
        @param volume: volume to be pre-processed
        @type volume: L{pytom_volume.vol}
        @param bypassFlag: Set if only bandpassFilter needed. False otherwise and all routines will be processed.
        @param downscale: not used anymore
        @param particle: particle Volume to be subtracted from input volume
        @type particle: L{pytom_volume.vol}
        @return: Returns modified volume
        @author: Thomas Hrabe  
        """
        
        from pytom_volume import vol
        
        if self._bandpassOn:

            from pytom.basic.filter import bandpassFilter

            # if frequencies specified in Nyquist, 0.5 being highest
            # fixed wrong adjustment of frequencies upon binning - FF
            if self._highestFrequency < 1:
                highestFrequency = self._highestFrequency*volume.sizeX()
                lowestFrequency = self._lowestFrequency*volume.sizeX()
            else:
                highestFrequency = self._highestFrequency
                lowestFrequency = self._lowestFrequency
            
            v = bandpassFilter(volume=volume, lowestFrequency=lowestFrequency, highestFrequency=highestFrequency,
                               bpf=0, smooth=self._bandpassSmooth)
            volume = v[0]

        if self._prerotateOn and (not bypassFlag):
            
            from pytom_volume import rotate
        
            rot = vol(volume.sizeX(),volume.sizeY(),volume.sizeZ())
            rotation = self.prerotate
            rotate( volume, rot, rotation[0], rotation[1], rotation[2])
            volume = rot 

        if self._weightingOn and (not bypassFlag):
            
            from pytom_volume import read
            from pytom_freqweight import weight
            from pytom.basic.fourier import fft, ifft
            
            wedgeSum = read(self._weightingFile)
            
            fVolume     = fft(volume)
            weighting   = weight(wedgeSum)
            weighting.apply(fVolume)     
            volume      = ifft(fVolume)
        
        if self._substractParticle and particle.__class__ == vol:
            volume -= particle            

        if self._taper >0:
            from pytom.tools.macros import volumesSameSize
            if self._taperMask is None or not volumesSameSize(volume,self._taperMask):
                from pytom.basic.functions import taper_edges
                volume, self._taperMask = taper_edges(volume, self._taper)
            else:
                volume = volume * self._taperMask
                
        return volume
    
    def noPrefiltering(self):
        """
        noPrefiltering : Disables prefiltering of particle and reference
        @author: Thomas Hrabe
        """
        self._bandpassOn = False
    
    def noWeighting(self):
        """
        noWeighting : Disables weighting of particle
        @author: Thomas Hrabe
        """
        self._weightingOn = False
    
    def noRotation(self):
        """
        noRotation : Disables prerotation of particle
        @author: Thomas Hrabe
        """    
        self._prerotateOn = False
    
    
    def setRotation(self,rotation=-1):
        """
        setRotation:
        @author: Thomas Hrabe
        """
        self._prerotateOn = not (isinstance(rotation, (int, long)))
        self._prerotate = rotation
        
    def setWeighting(self,weightingFile=''):
        """
        setWeighting
        @author: Thomas Hrabe
        """
        self._weightingOn = len(weightingFile) > 0
        self._weightingFile = weightingFile
    
    def setBandpass(self,lowestFrequency = None,highestFrequency = None,smooth=-1):
        """
        setBandpass
        @author: Thomas Hrabe
        """
        
        if not lowestFrequency and not highestFrequency:
            self._bandpassOn = False
            self._lowestFrequency = -1
            self._highestFrequency = -1
            self._bandpassFilter = None
            self._bandpassSmooth = -1
        elif lowestFrequency > highestFrequency:
            print 'Lowest  frequency: ', lowestFrequency
            print 'Highest frequency: ', highestFrequency
            raise RuntimeError('Preprocessing object: lowest frequency > highest frequency. Abort!')
        else:
            self._bandpassOn = 0.0 <= lowestFrequency < highestFrequency 
            self._lowestFrequency = lowestFrequency
            
            self._highestFrequency = highestFrequency
            self._bandpassFilter = None
            self._bandpassSmooth = smooth
    
    
    def fromXML(self,xmlObj):
        """
        fromXML : Assigns values to job attributes from XML object
        @param xmlObj: A xml object  
        @author: Thomas Hrabe 
        """        
        filter = xmlObj.xpath('Bandpass')
        
        if len(filter) > 0:
            filter = filter[0]
            self._bandpassOn = True
            self._lowestFrequency = float(filter.get('LowestFrequency'))
            self._highestFrequency = float(filter.get('HighestFrequency'))
            self._bandpassSmooth = float(filter.get('Smooth'))
        else:
            self._bandpassOn = False
    
        prerotation = xmlObj.xpath('Prerotation')
        
        if len(prerotation) > 0:
            prerotation = prerotation[0]
            
            self._prerotateOn = True
            
            z1 = float(prerotation.get('Phi'))
            z2 = float(prerotation.get('Psi'))
            x = float(prerotation.get('Theta'))
            
            self._prerotate = [z1,z2,x]
            
        else:
            self._prerotateOn = False
        
        weighting = xmlObj.xpath('WedgeWeight')
        
        if len(weighting) > 0:
            
            weighting = weighting[0]
            
            self._weightingOn = True
            
            self._weightingFile = weighting.get('File')
    
        else: 
            self._weightingOn = False
            
    def toXML(self,doc=-1):    
        """
        toXML : Compiles a XML file from preprocessing object
        @param doc: Optional XML document 
        @return: Preprocessing description in XML
        @author: Thomas Hrabe
        """ 
        from lxml import etree
        
        preObj = etree.Element("Preprocessing")
        
        if self._bandpassOn:
            
            filter = etree.Element('Bandpass',LowestFrequency = str(float(self._lowestFrequency)),HighestFrequency = str(float(self._highestFrequency)),Smooth = str(float(self._bandpassSmooth)))
            
            preObj.append(filter)
            
        if self._prerotateOn:
            
            prerotate = etree.Element('Prerotation')
            
            prerotate.set('Phi',str(self._prerotate[0]))
            prerotate.set('Psi',str(self._prerotate[1]))
            prerotate.set('Theta',str(self._prerotate[2]))
        
            preObj.append(prerotate)
        
        if self._weightingOn:
            
            weighting = etree.Element('WedgeWeight', File = self._weightingFile)
            
            preObj.append(weighting)
        
        return preObj
    
    
class Prerotation(PyTomClass):    
    """
    Prerotation: Performs first alignment of particles on a particleList for membrane proteins. We assume the particles were selected from a viron / nucleus (sphere like surface) 
    """
    def __init__(self,particleList = None,centerVector = None,symmetry = None):
        """
        @param particleList: The particle List 
        @param centerVector: Center of viron (sphere) in origin of particles 
        @param symmetry: We might already know some symmetry properties of the particles
        """
        from pytom.basic.structures import PointSymmetry
        
        self._particleList = particleList or []
        self._centerVector = centerVector or []
        self._symmetry = symmetry or PointSymmetry(1)
        
    def getParticleList(self):
        return self._particleList

    def calculatePolarCoordinates(self,randomizeZ1=False):
        """
	@deprecated: estimateOrienationOnSphere
        """
        self.estimateOrienationOnSphere( randomizeZ1)

    def estimateOrienationOnSphere(self, randomizeZ1=False):
        """
        estimateOrienationOnSphere: Determines prerotation info for each member in particle list. \
        Angles in rotation are set accordingly.
        @param randomizeZ1: Randomize Z1 rotation (default is false)    
        @type randomizeZ1: L{bool}
        """
        from pytom.angles.angleFnc import z2XFromCenterVector
        
        numberParticles = len(self._particleList);
        
        if randomizeZ1:
            import random
            random.seed()
        
        for i in xrange(numberParticles):
            particle = self._particleList[i]
            origin = particle.getPickPosition()
            [z2,x] = z2XFromCenterVector(origin.toVector(),self._centerVector)
            rotation = particle.getRotation()
            
            if randomizeZ1:
                rotation.setZ1(random.randint(0,360))
                
            rotation.setZ2(z2)
            rotation.setX(x)
            particle.setRotation(rotation)
            self._particleList[i] = particle
        
    def computeAverage(self,averageFilename):
        """
        computeAverage: Computes the average through L{pytom.alignment.aligner.ExMax}.expectation 
        @param averageFilename: The resulting file
        """
        from pytom.alignment.structures import ExpectationJob
        from pytom.alignment.ExMaxAlignment import ExMaxWorker
        
        pl = self._symmetry.apply(self._particleList)
        
        job = ExpectationJob(pl,averageFilename)
        
        worker = ExMaxWorker()
        worker.fromJob(job)
        
        worker.run()
