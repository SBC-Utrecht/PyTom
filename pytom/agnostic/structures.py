from pytom.gpu.initialize import xp, device
from pytom.basic.structures import Rotation as RotationPytomC
from pytom.agnostic.io import read
from pytom.agnostic.transform import fftshift, fourier_reduced2full
from pytom.agnostic.filter import applyFourierFilter
import scipy.optimize


class PyTomClassError(Exception):

    def __init__(self, value):
        self._value = value

    def __str__(self):
        print(self._value)


class PyTomClass(object):
    """
    PyTomClass : A virtual class.
    G{classtree PyTomClass}
    """

    def fromXML(self, xmlObj):
        """
        fromXML: Initializes child object from a XML element. You MUST overload this function in a child class.
        @param xmlObj: A xml object
        @type xmlObj: L{lxml.etree._Element}
        """
        print('fromXML: This is an abstract method!')
        assert False

    def toXML(self):
        """
        toXML: Serializes child object to a XML element. You MUST overload this function in a child class.
        """
        print('toXML: This is an abstract method!')
        assert False
        # raise RuntimeError("Function should not be called")

    def fromStr(self, string):
        """
        fromStr : Creates object from XML String. Do NOT overload this function in any child!
        @param string: The XML String
        @type string: str
        @author: Thomas Hrabe
        """
        from lxml import etree

        root = etree.fromstring(string)

        self.fromXML(root)

    def __str__(self):
        """
        __str__ : Prints object as XML String. Do NOT overload this function in any child!
        @author: Thomas Hrabe
        """
        from lxml import etree

        tree = self.toXML()
        self._xmlString = etree.tostring(tree, pretty_print=True).decode("utf-8")[:-1]

        return self._xmlString

    def fromXMLFile(self, filename):
        """
        fromXMLFile: Configures object according to file. Do NOT overload this function in any child!
        @param filename: Absolute / relative path to file
        @type filename: L{str}
        @author: Thomas Hrabe
        """
        from pytom.tools.files import readStringFile

        self._XMLfilename = filename
        lines = readStringFile(filename)

        self.fromStr(lines)

    def toXMLFile(self, filename):
        """
        toXMLFile : Dumps object to XML file. Do NOT overload this function in any child!
        @param filename: Absolute / relative path to file
        @type filename: L{str}
        @author: Thomas Hrabe
        """
        import pytom
        self._XMLfilename = filename
        versionString = '<!-- PyTom Version: ' + pytom.__version__ + ' -->\n'

        file = open(filename, "w")

        file.write(versionString + str(self))

        file.close()

    def __eq__(self, otherObject):
        """
        __eq__: Checks for object equality. Do NOT overload this function in any child!
        """
        return self.__str__() == otherObject.__str__()

    def __ne__(self, otherObject):
        """
        __ne__: Checks for object inequality. Do NOT overload this function in any child!
        """
        return (not self.__eq__(otherObject))

    def copy(self):
        """
        copy: Copies values of self to new object
        """
        import copy

        return copy.deepcopy(self)

    def xpath(self, query):
        """
        xpath: Wrapper for xml xpath function. Processes query on this objects XML. Do NOT overload this function in any child!
        @param query: This query will be processed
        @return: xml object - query result
        @author: Thomas Hrabe
        """
        xml = self.toXML()

        return xml.xpath(query)

    def xsltTransform(self, xsltStringIO):
        """
        xsltTransform: Performs a arbitrary xslt transformation on self.toXML(). Do NOT overload this function in any child!
        @param xsltStringIO: The transformation code as string
        @type xsltStringIO: StringIO object
        """

        from lxml import etree
        xslt_doc = etree.parse(xsltStringIO)
        transform = etree.XSLT(xslt_doc)

        selfXML = self.toXML()

        xsltResult = transform(selfXML)

        return xsltResult

    def check(self):
        """
        check: Performs logical check on self to make sure that job settings are correct. Will display error message if not.
        """
        print('check: This is an abstract method!')
        assert False

    def toHTMLFile(self, filename, xsltFile=None, xsltTransform=None):
        """
        toHTMLFile: Transforms object to HTML file according to xslt transformation. Throws exception if xsltFile and xsltTransform are None
        @param filename: Destination of HTML object
        @param xsltFile: File with the xslt transformation. Default is None. If None, you must specify xsltTransform
        @param xsltTransform: XSLT object. If None, then xsltFile must be set.
        """
        from lxml import etree
        if xsltFile == None and xsltTransform == None:
            raise RuntimeError('You did not provide a valid XSLT transformation.')

        if xsltFile:
            from pytom.tools.files import readStringFile
            import io
            xsltString = readStringFile(xsltFile)
            xsltTransform = io.StringIO(xsltString)
        else:
            raise RuntimeError('You did not provide a valid XSLT transformation.')

        selfHTML = self.xsltTransform(xsltTransform)

        htmlString = etree.tostring(selfHTML, pretty_print=True)

        file = open(filename, 'w')

        file.write(htmlString)

        file.close()


class Preprocessing(PyTomClass):
    """
    Preprocessing: Defines procedures used for preprocessing before each alignment
    step such as - bandpass filtering, - prerotation of volume, - rotation of weighting.
    The highest frequency in the bandpass filter will never be lower than 0.1, even
    when specified lower.
    Check code in apply function.
    """

    def __init__(self, lowestFrequency=None, highestFrequency=None, smooth=0,
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

    def setHighestFrequency(self, highestFrequency):
        self._highestFrequency = highestFrequency

    def setLowestFrequency(self, lowestFrequency):
        self._lowestFrequency = lowestFrequency

    def setSubstractParticle(self, on):
        self._substractParticle = on

    def apply(self, volume, bypassFlag=False, downscale=1, particle=None):
        """
        apply: Performs preprocessing of volume and reference
        @param volume: volume to be pre-processed
        @type volume: L{pytom.lib.pytom_volume.vol}
        @param bypassFlag: Set if only bandpassFilter needed. False otherwise and all routines will be processed.
        @param downscale: not used anymore
        @param particle: particle Volume to be subtracted from input volume
        @type particle: L{pytom.lib.pytom_volume.vol}
        @return: Returns modified volume
        @author: Thomas Hrabe
        """

        if self._bandpassOn:

            from pytom.agnostic.filter import bandpass as bandpassFilter

            # if frequencies specified in Nyquist, 0.5 being highest
            # fixed wrong adjustment of frequencies upon binning - FF
            if self._highestFrequency < 1:
                highestFrequency = self._highestFrequency * volume.shape[0]
                lowestFrequency = self._lowestFrequency * volume.shape[0]
            else:
                highestFrequency = self._highestFrequency
                lowestFrequency = self._lowestFrequency

            v = bandpassFilter(volume, low=lowestFrequency, high=highestFrequency,
                               smooth=self._bandpassSmooth)
            volume = v

        if self._prerotateOn and (not bypassFlag):
            from pytom.voltools import transform

            rot = xp.zeros_like(volume)
            rotation = [self.prerotate.getZ1(), self._prerotate.getX(), self._prerotate.getZ2()]
            transform(volume, output=rot, rotation=rotation, rotation_order='rzxz', device=device, interpolation='filt_bspline')
            volume = rot

        if self._weightingOn and (not bypassFlag):
            from pytom.agnostic.io import read
            from pytom.agnostic.filter import applyFourierFilter

            # TODO What exactly is read here? I interpreted this as a volume
            wedgeSum = read(self._weightingFile)

            applyFourierFilter(volume, wedgeSum)

        if self._substractParticle and particle.__class__ == vol:
            volume -= particle

        if self._taper > 0:
            from pytom.agnostic.tools import volumesSameSize
            if self._taperMask is None or not volumesSameSize(volume, self._taperMask):
                from pytom.agnostic.tools import taper_edges
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

    def setRotation(self, rotation=-1):
        """
        setRotation:
        @author: Thomas Hrabe
        """
        self._prerotateOn = not (isinstance(rotation, int))
        self._prerotate = rotation

    def setWeighting(self, weightingFile=''):
        """
        setWeighting
        @author: Thomas Hrabe
        """
        self._weightingOn = len(weightingFile) > 0
        self._weightingFile = weightingFile

    def setBandpass(self, lowestFrequency=None, highestFrequency=None, smooth=-1):
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
            print('Lowest  frequency: ', lowestFrequency)
            print('Highest frequency: ', highestFrequency)
            raise RuntimeError('Preprocessing object: lowest frequency > highest frequency. Abort!')
        else:
            self._bandpassOn = 0.0 <= lowestFrequency < highestFrequency
            self._lowestFrequency = lowestFrequency

            self._highestFrequency = highestFrequency
            self._bandpassFilter = None
            self._bandpassSmooth = smooth

    def fromXML(self, xmlObj):
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

            self._prerotate = [z1, z2, x]

        else:
            self._prerotateOn = False

        weighting = xmlObj.xpath('WedgeWeight')

        if len(weighting) > 0:

            weighting = weighting[0]

            self._weightingOn = True

            self._weightingFile = weighting.get('File')

        else:
            self._weightingOn = False

    def toXML(self, doc=-1):
        """
        toXML : Compiles a XML file from preprocessing object
        @param doc: Optional XML document
        @return: Preprocessing description in XML
        @author: Thomas Hrabe
        """
        from lxml import etree

        preObj = etree.Element("Preprocessing")

        if self._bandpassOn:
            filter = etree.Element('Bandpass', LowestFrequency=str(float(self._lowestFrequency)),
                                   HighestFrequency=str(float(self._highestFrequency)),
                                   Smooth=str(float(self._bandpassSmooth)))

            preObj.append(filter)

        if self._prerotateOn:
            prerotate = etree.Element('Prerotation')

            prerotate.set('Phi', str(self._prerotate[0]))
            prerotate.set('Psi', str(self._prerotate[1]))
            prerotate.set('Theta', str(self._prerotate[2]))

            preObj.append(prerotate)

        if self._weightingOn:
            weighting = etree.Element('WedgeWeight', File=self._weightingFile)

            preObj.append(weighting)

        return preObj


class Mask(PyTomClass):
    """
    Mask: A mask object. Used to control whether mask is spherical or not.
    If not, mask will return rotated, asymetrical masks for alignment and picking.
    """

    def __init__(self, filename=None, isSphere=True, binning=1):
        """
        Mask(filename=None,isSphere=True,binning =1 )
        @param filename: filename of mask volume (NOT XML, but .em or .mrc/,ccp4)
        @type filename: string
        @param isSphere: is it a sphere?
        @type isSphere: bool
        @param binning: binning factor
        @type binning: int
        """
        self._filename = filename or ''
        self._volume = None
        self._isSphere = isSphere
        self._binning = binning

    def isSphere(self):
        """
        isSphere: Returns true when this mask is a sphere.
        """
        return self._isSphere

    def getFilename(self):
        return self._filename;

    def setBinning(self, binning):
        """
        set binning for mask usage
        """
        self._binning = binning

    def setDownscale(self, scaleFactor):
        """to be deprecated
        @deprecated: Use setBinning instead!
        """
        self.setBinning(scaleFactor)

    def getVolume(self, rotation=None, bufferedRead=False):
        """
        getVolume: Returns this mask's volume. The volume returned is rotated only
        if rotation is set and self._isSphere == False.

        @param rotation: L{pytom.basic.structures.Rotation}
        @deprecated: bufferedRead: Not supported anymore
        """

        if self._volume is None:
            from pytom.tools.files import checkFileExists

            if not checkFileExists(self._filename):
                raise IOError('Could not find mask named: ' + self._filename)

            from pytom.agnostic.io import read
            from pytom.tools.files import checkFileExists
            from pytom.agnostic.transform import resize

            if not checkFileExists(self._filename):
                raise IOError('Particle ' + self._filename + ' does not exist!')

            try:
                self._volume = read(self._filename)
                if self._binning != 1:
                    self._volume, volumef = resize(volume=self._volume, factor=1. / self._binning,
                                                   interpolation='Fourier')
            except RuntimeError:
                raise RuntimeError('Error reading file ' + self._volume)

            # self._volume = read(self._filename,0,0,0,0,0,0,0,0,0,self._binning,
            # self._binning,self._binning)

        if rotation and not (self._isSphere):
            from pytom.voltools import transform
            from pytom.agnostic.tools import determineRotationCenter

            if rotation.__class__ == list:
                from pytom.agnostic.structures import Rotation
                rotation = Rotation(rotation)

            rotationCenter = determineRotationCenter(self._filename, self._binning)

            maskRot = xp.zeros_like(self._volume.shape)
            transform(self._volume, output=maskRot, rotation=[rotation.getZ1(), rotation.getX(),rotation.getZ2()],
                      center=[rotationCenter[0], rotationCenter[1], rotationCenter[2]], device=device, interpolation='filt_bspline')

            return maskRot
        else:
            return self._volume

    def toXML(self):
        """
        toXML : Compiles a XML object from job object
        @rtype: L{lxml.etree._Element}
        @return: XML Object
        @author: Thomas Hrabe
        """

        from lxml import etree

        maskElement = etree.Element("Mask", isSphere=self._isSphere.__str__(),
                                    Filename=self._filename, Binning=str(int(self._binning)))

        return maskElement

    def fromXML(self, xmlObj=-1):
        """
        fromXML : Assigns values to job attributes from XML object
        @param xmlObj: A xml object
        @type xmlObj: L{lxml.etree._Element}
        @author: Thomas Hrabe
        """
        from lxml.etree import _Element

        if xmlObj.__class__ != _Element:
            raise TypeError('Mask: You must provide a valid Mask XML object.')

        self._filename = xmlObj.get('Filename')
        self._isSphere = xmlObj.get('isSphere') == 'True'
        self._binning = xmlObj.get('Binning')

        if self._binning is None:
            self._binning = xmlObj.get('Downscale')

        if self._binning is None:
            self._binning = 1

        self._binning = int(self._binning)
        del (self._volume)
        self._volume = None

    def check(self):

        from pytom.tools.files import checkFileExists

        if not checkFileExists(self._filename):
            raise IOError('Could not find mask file: ' + str(self._filename))

    def convert2numpy(self):
        return self


class Reference(PyTomClass):
    """
    Reference: Stores information about the current reference.
    """

    def __init__(self, referenceFile='', generatedByParticleList=None):
        """
        @param referenceFile: path to reference
        @type referenceFile: L{str}
        @param generatedByParticleList: ParticleList from which this reference was generated
        @type generatedByParticleList: L{pytom.agnostic.structures.ParticleList}
        """

        from pytom.tools.files import checkFileExists

        self._referenceFile = referenceFile

        if len(self._referenceFile) > 0:
            preWedgeName = self._referenceFile[:len(self._referenceFile) - 3] + '-PreWedge.em'
            wedgeSumName = self._referenceFile[:len(self._referenceFile) - 3] + '-WedgeSumUnscaled.em'
            filterName = self._referenceFile[:len(self._referenceFile) - 3] + '-Filter.dat'
        else:
            preWedgeName = ''
            wedgeSumName = ''
            filterName = ''

        if checkFileExists(preWedgeName):
            self._preWedgeFile = preWedgeName
        else:
            self._preWedgeFile = ''

        if checkFileExists(wedgeSumName):
            self._referenceWeighting = wedgeSumName
        else:
            self._referenceWeighting = ''

        if checkFileExists(filterName):
            self.setFSCFilter(filterName)
        else:
            self.setFSCFilter('')

        if generatedByParticleList:
            self._generatedByParticleList = generatedByParticleList
        else:
            from pytom.agnostic.structures import ParticleList
            self._generatedByParticleList = ParticleList('/')

    def hasGeneratedByInfo(self):
        """
        hasGeneratedByInfo
        """
        return len(self._generatedByParticleList) > 0

    def getReferenceFilename(self):
        """
        @deprecated: Use getFilename instead!
        """
        return self._referenceFile

    def getFilename(self):
        return self._referenceFile

    def setFilename(self, filename):
        """
        @param filename: name of reference file
        @type filename: L{str}
        """
        self._referenceFile = filename

    def getWeightingFilename(self):
        return self._referenceWeighting

    def getPreWedgeFilename(self):
        return self._preWedgeFile

    def getGeneratedByParticleList(self):
        return self._generatedByParticleList

    def hasWeighting(self):
        return self._referenceWeighting != ''

    def hasPreWedge(self):
        return self._preWedgeFile != ''

    def wasGeneratedBy(self, particle):
        """
        wasGeneratedBy: Determines whether this reference was generated by a certain particle
        @param particle: The particle
        @return: True if reference was generated by particle, False if not
        """
        try:
            if particle.__class__ == str:
                p = self._generatedByParticleList.getParticleByFilename(particle)
            else:
                p = self._generatedByParticleList.getParticleByFilename(particle.getFilename())

            return True
        except Exception:
            return False

    def updateReference(self, referenceFile, generatedByParticleList=None):
        """
        updateReference: updates current object
        @param referenceFile:
        @param generatedByParticleList: ParticleList from which this reference was generated
        @type generatedByParticleList: L{pytom.agnostic.structures.ParticleList}
        """
        from pytom.tools.files import checkFileExists

        if not checkFileExists(referenceFile):
            raise RuntimeError('Reference file does not exist!')

        if len(self._referenceFile) > 0:
            preWedgeName = self._referenceFile[:len(self._referenceFile) - 3] + '-PreWedge.em'
            wedgeSumName = self._referenceFile[:len(self._referenceFile) - 3] + '-WedgeSum.em'
        else:
            preWedgeName = ''
            wedgeSumName = ''

        if checkFileExists(preWedgeName):
            self._preWedgeFile = preWedgeName
        else:
            self._preWedgeFile = ''

        if checkFileExists(wedgeSumName):
            self._referenceWeighting = wedgeSumName
        else:
            self._referenceWeighting = ''

        if generatedByParticleList:
            self._generatedByParticleList = generatedByParticleList
        else:
            from pytom.agnostic.structures import ParticleList
            self._generatedByParticleList = ParticleList('/')

    def getVolume(self, binning=1):
        """
        getVolume:
        @param binning: binning factor
        @type binning: int
        @return: Reference volume
        """

        from pytom.agnostic.io import read
        from pytom.tools.files import checkFileExists
        from pytom.agnostic.transform import resize

        if not checkFileExists(self._referenceFile):
            raise IOError('Particle ' + self._referenceFile + ' does not exist!')

        try:
            reference = read(self._referenceFile)
            if binning != 1:
                reference, referencef = resize(volume=reference, factor=1. / binning, interpolation='Fourier')
        except RuntimeError:
            raise RuntimeError('Error reading file ' + self._referenceFile)

        return reference

    def getWeighting(self):
        """
        getVolume:
        @return: Sum of wedges
        """

        from pytom.agnostic.io import read
        from pytom.tools.files import checkFileExists

        if self._referenceWeighting == '':
            return self._referenceWeighting

        if not checkFileExists(self._referenceWeighting):
            raise IOError('Reference ' + self._referenceWeighting + ' does not exist!')

        return read(self._referenceWeighting)

    def getPreWedge(self):
        """
        getPreWedge:
        @return: Reference before wedge weighting
        """

        from pytom.agnostic.io import read
        from pytom.tools.files import checkFileExists

        if self._preWedgeFile == '':
            return self._preWedgeFile

        if not checkFileExists(self._preWedgeFile):
            raise IOError('PreWedge ' + self._preWedgeFile + ' does not exist!')

        return read(self._preWedgeFile)

    def subtractParticle(self, particle, binning=1, verbose=False):
        """
        subtractParticle: Subtract a particle from the current reference
        @param particle:
        @param binning: Bin result
        @param verbose: talkative
        @type verbose: bool
        @return: [newReferenceVolume,newSumOfWedges]
        @rtype: [L{pytom.lib.pytom_volume.vol},L{pytom.lib.pytom_volume.vol}]
        @change: more accurate binning now in Fourier space - FF
        """
        # if flag is set and both files exist, do subtract particle from reference

        from pytom.agnostic.io import read
        from pytom.gpu.gpuFunctions import applyFourierFilter as convolute
        from pytom.agnostic.filter import rotateWeighting
        from pytom.agnostic.tools import invert_WedgeSum
        from pytom.agnostic.filter import filter_volume_by_profile
        from pytom.basic.resolution import read_fscFromAscii

        # find particle with rotation and shift settings in the generatedByParticleList
        particle = self._generatedByParticleList.getParticleByFilename(particle.getFilename())

        # read original data
        preWedge = read(self.getPreWedgeFilename())# average without wedge weighting
        wedgeSum = read(self.getWeightingFilename())  # sum of wedges

        particleWedge = particle.getWedge()

        particleTransformed = particle.getTransformedVolume()
        if self.getFSCFilter() != '':
            fsc_fil = read_fscFromAscii(self.getFSCFilter())
            particleTransformed = filter_volume_by_profile(volume=particleTransformed, profile=fsc_fil)
        ptx,pty,ptz = particleTransformed.shape
        particleRotation = particle.getRotation()

        rotinvert = particleRotation.invert()
        if analytWedge:
            # > buggy version
            wedge = particleWedge.returnWedgeVolume(ptx, pty, ptz, False, rotinvert)
            # < buggy version
        else:
            wedge = rotateWeighting(particleWedge.returnWedgeVolume(ptx, pty, ptz, False),
                                    [rotinvert[0], rotinvert[2],rotinvert[1]], mask=None)


        # substract particle and wedge, consistent with pytom.alignment.alignmentFunctions.average
        newRefPreWedge = preWedge - particleTransformed
        newWedgeSum = wedgeSum - wedge
        lowval = lowlimit = .05 * newWedgeSum[0, 0, 0]
        invert_WedgeSum(invol=newWedgeSum, r_max=ptx// 2 - 2.,lowlimit=lowlimit, lowval=lowval)

        # weight new average
        referenceBig = convolute(newRefPreWedge, newWedgeSum)

        # bin to size
        if binning > 1:
            from pytom.agnostic.transform import resize  # , scale as scaleVolume
            from pytom.agnostic.io import readSubvolumeFromFourierspaceFile

            reference = resize(volume=referenceBig, factor=1. / float(binning), interpolation='Fourier')[0]
            # reference = scaleVolume(referenceBig,1/float(binning))


            rx,ry,rz = reference.shape

            newWedgeSum = readSubvolumeFromFourierspaceFile(newWedgeSum, rx, ry, rz)

        else:
            reference = referenceBig

        return [reference, newWedgeSum]

    def fromXML(self, xmlObj):
        """
        fromXML : Assigns values to job attributes from XML object
        @param xmlObj: A xml object
        @type xmlObj: L{lxml.etree._Element}
        @author: Thomas Hrabe
        @change: option to read only ParticleListFilename, FF
        """
        from lxml.etree import _Element

        if xmlObj.__class__ != _Element:
            raise TypeError('You must provide a valid XML-Reference object.')

        self._referenceFile = xmlObj.get('File')

        from pytom.agnostic.structures import ParticleList
        self._generatedByParticleList = ParticleList('/')

        generatedByXML = xmlObj.xpath('ParticleList')

        if len(generatedByXML) > 0:
            # onlyParticleListFilename==False
            if len(generatedByXML[0].xpath('Particle')) > 0:
                self._generatedByParticleList.fromXML(generatedByXML[0])
            else:
                particleListFilename = generatedByXML[0].get('Filename')
                if not particleListFilename:
                    raise IOError('Reference XML should either have Particles or Filename for ParticleList')
                self._generatedByParticleList.fromXMLFile(filename=particleListFilename)

        try:
            self._referenceWeighting = xmlObj.get('Weighting')
            self._preWedgeFile = xmlObj.get('PreWedge')

            if not self._referenceWeighting or not self._preWedgeFile:
                raise Exception()

        except:
            from pytom.tools.files import checkFileExists

            if len(self._referenceFile) > 0:
                preWedgeName = self._referenceFile[:len(self._referenceFile) - 3] + '-PreWedge.em'
                wedgeSumName = self._referenceFile[:len(self._referenceFile) - 3] + '-WedgeSum.em'
            else:
                preWedgeName = ''
                wedgeSumName = ''

            if checkFileExists(preWedgeName):
                self._preWedgeFile = preWedgeName
            else:
                self._preWedgeFile = ''

            if checkFileExists(wedgeSumName):
                self._referenceWeighting = wedgeSumName
            else:
                self._referenceWeighting = ''

    def setFSCFilter(self, fscFilterName=''):
        """
        set FSC-derived filter
        """
        self._fscFilterName = fscFilterName

    def getFSCFilter(self):
        return self._fscFilterName

    def toXML(self, onlyParticleListFilename=False):
        """
        toXML : Compiles a XML file from job object
        @author: Thomas Hrabe
        @param onlyParticleListFilename: print only particleList filename rather than all particles
        @type onlyParticleListFilename: L{bool}
        """
        from lxml import etree
        referenceElement = etree.Element('Reference', File=self._referenceFile, Weighting=self._referenceWeighting,
                                         PreWedge=self._preWedgeFile, FSCFilter=self._fscFilterName)

        if self._generatedByParticleList:
            if onlyParticleListFilename == False:
                referenceElement.append(self._generatedByParticleList.toXML())
            else:
                particleListElement = etree.Element('ParticleList', Filename=self._generatedByParticleList._XMLfilename)
                referenceElement.append(particleListElement)

        assert referenceElement.__class__ == etree._Element

        return referenceElement

    def check(self):

        from pytom.tools.files import checkFileExists

        if not checkFileExists(self._referenceFile):
            raise IOError('Could not find reference file: ' + str(self._referenceFile))


class Wedge(PyTomClass):
    """
    Wedge: used as an dummy class to distinguish between single tilt axis wedge and double tilt axis wedge in fromXML
    """

    def __init__(self, wedge_angles=[0.0, 0.0], cutoffRadius=0.0, tiltAxis='Y', smooth=0.0, wedge_3d_ctf_file='',
                 ctf_max_resolution=0.):
        """
        __init__: This constructor is compatible to L{pytom.agnostic.structures.SingleTiltWedge} and L{pytom.agnostic.structures.DoubleTiltWedge}.
        """
        if wedge_3d_ctf_file != '':
            # we need to initialize the 3d
            self._wedgeObject = Wedge3dCTF(filename=wedge_3d_ctf_file, ctf_max_resolution=ctf_max_resolution)
            self._type = 'Wedge3dCTF'
        else:
            try:

                if wedge_angles.__class__ == list and wedge_angles[0].__class__ == list and len(wedge_angles[0]) == 2 and \
                        wedge_angles[1].__class__ == list and len(wedge_angles[1]) == 2:
                    self._wedgeObject = DoubleTiltWedge(wedge_angles=wedge_angles, tiltAxis1='Y', rotation12=tiltAxis,
                                                        cutoffRadius=cutoffRadius, smooth=smooth)
                    self._type = 'DoubleTiltWedge'
                else:

                    raise RuntimeError('Do SingleTiltWedge')
            except:
                if wedge_angles.__class__ == list and wedge_angles[0].__class__ == list and wedge_angles[1].__class__ == list:
                    raise TypeError('Wrong parameters for SingleTiltWedge object error thrown by Wedge object!')
                self._wedgeObject = SingleTiltWedge(wedge_angles, cutoffRadius=cutoffRadius, tiltAxis=tiltAxis,
                                                    smooth=smooth)
                self._type = 'SingleTiltWedge'

    def returnWedgeFilter(self, wedgeSizeX=None, wedgeSizeY=None, wedgeSizeZ=None, rotation=None):
        """
        returnWedgeFilter : Returns a wedge filter specified by this object
        @param wedgeSizeX: volume size for x (original size)
        @param wedgeSizeY: volume size for y (original size)
        @param wedgeSizeZ: volume size for z (original size)
        @rtype: L{pytom.lib.pytom_freqweight.weight}
        @return: Weighting object. Remember, the wedge will be cutoff at size_x/2 if no cutoff provided in constructor or cutoff == 0!
        @author: Thomas Hrabe
        """
        if self._type == 'Wedge3dCTF':
            print('ERROR: wedge filter object not available for 3d ctf. exiting...')
            raise PyTomClassError
        else:
            return self._wedgeObject.returnWedgeFilter(wedgeSizeX, wedgeSizeY, wedgeSizeZ, rotation)

    def returnWedgeVolume(self, wedgeSizeX=None, wedgeSizeY=None, wedgeSizeZ=None, humanUnderstandable=False,
                          rotation=None):
        """
        returnWedgeVolume: Returns a wedge volume for later processing.
        @param wedgeSizeX: volume size for x (size of original volume in real space)
        @param wedgeSizeY: volume size for y (size of original volume in real space)
        @param wedgeSizeZ: volume size for z (size of original volume in real space)
        @param humanUnderstandable: if True (default is False), the volume will be transformed from reducedComplex to full and shifted afterwards
        @param rotation: rotation of wedge

        @rtype: L{pytom.lib.pytom_volume.vol}
        @author: Thomas Hrabe
        """

        return self._wedgeObject.returnWedgeVolume(wedgeSizeX, wedgeSizeY, wedgeSizeZ, humanUnderstandable, rotation)

    def getWedgeAngle(self):
        """
        get angles of wedge
        """
        if self._type == 'Wedge3dCTF':
            print('ERROR: wedge filter object not available for 3d ctf. exiting...')
            raise PyTomClassError
        else:
            return self._wedgeObject.getWedgeAngle()

    def getType(self):
        return self._type

    def apply(self, volume, rotation=None):
        """
        apply: Applies this wedge to a given volume
        @param volume: The volume to be filtered
        @type volume: L{pytom.lib.pytom_volume.vol} or L{pytom.lib.pytom_volume.vol_comp}
        @return: The filtered volume
        @rtype: L{pytom.lib.pytom_volume.vol}
        @author: Thomas Hrabe
        """

        return self._wedgeObject.apply(volume, rotation)

    def getWedgeObject(self):
        """
        getWedgeObject: Returns the underlying wedge object used
        """
        return self._wedgeObject

    def toSphericalFunc(self, b, radius):
        """Convert the wedge from real space to a spherical function.
        This function serves as an interface to FRM.
        @param b: Bandwidth of the spherical function.
        @return: a spherical function in numpy.array
        """
        return self._wedgeObject.toSphericalFunc(b, radius)

    def toXML(self):

        from lxml import etree

        wedgeElement = etree.Element('Wedge', Type=str(self._type))
        wedgeElement.append(self._wedgeObject.toXML())

        return wedgeElement

    def fromXML(self, xmlObj):
        """
        @param xmlObj: A xml object
        @type xmlObj: L{lxml.etree._Element}
        """
        if not xmlObj.tag in ['Wedge', 'WedgeInfo', 'SingleTiltWedge', 'DoubleTiltWedge', 'GeneralWedge', 'Wedge3dCTF']:
            raise TypeError('You must provide a XML-Wedge object.')

        wedgeType = []
        if xmlObj.tag == 'Wedge':
            wedgeType = xmlObj.get('Type')
            xmlObj = xmlObj[0]  # get the actual wedge inside
        else:
            wedgeType = xmlObj.tag

        if wedgeType == 'SingleTiltWedge':
            # is it a SingleTiltWedge?
            self._wedgeObject = SingleTiltWedge()
            self._wedgeObject.fromXML(xmlObj)

        elif wedgeType == 'WedgeInfo':
            # is it a old type wedge?
            self._wedgeObject = SingleTiltWedge()
            wedgeType = 'SingleTiltWedge'
            self._wedgeObject.fromXML(xmlObj)

        elif wedgeType == 'DoubleTiltWedge':
            # is it a DoubleTiltWedge?
            self._wedgeObject = DoubleTiltWedge()
            self._wedgeObject.fromXML(xmlObj)

        elif wedgeType == 'GeneralWedge':
            self._wedgeObject = GeneralWedge()
            self._wedgeObject.fromXML(xmlObj)

        elif wedgeType == 'Wedge3dCTF':
            self._wedgeObject = Wedge3dCTF()
            self._wedgeObject.fromXML(xmlObj)

        else:
            raise TypeError('The xml object provided does not contain wedge information!')

        self._type = wedgeType

    def setWedgeAngles(self, wedgeangles):
        if self._type == 'Wedge3dCTF':
            print('ERROR: wedge filter object not available for 3d ctf. exiting...')
            raise PyTomClassError
        else:
            self._wedgeObject.setWedgeAngles(wedgeangles)

    def getFilename(self):
        if self._type != 'Wedge3dCTF':
            print('filename only implemented for 3d ctf wedge')
        else:
            return self._wedgeObject.getFilename()

    def setFilename(self, filename):
        if self._type != 'Wedge3dCTF':
            print('filename only implemented for 3d ctf wedge')
        else:
            return self._wedgeObject.setFilename(filename)


class SingleTiltWedge(PyTomClass):
    """
    SingleTiltWedge : Saves all Wedge info such as wedge angle, and wedge \
    rotation. Supports asymmetric wedges.
    @author: Thomas Hrabe
    """

    def __init__(self, wedge_angle=0.0, rotation=None, cutoffRadius=0.0, tiltAxis='Y', smooth=0.0):
        """
        @param wedge_angle: The wedge angle. In wedge halfs.
        @type wedge_angle: either float (for symmetric wedge) or [float,float]  for \
            asymmetric wedge.
        @param rotation: deprecated!!! Only leave it here for compatibility reason.
        @deprecated: rotation
        @param cutoffRadius: Maximum frequency (in pixels) up to where wedge will \
            be generated. If omitted / 0, filter is fixed to size/2.
        @param tiltAxis: Default will be Y. You can choose between the strings "X",\
            "Y" or "custom". You can also specify a Rotation object
        @type tiltAxis: str or L{pytom.agnostic.structures.Rotation}
        @param smooth: Smoothing size of wedge at the edges in degrees. Default is 0.
        """
        if wedge_angle.__class__ == list:
            assert wedge_angle[0] >= 0
            assert wedge_angle[1] >= 0
            self._wedge_angle1 = wedge_angle[0]
            self._wedge_angle2 = wedge_angle[1]
        else:
            assert wedge_angle >= 0
            self._wedge_angle1 = wedge_angle
            self._wedge_angle2 = wedge_angle

        if rotation:
            pass  # print("average: Warning - input rotation will not be used because deprecated!")
        self._cutoffRadius = cutoffRadius
        self._tiltAxisRotation = Rotation(0.0, 0.0, 0.0)

        if tiltAxis.__class__ == str:
            if tiltAxis == 'Y':
                self._tiltAxisRotation = Rotation(0.0, 0.0, 0.0)
            elif tiltAxis == 'X':
                self._tiltAxisRotation = Rotation(90.0, 0.0, 0.0)
            elif tiltAxis == 'Z':
                self._tiltAxisRotation = Rotation(0.0, 0.0, 90.0)
        elif tiltAxis.__class__ == Rotation:
            self._tiltAxisRotation = tiltAxis
        else:
            raise TypeError(
                'TiltAxis parameter for SingleTiltWedge must either be string = "X" , "Y" , "Z" or a Rotation object!')

        self._smooth = smooth

        # for FRM
        self._wedge_vol = None  # cache for storing the wedge volume
        self._bw = None
        self._sf = None

    def setTiltAxisRotation(self, rotation):
        self._tiltAxisRotation = rotation

    def getWedgeAngle(self):
        """
        getWedgeAngle : Getter method for wedge_angle
        @return: self.wedge_angle in openingAngle/2. If its an asymmetric wedge, a list [angle1,angle2] will be returned
        @author: Thomas Hrabe
        """
        if self._wedge_angle1 == self._wedge_angle2:
            return self._wedge_angle1
        else:
            return [self._wedge_angle1, self._wedge_angle2]

    def getTiltAxisRotation(self):
        from pytom.agnostic.structures import Rotation

        return Rotation(self._tiltAxisRotation[0], self._tiltAxisRotation[1], self._tiltAxisRotation[2])

    def getTiltAxis(self):
        """
        @deprecated: Use getTiltAxisRotation instead
        """
        return self.getTiltAxisRotation()

    def returnWedgeFilter(self, wedgeSizeX, wedgeSizeY, wedgeSizeZ, rotation=None):
        """
        returnWedgeFilter : Returns a wedge filter specified by this object
        @param wedgeSizeX: volume size for x (original size)
        @param wedgeSizeY: volume size for y (original size)
        @param wedgeSizeZ: volume size for z (original size)
        @param rotation: Apply rotation to the wedge
        @type rotation: L{pytom.agnostic.structures.Rotation}
        @rtype: L{pytom.lib.pytom_freqweight.weight}
        @return: Weighting object. Remember, the wedge will be cutoff at size_x/2 if no cutoff provided in constructor or cutoff == 0!
        @author: Thomas Hrabe
        """
        # This import from own file ....
        if self._cutoffRadius == 0.0:
            cut = wedgeSizeX // 2
        else:
            cut = self._cutoffRadius
        weightObject = Weight(self._wedge_angle1, self._wedge_angle2, cut, wedgeSizeX, wedgeSizeY, wedgeSizeZ,
                              self._smooth)

        if (not isinstance(rotation, Rotation)) and (not isinstance(rotation, RotationPytomC)):
            rotation = Rotation()

        wedgeRotation = rotation * Rotation(self._tiltAxisRotation[0], self._tiltAxisRotation[1],
                                            self._tiltAxisRotation[2])
        weightObject.rotate(wedgeRotation.getZ1(), wedgeRotation.getZ2(), wedgeRotation.getX())

        return weightObject

    def returnWedgeVolume(self, wedgeSizeX, wedgeSizeY, wedgeSizeZ, humanUnderstandable=False, rotation=None):
        """
        returnWedgeVolume: Returns a wedge volume for later processing.
        @param wedgeSizeX: volume size for x (size of original volume in real space)
        @param wedgeSizeY: volume size for y (size of original volume in real space)
        @param wedgeSizeZ: volume size for z (size of original volume in real space)
        @param humanUnderstandable: if True (default is False), the volume will \
        be transformed from reducedComplex to full and shifted afterwards
        @param rotation: rotation of wedge

        @rtype: L{pytom.lib.pytom_volume.vol}
        @author: Thomas Hrabe
        """
        wedgeVolume = None

        if abs(self._wedge_angle1) < 1E-4 and abs(self._wedge_angle2) < 1E-4:

            if not humanUnderstandable:
                wedgeVolume = xp.ones((wedgeSizeX, wedgeSizeY, wedgeSizeZ // 2 + 1), dtype=xp.float32)
            else:
                wedgeVolume = xp.ones((wedgeSizeX, wedgeSizeY, wedgeSizeZ),dtype=xp.int32)

        else:
            # add a custom rotation to the previously specified one
            wedgeFilter = self.returnWedgeFilter(wedgeSizeX, wedgeSizeY, wedgeSizeZ, rotation)

            # return only the reduced complex part.
            wedgeVolume = wedgeFilter.getWeightVolume(True)

            if humanUnderstandable:
                # applies hermitian symmetry to full and ftshift so that even humans do understand
                from pytom.agnostic.transform import fourier_reduced2full as reducedToFull
                wedgeVolume = reducedToFull(wedgeVolume, (wedgeSizeZ % 2))
                wedgeVolume = xp.fft.fftshift(wedgeVolume)  # fftshift to shift corner to center of spectrum

        return wedgeVolume

    def apply(self, volume, rotation=None):
        """
        apply: Applies this wedge to a given volume
        @param volume: The volume to be filtered
        @type volume: L{pytom.lib.pytom_volume.vol} or L{pytom.lib.pytom_volume.vol_comp}
        @param rotation: rotate the wedge
        @type rotation: L{pytom.agnostic.structures.Rotation}
        @return: The filtered volume
        @rtype: L{pytom.lib.pytom_volume.vol}
        @author: Thomas Hrabe
        """
        if not volume.__class__ == xp.array((1)).__class__:
            raise TypeError('SingleTiltWedge: You must provide a xp.array here!')

        if self._wedge_angle1 > 0 or self._wedge_angle2 > 0:
            from pytom.agnostic.filter import applyFourierFilter as filter

            wedgeFilter = self.returnWedgeVolume(volume.shape[0], volume.shape[1], volume.shape[2], rotation)
            result = filter(volume, wedgeFilter)

            return result
        else:
            return volume

    def toSphericalFunc(self, b, radius=None):
        """Convert the wedge from real space to a spherical function.
        This function serves as an interface to FRM.
        @param b: Bandwidth of the spherical function.
        @param radius: Radius in frequency. Not used for SingleTiltWedge.
        @return: a spherical function in numpy.array
        """
        raise NotImplementedError

    def fromXML(self, xmlObj):
        """
        fromXML: Assigns values to result attributes from XML object
        @param xmlObj: A xml object
        @type xmlObj: L{lxml.etree._Element}
        @author: Thomas Hrabe
        """

        from lxml.etree import _Element

        if xmlObj.__class__ != _Element:
            raise TypeError('Is not a lxml.etree._Element! You must provide a valid XML-Wedge object.')

        if (xmlObj.get('Angle1') == None):
            self._wedge_angle1 = float(xmlObj.get('Angle'))
            self._wedge_angle2 = float(xmlObj.get('Angle'))
        else:
            self._wedge_angle1 = float(xmlObj.get('Angle1'))
            self._wedge_angle2 = float(xmlObj.get('Angle2'))

        self._cutoffRadius = float(xmlObj.get('CutoffRadius'))

        if (xmlObj.get('Smooth') != None):
            self._smooth = float(xmlObj.get('Smooth'))

        rotation = xmlObj.xpath('TiltAxisRotation')

        if len(rotation) > 0:
            rotation = rotation[0]
            self._tiltAxisRotation = [.0, .0, .0]
            self._tiltAxisRotation[0] = float(rotation.get('Z1'))
            self._tiltAxisRotation[1] = float(rotation.get('Z2'))
            self._tiltAxisRotation[2] = float(rotation.get('X'))
        else:
            self._tiltAxisRotation = [.0, .0, .0]

    def toXML(self):
        """
        toXML : Compiles a XML file from result object
        rtype : L{lxml.etree._Element}
        @author: Thomas Hrabe
        """
        from lxml import etree

        wedgeElement = etree.Element('SingleTiltWedge', Angle1=str(self._wedge_angle1),
                                     Angle2=str(self._wedge_angle2), CutoffRadius=str(self._cutoffRadius),
                                     Smooth=str(self._smooth))

        if ((not isinstance(self._tiltAxisRotation, int)) or
                (not self._tiltAxisRotation.__class__ == float)):
            rotationElement = etree.Element('TiltAxisRotation')

            rotationElement.set('Z1', str(float(self._tiltAxisRotation[0])))
            rotationElement.set('Z2', str(float(self._tiltAxisRotation[1])))
            rotationElement.set('X', str(float(self._tiltAxisRotation[2])))

            wedgeElement.append(rotationElement)

        return wedgeElement

    def setWedgeAngles(self, wedgeangles=(30, 30)):
        if len(wedgeangles) == 1:
            wedgeangles = [wedgeangles[0], wedgeangles[1]]
        self._wedge_angle1, self._wedge_angle2 = wedgeangles


class WedgeInfo(SingleTiltWedge):
    """
    WedgeInfo:
    @deprecated: Used for backward compatibility only..
    """


class DoubleTiltWedge(SingleTiltWedge):
    """
    DoubleTiltWedge: Represents a wedge determined for a double tilt series of projections
    """

    def __init__(self, wedge_angles=[[0, 0], [0, 0]], tiltAxis1='Y', rotation12=[90, 0, 0], cutoffRadius=0.0,
                 smooth=0.0):
        """Initialize a double tilt wedge
        @param wedge_angles: List of the tilt parameters [[tilt1.1,tilt1.2],[tilt2.1,tilt2.2]]. \
The missing region. All should be positive degrees.
        @param tiltAxis1: Specify tilt axis of first tilt here. The first tilt axis \
will be Y, the second might differ. Unless specified otherwise, it will \
be set to X -> the double tilt wedge is perfectly rotated by 90 degrees.
        @type tiltAxis1: either string (default: 'Y'), or a rotation object \
L{pytom.agnostic.structures.Rotation}
        @param rotation12: Rotation of 2nd tilt series with respect to 1st
        @type rotation12: rotation object L{pytom.agnostic.structures.Rotation} or 3-dim list
        @param cutoffRadius: Maximum frequency (in pixels) up to where wedge will \
be generated. If omitted / 0, filter is fixed to size/2.
        @param smooth: Smoothing size of wedge at the edges in degrees. Default is 0.
        """
        from pytom.agnostic.structures import SingleTiltWedge
        tilt1 = wedge_angles[0]
        tilt2 = wedge_angles[1]

        if tilt1.__class__ != list or tilt2.__class__ != list:
            raise TypeError('These are the wrong parameters for double tilt wedge!')

        if tiltAxis1.__class__ == str:
            if tiltAxis1 == 'Y':
                tiltAxis1 = Rotation(0, 0, 0)
        if tiltAxis1.__class__ == list:
            if len(tiltAxis1 != 3):
                raise TypeError('DoubleTiltWedge: tiltAxis1 must be 3-dim if chosen as array')
            else:
                tiltAxis1 = Rotation(tiltAxis1[0], tiltAxis1[1], tiltAxis1[2])

        if rotation12.__class__ == str:
            if rotation12 in ['X', 'Y']:
                # catch this error when somebody forgets to specify correctly and force wedges to be perpendicular.
                # happens only when double tilt wedge is specified incorrectly from wedge object
                rotation12 = Rotation(90, 0, 0)
            else:
                raise TypeError('DoubleTiltWedge: rotation12 must be 3-dim if chosen as array')

        if rotation12.__class__ == list:
            if len(rotation12) != 3:
                raise TypeError('DoubleTiltWedge: rotation12 must be 3-dim if chosen as array')
            else:
                rotation12 = Rotation(rotation12[0], rotation12[1], rotation12[2])

        tiltAxis2 = (rotation12 * tiltAxis1)

        self._wedge1 = SingleTiltWedge(tilt1, None, cutoffRadius, tiltAxis=tiltAxis1, smooth=smooth)
        self._wedge2 = SingleTiltWedge(tilt2, None, cutoffRadius, tiltAxis=tiltAxis2, smooth=smooth)

        # for FRM
        self._wedge_vol = None  # cache for storing the wedge volume
        self._bw = None
        self._sf = None

    def returnWedgeVolume(self, wedgeSizeX, wedgeSizeY, wedgeSizeZ, humanUnderstandable=False, rotation=None):
        """
        returnWedgeVolume: Returns a wedge volume for later processing.
        @param wedgeSizeX: volume size for x (size of original volume in real space)
        @param wedgeSizeY: volume size for y (size of original volume in real space)
        @param wedgeSizeZ: volume size for z (size of original volume in real space)
        @param humanUnderstandable: if True (default is False), the volume will \
            be transformed from reducedComplex to full and shifted afterwards
        @param rotation: rotation of 2nd wedge with respect to 1st
        @return: average of both wedges
        @rtype: L{pytom.lib.pytom_volume.vol}
        @author: Thomas Hrabe
        """

        w1 = self._wedge1.returnWedgeVolume(wedgeSizeX, wedgeSizeY, wedgeSizeZ, humanUnderstandable, rotation)
        w2 = self._wedge2.returnWedgeVolume(wedgeSizeX, wedgeSizeY, wedgeSizeZ, humanUnderstandable, rotation)

        # combine into one
        w = w1 + w2

        w[w>1] = 1



        return w

    def toSphericalFunc(self, b, radius=None):
        """Convert the wedge from real space to a spherical function.
        This function serves as an interface to FRM.
        @param b: Bandwidth of the spherical function.
        @param radius: Radius in frequency. Not used for DoubleTiltWedge.
        @return: a spherical function in numpy.array
        """
        raise NotImplementedError

    def fromXML(self, xmlObj):
        """
        @param xmlObj: A xml object
        @type xmlObj: L{lxml.etree._Element}
        """
        from lxml.etree import _Element

        if xmlObj.__class__ != _Element:
            raise TypeError('Is not a lxml.etree._Element! You must provide a valid XML-Wedge object.')

        try:
            wedges = xmlObj.xpath('SingleTiltWedge')

            if len(wedges) != 2:
                wedges = xmlObj.xpath('DoubleTiltWedge')[0].xpath('SingleTiltWedge')

            self._wedge1 = SingleTiltWedge()
            self._wedge2 = SingleTiltWedge()

            self._wedge1.fromXML(wedges[0])
            self._wedge2.fromXML(wedges[1])
        except:
            raise TypeError("Invalid DoubleTiltWedge format!")

    def toXML(self):
        from lxml import etree

        wedgeElement = etree.Element('DoubleTiltWedge')

        wedgeElement.append(self._wedge1.toXML())
        wedgeElement.append(self._wedge2.toXML())

        return wedgeElement

    def apply(self, volume, rotation=None):
        """
        apply: Applies this wedge to a given volume
        @param volume: The volume to be filtered
        @type volume: L{pytom.lib.pytom_volume.vol} or L{pytom.lib.pytom_volume.vol_comp}
        @param rotation: rotate the wedge
        @type rotation: L{pytom.agnostic.structures.Rotation}
        @return: The filtered volume
        @rtype: L{pytom.lib.pytom_volume.vol}
        @author: Thomas Hrabe
        """
        from pytom.agnostic.filter import applyFourierFilter

        if not volume.__class__ == xp.array:
            raise TypeError('You must provide a pytom.lib.pytom_volume.vol here!')

        wedgeVolume = self.returnWedgeVolume(volume.size_x(), volume.size_y(), volume.size_z(),
                                             humanUnderstandable=False, rotation=rotation)

        result = applyFourierFilter(volume, wedgeVolume)

        return result


class Wedge3dCTF(PyTomClass):
    """
    This is a class to allow pytom running with 3d ctf volumes from warp or relion as the missing wedge.
    Class handles reading of the wedge file from disk. Before being able to use the 3d ctf from warp it
    first needs to be converted to have the z axis as the reduced fourier space dimension.
    """

    def __init__(self, filename='', ctf_max_resolution=0.):
        """
        Class is mainly an io wrapper for the 3d ctf volumes so only needs a filename.
        @param filename: path to .mrc/.em file
        """
        self._filename = filename
        self._ctf_max_resolution = ctf_max_resolution

        # for FRM
        self._wedge_vol = None  # cache for storing the wedge volume
        self._bw = None
        self._sf = None

    def getFilename(self):
        return self._filename

    def setFilename(self, filename):
        self._filename = filename

    def get_ctf_max_resolution(self):
        return self._ctf_max_resolution

    def set_ctf_max_resolution(self, ctf_max_resolution):
        self._ctf_max_resolution = ctf_max_resolution

    def returnWedgeVolume(self, wedgeSizeX=None, wedgeSizeY=None, wedgeSizeZ=None, humanUnderstandable=False,
                          rotation=None):
        """parameters here are not needed but for compat with the wedge class"""
        if wedgeSizeX is not None or wedgeSizeY is not None or wedgeSizeZ is not None:
            print('WARNING! Wedge size specification wont have effect for 3d ctf wedge!')
        if rotation is not None:
            print('WARNING! Wedge rotation wont have effect for 3d ctf wedge!')
        wedge = read(self._filename)
        if not humanUnderstandable:
            return wedge
        else:
            return fftshift(fourier_reduced2full(wedge))

    def apply(self, volume, rotation=None):

        if rotation is not None:
            print('Warning: rotation will not be applied for 3d ctf wedge')

        if volume.__class__ is not xp.ndarray:
            raise TypeError('Wedge3dCTF: You must provide a numpy/cupy array here!')

        wedge = self.returnWedgeVolume()

        return applyFourierFilter(volume, wedge)

    def toSphericalFunc(self, b, radius=None):
        """Convert the wedge from real space to a spherical function.
        This function serves as an interface to FRM.
        @param b: Bandwidth of the spherical function.
        @param radius: Radius in frequency. Not used for SingleTiltWedge.
        @return: a spherical function in numpy.array
        """
        raise NotImplementedError

    def fromXML(self, xmlObj):
        """
        fromXML: Assigns values to result attributes from XML object
        @param xmlObj: A xml object
        @type xmlObj: L{lxml.etree._Element}
        @author: Thomas Hrabe
        """

        from lxml.etree import _Element

        if xmlObj.__class__ != _Element:
            raise TypeError('Is not a lxml.etree._Element! You must provide a valid XML-Wedge object.')

        self._filename = xmlObj.get('Filename')
        self._ctf_max_resolution = float(xmlObj.get('CtfMaxResolution'))

    def toXML(self):
        """
        toXML : Compiles a XML file from result object
        rtype : L{lxml.etree._Element}
        @author: Thomas Hrabe
        """
        from lxml import etree

        wedgeElement = etree.Element('Wedge3dCTF', Filename=self._filename,
                                     CtfMaxResolution=str(self._ctf_max_resolution))

        return wedgeElement


class GeneralWedge(PyTomClass):
    """A general weighting in the Fourier space. Maybe should not call it wedge anymore. Still use it for consistency.
    """

    def __init__(self, filename=None):
        """
        @param filename: name of wedge volume (.em or .mrc file)
        @type filename: L{str}
        """
        self._wedge_filename = filename
        self._weight_vol = None

    def _read_weight(self):
        if self._weight_vol is None:
            from pytom.agnostic.io import read
            w = read(self._wedge_filename)

            self._weight_vol = w

        return

    def returnWedgeFilter(self, wedgeSizeX=None, wedgeSizeY=None, wedgeSizeZ=None, rotation=None):
        from pytom.agnostic.structures import Weight as weight
        from pytom.agnostic.structures import Rotation

        self._read_weight()
        w = weight(self._weight_vol)  # if you initialize weight obj this way, it assumes the input volume is full and zero freq is at the center!

        if (not isinstance(rotation, Rotation)) and (not isinstance(rotation, RotationPytomC)):
            w.rotate(rotation.getZ1(), rotation.getZ2(), rotation.getX())

        return w

    def returnWedgeVolume(self, wedgeSizeX=None, wedgeSizeY=None, wedgeSizeZ=None, humanUnderstandable=False,
                          rotation=None):
        wedgeFilter = self.returnWedgeFilter(None, None, None, rotation)

        # return only the reduced complex part
        wedgeVolume = wedgeFilter.getWeightVolume(True)

        if humanUnderstandable:
            # applies hermitian symmetry to full and ftshift so that even humans do understand

            from pytom.agnostic.transform import fourier_reduced2full as reducedToFull
            wedgeVolume = reducedToFull(wedgeVolume, wedgeSizeX % 1)

            xp.fft.fftshift(wedgeVolume, True)

        return wedgeVolume

    def apply(self, volume, rotation=None):

        if not volume.__class__ == xp.array:
            raise TypeError('You must provide a xp.ndarray here!')

        from pytom.agnostic.filter import applyFourierFilter as filter
        wedgeFilter = self.returnWedgeFilter(None, None, None, rotation)
        result = filter(volume, wedgeFilter)
        return result

    def toSphericalFunc(self, b, radius):
        """Convert the wedge from real space to a spherical function.
        This function serves as an interface to FRM.
        @param b: Bandwidth of the spherical function.
        @param radius: Radius in frequency.
        @return: a spherical function in numpy.array
        """
        assert (b <= 128)

        # read the weight volume
        self._read_weight()

        size_x = self._weight_vol.size_x()
        size_y = self._weight_vol.size_y()
        size_z = self._weight_vol.size_z()

        # start sampling
        import numpy as np
        from math import pi, sin, cos
        from pytom.lib.pytom_numpy import vol2npy
        from scipy.ndimage.interpolation import map_coordinates
        v = vol2npy(self._weight_vol)

        x_ind = []
        y_ind = []
        z_ind = []
        for j in range(2 * b):
            for k in range(2 * b):
                the = pi * (2 * j + 1) / (4 * b)  # (0,pi)
                phi = pi * k / b  # [0,2*pi)

                x = cos(phi) * sin(the) * radius + size_x / 2
                y = sin(phi) * sin(the) * radius + size_y / 2
                z = cos(the) * radius + size_z / 2

                x_ind.append(x)
                y_ind.append(y)
                z_ind.append(z)

        # do interpolation
        inds = np.array([x_ind, y_ind, z_ind])
        res = map_coordinates(v, inds, order=3)

        # set the threshold
        res[res > 0.5] = 1
        res[res <= 0.5] = 0

        return xp.array(res)

    def fromXML(self, xmlObj):
        """
        @param xmlObj: A xml object
        @type xmlObj: L{lxml.etree._Element}
        """
        from lxml.etree import _Element

        if xmlObj.__class__ != _Element:
            raise TypeError('Is not a lxml.etree._Element! You must provide a valid XML-Wedge object.')

        if xmlObj.tag == 'GeneralWedge':
            element = xmlObj
        else:
            element = xmlObj.xpath('GeneralWedge')
            element = element[0]

        self._wedge_filename = element.get('Filename')
        self._weight_vol = None

        from pytom.tools.files import checkFileExists
        if not checkFileExists(self._wedge_filename):
            raise IOError('File ' + self._wedge_filename + ' does not exist!')

    def toXML(self):
        from lxml import etree
        wedgeElement = etree.Element('GeneralWedge', Filename=self._wedge_filename)
        return wedgeElement


class Particle(PyTomClass):
    """
    Particle: Stores information about individual particle.
    """

    def __init__(self, filename='', rotation=None, shift=None, wedge=None, className=0,
                 pickPosition=None, score=None, infoGUI=None):
        """
        @param filename: name of particle volume file (.em or .mrc file)
        @type filename: L{str}
        @param rotation: A pre rotation of the particle? (optional)
        @type rotation: L{pytom.agnostic.structures.Rotation}
        @param shift: A pre shift of particle? (optional)
        @type shift: L{pytom.agnostic.structures.Shift}
        @param wedge: Wedge information
        @type wedge: L{pytom.agnostic.structures.Wedge}
        @param className: Name of class this particle belongs to (numericals will be converted to string)
        @type className: int / float or string?
        @param pickPosition: Position of particle in original volume
        @type pickPosition: L{pytom.agnostic.structures.PickPosition}
        @param score: XCF Score of this particle, if any
        @type score: L{pytom.score.score.Score}
        """
        self._filename = filename
        from pytom.agnostic.structures import Rotation, Shift

        if not rotation:
            self._rotation = Rotation(0.0, 0.0, 0.0)
        elif rotation.__class__ == list:
            self._rotation = Rotation(rotation[0], rotation[1], rotation[2])
        elif rotation.__class__ == Rotation:
            self._rotation = rotation
        else:
            raise TypeError('Unknown type for rotation parameter!')

        # if not shift:
        if shift == None:
            self._shift = Shift(0.0, 0.0, 0.0)
        elif shift.__class__ == list:
            self._shift = Shift(shift[0], shift[1], shift[2])
        elif shift.__class__ == Shift:
            self._shift = shift
        else:
            raise TypeError('Unknown type for shift parameter!')

        if not infoGUI:
            self._infoGUI = InformationGUI('')
        elif infoGUI.__class__ == list:
            self._infoGUI = InformationGUI(infoGUI[0])
        elif infoGUI.__class__ == InformationGUI:
            self._infoGUI = infoGUI
        else:
            raise TypeError('Unknown type for InformationGUI parameter')

        if not pickPosition:
            from pytom.agnostic.structures import PickPosition
            self._pickPosition = PickPosition()
        else:
            self._pickPosition = pickPosition

        if wedge == []:
            wedge = None

        self._wedge = wedge or Wedge()

        self._className = str(className) or None

        self._score = score
        self._xmlCounter = 0
        self._scoreValue = -1.

    def getRotation(self):
        return self._rotation

    def getShift(self):
        return self._shift

    def getWedge(self):
        """
        get wedge object of particle
        @return: wedge
        @rtype: L{pytom.agnostic.structures.Wedge}
        """
        return self._wedge

    def getFilename(self):
        return self._filename

    def getClassName(self):
        """
        @deprecated: Use getClass instead
        """
        return self.getClass()

    def getClass(self):
        """
        @return: class of particle
        @rtype: str
        """
        if self._className != None:
            return str(self._className)
        else:
            return None

    def getPickPosition(self):
        return self._pickPosition

    def getScore(self):
        return self._score

    def getInfoGUI(self):
        return self._infoGUI

    def getVolume(self, binning=1):
        """
        read Volume from disk. If specified volume will be resized in Fourier space
        @param binning: binning factor (e.g., 2 makes it 2 times smaller in each dim)
        @type binning: int or float
        @return: volume
        @rtype: L{pytom.lib.pytom_volume.vol}
        @change: FF
        """
        from pytom.agnostic.io import read
        from pytom.tools.files import checkFileExists
        from pytom.agnostic.transform import resize

        if not checkFileExists(self._filename):
            raise IOError('Particle ' + self._filename + ' does not exist!')

        try:
            # volume = read(self._filename, 0,0,0,0,0,0,0,0,0, binning, binning, binning)
            volume = read(self._filename)#, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1)
            if binning != 1:
                volume, volumef = resize(volume=volume, factor=1. / binning, interpolation='Fourier')
        except RuntimeError:
            raise RuntimeError('Error reading file ' + self._filename)

        return volume

    def setPickPosition(self, pickPosition):
        """
        setPickPosition
	    @param pickPosition: position of particle in tomogram
	    @type pickPosition: L{PickPosition}
        """

        if pickPosition.__class__ == list and len(pickPosition) == 3:
            pickPosition = PickPosition(pickPosition[0], pickPosition[1], pickPosition[2])
        elif pickPosition.__class__ != PickPosition:
            raise TypeError('You must provide a PickPosition object here!')
        else:
            pass

        self._pickPosition = pickPosition

    def setShift(self, shift):
        """
        set shift of particle
        @param shift: Shift
        @type shift: L{pytom.agnostic.structures.Shift}
        """
        if not shift.__class__ == Shift:
            raise TypeError('Shift argument must be of class shift!')

        self._shift = shift

    def setFilename(self, newFilename=None):
        """
        set filename of particle
        @param newFilename: name of file
        @type newFilename: L{str}
        """
        self._filename = newFilename or self._filename

    def setWedge(self, wedge):
        """
        set wedge in particle
        @param wedge: Wedge
        @type wedge: L{pytom.agnostic.structures.Wedge}
        """
        self._wedge = wedge

    def setRotation(self, rotation):
        from pytom.agnostic.structures import Rotation

        if not rotation.__class__ == Rotation:
            raise TypeError('rotation parameter must be of type Rotation!')

        self._rotation = Rotation(rotation)

    def setScore(self, score):
        """
	    set Score type of particle
	    @param score: score type
	    @type score: L{pytom.score.score.Score}
        """
        from pytom.basic.score import Score

        if ((not score.__class__ == Score) and (not score.__class__.__bases__[0] == Score)):
            raise TypeError('score parameter must be of type Score')

        self._score = score

    def setScoreValue(self, scoreValue):
        """
        set score value (correlation coefficient) of a particle
        @param scoreValue: value of score
        @type scoreValue: L{float}
        """
        if not hasattr(self, 'self._score'):
            from pytom.basic.score import Score
            # set dummy score for particles
            dumscore = Score()
            dumscore.ctor()
            self.setScore(score=dumscore)
        self._score.setValue(value=scoreValue)
        self._scoreValue = scoreValue

    def getScoreValue(self):
        """
        get score value
        @return: scoreValue
        @rtype: L{float}
        """
        if not hasattr(self, 'self._score'):
            return self._scoreValue
        else:
            return self._score.getValue()

    def setClass(self, className):
        """
	    set class of particle
	    @param className: Name of class
	    @type className: L{str}

	    """
        self._className = str(className)

    def toXML(self):
        """
        toXML : Compiles a XML file from result object
        rtype : L{lxml.etree._Element}
        @author: Thomas Hrabe
        """
        from lxml import etree

        particle_element = etree.Element('Particle', Filename=self._filename)
        particle_element.append(self._rotation.toXML())
        particle_element.append(self._shift.toXML())
        particle_element.append(self._pickPosition.toXML())
        particle_element.append(self._wedge.toXML())
        particle_element.append(self._infoGUI.toXML())

        if self._score != None:
            particle_element.append(self._score.toXML())

        if not self._className == None:
            class_element = etree.Element('Class', Name=str(self._className))
            particle_element.append(class_element)

        return particle_element

    def fromXML(self, xmlObj):
        """
        @param xmlObj: A xml object
        @type xmlObj: L{lxml.etree._Element}
        @author: Thomas Hrabe
        """

        from lxml.etree import _Element
        from pytom.agnostic.structures import Rotation, Shift, PickPosition, Wedge

        if xmlObj.__class__ != _Element:
            from pytom.basic.exceptions import ParameterError
            raise ParameterError('Is not a lxml.etree._Element! You must provide a valid XML-Particle object.')

        if xmlObj.tag == 'Particle':
            particle_element = xmlObj
        else:
            particle_element = xmlObj.xpath('Particle')
            particle_element = particle_element[0]

        self._filename = particle_element.get('Filename')

        rotation_element = particle_element.xpath('Rotation')
        if len(rotation_element) > 0:
            rotation_element = rotation_element[0]
            self._rotation = Rotation([0.0, 0.0, 0.0])
            self._rotation.fromXML(rotation_element)
        else:
            self._rotation = Rotation([0.0, 0.0, 0.0])

        shift_element = particle_element.xpath('Shift')
        if len(shift_element) > 0:
            shift_element = shift_element[0]
            self._shift = Shift([0.0, 0.0, 0.0])
            self._shift.fromXML(shift_element)
        else:
            self._shift = Shift([0.0, 0.0, 0.0])

        source_element = particle_element.xpath('InfoGUI')
        if len(source_element) > 0:
            source_element = source_element[0]
            self._infoGUI = InformationGUI([''])
            self._infoGUI.fromXML(source_element)
        else:
            self._infoGUI = InformationGUI([''])

        position_element = particle_element.xpath('PickPosition')
        if len(position_element) > 0:
            position_element = position_element[0]
            self._pickPosition = PickPosition()
            self._pickPosition.fromXML(position_element)
        else:
            self._pickPosition = PickPosition()

        wedgeXML = particle_element.xpath('Wedge')
        if len(wedgeXML) == 0:
            wedgeXML = particle_element.xpath('SingleTiltWedge')
            if len(wedgeXML) == 0:
                wedgeXML = particle_element.xpath('WedgeInfo')
            if len(wedgeXML) == 0:
                wedgeXML = particle_element.xpath('DoubleTiltWedge')
            if len(wedgeXML) == 0:
                wedgeXML = particle_element.xpath('GeneralWedge')
        self._wedge = Wedge()
        if len(wedgeXML) > 0:
            self._wedge.fromXML(wedgeXML[0])

        class_element = particle_element.xpath('Class')
        if len(class_element) > 0 and class_element[0] != 'None':
            class_element = class_element[0]
            self._className = str(class_element.get('Name'))
        else:
            self._className = None

        score_element = particle_element.xpath('Score')
        if len(score_element) > 0:
            from pytom.basic.score import fromXML
            self._score = fromXML(score_element[0])

    def copy(self):
        """
        copy: Copies self to new Particle object
        """
        p = Particle('')
        pXML = self.toXML()

        p.fromXML(pXML)
        return p

    def check(self):

        from pytom.tools.files import checkFileExists

        if not checkFileExists(self._filename):
            raise IOError('Can not find particle at path: ' + str(self._filename))

    def getTransformedVolume(self, binning=1):
        """
        getTransformedVolume: Returns particle volume with applied inverse rotation and inverse shift
        @param binning: binning factor
        @type binning: C{int}
        @rtype: L{pytom.lib.pytom_volume.vol}
        """
        from pytom.agnostic.structures import Shift, Rotation
        from pytom.voltools import transform

        volume = self.getVolume(binning)

        shift = self.getShift()
        rotation = self.getRotation()

        if rotation != Rotation(0, 0, 0) or shift != Shift(0, 0, 0):

            volumeTransformed = xp.zeros((volume.shape[0], volume.shape[1], volume.shape[2]),dtype=xp.float32)

            transform(volume, output=volumeTransformed, rotation=[-rotation[1], -rotation[0], -rotation[2]],
                      center=[int(volume.shape[0] // 2), int(volume.shape[0] // 2), int(volume.shape[0] // 2)],
                      translation=[-shift[0] / binning, -shift[1] / binning, -shift[2] / binning], device=device,
                      interpolation='filt_bspline')

            return volumeTransformed

        else:
            return volume

    def setPickPositionCenter(self, x, y, z):
        """
        setPickPositionCenter: Assuming the pick position is related to the upper left pixel it will transform \
        the pick position coordinates to the new center.
        """
        self._pickPosition + [x, y, z]

    def __cmp__(self, otherParticle):

        if not otherParticle.__class__ == Particle:
            raise TypeError('Can only compare a Particle to another Particle')

        return cmp(self.getScore(), otherParticle.getScore())


def keyFunctionForParticleSortingByScore(particle):
    return particle.getScore()


def keyFunctionForParticleSortingByClassLabel(particle):
    return str(particle.getClass())


class ParticleList(PyTomClass):
    """
    ParticleList: Stores a list of particle objects L{pytom.agnostic.structures.Particle}
    """

    def __init__(self, directory=None, pl=None):
        """
        __init__:
        @param directory: Source directory of particle files
        @type directory: C{str}
        @param pl: A particle list.
        @type pl: L{pytom.agnostic.structures.ParticleList}
        """
        if directory and directory.__class__ == str and len(directory) > 0:
            if directory[len(directory) - 1] == '/':
                self._directory = directory
            else:
                self._directory = directory + '/'
        else:
            self._directory = None
        self._particleList = pl or []
        self._XMLfilename = ''

    def getDirectory(self):
        return self._directory

    def setDirectory(self, directory):
        """
        @param directory: Source directory of particle files
        @type directory: L{str}
        """
        from pytom.tools.files import checkDirExists
        if not checkDirExists(directory):
            raise RuntimeError('ParticleList.setDirectory: directory does not exist!')
        self._directory = directory

    def loadDirectory(self, av3Indexing=False, av3Prefix=''):
        """
        readDirectory: Appends all em/mrc/ccp4 files in self._directory to self._particleList (as unaligned particles).
        Make sure you do not have any other em,mrc,ccp4 data in that folder as it will also be considered as particles.
        @param av3Indexing: use running index when files were saved with av3. If you use this branch, make sure all!!! \
        em files in the directory have the same prefix!
        @param av3Prefix: Set file prefix, including _ !
        """
        import dircache
        from pytom.agnostic.structures import Particle
        from pytom.tools.files import checkDirExists

        if not checkDirExists(self._directory):
            raise IOError('ParticleList.loadDirectory: Directory ' + self._directory + ' does not exist!')

        pl = dircache.listdir(self._directory)
        numberParticles = 0

        for pli in pl:
            if pli.find('.em') > 0:
                if not av3Indexing:
                    self._particleList.append(Particle(self._directory + pli))
                numberParticles = numberParticles + 1
            if pli.find('.mrc') > 0:
                if not av3Indexing:
                    self._particleList.append(Particle(self._directory + pli))
                numberParticles = numberParticles + 1
            if pli.find('.ccp4') > 0:
                if not av3Indexing:
                    self._particleList.append(Particle(self._directory + pli))
                numberParticles = numberParticles + 1
        if av3Indexing:
            for i in range(numberParticles):
                self._particleList.append(Particle(self._directory + av3Prefix + str(i + 1) + '.em'));

    def __len__(self):
        return len(self._particleList)

    def append(self, particle):
        """
        append: Appends particle to self.particleList
        @param particle: particle
        @type particle: L{pytom.agnostic.structures.Particle}
        """
        from pytom.agnostic.structures import Particle

        assert particle.__class__ == Particle
        self._particleList.append(particle)

    def toXML(self):
        """
        toXML : Compiles a XML file from result object
        rtype : L{lxml.etree._Element}
        @author: Thomas Hrabe
        """
        from lxml import etree

        directory_element = etree.Element('ParticleList')

        if len(self._particleList) > 0:
            for p in self._particleList:
                directory_element.append(p.toXML())

        return directory_element

    def toShortXML(self):
        """
        write to short XML containing only filename of particle list rather than all particles
        """
        from lxml import etree
        if type(self._directory) != str:
            particleList_element = etree.Element('ParticleList', Path='./')
        else:
            particleList_element = etree.Element('ParticleList', Path=self._directory)
        return particleList_element

    def toHTMLFile(self, filename):
        """
        toHTMLFile: Overrides parent method and stores ParticleList to HMTL
        @param filename: HTML filename
        @type filename: C{str}
        """
        from pytom.tools.files import getPytomPath

        super(self.__class__, self).toHTMLFile(filename, getPytomPath() + '/xslt/ParticleList.xsl')

    def fromXML(self, xmlObj):
        """
        @param xmlObj: A xml object
        @type xmlObj: L{lxml.etree._Element}
        @author: Thomas Hrabe
        """
        from lxml.etree import _Element

        if xmlObj.__class__ != _Element:
            from pytom.basic.exceptions import ParameterError
            raise ParameterError('Is not a lxml.etree._Element! You must provide a valid XML object.')

        if xmlObj.tag == 'ParticleList':
            directory_element = xmlObj
        else:
            directory_element = xmlObj.xpath('ParticleList')
            directory_element = directory_element[0]

        self._particleList = []
        particles = directory_element.xpath('Particle')

        if len(particles) > 0:
            from pytom.agnostic.structures import Particle
            for p in particles:
                pp = Particle('')
                pp.fromXML(p)
                self._particleList.append(pp)

    def __getitem__(self, key):
        """
        __getitem__: Retreive particle at position definded by key
        @param key:
        @type key: int or string
        @rtype: L{pytom.basic.Particle} if key is a int or string. L{pytom.basic.ParticleList} if key is a slice initialised with [start:stop]. Step is ignored completely
        """
        if isinstance(key, int):

            if key < len(self):
                return self._particleList[key]
            else:
                raise IndexError('Index out of range.')

        elif key.__class__ == str:
            return self.getParticleByFilename(key)

        elif key.__class__ == slice:

            start = key.start
            stop = key.stop
            step = key.step

            if not start:
                start = 0

            if not stop:
                stop = len(self)

            if not step:
                step = 1

            if stop >= 922337203685477580:
                stop = len(self._particleList)

            newParticleList = ParticleList(self.getDirectory())

            for i in range(start, stop, step):
                newParticleList.append(self._particleList[i])

            return newParticleList

        else:
            assert False

    def __setitem__(self, key, value):
        """
        __setitem__: Enables syntax such as particleList[key] = value
        @param key: The index
        @type key: int
        @param value: A particle
        @type value: L{pytom.agnostic.structures.Particle}
        """
        if isinstance(key, int):
            if key < len(self._particleList):

                from pytom.agnostic.structures import Particle

                if value.__class__ == Particle:
                    self._particleList[key] = value
                else:
                    raise TypeError('Can only append Particles to ParticleList')
            else:
                raise IndexError('Index out of range.')
        elif key.__class__ == slice:

            start = key.start
            stop = key.stop
            step = key.step

            if not start:
                start = 0

            if not stop:
                stop = len(self)

            if not step:
                step = 1

            if not value.__class__ == ParticleList:
                raise TypeError('You must provide a ParticleList here!')

            if not ((stop - start) == len(value)):
                raise RuntimeError('Range and ParticleList length must match!')

            for i in range(start, stop, step):
                self[i] = value[i - start]

        else:
            assert False

    def fromAlignmentList(self, filename):
        """
        @param filename: filename of 'AlignmentList' (whatever this is ...)
        @type filename: L{str}
        """
        from pytom.tools.files import checkFileExists
        from pytom.alignment.structures import AlignmentList

        if not checkFileExists(filename):
            raise RuntimeError('AlignmentList ' + filename + ' does not exist!')
        al = AlignmentList()
        al.fromXMLFile(filename)
        pl = al.toParticleList()
        for p in pl:
            self.append(p)

    def fromMOTL(self, motlFile, sourceDirectory, particlePrefix, originVolume='', particleStartOffset=0, score=None):
        """
        fromMOTL: Loads a MOTL.em into this ParticleList
        @param motlFile: Path to MOTL.em
        @type motlFile: C{str}
        @param sourceDirectory: Directory where all particles are stored
        @param particlePrefix: Prefix of each particle in sourceDirectory (prefix == strings before the index, including _ )
        @param originVolume: Name / Path to volume where particles come from ('' by default)
        @param particleStartOffset: Offset in particle indexing. Motls generated with MATLAB will most likely start with a 1 as first particle (not with 0 as in C or Python)
        @type particleStartOffset: C{int}
        """
        from pytom.agnostic.io import read
        from pytom.tools.files import checkFileExists

        from pytom.agnostic.structures import Rotation, Shift, PickPosition, Particle
        # if not score:
        if score == None:
            from pytom.basic.score import xcfScore as score

        checkFileExists(motlFile)

        motifList = read(motlFile)

        numberParticles = motifList.size_y()

        for i in range(numberParticles):

            filename = sourceDirectory + particlePrefix + str(int(motifList(3, i, 0)) + particleStartOffset) + '.em'

            if not checkFileExists(filename):
                print('File ', filename, 'does not exist! ', i, ' (Matlab index)')

            currentClass = str(motifList(19, i, 0))

            origin = PickPosition(motifList(7, i, 0) - 1, motifList(8, i, 0) - 1, motifList(9, i, 0) - 1, originVolume)

            shift = Shift(motifList(10, i, 0), motifList(11, i, 0), motifList(12, i, 0))

            rotation = Rotation(motifList(16, i, 0), motifList(17, i, 0), motifList(18, i, 0))

            s = score()
            s.setValue(motifList(0, i, 0))

            particle = Particle(filename, rotation=rotation, shift=shift, className=currentClass, pickPosition=origin,
                                score=s)

            self.append(particle)

    def toMOTL(self, filename):
        """
        @warning: Resulting MOTL - Particle Index (4th entry) is increment 1:numberParticles. Does not neccessarily match to the particleFilename! Check L{copyFiles} for that.
        @param filename: The filename
        """
        from pytom.lib.pytom_volume import vol

        numberParticles = len(self)

        motl = vol(20, numberParticles, 1)
        motl.setAll(0)

        for i, particle in enumerate(self._particleList):

            # set the CCC
            score = particle.getScore()
            if score != None:
                motl(float(score.getValue()), 0, i, 0)

            origin = particle.getPickPosition()

            motl(i, 3, i, 0)

            motl(origin.getX() + 1, 7, i, 0)
            motl(origin.getY() + 1, 8, i, 0)
            motl(origin.getZ() + 1, 9, i, 0)

            shift = particle.getShift()

            motl(shift.getX(), 10, i, 0)
            motl(shift.getY(), 11, i, 0)
            motl(shift.getZ(), 12, i, 0)

            rotation = particle.getRotation()

            motl(rotation.getZ1(), 16, i, 0)
            motl(rotation.getZ2(), 17, i, 0)
            motl(rotation.getX(), 18, i, 0)

            if particle.getClassName() and particle.getClassName() != 'None':
                motl(float(particle.getClassName()), 19, i, 0)

        motl.write(filename)

    def setClassAllParticles(self, className):
        """
        setClassAllParticles: Updates class membership of all particles
        @param className: The new class name
        """

        for i, particle in enumerate(self._particleList):
            particle.setClass(className)
            self._particleList[i] = particle

    def getClassList(self):
        """
        Get the class list of all the particles.
        """
        res = []
        for p in self._particleList:
            res.append(p.getClass())

        return res

    def particlesFromClass(self, className):
        """
        particlesFromClass: Selects all members of className
        @param className: Name of particle class
        @type className: str
        @return: ParticleList
        """
        from pytom.agnostic.structures import ParticleList as ParticleListConstructor

        query = '/ParticleList/Particle[Class/@Name="' + str(className) + '"]'

        particles = self.xpath(query)

        newParticleList = ParticleListConstructor(self._directory, [])

        for particle in particles:
            p = Particle('')
            p.fromXML(particle)
            newParticleList.append(p)

        return newParticleList

    def splitByClass(self, verbose=False):
        """
        splitByClass: Splits current self.particleList and returns N ParticleLists if there are N classes.
        This function will resort the current particle list.
        @return: List of L{pytom.agnostic.structures.ParticleList}
        """
        # sort list to ascending order
        self.sortByClassLabel()

        classes = self.xpath('/ParticleList/Particle/Class/@Name')

        classNames = []

        for i, aClass in enumerate(classes):
            if not str(aClass) in classNames:
                classNames.append(str(aClass))

        if verbose:
            print(classNames)

        particleListList = []

        for className in classNames:

            particleList = self.particlesFromClass(className)

            if verbose:
                print(particleList)

            particleListList.append(particleList)

        return particleListList

    def getParticlesByString(self, string):
        """
        getParticlesByString: Will get all particles whose filename property contain string
        @param string: The substring searched
        @type string: str
        @return: A new particle list
        """

        pxml = self.xpath('/ParticleList/Particle[contains(@Filename,"' + string + '")]')

        if len(pxml) == 0:
            raise RuntimeError('Particles with property ' + string + ' not found in this list!')

        plnew = ParticleList()

        for px in pxml:
            p = Particle()
            p.fromXML(px)
            plnew.append(p)

        return plnew

    def getParticleByFilename(self, filename):
        """
        getParticleByFilename:
        @param filename: The filename
        @returns: A particle if successfull or
        """

        for p in self._particleList:
            if p.getFilename() == filename:
                return p

        raise KeyError('Particle not found in list!')

    def __add__(self, particleList):
        """
        __add__ : Concatenates two ParticleLists
        @param particleList:
        """

        from pytom.agnostic.structures import ParticleList

        if not particleList.__class__ == ParticleList:
            raise TypeError('ParticleList: Can not concatenate this particleList to a non ParticleList object!')

        for i in range(len(particleList)):
            self._particleList.append(particleList[i])

        return self

    def setWedgeAllParticles(self, wedge):
        """
        setWedgeAllParticles: Sets wedge properties to all particles
        @param wedge: Wedge
        @type wedge: L{pytom.agnostic.structures.Wedge}
        """
        for i in range(len(self)):
            p = self._particleList[i]
            p.setWedge(wedge)
            self._particleList[i] = p

    def copy(self):
        pl = ParticleList()
        plXML = self.toXML()
        pl.fromXML(plXML)
        return pl

    def copyFiles(self, destinationDirectory, particlePrefix):
        """
        copyFiles: Copies particles in this list to a directory and renames each particle to particlePrefix_particleIndexInList.em. The index is designed to be consistent with L{toMOTL}
        @param destinationDirectory: The directory
        @type destinationDirectory: str
        @param particlePrefix: Prefix of particle name
        @type particlePrefix: str
        """
        from pytom.tools.files import checkDirExists

        if not checkDirExists(destinationDirectory):
            raise IOError('Destination directory does not exist.' + destinationDirectory)

        for i, particle in enumerate(self._particleList):
            particleVolume = particle.getVolume()
            particleVolume.write(destinationDirectory + particlePrefix + '_' + str(i + 1) + '.em')

    def sortParticlesByIndex(self, particlePrefix=None, directoryName=None, startIndex=1, fileType='.em'):
        """
        sortParticlesByIndex: Will sort this particleList so that particles are \
        in av3 index order. Must be processed before this is exported as motl

        @param particlePrefix: ParticlePrefix
        @param directoryName: Source directory where files are stored
        @param startIndex: starting index of the particle, default is 1
        @param fileType: File type of extracted files [.em , .mrc, .ccp4]. .em is default. Obsolete.
        @return: A new particle list
        """
        import os

        if not particlePrefix:
            print('Warning: Particle Prefix is empty. Will return original particle list.')
            return self

        if directoryName.__class__ != str:
            directoryName = self._directory

        if directoryName[-1] != os.path.sep:
            directoryName += os.path.sep

        numberParticles = len(self)

        newParticleList = ParticleList(self._directory)

        for index in range(startIndex, numberParticles + startIndex):
            try:
                particle = self.getParticleByFilename(directoryName + particlePrefix + str(index) + '.em')
            except:
                continue

            newParticleList.append(particle)

        return newParticleList

    def check(self):
        """
        check: Performs logical check on self to make sure that job settings are correct. Will display error message if not.
        """
        for particle in self._particleList:
            particle.check()

    def setClassFromLocalizationClassList(self, filename):
        from pytom.tools.files import readStringFile
        lines = readStringFile(filename)
        class_labels = [int(i) for i in lines.split('\n')[:-1]]

        for index in range(len(self)):
            particle = self[index]
            particle.setClass(class_labels[index])
            self[index] = particle

    def exportClassList(self, filename):
        try:
            f = open(filename, 'w')
            for p in self._particleList:
                f.write(str(p.getClass()) + '\n')
        finally:
            f.close()

    def setPickPositionFromParticleList(self, pl):
        assert len(self._particleList) == len(pl)
        for i in range(len(pl)):
            self._particleList[i].setPickPosition(pl[i].getPickPosition())

    def fromLocalizationResult(self, filename):
        from pytom.localization.structures import readParticleFile
        particles = readParticleFile(filename)

        for i in range(len(particles)):
            particle = particles[i]
            p = particle.toParticle()

            self.append(p)

    def averageScore(self):
        """
        @return: average score value
        @return: average score of particles stored in list
        @rtype: L{float}
        """
        scores = 0.0
        for particle in self:
            score = particle.getScore()
            if score.__class__ == None:
                print('Warning: score was not set for particleList.averageScore!')
            else:
                scores = scores + score.getValue()
        return scores / float(len(self))

    def sumOfScores(self):
        """
        @deprecated: replaced by averageScore
        """
        return self.averageScore()

    def scoresAsVector(self):
        """
        list of scores
        @return: scores
        @rtype: L{list}
        """
        scores = []
        for particle in self:
            score = particle.getScore()
            scores.append(score.getValue())
        return scores

    def classMoves(self, newParticleList, verbose=False):
        """
        classMoves:
        @deprecated: use classDifference instead!
        """
        return self.classDifference(newParticleList)

    def classDifference(self, newParticleList, verbose=False):
        """
        classDifference: Determines the difference of this (Ground Truth) particleList compared to the newParticleList
        @return: Values determined by L{pytom.classification.classificationResults.assessClassification}
        """
        from pytom.classification.classificationResults import assessClassification

        return assessClassification(newParticleList, self, verbose)

    def setCenterInPickVolume(self, x, y, z):
        """
        setCenterInPickVolume: Assuming the pick position is related to the \
        upper left pixel it will transform the pick position coordinates to \
        the new center. Transforms all particles in list.

        @param x: x-coordinate of center
        @type x: L{float}
        @param y: y-coordinate
        @type y: L{float}
        @param z: z-coordinate
        @type z: L{float}

        @lastchange: 10.12.2012: changed name of function, FF
        """

        for i in range(len(self)):
            particle = self[i]
            particle.setPickPositionCenter(x, y, z)
            self[i] = particle

    def average(self, averageFileName, progressBar=False, createInfoVolumes=False,
                _mpiParallel=False, weighting=False, gpuIDs=[]):
        """
        average: Calculates average of ParticleList, stored as file
        @param averageFileName: Filename of average result such as average.em
        @param progressBar: Print ProgressBar during averaging. Default is False!
        @param createInfoVolumes: Create info data (wedge sum, inverted density) too? False by default.
        @param _mpiParallel: run the average in parallel
        @param weighting: create the average using weighting scheme
        @return: Reference object for resulting average
        @rtype: L{pytom.agnostic.structures.Reference}
        """

        gpuIDs = [] if gpuIDs is None else gpuIDs
        if _mpiParallel:
            from pytom.alignment.alignmentFunctions import _disrtibuteAverageMPI

            return _disrtibuteAverageMPI(self, averageFileName, progressBar, False, createInfoVolumes)
        elif gpuIDs:
            from pytom.alignment.alignmentFunctions import averageGPU as average
            return average(self, averageFileName, progressBar, False, createInfoVolumes, weighting)
        else:
            from pytom.alignment.alignmentFunctions import average

            return average(self, averageFileName, progressBar, False, createInfoVolumes, weighting)

    def averageClasses(self, directory, progressBar=False):
        """
        averageClasses: Will create averages of all particle classes
        @param directory: Directory where averages will be stored into
        @type directory: str
        @param progressBar: Print ProgressBar during averaging. Default is False!
        @type progressBar: bool
        """
        from pytom.tools.files import checkDirExists
        import pytom.lib.pytom_mpi as pytom_mpi

        if not checkDirExists(directory):
            raise RuntimeError('Directory specified does not exist!')

        pls = self.splitByClass()

        for i in range(len(pls)):
            particle = pls[i][0]
            className = particle.getClassName()
            pls[i].average(directory + '/class_' + str(className) + '.em', progressBar, False,
                           pytom_mpi.isInitialised())

    def splitOddEven(self, verbose=False):
        """
        return particle lists for odd and even indices

        @param verbose: Verbose mode. Default -> False
        @type verbose: L{bool}
        @return: even, odd
        @rtype: L{ParticleList}, L{ParticleList}
        @author: FF
        """
        even = ParticleList('/')
        odd = ParticleList('/')
        if verbose:
            print('Number of particles in particle list: ', len(self))
        for particleCounter in range(len(self)):
            if particleCounter % 2 == 0:
                even.append(self[particleCounter])
            else:
                odd.append(self[particleCounter])
        return odd, even

    def updateFromOddEven(self, odd, even, verbose=False):
        """
        update particle list from Odd and Even sublists
        @type odd: L{ParticleList}
        @type even: L{ParticleList}
        """
        assert type(odd) == ParticleList
        assert type(even) == ParticleList
        assert len(even) + len(odd) == len(self)
        if verbose:
            print('Number of particles in particle list: ', len(self), ", odd:", len(odd), ", even:", len(even))
        iodd = 0
        ieven = 0
        for ii in range(len(self)):
            if ii % 2 == 0:
                self._particleList[ii] = even[ieven]
                ieven = ieven + 1
            else:
                self._particleList[ii] = odd[iodd]
                iodd = iodd + 1

    def determine_resolution(self, criterion=0.5, number_bands=None, mask=None, verbose=False, plot='',
                            keepHalfsetAverages=False, halfsetPrefix='', parallel=True, randomize=0.8):
        """
        determine_resolution
        @param criterion: The resolution criterion
        @param number_bands: Will use cubesize_x / 2 as default if not specified
        @param mask: A mask used for specifying location of particle. Can be None
        @param verbose: Verbose mode. Default -> False
        @param plot: Plot FSC curve to disk? Provide svg or png filename here. Default is '' -> no plot!
        @param keepHalfsetAverages: Delete even / odd averages. Default is false -> averages will not be kept
        @param halfsetPrefix: Prefix for half set files. Default is ''
        @param parallel: If True (default), this function will enable parallel averaging if possible.
        @param randomize: if you want to correct your resolution using phases randomization set this variable to a value between 0 and 1.
        @type randomize: float
        @return: [Resolution in Nyquist , resolution in band, number_bands]
        @todo: Change return type to L{pytom.agnostic.structures.resolution} to make it clearer.
        """
        if len(self) < 2:
            raise RuntimeError('ParticleList must have at least 2 elements to determine resolution!')

        from pytom.agnostic.io import read
        from pytom.basic.correlation import fsc, determine_resolution
        import pytom.lib.pytom_mpi as pytom_mpi
        from pytom.lib.pytom_numpy import vol2npy

        import os

        even = ParticleList('/')
        odd = ParticleList('/')

        if verbose:
            print('Number of particles in particle list: ', len(self))

        for particleCounter in range(len(self)):

            if particleCounter % 2 == 0:
                even.append(self[particleCounter])
            else:
                odd.append(self[particleCounter])

        if verbose:
            print('Averaging even:')

        even.average(halfsetPrefix + 'even.em', verbose, _mpiParallel=(pytom_mpi.isInitialised() and parallel))

        if verbose:
            print('Averaging odd:')

        odd.average(halfsetPrefix + 'odd.em', verbose, _mpiParallel=(pytom_mpi.isInitialised() and parallel))

        oddVolume = read(halfsetPrefix + 'odd.em')
        evenVolume = read(halfsetPrefix + 'even.em')

        if not keepHalfsetAverages:
            os.system('rm ' + halfsetPrefix + 'odd.em')
            os.system('rm ' + halfsetPrefix + 'odd-PreWedge.em')
            os.system('rm ' + halfsetPrefix + 'odd-WedgeSumUnscaled.em')
            os.system('rm ' + halfsetPrefix + 'even.em')
            os.system('rm ' + halfsetPrefix + 'even-PreWedge.em')
            os.system('rm ' + halfsetPrefix + 'even-WedgeSumUnscaled.em')

        if not number_bands:
            number_bands = oddVolume.size_x() / 2

        if verbose:
            print('Using ', number_bands, ' shells for FSC')

        # oddVolume = vol2npy(oddVolume).copy()
        # evenVolume = vol2npy(evenVolume).copy()

        calc_fsc = fsc(oddVolume, evenVolume, number_bands, mask, verbose)
        #TODO: After this the variable f seems to appear (maybe original name of (calc_)fsc)
        #      This should always have errored out, is this dead code?
        if randomize is None:
            for (ii, fscel) in enumerate(f):
                f[ii] = 2. * fscel / (1. + fscel)
            r = determine_resolution(f, fscCriterion, verbose)
        else:
            randomizationFrequency = np.floor(determine_resolution(np.array(f), float(randomize), verbose)[1])
            oddVolumeRandomizedPhase = correlation.randomize_phase_beyond_freq(vol2npy(v1), randomizationFrequency)
            evenVolumeRandomizedPhase = correlation.randomize_phase_beyond_freq(vol2npy(v2), randomizationFrequency)
            write('randOdd.mrc', oddVolumeRandomizedPhase)
            write('randEven.mrc', evenVolumeRandomizedPhase)
            oddVolumeRandomizedPhase = read('randOdd.mrc')
            evenVolumeRandomizedPhase = read('randEven.mrc')
            fsc2 = fsc(oddVolumeRandomizedPhase, evenVolumeRandomizedPhase, number_bands, mask, verbose)
            fsc_true = list(correlation.calc_fsc_true(np.array(f), np.array(fsc2)))
            for (ii, fscel) in enumerate(fsc_true):
                fsc_true[ii] = 2. * fscel / (1. + fscel)
            r = determine_resolution(fsc_true, fscCriterion, verbose)
        # randomizationFrequency = np.floor(determine_resolution(fsc, 0.8, verbose)[1])

        # oddVolumeRandomizedPhase = randomize_phase_beyond_freq(oddVolume, randomizationFrequency)
        # evenVolumeRandomizedPhase = randomize_phase_beyond_freq(oddVolume, randomizationFrequency)
        # fsc2 = fsc(oddVolumeRandomizedPhase, evenVolumeRandomizedPhase, number_bands, mask, verbose)
        if verbose:
            print('FSC list:')
            print(calc_fsc)
            print('FSC_Random:\n', fsc2)
            print('FSC_true:\n', fsc_true)

        if not plot == '':
            try:
                from pytom.plotting.plot import plotFSC
                plotFSC(calc_fsc, plot)
            except:
                pass

        return r

    def particleGallery(self, destinationFolder, transform=True, applyClassColorLabeling=False, highest_frequency=None):
        """
        particleGallery: will create a particle gallery in destination folder
        @param destinationFolder: Destination folder where slices will be written
        @param transform: Apply transformation (according to alignment parameters) to each particle? Default is True!
        @param applyClassColorLabeling: Colorlabel particles according to their class? Default is false!
        """
        from pytom.agnostic.filter import bandpass as lowpassFilter
        from pytom.tools.toImage import volumeToPNG
        from pytom.tools.files import checkDirExists, dump2TextFile

        if not checkDirExists(destinationFolder):
            raise IOError('Destination folder does not exist!')

        html = '<html><title>Particle Gallery</title><body><table border=0>\n'
        html += '<tr><td>Filename</td><td>X Slice</td><td>Y Slice</td><td>Z Slice</td></tr>\n'

        if applyClassColorLabeling:
            self.sortByClassLabel()
            colors = ['#4169E1', '#FA8072']

            currentClassLabel = self[0].getClass()
            currentColor = 0

        for particleIndex in range(len(self)):

            particle = self[particleIndex]

            if transform:
                particleVolume = particle.getTransformedVolume()
            else:
                particleVolume = particle.getVolume()

            if highest_frequency:
                particleVolume = lowpassFilter(particleVolume, high=highest_frequency, sigma=float(highest_frequency) / 10)

            if applyClassColorLabeling:
                if particle.getClass() != currentClassLabel:
                    currentColor = (currentColor + 1) % len(colors)
                    currentClassLabel = particle.getClass()

                html += '<tr><td><font color="' + colors[
                    currentColor] + '">' + particle.getFilename() + '</font></td>\n'
            else:
                html += '<tr><td>' + particle.getFilename() + '</td>\n'

            volumeToPNG(particleVolume, destinationFolder + '/' + str(particleIndex) + '-x.png', projectionAxis='x')
            html += '<td><img src="' + './' + str(particleIndex) + '-x.png"/></td>\n'

            volumeToPNG(particleVolume, destinationFolder + '/' + str(particleIndex) + '-y.png', projectionAxis='y')
            html += '<td><img src="' + './' + str(particleIndex) + '-y.png"/></td>\n'

            volumeToPNG(particleVolume, destinationFolder + '/' + str(particleIndex) + '-z.png', projectionAxis='z')
            html += '<td><img src="' + './' + str(particleIndex) + '-z.png"/></td></tr>\n'

        html += '</body></html>'
        dump2TextFile(destinationFolder + '/gallery.html', html, False)

    def splitNSublists(self, n):
        """
        splitNSublists: Will split list into n+1 equally sized sublists.
        @param n:
        @return: List of sublists.
        @author: Thomas Hrabe
        """
        if n >= len(self) or n == 0 or n == 1:
            # do nothing
            return [self]

        from math import floor
        stepSize = int(floor(len(self) / float(n)))

        lastIndex = 0
        lists = []

        for i in range(stepSize, len(self), stepSize):
            lists.append(self[lastIndex:i])
            lastIndex = i
            if i + stepSize >= len(self):
                lists.append(self[i:len(self)])

        return lists

    def setClassesFromList(self, otherParticleList):
        """
        setClassesFromList: Assign class membership to each particle according to other particle list
        @param otherParticleList: particle list
        @type otherParticleList: L{ParticleList}
        """
        for particle in self._particleList:
            otherParticle = otherParticleList.getParticleByFilename(particle.getFilename())

            particle.setClass(otherParticle.getClassName())

    def sortByScore(self, sortType='descending'):
        """
        sortByScore:
        @param sortType: can be either 'descending' (default) or ascending
        """
        self._particleList = sorted(self._particleList, key=keyFunctionForParticleSortingByScore,
                                    reverse=sortType == 'descending')

    def sortByClassLabel(self, sortType='ascending'):
        """
        sortByClassLabel:
        @param sortType: can be either 'ascending' (default) or 'descending'
        @type sortType: str
        """
        self._particleList = sorted(self._particleList, key=keyFunctionForParticleSortingByClassLabel,
                                    reverse=(sortType == 'descending'))

    def createDeNovoReference(self, name, num=0):
        """
        createDeNovoReference: Use the particle list to create a template for reference-free alignment.\
        Will randomly rotate all particles and average them.
        @param name: the created reference name. It is automatically written to the disk.
        @param num: use how many particles for creating the template. (Default is 0: use all!)
        @return: The randomized particle list
        """
        import numpy
        total_num = len(self)
        assert num >= 0 and num <= total_num

        if num == 0 or num == total_num:
            random_index = list(range(0, total_num))
        else:
            random_index = numpy.random.randint(0, total_num, num)

        random_pl = ParticleList()
        for i in random_index:
            p = self._particleList[i]
            # randomize the rotation
            random_rotation = Rotation(numpy.random.rand() * 360, numpy.random.rand() * 360, numpy.random.rand() * 180)
            p.setRotation(random_rotation)
            random_pl.append(p)

        # create the average and write it to the disk
        random_pl.average(name, True)

        return random_pl

    def getVarianceMap(self, average, verbose=True):
        """
        Calculate the variance map on the aligned particle list.
        @param average: average of the aligned particle list
        @type average: L{pytom.lib.pytom_volume.vol}
        @param verbose: verbose mode
        @type verbose: L{boolean}
        @return: 3D variance map
        @rtype: L{pytom.lib.pytom_volume.vol}
        """
        from pytom.lib.pytom_volume import variance
        from pytom.tools.ProgressBar import FixedProgBar
        from pytom.agnostic.filter import applyFourierFilter

        if verbose:
            progressBar = FixedProgBar(0, len(self._particleList), '')
            progressBar.update(0)
            num_processed = 0

        sum_w = None
        sum_var = None

        for p in self._particleList:
            v = p.getTransformedVolume()
            # rotate the wedge
            w = p.getWedge()
            if sum_w is None:
                sum_w = w.returnWedgeVolume(v.shape[0], v.shape[1], v.shape[2], False, p.getRotation().invert())
            else:
                sum_w += w.returnWedgeVolume(v.shape[0], v.shape[1], v.shape[2], False, p.getRotation().invert())

            # calculate the variance
            wa = w.apply(average, p.getRotation().invert())
            # rescale the vol
            vv = v / xp.sqrt(variance(v, False)) * xp.sqrt(variance(wa, False))
            vv.shiftscale(wa.mean() - vv.mean(), 1)

            var = wa - vv
            var = var**2

            if sum_var is None:
                sum_var = var
            else:
                sum_var += var

            if verbose:
                num_processed += 1
                progressBar.update(num_processed)

        # take care the wedge weighting
        result = applyFourierFilter(sum_var, 1/sum_w)

        return result

    def getSTDMap(self, average=None, mask=None, verbose=True):
        """
        Calculate the standard deviation map using getVarianceMap function.

        @param average: average of the aligned particle list
        @type average: L{pytom.lib.pytom_volume.vol}
        @param verbose: verbose mode
        @type verbose: L{boolean}
        @return: 3D standard deviation map
        @rtype: L{pytom.lib.pytom_volume.vol}
        """
        if average is None:
            from pytom.alignment.alignmentFunctions import average2
            average, fsc = average2(self, False, False, False, mask, verbose=verbose)

        from pytom.lib.pytom_volume import abs, power
        vm = self.getVarianceMap(average, verbose)

        std_map = abs(vm)  # sometimes it can have negative values
        power(std_map, 0.5)

        if mask:
            map = std_map * mask  # get rid of the thing outside the mask
        else:
            map = std_map

        return map

    def getCVMap(self, average=None, std_map=None, threshold=None, negative_density=True, verbose=True):
        """Calculate the coefficient of variance map.
        """
        from pytom.lib.pytom_volume import vol, variance, mean, abs
        if average is None:
            from pytom.alignment.alignmentFunctions import average2
            average, fsc = average2(self, False, False, False, verbose=verbose)

        if std_map is None:
            std_map = self.getSTDMap(average, verbose=verbose)

        if threshold is None:
            if negative_density:
                threshold = mean(average) - 3 * variance(average, False) ** 0.5
            else:
                threshold = mean(average) + 3 * variance(average, False) ** 0.5

        avg_copy = vol(average)
        #        if negative_density:
        #            limit(avg_copy, 0,0, threshold,0, False, True)
        #        else:
        #            limit(avg_copy, threshold,0, 0,0, True, False)

        cv_map = std_map / abs(avg_copy)

        return cv_map

    def getRandomPortion(self, portion=0.5):
        """Get a random portion of particles out.

        @param portion: percentage or int.
        """
        if portion < 1:
            n = int(len(self) * float(portion))
        else:
            n = portion

        if n < 0 or n > len(self):
            raise ValueError("Invalid argument.")

        import random
        sample = list(range(len(self)))
        random_sample = [sample[i] for i in sorted(random.sample(range(len(sample)), n))]

        res = ParticleList()
        for s in random_sample:
            res.append(self[s])

        return res

    def loadCoordinateFile(self, filename, name_prefix=None, wedge_angle=None, infoGUI=None):
        """
        Initialize the particle list using the given coordinate file.
        The coordinate file simply contains three columns of X, Y and Z separated
        by tabs (compatible with EMAN2).

        @param filename: Coordinate file name.
        @param name_prefix: Particle name prefix
        @param wedge_angle: angle(s) specifying single axis tilt
        @type wedge_angle: float or 2-dim list of floats
        """
        if not name_prefix:
            name_prefix = './particle_'

        try:
            self._particleList
        except:
            self._particleList = []
        if wedge_angle:
            wedge = SingleTiltWedge(wedge_angle=wedge_angle)
        try:
            f = open(filename, 'r')
            ff = [line for line in f.readlines()]
            i = 0

            for line in ff:
                if '#' in line or not line: continue
                try:
                    x, y, z = [float(n) for n in line.split('\t')]
                except:
                    x, y, z = [float(n) for n in line.split()]
                p = Particle(name_prefix + str(i) + '.em', rotation=None, shift=None,
                             wedge=None, className=0, pickPosition=PickPosition(x, y, z),
                             score=None, infoGUI=infoGUI)
                if wedge_angle:
                    p.setWedge(wedge)
                self._particleList.append(p)
                i += 1
        finally:
            f.close()

    def loadCoordinateFileHeader(self, filename):
        header = [line.strip('#').split() for line in open(filename, 'r').readlines() if '#' in line]
        headerInfo = ['', 1, 1]

        for l in header:
            if not l or len(l) < 2:
                continue
            if l[0] == 'TOMONAME':
                headerInfo[0] = l[1]
            elif l[0] == 'MARKERINDEX':
                headerInfo[1] = l[1]
            elif l[0] == 'BINNINGFACTOR':
                headerInfo[2] = int(l[1])

        return headerInfo

    def listFromBestScorePercentage(self, percentage):
        """
        listFromBestParticles: Returns a particle list from all N% best scored particles
        @param percentage: The percentage, set 0 - 100 (note - 0 will return an empty list!)
        @type percentage: float
        @return: A new particle List
        """
        self.sortByScore(sortType='descending')
        numberParticles = int(len(self) * percentage / 100.0)
        return self[:numberParticles]

    def updateFromPeaks(self, peaks):
        """
        update shifts and rotations in particle list according to values in peaks
        @param peaks: list of peaks
        @type peaks: list of L{pytom.alignment.structures.Peak}
        @author: FF
        """
        from pytom.alignment.structures import Peak
        assert type(peaks[0]) == Peak
        assert len(self) == len(peaks)
        for (ii, part) in enumerate(self._particleList):
            part.setRotation(rotation=peaks[ii].getRotation())
            part.setShift(shift=peaks[ii].getShift())
            part.setScoreValue(scoreValue=peaks[ii].getScoreValue())

    def setFileName(self, filename):
        """
        set filename of particle list
        @type filename: L{str}
        """
        self._XMLfilename = filename

    def getFileName(self):
        """
        set filename of particle list
        """
        return self._XMLfilename

    def pickle(self):
        """Make this particle list picklable.
        """
        from pytom.alignment.FRMAlignment import FRMScore

        for p in self._particleList:
            if p.getScore().__class__ != FRMScore:
                p.setScore(FRMScore(p.getScore().getValue()))

    def addRotation(self, rot):
        """
        add a rotation to a particle list, i.e., all particles causing effective rotation of average by the inverse
        @param rot: rotation to be applied IN ADDITION to the one stored for each particle
        @type rot: L{pytom.agnostic.structures.Rotation}
        @author: FF
        """
        assert isinstance(rot, Rotation), "addRotation: rot must be Rotation!"
        for particle in self._particleList:
            trot = particle.getRotation()
            particle.setRotation(rotation=trot * rot)

    def addShift(self, translation):
        """
        add a translation to a particle list, i.e., all particles causing effective translation of average by the
        inverse. In pseudocode: p.newTranslation = p.currentTranslation + p.currentRotation.rotate(newTranslation)
        @param translation: translation to be applied IN ADDITION to the one stored for each particle
        @type translation: L{pytom.agnostic.structures.Shift}
        @author: FF
        """
        assert isinstance(translation, Shift), "addShift: translation must be Shift!"
        for particle in self._particleList:
            ttrans = particle.getShift()
            trot = particle.getRotation()
            particle.setShift(shift=ttrans + translation.rotate(rot=trot))


class InformationGUI(PyTomClass):
    """
    InformationGUI: Stores information used by GUI: project directory.
    """

    def __init__(self, tomogramName='', refMarkIndex=0, binningFactor=1):
        self._tomoname = tomogramName
        self._refMarkIndex = int(refMarkIndex)
        self._binningFactor = int(binningFactor)

    def getTomoName(self):
        return self._tomoname

    def setTomoName(self, tomoname):
        self._tomoname = str(tomoname)

    def getRefMarkIndex(self):
        return self._refMarkIndex

    def setRefMarkIndex(self, index):
        self._refMarkIndex = int(index)

    def getBinningFactor(self):
        return self._binningFactor

    def setBinningFactor(self, factor):
        self._binningFactor = int(factor)

    def toXML(self):
        from lxml import etree

        info_element = etree.Element('InfoTomogram', TomoName=str(self._tomoname),
                                     RefMarkIndex=str(self._refMarkIndex), BinningFactor=str(self._binningFactor))
        return info_element

    def fromXML(self, xmlObj):
        from lxml.etree import _Element

        if xmlObj.__class__ != _Element:
            raise TypeError('Is not a lxml.etree._Element! You must provide a valid XMLobject.')

        if xmlObj.tag == 'InfoGUI':
            info_element = xmlObj
        else:
            TypeError('InfoTomogram: You must provide a valid InformationGUI XML object.')

        self._tomoname = str(info_element.get('TomoName'))
        self._refMarkIndex = int(info_element.get('RefMarkIndex'))
        self._binningFactor = int(info_element.get('BinningFactor'))

        if self._binningFactor < 1:
            self._binningFactor = 1


class Rotation(PyTomClass):
    """
    Rotation: Stores a 3D Euler Rotation according to multiple supported \
    paradigms. supported are: (ZXZ = default, ZYZ, YZY, XYZ, ZYX, ...).
    In order to use the pytom transformations, stick to getZ1,getZ2,getX \
    methods. They will always return ZXZ angles
    """

    def __init__(self, z1=None, z2=None, x=None, paradigm='ZXZ'):

        if z1.__class__ == list or z1.__class__ == Rotation:
            self._z1 = float(z1[0])
            self._z2 = float(z1[1])
            self._x = float(z1[2])
        elif z1 == None and z2 == None and x == None:
            self._z1 = 0.0
            self._z2 = 0.0
            self._x = 0.0
        else:
            self._z1 = float(z1)
            self._z2 = float(z2)
            self._x = float(x)

        self._checkParadigm(paradigm)

        self._paradigm = paradigm

    def _checkParadigm(self, paradigm):
        """
        _checkParadigm
        @param paradigm: specification of rotation 'ZXZ', 'ZYZ', 'YZY', 'XYZ', 'ZYX'
        @type paradigm: str
        """
        from pytom.basic.exceptions import ParameterError
        if (len(paradigm) != 3) or not (paradigm in ['ZXZ', 'ZYZ', 'YZY', 'XYZ', 'ZYX']):
            raise ParameterError('Rotation paradigm can be ZXZ, ZYZ, YZY, XYZ, ZYX. \
	    The one you specified is wrong.')

    def getParadigm(self):
        """
        getParadigm
        """
        return self._paradigm

    def setParadigm(self, paradigm):
        """
        setParadigm
        @param paradigm: The rotation paradigm - specification of rotation 'ZXZ', 'ZYZ', 'YZY', 'XYZ', 'ZYX'
        @type paradigm: str
        """
        self._checkParadigm(paradigm)
        self._paradigm = paradigm

    def getPhi(self):
        """
        getPhi: First rotation in current paradigm
        """
        return self._z1

    def getPsi(self):
        """
        getPsi: Third rotation in current paradigm
        """
        return self._z2

    def getTheta(self):
        """
        getTheta: Second rotation in current paradigm
        """
        return self._x

    def setPhi(self, phi):
        """
        setPhi: Set first rotation in current paradigm
        """
        self._z1 = phi

    def setPsi(self, psi):
        """
        setPsi: Set third rotation in current paradigm
        """
        self._z2 = psi

    def setTheta(self, theta):
        """
        setTheta: Set second rotation in current paradigm
        """
        self._x = theta

    def getZ1(self):
        return self._z1

    def getZ2(self):
        return self._z2

    def getX(self):
        return self._x

    def setZ1(self, z1):
        """
        set z1 rotation
        @param z1: z1 rotation
        @type z1: float
        """
        self._z1 = z1

    def setZ2(self, z2):
        """
        set z2 rotation
        @param z2: z2 rotation
        @type z2: float
        """
        self._z2 = z2

    def setX(self, x):
        """
        set x rotation
        @param x: x rotation
        @type x: float
        """
        self._x = x

    def toXML(self):
        """
        """
        from lxml import etree

        rotation_element = etree.Element('Rotation', Z1=str(self._z1), Z2=str(self._z2), X=str(self._x),
                                         Paradigm=self._paradigm)
        return rotation_element

    def fromXML(self, xmlObj):
        from lxml.etree import _Element

        if xmlObj.__class__ != _Element:
            raise TypeError('Is not a lxml.etree._Element! You must provide a valid XMLobject.')

        if xmlObj.tag == 'Rotation':
            rotation_element = xmlObj
        else:
            TypeError('Rotation: You must provide a valid Rotation XML object.')

        self._z1 = float(rotation_element.get('Z1'))
        self._z2 = float(rotation_element.get('Z2'))
        self._x = float(rotation_element.get('X'))

        self._paradigm = rotation_element.get('Paradigm')
        if self._paradigm == None:
            self._paradigm = 'ZXZ'

    def toList(self):
        """
        toList:
        @deprecated: Yes, use toVector instead!
        @return: List [z1,z2,x]
        """
        return [self._z1, self._z2, self._x]

    def toVector(self, convention='zzx'):
        """
        toVector:
        @param convention: how do you want your angles returned. default zzx
        @type z2: string
        @return: List [z1,z2,x]
        """

        dictAngles = {'z1': self._z1, 'z2': self._z2, 'x1': self._x}
        num_repeat = {'z': 1, 'x': 1, 'y': 1}
        outname = []
        for angle in convention:
            q = angle + str(num_repeat[angle])
            num_repeat[angle] += 1
            outname.append(dictAngles[q])

        return outname

    def toQuaternion(self):
        """
        toQuaternion:
        @return: This rotation L{pytom.angles.quaternions.Quaternion} object
        """
        from pytom.angles.quaternions import Quaternion

        return Quaternion(self._z1, self._z2, self._x)

    def __getitem__(self, key):
        """
        __getitem__:Returns rotation elements. 0th is z1, 1st is z2, 2nd is x
        """
        if key < 0:
            raise IndexError('Index out of range for Rotation class!')
        if key == 0:
            return self._z1
        if key == 1:
            return self._z2
        if key == 2:
            return self._x
        if key > 2:
            raise IndexError('Index out of range for Rotation class!')

    def __setitem__(self, key, value):
        """
        __setitem__: Enables syntax such as Rotation[key] = value
        @param key: The index
        @type key: int
        @param value: The valueset
        @type value: L{float}
        """
        if isinstance(key, int):
            if key < 0:
                raise IndexError('Index out of range for Rotation class!')
            if key == 0:
                self._z1 = value
            elif key == 1:
                self._z2 = value
            elif key == 2:
                self._x = value
            else:
                raise IndexError('Index out of range for Rotation class!')
        else:
            raise KeyError('Rotation class does not support slices!')

    def copy(self):
        return Rotation(self._z1, self._z2, self._x)

    def __add__(self, firstRotation):
        """
        __add__ : Adds this as a second Rotation to another Rotation object ->   R_self(R_firstRotation(X))
        @param firstRotation: L{Rotation}
        @todo: Do in Quaternion space
        """
        if self._paradigm == 'ZXZ':
            from pytom.angles.angleFnc import matToZXZ as matToR

        assert firstRotation.__class__ == Rotation
        m1 = firstRotation.toMatrix()
        m2 = self.toMatrix()
        m3 = m2 * m1
        r = matToR(m3)
        return Rotation(z1=r[0], z2=r[1], x=r[2])

    def __mul__(self, firstRotation):
        """
        __mul__: Does the same as self.__add__(firstRotation)
        """
        return self + firstRotation

    def invert(self):
        """
        invert: Will return this object as inverted rotation.
        @return: inverted Rotation
        """
        return Rotation(-self.getZ2(), -self.getZ1(), -self.getX())

    def toMatrix(self, fourByfour=False):
        """
        toMatrix: Returns a rotation matrix of this object
        @param fourByfour: return 4x4 matrix (default is False, return 3x3)
        @rtype: 3x3 list
        """
        if self._paradigm == 'ZXZ':
            from pytom.angles.angleFnc import zxzToMat as toMat
            res = toMat(*self.toVector())
        else:
            raise RuntimeError("Rotation: Paradigm other than ZXZ not supported yet!")

        if fourByfour is False:
            return res
        else:
            from pytom.tools.maths import Matrix
            mtx = Matrix(4, 4)
            mtx[0, 0] = res[0, 0]
            mtx[0, 1] = res[0, 1]
            mtx[0, 2] = res[0, 2]
            mtx[1, 0] = res[1, 0]
            mtx[1, 1] = res[1, 1]
            mtx[1, 2] = res[1, 2]
            mtx[2, 0] = res[2, 0]
            mtx[2, 1] = res[2, 1]
            mtx[2, 2] = res[2, 2]
            mtx[3, 3] = 1

            return mtx

    def __eq__(self, otherRotation):
        """
        __eq__: Compares two rotations
        """

        m1 = self.toMatrix()
        m2 = otherRotation.toMatrix()

        return m1 == m2

    def convert2pytomc(self):
        from pytom.basic.structures import Rotation
        return Rotation(*self.toVector())


class Shift(PyTomClass):
    """
    Shift
    """

    def __init__(self, x=None, y=0.0, z=0.0):
        """
        @param x: x-component of vector or vector
        @type x: float or 3-dim list
        @param y: y-component of vector or vector
        @type y: float
        @param z: z-component of vector or vector
        @type z: float
        """
        if x.__class__ == Shift or x.__class__ == list:
            self._x = float(x[0])
            self._y = float(x[1])
            self._z = float(x[2])
        elif x is not None:
            self._x = float(x)
            self._y = float(y)
            self._z = float(z)
        else:
            self._x = 0.0
            self._y = 0.0
            self._z = 0.0

    def __getitem__(self, key):
        """
        __getitem__:Returns rotation elements. 0th is x, 1st is y, 2nd is z
        """
        if key < 0:
            raise IndexError('Index out of range for Rotation class!')
        if key == 0:
            return self._x
        if key == 1:
            return self._y
        if key == 2:
            return self._z
        if key > 2:
            raise IndexError('Index out of range for Rotation class!')

    def __setitem__(self, key, value):
        """
        __setitem__: Enables syntax such as Shift[key] = value
        @param key: The index
        @type key: int
        @param value: A value
        @type value: L{float}
        """
        if isinstance(key, int):
            if key < 0:
                raise IndexError('Index out of range for Rotation class!')
            if key == 0:
                self._x = value
            elif key == 1:
                self._y = value
            elif key == 2:
                self._z = value
            else:
                raise IndexError('Index out of range for Rotation class!')
        else:
            raise KeyError('Rotation class does not support slices!')

    def getX(self):
        return self._x

    def getY(self):
        return self._y

    def getZ(self):
        return self._z

    def setX(self, x):
        assert x.__class__ == float or isinstance(x, int)
        self._x = x

    def setY(self, y):
        assert y.__class__ == float or isinstance(y, int)
        self._y = y

    def setZ(self, z):
        assert z.__class__ == float or isinstance(z, int)
        self._z = z

    def toXML(self):
        """
        """
        from lxml import etree

        shift_element = etree.Element('Shift', X=str(self._x), Y=str(self._y), Z=str(self._z))

        return shift_element

    def fromXML(self, xmlObj):
        """

        """
        from lxml.etree import _Element

        if xmlObj.__class__ != _Element:
            raise TypeError('Shift: Is not a lxml.etree._Element! You must provide a valid XMLobject.')

        if xmlObj.tag == 'Shift':
            shift_element = xmlObj
        else:
            TypeError('Shift: You must provide a valid Shift object.')

        self._x = float(shift_element.get('X'))
        self._y = float(shift_element.get('Y'))
        self._z = float(shift_element.get('Z'))

    def toVector(self):
        """
        toVector:
        @return: [x,y,z]
        """
        return [float(self._x), float(self._y), float(self._z)]

    def toMatrix(self):
        from pytom.tools.maths import Matrix
        mtx = Matrix(4, 4)
        mtx[0, 0] = 1
        mtx[1, 1] = 1
        mtx[2, 2] = 1
        mtx[3, 3] = 1
        mtx[0, 3] = float(self._x)
        mtx[1, 3] = float(self._y)
        mtx[2, 3] = float(self._z)

        return mtx

    def __getitem__(self, key):
        """
        __getitem__:Returns shift elements. 0th is x, 1st is y, 2nd is z
        """
        if key == 0:
            return self._x
        if key == 1:
            return self._y
        if key == 2:
            return self._z

    def getLength(self):
        """
        getLength:
        @deprecated: Use len(Object) instead!
        """
        print('Shift.getLength is deprecated! Use len(Shift) instead!')
        return len(self)

    def __len__(self):
        from pytom.agnostic.tools import euclidian_distance

        v1 = [0.0, 0.0, 0.0]
        v2 = self.toVector()

        return euclidian_distance(v1, v2)

    def scale(self, scalefactor):
        """
        Scale the shift by a factor
        @param scalefactor: 'stretch' of shift
        @type scalefactor: L{float}
        """
        self._x = self._x * float(scalefactor)
        self._y = self._y * float(scalefactor)
        self._z = self._z * float(scalefactor)

    def invert(self):
        """
        invert: Will return this object as inverted shift.
        @return: inverted shift (note: shift itself remains unaltered)
        @rtype: L{pytom.agnostic.structures.Shift}
        """
        return Shift(-self.getX(), -self.getY(), -self.getZ())

    def addVector(self, vector):
        if len(vector) == 3:
            self._x += vector[0]
            self._y += vector[1]
            self._z += vector[2]

    def rotate(self, rot):
        """
        rotate shift by rotation
        @param rot: rotation
        @type rot: L{pytom.agnostic.structures.Rotation}
        @return: rotated shift (note: shift itself remains unaltered)
        @rtype: L{pytom.agnostic.structures.Shift}
        """
        assert isinstance(rot, Rotation), "rot must be of type Rotation"
        m = rot.toMatrix(fourByfour=True) * self.toMatrix()
        return Shift(x=m.getColumn(3)[0], y=m.getColumn(3)[1], z=m.getColumn(3)[2])

    def __mul__(self, otherShift):
        """
        __mul__: Overloaded * operator
        @param otherShift: The other shift (or scalar) to cross product to the current shift. If scalar apply element wise multiplication.
        @type otherShift: Shift
        @return: A new Shift object with the current result
        """

        if not otherShift.__class__ == self.__class__ and not otherShift.__class__ in [int, float]:
            raise NotImplementedError("Shift mult: Add partner must be another Shift or int, float!")

        if otherShift.__class__ in [int, float]:
            newShift = Shift(self[0], self[1], self[2])
            newShift.scale(otherShift)
            return newShift

        x = self[1] * otherShift[2] - self[2] * otherShift[1]
        y = self[2] * otherShift[0] - self[0] * otherShift[2]
        z = self[0] * otherShift[1] - self[1] * otherShift[0]

        return Shift(x, y, z)

    def __rmul__(self, otherShift):
        return self * otherShift

    def __add__(self, otherShift):
        """
        __add__: Overloaded + operator
        @param otherShift: The other shift (or scalar) to add to the current shift
        @type otherShift: Shift
        @return: A new Shift object with the current result
        """
        if not otherShift.__class__ == self.__class__ and not otherShift.__class__ in [int, float]:
            raise NotImplementedError("Shift add: Add partner must be another Shift or int,float!")

        if otherShift.__class__ in [int, float]:
            return Shift(self[0] + otherShift, self[1] + otherShift, self[2] + otherShift)

        return Shift(self[0] + otherShift[0], self[1] + otherShift[1], self[2] + otherShift[2])

    def __radd__(self, otherShift):
        return self + otherShift

    def __sub__(self, otherShift):
        """
        __sub__: Overloaded - operator
        @param otherShift: The other shift (or scalar) to subtract to the current shift
        @type otherShift: Shift
        @return: A new Shift object with the current result
        """
        if not otherShift.__class__ == self.__class__ and not otherShift.__class__ in [int, float]:
            raise NotImplementedError("Shift sub: Add partner must be another Shift, too!")

        return self + (-1 * otherShift)

    def convert2pytomc(self):
        from pytom.basic.structures import Shift
        return Shift(*self.toVector())


class PickPosition(PyTomClass):
    """
    PickPosition: Specifies the position of a particle within a tomogram.
    """

    def __init__(self, x=0.0, y=0.0, z=0.0, originFilename=''):

        if x.__class__ == list:
            self._x = x[0]
            self._y = x[1]
            self._z = x[2]
        else:
            self._x = x
            self._y = y
            self._z = z

        self._originFilename = originFilename

    def getX(self):
        return self._x

    def getY(self):
        return self._y

    def getZ(self):
        return self._z

    def setX(self, x):
        self._x = x

    def setY(self, y):
        self._y = y

    def setZ(self, z):
        self._z = z

    def scale(self, scalefactor):
        """
        Scale the shift by a factor
        @param scalefactor: 'stretch' of shift
        @type scalefactor: L{float}
        """
        self._x = self._x * float(scalefactor)
        self._y = self._y * float(scalefactor)
        self._z = self._z * float(scalefactor)

    def getOriginFilename(self):
        return self._originFilename

    def toXML(self):
        """
        """
        from lxml import etree

        shift_element = etree.Element('PickPosition', X=str(self._x), Y=str(self._y), Z=str(self._z),
                                      Origin=str(self._originFilename))

        return shift_element

    def fromXML(self, xmlObj):
        from lxml.etree import _Element

        if xmlObj.__class__ != _Element:
            raise TypeError('Is not a lxml.etree._Element! You must provide a valid XMLobject.')

        if xmlObj.tag == 'PickPosition':
            shift_element = xmlObj
        else:
            TypeError('PickPosition: You must provide a valid PickPosition XML object.')

        self._x = float(shift_element.get('X'))
        self._y = float(shift_element.get('Y'))
        self._z = float(shift_element.get('Z'))
        self._originFilename = shift_element.get('Origin')

    def toVector(self):
        """
        toVector:
        @return: [x,y,z]
        """
        return [self._x, self._y, self._z]

    def __add__(self, vector):
        """
        __add__: Add a [x,y,z] list (vector) to the current position
        """
        if not (vector.__class__ == list and len(vector) == 3):
            RuntimeError('PickPosition: You can only add a [x,y,z] list to PickPosition!')

        self._x += vector[0]
        self._y += vector[1]
        self._z += vector[2]

    def __sub__(self, vector):

        from pytom.tools.maths import scale

        self.__add__(scale(vector, -1))


class Symmetry(PyTomClass):

    #     def __init__(self,nfold=1,z2=0,x=0, search_ang=0, search_axis_z2=0, search_axis_x=0):
    #         """
    #         __init__ : Implementation for backward compatibility only. Overwrites self to PointSymmetry object
    #         """
    #
    #         self = PointSymmetry(nfold=1,z2=0,x=0, search_ang=0, search_axis_z2=0, search_axis_x=0)

    def fromXML(self, xmlObj):
        """
        fromXML:
        @param xmlObj: A xml object
        @type xmlObj: L{lxml.etree._Element}
        @author: Thomas Hrabe
        """
        from lxml.etree import _Element

        if xmlObj.__class__ != _Element:
            raise TypeError('Is not a lxml.etree._Element! You must provide a valid XMLobject.')

        if xmlObj.tag == 'Symmetry':
            symmetryElement = xmlObj
        else:
            TypeError('SymetryObject: you must provide a valid Symmetry XML object.')

        symmetryType = xmlObj.get("Type")
        symmetry = None

        if symmetryType is None or symmetryType == 'PointSymmetry':
            symmetry = PointSymmetry()
            symmetry.fromXML(xmlObj)
        elif symmetryType == 'HelicalSymmetry':
            symmetry = HelicalSymmetry()
            symmetry.fromXML(xmlObj)

        return symmetry

    def fromStr(self, xmlString):
        if self.__class__ == Symmetry:
            from lxml import etree
            root = etree.fromstring(xmlString)
            symmetryObject = self.fromXML(root)

            return symmetryObject
        else:
            raise RuntimeError(
                'Children of Symmetry may not use fromStr method. Use the parent\'s Symmetry.fromStr instead!')

    def copy(self):
        sym = Symmetry()
        symXML = self.toXML()
        return sym.fromXML(symXML)

    def apply(self, particleList):

        raise RuntimeError('Symmetry - apply : this is a virtual method')

    def applyToParticle(self, particleList):

        raise RuntimeError('Symmetry - applyToParticle : this is a virtual method')

    def setXTilt(self, x):
        self._x = x

    def setZ2Tilt(self, z2):
        self._z2 = z2

    def getNFold(self):
        """
        getNFold
        """
        return self._nfold

    def isOneFold(self):
        """
        isOneFold: Returns True if self._numberSymmetries == 1
        @rtype: Bool
        """
        return self._nfold == 1

    def getZ2(self):
        return self._z2

    def getX(self):
        return self._x


class PointSymmetry(Symmetry):
    """
    PointSymmetry: Stores all neccessary information about the symmetry of an object
    """

    def __init__(self, nfold=1, z2=0, x=0, search_ang=0, search_axis_z2=0, search_axis_x=0):
        """
        __init__:
        @param nfold: symmetry number around a given axis (Default is the Z axis)
        @param z2: Phi of symmetry axis
        @param x: rotation around x-axis (the) of symmetry axis
        @param search_ang: this and the next two paramters are introduced to tackle the situation that after the symmetry rotation we still have to search for the best match
        @param search_axis_z2: search axis (z2)
        @param search_axis_x: search axis (x)
        """
        self._nfold = int(nfold)
        self._z2 = float(z2)
        self._x = float(x)

        self._search_ang = int(search_ang)
        self._search_axis_z2 = float(search_axis_z2)
        self._search_axis_x = float(search_axis_x)

    def getAngleList(self):
        """
        getAngleList
        """
        from pytom.angles.localSampling import LocalSampling

        symmetryIncrement = 360 / self._nfold

        return LocalSampling(0, symmetryIncrement, 0, self._z2, self._x)

    def toXML(self):
        """
        toXML : Compiles a XML file from result object
        rtype : L{lxml.etree._Element}
        @author: Thomas Hrabe
        """
        from lxml import etree

        return etree.Element('Symmetry', Type="PointSymmetry", NFold=str(self._nfold), Z2=str(self._z2), X=str(self._x),
                             SearchAngle=str(self._search_ang), SearchAxisZ2=str(self._search_axis_z2),
                             SearchAxisX=str(self._search_axis_x))

    def fromXML(self, xmlObj):
        """
        fromXML:
        @param xmlObj: A xml object
        @type xmlObj: L{lxml.etree._Element}
        @author: Thomas Hrabe
        """
        if xmlObj.tag == 'Symmetry':
            symmetryElement = xmlObj
        try:
            self._nfold = int(symmetryElement.get('NFold'))
            self._z2 = float(symmetryElement.get('Z2'))
            self._x = float(symmetryElement.get('X'))
        except TypeError:
            self._nfold = int(symmetryElement.get('NumberSymmetries'))
            self._z2 = float(symmetryElement.get('Psi'))
            self._x = float(symmetryElement.get('Theta'))

        try:
            searchAngle = symmetryElement.get('SearchAngle')
            if ',' in searchAngle:  # it is a list specifying the search angle range
                self._search_ang = [int(i) for i in searchAngle[1:-1].split(',')]
            else:
                self._search_ang = int(searchAngle)

            self._search_axis_z2 = float(symmetryElement.get('SearchAxisZ2'))
            self._search_axis_x = float(symmetryElement.get('SearchAxisX'))
        except TypeError:
            self._search_ang = 0
            self._search_axis_z2 = 0
            self._search_axis_x = 0

    def apply(self, particleList):
        """
        apply: Applies this symmetry to a particle List
        @param particleList: The particle List
        @return: The symmetrized particle list
        """
        from pytom.tools.macros import frange

        if self.isOneFold():
            return particleList

        from pytom.agnostic.structures import ParticleList, Rotation
        newList = ParticleList(particleList.getDirectory())

        for i in range(len(particleList)):
            particle = particleList[i]

            for inplaneAngle in frange(0, 360, 360.0 / float(self._nfold)):
                # get a new object -> make sure we do not modify any property of the object
                p2 = particle.copy()
                particleRotation = p2.getRotation()

                # read below line from right to left!
                # multiply inverted particle Rotation with inverted symmetry axis
                # add inplane rotation for symmetry
                # orient back to alignment axis
                newRotation = Rotation(0, self._z2, self._x) * (Rotation(inplaneAngle, 0, 0) * (
                            Rotation(-self._z2, 0, -self._x) * particleRotation.invert()))

                # finally, we must invert the new rotation to put it into the correct start position for averaging!
                p2.setRotation(newRotation.invert())

                newList.append(p2)

        return newList

    def calculate_rotation_cc(self, volume, mask=None, axis='Z', angularIncrement=1, plot=False, verbose=False):
        """
        calculate_rotation_cc: For backward compatibility
        @deprecated: Use calculateCorrelationProfile instead
        """
        return self.calculateCorrelationProfile(self, volume, mask, axis, angularIncrement, plo, verbose)

    def calculateCorrelationProfile(self, volume, mask=None, axis='Z', angularIncrement=1, plot=False, verbose=False,
                                    interpolation='filt_bspline', rot_order='zxz'):
        """
        calculateCorrelationProfile: Will rotate volume around defined axis and correlate with original volume to determine symmetry.
        @param volume: The volume
        @param mask: A mask used for correlation & normalizaion. If none, mask of 1s will be used.
        @param axis: Symmetry axis. Z is default. Specify alternatives through x,z rotation like [10,20] for a tilded z axis by x=10 and z=20
        @param angularIncrement: Angular increment. 1 is default
        @param plot: Plot outcome? Default is False
        @param verbose:
        @return: Will return a list of correlation coefficients
        @author: Thomas Hrabe
        """

        from pytom.angles.localSampling import LocalSampling
        from pytom.agnostic.correlation import nxcc
        from pytom.agnostic.structures import Mask
        from pytom.voltools import StaticVolume

        volumeTexture = StaticVolume(volume, device=device, interpolation=interpolation)

        if axis == 'Z' or axis.__class__ != list:
            rotations = LocalSampling(0, angularIncrement)
        elif axis.__class__ == list and len(axis) == 2:
            rotations = LocalSampling(0, angularIncrement, startZ2=0, startX=0)
        else:
            raise ParameterError("Symmetry object: axis must either be 'Z' or a list of two tilt angles")

        if not mask:
            mask = vol(volume.size_x(), volume.size_y(), volume.size_z())
            mask.setAll(1)
        elif mask.__class__ == Mask:
            mask = mask.getVolume()

        ccList = []

        rotatedVolume = xp.zeros((volume.size_x(), volume.size_y(), volume.size_z()), dtype=xp.float32)

        currentRotation = rotations.nextRotation()

        while currentRotation != [None, None, None]:
            rotation = [currentRotation[0], currentRotation[2], currentRotation[1]]
            volumeTexture.transform(volume, output=rotatedVolume, rotation=rotation, rotation_order=rot_order)
            # rotatedVolume.write('rot'+str(currentRotation[0])+'.em')
            cc = nxcc(volume, rotatedVolume, mask)
            if verbose:
                print('Rotation ', currentRotation, ' Value', cc)

            ccList.append(cc)

            currentRotation = rotations.nextRotation()

        if plot:
            import matplotlib
            try:
                matplotlib.use('Qt5Agg')
            except:
                pass
            from matplotlib import pyplot
            pyplot.plot(ccList)
            pyplot.show()

        return ccList

    def applyToParticle(self, volume, interpolation='filt_bspline', rot_order='zxz'):
        """
        applyToParticle: symmetrize particle
        @param volume: The volume
        """
        from pytom.agnostic.structures import Rotation
        from pytom.voltools import StaticVolume

        volumeTexture = StaticVolume(volume, device=device, interpolation=interpolation)

        size_x,size_y,size_z = volume.shape

        symVolume = xp.zeros((size_x, size_y, size_z), dtype=xp.float32)
        rotVolume = xp.zeros((size_x, size_y, size_z), dtype=xp.float32)


        symmetryAngle = 360 / float(self._nfold)

        for i in range(self._nfold):
            if i == 0:
                symVolume = volume
                continue

            finalRotation = Rotation(0, self._z2, self._x) * Rotation(symmetryAngle * i, 0, 0) * Rotation(-self._z2, 0,
                                                                                                          -self._x)
            rotVolume *= 0
            # rotation = finalRotation.toVector(rot_order)
            # volumeTexture.transform(volume, output=rotVolume, rotation=rotation, rotation_order=rot_order)#finalRotation.getZ1(), finalRotation.getZ2(), finalRotation.getX())

            # search for the best match after the symmetry rotation
            if self._search_ang != 0:
                from pytom.agnostic.correlation import xcc
                tmp = xp.zeros((size_x, size_y, size_z), dtype=xp.float32)
                max_cc = None
                best_ang = None
                if isinstance(self._search_ang, int):
                    best_ang = self._search_ang
                else:
                    for ang in range(self._search_ang[0], self._search_ang[1]):
                        rot = Rotation(0, self._search_axis_z2, self._search_axis_x) * Rotation(ang, 0, 0) * Rotation(
                            -self._search_axis_z2, 0, -self._search_axis_x)

                        nrot = rot * finalRotation
                        rotation = nrot.toVector(rot_order)
                        volumeTexture.transform(volume, output=tmp, rotation=rotation, rotation_order=rot_order)
                        cc = xcc(volume, tmp)
                        if max_cc is None:
                            max_cc = cc
                            best_ang = ang
                        elif cc > max_cc:
                            max_cc = cc
                            best_ang = ang
                        else:
                            pass

                rot = Rotation(0, self._search_axis_z2, self._search_axis_x) * Rotation(best_ang, 0, 0) * Rotation(
                    -self._search_axis_z2, 0, -self._search_axis_x)

                finalRotation = rot * finalRotation

            rotation = finalRotation.toVector(rot_order)
            volumeTexture.transform(volume, output=rotVolume, rotation=rotation,
                                    rotation_order=rot_order)  # finalRotation.getZ1(), finalRotation.getZ2(), finalRotation.getX())

            symVolume = symVolume + rotVolume

        return symVolume


# TODO remove in next cleaning
# class MultiSymmetries(PyTomClass):
#     """Multiple symmetries support.
#     """
#
#     def __init__(self, symmetries=None):
#         if symmetries is not None:
#             if symmetries.__class__ != list:
#                 raise TypeError('The symmetries should be a list containing all the symmetries object!')
#             for s in symmetries:
#                 if s.__class__ != Symmetry:
#                     raise TypeError('Each of the object in the list should be of type Symmetry!')
#             self.symmetries = symmetries
#         else:
#             self.symmetries = []
#
#     def toXML(self):
#         from lxml import etree
#
#         jobElement = etree.Element("MultiSymmetries")
#         for s in self.symmetries:
#             jobElement.append(s.toXML())
#
#         return jobElement
#
#     def fromXML(self, xmlObj):
#         from lxml.etree import _Element
#         if xmlObj.__class__ != _Element:
#             raise TypeError('Is not a lxml.etree._Element! You must provide a valid XMLobject.')
#         if xmlObj.tag == 'MultiSymmetries':
#             element = xmlObj
#         else:
#             TypeError('MultiSymmetries: You must provide a valid MultiSymmetries XML object.')
#
#         symmetries_elements = element.xpath('Symmetry')
#         self.symmetries = []
#         for s in symmetries_elements:
#             sym = Symmetry()
#             sym = sym.fromXML(s)
#             self.symmetries.append(sym)
#
#     def applyToParticle(self, volume):
#         symVolume = volume
#         for sym in self.symmetries:
#             symVolume = sym.applyToParticle(symVolume)
#
#         return symVolume


class HelicalSymmetry(Symmetry):

    def __init__(self, nfold=None, symmetryAngle=None, repeatShift=None, isRightSymmetry=True, x=0, z2=0):

        self._nfold = nfold
        self._symmetryAngle = symmetryAngle
        self._repeatShift = repeatShift
        self._isRightSymmetry = isRightSymmetry
        self._x = x
        self._z2 = z2

    def apply(self, particleList):
        """
        apply: Applies this symmetry to a particle List
        @param particleList: The particle List
        @return: The symmetrized particle list
        """

        if self.isOneFold():
            return particleList

        from pytom.agnostic.structures import ParticleList, Rotation
        newList = ParticleList(particleList.getDirectory())

        factor = 1 if self._isRightSymmetry else -1

        for i in range(len(particleList)):
            particle = particleList[i]

            for i in range(int(self._nfold)):
                # get a new object -> make sure we do not modify any property of the object
                p2 = particle.copy()
                particleRotation = p2.getRotation()
                particleShift = p2.getShift()

                # read below line from right to left!
                # multiply inverted particle Rotation to register it on common rotation axis
                # add symmetry rotation
                newRotation = Rotation(self._symmetryAngle * i, 0, 0) * particleRotation.invert()
                newShift = Shift(0, 0, factor * self._repeatShift * i) + particleShift.invert()

                # finally, we must invert the new rotation to put it into the correct start position
                p2.setRotation(newRotation.invert())
                p2.setShift(newShift)

                newList.append(p2)

                p3 = particle.copy()
                particleRotation = p3.getRotation()
                particleShift = p3.getShift()

                # read below line from right to left!
                # multiply inverted particle Rotation to register it on common rotation axis
                # orient back to alignment axis
                newRotation = Rotation(-1 * self._symmetryAngle * i, 0, 0) * particleRotation.invert()
                newShift = Shift(0, 0, -1 * factor * self._repeatShift * i) + particleShift.invert()
                # finally, we must invert the new rotation to put it into the correct start position
                p3.setRotation(newRotation.invert())
                p3.setShift(newShift)

                newList.append(p3)

        return newList

    def applyToPatricle(self, particle, interpolation='filt_bspline'):
        from pytom.voltools import StaticVolume

        particleTexture = StaticVolume(particle, device=device, interpolation=interpolation)

        result = xp.zeros(particle.shape,dtype=xp.float32)
        temp = xp.zeros_like(result)

        factor = 1 if self._isRightSymmetry else -1

        for i in range(self._nfold):
            particleTexture.transform(rotation=[self._symmetryAngle * i, 0, 0], rotation_order='zxz', output=temp,
                                      translation2=[0, 0, factor * self._repeatShift * i])
            result += temp
            temp *= 0
            particleTexture.transform(rotation=[-self._symmetryAngle * i, 0, 0], rotation_order='zxz', output=temp,
                                      translation2=[0, 0, -1 * factor * self._repeatShift * i])
            result += temp

        return result

    def toXML(self):
        """
        toXML : Compiles a XML file from result object
        rtype : L{lxml.etree._Element}
        @author: Thomas Hrabe
        """
        from lxml import etree
        return etree.Element('Symmetry', Type="HelicalSymmetry", NFold=str(self._nfold),
                             RepeatShift=str(self._repeatShift), Z2=str(self._z2), X=str(self._x),
                             SymmetryAngle=str(self._symmetryAngle), IsRightSymmetry=str(self._isRightSymmetry))

    def fromXML(self, xmlObj):
        """
        fromXML:
        @param xmlObj: A xml object
        @type xmlObj: L{lxml.etree._Element}
        @author: Thomas Hrabe
        """
        from lxml.etree import _Element

        if xmlObj.__class__ != _Element:
            raise RuntimeError('Is not a lxml.etree._Element! You must provide a valid XMLobject.')

        if xmlObj.tag == 'Symmetry' and xmlObj.get('Type') == 'HelicalSymmetry':
            symmetry_element = xmlObj
        else:
            raise RuntimeError('Is not a HelicalSymmetry object! You must provide a valid HelicalSymmetry.')

        self._nfold = float(symmetry_element.get('NFold'))
        self._symmetryAngle = float(symmetry_element.get('SymmetryAngle'))
        self._repeatShift = float(symmetry_element.get('RepeatShift'))
        self._z2 = float(symmetry_element.get('Z2'))
        self._x = float(symmetry_element.get('X'))
        self._isRightSymmetry = symmetry_element.get('IsRightSymmetry') == 'True'


class Weight():
    def __init__(self, wedge_angle1=0, wedge_angle2=0, cutOffRadius=0, size_x=0, size_y=0, size_z=0, smooth=0, rotation=None):
        self.wedge_angle1 = wedge_angle1
        self.wedge_angle2 = wedge_angle2
        self.cutOffRadius = cutOffRadius
        self.size_x = size_x
        self.size_y = size_y
        self.size_z = size_z
        self.smooth=smooth
        self.rotation = rotation

    def rotate(self, phi, psi, theta):
        self.rotation = [phi, theta, psi]

    def apply(self, volume):
        from pytom.agnostic.filter import applyFourierFilter

        wedge = self.getWeightVolume()

        return applyFourierFilter(volume, wedge)

    def getWeightVolume(self, reducedComplex=True):
        from pytom.agnostic.filter import create_wedge
        from pytom.agnostic.transform import fourier_reduced2full

        wedge = create_wedge(self.wedge_angle1, self.wedge_angle2, self.cutOffRadius, self.size_x, self.size_y, self.size_z, self.smooth, self.rotation)

        if reducedComplex == False:
            wedge = fourier_reduced2full(wedge)

        return wedge


class Alignment:
    def __init__(self, vol1, vol2, score, mask=None, iniRot=None, iniTrans=None,
                 opti='fmin_powell', interpolation='linear', verbose=False):
        """
        alignment of a particle against a reference

        @param vol1: (constant) volume
        @type vol1: L{pytom.lib.pytom_volume.vol}
        @param vol2: volume that is matched to reference
        @type vol2: L{pytom.lib.pytom_volume.vol}
        @param score: score for alignment - e.g., pytom.basic.correlation.nxcc
        @type score: L{pytom.basic.correlation}
        @param mask: mask correlation is constrained on
        @type mask: L{pytom.lib.pytom_volume.vol}
        @param iniRot: initial rotation of vol2
        @type iniRot: L{pytom.basic.Rotation}
        @param iniTrans: initial translation of vol2
        @type iniTrans: L{pytom.basic.Shift}
        @param opti: optimizer ('fmin_powell', 'fmin', 'fmin_cg', 'fmin_slsqp', 'fmin_bfgs')
        @param interpolation: interpolation type - 'linear' (default) or 'spline'
        @type interpolation: str
        @type opti: L{str}

        @author: FF
        """
        from pytom.agnostic.normalise import normaliseUnderMask, mean0std1
        from pytom.agnostic.tools import volumesSameSize
        from pytom.voltools import StaticVolume

        assert isinstance(interpolation, str), "interpolation must be of type str"
        self.interpolation = interpolation

        self.verbose = verbose
        if not volumesSameSize(vol1, vol2):
            raise RuntimeError('Vol1 and vol2 must have same size!')

        # normalize constant volume
        if not mask is None:
            (v, p) = normaliseUnderMask(vol1, mask)
        else:
            v = mean0std1(vol1, True)

        self.vol1 = v
        self.vol2 = StaticVolume(xp.array(vol2), interpolation=self.interpolation, device=device)
        self.rotvol2 = xp.zeros_like(self.vol1,dtype=xp.float32)
        self.mask = mask

        if not iniRot:
            iniRot = Rotation()
        if not iniTrans:
            iniTrans = Shift()
        self.rot_trans = self.transRot2vector(rot=iniRot, trans=iniTrans)

        self.score = score
        self.val = -100000.
        self.centX = int(self.vol1.shape[0] // 2)
        self.centY = int(self.vol1.shape[1] // 2)
        self.centZ = int(self.vol1.shape[2] // 2)
        self.binning = 1

        # set optimizer
        self.opti = opti
        if opti == 'fmin':
            self.optimizer = scipy.optimize.fmin
        elif opti == 'fmin_slsqp':
            self.optimizer = scipy.optimize.fmin_slsqp
        elif opti == 'fmin_cg':
            self.optimizer = scipy.optimize.fmin_cg
        elif opti == 'fmin_bfgs':
            self.optimizer = scipy.optimize.fmin_bfgs
        elif opti == 'fmin_powell':
            self.optimizer = scipy.optimize.fmin_powell
        else:
            raise TypeError('opti must be of type str')

    def transRot2vector(self, rot, trans):
        """
        convert rotation and translation to 6-dimensional vector
        @param rot: rotation of vol2
        @type rot: L{pytom.basic.Rotation}
        @param trans: translation of vol2
        @type trans: L{pytom.basic.Shift}
        @return: rotation and translation vector (6-dim: z1,z2,x \
            angles, x,y,z, translations)
        @rtype: L{list}
        @author: FF

        """
        self.rot_trans = 6 * [0.]
        for ii in range(0, 3):
            self.rot_trans[ii] = rot[ii]
        for ii in range(3, 6):
            self.rot_trans[ii] = trans[ii - 3]
        return self.rot_trans

    def vector2transRot(self, rot_trans):
        """
        convert 6-dimensional vector to rotation and translation

        @param rot_trans: rotation and translation vector (6-dim: z1,z2,x \
            angles, x,y,z, translations)
        @type rot_trans: L{list}
        @return: rotation of vol2, translation of vol2
        @rtype: L{pytom.basic.Rotation}, L{pytom.basic.Shift}
        @author: FF
        """
        rot = Rotation()
        trans = Shift()
        for ii in range(0, 3):
            rot[ii] = self.rot_trans[ii]
        for ii in range(3, 6):
            trans[ii - 3] = self.rot_trans[ii]
        return rot, trans

    def set_rot(self, rot):
        """
        set initial rotation

        @param rot: rotation of vol2
        @type rot: L{pytom.basic.Rotation}
        """
        self.rot_trans[0] = rot[0]
        self.rot_trans[1] = rot[1]
        self.rot_trans[2] = rot[2]

    def set_trans(self, trans):
        """
        set initial translation

        @param trans: translation of vol2
        @type trans: L{pytom.basic.Shift}
        """
        self.rot_trans[3] = trans[0]
        self.rot_trans[4] = trans[1]
        self.rot_trans[5] = trans[2]

    def set_searchVol(self, vol1):
        """
        set search volume (vol1 internally)

        @param vol1: search volume
        @type vol1: L{pytom.lib.pytom_volume.vol}

        """
        from pytom.agnostic.normalise import normaliseUnderMask, mean0std1

        if self.mask:
            (self.vol1, p) = normaliseUnderMask(vol1, self.mask)
        else:
            self.vol1 = mean0std1(vol1, True)

    def evalScore(self, rot_trans):
        """
        evaluate score for given rotation and translation - NEGATIVE cc because\
        optimizer MINIMIZES

        @param rot_trans: rotation and translation vector (6-dim: z1,z2,x \
            angles, x,y,z, translations)
        @type rot_trans: L{list}
        @author: FF
        """
        self.vol2.transform(output=self.rotvol2, rotation=(rot_trans[0], rot_trans[2], rot_trans[1]),
                            rotation_order='rzxz', center=(self.centX, self.centY, self.centZ), translation=(
                rot_trans[3] / self.binning, rot_trans[4] / self.binning, rot_trans[5] / self.binning))

        # score is a fucntion passed to the Alignment structure
        self.val = float(-1. * (self.score(volume=self.vol1,
                                           template=self.rotvol2,
                                           mask=self.mask,
                                           volume_is_normalized=True)))

        return self.val

    def localOpti(self, iniRot=None, iniTrans=None):
        """
        @param iniRot: initial rotation of vol2
        @type iniRot: L{pytom.basic.Rotation}
        @param iniTrans: initial translation of vol2
        @type iniTrans: L{pytom.basic.Shift}
        @return: opti_score, opti_rot, opti_trans
        @rtype: L{float}, L{pytom.basic.structures.Rotation}, L{pytom.basic.structures.Shift}
        @author: FF
        """

        if not (type(iniRot) == type(None)):
            self.set_rot(rot=iniRot)
        if not (type(iniTrans) == type(None)):
            self.set_trans(trans=iniTrans)

        if self.verbose:
            # alignment score before optimization
            print("CC before optimization %1.3f" % (-1. * self.evalScore(self.rot_trans)))

        # optimize scoring function
        maxiter = 20
        if (self.opti == 'fmin') or (self.opti == 'fmin_powell'):
            rot_trans = self.optimizer(self.evalScore, self.rot_trans,
                                       xtol=0.05, ftol=0.001, maxiter=maxiter, maxfun=maxiter * 20)
            # xtol=0.0001, ftol=0.0001, maxiter=maxiter, maxfun=None)
        elif self.opti == 'fmin_cg':
            rot_trans = self.optimizer(self.evalScore, self.rot_trans,
                                       gtol=0.0000001,
                                       maxiter=maxiter)
        elif self.opti == 'fmin_slsqp':
            rot_trans = self.optimizer(self.evalScore, self.rot_trans,
                                       iter=maxiter, acc=1e-03)
        elif self.opti == 'fmin_bfgs':
            rot_trans = self.optimizer(self.evalScore, self.rot_trans,
                                       maxiter=maxiter, epsilon=1e-06)
        elif self.opti == 'leastsq':
            rot_trans, success = self.optimizer(self.evalScore, self.rot_trans,
                                                maxfev=maxiter, epsfcn=0.0, factor=10)
        self.rot_trans = rot_trans

        # alignment score before optimization
        finscore = self.evalScore(self.rot_trans)
        rot, trans = self.vector2transRot(rot_trans)
        if self.verbose:
            print("CC after optimization %1.3f" % (-1. * finscore))
            print("rot_trans = ", rot_trans)
            print(rot, trans)

        return -1. * finscore, rot, trans
