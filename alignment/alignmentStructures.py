from pytom.basic.structures import PyTomClass
from pytom.gpu.initialize import xp, device
from pytom.angles.localSampling import LocalSampling


class SamplingParameters(PyTomClass):
    """
    class to store parameters for sampling
    """

    def __init__(self, rotations=LocalSampling(shells=3, increment=3., z1Start=0., z2Start=0., xStart=0.),
                 binning=1, adaptive_res=0.1, sample_info=None):
        """
        @param rotations: rotations for orientation sampling
        @type rotations: L{pytom.angles.angle.AngleObject}
        @param binning: binning FACTOR (1=no binning, 2=2x2x2 voxel-> 1 voxel, etc.)
        @type binning: L{int}
        @param adaptive_res: adaptive resolution offset
        @type adaptive_res: L{float}
        @param sample_info: info about sample
        @type sample_info: L{pytom.basic.structures.SampleInformation}
        @author: FF
        """
        from pytom.angles.angle import AngleObject
        from pytom.basic.structures import SampleInformation

        if rotations == None:
            rotations = LocalSampling(shells=3, increment=3., z1Start=0., z2Start=0., xStart=0.)
        assert type(rotations) == LocalSampling or type(
            rotations) == AngleObject, "SamplingParameters: rotations must be of type LocalSampling"
        self.rotations = rotations
        assert type(binning) == int, "SamplingParameters: binning must be of type int"
        self.binning = binning
        assert type(adaptive_res) == float, "SamplingParameters: adaptive_res must be of type float"
        self.adaptive_res = adaptive_res
        if sample_info == None:
            sample_info = SampleInformation()
        assert type(
            sample_info) == SampleInformation, "SamplingParameters: sample_info must be of type SampleInformation"
        self.sampleInformation = sample_info

    def toXML(self):
        """
        generate XML of Sampling Parameters
        @return: xml object of SamplingParameters
        @rtype: L{etree.Element}
        @author: FF
        """
        from lxml import etree

        xmlObj = etree.Element("SamplingParameters")
        xmlObj.append(self.rotations.toXML())
        xmlObj.append(self.sampleInformation.toXML())
        xmlObj.set("Binning", str(self.binning))
        if self.adaptive_res is False:
            xmlObj.set("AdaptiveResolution", '0')
        else:
            xmlObj.set("AdaptiveResolution", str(self.adaptive_res))
        return xmlObj

    def fromXML(self, xmlObj):
        """
        @param xmlObj: xml object
        @type xmlObj: L{lxml.etree._Element}
        @author: FF
        """
        from lxml.etree import _Element
        from pytom.angles.angle import AngleObject
        from pytom.basic.structures import SampleInformation

        if xmlObj.__class__ != _Element:
            raise Exception('You must provide a valid XML object.')
        if xmlObj.tag == "SamplingParameters":
            samplingDescription = xmlObj
        else:
            samplingDescription = xmlObj.xpath('SamplingParameters')
            if len(samplingDescription) == 0:
                raise Exception("This XML is not a SamplingParameters.")
            samplingDescription = samplingDescription[0]

        self.binning = int(samplingDescription.get('Binning'))
        self.adaptive_res = float(samplingDescription.get('AdaptiveResolution'))

        rot = AngleObject()
        self.rotations = rot.fromXML(samplingDescription.xpath('Angles')[0])

        try:
            si = samplingDescription.xpath('SampleInformation')[0]
            self.sampleInformation = SampleInformation()
            self.sampleInformation.fromXML(si)
        except:
            self.sampleInformation = SampleInformation()


from pytom.score.score import FLCFScore


class ScoringParameters(PyTomClass):
    """
    class to store parameters for subtomo alignment score
    @author: FF
    """

    def __init__(self, score=FLCFScore(), ref=None, weighting=False, compoundWedge=False, mask=None,
                 preprocessing=None, fsc_criterion=0.143, symmetries=None):
        """
        @param mask: mask
        @type mask: L{pytom.basic.structures.Mask}
        @param preprocessing: pre-processing parameters (bandpass filter)
        @type preprocessing: L{pytom.alignment.preprocessing.Preprocessing}
        @param weighting: weighting of particles in average by exp(cc)
        @type weighting: bool
        @param compoundWedge: compoundWedge weighting of reference - i.e. used unweighted average for alignment but \
                                rather weight particle accordingly
        @type compoundWedge: bool
        @param symmetries: apply point symmetry for average
        @type symmetries: ???
        @param fsc_criterion: FSC criterion
        @type fsc_criterion: float

        @author: FF
        """
        if 'cpu' in device:
            from pytom.basic.structures import Mask, Reference, Symmetry
            from pytom.alignment.preprocessing import Preprocessing
            from pytom.score.score import Score, FLCFScore
        else:
            from pytom.tompy.structures import Mask, Reference, Symmetry, Preprocessing
            from pytom.score.score import Score, FLCFScore

        assert type(score) == Score or type(
            score) == FLCFScore, "ScoringParameters: input score not of pytom type Score"
        self.score = score
        assert (type(ref) == Reference) or (type(ref) == type(None)), \
            "ScoringParameters: input ref not of pytom type Reference or None"
        self.reference = ref
        assert (type(mask) == Mask) or (
                    type(mask) == type(None)), "ScoringParameters: input mask not of pytom type Mask"
        self.mask = mask
        assert (type(preprocessing) == Preprocessing) or (type(preprocessing) == type(None)), \
            "ScoringParameters: input preprocessing not of pytom type Preprocessing"
        self.preprocessing = preprocessing
        assert type(fsc_criterion) == float, "ScoringParameters: input fsc_criterion not of type float"
        self.fsc_criterion = fsc_criterion
        assert type(symmetries) == Symmetry or type(symmetries) == type(None)
        self.symmetries = symmetries
        self.weighting = weighting
        self.compoundWedge = compoundWedge

    def toXML(self):
        """
        generate xml node of Scoring Parameters
        @return: xml object of ScoringParameters
        @rtype: L{etree.Element}
        @author: FF
        """
        from lxml import etree

        xmlObj = etree.Element("ScoringParameters")
        xmlObj.append(self.score.toXML())
        if self.reference:
            xmlObj.append(self.reference.toXML(onlyParticleListFilename=True))
        if self.mask:
            xmlObj.append(self.mask.toXML())
        if self.preprocessing:
            xmlObj.append(self.preprocessing.toXML())
        if self.symmetries != None:
            xmlObj.append(self.symmetries.toXML())
        xmlObj.set("Weighting", str(self.weighting))
        xmlObj.set("CompoundWedge", str(self.compoundWedge))
        return xmlObj

    def fromXML(self, xmlObj):
        """
        read ScoringParameters from xml object
        @param xmlObj: xml object
        @type xmlObj: L{lxml.etree._Element}
        @author: FF
        """
        from lxml.etree import _Element
        from pytom.score.score import fromXML as fromXMLScore
        from pytom.basic.structures import Mask, Reference, MultiSymmetries
        from pytom.alignment.preprocessing import Preprocessing

        if xmlObj.__class__ != _Element:
            raise Exception('You must provide a valid XML object.')
        if xmlObj.tag == "ScoringParameters":
            scoringDescription = xmlObj
        else:
            scoringDescription = xmlObj.xpath('SamplingParameters')
            if len(scoringDescription) == 0:
                raise Exception("This XML is not a SamplingParameters.")
            scoringDescription = scoringDescription[0]

        if scoringDescription.get('Weighting'):
            tmp = str(scoringDescription.get('Weighting'))
            if tmp.lower() == 'true':
                self.weighting = True
            else:
                self.weighting = False
        else:
            self.weighting = False

        if scoringDescription.get('CompoundWedge'):
            tmp = str(scoringDescription.get('CompoundWedge'))
            if tmp.lower() == 'true':
                self.compoundWedge = True
            else:
                self.CompoundWedge = False
        else:
            self.CompoundWedge = False

        if scoringDescription.xpath('Score'):
            self.score = fromXMLScore(scoringDescription.xpath('Score')[0])
        else:
            self.score = FLCFScore()

        self.reference = Reference()
        self.reference.fromXML(xmlObj=scoringDescription.xpath('Reference')[0])

        self.mask = Mask()
        if scoringDescription.xpath('Mask'):
            self.mask.fromXML(xmlObj=scoringDescription.xpath('Mask')[0])

        self.preprocessing = Preprocessing()
        if scoringDescription.xpath('Preprocessing'):
            self.preprocessing.fromXML(xmlObj=scoringDescription.xpath('Preprocessing')[0])

        if scoringDescription.xpath('MultiSymmetries'):
            self.symmetries = MultiSymmetries()
            self.symmetries.fromXML(xmlObj=scoringDescription.xpath('MultiSymmetries')[0])
        else:
            self.symmetries = None


class GLocalSamplingJob(PyTomClass):
    """
    Gold standard Local Sampling Job
    """

    def __init__(self, pl=None, ref=None, mask=None, sample_info=None, rotations=None,
                 preprocessing=None, dest='.', max_iter=10, score=FLCFScore(), weighting=False, compoundWedge=False,
                 binning=1, symmetries=None, adaptive_res=0.1, fsc_criterion=0.143, gpuIDs=None):
        """
        @param pl: particle list
        @type pl: L{pytom.basic.structures.ParticleList}
        @param ref: reference
        @type ref: L{pytom.basic.structures.Reference}
        @param mask: mask
        @type mask: L{pytom.basic.structures.Mask}
        @param sample_info: info about sample ?!
        @param rotations: rotations for orientation sampling
        @type rotations: L{pytom.angles.angle.AngleObject}
        @param preprocessing: pre-processing parameters (mostly bandpass filter)
        @type preprocessing: L{pytom.alignment.preprocessing.Preprocessing}
        @param dest: destination folder
        @type dest: L{str}
        @param max_iter: maximum iteration of quasi ExMax alignment
        @type max_iter: L{int}
        @param score: Scoring object
        @type score: L{pytom.score.score.Score}
        @param weighting: Compound wedge weighting of average
        @type weighting: L{bool}
        @param compoundWedge: compound wedge weighting of reference for alignment
        @type compoundWedge: bool
        @param binning: binning FACTOR (1=no binning, 2=2x2x2 voxel-> 1 voxel, etc.)
        @type binning: L{int}
        @param symmetries: apply point symmetry for average
        @type symmetries: ???
        @param adaptive_res: adaptive resolution offset
        @type adaptive_res: L{float}
        @param fsc_criterion: FSC criterion
        @type fsc_criterion: L{float}

        @author: FF
        """
        from pytom.basic.structures import ParticleList, Reference, Mask, MultiSymmetries
        from pytom.score.score import Score
        from pytom.angles.angle import AngleObject
        from pytom.alignment.preprocessing import Preprocessing

        assert (type(pl) == ParticleList) or (type(pl) == type(None)), \
            "GLocalSamplingJob: pl must be particleList or None"
        self.particleList = pl
        assert type(dest) == str, "GLocalSamplingJob: dest must be a string!"
        self.destination = dest
        assert type(max_iter) == int, "GLocalSamplingJob: max_iter must be an integer!"
        self.max_iter = max_iter

        # set scoring parameters
        assert (type(ref) == Reference) or (type(ref) == type(None)), \
            "GLocalSamplingJob: ref must be Reference or None"
        if score is None:
            score = FLCFScore()
        assert type(score) == Score or type(score) == FLCFScore, "GLocalSamplingJob: score is of type Score"
        assert (type(mask) == Mask) or (type(mask) == type(None)), \
            "GLocalSamplingJob: mask is of type Mask"
        if preprocessing is None:
            self.preprocessing = Preprocessing()
        else:
            self.preprocessing = preprocessing
        assert type(self.preprocessing) == Preprocessing, \
            "GLocalSamplingJob: preprocessing is of type Preprocessing"
        assert type(fsc_criterion) == float, "GLocalSamplingJob: fsc_criterion is a float"
        assert type(weighting) == bool, "GLocalSamplingJob: weighting is bool"
        assert (type(symmetries) == MultiSymmetries) or (type(symmetries) == type(None))
        self.scoringParameters = ScoringParameters(score=score, ref=ref, weighting=weighting,
                                                   compoundWedge=compoundWedge, mask=mask,
                                                   preprocessing=preprocessing, fsc_criterion=fsc_criterion,
                                                   symmetries=symmetries)
        if rotations is None:
            from pytom.angles.localSampling import LocalSampling
            rotations = LocalSampling(shells=3, increment=3., z1Start=0., z2Start=0., xStart=0.)
        assert type(rotations) == LocalSampling or type(
            rotations) == AngleObject, "GLocalSamplingJob: rotations is LocalSampling"
        assert type(binning) == int, "GLocalSamplingJob: binning is an int"
        assert type(adaptive_res) == float, "GLocalSamplingJob: adaptive_res is of type float"
        self.samplingParameters = SamplingParameters(rotations=rotations,
                                                     binning=binning, adaptive_res=adaptive_res,
                                                     sample_info=sample_info)
        self.gpu = gpuIDs

        self.gpu = gpuIDs

    def fromXML(self, xmlObj):
        """
        @param xmlObj: xml object
        @type xmlObj: L{lxml.etree.Element}
        @author: FF
        """
        from lxml.etree import _Element
        from pytom.basic.structures import ParticleList

        if xmlObj.__class__ != _Element:
            raise Exception('You must provide a valid XML object.')
        if xmlObj.tag == "GLocalSamplingJob":
            jobDescription = xmlObj
        else:
            jobDescription = xmlObj.xpath('GLocalSamplingJob')
            if len(jobDescription) == 0:
                raise Exception("This XML is not a GLocalSamplingJob.")
            jobDescription = jobDescription[0]

        pl = ParticleList('.')
        particleList_element = jobDescription.xpath('ParticleList')[0]
        filename = particleList_element.get('Filename')
        pl.fromXMLFile(filename=filename)

        self.destination = jobDescription.get('Destination')
        self.max_iter = int(jobDescription.get('MaxIterations'))

        try:
            samplingDescription = jobDescription.xpath('SamplingParameters')[0]
        except IndexError:
            raise Exception("GLocalSamplingJob must contain SamplingParameters")
        self.samplingParameters = SamplingParameters()
        self.samplingParameters.fromXML(xmlObj=samplingDescription)

        try:
            scoringDescription = jobDescription.xpath('ScoringParameters')[0]
        except IndexError:
            raise Exception("GLocalSamplingJob must contain ScoringParameters")
        self.scoringParameters = ScoringParameters()
        self.scoringParameters.fromXML(xmlObj=scoringDescription)

    def toXML(self):
        """
        write Job to XML file

        @author: FF
        """
        from lxml import etree

        jobElement = etree.Element("GLocalSamplingJob")
        if self.particleList:
            particleListElement = etree.Element("ParticleList")
            particleListElement.set("Filename", self.particleList._XMLfilename)
            jobElement.append(particleListElement)
            # if not shortParticleList:
            #    jobElement.append(self.particleList.toXML())
            # else:
            #    jobElement.append(self.particleList.toShortXML())
        jobElement.set("Destination", self.destination)
        jobElement.set("MaxIterations", str(self.max_iter))

        jobElement.append(self.scoringParameters.toXML())
        jobElement.append(self.samplingParameters.toXML())
        return jobElement

    def prettyPrintShortXML(self):
        """
        print job with short XML for particle list
        @author: FF
        """
        from lxml import etree
        return etree.tostring(self.toXML(), pretty_print=True)

    def getParticleList(self):
        """
        get particle list
        @return: particle list
        @rtype: L{pytom.basic.structure.ParticleList}
        @author: FF
        """
        return self.particleList

    @property
    def getParticleListtoXML(self):
        """
        @author: FF
        """
        return self.particleList.toXML()

    def check(self):
        """
        @author: FF

        """
        from pytom.tools.files import checkDirExists
        self.particleList.check()
        self.reference.check()
        self.mask.check()
        if not checkDirExists(self.destination):
            raise RuntimeError('Destination path not found! ' + self.destination)