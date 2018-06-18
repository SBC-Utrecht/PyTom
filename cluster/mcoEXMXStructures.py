'''
Created on Dec 22, 2010

@author: hrabe
'''

from pytom.basic.structures import PyTomClass
from pytom.alignment.ExMaxAlignment import ExMaxWorker

class ClusterSwap(PyTomClass):
    """
    ClusterSwap: Saves conditions under which a class swap of a particle happened.
    """
    
    def __init__(self,particle=None,oldClusterName=None,oldScore=None,newClusterName=None,newScore=None):
        
        self._particle = particle
        self._oldClusterName = oldClusterName
        self._oldScore = oldScore
        self._newClusterName = newClusterName
        self._newScore = newScore
    
    def toXML(self):    
        
        from lxml import etree
        
        swapElement = etree.Element('ClassSwap',OldClusterName = str(self._oldClusterName),NewClusterName=str(self._newClusterName))
        
        swapElement.append(self._particle.toXML())
        
        os = self._oldScore.toXML()
        os.set('OriginalType',str(os.tag))
        os.tag = 'OldScore'
        swapElement.append(os)
        
        ns = self._newScore.toXML()
        ns.set('OriginalType',str(ns.tag))
        ns.tag = 'NewScore'
        swapElement.append(ns)
    
        return swapElement
    
    def fromXML(self,xmlObj):
        
        from lxml.etree import _Element
        from pytom.basic.structures import Particle
        from pytom.score.score import fromXML as scoreFromXML
        
        if xmlObj.__class__ != _Element :
            raise Exception('Is not a lxml.etree._Element! You must provide a valid XMLobject.')
        
        
        self._oldClusterName = str(xmlObj.get('OldClusterName')) 
        self._newClusterName = str(xmlObj.get('NewClusterName'))
        
        oldScore = xmlObj.xpath('OldScore')[0]
        oldScore.tag = oldScore.get('OriginalType')
        self._oldScore = scoreFromXML(oldScore)
        
        newScore = xmlObj.xpath('NewScore')[0]
        newScore.tag = newScore.get('OriginalType')
        self._newScore = scoreFromXML(newScore)
        
        particle = xmlObj.xpath('Particle')[0]
        self._particle = Particle('')
        self._particle.fromXML(particle)

class SwapList(PyTomClass):      
  
    def __init__(self,swaps = None):
        
        if swaps != None and len(swaps) > 0:
            for swap in swaps:
                if not swap.__class__ == ClusterSwap:
                    raise TypeError('All list elements must conform to the ClusterSwap class')
        
            self._swapList = swaps
        else:
            self._swapList = []
            
            
    def append(self,swap):
    
        self._swapList.append(swap)
       
    def toXML(self):
        from lxml import etree
        
        listElement = etree.Element('ClusterSwapList')
        
        for swap in self._swapList:
            listElement.append(swap.toXML())
        
        return listElement
  
    def fromXML(self,xmlObj):
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element :
            raise Exception('Is not a lxml.etree._Element! You must provide a valid XMLobject.')
        
        swaps = xmlObj.xpath('ClassSwap')
        
        self._swapList = []
        
        for swap in swaps:
            s = ClusterSwap()
            s.fromXML(swap)
            self._swapList.append(s)
        
        
        
class MCOEXMXJob(PyTomClass):
    
    """
    MCOEXMXJob: Job description for deterministic classification K-Means style.
    """
    
    def __init__(self,particleList=None,numberIterations=None,destinationDirectory=None,mask=None,score=None,preprocessing=None,wedgeInfo=None,binning=None,sampleInformation=None,numberClasses=None,endThreshold=None,symmetry = None,doAlignment = False,frmBandwidth = None):
        """
        __init__:
        @param particleList: The particles
        @type particleList: L{pytom.basic.structures.ParticleList}
        @param numberIterations: Number of kmeans rounds
        @type numberIterations: int
        @param destinationDirectory: Where will results be stored? 
        @type destinationDirectory: str
        @param mask: Score mask
        @type mask: L{pytom.basic.structures.Mask}
        @param score: The score
        @type score: Any child of L{pytom.alignment.score.Score}
        @param preprocessing: Bandpass and other things here 
        @type preprocessing: L{pytom.alignment.preprocessing.Preprocessing}
        @param wedgeInfo: 
        @type wedgeInfo: L{pytom.basic.structures.WedgeInfo}
        @param binning: Binning factor. 
        @type binning:
        @param sampleInformation:   
        @type sampleInformation:
        @param numberClasses: Optional number classes. Particle List will be radomised if no predifined class membership is provided.
        @type numberClasses: int
        @param endThreshold: End threshold of how many particles at least should jump to continue (value as fraction between 0 and 1)  
        @type endThreshold:  float
        @param symmetry: Set symmetry for particles, default is none
        @type symmetry: Any child of L{pytom.basic.structures.Symmetry} 
        @param doAlignment: Include FRM alignment in classification. Default is False.
        @type doAlignment: boolean 
        @param frmBandwidth: The bandwidth of the spherical harmonics - lowestBand used, highestBand used 
        @type frmBandwidth: [lowestBand,highestBand]
        """
        from pytom.basic.structures import PointSymmetry,HelicalSymmetry,Reference
        from pytom.angles.angleList import OneAngleList
        from pytom.alignment.ExMaxAlignment import ExMaxJob
        
        #self._particleList = particleList
        if numberIterations:
            self._numberIterations = int(numberIterations)
        else:
            self._numberIterations = 0
        if destinationDirectory:
            self._destinationDirectory = str(destinationDirectory)
        else:
            self._destinationDirectory = None
            
        self._numberClasses = numberClasses
        self._endThreshold = endThreshold
        
        rotations = OneAngleList()
        rotations.append(0,0,0)
        
        if not symmetry.__class__ in [PointSymmetry,HelicalSymmetry]:
            symmetry = PointSymmetry(1)
        
        
        self._exMaxJob = ExMaxJob(particleList,self._destinationDirectory,Reference(''),score,rotations,mask,
                                                          symmetry,0,1,preprocessing,
                                                          -1.0,binning,sampleInformation,0.5,adaptiveResolution=False)
        
        self._doAlignment = doAlignment
        self._frmBandwidth = frmBandwidth
    
    def getAdaptiveResolution(self):
        return self._exMaxJob.getAdaptiveResolution()
    
    def setAdaptiveResolution(self,adaptiveResolution):
        self._exMaxJob.setAdaptiveResolution(adaptiveResolution)
    
    def getFRMBandwidth(self):
        return self._frmBandwidth
    
    def setFRMBandwidth(self,value):
        assert value.__class__ == list and len(list) == 2
        
        self._frmBandwidth = value
        
    def getDoAlignment(self):
        return self._doAlignment
    
    def setDoAlignment(self,value):
        assert value.__class__ == bool
        
        self._doAlignment = value
        
    def getParticleList(self):
        #return self._particleList
        return self._exMaxJob.getParticleList()
    
    def getDestinationDirectory(self):
        return self._destinationDirectory
    
    def getNumberIterations(self):
        return self._numberIterations
    
    def getExMaxJob(self):
        return self._exMaxJob
    
    def getNumberClasses(self):
        return self._numberClasses
    
    def getEndThreshold(self):
        return self._endThreshold
    
    def setDestinationDirectory(self,directory):
        from pytom.tools.files import checkDirExists
        
        if not checkDirExists(directory):
            from os import mkdir
            mkdir(directory)
        
        self._destinationDirectory = directory
    
    def setParticleList(self,particleList):
        #self._particleList = particleList
        self._exMaxJob.setParticleList(particleList)
        
    def setNumberIterations(self,numberIterations):
        self._numberIterations = numberIterations
        
    def toXML(self):
        
        from lxml import etree
        
        job_element = etree.Element('MCOEXMXJob')
        job_element.append(self._exMaxJob.toXML())
        
        #job_element.append(self._particleList.toXML())
        job_element.set('NumberIterations',str(self._numberIterations))
        job_element.set('DestinationDirectory',self._destinationDirectory)
        job_element.set('NumberClasses',str(self._numberClasses))
        job_element.set('EndThreshold',str(self._endThreshold))
        job_element.set('FRMAlignment',str(self._doAlignment))
        job_element.set('FRMBandwidth',str(self._frmBandwidth))
        
        return job_element
        
    def fromXML(self,xmlObj):
        
        from lxml.etree import _Element;
        
        if xmlObj.__class__ != _Element :
            raise Exception('Is not a lxml.etree._Element! You must provide a valid XMLobject.');
        
        from pytom.basic.structures import ParticleList
        from pytom.alignment.ExMaxAlignment import ExMaxJob
               
        self._numberIterations = int(float(xmlObj.get('NumberIterations')))
        self._destinationDirectory = str(xmlObj.get('DestinationDirectory'))
        self._numberClasses = int(float(xmlObj.get('NumberClasses')))
        self._endThreshold = float(xmlObj.get('EndThreshold'))
        
        
        if xmlObj.get('FRMAlignment') == []:
            self._doAlignment = False
        else:
            self._doAlignment = xmlObj.get('FRMAlignment') == 'True'
        
        
        if not xmlObj.get('FRMBandwidth') or xmlObj.get('FRMBandwidth') == [] or xmlObj.get('FRMBandwidth') == 'None':
            self._frmBandwidth = None
        else:
            self._frmBandwidth = [int(i) for i in xmlObj.get('FRMBandwidth')[1:-1].split(',')]
            
        ex = xmlObj.xpath('ExpectationMaximisationJob')
        self._exMaxJob = ExMaxJob(0,0,0,0,0,0,0,0,0,0,0,0)
        self._exMaxJob.fromXML(ex[0])


class MCOEXMXWorker(ExMaxWorker):        
    """
    MCOEXMXWorker: Will perform all operations neccessary to compare two particles. 
    I.e, determines the class label of a particle according to a class distribution for a certain number of rotations.
    It is driven mainly by job objects defining the next operation. Accepted job structures are L{pytom.alignment.structures.MaximisationJob} or L{pytom.alignment.structures.ExpectationJob} 
    @author: Thomas Hrabe
    """
    
    def getDoAlignment(self):
        return self._doAlignment
    
    def setDoAlignment(self,value):
        assert value.__class__ == bool
        
        self._doAlignment = value
    
    def getFRMBandwidth(self):
        return self._frmBandwidth
    
    def setFRMBandwidth(self,value):
        assert value.__class__ == list and len(list) == 2
        
        self._frmBandwidth = value
    
    def _expectation(self):
        """
        _expectation : Compares two particles. Does not do any alignemnt unless job enables FRMAlignment, this function overloads ExMaxWorker._expectation.
        @rtype: List of two objects
        @return: [endShift,endRotation] contains the final shift and final rotation.        
        @author: Thomas Hrabe
        """
        from pytom_volume import vol,shift,rotate,read
        from pytom.alignment.structures import Peak,MaximisationResult
        from pytom.basic.structures import Particle,Reference,Shift,Rotation,Wedge
        from pytom.tools.macros import volumesSameSize
        from pytom.alignment.alignmentFunctions import compareTwoVolumes
        from pytom.angles.angleList import OneAngleList
        from pytom.alignment.alignmentFunctions import FRMAlignmentWrapper    
        
        try:
            from frm import frm_align_vol2
        except ImportError:
            self._doAlignment = False
            
        if self._particle.__class__ == Particle:
            from pytom_volume import read
            particleFile    = self._particle.getFilename()
            particle        = read(particleFile,0,0,0,0,0,0,0,0,0,self._binning,self._binning,self._binning)
        else:
            raise TypeError('Provided particle must be of type pytom.basic.structures.Particle')
            
        if self._scoreObject.__class__ == str:
            from pytom.score.score import Score
            string = self._scoreObject
            self._scoreObject = Score() 
            self._scoreObject.fromStr(string)
            
        if self._reference.wasGeneratedBy(self._particle) and self._reference.hasPreWedge() and self._reference.hasWeighting():
            #substract particle from reference if information is available
            [reference, self._referenceWeighting] = self._reference.subtractParticle(self._particle,self._binning)
            referenceObject = self._reference
        else:
            if self._reference.__class__ == Reference:
                from pytom_volume import read
                referenceObject = self._reference
                referenceFile   = self._reference.getReferenceFilename()
                reference = read(referenceFile,0,0,0,0,0,0,0,0,0,self._binning,self._binning,self._binning)
            
            if self._referenceWeighting == str and len(self._referenceWeighting) >0:
                self._referenceWeightingFile        = self._referenceWeighting
                
                if binning == 1:
                    from pytom_volume import read
                    self._referenceWeighting        = read(self._referenceWeighting)
                else:
                    from pytom.basic.files import readSubvolumeFromFourierspaceFile
                    self._referenceWeighting        = readSubvolumeFromFourierspaceFile(self._referenceWeightingFile,reference.sizeX(),reference.sizeY(),reference.sizeZ())
        
        
        self._mask.setBinning(self._binning)
        
        pScore = self._particle.getScore()
        
        if not pScore:
            pScore = self._scoreObject


        if self._doAlignment:
            #process FRM alignment for each particle and class 
            peakPrior = self._scoreObject.getPeakPrior()
            
            newShift, newRotation, score = FRMAlignmentWrapper(particle,self._particle.getWedge() , reference, Wedge([0,0]),self._frmBandwidth, int(self._preprocessing.getHighestFrequency()), self._mask, self._scoreObject.getPeakPrior())

            self._scoreObject.setValue(score)

            result = MaximisationResult(self._particle,referenceObject,self._scoreObject,newShift,newRotation,OneAngleList(self._particle.getRotation()))
        else:
            #process classification without alignment only
            peakV = compareTwoVolumes(particle,reference,self._referenceWeighting,self._particle.getWedge() ,OneAngleList(self._particle.getRotation()),self._particle.getShift(),self._scoreObject,self._mask,self._preprocessing,self._binning)

            self._scoreObject.setValue(peakV)
        
            result = MaximisationResult(self._particle,referenceObject,self._scoreObject,self._particle.getShift(),self._particle.getRotation(),OneAngleList(self._particle.getRotation()))
            
        return result

class KMDataSet:
    """
    KMDataSet: K-means data set
    """
    def __init__(self, data=[]):
        """
        @param data: list of data, must have the same length
        @type data: list
        """
        self.data = data
        self.refresh()
        
        if not self.verify():
            raise Exception('Not valid data set!')
        
        self.clusters=[]
    
    def refresh(self):
        self.num = len(self.data)
        if self.num > 0:
            self.dim = len(self.data[0])
        else:
            self.dim = 0
    
    def verify(self):
        """
        verify: verify whether the input data set is valid or not
        """
        for d in self.data:
            if len(d) != self.dim:
                return False
        
        return True
    
    def getData(self, n):
        """
        getData: get the data according to the index
        @param n: data index
        @type n: integer
        @return: data
        """
        if n < 0 or n >= self.num:
            return None
        else:
            return self.data[n]
    
    def initializeClusters(self, k, centroids=None):
        """
        initializeClusters: initialize the clusters
        @param k: num of clusters
        @type k: integer
        """
        if not centroids:
            # not given, randomize the centroids firstly
            from random import uniform
            self.min = []
            self.max = []
            
            for i in xrange(self.dim):
                tmp = [d[i] for d in self.data]
                self.min.append(min(tmp))
                self.max.append(max(tmp))
            
            centroids = []
            for i in xrange(k):
                c = []
                for j in xrange(self.dim):
                    c.append(uniform(self.min[j], self.max[j]))
                centroids.append(c)
        
        for c in centroids:
            self.clusters.append(Cluster(self, c))

    def findSrcCluster(self, i):
        """
        findSrcCluster: find the source cluster of a specific data
        @param i: the index of the data
        @type i: integer
        @return: the cluster
        @rtype: {pytom.cluster.kmeans.Cluster}
        """
        for c in self.clusters:
            if i in c.members:
                return c
        else:
            return None

    def findDestCluster(self, i, proximity_fnc):
        """
        findDestCluster: find the destination cluster where the data should go according to the proximity function
        @param i: the index of the data
        @type i: integer
        @param proximity_fnc: the proximity function
        @type proximity_fnc: function
        @return: the cluster
        @rtype: {pytom.cluster.kmeans.Cluster}
        """
        # find the cloest centroid
        d = self.getData(i)
        dists = []
        for c in self.clusters:
            dists.append(proximity_fnc(d, c.centroid))
            
        return self.clusters[dists.index(min(dists))]
    
    def assignCluster(self, i, proximity_fnc):
        """
        assignCluster: assign a data to a specific cluster according to the proximity function with the cluster centroids
        @param i: the index of the cluster
        @type i:  i: integer
        @param proximity_fnc: the proximity function
        @type proximity_fnc: function
        @return: a boolean value indicating whether the data's cluster membership has been changed during this procedure or not
        @rtype: boolean
        """
        src = self.findSrcCluster(i)
        dest = self.findDestCluster(i, proximity_fnc)
        if src != dest:
            if not src:
                dest.add(i)
            else:
                src.remove(i)
                dest.add(i)
            return True
        else:
            return False
    
    def updateClusters(self):
        """
        updateClusters: update the clusters' centroids
        """
        for c in self.clusters:
            c.updateCentroid()
    
    def findEmptyClusters(self):
        res = []
        for i in xrange(len(self.clusters)):
            if self.clusters[i].isEmpty():
                res.append(i)
        return res
    
    def removeCluster(self, n):
        self.clusters.pop(n)
    
    def addCluster(self, c):
        self.clusters.append(c)
    
    def sortClusters(self, proximity_fnc):
        def cmp(x,y):
            if x.getSSE(proximity_fnc)-y.getSSE(proximity_fnc) < 0:
                return -1
            else:
                return 1
        
        self.clusters.sort(cmp)
    
    def splitCluster(self, n):
        """
        splitCluster: split the cluster into two sub-clusters with lower SSEs
        """
        from pytom.cluster.kmeans import clusterKMeans
        c = self.clusters[n]
        data = c.retrieveAll()
        res = clusterKMeans(data, 2, scale=False, plot=False) # do not scale
        a = Cluster(self, res.clusters[0].centroid)
        b = Cluster(self, res.clusters[1].centroid)
        for i in res.clusters[0].members:
            a.add(c.members[i])
        for i in res.clusters[1].members:
            b.add(c.members[i])
        
        self.removeCluster(n)
        self.clusters.extend([a,b])
        self.updateClusters()
    
    def getTotalSSE(self, proximity_fnc):
        res = 0.
        for c in self.clusters:
            res = res + c.getSSE(proximity_fnc)
        return res
        
    def getAllDataIndexWithoutCluster(self):
        res = range(self.num)
        for c in self.clusters:
            for n in c.members:
                res.remove(n)
        return res
    
class Cluster():
    """
    Cluster: the cluster associated with the K-means data set
    """
    def __init__(self, dataset, centroid, name=None):
        """
        @param dataset: k-means data set
        @type dataset: {pytom.cluster.kmeans.KMDataSet}
        @param centroid: the centroid of this cluster
        @type centroid: list
        """
        self.dataset = dataset
        self.centroid = centroid
        self.members = []
        self.name = name
    
    def add(self, n):
        """
        add: add a new data into this cluster
        @param n: data index
        @type n: integer
        """
        self.members.append(n)
        
    def remove(self, n):
        """
        remove: remove a data out of this cluster
        @param n: data index
        @type n: integer
        """
        if n in self.members:
            self.members.remove(n)
    
    def retrieve(self, n):
        """
        retrieve: retrieve the data contained in this cluster. If the data is not contained in this cluster, return None.
        @param n: data index
        @type n: integer
        @return: data
        @rtype: list
        """
        if n in self.members:
            return self.dataset.getData(n)
        else:
            return None
    
    def retrieveAll(self):
        """
        retrieveAll: retrieve all the data in this cluster
        @return: list of data
        @rtype: list of list
        """
        res = []
        self.members.sort()
        for i in self.members:
            res.append(self.dataset.getData(i))
        return res
    
    def updateCentroid(self):
        """
        updateCentroid: update the centroid of the cluster
        """
        if self.members == []:
            return
        dim = self.dataset.dim
        sum = [0. for i in xrange(dim)]
        for n in self.members:
            d = self.dataset.getData(n)
            sum = [sum[i]+d[i] for i in xrange(dim)]
        self.centroid = [i/len(self.members) for i in sum]
    
    def getSSE(self, proximity_fnc):
        """
        getSSE: get the sum of squared error of this cluster
        @param proximity_fnc: proximity function
        @type proximity_fnc: function
        @return: the sum of squared error
        @rtype: float
        """
        sse = 0.
        self.radius = 0.
        for n in self.members:
            d = self.dataset.getData(n)
            dist = proximity_fnc(d, self.centroid)
            sse = sse + dist
            if dist > self.radius:
                self.radius = dist
        return sse
    
    def isEmpty(self):
        if len(self.members) == 0:
            return True
        else:
            return False
