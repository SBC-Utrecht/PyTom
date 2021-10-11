'''
Created on May 19, 2010

@author: chen
'''

from pytom.basic.structures import PyTomClass

class PeakJob(PyTomClass):
    """
    PeakJob: stores all the infos needed for calculation of the peak score
    """
    
    def __init__(self, volume='', reference='', mask='', wedge='', rotations='', score='', jobID=0, members=1, dstDir='./', bandpass=None):
        """
        @param volume: target volume
        @type volume: L{pytom.localization.structures.Volume}
        @param reference: reference volume
        @type reference: L{pytom.basic.structures.Reference}
        @param mask: mask volume
        @type mask: L{pytom.basic.structures.Mask}
        @param wedge: wedge information
        @type wedge: L{pytom.basic.structures.WedgeInfo}
        @param rotations: rotation list
        @type rotations: L{pytom.angles.angle}
        @param score: score function
        @type score: L{pytom.score.score}
        @param jobID: job identification
        @type jobID: integer
        @param members: how many members are there available to accomplish this job (1 means only itself)
        @type members: integer
        @param dstDir: destination directory where the result is written to
        @type dstDir: string
        @param bandpass: bandpass object that will be applied to the reference
        @type bandpass: L{pytom.basic.structure.BandPassFilter}
        """
        self.volume = volume
        self.reference = reference
        self.mask = mask
        self.wedge = wedge
        self.rotations = rotations
        self.score = score
        self.jobID = jobID
        self.members = members
        if dstDir[-1] == '/':
            self.dstDir = dstDir
        else:
            self.dstDir = dstDir + '/'
        self.bandpass = bandpass
            
    def copy(self, fromJob):
        self.volume = fromJob.volume
        self.reference = fromJob.reference
        self.mask = fromJob.mask
        self.wedge = fromJob.wedge
        self.rotations = fromJob.rotations
        self.score = fromJob.score
        self.jobID = fromJob.jobID
        self.members = fromJob.members
        self.dstDir = fromJob.dstDir
        self.bandpass = fromJob.bandpass
        
    def fromXML(self, xmlObj):
        """
        fromXML : Assigns values to job attributes from XML object
        @param xmlObj: A xml object  
        @author: chen 
        """
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element:
            raise Exception('You must provide a valid XML object.')
        
        if xmlObj.tag == "JobDescription":
            jobDescription = xmlObj
        else:  
            jobDescription = xmlObj.xpath('JobDescription')
            
            if len(jobDescription) == 0:
                raise Exception("This XML is not an JobDescription.")
            
            jobDescription = jobDescription[0]
        
        id = jobDescription.get('ID')
        if id != None and id != 'None':
            self.jobID = int(id)
            
        members = jobDescription.get('Members')
        if members != None and members != 'None':
            self.members = int(members)
        
        dstDir = jobDescription.get('Destination')
        if dstDir != None:
            if dstDir[-1] == '/':
                self.dstDir = dstDir
            else:
                self.dstDir = dstDir + '/'

        from pytom.basic.score import fromXML as fromXMLScore
        from pytom.basic.structures import Mask,Reference,Wedge
#        from pytom.angles.angleList import AngleList
        from pytom.localization.structures import Volume
        
        e = jobDescription.xpath('Volume')[0]
        v = Volume()
        v.fromXML(e)
        self.volume = v
        
        ref = jobDescription.xpath('Reference')[0]
        self.reference = Reference('')
        self.reference.fromXML(ref)
        
        wedgeXML = jobDescription.xpath('Wedge')
        
        if len(wedgeXML) == 0:
            wedgeXML = jobDescription.xpath('SingleTiltWedge')
            
        if len(wedgeXML) == 0:
            wedgeXML = jobDescription.xpath('WedgeInfo')
        
        if len(wedgeXML) == 0:
            wedgeXML = jobDescription.xpath('DoubleTiltWedge')
        
        assert len(wedgeXML) > 0
                
        self.wedge = Wedge()
        self.wedge.fromXML(wedgeXML[0])
        
        
        mask = jobDescription.xpath('Mask')[0]
        self.mask = Mask('')
        self.mask.fromXML(mask)
        
        score = jobDescription.xpath('Score')
        self.score = fromXMLScore(score[0])
        
        rot = jobDescription.xpath('Angles')[0]
        from pytom.angles.angle import AngleObject
        ang = AngleObject()
        self.rotations = ang.fromXML(rot)
#        self.rotations = AngleList()
#        self.rotations.fromXML(rot)
        
        bp = jobDescription.xpath('BandPassFilter')
        if bp != []:
            bp = bp[0]
            from pytom.basic.structures import BandPassFilter
            self.bandpass = BandPassFilter(0,0,0)
            self.bandpass.fromXML(bp)
        else:
            self.bandpass = None
        
    def toXML(self):        
        """
        toXML : Compiles a XML file from job object
        @author: chen
        """    
        
        from lxml import etree
             
        jobElement = etree.Element("JobDescription", ID=str(self.jobID), Members=str(self.members), Destination=str(self.dstDir))
        
        jobElement.append(self.volume.toXML())
        jobElement.append(self.reference.toXML())
        jobElement.append(self.mask.toXML())
        jobElement.append(self.wedge.toXML())
        jobElement.append(self.rotations.toXML())
        jobElement.append(self.score.toXML())
        if self.bandpass:
            jobElement.append(self.bandpass.toXML())
        
        return jobElement
    
    def check(self): 
        """
        check: Performs check whether all settings are valid. Paths and Files exist
        @author: chen 
        """
        
        from pytom.tools.files import checkFileExists,checkDirExists
        
        returnValue = checkFileExists(self.volume.getFilename())
        if not returnValue:
            raise IOError('File: ' + str(self.volume) + ' not found!')
            
        returnValue = checkFileExists(self.reference.getReferenceFilename())
        if not returnValue:
            raise IOError('File: ' + str(self.reference) + ' not found!')
        
        returnValue = checkFileExists(self.mask.getFilename())
        if not returnValue:
            raise IOError('File: ' + str(self.mask) + ' not found!')
        
        returnValue = checkDirExists(self.dstDir[:-1])
        if not returnValue:
            raise IOError('Directory: ' + str(self.dstDir) + ' not found!')
        
        return returnValue

    def send(self, source, destination):
        """
        send: Send the job-relevant message from source to destination
        @param source: source machine id gained from pytom_mpi
        @type source: int
        @param destination: destination machine id
        @type destination: int
        @author: chen
        """
        
        from pytom.localization.peak_job_msg import PeakJobMsg
        
#        self.check()
        msg = PeakJobMsg(str(source), str(destination))
        msg.setJob(self)
        
        import pytom_mpi
        print(f'destination: {destination}\ntype: {type(destination)}')
        pytom_mpi.send(str(msg), int(destination))


class PeakResult(PyTomClass):
    """
    PeakResult: stores the result volume information of the calculation of the peak from the workers
    """
    
    def __init__(self, result='', orient='', jobID=0):
        """
        @param result: result volume
        @type result: L{pytom.localization.structures.Volume}
        @param orient: volume storing the orientation information
        @type orient: L{pytom.localization.structures.Orientation}
        @param jobID: identify the result corresponds to which job
        @type jobID: integer
        """
        self.result = result
        self.orient = orient
        self.jobID = jobID
        
    def fromXML(self, xmlObj):
        """
        fromXML : Assigns values to result attributes from XML object
        @param xmlObj: A xml object  
        @author: chen
        """
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element:
            raise Exception('Is not a lxml.etree._Element! You must provide a valid XML object.'); 
        
        if xmlObj.tag == "Result":
            result = xmlObj
        else:  
            result = xmlObj.xpath('Result')
            
            if len(result) == 0:
                raise Exception("This XML is not an PeakResult. No Result provided.")
            
            result = result[0];
            
        self.jobID = int(result.get('ID'))
        
        from pytom.localization.structures import Volume, Orientation
        e = result.xpath('Volume')[0]
        p = Volume()
        p.fromXML(e)
        self.result = p
        
        e = result.xpath('Orientation')[0]
        o = Orientation()
        o.fromXML(e)
        self.orient = o
        
    def toXML(self):
        """
        toXML : Compiles a XML from result object
        @author: chen
        """    
        from lxml import etree;
        
        resultElement = etree.Element("Result", ID=str(self.jobID));
        
        resultElement.append(self.result.toXML());
        resultElement.append(self.orient.toXML());
        
        return resultElement;
    
    def check(self): 
        """
        check: Performs check whether all settings are valid. Paths and Files exist
        @author: chen 
        """
        from pytom.tools.files import checkFileExists
        
        returnValue = checkFileExists(self.result.getFilename())
        if not returnValue:
            raise IOError(str(self.result) + ' not found!')
        
        returnValue = checkFileExists(self.orient.getFilename())
        if not returnValue:
            raise IOError(str(self.orient) + ' not found!')
         
        return returnValue
    
    def send(self, source, destination):
        
        from pytom.localization.peak_job_msg import PeakResultMsg
        
        msg = PeakResultMsg(str(source), str(destination))
        msg.setResult(self)
        
        import pytom_mpi
        print(f'destination: {destination}\ntype: {source}')
        pytom_mpi.send(str(msg), int(destination))


class JobInfo():
    """
    JobInfo: Class for storing the job information
    @param jobID: current job id
    @type jobID: integer
    @param originalJobID: the father job id
    @type originalJobID: integer
    @param splitType: split type of the job
    @type splitType: "Ang" or "Vol"
    """
    def __init__(self, jobID=-1, originalJobID=-1, splitType=None):
        self.jobID = jobID
        self.originalJobID = originalJobID
        self.splitType = splitType


if __name__ == '__main__':
    from pytom.basic.structures import Mask, Reference, WedgeInfo
    from pytom.localization.structures import Volume
    from pytom.basic.score import FLCFScore
    from pytom.localization.peak_job import PeakJob
    
    v = Volume('/fs/pool/pool-foerster/apps/src/molmatch/test/testvol.em')
    ref = Reference('/fs/pool/pool-foerster/apps/src/molmatch/test/templ.em')
    m = Mask('/fs/pool/pool-foerster/apps/src/molmatch/test/mask_15.em', True)

    w = WedgeInfo(0)
    s = FLCFScore()
    
    from pytom.angles.globalSampling import GlobalSampling
    r = GlobalSampling('/fs/home/ychen/develop/pytom/trunk/pytom/pytomc/libs/libtomc/data/common_data/angles_90_26.em')
    
    job = PeakJob(v, ref, m, w, r, s)
    job.toXMLFile('JobInfo.xml')
    job.toHTMLFile('JobInfo.html')
        
        
        
