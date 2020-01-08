'''
Created on May 28, 2010

@author: chen
'''
import unittest

class pytom_LocalTest(unittest.TestCase):
    """
    pytom_LocalTest: Test case for package localization
    """
    def setUp(self):
        """set up"""
        self.testfilename = './testData/emd_1480.map.em_bin_4.em'

    # Module Structures---------------------------------
    def Volume_Test(self):
        from pytom.localization.structures import Volume
        
        filename = self.testfilename
        a = Volume(filename)
        
        self.assertTrue( a.getFilename() == filename, 
            msg='filename not stored in volume structure')
        
        xmlObj = a.toXML()
        b = Volume()
        b.fromXML(xmlObj)
        
        self.assertTrue( b.getFilename() == filename, 
            msg='filename not updated in volume structure')
        
    def Orientation_Test(self):
        from pytom.localization.structures import Orientation
        
        a = Orientation(self.testfilename)
        
        self.assertTrue( a.getFilename() == self.testfilename,
            msg='filename not stored in Orientation structure')
        
        xmlObj = a.toXML()
        b = Orientation()
        b.fromXML(xmlObj)
        
        self.assertTrue( b.getFilename() == self.testfilename,
            msg='filename not stored in copied Orientation structure')
        
    def FoundParticle_Test(self):
        from pytom.localization.structures import FoundParticle
        from pytom.basic.structures import PickPosition, Rotation
        from pytom.score.score import FLCFScore
        
        r = Rotation([1,2,3])
        p = PickPosition([4,5,6], originFilename='originalFilename')
        s = FLCFScore(); s.setValue(0.2)
        
        a = FoundParticle(p, r, s, 'particleFilename')
        xmlObj = a.toXML()
        b = FoundParticle()
        b.fromXML(xmlObj)
        
        self.assertTrue( b.pos.toVector() == [4,5,6],
            msg='position not transferred from xml')
        self.assertTrue( b.pos.getOriginFilename() == 'originalFilename',
            msg='filename not transferred from xml')
        self.assertTrue( b.orient.toList() == [1,2,3],
            msg='orientation not transferred from xml')
        self.assertTrue( float(b.score.getValue()) == 0.2,
            msg='score not transferred from xml')
        self.assertTrue( b.filename == 'particleFilename',
            msg='particle filename not transferred from xml')
        
    # Module PeakJob------------------------------------
    def PeakJob_Test(self):
        from pytom.localization.peak_job import PeakJob
        from pytom.basic.structures import Mask, Reference, WedgeInfo
        from pytom.localization.structures import Volume
        from pytom.score.score import FLCFScore
        from pytom.angles.angleList import AngleList
        
        v = Volume( self.testfilename)
        ref = Reference(self.testfilename)
        m = Mask(self.testfilename)
        w = WedgeInfo(30)
        s = FLCFScore()
    
        rot = AngleList([[1,1,1],[2,2,2],[3,3,3]])
        
        a = PeakJob(v, ref, m, w, rot, s, 1)
        xmlObj = a.toXML()
        b = PeakJob()
        b.fromXML(xmlObj)
        
        self.assertTrue( b.volume.getFilename() == a.volume.getFilename(),
            msg='')
        self.assertTrue( b.jobID == a.jobID,
            msg='')
        self.assertTrue( b.reference.getReferenceFilename() == a.reference.getReferenceFilename(),
            msg='')
        self.assertTrue( b.mask.getFilename() == a.mask.getFilename(),
            msg='')
        self.assertTrue( b.wedge.getWedgeAngle() == a.wedge.getWedgeAngle(),
            msg='')
        self.assertTrue( b.score.getScoreFunc() == a.score.getScoreFunc(),
            msg='')
        
    def PeakResult_Test(self):
        from pytom.localization.peak_job import PeakResult
        from pytom.localization.structures import Volume, Orientation
        
        v = Volume( self.testfilename)
        o = Orientation(self.testfilename)
        a = PeakResult(v, o, 2)
        
        xmlObj = a.toXML()
        b = PeakResult()
        b.fromXML(xmlObj)
        
        self.assertTrue( b.jobID == a.jobID,
            msg='')
        self.assertTrue( b.result.getFilename() == a.result.getFilename(),
            msg='')
        self.assertTrue( b.orient.getFilename() == a.orient.getFilename(),
            msg='')
        
    # Module PeakJobMsg------------------------------------
    def PeakJobMsg_Test(self):
        from pytom.localization.peak_job_msg import PeakJobMsg
        from pytom.localization.peak_job import PeakJob
        from pytom.basic.structures import Mask, Reference, WedgeInfo
        from pytom.localization.structures import Volume
        from pytom.score.score import FLCFScore
        from pytom.angles.angleList import AngleList
        
        v = Volume( self.testfilename)
        ref = Reference(self.testfilename)
        m = Mask(self.testfilename)
        w = WedgeInfo(30)
        s = FLCFScore()
    
        rot = AngleList([[1,1,1],[2,2,2],[3,3,3]])
        
        j = PeakJob(v, ref, m, w, rot, s)
        
        a = PeakJobMsg(str(0), str(1))
        a.setJob(j)
        xmlObj = a.toXML()
        b = PeakJobMsg()
        b.fromXML(xmlObj)
        
        self.assertTrue( b.getSender() == a.getSender(),
            msg='')
        self.assertTrue( b.getRecipient() == a.getRecipient(),
            msg='')
        
    def PeakResultMsg_Test(self):
        from pytom.localization.peak_job_msg import PeakResultMsg
        from pytom.localization.peak_job import PeakResult
        from pytom.localization.structures import Volume, Orientation
        
        v = Volume( self.testfilename)
        o = Orientation(self.testfilename)
        r = PeakResult(v, o)
        
        a = PeakResultMsg(str(1), str(0))
        a.setResult(r)
        xmlObj = a.toXML()
        b = PeakResultMsg()
        b.fromXML(xmlObj)
        
        self.assertTrue( b.getSender() == a.getSender(),
            msg='')
        self.assertTrue( b.getRecipient() == a.getRecipient(),
            msg='')
        
    
    #def BackwardCompatibility_Test(self):
    #    from pytom.localization.peak_job import PeakJob
    #    job = PeakJob()
    #    job.fromXMLFile('./testData/xmlFiles/localizationJob.xml')

        
    
        
def local_TestSuite():
    """
    local_TestSuite: Create a test suite for package localization
    """
    from unittest import TestSuite
    suite = TestSuite()
    
    suite.addTest(pytom_LocalTest("Volume_Test"))
    suite.addTest(pytom_LocalTest("Orientation_Test"))
    suite.addTest(pytom_LocalTest("FoundParticle_Test"))
    suite.addTest(pytom_LocalTest("PeakJob_Test"))
    suite.addTest(pytom_LocalTest("PeakResult_Test"))
    suite.addTest(pytom_LocalTest("PeakJobMsg_Test"))
    suite.addTest(pytom_LocalTest("PeakResultMsg_Test"))
    
    return suite

if __name__ == '__main__':
    
    suite = local_TestSuite()
    
    from unittest import TextTestRunner
    runner = TextTestRunner()
    runner.run(suite)
