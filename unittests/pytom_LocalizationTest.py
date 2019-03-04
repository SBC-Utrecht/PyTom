'''
Created on May 28, 2010

@author: chen
'''
import unittest

class pytom_LocalTest(unittest.TestCase):
    """
    pytom_LocalTest: Test case for package localization
    """
    # Module Structures---------------------------------
    def Volume_Test(self):
        from pytom.localization.structures import Volume
        
        filename = './testData/ribo.em'
        a = Volume(filename)
        
        assert a.getFilename() == filename
        
        xmlObj = a.toXML()
        b = Volume()
        b.fromXML(xmlObj)
        
        assert b.getFilename() == filename
        
    def Orientation_Test(self):
        from pytom.localization.structures import Orientation
        
        filename = './testData/ribo.em'
        a = Orientation(filename)
        
        assert a.getFilename() == filename
        
        xmlObj = a.toXML()
        b = Orientation()
        b.fromXML(xmlObj)
        
        assert b.getFilename() == filename
        
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
        
        assert b.pos.toVector() == [4,5,6]
        assert b.pos.getOriginFilename() == 'originalFilename'
        assert b.orient.toList() == [1,2,3]
        assert float(b.score.getValue()) == 0.2
        assert b.filename == 'particleFilename'
        
    # Module PeakJob------------------------------------
    def PeakJob_Test(self):
        from pytom.localization.peak_job import PeakJob
        
        from pytom.basic.structures import Mask, Reference, WedgeInfo
        from pytom.localization.structures import Volume
        from pytom.score.score import FLCFScore
        
        v = Volume('./testData/ribo.em')
        ref = Reference('./testData/ribo.em')
        m = Mask('./testData/ribo.em')
        w = WedgeInfo(30)
        s = FLCFScore()
    
        from pytom.angles.angleList import AngleList
        rot = AngleList([[1,1,1],[2,2,2],[3,3,3]])
        
        a = PeakJob(v, ref, m, w, rot, s, 1)
        xmlObj = a.toXML()
        b = PeakJob()
        b.fromXML(xmlObj)
        
        assert b.volume.getFilename() == a.volume.getFilename()
        assert b.jobID == a.jobID
        assert b.reference.getReferenceFilename() == a.reference.getReferenceFilename()
        assert b.mask.getFilename() == a.mask.getFilename()
        assert b.wedge.getWedgeAngle() == a.wedge.getWedgeAngle()
        assert b.score.getScoreFunc() == a.score.getScoreFunc()
        
    def PeakResult_Test(self):
        from pytom.localization.peak_job import PeakResult
        from pytom.localization.structures import Volume, Orientation
        
        v = Volume('./testData/ribo.em')
        o = Orientation('./testData/ribo.em')
        a = PeakResult(v, o, 2)
        
        xmlObj = a.toXML()
        b = PeakResult()
        b.fromXML(xmlObj)
        
        assert b.jobID == a.jobID
        assert b.result.getFilename() == a.result.getFilename()
        assert b.orient.getFilename() == a.orient.getFilename()
        
    # Module PeakJobMsg------------------------------------
    def PeakJobMsg_Test(self):
        from pytom.localization.peak_job_msg import PeakJobMsg
        
        from pytom.localization.peak_job import PeakJob
        from pytom.basic.structures import Mask, Reference, WedgeInfo
        from pytom.localization.structures import Volume
        from pytom.score.score import FLCFScore
        
        v = Volume('./testData/ribo.em')
        ref = Reference('./testData/ribo.em')
        m = Mask('./testData/ribo.em')
        w = WedgeInfo(30)
        s = FLCFScore()
    
        from pytom.angles.angleList import AngleList
        rot = AngleList([[1,1,1],[2,2,2],[3,3,3]])
        
        j = PeakJob(v, ref, m, w, rot, s)
        
        a = PeakJobMsg(str(0), str(1))
        a.setJob(j)
        xmlObj = a.toXML()
        b = PeakJobMsg()
        b.fromXML(xmlObj)
        
        assert b.getSender() == a.getSender()
        assert b.getRecipient() == a.getRecipient()
        
    def PeakResultMsg_Test(self):
        from pytom.localization.peak_job_msg import PeakResultMsg
        
        from pytom.localization.peak_job import PeakResult
        from pytom.localization.structures import Volume, Orientation
        
        v = Volume('./testData/ribo.em')
        o = Orientation('./testData/ribo.em')
        r = PeakResult(v, o)
        
        a = PeakResultMsg(str(1), str(0))
        a.setResult(r)
        xmlObj = a.toXML()
        b = PeakResultMsg()
        b.fromXML(xmlObj)
        
        assert b.getSender() == a.getSender()
        assert b.getRecipient() == a.getRecipient()
        
    
    def BackwardCompatibility_Test(self):
        from pytom.localization.peak_job import PeakJob
        job = PeakJob()
        job.fromXMLFile('./testData/xmlFiles/localizationJob.xml')

        
    
        
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
