import unittest
from pytom_volume import read, vol
from pytom_volume import transformSpline
from pytom.basic.structures import Rotation, Shift
from pytom.basic.correlation import nxcc
from pytom.alignment.localOptimization import Alignment

class pytom_LocalOptimizationTest(unittest.TestCase):
    
    def localOpti_Test(self):
        """
        test that local optimization of CC works
        """
        from pytom.angles.angleFnc import differenceAngleOfTwoRotations
        tmp  = read("testData/emd_1480.map.em_bin_4.em")
        vol2 = read("testData/emd_1480.map.em_bin_4.em")
        
        mask = vol(tmp.sizeX(),tmp.sizeX(),tmp.sizeX())
        mask.setAll(1.)
        trans = Shift(1.2, -3.4, 2.2)
        rot   = Rotation(10.,20.,30.)
        iniTrans = Shift(0., 0., 0.)
        iniRot   = Rotation(5.,25.,20.)
        
        vol1 = vol(tmp.sizeX(),tmp.sizeY(),tmp.sizeZ())
        transformSpline( tmp, vol1, rot[0], rot[1], rot[2],
                    int(tmp.sizeX()//2),int(tmp.sizeY()//2),int(tmp.sizeY()//2),
                    0,0,0,
                    trans[0],trans[1],trans[2])
        # for comparison scores
        alignment = Alignment(vol1=vol1, vol2=vol2, score=nxcc, mask=mask,
                              iniRot=rot, iniTrans=trans, opti='fmin_powell', verbose=0)
        finscore, optiRot, optiTrans = alignment.localOpti( iniRot=iniRot, iniTrans=iniTrans)
        cc_opti = alignment.cc(rot=optiRot, trans=optiTrans)

        cc_true = alignment.cc(rot=rot, trans=trans)
        self.assertTrue( cc_opti-cc_true < 0.001, 'Local Opti test broken :(')
        self.assertTrue(differenceAngleOfTwoRotations(rotation1=rot, rotation2=optiRot) < 0.1,
                        'difference between determined angles too large :(')
        #print optiRot
        #print differenceAngleOfTwoRotations(rotation1=rot, rotation2=optiRot)
        # test consistency of determined rotation and shift
        vol3 = vol(tmp.sizeX(),tmp.sizeY(),tmp.sizeZ())
        transformSpline( tmp, vol3, optiRot[0], optiRot[1], optiRot[2],
                    int(tmp.sizeX()/2),int(tmp.sizeY()/2),int(tmp.sizeY()/2),
                    0,0,0, optiTrans[0], optiTrans[1], optiTrans[2])
        alignment2 = Alignment( vol1=vol1, vol2=vol3, score=nxcc, mask=mask,
                iniRot=rot, iniTrans=trans, opti='fmin_powell', verbose=0)
        finscore2, optiRot2, optiTrans2 = alignment2.localOpti( iniRot=Rotation(0,0,0),
               iniTrans=Shift(0,0,0))
        self.assertTrue(differenceAngleOfTwoRotations(rotation1=Rotation(0,0,0), rotation2=optiRot2) < 0.2,
                        'difference between determined angles too large :(')

    def runTest(self):
        self.localOpti_Test()
        
if __name__ == '__main__':
    unittest.main()
