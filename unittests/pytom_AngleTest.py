import unittest

class pytom_AngleTest(unittest.TestCase):
    
    def EulerAngleList_oldRotations_Test(self):
        from pytom.angles.angleList import EulerAngleList
        angle1 = EulerAngleList(0,0,0,1,1,1,10,10,10)
        angle1.focusRotation([5,5,5])
    
    def EulerAngleList_numberRotations_Test(self):
        from pytom.angles.angleList import EulerAngleList
        ang2 = EulerAngleList(0,0,0,1,1,1,2,2,2)
        assert len(ang2) == 7
    
    def LocalSampling_numberRotations_Test(self):
        from pytom.angles.localSampling import LocalSampling
        ang = LocalSampling(6,3,0,0,0)
        assert ang.numberRotations() == 521
    
    def LocalSampling_reset_Test(self):
        from pytom.angles.localSampling import LocalSampling
        ang = LocalSampling(2,3,0,0,0)
        rotations = ang.getAllRotations()
        ang.reset()
        rotations2 = ang.getAllRotations()
        assert rotations == rotations2
    
    def AngleEMList_focusRotation_Test(self):
        from pytom.angles.globalSampling import GlobalSampling
        from pytom_volume import vol
        import os
        
        v = vol(3,2,1)
        v.setAll(1)
        v.write('./angs.em')
        ang = GlobalSampling('angs.em',True)
        os.system('rm angs.em')
        ang2 = ang.focusRotation([1,1,1])
        assert ang2.__class__ == GlobalSampling
    
    def LocalSampling_focusRotation_Test(self):
        from pytom.angles.localSampling import LocalSampling

        ang = LocalSampling(1,3,4,50,23)
        ang2 = ang.focusRotation([1,1,1],0)
        assert ang2.__class__ == LocalSampling
    
    def LocalSampling_init_Test(self,a=1.0,b=3.0, c=123.370865619, d=134.740207266, e=95.2352127297):
        from pytom.angles.localSampling import LocalSampling
        ang = LocalSampling(a,b,c,d,e)
        ang2 = LocalSampling()
        ang2.fromStr('<Angles StartX="95.2352127297" StartZ1="123.370865619" StartZ2="134.740207266" NumberShells="1" AngleIncrement="3.0" ShellIncrement="3.0" Type="LocalSampling"/>')
        assert ang.numberRotations() == ang2.numberRotations()
    
    def AV3Sampling_Test(self):
        from pytom.angles.localSampling import LocalSampling as AV3Sampling
        a = AV3Sampling()
        a.fromStr('<Angles Type="Equidistant"><Parameters Increment="3.0" Shells="3.0" Phi_old="0" Psi_old="0" Theta_old="0" ShellsParameter="1" IncrementParameter="1"/></Angles>')
        assert str(a) == '<Angles Type="AV3Sampling" Increment="3.0" Shells="3.0" Phi_old="0.0" Psi_old="0.0" Theta_old="0.0" ShellsParameter="1" IncrementParameter="1"/>'
    
        from pytom.angles.angle import AngleObject
        a = AngleObject()
        a.fromStr('<Angles Type="Equidistant"><Parameters Increment="3.0" Shells="3.0" Phi_old="0" Psi_old="0" Theta_old="0" ShellsParameter="1" IncrementParameter="1"/></Angles>')

        assert str(a) == '<Angles Type="AV3Sampling" Increment="3.0" Shells="3.0" Phi_old="0.0" Psi_old="0.0" Theta_old="0.0" ShellsParameter="1" IncrementParameter="1"/>'
      
    def GlobalLocalCombined_Test(self):
        from pytom.angles.combined import GlobalLocalCombined
        from pytom.angles.localSampling import LocalSampling
        from pytom.angles.globalSampling import GlobalSampling
        
        ls = LocalSampling(1.0,3.0,123.370865619,134.740207266,95.2352127297)
        gs = GlobalSampling('angles_07_45123.em')  
        gl = GlobalLocalCombined(gs,ls)
    
    def conversion_Test(self,phi1=0,psi1=0,the1=0,phi2=90,psi2=90,the2=90):
        from pytom.angles.angleFnc import zxzToMat
        
        m = zxzToMat(phi1,psi1,the1,False)
        assert m.sum() == 3.0 and m[0][0] == 1.0 and m[1][1]== 1.0 and m[2][2]== 1.0
        m = zxzToMat(phi2,psi2,the2,False)
    
        from pytom.angles.angleFnc import matToZXZ
        a = matToZXZ(m)
        assert a[0] == phi2 and a[1] == psi2 and a[2] == the2

    def angleFromResolution_Test(self,resolution = 30, particleDiameter = 250,value = 6.8754935415698784):
        from pytom.basic.resolution import angleFromResolution
        assert value == angleFromResolution(resolution,particleDiameter)
        
    
    def pointRotateZXZTest(self,vector=[0,0,1],phi=90,psi=90,theta=90,result=[1,0,0]):
        from pytom.angles.angleFnc import pointRotateZXZ
        r = pointRotateZXZ(vector,phi,psi,theta)
        assert r[0] == result[0] and r[1] == result[1] and r[2] == result[2]
        
    def pointRotateZXZ_Tests(self):
        self.pointRotateZXZTest();
        self.pointRotateZXZTest([0,0,1],0,0,0,[0,0,1]);
        self.pointRotateZXZTest([0,0,1],0,0,90,[0,-1,0]);
        self.pointRotateZXZTest([0,0,1],0,0,-90,[0,1,0]);
        self.pointRotateZXZTest([0,0,1],10,0,0,[0,0,1]);
        self.pointRotateZXZTest([0,0,1],0,10,0,[0,0,1]);
        self.pointRotateZXZTest([0,0,1],10,10,0,[0,0,1]);
        self.pointRotateZXZTest([0,0,1],10,10,10,[0.030153689906001091, -0.17101006209850311, 0.98480772972106934]);
        
    def EquidistantList_Test(self):
        from pytom.angles.localSampling import LocalSampling
        
        # ang = LocalSampling(0, 2)
        # assert len(ang) == 180
        ang = LocalSampling(1, 2)
        ang2 = ang.focusRotation(ang.nextRotation())

        assert len(ang2) == 1332

           
    def zxzToMatToZXZ_T(self,z1,z2,x):
        from pytom.angles.angleFnc import zxzToMat,matToZXZ
        m = zxzToMat(z1,z2,x)
        r = matToZXZ(m)
        print(z1,z2,x,r.getZ1(),r.getZ2(),r.getX())
        assert abs(r.getZ1() - z1) < 0.00001 and abs(r.getZ2() - z2) < 0.00001 and abs(r.getX() - x) < 0.00001
    
    def zxzToMatToZXZ_Test(self):
        self.zxzToMatToZXZ_T(0, 0, 0)
        self.zxzToMatToZXZ_T(0, 290, 0)
        self.zxzToMatToZXZ_T(0, 0, 10)
        self.zxzToMatToZXZ_T(0, 10, 20)
        self.zxzToMatToZXZ_T(20, 0, 10)
        self.zxzToMatToZXZ_T(10, 20, 30)
        self.zxzToMatToZXZ_T(30, 10, 20)
        self.zxzToMatToZXZ_T(20, 30, 10)
        
    def axisAngleToMatToAxisAngle_T(self,axis,angle):
        from pytom.angles.angleFnc import zxzToAxisAngle,axisAngleToZXZ
        from pytom.tools.maths import epsilon
        r = axisAngleToZXZ(axis,angle)
        a = zxzToAxisAngle(r.getZ1(),r.getZ2(),r.getX())
        assert abs(angle-a[0]) < epsilon and abs(axis[0] - a[1][0]) < epsilon and abs(axis[1]- a[1][1]) < epsilon and abs(axis[2] - a[1][2]) < epsilon
        
        
    def axisAngleToMatToAxisAngle_Test(self):
        self.axisAngleToMatToAxisAngle_T([1,0,0],0)
        self.axisAngleToMatToAxisAngle_T([1,0,0],180)
        self.axisAngleToMatToAxisAngle_T([0,1,0],180)
        self.axisAngleToMatToAxisAngle_T([0,0,1],180)
        self.axisAngleToMatToAxisAngle_T([0,0,1],10)
        self.axisAngleToMatToAxisAngle_T([0,1,0],10)
        self.axisAngleToMatToAxisAngle_T([1,0,0],10)
        
    
    def indexAccess_T(self,angleObject=None,rotation=None,key=None):
        if not angleObject:
            from pytom.angles.angleList import AngleList 
            rl = [[0,0,0],[1,1,1],[2,2,2],[3,3,3]]
            angleObject = AngleList(rl)
        if not rotation:    
            rotation = [1,1,1]
        if not key:
            key = 1
        r = angleObject[key]
        assert r == rotation


    def indexAccess_Test(self):
        self.indexAccess_T()

        
    def AngleList_sliceAccess_Test(self):
        from pytom.angles.angleList import AngleList
        rl = [[0,0,0],[1,1,1],[2,2,2],[3,3,3]]
        angleObject = AngleList(rl)
        rotations = [[1,1,1],[2,2,2]]
        res = angleObject[1:3]
        assert res == rotations
    
    
    def LocalSamplingNumberAngles_Test(self):
        from pytom.angles.angle import AngleObject
        ang = AngleObject()
        ang.fromStr('<Angles Type="Equidistant"><Parameters Increment="3.0" \
             Shells="3.0" Phi_old="0" Psi_old="0" Theta_old="0" \
             ShellsParameter="1" IncrementParameter="1"/></Angles>')
        assert len(ang) == 148

    def matrixMult_Test(self):
        from pytom.basic.structures import Rotation
        from pytom.angles.angleFnc import differenceAngleOfTwoRotations

        rot1 = Rotation(z1=0.12480606462, z2=-0.00132220288025, x=0.106087546687)
        rot2 = Rotation(0,0,0)
        self.assertTrue(differenceAngleOfTwoRotations(rotation1=rot1, rotation2=rot2) < .2,
                        'difference between determined angles too large :(')

    
    def runTest(self):
        self.matrixMult_Test()
        self.LocalSamplingNumberAngles_Test()
        self.axisAngleToMatToAxisAngle_Test()
        
if __name__ == '__main__':
    unittest.main()
