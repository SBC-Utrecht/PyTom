import unittest;

class pytom_AngleTest(unittest.TestCase):
    
    def test_EulerAngleList_oldRotations(self):
        from pytom.angles.angleList import EulerAngleList
        angle1 = EulerAngleList(0,0,0,1,1,1,10,10,10)
        angle1.focusRotation([5,5,5])
    
    def test_EulerAngleList_numberRotations(self):
        from pytom.angles.angleList import EulerAngleList
        ang2 = EulerAngleList(0,0,0,1,1,1,2,2,2)
        self.assertTrue( len(ang2) == 7)
    
    def test_LocalSampling_numberRotations(self):
        from pytom.angles.localSampling import LocalSampling
        ang = LocalSampling(6,3,0,0,0)
        self.assertTrue( ang.numberRotations() == 521)
    
    def test_LocalSampling_reset(self):
        from pytom.angles.localSampling import LocalSampling
        ang = LocalSampling(2,3,0,0,0)
        rotations = ang.getAllRotations()
        ang.reset()
        rotations2 = ang.getAllRotations()
        self.assertTrue( rotations == rotations2)
    
    def test_AngleEMList_focusRotation(self):
        from pytom.angles.globalSampling import GlobalSampling
        from pytom_volume import vol
        import os
        
        v = vol(3,2,1)
        v.setAll(1)
        v.write('./angs.em')
        ang = GlobalSampling('angs.em',True)
        os.system('rm angs.em')
        ang2 = ang.focusRotation([1,1,1])
        self.assertTrue( ang2.__class__ == GlobalSampling)
    
    def test_LocalSampling_focusRotation(self):
        from pytom.angles.localSampling import LocalSampling

        ang = LocalSampling(1,3,4,50,23)
        ang2 = ang.focusRotation([1,1,1],0)
        self.assertTrue( ang2.__class__ == LocalSampling)
    
    def test_LocalSampling_init(self,a=1.0,b=3.0, c=123.370865619, d=134.740207266, e=95.2352127297):
        from pytom.angles.localSampling import LocalSampling
        ang = LocalSampling(a,b,c,d,e)
        ang2 = LocalSampling()
        ang2.fromStr('<Angles StartX="95.2352127297" StartZ1="123.370865619" StartZ2="134.740207266" NumberShells="1" AngleIncrement="3.0" ShellIncrement="3.0" Type="LocalSampling"/>')
        self.assertTrue( ang.numberRotations() == ang2.numberRotations())
    
    def test_AV3Sampling(self):
        from pytom.angles.localSampling import LocalSampling as AV3Sampling
        a = AV3Sampling()
        a.fromStr('<Angles Type="Equidistant"><Parameters Increment="3.0" Shells="3.0" Phi_old="0" Psi_old="0" Theta_old="0" ShellsParameter="1" IncrementParameter="1"/></Angles>')
        self.assertTrue( str(a) == '<Angles Type="AV3Sampling" Increment="3.0" Shells="3.0" Phi_old="0.0" Psi_old="0.0" Theta_old="0.0" ShellsParameter="1" IncrementParameter="1"/>')
    
        from pytom.angles.angle import AngleObject
        a = AngleObject()
        a.fromStr('<Angles Type="Equidistant"><Parameters Increment="3.0" Shells="3.0" Phi_old="0" Psi_old="0" Theta_old="0" ShellsParameter="1" IncrementParameter="1"/></Angles>')
        #self.assertTrue( str(a) == '<Angles Type="AV3Sampling" Increment="3.0" Shells="3.0" Phi_old="0.0" Psi_old="0.0" Theta_old="0.0" ShellsParameter="1" IncrementParameter="1"/>')
      
    def test_GlobalLocalCombined(self):
        from pytom.angles.combined import GlobalLocalCombined
        from pytom.angles.localSampling import LocalSampling
        from pytom.angles.globalSampling import GlobalSampling
        
        ls = LocalSampling(1.0,3.0,123.370865619,134.740207266,95.2352127297)
        gs = GlobalSampling('angles_07_45123.em')  
        gl = GlobalLocalCombined(gs,ls)
    
    def test_conversion(self,phi1=0,psi1=0,the1=0,phi2=90,psi2=90,the2=90):
        from pytom.angles.angleFnc import zxzToMat
        
        m = zxzToMat(phi1,psi1,the1,False)
        self.assertTrue( m.trace() == 3.0, msg='trace of matrix is wrong')
        self.assertTrue( m.getRow(0)[0] == 1.0 and m.getRow(1)[1]== 1.0 and m.getRow(2)[2]== 1.0)
        m = zxzToMat(phi2,psi2,the2,False)
    
        from pytom.angles.angleFnc import matToZXZ
        a = matToZXZ(m)
        self.assertTrue( a[0] == phi2 and a[1] == psi2 and a[2] == the2)

    def test_angleFromResolution(self,resolution = 30, particleDiameter = 250,value = 6.8754935415698784):
        from pytom.basic.resolution import angleFromResolution
        self.assertTrue( value == angleFromResolution(resolution,particleDiameter))
        
    
    def pointRotateZXZTest(self,vector=[0,0,1],phi=90,psi=90,theta=90,result=[1,0,0]):
        from pytom.angles.angleFnc import pointRotateZXZ
        r = pointRotateZXZ(vector,phi,psi,theta)
        self.assertTrue( r[0] == result[0] and r[1] == result[1] and r[2] == result[2])
        
    def test_pointRotateZXZ(self):
        self.pointRotateZXZTest();
        self.pointRotateZXZTest([0,0,1],0,0,0,[0,0,1]);
        self.pointRotateZXZTest([0,0,1],0,0,90,[0,-1,0]);
        self.pointRotateZXZTest([0,0,1],0,0,-90,[0,1,0]);
        self.pointRotateZXZTest([0,0,1],10,0,0,[0,0,1]);
        self.pointRotateZXZTest([0,0,1],0,10,0,[0,0,1]);
        self.pointRotateZXZTest([0,0,1],10,10,0,[0,0,1]);
        self.pointRotateZXZTest([0,0,1],10,10,10,[0.030153689906001091, -0.17101006209850311, 0.98480772972106934]);
        
    def test_LocalSamling(self):
        """
        test local sampling used for GLocal etc
        """
        from pytom.angles.localSampling import LocalSampling
        
        ang = LocalSampling(1, 2)
        ang2 = ang.nextRotation()
        self.assertAlmostEqual(first=358., second=ang2[1], places=3, msg='next angle is crap')
        angs=LocalSampling(shells=3, increment=3, z1Start=30.0, z2Start=40.0, xStart=70.0)
        #ang2 = ang.nextRotation()
        prevang = None
        mindist = 1000.
        maxdist = -1.
        for ii in range(0,len(angs)-1):
            (angdist, prevang) = self.compare_angles( angs, prevang)
            if angdist > maxdist:
                maxdist = angdist
            if angdist < mindist:
                mindist = angdist
        self.assertTrue( mindist < 0.1, 'LocalSampling does not get within 0.1 of starting angle')
        self.assertTrue( 3*3. < maxdist < 11. , 'LocalSampling max not in interval 9 ... 11.')

    def compare_angles(self, ang, prevang):
        from pytom.tools.maths import rotation_distance
        ang2 = ang.nextRotation()
        angdist = ang.distanceFunction(rotation=ang2)
        if prevang:
            neighbor_dist =  rotation_distance(ang1=prevang, ang2=ang2)
            #print(neighbor_dist)
        return angdist, ang2

    def zxzToMatToZXZ_T(self,z1,z2,x):
        from pytom.angles.angleFnc import zxzToMat,matToZXZ
        m = zxzToMat(z1=z1,z2=z2,x=x)
        r = matToZXZ(m)
        self.assertAlmostEqual(first=r.getZ1(), second=z1, places=3, msg='numerical issue in z1')
        self.assertAlmostEqual(first=r.getZ2(), second=z2, places=3, msg='numerical issue in z2')
        self.assertAlmostEqual(first=r.getX(), second=x, places=3, msg='numerical issue in x')
    
    def test_zxzToMatToZXZ(self):
        self.zxzToMatToZXZ_T(0, 0, 0)
        self.zxzToMatToZXZ_T(10, 0, 0.5)
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
        #self.assertTrue( abs(angle-a[0]) < epsilon and abs(axis[0] - a[1][0]) < epsilon and abs(axis[1]- a[1][1]) < epsilon and abs(axis[2] - a[1][2]) < epsilon)
        self.assertAlmostEqual( first=angle, second=a[0], places=2, msg='matrix conversion 1'  )
        self.assertAlmostEqual( first=axis[0], second=a[1][0], places=2, msg='xxx')
        self.assertAlmostEqual( first=axis[1], second=a[1][1], places=2, msg='xxx') 
        self.assertAlmostEqual( first=axis[2], second=a[1][2], places=2, msg='xxx') 
        
    def test_axisAngleToMatToAxisAngle(self):
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
        self.assertTrue( r == rotation)

    def test_indexAccess(self):
        self.indexAccess_T()
        
    def test_AngleList_sliceAccess(self):
        from pytom.angles.angleList import AngleList
        rl = [[0,0,0],[1,1,1],[2,2,2],[3,3,3]]
        angleObject = AngleList(rl)
        rotations = [[1,1,1],[2,2,2]]
        res = angleObject[1:3]
        self.assertTrue( res == rotations)
    
    def test_LocalSamplingNumberAngles(self):
        """
        test I/O of localSampling
        """
        from pytom.angles.localSampling import LocalSampling
        ang = LocalSampling()
        #print(ang)
        ang.fromStr('<Angles Type="Equidistant"><Parameters Increment="3.0" \
             Shells="3.0" Phi_old="0" Psi_old="0" Theta_old="0" \
             ShellsParameter="1" IncrementParameter="1"/></Angles>')
        self.assertTrue( len(ang) == 148, 
            'Number of angles for LocalSampling different from expevctation')

    def test_matrixMult(self):
        from pytom.basic.structures import Rotation
        from pytom.angles.angleFnc import differenceAngleOfTwoRotations

        rot1 = Rotation(z1=0.12480606462, z2=-0.00132220288025, x=0.106087546687)
        rot2 = Rotation(0,0,0)
        self.assertTrue(differenceAngleOfTwoRotations(rotation1=rot1, rotation2=rot2) < .2,
                        'difference between determined angles too large :(')

        
if __name__ == '__main__':
    unittest.main()
