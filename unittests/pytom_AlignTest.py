import unittest

class pytom_AlignTest(unittest.TestCase):
    
    def ExpectationMaximisationJobCheck_T(self,particleList,mask,reference,destination):
        from pytom.alignment.structures import ExpectationMaximisationJob
        emj = ExpectationMaximisationJob(particleList,destination,reference,0,0,mask,0,0,0,0,0,0.1,1,None)
        emj.check()


    def ExpectationMaximisationJobCheck_Test(self):
        from pytom.basic.structures import Mask,Reference,ParticleList,Particle
        from pytom_volume import vol
        import os
        
        pV = vol(10,10,10)
        pV2 = vol(10,10,10)
        
        pV.write('pV.em')
        pV2.write('pV2.em')
        
        pl = ParticleList('./')
        
        pl.append(Particle('pV.em'))
        pl.append(Particle('pV2.em'))
        
        mask = Mask('pV.em')
    
        reference = Reference('pV2.em')
    
        self.ExpectationMaximisationJobCheck_T(pl, mask, reference, './')
    
        pl.append(Particle('pV3.em'))
        try:
            self.ExpectationMaximisationJobCheck_T(pl, mask, reference, './')
            assert False
        except:
            pass
        
        os.system('rm pV.em')
        os.system('rm pV2.em')


    def Peak_Test(self):
        from pytom.alignment.structures import Peak
        from pytom.basic.structures import Rotation,Shift
        p1 = Peak(1,Rotation(0,0,0),Shift(0,0,0))
        p2 = Peak(0,Rotation(1,1,1),Shift(1,1,1))
        
        assert p2 < p1
        assert not p1 < p2 
    
    def bestAlignment_Test(self):
        from pytom.basic.structures import Rotation,Shift
        self.bestAlignment_T()
        self.bestAlignment_T(rotation = Rotation(12,45,3),shiftO=Shift(0,0,0))
        self.bestAlignment_T(rotation = Rotation(12,45,3),shiftO=Shift(-4,3,2))
    
    def bestAlignment_T(self,particle=None,reference=None,referenceWeighting=None,wedgeInfo=None,rotation=None,shiftO=None,scoreObject=None,mask=None,preprocessing=None):
        
        from pytom.alignment.alignmentFunctions import bestAlignment
        from pytom.angles.localSampling import LocalSampling
        from pytom.basic import filter
        from pytom.score.score import nxcfScore
        from pytom.basic.structures import WedgeInfo,Rotation,Shift
        from pytom.alignment.preprocessing import Preprocessing
        from pytom_volume import vol, initSphere, rotate, shift
        from pytom.basic.structures import Mask
        
        v = vol(32,32,32)
        v.setAll(0)
        vMod = vol(32,32,32)
        vRot = vol(32,32,32)
        v.setV(1,10,10,10)
        v.setV(1,20,20,20)
        v.setV(1,15,15,15)
        v.setV(1,7,21,7)
        
        pre = preprocessing or Preprocessing()
        scr = scoreObject or nxcfScore()
        
        wedgeInfo = wedgeInfo or WedgeInfo(30.0,[10.0,20.0,30.0],0.0)
        referenceWeighting = referenceWeighting or ''
        
        rotation = rotation or Rotation(10,20,30)
        shiftO = shiftO or Shift(1,-3,5)
        
        if not mask:
            m = vol(32,32,32)
            initSphere(m,m.sizeX()//2-3,1,0, 32//2, 32//2, 32//2)
            m.write('testMask.em')
            mask = Mask('testMask.em')
        
        rotate(v,vRot,rotation.getPhi(),rotation.getPsi(),rotation.getTheta())
        shift(vRot,vMod,shiftO.getX(),shiftO.getY(),shiftO.getZ())
        vMod = vMod * mask.getVolume()
        
        ang = LocalSampling(1,90,rotation.getZ1(),rotation.getZ2(),rotation.getX())
        
        peak = bestAlignment(vMod,v,referenceWeighting, wedgeInfo, ang, scr, mask, pre)
        
        import os
        os.system('rm testMask.em')
        print(peak.getShift(), shiftO)
        assert peak.getRotation() == rotation
        assert peak.getShift() == shiftO
        
   
    def Aligner_Test(self):
        self.Aligner_T()
    
    def Aligner_T(self,angle = [15,15,15],shift = [5,3,2],wedge = 30,snr = 0.1):
        from pytom.simulation import EMSimulation
        from pytom.alignment.ExMaxAlignment import ExMax
        from pytom.angles.angleList import EulerAngleList
        from pytom.score.score import FLCFScore as score
        from pytom.basic.structures import WedgeInfo
        import pytom_volume
        
        a = pytom_volume.read('/fs/home/hrabe/data/volumes/alignmentTest/ribOrig.em')
        b = EMSimulation.simpleSimulation(a,angle,shift,wedge,snr)
        
        ang = EulerAngleList(0,0,0,10,10,10,30,30,30)
        
        al = ExMax(a,b,score(),WedgeInfo(wedge),ang)
        result = al.run(2)
        
        if not (shift[0] == result.shift[0]*-1) and (shift[1] == result.shift[1]*-1) and (shift[2] == result.shift[2]*-1):
            print(shift,'!=',result.shift)
            
        if not angle == result.rotation:
            print(angle,'!=',result.rotation) 
        
        assert (shift[0] == result.shift[0]*-1) and (shift[1] == result.shift[1]*-1) and (shift[2] == result.shift[2]*-1)    
        assert angle == result.rotation

    def Expectation_Test(self):
        
        from pytom_volume import vol,rotate,initSphere,peak,read,abs,peak
        from pytom.basic.fourier import iftshift
        from pytom.basic.structures import WedgeInfo,Symmetry
        from pytom.alignment.structures import AlignmentList,MaximisationResult
        from pytom.alignment.manager import Manager
        from pytom.score.score import xcfScore
        from pytom.basic.filter import lowpassFilter
        import os
        
        particle = vol(64,64,64)
        mask = vol(64,64,64)
        particle.setAll(0)
        particle.setV(1,32,32,32)
        
        initSphere(mask,25,6,0)
        
        wedge = WedgeInfo(30,[0,0,0],0)
        
        p1 = wedge.apply(particle)
        p2 = wedge.apply(particle)
        p3 = wedge.apply(particle)
        p4 = wedge.apply(particle)
        p5 = wedge.apply(particle)
        
        r1 = [0,0,0]
        
        r2 = [0,0,90]
        r3 = [0,90,0]
        r4 = [90,0,0]
        r5 = [-100,60,45]
        
        rr = vol(64,64,64)
        
        rotate(p1,rr,r1[0],r1[1],r1[2])
        p1 = rr * mask
        
        rotate(p2,rr,r2[0],r2[1],r2[2])
        p2 = rr * mask
        
        rotate(p3,rr,r3[0],r3[1],r3[2])
        p3 = rr * mask
        
        rotate(p4,rr,r4[0],r4[1],r4[2])
        p4 = rr * mask
        
        
        rotate(p5,rr,r5[0],r5[1],r5[2])
        p5 = rr * mask
        
        p1.write('./p1.em')
        p2.write('./p2.em')
        p3.write('./p3.em')
        p4.write('./p4.em')
        p5.write('./p5.em')
        
        m1 = MaximisationResult('./p1.em','none',xcfScore(),[0,0,0],r1,wedge,1)
        m2 = MaximisationResult('./p2.em','none',xcfScore(),[0,0,0],r2,wedge,1)
        m3 = MaximisationResult('./p3.em','none',xcfScore(),[0,0,0],r3,wedge,1)
        m4 = MaximisationResult('./p4.em','none',xcfScore(),[0,0,0],r4,wedge,1)
        m5 = MaximisationResult('./p5.em','none',xcfScore(),[0,0,0],r5,wedge,1)
        
        al = AlignmentList()
        al.append(m1)
        al.append(m2)
        al.append(m3)
        al.append(m4)
        al.append(m5)
        
        
        ma = Manager()
        ma.alignmentList = al
        ma.destination = './'
        """
        this is the real test! everything before was preparing testing data
        """
             
        ma.sequentialExpectation('res.em')
        
        
        
        res = read('./res.em')
        
        lowParticle = lowpassFilter(particle,0.5)
        dif = res - lowParticle[0]
        abs(dif)
        p = peak(dif)
        pr = peak(res)
        

        
        assert res.getV(pr[0],pr[1],pr[2]) > 0.45 and pr == [32,32,32]
        """
        in case assert fails, check these remaining volumes for errors 
        """
        os.system('rm ./p1.em')
        os.system('rm ./p2.em')
        os.system('rm ./p3.em')
        os.system('rm ./p4.em')
        os.system('rm ./p5.em')
        os.system('rm ./res.em')
        os.system('rm ./res-INV.em')
        os.system('rm ./res-WedgeSum.em')
        os.system('rm ./res-WedgeSumHuman.em')
        os.system('rm ./res-PreWedge.em')
        
    def Symmetry_Test(self):
        self.Symmetry_T()
        self.Symmetry_T(40,40,40,-3,3)
        self.Symmetry_T(40,40,40,0,0,0,0,30)
        self.Symmetry_T(40,40,40,0,0,0,20,0)
        self.Symmetry_T(40,40,40,0,0,0,20,30)
        self.Symmetry_T(40,40,40,1,5,0,20,30)
        self.Symmetry_T(40,40,40,-1,5,0,78,-35)
        self.Symmetry_T(40,40,40,1,-5,0,-78,35)
        
    def Symmetry_T(self,x=40,y=40,z=40,sx=2,sy=5,sz=0,psi=0,theta=0):
        from pytom_volume import read,vol,rotate,shift
        from pytom.angles.angleFnc import pointRotateZXZ
        from pytom.basic.structures import WedgeInfo,Symmetry
        from pytom.alignment.structures import MaximisationResult,ExpectationJob
        from pytom.alignment.alignmentFunctions import cloneMaximisationResult
        from pytom.alignment.aligner import ExMax
        import os
        
        v = vol(60,60,60)
        v.setAll(0)
        
        v.setV(1,x,y,z)
        
        rot = vol(60,60,60)
        rot.setAll(0)
        rotate(v,rot,0,psi,theta)
        
        shift(rot,v,sx,sy,sz)

        v.write('./v.em')
        
        mr = MaximisationResult('./v.em',shift=[sx,sy,sz],rotation=[0,psi,theta],wedge=WedgeInfo(0,[0,0,0],0))
        
        mrs = cloneMaximisationResult(mr,Symmetry(3))
        
        ej = ExpectationJob(None,'./avg.em')
        
        ej.appendMaximisationResult(mr)
        
        for m in mrs:
            ej.appendMaximisationResult(m)
            #print m
        
        exMax = ExMax()
        
        exMax.fromJob(ej)
        exMax.run()
        
        rv = read('./avg.em')
        
        centerX = rv.sizeX()/2 
        centerY = rv.sizeY()/2 
        centerZ = rv.sizeZ()/2 
        
        v = [0 for _ in xrange(3)]
        v[0] = x - centerX
        v[1] = y - centerY
        v[2] = z - centerZ
        
        v = pointRotateZXZ(v,0,psi,theta)
        
        v1 = pointRotateZXZ(v,-psi,0,-theta)
        assert round(centerX + v1[0]) == x and round(centerY + v1[1]) == y and round(centerZ + v1[2]) == z
        assert 0 < (rv.getV(int(round(centerX + v1[0])),int(round(centerY + v1[1])),int(round(centerZ + v1[2])))) 
        
        #-----------------------------------------
        v2 = pointRotateZXZ(v,-psi,-120,-theta)
        assert 0 < (rv.getV(int(round(centerX + v2[0])),int(round(centerY + v2[1])),int(round(centerZ + v2[2]))))
        
        #-----------------------------------------
        v3 = pointRotateZXZ(v,-psi,-240,-theta)
        assert 0 < (rv.getV(int(round(centerX + v3[0])),int(round(centerY + v3[1])),int(round(centerZ + v3[2]))))
        
        os.system('rm ./v.em')
        os.system('rm ./avg.em')
        os.system('rm ./avg-INV.em')
        os.system('rm ./avg-WedgeSum.em')
        os.system('rm ./avg-WedgeSumHuman.em')
        os.system('rm ./avg-PreWedge.em')
        
    
    def BackwardCompatibility_Test(self):
        
        from pytom.alignment.structures import ExpectationMaximisationJob
        job = ExpectationMaximisationJob(0,0,0,0,0,0,0,0,0,0,0,0)
        job.fromXMLFile('./testData/xmlFiles/alignmentJob.xml')
        
        job2 = ExpectationMaximisationJob(0,0,0,0,0,0,0,0,0,0,0,0)
        job2.fromXMLFile('./testData/xmlFiles/alignmentJob2.xml')

    def runTest(self):
        self.bestAlignment_Test()
        
if __name__ == '__main__':
    unittest.main()

