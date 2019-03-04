'''
Created on Aug 31, 2009

@author: hrabe
'''

import unittest

class pytom_StructuresTest(unittest.TestCase):
    
    def Wedge_getWedgeVolume_Test(self):
        
        from pytom.basic.structures import Wedge,Rotation
        from pytom_volume import vol
    
        wi = Wedge(30,0,Rotation(10,20,30))
        f1 = wi.returnWedgeFilter(60,60,60)
        f1.rotate(10,20,30)
        wf1 = f1.getWeightVolume(True)
        
        assert wf1.__class__ == vol


    def Wedge_apply_Test(self):
        
        from pytom_volume import vol,rotate,initSphere

        particle = vol(64,64,64)
        particle.setAll(0)
        particle.setV(1,33,33,33)
        
        from pytom.basic.structures import Wedge
        
        wedge1 = Wedge(30,0)
        
        p1 = wedge1.apply(particle)
        
        from pytom.basic.filter import wedgeFilter
        
        p2 = wedgeFilter(particle,30,32,True)
        
        assert p1.equalsTo(p2[0])


    def DoubleTiltWedge_Test(self):
        from pytom.basic.structures import DoubleTiltWedge,Rotation

        w1 = DoubleTiltWedge([[30,30],[30,30]]) 
        w1.returnWedgeVolume(100,100,100,humanUnderstandable = False)
        w1.returnWedgeVolume(100,100,100,humanUnderstandable = True)

        w2 = DoubleTiltWedge([[30,30],[30,30]],rotation12 = Rotation(10,0,20))
        w2.returnWedgeVolume(100,100,100,humanUnderstandable = False)
        w2.returnWedgeVolume(100,100,100,humanUnderstandable = True)


    def ParticleList_getParticleByFilename_Test(self):
        from pytom.basic.structures import Particle,ParticleList
        import os
        pl = ParticleList(' ',[])
        pl2 = ParticleList(' ',[])
        p1 = Particle('spikesBin2_NORMED/p1_1.em')
        p2 = Particle('p1_1.em')
        
        pl2.append(p1)
        pl2.append(p2)
        
        pl2.toXMLFile('test.xml')
        pl.fromXMLFile('test.xml')
        
        p= pl2.getParticleByFilename('p1_1.em')
        assert p.__str__() == p2.__str__()
        
        p = pl2.getParticleByFilename('spikesBin2_NORMED/p1_1.em')
        assert p.__str__() == p1.__str__()
        
        p = pl.getParticleByFilename('p1_1.em')
        assert p.__str__() == p2.__str__()
        
        p = pl.getParticleByFilename('spikesBin2_NORMED/p1_1.em')
        assert p.__str__() == p1.__str__()
        
        os.system('rm test.xml')
        
    def ParticleList_sliceTest(self):
        from pytom.basic.structures import Particle,ParticleList
        pl = ParticleList('/',[])
        pl2 = ParticleList('/',[])
        p1 = Particle('1')
        p2 = Particle('2')
        p3 = Particle('3')
        
        pl.append(p1)
        pl.append(p2)
        pl.append(p3)
        
        pl2.append(p1)
        
        pl3 = pl[0:1]
        assert pl2 == pl3
        
        pl2.append(p2)
        pl3 = pl[0:2]
        
        assert pl2 == pl3
        
        pl3 = pl[0:]
    
        assert pl == pl3
    
    def ParticleList_assignTest(self):
        from pytom.basic.structures import Particle,ParticleList
        
        pl = ParticleList('/',[])
        pl2 = ParticleList('/',[])
        p1 = Particle('1')
        p2 = Particle('2')
        p3 = Particle('3')
        
        pl.append(p1)
        pl.append(p2)
        pl.append(p3)
        
        pl2.append(p3)
        pl2[0] = p1
        
        assert pl2[0] == p1
        
        pl2.append(p2)
        pl2.append(p3)
        
        pl[1] = p3
        
        pl2[1:2] = pl[1:2] 
    
        assert pl2[1] == pl[1] and pl2[2] == pl[2]
        
        
    
    def ParticleList_addTest(self):
        from pytom.basic.structures import Particle,ParticleList
        pl1 = ParticleList('/',[])
        pl2 = ParticleList('/',[])
        pl3 = ParticleList('/',[])
        p1 = Particle('1')
        p2 = Particle('2')
        
        pl1.append(p1)
        pl2.append(p2)
        
        pl3.append(p1)
        pl3.append(p2)
        
        pl4 = pl1 + pl2
        
        assert pl3 == pl4
    
    def ParticleList_checkTest(self):
        from pytom.basic.structures import ParticleList,Particle
        from pytom_volume import vol
        import os
        
        pV = vol(10,10,10)
        pV2 = vol(10,10,10)
        
        pV.write('pV.em')
        pV2.write('pV2.em')
        
        pl = ParticleList('./')
        
        pl.append(Particle('pV.em'))
        pl.append(Particle('pV2.em'))
        
        pl.check()
        
        pl.append(Particle('pV3.em'))
        
        try:
            pl.check()
            assert False
        except:
            pass
        
        os.system('rm pV.em')
        os.system('rm pV2.em')
    
    
    def Mask_checkTest(self):
        from pytom.basic.structures import Mask
        from pytom_volume import vol
        import os
        
        pV = vol(10,10,10) 
        pV.write('m.em')
        
        m = Mask('m.em')
        
        m.check()
        
        try:
            m = Mask('m.em')
            m.check()
            assert False
        except:
            pass
        
        os.system('rm m.em')
    
    def Reference_checkTest(self):
        from pytom.basic.structures import Reference
        from pytom_volume import vol
        import os
        
        pV = vol(10,10,10) 
        pV.write('r.em')
        
        r = Reference('r.em')
        
        r.check()
        
        try:
            r = Reference('r.em')
            r.check()
            assert False
        except:
            pass
        
        os.system('rm r.em')
      
    def ReferenceList_accessTest(self):
        
        from pytom.basic.structures import Reference,ReferenceList
        
        ref1 = Reference('a')
        ref2 = Reference('b')
        
        refList = ReferenceList()
        refList.append(ref1)
        refList.append(ref2)
        
        ref3 = refList[0]
        
        assert ref1 == ref3
        
        ref4 = refList['b']
        
        assert ref2 == ref4
        
        
    def AlignmentList_sortByParticleListTest(self):
        from pytom.alignment.structures import AlignmentList,MaximisationResult
        from pytom.basic.structures import Particle,ParticleList
        from pytom.basic.structures import Reference
        from pytom.angles.localSampling import LocalSampling
        from pytom.score.score import xcfScore
        
        p1 = Particle('a')
        p2 = Particle('b')
        p3 = Particle('c')
        
        pl = ParticleList('/')
        pl.append(p1)
        pl.append(p2)
        pl.append(p3)
        
        ang = LocalSampling(6.0,3.0,1.0,0.0,1.0)
        
        r1 = MaximisationResult('a',Reference('b','x'),xcfScore(),[0.0,0.0,0.0],[0.0,0.0,0.0],ang)
        r2 = MaximisationResult('b',Reference('b','x'),xcfScore(),[0.0,0.0,0.0],[0.0,0.0,0.0],ang)
        r3 = MaximisationResult('c',Reference('b','x'),xcfScore(),[0.0,0.0,0.0],[0.0,0.0,0.0],ang)
        
        al = AlignmentList()
        al.append(r1)
        al.append(r3)
        al.append(r2)
    
        al.sortByParticleList(pl)
        pl2 = al.toParticleList()
        
        assert pl[0].getFilename() == pl2[0].getFilename()
        assert pl[1].getFilename() == pl2[1].getFilename()
        assert pl[2].getFilename() == pl2[2].getFilename()
         
        
        
    def Shift_OperatorTest(self):
        
        from pytom.basic.structures import Shift
        
        s1 = Shift(0,0,0)
        s2 = Shift(1,1,1)
        
        s3 = s1 + 1
        assert s2 == s3
        
        s3 = 1 + s1
        assert s2 == s3
        
        s3 = s2 - 1
        assert s1 == s3

        s4 = Shift(5,5,5)
        s3 = s2 * 5
        assert s3 == s4
        
        s3 = 5 * s2
        assert s3 == s4
        
        s5 = Shift(1,1,1)
        s6 = Shift(1,0,0)
        s7 = Shift(0,1,-1)
        
        s3 = s5 * s6
        assert s3 == s7
        
        
        
        
        