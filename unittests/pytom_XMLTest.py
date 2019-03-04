# -*- coding: iso-8859-1 -*-
import unittest

def ScoreFromXML(score):
    
    from pytom.score.score import fromXML
    
    score2 = fromXML(score.toXML())
    
    assert score2.__str__() == score.__str__()

class pytom_XMLTest(unittest.TestCase):
    
    def ExpectationMaximisationJobXML_Test(self):
        
        from pytom.angles.angleList import EulerAngleList
        from pytom.score.score import FLCFScore as score
        from pytom.alignment.structures import ExpectationMaximisationJob as EM
        from pytom.alignment.preprocessing import Preprocessing
        from pytom.basic.structures import Reference,PointSymmetry,SampleInformation
        
        lowBand = 0
        highBand = 0.5
        smooth = 4
        sourcePath = 'sourcePath'
        destinationPath = 'destinationPath'
        reference = Reference('reference')
        sampleInfo = SampleInformation(1,1)
        
        sc = score()
        
        ang = EulerAngleList(0,359,359,10,1,1,360,0,0)
        
        pre = Preprocessing(lowBand,highBand,smooth)
        
        destinationPath = destinationPath + str(lowBand)+'-'+str(highBand)+'-sm'+str(smooth)+'/'
        
        mask='/fs/home/hrabe/projects/alignmentScore/simulations/mask/mask60.em'
        #mask=''
        emjob = EM(sourcePath,destinationPath,reference,sc,ang,mask,PointSymmetry(),0,10,pre,2.0,1,sampleInfo)
        
        jstr = str(emjob)
        
        
        emjob2 = EM('zzz','sasdf',reference,sc,ang,mask,PointSymmetry(),0,10,pre,0.2,1,sampleInfo)
        
        emjob2.fromStr(jstr)
        
        #print(emjob,emjob2)
        
        assert emjob == emjob2
        
        
    def MaximisationJobXML_Test(self):

        from pytom.alignment.structures import MaximisationJob
        from pytom.alignment.preprocessing import Preprocessing
        from pytom.basic.structures import Mask,Reference
        import pytom.angles.angleList
        import pytom.score.score
        
        p = Preprocessing(0.1,0.2,5,[1.2,2.3,4.5])
        
        ang = pytom.angles.angleList.EulerAngleList(0,0,0,10,10,10,30,30,30)
        sc = pytom.score.score.nxcfScore()
        j = MaximisationJob('a',Reference('b'),sc,ang,Mask('f'),1,p)
        
        jj = MaximisationJob()
        s = str(j)
        jj.fromStr(s)
        
        assert j == jj    
        
        
    def MaximisationResultXML_Test(self):
        
        import pytom.parallel.alignmentMessages
        import pytom.alignment.structures
        import pytom.score.score
        from pytom.basic.structures import Reference
        from pytom.angles.localSampling import LocalSampling
        
        sc = pytom.score.score.nxcfScore()
        sc.setValue(0.324)
        
        ang = LocalSampling(6.0,3.0,1.0,0.0,1.0)
        
        res = pytom.alignment.structures.MaximisationResult('a',Reference('b'),sc,[0.0,0.0,0.0],[0.0,0.0,0.0],ang)
        rm = pytom.parallel.alignmentMessages.MaximisationResultMsg('','')
        rm.setResult(res)
        
        xmlStr = str(rm)
        
        rm2=pytom.parallel.alignmentMessages.MaximisationResultMsg('','')
        rm2.fromStr(xmlStr)

        assert str(rm) == str(rm2)
        
    def ExpectationJobXML_Test(self):
        
        from pytom.alignment.structures import MaximisationResult
        from pytom.alignment.structures import ExpectationJob
        from pytom.score.score import xcfScore
        from pytom.basic.structures import Reference
        
        sc = xcfScore()
        sc.setValue(0.23)
        sc2 = xcfScore()
        sc2.setValue(0.54)
        
        m = MaximisationResult('a',Reference('b'),sc,[1.0,1.0,1.0],[10.0,20.0,30.0])
        n = MaximisationResult('c',Reference('b'),sc2,[0.0,1.0,2.0],[30.0,40.0,50.0])
        
        e = ExpectationJob(particleList = None , newAverageName = '')
        
        e.appendMaximisationResult(m)
        e.appendMaximisationResult(n)
        
        eStr = e.__str__()
        
        e2 = ExpectationJob(particleList = None , newAverageName = '')
        e2.fromStr(eStr)
        
        e2Str = e2.__str__()
        
        assert eStr == e2Str
        
    def ExpectationJobMsgXML_Test(self):
        
        from pytom.alignment.structures import MaximisationResult
        from pytom.alignment.structures import ExpectationJob
        from pytom.basic.structures import Reference
        from pytom.parallel.alignmentMessages import ExpectationJobMsg
        from pytom.score.score import xcfScore
        
        sc = xcfScore()
        sc.setValue(0.23)
        sc2 = xcfScore()
        sc2.setValue(0.54)
        
        m = MaximisationResult('a',Reference('b'),sc,[1.0,1.0,1.0],[10.0,20.0,30.0])
        n = MaximisationResult('c',Reference('b'),sc2,[0.0,1.0,2.0],[30.0,40.0,50.0])
        
        e = ExpectationJob(None,'')
        
        e.appendMaximisationResult(m)
        e.appendMaximisationResult(n)
        
        msg = ExpectationJobMsg()
        msg.setJob(e)
        
        msgStr = msg.__str__()
        msg2 = ExpectationJobMsg()
        
        msg2.fromStr(msgStr)
        msg2Str=msg2.__str__()
        
        
        assert msgStr == msg2Str
    
    def ExpectationResultXML_Test(self):
        
        from pytom.alignment.structures import ExpectationResult
        
        ex = ExpectationResult('a')
        ex2 = ExpectationResult('b')
        
        s = ex.__str__()
        
        ex2.fromStr(s)
        
        assert ex.__str__() == ex2.__str__()

    
    def ExpectationResultMsgXML_Test(self):
        
        from pytom.alignment.structures import ExpectationResult
        from pytom.parallel.alignmentMessages import ExpectationResultMsg
        
        ex = ExpectationResult('a')
        ex2 = ExpectationResult('b')
        em = ExpectationResultMsg()
        em.setResult(ex)
        
        em2 = ExpectationResultMsg()
        em2.setResult(ex2)
        
        s = em.__str__()
        
        em2.fromStr(s)
        
        assert em.__str__() == em2.__str__()
        
        
    def EulerAngleListXML_Test(self):
        
        from pytom.angles.angleList import EulerAngleList
        
        ang = EulerAngleList(0,0,0,10,10,10,350,350,170)
        
        sang = ang.__str__()
        
        ang2 = EulerAngleList(0,0,0,10,10,10,10,10,10)
        
        ang2.fromStr(sang)
        
        assert sang == ang2.__str__()
        
        
    def AngleEquidistantXML_Test(self):
        
        from pytom.angles.localSampling import LocalSampling
        
        ang = LocalSampling(10,3,0.0,0.0,0.0)
        
        sang = ang.__str__()
        ang2 = LocalSampling()    
        ang2.fromStr(sang)
    
        assert sang == ang2.__str__()

        
    def AngleListFromEMXML_Test(self):
        
        from pytom.angles.globalSampling import GlobalSampling
        ang = GlobalSampling('angles_90_26.em',[[1.0,2.0,3.0],[4.0,5.0,6.0]])
        ang2 = GlobalSampling('angles_90_2.em')
        
        sang = ang.__str__()        
        
        ang2.fromStr(sang)
        assert ang2.__str__() == sang

        
    def PreprocessingXML_Test(self):

        
        from pytom.alignment.preprocessing import Preprocessing
        
        p = Preprocessing(0.1,0.2,4,[1.2,2.3,4.5],'File')
        
        pp = Preprocessing()
        
        x = p.toXML()
        
        pp.fromXML(x)
        
        assert p == pp
            

        
    def WedgeInfoXML_Test(self):
        
        from pytom.basic.structures import WedgeInfo
        
        wi = WedgeInfo(30.0,[10.1,20.2,30.3],0.0,'X',5.0)
        wi2 = WedgeInfo(0)
        
        wistr = str(wi)
        
        wi2.fromStr(wistr)
        assert wi == wi2
        
        wedge = wi.returnWedgeFilter(30,30,30)
        
        from pytom_freqweight import weight
        
        assert wedge.__class__ == weight
        
        wedge = wi.returnWedgeVolume(30,30,30,False)
        
        from pytom_volume import vol
        
        assert wedge.__class__ == vol
        assert wedge.sizeX() == 30 and wedge.sizeY() == 30 and wedge.sizeZ() == 16
        
        wedge = wi.returnWedgeVolume(30,30,30,True)
        assert wedge.sizeX() == 30 and wedge.sizeY() == 30 and wedge.sizeZ() == 30
    
    def AsymmetricWedgeInfoXML_Test(self):
        
        from pytom.basic.structures import WedgeInfo
        
        wi = WedgeInfo([20.0,30.0],[10.1,20.2,30.3],0.0,'X')
        wi2 = WedgeInfo(0)
        
        wistr = wi.__str__()
        
        wi2.fromStr(wistr)
        assert wi == wi2
        
        wedge = wi.returnWedgeFilter(30,30,30)
        
        from pytom_freqweight import weight
        
        assert wedge.__class__ == weight
        
        wedge = wi.returnWedgeVolume(30,30,30,False)
        
        from pytom_volume import vol
        
        assert wedge.__class__ == vol
        assert wedge.sizeX() == 30 and wedge.sizeY() == 30 and wedge.sizeZ() == 16
        
        wedge = wi.returnWedgeVolume(30,30,30,True)
        assert wedge.sizeX() == 30 and wedge.sizeY() == 30 and wedge.sizeZ() == 30
    
        
    def ReferenceXML_Test(self):

        from pytom.basic.structures import Reference
        
        r1 = Reference('a')
        r2 = Reference('b')
        
        rstr = str(r1)
        
        r2.fromStr(rstr)
        
        assert rstr == str(r2)
        
        
    def StatusMessageXML_Test(self):

        from pytom.parallel.messages import StatusMessage
        
        s = StatusMessage('a','b')
        s.setStatus('Message')
        ss = StatusMessage('c','d')
        s.setStatus('Wurst')
        
        sstr = s.__str__()
        
        ss.fromStr(sstr)
        
        assert ss.__str__() == sstr

        
    def AlignmentListXML_Test(self):
        
        from pytom.alignment.structures import MaximisationResult,AlignmentList
        from pytom.basic.structures import Reference
        from pytom.score.score import xcfScore
        
        sc = xcfScore()
        sc.setValue(0.23)
        sc2 = xcfScore()
        sc2.setValue(0.54)
        
        m = MaximisationResult('a',Reference('b'),sc,[1.0,1.0,1.0],[10.0,20.0,30.0])
        n = MaximisationResult('c',Reference('b'),sc2,[0.0,1.0,2.0],[30.0,40.0,50.0])
        
        al = AlignmentList() 
        al2 = AlignmentList()
        
        al.append(m)
        al.append(n)
        s1 = str(al)
        
        al2.fromStr(s1)
   
        
        assert al == al2
        
        
        
    def ParticleXML_Test(self,wa=30.0,pn='a',pr=[3.2,4.3,2.1],ps=[0.5,-2.1,1.5],className='A'):
        
        from pytom.basic.structures import Particle,PickPosition,Wedge
        wi = Wedge(wa)
        p = Particle(pn,pr,ps,wi,className,PickPosition())    
        
        p2 = Particle(pn+'b')
        
        pstr = str(p)
        
        p2.fromStr(pstr)

        assert pstr == str(p2)
    
    def MaskXML_Test(self,filename='./testData/ribo.em',isSphere=True,binning=1):
        from pytom.basic.structures import Mask
        
        m1 = Mask(filename,isSphere,binning)
        m2 = Mask('')
        
        m2.fromStr(str(m1))
        
        assert m1 == m2
    
    def ParticleListXML_Test(self,dirName='dir'):
        
        from pytom.basic.structures import ParticleList,Particle,Wedge
        wi = Wedge([30.0,30.0])
        p = Particle('a',[3.2,4.3,2.1],[0.5,-2.1,1.5],30)
        p2 = Particle('b',[3.2,4.3,2.1],[0.5,-2.1,1.4],20)
        
        fl = ParticleList()
        
        fl.append(p)
        fl.append(p2)
        fl.setWedgeAllParticles(wi)
        
        fl2 = ParticleList()
        
        flstr = str(fl)
        
        fl2.fromStr(flstr)

        assert flstr == str(fl2)
        
    def GrowingAverageJobXML_Test(self,particleList=[],angleObject=[],wedgeObject=[],maskFile=[],scoreObject=[],startParticleNumber=0):      
        
        if particleList == []:
            from pytom.basic.structures import ParticleList
            particleList = ParticleList()
        
        if angleObject == []:
            from pytom.angles.localSampling import LocalSampling
            angleObject = LocalSampling()
        
        if wedgeObject == []:
            from pytom.basic.structures import Wedge
            wedgeObject = Wedge()
        
        if maskFile == []:
            maskFile = 'a'
            
        if scoreObject == []:
            from pytom.score.score import xcfScore
            scoreObject = xcfScore()
            
        from pytom.alignment.structures import GrowingAverageJob
        
                            
        j1 = GrowingAverageJob(particleList,angleObject,maskFile,scoreObject,startParticleNumber)
        j2 = GrowingAverageJob(particleList,angleObject,maskFile,scoreObject,1)
        
        j1str = str(j1)
        
        j2.fromStr(j1str)
        
        assert j1 == j2
        
        
    def ScoresFromXML_Test(self):
        
        from pytom.score.score import xcfScore,nxcfScore,SOCScore,FLCFScore,RScore,FSCScore
        
        ScoreFromXML(xcfScore())
        ScoreFromXML(nxcfScore())
        ScoreFromXML(FLCFScore())
        ScoreFromXML(SOCScore())
        ScoreFromXML(RScore())
        ScoreFromXML(FSCScore())
        
    def SymmetryXML_Test(self):
        
        from pytom.basic.structures import Symmetry,PointSymmetry,HelicalSymmetry
        
        sym1 = PointSymmetry(2,30,40)
        sym1Str = str(sym1)
        
        symFactory = Symmetry()
        sym2 = symFactory.fromStr(sym1Str)
        
        assert sym2 == sym1
        
        sym3 = HelicalSymmetry(1,10,8)
        sym3Str = str(sym3)
        
        sym4 = symFactory.fromStr(sym3Str)
        assert sym3 == sym4
        
    def PeakPriorXML_Test(self):
        
        from pytom.score.score import nxcfScore,PeakPrior
        
        d = PeakPrior('',1.0,2.0)
        d2 = PeakPrior()
        
        d2.fromStr(str(d))
        
        assert d == d2    
    
    def EquidistantList_T(self,value = None):
        
        from pytom.angles.localSampling import LocalSampling
        
        ang  = LocalSampling()
        ang2 = LocalSampling()
        
        if value.__class__ == str:
            ang.fromStr(value)
        elif value:
            ang.fromXML(value)
        
        angXML = ang.toXML()
        ang2.fromXML(angXML)
        
        assert ang == ang2
        assert ang.numberRotations() == ang2.numberRotations()
        
    def EquidistantListXML_Test(self):
        #self.EquidistantList_T()
        self.EquidistantList_T('<Angles NumberShells="1" StartZ2="-40.4622274921" StartZ1="8.0" StartX="106.220088487" AngleIncrement="1" Type="LocalSampling"/>')
    
    def RestrictedInplaneEuidistantListXML_Test(self):
        from pytom.angles.localSampling import RestrictedInplaneLocalSampling


        ang = RestrictedInplaneLocalSampling(1,10,10,20,10,None,90)
        ang2 = RestrictedInplaneLocalSampling()
        ang2.fromStr(str(ang))
    
        ang = RestrictedInplaneLocalSampling(1,10,4.5,3.2,-84,12,34)
        ang2 = RestrictedInplaneLocalSampling()
        ang2.fromStr(str(ang))
        
    def CorrelationMatrixJobXML_Test(self):
        from pytom.cluster.correlationMatrixStructures import CorrelationMatrixJob
        from pytom.basic.structures import ParticleList,Mask
        
        job = CorrelationMatrixJob(ParticleList('/'),Mask(''),'c',True,3.0,0.1,0.2)
        
        job2 = CorrelationMatrixJob()
        job2.fromStr(str(job))
        
        #print job
        #print job2
        
        assert job == job2
        
    def ReferenceListXML_Test(self):
        
        from pytom.basic.structures import Reference,ReferenceList
        
        ref1 = Reference('a')
        ref2 = Reference('b')
        
        refList = ReferenceList()
        refList.append(ref1)
        refList.append(ref2)
        
        refList2 = ReferenceList()
        
        refListStr = str(refList)
        refList2.fromStr(refListStr)
        
        assert refList == refList2
        
        
    def CorrelationVectorJobXML_Test(self):
        
        from pytom.cluster.correlationMatrixStructures import CorrelationVectorJob
        from pytom.basic.structures import Particle,ParticleList,Mask
        
        job = CorrelationVectorJob(Particle(''),ParticleList('/'),Mask(''),1,True,2.0,0.1,0.2)
        job2 = CorrelationVectorJob()
        
        job2.fromStr(str(job))
        
        #print job
        #print job2
        
        assert job == job2
        
        
    def ShiftXML_T(self,x,y,z):
        from pytom.basic.structures import Shift 
        
        s1 = Shift(x,y,z)
        s2 = Shift(0,0,0)
        
        s2.fromStr(str(s1))
        
        assert s1 == s2
        
    def ShiftXML_Test(self):
        
        self.ShiftXML_T(1.0, 1.0, 1.0)
        self.ShiftXML_T(1, 1, 1)
        self.ShiftXML_T(0.1, 0.2, 0.3)
        
    def MultiRefAlignJobXML_Test(self): 
        
        from pytom.alignment.MultiRefStructures import MultiRefAlignJob
        from pytom.basic.structures import Reference,ReferenceList,ParticleList,Particle,WedgeInfo,Mask,Symmetry,SampleInformation
        from pytom.score.score import FLCFScore as score
        from pytom.alignment.preprocessing import Preprocessing
        from pytom.angles.fromFile import AngleListFromEM
        
        wi = WedgeInfo(30.0)
        
        p = Particle('a',[3.2,4.3,2.1],[0.5,-2.1,1.5],30)
        p2 = Particle('b',[3.2,4.3,2.1],[0.5,-2.1,1.4],20)
        
        pl = ParticleList('/',[])
        
        ref1 = Reference('a')
        ref2 = Reference('b')
        
        refList = ReferenceList()
        refList.append(ref1)
        refList.append(ref2)
        
        rotations= AngleListFromEM('angles_90_26.em',refinementParameters=[2,10])
        
        j1 = MultiRefAlignJob(pl,'',refList,score(),rotations,Mask(''),wi,Symmetry(1.0),1,Preprocessing(),0.1,1,SampleInformation(4.7,250.0),0.5,numberClassificationIterations=10)
        j2 = MultiRefAlignJob()
        
        j2.fromStr(str(j1))

        assert j1 == j2
        
        
#    def SimulatedAnnealingAlignJobXML_Test(self):
#        from pytom.alignment.SimulatedAnnealingStructures import SimulatedAnnealingJob
#        from pytom.basic.structures import Reference,ReferenceList,ParticleList,Particle,WedgeInfo,Mask,Symmetry,SampleInformation
#        from pytom.score.score import FLCFScore as score
#        from pytom.alignment.preprocessing import Preprocessing
#        from pytom.angles.fromFile import AngleListFromEM
#        
#        wi = WedgeInfo(30.0)
#        
#        p = Particle('a',[3.2,4.3,2.1],[0.5,-2.1,1.5],30)
#        p2 = Particle('b',[3.2,4.3,2.1],[0.5,-2.1,1.4],20)
#        
#        pl = ParticleList('/',[])
#        
#        ref1 = Reference('a')
#        ref2 = Reference('b')
#        
#        refList = ReferenceList()
#        refList.append(ref1)
#        refList.append(ref2)
#        
#        rotations= AngleListFromEM('angles_90_26.em',refinementParameters=[2,10])
#        
#        j1 = SimulatedAnnealingJob(pl,'',refList,score(),rotations,Mask(''),wi,Symmetry(1.0),1,Preprocessing(),0.1,1,SampleInformation(4.7,250.0),0.5,10,1.0)
#        j2 = SimulatedAnnealingJob()
#        
#        j2.fromStr(str(j1))
#        
#        assert j1 == j2
        
    def MCOEXMXJobXML_Test(self):
        from pytom.cluster.mcoEXMXStructures import MCOEXMXJob
        from pytom.basic.structures import ParticleList,WedgeInfo,Mask,SampleInformation
        from pytom.score.score import FLCFScore as score
        from pytom.alignment.preprocessing import Preprocessing
        
        job = MCOEXMXJob(ParticleList('/'),10,'a',Mask(''),score(),Preprocessing(),WedgeInfo(30.0),1,SampleInformation(4.7,250.0),5,0.01)
        job2 = MCOEXMXJob(ParticleList('/'),10,'a',Mask(''),score(),Preprocessing(),WedgeInfo(30.0),1,SampleInformation(4.7,250.0),1,0.01)
        
        assert not (job == job2)
        job2.fromStr(str(job))

        assert job == job2
        
        
    def AnnealingTemperature_Test(self):
        from pytom.cluster.mcoACStructures import AnnealingTemperature,SigmaTemperature,temperatureFromXML
        
        temp1 = AnnealingTemperature(10,10)
        temp2 = SigmaTemperature(3,1)
        
        temp3 = temperatureFromXML(temp1.toXML())
        
        assert temp1 == temp3
        
        temp4 = temperatureFromXML(temp2.toXML())
        assert temp2 == temp4
        
    def AnnealingJob_Test(self):
        from pytom.cluster.mcoACStructures import MCOACJob,AnnealingTemperature,MetropolisCriterion
        from pytom.basic.structures import ParticleList,WedgeInfo,Mask,SampleInformation
        from pytom.score.score import FLCFScore as score
        from pytom.alignment.preprocessing import Preprocessing
        job = MCOACJob(ParticleList('/'),'a',Mask(''),score(),Preprocessing(),WedgeInfo(30.0),1,SampleInformation(4.7,250.0),5,AnnealingTemperature(10,3),MetropolisCriterion(),0.01)
        job2 = MCOACJob(ParticleList('/'),'a',Mask(''),score(),Preprocessing(),WedgeInfo(30.0),1,SampleInformation(4.7,250.0),1,AnnealingTemperature(1,1),MetropolisCriterion(),0.1) 
        
        assert not (job == job2)
        
        job2.fromStr(str(job))

        assert job == job2
        
        
    def ClusterSwap_Test(self):
        from pytom.cluster.mcoEXMXStructures import ClusterSwap
        from pytom.basic.structures import Particle
        from pytom.score.score import FLCFScore as score
        
        
        s1 = ClusterSwap(Particle('a'),'1',score(),'2',score())
        s2 = ClusterSwap(Particle('a'),'1',score(),'1',score())
        
        s2.fromStr(str(s1))
        
        assert s1 == s2
        
    def SwapList_Test(self):
    
        from pytom.cluster.mcoEXMXStructures import ClusterSwap,SwapList
        from pytom.basic.structures import Particle
        from pytom.score.score import FLCFScore as score
        
        
        s1 = ClusterSwap(Particle('a'),'1',score(),'2',score())
        s2 = ClusterSwap(Particle('a'),'1',score(),'1',score())
        
        sl1 = SwapList([s1,s2])
        sl2 = SwapList()
        
        sl2.append(s1)
        
        sl2.fromStr(str(sl1))
        
        assert sl1 == sl2
        
    def ProjectionXML_Test(self):        
        from pytom.reconstruction.reconstructionStructures import Projection
        
        p1 = Projection('p1', 20)
        p2 = Projection('p2', 30)
        
        assert not (p1 == p2)
        p1.fromStr(str(p2))
        
        assert p1 == p2
               
        
    def ProjectionListXML_Test(self):
        from pytom.reconstruction.reconstructionStructures import Projection, ProjectionList
        
        p1 = Projection('p1', 20)
        p2 = Projection('p2', 30)
        p3 = Projection('p3', 40)
               
        pl1 = ProjectionList([p1, p2])
        pl2 = ProjectionList([p1, p3])
        
        assert not (pl1 == pl2)
        p3.fromStr(str(p2))
        
        assert pl1 == pl2
        
    def AnnealingCriterionXML_Test(self):
        from pytom.cluster.mcoACStructures import MetropolisCriterion
        
        m = MetropolisCriterion()
        m2 = MetropolisCriterion(True)
        
        m2.fromStr(str(m))
        
        assert m == m2
        
        
    def SingleTiltWedgeXML_Test(self):
        from pytom.basic.structures import Wedge,Particle
        
        w2 = Wedge()
        w = Wedge([30.0,30.0])
        w2.fromStr(str(w))
        assert w == w2
        
        p = Particle('a')
        p.setWedge(w2)
        
        p2 = Particle('b')
        
        p2.fromStr(str(p))
        
        assert p == p2
    
    def DoubleTiltWedgeXML_Test(self):
        from pytom.basic.structures import Wedge,Particle
        
        w2 = Wedge()
        w = Wedge([[30.0,30.0],[30.0,30.0]])
        
        w2.fromStr(str(w))
        
        assert w == w2
        
        p = Particle('a')
        p.setWedge(w2)
        
        p2 = Particle('b')
        
        p2.fromStr(str(p))
        
        assert p == p2


