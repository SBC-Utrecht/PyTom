import unittest

class pytom_ScoreTest(unittest.TestCase):
    
    
    def xcfScore_Test(self,v=0,wedge=30.0):
        
        from pytom.score.score import xcfScore as score
        from pytom_volume import vol,initSphere,peak
        from pytom.simulation.EMSimulation import simpleSimulation
        from pytom.basic.structures import WedgeInfo
        
        if v.__class__ == int:
            v = vol(32,32,32);
            initSphere(v,10,2,0,15,15,15)
            
        wi = WedgeInfo(wedge,[10.0,20.0,30.0],0.0)
        
        s = simpleSimulation(v,[10,10,10],[1,2,3],wi,0.1)
        
        sc = score()
        
        c = sc.scoringCoefficient(s,v)
        cf = sc.scoringFunction(s,v)
        
        p= peak(cf)
        print(c, cf.getV(p[0],p[1],p[2]) )
        assert c == (cf.getV(p[0],p[1],p[2]))
        

    def nxcfScore_Test(self,v=0,wedge=30.0):
        
        from pytom.score.score import nxcfScore as score
        from pytom_volume import vol,initSphere,peak
        from pytom.simulation.EMSimulation import simpleSimulation
        from pytom.basic.structures import WedgeInfo
        
        if v.__class__ == int:
            v = vol(32,32,32)
            initSphere(v,10,2,0,15,15,15)
        
            
        wi = WedgeInfo(wedge,[10.0,20.0,30.0],0.0)
        
        s = simpleSimulation(v,[10,10,10],[1,2,3],wi,0.1)
        
        sc = score()
        
        c = sc.scoringCoefficient(s,v)
        cf = sc.scoringFunction(s,v)
        
        p = peak(cf)
        print(c, cf.getV(p[0],p[1],p[2]))
        assert c == cf.getV(p[0],p[1],p[2])
        
    def flcfScore_Test(self,v=0,wedge=30.0):
        
        from pytom.score.score import FLCFScore as score
        from pytom_volume import vol,initSphere,peak
        from pytom.simulation.EMSimulation import simpleSimulation
        from pytom.basic.structures import WedgeInfo
        
        if v.__class__ == int:
            v = vol(32,32,32)
            initSphere(v,10,2,0,15,15,15)
            
        wi = WedgeInfo(wedge,[10.0,20.0,30.0],0.0)
        
        s = simpleSimulation(v,[10,10,10],[1,2,3],wi,0.1)
        
        sc = score()
        c  = sc.scoringCoefficient(s,v)
        cf = sc.scoringFunction(s,v)
        
        p= peak(cf)
        print(c, cf.getV(p[0], p[1], p[2]))
        assert c == cf.getV(p[0],p[1],p[2])
        
    def socScore_Test(self,v=0,wedge=30.0):
        
        from pytom.score.score import SOCScore as score
        from pytom_volume import vol,initSphere,peak
        from pytom.simulation.EMSimulation import simpleSimulation
        from pytom.basic.structures import WedgeInfo
        
        if v.__class__ == int:
            v = vol(32,32,32);
            initSphere(v,10,2,0,15,15,15)
            
        wi = WedgeInfo(wedge,[10.0,20.0,30.0],0.0)
                
        s = simpleSimulation(v,[10,10,10],[1,2,3],wi,0.1)
        
        sc = score()
        
        c = sc.scoringCoefficient(s,v)
        cf = sc.scoringFunction(s,v)
        
        p= peak(cf)
        print(c, cf.getV(p[0], p[1], p[2]))
        assert c == cf.getV(p[0],p[1],p[2])    
    
    def RScore_Test(self):
        """
        """
        