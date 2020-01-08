import unittest

class pytom_ScoreTest(unittest.TestCase):
    
    def setUp(self):
        """set up"""
        from pytom_volume import vol, initSphere
        from pytom.basic.structures import WedgeInfo
        from pytom.simulation.EMSimulation import simpleSimulation

        self.wedge = 30.
        self.shift = [1, 2, 3]
        self.rotation = [0, 0, 0]
        # create sphere
        self.v = vol(32,32,32)
        initSphere(self.v,10,2,0,15,15,15)
        self.wi = WedgeInfo(wedgeAngle=self.wedge, rotation=[10.0,20.0,30.0], 
              cutoffRadius=0.0)
        self.s = simpleSimulation( volume=self.v, rotation=self.rotation, 
              shiftV=self.shift, wedgeInfo=self.wi, SNR=0.1)

    def test_xcfScore(self):
        from pytom.score.score import xcfScore as score
        from pytom_volume import peak

        sc = score()
        # consistency of scoring coefficient and scoring function - difference 
        # due to sub-pixel accuracy for score
        c  = sc.scoringCoefficient( self.s, self.v)
        cf = sc.scoringFunction( self.s, self.v)
        p= peak(cf)
        self.assertAlmostEqual( first=c, second=cf.getV(16,16,16), places=2, 
            msg='Scoring coefficient and scoring function XCF inconsistent')
        self.assertLess( c, cf.getV(p[0],p[1],p[2]), 
            'Scoring coefficient and scoring function XCF inconsistent')

    def test_nxcfScore(self):
        """
        test nxcf score
        """
        from pytom.score.score import nxcfScore as score
        from pytom_volume import peak
        sc = score()
        # check auto-correlation coefficient
        c  = sc.scoringCoefficient( self.s, self.s)
        self.assertAlmostEqual( first=c,  second=1., places=4, 
             msg='XCFScore: Auto-correlation not == 1')
        # consistency of scoring coefficient and scoring function - difference 
        # due to sub-pixel accuracy for score
        c  = sc.scoringCoefficient( self.s, self.v)
        cf = sc.scoringFunction( self.s, self.v)
        p= peak(cf)
        self.assertAlmostEqual( first=c, second=cf.getV(16,16,16), places=2, 
            msg='Scoring coefficient and scoring function NXCF inconsistent')
        self.assertLess( c, cf.getV(p[0],p[1],p[2]), 
            'Scoring coefficient and scoring function NXCF inconsistent')
        
        
    def test_flcfScore(self):
        """
        test FLCF score
        """
        from pytom.score.score import FLCFScore as score
        from pytom_volume import peak
        
        sc = score()
        # check auto-correlation coefficient
        c  = sc.scoringCoefficient( self.s, self.s)
        self.assertAlmostEqual( first=c,  second=1., places=5, 
             msg='SocScore: Auto-correlation not == 1')
        # consistency of scoring coefficient and scoring function - difference 
        # due to sub-pixel accuracy for score
        c  = sc.scoringCoefficient( self.s, self.v)
        cf = sc.scoringFunction( self.s, self.v)
        p= peak(cf)
        self.assertAlmostEqual( first=c, second=cf.getV(p[0],p[1],p[2]), places=2, 
            msg='Scoring coefficient and scoring function FLCF inconsistent')
        
    def test_socScore(self):
        """
        second order correlation score
        """
        from pytom.score.score import SOCScore as score
        from pytom_volume import peak

        sc = score()
        # check auto-correlation coefficient
        c  = sc.scoringCoefficient( self.s, self.s)
        self.assertAlmostEqual( first=c,  second=1., places=5, 
             msg='SocScore: Auto-correlation not == 1')
        # consistency of scoring coefficient and scoring function - difference 
        # due to sub-pixel accuracy for score
        c  = sc.scoringCoefficient( self.s, self.v)
        cf = sc.scoringFunction( self.s, self.v)
        p= peak(cf)
        self.assertAlmostEqual( first=c, second=cf.getV(p[0],p[1],p[2]), places=2, 
            msg='Scoring coefficient and scoring function SOC inconsistent')
    
    def RScore_Test(self):
        """
        """
        
    def runTest(self):
        self.test_xcfScore()
        self.test_nxcfScore()
        self.test_flcfScore()
        self.test_socScore()

if __name__ == '__main__':
    unittest.main()
