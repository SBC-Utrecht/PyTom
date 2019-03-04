import unittest;

class pytom_SimulationTest(unittest.TestCase):

    def gaussianNoise_Test(self):
        eps = 0.1;
        from pytom_volume import vol,mean,variance,gaussianNoise;
        from math import sqrt;
        
        m = 0;
        s = 1;
        
        v = vol(10,10,10);
        v.setAll(0);
        gaussianNoise(v,float(m),float(s));
        
        mv = mean(v);
        sv = sqrt(variance(v,False));
        
        assert abs(0+mv) <=abs(m-eps);
        assert sv <=abs(s+eps);
        