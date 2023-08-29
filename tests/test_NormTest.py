import unittest;
epsilon = 0.0001;


class pytom_NormTest(unittest.TestCase):
    
    def mean0std1_Test(self):

        import pytom.lib.pytom_volume as pytom_volume
        from pytom.basic.normalise import mean0std1
        from helper_functions import create_RandomParticleList, installdir

        v = pytom_volume.read(f'{installdir}/../tests/testData/ribo.em')
        
        m = mean0std1(v,True)
        
        me = pytom_volume.mean(m)
        var = pytom_volume.variance(m,False)
        
        assert epsilon >= me >= -epsilon #mean value test 
        assert 1+epsilon >= var >= 1-epsilon

    def zeroVol_Test(self):
        import pytom.lib.pytom_volume as pytom_volume
        from pytom.basic.normalise import normaliseUnderMask
        from pytom.basic.correlation import stdValueUnderMask, meanValueUnderMask

        vol = pytom_volume.vol(64,64,64)
        vol.setAll(.0)
        vol.setV(1., 32, 32, 32)
        mask = pytom_volume.vol(64,64,64)
        pytom_volume.initSphere(mask,13,0,0,33,33,33)
        normalised, p = normaliseUnderMask(vol, mask)
        std = stdValueUnderMask(normalised, mask, meanValueUnderMask(normalised, mask, p), p)
        self.assertTrue(1 + epsilon >= std >= 1 - epsilon)

    def runTest(self):
        self.mean0std1_Test()
        self.zeroVol_Test()


if __name__ == '__main__':
    unittest.main()

