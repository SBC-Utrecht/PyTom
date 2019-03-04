import unittest;
epsilon = 0.0001;


class pytom_NormTest(unittest.TestCase):
    
    def mean0std1_Test(self):

        import pytom_volume
        from pytom.basic.normalise import mean0std1
        v = pytom_volume.read('./testData/ribo.em')
        
        m = mean0std1(v,True)
        
        me = pytom_volume.mean(m)
        var = pytom_volume.variance(m,False)
        
        assert epsilon >= me >= -epsilon #mean value test 
        assert 1+epsilon >= var >= 1-epsilon
