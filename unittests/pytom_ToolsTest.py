'''
Created on May 17, 2010

@author: hrabe
'''


import unittest;

class pytom_ToolsTest(unittest.TestCase):
    
    def VolumeStorage_Test(self):
        from pytom.tools.memory import read
        from pytom_volume import vol,gaussianNoise
        import os
        
        noise = vol(10,10,10)
        gaussianNoise(noise,0,1)
        
        noise.write('./test.em')
        
        noise2 = read('./test.em',0,0,0,0,0,0,0,0,0,0,0,0)
        
        assert noise.equalsTo(noise2)
        noise3 = read('./test.em',0,0,0,0,0,0,0,0,0,0,0,0)
        
        
        assert noise2.equalsTo(noise3)
        
        os.system('rm ./test.em')
        
        
        
    