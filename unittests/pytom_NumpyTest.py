'''
Created on May 4, 2010

@author: hrabe
'''
import unittest;



class pytom_NumpyTest(unittest.TestCase):
    
    
    def numpy_T(self,v=None):
        
        import pytom_numpy
        
        if not v:
            from pytom_volume import vol,gaussianNoise
            
            v = vol(10,10,10)
            gaussianNoise(v,0,1)
            
        nv = pytom_numpy.vol2npy(v)
        #print 'hier'
        #print nv.shape,nv.strides,nv.dtype
        #print v.strideX(),v.strideY(),v.strideZ()
        vv = pytom_numpy.npy2vol(nv,3)

        assert v.equalsTo(vv)
        
            
    def numpy_Test(self):
        from pytom_volume import vol,gaussianNoise
        self.numpy_T()
        
        """
        v = vol(512,512,128)
        gaussianNoise(v,0,1)
        print 'in'
        self.numpy_T(v)
        print 'out'
        """