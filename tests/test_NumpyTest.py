'''
Created on May 4, 2010

@author: hrabe
'''
import unittest
import numpy as np
from pytom_volume import vol, gaussianNoise
import pytom_numpy


class pytom_NumpyTest(unittest.TestCase):
    def setUp(self):
        from functools import reduce
        size = (10, 11, 12)
        self.N = reduce(lambda x, y: x * y, size)
        self.random_vol = vol(*size)
        gaussianNoise(self.random_vol, 0, 1)
        self.random_npy = np.asfortranarray(np.random.normal(0, 1, size).astype(np.float32))
    
    def forward(self):
        nv = pytom_numpy.vol2npy(self.random_vol).copy(order='F')
        pv = pytom_numpy.npy2vol(nv, len(nv.shape))
        self.assertTrue(self.random_vol.equalsTo(pv), msg='numpy issue :(')

    def backward(self):
        copy = self.random_npy.copy(order='F')
        pv = pytom_numpy.npy2vol(copy, len(copy.shape))
        nv = pytom_numpy.vol2npy(pv)
        print((self.random_npy != nv).sum())
        self.assertTrue((self.random_npy != nv).sum() == 0, msg='numpy issue :(')
            
    def test_conversion(self):
        self.forward()
        self.backward()


if __name__ == '__main__':
    unittest.main()
