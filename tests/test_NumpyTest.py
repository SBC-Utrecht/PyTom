"""
Created on Jan 12, 2022

@author: Marten Chaillet
"""
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
        self.assertTrue((self.random_npy != nv).sum() == 0, msg='numpy issue :(')

    def known_issue(self):
        from pytom.basic.structures import Wedge
        from pytom.agnostic.filter import create_wedge
        w1, w2 = 30., 30.
        w = Wedge([w1, w2])
        sx, sy, sz = 12, 12, 12
        smooth = w._wedgeObject._smooth

        # no issues
        temp = w.returnWedgeVolume(sx, sy, sz)
        wedge_pytom = pytom_numpy.vol2npy(temp).copy().astype(np.complex64)
        wedge_numpy = create_wedge(w1, w2, sx // 2, sx, sy, sz, smooth).astype(np.complex64)
        self.assertTrue((wedge_pytom != wedge_numpy).sum() == 0, msg='somehting is very wrong with conversion')

        # here there is an issue
        # this is just a known location where I know there is a mismatch of 8 (McHaillet)
        wedge_pytom = pytom_numpy.vol2npy(w.returnWedgeVolume(sx, sy, sz)).copy().astype(np.complex64)
        wedge_numpy = create_wedge(w1, w2, sx // 2, sx, sy, sz, smooth).astype(np.complex64)
        self.assertTrue((wedge_pytom != wedge_numpy).sum() == 8, msg='this is a known issue')

    def test_conversion(self):
        self.forward()
        self.backward()
        self.known_issue()


if __name__ == '__main__':
    unittest.main()
