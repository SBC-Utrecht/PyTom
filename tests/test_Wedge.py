"""
Created on Jan 12, 2022

@author: Marten Chaillet
"""
import unittest
import numpy as np


class pytom_NumpyTest(unittest.TestCase):
    def setUp(self):
        self.SX, self.SY, self.SZ = 100, 100, 90
        self.SZ_uneven = 91
        self.SX_uneven = 101

    def wedge_pytom_vs_numpy(self):
        from pytom.basic.structures import Wedge
        from pytom.agnostic.filter import create_wedge
        from pytom_numpy import vol2npy
        w1, w2 = 30., 30.
        w = Wedge([w1, w2])
        smooth = w._wedgeObject._smooth

        # here there is no problem
        SX, SY, SZ = self.SX, self.SY, self.SZ
        temp = w.returnWedgeVolume(SX, SY, SZ)
        wedge_pytom = vol2npy(temp).copy().astype(np.float32)
        wedge_numpy = create_wedge(w1, w2, SX // 2, SX, SY, SZ, smooth).astype(np.float32)
        self.assertTrue((wedge_pytom != wedge_numpy).sum() == 0, "now there is a big issue")

        # here there is again a problem with uneven z height
        SX, SY, SZ = self.SX_uneven, self.SY, self.SZ
        temp = w.returnWedgeVolume(SX, SY, SZ)
        wedge_pytom = vol2npy(temp).copy().astype(np.float32)
        wedge_numpy = create_wedge(w1, w2, SX // 2, SX, SY, SZ, smooth).astype(np.float32)
        self.assertTrue((wedge_pytom != wedge_numpy).sum() == 0, "same, this was also working ")

        #
        SX, SY, SZ = self.SX, self.SY, self.SZ_uneven
        temp = w.returnWedgeVolume(SX, SY, SZ)
        wedge_pytom = vol2npy(temp).copy().astype(np.float32)
        wedge_numpy = create_wedge(w1, w2, SX // 2, SX, SY, SZ, smooth).astype(np.float32)
        self.assertTrue((wedge_pytom != wedge_numpy).sum() == 7843, "something changed with this know issue")

    def test_conversion(self):
        self.wedge_pytom_vs_numpy()


if __name__ == '__main__':
    unittest.main()
