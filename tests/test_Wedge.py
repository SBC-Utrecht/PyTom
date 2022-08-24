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
        self.SY_uneven = 99
        self.smooth = 3

    def wedge_pytom_vs_numpy(self, w1=30., w2=30.):
        from pytom.basic.structures import Wedge
        from pytom.agnostic.filter import create_wedge
        from pytom_numpy import vol2npy
        w = Wedge([w1, w2])
        smooth = w._wedgeObject._smooth

        # here there is no problem
        SX, SY, SZ = self.SX, self.SY, self.SZ
        temp = w.returnWedgeVolume(SX, SY, SZ)
        wedge_pytom = vol2npy(temp).copy().astype(np.float32)
        wedge_numpy = create_wedge(w1, w2, SX // 2, SX, SY, SZ, smooth).astype(np.float32)
        self.assertTrue((wedge_pytom != wedge_numpy).sum() == 0, "now there is a big issue")

        # test uneven X
        SX, SY, SZ = self.SX_uneven, self.SY, self.SZ
        temp = w.returnWedgeVolume(SX, SY, SZ)
        wedge_pytom = vol2npy(temp).copy().astype(np.float32)
        wedge_numpy = create_wedge(w1, w2, SX // 2, SX, SY, SZ, smooth).astype(np.float32)
        self.assertTrue((wedge_pytom != wedge_numpy).sum() == 0, "failing with uneven X")

        # test uneven Y
        SX, SY, SZ = self.SX, self.SY_uneven, self.SZ
        temp = w.returnWedgeVolume(SX, SY, SZ)
        wedge_pytom = vol2npy(temp).copy().astype(np.float32)
        wedge_numpy = create_wedge(w1, w2, SX // 2, SX, SY, SZ, smooth).astype(np.float32)
        self.assertTrue((wedge_pytom != wedge_numpy).sum() == 0, "failing with uneven Y")

        # test uneven Z
        run = False
        if run:
            SX, SY, SZ = self.SX, self.SY, self.SZ_uneven
            temp = w.returnWedgeVolume(SX, SY, SZ)
            wedge_pytom = vol2npy(temp).copy().astype(np.float32)
            wedge_numpy = create_wedge(w1, w2, SX // 2, SX, SY, SZ, smooth).astype(np.float32)
            # self.assertTrue((wedge_pytom != wedge_numpy).sum() == 7843, "failing with uneven Z")
            self.assertTrue((wedge_pytom != wedge_numpy).sum() == 0, "failing with uneven Z")

    def test_conversion(self):
        wedges = [(30., 30.), (10., 80.), (34.23, 51.02)]
        for w in wedges:
            print(w)
            self.wedge_pytom_vs_numpy(w1=w[0], w2=w[1])


if __name__ == '__main__':
    unittest.main()
