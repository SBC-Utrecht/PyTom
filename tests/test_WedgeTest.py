"""
Created on Jan 12, 2022

@author: Marten Chaillet
"""
import unittest
import numpy as np
from pytom.basic.structures import Wedge, Rotation
from pytom.agnostic.correlation import nxcc
from pytom.agnostic.structures import Wedge as WedgeNp
from pytom.lib.pytom_numpy import vol2npy


class pytom_NumpyTest(unittest.TestCase):
    def setUp(self):
        self.SX, self.SY, self.SZ = 22, 22, 20
        self.SZ_uneven = 19
        self.SX_uneven = 23
        self.SY_uneven = 19
        self.smooth = 2
        self.rotation = Rotation(23, 9., 51.222)

    def wedge_pytom_vs_numpy_red(self, w1=30., w2=30., known_error=0):
        w = Wedge([w1, w2])
        w_np = WedgeNp([w1, w2])

        # here there is no problem
        SX, SY, SZ = self.SX, self.SY, self.SZ
        temp = w.returnWedgeVolume(SX, SY, SZ)
        wedge_pytom = vol2npy(temp).copy().astype(np.float32)
        wedge_numpy = w_np.returnWedgeVolume(SX, SY, SZ)
        self.assertTrue((wedge_pytom != wedge_numpy).sum() == 0, "now there is a big issue")

        # test uneven X
        SX, SY, SZ = self.SX_uneven, self.SY, self.SZ
        temp = w.returnWedgeVolume(SX, SY, SZ)
        wedge_pytom = vol2npy(temp).copy().astype(np.float32)
        wedge_numpy = w_np.returnWedgeVolume(SX, SY, SZ)
        self.assertTrue((wedge_pytom != wedge_numpy).sum() == known_error[0], "failing with uneven X")

        # test uneven Y
        SX, SY, SZ = self.SX, self.SY_uneven, self.SZ
        temp = w.returnWedgeVolume(SX, SY, SZ)
        wedge_pytom = vol2npy(temp).copy().astype(np.float32)
        wedge_numpy = w_np.returnWedgeVolume(SX, SY, SZ)
        self.assertTrue((wedge_pytom != wedge_numpy).sum() == known_error[1], "failing with uneven Y")

        # test uneven Z
        SX, SY, SZ = self.SX, self.SY, self.SZ_uneven
        temp = w.returnWedgeVolume(SX, SY, SZ)
        wedge_pytom = vol2npy(temp).copy().astype(np.float32)
        wedge_numpy = w_np.returnWedgeVolume(SX, SY, SZ)
        self.assertTrue((wedge_pytom != wedge_numpy).sum() == known_error[2], "failing with uneven Z")

        # FOR DEVELOPERS: set to true if you want to plot the wedges
        # shows the known issue for uneven z-height for pytom_volume
        run = False
        if run and w1 == w2:
            import matplotlib
            try:
                matplotlib.use('Qt5Agg')
            except:
                pass
            import matplotlib.pyplot as plt

            SX, SY, SZ = self.SX_uneven, self.SY, self.SZ
            temp = w.returnWedgeVolume(SX, SY, SZ)
            wedge_pytom = vol2npy(temp).copy().astype(np.float32)
            wedge_numpy = w_np.returnWedgeVolume(SX, SY, SZ)
            fig, ax = plt.subplots(1, 3)
            ax[0].imshow(wedge_numpy[:, 0, :])
            ax[0].set_title('agnostic wedge')
            ax[1].imshow(wedge_pytom[:, 0, :])
            ax[1].set_title('pytom_volume wedge')
            ax[2].imshow(wedge_pytom[:, 0, :] - wedge_numpy[:, 0, :])
            ax[2].set_title('difference')
            plt.show()

            SX, SY, SZ = self.SX, self.SY, self.SZ_uneven
            temp = w.returnWedgeVolume(SX, SY, SZ)
            wedge_pytom = vol2npy(temp).copy().astype(np.float32)
            wedge_numpy = w_np.returnWedgeVolume(SX, SY, SZ)
            fig, ax = plt.subplots(1, 3)
            ax[0].imshow(wedge_numpy[:, 0, :])
            ax[0].set_title('agnostic wedge')
            ax[1].imshow(wedge_pytom[:, 0, :])
            ax[1].set_title('pytom_volume wedge')
            ax[2].imshow(wedge_pytom[:, 0, :] - wedge_numpy[:, 0, :])
            ax[2].set_title('difference')
            plt.show()

    def wedge_pytom_vs_numpy_full(self, w1=30., w2=30., known_error=(0,0)):
        w = Wedge([w1, w2])
        w_np = w.convert2numpy()

        # here there is no problem
        SX, SY, SZ = self.SX, self.SY, self.SZ
        temp = w.returnWedgeVolume(SX, SY, SZ, humanUnderstandable=True)
        wedge_pytom = vol2npy(temp).copy().astype(np.float32)
        wedge_numpy = w_np.returnWedgeVolume(SX, SY, SZ, humanUnderstandable=True)
        self.assertTrue((wedge_pytom != wedge_numpy).sum() == 0, "now there is a big issue")

        # test uneven X
        SX, SY, SZ = self.SX_uneven, self.SY, self.SZ
        temp = w.returnWedgeVolume(SX, SY, SZ, humanUnderstandable=True)
        wedge_pytom = vol2npy(temp).copy().astype(np.float32)
        wedge_numpy = w_np.returnWedgeVolume(SX, SY, SZ, humanUnderstandable=True)
        self.assertTrue((wedge_pytom != wedge_numpy).sum() == known_error[0], "failing with uneven X")

        # test uneven Y
        SX, SY, SZ = self.SX, self.SY_uneven, self.SZ
        temp = w.returnWedgeVolume(SX, SY, SZ, humanUnderstandable=True)
        wedge_pytom = vol2npy(temp).copy().astype(np.float32)
        wedge_numpy = w_np.returnWedgeVolume(SX, SY, SZ, humanUnderstandable=True)
        self.assertTrue((wedge_pytom != wedge_numpy).sum() == known_error[1], "failing with uneven Y")

    def smoothed_wedge(self, w1, w2):
        w = Wedge([w1, w2], smooth=self.smooth)
        w_np = w.convert2numpy()

        SX, SY, SZ = self.SX, self.SY, self.SZ
        temp = w.returnWedgeVolume(SX, SY, SZ)
        wedge_pytom = vol2npy(temp).copy().astype(np.float32)
        wedge_numpy = w_np.returnWedgeVolume(SX, SY, SZ)
        self.assertTrue(nxcc(wedge_pytom, wedge_numpy) > 0.7, "now there is a big issue")

        # test uneven X
        SX, SY, SZ = self.SX_uneven, self.SY, self.SZ
        temp = w.returnWedgeVolume(SX, SY, SZ)
        wedge_pytom = vol2npy(temp).copy().astype(np.float32)
        wedge_numpy = w_np.returnWedgeVolume(SX, SY, SZ)
        self.assertTrue(nxcc(wedge_pytom, wedge_numpy) > 0.7, "failing with uneven X")

        # test uneven Y
        SX, SY, SZ = self.SX, self.SY_uneven, self.SZ
        temp = w.returnWedgeVolume(SX, SY, SZ)
        wedge_pytom = vol2npy(temp).copy().astype(np.float32)
        wedge_numpy = w_np.returnWedgeVolume(SX, SY, SZ)
        self.assertTrue(nxcc(wedge_pytom, wedge_numpy) > 0.7, "failing with uneven Y")

    def rotated_wedge(self, w1, w2):
        """
        THIS TEST IS SKIPPED

        Both pytom_volume and agnostic wedge produce incorrect results with the rotation parameter.
        Do not use this option in your functions. Generate the wedge without rotation, and then rotate
        it with some interpolation function.
        """
        # from pytom.agnostic.io import write
        w = Wedge([w1, w2], smooth=self.smooth)
        w_np = w.convert2numpy()

        SX, SY, SZ = self.SX, self.SY, self.SZ
        temp = w.returnWedgeVolume(SX, SY, SZ, rotation=self.rotation)
        wedge_pytom = vol2npy(temp).copy().astype(np.float32)
        wedge_numpy = w_np.returnWedgeVolume(SX, SY, SZ, rotation=self.rotation)
        # write('wedge_pytom.mrc', np.fft.ifftshift(wedge_pytom, axes=(0, 1)))
        # write('wedge_numpy.mrc', np.fft.ifftshift(wedge_numpy, axes=(0, 1)))
        print(nxcc(wedge_pytom, wedge_numpy))
        self.assertTrue(nxcc(wedge_pytom, wedge_numpy) > 0.90, "now there is a big issue")

        # test uneven X
        SX, SY, SZ = self.SX_uneven, self.SY, self.SZ
        temp = w.returnWedgeVolume(SX, SY, SZ, rotation=self.rotation)
        wedge_pytom = vol2npy(temp).copy().astype(np.float32)
        wedge_numpy = w_np.returnWedgeVolume(SX, SY, SZ, rotation=self.rotation)
        print(nxcc(wedge_pytom, wedge_numpy))
        self.assertTrue(nxcc(wedge_pytom, wedge_numpy) > 0.90, "failing with uneven X")

        # test uneven Y
        SX, SY, SZ = self.SX, self.SY_uneven, self.SZ
        temp = w.returnWedgeVolume(SX, SY, SZ, rotation=self.rotation)
        wedge_pytom = vol2npy(temp).copy().astype(np.float32)
        wedge_numpy = w_np.returnWedgeVolume(SX, SY, SZ, rotation=self.rotation)
        print(nxcc(wedge_pytom, wedge_numpy))
        self.assertTrue(nxcc(wedge_pytom, wedge_numpy) > 0.90, "failing with uneven Y")

        # FOR DEVELOPERS: set to true if you want to plot the wedges
        # shows the known issue for uneven z-height for pytom_volume
        run = True
        if run and w1 == w2:
            import matplotlib
            try:
                matplotlib.use('Qt5Agg')
            except:
                pass
            import matplotlib.pyplot as plt

            SX, SY, SZ = self.SX, self.SY, self.SZ
            temp = w.returnWedgeVolume(SX, SY, SZ, rotation=self.rotation)
            wedge_pytom = vol2npy(temp).copy().astype(np.float32)
            wedge_numpy = w_np.returnWedgeVolume(SX, SY, SZ, rotation=self.rotation)
            fig, ax = plt.subplots(1, 3)
            ax[0].imshow(wedge_numpy[:, 0, :])
            ax[0].set_title('agnostic wedge')
            ax[1].imshow(wedge_pytom[:, 0, :])
            ax[1].set_title('pytom_volume wedge')
            ax[2].imshow(wedge_pytom[:, 0, :] - wedge_numpy[:, 0, :])
            ax[2].set_title('difference')
            plt.show()

    def test_conversion(self):
        # (wedge 1, wedge 2, known error)
        wedges = [(30., 30., 365), (5., 85., 317), (34.23, 51.02, 375)]
        error_red = [(630, 174, 365), (394, 132, 317), (514, 138, 375)]  # number of pixels that have different values
        error_full = [(1090, 312), (668, 212), (920, 240)]
        for w, errred, errfull in zip(wedges, error_red, error_full):
            # ======== reduced test
            self.wedge_pytom_vs_numpy_red(w1=w[0], w2=w[1], known_error=errred)
            # ======== full wedge test
            self.wedge_pytom_vs_numpy_full(w1=w[0], w2=w[1], known_error=errfull)
            # ======== smoothed wedge test
            self.smoothed_wedge(w1=w[0], w2=w[1])
            # self.rotated_wedge(w1=w[0], w2=w[1])


if __name__ == '__main__':
    unittest.main()
