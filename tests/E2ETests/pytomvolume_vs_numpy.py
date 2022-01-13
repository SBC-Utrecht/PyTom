"""
Test CTF with astigmatism and rectangular images shapes.

@author: Marten Chaillet
"""
import unittest
import numpy as np
from pytom.agnostic.correlation import xcc, nxcc


class NumericalTest(unittest.TestCase):
    def setUp(self):
        from pytom.agnostic.io import read
        from pytom.agnostic.tools import paste_in_center
        from pytom_numpy import npy2vol
        self.wedge_angles = [30., 30.]
        self.dims = (100, 100, 71)
        self.ndims = 3
        self.cutoff = self.dims[0] // 2

        # initiate random values for testing
        self.random_np = np.random.rand(*self.dims).astype(np.float32)
        self.random_copy = self.random_np.copy(order='F')
        self.random_vol = npy2vol(self.random_copy, self.ndims)

        # initiate mask for stdv
        self.mask = np.zeros_like(self.random_np, dtype=np.float32)
        self.mask_np = paste_in_center(read('../testData/human_ribo_mask_32_8_5.mrc').astype(np.float32), self.mask)
        self.mask_copy = self.mask_np.copy(order='F')
        self.mask_vol = npy2vol(self.mask_copy, self.ndims)

        # rotation reference
        self.reference_np = read('../testData/human_ribo_21A.em').astype(np.float32)
        self.reference_copy = self.reference_np.copy(order='F')
        self.reference_vol = npy2vol(self.reference_copy, self.ndims)

        # empty place to put the rotation
        self.template_np = np.zeros_like(self.reference_np, dtype=np.float32)
        self.template_copy = self.template_np.copy(order='F')
        self.template_vol = npy2vol(self.template_copy, self.ndims)

        u = np.random.uniform(0.0, 1.0, (2,))
        theta = np.arccos(2 * u[0] - 1)
        phi = 2 * np.pi * u[1]
        psi = np.random.uniform(0.0, 2 * np.pi)
        self.rotation_angles = (theta, phi, psi)

    def wedgePyTomVol(self):
        from pytom.basic.structures import Wedge
        from pytom_numpy import vol2npy
        w = Wedge(wedgeAngles=self.wedge_angles, cutoffRadius=self.cutoff)
        wedge = w.returnWedgeVolume(*self.dims)
        wedge_np = vol2npy(wedge).copy()
        return wedge_np.astype(np.float32)

    def wedgeNpCp(self):
        from pytom.agnostic.filter import create_wedge
        return create_wedge(*self.wedge_angles, self.cutoff, *self.dims).astype(np.float32)
    # Run each test for numpy and cupy separetely!

    def fftPyTomVol(self):
        from pytom_volume import reducedToFull
        from pytom_numpy import vol2npy
        from pytom.basic.fourier import fft
        ft = reducedToFull(fft(self.random_vol))
        ft_np = vol2npy(ft).copy()
        return ft_np.astype(np.complex64)

    def fftNpCp(self):
        return np.fft.fftn(self.random_np).astype(np.complex64)

    def rotatePyTomVol(self):
        from pytom_volume import rotateSpline as rotate  # for more accuracy
        from pytom_numpy import vol2npy
        rotate(self.reference_vol, self.template_vol, *self.rotation_angles)
        rotated_np = vol2npy(self.template_vol).copy()
        return rotated_np.astype(np.float32)

    def rotateVoltools(self):
        import pytom.voltools as vt
        self.texture = vt.StaticVolume(self.reference_np, interpolation='filt_bspline', device='cpu')
        self.texture.transform(rotation=self.rotation_angles, rotation_order='rzxz', output=self.template_np)
        return self.template_np

    def calcStdvPyTom(self):
        from pytom.basic.correlation import meanUnderMask, stdUnderMask
        from pytom_numpy import vol2npy
        from pytom_volume import sum
        p = sum(self.mask_vol)
        meanV = meanUnderMask(self.random_vol, self.mask_vol, p)
        stdV = stdUnderMask(self.random_vol, self.mask_vol, p, meanV)
        stdV_np = vol2npy(stdV).copy()
        return stdV_np.astype(np.float32)

    def calcStdvNpCp(self):
        from pytom.agnostic.normalise import meanVolUnderMask, stdVolUnderMask
        meanV = meanVolUnderMask(self.random_np, self.mask_np)
        stdV = stdVolUnderMask(self.random_np, self.mask_np, meanV)
        return stdV

    def rmsd(self, pytomvol, npcp):
        return np.sqrt(((pytomvol - npcp) ** 2).sum() / npcp.size)

    def wedge(self):
        w_pytomvol = self.wedgePyTomVol()
        w_npcp = self.wedgeNpCp()
        if np.all(w_pytomvol == w_npcp):
            print('pytom_volume wedge is identical to numpy/cupy wedge')
        ccc = nxcc(w_pytomvol, w_npcp)
        print('wedge ccc: ', ccc)
        print('wedge rmsd: ', self.rmsd(w_pytomvol, w_npcp))
        # self.display_diff_vol(np.abs(w_pytomvol - w_npcp))
        self.assertAlmostEqual(ccc, 1.0, places=4, msg='correlation not sufficient between pytom_volume wedge and '
                                                       'numpy/cupy wedge')

    def fft(self):
        ft_pytomvol = self.fftPyTomVol()
        ft_npcp = self.fftNpCp()
        # self.display_fft(ft_pytomvol, ft_npcp)
        ccc = nxcc(np.abs(ft_pytomvol), np.abs(ft_npcp))  # correlate abs, because nxcc does not accept complex
        print('fft ccc: ', ccc)
        print('fft rmsd real: ', self.rmsd(ft_pytomvol.real, ft_npcp.real))
        print('fft rmsd imag: ', self.rmsd(ft_pytomvol.imag, ft_npcp.imag))
        # self.display_diff_fft(np.abs(ft_pytomvol.real - ft_npcp.real), np.abs(ft_pytomvol.imag - ft_npcp.imag))
        self.assertAlmostEqual(ccc, 1.0, places=4, msg='correlation not sufficient between pytom_volume fft and '
                                                       'numpy/cupy fft')

    def rotation(self):
        rt_pytomvol = self.rotatePyTomVol()
        rt_npcp = self.rotateVoltools()
        # self.display_vol(rt_pytomvol, rt_npcp)
        ccc = nxcc(rt_pytomvol, rt_npcp)
        print('rotation ccc: ', ccc)
        print('rotation rmsd: ', self.rmsd(rt_pytomvol, rt_npcp))
        # self.display_diff_vol(np.abs(rt_pytomvol - rt_npcp))
        self.assertAlmostEqual(ccc, 1.0, places=1, msg='correlation not sufficient between pytom_volume rotateSpline '
                                                       'and voltools cpu rotate')

    def calcStdv(self):
        std_pytomvol = self.calcStdvPyTom()
        std_npcp = self.calcStdvNpCp()
        # self.display_vol(std_pytomvol, std_npcp)
        ccc = nxcc(std_pytomvol, std_npcp)
        print('stdv ccc: ', ccc)
        print('stdv rmsd: ', self.rmsd(std_pytomvol, std_npcp))
        # self.display_diff_vol(np.abs(std_pytomvol - std_npcp))
        self.assertAlmostEqual(ccc, 1.0, places=4, msg='correlation not sufficient between pytom_volume calcStdv '
                                                       'and voltools cpu calcStdv')

    def display_fft(self, v1, v2):
        import matplotlib
        matplotlib.use('Qt5Agg')
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 2)
        ax[0].imshow(np.log(np.abs(v1[..., 0])))
        ax[1].imshow(np.log(np.abs(v2[..., 0])))
        ax[0].set_title('pytom_volume')
        ax[1].set_title('numpy/cupy')
        plt.show()

    def display_vol(self, v1, v2):
        import matplotlib
        matplotlib.use('Qt5Agg')
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 2)
        ax[0].imshow(v1[..., v1.shape[2] // 2])
        ax[1].imshow(v2[..., v2.shape[2] // 2])
        ax[0].set_title('pytom_volume')
        ax[1].set_title('numpy/cupy')
        plt.show()

    def display_diff_vol(self, v1):
        import matplotlib
        matplotlib.use('Qt5Agg')
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        ax.imshow(v1[...,v1.shape[2] // 2])
        plt.show()

    def display_diff_fft(self, real, imag):
        import matplotlib
        matplotlib.use('Qt5Agg')
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 2)
        ax[0].imshow(real[..., 0])
        ax[1].imshow(imag[..., 0])
        ax[0].set_title('diff real part')
        ax[1].set_title('diff imag part')
        plt.show()

    def test_stability(self):
        self.wedge()
        self.rotation()
        self.fft()
        self.calcStdv()


if __name__ == '__main__':
    unittest.main()
