"""
Test CTF with astigmatism and rectangular images shapes.

@author: Marten Chaillet
"""
import unittest
import numpy as np

class CTFTest(unittest.TestCase):
    def setUp(self):
        """Initialize simulation parameters"""
        from pytom.simulation.microscope import convert_defocusU_defocusV_to_defocus_astigmatism as convert
        from pytom.simulation.potential import create_gold_marker
        import pytom.simulation.physics as phys
        import os

        self.voltage = 300e3
        self.pixel_size_angstrom = 2
        real, imag = create_gold_marker(self.pixel_size_angstrom, phys.V_WATER, oversampling=1, solvent_factor=1.0,
                                       imaginary=True, voltage=self.voltage)
        self.gold_bead = real + 1j * imag
        self.gold_size = self.gold_bead.shape

        self.x, self.y, self.z = 1000, 1000, 300
        self.potential = np.zeros((self.x, self.y, self.z), dtype=self.gold_bead.dtype)
        self.potential[self.x//2:self.x//2 + self.gold_size[0], self.y//2:self.y//2+self.gold_size[1],
                                self.z//2:self.z//2+self.gold_size[2]] = self.gold_bead

        self.defocusU, self.defocusV, self.ast_angle_deg = 5e-6, 7e-6, -69
        self.defocus, self.astigmatism = convert(self.defocusU, self.defocusV)
        self.Cs = 2.7e-3
        self.amplitude_contrast = 0.07
        self.sigma_decay = 0.4

        # create temporary dir for storing simulation data
        if not os.path.exists('temp_simulation'):
            os.mkdir('temp_simulation')

        # specific defocus and msdz, but otherwise default parameters for ctf function
        self.param_sim = {
            'save_path':            './temp_simulation',
            'angles':               [0],
            'nodes':                1,
            'pixel_size':           self.pixel_size_angstrom*1e-10,
            'binning':              2,
            'dose':                 100,
            'voltage':              300e3,
            'defocus':              self.defocus,
            'astigmatism':          self.astigmatism,
            'ast_angle':            self.ast_angle_deg,
            'msdz':                 5e-9,
            'camera_type':          'K2SUMMIT',
            'camera_folder':        '../simulation/detectors',
        }

    def tearDown(self):
        """Remove all the files gathered during simulation"""
        from os import path
        dir = self.param_sim['save_path']
        self.remove_file(path.join(dir, 'projections.mrc'))
        self.remove_file(path.join(dir, 'noisefree_projections.mrc'))
        self.remove_file(path.join(dir, 'alignment_simulated.txt'))
        # self.remove_file(path.join(dir, 'reconstruction.mrc'))
        self.remove_tree(dir)

    def remove_tree(self, foldername):
        """Assert folder exists, then remove its content and itself"""
        from os import path
        from shutil import rmtree

        foldercheck = path.exists(foldername)
        if not foldercheck:
            print(foldername + " does not exist")
        self.assertTrue(foldercheck, msg="folder " + foldername + " does not exist")
        if foldercheck:
            rmtree(foldername)

    def remove_file(self, filename):
        """Assert that file exists, then remove it"""
        from os import remove
        from os import path

        filecheck = path.exists(filename)
        if not filecheck:
            print(filename + " does not exist")
        self.assertTrue(filecheck, msg="file " + filename + " does not exist")
        if filecheck:
            remove(filename)

    def simulateProjection(self, c=''):
        """Run the simulation, output here will be written to some temp storage"""
        from pytom.simulation.MicrographModeller import generate_tilt_series_cpu
        from pytom.tompy.io import read
        import os

        if not os.path.exists(self.param_sim['save_path'] + c):
            os.mkdir(self.param_sim['save_path'] + c)

        generate_tilt_series_cpu(self.param_sim['save_path'] + c,
                                 self.param_sim['angles'],
                                 nodes=self.param_sim['nodes'],
                                 pixel_size=self.param_sim['pixel_size'],
                                 binning=self.param_sim['binning'],
                                 dose=self.param_sim['dose'],
                                 voltage=self.param_sim['voltage'],
                                 defocus=self.param_sim['defocus'],
                                 astigmatism=self.astigmatism,
                                 astigmatism_angle=self.ast_angle_deg,
                                 msdz=self.param_sim['msdz'],
                                 camera_type=self.param_sim['camera_type'],
                                 camera_folder=self.param_sim['camera_folder'],
                                 grandcell=self.potential)

        return read(os.path.join(self.param_sim['save_path'] + c, 'projections.mrc')).squeeze()

    def basicCTF(self):
        """Convolute gold marker with a ctf, then use a wiener filter on this marker. CCC is calculated between
        original gold bead and the wiener filtered one."""
        from pytom.tompy.correlation import nxcc
        from pytom.simulation.microscope import create_ctf

        # import matplotlib
        # matplotlib.use('Qt5Agg')
        # import matplotlib.pyplot as plt

        # create an image of a point
        point = self.potential[:, :, self.z//2 + self.gold_size[2] // 2].real
        point = point[2:,101:]  # crop with annoying lengths to ensure the ctf function still works

        # get a ctf
        ctf = create_ctf(point.shape,
                         self.param_sim['pixel_size'],
                         self.defocusU,
                         self.amplitude_contrast,
                         self.param_sim['voltage'],
                         self.Cs,
                         sigma_decay=self.sigma_decay,
                         defocusV=self.defocusV,
                         defocus_angle_deg=self.ast_angle_deg)

        point_conv = np.fft.ifftn(np.fft.fftn(point) * np.fft.ifftshift(ctf)).real

        ctf_corr = create_ctf(point.shape,
                         self.param_sim['pixel_size'],
                         self.defocusU,
                         self.amplitude_contrast,
                         self.param_sim['voltage'],
                         self.Cs,
                         sigma_decay=0,
                         defocusV=self.defocusV,
                         defocus_angle_deg=self.ast_angle_deg)

        ctf_corr = np.fft.ifftshift(ctf_corr)

        point_corrected = np.fft.ifftn(np.conjugate(ctf_corr) * np.fft.fftn(
            point_conv) / (ctf_corr ** 2 + 0.1)).real

        # fig, ax = plt.subplots(1, 3)
        # ax[0].imshow(point)
        # ax[1].imshow(point_conv)
        # ax[2].imshow(point_corrected)
        # plt.show()

        ccc = nxcc(point, point_corrected)

        # test if correlation is within the limits
        # print('cross correlation of imaged point source with a ctf = ', ccc)
        return ccc

    def complexCTF(self):
        """Gold marker is simulated with multislice method. Afterwards it is wiener filtered and cross-correlated
        with the original gold bead under a mask."""
        from pytom.tompy.correlation import nxcc
        from pytom.simulation.microscope import create_ctf
        from pytom.tompy.tools import create_circle

        # import matplotlib
        # matplotlib.use('Qt5Agg')
        # import matplotlib.pyplot as plt

        # generate a projection of a point with multislice
        projection = self.simulateProjection()

        ctf = create_ctf(projection.shape,
                         self.param_sim['pixel_size'],
                         self.defocusU,
                         self.amplitude_contrast,
                         self.param_sim['voltage'],
                         self.Cs,
                         sigma_decay=self.sigma_decay,
                         defocusV=self.defocusV,
                         defocus_angle_deg=self.ast_angle_deg)

        ctf = np.fft.ifftshift(ctf)

        projection_corrected = np.fft.ifftn(np.conjugate(ctf) * np.fft.fftn(projection) / (ctf ** 2 + 0.1)).real

        # create an image of a point
        gold_slice = self.potential[:, :, self.z // 2 + self.gold_size[2] // 2].real

        # fig, ax = plt.subplots(1, 3)
        # ax[0].imshow(gold_slice)
        # ax[1].imshow(projection)
        # ax[2].imshow(projection_corrected)
        # plt.show()

        # add a mask for the correlation because of the noise
        mask = create_circle(projection.shape, radius=0.9*(self.gold_size[0]//2),
                             sigma=0.1*(self.gold_size[0]//2), center=[self.x//2+self.gold_size[0]//2,
                                                                  self.y//2+self.gold_size[1]//2])
        # calculate correlation
        ccc = nxcc(gold_slice, projection_corrected, mask=mask)

        # test if correlation is within the limits
        # print('cross correlation of imaged point source with a ctf = ', ccc)

        return ccc

    def test_CTFs(self):
        ccc = self.basicCTF()
        self.assertGreater(ccc, 0.7, msg='correlation not sufficient between original gold marker and wiener filtered '
                                         'gold marker')
        ccc = self.complexCTF()
        self.assertGreater(ccc, 0.8, msg='correlation not sufficient between original gold marker and wiener filtered '
                                         'simulation')


if __name__ == '__main__':
    unittest.main()
