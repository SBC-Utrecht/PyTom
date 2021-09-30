"""
Test CTF with astigmatism and rectangular images shapes.
todo simulation should also have a unittest for potential generation, Fourier shell scaling, and grandmodel generation.

@author: Marten Chaillet
"""
import unittest
import numpy as np

class CTFTest(unittest.TestCase):
    def setUp(self):
        """Initialize simulation parameters"""
        from pytom.simulation.microscope import convert_defocusU_defocusV_to_defocus_astigmatism as convert
        import os

        self.x, self.y, self.z = 1000, 1000, 33
        self.potential = np.zeros((self.x, self.y, self.z))
        self.potential[self.x//2, self.y//2, self.z//2] = 1

        self.defocusU, self.defocusV, self.ast_angle_deg = 3e-6, 5e-6, -69
        self.defocus, self.astigmatism = convert(self.defocusU, self.defocusV)
        self.Cs = 2.7e-3
        self.amplitude_contrast = 0
        self.sigma_decay = 0.3

        # create temporary dir for storing simulation data
        if not os.path.exists('temp_simulation'):
            os.mkdir('temp_simulation')

        # specific defocus and msdz, but otherwise default parameters for ctf function
        self.param_sim = {
            'save_path':            './temp_simulation',
            'angles':               [0],
            'nodes':                1,  # todo change to multiple if possible ??
            'pixel_size':           1e-10,
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

        # TODO remove redundant ctf functions once tested
        # TODO fix astigmatism passing to simulation

        # self.param_rec = {
        #     'save_path':            './temp_simulation',
        #     'weighting':            1,
        #     'reconstruction_bin':   1
        # }

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
                                 msdz=self.param_sim['msdz'],
                                 camera_type=self.param_sim['camera_type'],
                                 camera_folder=self.param_sim['camera_folder'],
                                 grandcell=self.potential)

        return read(os.path.join(self.param_sim['save_path'] + c, 'noisefree_projections.mrc')).squeeze()

    def test_CTF(self):
        """Run a multislice simulation of a point source. Compare its CTF with a ctf from the create_ctf() function."""
        from pytom.tompy.correlation import nxcc
        from pytom.simulation.microscope import create_ctf

        import matplotlib
        matplotlib.use('Qt5Agg')
        import matplotlib.pyplot as plt

        # generate a projection of a point with multislice
        projection = self.simulateProjection()

        # create an image of a point
        point = self.potential[:, :, self.z//2]

        point_ft = np.fft.fftn(point)

        # get a ctf
        ctf = create_ctf(projection.shape,
                         self.param_sim['pixel_size'] * self.param_sim['binning'],
                         self.defocusU,
                         self.amplitude_contrast,
                         self.param_sim['voltage'],
                         self.Cs,
                         sigma_decay=self.sigma_decay,
                         defocusV=self.defocusV,
                         defocus_angle_deg=self.ast_angle_deg)

        ctf_abs = np.abs(ctf)

        point_conv = np.fft.fftn(point_ft * np.fft.ifftshift(ctf)).real
        point_power_spectrum = np.abs(np.fft.fftshift(np.fft.fftn(point_conv)))

        # fig, ax = plt.subplots(1, 3)
        # ax[0].imshow(point_conv)
        # ax[1].imshow(point_power_spectrum)
        # ax[2].imshow(ctf_abs)
        # plt.show()

        ccc1 = nxcc(point_power_spectrum, ctf_abs)

        # test if correlation is within the limits
        print('cross correlation of imaged point source with a ctf = ', ccc1)
        self.assertGreater(ccc1, 0.99, msg='correlation not sufficient between point convoluted with ctf and ctf')

        projection_power_spectrum = np.abs(np.fft.fftshift(np.fft.fftn(projection)))

        fig, ax = plt.subplots(1, 3)
        ax[0].imshow(projection)
        ax[1].imshow(projection_power_spectrum)
        ax[2].imshow(ctf_abs)
        plt.show()

        # calculate correlation
        ccc2 = nxcc(projection_power_spectrum, ctf_abs)

        # test if correlation is within the limits
        print('cross correlation of imaged point source with a ctf = ', ccc2)
        # self.assertGreater(ccc2, 0.80, msg='correlation not sufficient between ctf of multislice imaged point and ctf')


if __name__ == '__main__':
    unittest.main()
