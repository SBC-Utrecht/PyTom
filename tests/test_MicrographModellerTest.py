"""
Test MicrographModeller tomogram simulation.
todo simulation should also have a unittest for potential generation, Fourier shell scaling, and grandmodel generation.

@author: Marten Chaillet
"""
import unittest
import os
import numpy as np

# TODO: this shouldn't happen, we shouldn't force interactive backend outside of pytomGUI
@unittest.skipIf(os.environ.get('AM_I_IN_A_DOCKER_CONTAINER', False),
                 "The tests below call matplotlib and force to use qt5agg which is unavailable on default docker")
class MicrographModellerTest(unittest.TestCase):
    def setUp(self):
        """Initialize simulation parameters"""
        from pytom.simulation.potential import iasa_integration
        import os

        self.param_pot = {
            'pdb':                  './testData/3j9m.cif',
            'voxel_size':           5,
            'oversampling':         2,
            'solvent_masking':      True,
            'absorption_contrast':  True,
            'voltage':              300e3
        }

        self.potential = iasa_integration(self.param_pot['pdb'],
                                          voxel_size=self.param_pot['voxel_size'],
                                          oversampling=self.param_pot['oversampling'],
                                          solvent_masking=self.param_pot['solvent_masking'],
                                          absorption_contrast=self.param_pot['absorption_contrast'],
                                          voltage=self.param_pot['voltage'])

        if self.potential.shape[0] % 2:
            self.potential = np.pad(self.potential, pad_width=(0, 1), mode='constant', constant_values=0)

        # create temporary dir for storing simulation data
        if not os.path.exists('temp_simulation'):
            os.mkdir('temp_simulation')

        # specific defocus and msdz, but otherwise default parameters for ctf function
        self.param_sim = {
            'save_path':            './temp_simulation',
            'angles':               list(range(-60, 60 + 3, 3)),
            'nodes':                1,  # todo change to multiple if possible ??
            'pixel_size':           5e-10,
            'binning':              2,
            'dose':                 80,
            'voltage':              300e3,
            'defocus':              3e-6,
            'msdz':                 5e-9,
            'camera_type':          'K2SUMMIT',
            'camera_folder':        '../pytom/simulation/detectors',
        }

        self.param_rec = {
            'save_path':            './temp_simulation',
            'weighting':            1,
            'reconstruction_bin':   1
        }

    def tearDown(self):
        """Remove all the files gathered during simulation"""
        from os import path
        dir = self.param_sim['save_path']
        self.remove_tree(path.join(dir, 'projections'))
        self.remove_file(path.join(dir, 'projections.mrc'))
        self.remove_file(path.join(dir, 'noisefree_projections.mrc'))
        self.remove_file(path.join(dir, 'alignment_simulated.txt'))
        self.remove_file(path.join(dir, 'reconstruction.mrc'))
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

    def simulateTomogram(self, c=''):
        """Run the simulation, output here will be written to some temp storage"""
        from pytom.simulation.MicrographModeller import generate_tilt_series_cpu, reconstruct_tomogram
        from pytom.agnostic.io import read
        from os import path
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

        reconstruct_tomogram(self.param_rec['save_path'] + c,
                             weighting=self.param_rec['weighting'],
                             reconstruction_bin=self.param_rec['reconstruction_bin'])

        return read(path.join(self.param_rec['save_path'] + c, 'reconstruction.mrc'))

    def test_Simulation(self):
        """Run two simulations and test their correlation. Both will have a different realization of noise and will
        slightly differ."""
        from pytom.agnostic.correlation import nxcc
        from pytom.agnostic.tools import create_sphere
        from pytom.simulation.support import reduce_resolution

        # generate two different realization of tomogram noise
        tomo_1 = self.simulateTomogram()
        tomo_2 = self.simulateTomogram()

        # mask for correlation
        r = int(tomo_1.shape[0] / 2 * 0.8)
        mask = create_sphere(tomo_1.shape, radius=r, sigma=r / 20., num_sigma=2)

        # calculate cross-correlation coefficient of the two tomograms
        spacing = self.param_sim['pixel_size'] * 1e10
        cc = nxcc(reduce_resolution(tomo_1, spacing, 2 * spacing * 8),
                  reduce_resolution(tomo_2, spacing, 2 * spacing * 8),
                  mask=mask)

        print('normalized cross correlation of two simulations of identical volume after binning both subtomograms 8 '
              'times = ', cc)
        self.assertGreater(cc, 0.75, msg='correlation is not sufficient between simulations')


if __name__ == '__main__':
    unittest.main()

