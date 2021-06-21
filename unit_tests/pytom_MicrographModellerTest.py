"""
Test MicrographModeller tomogram simulation.
todo simulation should also have a unittest for potential generation, Fourier shell scaling, and grandmodel generation
"""
import unittest


class MicrographModellerTest(unittest.TestCase):
    def setUp(self):
        """Initialize simulation parameters"""
        from pytom.simulation.potential import iasa_integration

        self.param_pot = {
            'pdb':                  './testData/3eoj.pdb',  #todo where to save?
            'voxel_size':           2.5,
            'oversampling':         2,
            'absorption_contrast':  True,
            'voltage':              300e3
        }

        real, imag = iasa_integration(self.param_pot['pdb'],
                                      voxel_size=self.param_pot['voxel_size'],
                                      oversampling=self.param_pot['oversampling'],
                                      absorption_contrast=self.param_pot['absorption_contrast'],
                                      voltage=self.param_pot['voltage'])
        self.potential = real + 1j * imag

        # specific defocus and msdz, but otherwise default parameters for ctf function
        self.param_sim = {
            'save_path':            './testData',  #todo where to save?
            'angles':               list(range(-60, 60 + 3, 3)),
            'nodes':                1,  #todo change to multiple if possible?
            'pixel_size':           5e-10,
            'binning':              2,
            'dose':                 80,
            'voltage':              300e3,
            'defocus':              3e-6,
            'msdz':                 5e-9,
            'camera_type':          'K2SUMMIT',
            'camera_folder':        'pytom/data/detectors',
        }

        self.param_rec = {
            'save_path':            './testData',
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

    def remove_tree(self, foldername):
        """Assert folder exists, then remove its content and itself"""
        from os import path
        from shutil import rmtree

        foldercheck = path.exists(foldername)
        if not foldercheck:
            print(foldername + " does not exist")
        self.assertTrue(foldercheck, msg="file " + foldername + " does not exist")
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

    def simulateTomogram(self):
        """Run the simulation, output here will be written to some temp storage"""
        from pytom.simulation.MicrographModeller import generate_tilt_series_cpu, reconstruct_tomogram
        from pytom.tompy.io import read
        from os import path

        generate_tilt_series_cpu(self.param_sim['save_path'],
                                 self.param_sim['angles'],
                                 nodes=self.param_sim['nodes'],
                                 pixel_size=self.param_sim['pixel_size'],
                                 binning=self.param_sim['binning'],
                                 dose=self.param_sim['dose'],
                                 voltage=self.param_sim['voltage'],
                                 defocus=self.param_sim['defocus'],
                                 msdz=self.param_sim['msdz'],
                                 camera_type=self.param_sim['camera_type'],
                                 camera_folder=self.param_sim['camera_folder'])

        reconstruct_tomogram(self.param_rec['save_path'],
                             weighting=self.param_rec['weighting'],
                             reconstruction_bin=self.param_rec['reconstruction_bin'])

        return read(path.join(self.param_rec['save_path'], 'reconstruction.mrc'))

    def test_Simulation(self):
        """Run two simulations and test their correlation. Both will have a different realization of noise and will
        slightly differ."""
        from pytom.tompy.correlation import xcc

        # generate two different realization of tomogram noise
        tomo_1 = self.simulateTomogram()
        tomo_2 = self.simulateTomogram()

        # calculate cross-correlation coefficient of the two tomograms
        cc = xcc(tomo_1, tomo_2)

        # self.assertAlmostEqual(first=1.0, second=cc, places=1)
        self.assertGreater(cc, 0.95, msg='correlation is not sufficient between simulations')


if __name__ == '__main__':
    unittest.main()

