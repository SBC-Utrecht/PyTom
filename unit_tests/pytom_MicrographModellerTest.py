"""
Test MicrographModeller tomogram simulation.
"""

import unittest

class MicrographModellerTest(unittest.TestCase):

    def setUp(self):
        """Initialize simulation parameters"""

        self.pdb = './testData/3eoj.pdb'
        self.voxel_size = 2.5
        self.oversampling = 2
        self.absorption_contrast = True
        self.voltage = 300e3

        self.detector = 'K2SUMMIT'
        self.detector_data = 'pytom/data/detectors'

        self.tilt_angles = list(range(-60, 60, 3))


if __name__ == '__main__':
    unittest.main()