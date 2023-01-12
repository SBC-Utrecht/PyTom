import unittest
import numpy as np
from pytom.voltools.utils import transform_matrix


class pytom_ConversionTest(unittest.TestCase):
    def setUp(self):
        # params from xf file
        self.xf_matrix = np.array([[0.9990135, -0.0829354], [0.0829354,  0.9990134]])
        self.xf_shift_x = 22.127
        self.xf_shift_y = -24.819

        # params pytom
        self.scale = 1.0024501
        self.in_plane_rotation = 4.7456585
        self.shift_x = 24.16354487
        self.shift_y = -22.95940421

    def forward(self):
        # imod (.xf) to pytom (alignmentResults.txt)
        shift_x, shift_y = np.dot(self.xf_matrix, np.array([self.xf_shift_x, self.xf_shift_y]))
        scale = np.sqrt(self.xf_matrix[0, 0] ** 2 + self.xf_matrix[1, 0] ** 2)
        in_plane_rotation = np.rad2deg(np.arctan2(self.xf_matrix[1, 0], self.xf_matrix[0, 0]))
        self.assertAlmostEqual(self.shift_x, shift_x, 2)
        self.assertAlmostEqual(self.shift_y, shift_y, 2)
        self.assertAlmostEqual(self.scale, scale, 2)
        self.assertAlmostEqual(self.in_plane_rotation, in_plane_rotation, 2)

    def backward(self):
        # pytom (alignmentResults.txt) to imod (.xf)
        m = transform_matrix(rotation=(self.in_plane_rotation, 0, 0),
                             rotation_order='rzxz',
                             scale=(self.scale,) * 3)[:2, :2]
        # shift_x, shift_y = np.dot(m, np.array([self.shift_x, self.shift_y]))
        # self.assertAlmostEqual(self.xf_shift_x, shift_x, 2)
        # self.assertAlmostEqual(self.xf_shift_y, shift_y, 2)
        self.assertAlmostEqual(np.abs(self.xf_matrix - m.T).sum(), 0, 2)

    def runTest(self):
        self.forward()
        self.backward()


if __name__ == '__main__':
    unittest.main()
