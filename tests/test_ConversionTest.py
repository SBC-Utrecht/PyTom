import unittest
import numpy as np
import os
import pathlib
from pytom.voltools.utils import transform_matrix
from pytom.agnostic.io import read
from pytom.basic.files import pl2star


# test data folder and test files
test_data = pathlib.Path('testData')
subtomo_star = test_data.joinpath('subtomo_format.star')
pickpos_xml = test_data.joinpath('pickpos_format.xml')
# output files from conversion should be named
subtomo_xml = 'subtomo_format.xml'
subtomo_copy_star = 'subtomo_format_copy.star'
pickpos_star = 'pickpos_format.star'


class pytom_ConversionTest(unittest.TestCase):
    @classmethod
    def tearDownClass(cls):
        test_data.joinpath(subtomo_xml).unlink()
        test_data.joinpath(subtomo_copy_star).unlink()
        test_data.joinpath(pickpos_star).unlink()

    def imod_forward(self):
        # imod (.xf) to pytom (alignmentResults.txt)
        shift_x, shift_y = np.dot(self.xf_matrix, np.array([self.xf_shift_x, self.xf_shift_y]))
        scale = np.sqrt(self.xf_matrix[0, 0] ** 2 + self.xf_matrix[1, 0] ** 2)
        in_plane_rotation = np.rad2deg(np.arctan2(self.xf_matrix[1, 0], self.xf_matrix[0, 0]))
        self.assertAlmostEqual(self.shift_x, shift_x, 2)
        self.assertAlmostEqual(self.shift_y, shift_y, 2)
        self.assertAlmostEqual(self.scale, scale, 2)
        self.assertAlmostEqual(self.in_plane_rotation, in_plane_rotation, 2)

    def imod_backward(self):
        # pytom (alignmentResults.txt) to imod (.xf)
        m = transform_matrix(rotation=(self.in_plane_rotation, 0, 0),
                             rotation_order='rzxz',
                             scale=(self.scale,) * 3)[:2, :2]
        # shift_x, shift_y = np.dot(m, np.array([self.shift_x, self.shift_y]))
        # self.assertAlmostEqual(self.xf_shift_x, shift_x, 2)
        # self.assertAlmostEqual(self.xf_shift_y, shift_y, 2)
        self.assertAlmostEqual(np.abs(self.xf_matrix - m.T).sum(), 0, 2)

    def test_imod_parameter_conversion(self):
        # params from xf file
        self.xf_matrix = np.array([[0.9990135, -0.0829354], [0.0829354, 0.9990134]])
        self.xf_shift_x = 22.127
        self.xf_shift_y = -24.819

        # params pytom
        self.scale = 1.0024501
        self.in_plane_rotation = 4.7456585
        self.shift_x = 24.16354487
        self.shift_y = -22.95940421

        self.imod_forward()
        self.imod_backward()

    def test_subtomo_format_conversion(self):
        # do forward and backward conversion of the subtomo star format
        os.system(f'convert.py -f {subtomo_star} -t testData --outname subtomo_format -o xml')
        os.system(f'convert.py -f {test_data.joinpath(subtomo_xml)} -t testData --outname subtomo_format_copy '
                  f'-o star --rln-voltage 200 --rln-spherical-aberration 2.7 --pixelSize 6.896')

        # read them in and test whether the copy has still all the same header as the original
        star_original = read(subtomo_star)
        star_copy = read(test_data.joinpath(subtomo_copy_star))
        original_in_copy = [col in star_copy.dtype.names for col in star_original.dtype.names]
        self.assertTrue(all(original_in_copy), msg='after back and forth conversion not all header are present')

        # check whether the function raises the right exception when voltage or Cs is missing
        with self.assertRaises(ValueError):
            pl2star(test_data.joinpath(subtomo_xml), test_data)
        with self.assertRaises(ValueError):
            pl2star(test_data.joinpath(subtomo_xml), test_data, rln_spherical_aberration=2.7)
        with self.assertRaises(ValueError):
            pl2star(test_data.joinpath(subtomo_xml), test_data, rln_voltage=200)

    def test_pickpos_conversion(self):
        # convert forward
        os.system(f'convert.py -f {pickpos_xml} -t testData --outname pickpos_format -o star '
                  f'--pixelSize 20')

        # test if all the headers are present for the pickpos format that can be read by Warp
        star_pickpos = read(test_data.joinpath(pickpos_star))
        required_headers = ['CoordinateX', 'CoordinateY', 'CoordinateZ',
                            'MicrographName','Magnification', 'DetectorPixelSize',
                            'GroupNumber', 'AngleRot', 'AngleTilt', 'AnglePsi']
        headers_present = [col in star_pickpos.dtype.names for col in required_headers]
        self.assertTrue(all(headers_present), msg='pick position star file for Warp does not have all the required '
                                                  'data')


if __name__ == '__main__':
    unittest.main()
