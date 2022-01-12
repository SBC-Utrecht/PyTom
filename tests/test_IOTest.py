"""
Created on Jul 30, 2012

@author: hrabe
"""

import unittest
import numpy as np
import os


class pytom_IOTest(unittest.TestCase):
    def setUp(self):
        self.fnames = []
        self.epsilon = 1E-4
        self.outfolder = 'IOTest'
        self.rotation_angles = [10.5, 23.5, -32.5]
        self.sx, self.sy, self.sz=10,11,12
        self.tilt_angle = 132.23
        self.pixel_size = 2.62
        if not os.path.exists(self.outfolder):
            os.mkdir(self.outfolder)

    def tearDown(self):
        from tests.helper_functions import remove_tree
        remove_tree(self.outfolder)

    def read(self):
        """
        test general read function
        """
        from pytom.basic.files import read

        v = read(f'./testData/ribo.em')
        self.assertTrue( v.sizeX() == 100 and v.sizeY() == 100 and v.sizeZ() == 100)
        
        v = read(f'./testData/ribo.em',binning=[4,4,4])
        assert v.sizeX() == 25 and v.sizeY() == 25 and v.sizeZ() == 25
    
    def readem(self):
        """
        test readem
        """
        from pytom.basic.files import read_em

        int2_file = f'./testData/int2.em'
        (image, header) = read_em(int2_file)
        inibytes = header.get_1st4bytes()
        self.assertTrue(inibytes[3] == 5,
                        'Datatype number is wrong! should have been changed to 5 from 2 in initial file')
        dims = header.get_dim()

        #now binning
        (image, header) = read_em(int2_file, binning=[4, 4, 1])
        bdims = header.get_dim()
        self.assertTrue(bdims[0]*4 == dims[0],
                        'Dimensions of binned file are wrong!')

    def write_read_MRC(self):
        self.data_write_read('mrc')

    def write_read_EM(self):
        from pytom.agnostic.io import write, read

        fname = f'{self.outfolder}/dummy_reading.em'
        self.fnames.append(fname)
        data = np.ones((self.sx, self.sy, self.sz))

        write(fname, data)
        self.assertTrue(os.path.exists(fname))

        data2 = read(fname)
        self.assertTrue((data != data2).sum() == 0)

    def write_read_STAR(self):
        from pytom.agnostic.io import write, read
        from pytom.basic.datatypes import DATATYPE_RELION31_EXTENDED, fmtR31EXTENDED, headerRelion31EXTENDEDSubtomo

        fname = f'testData/example_starfile.star'

        data2 = read(fname, dtype=DATATYPE_RELION31_EXTENDED)
        oname = f'{self.outfolder}/dummy_star.star'
        self.fnames.append(oname)
        write(oname, data2, fmt=fmtR31EXTENDED, header=headerRelion31EXTENDEDSubtomo)
        self.assertTrue(os.path.exists(oname))

        data2 = read(fname)
        oname = f'{self.outfolder}/dummy_star2.star'
        self.fnames.append(oname)
        write(oname, data2)
        self.assertTrue(os.path.exists(oname))

        fname = 'testData/relion_with_optics_group.star'
        data3 = read(fname)
        data4 = read(fname, read_optics_group=True)

        self.assertTrue(len(data3) == 60,
                        f'Error reading data from {fname}: epected {78} lines, found {len(data3)} lines')
        self.assertTrue(len([data4]) == 1,
                        f'Error reading data from {fname}: epected { 1} lines, found {len([data4])} lines')

    def write_read_REC(self):
        self.data_write_read('rec')

    def write_read_ST(self):
        self.data_write_read('st')

    def read_LOG(self):
        from pytom.agnostic.io import read
        fname = f'testData/taSolution.log'
        a = read(fname)

    def write_read_TXT(self):
        from pytom.agnostic.io import write, read
        from pytom.basic.datatypes import DATATYPE_METAFILE, FMT_METAFILE, HEADER_METAFILE

        fname = f'{self.outfolder}/data.meta'

        data = np.zeros((10), dtype=DATATYPE_METAFILE)
        data['TiltAngle'] = range(1, 11)

        for i in range(10):
            data['FileName'][i] = f'sorted_{i:02d}.mrc'

        write(fname, data, fmt=FMT_METAFILE, header=HEADER_METAFILE)
        self.assertTrue(os.path.exists(fname), f'writing {fname} failed')

        data = read(fname, dtype=DATATYPE_METAFILE)
        self.assertTrue(data['TiltAngle'][9] == 10, f'reading {fname} failed')

    def read_size(self):
        from pytom.agnostic.io import read_size

        for d, sized in ( ('x', self.sx), ('y', self.sy), ('z', self.sz)):
            x = read_size(self.fnames[0], d)
            self.assertTrue(x==sized, f'size in {d}-dimension is off: found {x}, expected {sized}')

        x,y,z = read_size(self.fnames[0])
        self.assertTrue(x==self.sx and y==self.sy and z==self.sz,
                        f'Reading sizes failed:\n\texpected {self.sx}.{self.sy},{self.sz}\n\tfound {x},{y},{z}')

    def read_tilt_angle(self):
        from pytom.agnostic.io import read_tilt_angle

        a = read_tilt_angle(self.fnames[1])
        self.assertTrue(np.abs(a-self.tilt_angle) < self.epsilon,
                        f'Tilt Angle is not correct: found {a}, expected {self.tilt_angle}')

    def write_tilt_angle(self):
        from pytom.agnostic.io import write_tilt_angle, read_tilt_angle

        fname = self.fnames[1]
        self.rotation_angles[0] = ta = -23.123
        write_tilt_angle(fname, self.rotation_angles[0])

        a = read_tilt_angle(fname)
        self.assertTrue(np.abs(a-ta) < self.epsilon,
                        f'Tilt Angle is not correct: found {a}, expected {ta}')

    def read_header(self):
        from pytom.agnostic.io import read_header
        header = read_header(self.fnames[1])
        assert header

    def write_rotation_angle(self):
        from pytom.agnostic.io import write_rotation_angles
        fname = self.fnames[1]
        tas = [-23, 234.12, 121.22]
        write_rotation_angles(fname,z1=tas[0],x=tas[1],z2=tas[2])

        self.read_rotation_angle(tas=tas)

    def read_rotation_angle(self, tas=None):
        from pytom.agnostic.io import read_rotation_angles

        tas = self.rotation_angles if tas is None else tas

        a = read_rotation_angles(self.fnames[1])
        for n, angle in enumerate(a):
            self.assertTrue(np.abs(angle - tas[n]) < self.epsilon,
                        f'Rotation Angle {n} is not correct: found {angle}, expected {tas[n]}')

    def read_pixelsize(self):
        from pytom.agnostic.io import read_pixelsize

        fname = self.fnames[1]

        a = read_pixelsize(fname)[0]
        self.assertTrue(np.abs(a - self.pixel_size) < self.epsilon,
                    f'Tilt Angle is not correct: found {a}, expected {self.pixel_size}')

    def data_write_read(self, extension='mrc'):
        from pytom.agnostic.io import write, read

        fname = f'{self.outfolder}/dummy_reading.{extension.lower()}'
        self.fnames.append(fname)
        data = np.random.random((self.sx,self.sy,self.sz)).astype(np.float32)

        write(fname, data, pixel_size=self.pixel_size, tilt_angle=self.tilt_angle)
        self.assertTrue(os.path.exists(fname), f'writing {fname} failed')

        data2 = read(fname)

        write(f'{self.outfolder}/testing.mrc', data2 -data)
        write(f'{self.outfolder}/testing2.mrc', data2 )

        self.assertTrue(np.abs(data - data2).sum() < self.epsilon, f'reading {fname} failed: input and output files are different')

    def runTest(self):
        self.read()
        self.readem()
        self.write_read_EM()
        self.write_read_MRC()
        self.write_read_REC()
        self.write_read_ST()
        self.write_read_TXT()
        self.read_LOG()
        self.write_read_STAR()
        self.read_size()
        self.read_tilt_angle()
        self.write_tilt_angle()
        self.write_rotation_angle()
        self.read_header()
        self.read_pixelsize()

        
if __name__ == '__main__':
    unittest.main()
