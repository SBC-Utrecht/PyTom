import os
import pathlib
import shutil
import numpy as np
import unittest
from pytom.agnostic.correlation import flcf, nxcc
from pytom.agnostic.io import read
from pytom.agnostic.tools import create_circle
from pytom.basic.datatypes import DATATYPE_ALIGNMENT_RESULTS_RO, fmtAlignmentResultsRo as fmt, \
    HEADER_ALIGNMENT_RESULTS_RO as header
from pytom.basic.files import loadtxt, savetxt
from pytom.basic.transformations import general_transform2d
from pytom.lib.pytom_volume import read as readC
from pytom.lib.pytom_numpy import vol2npy
from pytom.gpu.initialize import xp

only_run_cpu = False
try:
    from pytom.gpu.initialize import initialize_gpu
    initialize_gpu(0)
except ImportError:
    only_run_cpu = True
if only_run_cpu:
    from numpy.fft import fftshift, fftn, ifftn
else:
    from cupy.fft import fftshift, fftn, ifftn

imod_reconstruction_folder = pathlib.Path('ImodReconstruction')
imod_aligned_images_folder = pathlib.Path('ImodAlignedImages')
raw_tilt_images_folder = pathlib.Path('RawTiltImages')


def submit(cmd):
    os.system(cmd)


def gaussian_filter(data, high, sigma):
    a = create_circle(data.shape, high, sigma)
    return ifftn(fftn(data) * fftshift(a)).real


class pytom_AlignmentReconstructionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not imod_reconstruction_folder.exists():
            os.system('tar -zxvf ../testData/imod_reconstruction.tar.gz -C ./ ')
        imod_reconstruction_folder.mkdir(exist_ok=True)
        imod_aligned_images_folder.mkdir(exist_ok=True)
        raw_tilt_images_folder.mkdir(exist_ok=True)

        # convert imod files for pytom reading
        submit(f'mrcs2mrc.py -f {imod_reconstruction_folder.joinpath("tomo1_bin6.mrc.st")} '
               f'-t {raw_tilt_images_folder} --prefix sorted ')
        # Import IMOD ali parameters into alignmentResults.txt
        cmd = f'convert.py -f {imod_reconstruction_folder.joinpath("tomo1_bin6.mrc.xf")} -o txt ' \
              f'--tlt-file {imod_reconstruction_folder.joinpath("tomo1.mrc.tlt")} -t {raw_tilt_images_folder} ' \
              f'--sortedFolder {raw_tilt_images_folder.absolute()} '
        submit(cmd)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(imod_reconstruction_folder)
        shutil.rmtree(imod_aligned_images_folder)
        shutil.rmtree(raw_tilt_images_folder)

    @staticmethod
    def run_imod():
        # created aligned tilt stack, tomogram, and unpack .ali stack to alignment folder
        submit(f'cd {imod_reconstruction_folder}; submfg newst.com ')
        submit(f'cd {imod_reconstruction_folder}; submfg newst_bin.com ')
        submit(f'cd {imod_reconstruction_folder}; submfg tilt.com ')
        submit(f'cd {imod_reconstruction_folder}; trimvol -yz tomo1_bin12.mrc.rec tomogram_IMOD.mrc ')
        submit(f'mrcs2mrc.py -f {imod_reconstruction_folder.joinpath("tomo1_bin12.mrc.ali")} '
               f'-t {imod_aligned_images_folder} --prefix sorted ')
        submit(f'mrcs2mrc.py -f {imod_reconstruction_folder.joinpath("tomo1_bin6.mrc.ali")} '
               f'-t {imod_aligned_images_folder} --prefix unbinned_sorted ')

    @staticmethod
    def create_alignment_file():
        # get angles from taSolution.log file
        tilt_angles = np.loadtxt(imod_reconstruction_folder.joinpath('tomo1.mrc.tlt'))

        # create dummy alignmentResults.txt file: scaling =1, other ali params = 0, filename need update
        alignment_data = np.zeros((len(tilt_angles)), dtype=DATATYPE_ALIGNMENT_RESULTS_RO)
        alignment_data['TiltAngle'] = tilt_angles
        alignment_data['Magnification'] = 1
        alignment_data['OperationOrder'] = 'TSR'

        # iterdir() does not sort the files, so need to explicitly pass file list to sorted()
        for n, file_name in enumerate(sorted([f for f in imod_aligned_images_folder.iterdir() if
                                              f.name.startswith('sorted')])):
            alignment_data['FileName'][n] = file_name.name

        savetxt(imod_aligned_images_folder.joinpath('alignmentResults.txt'), alignment_data, fmt=fmt, header=header)

    def setUp(self):
        # Image 33 is chosen as it has the highest shift and scaling of all images in the tilt series
        self.imageID = 33  # 25 or 33

        # WBP command str. Needs formatting in tests
        self.WBPFormat = """cd {}
                reconstructWB.py \\
                    --alignResultFile {} \\
                    --tomogram {} \\
                    --projBinning {} \\
                    --applyWeighting -1  \\
                    --size 300,300,150
                """

        self.WBPFormat2 = """cd {}
                reconstructWB.py \\
                    --alignResultFile {} \\
                    --tomogram {} \\
                    --projBinning {} \\
                    --applyWeighting -1  \\
                    --size 300,300,150 \\
                    {}
                """

        # we need the pytom cpu reconstruction for each test, but we dont want to reconstruct it multiple times
        # via setUp() => only reconstruct if it does not yet exist
        outfile = 'tomogram_PYTOM.mrc'
        if not raw_tilt_images_folder.joinpath(outfile).exists():
            submit(self.WBPFormat.format(raw_tilt_images_folder, 'alignmentResults.txt', outfile, 2))

    def test_pytom_vs_imod_reconstruction(self):
        if shutil.which('submfg') is None:
            raise unittest.SkipTest('Cannot run unittest. IMOD not installed.')

        self.run_imod()
        self.create_alignment_file()

        self.check_image_alignment()

        # check the shift between the two
        # expected shift is -1 in x
        self.check_tomogram_shift(imod_reconstruction_folder.joinpath('tomogram_IMOD.mrc'),
                                  raw_tilt_images_folder.joinpath('tomogram_PYTOM.mrc'), expected=[0, 0, 0])

    def check_tomogram_shift(self, fname1, fname2, expected=[0, 0, 0]):
        """
        @expected: is the expected shift between the two tomograms
        """
        a = read(fname1)
        b = read(fname2)

        c = flcf(a, b)

        loc = np.unravel_index(c.argmax(), c.shape)
        shift = (loc[0]-b.shape[0]//2, loc[1]-b.shape[1]//2, loc[2]-b.shape[2]//2)

        self.assertTrue(shift[0] == expected[0] and shift[1] == expected[1] and shift[2] == expected[2],
                        f'shifts between tomograms not zero, namely: {shift[0]} {shift[1]} {shift[2]}')

    def check_image_alignment(self):
        converted_alignment = loadtxt(raw_tilt_images_folder.joinpath('alignmentResults.txt'),
                                      dtype=DATATYPE_ALIGNMENT_RESULTS_RO)

        data = readC(str(converted_alignment['FileName'][self.imageID]))
        indata = readC(str(imod_aligned_images_folder.joinpath(f'unbinned_sorted_{self.imageID:02d}.mrc')))

        offset = 100
        indata = vol2npy(indata).copy()
        # force cast to cupy if needed
        indata = xp.array(indata)
        indata = indata[offset:-offset, offset:-offset]
        indata = gaussian_filter(indata, 30, 10)
        indata = (indata - indata.mean()) / indata.std()

        # get the alignment parameters
        # order = str2list(a['OperationOrder'][0])
        # assert order == [1, 0, 2], 'operation order not as expected from an imod alignment'
        shift = [converted_alignment['AlignmentTransX'][self.imageID],
                 converted_alignment['AlignmentTransY'][self.imageID]]
        mag = converted_alignment['Magnification'][self.imageID]
        rot = converted_alignment['InPlaneRotation'][self.imageID]

        """
        imod docs say the following:
        from coordinate Xi and Yi in the input image to Xo and Yo in output.
        Xc and Yc is center of input image. nxb nyb is the shape of the output image.
        xt and yt are the translations. And a represents the 2x2 matrix for rotation, scaling, and skewing.
        Xo = a11(Xi - Xc) + a12(Yi - Yc) + nxb/2. + xt  
        Yo = a21(Xi - Xc) + a22(Yi - Yc) + nyb/2. + yt
        """
        # imod center is (s - 1) / 2
        center = ((data.sizeX() - 1) / 2., (data.sizeY() - 1) / 2., 0)
        out = general_transform2d(data, shift=shift, scale=float(mag), rot=float(rot), order=[1, 0, 2], crop=True,
                                  center=center)
        res = vol2npy(out).copy()
        # Force to cupy if needed
        res = xp.array(res)

        # crop volumes
        res = res[offset:-offset, offset:-offset]

        # apply gaussian to pytom image
        res = gaussian_filter(res, 30, 10)
        res = (res - res.mean()) / res.std()

        # calculate xcorr between the two
        ccc = nxcc(res, indata)

        self.assertTrue(ccc > 0.99, f'cross-correlation not good enough: {ccc:.3f}')

    def test_alignment_reconstruction(self):
        # make sure that a tomogram was reconstructed with the proper dimensions
        tomo = read(raw_tilt_images_folder.joinpath('tomogram_PYTOM.mrc'))
        self.assertTrue(tomo.shape == (300, 300, 150), msg='reconstructed tomogram does not have expected dimensions')

    @unittest.skipIf(only_run_cpu,
                     "The test below uses a GPU and cannot be executed in a docker environment")
    def test_compare_cpu_gpu_reconstruction(self):
        # always run on gpu 0
        outfile_gpu = 'tomogram_PYTOM_gpu.mrc'
        submit(self.WBPFormat2.format(raw_tilt_images_folder, 'alignmentResults.txt', outfile_gpu, 2, '-g 0'))
        self.check_tomogram_shift(raw_tilt_images_folder.joinpath('tomogram_PYTOM.mrc'),
                                  raw_tilt_images_folder.joinpath(outfile_gpu), expected=[0, 0, 0])


if __name__ == '__main__':
    unittest.main()
