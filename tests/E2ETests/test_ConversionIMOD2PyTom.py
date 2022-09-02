import os
import numpy as np
import unittest
from pytom.basic.datatypes import DATATYPE_ALIGNMENT_RESULTS_RO, fmtAlignmentResultsRo as fmt, \
    HEADER_ALIGNMENT_RESULTS_RO as header, DATATYPE_TASOLUTION as dtype_ta
from pytom.basic.files import loadtxt, savetxt
from pytom.agnostic.correlation import FLCF, nxcc
from pytom.agnostic.io import read, write
from pytom.agnostic.tools import create_circle, convert_operation_order_list2str as list2str, \
    convert_operation_order_str2list as str2list
from numpy.fft import fftshift, fftn, ifftn
from pytom_volume import read as readC, vol, pasteCenter
from pytom_numpy import vol2npy
from pytom.basic.transformations import general_transform2d, resize


class pytom_MyFunctionTest(unittest.TestCase):
    def create_folder(self, folder):
        if not os.path.exists(folder): os.mkdir(folder)

    def submit(self, cmd):
        os.system(cmd)

    def setUp(self):
        # tomogram id
        self.id = 1

        # Image 33 is chosen as it has the highest shift and scaling of all images in the tilt series
        self.imageID = 33  # 25 or 33

        func = self.submit

        # Clone data repo if not existing
        if not os.path.exists('PyTomUnitTests'):
            gitcmd = 'git clone https://github.com/McHaillet/PyTomUnitTests.git'
            func(gitcmd)

        # main folders
        self.datadir = 'PyTomUnitTests/ImodReconstruction'
        self.alidir = 'PyTomUnitTests/ImodAlignedImages'
        self.rawdir = 'PyTomUnitTests/RawTiltImages'

        for f in (self.rawdir, self.alidir, self.datadir):
            self.create_folder(f)

        # Check if possible to run IMOD
        if len(os.popen('which submfg').read()) < 5:
            raise Exception('Cannot run unittest. IMOD not installed.')

        # created aligned tilt stack, tomogram, and unpack ali stack to Ali folder
        func(f'cd {self.datadir}; submfg newst.com ')
        func(f'cd {self.datadir}; submfg newst_bin1.com ')
        func(f'cd {self.datadir}; submfg tilt.com ')
        func(f'cd {self.datadir}; trimvol -yz tomo{self.id}.mrc_full.rec tomogram_IMOD.mrc ')
        func(f'mrcs2mrc.py -f {self.datadir}/tomo{self.id}.mrc.ali -t {self.alidir} --prefix sorted ')
        func(f'mrcs2mrc.py -f {self.datadir}/tomo{self.id}.mrc.bin1.ali -t {self.alidir} --prefix unbinned_sorted ')
        func(f'mrcs2mrc.py -f {self.datadir}/tomo{self.id}.mrc.st -t {self.rawdir} --prefix sorted ')

        # Markerfile is not used but must exist
        func(f'touch {self.rawdir}/markerfile.txt')
        func(f'touch {self.alidir}/markerfile.txt')

        # get angles from ta file
        ta = loadtxt(os.path.join(self.datadir, 'taSolution.log'), dtype=dtype_ta, skip_header=3)
        angles = ta['Tilt']

        # create dummy alignmentResults.txt file: scaling =1, other ali params = 0, filename need update
        a = np.zeros((35), dtype=DATATYPE_ALIGNMENT_RESULTS_RO)
        a['TiltAngle'] = angles
        a['Magnification'] = 1
        a['OperationOrder'] = 'TSR'  # 'TSR'?

        for n, fname in enumerate([fname for fname in os.listdir(self.alidir) if fname.startswith('sorted')]):
            a['FileName'][n] = fname

        savetxt(f'{self.alidir}/alignmentResults.txt', a, fmt=fmt, header=header)

        # Import IMOD ali parameters into alignmentResults.txt
        cmd = f'convert.py -f {self.datadir}/taSolution.log -o txt ' \
              f'--alignxf {self.datadir}/tomo{self.id}.mrc.xf -t {self.rawdir} ' \
              f'--sortedFolder {os.getcwd()}/{self.rawdir}'
        func(cmd)

        # WBP command str. Needs formatting in tests
        self.WBPFormat = """cd {}
        reconstructWB.py \\
            --alignResultFile {} \\
            --tomogram {} \\
            --projBinning {} \\
            --applyWeighting -1  \\
            --size 464,480,464
        """

        self.WBPFormat2 = """cd {}
        reconstructWB.py \\
            --alignResultFile {} \\
            --tomogram {} \\
            --projBinning {} \\
            --applyWeighting -1  \\
            --size 464,480,464 \\
            --projectionDirectory {} \\
            --projectionPrefix {}
        """

    def checkReconfromIMODali(self):
        # make a reconstruction from the imod aligned and prebinned tilts
        # the alignmentresults only contain the tilt angles
        outfile = 'tomogram_IMOD_ALI.mrc'
        cmd = self.WBPFormat.format(self.alidir, 'alignmentResults.txt', outfile, 1)
        self.submit(cmd)

        # check the shift between the two
        self.check_shift(f'{self.datadir}/tomogram_IMOD.mrc', f'{self.alidir}/{outfile}')

    def check_shift(self, fname1, fname2, exp=[0,0,0]):
        """
        @exp: is the expexted shift between the two tomos
        """
        a = read(fname1)
        b = read(fname2)

        c = FLCF(a, b)

        loc = np.unravel_index(c.argmax(), c.shape)
        shift = (loc[0]-b.shape[0]//2, loc[1]-b.shape[1]//2, loc[2]-b.shape[2]//2)

        self.assertTrue(shift[0] == exp[0] and shift[1] == exp[1] and shift[2] == exp[2],
                        f'shifts between tomograms not zero, namely: {shift[0]} {shift[1]} {shift[2]}')

    def checkAlignmentUnbinned(self, plot=False):
        from pytom.agnostic.transform import resize

        fname  = f'{self.rawdir}/alignmentResults.txt'
        a = loadtxt(fname, dtype=DATATYPE_ALIGNMENT_RESULTS_RO)

        # square = True
        # before = False
        offset = 500

        r = self.imageID
        fname = a['FileName'][r]
        raw = readC(fname)
        data = vol(3838, 3838, 1)
        pasteCenter(raw, data)

        mm = f'{self.alidir}/unbinned_sorted_{r:02d}.mrc'
        raw = readC(mm)
        indata = vol(3838, 3838, 1)
        pasteCenter(raw, indata)  # square the imod image
        indata = vol2npy(indata).copy()

        # get the alignment parameters
        order = str2list(a['OperationOrder'][0])
        assert order == [1, 0, 2], 'operation order not as expected from an imod alignment'
        tx, ty = a['AlignmentTransX'][r], a['AlignmentTransY'][r]
        shift = [tx, ty]
        mag = a['Magnification'][r]
        rot = a['InPlaneRotation'][r]

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
        out = general_transform2d(data, shift=shift, scale=float(mag), rot=float(rot), order=order, crop=True,
                                  center=center)
        # (s // 2 - 1) gave best so far
        res = vol2npy(out).copy()

        # crop volumes
        res = res[offset:-offset, offset:-offset]
        indata = indata[offset:-offset, offset:-offset]

        # bin them
        res = resize((res - res.mean()) / res.std(), 1 / 8)
        indata = resize((indata - indata.mean()) / indata.std(), 1 / 8)

        # apply gaussian to pytom image
        res = self.gaussian_filter(res, 30, 240)

        # calculate xcorr between the two
        ccc = nxcc(res, indata)

        if plot:
            from pytom.plotting.plottingFunctions import subplots, show
            fig, ax = subplots(1, 3, figsize=(30, 15))

            ax[0].imshow(res)
            ax[1].imshow(indata)
            ax[2].imshow((res-indata))

            xx,yy = 20 , 30
            ax[0].text(xx, yy, "PyTom Aligned", color="orange")
            ax[1].text(xx, yy, "IMOD Aligned", color="orange")
            ax[2].text(xx, yy, "Abs Diff", color="orange")
            ax[0].set_title(f'Rotation Order: {list2str(order)} ({"".join(map(str, order))})\nNXCC: {ccc:.3f}')

            fig.tight_layout()
            show()

        self.assertTrue(ccc > 0.98, f'cross-correlation not good enough: {ccc:.3f}')

    def checkAlignmentPyTom(self, plot=False):
        fname  = f'{self.rawdir}/alignmentResults.txt'
        a = loadtxt(fname, dtype=DATATYPE_ALIGNMENT_RESULTS_RO)

        square = True
        before = False

        r = self.imageID
        fname = a['FileName'][r]

        if square:

            raw = readC(fname)
            if before:
                dd = raw.numelem()
                raw = resize(raw, 1/8)[0]
                raw /= dd/raw.numelem()
                bin = 8
                data = vol(480,480,1)
            else:
                data = vol(3838, 3838, 1)
                bin =1

            pasteCenter(raw, data)
        else:

            data = readC(fname)

            if before:
                bin = 8
                data = resize(data, 1/8)[0]
            else:
                bin=1

        mm = f'{self.alidir}/sorted_{r:02d}.mrc'
        indata = readC(mm)
        indata = vol2npy(indata).copy()
        off = 0
        n, best = 0, -1
        if plot:
            from pytom.plotting.plottingFunctions import subplots, show
            fig, ax = subplots(3,6,figsize=(30, 15))
        for order in ((1, 2, 0), (0, 2, 1), (0, 1, 2), (2, 0, 1), (1, 0, 2), (2, 1, 0)):
            tx, ty = a['AlignmentTransX'][r]/bin, a['AlignmentTransY'][r]/bin

            shift = [tx, ty]
            mag = a['Magnification'][r]
            rot = a['InPlaneRotation'][r]

            """
            imod docs say the following:
            from coordinate Xi and Yi in the input image to Xo and Yo in output.
            Xc and Yc is center of input image. nxb nyb is the shape of the output image.
            xt and yt are the translations. And a represents the 2x2 matrix for rotation, scaling, and skewing.
            Xo = a11(Xi - Xc) + a12(Yi - Yc) + nxb/2. + xt  
            Yo = a21(Xi - Xc) + a22(Yi - Yc) + nyb/2. + yt
            """

            center = ((data.sizeX() - 1) / 2., (data.sizeY() - 1) / 2., 0)
            out = general_transform2d(data, shift=shift, scale=float(mag), rot=float(rot), order=order, crop=True,
                                      center=center)

            if not before:
                out = resize(out, 1 / 8)[0]  # could try to resize with different center
                out /= data.numelem() / out.numelem()
                res = vol2npy(out).copy()
                res[res<20] = 621
            else:
                res = vol2npy(out).copy()

            res = self.gaussian_filter(res, 30, 240)

            if square:
                res = res[8:-8,:]

            t, e = 32, 400
            ccc = nxcc(res[t+off:e,t:e], indata[t:e-off, t:e])

            if plot:
                ax[0][n].imshow(res[t+off:e,t:e])
                ax[1][n].imshow(indata[t+off:e,t:e])
                ax[2][n].imshow((res-indata)[t+off:e,t:e])
                if order == (2, 1, 0):
                    print(list2str(order))
                    write(f'imod_{list2str(order)}.mrc', indata[t+off:e,t:e])
                    write(f'pytom_{list2str(order)}.mrc', res[t+off:e,t:e])

                xx,yy = 20,30
                ax[0][n].text(xx, yy, "PyTom Aligned", color="orange")
                ax[1][n].text(xx, yy, "IMOD Aligned", color="orange")
                ax[2][n].text(xx, yy, "Abs Diff", color="orange")
                ax[0][n].set_title(f'Rotation Order: {list2str(order)} ({"".join(map(str, order))})\nNXCC: {ccc:.3f}')

            if best < ccc:
                best = ccc
                bestOrder = order
            n += 1

        if plot:
            for i in range(18):
                x,y=i//6, i%6
                ax[x][y].set_xticks([])
                ax[x][y].set_yticks([])

            fig.tight_layout()
            show()

        self.assertTrue(best > 0.8, f'cross-correlation not good enough: {best:.3f}')
        self.assertTrue(bestOrder in ((2,0,1), (1,0,2), (2, 1, 0)), f'best order not expexted: {order}, '
                                                                    f'{list2str(order)}')

    def gaussian_filter(self, data, high, sigma):
        a = create_circle(data.shape, high, sigma)
        return ifftn(fftn(data) * fftshift(a)).real

    def checkPyTomMimickingImod(self):
        # pytom alignment after converting .log and .xf file
        # alignmentresults contains all the jazz
        outfile1 = 'tomogram_PyTom.mrc'
        cmdWBP = self.WBPFormat.format(self.rawdir, 'alignmentResults.txt', outfile1, 8)
        self.submit(cmdWBP)

        # reconstruction from IMOD unbinned, aligned tilts
        # alignmentresults only contains the tilt angles
        outfile2 = 'tomogram_IMOD_ALI_bin1.mrc'
        cmdWBP = self.WBPFormat2.format(self.alidir, 'alignmentResults.txt', outfile2, 8, '.', 'unbinned_sorted')
        self.submit(cmdWBP)

        # check the shift between these two
        self.check_shift(f'{self.alidir}/{outfile2}', f'{self.rawdir}/{outfile1}', exp=[1, 0, 0])
        # expected shift here is not fully zero. pytom will bin the image before alignment
        # which apparently changes the alignment center slightly.

    def runTest(self):
        self.checkReconfromIMODali()
        self.checkAlignmentPyTom(plot=False)
        self.checkAlignmentUnbinned(plot=False)
        self.checkPyTomMimickingImod()


if __name__ == '__main__':
    unittest.main()
