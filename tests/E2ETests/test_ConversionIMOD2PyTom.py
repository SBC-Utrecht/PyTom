import os
import sys
from pytom.agnostic.io import read, write
import numpy as np
import unittest
from numpy import sqrt
from os import system as exe
from os.path import join as merge
from pytom.basic.structures import Particle, ParticleList
from os.path import exists
from pytom.basic.datatypes import DATATYPE_ALIGNMENT_RESULTS_RO, fmtAlignmentResultsRo as fmt, HEADER_ALIGNMENT_RESULTS_RO as header
from pytom.basic.files import loadtxt, savetxt
from pytom.agnostic.correlation import FLCF
from pytom.agnostic.io import read

class pytom_MyFunctionTest(unittest.TestCase):

    def create_folder(self, folder):
        if not exists(folder): os.mkdir(folder)

    def submit(self, cmd):
        print(cmd)
        os.system(cmd)

    def setUp(self):

        # tomogram id
        self.id=1

        # Image 33 is chosen as it has the highest shift and scaling of all images in the tilt series
        self.imageID = 33

        func = self.submit

        # Clone data repo if not existing
        if not exists('PyTomUnitTests'):
            gitcmd = 'git clone https://github.com/gijsschot/PyTomUnitTests.git'
            os.system(gitcmd)

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
        func(f'cd {self.datadir}; submfg tilt.com ')
        func(f'cd {self.datadir}; trimvol -yz tomo{self.id}.mrc_full.rec tomogram_IMOD.mrc ')
        func(f'mrcs2mrc.py -f {self.datadir}/tomo{self.id}.mrc.ali -t {self.alidir} --prefix sorted ')
        func(f'mrcs2mrc.py -f {self.datadir}/tomo{self.id}.mrc.st -t {self.rawdir} --prefix sorted ')

        # Markerfile is not used but must exist
        func(f'touch {self.rawdir}/markerfile.txt')
        func(f'touch {self.alidir}/markerfile.txt')

        # create dummy alignmentResults.txt file: scaling =1, other ali params = 0, filename need update
        a = np.zeros((35), dtype=DATATYPE_ALIGNMENT_RESULTS_RO)
        angles = loadtxt(f'{self.datadir}/tomo{self.id}.mrc.tlt')
        a['TiltAngle'] = angles
        a['Magnification'] = 1
        a['OperationOrder'] = 'RTS'

        for n, fname in enumerate([fname for fname in os.listdir(self.alidir) if fname.startswith('sorted')]):
            a['FileName'][n] = fname

        savetxt(f'{self.alidir}/alignmentResults.txt', a, fmt=fmt, header=header)

        # Import IMOD ali parameters into alignmentResults.txt
        cmd = f'convert.py -f {self.datadir}/taSolution.log -o txt --alignxf {self.datadir}/tomo{self.id}.mrc.xf -t {self.rawdir} --sortedFolder {os.getcwd()}/{self.rawdir}'
        func(cmd)

        # WBP command str. Needs formatting in tests
        self.WBPFormat = """cd {}

reconstructTomogram.py \
    --tiltSeriesName ./sorted \
    --firstIndex 0 \
    --lastIndex 34 \
    --referenceIndex 17 \
    --markerFile markerfile.txt \
    --referenceMarkerIndex 0 \
    --expectedRotationAngle 0 \
    --projectionTargets temp_files_unweighted/sorted_aligned \
    --projectionBinning {} \
    --lowpassFilter 0.9  \
    --weightingType -1  \
    --tomogramFile {} \
    --tiltSeriesFormat mrc \
    --fileType mrc  \
    --tomogramSizeX 464  \
    --tomogramSizeY 480 \
    --tomogramSizeZ 464 \
    --reconstructionCenterX 0 \
    --reconstructionCenterY 0 \
    --reconstructionCenterZ 0 \
    --specimenAngle 0 
"""

    def checkReconfromIMODali(self):
        outfile = 'tomogram_IMOD_ALI.mrc'
        cmd = self.WBPFormat.format(self.alidir, 1, outfile)
        os.system(cmd)
        self.check_shift(f'{self.datadir}/tomogram_IMOD.mrc', f'{self.alidir}/{outfile}')

    def check_shift(self, fname1, fname2):
        a = read(fname1)
        b = read(fname2)

        c = FLCF(a,b)

        loc = np.unravel_index(c.argmax(), c.shape)
        shift = (loc[0]-b.shape[0]//2, loc[1]-b.shape[1]//2, loc[2]-b.shape[2]//2)

        print(shift)
        self.assertTrue(np.abs(shift[0]) < 2 and np.abs(shift[1]) < 2 and np.abs(shift[1]) < 2, f'shifts between tomograms: {shift[0]} {shift[1]} {shift[2]}')

    def checkAlignmentPyTom(self):
        from pytom_volume import read as readC, vol, pasteCenter
        from pytom_numpy import vol2npy
        from pytom.agnostic.correlation import nxcc
        from pytom.basic.transformations import general_transform2d, resize
        from pytom.plotting.plottingFunctions import subplots, show
        from pytom.agnostic.normalise import mean0std1
        from pytom.agnostic.tools import convert_operation_order_list2str as list2str

        fname  = f'{self.rawdir}/alignmentResults.txt'
        a = loadtxt(fname, dtype=DATATYPE_ALIGNMENT_RESULTS_RO)

        square = True
        before = False
        plot = True

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
        indata= readC(mm)

        indata = vol2npy(indata).copy()
        off = 0
        n, best  = 0, -1
        if plot:
            fig,ax = subplots(3,6,figsize=(30,15))
        for order in ((1,2,0), (0,2,1), (0,1,2), (2,0,1), (1,0,2), (2,1,0)):
            tx, ty = a['AlignmentTransX'][r]/bin, a['AlignmentTransY'][r]/bin

            shift = [tx,ty] #[-a['AlignmentTransX'][r], -a['AlignmentTransY'][r]]
            mag = a['Magnification'][r]
            rot = a['InPlaneRotation'][r]

            out = general_transform2d(data, shift=shift, scale=float(mag), rot=float(rot), order=order, crop=True,
                                      center=(data.sizeX()/2-1, data.sizeY()/2-1, 0))

            if not before:
                out = resize(out, 1 / 8)[0]
                out /= data.numelem() / out.numelem()
                res = vol2npy(out).copy()
                res[res<20] = 621
            else:
                res = vol2npy(out).copy()

            res = self.gaussian_filter(res, 30, 240)

            if square:
                res = res[8:-8,:]

            print(res.shape, indata.shape)

            t, e = 32, 400
            ccc = nxcc(res[t+off:e,t:e], indata[t:e-off, t:e])


            if plot:
                ax[0][n].imshow(res[t+off:e,t:e])
                ax[1][n].imshow(indata[t+off:e,t:e])
                ax[2][n].imshow(np.abs(res-indata)[t+off:e,t:e])

                xx,yy = 20,30
                ax[0][n].text(xx, yy, "PyTom Aligned", color="orange")
                ax[1][n].text(xx, yy, "IMOD Aligned", color="orange")
                ax[2][n].text(xx, yy, "Abs Diff", color="orange")
                ax[0][n].set_title(f'Rotation Order: {list2str(order)} ({"".join(map(str, order))})\nNXCC: {ccc:.3f}')

            write(f'sorted_aligned_{order[0]}{order[1]}{order[2]}.mrc', res)




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

        self.assertTrue(best > 0.7, f'cross-correlation not good enough: {best:.3f}')
        self.assertTrue(bestOrder in ((0,2,1), (1,2,0), (0,1,2)), order)

    def gaussian_filter(self, data, high, sigma):
        from pytom.agnostic.tools import create_sphere, create_circle
        from numpy.fft import fftshift, fftn, ifftn
        from pytom.plotting.plottingFunctions import subplots, show
        a = create_circle(data.shape, high, sigma)

        return ifftn(fftn(data) * fftshift(a)).real

    def lowpass(self, image, lowpassFilter, imdim):
        from pytom.basic.filter import filter
        import pytom_freqweight
        print('aa')

        lpf = pytom_freqweight.weight(0.0, lowpassFilter * imdim // 2, imdim, imdim // 2 + 1, 1,
                                  lowpassFilter / 5. * imdim)
        filtered = filter(volume=image, filterObject=lpf, fourierOnly=False)
        image = filtered[0]
        return image

    def checkPyTomMimickingImod(self):
        outfile = 'tomogram_PyTom.mrc'
        cmdWBP = self.WBPFormat.format(self.rawdir, 8, outfile)
        os.system(cmdWBP)
        self.check_shift(f'{self.datadir}/tomogram_IMOD.mrc', f'{self.rawdir}/{outfile}')


    def runTest(self):
        self.checkReconfromIMODali()
        self.checkAlignmentPyTom()
        self.checkPyTomMimickingImod()

if __name__ == '__main__':
    unittest.main()
