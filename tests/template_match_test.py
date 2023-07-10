import unittest
import numpy as np
import os


only_run_cpu = False

if os.environ.get('AM_I_IN_A_DOCKER_CONTAINER', False):
    only_run_cpu = True
else:
    os.environ['PYTOM_GPU'] = '0'

    from pytom.gpu.initialize import device

    if 'gpu' not in device:
        print('WARNING! GPU support could not be loaded for template matching unittest, only running CPU TM test.')
        only_run_cpu = True


# Remaining imports need to be done after setting GPU env variable:
# * otherwise numpy/cupy backend and the device are globally set before the GPU env variable is specified
from pytom.localization.extractPeaks import templateMatchingGPU, extractPeaks
from pytom.agnostic.io import read
from pytom.basic.files import read_em
from pytom.agnostic.structures import Wedge
from pytom.agnostic.tools import create_sphere
from pytom.angles.globalSampling import GlobalSampling
from pytom.lib.pytom_volume import vol, pasteCenter, initSphere
from pytom.lib.pytom_numpy import vol2npy


class PytomTempMatchTest(unittest.TestCase):
    """
    pytom_LocalTest: Test case for package localization
    """
    def setUp(self):
        """set up"""
        self.testfilename = f'./testData/emd_1480.map.em_bin_4.em'

    def test_localization_cpu(self):
        # read the template as pytom volume
        template, _ = read_em(self.testfilename)
        sx, sy, sz = template.size_x(), template.size_y(), template.size_z()

        # make the mask the same as the template
        mask = vol(sx, sy, sz)
        initSphere(mask, sx//4, 0, 0, sx//2, sy//2, sz//2)

        # create mock tomogram
        volume = vol(int(sx * 1.5), int(sy * 2), sz)
        volume.setAll(.0)
        pasteCenter(template, volume)

        # set missing wedge
        wedge = Wedge([30., 30.])

        # load angles for orientational search
        rotations = GlobalSampling('angles_90_26.em')

        # start cpu template matching
        scores, angles, _, _ = extractPeaks(volume, template, rotations, mask=mask, maskIsSphere=True, wedgeInfo=wedge,
                                            nodeName='0', verboseMode=False, moreInfo=False)

        scores_np = vol2npy(scores).copy()
        angles_np = vol2npy(angles).copy()

        ind = np.unravel_index(scores_np.argmax(), scores_np.shape)
        angle = rotations[int(angles_np[ind])]

        self.assertTrue(np.array(angle).sum() == .0,
                        msg='angle deviates too much from expected')

    @unittest.skipIf(only_run_cpu,
                     "The test below uses a GPU and cannot be executed in a docker environment")
    def test_localization_gpu(self):
        # template
        template = read(self.testfilename, keepnumpy=True)

        # make the mask the same as the template
        mask = create_sphere(template.shape, template.shape[0] // 4)

        # create a mock tomogram
        sx, sy, sz = template.shape
        volume = np.zeros((int(sx * 1.5), int(sy * 2), sz))  # create a mock tomogram
        volume[0:sx, 0:sy, :] = template.copy()

        # set missing wedge
        wedge = Wedge([30., 30.])

        # load angles for orientational search
        rotations = GlobalSampling('angles_90_26.em')

        # start gpu template matching
        scores, angles, _, _ = templateMatchingGPU(volume, template, rotations, mask=mask, maskIsSphere=True,
                                                   wedgeInfo=wedge, gpuID=0)

        ind = np.unravel_index(scores.argmax(), scores.shape)
        angle = rotations[int(angles[ind])]

        self.assertTrue(np.array(angle).sum() == .0,
                        msg='angle deviates too much from expected')


if __name__ == '__main__':
    unittest.main()

