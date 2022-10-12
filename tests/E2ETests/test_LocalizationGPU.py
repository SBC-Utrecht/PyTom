from pytom.gpu.gpuStructures import TemplateMatchingGPU
from pytom.agnostic.io import read
from pytom.agnostic.structures import Wedge
from pytom.angles.globalSampling import GlobalSampling
import unittest
import numpy as np
import time


class pytom_LocalTest(unittest.TestCase):
    """
    pytom_LocalTest: Test case for package localization
    """
    def setUp(self):
        """set up"""
        self.testfilename = f'../testData/emd_1480.map.em_bin_4.em'

    def run_localization(self, gpu_id):
        template = read(self.testfilename)
        mask = template.copy()
        sx, sy, sz = template.shape
        # volume = np.zeros((int(sx * 1.5), int(sy * 2), sz))
        # volume[0:sx, 0:sy, :] = template.copy()
        volume = template.copy()
        SX, SY, SZ = volume.shape
        wedge = Wedge([30., 30.])
        angles = GlobalSampling('angles_38.53_256.em')[:]

        # convolve the search volume with the wedge
        wedge_volume = wedge.returnWedgeVolume(SX, SY, SZ).astype(np.float32)
        volume = np.real(np.fft.irfftn(np.fft.rfftn(volume) * wedge_volume))

        # wedge for template
        wedge = wedge.returnWedgeVolume(sx, sy, sz).astype(np.float32)

        input = (volume, template, mask, wedge, angles, volume.shape)

        tm_process = TemplateMatchingGPU(0, gpu_id, input=input)
        tm_process.start()

        sleep_time, max_sleep_time = 0, 3600
        while tm_process.is_alive() and sleep_time < max_sleep_time:
            time.sleep(1)
            sleep_time += 1

        self.assertTrue(tm_process.completed, msg='something went wrong with TM, process did not complete')

    def test_localization(self):
        self.run_localization(0)


if __name__ == '__main__':
    unittest.main()

