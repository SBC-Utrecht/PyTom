import unittest, os
import numpy as np
from pytom_volume import vol, transform, transformCubic, transformSpline, transformFourierSpline, variance
import pytom.voltools as vt


class pytom_InterpolationTest(unittest.TestCase):
    def setUp(self):
        """set up"""
        # a point somewhere
        self.dims = (128, 128, 128)
        self.center = tuple([(d - 1) / 2 for d in self.dims])
        self.x, self.y, self.z = self.dims
        self.order = 'rzxz'
        self.forward = (31.7, 0, 0)  # ZXZ rotation
        self.backward = (-31.7, 0, 0)
        self.vt_interpolation = ['linear', 'bspline', 'bspline_simple',
                                 'filt_bspline', 'filt_bspline_simple']
        self.pt_interpolation = {'linear': transform, 'cubic': transformCubic,
                                 'spline': transformSpline}  #, 'fourier_spline': transformFourierSpline}

    def test_pytom_volume(self):
        box = vol(self.x, self.y, self.z)
        box.setAll(0.)
        box.setV(1, 10, 100, 27)

        box_rot = vol(self.x, self.y, self.z)
        box_org = vol(self.x, self.y, self.z)

        for name, func in self.pt_interpolation.items():
            box_rot.setAll(0.)
            box_org.setAll(0.)

            func(box, box_rot, self.forward[0], self.forward[1], self.forward[2],
                 self.center[0], self.center[1], self.center[2], 0, 0, 0, 0, 0, 0)
            func(box_rot, box_org, self.backward[0], self.backward[1], self.backward[2],
                 self.center[0], self.center[1], self.center[2], 0, 0, 0, 0, 0, 0)

            diff = box - box_org
            print(f'variance for pt {name} interpolation: {variance(diff, False)}')
            self.assertTrue(expr=variance(diff, False) < 1e-6,
                            msg=f"point box after forward backward rotation differs too much for pt {name}")

    def test_voltools_cpu(self):
        box = np.zeros(self.dims, dtype=np.float32)
        box[10, 100, 27] = 1.

        box_rot = np.zeros_like(box)
        box_org = np.zeros_like(box)

        for name in self.vt_interpolation:
            vt.transform(box, rotation=self.forward, rotation_order=self.order, center=self.center, output=box_rot,
                         device='cpu', interpolation=name)
            vt.transform(box_rot, rotation=self.backward, rotation_order=self.order, center=self.center, output=box_org,
                         device='cpu', interpolation=name)
            print(f'variance for vt {name} interpolation on cpu: {(box - box_org).var()}')
            self.assertTrue(expr=(box - box_org).var() < 1e-6,
                            msg=f"point box after forward backward rotation differs too much for vt cpu {name}")

    @unittest.skipIf(os.environ.get('AM_I_IN_A_DOCKER_CONTAINER', False),
                     "The test below uses a GPU and cannot execute in a docker environment")
    def test_voltools_gpu(self):
        import cupy as cp
        box = cp.zeros(self.dims, dtype=cp.float32)
        box[10, 100, 27] = 1.

        box_rot = cp.zeros_like(box)
        box_org = cp.zeros_like(box)

        for name in self.vt_interpolation:
            vt.transform(box, rotation=self.forward, rotation_order=self.order, center=self.center, output=box_rot,
                         device='gpu:0', interpolation=name)
            vt.transform(box_rot, rotation=self.backward, rotation_order=self.order, center=self.center, output=box_org,
                         device='gpu:0', interpolation=name)
            print(f'variance for vt {name} interpolation on gpu: {(box - box_org).var()}')
            self.assertTrue(expr=(box - box_org).var() < 1e-6,
                            msg=f"point box after forward backward rotation differs too much for vt gpu {name}")


if __name__ == '__main__':
    unittest.main()
