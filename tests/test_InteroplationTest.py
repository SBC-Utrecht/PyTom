import unittest, os, numba
import numpy as np
import pytom.voltools as vt
import time
from pytom_volume import vol, transform, transformCubic, transformSpline, transformFourierSpline, variance
from pytom.agnostic.interpolation import fill_values_real_spline, fill_values_real_spline_parallel


class pytom_InterpolationTest(unittest.TestCase):
    def setUp(self):
        """set up"""
        self.eps = 1e-5

        # a point somewhere
        box_size = 64
        self.dims = (box_size, ) * 3
        self.center = tuple([(d - 1) / 2 for d in self.dims])
        self.point = (int(10. / 128. * box_size), int(100. / 128. * box_size), int(27. / 128. * box_size))
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
        box.setV(1, *self.point)

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
            self.assertTrue(expr=variance(diff, False) < self.eps,
                            msg=f"point box after forward backward rotation differs too much for pt {name}")

    def test_voltools_cpu(self):
        box = np.zeros(self.dims, dtype=np.float32)
        box[self.point[0], self.point[1], self.point[2]] = 1.

        box_rot = np.zeros_like(box)
        box_org = np.zeros_like(box)

        for name in self.vt_interpolation:
            vt.transform(box, rotation=self.forward, rotation_order=self.order, center=self.center, output=box_rot,
                         device='cpu', interpolation=name)
            vt.transform(box_rot, rotation=self.backward, rotation_order=self.order, center=self.center, output=box_org,
                         device='cpu', interpolation=name)
            print(f'variance for vt {name} interpolation on cpu: {(box - box_org).var()}')
            self.assertTrue(expr=(box - box_org).var() < self.eps,
                            msg=f"point box after forward backward rotation differs too much for vt cpu {name}")

    @unittest.skipIf(os.environ.get('AM_I_IN_A_DOCKER_CONTAINER', False),
                     "The test below uses a GPU and cannot execute in a docker environment")
    def test_voltools_gpu(self):
        try:
            import cupy as cp
        except ImportError:
            raise unittest.SkipTest("No working cupy install found.")
        box = cp.zeros(self.dims, dtype=cp.float32)
        box[self.point[0], self.point[1], self.point[2]] = 1.

        box_rot = cp.zeros_like(box)
        box_org = cp.zeros_like(box)

        for name in self.vt_interpolation:
            vt.transform(box, rotation=self.forward, rotation_order=self.order, center=self.center, output=box_rot,
                         device='gpu:0', interpolation=name)
            vt.transform(box_rot, rotation=self.backward, rotation_order=self.order, center=self.center, output=box_org,
                         device='gpu:0', interpolation=name)
            print(f'variance for vt {name} interpolation on gpu: {(box - box_org).var()}')
            self.assertTrue(expr=(box - box_org).var() < self.eps,
                            msg=f"point box after forward backward rotation differs too much for vt gpu {name}")

    def test_numba_interpolation(self):
        box = np.zeros(self.dims, dtype=np.float32)
        box[self.point[0], self.point[1], self.point[2]] = 1.
        box_rot = np.zeros_like(box)
        box_org = np.zeros_like(box)

        mtx_forward = vt.utils.transform_matrix(rotation=self.forward, rotation_order=self.order, center=self.center)
        mtx_backward = vt.utils.transform_matrix(rotation=self.backward, rotation_order=self.order, center=self.center)

        fill_values_real_spline(box, box_rot, mtx_forward)
        fill_values_real_spline(box_rot, box_org, mtx_backward)

        # ============= time the code
        # t1 = time.time()
        # for i in range(10):
        #     box_rot *= 0
        #     box_org *= 0
        #     fill_values_real_spline(box, box_rot, mtx_forward)
        #     fill_values_real_spline(box_rot, box_org, mtx_backward)
        # print('execution numba 1 thread: ', (time.time() - t1) / 10)

        print(f'variance for spline interpolation with numba: {(box - box_org).var()}')
        self.assertTrue(expr=(box - box_org).var() < self.eps,
                        msg=f"point box after forward backward rotation differs too much for numba spline")

    @unittest.skipIf(os.environ.get('AM_I_IN_A_DOCKER_CONTAINER', False),
                     "The test below uses multiple cores and cannot execute in a docker environment")
    def test_numba_interpolation_parallel(self):
        print('before parallel: ', numba.get_num_threads())
        numba.set_num_threads(4)
        print('after setting for parallel: ', numba.get_num_threads())

        box = np.zeros(self.dims, dtype=np.float32)
        box[self.point[0], self.point[1], self.point[2]] = 1.
        box_rot = np.zeros_like(box)
        box_org = np.zeros_like(box)

        mtx_forward = vt.utils.transform_matrix(rotation=self.forward, rotation_order=self.order, center=self.center)
        mtx_backward = vt.utils.transform_matrix(rotation=self.backward, rotation_order=self.order, center=self.center)

        fill_values_real_spline_parallel(box, box_rot, mtx_forward)
        fill_values_real_spline_parallel(box_rot, box_org, mtx_backward)

        # parallel function only becomes feasible for larger array dimensions than specified here
        t1 = time.time()
        for i in range(10):
            box_rot *= 0
            box_org *= 0
            fill_values_real_spline_parallel(box, box_rot, mtx_forward)
            fill_values_real_spline_parallel(box_rot, box_org, mtx_backward)
        print('execution numba 4 threads: ', (time.time() - t1) / 10)

        print(f'variance for spline interpolation with numba parallel: {(box - box_org).var()}')
        self.assertTrue(expr=(box - box_org).var() < self.eps,
                        msg=f"point box after forward backward rotation differs too much for numba spline parallel")


if __name__ == '__main__':
    unittest.main()
