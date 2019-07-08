import numpy as np
from pycuda import autoinit
from volume import Volume

from utils.matrices import scale_matrix, shear_matrix, rotation_matrix, translation_matrix
from utils.kernels import get_transform_kernel, gpuarray_to_texture, VoltoolsElementwiseKernel, get_correlation_kernels
from pycuda import gpuarray as gu
from pycuda import driver
from pycuda.cumath import fabs
from typing import Union, Tuple

from pycuda.reduction import ReductionKernel
from pycuda.elementwise import ElementwiseKernel

from skcuda.fft import fft, ifft, Plan
from skcuda.linalg import conj
# from skcuda.linalg import sqrt as cusqrt

cross_correlatiob_mod = driver.SourceModule("""
    _global_ void update_scores_angles(float32 *scores, int *angles, float32 *ccc_map, int angleID)
    {
        int idx = threadIdx.x + threadIdx.y*10 + threadIdx.z*100;
           
        if scores[idx] < ccc_map[idx] {
            scores[idx] = ccc_map[idx];
            angles[idx] = angleID;
        }
          
    }
    __global__ void insert_volume(float32 *padded_volume, float32 *template, int sizePad, int sizeTem)
    {
        int idx = threadIdx.x + threadIdx.y * threadDimx.x + blockSize.x*threadDim.x*blockIdx.x
        int pad_idx = (1 + 2 * idx) * sizePad
        padded_volume[pad_idx] = template[idx]
""")


update_scores_angles = cross_correlatiob_mod.get_function('update_scores_angles')


class TemplateMatchingPlan():
    def __init__(self, volume, template, gpu):
        self.gpu = gpu
        volume_gpu = gu.to_gpu(volume)
        self.fwd_plan = Plan(volume.shape, volume.dtype, np.complex64)
        self.volume_fft = gu.zeros_like(volume_gpu, dtype=np.complex64)
        fft(volume_gpu, self.volume_fft, self.fwd_plan)
        self.template_fft = gu.zeros_like(volume_gpu, dtype=np.complex64)
        self.ccc_map = gu.zeros_like(volume_gpu, dtype=np.float32)
        self.norm_volume = gu.prod(volume_gpu.shape)
        #self.scores = gu.zeros_like(volume_gpu, dtype=np.float32)
        #self.angles = gu.zeros_like(volume_gpu, dtype=np.float32)
        self.padded_volume = gu.zeros_like(volume_gpu, dtype=np.float32)
        del volume_gpu
        self.inv_plan = Plan(volume.shape, np.complex64, volume.dtype)
        self.template = Volume(template)


def prepare_template_matching(volume, template, gpu=True):
    plan = TemplateMatchingPlan(volume, template, gpu=gpu)
    return plan

def pad_volume(plan):
    dx,dy,dz = plan.template.d_data.shape
    px,py,pz = plan.padded_volume.shape

    plan.padded_volume = 0.


def cross_correlate(plan, normalize=True):#volume, template, volume_fft, plan, inv_plan, template_fft, ccc_map, norm_volume, normalize=True,):
    norm_template = np.sum(plan.template.d_data)
    fft(plan.template.d_data, plan.template_fft, plan.fwd_plan)
    conj(plan.template_fft, overwrite=True)
    volume_fft = plan.volume_fft * plan.template_fft
    ifft(volume_fft, plan.ccc_map, plan.inv_plan)
    if normalize:
        plan.ccc_map = plan.ccc_map / plan.norm_volume / norm_template

def template_matching_gpu(volume: np.ndarray, template: np.ndarray, angle_list: Union[tuple, list],
                          return_cpu: bool, normalize: bool) -> Union[np.ndarray, gu.ndarray]:

    plan = prepare_template_matching(volume, template, gpu=True)

    for n, angs in enumerate(angle_list):
        rotate(plan.template, rotation=angs, rotation_order='szxz', return_cpu=False)
        pad_volume(plan, size=volume.shape[0])
        cross_correlate(plan)
        update_scores_angles(plan.scores, plan.angles, ccc_map, np.float32(n), block=(10,10,10))

    if return_cpu:
        return plan.scores.get(), plan.angles.get()
    else:
        return plan.scores, plan.angles

def cpu_correlation(vol1: np.ndarray, vol2: np.ndarray) -> float:
    Xm = np.mean(vol1)
    ym = np.mean(vol2)
    r_num_1 = np.sum((vol1 - Xm) * (vol2 - ym))
    r_den_1 = np.sqrt(np.sum((vol1 - Xm) ** 2) * np.sum((vol2 - ym) ** 2))
    r = r_num_1 / r_den_1
    return r

def correlation(vol1: Union[Volume], vol2: Union[Volume]) -> float:

    if vol1.dtype != vol2.dtype:
        raise ValueError(f'Volumes dtype must be of same type (vol1: {vol1.dtype}; vol2: {vol2.dtype}).')

    # precalculate means
    v1mean = vol1.mean()
    v2mean = vol2.mean()

    # kernels
    num = vol1.kernel('cor_num')
    den = vol1.kernel('cor_den')

    r_num = num(vol1.d_data, vol2.d_data, v1mean, v2mean).get()
    r_den = np.sqrt((den(vol1.d_data, v1mean) * den(vol2.d_data, v2mean)).get())

    return r_num / r_den


if __name__ == '__main__':
    from transforms import *

    np.random.seed(1)
    template = np.zeros((32, 32, 32), dtype=np.float32)
    template[15:, 15:, 15:] = 100
    volume = np.random.rand(128, 128, 128)
    volume = volume.astype(np.float32)


    ccc_map = template_matching_gpu(volume, template, (0, 0, 0), return_cpu=True)

    print('max', ccc_map.max())
    print('min', ccc_map.min())

    norm_ccc_map = ccc_map / np.prod(template.shape) / template.sum()
    print('n max', norm_ccc_map.max())# / template.sum())
    print('n min', norm_ccc_map.min())