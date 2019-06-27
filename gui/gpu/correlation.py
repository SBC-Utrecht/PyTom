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
""")


update_scores_angles = cross_correlatiob_mod.get_function('update_scores_angles')


def cross_correlation(volume, template, angle_list=[], return_cpu=False, normalize=True):
    # todo: support for float64/complex128
    # todo: check sizes
    # todo: pad the template
    # todo: add support for gpuarray and numpy

    plan = Plan(volume.shape, volume.dtype, np.complex64)
    inv_plan = Plan(volume.shape, np.complex64, volume.dtype)
    volume_fft = gu.zeros_like(volume.d_data, dtype=np.complex64)
    template_fft = gu.zeros_like(template.d_data, dtype=np.complex64)
    ccc_map = gu.zeros_like(volume.d_data, dtype=np.float32)
    norm_volume = gu.prod(volume.d_data.shape)

    multiple = 1
    if type(angle_list) == tuple and len(angle_list) == 3:
        angle_list = [angle_list]
        multiple = 0
    if multiple:
        scores = gu.zeros_like(volume.d_data, dtype=np.float32)
        angles = gu.zeros_like(volume.d_data, dtype=np.int)

    cnt = 0
    for angs in angle_list:
        rotate(template, rotation=angs, rotation_order='szxz', return_cpu=False)
        norm_template = np.sum(template.d_data)
        fft(volume.d_data, volume_fft, plan)
        fft(template.d_data, template_fft, plan)
        conj(template_fft, overwrite = True)
        fft_ccmap = volume_fft * template_fft
        ifft(fft_ccmap, ccc_map, inv_plan)
        if normalize:
            ccc_map = ccc_map / norm_volume / norm_template
        if multiple:
            update_scores_angles(scores, angles, ccc_map, cnt, block=(10,10,10))
        cnt += 1

    if return_cpu:
        if multiple:
            return scores.get(), angles.get()
        else:
            return ccc_map.get()

    else:
        if multiple:
            return scores, angles
        else:
            return ccc_map

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
    template = np.zeros((128, 128, 128), dtype=np.float32)
    template[32:64, 32:64, 32:64] = 100
    # template[64, 64, 64] = 100


    volume = np.random.rand(128, 128, 128)
    # volume = np.random.rand(128, 128, 128)
    volume = volume.astype(np.float32)
    # volume[64, 64, 64] = 0

    vol = Volume(volume)
    tem = Volume(template)

    ccc_map = cross_correlation(vol, tem, (0, 0, 0), return_cpu=True)

    print('max', ccc_map.max())# / template.sum())
    print('min', ccc_map.min())

    norm_ccc_map = ccc_map / np.prod(template.shape) / template.sum()
    print('n max', norm_ccc_map.max())# / template.sum())
    print('n min', norm_ccc_map.min())