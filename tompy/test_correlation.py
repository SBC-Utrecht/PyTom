import numpy as np
#from voltools import Volume
from pycuda import autoinit
from pycuda import gpuarray as gu
from pycuda import driver
from pycuda.cumath import fabs
from pycuda.tools import DeviceData

from typing import Union, Tuple

from pycuda.reduction import ReductionKernel
from pycuda.elementwise import ElementwiseKernel
from pycuda.compiler import SourceModule

from skcuda.fft import fft, ifft, Plan
from skcuda.linalg import conj

max_threads = DeviceData().max_threads
                                                                                
cc_mod = SourceModule("""
__global__ void update_scores_angles(float *scores, float *angles, float *ccc_map, float angleId, int num_elements, int dimx)
{ 
   const int idx = threadIdx.x + threadIdx.y*blockDim.x + blockIdx.x*blockDim.x*blockDim.y;

   
   for (int i=0; i < dimx; i++) {
       if (idx+i < num_elements) {
           if (scores[idx+i] < ccc_map[idx+i]) {
               scores[idx+i] = ccc_map[idx+i];
               angles[idx] = angleId;
           }
       }       
   }
}
""")

kernel = ElementwiseKernel( "float *scores, float *angles, float *ccc_map, float angleId, int num_elements", 
                            """
    if (i < num_elements) {
        if (scores[i] < ccc_map[i]) {
            scores[i] = ccc_map[i];
            angles[i] = angleId;
        }
    }""", 'update_scores_angles2' )

update_scores_angles = cc_mod.get_function('update_scores_angles')


class TemplateMatchingPlan():
    def __init__(self, volume, template, mask, gpu):
        self.gpu = gpu
        self.volume = gu.to_gpu(volume)
        self.template = Volume(template)
        self.mask = gu.to_gpu(mask)
        self.fwd_plan = Plan(volume.shape, volume.dtype, xp.complex64)
        self.inv_plan = Plan(volume.shape, xp.complex64, volume.dtype)
        self.volume_fft = gu.zeros_like(self.volume, dtype=xp.complex64)
        self.template_fft = gu.zeros_like(self.template.d_data, dtype=xp.complex64)
        self.ccc_map = gu.zeros_like(self.volume, dtype=xp.float32)
        self.norm_volume = xp.prod(volume.shape)
        self.scores = gu.zeros_like(self.volume, dtype=xp.float32)
        self.angles = gu.zeros_like(self.volume, dtype=xp.float32)
    pass

def prepare_template_matching(volume, template, mask, gpu=True):
    plan = TemplateMatchingPlan(volume, template, gpu=gpu)
    return plan

def cross_correlate(plan, normalize=True):
    #norm_template = xp.sum(plan.template.d_data)
    fft(plan.volume, plan.volume_fft, plan.fwd_plan)
    fft(plan.template.d_data, plan.template_fft, plan.fwd_plan)
    conj(plan.template_fft, overwrite=True)
    volume_fft = plan.volume_fft * plan.template_fft
    ifft(volume_fft, plan.ccc_map, plan.inv_plan)
    #if normalize:
    #    plan.ccc_map = plan.ccc_map / plan.norm_volume / norm_template

def rotate(volume, angles):
    volume.transform(rotation=angles, rotation_order='szyx', rotation_units='deg', around_center=True)

def template_matching_gpu(volume, template, mask=None, angle_list=[], return_cpu=False, normalize=True):
    n_threads = min(int(volume.shape[0]), int(xp.floor(xp.sqrt(max_threads))) )
    n_blocks  = int(xp.ceil(xp.prod(volume.shape[1:]) / n_threads**2))

    plan = prepare_template_matching(volume, template, gpu=True)
    num_vox = xp.int32(xp.prod(volume.shape))
    dimx = xp.int32(volume.shape[0])
    for angleId, angs in enumerate(angle_list):
        rotate(plan.template, angs)#, rotation_order='szxz', return_cpu=False)
        cross_correlate(plan)
        #kernel(plan.scores, plan.angles, plan.ccc_map, xp.float32(angleId), xp.int32(a))
        update_scores_angles(plan.scores, plan.angles, plan.ccc_map, xp.float32(angleId), num_vox, dimx, grid=(n_blocks,1), block=(n_threads, 1, 1))

    if return_cpu:
        return plan.scores.get(), plan.angles.get()
    else:
        return plan.scores, plan.angles



if __name__=='__main__': 
    import sys
    from scipy.ndimage import rotate as ROTATE
    num_angles, size = map(int, sys.argv[1:3])
    
    template = xp.zeros((size, size, size), dtype=xp.float32)

    temp = xp.zeros((64,64,64))
    temp[16:-16,16:-16,16:-16] = 1.
    
    for i in range(3):
        tr = ROTATE(temp,30*i,axes=(2,0),reshape=False,order=5)
        ii = i+1
        x,y,z = ii*75, ii*75, ii*75
        template[x-32:x+32, x-32:x+32, x-32:x+32] = tr

    volume = xp.random.rand(size, size, size)
    volume = volume.astype(xp.float32)
    print(volume.shape)
    start = driver.Event()
    end = driver.Event()

    start.record()

    scores, angles = template_matching_gpu(volume, template, [[0,0,0],]*num_angles, return_cpu=True)

    end.record()
    end.synchronize()

    print('exec time (s):', start.time_till(end)*1e-3)

    norm_ccc_map = scores / xp.prod(template.shape) / template.sum()

