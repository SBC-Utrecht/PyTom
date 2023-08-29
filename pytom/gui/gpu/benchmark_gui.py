import numpy as np
from voltools import Volume
from pycuda import autoinit
from pycuda import gpuarray as gu
from pycuda import driver
from pycuda.cumath import fabs
from pycuda.tools import DeviceData
from pytom.agnostic.tools import paste_in_center

from typing import Union, Tuple

from pycuda.reduction import ReductionKernel
from pycuda.elementwise import ElementwiseKernel
from pycuda.compiler import SourceModule
from skcuda.misc import sum, max, min
from skcuda.fft import fft, ifft, Plan
from skcuda.linalg import conj

max_threads = DeviceData().max_threads
                                                                                
cc_mod = SourceModule("""
__global__ void update_scores_angles(float *scores, float *angles, float *ccc_map, float angleId, int num_elements, int dimx)
{ 
   const int idx = (threadIdx.x + blockIdx.x*dimx)*dimx;

   
   for (int i=0; i < dimx; i++) {
       if (idx +i < num_elements){
           if (scores[idx+i] < ccc_map[idx+i]) {
               scores[idx+i] = ccc_map[idx+i];
               angles[idx+i] = angleId;
            }
       }       
   }
   
}

__global__ void pasteCenter( float *paddedVolume, int spx, int spy, int spz, float *originalVolume, int sox, int soy, int soz)
{
    
	int indX = blockIdx.x * blockDim.x + threadIdx.x;
	int strideX = blockDim.x * gridDim.x;

	int indY = blockIdx.y * blockDim.y + threadIdx.y;
	int strideY = blockDim.y * gridDim.y;

	int indZ = blockIdx.z * blockDim.z + threadIdx.z;
	int strideZ = blockDim.z * gridDim.z;

    int sizePaddingX = (spx-sox)/2;
    int sizePaddingY = (spy-soy)/2;
    int sizePaddingZ = (spz-soz)/2;
        

    for(int i = indX; i < sox; i += strideX)
    {
        for(int j = indY; j < soy; j += strideY)
        {
            for(int k = indZ; k < soz; k += strideZ)
            {
                paddedVolume[(i+sizePaddingX)*spy*spz + (j+sizePaddingY)*spz + k + sizePaddingZ] = originalVolume[i*soy*soz + j*soz + k];
            }
        }
    }
    
}


""")

"""
    const int j = blockIdx.x * blockDim.x + threadIdx.x + paddingSize
    const ini ii= 
    
    
    const int idx =  threadIdx.x + * sVol + threadIdx.y * sVol * sVol;
    int sizePadding = (sPad-sVol)/2;
    
    int pidx = sizePadding + (threadIdx.x+sizePadding)*sPad + (threadIdx.y+sizePadding)*sPad*sPad;
    if (idx < sVol*sVol*sVol){
        for (int i=0; i < sVol; i++) {
            if (volume[idx+i){
            padded_volume[pidx+i] = volume[idx+i];
        }
    }
    __syncthreads();
}

"""


linearAdd = ElementwiseKernel(
        "double *d_t, double *d_m, double a, double b",
        "d_t[i] = ((d_t[i] - a) * b) * d_m[i]",
        "linear_combination")



update_scores_angles = cc_mod.get_function('update_scores_angles')
paste_in_center_gpu = cc_mod.get_function("pasteCenter")

class TemplateMatchingPlan():
    def __init__(self, volume, template, mask, wedge, std_v, gpu=True):
        self.volume = gu.to_gpu(volume)

        self.template = Volume(template)
        self.templatePadded = gu.zeros_like(self.volume, dtype=np.float32)

        self.mask = Volume(mask)
        self.maskPadded = gu.zeros_like(self.volume, dtype=np.float32)
        self.sOrg = mask.shape
        self.sPad = volume.shape
        print(self.sPad, self.sOrg)
        rotate(self.mask, [0,0,0], self.maskPadded, self.sPad, self.sOrg)
        #paste_in_center_gpu(self.template.d_data, self.templatePadded, np.int32(self.sPad), np.int32(self.maskSize), block=(10, 10, 10), grid=(8,1,1))
        #rotate(self.template, [0, 0, 0], self.templatePadded, self.sPad, self.maskSize)
        print(volume.shape, std_v.shape, wedge.shape)
        self.wedge = gu.to_gpu(wedge)
        self.std_v = gu.to_gpu(std_v)

        self.fwd_plan = Plan(volume.shape, volume.dtype, np.complex64)
        self.inv_plan = Plan(volume.shape, np.complex64, volume.dtype)

        self.volume_fft = gu.zeros_like(self.volume, dtype=np.complex64)
        self.template_fft = gu.zeros_like(self.volume, dtype=np.complex64)

        self.ccc_map = gu.zeros_like(self.volume, dtype=np.float32)
        self.norm_volume = np.prod(volume.shape)

        self.scores = gu.ones_like(self.volume, dtype=np.float32)*-1000
        self.angles = gu.ones_like(self.volume, dtype=np.float32)*-1000
        self.p = sum(self.mask.d_data)



    pass

def prepare_template_matching(volume, template, mask, wedge, std_v, gpu=True):
    plan = TemplateMatchingPlan(volume, template, mask, wedge, std_v, gpu=gpu)
    return plan

def cross_correlate(plan, normalize=True):
    #norm_template = sum(plan.mask.d_data)

    fft(plan.volume, plan.volume_fft, plan.fwd_plan)
    fft(plan.templatePadded, plan.template_fft, plan.fwd_plan)

    conj(plan.template_fft, overwrite=True)
    volume_fft = plan.volume_fft * plan.template_fft
    ifft(volume_fft, plan.ccc_map, plan.inv_plan, scale=True)

    plan.ccc_map /= np.float32(plan.p.get())*plan.std_v



def rotate(volume, angles, out, sPad, sOrg):
    spx, spy, spz = sPad
    sox, soy, soz = sOrg
    volume.transform(rotation=angles, rotation_order='szyx', rotation_units='deg', around_center=True)
    paste_in_center_gpu(out, np.int32(spx), np.int32(spy), np.int32(spz), volume.d_data,
                        np.int32(sox), np.int32(soy), np.int32(soz), block=(10, 10, 10), grid=(8,1,1))


def calc_std_v(plan):
    fft(plan.volume, plan.volume_fft, plan.fwd_plan)
    fft(plan.maskPadded, plan.template_fft, plan.fwd_plan)
    conj(plan.template_fft, overwrite=True)
    volume_fft = plan.volume_fft * plan.template_fft
    ifft(volume_fft, plan.std_v, plan.inv_plan, scale=True)
    plan.std_v = plan.std_v / ( plan.p)
    

def meanUnderMask(volume, mask=None, p=1, gpu=False):
    """
    meanValueUnderMask: Determines the mean value under a mask
    @param volume: The volume
    @type volume:  L{pytom.lib.pytom_volume.vol}
    @param mask:  The mask
    @type mask:  L{pytom.lib.pytom_volume.vol}
    @param p: precomputed number of voxels in mask
    @type p: float
    @return: A value (scalar)
    @rtype: single
    @change: support None as mask, FF 08.07.2014
    """

    return sum((volume*mask)) / p

def stdUnderMask(volume, mask, meanValue, p=None, gpu=False):
    """
    stdValueUnderMask: Determines the std value under a mask

    @param volume: input volume
    @type volume:  L{pytom.lib.pytom_volume.vol}
    @param mask: mask
    @type mask:  L{pytom.lib.pytom_volume.vol}
    @param p: non zero value numbers in the mask
    @type p: L{float} or L{int}
    @return: A value
    @rtype: L{float}
    @change: support None as mask, FF 08.07.2014
    """
    res = (meanUnderMask(volume**2, mask, p=p) - meanValue**2)**0.5
    return res

def update_temp(plan):
    meanT = meanUnderMask(plan.template.d_data, plan.mask.d_data, p=plan.p)
    stdT =  1. / stdUnderMask(plan.template.d_data, plan.mask.d_data, meanT, p=plan.p)
    linearAdd(plan.templatePadded, plan.maskPadded, np.float32(meanT.get()), np.float32(stdT.get()))

def template_matching_gpu(volume, template, mask, wedge, std_v, angle_list=[], isSphere=True, return_cpu=False, normalize=True):

    plan = prepare_template_matching(volume, template, mask, wedge, std_v)
    plan.num_vox = np.int32(np.prod(volume.shape))
    dimx = np.int32(volume.shape[0])
    dimy = int(volume.shape[1])
    dimz = int(volume.shape[2])
    print(dimx, dimy, dimz)

    for angleId, angs in enumerate(angle_list):
        rotate(plan.template, angs, plan.templatePadded, plan.sPad, plan.sOrg)#, rotation_order='szxz', return_cpu=False)
        #multiply_wedge(plan.template)
        #if isSphere:
        #rotate(plan.mask, angs)
        #    calc_std_v(plan)

        update_temp(plan)
        cross_correlate(plan)
        update_scores_angles(plan.scores, plan.angles, plan.ccc_map, np.float32(angleId), plan.num_vox, dimx, grid=(dimy, 1),
                             block=(dimz, 1, 1))



    if return_cpu:
        return plan.scores.get(), plan.angles.get(), plan
    else:
        return plan.scores, plan.angles


def plot_central_sections(vol):
    import matplotlib
    try:
        matplotlib.use('Qt5Agg')
    except:
        pass
    from pylab import subplots, show

    fig,ax = subplots(1,3,figsize=(15,5))
    dx,dy,dz = vol.shape
    ax[0].imshow(vol[dx//2,:,:])
    ax[1].imshow(vol[:,dy//2,:])
    ax[2].imshow(vol.sum(axis=2))
    show()

if __name__=='__main__':


    import sys
    from scipy.ndimage import rotate as ROTATE
    from pytom.lib.pytom_freqweight import weight
    from pytom.lib.pytom_numpy import vol2npy
    import pytom.lib.pytom_volume as pytom_volume

    from pytom.agnostic.correlation import meanVolUnderMask, stdVolUnderMask
    import pytom.basic.correlation as corr
    from pytom.basic.files import read as readd
    from pytom.basic.files import write_em
    from pytom.lib.pytom_numpy import vol2npy, npy2vol
    from pytom.lib.pytom_volume import pasteCenter, vol
    import numpy as np

    num_angles, size = map(int, sys.argv[1:3])
    size2=int(sys.argv[3])
    csize =int(sys.argv[4])

    try:
        start, end = map(int, sys.argv[5:7])
    except:
        start, end = 0, csize


    wedge_angle = 30
    wedgeFilter = weight(wedge_angle,0,end-start, size, size)
    wedgeVolume = wedgeFilter.getWeightVolume(True)

    filterVolume = pytom_volume.reducedToFull(wedgeVolume)
    wedgeV = vol2npy(filterVolume).copy()



    from pytom.agnostic.io import read
    import mrcfile

    # NDARRAYS
    voluNDA = mrcfile.open('tomo.mrc', permissive=True).data.copy()
    tempNDA = read('template.em')
    maskNDA = read("mask.em")

    sox, soy,soz = tempNDA.shape
    spx, spy,spz = voluNDA.shape

    #GPUARRAYS
    voluGPU = gu.to_gpu(voluNDA.astype(np.float32))
    tempGPU = gu.to_gpu(tempNDA.astype(np.float32))
    maskGPU = gu.to_gpu(maskNDA)

    #Pytom C arrays
    voluVOL = readd('tomo.mrc')
    tempVOL = readd('template.em')
    maskVOL = readd('mask.em')


    # PASTE TEMPLATE IN CENTER
    tempCentNDA = np.zeros_like(voluNDA)
    tempCentNDA = paste_in_center(tempNDA, tempCentNDA)

    tempCentGPU = gu.zeros_like(voluGPU)
    paste_in_center_gpu(tempCentGPU, np.int32(spx), np.int32(spy), np.int32(spz), tempGPU,
                        np.int32(sox), np.int32(soy), np.int32(soz), block=(10, 10, 10), grid=(8,1,1))
    TCG = tempCentGPU.get()

    tempCentVOL = vol(voluVOL.size_x(), voluVOL.size_y(), voluVOL.size_z())
    tempCentVOL.setAll(0)
    pasteCenter(tempVOL, tempCentVOL)
    TCV = vol2npy(tempCentVOL)

    print('\n center template results')
    for eval in (tempCentNDA, TCG, TCV):
        print((eval-TCV).mean(), eval.min(), eval.max())


    #PASTE MASK IN CENTER
    maskCentNDA = np.zeros_like(voluNDA)
    maskCentNDA = paste_in_center(maskNDA, maskCentNDA)

    maskCentGPU = gu.zeros_like(voluGPU)
    paste_in_center_gpu(maskCentGPU, np.int32(spx), np.int32(spy), np.int32(spz), maskGPU,
                        np.int32(sox), np.int32(soy), np.int32(soz), block=(10, 10, 10), grid=(8,1,1))
    MCG = maskCentGPU.get()

    maskCentVOL = vol(voluVOL.size_x(), voluVOL.size_y(), voluVOL.size_z())
    maskCentVOL.setAll(0)
    pasteCenter(maskVOL, maskCentVOL)
    MCV = vol2npy(maskCentVOL)

    print('\n center mask results (nd, gpu, vol)')
    for eval in (maskCentNDA, MCG, MCV):
        print(eval.mean(), eval.min(), eval.max())


    #CALCULATE STDV
    import pytom.agnostic.correlation as corrNDA
    import pytom.agnostic.testGPUcorr as corrGPU
    import pytom.basic.correlation as corrVOL


    meanNDA = corrNDA.meanVolUnderMask(voluNDA, maskCentNDA, xp=np)
    stdNDA = corrNDA.meanVolUnderMask(voluNDA, maskCentNDA, meanNDA)

    meanVOL = corrVOL.meanUnderMask(voluVOL, maskCentVOL, maskNDA.sum())
    stdVOL  = corrVOL.stdUnderMask(voluVOL, maskCentVOL, maskNDA.sum(), meanVOL)
    write_em('stdV.em', stdVol)
    std_v = vol2npy(stdVOL).copy()
    std_v = std_v.transpose(2,1,0).copy()

    plan = prepare_template_matching(voluNDA, tempNDA, maskNDA, wedge, std_v)



    meanVal = corrVOL.meanValueUnderMask(tempVOL, maskVOL, mask.sum())
    stdVal = corrVOL.stdValueUnderMask(tempVOL, maskVOL, meanVal, mask.sum())

    templateVol = ((tempVOL - meanVal)/stdVal)*maskVOL
    tvol = vol2npy(templateVol)


    meanTGPU = meanUnderMask(tempGPU, maskGPU, p=p)
    stdTGPU  = stdUnderMask(tempGPU, maskGPU, meanTGPU, p=p)
    stdTGPU2 = stdUnderMask(tempGPU, maskGPU, np.float32(meanVal), p=p)

    linearAdd(tempGPU, maskGPU, np.float64(meanTGPU.get()), np.float64(1./stdTGPU.get()))
    tgpu = tempGPU.get()
    print(tvol.mean(), tvol.max(), tvol.min(), tvol.std())
    print(tgpu.mean(), tgpu.max(), tgpu.min(), tgpu.std())
    print()

    print(stdVal, meanVal, stdTGPU, meanTGPU, stdTGPU2)

    import sys
    sys.exit()
    #mask = paste_in_center(mask, np.zeros_like(volume))


    start = driver.Event()
    end = driver.Event()

    start.record()
 
    #scores, angles, plan = template_matching_gpu(volume, temp, mask, wedgeV, np.fft.fftshift(std_v).astype(np.float32), [[0,0,0],]*num_angles, return_cpu=True)

    end.record()
    end.synchronize()
    print('exec time (s):', start.time_till(end) * 1e-3, 7000*start.time_till(end) * 1e-3/num_angles/60)

    map = ((plan.templatePadded.get()))
    print(map.max(), map.min(), map.mean())
    mrcfile.new('templatePadded.mrc',map.astype(np.float32),overwrite=True)
    mrcfile.new('in.mrc',temp*mask,overwrite=True)
