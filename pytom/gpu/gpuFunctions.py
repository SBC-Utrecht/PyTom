from pytom.gpu.initialize import device, xp

def applyFourierFilter(particle, filter):
    return xp.fft.irfftn(xp.fft.rfftn(particle) * filter)

def applyFourierFilterFull(particle, filter):
    return xp.fft.ifftn(xp.fft.fftn(particle) * filter)

def applyRotatedStaticFourierFilter(particle, filter, rotation, rotation_order='rzxz'):
    return applyFourierFilter(particle, filter.transform(rotation=rotation, rotation_order=rotation_order))

def applyRotatedFilter(particle, filter, rotation, rotation_order='rzxz'):
    return applyFourierFilter(particle, transform(filter, rotation=rotation, rotation_order=rotation_order, device=device))

def applyFourierFilters(particle, filters):
    '''Apply a series of filters in fourier space to a real space object
    @param particle: 2D/3D image or volume
    @type particle: ndarray
    @param filters: list of filters
    @type filters: list of cupy or numpy ndarrays
    @return: 2D/3D filtered image/volume
    @rtype: cupy or numpy ndarray'''
    ft_particle = xp.fft.rfftn(particle)
    for filter in filters:
        ft_particle *= filter
    return xp.fft.irfftn(ft_particle)

def cross_correlation(volume_fft, template, norm, p):
    return xp.fft.fftshift((xp.fft.irfftn(volume_fft * xp.fft.rfftn(template).conj())) /(norm * p)).real

def find_coords_max_ccmap(volume):
    return xp.unravel_index(volume.argmax(),volume.shape)


def subPixelPeak(inp, k=1, ignore_border=1):
    from pytom.voltools import transform
    is2d = len(inp.shape.squeeze()) == 2
    if is2d:
        inp2 = inp.squeeze()[ignore_border:-ignore_border, ignore_border:-ignore_border]
        scale = (k,k,1)
        out = xp.array(inp.squeeze(), dtype=xp.float32)
        x, y = xp.array(xp.unravel_index(inp2.argmax(), inp2.shape)) + ignore_border
        out = xp.expand_dims(out, 2)
        translation = [inp.shape[0] // 2 - x, inp.shape[1] // 2 - y, 0]
    else:
        inp2 = inp[ignore_border:-ignore_border, ignore_border:-ignore_border, ignore_border:-ignore_border]
        scale = (k,k,k)
        out = xp.array(inp, dtype=xp.float32)
        x, y, z = xp.array(xp.unravel_index(inp2.argmax(), inp2.shape)) + ignore_border
        translation = [inp.shape[0] // 2 - x, inp.shape[1] // 2 - y, inp.shape[2]//2 - z]

    zoomed = xp.zeros_like(out, dtype=xp.float32)
    transform(out, output=zoomed, scale=scale, translation=translation, device='gpu:0', interpolation='filt_bspline')

    shifts = [x2,y2,z2] = xp.unravel_index(zoomed.argmax(), zoomed.shape)

    return [zoomed[x2,y2,z2], shifts[3-is2d]]

def subPixelShifts(volume, scale=[1,1,1], cutout=20, shifts=True):
    from pytom.voltools import transform
    import numpy as np

    volC = xp.array([volume.shape[0]//2, volume.shape[1]//2, volume.shape[2]//2])
    cx,cy,cz = center = find_coords_max_ccmap(volume)

    x,y,z = volume.shape
    if cx-cutout//2 < 0 or  cy-cutout//2 < 0 or  cz-cutout//2 < 0 or  cx+cutout//2 >= x or  cy+cutout//2 >=y or  cz+cutout//2 >= z:
        return [volume.max(), xp.array([cx-x//2, cy-y//2, cz-z//2])]


    cutvol = volume[cx-cutout//2:cx+cutout//2, cy-cutout//2:cy+cutout//2, cz-cutout//2:cz+cutout//2]
    nx,ny,nz = [int(np.ceil(cutout / s)) for s in scale]
    nC = xp.array([nx//2, ny//2, nz//2])
    centered = xp.zeros((nx,ny,nz),dtype=xp.float32)
    centered[nx//2-cutout//2:nx//2+cutout//2, ny//2-cutout//2:ny//2+cutout//2, nz//2-cutout//2:nz//2+cutout//2] = cutvol
    interpolated_peak = xp.zeros_like(centered)
    transform(centered, scale=scale, device=device, interpolation='filt_bspline', output=interpolated_peak)
    interpol_center = find_coords_max_ccmap(interpolated_peak)
    shifts = xp.array(center) - volC + (xp.array(interpol_center) - xp.array(nC)) * xp.array(scale)

    return [interpolated_peak.max(), shifts]

def subPixelShiftsImage(volume, scale=[1,1,1], cutout=20, shifts=True):
    from pytom.voltools import transform
    import numpy as np


    if shifts: volC = xp.array([volume.shape[0]//2, volume.shape[1]//2, volume.shape[2]//2])
    else: volC = xp.array([0,0,0],dtype=xp.int32)



    cx,cy,cz = center = find_coords_max_ccmap(volume)
    cutvol = volume[cx-cutout//2:cx+cutout//2, cy-cutout//2:cy+cutout//2, :]
    nx,ny,nz = [int(np.ceil(cutout / s)) for s in scale]
    nC = xp.array([nx//2, ny//2, 0])
    centered = xp.zeros((nx,ny,1),dtype=xp.float32)
    centered[nx//2-cutout//2:nx//2+cutout//2, ny//2-cutout//2:ny//2+cutout//2, :] = cutvol

    interpolated_peak = transform(centered, scale=scale, device=device, interpolation='filt_bspline')


    interpol_center = find_coords_max_ccmap(interpolated_peak)
    shifts = xp.array(center) - volC + (xp.array(interpol_center) - xp.array(nC)) * xp.array(scale)

    return [interpolated_peak.max(), shifts]



def add_transformed_particle_to_sum(particle, sum, rotation=None, rotation_order='rzxz', translation=None, scale=None):
    from pytom.voltools import transform
    rot = transform(particle, rotation=rotation, translation=translation, scale=scale, rotation_order=rotation_order,
                    device=device, interpolation='filt_bspline')
    sum += xp.array(rot)
    return sum

def meanVolUnderMaskPlanned(volume, mask, ifftnP, fftnP, plan):
    """
    meanUnderMask: calculate the mean volume under the given mask (Both should have the same size)
    @param volume: input volume
    @type volume:  L{numpy.ndarray} L{cupy.ndarray}
    @param mask: mask
    @type mask:  L{numpy.ndarray} L{cupy.ndarray}
    @param p: non zero value numbers in the mask
    @type p: L{int} or L{float}
    @param gpu: Boolean that indicates if the input volumes stored on a gpu.
    @type gpu: L{bool}
    @return: the calculated mean volume under mask
    @rtype:  L{numpy.ndarray} L{cupy.ndarray}
    @author: Gijs van der Schot
    """
    volume_fft = fftnP(volume.astype(xp.complex64), plan=plan)
    mask_fft = fftnP(mask.astype(xp.complex64), plan=plan)
    res = xp.fft.fftshift(ifftnP(volume_fft * xp.conj(mask_fft), plan=plan)) / mask.sum()
    return res.real

def stdVolUnderMaskPlanned(volume, mask, meanV, ifftnP, fftnP, plan):
    """
    stdUnderMask: calculate the std volume under the given mask
    @param volume: input volume
    @type volume:  L{numpy.ndarray} L{cupy.ndarray}
    @param mask: mask
    @type mask:  L{numpy.ndarray} L{cupy.ndarray}
    @param p: non zero value numbers in the mask
    @type p: L{int}
    @param meanV: mean volume under mask, which should already been caculated
    @type meanV:  L{numpy.ndarray} L{cupy.ndarray}
    @return: the calculated std volume under mask
    @rtype:  L{numpy.ndarray} L{cupy.ndarray}
    @author: GvdS
    """

    meanV2 = meanV * meanV
    vol2 = volume * volume
    var = meanVolUnderMaskPlanned(vol2, mask, ifftnP, fftnP, plan) - meanV2
    var[var<1E-09] = 1

    return var**0.5

if 'gpu' in device:
    argmax = xp.RawKernel(r'''

    extern "C"  __device__ void warpReduce(volatile float* sdata, volatile int* maxid, int tid, int blockSize) {
        if (blockSize >= 64 && sdata[tid] < sdata[tid + 32]){ sdata[tid] = sdata[tid + 32]; maxid[tid] = maxid[tid+32];} 
        if (blockSize >= 32 && sdata[tid] < sdata[tid + 16]){ sdata[tid] = sdata[tid + 16]; maxid[tid] = maxid[tid+16];}
        if (blockSize >= 16 && sdata[tid] < sdata[tid +  8]){ sdata[tid] = sdata[tid +  8]; maxid[tid] = maxid[tid+ 8];}
        if (blockSize >=  8 && sdata[tid] < sdata[tid +  4]){ sdata[tid] = sdata[tid +  4]; maxid[tid] = maxid[tid+ 4];}
        if (blockSize >=  4 && sdata[tid] < sdata[tid +  2]){ sdata[tid] = sdata[tid +  2]; maxid[tid] = maxid[tid+ 2];}
        if (blockSize >=  2 && sdata[tid] < sdata[tid +  1]){ sdata[tid] = sdata[tid +  1]; maxid[tid] = maxid[tid+ 1];}
    }

    extern "C" __global__ 
    void argmax(float *g_idata, float *g_odata, int *g_mdata, int n) {
        __shared__ float sdata[1024]; 
        __shared__ int maxid[1024];
         
        for (int i=0; i < 1024; i++){ 
            sdata[i] = 0;
            maxid[i] = 0;}
            
        __syncthreads();                                                                                                                              
        
        int blockSize = blockDim.x;                                                                                                                                                   
        unsigned int tid = threadIdx.x;                                                                                                                                        
        int i = blockIdx.x*(blockSize)*2 + tid;                                                                                                                       
        int gridSize = blockSize*gridDim.x*2;                                                                                                                         
        
        while (i < n) {
            if (sdata[tid] < g_idata[i]){
                sdata[tid] = g_idata[i];
                maxid[tid] = i;
            }
            if (sdata[tid] < g_idata[i+blockSize]){
                sdata[tid] = g_idata[i+blockSize];
                maxid[tid] = i + blockSize; 
            }
            i += gridSize; 
        };
         __syncthreads();                                                                                                                                                       
        
        if (blockSize >= 1024){ if (tid < 512 && sdata[tid] < sdata[tid + 512]){ sdata[tid] = sdata[tid + 512]; maxid[tid] = maxid[tid+512]; } __syncthreads(); }
        if (blockSize >=  512){ if (tid < 256 && sdata[tid] < sdata[tid + 256]){ sdata[tid] = sdata[tid + 256]; maxid[tid] = maxid[tid+256]; } __syncthreads(); }
        if (blockSize >=  256){ if (tid < 128 && sdata[tid] < sdata[tid + 128]){ sdata[tid] = sdata[tid + 128]; maxid[tid] = maxid[tid+128]; } __syncthreads(); }
        if (blockSize >=  128){ if (tid <  64 && sdata[tid] < sdata[tid +  64]){ sdata[tid] = sdata[tid +  64]; maxid[tid] = maxid[tid+ 64]; } __syncthreads(); }
        if (tid < 32){ warpReduce(sdata, maxid, tid, blockSize);}                                                                                                                                                                      
        if (tid == 0) {g_odata[blockIdx.x] = sdata[0]; g_mdata[blockIdx.x] = maxid[0];} 
        

    }''', 'argmax')
else:
    argmax = xp.argmax()
