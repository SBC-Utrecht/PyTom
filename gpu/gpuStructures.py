import threading
import importlib
from pytom.tompy.io import read, write
from pytom.gpu.initialize import device, xp

class TemplateMatchingPlan():
    def __init__(self, volume, template, mask, wedge, cp, vt, calc_stdV, pad, get_fft_plan, deviceid):

        import pytom.voltools as vt
        from pytom.tompy.filter import applyFourierFilter

        self.volume = cp.asarray(volume,dtype=cp.float32,order='C')


        self.sPad = volume.shape

        self.mask = cp.asarray(mask,dtype=cp.float32, order='C')
        self.maskPadded = cp.zeros_like(self.volume).astype(cp.float32)
        self.sOrg = mask.shape
        pad(self.mask, self.maskPadded, self.sPad, self.sOrg)

        self.templateOrig = cp.asarray(template, dtype=cp.float32,order='C')
        self.template = cp.asarray(template, dtype=cp.float32,order='C')
        self.texture = vt.StaticVolume(self.template, interpolation='filt_bspline', device=device)
        self.templatePadded = cp.zeros_like(self.volume)

        self.wedge = cp.asarray(wedge,order='C')
        self.volume_fft = cp.fft.rfftn(self.volume) #* cp.array(wedgeVolume,dtype=cp.float32)

        self.ccc_map = cp.zeros_like(self.volume)

        self.scores = cp.ones_like(self.volume)*-1000
        self.angles = cp.ones_like(self.volume)*-1000
        self.p = self.mask.sum()

        if not get_fft_plan is None:
            self.volume_fft2 = cp.fft.fftn(self.volume)
            self.fftplan = get_fft_plan(self.volume.astype(cp.complex64))
        else:
            self.fftplan = None
            self.volume_fft2 = None

        calc_stdV(self)
        cp.cuda.stream.get_current_stream().synchronize()


class TemplateMatchingGPU(threading.Thread):
    def __init__(self, jobid, deviceid, input):
        threading.Thread.__init__(self)

        import numpy as np
        from pytom.tompy.tools import paste_in_center
        from pytom.tompy.transform import rotate3d
        from pytom_volume import rotateSpline as rotate
        import cupy as cp

        cp.cuda.Device(deviceid).use()

        from cupy import sqrt, float32
        from cupy.fft import fftshift, rfftn, irfftn, ifftn, fftn
        import voltools as vt
        from pytom.tompy.io import write
        from cupyx.scipy.ndimage import map_coordinates


        if 1:
            from cupyx.scipy.fftpack.fft import get_fft_plan
            from cupyx.scipy.fftpack.fft import fftn as fftnP
            from cupyx.scipy.fftpack.fft import ifftn as ifftnP
            self.fftnP = fftnP
            self.ifftnP = ifftnP
        else:
            get_fft_plan = None

        self.cp = cp
        self.map_coordinates = map_coordinates
        self.Device = cp.cuda.Device
        self.jobid = jobid
        self.deviceid = deviceid
        self.active = True
        self.input = input

        self.fftn = fftn
        self.ifftn = ifftn
        self.fftshift = fftshift
        self.rfftn = rfftn
        self.irfftn = irfftn
        self.sqrt = sqrt
        self.float32 = float32


        self.update_scores_angles = cp.RawKernel(r"""
        extern "C"  __global__ void update_scores_angles(float *scores, float *angles, float *ccc_map, float angleId, int num_elements, int dimx)
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
           __syncthreads();

        }
        """, 'update_scores_angles')


        self.updateResFromIdx = cp.ElementwiseKernel(
            'float32 scores, float32 angles, float32 ccc_map, float32 angleId',
            'float32 out, float32 out2',
            'if (scores < ccc_map) {out = ccc_map; out2 = angleId;}',
            'updateResFromIdx')

        self.plan = TemplateMatchingPlan(input[0], input[1], input[2], input[3], cp, vt, self.calc_stdV, self.pad, get_fft_plan, deviceid)

        print("Initialized job_{:03d} on device {:d}".format(self.jobid, self.deviceid))

    def run(self):
        print("RUN")
        if 1:
            self.template_matching_gpu(self.input[4], self.input[5])
            self.completed = True
        else:
            self.completed = False
        self.active = False

    def template_matching_gpu(self, angle_list, dims, isSphere=True, verbose=True):

        self.Device(self.deviceid).use()

        import voltools as vt
        from pytom.tompy.io import write
        sx, sy, sz = self.plan.template.shape
        SX, SY, SZ = self.plan.templatePadded.shape
        CX, CY, CZ = SX // 2, SY // 2, SZ // 2
        cx, cy, cz = sx // 2, sy // 2, sz // 2
        mx, my, mz = sx %  2, sy %  2, sz %  2

        for angleId, angles in enumerate(angle_list):

            # Rotate
            #self.plan.template = self.rotate3d(self.plan.templateOrig, phi=phi,the=the,psi=psi)
            self.plan.texture.transform(rotation=(angles[0], angles[2], angles[1]), rotation_order='rzxz', output=self.plan.template)

            # Add wedge
            #print(self.rfftn(self.plan.template).shape)
            self.plan.template = self.irfftn(self.rfftn(self.plan.template) * self.plan.wedge)

            # Normalize template
            meanT = self.meanUnderMask(self.plan.template, self.plan.mask, p=self.plan.p)
            stdT = self.stdUnderMask(self.plan.template, self.plan.mask, meanT, p=self.plan.p)
            self.plan.template = ((self.plan.template - meanT) / stdT) * self.plan.mask
            print(meanT, stdT, self.plan.p)
            # Paste in center
            self.plan.templatePadded[CX-cx:CX+cx+mx, CY-cy:CY+cy+my,CZ-cz:CZ+cz+mz] = self.plan.template
            self.plan.wedge.shape

            # Cross-correlate and normalize by stdV
            self.plan.ccc_map = self.normalized_cross_correlation(self.plan.volume_fft2, self.plan.templatePadded, self.plan.stdV, self.plan.p, plan=self.plan.fftplan)

            # Update the scores and angles
            self.updateResFromIdx(self.plan.scores, self.plan.angles, self.plan.ccc_map, angleId, self.plan.scores, self.plan.angles)

            #self.cp.cuda.stream.get_current_stream().synchronize()

    def is_alive(self):
        return self.active

    def calc_stdV(self, plan):
        stdV = self.meanVolUnderMask2(plan.volume**2, plan) - self.meanVolUnderMask2(plan.volume, plan)**2
        stdV[stdV < self.float32(1e-09)] = 1
        plan.stdV = self.sqrt(stdV)

    def meanVolUnderMask2(self, volume, plan):
        res = self.fftshift(self.ifftn(self.fftn(volume) * self.fftn(plan.maskPadded).conj())) / plan.mask.sum()
        return res.real

    def meanUnderMask(self,volume, mask=None, p=1, gpu=False):
        """
        meanValueUnderMask: Determines the mean value under a mask
        @param volume: The volume
        @type volume:  L{pytom_volume.vol}
        @param mask:  The mask
        @type mask:  L{pytom_volume.vol}
        @param p: precomputed number of voxels in mask
        @type p: float
        @return: A value (scalar)
        @rtype: single
        @change: support None as mask, FF 08.07.2014
        """
        return (volume * mask).sum() / p

    def stdUnderMask(self, volume, mask, meanValue, p=None, gpu=False):
        """
        stdValueUnderMask: Determines the std value under a mask

        @param volume: input volume
        @type volume:  L{pytom_volume.vol}
        @param mask: mask
        @type mask:  L{pytom_volume.vol}
        @param p: non zero value numbers in the mask
        @type p: L{float} or L{int}
        @return: A value
        @rtype: L{float}
        @change: support None as mask, FF 08.07.2014
        """
        return self.sqrt(self.meanUnderMask(volume**2, mask, p=p) - meanValue**2)

    def normalized_cross_correlation(self, volume_fft, template, norm, p=1, plan=None):
        return xp.real(self.fftshift(self.ifftnP(volume_fft * xp.conj(self.fftnP(template, plan=plan)), plan=plan)) / (p * norm))
        #print(volume_fft.shape)

        print('cc')
        cc = xp.conj(xp.fft.fftn(template))
        print('ccc')
        aa = xp.fft.fftn(template)
        res = xp.abs(xp.fft.fftshift(xp.fft.ifftn(aa * cc )) )
        #res = xp.abs(xp.fft.fftshift(xp.fft.ifftn(volume_fft * xp.fft.fftn(template).conj() )) / (p * norm))
        return res


    def pad(self, volume, out, sPad, sOrg):
        SX, SY, SZ = sPad
        sx, sy, sz = sOrg
        print(sPad, sOrg)
        out[SX // 2 - sx // 2:SX // 2 + sx // 2 + sx % 2, SY // 2 - sy // 2:SY // 2 + sy // 2 + sy % 2,
            SZ // 2 - sz // 2:SZ // 2 + sz // 2 + sz % 2] = volume

    def rotate3d(self, data, phi=0, psi=0, the=0, center=None, order=1, output=None):
        """Rotate a 3D data using ZXZ convention (phi: z1, the: x, psi: z2).

        @param data: data to be rotated.
        @param phi: 1st rotate around Z axis, in degree.
        @param psi: 3rd rotate around Z axis, in degree.
        @param the: 2nd rotate around X axis, in degree.
        @param center: rotation center.

        @return: the data after rotation.
        """
        # Figure out the rotation center
        if center is None:
            cx = data.shape[0] / 2
            cy = data.shape[1] / 2
            cz = data.shape[2] / 2
        else:
            assert len(center) == 3
            (cx, cy, cz) = center

        # Transfer the angle to Euclidean
        phi = -float(phi) * self.cp.pi / 180.0
        the = -float(the) * self.cp.pi / 180.0
        psi = -float(psi) * self.cp.pi / 180.0
        sin_alpha = self.cp.sin(phi)
        cos_alpha = self.cp.cos(phi)
        sin_beta = self.cp.sin(the)
        cos_beta = self.cp.cos(the)
        sin_gamma = self.cp.sin(psi)
        cos_gamma = self.cp.cos(psi)



        # Calculate inverse rotation matrix
        Inv_R = self.cp.zeros((3, 3), dtype=self.cp.float32)

        Inv_R[0][0] = cos_alpha * cos_gamma - cos_beta * sin_alpha * sin_gamma
        Inv_R[0][1] = -cos_alpha * sin_gamma - cos_beta * sin_alpha * cos_gamma
        Inv_R[0][2] = sin_beta * sin_alpha

        Inv_R[1][0] = sin_alpha * cos_gamma + cos_beta * cos_alpha * sin_gamma
        Inv_R[1][1] = -sin_alpha * sin_gamma + cos_beta * cos_alpha * cos_gamma
        Inv_R[1][2] = -sin_beta * cos_alpha

        Inv_R[2][0] = sin_beta * sin_gamma
        Inv_R[2][1] = sin_beta * cos_gamma
        Inv_R[2][2] = cos_beta



        grid = self.cp.mgrid[-cx:data.shape[0]-cx, -cy:data.shape[1]-cy, -cz:data.shape[2]-cz]
        temp = grid.reshape((3, grid.size // 3))
        temp = self.cp.dot(Inv_R, temp)
        grid = self.cp.reshape(temp, grid.shape)
        grid[0] += cx
        grid[1] += cy
        grid[2] += cz

        # Interpolation
        #self.map_coordinates(data, grid, output=output)
        return self.map_coordinates(data, grid.reshape(len(grid), -1), order=order).reshape(grid.shape[1:])


class GLocalAlignmentPlan():
    def __init__(self, particle, reference, mask, wedge, maskIsSphere=True, cp=xp, device='cpu', interpolation='filt_bspline', taper=0):
        id = int(device.split(":")[1])

        cp.cuda.Device(id).use()
        self.cp = cp
        self.device = device
        print(f'start init on: {device}')
        from pytom.voltools import StaticVolume
        from pytom.tompy.tools import taper_edges
        from pytom.tompy.transform import fourier_reduced2full
        from cupyx.scipy.fftpack.fft import get_fft_plan
        from cupyx.scipy.fftpack.fft import fftn as fftnP
        from cupyx.scipy.fftpack.fft import ifftn as ifftnP
        from pytom.tompy.io import read_size

        # Allocate volume and volume_fft
        self.volume        = cp.zeros(read_size(particle.getFilename()), dtype=cp.float32)
        self.volume_fft    = cp.zeros_like(self.volume, dtype=cp.complex64)

        # Allocate planned-fftn-related functions and structure
        self.ifftnP        = ifftnP
        self.fftnP         = fftnP
        self.fftplan       = get_fft_plan(self.volume.astype(cp.complex64))

        # Allocate mask-related objects
        mask               = mask.getVolume().get()
        self.mask          = cp.array(mask, dtype=cp.float32)
        self.rotatedMask   = cp.array(mask, dtype=cp.float32)
        self.maskTex       = StaticVolume(self.mask, device=device, interpolation=interpolation)
        self.mask_fft      = self.fftnP(self.rotatedMask.astype(cp.complex64),plan=self.fftplan)
        self.p             = self.rotatedMask.sum()

        # Allocate wedge-related objects
        self.wedge         = wedge
        self.wedgeAngles   = wedge.getWedgeAngle()
        wedgePart          = wedge.returnWedgeVolume(*self.volume.shape,humanUnderstandable=True).get()
        self.rotatedWedge  = cp.array(wedgePart, dtype=cp.float32)
        self.wedgePart     = cp.fft.fftshift(wedgePart)
        self.wedgeTex      = StaticVolume(self.rotatedWedge.copy(), device=self.device, interpolation=interpolation)
        del wedgePart

        # Allocate reference related objects
        reference          = cp.array((reference.getVolume()).get(), dtype=cp.float32) * self.mask
        self.referenceTex  = StaticVolume(reference, device=device, interpolation=interpolation)
        self.rotatedRef    = reference.copy()
        self.simulatedVolume = cp.zeros_like(self.volume, dtype=cp.float32)
        del reference

        # Allocate taper-related arrays
        dummy, taperMask   = taper_edges(self.volume, self.volume.shape[0]/10)
        self.taperMask          = cp.array(taperMask.get(), dtype=cp.float32)
        self.masked_taper  = self.rotatedMask * self.taperMask
        del dummy, taperMask

        # Allocate required volumes
        self.stdV          = cp.zeros_like(self.volume, dtype=cp.float32)
        self.ccc_map       = cp.zeros_like(self.volume, dtype=cp.float32)
        self.sumParticles  = cp.zeros_like(self.volume, dtype=cp.float32)
        self.sumWeights    = cp.zeros_like(self.wedgePart[:,:,:self.volume.shape[2]//2+1], dtype=cp.float32)

        # Allocate arrays used for normalisation of reference
        self.num_threads   = 512
        self.nblocks       = int(self.cp.ceil(self.simulatedVolume.size / self.num_threads / 2))
        self.fast_sum_mean = self.cp.zeros((4 * self.nblocks * 2), dtype=self.cp.float32)
        self.fast_sum_stdv = self.cp.zeros((4 * self.nblocks * 2), dtype=self.cp.float32)

        # Allocate arrays used in subPixelMas3D
        b                  = max(0, int(self.volume.shape[0] // 2 - 4 / 0.1))
        zx, zy, zz         = self.volume.shape
        croppedForZoom     = self.volume[b:zx - b, b:zy - b, b:zz - b]
        nblocks            = int(self.cp.ceil(croppedForZoom.size / 1024 / 2))
        self.zoomed        = self.cp.zeros_like(croppedForZoom, dtype=cp.float32)
        self.fast_sum      = self.cp.zeros((nblocks), dtype=cp.float32)
        self.max_id        = self.cp.zeros((nblocks), dtype=cp.int32)
        del croppedForZoom

        # Kernels
        self.normalize     = cp.ElementwiseKernel( 'T ref, T mask, raw T mean, raw T stdV ', 'T z', 'z = ((ref - mean[i*0]) / stdV[i*0]) * mask', 'norm2')
        self.sumMeanStdv   = cp.RawKernel(r'''
    __device__ void warpReduce(volatile float* sdata, volatile float* stdv, int tid, int blockSize) {
        if (blockSize >= 64) {sdata[tid] += sdata[tid + 32]; stdv[tid] += stdv[tid+32];}
        if (blockSize >= 32) {sdata[tid] += sdata[tid + 16]; stdv[tid] += stdv[tid+16];}
        if (blockSize >= 16) {sdata[tid] += sdata[tid +  8]; stdv[tid] += stdv[tid+ 8];}
        if (blockSize >=  8) {sdata[tid] += sdata[tid +  4]; stdv[tid] += stdv[tid+ 4];}
        if (blockSize >=  4) {sdata[tid] += sdata[tid +  2]; stdv[tid] += stdv[tid+ 2];}
        if (blockSize >=  2) {sdata[tid] += sdata[tid +  1]; stdv[tid] += stdv[tid+ 1];} }



    extern "C" __global__ 
    void sumMeanStdv(float *g_idata, float *mask, float *g_mean, float * g_stdv,  int n) {
                                                                                                                                                         
        __shared__ float mean[1024];
        __shared__ float stdv[1024];
        int blockSize = blockDim.x; 
        unsigned int tid = threadIdx.x;                                                                                                                                        
        int i = blockIdx.x*(blockSize)*2 + tid;                                                                                                                       
        int gridSize = blockSize*gridDim.x*2;                                                                                                                         

        mean[tid] = 0.;                                                                                                                                                       
        stdv[tid] = 0.;
        
        while (i < n) {
            mean[tid] += g_idata[i] * mask[i] + g_idata[i + blockSize] * mask[i+blockSize];  
            stdv[tid] += g_idata[i] * g_idata[i] * mask[i] + g_idata[i + blockSize] * g_idata[i + blockSize] * mask[i+blockSize];
            i += gridSize;}
         __syncthreads();                                                                                                                                                       

        if (blockSize >= 1024){ if (tid < 512) { mean[tid] += mean[tid + 512]; stdv[tid] += stdv[tid+512];} __syncthreads(); }                                                                                                    
        if (blockSize >= 512) { if (tid < 256) { mean[tid] += mean[tid + 256]; stdv[tid] += stdv[tid+256];} __syncthreads(); }                                                                          
        if (blockSize >= 256) { if (tid < 128) { mean[tid] += mean[tid + 128]; stdv[tid] += stdv[tid+128];} __syncthreads(); }                                                                          
        if (blockSize >= 128) { if (tid <  64) { mean[tid] += mean[tid +  64]; stdv[tid] += stdv[tid+ 64];} __syncthreads(); }                                                                          
        if (tid < 32){ warpReduce(mean, stdv, tid, blockSize);}                                                                                                                                                                      
        if (tid == 0) {g_mean[blockIdx.x] = mean[0]; g_stdv[blockIdx.x] = stdv[0];} 
       

    }''', 'sumMeanStdv')
        self.argmax        = self.cp.RawKernel(r'''

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
                    /* 
                    for (int i=0; i < 1024; i++){ 
                        sdata[i] = 0;
                        maxid[i] = 0;}

                    __syncthreads();                                                                                                                              
                    */
                    int blockSize = blockDim.x;                                                                                                                                                   
                    unsigned int tid = threadIdx.x;                                                                                                                                        
                    int i = blockIdx.x*(blockSize)*2 + tid;                                                                                                                       
                    int gridSize = blockSize*gridDim.x*2;                                                                                                                         

                    while (i < n) {
                        //if (sdata[tid] < g_idata[i]){
                            sdata[tid] = g_idata[i];
                            maxid[tid] = i;
                        //}
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

        print(f'init completed on: {device}')

    def clean(self):
        del self.masked_taper

        del self.fftplan
        del self.ifftnP
        del self.fftnP
        del self.cp

        del self.stdV
        del self.ccc_map

        del self.volume_fft
        del self.volume

        del self.wedge
        del self.wedgePart
        del self.wedgeTex
        del self.rotatedWedge
        del self.wedgeAngles

        del self.referenceTex
        del self.rotatedRef
        del self.simulatedVolume

        del self.mask
        del self.maskTex
        del self.rotatedMask
        del self.mask_fft

        del self.normalize
        del self.sumMeanStdv
        del self.argmax

        del self.max_id
        del self.fast_sum
        del self.fast_sum_mean
        del self.fast_sum_stdv
        del self.zoomed
        del self.nblocks
        del self.num_threads



    # def meanUnderMask(self, volume, mask, p, xp, num_threads=1024):
    #
    #     nblocks = int(xp.ceil(volume.size/num_threads/2))
    #     self.fast_sum = xp.zeros((4*nblocks*2),dtype=xp.float32)
    #     self.reduce6((nblocks, 1,), (num_threads, 1, 1), (volume, mask, self.fast_sum, volume.size), shared_mem=4*num_threads)
    #     meanT = self.fast_sum.sum() / p
    #     #
    #     # print(meanT, (volume * mask).sum()/p)
    #     return meanT
    #
    # def stdUnderMask(self, volume, mask, meanV, p, xp, num_threads=512):
    #
    #     nblocks = int(xp.ceil(volume.size / num_threads/2))
    #     self.fast_sum = xp.zeros((4*nblocks*2), dtype=xp.float32)
    #     self.reduce7((nblocks, 1,), (num_threads, 1, 1), (volume, mask, self.fast_sum, volume.size), shared_mem=4*num_threads)
    #     dd = self.fast_sum.sum() / p
    #     var = dd - meanV*meanV
    #     # print('variance', abs(var - (volume*volume*mask).sum()/p - meanV*meanV), meanV)
    #
    #     return xp.sqrt(var)

    def normalizeVolume(self):
        self.fast_sum_stdv *= 0
        self.fast_sum_mean *= 0
        self.sumMeanStdv((self.nblocks, 1,), (self.num_threads, 1, 1), (self.simulatedVolume, self.mask, self.fast_sum_mean, self.fast_sum_stdv, self.simulatedVolume.size),
                         shared_mem=8 * self.num_threads)
        meanT = self.fast_sum_mean.sum() / self.p
        stdT = self.cp.sqrt(self.fast_sum_stdv.sum()/self.p - meanT*meanT)

        self.simulatedVolume = ((self.simulatedVolume - meanT)/stdT) * self.mask

    def cross_correlation(self):
        self.simulatedVolume = self.fftnP(self.simulatedVolume.astype(self.cp.complex64), plan=self.fftplan)
        self.ccc_map = ( self.cp.fft.ifftshift(self.ifftnP(self.volume_fft * self.simulatedVolume.conj(), plan=self.fftplan))).real
        self.ccc_map *= self.stdV

    def wedgeRotatedRef(self):
        self.simulatedVolume = self.ifftnP(self.fftnP(self.rotatedRef.astype(self.cp.complex64), plan=self.fftplan) * self.wedgePart,
                                  plan=self.fftplan).real.astype(self.cp.float32) * self.masked_taper

    def wedgeParticle(self):
        self.volume_fft = self.fftnP(self.volume.astype(self.cp.complex64),plan=self.fftplan) * self.wedgePart
        self.volume_fft = self.volume_fft.astype(self.cp.complex64)

    def subPixelMax3D(self, ignore_border=1, k=0.1, profile=False):

        from pytom.voltools import transform


        ox, oy, oz = self.ccc_map.shape
        ib = ignore_border
        cropped_volume = self.ccc_map[ib:ox - ib, ib:oy - ib, ib:oz - ib].astype(self.cp.float32)

        if profile:
            stream = self.cp.cuda.Stream.null
            t_start = stream.record()

        # x,y,z = self.cp.array(maxIndex(cropped_volume)) + ignore_border
        x, y, z = self.cp.array(self.cp.unravel_index(cropped_volume.argmax(), cropped_volume.shape)) + ignore_border

        dx, dy, dz = self.ccc_map.shape
        translation = [dx // 2 - x, dy // 2 - y, dz // 2 - z]

        if profile:
            t_end = stream.record()
            t_end.synchronize()

            time_took = self.cp.cuda.get_elapsed_time(t_start, t_end)
            print(f'initial find max time: \t{time_took:.3f}ms')
            t_start = stream.record()

        b = border = max(0, int(self.ccc_map.shape[0] // 2 - 4 / k))
        zx, zy, zz = self.ccc_map.shape
        out = self.ccc_map[b:zx - b, b:zy - b, b:zz - b].astype(self.cp.float32)
        self.zoomed *= 0
        transform(out, output=self.zoomed, scale=(k, k, k), device=self.device, translation=translation,
                  interpolation='filt_bspline')

        if profile:
            t_end = stream.record()
            t_end.synchronize()

            time_took = self.cp.cuda.get_elapsed_time(t_start, t_end)
            print(f'transform finished in \t{time_took:.3f}ms')
            t_start = stream.record()
        num_threads = 1024
        nblocks = int(self.cp.ceil(self.zoomed.size / num_threads / 2))
        self.max_id *= 0
        self.fast_sum *= 0
        self.argmax((nblocks, 1,), (num_threads, 1, 1), (self.zoomed, self.fast_sum, self.max_id, self.zoomed.size), shared_mem=8 * num_threads)
        x2, y2, z2 = self.cp.unravel_index(self.max_id[self.fast_sum.argmax()], self.zoomed.shape)

        peakValue = self.zoomed[x2][y2][z2]
        peakShift = [x + (x2 - self.zoomed.shape[0] // 2) * k - zx // 2,
                     y + (y2 - self.zoomed.shape[1] // 2) * k - zy // 2,
                     z + (z2 - self.zoomed.shape[2] // 2) * k - zz // 2]

        if profile:
            t_end = stream.record()
            t_end.synchronize()

            time_took = self.cp.cuda.get_elapsed_time(t_start, t_end)
            print(f'argmax finished in \t{time_took:.3f}ms')
            t_start = stream.record()
            print()

        return [peakValue, peakShift]

    def maxIndex(self, volume, num_threads=1024):
        nblocks = int(self.cp.ceil(volume.size / num_threads / 2))
        fast_sum = -1000000 * self.cp.ones((nblocks), dtype=self.cp.float32)
        max_id = self.cp.zeros((nblocks), dtype=self.cp.int32)
        self.argmax((nblocks, 1,), (num_threads, 1, 1), (volume, fast_sum, max_id, volume.size),
                    shared_mem=16 * num_threads)
        mm = min(max_id[fast_sum.argmax()], volume.size - 1)
        indices = self.cp.unravel_index(mm, volume.shape)
        return indices

    def subPixelMax33D(self, volume, k=.01, ignore_border=50, interpolation='filt_bspline', profile=True):
        """
        Function to find the highest point in a 3D array, with subpixel accuracy using cubic spline interpolation.

        @param inp: A 3D numpy/cupy array containing the data points.
        @type inp: numpy/cupy array 3D
        @param k: The interpolation factor used in the spline interpolation, k < 1 is zoomed in, k>1 zoom out.
        @type k: float
        @return: A list of maximal value in the interpolated volume and a list of  x position, the y position and
            the value.
        @returntype: list
        """

        from pytom.voltools import transform
        from pytom.tompy.io import write

        ox, oy, oz = volume.shape
        ib = ignore_border
        cropped_volume = volume[ib:ox - ib, ib:oy - ib, ib:oz - ib].astype(self.cp.float32)

        if profile:
            stream = self.cp.cuda.Stream.null
            t_start = stream.record()

        # x,y,z = self.cp.array(maxIndex(cropped_volume, Argmax)) + ignore_border
        x, y, z = self.cp.array(self.cp.unravel_index(cropped_volume.argmax(), cropped_volume.shape)) + ignore_border

        if profile:
            t_end = stream.record()
            t_end.synchronize()

            time_took = self.cp.cuda.get_elapsed_time(t_start, t_end)
            print(f'initial find max time: \t{time_took:.3f}ms')
            t_start = stream.record()

        b = border = max(0, int(volume.shape[0] // 2 - 4 / k))
        zx, zy, zz = volume.shape

        croppedForZoom = volume[b:zx - b, b:zy - b, b:zz - b]
        out = self.cp.array(croppedForZoom, dtype=self.cp.float32)

        dx, dy, dz = volume.shape
        translation = [dx // 2 - x, dy // 2 - y, dz // 2 - z]
        zoomed = self.cp.zeros_like(out, dtype=self.cp.float32)

        transform(out, output=zoomed, scale=(k, k, k), device=self.device, translation=translation,
                  interpolation=interpolation)
        if profile:
            t_end = stream.record()
            t_end.synchronize()

            time_took = self.cp.cuda.get_elapsed_time(t_start, t_end)
            print(f'transform finished in \t{time_took:.3f}ms')
            t_start = stream.record()

        # x2, y2, z2 = self.maxIndex(zoomed)
        x2,y2,z2 = self.cp.unravel_index(zoomed.argmax(), zoomed.shape)

        peakValue = zoomed[x2][y2][z2]
        peakShift = [x + (x2 - zoomed.shape[0] // 2) * k - volume.shape[0] // 2,
                     y + (y2 - zoomed.shape[1] // 2) * k - volume.shape[1] // 2,
                     z + (z2 - zoomed.shape[2] // 2) * k - volume.shape[2] // 2]
        if profile:
            t_end = stream.record()

        return [peakValue, peakShift]

    def calc_stdV(self):
        meanV = (self.cp.fft.fftshift(self.ifftnP(self.volume_fft * self.cp.conj(self.mask_fft), plan=self.fftplan))/self.p).real
        self.stdV = 1 / (self.stdVolUnderMaskPlanned(self.volume, meanV) * self.p)
        del meanV
    def meanVolUnderMaskPlanned(self, volume):
        """
        meanUnderMask: calculate the mean volume under the given mask (Both should have the same size)
        @param volume: input volume
        @type volume:  L{numpy.ndarray} L{cupy.ndarray}
        @param mask: mask
        @type mask:  L{numpy.ndarray} L{cupy.ndarray}
        @return: the calculated mean volume under mask
        @rtype:  L{numpy.ndarray} L{cupy.ndarray}
        @author: Gijs van der Schot
        """
        volume_fft = self.fftnP(volume.astype(self.cp.complex64), plan=self.fftplan)
        mask_fft = self.fftnP(self.rotatedMask.astype(self.cp.complex64), plan=self.fftplan)
        res = self.cp.fft.fftshift(self.ifftnP(volume_fft * self.cp.conj(mask_fft), plan=self.fftplan)) / self.p

        return res.real

    def stdVolUnderMaskPlanned(self, volume, meanV):
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
        var = self.meanVolUnderMaskPlanned(vol2) - meanV2
        var[var < 1E-09] = 1

        return self.cp.sqrt(var)

    def addParticleAndWedgeToSum(self, particle, bestPeak, centerCoordinates, rotation_order='rzxz'):
        from pytom.voltools import transform

        bestRotation = [-1 * bestPeak.getRotation()[1], -1 * bestPeak.getRotation()[2], -1 * bestPeak.getRotation()[0]]
        bestShift = list((self.cp.array(bestPeak.getShift().toVector()) * -1).get())

        # Add rotated particle to sum
        self.rotatedRef *= 0
        transform(particle, rotation=bestRotation, translation=bestShift, output=self.rotatedRef,
                        rotation_order=rotation_order, device=self.device, interpolation='filt_bspline')
        self.sumParticles += self.rotatedRef

        # Add rotated wedge to sum
        self.rotatedWedge *= 0
        self.wedgeTex.transform(rotation=bestRotation, center=centerCoordinates, rotation_order=rotation_order,
                                output=self.rotatedWedge)
        self.sumWeights += self.cp.fft.fftshift(self.rotatedWedge)[:,:,:self.rotatedWedge.shape[2]//2+1]

    def updateWedge(self, wedge, interpolation='filt_bspline'):
        from pytom.voltools import StaticVolume

        wedgeAngles = wedge.getWedgeAngle()
        if type(wedgeAngles) == float: wedgeAngles = [wedgeAngles, wedgeAngles]
        if type(self.wedgeAngles) == float: self.wedgeAngles = [self.wedgeAngles, self.wedgeAngles]

        if wedgeAngles[0] == self.wedgeAngles[0] and wedgeAngles[1] == self.wedgeAngles[1]:
            return

        else:
            self.wedgeAngles  = wedgeAngles
            self.wedge        = wedge
            wedgePart         = wedge.returnWedgeVolume(*self.volume.shape, humanUnderstandable=True).astype(self.cp.float32).get()
            self.rotatedWedge = self.cp.array(wedgePart, dtype=self.cp.float32)
            self.wedgePart    = self.cp.fft.fftshift(wedgePart)
            self.wedgeTex     = StaticVolume(self.rotatedWedge.copy(), device=self.device, interpolation=interpolation)

def rotate3d(data, phi=0, psi=0, the=0, center=None, order=1, output=None):
    """Rotate a 3D data using ZXZ convention (phi: z1, the: x, psi: z2).

    @param data: data to be rotated.
    @param phi: 1st rotate around Z axis, in degree.
    @param psi: 3rd rotate around Z axis, in degree.
    @param the: 2nd rotate around X axis, in degree.
    @param center: rotation center.

    @return: the data after rotation.
    """
    from cupyx.scipy.ndimage import map_coordinates
    import cupy as cp

    # Figure out the rotation center
    if center is None:
        cx = data.shape[0] / 2
        cy = data.shape[1] / 2
        cz = data.shape[2] / 2
    else:
        assert len(center) == 3
        (cx, cy, cz) = center

    # Transfer the angle to Euclidean
    phi = -float(phi) * cp.pi / 180.0
    the = -float(the) * cp.pi / 180.0
    psi = -float(psi) * cp.pi / 180.0
    sin_alpha = cp.sin(phi)
    cos_alpha = cp.cos(phi)
    sin_beta = cp.sin(the)
    cos_beta = cp.cos(the)
    sin_gamma = cp.sin(psi)
    cos_gamma = cp.cos(psi)



    # Calculate inverse rotation matrix
    Inv_R = cp.zeros((3, 3), dtype=cp.float32)

    Inv_R[0][0] = cos_alpha * cos_gamma - cos_beta * sin_alpha * sin_gamma
    Inv_R[0][1] = -cos_alpha * sin_gamma - cos_beta * sin_alpha * cos_gamma
    Inv_R[0][2] = sin_beta * sin_alpha

    Inv_R[1][0] = sin_alpha * cos_gamma + cos_beta * cos_alpha * sin_gamma
    Inv_R[1][1] = -sin_alpha * sin_gamma + cos_beta * cos_alpha * cos_gamma
    Inv_R[1][2] = -sin_beta * cos_alpha

    Inv_R[2][0] = sin_beta * sin_gamma
    Inv_R[2][1] = sin_beta * cos_gamma
    Inv_R[2][2] = cos_beta



    grid = cp.mgrid[-cx:data.shape[0]-cx, -cy:data.shape[1]-cy, -cz:data.shape[2]-cz]
    temp = grid.reshape((3, grid.size // 3))
    temp = cp.dot(Inv_R, temp)
    grid = cp.reshape(temp, grid.shape)
    grid[0] += cx
    grid[1] += cy
    grid[2] += cz
    dataout = cp.zeros_like(data)

    # Interpolation

    dataout = map_coordinates(data, grid.reshape(len(grid), -1), order=order).reshape(grid.shape[1:])

    return dataout

class GLocalAlignmentGPU(threading.Thread):
    def __init__(self, jobid, deviceid, input):
        threading.Thread.__init__(self)

        import numpy as np
        from pytom.tompy.tools import paste_in_center
        from pytom.tompy.transform import rotate3d
        from pytom_volume import rotateSpline as rotate
        import cupy as cp

        cp.cuda.Device(deviceid).use()

        from cupy import sqrt, float32
        from cupy.fft import fftshift, rfftn, irfftn, ifftn, fftn
        import voltools as vt
        from pytom.tompy.io import write
        from cupyx.scipy.ndimage import map_coordinates


        if 1:
            from cupyx.scipy.fftpack.fft import get_fft_plan
            from cupyx.scipy.fftpack.fft import fftn as fftnP
            from cupyx.scipy.fftpack.fft import ifftn as ifftnP
            self.fftnP = fftnP
            self.ifftnP = ifftnP
        else:
            get_fft_plan = None

        self.cp = cp
        self.map_coordinates = map_coordinates
        self.Device = cp.cuda.Device
        self.jobid = jobid
        self.deviceid = deviceid
        self.active = True
        self.input = input

        self.fftn = fftn
        self.ifftn = ifftn
        self.fftshift = fftshift
        self.rfftn = rfftn
        self.irfftn = irfftn
        self.sqrt = sqrt
        self.float32 = float32


        self.update_scores_angles = cp.RawKernel(r"""
        extern "C"  __global__ void update_scores_angles(float *scores, float *angles, float *ccc_map, float angleId, int num_elements, int dimx)
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
           __syncthreads();

        }
        """, 'update_scores_angles')


        self.updateResFromIdx = cp.ElementwiseKernel(
            'float32 scores, float32 angles, float32 ccc_map, float32 angleId',
            'float32 out, float32 out2',
            'if (scores < ccc_map) {out = ccc_map; out2 = angleId;}',
            'updateResFromIdx')

        self.plan = TemplateMatchingPlan(input[0], input[1], input[2], input[3], cp, vt, self.calc_stdV, self.pad, get_fft_plan, deviceid)

        print("Initialized job_{:03d} on device {:d}".format(self.jobid, self.deviceid))

    def run(self):
        print("RUN")
        if 1:
            self.glocal_alignment_gpu(self.input[4], self.input[5])
            self.completed = True
        else:
            self.completed = False
        self.active = False

    def glocal_alignment_gpu(self, angle_list, dims, isSphere=True, verbose=True):

        self.Device(self.deviceid).use()

        import pytom.voltools as vt
        from pytom.tompy.io import write
        sx, sy, sz = self.plan.template.shape
        cx, cy, cz = sx // 2, sy // 2, sz // 2
        mx, my, mz = sx %  2, sy %  2, sz %  2

        for angleId, angles in enumerate(angle_list):

            # Rotate
            #self.plan.template = self.rotate3d(self.plan.templateOrig, phi=phi,the=the,psi=psi)
            self.plan.texture.transform(rotation=(angles[0], angles[2], angles[1]), rotation_order='rzxz', output=self.plan.template)

            # Add wedge
            #print(self.rfftn(self.plan.template).shape)
            self.plan.template = self.irfftn(self.rfftn(self.plan.template) * self.plan.wedge)

            # Normalize template
            meanT = self.meanUnderMask(self.plan.template, self.plan.mask, p=self.plan.p)
            stdT = self.stdUnderMask(self.plan.template, self.plan.mask, meanT, p=self.plan.p)
            self.plan.template = ((self.plan.template - meanT) / stdT) * self.plan.mask

            # write('template_gpu.em', self.plan.templatePadded)

            # Cross-correlate and normalize by stdV
            self.plan.ccc_map = self.normalized_cross_correlation(self.plan.volume_fft2, self.plan.template, self.plan.stdV, self.plan.p, plan=self.plan.fftplan)

            # Update the scores and angles
            self.updateResFromIdx(self.plan.scores, self.plan.angles, self.plan.ccc_map, angleId, self.plan.scores, self.plan.angles)

            #self.cp.cuda.stream.get_current_stream().synchronize()

    def is_alive(self):
        return self.active

    def calc_stdV(self, plan):
        stdV = self.meanVolUnderMask2(plan.volume**2, plan) - self.meanVolUnderMask2(plan.volume, plan)**2
        stdV[stdV < self.float32(1e-09)] = 1
        plan.stdV = self.sqrt(stdV)

    def meanVolUnderMask2(self, volume, plan):
        res = self.fftshift(self.ifftn(self.fftn(volume) * self.fftn(plan.maskPadded).conj())) / plan.mask.sum()
        return res.real

    def meanUnderMask(self,volume, mask=None, p=1, gpu=False):
        """
        meanValueUnderMask: Determines the mean value under a mask
        @param volume: The volume
        @type volume:  L{pytom_volume.vol}
        @param mask:  The mask
        @type mask:  L{pytom_volume.vol}
        @param p: precomputed number of voxels in mask
        @type p: float
        @return: A value (scalar)
        @rtype: single
        @change: support None as mask, FF 08.07.2014
        """
        return (volume * mask).sum() / p

    def stdUnderMask(self, volume, mask, meanValue, p=None, gpu=False):
        """
        stdValueUnderMask: Determines the std value under a mask

        @param volume: input volume
        @type volume:  L{pytom_volume.vol}
        @param mask: mask
        @type mask:  L{pytom_volume.vol}
        @param p: non zero value numbers in the mask
        @type p: L{float} or L{int}
        @return: A value
        @rtype: L{float}
        @change: support None as mask, FF 08.07.2014
        """
        return self.sqrt(self.meanUnderMask(volume**2, mask, p=p) - meanValue**2)

    def normalized_cross_correlation(self, volume_fft, template, norm, p=1, plan=None):
        return self.cp.real(self.fftshift(self.ifftnP(volume_fft * self.cp.conj(self.fftnP(template, plan=plan)), plan=plan)) / (p * norm))
        #print(volume_fft.shape)

        print('cc')
        cc = self.cp.conj(self.cp.fft.fftn(template))
        print('ccc')
        aa = self.cp.fft.fftn(template)
        res = self.cp.abs(self.cp.fft.fftshift(self.cp.fft.ifftn(aa * cc )) )
        #res = self.cp.abs(self.cp.fft.fftshift(self.cp.fft.ifftn(volume_fft * xp.fft.fftn(template).conj() )) / (p * norm))
        return res


    def pad(self, volume, out, sPad, sOrg):
        SX, SY, SZ = sPad
        sx, sy, sz = sOrg
        print(sPad, sOrg)
        out[SX // 2 - sx // 2:SX // 2 + sx // 2 + sx % 2, SY // 2 - sy // 2:SY // 2 + sy // 2 + sy % 2,
            SZ // 2 - sz // 2:SZ // 2 + sz // 2 + sz % 2] = volume

    def rotate3d(self, data, phi=0, psi=0, the=0, center=None, order=1, output=None):
        """Rotate a 3D data using ZXZ convention (phi: z1, the: x, psi: z2).

        @param data: data to be rotated.
        @param phi: 1st rotate around Z axis, in degree.
        @param psi: 3rd rotate around Z axis, in degree.
        @param the: 2nd rotate around X axis, in degree.
        @param center: rotation center.

        @return: the data after rotation.
        """
        # Figure out the rotation center
        if center is None:
            cx = data.shape[0] / 2
            cy = data.shape[1] / 2
            cz = data.shape[2] / 2
        else:
            assert len(center) == 3
            (cx, cy, cz) = center

        # Transfer the angle to Euclidean
        phi = -float(phi) * self.cp.pi / 180.0
        the = -float(the) * self.cp.pi / 180.0
        psi = -float(psi) * self.cp.pi / 180.0
        sin_alpha = self.cp.sin(phi)
        cos_alpha = self.cp.cos(phi)
        sin_beta = self.cp.sin(the)
        cos_beta = self.cp.cos(the)
        sin_gamma = self.cp.sin(psi)
        cos_gamma = self.cp.cos(psi)



        # Calculate inverse rotation matrix
        Inv_R = self.cp.zeros((3, 3), dtype=self.cp.float32)

        Inv_R[0][0] = cos_alpha * cos_gamma - cos_beta * sin_alpha * sin_gamma
        Inv_R[0][1] = -cos_alpha * sin_gamma - cos_beta * sin_alpha * cos_gamma
        Inv_R[0][2] = sin_beta * sin_alpha

        Inv_R[1][0] = sin_alpha * cos_gamma + cos_beta * cos_alpha * sin_gamma
        Inv_R[1][1] = -sin_alpha * sin_gamma + cos_beta * cos_alpha * cos_gamma
        Inv_R[1][2] = -sin_beta * cos_alpha

        Inv_R[2][0] = sin_beta * sin_gamma
        Inv_R[2][1] = sin_beta * cos_gamma
        Inv_R[2][2] = cos_beta



        grid = self.cp.mgrid[-cx:data.shape[0]-cx, -cy:data.shape[1]-cy, -cz:data.shape[2]-cz]
        temp = grid.reshape((3, grid.size // 3))
        temp = self.cp.dot(Inv_R, temp)
        grid = self.cp.reshape(temp, grid.shape)
        grid[0] += cx
        grid[1] += cy
        grid[2] += cz

        # Interpolation
        #self.map_coordinates(data, grid, output=output)
        return self.map_coordinates(data, grid.reshape(len(grid), -1), order=order).reshape(grid.shape[1:])
