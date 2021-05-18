import threading
import importlib
from pytom.tompy.io import read, write
from pytom.gpu.initialize import device, xp

class TemplateMatchingPlan():
    def __init__(self, volume, template, mask, wedge, cp, vt, calc_stdV, pad, get_fft_plan, deviceid):
        from pytom.tompy.io import read
        from pytom.basic.files import read as readC
        import pytom.voltools as vt

        cp.cuda.Device(deviceid).use()

        self.volume = cp.asarray(volume,dtype=cp.float32,order='C')

        self.sPad = volume.shape

        self.mask = cp.asarray(mask, dtype=cp.float32, order='C')
        self.maskPadded = cp.zeros_like(self.volume).astype(cp.float32)
        self.sOrg = mask.shape
        pad(self.mask, self.maskPadded, self.sPad, self.sOrg)

        self.templateOrig = cp.asarray(template, dtype=cp.float32,order='C')

        write('te.em', template)
        self.templateVol = readC('te.em')
        self.where = 'gpu'
        self.template = cp.asarray(template, dtype=cp.float32,order='C')
        self.texture = vt.StaticVolume(self.template, interpolation='filt_bspline', device=f'gpu:{deviceid}')
        self.templatePadded = cp.zeros_like(self.volume)

        self.wedge = cp.asarray(wedge,order='C')
        #self.volume_fft = cp.fft.rfftn(self.volume) #* cp.array(wedgeVolume,dtype=cp.float32)
        calc_stdV(self)
        cp.cuda.stream.get_current_stream().synchronize()
        # del self.volume_fft

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

        cp.cuda.stream.get_current_stream().synchronize()


class TemplateMatchingGPU(threading.Thread):
    def __init__(self, jobid, deviceid, input):
        threading.Thread.__init__(self)

        import cupy as cp
        cp.cuda.Device(deviceid).use()

        from cupy import sqrt, float32
        from cupy.fft import fftshift, rfftn, irfftn, ifftn, fftn
        import pytom.voltools as vt
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
        self.Device(self.deviceid).use()
        if 1:
            self.template_matching_gpu(self.input[4], self.input[5])
            self.completed = True
        else:
            self.completed = False
        self.active = False

    def template_matching_gpu(self, angle_list, dims, isSphere=True, verbose=True):
        from pytom_numpy import vol2npy
        from pytom_volume import rotateSpline as rotate, vol
        import time
        ref = vol(self.plan.templateVol.sizeX(), self.plan.templateVol.sizeY(), self.plan.templateVol.sizeZ())
        self.Device(self.deviceid).use()

        sx, sy, sz = self.plan.template.shape
        SX, SY, SZ = self.plan.templatePadded.shape
        CX, CY, CZ = SX // 2, SY // 2, SZ // 2
        cx, cy, cz = sx // 2, sy // 2, sz // 2
        mx, my, mz = sx %  2, sy %  2, sz %  2

        for angleId, angles in enumerate(angle_list):
            # Rotate
            #self.plan.template = self.rotate3d(self.plan.templateOrig, phi=phi,the=the,psi=psi)

            if self.plan.where  == 'gpu':
                self.plan.texture.transform(rotation=(angles[0], angles[2], angles[1]), rotation_order='rzxz', output=self.plan.template)
                #self.plan.template *= self.plan.mask
            else:
                #ref.setAll(0.)
                rotate(self.plan.templateVol, ref, angles[0], angles[1], angles[2])
                self.plan.template = xp.array(vol2npy(ref).copy(),dtype=xp.float32)

            # Add wedge
            self.plan.template = self.irfftn(self.rfftn(self.plan.template) * self.plan.wedge, s=self.plan.template.shape)

            # Normalize template
            meanT = self.meanUnderMask(self.plan.template, self.plan.mask, p=self.plan.p)
            stdT = self.stdUnderMask(self.plan.template, self.plan.mask, meanT, p=self.plan.p)

            self.plan.template = ((self.plan.template - meanT) / stdT) * self.plan.mask

            # Paste in center
            self.plan.templatePadded[CX-cx:CX+cx+mx, CY-cy:CY+cy+my,CZ-cz:CZ+cz+mz] = self.plan.template

            # write('volume_gpu.em', self.ifftnP(self.plan.volume_fft2, plan=self.plan.fftplan).real)
            # write('template_gpu.em', self.plan.templatePadded)
            # write('m_gpu.em', self.plan.maskPadded)
            # write('stdV_gpu.em', self.plan.stdV)

            # Cross-correlate and normalize by stdV
            self.plan.ccc_map = self.normalized_cross_correlation(self.plan.volume_fft2, self.plan.templatePadded, self.plan.stdV, self.plan.p, plan=self.plan.fftplan)

            # Update the scores and angles
            self.updateResFromIdx(self.plan.scores, self.plan.angles, self.plan.ccc_map, angleId, self.plan.scores, self.plan.angles)
            #self.cp.cuda.stream.get_current_stream().synchronize()

            if angleId % 100 == 0: print(angleId, angles, self.plan.scores.max(), angle_list[int(self.plan.angles[self.cp.unravel_index(self.plan.scores.argmax(), self.plan.scores.shape)].get())])

    def is_alive(self):
        return self.active

    def calc_stdV(self, plan):
        self.Device(self.deviceid).use()
        stdV = self.meanVolUnderMask2(plan.volume**2, plan) - self.meanVolUnderMask2(plan.volume, plan)**2
        stdV[stdV < self.float32(1e-09)] = 1
        plan.stdV = self.sqrt(stdV)

    def meanVolUnderMask2(self, volume, plan):
        self.Device(self.deviceid).use()
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
        return self.fftshift(self.ifftnP(volume_fft * self.fftnP(template, plan=plan).conj(), plan=plan)).real / (p * norm)
        # #print(volume_fft.shape)
        #
        # print('cc')
        # cc = xp.conj(xp.fft.fftn(template))
        # print('ccc')
        # aa = xp.fft.fftn(template)
        # res = xp.abs(xp.fft.fftshift(xp.fft.ifftn(aa * cc )) )
        # #res = xp.abs(xp.fft.fftshift(xp.fft.ifftn(volume_fft * xp.fft.fftn(template).conj() )) / (p * norm))
        # return res


    def pad(self, volume, out, sPad, sOrg):
        SX, SY, SZ = sPad
        sx, sy, sz = sOrg
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
        from pytom.tompy.io import read_size, write
        from pytom.gpu.kernels import argmax_text, meanStdv_text
        from pytom.tompy.filter import bandpass

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
        self.wedgePart     = cp.fft.fftshift(wedgePart).astype(cp.float32)
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
        self.num_threads   = 1024
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
        self.normalize     = cp.ElementwiseKernel('T ref, T mask, raw T mean, raw T stdV ', 'T z', 'z = ((ref - mean[i*0]) / stdV[i*0]) * mask', 'norm2')
        self.sumMeanStdv   = cp.RawKernel(meanStdv_text, 'sumMeanStdv')
        self.argmax        = self.cp.RawKernel(argmax_text, 'argmax')

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

    def subPixelMax3D(self, ignore_border=1, k=0.1, profile=False, check=True):

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
        if b> 0:
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
            num_threads = self.num_threads
            nblocks = int(self.cp.ceil(self.zoomed.size / num_threads / 2 ))
            self.max_id *= 0
            self.fast_sum *= 0

            self.argmax((nblocks, 1,), (num_threads, 1, 1), (self.zoomed, self.fast_sum, self.max_id, self.zoomed.size), shared_mem=8 * num_threads)

            # if check:
            #     zzz = self.zoomed.argmax()
            #     a = self.zoomed.max()
            #
            #     if self.max_id[self.fast_sum.argmax()] != zzz:
            #         print(a, self.fast_sum.max())
            #         print(self.cp.unravel_index(self.max_id[self.fast_sum.argmax()], self.zoomed.shape))
            #         print(self.cp.unravel_index(zzz, self.zoomed.shape))

            x2, y2, z2 = self.cp.unravel_index(self.max_id[self.fast_sum.argmax()], self.zoomed.shape)


        else:
            out = self.ccc_map.copy().astype(self.cp.float32)
            self.zoomed *= 0
            transform(out, output=self.zoomed, scale=(k, k, k), device=self.device, translation=translation,
                      interpolation='filt_bspline')
            zx, zy, zz = self.ccc_map.shape
            x2, y2, z2 = self.cp.unravel_index(self.zoomed.argmax(), self.zoomed.shape)

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
        import pytom.voltools as vt
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


class CCCPlan():
    def __init__(self, particleList, maskname, freq, cp=xp, device='cpu', interpolation='filt_bspline',
                 max_num_part=65, profile=True, binning=1):
        id = int(device.split(":")[1])

        cp.cuda.Device(id).use()
        self.cp = cp
        self.device = device
        print(f'start init on: {device}')
        from pytom.voltools import transform, StaticVolume
        from pytom.tompy.tools import taper_edges
        from pytom.tompy.transform import fourier_reduced2full, resize
        from cupyx.scipy.fftpack.fft import get_fft_plan
        from cupyx.scipy.fftpack.fft import fftn as fftnP
        from cupyx.scipy.fftpack.fft import ifftn as ifftnP
        from pytom.tompy.io import read_size, read, write
        from pytom.tompy.tools import create_sphere

        maskFull = read(maskname)

        if binning != 1:
            mask = resize(maskFull, 1. / binning)[0] # This is not a good idea
            mask = maskFull[::binning, ::binning, ::binning]
        else:
            mask = maskFull

        wedgeVol = particleList[0].getWedge().convert2numpy().getWedgeObject().returnWedgeVolume(*mask.shape)
        wedge = fourier_reduced2full(wedgeVol, isodd=mask.shape[-1] % 2)
        lpf = xp.fft.fftshift(create_sphere(mask.shape, freq))

        print(mask.shape)
        self.volume = cp.zeros(mask.shape, dtype=cp.float32)
        self.vg = cp.zeros_like(self.volume, dtype=cp.float32)
        self.vf = cp.zeros_like(self.volume, dtype=cp.float32)

        self.binning = binning
        self.max_particles_on_card = max_num_part
        self.particleList = particleList
        self.formatting = ['', ] * self.max_particles_on_card * 2

        # General functions
        self.transform = transform
        self.ifftnP = ifftnP
        self.fftnP = fftnP
        self.fftplan = get_fft_plan(self.volume.astype(cp.complex64))
        self.ifftplan = get_fft_plan(self.volume.astype(cp.complex64))
        self.fftplan2 = get_fft_plan(maskFull.astype(cp.complex64))
        self.ifftplan2 = get_fft_plan(self.volume.astype(cp.complex64))

        self.read = read
        self.write = write
        self.xp = cp
        self.profile = profile
        self.stream = cp.cuda.Stream.null
        self.staticVolume = StaticVolume
        self.device = device
        self.interpolation = interpolation
        self.read_time = 0

        # Mask info
        self.mask = self.cp.array(mask, dtype=xp.float32)
        self.p = mask.sum()

        # Wedge
        self.wedge = StaticVolume(cp.fft.fftshift(wedge), device=device, interpolation=interpolation)
        self.wedgeZero = self.cp.array(wedge, dtype=self.cp.float32)
        self.rot_wf = self.cp.zeros_like(self.volume)
        self.rot_wg = self.cp.zeros_like(self.volume)
        self.cwf, self.cwg = [-9999, -9999, -9999], [-9999, -9999, -9999]

        # Low Pass Filter
        self.lpf = self.cp.array(lpf, dtype=self.cp.float32)

        # Allocate arrays used for normalisation of reference
        self.num_threads = 1024
        rr = 1
        self.nblocks = int(self.cp.ceil(self.volume.size / self.num_threads / 2))
        self.fast_sum_mean = self.cp.zeros((self.nblocks * rr), dtype=self.cp.float32)
        self.fast_sum_stdv = self.cp.zeros((self.nblocks * rr), dtype=self.cp.float32)

        # Kernels

        self.normalize = cp.RawKernel(r'''
extern "C" __global__
void normalize( float* volume, const float* mask, float mean, float stdV, int n ){
     int tid = blockDim.x * blockIdx.x + threadIdx.x;
     
     if (tid < n) { volume[tid] = ((volume[tid] - mean)/ stdV) * mask[tid];}
}
''', 'normalize')

        self.sumMeanStdv = cp.RawKernel(r'''
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
            mean[tid] += g_idata[i] * mask[i] ; 
            stdv[tid] += g_idata[i] * g_idata[i] * mask[i];
            if (i + blockSize < n){
                mean[tid] += g_idata[i + blockSize] * mask[i+blockSize];
                stdv[tid] += g_idata[i + blockSize] * g_idata[i + blockSize] * mask[i+blockSize];}
            i += gridSize;}
         __syncthreads();                                                                                                                                                       

        if (blockSize >= 1024){ if (tid < 512) { mean[tid] += mean[tid + 512]; stdv[tid] += stdv[tid+512];} __syncthreads(); }                                                                                                    
        if (blockSize >= 512) { if (tid < 256) { mean[tid] += mean[tid + 256]; stdv[tid] += stdv[tid+256];} __syncthreads(); }                                                                          
        if (blockSize >= 256) { if (tid < 128) { mean[tid] += mean[tid + 128]; stdv[tid] += stdv[tid+128];} __syncthreads(); }                                                                          
        if (blockSize >= 128) { if (tid <  64) { mean[tid] += mean[tid +  64]; stdv[tid] += stdv[tid+ 64];} __syncthreads(); }                                                                          
        if (tid < 32){ warpReduce(mean, stdv, tid, blockSize);}                                                                                                                                                                      
        if (tid == 0) {g_mean[blockIdx.x] = mean[0]; g_stdv[blockIdx.x] = stdv[0];} 
       __syncthreads(); 

    }''', 'sumMeanStdv')

        self.sumKernel = cp.RawKernel(r'''
            __device__ void warpReduceSum(volatile float* sdata, int tid, int blockSize) {
                if (blockSize >= 64) {sdata[tid] += sdata[tid + 32];}
                if (blockSize >= 32) {sdata[tid] += sdata[tid + 16];}
                if (blockSize >= 16) {sdata[tid] += sdata[tid +  8];}
                if (blockSize >=  8) {sdata[tid] += sdata[tid +  4];}
                if (blockSize >=  4) {sdata[tid] += sdata[tid +  2];}
                if (blockSize >=  2) {sdata[tid] += sdata[tid +  1];} }



            extern "C" __global__ 
            void sumKernel(float *g_idata, float *g_i2data, float *g_mdata, float *g_mean,  int n) {

                __shared__ float mean[1024];
                int blockSize = blockDim.x; 
                unsigned int tid = threadIdx.x;                                                                                                                                        
                int i = blockIdx.x*(blockSize)*2 + tid;                                                                                                                       
                int gridSize = blockSize*gridDim.x*2;                                                                                                                         

                mean[tid] = 0.;                                                                                                                                                       

                while (i < n) {
                    mean[tid] += g_idata[i] * g_i2data[i] * g_mdata[i]; 
                    if (i + blockSize < n){mean[tid] += g_idata[i + blockSize] * g_i2data[i + blockSize] * g_mdata[i + blockSize];}
                    i += gridSize;}
                 __syncthreads();                                                                                                                                                       

                if (blockSize >= 1024){ if (tid < 512) { mean[tid] += mean[tid + 512];} __syncthreads(); }                                                                                                    
                if (blockSize >= 512) { if (tid < 256) { mean[tid] += mean[tid + 256];} __syncthreads(); }                                                                          
                if (blockSize >= 256) { if (tid < 128) { mean[tid] += mean[tid + 128];} __syncthreads(); }                                                                          
                if (blockSize >= 128) { if (tid <  64) { mean[tid] += mean[tid +  64];} __syncthreads(); }                                                                          
                if (tid < 32){ warpReduceSum(mean, tid, blockSize);}                                                                                                                                                                      
                if (tid == 0) {g_mean[blockIdx.x] = mean[0];} 
               __syncthreads(); 

            }''', 'sumKernel')


        # Results
        self.results = cp.zeros((len(self.particleList), len(self.particleList)), dtype=cp.float32)
        self.currently_loaded = cp.ones((self.max_particles_on_card * 2), dtype=cp.int32) * -1
        self.needed = cp.ones_like(self.currently_loaded) * -1

        print(f'init completed on: {device}')

    def normalizeVolume(self, volume):
        self.fast_sum_stdv *= 0.
        self.fast_sum_mean *= 0.
        self.sumMeanStdv((self.nblocks, 1,), (self.num_threads, 1, 1),
                         (volume, self.mask, self.fast_sum_mean, self.fast_sum_stdv, self.simulatedVolume.size),
                         shared_mem=8 * self.num_threads)
        meanT = self.fast_sum_mean.sum() / self.p
        stdT = self.cp.sqrt(self.fast_sum_stdv.sum() / self.p - meanT * meanT)

        return ((self.simulatedVolume - meanT) / stdT) * self.mask

    def ccc_go(self, pairs):
        from pytom.tompy.filter import create_wedge
        from pytom.tompy.transform import fourier_reduced2full as r2f
        ttt = 0
        vg, vf = 0, 0
        for pair in pairs:
            if self.profile:
                tall = self.stream.record()
            tt = 0
            # print(int(self.currently_loaded[pair[0]]), int(self.currently_loaded[pair[1]]))
            g = self.particleList[int(self.currently_loaded[pair[0]])]
            f = self.particleList[int(self.currently_loaded[pair[1]])]

            shift_f = f.getShift()
            shift_f.scale(-1)
            shift_g = g.getShift()
            shift_g.scale(-1)

            wf_rotation_zzx = f.getRotation().invert().toVector()
            wf_rotation = wf = [wf_rotation_zzx[0], wf_rotation_zzx[2], wf_rotation_zzx[1]]
            wg_rotation_zzx = g.getRotation().invert().toVector()
            wg_rotation = wg = [wg_rotation_zzx[0], wg_rotation_zzx[2], wg_rotation_zzx[1]]

            if self.cwf == f.getFilename():
                dowf = False
            else:
                dowf = True
                self.cwf = f.getFilename()

            if self.cwg == g.getFilename():
                dowg = False
            else:
                dowg = True
                self.cwg = g.getFilename()
            if self.profile:
                t_start = self.stream.record()
            # self.transform(self.particles[int(pair[0]),:,:,:], output=self.vf, interpolation='filt_bspline', rotation=wf_rotation, device=self.device)
            # self.transform(self.particles[int(pair[1]),:,:,:], output=self.vg, interpolation='filt_bspline', rotation=wg_rotation, device=self.device)

            if dowg:
                self.formatting[pair[0]].transform(output=self.vg, rotation=wg_rotation, rotation_order='rzxz',
                                                   translation2=shift_g.toVector())
                self.fvg = self.fftnP(self.vg.astype(self.cp.complex64), plan=self.fftplan)

                self.wedge.transform(output=self.rot_wg, rotation=wg_rotation, rotation_order='rzxz')
                self.rot_wg = self.xp.fft.fftshift(self.rot_wg)

            if dowf:
                self.formatting[pair[1]].transform(output=self.vf, rotation=wf_rotation, rotation_order='rzxz',
                                                   translation2=shift_f.toVector())
                self.fvf = self.fftnP(self.vf.astype(self.cp.complex64), plan=self.fftplan)

                self.wedge.transform(output=self.rot_wf, rotation=wf_rotation, rotation_order='rzxz')
                self.rot_wf = self.xp.fft.fftshift(self.rot_wf)

            #
            # self.rot_wf *= 0.
            # self.rot_wg *= 0.
            if self.profile:
                t_end = self.stream.record()
                t_end.synchronize()

                time_took = self.cp.cuda.get_elapsed_time(t_start, t_end)
                print(f'rotation parts: \t{time_took:.3f}ms')
                t_start = self.stream.record()

            # self.rot_wg = self.xp.fft.fftshift(r2f(create_wedge(30,30,30,60,60,60, rotation=wg_rotation)))
            # self.rot_wf = self.xp.fft.fftshift(r2f(create_wedge(30,30,30,60,60,60, rotation=wf_rotation)))
            #

            # self.write('wg.mrc', self.xp.fft.fftshift(self.rot_wg))
            # self.write('wf.mrc', self.rot_wf)
            # raise Exception('Test')

            # write(f'GPU/vol_{self.currently_loaded[pair[0]]}_{self.currently_loaded[pair[1]]}.mrc', self.ifftnP(self.fftnP(self.vf.astype(self.cp.complex64), plan=self.fftplan) * self.xp.fft.fftshift(self.rot_wg*self.rot_wf), plan=self.fftplan).real.astype(self.xp.float32))
            # # write(f'GPU/vol_{self.currently_loaded[pair[1]]}_{self.currently_loaded[pair[0]]}.mrc', self.ifftnP(self.fftnP(self.vg.astype(self.cp.complex64), plan=self.fftplan) * self.xp.fft.fftshift(self.rot_wf*self.rot_wg), plan=self.fftplan).real.astype(self.xp.float32))
            # a = self.ifftnP(self.fftnP(self.vf.astype(self.cp.complex64), plan=self.fftplan) * self.xp.fft.fftshift(self.rot_wg), plan=self.fftplan).real.astype(self.xp.float32)
            # b = self.ifftnP(self.fftnP(self.vg.astype(self.cp.complex64), plan=self.fftplan) * self.xp.fft.fftshift(self.rot_wf), plan=self.fftplan).real.astype(self.xp.float32)
            # #
            # #
            # write(f'GPU/vol_{self.currently_loaded[pair[0]]}_{self.currently_loaded[pair[1]]}.mrc', a)
            # write(f'GPU/vol_{self.currently_loaded[pair[1]]}_{self.currently_loaded[pair[0]]}.mrc', b)

            w = self.rot_wg * self.rot_wf
            self.wvf = self.ifftnP(self.fvf * w, plan=self.ifftplan).real.astype(self.xp.float32)
            self.wvg = self.ifftnP(self.fvg * w, plan=self.ifftplan2).real.astype(self.xp.float32)

            # write(f'GPU/vol_{self.currently_loaded[pair[0]]}_{self.currently_loaded[pair[1]]}.mrc', self.wvf)
            # write(f'GPU/vol_{self.currently_loaded[pair[1]]}_{self.currently_loaded[pair[0]]}.mrc', self.wvg)

            if self.profile:
                t_end = self.stream.record()
                t_end.synchronize()

                time_took = self.cp.cuda.get_elapsed_time(t_start, t_end)
                print(f'apply wedge: \t{time_took:.3f}ms')
                t_start = self.stream.record()

            score = self.nxcc(self.wvf, self.wvg)

            if self.profile:
                t_end = self.stream.record()
                t_end.synchronize()

                time_took = self.cp.cuda.get_elapsed_time(t_start, t_end)
                print(f'nxcc calculate: \t{time_took:.3f}ms')

            self.results[self.currently_loaded[pair[0]], self.currently_loaded[pair[1]]] += self.xp.float32(score)
            if self.profile:
                t_end = self.stream.record()
                t_end.synchronize()

                time_took = self.cp.cuda.get_elapsed_time(tall, t_end)
                print(f'total time: \t{time_took:.3f}ms\n')

    def nxcc(self, volume, template):
        """
        nxcc: Calculates the normalized cross correlation coefficient in real space
        @param volume: A volume
        @type volume:  L{pytom_volume.vol}
        @param template: A template that is searched in volume. Must be of same size as volume.
        @type template:  L{pytom_volume.vol}
        @param mask: mask to constrain correlation
        @type mask: L{pytom_volume.vol}
        @param volumeIsNormalized: speed up if volume is already normalized
        @type volumeIsNormalized: L{bool}
        @return: A value between -1 and 1
        @raise exception: Raises a runtime error if volume and template have a different size.
        @author: Thomas Hrabe
        @change: flag for pre-normalized volume, FF
        """

        from pytom.tompy.macros import volumesSameSize
        from pytom.tompy.io import write

        if not volumesSameSize(volume, template):
            raise RuntimeError('Volume and template must have same size!')


        volume = self.normaliseUnderMask(volume)
        template = self.normaliseUnderMask(template)


        self.fast_sum_mean *= 0.
        self.sumKernel((self.nblocks, 1, 1), (self.num_threads,),
                       (volume, template, self.mask, self.fast_sum_mean, volume.size),
                       shared_mem=16 * self.num_threads)

        ncc = self.fast_sum_mean.sum()
        ncc /= self.p

        return float(ncc)

    def normaliseUnderMask(self, volume):

        self.fast_sum_stdv *= 0.
        self.fast_sum_mean *= 0.

        self.sumMeanStdv((self.nblocks, 1, 1), (self.num_threads,),
                         (volume, self.mask, self.fast_sum_mean, self.fast_sum_stdv, volume.size),
                         shared_mem=16 * self.num_threads)
        # if self.profile:
        #     t_end = self.stream.record()
        #     t_end.synchronize()
        #
        #     time_took = self.cp.cuda.get_elapsed_time(t_start, t_end)
        #     print(f'nxcc norm s: \t{time_took:.3f}ms')
        #     t_start = self.stream.record()

        # print(meanT[meanT.argmax()-4:meanT.argmax()+2])

        meanT = self.fast_sum_mean.sum() / self.p
        stdT = self.cp.sqrt(self.fast_sum_stdv.sum() / self.p - meanT * meanT)


        return ((volume - meanT)/ stdT) * self.mask


        self.normalize((int(self.cp.ceil(volume.size / self.num_threads)),1,1),(self.num_threads,1,1), (volume, self.mask, meanT, stdT, volume.size))

        return volume

        print(volume.sum())

    def store_particle(self, id, filename, profile=True, binning=1):
        from pytom.tompy.transform import resize, fourier_reduced2full, resizeFourier, fourier_full2reduced
        if self.profile:
            t_start = self.stream.record()

        if self.binning != 1:
            volume = read(filename)
            volumef = fourier_full2reduced(self.fftnP(volume, plan=self.fftplan2))
            volumef = resizeFourier(volumef, 1. / self.binning, (volumef.shape[0] % 2) == 1)
            # volume,volumef = resize(read(filename), 1./self.binning)
            p = self.ifftnP(fourier_reduced2full(volumef, isodd=(volumef.shape[0] % 2) == 1) * self.lpf,
                            plan=self.ifftplan).real.astype(self.xp.float32)
        else:
            p = self.ifftnP(self.fftnP(read(filename), plan=self.fftplan) * self.lpf, plan=self.ifftplan).real.astype(
                self.xp.float32)

        self.formatting[id] = self.staticVolume(p, device=self.device, interpolation=self.interpolation)

        if self.profile:
            t_end = self.stream.record()
            t_end.synchronize()

            time_took = self.cp.cuda.get_elapsed_time(t_start, t_end)
            print(f'read/bin particle: \t{time_took:.3f}ms')


class NonUniformFFTPlan():
    def __init__(self):
        from pytom.gpu.kernels import spmvh_kernel_text, spmv_kernel_text, sum_weighted_norm_complex_array_text, sum_norm_complex_array_text, cTensorCopy_text, cTensorMultiply_text, sum_text

        #kernels
        self.spmvh = xp.RawKernel(spmvh_kernel_text, 'spmvh')
        self.spmv = xp.RawKernel(spmv_kernel_text, 'spmv')
        self.sum_weighted_norm_complex_array = xp.RawKernel(sum_weighted_norm_complex_array_text, 'sum_weighted_norm_complex_array')
        self.sum_norm_complex_array = xp.RawKernel(sum_norm_complex_array_text, 'sum_norm_complex_array')
        self.sum = xp.RawKernel(sum_text, 'sum')
        self.cTensorCopy = xp.RawKernel(cTensorCopy_text, 'cTensorCopy')
        self.cTensorMultiply = xp.RawKernel(cTensorMultiply_text, 'cTensorMultiply')

        self.num_threads = 1024

class NonUniformFFT():
    def __init__(self, API='cuda', device_number=None,
                 verbosity=0):
        """
        Constructor.

        :param API: The API for the heterogeneous system. API='cuda'
                    or API='ocl'
        :param platform_number: The number of the platform found by the API.
        :param device_number: The number of the device found on the platform.
        :param verbosity: Defines the verbosity level, default value is 0
        :type API: string
        :type platform_number: integer
        :type device_number: integer
        :type verbosity: integer
        :returns: 0
        :rtype: int, float

        :Example:

        >>> from pynufft import NUFFT_hsa
        >>> NufftObj = NUFFT_hsa(API='cuda', platform_number=0,
                                         device_number=0, verbosity=0)
        """

        import numpy
        self.dtype = numpy.complex64
        self.verbosity = verbosity

        if self.verbosity > 0:
            print('The choosen API by the user is ', API)

        device_number = 0 if device_number is None else device_number

        xp.cuda.Device(device_number).use()
        self.wavefront = 32
        self.plan = NonUniformFFTPlan()

    def set_wavefront(self, wf):
        self.wavefront = int(wf)  # api.DeviceParameters(device).warp_size
        if self.verbosity > 0:
            print('Wavefront of OpenCL (as wrap of CUDA) = ', self.wavefront)


    def planH(self, om, Nd, Kd, Jd, ft_axes=None, batch=None, radix=None):
        """
        Design the multi-coil or single-coil memory reduced interpolator.


        :param om: The M off-grid locations in the frequency domain.
                   Normalized between [-pi, pi]
        :param Nd: The matrix size of equispaced image.
                    Example: Nd=(256, 256) for a 2D image;
                             Nd = (128, 128, 128) for a 3D image
        :param Kd: The matrix size of the oversampled frequency grid.
                   Example: Kd=(512,512) for 2D image;
                   Kd = (256,256,256) for a 3D image
        :param Jd: The interpolator size.
                   Example: Jd=(6,6) for 2D image;
                   Jd = (6,6,6) for a 3D image
        :param ft_axes: The dimensions to be transformed by FFT.
                   Example: ft_axes = (0, 1) for 2D,
                   ft_axes = (0, 1, 2) for 3D;
                   ft_axes = None for all dimensions.
        :param batch: Batch NUFFT.
                    If provided, the shape is Nd + (batch, ).
                    The last axis is the number of parallel coils.
                    batch = None for single coil.
        :param radix: ????.
                    If provided, the shape is Nd + (batch, ).
                    The last axis is the number of parallel coils.
                    batch = None for single coil.
        :type om: numpy.float array, matrix size = (M, ndims)
        :type Nd: tuple, ndims integer elements.
        :type Kd: tuple, ndims integer elements.
        :type Jd: tuple, ndims integer elements.
        :type ft_axes: tuple, selected axes to be transformed.
        :type batch: int or None
        :returns: 0
        :rtype: int, float
        """
        from pynufft import helper
        import numpy

        self.ndims = len(Nd)  # dimension
        if ft_axes is None:
            ft_axes = range(0, self.ndims)
        self.ft_axes = ft_axes

        self.st = helper.plan(om, Nd, Kd, Jd, ft_axes=ft_axes,
                              format='pELL', radix=radix)
        if batch is None:
            self.parallel_flag = 0
        else:
            self.parallel_flag = 1

        if batch is None:
            self.batch = numpy.uint32(1)
        else:
            self.batch = numpy.uint32(batch)

        self.Nd = self.st['Nd']  # backup
        self.Kd = self.st['Kd']

        if self.batch == 1 and (self.parallel_flag == 0):
            self.multi_Nd = self.Nd
            self.multi_Kd = self.Kd
            self.multi_M = (self.st['M'],)
            # Broadcasting the sense and scaling factor (Roll-off)
            # self.sense2 = self.sense*numpy.reshape(self.sn, self.Nd + (1, ))
        else:  # self.batch is 0:
            self.multi_Nd = self.Nd + (self.batch,)
            self.multi_Kd = self.Kd + (self.batch,)
            self.multi_M = (self.st['M'],) + (self.batch,)

        self.invbatch = 1.0 / self.batch
        self.Kdprod = numpy.uint32(numpy.prod(self.st['Kd']))
        self.Jdprod = numpy.uint32(numpy.prod(self.st['Jd']))
        self.Ndprod = numpy.uint32(numpy.prod(self.st['Nd']))

        self.Nd_elements, self.invNd_elements = helper.strides_divide_itemsize(self.st['Nd'])

        # only return the Kd_elements
        self.Kd_elements = helper.strides_divide_itemsize(self.st['Kd'])[0]
        self.NdCPUorder, self.KdCPUorder, self.nelem = helper.preindex_copy(self.st['Nd'], self.st['Kd'])
        self.offload()

        return 0

    def offload(self):  # API, platform_number=0, device_number=0):
        """
        self.offload():

        Off-load NUFFT to the opencl or cuda device(s)

        :param API: define the device type, which can be 'cuda' or 'ocl'
        :param platform_number: define which platform to be used.
                                The default platform_number = 0.
        :param device_number: define which device to be used.
                            The default device_number = 0.
        :type API: string
        :type platform_number: int
        :type device_number: int
        :return: self: instance
        """
        import numpy
        self.pELL = {}  # dictionary

        self.pELL['nRow'] = numpy.uint32(self.st['pELL'].nRow)
        self.pELL['prodJd'] = numpy.uint32(self.st['pELL'].prodJd)
        self.pELL['sumJd'] = numpy.uint32(self.st['pELL'].sumJd)
        self.pELL['dim'] = numpy.uint32(self.st['pELL'].dim)
        self.pELL['Jd'] = self.to_device(self.st['pELL'].Jd.astype(numpy.uint32))
        self.pELL['meshindex'] = self.to_device(self.st['pELL'].meshindex.astype(numpy.uint32))
        self.pELL['kindx'] = self.to_device(self.st['pELL'].kindx.astype(numpy.uint32))
        self.pELL['udata'] = self.to_device(self.st['pELL'].udata.astype(self.dtype))

        self.volume = {}

        self.volume['Nd_elements'] = self.to_device(numpy.asarray(self.Nd_elements, dtype=numpy.uint32))
        self.volume['Kd_elements'] = self.to_device(numpy.asarray(self.Kd_elements, dtype=numpy.uint32))
        self.volume['invNd_elements'] = self.to_device(self.invNd_elements.astype(numpy.float32))
        self.volume['Nd'] = self.to_device(numpy.asarray(self.st['Nd'], dtype=numpy.uint32))
        self.volume['NdGPUorder'] = self.to_device(self.NdCPUorder)
        self.volume['KdGPUorder'] = self.to_device(self.KdCPUorder)
        self.volume['gpu_coil_profile'] = numpy.ones(self.multi_Nd, dtype=self.dtype)

        Nd = self.st['Nd']
        self.tSN = {}
        self.tSN['Td_elements'] = self.to_device(numpy.asarray(self.st['tSN'].Td_elements, dtype=numpy.uint32))
        self.tSN['invTd_elements'] = self.to_device(self.st['tSN'].invTd_elements.astype(numpy.float32))
        self.tSN['Td'] = self.to_device(numpy.asarray(self.st['tSN'].Td, dtype=numpy.uint32))
        self.tSN['Tdims'] = self.st['tSN'].Tdims
        self.tSN['tensor_sn'] = self.to_device(self.st['tSN'].tensor_sn.astype(numpy.float32))

        self.Ndprod = numpy.int32(numpy.prod(self.st['Nd']))
        self.Kdprod = numpy.int32(numpy.prod(self.st['Kd']))
        self.M = numpy.int32(self.st['M'])

        # if self.batch > 1:  # batch mode
        #     self.fft = reikna.fft.FFT(numpy.empty(self.st['Kd'] + (self.batch,), dtype=self.dtype),self.ft_axes).compile(self.thr, fast_math=False)
        # else:  # elf.Reps == 1 Batch mode is wrong for
        #     self.fft = reikna.fft.FFT(numpy.empty(self.st['Kd'], dtype=self.dtype),self.ft_axes).compile(self.thr, fast_math=False)

        self.zero_scalar = self.dtype(0.0 + 0.0j)
        del self.st['pELL']
        if self.verbosity > 0:
            print('End of offload')

    def reset_sense(self):
        self.volume['gpu_coil_profile'] *= 0.

    def set_sense(self, coil_profile):
        if coil_profile.shape != self.multi_Nd:
            print('The shape of coil_profile is ', coil_profile.shape)
            print('But it should be', self.Nd + (self.batch,))
            raise ValueError
        else:
            self.volume['gpu_coil_profile'] = self.to_device(
                coil_profile.astype(self.dtype))
            if self.verbosity > 0:
                print('Successfully loading coil sensitivities!')

        # if coil_profile.shape == self.Nd + (self.batch, ):

    def to_device(self, image, shape=None):
        import numpy
        print(str(image.dtype))
        ddict = {'complex64' : xp.complex64, 'float32': xp.float32, 'uint32': xp.uint32, 'int32': xp.int32, 'complex128': xp.complex128}
        g_image = xp.array(image, dtype=ddict[str(image.dtype)])
        return g_image

    def x2xx(self, x):
        # xx = self.thr.array(xx.shape, dtype = self.dtype)
        # self.thr.copy_array(z, dest=xx, )
        # size = int(xx.nbytes/xx.dtype.itemsize))
        # Hack of error in cuda backends; 8 is the byte of numpy.complex64
        # size = int(xx.nbytes/8)
        import numpy
        xx = x.copy()
        nblocks = int(xp.ceil(self.Ndprod / self.plan.num_threads))
        nthreads = min(self.plan.num_threads, self.Ndprod)

        self.plan.cTensorMultiply((nblocks, 1, 1), (nthreads,1 ,1),
                                  (numpy.uint32(self.batch), numpy.uint32(self.tSN['Tdims']), self.tSN['Td'], self.tSN['Td_elements'],
                                   self.tSN['invTd_elements'], self.tSN['tensor_sn'], xx, numpy.uint32(0)))
                                 # local_size=None,
                                 # global_size=int(self.batch * self.Ndprod))

        return xx

    def q2y(self, q):
        """
        Private: interpolation by the Sparse Matrix-Vector Multiplication
        """

        nblocks = int(xp.ceil(self.pELL['nRow'] * self.wavefront / self.plan.num_threads))
        nthreads = min(self.plan.num_threads, xp.ceil(self.pELL['nRow'] * self.wavefront))

        y = xp.zeros(self.multi_M, dtype=xp.complex64)
        print(q.max(), q.dtype)
        self.plan.spmv((nblocks, 1, 1),(nthreads,1,1), (self.batch, self.pELL['nRow'],self.pELL['prodJd'],self.pELL['sumJd'], self.pELL['dim'],
                                 self.pELL['Jd'], self.pELL['meshindex'],self.pELL['kindx'], self.pELL['udata'], q, y),shared_mem=8 * nthreads)
                                 # local_size=int(self.wavefront), global_size=int(self.pELL['nRow'] * self.batch * self.wavefront))


        return y

    def y2q(self, y):
        """
        Private: gridding by the Sparse Matrix-Vector Multiplication
        However, serial atomic add is far too slow and inaccurate.
        """
        import matplotlib
        matplotlib.use('Qt5Agg')
        from pylab import imshow, show, plot

        k   = xp.zeros(self.multi_Kd, dtype=xp.complex64)
        res = xp.zeros(self.multi_Kd, dtype=xp.complex64) # array which saves the residue of two sum

        nblocks = int(xp.ceil(self.pELL['nRow'] / self.plan.num_threads))
        nthreads = min(self.pELL['nRow'], self.plan.num_threads)


        # print(self.pELL['udata'].shape)
        # for kk in range(4):
        #     plot(abs(self.pELL['udata'][kk]).get())
        #     show()



        self.plan.spmvh((nblocks,1,1),(nthreads, 1, 1),
                                   ((self.batch), (self.pELL['nRow']), (self.pELL['prodJd']), (self.pELL['sumJd']),
                                    (self.pELL['dim']),self.pELL['Jd'], self.pELL['meshindex'], self.pELL['kindx'],
                                    self.pELL['udata'], k, res, y))

        #     local_size=None,
        #     global_size=int(self.pELL['nRow'] * self.batch)  # *
        #     #                                             int(self.pELL['prodJd']) * int(self.batch))
        # )

        # imshow(abs(y.get())**2,norm=matplotlib.colors.LogNorm())
        # show()
        rr = k + res
        print(rr.imag.max())
        return rr

    def q2xx(self, k):
        """
        Private: the inverse FFT and image cropping (which is the reverse of
        _xx2k() method)
        """
        import numpy
        print('k', k.shape)
        x = xp.fft.fftshift( xp.fft.ifftn(xp.fft.fftshift(k)))
        xx = xp.zeros(self.multi_Nd, dtype=self.dtype)

        nblocks = int(xp.ceil(self.Ndprod / self.plan.num_threads))

        self.plan.cTensorCopy((nblocks, 1, 1),(self.plan.num_threads,1,1),(self.batch, numpy.uint32(self.ndims), self.volume['Nd_elements'],
                                     self.volume['Kd_elements'], self.volume['invNd_elements'], x, xx, numpy.int32(-1)))
            # local_size=None,
            # global_size=int(self.Ndprod))

        return xx


    def release(self):
        del self.volume
        del self.plan
        del self.pELL


    def sum_weighted_norm_complex_array(self, volume, weights):
        nblocks = int(xp.ceil(volume.size / self.plan.num_threads / 2))
        fast_sum_mean = xp.zeros((nblocks), dtype=xp.float32)
        print(fast_sum_mean.sum(), weights.shape, volume.shape)
        self.plan.sum_weighted_norm_complex_array((nblocks, 1,), (self.plan.num_threads, 1, 1),
                                     (volume, weights, fast_sum_mean, volume.size),
                                     shared_mem=8 * self.plan.num_threads)

        if fast_sum_mean.size > self.plan.num_threads:
            nblocks2 = int(xp.ceil(fast_sum_mean.size / self.plan.num_threads / 2))
            fast_sum = xp.zeros((nblocks), dtype=xp.float32)
            self.plan.sum((nblocks2, 1,), (self.plan.num_threads, 1, 1), (fast_sum_mean, fast_sum, fast_sum_mean.size),
                          shared_mem=8*self.plan.num_threads)
            return fast_sum.sum()
        else:
            return fast_sum_mean.sum()

    def sum_norm_complex_array(self, volume):
        nblocks = int(xp.ceil(volume.size / self.plan.num_threads / 2))
        fast_sum_mean = xp.zeros((nblocks), dtype=xp.float32)
        self.plan.sum_norm_complex_array((nblocks, 1,), (self.plan.num_threads, 1, 1),
                                                  (volume, fast_sum_mean, volume.size),
                                                  shared_mem=8 * self.plan.num_threads)

        if fast_sum_mean.size > self.plan.num_threads:
            nblocks2 = int(xp.ceil(fast_sum_mean.size / self.plan.num_threads / 2))
            fast_sum = xp.zeros((nblocks), dtype=xp.float32)
            self.plan.sum((nblocks2, 1,), (self.plan.num_threads, 1, 1), (fast_sum_mean, fast_sum, fast_sum_mean.size),
                          shared_mem=8 * self.plan.num_threads)
            return fast_sum.sum()
        else:
            return fast_sum_mean.sum()

    def solve(self, gy, *args, **kwargs):

        """
        The solver of NUFFT_hsa

        :param gy: data, reikna array, (M,) size
        :param solver: could be 'cg', 'L1TVOLS', 'L1TVLAD'
        :param maxiter: the number of iterations
        :type gy: numpy/cupy array, dtype = numpy.complex64 or cupy.complex64
        :type maxiter: int
        :type damping: numpy.array with size Kd of dtype complex64
        :type weights: numpy.array with size M of dtype complex64
        :return: numpy array with size Nd
        """

        gy = self.to_device(gy)
        q_hat_iter = self.y2q(gy)

        damping = xp.ones_like(q_hat_iter,dtype=xp.float32) if not 'damping' in kwargs.keys() else kwargs['damping']
        weights = xp.ones_like(gy,dtype=xp.float32) if not 'weights' in kwargs.keys() else kwargs['weights']

        epsilon = 1E-6

        r_iter = self.q2y(q_hat_iter)

        r_iter -= gy

        z_hat_iter = self.y2q(weights * r_iter)
        p_hat_iter = z_hat_iter.copy()

        rsold = self.sum_weighted_norm_complex_array(r_iter, weights)
        zsold = self.sum_weighted_norm_complex_array(z_hat_iter, damping)

        for pp in range(0, kwargs['maxiter']):
            q_hat = damping * p_hat_iter
            y_iter = self.q2y(q_hat)
            rsnew = self.sum_weighted_norm_complex_array(y_iter, weights)
            alpha = rsold / rsnew
            rsold = rsnew
            q_hat_iter += alpha * damping * p_hat_iter
            r_iter -= alpha * weights * y_iter
            z_hat_iter = self.y2q(weights * r_iter)
            zsnew = self.sum_weighted_norm_complex_array(z_hat_iter, damping)
            beta = zsnew / zsold
            zsold = zsnew
            p_hat_iter = beta * p_hat_iter + z_hat_iter

            if rsnew < epsilon:
                break

        x2 = xp.fft.fftshift(self.q2xx(q_hat_iter))  # or is it f_hat?
        print(x2.shape, x2.dtype, self.Ndprod)
        try:
            x2 /= self.volume['SnGPUArray']
        except:
            import numpy
            nblocks = int(numpy.ceil(self.Ndprod / self.plan.num_threads))

            self.plan.cTensorMultiply((nblocks,1,1), (self.plan.num_threads,1,1), (numpy.uint32(self.batch),
                                      numpy.uint32(self.tSN['Tdims']),
                                      self.tSN['Td'],
                                      self.tSN['Td_elements'],
                                      self.tSN['invTd_elements'],
                                      self.tSN['tensor_sn'],
                                      x2,
                                      numpy.uint32(1)))

        return x2



