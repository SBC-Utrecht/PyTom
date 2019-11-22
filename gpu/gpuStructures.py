import threading
import importlib
from pytom.tompy.io import read, write

class TemplateMatchingPlan():
    def __init__(self, volume, template, mask, wedge, cp, vt, calc_stdV, pad, get_fft_plan, deviceid):
        print(deviceid)
        cp.cuda.Device(deviceid).use()
        import voltools
        importlib.reload(voltools)

        self.volume = cp.asarray(volume,dtype=cp.float32,order='C')
        self.volume_fft = cp.fft.rfftn(self.volume)
        self.sPad = volume.shape

        self.mask = cp.asarray(mask,dtype=cp.float32, order='C')
        self.maskPadded = cp.zeros_like(self.volume).astype(cp.float32)
        self.sOrg = mask.shape
        pad(self.mask, self.maskPadded, self.sPad, self.sOrg)

        self.templateOrig = cp.asarray(template, dtype=cp.float32,order='C')
        self.template = cp.asarray(template, dtype=cp.float32,order='C')
        self.texture = voltools.StaticVolume(self.template, interpolation=vt.Interpolations.FILT_BSPLINE, deviceid=deviceid)
        self.templatePadded = cp.zeros_like(self.volume)

        self.wedge = cp.asarray(wedge,order='C')

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

        sx, sy, sz = self.plan.template.shape
        SX, SY, SZ = self.plan.templatePadded.shape
        CX, CY, CZ = SX // 2, SY // 2, SZ // 2
        cx, cy, cz = sx // 2, sy // 2, sz // 2
        mx, my, mz = sx %  2, sy %  2, sz %  2

        for angleId, angles in enumerate(angle_list):

            # Rotate
            self.plan.template = self.rotate3d(self.plan.templateOrig, phi=phi,the=the,psi=psi)
            #self.plan.texture.transform(rotation=(angles[0], angles[2], angles[1]), rotation_order='rzxz', output=self.plan.template)

            # Add wedge
            self.plan.template = self.irfftn(self.rfftn(self.plan.template) * self.plan.wedge)

            # Normalize template
            meanT = self.meanUnderMask(self.plan.template, self.plan.mask, p=self.plan.p)
            stdT = self.stdUnderMask(self.plan.template, self.plan.mask, meanT, p=self.plan.p)
            self.plan.template = ((self.plan.template - meanT) / stdT) * self.plan.mask

            # Paste in center
            self.plan.templatePadded[CX-cx:CX+cx+mx, CY-cy:CY+cy+my,CZ-cz:CZ+cz+mz] = self.plan.template

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
        return self.fftshift(self.ifftnP(volume_fft * self.fftnP(template, plan=plan).conj(), plan=plan).real) / (p * norm)
        #return self.fftshift(self.irfftn(volume_fft * self.rfftn(template).conj() )) / (p * norm)

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
    def __init__(self, particle, reference, mask, wedge, maskIsSphere=True):
        import pycuda.gpuarray as gu
        from voltools.volume import Volume
        self.particle = gu.to_gpu(particle)
        self.template = Volume(reference)
        self.wedge = Volume(wedge)
        self.mask = Volume(mask)
        self.mask.d_data = gu.to_gpu(mask)
        self.fwd_plan = Plan(particle.shape, volume.dtype, np.complex64)
        self.inv_plan = Plan(particle.shape, np.complex64, volume.dtype)
        self.volume_fft = gu.zeros_like(self.particle, dtype=np.complex64)
        self.template_fft = gu.zeros_like(self.reference.d_data, dtype=np.complex64)
        self.ccc_map = gu.zeros_like(self.volume, dtype=np.float32)


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
    print(data.shape, grid.shape,dataout.shape)
    dataout = map_coordinates(data, grid.reshape(len(grid), -1), order=order).reshape(grid.shape[1:])

    return dataout
