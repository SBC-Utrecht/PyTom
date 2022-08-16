from pytom.gpu.initialize import xp, device
from pytom.voltools import transform
import numpy

def backProjectGPU2(projections, vol_bp, vol_phi, proj_angles, recPosVol=[0,0,0], vol_offsetProjections=[0,0,0], interpolation='filt_bspline'):

    assert len(proj_angles) == projections.shape[2] #'Number of angles and projections should match'

    for n in range(projections.shape[2]):
        vol_bp += back_projection(projections[:,:,n], proj_angles[n], interpolation)

    return vol_bp

def back_projection(image, angle, interpolation=None):
    '''
    back_projection: Returns 3D volume with values projected under an angle
    @param image: 2d projection image
    @type image: 2d numpy/cupy array of floats
    @param angle: Tilt angle of projection
    @type angle: float32
    @param interpolation: interpolation type using in rotation. filtered bspline is recommended.
    @return:
    @author Gijs van der Schot
    '''

    dimx, dimy = image.shape
    dims = [dimx, dimy, dimx]

    ndim = (dimx*3)//2

    out = xp.dstack([image] * (ndim))

    bp = xp.zeros_like(out)
    transform(out, output=bp, rotation=(0, angle, 0), rotation_order='sxyz', interpolation=interpolation,
              center=numpy.array([dims[0] // 2, dims[1] // 2, ndim//2]), device=device)

    return bp[:,:,ndim//2-dimx//2:ndim//2+dimx//2].copy()

def backProjectGPU(projections, reconstruction, vol_phi, proj_angles, recPosVol=None, vol_offsetProjections=None,
                   interpolation=''):
    from pytom.gpu.initialize import xp, device
    from pytom.gpu.kernels import reconstruction_wbp_text

    cos = xp.cos
    sin = xp.sin

    # preparation
    theta_angles = xp.deg2rad(xp.array(proj_angles))

    dims = xp.array(reconstruction.shape, dtype=xp.int32)

    assert len(theta_angles) == projections.shape[2]  #'Number of angles and projections should match'

    center_recon = xp.zeros((3),dtype=xp.int32)
    center_recon[0] = (dims[0] // 2 + 1) - recPosVol[0,0]
    center_recon[1] = (dims[1] // 2 + 1) - recPosVol[0,1]
    center_recon[2] = (dims[2] // 2 + 1) - recPosVol[0,2]

    dims_proj = xp.array(projections.shape, dtype=xp.int32)

    nthreads = 1024
    nblocks = int(xp.ceil(reconstruction.size / nthreads).get())

    reconstruction_wbp = xp.RawKernel(reconstruction_wbp_text, 'reconstruction_wbp')

    center_proj = xp.zeros_like(vol_offsetProjections, dtype=xp.int32)
    center_proj[:, 0] = (dims_proj[0] // 2 + 1) + vol_offsetProjections[:, 0]
    center_proj[:, 1] = (dims_proj[1] // 2 + 1) + vol_offsetProjections[:, 1]

    for n in range(projections.shape[2]):
        # get projection
        src = xp.array(projections[:, :, n], dtype=xp.float32)

        Z1 = Z2 = 0.0
        Y = theta_angles[n]

        tr11 = cos(Y)*cos(Z1)*cos(Z2)-sin(Z1)*sin(Z2)
        tr21 = cos(Y)*sin(Z1)*cos(Z2)+cos(Z1)*sin(Z2)
        tr31 = -sin(Y)*cos(Z2)
        tr12 = -cos(Y)*cos(Z1)*sin(Z2)-sin(Z1)*cos(Z2)
        tr22 = -cos(Y)*sin(Z1)*sin(Z2)+cos(Z1)*cos(Z2)
        tr32 = sin(Y)*sin(Z2)
        tr = xp.array([tr11, tr21, tr31, tr12, tr22, tr32],dtype=xp.float32)

        reconstruction_wbp((nblocks,1,1,), (nthreads,1,1), (src, center_proj[n,:], dims_proj, reconstruction, center_recon, dims, tr, reconstruction.size))


def exactFilter(tilt_angles, tiltAngle, sX, sY, sliceWidth, arr=[]):
    """
    exactFilter: Generates the exact weighting function required for weighted backprojection - y-axis is tilt axis
    Reference : Optik, Exact filters for general geometry three dimensional reconstuction, vol.73,146,1986.
    @param tilt_angles: list of all the tilt angles in one tilt series
    @param titlAngle: tilt angle for which the exact weighting function is calculated
    @param sizeX: size of weighted image in X
    @param sizeY: size of weighted image in Y

    @return: filter volume

    """

    from cupy import array, matrix, sin, pi, arange, float32, column_stack, argmin, clip, ones, ceil

    # Using Friedel Symmetry in Fourier space.
    # sY = sY // 2 + 1

    # Calculate the relative angles in radians.
    diffAngles = (array(tilt_angles) - tiltAngle) * pi / 180.

    # Closest angle to tiltAngle (but not tiltAngle) sets the maximal frequency of overlap (Crowther's frequency).
    # Weights only need to be calculated up to this frequency.
    sampling = min(abs(diffAngles)[abs(diffAngles) > 0.001])
    crowtherFreq = min(sX // 2, int(ceil(1 / sin(sampling))))
    arrCrowther = matrix(abs(arange(-crowtherFreq, min(sX // 2, crowtherFreq + 1))))

    # Calculate weights
    wfuncCrowther = 1. / (clip(1 - array(matrix(abs(sin(diffAngles))).T * arrCrowther) ** 2, 0, 2)).sum(axis=0)

    # Create full with weightFunc
    wfunc = ones((sX, sY, 1), dtype=float32)
    wfunc[sX // 2 - crowtherFreq:sX // 2 + min(sX // 2, crowtherFreq + 1), :, 0] = column_stack(
        ([(wfuncCrowther), ] * (sY))).astype(float32)
    return wfunc

def fourierReconstructGPU(vol_img, vol_bp, vol_phi, vol_the, recPosVol=[0,0,0], vol_offsetProjections=[0,0,0]):
    from pytom_numpy import vol2npy
    import cupy as cp
    from pytom.voltools import transform
    from cupy.fft import fftshift, ifftn, fftn
    import numpy
    from cupyx.scipy.ndimage import rotate
    interpolation = 'filt_bspline'
    projections = vol2npy(vol_img)
    projections = cp.array(projections)
    proj_angles = vol2npy(vol_the)[0, 0, :]
    dims = (projections.shape[0], projections.shape[0], projections.shape[0])
    tot = cp.zeros(dims, cp.complex64)

    
    for n in range(projections.shape[2]):
        org = cp.zeros(dims, dtype=cp.complex64)
        rot = cp.zeros(dims, dtype=cp.complex64)
        ff = fftshift(fftn(fftshift(projections[:,:,n])))

        org[:, :, dims[0]//2] = ff
        #print(n,ff.shape, org.shape)
        if 0:
            rot.real = rotate(org.real, proj_angles[n], axes=(2, 0), reshape=False)
            rot.imag = rotate(org.imag, proj_angles[n], axes=(2, 0), reshape=False)
        else:
            rot.real = transform(org.real, rotation=(0, proj_angles[n],0), rotation_order='sxyz',interpolation=interpolation)
            rot.imag = transform(org.imag, rotation=(0, proj_angles[n],0), rotation_order='sxyz',interpolation=interpolation)
        tot += rot
        del org, rot
    
    if dims[0] % 2:
        t = cp.zeros_like(tot)
        t[1:, 1:, 1:] = tot[:-1, :-1, :-1]
    else:
        t = tot

    return fftshift((ifftn(fftshift(t)).real)).get()
