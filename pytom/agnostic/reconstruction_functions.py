from pytom.gpu.initialize import xp, device
from pytom.voltools import transform
import numpy


def backProjectGPU(projections, reconstruction, vol_phi, proj_angles, recPosVol=None, vol_offsetProjections=None,
                   interpolation=''):
    """
    Docstring generated with chatGPT!

    ===> This is the stand-up comedian version:

    "Alright folks, are you ready to backproject like a pro? I present to you: the backProjectGPU function. It's like
    taking a trip back in time, but instead of dinosaurs, you got 3D volumes.

    Here's how it works. You give it some projections, a reconstruction and some angles, and it'll give you back a 3D
    volume that's been backprojected using the GPU. It's like magic, but with less rabbits and top hats.

    Now, I know what you're thinking: "But wait, what's a vol_phi?" Well, it's just a fancy way of saying the angles
    of the projections in radians. But don't worry, you don't need to know what radians are, just trust me on this one.

    And if you're feeling adventurous, you can even give it a recPosVol and vol_offsetProjections. But just like with
    radians, you don't really need to know what they are, just trust the process.

    So, give it a try, backproject like a pro, and let me know how it goes! And remember, always backproject
    responsibly."

    ===> This is the proper version:

    backProjectGPU: Returns 3D volume by backprojecting the projections using the GPU
    @param projections: 3D array with projections
    @type projections: cupy array of floats
    @param reconstruction: 3D array for the reconstruction volume
    @type reconstruction: cupy array of floats
    @param vol_phi: 1D array of tilt angles (y-axis rotation) of projections in degrees
    @type vol_phi: cupy array of floats
    @param proj_angles: 1D array of skew angle (x-axis rotation) of projections in degrees
    @type proj_angles: cupy array of floats
    @param recPosVol: position of the volume in the reconstruction space
    @type recPosVol: list of integers
    @param vol_offsetProjections: projection offset in the projection space
    @type vol_offsetProjections: list of integers
    @param interpolation: interpolation type using in rotation. filtered bspline is recommended.
    @type interpolation: str
    @return: reconstruction
    @author Gijs van der Schot, Marten Chaillet
    """
    from pytom.gpu.initialize import xp, device
    from pytom.gpu.kernels import reconstruction_wbp_text

    cos = xp.cos
    sin = xp.sin

    # preparation
    theta_angles = xp.deg2rad(xp.array(proj_angles))

    dims = xp.array(reconstruction.shape, dtype=xp.int32)

    assert len(theta_angles) == projections.shape[2]  #'Number of angles and projections should match'

    center_recon = xp.zeros((3), dtype=xp.int32)
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

        reconstruction_wbp((nblocks,1,1,), (nthreads,1,1), (src, center_proj[n,:], dims_proj, reconstruction,
                                                            center_recon, dims, tr, reconstruction.size))


def exactFilter(tilt_angles, tiltAngle, sX, sY, sliceWidth, arr=[]):
    """
    exactFilter: Generates the exact weighting function required for weighted backprojection - y-axis is tilt axis
    Reference : Optik, Exact filters for general geometry three dimensional reconstuction, vol.73,146,1986.
    @param tilt_angles: List of tilt angles
    @type tilt_angles: list  of floats
    @param tiltAngle: Tilt angle of projection
    @type tiltAngle: float
    @param sX: size of the x-axis
    @type sX: int
    @param sY: size of the y-axis
    @type sY: int
    @param sliceWidth: width slice through the object
    @type sliceWidth: int
    @param arr: array
    @type arr: list
    @return: 2D weight array of floats
    @author Gijs van der Schot
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
    """
    fourierReconstructGPU: Reconstructs a 3D volume from a set of 2D projections, using Fourier reconstruction method.
    @param vol_img: 3D volume of projections
    @type vol_img: cupy array of floats
    @param vol_bp: 3D volume of back-projections
    @type vol_bp: cupy array of floats
    @param vol_phi: 1D array of angles in degrees
    @type vol_phi: cupy array of floats
    @param vol_the: 1D array of theta angles in degrees
    @type vol_the: cupy array of floats
    @param recPosVol: Position of the reconstruction in the 3D volume, default is [0,0,0]
    @type recPosVol: List of integers
    @param vol_offsetProjections: Offset of the projections in the 3D volume, default is [0,0,0]
    @type vol_offsetProjections: List of integers
    @return: 3D volume of reconstructed values
    """
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
        ff = fftshift(fftn(fftshift(projections[:, :, n])))
        org[:, :, dims[0]//2] = ff
        rot.real = transform(org.real, rotation=(0, proj_angles[n],0), rotation_order='sxyz',
                             interpolation=interpolation)
        rot.imag = transform(org.imag, rotation=(0, proj_angles[n],0), rotation_order='sxyz',
                             interpolation=interpolation)
        tot += rot
        del org, rot
    
    if dims[0] % 2:
        t = cp.zeros_like(tot)
        t[1:, 1:, 1:] = tot[:-1, :-1, :-1]
    else:
        t = tot

    return fftshift((ifftn(fftshift(t)).real)).get()
