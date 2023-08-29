from pytom.gpu.initialize import xp
from pytom.gpu.kernels import reconstruction_wbp_text


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
