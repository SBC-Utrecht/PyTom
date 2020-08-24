from pytom.gpu.initialize import xp, device
from pytom.voltools import transform
import numpy

def backProjectGPU(projections, vol_bp, vol_phi, proj_angles, recPosVol=[0,0,0], vol_offsetProjections=[0,0,0], interpolation='filt_bspline'):

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


def backProjectGPU2(vol_img, vol_bp, vol_phi, vol_the, recPosVol=[0,0,0], vol_offsetProjections=[0,0,0]):
    from pytom_numpy import vol2npy
    import cupy as cp
    from cupy import cos, sin
    import voltools

    projections = vol2npy(vol_img)
    proj_angles = vol2npy(vol_the)[0,0,:]

    dims = (vol_bp.sizeX(), vol_bp.sizeY(), vol_bp.sizeZ())
    # TODO add offsets

    # preparation
    reconstruction = cp.zeros(dims, dtype=cp.float32)
    theta_angles = cp.deg2rad(cp.array(proj_angles))

    # theta_angles_rad = (theta_angles)
    assert len(theta_angles) == projections.shape[2] #'Number of angles and projections should match'
    # create coordinates and stack them
    x, y, z = cp.meshgrid(cp.arange(0, dims[0]), cp.arange(0, dims[1]), cp.arange(0, dims[2]), indexing='ij')
    # backproject each projection into reconstruction

    for n in range(projections.shape[2]):
        # get projection
        src = cp.array(projections[:,:,n])

        # previous
        Z1 = Z2 = 0.0
        Y = theta_angles[n]
        x_offset = y_offset = z_offset = 0
        x_proj_offset = y_proj_offset = 0
        tr11 = cos(Y)*cos(Z1)*cos(Z2)-sin(Z1)*sin(Z2)
        # tr21 = cos(Y)*sin(Z1)*cos(Z2)+cos(Z1)*sin(Z2)
        tr31 = -sin(Y)*cos(Z2)
        # tr12 = -cos(Y)*cos(Z1)*sin(Z2)-sin(Z1)*cos(Z2)
        tr22 = -cos(Y)*sin(Z1)*sin(Z2)+cos(Z1)*cos(Z2)
        # tr32 = sin(Y)*sin(Z2)
        reconstruction_center_x = dims[0] // 2 - x_offset
        reconstruction_center_y = dims[1] // 2 - y_offset
        reconstruction_center_z = dims[2] // 2 - z_offset
        projection_center_x = dims[0] // 2 + x_proj_offset
        projection_center_y = dims[1] // 2 + y_proj_offset
        # interpolation_position_x = tr11 * (x - reconstruction_center_x) + tr21 * (y - reconstruction_center_y) + tr31 * (z - reconstruction_center_z) + projection_center_x
        # interpolation_position_y = tr12 * (x - reconstruction_center_x) + tr22 * (y - reconstruction_center_y) + tr32 * (z - reconstruction_center_z) + projection_center_y
        interpolation_position_x = tr11 * (x - reconstruction_center_x) + tr31 * (z - reconstruction_center_z) + projection_center_x
        interpolation_position_y = tr22 * (y - reconstruction_center_y) + projection_center_y
        proj_position_x = interpolation_position_x.astype(int)
        proj_position_y = interpolation_position_y.astype(int)
        interpolation_offset_x = (interpolation_position_x - proj_position_x)
        interpolation_offset_y = (interpolation_position_y - proj_position_y)
        m0, m1, m2 = cp.where((proj_position_x >= 1) & (proj_position_x < dims[0]) & (proj_position_y >= 1) & (proj_position_y < dims[1]))
        value_x1 = src[proj_position_x[m0, m1, m2]-1, proj_position_y[m0, m1, m2]-1] + \
                   (src[proj_position_x[m0, m1, m2], proj_position_y[m0, m1, m2]-1] - src[proj_position_x[m0, m1, m2]-1, proj_position_y[m0, m1, m2]-1]) * interpolation_offset_x[m0, m1, m2]
        value_x2 = src[proj_position_x[m0, m1, m2]-1, proj_position_y[m0, m1, m2]] + \
                   (src[proj_position_x[m0, m1, m2], proj_position_y[m0, m1, m2]] - src[proj_position_x[m0, m1, m2]-1, proj_position_y[m0, m1, m2]]) * interpolation_offset_x[m0, m1, m2]
        value_y = value_x1 + (value_x2 - value_x1) * interpolation_offset_y[m0, m1, m2]
        reconstruction[m0, m1, m2] += value_y
        del value_y
        del value_x1
        del value_x2
        del m0, m1, m2
        # import napari
        # with napari.gui_qt():
        #     viewer = napari.Viewer()
        #     viewer.add_image(reconstruction.get(), name='recon', interpolation='bilinear')
        #     viewer.add_image(np.log10(np.abs(np.fft.fftshift(np.fft.fftn(reconstruction.get())))), name='fft', interpolation='bilinear')
    rec = reconstruction.get()
    return rec

    for i in range(dims[0]):
        for j in range(dims[1]):
            for k in range(dims[2]):
                vol_bp.setV(float(rec[i][j][k]), int(i), int(j), int(k))


    return vol_bp

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
    import voltools
    from cupy.fft import fftshift, ifftn, fftn
    import numpy
    from cupyx.scipy.ndimage import rotate
    interpolation = voltools.Interpolations.FILT_BSPLINE
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
            rot.real =voltools.transform(org.real, rotation=(0, proj_angles[n],0), rotation_order='sxyz',interpolation=interpolation)
            rot.imag =voltools.transform(org.imag, rotation=(0, proj_angles[n],0), rotation_order='sxyz',interpolation=interpolation)
        tot += rot
        del org, rot
    
    if dims[0] % 2:
        t = cp.zeros_like(tot)
        t[1:, 1:, 1:] = tot[:-1, :-1, :-1]
    else:
        t = tot

    return fftshift((ifftn(fftshift(t)).real)).get()
