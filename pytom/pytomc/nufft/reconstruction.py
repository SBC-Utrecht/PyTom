import swig_nufft
import numpy as np

def fourier_3d_iter_reconstruct(projections, tilt_angles, iteration=10):
    """Do the reconstruction on a projection list using direct Fourier inversion (3D iterative).
    @param projections: projections
    @param tilt_angles: tilt angles of projections
    """
    from pytom.basic.structures import Rotation

    if len(projections) != len(tilt_angles):
        raise Exception("Length of projections and tilt angles not consistent!")

    freal = None # Fourier coefficients
    fimag = None
    kx = None # nodes positions
    ky = None
    kz = None
    weights = None # weights

    for i in range(len(tilt_angles)):
        proj = projections[i]
        NX, NY, NZ = proj.shape
        assert NZ == 1 and NX == NY # for now only supports the square
        proj = proj.reshape((NX, NY)) # make it 2D

        if i == 0:
            # nodes in 0 degree
            x, y = np.meshgrid(np.array(np.linspace(-0.5, 0.5, NX, endpoint=False), dtype="float32"), np.array(np.linspace(-0.5, 0.5, NY, endpoint=False), dtype="float32"))
            x = x.reshape((1, x.size))
            y = y.reshape((1, y.size))
            z = np.zeros(x.size, dtype="float32")
            xyz = np.vstack((x,y,z))

        # do the 2D Fourier transform on each projection
        fproj = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(proj)))
        fproj_real = np.array(np.real(fproj), dtype="float32")
        fproj_imag = np.array(np.imag(fproj), dtype="float32")

        if freal is None:
            freal = fproj_real.reshape((fproj_real.size), order='F') # to be consistent with the order of kx,ky,kz
            fimag = fproj_imag.reshape((fproj_imag.size), order='F')
        else:
            freal = np.hstack((freal, fproj_real.reshape((fproj_real.size), order='F') ))
            fimag = np.hstack((fimag, fproj_imag.reshape((fproj_imag.size), order='F') ))

        # calculate the nodes positions
        ang = tilt_angles[i] # get the tilt angle of this project
        rot = Rotation([270, 90, ang]) # rotate around y axis
        m = rot.toMatrix()
        m = np.array([m.getRow(0), m.getRow(1), m.getRow(2)], dtype="float32")
        kk = np.dot(m, xyz)
        if kx is None:
            kx = kk[0]
            ky = kk[1]
            kz = kk[2]
        else:
            kx = np.hstack((kx, kk[0]))
            ky = np.hstack((ky, kk[1]))
            kz = np.hstack((kz, kk[2]))

        # compute the weights
        if i == 0:
            ang_interval = tilt_angles[1]-tilt_angles[0]
        elif i == len(tilt_angles)-1:
            ang_interval = tilt_angles[-1]-tilt_angles[-2]
        else:
            ang_interval = (tilt_angles[i+1] - tilt_angles[i-1])/2
        w = np.abs(ang_interval/360.*(np.array(list(range(NX//2, -NX//2, -1)), dtype="float32")*2))
        wei = np.tile(w, (NY, 1))
        if weights is None:
            weights = wei.reshape((wei.size))
        else:
            weights = np.hstack((weights, wei.reshape((wei.size))))

    # reconstruct
    M = freal.size
    Z = NX
    res_real = np.zeros(NX*NY*Z, dtype="float32")
    res_imag = np.zeros(NX*NY*Z, dtype="float32")
    weights = np.float32(weights)

    swig_nufft.fourier_3d_iter_reconstruct(freal, fimag, NX, Z, M, weights, kx, ky, kz, res_real, res_imag, iteration)

    # resize and return the real part
    res = res_real.reshape((NX,NY,Z))

    return res


def fourier_2d1d_iter_reconstruct(projections, tilt_angles, iteration=10, err=None):
    """Do the reconstruction on a projection list using direct Fourier inversion (2D+1D iterative).
    For noisy situation, the #iter should be small
    @param projections: projections
    @param tilt_angles: tilt angles of projections
    @param err: the percentage of allowed error
    """
    from pytom.basic.structures import Rotation

    if len(projections) != len(tilt_angles):
        raise Exception("Length of projections and tilt angles not consistent!")

    freal = None # Fourier coefficients
    fimag = None
    kx = None # nodes positions
    kz = None
    weights = None # weights

    for i in range(len(tilt_angles)):
        proj = projections[i]
        if len(proj.shape) == 3:
            NX, NY, NZ = proj.shape
            assert NZ == 1 and NX == NY # for now only supports the square
            proj = proj.reshape((NX, NY)) # make it 2D
        elif len(proj.shape) == 2:
            NX, NY = proj.shape
            NZ = 1
        else:
            raise Exception("Projection shape incorrect!")

        if i == 0:
            # nodes in 0 degree
            x, y = np.meshgrid(np.array(np.linspace(-0.5, 0.5, NX, endpoint=False), dtype="float32"), np.array(np.linspace(-0.5, 0.5, NY, endpoint=False), dtype="float32"))
            x = x.reshape((1, x.size))
            y = y.reshape((1, y.size))
            z = np.zeros(x.size, dtype="float32")
            xyz = np.vstack((x,y,z))

            freal = np.zeros((NX, NY, len(tilt_angles)), dtype="float32")
            fimag = np.zeros((NX, NY, len(tilt_angles)), dtype="float32")

        # do the 2D Fourier transform on each projection
        fproj = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(proj)))
        freal[:,:,i] = np.array(np.real(fproj), dtype="float32")
        fimag[:,:,i] = np.array(np.imag(fproj), dtype="float32")

        # calculate the nodes positions
        ang = tilt_angles[i] # get the tilt angle of this projection
        rot = Rotation([270, 90, ang]) # rotate around y axis
        m = rot.toMatrix()
        m = np.array([m.getRow(0), m.getRow(1), m.getRow(2)], dtype="float32")
        kk = np.dot(m, xyz)
        kk = np.array(kk, dtype="float32")
        if kx is None:
            kx = kk[0][:NX]
            kz = kk[2][:NX]
        else:
            kx = np.hstack((kx, kk[0][:NX]))
            kz = np.hstack((kz, kk[2][:NX]))

        # compute the weights
        if i == 0:
            ang_interval = tilt_angles[1]-tilt_angles[0]
        elif i == len(tilt_angles)-1:
            ang_interval = tilt_angles[-1]-tilt_angles[-2]
        else:
            ang_interval = (tilt_angles[i+1] - tilt_angles[i-1])/2
        w = np.abs(ang_interval/180.*np.pi*(np.array(list(range(NX//2, -NX//2, -1)), dtype="float32")*2)/(NX//2)**2)
        w[NX//2] = np.pi/(4*len(tilt_angles)*(NX//2)**2) # take care of the duplicates
        if weights is None:
            weights = w.reshape((w.size))
        else:
            weights = np.hstack((weights, w.reshape((w.size))))

    # change the Y and Z axis
    freal = np.swapaxes(freal, 1,2)
    freal = freal.reshape((freal.size), order="F")
    fimag = np.swapaxes(fimag, 1,2)
    fimag = fimag.reshape((fimag.size), order="F")

    # construct the damping factors
    damping = np.ones((NX, NX), dtype="float32")
    # radius = 25
    # sigma = 5
    # [xx, yy] = np.mgrid[0:NX, 0:NX]
    # rr = np.sqrt((xx-NX/2)**2+(yy-NX/2)**2)
    # ind = rr>radius
    # damping[ind] = np.exp(-((rr[ind] - radius)/sigma)**2/2)
    damping = damping.reshape((damping.size))

    # estimate the variance of the noise
    if err and err > 0:
        data = np.array(projections)
        var_all = np.var(data)
        var_err = err*var_all
        threshold = var_err
    else:
        threshold = -1

    # reconstruct
    M = freal.size
    Z = NX
    res_real = np.zeros(NX*NY*Z, dtype="float32")
    res_imag = np.zeros(NX*NY*Z, dtype="float32")
    weights = np.float32(weights)

    # weird!
    print(weights.shape, weights.dtype)
    swig_nufft.fourier_2d1d_iter_reconstruct(freal, fimag, NX, Z, M, weights, kz, kx, res_real, res_imag, iteration, threshold, damping)

    # resize and return the real part
    res = res_real.reshape((NX,NY,Z), order="F") # weird!
    res = np.swapaxes(res, 1,2) # weird!

    return res


def fourier_2d1d_gridding_reconstruct(projections, tilt_angles):
    """Do the reconstruction on a projection list using direct Fourier inversion (2D+1D gridding).
    @param projections: projections
    @param tilt_angles: tilt angles of projections
    """
    from pytom.basic.structures import Rotation

    if len(projections) != len(tilt_angles):
        raise Exception("Length of projections and tilt angles not consistent!")

    freal = None # Fourier coefficients
    fimag = None
    kx = None # nodes positions
    kz = None
    weights = None # weights

    for i in range(len(tilt_angles)):
        proj = projections[i]
        NX, NY, NZ = proj.shape
        assert NZ == 1 and NX == NY # for now only supports the square
        proj = proj.reshape((NX, NY)) # make it 2D

        if i == 0:
            # nodes in 0 degree
            x, y = np.meshgrid(np.array(np.linspace(-0.5, 0.5, NX, endpoint=False), dtype="float32"), np.array(np.linspace(-0.5, 0.5, NY, endpoint=False), dtype="float32"))
            x = x.reshape((1, x.size))
            y = y.reshape((1, y.size))
            z = np.zeros(x.size, dtype="float32")
            xyz = np.vstack((x,y,z))

            freal = np.zeros((NX, NY, len(tilt_angles)), dtype="float32")
            fimag = np.zeros((NX, NY, len(tilt_angles)), dtype="float32")

        # do the 2D Fourier transform on each projection
        fproj = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(proj)))
        freal[:,:,i] = np.array(np.real(fproj), dtype="float32")
        fimag[:,:,i] = np.array(np.imag(fproj), dtype="float32")

        # calculate the nodes positions
        ang = tilt_angles[i] # get the tilt angle of this project
        rot = Rotation([270, 90, ang]) # rotate around y axis
        m = rot.toMatrix()
        m = np.array([m.getRow(0), m.getRow(1), m.getRow(2)], dtype="float32")
        kk = np.dot(m, xyz)
        if kx is None:
            kx = kk[0][:NX]
            kz = kk[2][:NX]
        else:
            kx = np.hstack((kx, kk[0][:NX]))
            kz = np.hstack((kz, kk[2][:NX]))

        # compute the weights
        if i == 0:
            ang_interval = tilt_angles[1]-tilt_angles[0]
        elif i == len(tilt_angles)-1:
            ang_interval = tilt_angles[-1]-tilt_angles[-2]
        else:
            ang_interval = (tilt_angles[i+1] - tilt_angles[i-1])/2
        w = np.abs(ang_interval/180.*np.pi*(np.array(list(range(NX//2, -NX//2, -1)), dtype="float32")*2)/(NX//2)**2)
        w[NX//2] = np.pi/(4*len(tilt_angles)*(NX//2)**2) # take care of the duplicates
        if weights is None:
            weights = w.reshape((w.size))
        else:
            weights = np.hstack((weights, w.reshape((w.size))))

    # change the Y and Z axis
    freal = np.swapaxes(freal, 1,2)
    freal = freal.reshape((freal.size), order="F")
    fimag = np.swapaxes(fimag, 1,2)
    fimag = fimag.reshape((fimag.size), order="F")

    # reconstruct
    M = freal.size
    Z = NX
    res_real = np.zeros(NX*NY*Z, dtype="float32")
    res_imag = np.zeros(NX*NY*Z, dtype="float32")
    weights = np.float32(weights)

    # weird!
    swig_nufft.fourier_2d1d_gridding_reconstruct(freal, fimag, NX, Z, M, weights, kz, kx, res_real, res_imag)

    # resize and return the real part
    res = res_real.reshape((NX,NY,Z), order="F") # weird!
    res = np.swapaxes(res, 1,2) # weird!

    return res

