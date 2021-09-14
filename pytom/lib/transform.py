import swig_nufft
import numpy as np

def fourier_interpolate_vol(v, nodes):
    """Given a real 3D volume data in PyTom format, interpolate at specified locations in Fourier space.
    @param v: 3D volume data in PyTom format.
    @dtype v: Pytom_volume.vol
    @param nodes: a list of locations in Fourier space. Each node must be in the range [-0.5, 0.5).
    @dtype nodes: List. [[-0.5, -0.5, 0.3], ...]
    @return: interpolated complex values in a list.
    @rtype: List.
    """
    
    dim_x, dim_y, dim_z = v.shape

    # shift it
    v = np.fft.fftshift(v)

    # linearize it
    v = v.reshape((v.size))

    # linearize it
    nodes = np.array(nodes, dtype="double")
    nodes = nodes.reshape((nodes.size))

    # call the function
    res = nufft.fourier_interpolate_3d(v, dim_x, dim_y, dim_z, nodes, nodes.size/3)

    # convert to complex format
    res = res[::2] + 1j*res[1::2]

    return res

def fourier_rotate_vol(v, angle):
    """
    """
    if v.dtype != np.dtype('float32') or np.isfortran(v) is False:
        v = np.array(v, dtype="float32", order="F")

    # get the rotation matrix
    from pytom.basic.structures import Rotation
    rot = Rotation(angle)
    m = rot.toMatrix()
    m = m.transpose() # get the invert rotation matrix
    m = np.array([m.getRow(0), m.getRow(1), m.getRow(2)], dtype="float32")
    m = m.flatten()

    dim_x, dim_y, dim_z = v.shape

    # # Do the ifftshift first and then fftshift!
    # v = np.fft.fftshift(np.fft.ifftshift(v))

    # linearize it
    v = v.reshape((v.size))

    # allocate the memory for the result
    res = np.zeros(v.size*2, dtype="float32")

    # call the low level c function
    swig_nufft.fourier_rotate_vol(v, dim_x, dim_y, dim_z, m, res)

    res = res[::2] + 1j*res[1::2]
    res = res.reshape((dim_x, dim_y, dim_z), order='F')

    # inverse fft
    ans = np.fft.fftshift(np.real(np.fft.ifftn(np.fft.ifftshift(res))))

    return np.array(ans, dtype="float32", order="F")

def fourier_shift_vol(v, shift):
    pass
