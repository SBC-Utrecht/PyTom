import swig_nufft
import numpy as np

def fourier_interpolate_3d(v, dim_x, dim_y, dim_z, nodes, num_nodes):
	"""Wrapper to the low level c function. Should not be used directly.
	@param v: linearized 3D volume data in numpy format.
	@dtype v: numpy.array
	@param dim_x:
	"""
	if v.__class__ != np.ndarray or nodes.__class__ != np.ndarray:
		raise RuntimeError("Input args must have numpy.ndarray type!")
	
	# allocate the result memory
	res = np.zeros(num_nodes*2, dtype="float32")

	swig_nufft.fourier_interpolate_3d(v, dim_x, dim_y, dim_z, nodes, res)

	return res

def node2index(node, size):
	if node < -0.5 or node >= 0.5:
		raise RuntimeError("Node should be in [-0.5, 0.5) !")

	return int(node*size+size/2)

def index2node(index, size):
	if index < 0 or index > size:
		raise RuntimeError("Index should be in range [0, size] !")

	return float(index)/size - 0.5

def node2index_3d(node, size):
	if size.__class__ == int:
		size = [size, size, size]
	elif size.__class__ == list:
		pass
	else:
		raise RuntimeError("Size must be an int or list!")

	return [node2index(node[0], size[0]), node2index(node[1], size[1]), node2index(node[2], size[2])]

def index2node_3d(index, size):
	if size.__class__ == int:
		size = [size, size, size]
	elif size.__class__ == list:
		pass
	else:
		raise RuntimeError("Size must be an int or list!")

	return [index2node(index[0], size[0]), index2node(index[1], size[1]), index2node(index[2], size[2])]

def get_spherical_nodes(r, size, b):
    """Get the nodes in NFFT format
    """
    from math import pi, cos, sin

    res = []
    for j in xrange(2*b):
        for k in xrange(2*b):
            the = pi*(2*j+1)/(4*b) # (0,pi)
            phi = pi*k/b # [0,2*pi)
            res.append([r*cos(phi)*sin(the)/size, r*sin(phi)*sin(the)/size, r*cos(the)/size])

    return res