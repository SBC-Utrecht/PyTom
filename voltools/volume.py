import numpy as np
import numpy as cp
from typing import Tuple, Union
from .transforms import _get_transform_kernel, _bspline_prefilter, affine
from .utils import compute_elementwise_launch_dims, switch_to_device,\
    scale_matrix, shear_matrix, rotation_matrix, translation_matrix, transform_matrix, get_available_devices

try:
    import cupy as cp
except ImportError:
    pass

class StaticVolume:

    def __init__(self, data: np.ndarray, interpolation: str = 'linear', device: str = 'gpu'):

        if data.ndim != 3:
            raise ValueError('Expected a 3D array')
        if device not in get_available_devices():
            raise ValueError(f'Unknown device ({device}), must be one of {get_available_devices()}')
        switch_to_device(device)

        self.device = device
        self.interpolation = interpolation

        if device.startswith('gpu'):
            data = cp.array(data)
            self.shape = data.shape
            self.d_shape = cp.asarray(data.shape, dtype=cp.uint32)
            self.d_type = data.dtype
            self.affine_kernel = _get_transform_kernel(interpolation)

            # init texture
            ch = cp.cuda.texture.ChannelFormatDescriptor(32, 0, 0, 0, cp.cuda.runtime.cudaChannelFormatKindFloat)
            arr = cp.cuda.texture.CUDAarray(ch, *self.shape[::-1])
            self.res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=arr)
            self.tex = cp.cuda.texture.TextureDescriptor((cp.cuda.runtime.cudaAddressModeBorder,
                                                         cp.cuda.runtime.cudaAddressModeBorder,
                                                         cp.cuda.runtime.cudaAddressModeBorder),
                                                         cp.cuda.runtime.cudaFilterModeLinear,
                                                         cp.cuda.runtime.cudaReadModeElementType)
            self.tex_obj = cp.cuda.texture.TextureObject(self.res, self.tex)

            # prefilter if required and upload to texture
            if interpolation.startswith('filt_bspline'):
                data = _bspline_prefilter(data)
            arr.copy_from(data)

            # workgroup dims
            self.dim_grid, self.dim_blocks = compute_elementwise_launch_dims(self.shape)

            del data

        elif device == 'cpu':
            self.shape = data.shape
            self.data = data

    def affine(self, transform_m: np.ndarray, profile: bool = False, output: cp.ndarray = None) -> Union[np.ndarray, None]:

        if self.device.startswith('gpu'):

            if profile:
                stream = cp.cuda.Stream.null
                t_start = stream.record()

            # kernel setup
            xform = cp.asarray(transform_m)

            if output is None:
                output_vol = cp.zeros(self.shape, dtype=self.d_type)
            else:
                output_vol = output
            # launch
            self.affine_kernel(self.dim_grid, self.dim_blocks, (output_vol, self.tex_obj, xform, self.d_shape))

            if profile:
                t_end = stream.record()
                t_end.synchronize()

                time_took = cp.cuda.get_elapsed_time(t_start, t_end)
                print(f'transform finished in {time_took:.3f}ms')

            del xform
            if output is None:
                return output_vol
            else:
                return None

        elif self.device == 'cpu':
            return affine(self.data, transform_m, self.interpolation, profile, output, self.device)

    def transform(self, scale: Union[float, Tuple[float, float, float], np.ndarray] = None,
                  shear: Union[float, Tuple[float, float, float], np.ndarray] = None,
                  rotation: Union[Tuple[float, float, float], np.ndarray] = None,
                  axisrotation: Union[Tuple[float, float, float], np.ndarray] = None,
                  rotation_units: str = 'deg', rotation_order: str = 'rzxz',
                  translation: Union[Tuple[float, float, float], np.ndarray] = None,
                  translation2: Union[Tuple[float, float, float], np.ndarray] = None,
                  center: Union[Tuple[float, float, float], np.ndarray] = None,
                  profile: bool = False, output = None, matrix: np.ndarray = None) -> Union[np.ndarray, None]:

        if center is None:
            center = np.divide(self.shape, 2, dtype=np.float32)
        # passing just one float is uniform scaling
        if isinstance(scale, float):
            scale = (scale, scale, scale)
        if isinstance(shear, float):
            shear = (shear, shear, shear)

        if matrix is None:
            matrix = transform_matrix(scale, shear, rotation, axisrotation, rotation_units, rotation_order, translation,
                                      translation2, center)
        return self.affine(matrix, profile, output)

    def translate(self,
                  translation: Tuple[float, float, float],
                  profile: bool = False, output = None) -> Union[np.ndarray, None]:

        m = translation_matrix(translation)
        return self.affine(m, profile, output)

    def shear(self,
              coefficients: Union[float, Tuple[float, float, float]],
              profile: bool = False, output = None) -> Union[np.ndarray, None]:

        # passing just one float is uniform scaling
        if isinstance(coefficients, float):
            coefficients = (coefficients, coefficients, coefficients)

        m = shear_matrix(coefficients)
        return self.affine(m, profile, output)

    def scale(self,
              coefficients: Union[float, Tuple[float, float, float]],
              profile: bool = False, output = None) -> Union[np.ndarray, None]:

        # passing just one float is uniform scaling
        if isinstance(coefficients, float):
            coefficients = (coefficients, coefficients, coefficients)

        m = scale_matrix(coefficients)
        return self.affine(m, profile, output)

    def rotate(self,
               rotation: Tuple[float, float, float],
               rotation_units: str = 'deg',
               rotation_order: str = 'rzxz',
               profile: bool = False, output = None) -> Union[np.ndarray, None]:

        m = rotation_matrix(rotation=rotation, rotation_units=rotation_units, rotation_order=rotation_order)
        return self.affine(m, profile, output)
