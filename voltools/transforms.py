import time
import numpy as np
from scipy.ndimage import affine_transform
from typing import Union, Tuple
from pathlib import Path

from .utils import scale_matrix, shear_matrix, rotation_matrix, translation_matrix, transform_matrix, \
    compute_prefilter_workgroup_dims, compute_elementwise_launch_dims, get_available_devices, switch_to_device


_INTERPOLATIONS = {
    'linear': 'linearTex3D',
    'bspline': 'cubicTex3D',
    'bspline_simple': 'cubicTex3DSimple',
    'filt_bspline': 'cubicTex3D',
    'filt_bspline_simple': 'cubicTex3DSimple'
}
AVAILABLE_INTERPOLATIONS = list(_INTERPOLATIONS.keys())
AVAILABLE_DEVICES = get_available_devices()

if 'gpu' in AVAILABLE_DEVICES:
    import cupy as cp


def transform(volume: np.ndarray,
              scale: Union[float, Tuple[float, float, float], np.ndarray] = None,
              shear: Union[float, Tuple[float, float, float], np.ndarray] = None,
              rotation: Union[Tuple[float, float, float], np.ndarray] = None,
              rotation_units: str = 'deg', rotation_order: str = 'rzxz',
              translation: Union[Tuple[float, float, float], np.ndarray] = None,
              center: Union[Tuple[float, float, float], np.ndarray] = None,
              interpolation: str = 'linear',
              profile: bool = False,
              output = None, device: str = 'cpu'):

    if center is None:
        center = np.divide(volume.shape, 2, dtype=np.float32)

    # passing just one float is uniform scaling
    if isinstance(scale, float):
        scale = (scale, scale, scale)
    if isinstance(shear, float):
        shear = (shear, shear, shear)

    m = transform_matrix(scale, shear, rotation, rotation_units, rotation_order, translation, center)
    return affine(volume, m, interpolation, profile, output, device)


def translate(volume: np.ndarray,
              translation: Tuple[float, float, float],
              interpolation: str = 'linear',
              profile: bool = False,
              output = None, device: str = 'cpu'):

    m = translation_matrix(translation)
    return affine(volume, m, interpolation, profile, output, device)


def shear(volume: np.ndarray,
          coefficients: Union[float, Tuple[float, float, float]],
          interpolation: str = 'linear',
          profile: bool = False,
          output = None, device: str = 'cpu'):

    # passing just one float is uniform scaling
    if isinstance(coefficients, float):
        coefficients = (coefficients, coefficients, coefficients)

    m = shear_matrix(coefficients)
    return affine(volume, m, interpolation, profile, output, device)


def scale(volume: np.ndarray,
          coefficients: Union[float, Tuple[float, float, float]],
          interpolation: str = 'linear',
          profile: bool = False,
          output = None, device: str = 'cpu'):

    # passing just one float is uniform scaling
    if isinstance(coefficients, float):
        coefficients = (coefficients, coefficients, coefficients)

    m = scale_matrix(coefficients)
    return affine(volume, m, interpolation, profile, output, device)


def rotate(volume: np.ndarray,
           rotation: Tuple[float, float, float],
           rotation_units: str = 'deg',
           rotation_order: str = 'rzxz',
           interpolation: str = 'linear',
           profile: bool = False,
           output = None, device: str = 'cpu'):

    m = rotation_matrix(rotation=rotation, rotation_units=rotation_units, rotation_order=rotation_order)
    return affine(volume, m, interpolation, profile, output, device)


def affine(volume: np.ndarray,
           transform_m: np.ndarray,
           interpolation: str = 'linear',
           profile: bool = False,
           output = None,
           device: str = 'cpu'):

    if device not in AVAILABLE_DEVICES:
        raise ValueError(f'Unknown device ({device}), must be one of {AVAILABLE_DEVICES}')

    if device == 'cpu':

        if profile:
            t_start = time.time()

        # set parameters for scipy affine transform
        if interpolation == 'linear':
            order = 1
        else:
            order = 3

        if not interpolation.startswith('filt_bspline'):
            prefilter = False
        else:
            prefilter = True

        # run affine transformation
        output_vol = affine_transform(volume, transform_m, output=output, order=order, prefilter=prefilter)

        if profile:
            t_end = time.time()
            time_took = (t_end - t_start) * 1000
            print(f'transform finished in {time_took:.3f}ms')

        if output is not None:
            return output
        else:
            return output_vol

    elif device.startswith('gpu'):
        switch_to_device(device)

        if profile:
            stream = cp.cuda.Stream.null
            t_start = stream.record()

        volume = cp.asarray(volume)
        volume_shape = volume.shape

        # texture setup
        ch = cp.cuda.texture.ChannelFormatDescriptor(32, 0, 0, 0, cp.cuda.runtime.cudaChannelFormatKindFloat)
        arr = cp.cuda.texture.CUDAarray(ch, *volume_shape[::-1])  # CUDAArray: last dimension=fastest changing dimension
        res = cp.cuda.texture.ResourceDescriptor(cp.cuda.runtime.cudaResourceTypeArray, cuArr=arr)
        tex = cp.cuda.texture.TextureDescriptor((cp.cuda.runtime.cudaAddressModeBorder,
                                                 cp.cuda.runtime.cudaAddressModeBorder,
                                                 cp.cuda.runtime.cudaAddressModeBorder),
                                                cp.cuda.runtime.cudaFilterModeLinear,
                                                cp.cuda.runtime.cudaReadModeElementType)
        texobj = cp.cuda.texture.TextureObject(res, tex)

        # prefilter if required and upload to texture
        if interpolation.startswith('filt_bspline'):
            volume = _bspline_prefilter(volume)
            arr.copy_from(volume)
        else:
            arr.copy_from(volume)

        # kernel setup
        kernel = _get_transform_kernel(interpolation)
        dims = cp.asarray(volume_shape, dtype=cp.uint32)
        xform = cp.asarray(transform_m)
        dim_grid, dim_blocks = compute_elementwise_launch_dims(volume_shape)

        if output is None:
            volume.fill(0.0)  # reuse input array
        else:
            volume = output

        kernel(dim_grid, dim_blocks, (volume, texobj, xform, dims))

        if profile:
            t_end = stream.record()
            t_end.synchronize()

            time_took = cp.cuda.get_elapsed_time(t_start, t_end)
            print(f'transform finished in {time_took:.3f}ms')

        if output is None:
            del texobj, xform, dims
            return volume.get()
        else:
            del texobj, xform, dims
            return None

    else:
        raise ValueError(f'No instructions for {device}.')

def _get_transform_kernel(interpolation: str = 'linear'):

    if interpolation not in AVAILABLE_INTERPOLATIONS:
        raise ValueError(f'Interpolation must be one of {interpolation}')

    code = f'''
        #include "helper_math.h"
        #include "bspline.h"
        #include "helper_interpolation.h"
        
        extern "C" {{
            inline __device__ int get_z_idx(const int i, const uint4* const dims) {{
                return i / (dims[0].y * dims[0].z);
            }}
            inline __device__ int get_y_idx(const int i, const uint4* const dims) {{
                return (i % (dims[0].y * dims[0].z)) / dims[0].z;
            }}
            inline __device__ int get_x_idx(const int i, const uint4* const dims) {{
                return (i % (dims[0].y * dims[0].z)) % dims[0].z;
            }}

            __global__ void transform(float* const volume,
                                      cudaTextureObject_t texture,
                                      const float4* const xform,
                                      const uint4* const dims) {{
            
                unsigned tid = threadIdx.x;
                unsigned total_threads = gridDim.x*blockDim.x;
                unsigned cta_start = blockDim.x*blockIdx.x;
                unsigned i;
                unsigned n = dims[0].x * dims[0].y * dims[0].z;
                
                for (i = cta_start + tid; i < n; i += total_threads) {{
                    int z = get_x_idx(i, dims);
                    int y = get_y_idx(i, dims);
                    int x = get_z_idx(i, dims);
                    
                    float4 voxf = make_float4(((float)x), ((float)y), ((float)z), 1.0f);
                    
                    float3 ndx;
                    ndx.z = dot(voxf, xform[0]) + .5f;
                    ndx.y = dot(voxf, xform[1]) + .5f;
                    ndx.x = dot(voxf, xform[2]) + .5f;
                    
                    if (ndx.x < 0 || ndx.y < 0 || ndx.z < 0 || ndx.x >= dims[0].z || ndx.y >= dims[0].y || ndx.z >= dims[0].x) {{
                        continue;
                    }}
                    
                    volume[i] = {_INTERPOLATIONS[interpolation]}(texture, ndx);
                }}
            }}
        }}
    '''
    incl_path = str((Path(__file__).parent / 'kernels').absolute())
    kernel = cp.RawKernel(code=code, name='transform', options=('-I', incl_path))
    return kernel

def _bspline_prefilter(volume):

    code = f'''
        #include "helper_math.h"
        #include "bspline.h"
    '''
    incl_path = str((Path(__file__).parent / 'kernels').absolute())
    prefilter_x = cp.RawKernel(code=code, name='SamplesToCoefficients3DX', options=('-I', incl_path))
    prefilter_y = cp.RawKernel(code=code, name='SamplesToCoefficients3DY', options=('-I', incl_path))
    prefilter_z = cp.RawKernel(code=code, name='SamplesToCoefficients3DZ', options=('-I', incl_path))

    slice_stride = volume.strides[1]
    dim_grid, dim_block = compute_prefilter_workgroup_dims(volume.shape)
    dims = cp.asarray(volume.shape[::-1], dtype=cp.int32)

    prefilter_x(dim_grid[0], dim_block[0], (volume, slice_stride, dims))
    prefilter_y(dim_grid[1], dim_block[1], (volume, slice_stride, dims))
    prefilter_z(dim_grid[2], dim_block[2], (volume, slice_stride, dims))

    return volume
