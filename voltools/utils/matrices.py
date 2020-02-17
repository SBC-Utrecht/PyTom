import numpy as np
from typing import Tuple, Union

# Rotation routines are heavily based on Christoph Gohlike's transformations.py
# axis sequences for Euler angles
_NEXT_AXIS = [1, 2, 0, 1]
# map axes strings to/from tuples of inner axis, parity, repetition, frame
_AXES2TUPLE = {
    'sxyz': (0, 0, 0, 0), 'sxyx': (0, 0, 1, 0), 'sxzy': (0, 1, 0, 0),
    'sxzx': (0, 1, 1, 0), 'syzx': (1, 0, 0, 0), 'syzy': (1, 0, 1, 0),
    'syxz': (1, 1, 0, 0), 'syxy': (1, 1, 1, 0), 'szxy': (2, 0, 0, 0),
    'szxz': (2, 0, 1, 0), 'szyx': (2, 1, 0, 0), 'szyz': (2, 1, 1, 0),
    'rzyx': (0, 0, 0, 1), 'rxyx': (0, 0, 1, 1), 'ryzx': (0, 1, 0, 1),
    'rxzx': (0, 1, 1, 1), 'rxzy': (1, 0, 0, 1), 'ryzy': (1, 0, 1, 1),
    'rzxy': (1, 1, 0, 1), 'ryxy': (1, 1, 1, 1), 'ryxz': (2, 0, 0, 1),
    'rzxz': (2, 0, 1, 1), 'rxyz': (2, 1, 0, 1), 'rzyz': (2, 1, 1, 1)}
_TUPLE2AXES = dict((v, k) for k, v in _AXES2TUPLE.items())
AVAILABLE_ROTATIONS = list(_AXES2TUPLE.keys())
AVAILABLE_UNITS = ['rad', 'deg']


def translation_matrix(translation: Union[Tuple[float, float, float], np.ndarray],
                       dtype: np.dtype = np.float32) -> np.ndarray:
    m = np.identity(4, dtype=dtype)
    m[:3, 3] = translation[:3]
    m[:3, 3] = -1 * m[:3, 3]
    return m


def rotation_matrix(rotation: Union[Tuple[float, float, float], np.ndarray],
                    rotation_units: str = 'deg', rotation_order: str = 'rzxz',
                    dtype: np.dtype = np.float32) -> np.ndarray:

    if rotation_units not in AVAILABLE_UNITS:
        raise ValueError(f'Rotation units must be one of {AVAILABLE_UNITS}')

    if rotation_order not in AVAILABLE_ROTATIONS:
        raise ValueError(f'Rotation order must be one of {AVAILABLE_ROTATIONS}')

    # Units conversion
    if rotation_units == 'deg':
        ai, aj, ak = np.deg2rad(rotation)[:3]
    else:
        ai, aj, ak = rotation[:3]

    # CCW notation
    ai, aj, ak = -1 * ai, -1 * aj, -1 * ak

    # General rotation calculations
    firstaxis, parity, repetition, frame = _AXES2TUPLE[rotation_order]

    i = firstaxis
    j = _NEXT_AXIS[i+parity]
    k = _NEXT_AXIS[i-parity+1]

    if frame:
        ai, ak = ak, ai
    if parity:
        ai, aj, ak = -ai, -aj, -ak

    si, sj, sk = np.sin([ai, aj, ak])
    ci, cj, ck = np.cos([ai, aj, ak])
    cc, cs = ci*ck, ci*sk
    sc, ss = si*ck, si*sk

    # Matrix assembly
    m = np.identity(4, dtype=dtype)

    if repetition:
        m[i, i] = cj
        m[i, j] = sj*si
        m[i, k] = sj*ci
        m[j, i] = sj*sk
        m[j, j] = -cj*ss+cc
        m[j, k] = -cj*cs-sc
        m[k, i] = -sj*ck
        m[k, j] = cj*sc+cs
        m[k, k] = cj*cc-ss
    else:
        m[i, i] = cj*ck
        m[i, j] = sj*sc-cs
        m[i, k] = sj*cc+ss
        m[j, i] = cj*sk
        m[j, j] = sj*ss+cc
        m[j, k] = sj*cs-sc
        m[k, i] = -sj
        m[k, j] = cj*si
        m[k, k] = cj*ci

    return m


def shear_matrix(coefficients: Union[Tuple[float, float, float], np.ndarray],
                 dtype: np.dtype = np.float32) -> np.ndarray:
    m = np.identity(4, dtype)
    m[1, 2] = coefficients[2]
    m[0, 2] = coefficients[1]
    m[0, 1] = coefficients[0]
    return m


def scale_matrix(coefficients: Union[Tuple[float, float, float], np.ndarray],
                 dtype: np.dtype = np.float32) -> np.ndarray:
    m = np.identity(4, dtype)
    m[0, 0] = coefficients[0]
    m[1, 1] = coefficients[1]
    m[2, 2] = coefficients[2]
    return m


def transform_matrix(scale: Union[Tuple[float, float, float], np.ndarray] = None,
                     shear: Union[Tuple[float, float, float], np.ndarray] = None,
                     rotation: Union[Tuple[float, float, float], np.ndarray] = None,
                     rotation_units: str = 'deg', rotation_order: str = 'rzxz',
                     translation: Union[Tuple[float, float, float], np.ndarray] = None,
                     center: Union[Tuple[float, float, float], np.ndarray] = None,
                     dtype: np.dtype = np.float32) -> np.ndarray:
    """
    Returns composed transformation matrix according to the passed arguments.
    Transformation order: scale, shear, rotation, translation
    If center is passed, the transformation is done around the center (shape/2)
    and the transformation order is: pre-translation, scale, shear, rotation, post-translation, translation
    """

    m = np.identity(4, dtype=dtype)

    # Translation
    if translation is not None:
        m = np.dot(m, translation_matrix(translation, dtype))

    # Post-translation
    if center is not None:
        m = np.dot(m, translation_matrix(tuple(-1 * i for i in center), dtype))

    # Rotation
    if rotation is not None:
        m = np.dot(m, rotation_matrix(rotation, rotation_units, rotation_order, dtype))

    # Shear
    if shear is not None:
        m = np.dot(m, shear_matrix(shear, dtype))

    # Scale
    if scale is not None:
        m = np.dot(m, scale_matrix(scale, dtype))

    # Pre-translation
    if center is not None:
        m = np.dot(m, translation_matrix(center, dtype))

    # Keep it homo, geneous
    m /= m[3, 3]

    return m
