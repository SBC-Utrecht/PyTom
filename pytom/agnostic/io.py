#! /usr/bin/env python
# -*- coding: utf-8 -*-

from pytom.gpu.initialize import xp, device
import numpy as np

def read(filename, ndarray=True, order='F', keepnumpy=False, deviceID=None):
    """Read EM file. Now only support read the type float32 on little-endian machines.

    @param filename: file name to read.

    @return The data from the EM file in ndarray
    """
    assert filename
    assert filename.split('.')[-1].lower() in ('em', 'mrc')

    if filename.endswith('.em'):
        data = read_em(filename, order=order, keepnumpy=keepnumpy, deviceID=deviceID)
    elif filename.endswith('.mrc') or filename.endswith('.rec'):
        data = read_mrc(filename, order=order, keepnumpy=keepnumpy, deviceID=deviceID)
    elif filename.endswith('.txt'):
        from pytom.basic.datatypes import DATATYPE_METAFILE, DATATYPE_ALIGNMENT_RESULTS, DATATYPE_MARKER_RESULTS, DATATYPE_MARKERFILE
        from pytom.basic.files import loadtxt

        success =False
        data = None

        for dtype in (DATATYPE_METAFILE, DATATYPE_ALIGNMENT_RESULTS, DATATYPE_MARKER_RESULTS, DATATYPE_MARKERFILE):

            try:
                data = loadtxt(filename, dtype=dtype)
                success = True
            except Exception as e:
                pass
            if success:
                return data

        if not success: raise Exception(f'No valid interpreter found for {filename}')
    else:
        raise Exception('Invalid filetype. Please provide a *.em or a *.mrc file, or a pytom txt file.')

    if ndarray:
        return data
    else:
        return n2v(data)

def read_mrc(filename, order='F', keepnumpy=False, deviceID=None):
    import numpy as np
    f = open(filename, 'r')
    try:
        dt_header = np.dtype('int32')
        header = np.fromfile(f, dt_header, 256)
        x = header[0]
        y = header[1]
        z = header[2]

        if header[23]:
            extended_header = np.fromfile(f, np.float32, header[23]//4)

        # read the data type
        dt = header[3]
        default_type = False
        if dt == 0:  # byte
            dt_data = np.dtype('int8')
        elif dt == 1:  # short
            dt_data = np.dtype('int16')
        elif dt == 2:  # long
            dt_data = np.dtype('float32') # little-endian, float32
            default_type = True
        elif dt == 4:  # float32
            dt_data = np.dtype('complex32')
        elif dt == 6:  # float complex
            dt_data = np.dtype('uint16')
        else:
            raise Exception("Data type not supported yet!")

        v = np.fromfile(f, dt_data, x * y * z)
    finally:
        f.close()


    if keepnumpy:
        volume = np.array(v.reshape((x, y, z), order=order), dtype='float32').copy()  # fortran-order array
    else:
        if not deviceID == None:
            id = int(deviceID.split(":")[1])
            xp.cuda.Device(id).use()
        volume = xp.array(v.reshape((x, y, z), order=order), dtype=xp.float32).copy()  # fortran-order array

    return volume

def read_em(filename, order='F', keepnumpy=False, deviceID=None):
    """Read EM file. Now only support read the type float32 on little-endian machines.

    @param filename: file name to read.

    @return The data from the EM file in ndarray
    """
    f = open(filename, 'r')

    try:
        dt_header =np.dtype('int32')
        header = np.fromfile(f, dt_header, 128)
        x = header[1]
        y = header[2]
        z = header[3]

        # read the data type
        dt = int(hex(header[0])[2])
        default_type = False

        if dt == 1: # byte
            raise Exception("Data type not supported yet!")
        elif dt == 2: # short
            dt_data = np.dtype('<i2')
        elif dt == 4: # long
            dt_data = np.dtype('<i4')
        elif dt == 5: # float32
            dt_data = np.dtype('<f4') # little-endian, float32
            default_type = True
        elif dt == 8: # float complex
            raise Exception("Data type not supported yet!")
        elif dt == 9: # double
            dt_data = np.dtype('<f8')
        elif dt == 10: # double complex
            raise Exception("Data type not supported yet!")
        else:
            raise Exception("Data type not supported yet!")

        v = np.fromfile(f, dt_data, x*y*z)
    finally:
        f.close()

    if keepnumpy:
        volume = np.array(v.reshape((x, y, z), order=order), dtype='float32').copy() # fortran-order array
    else: # if the input data is not the default type, convert
        if not deviceID == None:
            id = int(deviceID.split(":")[1])
            xp.cuda.Device(id).use()
        volume = xp.array(v.reshape((x, y, z), order=order), dtype='float32').copy() # fortran-order array

    return volume

def write(filename, data, tilt_angle=0, pixel_size=1, order='F'):
    """Write EM or MRC file. Now only support written in type float32 on little-endian machines.

    @param filename: file name.
    @param data: data which should be written to file. Can be ndarray or pytom volume
    @param tilt_angle: if data is a projection, this is the tilt angle of the projection.
    @param ndarray: is the data a ndarray or not (if not it is assumed to be a pytom volume).
    @param data: data to write.
    """
    assert filename
    assert filename.split('.')[-1].lower() in ('em', 'mrc')


    try:
        data = data.get()
    except:
        pass

    if not type(data) == type(np.array((1))):
        from pytom_numpy import vol2npy
        try: data = vol2npy(data).copy()
        except: raise Exception('Invalid data type of data. Please provide an ndarray or a pytom volume object.')

    if data.dtype != np.dtype("float32"):
        data = data.astype(np.float32)

    try:
        data = data.get()
    except:
        pass

    if filename.endswith('.em'):
        write_em(filename, data, tilt_angle, order=order)
    elif filename.endswith('.mrc'):
        write_mrc(filename, data, tilt_angle, pixel_size=pixel_size, order=order)
    else:
        raise Exception('Unsupported file format, cannot write an {}-file'.format(filename.split('.')[-1]))

def binary_string(values, type):
    import numpy as np
    return np.array(values, type).tostring()

def write_mrc(filename, data, tilt_angle=0, pixel_size=1, inplanerot=0, magnification=1., dx=0., dy=0., current_tilt_angle=999,
              order='F'):
    import numpy as np
    try:
        data = data.get()
    except:
        pass
    if data.dtype != np.dtype('float32'): # if the origin data type is not float32, convert
        data = data.astype(np.float32)

    assert len(data.shape) < 4 and len(data.shape) > 0

    size = np.ones((3),dtype=np.int32)
    size[:len(data.shape)] = data.shape

    mode = 2
    cell_angles = (0, 0, 0)
    cell_size= pixel_size*size
    dmin, dmax, dmean = data.min(), data.max(), data.mean()

    if current_tilt_angle == 999:
        current_tilt_angle=tilt_angle

    extra = [0,0,0,20140] + [0]*13 + [dx, dy, tilt_angle, inplanerot, magnification, current_tilt_angle, 0, 0] # set tiltangle at 43th byte

    if np.little_endian:
        machst = 0x00004144
    else:
        machst = 0x11110000

    rms = data.std()

    strings = [
        binary_string(size, np.int32),  # nx, ny, nz
        binary_string(mode, np.int32),  # mode
        binary_string((0, 0, 0), np.int32),  # nxstart, nystart, nzstart
        binary_string(size, np.int32),  # mx, my, mz
        binary_string(cell_size, np.float32),  # cella
        binary_string(cell_angles, np.int32),  # cellb
        binary_string((1, 2, 3), np.int32),  # mapc, mapr, maps
        binary_string((dmin, dmax, dmean), np.float32),  # dmin, dmax, dmean
        binary_string(0, np.int32),  # ispg
        binary_string(0, np.int32),  # nsymbt
        binary_string(extra, np.float32),  # extra
        binary_string((0., 0., 0.), np.int32),  # origin
        binary_string(0, np.int32),  # map
        binary_string(machst, np.int32),  # machst
        binary_string(rms, np.float32),  # rms
        binary_string(0, np.int32),  # nlabl
        binary_string([0]*200, np.int32),
    ]

    header = b"".join(strings)

    f = open(filename, 'wb')
    try:
        f.write(header)
        f.write(data.tostring(order=order)) # fortran-order array
    finally:
        f.close()

def write_em(filename, data, tilt_angle=0,order='F'):
    """Write EM file. Now only support written in type float32 on little-endian machines.

    @param filename: file name.
    @param data: data to write.
    """
    if data.dtype != np.dtype('float32'): # if the origin data type is not float32, convert
        data = np.array(data, dtype='float32')
    
    header = np.zeros(128, dtype='int32')
    header[0] = 83886086 # '0x5000006', TODO: hard-coded, to be changed!
    header[24+18] = int(tilt_angle*1000) # set tilt angle

    if len(data.shape) == 3:
        header[1:4] = data.shape
    elif len(data.shape) == 2:
        header[1:3] = data.shape
        header[3] = 1
    elif len(data.shape) == 1:
        header[1] = data.shape[0]
        header[2:4] = 1
    else:
        raise Exception("Input data shape invalid!")

    f = open(filename, 'wb')
    try:
        f.write(header.tostring())
        f.write(data.tostring(order=order)) # fortran-order array
    finally:
        f.close()

def read_size(filename, dim=''):
    emfile = filename.endswith('.em')*1
    f = open(filename, 'r')
    try:
        dt_header = np.dtype('int32')
        header = np.fromfile(f, dt_header, 4)
        x = header[0+emfile]
        y = header[1+emfile]
        z = header[2+emfile]
    except:
        raise Exception("reading of MRC file failed")

    f.close()

    if dim=='x':
        return x
    elif dim=='y':
        return y
    elif dim=='z':
        return z
    else:
        return [x,y,z]

def read_header(filename):
    emfile = filename.endswith('.em') * 1
    f = open(filename, 'rb')
    try:
        if emfile:
            header_data = np.fromfile(f, np.dtype('int32'), 128)
        else:
            """
            header_data = []

            header_data += list(np.fromfile(f, np.dtype('int32'), 10))
            header_data += list(np.fromfile(f, np.dtype('float32'), 3))
            header_data += list(np.fromfile(f, np.dtype('int32'), 6))
            header_data += list(np.fromfile(f, np.dtype('float32'), 3))
            header_data += list(np.fromfile(f, np.dtype('int32'), 2))
            header_data += list(np.fromfile(f, np.dtype('float32'), 25))
            header_data += list(np.fromfile(f, np.dtype('int32'), 5))
            header_data += list(np.fromfile(f, np.dtype('float32'), 1))
            header_data += list(np.fromfile(f, np.dtype('int32'), 201))


            print(header_data)
            """
            header_data = np.fromfile(f, np.dtype('float32'), 256)
            #"""
    finally:
        f.close()

    return header_data

def read_tilt_angle(filename):
    """Reads the EM header only.
    @param filename: The em file
    """
    emfile = filename.endswith('.em') * 1
    f = open(filename, 'rb')
    try:
        if emfile:
            header_data = np.fromfile(f, np.dtype('int32'), 128)
            tilt_angle = header_data[42]/1000.
        else:
            header_data = np.fromfile(f, np.dtype('float32'), 256)
            tilt_angle = header_data[43]
    finally:
        f.close()

    return tilt_angle

def read_rotation_angles(filename):
    f = open(filename, 'rb')
    try:
        header_data = np.fromfile(f, np.dtype('float32'), 256)
        z1,x,z2 = header_data[43:46]
    finally:
        f.close()

    return [z1,x,z2]

def read_pixelsize(filename, dim=''):
    import sys
    import numpy as np

    f = open(filename, 'rb')
    try:
            header_data = []

            header_data += list(np.fromfile(f, np.dtype('int32'), 10))
            header_data += list(np.fromfile(f, np.dtype('float32'), 3))
            header_data += list(np.fromfile(f, np.dtype('int32'), 6))
            header_data += list(np.fromfile(f, np.dtype('float32'), 3))
            header_data += list(np.fromfile(f, np.dtype('int32'), 2))
            header_data += list(np.fromfile(f, np.dtype('float32'), 25))
            header_data += list(np.fromfile(f, np.dtype('int32'), 5))
            header_data += list(np.fromfile(f, np.dtype('float32'), 1))
            header_data += list(np.fromfile(f, np.dtype('int32'), 201))

    finally:
        f.close()

    if filename.endswith('.em'):
        sx,sy,sz = header_data[1:4]
    else:
        sx,sy,sz = header_data[:3]

    x,y,z = np.array(header_data[10:13])/np.array((sx,sy,sz), np.float32)
    if dim == 'x':
        return x
    elif dim == 'y':
        return y
    elif dim == 'z':
        return z
    else:
        return [x, y, z]


def write_rotation_angles(filename, z1=0, z2=0, x=0):
    assert filename.endswith('.mrc')
    header = read_header(filename)
    header[43:46] = [z1, x, z2]
    data = read(filename)
    f = open(filename, 'wb')
    try:
        f.write(header.tostring())
        f.write(data.tostring(order='F')) # fortran-order array
    finally:
        f.close()

def write_tilt_angle(filename, tilt_angle):
    emfile = filename.endswith('.em') * 1
    header = read_header(filename)
    if emfile:
        header[42] = int(round(tilt_angle*1000))
    else:
        header[43] = tilt_angle
    data = read(filename)

    f = open(filename, 'wb')
    try:
        f.write(tostring(header))
        f.write(data.tostring(order='F')) # fortran-order array
    finally:
        f.close()

def tostring(header):
    strings = []
    for value in header:

        strings.append(binary_string(value,type(value)))

    header = b"".join(strings)
    return header

def n2v(data):
    """Transfer Numpy array into a Pytom volume.
    Note the data will be shared between Numpy and Pytom!

    @param data: data to convert.
    """
    try:
        from pytom_volume import vol
        from pytom_numpy import npy2vol
    except:
        raise ImportError("Pytom library is not installed or set properly!")

    if data.dtype != np.dtype("float32"):
        data = np.array(data, dtype="float32")
    
    if len(data.shape) == 3:
        if np.isfortran(data):
            # Fortran order
            v = npy2vol(data, 3)
        else:
            vv = np.asfortranarray(data)
            v = npy2vol(vv, 3)
    elif len(data.shape) == 2:
        if np.isfortran(data):
            data = data.reshape(data.shape[0], data.shape[1], 1)
            v = npy2vol(data, 2)
        else:
            data = data.reshape(data.shape[0], data.shape[1], 1)
            vv = np.asfortranarray(data)
            v = npy2vol(vv, 2)
    else:
        raise Exception("Data shape invalid!")

    return v

def readSubvolumeFromFourierspaceFile(filename, sizeX, sizeY, sizeZ):
    """
    readSubvolumeFromFourierspaceFile: This function is required when data \
    (in real space) is read in binned mode and a related fourier space file
    like a wedge needs to be read alongside.
    Works only if fourier file is reduced complex without any shift applied.
    @param filename: The fourier space file name
    @param sizeX: X final size of subvolume if it was complete
    (what L{pytom.basic.structures.Wedge.returnWedgeVolume} with
    humanUnderstandable == True returns)
    @param sizeY: Y final size of subvolume if it was complete
    (what L{pytom.basic.structures.Wedge.returnWedgeVolume}
    with humanUnderstandable == True returns)
    @param sizeZ: Z final size of subvolume if it was complete
    (what L{pytom.basic.structures.Wedge.returnWedgeVolume}
    with humanUnderstandable == True returns)
    @return: A subvolume
    @author: Thomas Hrabe
    """
    from pytom.agnostic.io import read
    from pytom.basic.fourier import fourierSizeOperation
    [newX, newY, newZ] = fourierSizeOperation(sizeX, sizeY, sizeZ,
                                              reducedToFull=False)


    if filename.__class__ == str:
        originalVolume = read(filename)
    elif filename.__class__ == xp.array:
        # open a backdoor for this function to take volumes, but
        # this should be rather an exception -> not fully documented
        originalVolume = filename
    else:
        raise TypeError('Filename must be a string')

    originalVolume = xp.fft.fftshift(originalVolume,axes=(0,1))
    newVolume = originalVolume[sizeX//2-newX//2:sizeX//2+newX//2, sizeY//2-newY//2:sizeY//2+newY//2,:newZ]
    newVolume = xp.fft.fftshift(newVolume,axes=(0,1))


    return newVolume