
import numpy
from numpy import float32
from numpy import floor, ceil, newaxis, zeros, rot90

def read_mrc(mrc_fname,binning=[1,1,1]):
    a = MRC_Data(mrc_fname,'mrc')
    
    [dz,dy,dx] = a.data_size #420,420,420                                                                                   
    input_image = a.read_matrix((0,0,0),(dz,dy,dx),(1,1,1),None).copy()
    if input_image.shape[0] == 1:
        input_image = input_image[0]
    if binning[0] >1:
        input_image = downsample(input_image,binning[0])
    return input_image

def add_angle(mrc_fname, angle_tomo, single=True):
    from npy2mrc import write_mrc2000_grid_data, Array_Grid_Data
    if single:  mrc_fname, angle_tomo = [mrc_fname],[angle_tomo]
    for fname, ang in zip(mrc_fname, angle_tomo):
        a = MRC_Data(fname,'mrc')
        [dz,dy,dx] = a.data_size 
        write_mrc2000_grid_data( Array_Grid_Data(a.read_matrix((0,0,0),(dz,dy,dx),(1,1,1),None).copy()), fname, angle_tomo=angle)
    
def read_angle(mrc_fname, single=True,extra='user'):
    if single:
        mrc_fname = [mrc_fname]
    angles = []
    for fname in mrc_fname:
        angles.append( float(read_mrc_header(fname)[extra][0])/1000.)
    if single: angles = angles[0]
    return angles
    
def read_mrc_header(mrc_fname):
    a = MRC_Data(mrc_fname,'mrc')
    print(a)
    return a.header

def square_mrc(mrc_fname):
    i = read_mrc(mrc_fname)

    if len(i.shape)==3:
        o = i[i.shape[0],:min(i.shape[1:]),:min(i.shape[1:])]
    elif len(i.shape)==2:
        o = i[:min(i.shape),:min(i.shape)]
    else:
        raise Exception('no valid dimension mrcfile')
    
    convert_numpy_array3d_mrc(o[newaxis,:,:],mrc_fname)


def transpose_mrc_2d(mrc_fname):
    '''Transposes a mrc file'''
    i = read_mrc(mrc_fname)
    
    if len(i.shape)==2:
        convert_numpy_array3d_mrc((i.T)[newaxis,:,:],mrc_fname)
    elif len(i.shape)==3:
        dimz,dimy,dimx = i.shape
        o = zeros((dimz,dimx,dimy),dtype=float)
        for z in range(dimz):
            o[z,:,:] = i[z,:,:].T
        convert_numpy_array3d_mrc(o,mrc_fname)
    else:
        print('wrong dimensions')

def rot90_mrc_2d(mrc_fname,times=1):
    i = read_mrc(mrc_fname)
    o = rot90(i,times)
    convert_numpy_array3d_mrc((o)[newaxis,:,:],mrc_fname)

def downsample(img, factor):
    ory,orx = img.shape
    offsetx, offsety = crop_pix(img,factor)
    offsetx = int(offsetx)
    offsety = int(offsety)

    starty,endy,startx,endx = int(offsety//2), int( -(offsety//2+offsety%2)+((not offsety)*ory)),int(offsetx//2),\
                              int(-(offsetx//2+offsetx%2)+((not offsetx)*orx))
    img = img[starty:endy,startx:endx]

    image_size = img.shape[0]
    
    ds = (image_size//factor)
    ds2 =  img.shape[1]//factor
    return img.reshape(ds,image_size//ds,ds2,image_size//ds).mean(-1).mean(1)

def crop_pix(img,factor):
    ydim,xdim = img.shape
    offsetx,offsety = xdim%factor,ydim%factor
    #if ydim > xdim:
    #    offsety += ydim-xdim
    #if ydim < xdim:
    #    offsetx += xdim-ydim
    return offsetx, offsety



# -----------------------------------------------------------------------------
# Read MRC or CCP4 map file format electron microscope data.
# Byte swapping will be done if needed.
#

# -----------------------------------------------------------------------------
# file_type can be 'mrc' or 'ccp4' or 'imod'.
#

from numpy import float32

def allocate_array(size, value_type = float32, step = None, progress = None,
                   reverse_indices = True, zero_fill = False):

    if step is None:
        msize = size
    else:
        msize = [1+(sz-1)/st for sz,st in zip(size, step)]

    shape = list(msize)
    if reverse_indices:
        shape.reverse()

    if zero_fill:
        from numpy import zeros as alloc
    else:
        from numpy import empty as alloc

    try:
        m = alloc(shape, value_type)
    except ValueError:
        # numpy 1.0.3 sometimes gives ValueError, sometimes MemoryError
        report_memory_error(msize, value_type)
    except MemoryError:
        report_memory_error(msize, value_type)

    if progress:
        progress.array_size(msize, m.itemsize)

    return m



# -----------------------------------------------------------------------------
# Read an array from a binary file making at most one copy of array in memory.
#
def read_full_array(path, byte_offset, size, type, byte_swap,
                    progress = None, block_size = 2**20):

    a = allocate_array(size, type)
    
    file = open(path, 'rb')
    file.seek(byte_offset)

    if progress:
        progress.close_on_cancel(file)
        a_1d = a.ravel()
        n = len(a_1d)
        nf = float(n)
        for s in range(0,n,block_size):
            b = a_1d[s:s+block_size]
            file.readinto(b)
            progress.fraction(s/nf)
        progress.done()
    else:
        file.readinto(a)
        
    file.close()

    if byte_swap:
        a.byteswap(True)

    return a

# -----------------------------------------------------------------------------
# Return axes corresponding to cell_angles given in degrees.
#
def skew_axes(cell_angles):

  if tuple(cell_angles) == (90,90,90):
    # Use exact trig functions for this common case.
    ca = cb = cg = c1 = 0
    sg = c2 = 1
  else:
    # Convert to radians
    from math import pi, sin, cos, sqrt
    alpha, beta, gamma = [a * pi / 180 for a in cell_angles]
    cg = cos(gamma)
    sg = sin(gamma)
    cb = cos(beta)
    ca = cos(alpha)
    c1 = (ca - cb*cg)/sg
    c2 = sqrt(1 - cb*cb - c1*c1)

  axes = ((1, 0, 0), (cg, sg, 0), (cb, c1, c2))
  return axes


def rotation_from_axis_angle(axis, angle):

  from chimera import Xform, Vector
  r = Xform.rotation(Vector(*axis), angle)
  m = r.getOpenGLMatrix()
  rm = ((m[0], m[4], m[8]),
        (m[1], m[5], m[9]),
        (m[2], m[6], m[10]))
  return rm

def multiply_matrices(*mlist):

  if len(mlist) == 2:
    m1, m2 = mlist
    p = [[0,0,0,0],
         [0,0,0,0],
         [0,0,0,0]]
    for r in range(3):
      for c in range(4):
        p[r][c] = m1[r][0]*m2[0][c] + m1[r][1]*m2[1][c] + m1[r][2]*m2[2][c]
        p[r][3] += m1[r][3]
    p = tuple(map(tuple, p))
  else:
    p = multiply_matrices(*mlist[1:])
    p = multiply_matrices(mlist[0], p)
  return p

def invert_matrix(tf):

    from numpy import array, zeros, float, linalg
    tf = array(tf)
    r = tf[:,:3]
    t = tf[:,3]
    tfinv = zeros((3,4), float)
    rinv = tfinv[:,:3]
    tinv = tfinv[:,3]
    from numpy.linalg import inv as matrix_inverse
    from numpy import dot as matrix_multiply
    rinv[:,:] = matrix_inverse(r)
    tinv[:] = matrix_multiply(rinv, -t)
    return tfinv

def read_array(path, byte_offset, ijk_origin, ijk_size, ijk_step,
               full_size, type, byte_swap, progress = None):

    if (tuple(ijk_origin) == (0,0,0) and
        tuple(ijk_size) == tuple(full_size) and
        tuple(ijk_step) == (1,1,1)):
        m = read_full_array(path, byte_offset, full_size,
                            type, byte_swap, progress)
        return m

    matrix = allocate_array(ijk_size, type, ijk_step, progress)

    file = open(path, 'rb')

    if progress:
        progress.close_on_cancel(file)
        
    # Seek in file to read needed 1d slices.
    io, jo, ko = ijk_origin
    isize, jsize, ksize = ijk_size
    istep, jstep, kstep = ijk_step
    element_size = matrix.itemsize
    jbytes = full_size[0] * element_size
    kbytes = full_size[1] * jbytes
    ibytes = isize * element_size
    ioffset = io * element_size
    from numpy import fromstring
    for k in range(ko, ko+ksize, kstep):
      if progress:
        progress.plane((k-ko)/kstep)
      kbase = byte_offset + k * kbytes
      for j in range(jo, jo+jsize, jstep):
        offset = kbase + j * jbytes + ioffset
        file.seek(offset)
        data = file.read(ibytes)
        slice = fromstring(data, type)
        matrix[(k-ko)/kstep,(j-jo)/jstep,:] = slice[::istep]

    file.close()

    if byte_swap:
      matrix.byteswap(True)

    return matrix


    
# -----------------------------------------------------------------------------
# Apply scaling and skewing transformations.
#
def scale_and_skew(ijk, step, cell_angles):
  xa, ya, za = skew_axes(cell_angles)

  i, j, k = ijk
  d0, d1, d2 = step
  x,y,z = i*d0, j*d1, k*d2

  xyz = tuple(x*xa[a]+y*ya[a]+z*za[a] for a in (0,1,2))
  return xyz


def read_array(path, byte_offset, ijk_origin, ijk_size, ijk_step,
               full_size, type, byte_swap, progress = None):

    if (tuple(ijk_origin) == (0,0,0) and
        tuple(ijk_size) == tuple(full_size) and
        tuple(ijk_step) == (1,1,1)):
        m = read_full_array(path, byte_offset, full_size,
                            type, byte_swap, progress)
        return m

    matrix = allocate_array(ijk_size, type, ijk_step, progress)

    file = open(path, 'rb')

    if progress:
        progress.close_on_cancel(file)
        
    # Seek in file to read needed 1d slices.
    io, jo, ko = ijk_origin
    isize, jsize, ksize = ijk_size
    istep, jstep, kstep = ijk_step
    element_size = matrix.itemsize
    jbytes = full_size[0] * element_size
    kbytes = full_size[1] * jbytes
    ibytes = isize * element_size
    ioffset = io * element_size
    from numpy import fromstring
    for k in range(ko, ko+ksize, kstep):
      if progress:
        progress.plane((k-ko)/kstep)
      kbase = byte_offset + k * kbytes
      for j in range(jo, jo+jsize, jstep):
        offset = kbase + j * jbytes + ioffset
        file.seek(offset)
        data = file.read(ibytes)
        slice = fromstring(data, type)
        matrix[(k-ko)/kstep,(j-jo)/jstep,:] = slice[::istep]

    file.close()

    if byte_swap:
      matrix.byteswap(True)

    return matrix


    
class MRC_Data:

  def __init__(self, path, file_type):

    self.path = path

    import os.path
    self.name = os.path.basename(path)
    
    file = open(path, 'rb')

    file.seek(0,2)                              # go to end of file
    file_size = file.tell()
    file.seek(0,0)                              # go to beginning of file

    # Infer file byte order from column axis size nc.  Requires nc < 2**16
    # Was using mode value but 0 is allowed and does not determine byte order.
    self.swap_bytes = 0
    from numpy import int32,float32,float64
    nc = self.read_values(file, int32, 1)
    self.swap_bytes = not (nc > 0 and nc < 65536)
    file.seek(0,0)

    v = self.read_header_values(file, file_size, file_type)

    if v.get('imodStamp') == 1146047817:
      unsigned_8_bit = (v['imodFlags'] & 0x1 == 0)
    else:
      unsigned_8_bit = (file_type == 'imod' or v['type'] == 'mrc')
    self.element_type = self.value_type(v['mode'], unsigned_8_bit)

    self.check_header_values(v, file_size, file)
    self.header = v             # For dumpmrc.py standalone program.
    
    self.data_offset = file.tell()
    file.close()
    
    # Axes permutation.
    # Names c,r,s refer to fast, medium, slow file matrix axes.
    # Names i,j,k refer to x,y,z spatial axes.
    mapc, mapr, maps = v['mapc'], v['mapr'], v['maps']
    if (1 in (mapc, mapr, maps) and
        2 in (mapc, mapr, maps) and
        3 in (mapc, mapr, maps)):
      crs_to_ijk = (mapc-1,mapr-1,maps-1)
      ijk_to_crs = [None,None,None]
      for a in range(3):
        ijk_to_crs[crs_to_ijk[a]] = a
    else:
      crs_to_ijk = ijk_to_crs = (0, 1, 2)
    self.crs_to_ijk = crs_to_ijk
    self.ijk_to_crs = ijk_to_crs

    crs_size = v['nc'], v['nr'], v['ns']
    self.matrix_size = [int(s) for s in crs_size]
    self.data_size = [int(crs_size[a]) for a in ijk_to_crs]

    self.unit_cell_size = mx, my, mz = v['mx'], v['my'], v['mz']
    xlen, ylen, zlen = v['xlen'], v['ylen'], v['zlen']
    if mx > 0 and my > 0 and mz > 0 and xlen > 0 and ylen > 0 and zlen > 0:
      self.data_step = (xlen/mx, ylen/my, zlen/mz)
    else:
      self.data_step = (1.0, 1.0, 1.0)

    alpha, beta, gamma = (v['alpha'], v['beta'], v['gamma'])
    if not valid_cell_angles(alpha, beta, gamma, path):
      alpha = beta = gamma = 90
    self.cell_angles = (alpha, beta, gamma)

    from math import isnan
    if (v['type'] == 'mrc2000' and
        (v['zorigin'] != 0 or v['xorigin'] != 0 or v['yorigin'] != 0) and
        not (isnan(v['xorigin']) or isnan(v['yorigin']) or isnan(v['zorigin']))):
      #
      # This is a new MRC 2000 format file.  The xyz origin header parameters
      # are used instead of using ncstart, nrstart nsstart for new style files,
      # provided the xyz origin specified is not zero.  It turns out the
      # xorigin, yorigin, zorigin values are zero in alot of new files while
      # the ncstart, nrstart, nsstart give the correct (non-zero) origin. So in
      # cases where the xyz origin parameters and older nrstart, ncstart,
      # nsstart parameters specify different origins the one that is non-zero
      # is preferred.  And if both are non-zero, the newer xorigin, yorigin,
      # zorigin are used.
      #
      self.data_origin = (v['xorigin'], v['yorigin'], v['zorigin'])
    else:
      crs_start = v['ncstart'], v['nrstart'], v['nsstart']
      ijk_start = [crs_start[a] for a in ijk_to_crs]
      # Check if ijk_start values appear to be uninitialized.
      limit = 10*max(max(mx,my,mz), max(self.data_size))
      if [s for s in ijk_start if abs(s) > limit]:
        self.data_origin = (0., 0., 0.)
      else:
	#from griddata import scale_and_skew
        self.data_origin = scale_and_skew(ijk_start, self.data_step,
                                          self.cell_angles)

    r = ((1,0,0),(0,1,0),(0,0,1))
    for lbl in v['labels']:
      if str(lbl).startswith('Chimera rotation: '):
        ax,ay,az,angle = map(float, lbl.rstrip('\0').split()[2:])
        r = rotation_from_axis_angle((ax,ay,az), angle)
    self.rotation = r
    
    self.min_intensity = v['amin']
    self.max_intensity = v['amax']

  # ---------------------------------------------------------------------------
  # Format derived from C header file mrc.h.
  #
  def read_header_values(self, file, file_size, file_type):

    MRC_USER = 29
    CCP4_USER = 15
    MRC_NUM_LABELS = 10
    MRC_LABEL_SIZE = 80
    MRC_HEADER_LENGTH = 1024

    from numpy import int32, float32
    i32 = int32
    f32 = float32
    
    v = {}
    v['nc'], v['nr'], v['ns'] = self.read_values(file, i32, 3)
    v['mode'] = self.read_values(file, i32, 1)
    v['ncstart'], v['nrstart'], v['nsstart'] = self.read_values(file, i32, 3)
    v['mx'], v['my'], v['mz'] = self.read_values(file, i32, 3)
    v['xlen'], v['ylen'], v['zlen'] = self.read_values(file, f32, 3)
    v['alpha'], v['beta'], v['gamma'] = self.read_values(file, f32, 3)
    v['mapc'], v['mapr'], v['maps'] = self.read_values(file, i32, 3)
    v['amin'], v['amax'], v['amean'] = self.read_values(file, f32, 3)
    v['ispg'], v['nsymbt'] = self.read_values(file, i32, 2)
    if file_type == 'ccp4':
      v['lskflg'] = self.read_values(file, i32, 1)
      v['skwmat'] = self.read_values(file, f32, 9)
      v['skwtrn'] = self.read_values(file, f32, 3)
      v['user'] = self.read_values(file, i32, CCP4_USER)
      v['map'] = file.read(4)   # Should be 'MAP '.
      v['machst'] = self.read_values(file, i32, 1)
      v['rms'] = self.read_values(file, f32, 1)
      v['type'] = 'ccp4'
    else:
      # MRC file
      user = file.read(4*MRC_USER)
      if user[-4:] == 'MAP ':
        # New style MRC 2000 format file with xyz origin
        v['user'] = self.read_values_from_string(user, i32, MRC_USER)[:-4]
        xyz_origin = self.read_values_from_string(user[-16:-4], f32, 3)
        v['xorigin'], v['yorigin'], v['zorigin'] = xyz_origin
        v['imodStamp'] = self.read_values_from_string(user[56:60], i32, 1)
        v['imodFlags'] = self.read_values_from_string(user[60:64], i32, 1)
        v['machst'] = self.read_values(file, i32, 1)
        v['rms'] = self.read_values(file, f32, 1)
        v['type'] = 'mrc2000'
      else:
        # Old style MRC has xy origin instead of machst and rms.
        v['user'] = self.read_values_from_string(user, i32, MRC_USER)
        v['xorigin'], v['yorigin'] = self.read_values(file, f32, 2)
        v['type'] = 'mrc'

    v['nlabl'] = self.read_values(file, i32, 1)
    labels = []
    for i in range(MRC_NUM_LABELS):
      labels.append(file.read(MRC_LABEL_SIZE))
    v['labels'] = labels

    # Catch incorrect nsymbt value.
    if v['nsymbt'] < 0 or v['nsymbt'] + MRC_HEADER_LENGTH > file_size:
      raise Exception('MRC header value nsymbt (%d) is invalid'
                          % v['nsymbt'])
    v['symop'] = file.read(v['nsymbt'])

    return v

  # ---------------------------------------------------------------------------
  #
  def value_type(self, mode, unsigned_8_bit):

    MODE_char   = 0
    MODE_short  = 1
    MODE_float  = 2
    MODE_ushort  = 6            # Non-standard
    
    from numpy import uint8, int8, int16, uint16, float32, dtype
    if mode == MODE_char:
      if unsigned_8_bit:
        t = dtype(uint8)
      else:
        t = dtype(int8)        # CCP4 or MRC2000
    elif mode == MODE_short:
      t = dtype(int16)
    elif mode == MODE_ushort:
      t = dtype(uint16)
    elif mode == MODE_float:
      t = dtype(float32)
    else:
      raise SyntaxError('MRC data value type (%d) ' % mode +
                          'is not 8 or 16 bit integers or 32 bit floats')

    return t

  # ---------------------------------------------------------------------------
  #
  def check_header_values(self, v, file_size, file):

    if v['nc'] <= 0 or v['nr'] <= 0 or v['ns'] <= 0:
      raise SyntaxError('Bad MRC grid size (%d,%d,%d)'
                          % (v['nc'],v['nr'],v['ns']))

    esize = self.element_type.itemsize
    data_size = int(v['nc']) * int(v['nr']) * int(v['ns']) * esize
    header_end = file.tell()
    if header_end + data_size > file_size:
      if v['nsymbt'] and (header_end - v['nsymbt']) + data_size == file_size:
        # Sometimes header indicates symmetry operators are present but
        # they are not.  This error occurs in macromolecular structure database
        # entries emd_1042.map, emd_1048.map, emd_1089.map, ....
        # This work around code allows the incorrect files to be read.
        file.seek(-v['nsymbt'], 1)
        v['symop'] = ''
      else:
        msg = ('File size %d too small for grid size (%d,%d,%d)'
               % (file_size, v['nc'],v['nr'],v['ns']))
        if v['nsymbt']:
          msg += ' and %d bytes of symmetry operators' % (v['nsymbt'],)
        raise SyntaxError(msg)

  # ---------------------------------------------------------------------------
  #
  def read_values(self, file, etype, count):

    from numpy import array
    esize = array((), etype).itemsize
    string = file.read(esize * count)
    
    if len(string) < esize * count:
        raise SyntaxError('MRC file is truncated.  Failed reading %d values, type %s' % (count, etype.__name__))
        #return []
    values = self.read_values_from_string(string, etype, count)
    return values

  # ---------------------------------------------------------------------------
  #
  def read_values_from_string(self, string, etype, count):
  
    from numpy import fromstring
    values = fromstring(string, etype)
    if self.swap_bytes:
      values = values.byteswap()
    if count == 1:
      return values[0]
    return values

  # ---------------------------------------------------------------------------
  # Reads a submatrix from a the file.
  # Returns 3d numpy matrix with zyx index order.
  #
  def read_matrix(self, ijk_origin, ijk_size, ijk_step, progress):

    # ijk correspond to xyz.  crs refers to fast,medium,slow matrix file axes.
    crs_origin = [ijk_origin[a] for a in self.crs_to_ijk]
    crs_size = [ijk_size[a] for a in self.crs_to_ijk]
    crs_step = [ijk_step[a] for a in self.crs_to_ijk]

    #from readarray import read_array
    matrix = read_array(self.path, self.data_offset,
                        crs_origin, crs_size, crs_step,
                        self.matrix_size, self.element_type, self.swap_bytes,
                        progress)
    if not matrix is None:
      matrix = self.permute_matrix_to_xyz_axis_order(matrix)
    
    return matrix

  # ---------------------------------------------------------------------------
  #
  def permute_matrix_to_xyz_axis_order(self, matrix):
    
    if self.ijk_to_crs == (0,1,2):
      return matrix

    kji_to_src = [2-self.ijk_to_crs[2-a] for a in (0,1,2)]
    m = matrix.transpose(kji_to_src)

    return m

  # ---------------------------------------------------------------------------
  #
  def symmetry_matrices(self):

    h = self.header
    for name in ('symop', 'xlen', 'ylen', 'zlen', 'alpha', 'beta', 'gamma'):
      if not name in h:
        return []

    # Read space group symmetries in fractional coordinates
    s = h['symop']
    nsym = len(s)/80
    from Crystal.space_groups import parse_symop
    try:
      usyms = [parse_symop(s[80*i:80*(i+1)].replace(' ', ''))
               for i in range(nsym)]
    except:
      try:
        msg = 'Unable to parse symmetry operators of %s\n%s\n' % (self.name, s)
      except:
        # Garbage in sym operator data can't be interpreted as text.
        msg = 'Unable to parse symmetry operators of %s\n' % (self.name,)
      print(msg)
      return []

    a, b, c = float(h['xlen']), float(h['ylen']), float(h['zlen'])
    if a == 0 or b == 0 or c == 0:
      return []
    from math import pi
    alpha, beta, gamma = [d*pi/180 for d in self.cell_angles]

    # Convert symmetries to unit cell coordinates
    import Crystal
    u2r = Crystal.unit_cell_to_xyz_matrix(a, b, c, alpha, beta, gamma)

    from Matrix import invert_matrix, multiply_matrices
    r2u = invert_matrix(u2r)
    syms = [multiply_matrices(u2r, u2u, r2u) for u2u in usyms]
  
    return syms

# -----------------------------------------------------------------------------
#
def valid_cell_angles(alpha, beta, gamma, path):

  err = None
  
  for a in (alpha, beta, gamma):
    if a <= 0 or a >= 180:
      err = 'must be between 0 and 180'

  if alpha + beta + gamma >= 360 and err is None:
    err = 'sum must be less than 360'

  if max((alpha, beta, gamma)) >= 0.5 * (alpha + beta + gamma) and err is None:
    err = 'largest angle must be less than sum of other two'

  if err:
    from sys import stderr
    stderr.write('%s: invalid cell angles %.5g,%.5g,%.5g %s.\n'
                 % (path, alpha, beta, gamma, err))
    return False

  return True

def invert_matrix(tf):

    from numpy import array, zeros, float, linalg
    tf = array(tf)
    r = tf[:,:3]
    t = tf[:,3]
    tfinv = zeros((3,4), float)
    rinv = tfinv[:,:3]
    tinv = tfinv[:,3]
    from numpy.linalg import inv as matrix_inverse
    from numpy import dot as matrix_multiply
    rinv[:,:] = matrix_inverse(r)
    tinv[:] = matrix_multiply(rinv, -t)
    return tfinv

def skew_axes(cell_angles):

  if tuple(cell_angles) == (90,90,90):
    # Use exact trig functions for this common case.
    ca = cb = cg = c1 = 0
    sg = c2 = 1
  else:
    # Convert to radians
    from math import pi, sin, cos, sqrt
    alpha, beta, gamma = [a * pi / 180 for a in cell_angles]
    cg = cos(gamma)
    sg = sin(gamma)
    cb = cos(beta)
    ca = cos(alpha)
    c1 = (ca - cb*cg)/sg
    c2 = sqrt(1 - cb*cb - c1*c1)

  axes = ((1, 0, 0), (cg, sg, 0), (cb, c1, c2))
  return axes


# -----------------------------------------------------------------------------
# Maintain a cache of data objects using a limited amount of memory.
# The least recently accessed data is released first.
# -----------------------------------------------------------------------------

class Data_Cache:

    def __init__(self, size):

        self.size = size
        self.used = 0
        self.time = 1
        self.data = {}
        self.groups = {}

  # ---------------------------------------------------------------------------
  #
    def cache_data(self, key, value, size, description, groups = []):

        self.remove_key(key)
        d = Cached_Data(key, value, size, description,
                self.time_stamp(), groups)
        self.data[key] = d

        for g in groups:
            gtable = self.groups
            if not gtable.has_key(g):
                gtable[g] = []
            gtable[g].append(d)

        self.used = self.used + size
        self.reduce_use()

    # ---------------------------------------------------------------------------
    #
    def lookup_data(self, key):

        data = self.data
        if data.has_key(key):
            d = data[key]
            d.last_access = self.time_stamp()
            v = d.value
        else:
            v = None
        self.reduce_use()
        return v

    # ---------------------------------------------------------------------------
    #
    def remove_key(self, key):

        data = self.data
        if data.has_key(key):
            self.remove_data(data[key])
        self.reduce_use()

  # ---------------------------------------------------------------------------
  #
    def group_keys_and_data(self, group):

        groups = self.groups
        if not groups.has_key(group):
            return []

        kd = map(lambda d: (d.key, d.value), groups[group])
        return kd

      # ---------------------------------------------------------------------------
  #
    def resize(self, size):

        self.size = size
        self.reduce_use()

  # ---------------------------------------------------------------------------
  #
    def reduce_use(self):

        if self.used <= self.size:
            return

        data = self.data
        dlist = data.values()
        dlist.sort(lambda d1, d2: cmp(d1.last_access, d2.last_access))
        import sys
        for d in dlist:
            if sys.getrefcount(d.value) == 2:
                self.remove_data(d)
            if self.used <= self.size:
                break

  # ---------------------------------------------------------------------------
  #
    def remove_data(self, d):

        del self.data[d.key]
        self.used = self.used - d.size
        d.value = None

        for g in d.groups:
            dlist = self.groups[g]
            dlist.remove(d)
            if len(dlist) == 0:
                del self.groups[g]

  # ---------------------------------------------------------------------------
  #
    def time_stamp(self):
        t = self.time
        self.time = t + 1
        return t

# -----------------------------------------------------------------------------
#
class Cached_Data:

  def __init__(self, key, value, size, description, time_stamp, groups):

    self.key = key
    self.value = value
    self.size = size
    self.description = description
    self.last_access = time_stamp
    self.groups = groups


data_cache = Data_Cache(size = 0)


# -------------------------------------------- #
#         Grid_Data: datatype of chimera       #
# -------------------------------------------- #


class Grid_Data:

  def __init__(self, size,
               value_type = float32,
               origin = (0,0,0),
               step = (1,1,1),
               cell_angles = (90,90,90),
               rotation = ((1,0,0),(0,1,0),(0,0,1)),
               symmetries = (),
               name = '',
               path = '',       # Can be list of paths
               file_type = '',
               grid_id = '',
               default_color = (.7,.7,.7,1)):

    # Path, file_type and grid_id are for reloading data sets.
    self.path = path
    self.file_type = file_type  # 'mrc', 'spider', ....
    self.grid_id = grid_id      # String identifying grid in multi-grid files.
    
    if name == '':
      name = self.name_from_path(path)
    self.name = name
    
    self.size = tuple(size)

    from numpy import dtype
    if not isinstance(value_type, dtype):
      value_type = dtype(value_type)
    self.value_type = value_type        # numpy dtype.

    # Parameters defining how data matrix is positioned in space
    self.origin = tuple(origin)
    self.original_origin = self.origin
    self.step = tuple(step)
    self.original_step = self.step
    self.cell_angles = tuple(cell_angles)
    self.rotation = tuple(map(tuple, rotation))
    self.symmetries = symmetries
    self.ijk_to_xyz_transform = None
    self.xyz_to_ijk_transform = None

    self.rgba = default_color            # preferred color for displaying data

    global data_cache
    self.data_cache = data_cache

    self.writable = False
    self.change_callbacks = []

    self.update_transform()

  # ---------------------------------------------------------------------------
  #
  def set_path(self, path, format = None):

    if path != self.path:
      self.path = path
      self.name = self.name_from_path(path)
      self.call_callbacks('path changed')
      
    if format and format != self.file_type:
      self.file_type = format
      self.call_callbacks('file format changed')

  # ---------------------------------------------------------------------------
  #
  def name_from_path(self, path):

    from os.path import basename
    if isinstance(path, (list,tuple)):  p = path[0]
    else:                               p = path
    name = basename(p)
    return name

  # ---------------------------------------------------------------------------
  #
  def set_origin(self, origin):

    if tuple(origin) != self.origin:
      self.origin = tuple(origin)
      self.update_transform()

    # TODO: Update symmetries for origin, step, cell angle and rotation changes

  # ---------------------------------------------------------------------------
  #
  def set_step(self, step):

    if tuple(step) != self.step:
      self.step = tuple(step)
      self.update_transform()

  # ---------------------------------------------------------------------------
  #
  def set_cell_angles(self, cell_angles):

    if tuple(cell_angles) != self.cell_angles:
      self.cell_angles = tuple(cell_angles)
      self.update_transform()

  # ---------------------------------------------------------------------------
  #
  def set_rotation(self, rotation):

    r = tuple(map(tuple, rotation))
    if r != self.rotation:
      self.rotation = r
      self.update_transform()

  # ---------------------------------------------------------------------------
  # Compute 3 by 4 matrices encoding rotation and translation.
  #
  def update_transform(self):

    
    saxes = skew_axes(self.cell_angles)
    rsaxes = [apply_rotation(self.rotation, a) for a in saxes]
    tf, tf_inv = transformation_and_inverse(self.origin, self.step, rsaxes)
    if tf != self.ijk_to_xyz_transform or tf_inv != self.xyz_to_ijk_transform:
      self.ijk_to_xyz_transform = tf
      self.xyz_to_ijk_transform = tf_inv
      self.coordinates_changed()

  # ---------------------------------------------------------------------------
  # A matrix ijk corresponds to a point in xyz space.
  # This function maps the xyz point to the matrix index.
  # The returned matrix index need not be integers.
  #
  def xyz_to_ijk(self, xyz):

    return map_point(xyz, self.xyz_to_ijk_transform)

  # ---------------------------------------------------------------------------
  # A matrix ijk corresponds to a point in xyz space.
  # This function maps the matrix index to the xyz point.
  #
  def ijk_to_xyz(self, ijk):

    return map_point(ijk, self.ijk_to_xyz_transform)
    
  # ---------------------------------------------------------------------------
  # Spacings in xyz space of jk, ik, and ij planes.
  #
  def plane_spacings(self):

    spacings = map(lambda u: 1.0/norm(u[:3]), self.xyz_to_ijk_transform)
    return spacings
    
  # ---------------------------------------------------------------------------
  #
  def matrix(self, ijk_origin = (0,0,0), ijk_size = None,
             ijk_step = (1,1,1), progress = None, from_cache_only = False):

    if ijk_size == None:
      ijk_size = self.size

    m = self.cached_data(ijk_origin, ijk_size, ijk_step)
    if m is None and not from_cache_only:
      m = self.read_matrix(ijk_origin, ijk_size, ijk_step, progress)
      self.cache_data(m, ijk_origin, ijk_size, ijk_step)

    return m
    
  # ---------------------------------------------------------------------------
  # Must overide this function in derived class to return a 3 dimensional
  # NumPy matrix.  The returned matrix has size ijk_size and
  # element ijk is accessed as m[k,j,i].  It is an error if the requested
  # submatrix does not lie completely within the full data matrix.  It is
  # also an error for the size to be <= 0 in any dimension.  These invalid
  # inputs might throw an exception or might return garbage.  It is the
  # callers responsibility to make sure the arguments are valid.
  #
  def read_matrix(self, ijk_origin = (0,0,0), ijk_size = None,
                  ijk_step = (1,1,1), progress = None):

    raise NotImplementedError('Grid %s has no read_matrix() routine' % self.name)
  
  # ---------------------------------------------------------------------------
  # Convenience routine.
  #
  def matrix_slice(self, matrix, ijk_origin, ijk_size, ijk_step):

    i1, j1, k1 = ijk_origin
    i2, j2, k2 = map(lambda i, s: i+s, ijk_origin, ijk_size)
    istep, jstep, kstep = ijk_step
    m = matrix[k1:k2:kstep, j1:j2:jstep, i1:i2:istep]
    return m
    
  # ---------------------------------------------------------------------------
  # Deprecated.  Used before matrix() routine existed.
  #
  def full_matrix(self, progress = None):

    matrix = self.matrix()
    return matrix
    
  # ---------------------------------------------------------------------------
  # Deprecated.  Used before matrix() routine existed.
  #
  def submatrix(self, ijk_origin, ijk_size):

    return self.matrix(ijk_origin, ijk_size)

  # ---------------------------------------------------------------------------
  #
  def cached_data(self, origin, size, step):

    dcache = self.data_cache
    if dcache is None:
      return None

    key = (self, tuple(origin), tuple(size), tuple(step))
    m = dcache.lookup_data(key)
    if not m is None:
      return m

    # Look for a matrix containing the desired matrix
    group = self
    kd = dcache.group_keys_and_data(group)
    for k, matrix in kd:
      orig, sz, st = k[1:]
      if (step[0] < st[0] or step[1] < st[1] or step[2] < st[2] or
          step[0] % st[0] or step[1] % st[1] or step[2] % st[2]):
        continue        # Step sizes not compatible
      if (origin[0] < orig[0] or origin[1] < orig[1] or origin[2] < orig[2] or
          origin[0] + size[0] > orig[0] + sz[0] or
          origin[1] + size[1] > orig[1] + sz[1] or
          origin[2] + size[2] > orig[2] + sz[2]):
        continue        # Doesn't cover.
      dstep = map(lambda a,b: a/b, step, st)
      offset = map(lambda a,b: a-b, origin, orig)
      if offset[0] % st[0] or offset[1] % st[1] or offset[2] % st[2]:
        continue        # Offset stagger.
      moffset = map(lambda o,s: o / s, offset, st)
      msize = map(lambda s,t: (s+t-1) / t, size, st)
      m = matrix[moffset[2]:moffset[2]+msize[2]:dstep[2],
                 moffset[1]:moffset[1]+msize[1]:dstep[1],
                 moffset[0]:moffset[0]+msize[0]:dstep[0]]
      dcache.lookup_data(key)			# update access time
      return m

    return None

  # ---------------------------------------------------------------------------
  #
  def cache_data(self, m, origin, size, step):

    dcache = self.data_cache
    if dcache is None:
      return

    key = (self, tuple(origin), tuple(size), tuple(step))
    elements = reduce(lambda a,b: a*b, m.shape, 1)
    bytes = elements * m.itemsize
    groups = [self]
    descrip = self.data_description(origin, size, step)
    dcache.cache_data(key, m, bytes, descrip, groups)

  # ---------------------------------------------------------------------------
  #
  def data_description(self, origin, size, step):

    description = self.name

    if origin == (0,0,0):
      bounds = ' (%d,%d,%d)' % tuple(size)
    else:
      region = (origin[0], origin[0]+size[0]-1,
		origin[1], origin[1]+size[1]-1,
		origin[2], origin[2]+size[2]-1)
      bounds = ' (%d-%d,%d-%d,%d-%d)' % region
    description += bounds

    if step != (1,1,1):
      description += ' step (%d,%d,%d)' % tuple(step)

    return description

  # ---------------------------------------------------------------------------
  #
  def clear_cache(self):

    dcache = self.data_cache
    if dcache is None:
      return

    for k,d in dcache.group_keys_and_data(self):
      dcache.remove_key(k)

  # ---------------------------------------------------------------------------
  #
  def add_change_callback(self, cb):

    self.change_callbacks.append(cb)

  # ---------------------------------------------------------------------------
  #
  def remove_change_callback(self, cb):

    self.change_callbacks.remove(cb)

  # ---------------------------------------------------------------------------
  # Code has modified matrix elements, or the value type has changed.
  #
  def values_changed(self):

    self.call_callbacks('values changed')

  # ---------------------------------------------------------------------------
  # Mapping of array indices to xyz coordinates has changed.
  #
  def coordinates_changed(self):

    self.call_callbacks('coordinates changed')

  # ---------------------------------------------------------------------------
  #
  def call_callbacks(self, reason):
    
    for cb in self.change_callbacks:
      cb(reason)
    
# -----------------------------------------------------------------------------
# Return 3 by 4 matrix where first 3 columns give rotation and last column
# is translation.
#
def transformation_and_inverse(origin, step, axes):
  
  ox, oy, oz = origin
  d0, d1, d2 = step
  ax, ay, az = axes

  tf = ((d0*ax[0], d1*ay[0], d2*az[0], ox),
        (d0*ax[1], d1*ay[1], d2*az[1], oy),
        (d0*ax[2], d1*ay[2], d2*az[2], oz))

  tf_inv = invert_matrix(tf)
  
  # Replace array by tuples
  tf_inv = tuple(map(tuple, tf_inv))
  
  return tf, tf_inv

# -----------------------------------------------------------------------------
# Apply scaling and skewing transformations.
#
def scale_and_skew(ijk, step, cell_angles):

  
  xa, ya, za = skew_axes(cell_angles)

  i, j, k = ijk
  d0, d1, d2 = step
  x,y,z = i*d0, j*d1, k*d2

  xyz = tuple(x*xa[a]+y*ya[a]+z*za[a] for a in (0,1,2))
  return xyz
    
# -----------------------------------------------------------------------------
#
def apply_rotation(r, v):
  
  rv = [r[a][0]*v[0] + r[a][1]*v[1] + r[a][2]*v[2] for a in (0,1,2)]
  return tuple(rv)

# -----------------------------------------------------------------------------
#
def map_point(p, tf):

  tfp = [0,0,0]
  for r in range(3):
    tfr = tf[r]
    tfp[r] = tfr[0]*p[0] +tfr[1]*p[1] + tfr[2]*p[2] + tfr[3]
  tfp = tuple(tfp)
  return tfp

# -----------------------------------------------------------------------------
#
class Grid_Subregion(Grid_Data):

  def __init__(self, grid_data, ijk_min, ijk_max, ijk_step = (1,1,1)):

    self.full_data = grid_data

    ijk_min = map(lambda a,s: ((a+s-1)/s)*s, ijk_min, ijk_step)
    self.ijk_offset = ijk_min
    self.ijk_step = ijk_step

    size = map(lambda a,b,s: max(0,(b-a+s)/s), ijk_min, ijk_max, ijk_step)
    origin = grid_data.ijk_to_xyz(ijk_min)
    step = [ijk_step[a]*grid_data.step[a] for a in range(3)]

    Grid_Data.__init__(self, size, grid_data.value_type,
                       origin, step, grid_data.cell_angles,
                       grid_data.rotation, grid_data.symmetries,
                       name = grid_data.name + ' subregion')
    self.rgba = grid_data.rgba
    self.data_cache = None      # Caching done by underlying grid.
        
  # ---------------------------------------------------------------------------
  #
  def read_matrix(self, ijk_origin, ijk_size, ijk_step, progress):

    origin, step, size = self.full_region(ijk_origin, ijk_size, ijk_step)
    m = self.full_data.matrix(origin, size, step, progress)
    return m
        
  # ---------------------------------------------------------------------------
  #
  def cached_data(self, ijk_origin, ijk_size, ijk_step):

    origin, step, size = self.full_region(ijk_origin, ijk_size, ijk_step)
    m = self.full_data.cached_data(origin, size, step)
    return m
        
  # ---------------------------------------------------------------------------
  #
  def full_region(self, ijk_origin, ijk_size, ijk_step):

    origin = map(lambda i,s,o: i*s+o,ijk_origin, self.ijk_step, self.ijk_offset)
    size = map(lambda a,b: a*b, ijk_size, self.ijk_step)
    step = map(lambda a,b: a*b, ijk_step, self.ijk_step)
    return origin, step, size

  # ---------------------------------------------------------------------------
  #
  def clear_cache(self):

    self.full_data.clear_cache()

# -----------------------------------------------------------------------------
#
def norm(v):

  from math import sqrt
  d = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
  return d



# -------------------------------------------- #
#    Array_Grid_Data: wrapper of numpy.array   #
# -------------------------------------------- #

class Array_Grid_Data(Grid_Data):
  
  def __init__(self, array, origin = (0,0,0), step = (1,1,1),
               cell_angles = (90,90,90),
               rotation = ((1,0,0),(0,1,0),(0,0,1)),
               symmetries = (),
               name = ''):

      self.array = array
      
      path = ''
      file_type = ''
      component_name = ''

      grid_size = list(array.shape)
      grid_size.reverse()

      value_type = array.dtype

      Grid_Data.__init__(self, grid_size, value_type,
                         origin, step, cell_angles = cell_angles,
                         rotation = rotation, symmetries = symmetries,
                         name = name, path = path, file_type = file_type)

      self.writable = True
  
  # ---------------------------------------------------------------------------
  #
  def read_matrix(self, ijk_origin, ijk_size, ijk_step, progress):
  
    return self.cached_data(ijk_origin, ijk_size, ijk_step)

  # ---------------------------------------------------------------------------
  #
  def cached_data(self, ijk_origin, ijk_size, ijk_step):

      m = self.matrix_slice(self.array, ijk_origin, ijk_size, ijk_step)
      return m





# -----------------------------------------------------------------------------
# Write an MRC 2000 format file.
#
# Header contains four byte integer or float values:
#
# 1	NX	number of columns (fastest changing in map)	
# 2	NY	number of rows					
# 3	NZ	number of sections (slowest changing in map)	
# 4	MODE	data type :					
# 		0	image : signed 8-bit bytes range -128 to 127
# 		1	image : 16-bit halfwords		
# 		2	image : 32-bit reals			
# 		3	transform : complex 16-bit integers	
# 		4	transform : complex 32-bit reals	
# 5	NXSTART	number of first column in map			
# 6	NYSTART	number of first row in map			
# 7	NZSTART	number of first section in map			
# 8	MX	number of intervals along X			
# 9	MY	number of intervals along Y			
# 10	MZ	number of intervals along Z			
# 11-13	CELLA	cell dimensions in angstroms			
# 14-16	CELLB	cell angles in degrees				
# 17	MAP# axis corresp to cols (1,2,3 for X,Y,Z)		
# 18	MAPR	axis corresp to rows (1,2,3 for X,Y,Z)		
# 19	MAPS	axis corresp to sections (1,2,3 for X,Y,Z)	
# 20	DMIN	minimum density value				
# 21	DMAX	maximum density value				
# 22	DMEAN	mean density value				
# 23	ISPG	space group number 0 or 1 (default=0)		
# 24	NSYMBT	number of bytes used for symmetry data (0 or 80)
# 25-49   EXTRA	extra space used for anything			
# 50-52	ORIGIN  origin in X,Y,Z used for transforms		
# 53	MAP	character string 'MAP ' to identify file type	
# 54	MACHST	machine stamp					
# 55	RMS	rms deviation of map from mean density		
# 56	NLABL	number of labels being used			
# 57-256 LABEL(20,10) 10 80-character text labels		
#

# -----------------------------------------------------------------------------
#
def write_mrc2000_grid_data(grid_data, path, options = {}, progress = None,angle_tomo=0):

    mtype = grid_data.value_type.type
    type = closest_mrc2000_type(mtype)

    f = open(path, 'wb')
    if progress:
        progress.close_on_cancel(f)

    header = mrc2000_header(grid_data, type)
    f.write(header)

    stats = Matrix_Statistics()
    isz, jsz, ksz = grid_data.size
    for k in range(ksz):
        matrix = grid_data.matrix((0,0,k), (isz,jsz,1))
        if type != mtype:
            matrix = matrix.astype(type)
        f.write(matrix.tostring())
        stats.plane(matrix)
        if progress:
            progress.plane(k)

    # Put matrix statistics in header
    header = mrc2000_header(grid_data, type, stats,angle_tomo=angle_tomo)
    f.seek(0)
    f.write(header)

    f.close()

# -----------------------------------------------------------------------------
#
def mrc2000_header(grid_data, value_type, stats = None,angle_tomo=0):

    size = grid_data.size
    
    from numpy import float32, int16, int8, int32
    if value_type == float32:         mode = 2
    elif value_type == int16:         mode = 1
    elif value_type == int8:          mode = 0

    cell_size = map(lambda a,b: a*b, grid_data.step, size)

    if stats:
        dmin, dmax = stats.min, stats.max
        dmean, rms = stats.mean_and_rms(size)
    else:
        dmin = dmax = dmean = rms = 0

    from numpy import little_endian
    if little_endian:
        machst = 0x00004144
    else:
        machst = 0x11110000

    #from chimera.version import release
    from time import asctime
    ver_stamp = ' %s' %  asctime()
    labels = [ver_stamp[:80]]

    if grid_data.rotation != ((1,0,0),(0,1,0),(0,0,1)):
       
        axis, angle = rotation_axis_angle(grid_data.rotation)
        r = 'Chimera rotation: %12.8f %12.8f %12.8f %12.8f' % (axis + (angle,))
        labels.append(r)

    nlabl = len(labels)
    # Make ten 80 character labels.
    labels.extend(['']*(10-len(labels)))
    labels = [l + (80-len(l))*'\0' for l in labels]
    labelstr = ''.join(labels)

    cell_size = list(cell_size)
    extra = [0]*25
    if angle_tomo: extra = int(angle_tomo*1000)
    strings = [
        binary_string(size, int32),  # nx, ny, nz
        binary_string(mode, int32),  # mode
        binary_string((0,0,0), int32), # nxstart, nystart, nzstart
        binary_string(size, int32),  # mx, my, mz
        binary_string(cell_size, float32), # cella
        binary_string(grid_data.cell_angles, float32), # cellb
        binary_string((1,2,3), int32), # mapc, mapr, maps
        binary_string((dmin, dmax, dmean), float32), # dmin, dmax, dmean
        binary_string(0, int32), # ispg
        binary_string(0, int32), # nsymbt
        binary_string(extra, int32), # extra
        binary_string(grid_data.origin, float32), # origin
        binary_string('MAP ',str), # map
        binary_string(machst, int32), # machst
        binary_string(rms, float32), # rms
        binary_string(nlabl, int32), # nlabl
        binary_string(labelstr,str),
        ]

    
    header = b"".join(strings)
    return header
    
# -----------------------------------------------------------------------------
#
class Matrix_Statistics:

    def __init__(self):

        self.min = None
        self.max = None
        self.sum = 0.0
        self.sum2 = 0.0

    def plane(self, matrix):

        from numpy import ravel, minimum, maximum, add, multiply, array, float32
        matrix_1d = matrix.ravel()
        dmin = minimum.reduce(matrix_1d)
        if self.min == None or dmin < self.min:
            self.min = dmin
        dmax = maximum.reduce(matrix_1d)
        if self.max == None or dmax > self.max:
            self.max = dmax
        self.sum += add.reduce(matrix_1d)
        # TODO: Don't copy array to get standard deviation.
        # Avoid overflow when squaring integral types
        m2 = array(matrix_1d, float32)
        multiply(m2, m2, m2)
        self.sum2 += add.reduce(m2)

    def mean_and_rms(self, size):

        vol = float(size[0])*float(size[1])*float(size[2])
        mean = self.sum / vol
        sdiff = self.sum2 - self.sum*self.sum
        if sdiff > 0:
            from math import sqrt
            rms = sqrt(sdiff) / vol
        else:
            rms = 0
        return mean, rms

# -----------------------------------------------------------------------------
#
def binary_string(values, type):
    
    from numpy import array
    return array(values, type).tostring()

# -----------------------------------------------------------------------------
#
def closest_mrc2000_type(type):

    from numpy import float, float32, float64
    from numpy import int, int8, int16, int32, int64, character
    from numpy import uint, uint8, uint16, uint32, uint64
    if type in (float32, float64, float, int32, int64, int, uint, uint16, uint32, uint64):
        ctype = float32
    elif type in (int16, uint8):
        ctype = int16
    elif type in (int8, character):
        ctype = int8
    else:
        raise TypeError('Volume data has unrecognized type %s' % type)

    return ctype



def convert_numpy_array3d_mrc(data, writename,angle_tomo=0):
    ''' USAGE convert_spimage3d_mrc( data_in_numpy_array_float64, [outputname].mrc )
    This function converts a numpy array (float64) to mrc. Writen to convert spimage
    objects to MRC
    '''
    data = data.real

    if not writename[-4:] == '.mrc': writename += '.mrc'
    write_mrc2000_grid_data( Array_Grid_Data(data), writename, angle_tomo=angle_tomo)

if __name__=="__main__":
    import sys
    import numpy
    
    filename = sys.argv[1]
    writename = sys.argv[2]
    angle = float(sys.argv[3])
    original_data = numpy.zeros((7,41,12),dtype=int)#spimage.sp_image_read(filename,0).image.real
    
    convert_numpy_array3d_mrc(original_data,writename,angle_tomo=angle)

