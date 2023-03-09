'''
Created on Aug 6, 2012

@author: yuxiangchen
'''

import numpy as np

def create_bfactor_fnc(size, ps, bfactor, FSC=None, apply_range=None):
    if FSC is None:
        FSC = np.ones(size/2)
    
    # transfer the pixel size to angstrom
    x = (ps*size)/np.arange(1, size/2+1, 1) # starts from 1!
    if apply_range is None:
        apply_range_pixel = [1, size/2] # starts from 1!
    else:
        assert apply_range[0] > apply_range[1]
        apply_range_pixel = [size*ps/apply_range[0], size*ps/apply_range[1]]
    
    # create the FSC weight
    FSC_weight = np.sqrt((2*np.array(FSC))/(1+np.array(FSC)))
    
    # calculate the decay function
    decay = FSC_weight*np.exp(-bfactor/(np.power(x, 2)*4))
    
    # apply only in certain range
    decay[:apply_range_pixel[0]-1] = 0
    decay[apply_range_pixel[1]-1:] = 0
    
    return decay

def create_bfactor_restore_fnc(size, ps, bfactor, FSC=None, apply_range=None):
    if FSC is None:
        FSC = np.ones(size/2)
    
    # transfer the pixel size to angstrom
    x = (ps*size)/np.arange(1, size/2+1, 1) # starts from 1!
    if apply_range is None:
        apply_range_pixel = [1, size/2] # starts from 1!
    else:
        assert apply_range[0] > apply_range[1]
        apply_range_pixel = [size*ps/apply_range[0], size*ps/apply_range[1]]
    
    # create the FSC weight
    FSC_weight = np.sqrt((2*np.array(FSC))/(1+np.array(FSC)))
    
    # calculate the decay function
    decay = FSC_weight*np.exp(-bfactor/(np.power(x, 2)*4))
    
    restore = 1.0/decay
    
    # apply only in certain range
    restore[:apply_range_pixel[0]-1] = 0
    restore[apply_range_pixel[1]-1:] = 0
    
    return restore

def create_bfactor_vol(size, ps, bfactor, FSC=None, apply_range=None):
    """Create a B-factor volume in Frequency space.
    @param size: The size of the volume, assuming it is a cube
    @param ps: The pixel size in angstrom
    @param bfactor: B factor
    @param FSC: Fourier Shell Correlation
    @param apply_range: The apply range (also in angstrom) of the B factor correction
    """
    if FSC is None:
        FSC = np.ones(size/2)
    
    # transfer the pixel size to angstrom
    x = (ps*size)/np.arange(1, size/2+1, 1) # starts from 1!
    if apply_range is None:
        apply_range_pixel = [1, size/2] # starts from 1!
    else:
        assert apply_range[0] > apply_range[1]
        apply_range_pixel = [size*ps/apply_range[0], size*ps/apply_range[1]]
    
    # create the FSC weight
    FSC_weight = np.sqrt((2*np.array(FSC))/(1+np.array(FSC)))
    
    # calculate the decay function
    decay = FSC_weight*np.exp(-bfactor/(np.power(x, 2)*4))
    
    # make the decay volume
    v = sph2cart(decay, size)
    
    # transfer to the volume format and multiple with the mask
    from pytom.lib.pytom_volume import vol, initSphere
    from pytom.lib.pytom_numpy import npy2vol
    vv = npy2vol(np.array(v, dtype='float32', order='F'), 3)
    
    if apply_range_pixel[0] == 1:
        mask = vol(size, size, size)
        initSphere(mask, apply_range_pixel[1]-1, 0,0, 
	        size/2, size/2, size/2) # minus 1 to make it consistent with python
    else:
        mask1 = vol(size, size, size)
        mask2 = vol(size, size, size)
        initSphere(mask1, apply_range_pixel[0]-1, 0,0, size/2, size/2, size/2)
        initSphere(mask2, apply_range_pixel[1]-1, 0,0, size/2, size/2, size/2)
        mask = mask2 - mask1

    return vv*mask
    

def sph2cart(data, size):
    from scipy.interpolate import interp1d
    i, j, k = np.ogrid[-size/2:size/2, -size/2:size/2, -size/2:size/2]
    R = np.sqrt(i**2 + j**2 + k**2)
    f = interp1d(np.linspace(0, size/2-1, size/2), data, kind='cubic', bounds_error=False, fill_value=0)
    out = f(R)
    
    return out

def create_bfactor_restore_vol(size, ps, bfactor, FSC=None, apply_range=None):
    from pytom.lib.pytom_volume import vol
    v = create_bfactor_vol(size, ps, bfactor, FSC, apply_range)
    unit = vol(v)
    unit.setAll(1)
    restore = unit/v
    
    return restore

def bfactor_restore(v, ps, bfactor, FSC=None, apply_range=None):
    from pytom.basic.fourier import convolute
    kernel = create_bfactor_restore_vol(v.sizeX(), ps, bfactor, FSC, apply_range) # assuming the v is a cube!
    out = convolute(v, kernel, True)
    return out
