#!/usr/bin/env python

'''
Created on Sep 29, 2010

@author: chen
'''

def vol2sf(vol, r, b, center=None):
    """Transfer a volume into a serial of spherical functions. The volume will be decomposed into a serials of concentric shells.
    It will sample the volume on equiangular 2B*2B grid. (Trilinear interpolation) For more detail, see "http://www.cs.dartmouth.edu/~geelong/sphere/".
    
    Parameters
    ----------
    vol: Target volume
         pytom.lib.pytom_volume.vol

    r: The radius of shell you want to get from the volume (in voxel)
       Integer

    b: Bandwidth, which determines the sampling on the sphere
       Integer

    center: The center of the circle (By default is the center of the volume)
            list [x, y, z]
    
    Returns
    -------
    f(the_0, phi_0) ... f(the_0, phi_2B-1), f(the_2B-1, phi_0) ... f(the_2B-1, phi_2B-1)
    List
    """
    if r>=vol.sizeX()/2 or r>=vol.sizeY()/2 or r>=vol.sizeZ()/2:
        raise RuntimeError("Given radius is larger than the volume!")
    elif r <= 0:
        raise RuntimeError("Radius should be larger than the 0!")
    
    from pytom.lib.pytom_volume import volTOsf
    from pytom.lib.pytom_numpy import vol2npy
    if center:
        m_x, m_y, m_z = center
    else:
        # the middle point of the volume
        m_x = int(vol.sizeX()/2); m_y = int(vol.sizeY()/2); m_z = int(vol.sizeZ()/2)
    res = volTOsf(vol, r, b, m_x, m_y, m_z)
    
    # transfer the data format from volume to numpy array
    tmp = vol2npy(res) # somehow have to split the steps like this, otherwise wont work
    
    import numpy as np
    res = np.array(tmp) # copy the memory, otherwise the memory is shared with volume

#    from math import pi, sin, cos
#    res = [] # result list
#    
#    from pytom_volume import interpolateSpline
#    
#    if center:
#        m_x, m_y, m_z = center
#    else:
#        # the middle point of the volume
#        m_x = int(vol.sizeX()/2); m_y = int(vol.sizeY()/2); m_z = int(vol.sizeZ()/2)
#    
#    for j in xrange(2*b):
#        for k in xrange(2*b):
#            the = pi*(2*j+1)/(4*b) # (0,pi)
#            phi = pi*k/b # [0,2*pi)
#            
#            # trilinear interpolation
#            x = r*cos(phi)*sin(the); # _x = int(ceil(x)); x_ = int(floor(x)); dx = x-x_
#            y = r*sin(phi)*sin(the); # _y = int(ceil(y)); y_ = int(floor(y)); dy = y-y_
#            z = r*cos(the); # _z = int(ceil(z)); z_ = int(floor(z)); dz = z-z_
##            c_000 = vol.getV(m_x+x_,m_y+y_,m_z+z_)
##            c_100 = vol.getV(m_x+_x,m_y+y_,m_z+z_)
##            c_001 = vol.getV(m_x+x_,m_y+y_,m_z+_z)
##            c_101 = vol.getV(m_x+_x,m_y+y_,m_z+_z)
##            c_010 = vol.getV(m_x+x_,m_y+_y,m_z+z_)
##            c_110 = vol.getV(m_x+_x,m_y+_y,m_z+z_)
##            c_011 = vol.getV(m_x+x_,m_y+_y,m_z+_z)
##            c_111 = vol.getV(m_x+_x,m_y+_y,m_z+_z)
##            
##            c_00 = c_000*(1-dx)+c_100*dx
##            c_01 = c_001*(1-dx)+c_101*dx
##            c_10 = c_010*(1-dx)+c_110*dx
##            c_11 = c_011*(1-dx)+c_111*dx
##            c_0 = c_00*(1-dy)+c_10*dy
##            c_1 = c_01*(1-dy)+c_11*dy
##            c = c_0*(1-dz)+c_1*dz
#            
#            c = interpolateSpline(vol, m_x+x, m_y+y, m_z+z)
#            res.append(c)
    
    return res


def fvol2sf(vol, r, b):
    """Transfer a volume in Fourier space into a serial of spherical functions.
    """
    if r>=vol.sizeX()/2 or r>=vol.sizeY()/2 or r>=vol.sizeZ()/2:
        raise RuntimeError("Given radius is larger than the volume!")
    elif r <= 0:
        raise RuntimeError("Radius should be larger than the 0!")
    
    from pytom.lib.pytom_volume import fvolTOsf
    res = fvolTOsf(vol, r, b)
    
    return res


def vol2sf_mean(vol, b, max_radius=None, center=None):
    """Obsolete.
    """
    from math import pi, sin, cos
    from pytom.lib.pytom_volume import interpolate
    res = []
    
    if max_radius is None:
        max_radius = min([vol.sizeX()/2, vol.sizeY()/2, vol.sizeZ()/2])
    
    if center is not None:
        m_x, m_y, m_z = center
    else:
        m_x = int(vol.sizeX()/2); m_y = int(vol.sizeY()/2); m_z = int(vol.sizeZ()/2)
    
    for j in range(2*b):
        for k in range(2*b):
            the = pi*(2*j+1)/(4*b) # (0,pi)
            phi = pi*k/b # [0,2*pi)
            sum = 0.
            for r in range(max_radius):
                x = r*cos(phi)*sin(the)
                y = r*sin(phi)*sin(the)
                z = r*cos(the)
                sum += interpolate(vol, m_x+x, m_y+y, m_z+z)
            res.append(sum/max_radius)
            
    return res

def norm_sf(sf):
    """Normalize the spherical function to mean=0 and std=1.
    """
    from numpy import mean, std
    m = mean(sf); s = std(sf)
    res = [(i-m)/s for i in sf]
    return res

def sf2file(sf, filename, conv=True):
    """Save the spherical function as a file.
    The file format is is same as "http://www.cs.dartmouth.edu/~geelong/sphere/".

    Parameters
    ----------
    sf: Spherical function
        List

    filename: File name to save as
              String
    """
    f = open(filename, 'w')
    
    for v in sf:
        f.write("%f\n" % v)
        if conv is False:
            f.write("0\n") # imaginary part is always 0
    
    f.close()

if __name__ == '__main__':
    import sys, getopt
    usage = './scriptname -v volume -r radius -b bandwidth -c flag_for_convolution -o output_filename'
    
    if len(sys.argv) == 1:
        print(usage)
        sys.exit()        
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hcv:r:b:o:", ["help"])
    except getopt.GetoptError:
        print('Command not right. Exit!')
        sys.exit()
    
    conv = False
    for o,a in opts:
        if o in ("-h", "--help"):
            print(usage)
            sys.exit()
        if o in ("-v"):
            vol_name = a
        if o in ("-r"):
            radius = int(a)
        if o in ("-b"):
            bw = int(a)
        if o in ("-c"):
            conv = True
        if o in ("-o"):
            filename = a
    
    from pytom.lib.pytom_volume import read
    vol = read(vol_name)
    sf = vol2sf(vol, radius, bw)
    sf2file(sf, filename, conv)
