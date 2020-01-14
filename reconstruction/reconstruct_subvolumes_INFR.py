#!/usr/bin/env python

'''
Created on Jan 28, 2013

@author: yuxiangchen
'''

if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Reconstruct the subvolumes from projections using Fourier space method.',
                          authors='Yuxiang Chen',
                          options=[ScriptOption('-p', 'Particle list.', True, False),
                                   ScriptOption('-d', 'Unbinned projection directory. Note the projections should be aligned but not weighted!', True, False),
                                   ScriptOption('-s', 'Output particle size (int).', True, False),
                                   ScriptOption('-b', 'Binning factor of the particle list.', True, False),
                                   ScriptOption('-o', 'Cutting offset of the particle list.', True, False),
                                   ScriptOption('-i', 'Number of iterations to run.', True, True),
                                   ScriptOption(['-h', '--help'], 'Help.', False, True)])
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        pl_filename, proj_dir, vol_size, binning, offset, iter, b_help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()
    if b_help is True:
        print(helper)
        sys.exit()
    
    # parse the argument
    vol_size = int(vol_size)
    binning = int(binning)
    offset = [int(i) for i in offset.split(',')]
    if iter:
        iter = int(iter)
    else:
        iter = 10
    
    # force the user to specify even-sized volume
    assert vol_size % 2 == 0
    
    # load projection list
    from pytom.reconstruction.reconstructionStructures import ProjectionList
    projections = ProjectionList()
    projections.loadDirectory(proj_dir)
    projections.sort()
    dim_x = projections[0].getXSize()
    dim_y = projections[0].getYSize()
    dim_z = dim_x # make z dim the same as x!
    
    # load particle list
    from pytom.basic.structures import ParticleList
    pl = ParticleList()
    pl.fromXMLFile(pl_filename)
    
    # read all the projections in tompy
    import numpy as np
    from pytom.tompy.io import read, write
    proj = []
    tilt_angles = []
    for p in projections:
        proj.append(read(p.getFilename()))
        tilt_angles.append(p.getTiltAngle())
    
    # reconstruct each particles
    from math import cos, sin, pi, ceil
    from pytom.tompy.transform import cut_from_projection
    from nufft.reconstruction import fourier_2d1d_iter_reconstruct, fourier_2d1d_gridding_reconstruct
    from pytom.tools.ProgressBar import FixedProgBar
    prog = FixedProgBar(0, len(pl)-1, '')
    
    i = 0
    for p in pl:
        prog.update(i)
        i += 1
        
        # transfer the coordinate system
        x, y, z = p.getPickPosition().toVector()
        x = (x+offset[0])*binning
        y = (y+offset[1])*binning
        z = (z+offset[2])*binning
        
        # cut out corresponding parts from projections
        subregions = []
        l = (vol_size/2)*2**0.5
        cl = int(np.round(l))
        for img, ang in zip(proj, tilt_angles):
            # project the coordinate to 2D image
            yy = y # assume the rotation axis is around y
            xx = (cos(ang*pi/180)*(x-dim_x/2)-sin(ang*pi/180)*(z-dim_z/2)) + dim_x/2
            
            # cut out corresponding parts from projections
            # get a bigger region to hold all the info
#            subregion = np.zeros((cl*2, cl*2, 1))
#            
#            # calculate the region size of certain angle
#            s = l*cos((45-abs(ang))*pi/180)
#            cs = int(np.round(s))
#            
#            # cut the small patch out
#            patch = cut_from_projection(img, [xx,yy], [cs*2, vol_size])
#            
#            # fill in the subregion
#            subregion[cl-cs:cl+cs, cl-vol_size/2:cl+vol_size-vol_size/2] = patch
#            subregions.append(subregion)
#        
#        # reconstruct
#        v = fourier_2d1d_iter_reconstruct(subregions, tilt_angles, iter)
#        
#        # get the center part
#        v = v[cl-vol_size/2:cl+vol_size-vol_size/2, cl-vol_size/2:cl+vol_size-vol_size/2, cl-vol_size/2:cl+vol_size-vol_size/2]
        
        
        subregions = []
        for img, ang in zip(proj, tilt_angles):
            # project the coordinate to 2D image
            yy = y # assume the rotation axis is around y
            xx = (cos(ang*pi/180)*(x-dim_x/2)-sin(ang*pi/180)*(z-dim_z/2)) + dim_x/2
            
            # cut the small patch out
            patch = cut_from_projection(img, [xx,yy], [vol_size, vol_size])
            patch = patch - np.mean(patch)
            
            # fill in the subregion
            subregions.append(patch)
        
        # reconstruct
        v = fourier_2d1d_iter_reconstruct(subregions, tilt_angles, iter)
        
        
        # write to the disk
        write(p.getFilename(), v)
    
    