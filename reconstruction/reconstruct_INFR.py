#!/usr/bin/env python

'''
Created on Jun 6, 2013

@author: yuxiangchen
'''
import numpy
import os
from pytom.gui.guiFunctions import datatype, loadstar


if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Reconstruct the whole tomogram from projections using iterative NFFT method.',
                          authors='Yuxiang Chen',
                          options=[ScriptOption('-d', 'Unbinned projection directory. Note the projections should be aligned but not weighted!', True, False),
                                   ScriptOption('-o', 'Output filename.', True, False),
                                   ScriptOption('-i', 'Number of iterations to run.', True, True),
                                   ScriptOption('-m', 'Meta file including tiltangles.', True, False),
                                   ScriptOption(['-h', '--help'], 'Help.', False, True)])
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        proj_dir, output_filename, iter, metafile, b_help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()
    if b_help is True:
        print(helper)
        sys.exit()
    
    # parse the arguments
    if not iter:
        iter = 10
    else:
        iter = int(iter)

    if metafile and os.path.exists(metafile):
        metadata = loadstar(metafile, dtype=datatype)
        tiltAngles = metadata['TiltAngle']
    else:
        tiltAngles = []
    # start reconstruction
    from tompy.io import read, write
    from nufft.reconstruction import fourier_2d1d_iter_reconstruct
    from pytom.reconstruction.reconstructionStructures import ProjectionList
    projections = ProjectionList()
    projections.loadDirectory(proj_dir, tiltAngles=tiltAngles)
    projections.sort()
    
    projs = []
    tilt_angles = []
    for p in projections:
        print(p.getTiltAngle(), p.getFilename())
        projs.append(read(p.getFilename()))
        tilt_angles.append(p.getTiltAngle())
    
    v = fourier_2d1d_iter_reconstruct(projs, tilt_angles, iter)
    write(output_filename, v)

