#!/usr/bin/env pytom

"""
Created on Sep 02, 2019

@author: dschulte
"""
from pytom.gpu.initialize import device, xp
from pytom.tompy.tools import taper_edges


def cut_patch(projection, ang, pick_position, vol_size=200, binning=8, dimz=None, offset=[0,0,0], projection_name=None):
    from pytom.tompy.tools import taper_edges
    #from pytom.gpu.initialize import xp
    from pytom.voltools import transform
    from pytom.tompy.transform import rotate3d, rotate_axis
    from pytom.tompy.transform import cut_from_projection
    from pytom.tompy.io import read

    # Get the size of the original projection
    dim_x = projection.shape[0]
    dim_z = dim_x if dimz is None else dimz

    x, y, z = pick_position
    x = (x + offset[0]) * binning
    y = (y + offset[1]) * binning
    z = (z + offset[2]) * binning

    # Get coordinates of the paricle adjusted for the tilt angle
    yy = y  # assume the rotation axis is around y
    xx = (xp.cos(ang * xp.pi / 180) * (x - dim_x / 2) - xp.sin(ang * xp.pi / 180) * (z - dim_z / 2)) + dim_x / 2

    # Cut the small patch out

    patch = projection[max(0, int(xp.floor(xx))-vol_size//2):int(xp.floor(xx))+vol_size//2, int(yy)-vol_size//2:int(yy)+vol_size//2,:]
    shiftedPatch = xp.zeros_like(patch)
    transform(patch, output=shiftedPatch, translation=[xx%1,yy%1,0], device=device, interpolation='filt_bspline')
    return shiftedPatch.squeeze()


if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.tompy.mpi import MPI
    from time import gmtime, strftime
    from pytom.polishing.reconstruct_local_alignment import polish_particles
    from pytom.tompy.io import read_size
    from pytom.tompy.mpi import MPI
    from pytom.reconstruction.reconstructionStructures import ProjectionList
    from pytom.reconstruction.reconstructionFunctions import toProjectionStackFromAlignmentResultsFile, alignImageUsingAlignmentResultFile
    from pytom.basic.structures import ParticleList
    from pytom.tompy.reconstruction_functions import backProjectGPU
    from pytom.tompy.io import write, read
    from pytom_numpy import vol2npy
    import numpy
    import os
    from time import time, sleep
    from pytom.gui.guiFunctions import loadstar
    from pytom.basic.datatypes import DATATYPE_METAFILE

    helper = ScriptHelper(sys.argv[0].split('/')[-1],  # script name
                          description='Reconstruct a tomogram on a GPU based on a alignment Results alignment.',
                          authors='Gvds',
                          options=[ScriptOption(['-p', '--projectionDirectory'], 'Projection Directory', True, False ),
                                   ScriptOption(['-o', '--outputDirectory'], 'Results are written to this directory',
                                                True, False),
                                   ScriptOption(['-b', '--coordinateBinning'],
                                                'Binning factor of the particle list.', True, True),
                                   ScriptOption(['-g', '--gpuID'], 'Which gpu do you want to use?', True, True),
                                   ScriptOption(['-m', '--metaFile'], 'Metafile containing tiltangles.', True,True),
                                   ScriptOption(['-a', '--alignmentResultsFile'], 'Alignment ResultsFile', True, True),
                                   ScriptOption(['-s', '--sizeReconstruction'], 'Size Reconstruction in pixels. Three numbers separated by a comma.', True, True),
                                   ScriptOption(['-v', '--verbose'], 'print intermediate output', False, True)
                                   ])
    try:
        proj_dir, outdir, coordinateBinning, gpu, metafile, alignmentfile, size, verbose = ll = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        print(helper)
        sys.exit()

    coordinateBinning = int(coordinateBinning) if coordinateBinning else 1
    metadata = loadstar(metafile,dtype=DATATYPE_METAFILE)
    tilt_angles = metadata['TiltAngle']
    size = [464, 464, 464] if size is None else list(map(int,size.split(',')))
    patches = xp.zeros((size[0], size[1], len(tilt_angles)),dtype=xp.float32)

    images = []


    tt = time()

    cur = 0
    missed = 0
    for i in range(patches.shape[2]):
        temp_image = alignImageUsingAlignmentResultFile(alignmentfile, i, weighting=-1, circleFilter=True, binning=coordinateBinning)
        patches[:,:,i] = temp_image[:,:]
        del temp_image
    print(time() - tt)
    vol_bp = xp.zeros((size[0], size[1], size[2]), dtype=xp.float32)
    tt = time()
    s = 100
    bp = backProjectGPU(patches[:s,:s,:], vol_bp[:s,:s,:s], 0, tilt_angles)

    print(time()-tt)
    write(f'{outdir}/reconstruction.mrc', bp)#[ndim//2-vol_size//2:ndim//2+vol_size//2 + vol_size%2,
                                           #ndim//2-vol_size//2:ndim//2+vol_size//2 + vol_size%2,
                                           #ndim//2-vol_size//2:ndim//2+vol_size//2 + vol_size%2])

