#!/usr/bin/env pytom
from numpy import arange, meshgrid, float32, zeros, sqrt, exp
import sys
import mrcfile



def create_ellipse(size, mj, mn1, mn2, outname, smooth, cutoff_SD=3):
    X,Y,Z = meshgrid(arange(size/1), arange(size/1), arange(size/1))

    X -= size/2-0.5
    Y -= size/2-0.5
    Z -= size/2-0.5

    R = sqrt( (X/mj)**2 + (Y/mn1)**2 + (Z/mn2)**2)

    print(R.max(), R.min())

    out = zeros((size,size,size),dtype=float32)
    out[ R <= 1] = 1

    if smooth:
        R2 = R.copy()
        R2[R <= 1] = 1
        sphere = exp(-1 * ((R2-1)/smooth)**2)
        sphere[sphere <= exp(-cutoff_SD**2/2.)] = 0
        out = sphere

    if outname: 
        mrcfile.new(outname, out.astype(float32), overwrite=True)


if __name__ == '__main__':
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options

    options = [ScriptOption(['-e','--edgeSize'], 'Size of edge of box (pixels)', True, False),
               ScriptOption(['-m','--major'],  'Length of the major axis of the ellipse (pixels).', True, False),
               ScriptOption(['-n','--minor1'], 'Length of the first minor axis of the ellipse (pixels).',True, True),
               ScriptOption(['-l','--minor2'], 'Length of the smallest asis ellipse (pixels).',True, True),
               ScriptOption(['-o', '--outputFileName'], 'Name of output file.', True, False),
               ScriptOption(['-s', '--sigma'], 'SD of gaussian dropoff of smoothing kernel',True,True),
               ScriptOption(['-c', '--cuttoff'], 'Cutoff SD.', True, True),
               ScriptOption(['-h', '--help'], 'Help.', False, True)]


    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name                                                                                                                                             
                          description='Create 3D Ellipse at center of a cube of size (-s). \nOptional smoothing from the edge of the ellipse to zero using a gaussian dropoff with a cutoff defined by -c.',
                          authors='Gijs van der Schot',
                          options=options)
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        size, major, minor1, minor2, outname, sigma, cutoff, help = parse_script_options(sys.argv[1:], helper)
        
    except Exception as e:
        print(e)
        sys.exit()
    if help is True:
        print(helper)
        sys.exit()

    size = round(float(size))
    major = float(major)
    
    if not minor1: 
        minor1 = major
    else: 
        minor1 = int(minor1)
    
    if not minor2: 
        minor2 = minor1
    else: 
        minor2 = int(minor2)

    if cutoff: cutoff=float(cutoff)
    else: cutoff=3

    if sigma: sigma=float(sigma)
    else: sigma=0

    if not outname: 
        sys.exit()
    
    create_ellipse(size,major,minor1, minor2, outname, sigma, cutoff_SD=cutoff)
