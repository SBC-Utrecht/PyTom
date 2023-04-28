#!/usr/bin/env pytom
from pytom.agnostic.io import read, write
import sys, os
try:
    import matplotlib
    matplotlib.use('Qt5Agg')
    from pylab import *
except:
    import matplotlib
    from pylab import *

from skimage.morphology import *
from scipy.ndimage.morphology import binary_fill_holes
from scipy.ndimage import label
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage import median_filter

if __name__=='__main__':

    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.basic.structures import ParticleList

    options = [ScriptOption(['-f', '--fileName'], 'Filename of model.', True, False),
               ScriptOption(['-o', '--outputName'], 'Filename of mask file.', True, False),
               ScriptOption(['-b', '--binaryDilationCycles'], 'Number of binary dilation cycles.', True, False),
               ScriptOption(['-h', '--help'], 'Help.', False, True)]

    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name             
                          description='Extract tilt images from mrcstack, and creation of meta data file.',
                          authors='Gijs van der Schot',
                          options=options)

    if len(sys.argv) == 1:
        print(helper)
        sys.exit()

    try:
        filename, outname, num_cycles, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()

    if help is True:
        print(helper)
        sys.exit()


    if os.path.exists(filename):
        data = read(filename)
    else:
        print(helper)
        sys.exit()

    if num_cycles is None:
        num_cycles = 0
    else:
        num_cylces = int(num_cycles)

    
    data = read(filename)
    mask = zeros_like(data, dtype=int)
    mask[data < data.mean()-data.std()] = 1
    mask = remove_small_objects(mask.astype(bool))
    mask = binary_fill_holes(mask)

    l, n = label(mask)

    dx,dy,dz = data.shape
    mask = (l==l[dx//2,dy//2,dz//2])
    mask = binary_dilation(mask)
    mask = median_filter(mask, 6)
 
    for i in range(num_cycles+2):
        mask = binary_dilation(mask)

    mask = gaussian_filter(mask*100.,2)
    mask[mask>78] = 78.
    mask /= mask.max()

    write(outname, mask.astype(float32))

