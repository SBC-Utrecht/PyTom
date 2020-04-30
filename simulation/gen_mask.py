#!/usr/bin/env pytom
from pytom.tompy.io import read, write
import sys, os
import matplotlib

matplotlib.use('Qt5Agg')
from pylab import *
from skimage.morphology import *
from scipy.ndimage.morphology import binary_fill_holes
from scipy.ndimage import label
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage import median_filter


def gen_mask_fsc(data, num_cycles, outname=None, num_stds=1, smooth=2, maskD=None):
    '''This script generates a structured mask around the particle.
    @:param data: 3d array that constitutes the particle
    @:type data: ndarray of float32/64
    @:param num_cycles: number of binaray dilation cycles after thresholding.
    @:type num_cycles: int
    @:param outname: filename of output file. If not None, output file is written.
    @:type outname: str
    @:param num_stds: 3d ndarray
    @:type num_stds: float
    @:param smooth: 3d ndarray
    @:type smooth: float
    @:return return mask
    @:rtype ndarray of float32'''

    mask = zeros_like(data, dtype=int)
    print(data.std(), num_stds)
    mask[data > 0] = 1#data.mean() + float(num_stds) * data.std()] = 1
    mask = remove_small_objects(mask.astype(bool))
    mask = binary_fill_holes(mask)

    l, n = label(mask)



    part, total = 0, 0
    for i in range(1, n+1):
        if total < (l == i).sum():
            part = i
            total = (l == i).sum()

    mask = (l == part)
    if not maskD is None:
        mask *= maskD > 0

    mask = binary_dilation(mask)
    mask = median_filter(mask, 6)

    for i in range(num_cycles):
        mask = binary_dilation(mask)

    mask = gaussian_filter(mask * 100., smooth)
    mask[mask > 78] = 78.
    mask /= mask.max()



    if outname is None:
        return mask
    else:
        write(outname, mask.astype(float32))
        return mask


if __name__ == '__main__':

    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.basic.structures import ParticleList

    options = [ScriptOption(['-f', '--fileName'], 'Filename of model.', True, False),
               ScriptOption(['-o', '--outputName'], 'Filename of mask file.', True, False),
               ScriptOption(['-n', '--numStd'],
                            'The particle threshold is set to mean - stdev * numStd. Default value = 1', True, True),
               ScriptOption(['-s', '--smooth'],
                            'Smooth factor used to soften the edges of the mask (pixels). Default = 2', True, True),
               ScriptOption(['-c', '--numDilationCycles'], 'Number of binary dilation cycles. Default = 2', True, True),
               ScriptOption(['-m', '--mask'], 'Number of binary dilation cycles. Default = 2', True, True),
               ScriptOption(['-h', '--help'], 'Help.', False, True)]

    helper = ScriptHelper(sys.argv[0].split('/')[-1],  # script name
                          description='Extract tilt images from mrcstack, and creation of meta data file.',
                          authors='Gijs van der Schot',
                          options=options)

    if len(sys.argv) == 1:
        print(helper)
        sys.exit()

    try:
        filename, outname, numSTD, smooth, num_cycles, mask, help = parse_script_options(sys.argv[1:], helper)
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

    num_cycles = 2 if num_cycles is None else int(num_cycles)
    numStd = 1 if num_cycles is None else float(numSTD)
    smooth = 2 if smooth is None else float(smooth)

    if not mask is None:
        mask = read(mask)

    gen_mask_fsc(data, num_cycles, outname, numSTD, smooth, maskD=mask)
