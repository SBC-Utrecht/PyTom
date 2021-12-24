#!/usr/bin/env pytom

from numpy import arange, meshgrid, float32, zeros, sqrt, exp
import sys
import mrcfile


def create_ellipse(size, mj, mn1, mn2, smooth, cutoff_SD=3):
    assert all([size > 0, mj > 0, mn1 > 0, mn2 > 0, smooth > 0]), 'size, radii, or smooth are <= 0'

    X,Y,Z = meshgrid(arange(size/1), arange(size/1), arange(size/1))

    X -= size/2-0.5
    Y -= size/2-0.5
    Z -= size/2-0.5

    R = sqrt((X/mj)**2 + (Y/mn1)**2 + (Z/mn2)**2)

    out = zeros((size,size,size),dtype=float32)
    out[ R <= 1] = 1

    if smooth:
        R2 = R.copy()
        R2[R <= 1] = 1
        sigma = (smooth / ((mj + mn1 + mn2) / 3))
        sphere = exp(-1 * ((R2 - 1) / sigma) ** 2)
        sphere[sphere <= exp(-cutoff_SD**2/2.)] = 0
        out = sphere

    return out.astype(float32)


if __name__ == '__main__':
    from pytom.tools.script_helper import ScriptHelper2, ScriptOption2
    from pytom.tools.parse_script_options import parse_script_options2

    helper = ScriptHelper2(
        sys.argv[0].split('/')[-1],  # script name
        description='Create a spherical or ellipsoidal mask for matching or alignment !',
        authors='Gijs van der Schot, Marten Chaillet',
        options=[
            ScriptOption2(['-o', '--outputFileName'], 'Name of output file.', 'string', 'required'),
            ScriptOption2(['-b', '--boxSize'], 'Size of edge of box (in pixels)', 'int', 'required'),
            ScriptOption2(['-r', '--radius'], 'Length of the major axis of the ellipse (in pixels).', 'float',
                          'required'),
            ScriptOption2(['-m', '--minor'], 'Length of the first minor axis of the ellipse (in pixels), will be set '
                                             'equal to major axis if not provided.', 'float', 'optional'),
            ScriptOption2(['-n', '--minor2'], 'Length of the smallest asis ellipse (in pixels), will be set '
                                              'equal to major axis if not provided.', 'float', 'optional'),
            ScriptOption2(['-s', '--sigma'], 'SD of gaussian dropoff of smoothing kernel (in pixels)', 'float',
                          'optional', 0),
            ScriptOption2(['-c', '--cuttoff'], 'How many SDs (in pixels) will be sampled for the smoothing', 'int',
                          'optional', 3)])

    options = parse_script_options2(sys.argv[1:], helper)

    output_file_name, box_size, major, minor1, minor2, sigma, cutoff = options

    minor1 = major if minor1 is None else minor1
    minor2 = major if minor2 is None else minor2
    
    out = create_ellipse(box_size, major, minor1, minor2, sigma, cutoff_SD=cutoff)

    mrcfile.new(output_file_name, out.astype(float32), overwrite=True)
