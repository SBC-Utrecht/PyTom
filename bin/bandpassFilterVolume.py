#! /usr/bin/env pytom

if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper2, ScriptOption2
    from pytom.tools.parse_script_options import parse_script_options2
    from pytom.basic.structures import Wedge
    from pytom.localization.extractPeaks import templateMatchingGPU
    from pytom.tompy.io import read, write
    from pytom.angles.globalSampling import GlobalSampling
    from pytom.bin.extractCandidates import extractCandidatesWithoutJobFile
    from pytom.score.score import FLCFScore
    import numpy
    from pytom.basic.structures import ParticleList
    import os

    helper = ScriptHelper2(
        sys.argv[0].split('/')[-1],  # script name
        description='Apply Bandpass Filter to a volume.',
        authors='GvdS',
        options=[ScriptOption2(['-v', '--volume'], 'File Name of volume', 'file', 'required'),
                 ScriptOption2(['-o', '--outputFileName'], 'Filename of output volume', 'string', 'required'),
                 ScriptOption2(['-l', '--lowestFrequency'], 'Lowest frequency', 'float', 'optional', 0.),
                 ScriptOption2(['-h', '--highestFrequency'], 'Highest frequency', 'float', 'optional', 1),
                 ScriptOption2(['-s', '--sigma'], 'Width of gaussian edge', 'float', 'optional', 0.),
                 ScriptOption2(['-g', '--gpuID'], 'Index of used GPU', 'int', 'optional'),
                 ])

    options = parse_script_options2(sys.argv[1:], helper)

    fname, outname, low, high, sigma, gpu = options

    from pytom.tompy.io import read, write
    from pytom.tompy.filter import bandpass

    a= read(fname)

    high = a.shape[0]//2 if high <= low else high

    b = bandpass(a,low=low,high=high,sigma=sigma)

    write(outname, b)
