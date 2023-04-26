#!/usr/bin/env python

import mrcfile
import sys, os
from pytom.agnostic.io import read, write

if __name__=='__main__':
    
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.basic.structures import ParticleList

    options = [ScriptOption(['-f', '--fileName'], 'Filename of mrc stack.', True, False),
               ScriptOption(['-t', '--targetDir'], 'Folder in which the output files are saved', True, False),
               ScriptOption(['-h', '--help'], 'Help.', False, True)]

    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name             
                          description='Extract tilt images from mrcstack, and creation of meta data file.',
                          authors='Gijs van der Schot',
                          options=options)

    if len(sys.argv) == 1:
        print(helper)
        sys.exit()

    try:
        filename, outdir, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()

    if help is True:
        print(helper)
        sys.exit()


    if os.path.exists(filename):
        data = mrcfile.open(filename, permissive=True).data.copy()
        print(data.shape)
    else:
        print(helper)
        sys.exit()

    x,y,z = data.shape
    sliceNR = min(x,y,z)

    angles = range(-60,61,2)
    mdocname = os.path.join(outdir, os.path.basename(filename).replace('.mrc', '.mdoc'))
    mdocfile = open(mdocname, 'w')

    d = '''[ZValue = {}]
TiltAngle = {}
SubFramePath = X:\{}

'''

    for sliceId in range(sliceNR):
        tiltangle = angles[sliceId]
        outname = os.path.join(outdir, os.path.basename(filename).replace('.mrc', '_sorted_{:02d}.mrc'.format(sliceId)))
        #write(outname, data[:,:,sliceId], tilt_angle=tiltangle)
        mrcfile.new(outname, data[sliceId, :,:].astype('float32'), overwrite=True)
        mdocfile.write(d.format(sliceId, tiltangle, os.path.basename(outname)))

    mdocfile.close()
    
