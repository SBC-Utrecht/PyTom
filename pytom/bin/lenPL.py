#!/usr/bin/env pytom

import mrcfile
import sys, os

if __name__=='__main__':
    
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.basic.structures import ParticleList

    options = [ScriptOption(['-f', '--fileName'], 'Filename of mrc stack.', True, False),
               ScriptOption(['-h', '--help'], 'Help.', False, True)]

    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name             
                          description='Extract tilt images from mrcstack, and creation of meta data file.',
                          authors='Gijs van der Schot',
                          options=options)

    if len(sys.argv) == 1:
        print(helper)
        sys.exit()

    try:
        filename, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()

    if help is True:
        print(helper)
        sys.exit()

    numberParticles = os.popen(f'grep Shift {filename} | wc -l').read()[:-1]

    print(f'number of particles in {os.path.basename(filename)}: ', numberParticles)
