#!/usr/bin/env pytom


import os
import sys
import glob

if __name__ == '__main__':
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options

    options = [ScriptOption(['-d','--directory'], 'folder with particle Lists', True, False),
               ScriptOption(['-o','--outputName'], 'Output name of xml', True, False),
               ScriptOption(['-f','--outputFileNames'], 'Output name of xml', True, False),
               ScriptOption(['-h', '--help'], 'Help.', False, True)]



    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Convert coordinate list to particle list.',
                          authors='Friedrich Foerster',
                          options=options)
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        directory, outname, outfnames, help = parse_script_options(sys.argv[1:], helper)
        
    except Exception as e:
        print(e)
        sys.exit()

    fnames = []

    if directory:
        fnames = [line for line in glob.glob(os.path.join(directory, '*.xml')) if line.endswith('.xml')]

    if outfnames:
        fnames = outfnames.split(',' )

    prefix = '''<!-- PyTom Version: 0.971 -->
<ParticleList>
''' 

    out = "" 

    suffix = '''
</ParticleList>
'''

    for fname in fnames:
        infile = open(fname, 'r')
        a = infile.read()

        for line in a.split('\n'):
            if not 'ParticleList' in line and not '<!' in line :
                out += line+'\n' 
        infile.close()

    outfile = open(outname,'w' )
    outfile.write(prefix+out+suffix)
    outfile.close()
