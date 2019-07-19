#!/usr/bin/env pytom
'''
Created on Nov 8, 2013

@author: thrabe
'''

if __name__ == '__main__':
 # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options 
    from pytom.basic.files import pdb2em
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Compile a electron density from PDB file\n\
                          http://pytom.org/doc/pytom/files.html',
                          authors='Thomas Hrabe',
                          options=[ScriptOption(['-p','--pdbFile'], 'A PDB file', arg=True, optional=False),
                                   ScriptOption(['-c','--chain'], 'A Chain', arg=True, optional=True),
                                   ScriptOption(['-s','--pixelSize'], 'Pixel size of output volume (in Angstrom)', arg=True, optional=True),
                                   ScriptOption(['-v','--volumeSize'], 'Volume length (size) in all dimensions', arg=True, optional=True),
                                   ScriptOption(['-o','--outputVolumePath'], 'Path to output volume ', arg=True, optional=False),
                                   ScriptOption(['-i','--invertDensity'],'Set if density should be negative', arg=False, optional=False),
                                   ScriptOption(['-h', '--help'], 'Help.', arg=False, optional=True)])


    if len(sys.argv) == 1:
        print helper
        sys.exit()
    try:
        pdbFile, chain, pixelSize, cubeSize, volumePath ,densityNegative , helpme = parse_script_options(sys.argv[1:], helper)
    except:
        sys.exit()
        
    if helpme is True:
        print helper
        sys.exit()
        pass
    
    volume = pdb2em(pdbFile, float(pixelSize), int(cubeSize), chain = chain,densityNegative = densityNegative)
    
    volume.write(volumePath)