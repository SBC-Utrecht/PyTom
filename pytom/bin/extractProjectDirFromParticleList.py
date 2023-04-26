#!/usr/bin/env python

'''
Created on Feb 6, 2021

@author: GvdS
'''


def extractProjectDirFromPL(pl_name, output_dir):
    from pytom.basic.structures import ParticleList
    import os

    pl = ParticleList()
    pl.fromXMLFile(pl_name)

    outputdir = output_dir if output_dir else os.path.dirname(pl_name)
    outputdir = './' if not output_dir else outputdir


    pls = pl.splitByProjectDir()



    if output:
        for plist in pls:
            d = os.path.basename(plist[0].getInfoGUI().getProjectDir())
            if d and d[-1] == '/':
                d = d[:-1]
            plist.toXMLFile(os.path.join(outputdir, f'{d}_{os.path.basename(pl_name)}'))
    else:
        return pls

if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.basic.structures import Particle, ParticleList, SingleTiltWedge
    import os
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Extract certain classes from the particle list.',
                          authors='Yuxiang Chen',
                          options=[ScriptOption(['-p'], 'Particle List', True, False),
                                   ScriptOption(['-d'], 'Output directory of particle lists', True, True),
                                   ScriptOption(['-h', '--help'], 'Help.', False, True)])
    
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        pl_name, output, bHelp = parse_script_options(sys.argv[1:], helper)
    except:
        sys.exit()
    if bHelp is True:
        print(helper)
        sys.exit()

    output = './' if output is None else output
    
    extractProjectDirFromPL(pl_name, output)
    
