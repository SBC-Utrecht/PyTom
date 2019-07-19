#!/usr/bin/env pytom

'''
Created on Sep 12, 2012

@author: yuxiangchen
'''

if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.basic.structures import Particle, ParticleList, SingleTiltWedge
    import os
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Create a particle list from the EM files in a directory. Documentation is available at\n\
                          http://www.pytom.org/doc/pytom/genParticleList.html',
                          authors='Yuxiang Chen',
                          options=[ScriptOption(['-d'], 'Directory', True, False),
                                   ScriptOption(['-p'], 'Particle name prefix', True, True),
                                   ScriptOption(['-w'], 'Wedge Angle (degree)', True, False),
                                   ScriptOption(['-o'], 'Output particle list', True, False),
                                   ScriptOption(['-h', '--help'], 'Help.', False, True)])
    
    if len(sys.argv) == 1:
        print helper
        sys.exit()
    try:
        dir_name, name_prefix, wedge_angle, output, bHelp = parse_script_options(sys.argv[1:], helper)
    except:
        sys.exit()
    if bHelp is True:
        print helper
        sys.exit()
    
    if name_prefix is None:
        name_prefix = ''
    wedge_angle = float(wedge_angle)
    w = SingleTiltWedge(wedge_angle)
    
    pl = ParticleList()
    all_files = os.listdir(dir_name)
    for fname in all_files:
        p = None
        name_suffix = fname.split('.')[-1]
        if len(name_prefix):
            if name_prefix in fname and name_suffix in ['em', 'mrc']:
                p = Particle(dir_name+'/'+fname)
        else:
            if name_suffix in ['em', 'mrc']:
                p = Particle(dir_name+'/'+fname)
        if p is not None:
            pl.append(p)
    
    pl.setWedgeAllParticles(w)
    
    pl.toXMLFile(output)
    