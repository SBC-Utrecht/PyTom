#!/usr/bin/env python

'''
Created on Mar 20, 2012

@author: yuxiangchen
'''

if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Select the particle list for alignment according to the class list.',
                          authors='Yuxiang Chen',
                          options=[ScriptOption('-p', 'Particle list.', True, False),
                                   ScriptOption('-c', 'Class label file.', True, False),
                                   ScriptOption('-o', 'Output particle list.', True, False),
                                   ScriptOption('-t', 'True positive class.', True, True),
                                   ScriptOption(['-h', '--help'], 'Help.', False, True)])
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        pl_filename, class_label_filename, output, tp_label , help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()
    if help is True:
        print(helper)
        sys.exit()
    
    if tp_label is None:
        tp_label = 1
    else:
        tp_label = int(tp_label)
    
    from pytom.basic.structures import ParticleList
    pl = ParticleList(".")
    pl.fromXMLFile(pl_filename)
    
    pl.setClassFromLocalizationClassList(class_label_filename)
    new_pl = pl.particlesFromClass(tp_label)
    new_pl.toXMLFile(output)