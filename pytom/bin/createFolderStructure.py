#!/usr/bin/env pytom

'''
Created on Jan 17, 2020

@author: GvdS
'''
if __name__ == '__main__':
    # parse command line arguments
    import sys, os
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.gui.guiFunctions import create_project_filestructure

    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Create a folder with specific folder structure, as used by the GUI.',
                          authors='GvdS',
                          options=[ScriptOption(['-d'], 'Name of Folder.', True, False),
                                   ScriptOption(['-h', '--help'], 'Help.', False, True)])

    if len(sys.argv) == 1:
        print(helper)
        sys.exit()

    try:
        dir_name, bHelp = parse_script_options(sys.argv[1:], helper)
    except:
        print(helper)
        sys.exit()

    if bHelp == True:
        print(helper)
        sys.exit()


    if os.path.exists(dir_name):
        raise Exception('Folder already exists. Please give non-existent foldername.')

    create_project_filestructure(dir_name)
