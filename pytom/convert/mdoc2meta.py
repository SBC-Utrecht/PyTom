#!/usr/bin/env pytom
import matplotlib
import matplotlib.backends.backend_qt5agg

try:
    matplotlib.use('Qt5Agg')
except:
    pass
import sys
import os
import pickle, json
import numpy as np

if sys.version_info[0] < 3:
    print(sys.version_info[0])
    raise Exception("The GUI requires Python 3")

global pytompath
pytompath = os.path.dirname(os.popen('dirname `which pytom`').read()[:-1])

if not pytompath:
    print('Pytom package is not available. Please load, or install Pytom.')
    sys.exit()

def update_env_vars(pytompath):
    '''Make sure all pytom functionality can be imported from within the script. '''
    try:
        from pytom_volume import read
    except:
        update_vars = False
        for search in ('LD_LIBRARY_PATH','PATH','PYTHONPATH'):
            # Check if env vars include all paths set in paths.csh
            query_string = "cat {}/bin/paths.csh | grep 'setenv {}' | grep -v '${}'".format(pytompath, search,search)
            string = os.popen(query_string).read()[:-1].split()[-1]
            for new_lib in (string.split(':')):
                new_lib = new_lib.replace("'","")

                if not new_lib in os.environ[search].split(':'):
                    os.environ[search] += ':'+new_lib
                    update_vars = True
        #If any of the env vars are updated reopen this script.
        if update_vars:
            if len(sys.argv) < 2:
                pythonVersion = 'python{d[0]}.{d[1]}'.format( d=sys.version_info )
                path = os.popen('which {}'.format(pythonVersion)).read()[:-1]
                sys.argv = [path] + sys.argv

            os.execv(sys.argv[0],sys.argv)
            #os.execv('/cm/shared/apps/python3/3.7/bin/python3.7', sys.argv)
update_env_vars(pytompath)

from pylab import *
from pytom.gui.guiFunctions import datatype, headerText, units, fmt, createMetaDataFiles

def mdoc2meta(path, mdocfiles, target, mdoc_only=True):
    createMetaDataFiles(path, mdocfiles, target,mdoc_only=mdoc_only)


if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Convert mdoc file(s) to GUI meta file(s).',
                          authors='Gijs van der Schot',
                          options=[ScriptOption(['-f', '--file'], 'Filename', True, True),
                                   ScriptOption(['-d', '--directory'], 'A directory of files.', True, True),
                                   ScriptOption(['-t', '--targetPath'], 'Path to new file.', True, True),
                                   ScriptOption(['-h', '--help'], 'Help.', False, True)])
    
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        filename, directory, target, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()
    if help is True:
        print(helper)
        sys.exit()
    
    if filename:
        #convert only one file
        mdoc2meta(os.path.dirname(filename),[os.path.basename(filename)],target,mdoc_only=True)

    elif directory:

        fileList = [mdocfile for mdocfile in os.listdir(directory) if mdocfile.endswith('.mdoc')]
        mdoc2meta(directory, fileList, target, mdoc_only=True)
