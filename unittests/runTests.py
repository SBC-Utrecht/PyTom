import sys
import os
import pickle, json
import numpy as np

global pytompath
pytompath = os.path.dirname(os.popen('dirname `which pytom`').read()[:-1])
print(pytompath)

if not pytompath:
    print('Pytom package is not available. Please load, or install Pytom.')
    sys.exit()

def update_env_vars(pytompath):
    '''Make sure all pytom functionality can be imported from within the script. '''


    if 0:
        from pytom_volume import read
    else:
        update_vars = False
        for search in ('LD_LIBRARY_PATH', 'PATH', 'PYTHONPATH'):
            # Check if env vars include all paths set in paths.csh
            query_string = "cat {}/bin/paths.csh | grep 'setenv {}' | grep -v '${}'".format(pytompath, search,
                                                                                            search)
            string = os.popen(query_string).read()[:-1].split()[-1]
            for new_lib in (string.split(':')):
                new_lib = new_lib.replace("'", "")

                if not new_lib in os.environ[search].split(':'):
                    os.environ[search] += ':' + new_lib
                    update_vars = True


        # If any of the env vars are updated reopen this script.
        if update_vars:
            
            if len(sys.argv) < 2:
                python3 = os.popen('which python3').read()[:-1]
                sys.argv = [python3] + sys.argv
            print(sys.argv)
            os.execv(sys.argv[0], sys.argv)
            # os.execv('/cm/shared/apps/python3/3.7/bin/python3.7', sys.argv)

update_env_vars(pytompath)

from pytom.unittests import pytom_Test

pytom_Test.run()

