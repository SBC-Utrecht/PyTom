#!/cm/shared/apps/python3/3.7/bin/python3.7
import sys
import os
import pickle, json
import numpy as np

if sys.version_info[0] < 3:
    print(sys.version_info[0])
    raise Exception("The GUI requires Python 3")

global pytompath
pytompath = os.path.dirname(os.popen('dirname `which pytom`').read()[:-1])
if not pytompath: pytompath = '/Users/gijs/Documents/pytom_private'

if not pytompath:
    print('Pytom package is not available. Please load, or install Pytom.')
    sys.exit()

def update_env_vars(pytompath):
    '''Make sure all pytom functionality can be imported from within the script. '''
    if 0:
        from pytom_volume import read
    else:
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



from filter import exactFilter

import matplotlib
matplotlib.use('Qt5Agg')




def exactFilter(tilt_angles, tiltAngle, sX, sY, sliceWidth):
    """
    exactFilter: Generates the exact weighting function required for weighted backprojection - y-axis is tilt axis
    Reference : Optik, Exact filters for general geometry three dimensional reconstuction, vol.73,146,1986.
    @param tilt_angles: list of all the tilt angles in one tilt series
    @param titlAngle: tilt angle for which the exact weighting function is calculated
    @param sizeX: size of weighted image in X
    @param sizeY: size of weighted image in Y

    @return: filter volume

    """

    from pytom_volume import read, vol, vol_comp
    import numpy as np

    sizeY = (sY // 2) + 1
    IijFunc_sum = [0.0] * sX
    IijFunc = []
    weightFunc = vol(sX, sizeY, 1)
    weightFunc.setAll(0.0)


    for n in range(1, len(tilt_angles) + 1):
        theta = tilt_angles[n - 1]
        reltheta = abs(theta - tiltAngle)
        if reltheta == 0:
            continue
        ds = float(sliceWidth / sX)
        fijh = float(1. / (ds * np.sin(reltheta * np.pi / 180)))
        for ix in range(-sX // 2, sX // 2):
            freq = ix
            if 0 <= freq <= fijh:
                Iij = 1. - (freq / fijh)
                IijFunc.append(Iij)
            elif -fijh <= freq <= 0:
                Iij = 1. + (freq / fijh)
                IijFunc.append(Iij)
            else:
                Iij = 0.0
                IijFunc.append(Iij)
        IijFunc_sum = np.array(IijFunc_sum) + np.array(IijFunc)
        IijFunc = []
    InterFunc = 1. + (IijFunc_sum)
    wf = 1. / InterFunc
    for ix in range(0, sX):
        for iy in range(0, sizeY):
            weightFunc.setV(wf[ix], ix, iy, 0)

    return weightFunc


#import matplotlib
#matplotlib.use('Qt5Agg')


sX,sY = 1000,1000
tilt_angles= range(-60,60,2)
tiltAngle = 0
sliceWidth = 1000


from pytom_numpy import vol2npy
from copy import deepcopy


d = exactFilter(tilt_angles, tiltAngle, sX, sY, sliceWidth)

data = deepcopy( vol2npy(d) )

print( data.shape, data.sum() )

from pylab import *

plot(data)
show()
