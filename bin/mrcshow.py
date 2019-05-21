#!/usr/bin/env pytom

import sys
from pytom.gui.mrcOperations import read_mrc, downsample
from pytom_volume import read
from pytom_numpy import vol2npy
import matplotlib
matplotlib.use('Qt5Agg')
from pylab import imshow, show, subplots
import copy

if sys.argv[1].endswith('.em'):
    vol = read(sys.argv[1])
    data = copy.deepcopy(vol2npy(vol)).T
else:
    import mrcfile

    data = copy.deepcopy( mrcfile.open(sys.argv[1],permissive=True).data)

print( 'max: {} min: {} mean: {} std: {}'.format( data.max(), data.min(), data.mean(), data.std()))
print(data.shape)
from scipy.ndimage import gaussian_filter

dd = downsample(data,4)
dd[dd>dd.mean()+5*dd.std()] = dd.mean()

fig,ax = subplots(1,1,figsize=(10,10))
ax.imshow( dd, cmap='gray')
show()
#data.close()
