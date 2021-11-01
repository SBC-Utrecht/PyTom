#!/usr/bin/env pytom

from pytom.agnostic.io import read
import sys
import numpy as np

a = (read(sys.argv[1]).squeeze())
b = (read(sys.argv[2]).squeeze())


if len(sys.argv) > 3 and sys.argv[3] != 'plot' and int(sys.argv[3]) == 1:
    a = np.fft.fftshift(a,axes=(0,1))
    b = np.fft.fftshift(b, axes=(0,1))

try:
    mask = read('mask.mrc')
    a = a * mask
    b = b * mask
except Exception as e:
    print(e)
    mask = np.ones_like(a)
#a = np.angle(np.fft.fftshift(np.fft.fftn(a)))
#b = np.angle(np.fft.fftshift(np.fft.fftn(b)))

d = np.abs(a-b)

for e in [d]:
    print(e[mask > 1E-4].mean(), e[mask>1E-4].std(), e.min(), e.max())

try:
    import matplotlib
    matplotlib.use('Qt5Agg')
    from pylab import *
except:
    import matplotlib
    from pylab import *


d = d.squeeze()

if 'plot' in sys.argv:
    fig,ax = subplots(1,3,figsize=(15,5))
    try:
        if 1:
            ax[0].imshow(a[a.shape[0]//2,:,:])
            ax[1].imshow(b[a.shape[0]//2,:,:])
            ax[2].imshow(d[a.shape[0]//2,:,:])
        else:
            ax[0].imshow(a.sum(axis=2))
            ax[1].imshow(b.sum(axis=2))
            ax[2].imshow(d.sum(axis=2))
    except:
        ax[0].imshow(a[:,:])
        ax[1].imshow(b[:,:])
        ax[2].imshow(d[:,:])

    ax[0].set_title(sys.argv[1])
    ax[1].set_title(sys.argv[2])

    show()