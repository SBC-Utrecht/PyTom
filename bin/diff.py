#!/usr/bin/env pytom

from pytom.tompy.io import read
import sys
import numpy as np

a = read(sys.argv[1])
b = read(sys.argv[2])


a = np.angle(np.fft.fftshift(np.fft.fftn(a)))
b = np.angle(np.fft.fftshift(np.fft.fftn(b)))


d = np.abs(a-b)

print(d.mean(), d.std(), d.min(), d.max())

import matplotlib

matplotlib.use('Qt5Agg')
from pylab import imshow, show

d = d.squeeze()
print(d.shape)
imshow(d)
show()
