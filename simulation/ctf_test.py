from simulateProjections import calcCTF
import numpy as np

import matplotlib
matplotlib.use('Qt5Agg')
from pylab import *



ctf = calcCTF(3E-6, np.zeros((70,70)), 1E-9, voltage=300E3, Cs=2.7E-3, sigma_decay_ctf=0, amplitude_contrast=0.08)

# ax,fig

imshow(ctf.real)
show()
imshow(ctf.imag)
show()