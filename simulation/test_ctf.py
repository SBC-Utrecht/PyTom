import matplotlib
matplotlib.use('Qt5Agg')
from pylab import imshow, show
import numpy as np
from simulateProjections import calcCTF

pixel = 0.5 * 1E-9
defocus = 1E-6

ctf = calcCTF(np.zeros((128,128)), pixel, defocus)

imshow(- ctf.real - ctf.imag)
show()