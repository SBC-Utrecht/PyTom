from math import sqrt
from simulateProjections import *
import matplotlib
matplotlib.use('Qt5Agg')
from pylab import *
import numpy as xp

# V = 300e3
# h = 6.626068e-34
# e = 1.60217646e-19
# m = 9.10938188e-31
# c = 2.99792458e8
#
# print(h/sqrt(e*V*m*(e/m*V/c**2 + 2 )))
# # print(wavelength_eV2nm(V))
# print(6.625 * 10 ** (-34) / ((2 * 9.1 * 10 ** (-31)) * 1.6 * (10 ** (-19)) * V * 1) ** 0.5)


 # = create_ctf(2E-6, zeros((512,512)), 1E-9, voltage=300E3, Cs=2.7E-3, sigma=0.4)
#
# ctf = calcCTF(2E-6, zeros((512,512)), 1E-9, voltage=300E3, Cs=2.7E-3, sigma_decay_ctf=0, amplitude_contrast=0.08)
#
# imshow(np.abs(ctf))
# show()

dzprop = 2E-9
Lambda = wavelength_eV2m(300E3)
pixelSize = 1E-9
SIZE = 512

xwm = pixelSize * SIZE *100  # pixelsize for multislice * size sample
q_true_pix_m = 1 / xwm
q_m = rr(SIZE, SIZE).astype(xp.complex128) * q_true_pix_m

# FRESNEL PROPAGATOR
# P_x = pi * Lambda * (q_m**2) * dzprop
# P = xp.cos(P_x) + 1j * xp.sin(P_x)
P = xp.exp(-1j * pi * Lambda * (q_m**2) * dzprop)

imshow(abs((P))**2)
show()

print(P[256,256])
print(P[1,1])

im = xp.fft.ifftn((P))*512

imshow(abs(im))
show()

print(im[256,256])
print(im.max())