from pytom.gui.mrcOperations import *
from simulateProjections import *
import numpy
import glob

folder = '/data2/mchaillet/simulation/pytom'

m = read_mrc(f'{folder}/ice/H2O_model.mrc')

print(f'average water electron density = {m[m>0.1].mean():5.3f}')
print(f'average water electron density = {(m[m>0.1]*0.94).mean():5.3f}')
print(f'std water electron density = {(m[m>0.1]*0.94).std():5.3f}')

print(m.shape)

# factor = 10
#
# size = max(m.shape)
# m2 = numpy.zeros((size, size, size))
# dx, dy, dz = m.shape
# sx, ex = (size - dx) // 2, size - int(numpy.ceil((size - dx) / 2.))
# sy, ey = (size - dy) // 2, size - int(numpy.ceil((size - dy) / 2.))
# sz, ez = (size - dz) // 2, size - int(numpy.ceil((size - dz) / 2.))
# m2[sx:ex, sy:ey, sz:ez] = m
# print(m2.shape)
# m_r = crop(m2, factor)
# print(m_r.shape)
#
# convert_numpy_array3d_mrc(m_r, f'{folder}/ice/H2O_model_binned.mrc')

for i in glob.glob(f'{folder}/particles/*.mrc'):
    n = read_mrc(i)
    print(f'average protein electron density {i.split("/")[-1]:30s} = { n[n>0].mean():5.3f}')
    print(f'std protein electron density {i.split("/")[-1]:30s} = { n[n>0].std():5.3f}')

#===================  output
# average water electron density =  0.013135817
# average protein electron density =  0.072276145