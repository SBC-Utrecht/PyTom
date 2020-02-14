from pytom.gui.mrcOperations import *
import numpy

folder = '/data2/mchaillet/simulation/pytom/'

m = read_mrc(f'{folder}ice/H2O_model.mrc')
print('average water electron density = ', m.mean())

n = read_mrc(f'{folder}particles/test_3cf3_model.mrc')
print('average protein electron density = ', n.mean())

#===================  output
# average water electron density =  0.013135817
# average protein electron density =  0.072276145