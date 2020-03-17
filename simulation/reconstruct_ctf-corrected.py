from simulateProjections import *
# import matplotlib
# from pylab import *
# matplotlib.use('Qt5Agg')
import numpy as xp

modelID = 0
outputFolder = '/data2/mchaillet/simulation/pytom'
prefix = os.path.join(outputFolder, f'model_{modelID}/ctfCorrected/sorted_ctf_')
suffix = '.mrc'
vol = [512,512,512]

reconstruct_tomogram(prefix,suffix,0,40,vol,list(xp.arange(-60,63,3)),outputFolder,modelID,weighting=-1)