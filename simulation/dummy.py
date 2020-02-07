from math import sqrt
from simulateProjections import *
import matplotlib
matplotlib.use('Qt5Agg')
from pylab import *
#
# V = 200e3
# h = 6.626068e-34
# e = 1.60217646e-19
# m = 9.10938188e-31
# c = 2.99792458e8
#
# print(h/sqrt(e*V*m*(e/m*V/c**2 + 2 )))
# # print(wavelength_eV2nm(V))
# print(6.625 * 10 ** (-34) / ((2 * 9.1 * 10 ** (-31)) * 1.6 * (10 ** (-19)) * V * 1) ** 0.5)
#
#
# ctf, amp = create_ctf(-2, zeros((512,512)), 0.175E-9, voltage=200E3, Cs=2.7E-3, sigma=0.4)
#
# imshow(ctf)
# show()

import configparser
from pytom.gui.guiFunctions import loadstar, datatype


config = configparser.ConfigParser()
config.read_file(open('dummy.conf'))

print(config.sections())

if 'General' in config.sections():
    outputFolder    = config['General']['OutputFolder']
    modelID         = int(config['General']['ModelID'])

    # meta file
    metadata        = loadstar(config['General']['MetaFile'],dtype=datatype)
    angles          = metadata['TiltAngle']
    defocus         = metadata['DefocusU'][0]*1E-6
    voltage         = metadata['Voltage'][0]*1E3
    sphericalAberration = metadata['SphericalAberration'][0]*1E-3
    amplitudeContrast   = metadata['AmplitudeContrast'][0]
    pixelSize           = metadata['PixelSpacing'][0]*1E-10
else:
    raise Exception('No general parameters specified in the config file.')


if 'GenerateModel' in config.sections():
    generateModel       = True
    models              = eval(config['GenerateModel']['Models'])
    SIZE                = int(config['GenerateModel']['Size'])
    waterdensity        = float(config['GenerateModel']['WaterDensity'])
    numberOfParticles   = int(config['GenerateModel']['NumberOfParticles'])
else:
    generateModel = False

if 'GenerateProjections' in config.sections():
    generateProjections = True
    multislice      = bool(config['GenerateProjections']['MultiSlice'])
    msdz            = float(config['GenerateProjections']['MultiSliceSize'])*1E-9
    heightBox       = int(config['GenerateProjections']['HeightBox'])
else:
    generateProjections = False

if 'AddEffectsMicroscope' in config.sections():
    addEffectsMicroscope = True
    SNR             = float(config['AddEffectsMicroscope']['SNR'])
    sigmaDecayCTF   = float(config['AddEffectsMicroscope']['SigmaDecayCTF'])
else:
    addEffectsMicroscope = False

if 'ReconstructTomogram' in config.sections():
    reconstructTomogram = True
    prefix      = config['ReconstructTomogram']['Prefix']
    suffix      = config['ReconstructTomogram']['Suffix']
    start_idx   = int(config['ReconstructTomogram']['StartIdx'])
    end_idx     = int(config['ReconstructTomogram']['EndIdx'])
    weighting   = int(config['ReconstructTomogram']['Weighting'])
else:
    reconstructTomogram = True



