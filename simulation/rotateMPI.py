from pytom.tompy.mpi import MPI
from pytom.gui.guiFunctions import loadstar, datatype
from voltools import transform
import pytom.tompy.io
import configparser
import numpy as xp


def rotate_volume(filename, outname, angle):

    print('rotating volume to ', angle)
    volume = pytom.tompy.io.read_mrc(filename)
    # either 'linear' or 'filt_bspline'
    rotated_volume = transform(volume, rotation=(0, angle, 0), rotation_order='sxyz', interpolation='filt_bspline', device='cpu')
    pytom.tompy.io.write(outname, rotated_volume)


def mpi_rotation(angles, heightBox, outputFolder, modelID):

    grandcell = pytom.tompy.io.read_mrc(f'{outputFolder}/model_{modelID}/grandmodel_{modelID}.mrc')

    print(grandcell.shape)

    size = grandcell.shape[0]
    height = grandcell.shape[2]

    if grandcell.shape[2] >= heightBox:
        raise Exception('Your model is larger than than the box that we will rotate (heightBox parameter)')

    volume = xp.zeros((size, size, heightBox), dtype=xp.float32)

    offset = (heightBox - height) // 2

    if (heightBox - height) % 2:
        volume[:, :, offset + 1:-offset] = grandcell[:, :, :]
    else:
        volume[:, :, offset:-offset] = grandcell[:, :, :]

    dir = f'{outputFolder}/model_{modelID}/rotations'

    filename = f'{dir}/rotated_volume_0.mrc'
    pytom.tompy.io.write(filename, volume)

    params = []
    angles2 = angles.copy()
    angles2.remove(0.)
    for ang in angles2:
        outfilename = f'{dir}/rotated_volume_{int(ang)}.mrc'
        params.append((filename, outfilename, ang))

    mpi = MPI()
    mpi.begin()
    mpi.parfor(rotate_volume, params)
    mpi.end()



if __name__ == '__main__':

    config = configparser.ConfigParser()
    try:
        config.read_file(open('simulation.conf'))
    except Exception as e:
        print(e)
        raise Exception('Could not open config file.')

    outputFolder = config['General']['OutputFolder']
    modelID = int(config['General']['ModelID'])
    metadata = loadstar(config['General']['MetaFile'], dtype=datatype)
    angles = list(metadata['TiltAngle']) # specified in degrees
    heightBox = int(config['Rotation']['heightBox'])

    mpi_rotation(angles, heightBox, outputFolder, modelID)


