from pytom.tompy.mpi import MPI
from pytom.gui.guiFunctions import loadstar, datatype
from voltools import transform
from pytom.gui.mrcOperations import *
import configparser
import numpy as xp



def rotate_volume(filename, outname, angle):

    print('rotating volume to ', angle)
    volume = read_mrc(filename)
    rotated_volume = transform(volume, rotation=(0, angle, 0), rotation_order='szyx', interpolation='filt_bspline', device='cpu')
    convert_numpy_array3d_mrc(rotated_volume, outname)


def mpi_rotation(grandcell, angles, heightBox, outputFolder, modelID):

    size = grandcell.shape[1]

    if grandcell.shape[0] >= heightBox:
        raise Exception('Your model is larger than than the box that we will rotate (heightBox parameter)')

    volume = xp.zeros((heightBox, size, size), dtype=xp.float32)

    offset = (heightBox - grandcell.shape[0]) // 2

    if (heightBox - grandcell.shape[0]) % 2:
        volume[offset + 1:-offset, :, :] = grandcell[:, :, :]
    else:
        volume[offset:-offset, :, :] = grandcell[:, :, :]

    filename = f'{outputFolder}/model_{modelID}/rotated_volume_0.mrc'
    convert_numpy_array3d_mrc(volume, filename)

    params = []
    angles2 = angles.copy()
    angles2.remove(0.)
    for ang in angles2:
        outfilename = f'{outputFolder}/model_{modelID}/rotated_volume_{int(ang)}.mrc'
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

    grandcell = read_mrc(f'{outputFolder}/model_{modelID}/grandmodel_{modelID}.mrc')

    print(grandcell.shape)

    mpi_rotation(grandcell, angles, heightBox, outputFolder, modelID)


