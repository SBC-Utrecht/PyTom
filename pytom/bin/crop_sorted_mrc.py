#!/usr/bin/env pytom

import os
import sys
import mrcfile
from numpy import array, shape


if __name__=='__main__':
    folder = sys.argv[1]
    backup = os.path.join(folder,'backup')

    if not os.path.exists(backup): os.mkdir(backup)

    os.system('mv {}/sorted*.mrc {}'.format(folder, backup))
    os.system('cp {}/markerfile.em .bb_mark/{}_markerfile.em'.format(folder, folder.replace('/','_')))
    os.system('mv {}/markerfile.em {}'.format(folder, backup))

    a  = [line for line in os.listdir(backup) if line.startswith('sorted') and line.endswith('.mrc')]

    for mrc in a:

        a = mrcfile.open(os.path.join(backup,mrc),permissive=True)
        data = a.data[:,:]
        a.close()
        shapeF = data.shape
        size = min(shapeF)
        datao = data[shapeF[0]//2-size//2:shapeF[0]//2+size//2, shapeF[1]//2-size//2:shapeF[1]//2+size//2]

        mrcfile.new(os.path.join(folder,mrc), datao, overwrite=True)
    
    if os.path.exists('{}/markerfile.em'.format(backup)):
        from pytom.basic.files import read, write_em
        from pytom_numpy import vol2npy, npy2vol
        import copy

        volume = read('{}/markerfile.em'.format(backup))

        marker = copy.deepcopy(vol2npy(volume)).T

        marker[:,:,1] -= max(64, abs(shapeF[0]-shapeF[1])//2)
        marker[:,:,1][marker[:,:,1] < -1] = -1

        markerVolume = npy2vol(array(marker.T, dtype='float32', order='F'), 3)
        
        write_em('{}/markerfile.em'.format(folder), markerVolume)
