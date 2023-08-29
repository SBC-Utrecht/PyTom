#!/usr/bin/env pytom
"""
    little script to convert coordinate list to particle list xml file
    FF Jan 2013
"""
from pytom.basic.structures import ParticleList, Rotation
from copy import deepcopy
import random
import numpy
from pytom.lib.pytom_volume import read
from pytom.lib.pytom_numpy import vol2npy
from pathlib import Path
import shutil
import os


def convertCoords2PL(coordinate_files, particleList_file, subtomoPrefix=None, wedge_angles=None, angleList=False, projDir='',
        tomogram_name=None):
    pl = ParticleList()
    for n, coordinate_file in enumerate(coordinate_files):
        if wedge_angles is not None:
            wedge_angle = wedge_angles[2*n:2*(n+1)]
        else:
            wedge_angle = None
        if subtomoPrefix is not None:
            name_prefix = subtomoPrefix[n]
        else:
            name_prefix = None
        infoGUI = pl.loadCoordinateFileHeader(coordinate_file)
        l2 = len(pl)
        pl.loadCoordinateFile(filename=coordinate_file, name_prefix=name_prefix, wedge_angle=wedge_angle,
                              infoGUI=infoGUI, projDir=projDir)

        try:
            cc = 180./numpy.pi
            for i in range(len(pl)-l2):
                z1, z2, x = random.choice(angleList)
                pl[-i-1].setRotation(rotation=Rotation(z1=z1*cc, z2=z2*cc, x=x*cc, paradigm='ZXZ'))
        except:
            pass
    if tomogram_name is not None:
        for particle in pl:
            particle.getPickPosition().setOriginFilename(tomogram_name)
    pl.toXMLFile(particleList_file)

def entry_point():
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options

    options = [ScriptOption(['-p','--particleList'], 'Particle List', True, False),
               ScriptOption(['-c','--coords'], 'Coordinate List (ascii file from EMAN2)', True, False),
               ScriptOption(['-s','--subtomoPrefix'],'path and filename for subtomogram files (e.g., MyPath/particle_)',
                            True, True),
               ScriptOption(['-w','--wedge_angles'], 'missing wedge angle(s) [counter-clock, clock] or single angle',
                            True, True),
               ScriptOption(['-r', '--randomizeParticleOrientation'], 'Randomize the orientation of the particles.', False, True),
               ScriptOption(['-a', '--angleList'], 'Randomize the rotations of the particles, using the supplied angle list (em/mrc format).', True, True),
               ScriptOption(['-t', '--tomogram_name'], 'OriginFilename to store for each of the particles', True, True),
               ScriptOption(['-h', '--help'], 'Help.', False, True)]


    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Convert coordinate list to particle list.',
                          authors='Friedrich Foerster',
                          options=options)
    if len(sys.argv) == 1:
        print(helper)
        return
    try:
        plName, coordName, subtomoPrefix, w, r, angleList, tomogram_name, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        raise e
    if help is True:
        print(helper)
        return
    if w:
        if len(w.split(',')) > 1:
            wedge_angle = []
            for kk in w.split(','):
                wedge_angle.append(float(kk))
        else:
            wedge_angle = [float(w)]
    else:
        wedge_angle = None

    if r:
        if not angleList:
            pytompath = Path(shutil.which('pytom')).parents[1]
            angleList = pytompath / 'angles/angleLists/angles_18_3040.em'

    if angleList:
        if not os.path.exists(angleList):
            raise Exception('Angle List is not existing.')
            angleList = None
        elif not angleList.split('.')[-1] in ('em','mrc'):
            raise Exception('File format should be mrc or em.')
            angleList = None
        else:
            vol = read(angleList)
            angleList = deepcopy(vol2npy(vol)).T.astype(numpy.float32)
            rotation = random.choice(angleList)
            if len(rotation) != 3 or type(rotation[0]) != numpy.float32:
                raise Exception('AngleList should contain three floats per row.')
                angleList = None

    coordName = coordName.split(',')
    if subtomoPrefix is not None:
        subtomoPrefix = subtomoPrefix.split(',')
        assert len(coordName) == len(subtomoPrefix)
    if isinstance(wedge_angle, list):
        assert len(wedge_angle) == len(coordName)*2
    convertCoords2PL(coordinate_files=coordName, particleList_file=plName, subtomoPrefix=subtomoPrefix,
                     wedge_angles=wedge_angle,angleList=angleList, tomogram_name=tomogram_name)

if __name__ == '__main__':
    entry_point()
