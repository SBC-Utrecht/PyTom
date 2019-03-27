#!/usr/bin/env pytom
"""
    little script to convert coordinate list to particle list xml file
    FF Jan 2013
"""
from pytom.basic.structures import ParticleList, Rotation

from copy import deepcopy
import random
import numpy
from pytom_volume import read
from pytom_numpy import vol2npy
import os


def convertCoords2PL(coordinate_files, particleList_file, subtomoPrefix=None, wedgeAngles=None, angleList=False):
    pl = ParticleList()
    for n, coordinate_file in enumerate(coordinate_files):
        wedgeAngle = wedgeAngles[2*n:2*(n+1)]
        pl.loadCoordinateFile(filename=coordinate_file, name_prefix=subtomoPrefix[n], wedgeAngle=wedgeAngle)
        try:
            z1, z2, x = random.choice(angleList)
            pl[-1].setRotation(rotation=Rotation(z1=z1, z2=z2, x=x, paradigm='ZXZ'))
        except:
            pass
    pl.toXMLFile(particleList_file)


if __name__ == '__main__':
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options

    options = [ScriptOption(['-p','--particleList'], 'Particle List', True, False),
               ScriptOption(['-c','--coords'], 'Coordinate List (ascii file from EMAN2)', True, False),
               ScriptOption(['-s','--subtomoPrefix'],'path and filename for subtomogram files (e.g., MyPath/particle_)',
                            True, True),
               ScriptOption(['-w','--wedgeAngles'], 'missing wedge angle(s) [counter-clock, clock] or single angle',
                            True, True),
               ScriptOption(['-r', '--randomizeParticleOrientation'], 'Randomize the orientation of the particles.', False, True),
               ScriptOption(['-a', '--angleList'], 'Randomize the rotations of the particles, using the supplied angle list (em/mrc format).', True, True),

               ScriptOption(['-h', '--help'], 'Help.', False, True)]


    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Convert coordinate list to particle list.',
                          authors='Friedrich Foerster',
                          options=options)
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        plName, coordName, subtomoPrefix, w, r, angleList, help = parse_script_options(sys.argv[1:], helper)
        print(plName, coordName, subtomoPrefix, w, r, angleList)
    except Exception as e:
        print(e)
        sys.exit()
    if help is True:
        print(helper)
        sys.exit()
    if w:
        if len(w.split(',')) > 1:
            wedgeAngle = []
            for kk in w.split(','):
                wedgeAngle.append(float(kk))
        else:
            wedgeAngle = float(w)
    else:
        wedgeAngle = None

    if r:
        if not angleList:
            pytompath = os.path.dirname(os.popen('dirname `which pytom`').read()[:-1])
            angleList = os.path.join(pytompath, 'angles/angleLists/angles_18_3040.em')

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
    subtomoPrefix = subtomoPrefix.split(',')
    assert len(coordName) == len(subtomoPrefix)
    assert lem(wedgeAngle) == len(coordName)*2
    convertCoords2PL(coordinate_files=coordName, particleList_file=plName, subtomoPrefix=subtomoPrefix,
                     wedgeAngles=wedgeAngle,angleList=angleList)
