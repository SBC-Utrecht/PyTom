#!/usr/bin/env pytom
"""
    little script to create particle lists based on originating tomogram and closest marker. Supply particleList and
    log file from reconstruction with marker positions (as written by gui).
    GvdS May 2019
"""
from pytom.basic.structures import ParticleList, Rotation
from pytom.gui.guiFunctions import datatypeMR
from copy import deepcopy
import random
import numpy
from pytom_volume import read
from pytom_numpy import vol2npy
import os
import lxml.etree as et
from pytom.gui.guiFunctions import loadstar

def correct_header_markerfile(markerfile):
    a = open(markerfile, 'r')
    txt = a.read()
    a.close()

    nothashed = '''MarkerIndex 
OffsetX (px)
OffsetY (px)
OffsetZ (px)
PositionX (px)
PositionY (px)
PositionZ (px)'''

    hashed = '''# MarkerIndex 
# OffsetX (px)
# OffsetY (px)
# OffsetZ (px)
# PositionX (px)
# PositionY (px)
# PositionZ (px)
# '''

    b = open(markerfile, 'w')
    b.write(txt.replace(nothashed, hashed))
    b.close()

def determineShiftXY(marker1, marker2):
    print([marker1[name].mean()-marker2[name].mean() for name in ('PositionX', 'PositionY')])
    print([marker1[name].std()-marker2[name].std() for name in ('PositionX', 'PositionY')])
    return [marker1[name].mean()-marker2[name].mean() for name in ('PositionX', 'PositionY')]

def determine_closest_marker(x,y,z,markers):
    ''' Determines the closest markerpoint to particle at position x, y, z.
    Uses a markerset formatted in numpy data structure defined in datatypeMR (see guiFunctions).'''

    refmarkers = list(zip(markers['PositionX'], markers['PositionY'],markers['PositionZ']))

    markerIndex = 0
    dist = 10000
    for n, (rx, ry, rz) in enumerate(refmarkers):
        #print(rx, ry, rz, x, y, z)
        tempdist = numpy.sqrt((rx-x)**2+(ry-y)**2+(rz-z)**2)
        if tempdist < dist:
            markerIndex = n
            dist = tempdist

    return markerIndex

def extractParticleListsClosestToRefMarker(xmlfile, markerfile, binning_factor=8, directory='./', projDirTemplate=''):
    from pytom.basic.structures import PickPosition, ParticleList
    pL = ParticleList()
    pL.fromXMLFile(os.path.join(directory, xmlfile))

    dict_particle_lists = {}

    for particle in pL:
        tomogram = particle.getPickPosition().getOriginFilename().split('/')[-1].split('.')[0]
        if not tomogram:
            tomogram = particle.getSourceInfo().getTomoName().split('.')[0]
        if not tomogram	in dict_particle_lists.keys():
            dict_particle_lists[tomogram] = ParticleList()
        dict_particle_lists[tomogram].append(particle)


    try: 
        markers = loadstar(markerfile, dtype=datatypeMR)
    except:
        correct_header_markerfile(markerfile)
        markers = loadstar(markerfile, dtype=datatypeMR)

    xmlsCM = []

    for pl_key in dict_particle_lists.keys():
        outLists = {}
        for particle in dict_particle_lists[pl_key]:
            x, y, z = particle.getPickPosition().toVector()
            x *= binning_factor
            y *= binning_factor
            z *= binning_factor

            closestMarkerIndex = determine_closest_marker(x,y,z, markers)
            projectionDirectory = projDirTemplate.replace('_CLOSEST_', '_{:04d}_'.format(closestMarkerIndex))
            markerPositionFile = f'{projectionDirectory}/markerLocations_irefmark_{closestMarkerIndex}.txt'

            realignmarkers = loadstar(markerPositionFile, dtype=datatypeMR)

            if not closestMarkerIndex in outLists.keys():
                outLists[closestMarkerIndex] = ParticleList()

            ox,oy = determineShiftXY(markers, realignmarkers)

            oz = markers['OffsetZ'][closestMarkerIndex]
            originFname = particle.getPickPosition().getOriginFilename()

            print(x,y,z, ox,oy, oz)
            pp = PickPosition(x=(x-ox)/binning_factor,y=(y-oy)/binning_factor,z=((z-oz)/binning_factor), originFilename=originFname)
            particle.setPickPosition(pp)
            outLists[closestMarkerIndex].append(particle)

        for markerIndex in outLists.keys():
            outfname = '.tempCM_particleList_{}_refMarkerIndex_{}.xml'.format(pl_key, markerIndex)
            outfname = os.path.join(directory, outfname)
            outLists[markerIndex].toXMLFile(outfname)
            xmlsCM.append([markerIndex, outfname])


    return xmlsCM

if __name__ == '__main__':
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options

    options = [ScriptOption(['-p','--particleList'], 'Particle List', True, False),
               ScriptOption(['-l','--logfileReconstruction'], 'Particle List', True, False),
               ScriptOption(['-b','--binningFactor'], 'Binning Factor for reconstruction.', True, False),
               ScriptOption(['-h', '--help'], 'Help.', False, True)]


    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Convert particle list to n particle list based on tomoname.',
                          authors='Gijs van der Schot',
                          options=options)
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        plName, logfile, binningFactor, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()

    if help is True:
        print(helper)
        sys.exit()

    if not binningFactor:
        binningFactor = 8
    else:
        try:
            binningFactor = int(binningFactor)
        except:
            print('Invalid binningFactor. Exit.')
            sys.exit()

    if not plName or not os.path.exists(plName):
        print('particleList does not exist. Exit. ')
        sys.exit()
    if not logfile or not os.path.exists(logfile):
        print('logfile does not exist. Exit.')
        sys.exit()

    fnames = extractParticleListsClosestToRefMarker(plName, logfile)

    
