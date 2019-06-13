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


def remove_element(el):
    parent = el.getparent()
    if el.tail.strip():
        prev = el.getprevious()
        if prev:
            prev.tail = (prev.tail or '') + el.tail
        else:
            parent.text = (parent.text or '') + el.tail
    parent.remove(el)


def determine_closest_marker(x,y,z,markers):
    ''' Determines the closest markerpoint to particle at position x, y, z.
    Uses a markerset formatted in numpy data structure defined in datatypeMR (see guiFunctions).'''

    refmarkers = list(zip(markers['PositionX'], markers['PositionY'],markers['PositionZ']))

    markerIndex = 0
    dist = 10000
    for n, (rx, ry, rz) in enumerate(refmarkers):
        tempdist = numpy.sqrt((rx-x)**2+(ry-y)**2+(rz-z)**2)
        if tempdist < dist:
            markerIndex = n

    return markerIndex

def extractParticleListsClosestToRefMarker(xmlfile, markerfile, directory='./'):
    from pytom.basic.structures import PickPosition, ParticleList
    pL = ParticleList()
    pL.fromXMLFile(xmlfile)

    dict_particle_lists = {}

    for particle in pL:
        fname = particle.getFilename()
        tomogram = 'tomogram_' + fname.split('tomogram_')[1][:3]
        if not tomogram	in dict_particle_lists.keys():
            dict_particle_lists[tomogram] = ParticleList()
        dict_particle_lists[tomogram].append(particle)


    markers = numpy.loadtxt(markerfile, dtype=datatypeMR)

    xmlsCM = []

    for pl_key in dict_particle_lists.keys():
        outLists = {}
        for particle in dict_particle_lists[pl_key]:
            x, y, z = particle.getPickPosition().toVector()
            closestMarkerIndex = determine_closest_marker(x,y,z, markers)


            if not closestMarkerIndex in outLists.keys():
                outLists[closestMarkerIndex] = ParticleList()

            ox = markers['OffsetX'][closestMarkerIndex]
            oy = markers['OffsetY'][closestMarkerIndex]
            oz = markers['OffsetZ'][closestMarkerIndex]
            originFname = particle.getOriginFilename()
            particle.setPickPosition(PickPosition(x=x+ox,y=y+oy,z=z+oz,originFilename=originFname))
            outLists[closestMarkerIndex].append(particle)
        for markerIndex in outLists.keys():
            outfname = '.tempCM_particleList_{}_refMarkerIndex_{}.xml'.format(pl_key, markerIndex)
            outLists[markerIndex].write(os.path.join(directory, outfname), pretty_print=True)
            xmlsCM.append(outfname)

    return xmlsCM

if __name__ == '__main__':
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options

    options = [ScriptOption(['-p','--particleList'], 'Particle List', True, False),
               ScriptOption(['-l','--logfileReconstruction'], 'Particle List', True, False),
               ScriptOption(['-h', '--help'], 'Help.', False, True)]


    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Convert particle list to n particle list based on tomoname.',
                          authors='Gijs van der Schot',
                          options=options)
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        plName, logfile, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()

    if help is True:
        print(helper)
        sys.exit()

    if not plName or not os.path.exists(plName):
        print('particleList does not exist. Exit. ')
        sys.exit()
    if not logfile or not os.path.exists(logfile):
        print('logfile does not exist. Exit.')
        sys.exit()

    extractParticleListsClosestToRefMarker(plName, logfile)

