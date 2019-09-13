#!/usr/bin/env pytom
"""
    little script to create new particle list of particles that belong to selected classes
    GvdS May 2019
"""
from pytom.basic.structures import ParticleList, Rotation

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

def extractParticlesOfCLassFromXML(xmlfile, classIDS, fname=''):
    tree = et.parse(xmlfile)
    try:
        classIDS = list(map(int,classIDS))
    except:
        print('converting class IDs to int failed. Please submit a valid class ID string.')
        return
    for particle in tree.xpath("Particle"):
        remove = True

        classID = int(particle.xpath('Class')[0].get('Name'))
        
        if not classID in classIDS:
            remove_element(particle)

    if fname:
        try:
            tree.write(fname.replace('.xml', '_deselected.xml'), pretty_print=True)
        except:
            print('writing {} failed.'.format(fname))
    else:
        print('No file written.')



if __name__ == '__main__':
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options

    options = [ScriptOption(['-p','--particleList'], 'Particle List', True, False),
               ScriptOption(['-o','--outputFileName'], 'Filename of output Particle List ', True, False),
               ScriptOption(['-c','--classes'],'Index of class: if more than one class separate by comma',
                            True, False),
               ScriptOption(['-h', '--help'], 'Help.', False, True)]


    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Convert coordinate list to particle list.',
                          authors='Friedrich Foerster',
                          options=options)
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        plName, outName, classIDS, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()
    if help is True:
        print(helper)
        sys.exit()

    if not os.path.exists(plName):
        sys.exit()

    classIDS = classIDS.split(',')

    extractParticlesOfCLassFromXML(plName,classIDS,outName)
