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

def extractParticlesOfCLassFromXML(xmlfile):
    
    #classIDS = list(map(int,classIDS))
    excludeList = []
    tomogram = ''
    pL = [1]
    while len(pL):
        tree = et.parse(xmlfile)
        tomogram = ''
        for n , particle in enumerate( tree.xpath("Particle") ):
            remove = True
            origin = particle.xpath('PickPosition')[0].get('Origin')
            if not origin:
                origin = particle.xpath('InfoTomogram')[0].get('TomoName')

            if not tomogram and not origin in excludeList: 
                tomogram = origin
                
            if not tomogram == origin or origin in excludeList:
                remove_element(particle)

        excludeList.append(tomogram)
        
        try:
            tree.write("particleList_{}.xml".format(os.path.basename(tomogram)), pretty_print=True)
        except:
            print('writing {} failed.'.format(fname))
            print('No file written.')

        pL = tree.xpath('Particle')
        print(os.path.basename(tomogram), len(pL))

if __name__ == '__main__':
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options

    options = [ScriptOption(['-p','--particleList'], 'Particle List', True, False),
               ScriptOption(['-h', '--help'], 'Help.', False, True)]


    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Convert particle list to n particle list based on tomoname.',
                          authors='Gijs van der Schot',
                          options=options)
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        plName, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()
    if help is True:
        print(helper)
        sys.exit()

    if not os.path.exists(plName):
        sys.exit()

    
    extractParticlesOfCLassFromXML(plName)
