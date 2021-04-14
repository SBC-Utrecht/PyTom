#!/usr/bin/env pytom


import os
import sys
import glob

if __name__ == '__main__':
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.basic.structures import ParticleList
    from lxml import etree as et

    options = [ScriptOption(['-d','--directory'], 'folder with particle Lists', True, True),
               ScriptOption(['-o','--outputName'], 'Output name of xml', True, False),
               ScriptOption(['-f','--FileNames'], 'Series of xml filesnames joined by a comma.', True, True),
               ScriptOption(['-w', '--wedgeAngles'],'Up to two wedge angles per xml file joined by a comma.'+
                                                    'The Wedge angle is defined as the angle between 90 degrees and '+
                                                    'the first tilt angle.\n\n NB: only wedge parameters of particles '+
                                                    'from files for which there are two given wedgeangles are updated.'+
                                                    ' If you submit 4 particleList, 8 wedge angles have to be given.',
                            True, True),

               ScriptOption(['-h', '--help'], 'Help.', False, True)]



    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Combine several particle lists into one particle list.',
                          authors='Friedrich Foerster',
                          options=options)
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        directory, outname, XMLfnames, wedgeangles, help = parse_script_options(sys.argv[1:], helper)
        print(directory, outname, XMLfnames, wedgeangles)
    except Exception as e:
        print(e)
        sys.exit()


    if help is True:
        print(helper)
        sys.exit()


    fnames = []

    if wedgeangles:
        wedgeangles = wedgeangles.split(',')
    else:
        wedgeangle= []

    if directory:
        fnames = [line for line in glob.glob(os.path.join(directory, '*.xml')) if line.endswith('.xml')]

    if XMLfnames:
        print(XMLfnames)
        fnames = XMLfnames.split(',' )

    pl = ParticleList()

    if wedgeangles: wedgelen = len(wedgeangles)
    else: wedgelen = 0

    for n, xmlfile in enumerate(fnames):
        tempPL = ParticleList()
        tempPL.fromXMLFile(xmlfile)
        for particle in tempPL:
            if not wedgelen >  n + 1: continue
            w = particle.getWedge()
            w.setWedgeAngles(wedgeangles[n*2:n*2+2])
        pl += tempPL

    a = pl.toXMLFile(outname)
