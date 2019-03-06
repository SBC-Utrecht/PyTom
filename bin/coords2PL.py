#!/usr/bin/env pytom
"""
    little script to convert coordinate list to particle list xml file
    FF Jan 2013
"""
from pytom.basic.structures import ParticleList

def convertCoords2PL(coordinate_files, particleList_file, subtomoPrefix=None, wedgeAngles=None):
    pl = ParticleList()
    for n, coordinate_file in enumerate(coordinate_files):
        wedgeAngle = wedgeAngles[2*n:2*(n+1)]
        pl.loadCoordinateFile(filename=coordinate_file, name_prefix=subtomoPrefix[n], wedgeAngle=wedgeAngle)
    pl.toXMLFile(particleList_file)


if __name__ == '__main__':
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Convert coordinate list to particle list.',
                          authors='Friedrich Foerster',
                          options=[ScriptOption(['-p','--particleList'], 
                           'Particle List', True, False),
                                   ScriptOption(['-c','--coords'], 
                       'Coordinate List (ascii file from EMAN2)',
                       True, False),
                                   ScriptOption(['-s','--subtomoPrefix'], 
                       'path and filename for subtomogram files (e.g., MyPath/particle_)',
                       True, True),
                                   ScriptOption(['-w','--wedgeAngles'], 
                       'missing wedge angle(s) [counter-clock, clock] or single angle',
                       True, True),
                                   ScriptOption(['-h', '--help'], 'Help.', 
                       False, True)])
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        plName, coordName, subtomoPrefix, w, help = parse_script_options(sys.argv[1:], helper)
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

    coordName = coordName.split(',')
    subtomoPrefix = subtomoPrefix.split(',')
    assert len(coordName) == len(subtomoPrefix)
    assert lem(wedgeAngle) == len(coordName)*2
    convertCoords2PL(coordinate_files=coordName, particleList_file=plName, subtomoPrefix=subtomoPrefix,
                     wedgeAngles=wedgeAngle)

