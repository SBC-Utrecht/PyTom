#!/usr/bin/env pytom


import os
import sys
import glob
from pytom.basic.structures import ParticleList


def parseChimeraOutputFile(chimeraOutputFile, ref_vector=[0, 0, 1], convention='zxz'):
    import os

    assert os.path.exists(chimeraOutputFile)



    try:
        vec = os.popen(
            "cat " + chimeraOutputFile + " | grep rotation_axis | head -1 | awk '{print $5, $4, $3}'").read()[:-1]
        rotation_vector = [float(line) for line in vec.split(',') if line]
        rotation_vector[0] *= -1
        rotation_vector[1] *= -1
        rotation_vector[2] *= -1

        rot = os.popen("cat " + chimeraOutputFile + " | grep rotation_angle | head -1 | awk '{print $2}'").read()[:-1]
        rotation_angle = float(rot.replace(',', ''))
        r = R.from_rotvec((numpy.array(rotation_vector)))  # -numpy.array(ref_vector)))
        z1, x, z2 = r.as_euler("zxz", degrees=True)
    except:
        raise Exception('Parsing chimera file failed.')

    print(z2, x, z1, rotation_angle)
    return z1 - rotation_angle, x, z2

def updatePL(fnames, outnames, directory='', wedgeangles=[]):
    if type(fnames) == str:
        fnames = [fnames]
    if type(outnames) == str:
        outnames = [outnames]

    try: wedgelen = len(wedgeangles)
    except: wedgelen = 0

    for n, xmlfile in enumerate(fnames):
        tempPL = ParticleList()
        tempPL.fromXMLFile(xmlfile)

        for particle in tempPL:
            if directory:
                filename = os.path.join(directory, os.path.basename(particle.getFilename()))
                particle.setFilename(filename)

            if not wedgelen >  n + 1:
                continue
            w = particle.getWedge()
            w.setWedgeAngles(wedgeangles[n*2:n*2+2])


        tempPL.toXMLFile(outnames[n])


if __name__ == '__main__':
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.basic.structures import ParticleList

    options = [ScriptOption(['-o', '--outputName'], 'Output name of xml', True, False),
               ScriptOption(['-f', '--fileName'], 'particleList filesname.', True, False),
               ScriptOption(['-d', '--subtomoDirectory'], 'Update directory of subtomogram reconstructions. If "particle_" is included in name, only dirname before prefix is considered', True, True),
               ScriptOption(['-w', '--wedgeAngles'],'Wedge angles for all particles. if one angle is given, both angles are updated with this angle.', True, True),
               ScriptOption(['-r', '--rotatePL'], 'Rotate the particle list according to either chimera output file or three euler angles (separated by ,).', True,True),
               ScriptOption(['-s', '--shiftToPickPos'], 'move the shift to pick position. The parameter sypplied is the binning factor to go from shift to pick position.'),
               ScriptOption(['-h', '--help'], 'Help.', False, True)]

    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Update paramters in a particle list.',
                          authors='Gijs van der Schot',
                          options=options)

    if len(sys.argv) == 1:
        print(helper)
        sys.exit()

    try:
        outname, XMLfnames, prefix, wedgeangles, rotate, shiftBin, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()

    if help is True:
        print(helper)
        sys.exit()

    fnames = []

    if wedgeangles:
        wedgeangles = wedgeangles.split(',')
        if len(wedgeangles) %2:
            wedgeangles = wedgeangles+[wedgeangles[-1]]
    else:
        wedgeangles= []

    try:
        wedgeangles = list(map(float,wedgeangles))
    except Exception as e:
        print("Wedge Angle Error: ", e)
        sys.exit()

    if prefix: directory = prefix.split('/particle_')[0]
    else: directory=''
    if XMLfnames:
        print(XMLfnames)
        fnames = XMLfnames.split(',' )

    try:
        if rotate and os.path.exists(rotate):
            chimeraFName = rotate
            par
        elif rotate != None:
            z1,x,z2 = list(map(float,rotate.split(',')))
        else:
            z1,x,z2 = 0,0,0
    except Exception as e:
        print('rotateError: ', e)
        #sys.exit()

    try:
        binningFactorRecon = float(shiftBin)
        binningFactorSubtomo = 1
    except Exception as e:
        print('binning factor error: ', e)



    updatePL(fnames, outname, directory=directory, wedgeangles=wedgeangles)

    print('Success: Particle List updated!')
