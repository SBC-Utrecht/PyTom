#!/usr/bin/env pytom

import numpy
import random
import os
import sys
import glob
from pytom.basic.structures import ParticleList, Rotation
from pytom.basic.files import read
from pytom_numpy import vol2npy
import matplotlib
matplotlib.use('Qt5Agg')
from pylab import *
from scipy.spatial.transform import Rotation as R
from pytom.agnostic.io import read_size

def get_size(particleList, directory):
    tempPL = ParticleList()
    tempPL.fromXMLFile(particleList)



    tomoName = tempPL[0].getPickPosition().getOriginFilename() #if tempPL[0].getPickPosition().getOriginFilename() else tempPL[0].getSourceInfo().getTomoName()

    if not os.path.exists(tomoName):
        tomoName = os.path.join(directory, tomoName)
        if not os.path.exists(tomoName):
            return 'Failed'

    try:
        dimx, dimy, dimz = read_size(tomoName)
    except:
        print('Failed')
        return 'Failed'

    return [dimx,dimy,dimz]

def mirrorParticleList(particleList, outname, directory='./'):

    sizes = get_size(particleList, directory)

    if sizes == 'Failed':
        print('Mirroring particle coordinates did not succeed. Please ensure the paths to the origin tomogram are correct')
        return
    dimx,dimy,dimz = sizes

    tempPL = ParticleList()
    tempPL.fromXMLFile(particleList)
    for particle in tempPL:
        pp = particle.getPickPosition()
        pp.setX(dimx - pp.getX())
        pp.setY(dimy - pp.getY())
        pp.setZ(dimz - pp.getZ())
        shift = particle.getShift()
        shift.invert()

    tempPL.toXMLFile(outname)

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

    return z1 - rotation_angle, x, z2

def updatePL(fnames, outnames, directory='', suffix='', wedgeangles=[], multiplypickpos=1, multiplyshift=None,
             new_center=[], sizeSubtomo=64, move_shift=False, binSubtomo=1, binRecon=1, rotation=[],
             anglelist='', mirror=False,  tomogram_dir='./', convention='zxz', scalePP=None):
    if type(fnames) == str:
        fnames = [fnames]
    if type(outnames) == str:
        outnames = [outnames]

    try: wedgelen = len(wedgeangles)
    except: wedgelen = 0

    for n, xmlfile in enumerate(fnames):
        tempPL = ParticleList()
        tempPL.fromXMLFile(xmlfile)

        for pid, particle in enumerate(tempPL):

            if pid == 0:
                tomoName = tempPL[0].getPickPosition().getOriginFilename() #if tempPL[
                    # 0].getPickPosition().getOriginFilename() else tempPL[0].getSourceInfo().getTomoName()

                if os.path.exists(tomoName):
                    shape = read_size(tomoName)
                else:
                    scalePP = None
                    print('no rescaling of Pick Positions is done')


            # Update directory to particle
            if directory:
                filename = os.path.join(directory, os.path.basename(particle.getFilename()))
                particle.setFilename(filename)

            # Add suffix to directory name in which particle is stored
            # e.g. suffix = _bin3
            #  --> Subtomograms/tomogram_000/particle_1.em will become Subtomograms/tomogram_000_bin3/particle_1.em
            if suffix:
                filename = particle.getFilename()
                filename = os.path.join( os.path.dirname(filename) + suffix, os.path.basename(filename))
                particle.setFilename(filename)

            # Update wedge angles of angle1 and angle2
            if wedgelen >  n + 1:
                w = particle.getWedge()
                w.setWedgeAngles(wedgeangles[n*2:n*2+2])

            # Multiply pick position
            if abs(multiplypickpos - 1) > 1E-3:
                pp = particle.getPickPosition()
                pp.scale(multiplypickpos)

            if not scalePP is None:
                pp = particle.getPickPosition()
                pos = pp.toVector()
                posnew = [(pos[0] - shape[0] // 2) * scalePP + shape[0] // 2,
                          (pos[1] - shape[1] // 2) * scalePP + shape[1] // 2,
                          (pos[2] - shape[2] // 2) * scalePP + shape[2] // 2]

            # Shift is multiply by the respective binning factor.
            if not (multiplyshift is None):
                shift = particle.getShift()
                shift.scale(multiplyshift)


            # Randomize the angles of all particles in particle list.
            if type(anglelist) == type(numpy.array([])):
                cc = 180. / numpy.pi
                import random
                z1, z2, x = random.choice(anglelist)
                particle.setRotation(rotation=Rotation(z1=z1 * cc, z2=z2 * cc, x=x * cc, paradigm='ZXZ'))

            shift = particle.getShift()
            angles = particle.getRotation().toVector(convention=convention)
            rot_particleList = R.from_euler(convention, angles, degrees=True)

            if new_center:
                new_center_vector = numpy.array(new_center) - sizeSubtomo//2
                new_center_vector_rotated = rot_particleList.apply(new_center_vector)
                shift.addVector( new_center_vector_rotated)

            if move_shift == True:
                pp = particle.getPickPosition()
                shift.scale( binSubtomo / binRecon)
                ss = shift.toVector()
                pp.setX(pp.getX() + ss[0])
                pp.setY(pp.getY() + ss[1])
                pp.setZ(pp.getZ() + ss[2])
                particle.setPickPosition(pp)
                shift.scale(0.)
                #print(particle)


            # Combine rotations from particleList and rotation
            if rotation:
                rot_rotation = R.from_euler(convention, rotation, degrees=True)
                combined_rotation = rot_particleList * rot_rotation
                z1, x, z2 = combined_rotation.as_euler(convention, degrees=True)
                particle.setRotation(rotation=Rotation(z1=z1, z2=z2, x=x, paradigm='ZXZ'))

            if mirror:
                # Update shifts
                shift.scale(-1)
                particle.setShift(shift)

                # Update angles as well
                rotationT = particle.getRotation()
                rotationT.setZ1(-1*rotationT.getZ1())
                rotationT.setZ2(-1*rotationT.getZ2())
                rotationT.setX(-1*rotationT.getX())
                particle.setRotation(rotationT)

        tempPL.toXMLFile(outnames[n])
        if mirror:
            mirrorParticleList(outnames[n], outnames[n], directory=tomogram_dir)

if __name__ == '__main__':
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.basic.structures import ParticleList

    options = [ScriptOption(['-o', '--outputName'], 'Output name of xml', True, False),
               ScriptOption(['-f', '--fileName'], 'particleList filesname.', True, False),
               ScriptOption(['-s', '--suffix'], 'Suffix placed behind last dirname before particle_??.em.', True, True),
               ScriptOption(['-d', '--subtomoDirectory'], 'Update directory of subtomogram reconstructions. If "particle_" is included in name, only dirname before prefix is considered', True, True),
               ScriptOption(['-w', '--wedgeAngles'],'Wedge angles for all particles. if one angle is given, both angles are updated with this angle.', True, True),
               ScriptOption(['-a', '--rotatePL'], 'Rotate the particle list according to either chimera output file or three euler angles (separated by ,).', True,True),
               ScriptOption(['-i', '--shiftToPickPos'], 'move the shift to pick position. The parameter sypplied is the binning factor to go from shift to pick position.', True, True),
               ScriptOption(['-h', '--help'], 'Help.', False, True)]

    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Update paramters in a particle list.',
                          authors='Gijs van der Schot',
                          options=options)

    if len(sys.argv) == 1:
        print(helper)
        sys.exit()

    try:
        outname, XMLfnames, suffix, prefix, wedgeangles, rotate, shiftBin, help = parse_script_options(sys.argv[1:], helper)
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

    if suffix:
        if not suffix.startswith('_'): suffix = '_'+suffix
    else:
        suffix=''

    if prefix: directory = prefix.split('/particle_')[0]
    else: directory=''
    if XMLfnames:
        fnames = XMLfnames.split(',' )

    try:
        if rotate and os.path.exists(rotate):
            chimeraFName = rotate
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



    updatePL(fnames, outname, directory=directory, wedgeangles=wedgeangles, suffix=suffix)

    print('Success: Particle List updated!')

