#! /usr/bin/env pytom

def getBinningFactorAndReferenceMarker(volumeFileName):
    import os
    folder = os.path.dirname(os.popen(f'ls -alrt {volumeFileName}').read().split()[-1])

    binning, ref = 8,1

    try:
        if os.path.exists(os.path.join(folder, 'WBP_Reconstruction.sh')):
            data = open(os.path.join(folder, 'WBP_Reconstruction.sh'))
            dd = [line for line in data.readlines() if '--projectionBinning' in line or '--referenceMarkerIndex' in line]

            for a in dd:
                if '--projectionBinning' in a:
                    rr = a.split('--projectionBinning')[-1][1:].split(' ')
                    binning = int([i for i in rr if i][0])
                if '--referenceMarkerIndex' in a:
                    rr = a.split('--referenceMarkerIndex')[-1][1:].split(' ')
                    ref = int([i for i in rr if i][0])
    except Exception as e:
        print(e)
        pass

    return binning, ref

if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper2, ScriptOption2
    from pytom.tools.parse_script_options import parse_script_options2
    from pytom.basic.structures import Wedge
    from pytom.localization.extractPeaks import templateMatchingGPU
    from pytom.tompy.io import read, write
    from pytom.angles.globalSampling import GlobalSampling
    from pytom.bin.extractCandidates import extractCandidatesWithoutJobFile
    from pytom.score.score import FLCFScore
    import numpy
    from pytom.basic.structures import ParticleList
    import os

    helper = ScriptHelper2(
        sys.argv[0].split('/')[-1],  # script name
        description='Template matching and peak extraction on GPU.',
        authors='GvdS',
        options=[ScriptOption2(['--templateMatching'], 'Template matching', 'no arguments', 'optional'),
                 ScriptOption2(['-v', '--volume'], 'Filename of volume', 'file', 'optional'),
                 ScriptOption2(['-t', '--template'], 'Filename of template', 'file', 'optional'),
                 ScriptOption2(['-m', '--mask'], 'Filename of mask', 'file', 'optional'),
                 ScriptOption2(['-a', '--angles'], 'Filename of angle file defining the sampling', 'file', 'optional'),
                 ScriptOption2(['-w', '--wedgeAngles'], 'Missing wedge angles', 'float,float', 'optional',[30.,30.]),
                 ScriptOption2(['-s', '--scoreFunction'], 'Scoring function used', 'string', 'optional', 'FLCF'),
                 ScriptOption2(['-g', '--gpuID'], 'Index of used GPU', 'int', 'optional'),
                 ScriptOption2(['--startEndXYZ'], 'Region of interest in all dims. [startX, startY, startZ, endX, endY, endZ]', 'int,int,int,int,int,int', 'optional', [0,0,0,0,0,0]),

                 ScriptOption2(['--scoreVolume'], 'Filename of volume', 'file', 'optional', 'scores.mrc'),
                 ScriptOption2(['--rotationVolume'], 'Filename of volume', 'file', 'optional', 'angles.mrc'),

                 ScriptOption2(['--extractPeaks'], 'Extract peaks', 'no arguments', 'optional'),
                 ScriptOption2(['--particleList'], 'Filename of output particleList', 'string', 'optional', 'particleList.xml'),
                 ScriptOption2(['--particlePath'], 'Path of particles', 'string', 'optional', './'),
                 ScriptOption2(['-r', '--radiusParticle'], 'Radius of the particle', 'int', 'optional', 8),
                 ScriptOption2(['-n', '--numberOfParticles'], 'Number of particles', 'int', 'optional', 1),
                 ScriptOption2(['--minScore'], 'Minimum correlation score', 'float', 'optional', 0.1),
                 ScriptOption2(['--projectDir'], 'Project Directory (GUI)', 'directory', 'optional'),
                 ScriptOption2(['--NoOutputVolumes'], 'write intermediate output files', 'no arguments', 'optional')])

    options = parse_script_options2(sys.argv[1:], helper)

    templateMatch, volumeFileName, template, mask, angles, wedge_angles, scoreFunction, gpu, start_end, scores, rots, \
    extractPeaks, plFilename, particlePath, radius, num_particles, min_score, proj_dir, no_output = options

    start_end = numpy.array(start_end)

    if templateMatch:
        rotations = GlobalSampling(angles)
        volume = read(volumeFileName, keepnumpy=True)
        if start_end.sum():
            volume = volume[start_end[0]:start_end[3], start_end[1]:start_end[4], start_end[2]:start_end[5]]
            print(volume.shape)

        template = read(template, keepnumpy=True)
        mask = read(mask, keepnumpy=True)

        wedge = Wedge(wedge_angles)

        [scoreVolume, angleVolume, a,b] = templateMatchingGPU(volume, template, rotations, scoreFunction, mask,
                                                              wedgeInfo=wedge, gpuID=gpu)

        if not no_output:
            write(scores, scoreVolume)
            write(rots, angleVolume)

    if extractPeaks:

        binning, ref = getBinningFactorAndReferenceMarker(volumeFileName)

        res = extractCandidatesWithoutJobFile(volumeFileName, scores, rots, angles, radius, num_particles,
                                              min_score, scoringFunction=FLCFScore, offset=start_end[:3])

        if plFilename:

            pl = ParticleList()
            wedge = Wedge(*wedge_angles)

            if particlePath[-1] != os.sep:
                particlePath += os.sep

            for particle in res:
                newParticle = particle.toParticle()
                newParticle.setWedge(wedge)
                newParticle.setFilename(particlePath + newParticle.getFilename())
                if not proj_dir is None:
                    newParticle.getInfoGUI().setProjectDir(proj_dir)
                newParticle.getPickPosition().setBinningFactor(binning)
                newParticle.getPickPosition().setRefMarkerIndex(ref)
                pl.append(newParticle)

            pl.toXMLFile(plFilename)
