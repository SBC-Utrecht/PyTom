#!/usr/bin/env python

from multiprocessing import Process
import os
import sys


def tiltalignment_all_markers(start, end, procs, tiltSeriesName, firstIndex, lastIndex, refIndex, markerFileName,
                              targets, weightingType, tomogramFolder, fnames, projIndices=True,
                              expectedRotationAngle=0):
    '''Aligns tilt images, looping over markers in markerfile as reference. OPTIONAL: procs defines the number of processes running parallel.'''

    lines =  open(fnames).readlines()

    tomogram_names = [line.split()[0] for line in lines]
    refIndices = [line.split()[1] for line in lines]
    referenceMarkerIndex = [line.split()[3] for line in lines]
    numMarkers = [line.split()[2] for line in lines]
    firstAngles = [line.split()[4] for line in lines]
    lastAngles = [line.split()[5] for line in lines]
    firstIndices = [line.split()[6] for line in lines]
    lastIndices = [line.split()[7] for line in lines]
    weightingTypes = [line.split()[8] for line in lines]
    expectedRotationAngles = [line.split()[9] for line in lines]
    tiltSeriesNames = [line.split()[10] for line in lines]
    markerFileNames = [line.split()[11] for line in lines]
    fixMarkers = ['--fixMarkerPos' * ('True' in line.split()[12]) for line in lines]
    refMarkIdFlag = [f'--refMarkIdTomo {line.split()[13]}' * (line.split()[13] != '*') for line in lines]
    for t in tomogram_names:
        if not os.path.exists(os.path.join(t, 'alignment')):
            os.mkdir(os.path.join(t, 'alignment'))


    procs = []


    alignFolder  = {'--fixMarkerPos': 'LocalMarkerRefinement', '': 'GlobalAlignment'}

    for index in range(start, end):

        string = '{}/marker_{}_{},{}/{}/{}'
        outdir = string.format(targets, '__', firstAngles[index], lastAngles[index], alignFolder[fixMarkers[index]],
                               os.path.basename(tiltSeriesName))
        tiltSeriesName = tiltSeriesNames[index]
        if referenceMarkerIndex[index] == 'all':
            refmarks = range(int(numMarkers[index]))
        else:
            refmarks = [int(referenceMarkerIndex[index])]

        for refmarkindex in refmarks:
            logfile = '{}/logfile.alignment.{}.txt'.format(outdir.replace('____', '_{:04d}_'.format(int(refmarkindex))),
                                                           os.path.basename(tiltSeriesName))


            # Check if output folder exists. If not, create all missing folders and subfolders
            folders = [folder for folder in os.path.dirname(os.path.join(tomogram_names[index], logfile)).split('/') if folder]
            temp_f = ''
            for temp in folders:
                temp_f = os.path.join(temp_f, temp)
                if not os.path.exists(temp_f): os.mkdir(temp_f)

            cmd = '''cd {}; generateAlignedTiltImages.py \
    --tiltSeriesName {}  \
	--firstIndex {} \
	--lastIndex {} \
	--referenceIndex {} \
	--referenceMarkerIndex {} \
	--markerFile {} \
	--projectionTargets {} \
	--projectionBinning 1 \
	--lowpassFilter 0.9 \
	--weightingType {} \
	--expectedRotationAngle {} \
    --numberProcesses {} {} {} > {}'''

            cmd = cmd.format(tomogram_names[index], tiltSeriesName, firstIndices[index],
                             lastIndices[index], refIndices[index], refmarkindex, markerFileNames[index], outdir,
                             weightingTypes[index], expectedRotationAngles[index], 1, fixMarkers[index],
                             refMarkIdFlag[index], logfile)

            while len(procs) > 19.1:
                time.sleep(1)
                procs = [proc for proc in procs if proc.is_alive()]
            print(cmd)
            p = Process(target=os.system, args=([cmd]))
            procs.append(p)
            p.start()
            print('job for {} with refMarkIndex {} started.'.format(tomogram_names[index], refmarkindex))


def extract_subtomograms(start, end, fnames, folders, binning_tomogram=1, binning_subtomograms=1, offset='0,0,0',
                         size=80, weighting_type=0):
    for i in range(start, end):
        cmd = 'reconstructWB.py -p {} --projectionDirectory {} -s {} -b {} -o {} --applyWeighting {}'.format(fnames[i],
                                                                                                             folders[i],
                                                                                                             size,
                                                                                                             binning_tomogram,
                                                                                                             offset,
                                                                                                             weighting_type)
        print(cmd)
        p = Process(target=os.system, args=([cmd]))
        p.start()


if __name__ == '__main__':
    import sys
    # from pytom.reconstruction.TiltAlignmentStructures import TiltAlignmentParameters, TiltSeries, TiltAlignment
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.reconstruction.reconstructionFunctions import alignWeightReconstruct
    from pytom_volume import read
    from pytom_numpy import vol2npy
    import os
    from multiprocessing import Process
    import time

    options = [ScriptOption(['--start'], 'Starting Index tomogram folder', arg=True, optional=False),
               ScriptOption(['--end'], 'Ending index tomogram folder', arg=True, optional=False),
               ScriptOption(['--numberProcesses'], 'Number of processes spawned by each process', arg=True,
                            optional=True),
               ScriptOption(['--tiltSeriesName'],
                            'Name tilt series - either prefix of sequential tilt series files expected as "tiltSeriesName_index.em/mrc"',
                            arg=True, optional=False),
               ScriptOption(['--firstIndex'], 'Index of first projection.', arg=True, optional=True),
               ScriptOption(['--lastIndex'], 'Index of last projection.', arg=True, optional=True),
               ScriptOption(['--referenceIndex'], 'Index of reference projection used for alignment.', arg=True,
                            optional=True),
               ScriptOption(['--markerFile'],
                            'Name of EM markerfile or IMOD wimp File containing marker coordinates.', arg=True,
                            optional=True),
               ScriptOption(['--projectionTargets'],
                            'Relative or absolute path to the aligned projections that will be generated + file prefix. default: "align/myTilt"',
                            arg=True, optional=False),
               ScriptOption(['--weightingType'], 'Type of weighting (-1 default r-weighting, 0 no weighting)',
                            arg=True, optional=True),
               ScriptOption(['--tomogramFolder'], 'Folder in which tomogram_XXX is located', arg=True,
                            optional=True),
               ScriptOption(['--fnames'], 'File with tomogram names.', arg=True, optional=False),
               ScriptOption(['--projIndices'], 'Use projection indices.', arg=False, optional=True),
               ScriptOption(['--expectedRotationAngle'], 'Use projection indices.', arg=True, optional=True),
               ScriptOption(['--help'], 'Help function.', arg=False, optional=True)
               ]

    helper = ScriptHelper(sys.argv[0].split('/')[-1],
                          description='BATCHMODE: Align (and weight) projections taking each marker in markerfile as a reference, and save them in alignment/unweighted_unbinned_marker_XXX. \n',
                          authors='Gijs van der Schot',
                          options=options)

    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    if 1:
        start, end, procs, tiltSeriesName, firstProj, lastProj, referenceIndex, markerFileName, projectionTargets, \
        weightingType, tomogramfolder, fnames, projIndices, expectedRotationAngle, help = input = \
            parse_script_options(sys.argv[1:], helper)
    else:  # Exception as e:
        sys.exit()

    start = int(start)
    end = int(end)
    print(start, end)
    tiltalignment_all_markers(start, end, procs, tiltSeriesName, firstProj, lastProj, referenceIndex,
                              markerFileName,
                              projectionTargets, weightingType, tomogramfolder, fnames, projIndices=projIndices,
                              expectedRotationAngle=expectedRotationAngle)
