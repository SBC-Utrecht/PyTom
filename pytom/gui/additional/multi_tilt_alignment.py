import os


def tiltalignment_all_markers(start, end, tiltSeriesName, targets, fname):
    """
    Aligns tilt images, looping over markers in markerfile as reference.
    OPTIONAL: procs defines the number of processes running parallel.
    """
    from multiprocessing import Process
    import time

    lines = open(fname).readlines()

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

    alignFolder = {'--fixMarkerPos': 'LocalMarkerRefinement', '': 'GlobalAlignment'}

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
            folders = [folder for folder in os.path.dirname(os.path.join(tomogram_names[index],
                                                                         logfile)).split('/') if folder]
            temp_f = ''
            for temp in folders:
                temp_f = os.path.join(temp_f, temp)
                if not os.path.exists(temp_f): os.mkdir(temp_f)

            # generate alignment is not run with the write options, which means the actual images wont be aligned and
            # written, only the alignmentresults will be written to the projectionsTargets folder
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

            cmd = cmd.format(tomogram_names[index], pytompath, tiltSeriesName, firstIndices[index],
                             lastIndices[index], refIndices[index], refmarkindex, markerFileNames[index], outdir,
                             weightingTypes[index], expectedRotationAngles[index], 1, fixMarkers[index],
                             refMarkIdFlag[index], logfile)

            while len(procs) > 19:
                time.sleep(1)
                procs = [proc for proc in procs if proc.is_alive()]
            print(cmd)
            p = Process(target=os.system, args=([cmd]))
            procs.append(p)
            p.start()
            print('job for {} with refMarkIndex {} started.'.format(tomogram_names[index], refmarkindex))


if __name__ == '__main__':
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options

    options = [ScriptOption(['--start'], 'Starting Index tomogram folder', arg=True, optional=False),
               ScriptOption(['--end'], 'Ending index tomogram folder', arg=True, optional=False),
               ScriptOption(['--tiltSeriesName'],
                            'Location of tilt images (such as sorted)',
                            arg=True, optional=False),
               ScriptOption(['--projectionTargets'],
                            'Relative or absolute path to the aligned projections that will be generated + file '
                            'prefix. default: "align/myTilt"',
                            arg=True, optional=False),
               ScriptOption(['--fnames'], 'File with tomogram names.', arg=True, optional=False),
               ScriptOption(['--help'], 'Help function.', arg=False, optional=True)]

    helper = ScriptHelper(sys.argv[0].split('/')[-1],
                          description='BATCHMODE: Align (and weight) projections taking each marker in markerfile as '
                                      'a reference, and save them in alignment/unweighted_unbinned_marker_XXX. \n',
                          authors='Gijs van der Schot',
                          options=options)

    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        start, end, tiltSeriesName, projectionTargets, fname, help = \
            parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        sys.exit()

    global pytompath
    pytompath = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

    start = int(start)
    end = int(end)
    tiltalignment_all_markers(start, end, tiltSeriesName, projectionTargets, fname)
