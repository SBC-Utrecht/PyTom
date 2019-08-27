from multiprocessing import Process
import os
import sys 

global pytompath
pytompath = os.path.dirname(os.popen('dirname `which pytom`').read()[:-1])
#pytompath = '/data/gijsvds/pytom-develop/pytom_python3/pytom'
import time

def tiltalignment_all_markers(start, end, procs, tiltSeriesName, firstIndex, lastIndex, refIndex, markerFileName,
                              targets, weightingType, tomogramFolder, fnames, projIndices=True, expectedRotationAngle=0):
    '''Aligns tilt images, looping over markers in markerfile as reference. OPTIONAL: procs defines the number of processes running parallel.'''


    tomogram_names = [line.split()[0] for line in open(fnames).readlines()]
    refIndices = [line.split()[1] for line in open(fnames).readlines()]
    referenceMarkerIndex = [line.split()[3] for line in open(fnames).readlines()]
    numMarkers = [line.split()[2] for line in open(fnames).readlines()]
    firstAngles = [line.split()[4] for line in open(fnames).readlines()]
    lastAngles = [line.split()[5] for line in open(fnames).readlines()]
    firstIndices = [line.split()[6] for line in open(fnames).readlines()]
    lastIndices = [line.split()[7] for line in open(fnames).readlines()]
    weightingTypes = [line.split()[8] for line in open(fnames).readlines()]
    expectedRotationAngles = [line.split()[9] for line in open(fnames).readlines()]
    tiltSeriesNames = [line.split()[10] for line in open(fnames).readlines()]
    for t in tomogram_names:
        if not os.path.exists(os.path.join(t,'alignment')):
            os.mkdir(os.path.join(t,'alignment'))


    procs = []

    for index in range(start,end):

        string = '{}/marker_{}_{},{}'
        outdir = string.format(targets, '__', firstAngles[index], lastAngles[index])
        tiltSeriesName = tiltSeriesNames[index]
        if referenceMarkerIndex[index] == 'all':
            refmarks = range(int(numMarkers[index]))
        else:
            refmarks = [int(referenceMarkerIndex[index])]

        for refmarkindex in refmarks:
            logfile = '{}/logfile.alignment.{}.txt'.format(outdir.replace('____', '_{:04d}_'.format(int(refmarkindex))),
                                                           os.path.basename(tiltSeriesName))

            if not os.path.exists(os.path.dirname(os.path.join(tomogram_names[index], logfile))):
                os.mkdir(os.path.dirname(os.path.join(tomogram_names[index], logfile)))

            cmd = '''cd {}; pytom {}/gui/additional/generateAlignedTiltImages.py \
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
    --numberProcesses {} > {}'''

            cmd = cmd.format(tomogram_names[index], pytompath, tiltSeriesName, firstIndices[index],
                  lastIndices[index], refIndices[index], refmarkindex, markerFileName, outdir,
                  weightingTypes[index], expectedRotationAngles[index], 1, logfile)

            while len(procs) > 19.1:
                time.sleep(1)
                procs = [proc for proc in procs if proc.is_alive()]

            p = Process(target=os.system,args=([cmd]))
            procs.append(p)
            p.start()
            print('job for {} with refMarkIndex {} started.'.format(tomogram_names[index], refmarkindex))

def extract_subtomograms(start,end,fnames,folders,binning_tomogram=1,binning_subtomograms=1,offset='0,0,0',size=80,weighting_type=0):
    for i in range(start,end):

        cmd = 'reconstructWB.py -p {} --projectionDirectory {} -s {} -b {} -o {} --applyWeighting {}'.format(fnames[i],folders[i],size,binning_tomogram,offset,weighting_type)
        print(cmd)
        p = Process(target=os.system,args=([cmd]))
        p.start()


if __name__=='__main__':
    import sys
    #from pytom.reconstruction.TiltAlignmentStructures import TiltAlignmentParameters, TiltSeries, TiltAlignment                                                                           
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.gui.additional.reconstructionFunctions import alignWeightReconstruct
    from pytom_volume import read
    from pytom_numpy import vol2npy
    import os
    from multiprocessing import Process
    import time

    try:
        options=[ScriptOption(['--start'],'Starting Index tomogram folder', arg=True,optional=False),
                 ScriptOption(['--end'],'Ending index tomogram folder', arg=True,optional=False),
                 ScriptOption(['--numberProcesses'],'Number of processes spawned by each process', arg=True,optional=False),
                 ScriptOption(['--tiltSeriesName'], 'Name tilt series - either prefix of sequential tilt series files expected as "tiltSeriesName_index.em/mrc"', arg=True, optional=False),
                 ScriptOption(['--firstIndex'], 'Index of first projection.', arg=True, optional=True),
                 ScriptOption(['--lastIndex'], 'Index of last projection.', arg=True, optional=True),
                 ScriptOption(['--referenceIndex'], 'Index of reference projection used for alignment.', arg=True, optional=True),
                 ScriptOption(['--markerFile'], 'Name of EM markerfile or IMOD wimp File containing marker coordinates.', arg=True, optional=False),
                 ScriptOption(['--projectionTargets'], 'Relative or absolute path to the aligned projections that will be generated + file prefix. default: "align/myTilt"', arg=True, optional=False),
                 ScriptOption(['--weightingType'], 'Type of weighting (-1 default r-weighting, 0 no weighting)', arg=True,optional=False),
                 ScriptOption(['--tomogramFolder'], 'Folder in which tomogram_XXX is located', arg=True,optional=False),
                 ScriptOption(['--fnames'], 'File with tomogram names.', arg=True,optional=False),
                 ScriptOption(['--projIndices'], 'Use projection indices.', arg=False, optional=True),
                 ScriptOption(['--expectedRotationAngle'], 'Use projection indices.', arg=True, optional=True),
                 ScriptOption(['--help'],'Help function.',arg=False,optional=True)
        ]

        helper = ScriptHelper(sys.argv[0].split('/')[-1],
                          description='BATCHMODE: Align (and weight) projections taking each marker in markerfile as a reference, and save them in alignment/unweighted_unbinned_marker_XXX. \n',
                          authors='Gijs van der Schot',
                          options = options)

        if len(sys.argv) == 1:
            print(helper)
            sys.exit()
        try:
            start, end, procs, tiltSeriesName, firstProj, lastProj, referenceIndex, markerFileName, projectionTargets, \
            weightingType,tomogramfolder,fnames, projIndices, expectedRotationAngle, help = input = \
                parse_script_options(sys.argv[1:], helper)
        except:# Exception as e:
            sys.exit()


        start = int(start)
        end = int(end)

        tiltalignment_all_markers(start, end, procs, tiltSeriesName, firstProj, lastProj, referenceIndex,markerFileName,
                                  projectionTargets, weightingType, tomogramfolder, fnames, projIndices=projIndices,
                                  expectedRotationAngle=expectedRotationAngle)


    except:
        options=[ScriptOption(['--start'],'Starting Index tomogram folder', arg=True,optional=False),
                 ScriptOption(['--end'],'Ending index tomogram folder', arg=True,optional=False),
                 ScriptOption(['--fnames'], 'File with tomogram names.', arg=True,optional=False),
                 ScriptOption(['--binningTomogram'], 'Binning factor Tomogram from which particles are picked (1: no binning, 3: 3x3 pixels form one superpixel.', arg=True,optional=True),
                 ScriptOption(['--binningSubtomograms'], 'File with tomogram names.', arg=True,optional=True),
                 ScriptOption(['--offset'], 'Offset subtomograms', arg=True,optional=True),
                 ScriptOption(['--size'], 'Size of subtomograms (always cubic).', arg=True,optional=True),
                 ScriptOption(['--weightingType'], 'Weighting type: \n\t\t-1=WB,\n\t\t1=exact,\n\t\t0=none.', arg=True,optional=True),
                 ScriptOption(['--help'],'Help function.',arg=False,optional=True)
        ]

        helper = ScriptHelper(sys.argv[0].split('/')[-1],
                          description='BATCHMODE: Align (and weight) projections taking each marker in markerfile as a reference, and save them in alignment/unweighted_unbinned_marker_XXX. \n',
                          authors='Gijs van der Schot',
                          options = options)

        if len(sys.argv) == 1:
            print(helper)
            sys.exit()
    

        try:
            start, end, fnames,binning_tomogram,binning_subtomograms,offset,size,weighting_type,help = parse_script_options(sys.argv[1:], helper)
        except Exception as e:                                                                                                                                                                                    
            print(sys.version_info)
            sys.exit()

        if help is True:
            print(helper)
            sys.exit()

        start = int(start)
        end = int(end)

        fnames,folders = zip(*[lines.split() for lines in open(fnames).readlines()])

        if binning_tomogram == None:
            binning_tomogram = 1

        if binning_subtomograms == None:
            binning_subtomograms  = 1

        if offset == None:
            offset = '0,0,0'

        if size == None:
            size = 128

        if weighting_type == None:
            weighting_type = 0

        extract_subtomograms(start,end,fnames,folders,binning_tomogram,binning_subtomograms,offset,size,weighting_type)
