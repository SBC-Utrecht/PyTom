from multiprocessing import Process
import os
import sys 

global pytompath
pytompath = os.path.dirname(os.popen('dirname `which pytom`').read()[:-1])

def tiltalignment_all_markers(start, end, procs, tiltSeriesName, firstIndex, lastIndex, refIndex, markerFileName, targets, weightingType, tomogramFolder,fnames): 
    '''Aligns tilt images, looping over markers in markerfile as reference. OPTIONAL: procs defines the number of processes running parallel.'''


    tomogram_names = [line.split()[0] for line in open(fnames).readlines()]
    for t in tomogram_names:
        if not os.path.exists(os.path.join(t,'alignment')):
            os.mkdir(os.path.join(t,'alignment'))

    for index in range(start,end):
        cmd = '''cd {}; pytom {}/gui/additional/generateAlignedTiltImages.py \
        --tiltSeriesName {}  \
	--firstIndex {} \
	--lastIndex {} \
	--referenceIndex {} \
	--markerFile {} \
	--projectionTargets {} \
	--projectionBinning 1 \
	--lowpassFilter 0.9 \
	--weightingType {} \
        --numberProcesses {} > alignment/logfile.alignment.txt'''
        cmd = cmd.format(tomogram_names[index], pytompath, tiltSeriesName, firstIndex, 
                         lastIndex, refIndex, markerFileName, targets, weightingType, 
                         procs) 

        p = Process(target=os.system,args=([cmd]))
        p.start()


def extract_subtomograms(start,end,fnames,folders,binning_tomogram=1,binning_subtomograms=1,offset='0,0,0',size=80,weighting_type=0):
    for i in range(start,end):

        cmd = 'reconstructWB.py -p {} --projectionDirectory {} -s {} -b {} -o {} --applyWeighting {}'.format(fnames[i],folders[i],size,binning_tomogram,offset,weighting_type)
        print cmd
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
                 ScriptOption(['--firstIndex'], 'Index of first projection.', arg=True, optional=False),
                 ScriptOption(['--lastIndex'], 'Index of last projection.', arg=True, optional=False),
                 ScriptOption(['--referenceIndex'], 'Index of reference projection used for alignment.', arg=True, optional=False),
                 ScriptOption(['--markerFile'], 'Name of EM markerfile or IMOD wimp File containing marker coordinates.', arg=True, optional=False),
                 ScriptOption(['--projectionTargets'], 'Relative or absolute path to the aligned projections that will be generated + file prefix. default: "align/myTilt"', arg=True, optional=False),
                 ScriptOption(['--weightingType'], 'Type of weighting (-1 default r-weighting, 0 no weighting)', arg=True,optional=False),
                 ScriptOption(['--tomogramFolder'], 'Folder in which tomogram_XXX is located', arg=True,optional=False),
                 ScriptOption(['--fnames'], 'File with tomogram names.', arg=True,optional=False),
                 ScriptOption(['--help'],'Help function.',arg=False,optional=True)
        ]

        helper = ScriptHelper(sys.argv[0].split('/')[-1],
                          description='BATCHMODE: Align (and weight) projections taking each marker in markerfile as a reference, and save them in alignment/unweighted_unbinned_marker_XXX. \n',
                          authors='Gijs van der Schot',
                          options = options)

        if len(sys.argv) == 1:
            print helper
            sys.exit()
        if 1:
            start, end, procs, tiltSeriesName, firstProj, lastProj, referenceIndex, markerFileName, projectionTargets, weightingType,tomogramfolder,fnames,help = parse_script_options(sys.argv[1:], helper)
        else:# Exception as e:
            print sys.version_info
            print e
            sys.exit()


        start = int(start)
        end = int(end)
    
        tiltalignment_all_markers(start,end,procs,tiltSeriesName,firstProj,lastProj,referenceIndex,markerFileName,projectionTargets,weightingType,tomogramfolder,fnames)


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
            print helper
            sys.exit()
    

        try:
            start, end, fnames,binning_tomogram,binning_subtomograms,offset,size,weighting_type,help = parse_script_options(sys.argv[1:], helper)
        except Exception as e:                                                                                                                                                                                    
            print sys.version_info
            print e
            sys.exit()

        if help is True:
            print helper
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
