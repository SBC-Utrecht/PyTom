import os
import sys

from pytom.gui.additional.multi_tilt_alignment import extract_subtomograms
from pytom.gui.sbatch_commands import cmd_batch_subtomo_extraction
from pytom.gui.guiFunctions import write_text2file

global pytompath
pytompath = os.path.dirname(os.popen('dirname `which pytom`').read()[:-1])


def submit_slurm(projectfolder,fnames,binning_tomogram,binning_subtomograms,offset,size,weighting_type,deploy):

    fd = open(os.path.join(projectfolder,fnames))
    f = [line.split() for line in fd.readlines()]
    fd.close()
    number_processes = 20

    for i in range(0,len(f),number_processes):
        start =  i
        end = i+number_processes

        cmd = cmd_batch_subtomo_extraction.format( projectfolder, pytompath, start, end, fnames, weighting_type, binning_tomogram, binning_subtomograms, offset, size)

        ff = '{}/jobscripts/subtomo_extraction_{:03d}.job'.format(projectfolder, i)
        write_text2file(cmd,ff, 'w' )

        if deploy:
            os.system('sbatch {}'.format(ff))
        

if __name__=='__main__':
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.gui.additional.reconstructionFunctions import alignWeightReconstruct
    from pytom_volume import read
    from pytom_numpy import vol2npy
    import os
    from multiprocessing import Process
    import time
    import sys 

    options=[ScriptOption(['--projectFolder'], 'Working directory. Other paths are relative to this folder.', arg=True,optional=False),
             ScriptOption(['--fnames'], 'File with tomogram names.', arg=True,optional=False),
             ScriptOption(['--binningTomogram'], 'Binning factor Tomogram from which particles are picked (1: no binning, 3: 3x3 pixels form one superpixel.', arg=True,optional=True),
             ScriptOption(['--binningSubtomograms'], 'File with tomogram names.', arg=True,optional=True),
             ScriptOption(['--offset'], 'Offset subtomograms', arg=True,optional=True),
             ScriptOption(['--size'], 'Size of subtomograms (always cubic).', arg=True,optional=True),
             ScriptOption(['--weightingType'], 'Weighting type: \n\t\t-1=WB,\n\t\t1=exact,\n\t\t0=none.', arg=True,optional=True),
             ScriptOption(['--deploy'],'Submit jobs to slurm. If not supplied jobscripts are generated but not submitted.',arg=False,optional=True),
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
        project_folder, fnames, binning_tomogram, binning_subtomograms, offset, size, weighting_type, deploy, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:                                                                                                                                                                                    
        print sys.version_info
        print e
        sys.exit()

    if help is True:
        print helper
        sys.exit()

    if project_folder == None: project_folder = '.'
    if fnames == None: sys.exit()
    if binning_tomogram ==None: binning_tomogram = 1
    if binning_subtomograms ==None: binning_subtomograms = 1
    if offset == None: offset='0,0,0'
    if size == None: size = 128
    if weighting_type == None: weighting_type = 0
                                                                    
    submit_slurm(project_folder,fnames,binning_tomogram,binning_subtomograms,offset,size,weighting_type,deploy)
