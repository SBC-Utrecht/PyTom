import os
from pytom.gui.sbatch_commands import cmd_batch_alignment
from pytom.gui.guiFunctions import write_text2file

global pytompath
pytompath = os.path.dirname(os.popen('dirname `which pytom`').read()[:-1])


def batch_tilt_alignment( number_tomonames, fnames_tomograms='', projectfolder='.', num_procs=20, num_procs_per_proc=1, tiltseriesname='sorted/sorted_aligned',
                         markerfile='sorted/markerfile.em',targets='alignment', firstindex=1, lastindex=21, refindex=11, weightingtype=0, deploy=False):
    '''BATCHMODE: tilt alignment. Submits a number of sbatch jobs to slurm queueing system. Each job calculates the tilt aligment for each marker in a markerfile.  It divides the number or jobs with respect to the num_procs.'''

    for n in range(number_tomonames):

        if not n % num_procs == 0: 
            continue

        cmd = cmd_batch_alignment.format( projectfolder, pytompath, n, min(number_tomonames, num_procs+n),
                                          num_procs_per_proc, tiltseriesname, markerfile, targets, projectfolder,
                                          firstindex, lastindex, refindex, weightingtype, fnames_tomograms)

        write_text2file(cmd,'{}/jobscripts/alignment_{:03d}.job'.format(projectfolder, n), 'w' )

        if deploy:
            os.system('sbatch {}/jobscripts/alignment_{:03d}.job'.format(projectfolder, n))


if __name__=="__main__":
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.gui.additional.reconstructionFunctions import alignWeightReconstruct
    from pytom_volume import read
    from pytom_numpy import vol2npy
    import os
    from multiprocessing import Process
    import time
    import sys 

    options=[ScriptOption(['--fileTomoFolders'],'File with tomogram folder names (Relative to projectfolder)', arg=True,optional=False),
             ScriptOption(['--projectFolder'],'base folder where all tomograms are located in. From this folder jobs will be submitted.', arg=True,optional=False),
             ScriptOption(['--numberOfTomograms'],'Number of tomograms.', arg=True,optional=False),
             ScriptOption(['--numberProcesses'],'Number of occupied CPU cores per node.', arg=True,optional=True),
             ScriptOption(['--numberProcsPerTomogram'],'Number of processes per tomogram. Useful to set > 1 if number of markers is >10, and number of tomograms is small.', arg=True,optional=True),
             ScriptOption(['--tiltSeriesName'],'Path to tilt images including prefix (relative to projectfolder)', arg=True,optional=True),
             ScriptOption(['--firstIndex'],'First index.', arg=True,optional=True),
             ScriptOption(['--lastIndex'],'Last index.', arg=True,optional=True),
             ScriptOption(['--referenceIndex'],'Reference index.', arg=True,optional=True),
             ScriptOption(['--markerFileName'],'Markerfile name, relative to projectfolder.', arg=True,optional=True),
             ScriptOption(['--projectionTargets'],'Path to output folder. In this folder unbinned)unweighted_marker_XXX is created. (Relative to projectfolder)', arg=True,optional=True),
             ScriptOption(['--weightingType'],'Last index.', arg=True,optional=True),
             ScriptOption(['--deploy'],'submit jobs to slurm squeue.', arg=False,optional=True),
             ScriptOption(['-h', '--help'], 'Help.', False, True)    ]


#num_procs=20, num_procs_per_proc=1, tiltseriesname='sorted/sorted_aligned',
#markerfile='sorted/markerfile.em',targets='alignment', firstindex=1, lastindex=21, refindex=11, weightingtype=0, deploy=False
    helper = ScriptHelper(sys.argv[0].split('/')[-1],
                          description='BATCHSUBMISSION TILT ALIGNMENT: submits  \n',
                          authors='Gijs van der Schot',
                          options = options)

    if len(sys.argv) == 1:
        print helper
        sys.exit()
    
    fnames_tomograms, projectfolder, numTomograms, num_procs, procPerTomo, tiltSeriesName, firstIndex, lastIndex, referenceIndex, markerFileName, projectionTargets, weightingType,deploy,help = parse_script_options(sys.argv[1:], helper)

    if not projectfolder: projectfolder = '/home/gijsvds/ost_2/03_Tomographic_Reconstruction'
    if not fnames_tomograms:
        tomonames = sorted([t  for t in os.listdir(projectfolder) if len(t.split('_')[-1])==3])
        fnames_tomograms = '{}/names_tomograms.txt'.format(projectfolder)
        out = open(fnames_tomograms,'w')
        for i in tomonames:
            out.write(i+'\n')
        out.close()
        numTomograms=len(tomonames)
    if not num_procs: num_procs=20
    if not procPerTomo: procPerTomo=1
    if not tiltSeriesName: tiltSeriesName= 'sorted/sorted'
    if not markerFileName: markerFileName='sorted/markerfile.em'
    if not projectionTargets: projectionTargets='alignment'
    if not weightingType: weightingType = 0
    if not firstIndex: firstIndex=1
    if not lastIndex: lastIndex=21
    if not referenceIndex: referenceIndex=11
#( number_tomonames, fnames_tomograms='', projectfolder='.', num_procs=20, num_procs_per_proc=1, tiltseriesname='sorted/sorted_aligned',
#                         markerfile='sorted/markerfile.em',targets='alignment', firstindex=1, lastindex=21, refindex=11, weightingtype=0, deploy=False)
    print  parse_script_options(sys.argv[1:], helper)
    
    if help is True:
        print helper
        sys.exit()
        

    batch_tilt_alignment(int(numTomograms), projectfolder=projectfolder, fnames_tomograms=fnames_tomograms, num_procs=num_procs, num_procs_per_proc=procPerTomo,
                         tiltseriesname=tiltSeriesName,markerfile=markerFileName,targets=projectionTargets,weightingtype=weightingType,
                         firstindex=firstIndex,lastindex=lastIndex,refindex=referenceIndex,deploy=deploy)
