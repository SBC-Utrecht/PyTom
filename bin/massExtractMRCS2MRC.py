#!/usr/bin/env pytom

d = '''#!/usr/bin/bash
#SBATCH --time        3:00:00
#SBATCH -N 1
#SBATCH --partition fastq
#SBATCH --ntasks-per-node 1
#SBATCH --job-name    CTF_Batch_ID_{}                                                                       
#SBATCH --output      {}/LogFiles/%j-%x_PhaseFlipImod_{}.out
#SBATCH --oversubscribe

module load imod/4.10.25

cd {}

ctfphaseflip -inp ../sorted/{}_sorted.st -o ctfCorrected.st -an ../sorted/{}_sorted.tlt -defF resultsCTFPlotter.defocus -defT 200 -iW 15 -pi 0.1724 -am 0.08 -cs 2.7 -vo 200 -AxisAngle 270
'''

if __name__=='__main__':
    import sys
    import os
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.basic.structures import ParticleList

    options = [ScriptOption(['-f', '--fileName'], 'Filename of project folders.', True, False),
               ScriptOption(['-h', '--help'], 'Help.', False, True)]

    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name            \ 
                          description='Extract tilt images from mrcstack, and creation of meta data file.',
                          authors='Gijs van der Schot',
                          options=options)

    if len(sys.argv) == 1:
        print(helper)
        sys.exit()

    try:
        filename, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()

    if help is True:
        print(helper)
        sys.exit()


    
    todo = [line.split()[0] for line in open(filename).readlines()]

    n = 0


    outcommands = []


    for folder in todo:
        if os.path.exists(os.path.join(folder, 'logfile.js')):
            a = [[os.path.join(folder, '03_Tomographic_Reconstruction', f, 'ctf/'),f] for f in os.listdir(os.path.join(folder, '03_Tomographic_Reconstruction')) if f.startswith('tomogram') and os.path.exists(os.path.join(folder, '03_Tomographic_Reconstruction', f, 'ctf/resultsCTFPlotter.defocus'))]
        else:
            a = [[os.path.join(folder,'ctf/'), os.path.basename(folder)]]
            folder = os.path.dirname(os.path.dirname(folder))

        for n, (ff, tom) in enumerate(a):

            cmd = d.format(n%20, folder, tom, ff, tom, tom)
            fname = os.path.join(ff, 'ctfCorrected.st')
            outfolder = os.path.join(ff, 'sorted_ctf')
            origdir = os.path.join(os.path.dirname(os.path.dirname(ff)), 'sorted')
            os.system(f'mrcs2mrc.py -f {fname} -t {outfolder} -p sorted_ctf -o {origdir}')

