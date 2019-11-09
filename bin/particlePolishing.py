#!/usr/bin/env pytom

"""
Created on Sep 02, 2019

@author: dschulte
"""

if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.tompy.mpi import MPI
    from time import gmtime, strftime
    from pytom.reconstruction.reconstruct_local_alignment import polish_particles
    from pytom.tompy.io import read_size
    from pytom.tompy.mpi import MPI
    from pytom.reconstruction.reconstructionStructures import ProjectionList

    helper = ScriptHelper(sys.argv[0].split('/')[-1],  # script name
                          description='Reconstruct a local alignment of particles based on a global alignment.',
                          authors='Douwe Schulte',
                          options=[ScriptOption(['-p', '--particleList'],
                                                'particleList', 'has arguments', 'required'),
                                   ScriptOption(['--projectionDirectory'], '', True, False ),
                                   ScriptOption(['-t', '--template'],
                                                'The path to an averaged subtomogram, to use instead of many '
                                                'subtomograms', 'string', 'required'),
                                   ScriptOption(['--outputDirectory'], 'Results are written to this directory',
                                                'string', 'optional', './'),
                                   ScriptOption(['-f', '--FSCPath'],
                                                "The path to an FSC file (.dat) to use as a filter for the cutouts.",
                                                'has arguments', 'optional',''),
                                   ScriptOption(['-b', '--coordinateBinning'],
                                                'Binning factor of the particle list.', 'uint', 'required'),
                                   ScriptOption(['-o', '--recOffset'],
                                                'Cutting offset of the particle list (int x, int y, int z).',
                                                'int,int,int', 'required'),
                                   ScriptOption(['-m', '--maxParticleShift'],
                                                'Max shift of particles. Shifts larger than max shift will be set to 0',
                                                'uint', 'optional', 5),
                                   ScriptOption(['-g', '--createGraphics'],
                                                'Flag to turn on to create graphical reports of intermediate steps of '
                                                'the particle polishing', 'no arguments', 'optional', False),
                                   ScriptOption(['-n', '--numberOfParticles'],
                                                'To take a subset of the particlelist for debugging purposes',
                                                'uint', 'optional', -1),

                                   ScriptOption(['-v', '--verbose'], 'print intermediate output', 'no_arguments',
                                                'optional', False)
                                   ])

    particleList, proj_dir, template, fsc_path, outdir, coordinateBinning, recOffset, max_shift, \
    create_graphics, number_of_particles, verbose = parse_script_options(sys.argv[1:], helper)

    vol_size = read_size(template, 'x')

    if vol_size % 2 != 0:
        raise ValueError("The particle size has to be an even number.")

    peak_border = vol_size // 2 - max_shift
    if peak_border < 0:
        raise ValueError("The given shiftStandardDeviation and peakBorderSize result in a maximal shift "
                         "bigger than the given volume. Please use a bigger volume or a smaller maximal shift.")

    # Split particlelist based on filename
    # Map filename to right projection directory

    mpi = MPI()
    mpi.begin()

    projections = ProjectionList()
    projections.loadDirectory(proj_dir)
    projections.sort()

    # create the list of projectionfilenames and tiltangles
    proj = []
    tilt_angles = []

    for p in projections:
        proj.append(p.getFilename())
        tilt_angles.append(p.getTiltAngle())

    polish_particles(particleList, proj_dir, template, coordinateBinning, recOffset, proj, tilt_angles, mpi,
                     outputDirectory=outdir, peak_border=peak_border, fsc_path=fsc_path,
                     create_graphics=create_graphics, number_of_particles=number_of_particles, verbose=verbose)

    mpi.end()
