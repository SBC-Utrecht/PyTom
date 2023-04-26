#!/usr/bin/env python

"""
Created on Sep 02, 2019

@author: dschulte
"""

if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.agnostic.mpi import MPI
    from time import gmtime, strftime
    from pytom.polishing.reconstruct_local_alignment import polish_particles
    from pytom.agnostic.io import read_size
    from pytom.agnostic.mpi import MPI
    from pytom.reconstruction.reconstructionStructures import ProjectionList
    import sys

    helper = ScriptHelper(sys.argv[0].split('/')[-1],  # script name
                          description='Reconstruct a local alignment of particles based on a global alignment.',
                          authors='Douwe Schulte',
                          options=[ScriptOption(['-p', '--particleList'], 'particleList', True, False),
                                   ScriptOption(['--projectionDirectory'], 'Projection Directory', True, False ),
                                   ScriptOption(['-t', '--template'], 'The path to an averaged subtomogram, to use instead of many subtomograms', True, False),
                                   ScriptOption(['--outputDirectory'], 'Results are written to this directory',
                                                True, False),
                                   ScriptOption(['-f', '--FSCPath'],
                                                "The path to an FSC file (.dat) to use as a filter for the cutouts.",
                                                True, True),
                                   ScriptOption(['-b', '--coordinateBinning'],
                                                'Binning factor of the particle list.', True, True),
                                   ScriptOption(['-o', '--recOffset'],
                                                'Cutting offset of the particle list (int x, int y, int z).',
                                                True, True),
                                   ScriptOption(['-m', '--maxParticleShift'],
                                                'Max shift of particles. Shifts larger than max shift will be set to 0',
                                                True, True),
                                   ScriptOption(['-g', '--createGraphics'],
                                                'Flag to turn on to create graphical reports of intermediate steps of '
                                                'the particle polishing', False, True),
                                   ScriptOption(['-n', '--numberOfParticles'],
                                                'To take a subset of the particlelist for debugging purposes',
                                                True, True),
                                   ScriptOption(['--metaFile'], 'Metafile containing tiltangles.', True,True),
                                   ScriptOption(['-v', '--verbose'], 'print intermediate output', False, True)
                                   ])

    try:
        particleList, proj_dir, template, outdir, fsc_path, coordinateBinning, recOffset, max_shift, \
        create_graphics, number_of_particles, metafile, verbose = ll = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        print(helper)
        sys.exit()

    print(ll)
    vol_size = read_size(template, 'x')
    if vol_size % 2 != 0:
        raise ValueError("The particle size has to be an even number.")

    coordinateBinning = int(coordinateBinning) if coordinateBinning else 1
    recOffset = list(map(int, recOffset.split(','))) if recOffset else [0,0,0]
    max_shift = int(max_shift) if max_shift else 6
    create_graphics = True if create_graphics else False
    number_of_particles = int(number_of_particles) if number_of_particles else 0
    peak_border = vol_size // 2 - max_shift
    if peak_border < 0:
        raise ValueError("The given shiftStandardDeviation and peakBorderSize result in a maximal shift "
                         "bigger than the given volume. Please use a bigger volume or a smaller maximal shift.")

    # Split particlelist based on filename
    # Map filename to right projection directory

    mpi = MPI()
    mpi.begin()

    projections = ProjectionList()
    projections.loadDirectory(proj_dir,metafile=metafile)
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
