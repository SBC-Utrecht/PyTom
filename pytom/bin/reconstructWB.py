#!/usr/bin/env pytom
"""
- Created on Jul 19, 2011
- Updated to accomodate updates to the ProjectionList structure on August, 2022

@author: Marten Chaillet, hrabe
"""
import sys, os
from pytom.tools.script_helper import ScriptHelper2, ScriptOption2
from pytom.tools.parse_script_options import parse_script_options2
from pytom.basic.structures import ParticleList, PickPosition, PyTomClassError
from pytom.reconstruction.reconstructionStructures import ProjectionList
from pytom.agnostic.io import write
from multiprocessing import cpu_count


if __name__ == '__main__':
    # This is outdated:
    # this script is significantly linked to
    # pytom/frontend/serverpages/createReconstructionJob.py
    # any changes like parameter names must be changed in the other script, too!

    helper = ScriptHelper2(sys.argv[0].split('/')[-1],
          description='Reconstruct particles in a particle list. Documentation is available at \n\
              http://www.pytom.org/doc/pytom/resonstructWB',
          authors='Marten Chaillet, Thomas Hrabe, FF',
          options=[ScriptOption2(['-t', '--tomogram'],
                                 'Reconstruct a tomogram. Specify name of tomogam here. You do not need a particle '
                                 'list for that!',
                                 'string', 'optional'),
                   ScriptOption2(['-p', '--particleList'],
                                 'XML particle list.',
                                 'file', 'optional'),
                   ScriptOption2(['--projectionList'],
                                 'XML projection list.',
                                 'file', 'optional'),
                   ScriptOption2(['--projectionDirectory'],
                                 'Directory + prefix containing projections.',
                                 'directory', 'optional'),
                   ScriptOption2(['--projectionPrefix'],
                                 'Prefix of projection files in the directory in case they are multiple copies or '
                                 'variants of the same projection.',
                                 'string', 'optional'),
                   ScriptOption2(['-w', '--applyWeighting'],
                                 'If projections are not weighted, apply weighting before. If omited, no weighting.',
                                 'int', 'optional', 0),
                   ScriptOption2(['-s', '--size'],
                                 'Size of particle cube / tomogram. Can provide single integer, or 3 separate by '
                                 'commas, e.g. 464,464,150 ',
                                 'string', 'optional', '100'),
                   ScriptOption2(['-b', '--coordinateBinning'],
                                 'Binning factor of coordinates. If particle coordinates are determined in binned '
                                 'volume (with respect to projections) this binning factor needs to be specified.',
                                 'int', 'optional'),
                   ScriptOption2(['-o', '--recOffset'],
                                 'Cropping offset of the binned tomogram. Important to keep this consistent between '
                                 'tomogram and subtomogram reconstruction.',
                                 'int,int,int', 'optional', [0, 0, 0]),
                   ScriptOption2(['--projBinning'],
                                 'Bin projections BEFORE reconstruction. 1 is no binning, 2 will merge two voxels to '
                                 'one, 3 -> 1, 4 ->1 ...',
                                 'int', 'optional', 1),
                   ScriptOption2(['-m', '--metafile'],
                                 'Supply a metafile to get tiltangles.',
                                 'file', 'optional'),
                   ScriptOption2(['-a', '--alignResultFile'],
                                 'Supply an alignResultFile.',
                                 'file', 'optional'),
                   ScriptOption2(['--scaleFactorParticle'],
                                 'Scale particles by this factor.',
                                 'float', 'optional', 1.),
                   ScriptOption2(['--particlePolishResultFile'],
                                 'Supply a particle polish result file.',
                                 'file', 'optional'),
                   ScriptOption2(['--ctfCorrectionCenter'],
                                 'What is the z-height of the reference center. Option not implemented currently.',
                                 'float', 'optional'),
                   ScriptOption2(['--specimenAngle'],
                                 'Rotation around y axis of the specimen in the tomogram.',
                                 'float', 'optional', .0),
                   ScriptOption2(['--lowpassNy'],
                                 'Fraction of nyquist for low-pass filter that is applied to projections after '
                                 'binning.',
                                 'float', 'optional', 0.9),
                   ScriptOption2(['--tilt-range'],
                                 'Restrict reconstruction to the specified tilt range.',
                                 'int,int', 'optional'),
                   ScriptOption2(['-n', '--numProcesses'],
                                 'Number of processes to split the job over. In case of GPUs this will be overwritten '
                                 'by the number of GPUs.',
                                 'int', 'optional', 1),
                   ScriptOption2(['-g', '--gpuID'],
                                 'Which GPUs do you want to use? This can be a single gpu (0) or multiple (1,'
                                 '3). Multiple GPUs is only implemented for subtomogram reconstruction.',
                                 'string', 'optional')])
    # TODO add alignment origin option. that way imod alignment center can be passed to function.

    tomogram, particle_list_xml, projection_list_xml, projection_directory, projection_prefix, weighting, size, \
    coordinate_binning, rec_offset, projection_binning, metafile, align_result_file, scale_factor_particle, \
    particle_polish_file, ctf_center, specimen_angle, low_pass_ny, tilt_range, nprocs, gpuIDs \
        = parse_script_options2(sys.argv[1:], helper)

    # parse the gpuIDs and overwrite the nprocs if neccessary
    gpuIDs = None if gpuIDs is None else list(map(int, gpuIDs.split(',')))
    nprocs = len(gpuIDs) if not gpuIDs is None else nprocs
    if nprocs > cpu_count():  # reduce cpus if not available
        nprocs = cpu_count()
    size = list(map(int, size.split(',')))

    # check if tomogram output directory is viable
    if tomogram is not None:
        output_dir, tomogram_name = os.path.split(tomogram)
        if output_dir != '' and not os.path.exists(output_dir):
            print('Invalid output location for tomogram, exiting...')
            sys.exit(0)
        _, ext = os.path.splitext(tomogram_name)
        if ext == '':
            tomogram += '.mrc'
        elif ext not in ['.mrc', '.em']:
            print('Invalid output format for tomogram. Needs to be either .em or .mrc. Exiting...')
            sys.exit(0)
        # else everything is fine

    # initialize the projection list structure
    projection_list = ProjectionList()
    if projection_list_xml is not None:
        projection_list.fromXMLFile(projection_list_xml)
        print('projections loaded from xml')
    elif projection_directory is not None:

        try:  # to load from directory
            projection_list.loadDirectory(projection_directory, metafile=metafile,
                                          prefix=('' if projection_prefix is None else projection_prefix))
            print('projections loaded from directory')
        except PyTomClassError:  # Exception raised only if no projections are found in dir
            pass

        if os.path.exists(os.path.join(projection_directory, 'alignmentResults.txt')):
            align_result_file = os.path.join(projection_directory, 'alignmentResults.txt')

    if align_result_file is not None:
        projection_list.load_alignment(align_result_file)
        print('alignment parameters for projections loaded')

    # check if the files are okay and that they are some tilt angles loaded
    if not projection_list.ready_for_reconstruction():
        print('Something went wrong in loading the projection list and the alignment parameters. Please check your '
              'alignResultsFile.txt and the directory where the projections are stored. Exiting...')
        sys.exit(0)

    # start reconstruction of tomogram or subtomograms depending on script options
    if tomogram is not None:

        if len(size) == 1:
            size = size * 3

        vol = projection_list.reconstructVolume(dims=size, reconstructionPosition=rec_offset,
                                                binning=projection_binning, weighting=weighting,
                                                specimen_angle=specimen_angle, low_pass_ny_fraction=low_pass_ny,
                                                scale_factor=scale_factor_particle, tilt_range=tilt_range,
                                                num_procs=nprocs)
        write(tomogram, vol)

    elif particle_list_xml is not None:

        if len(size) > 1:
            print('subtomogram reconstruction uses same cube size for each dimension, selecting x dim as cube size')
        cube_size = size[0]

        # transform the cropping offset
        sx, sy = projection_list[0].getXSize(), projection_list[0].getYSize()

        # TODO is this still needed?
        # from pytom.basic.datatypes import DATATYPE_ALIGNMENT_RESULTS
        # from pytom.agnostic.io import read_size
        # from pytom.gui.guiFunctions import loadstar
        #
        # lar = loadstar(alignResultFile, dtype=DATATYPE_ALIGNMENT_RESULTS)
        # sx,sy,sz = read_size(lar['FileName'][0])
        #
        # if abs((lar['InPlaneRotation'][0] % 180) - 90) < 45:
        #     t = sx
        #     sx = sy
        #     sy = t
        #     print('Inverted shape of images due to rotation axis being close to horizontal, '
        #           'and the reconstruction having the rotation axis vertically.')

        rec_offset[0] = -sx/2 + rec_offset[0] * coordinate_binning
        rec_offset[1] = -sy/2 + rec_offset[1] * coordinate_binning
        rec_offset[2] = -sx/2 + rec_offset[2] * coordinate_binning

        # set particle list in order to reconstruct subtomograms
        particle_list = ParticleList()
        try:
            particle_list.fromXMLFile(particle_list_xml)
        except RuntimeError:
            print('Error reading particleList XML file! Exiting...')
            sys.exit()

        for particle in particle_list:
            pickPosition = particle.getPickPosition()
            x = (pickPosition.getX() * coordinate_binning + rec_offset[0])
            y = (pickPosition.getY() * coordinate_binning + rec_offset[1])
            z = (pickPosition.getZ() * coordinate_binning + rec_offset[2])

            # if abs(scaleFactorParticle-1) > 0.000001:
            #     x = x * float(scaleFactorParticle)
            #     y = y * float(scaleFactorParticle)
            #     z = z # * float(scaleFactorParticle)

            particle.setPickPosition(PickPosition(x=x, y=y, z=z))

        projection_list.reconstructVolumes(particles=particle_list, cube_size=cube_size, binning=projection_binning,
                                           weighting=weighting, low_pass_ny_fraction=low_pass_ny,
                                           post_scale=scale_factor_particle, num_procs=nprocs, ctfcenter=ctf_center,
                                           polishResultFile=particle_polish_file, tilt_range=tilt_range,
                                           show_progress_bar=True, gpuIDs=gpuIDs)

    else:
        print('No valid tomogram or particle list provided. Exiting...')
        sys.exit(0)
