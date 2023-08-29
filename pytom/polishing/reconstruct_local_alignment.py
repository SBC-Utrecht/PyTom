#!/usr/bin/env pytom
# -*- coding: utf-8 -*-

"""
Created on Sep 02, 2019
Updated on Sep 08, 2019 (GS)
@author: dschulte
"""
import matplotlib

try:
    matplotlib.use('Qt5Agg')
except:
    pass
import pylab as pp
from pytom.gpu.initialize import xp, device
from pylab import imshow, show, savefig, subplots

def polish_particles(particle_list_filename, projection_directory, averaged_subtomogram, binning, offset,
                     projections, tilt_angles, mpi, fsc_path='', peak_border=75, outputDirectory='./',
                     create_graphics=False, number_of_particles=0, verbose=False ):
    """
    To polish a particle list based on (an) initial subtomogram(s).

    :param particle_list_filename: the filename of the particlelist
    :type particle_list_filename: str
    :param projection_directory: the directory of the projections
    :type projection_directory: str
    :param averaged_subtomogram: to give a path to an averaged subtomogram to be used instead of subtomograms of all
               particles separately
    :type averaged_subtomogram: str
    :param binning: the binning factor used
    :type binning: int
    :param offset: the offset used (x, y, z)
    :type offset: list(int, int, int)

    :param projections: a list with filenames of projections
    :type projections: list(str)
    :param tilt_angles: the list of tiltangles used
    :type tilt_angles: list(int)
    :param mpi: the instance of mpi which can be used for multi threading
    :type mpi: mpi instance
    :param create_graphics: to create plots of major parts of the algorithm, mainly used for debugging
               and initial creation
    :type create_graphics: bool
    :param number_of_particles: to use a subset of the particles for the particle polishing
    :type number_of_particles: int
    :param skip_alignment: skips the alignment phase, does not do particle polishing
    :type skip_alignment: bool

    :return: nothing, it writes everything to disk
    :returntype: void
    """
    assert number_of_particles == -1 or number_of_particles > 0
    assert binning > 0
    assert vol_size > 0
    assert vol_size % 2 == 0
    assert isinstance(projections, list)
    assert isinstance(vol_size, int)
    assert isinstance(binning, int)
    assert isinstance(offset, list) and len(offset) == 3
    assert isinstance(offset[0], int) and isinstance(offset[1], int) and isinstance(offset[2], int)
    assert isinstance(tilt_angles, list)
    assert isinstance(particle_list_filename, str)
    assert isinstance(projection_directory, str)
    assert isinstance(create_graphics, bool)
    assert isinstance(averaged_subtomogram, str)
    assert isinstance(number_of_particles, int)
    assert isinstance(skip_alignment, bool)

    import numpy as np
    import os
    from pytom.agnostic.io import read_size, read
    from pytom.basic.datatypes import fmtLAR, headerLocalAlignmentResults, LOCAL_ALIGNMENT_RESULTS

    # load particle list
    from pytom.basic.structures import ParticleList

    particlelist = ParticleList()
    particlelist.fromXMLFile(particle_list_filename)
    particle_list_name = os.path.splitext(os.path.basename(str(particle_list_filename)))[0]
    if number_of_particles > 0:
        particlelist = particlelist[:number_of_particles]

    if verbose:
        print(len(particlelist))
        print("{:s}> Creating the input array".format(gettime()))

    dimz = read_size(particlelist[0].getPickPosition().getOriginFilename(), 'z') * binning
    vol_size = read_size(particlelist[0].getFilename(), 'x')
    input_to_processes = []



    # data = {}
    # for projectioname in projections:
    #     data[projectioname] = read(projectioname)
    #
    # template = read(averaged_subtomogram)
    input_to_processes = []
    for particle_number, particle in enumerate(particlelist):

        rot = (particle.getRotation().getZ1(), particle.getRotation().getX(), particle.getRotation().getZ2())

        # loop over tiltrange, take patch and cross correlate with reprojected subtomogram
        for img, ang in zip(projections, tilt_angles):
            input_to_processes.append([averaged_subtomogram, ang, offset, vol_size,
                                       particle.getPickPosition().toVector(), rot, particle.getFilename(),
                                       particle_number, binning, img, create_graphics, fsc_path, dimz, peak_border,
                                       ])

    if verbose:
        print(len(input_to_processes))
        print("{:s}> Created the input array".format(gettime()))


    if verbose: print("{:s}> Started on running the process".format(gettime()))

    lists = list(zip(*input_to_processes))
    if verbose: print(len(list(lists)))

    output = mpi.parfor(run_single_tilt_angle, list(zip(lists[0], lists[1],lists[2],lists[3],lists[4],lists[5],lists[6],lists[7],lists[8],lists[9],lists[10],lists[11],lists[12],lists[13]))) # Some problem internally 23/10/2019

    results_file = os.path.join(outputDirectory, f"resultsPolish_{particle_list_name}.txt")
    np.savetxt(results_file, np.array(output, dtype=LOCAL_ALIGNMENT_RESULTS), fmt=fmtLAR, header=headerLocalAlignmentResults)

    if verbose: print("{:s}> Ran the processes".format(gettime()))



def run_single_tilt_angle_unpack(inp):
    """
    unpack a list with arguments to "run_single_tilt_angle"
    @param inp: the arguments to "run_single_tilt_angle" in the same order in a single list (iterable)
    @type inp: list/something with indexing
    @return: the value from "run_single_tilt_angle"
    @returntype: list
    """
    return run_single_tilt_angle(*inp)

def cut_patch(projection, ang, pick_position, vol_size=200, binning=8, dimz=0, offset=[0,0,0]):
    from pytom.agnostic.transform import rotate3d, rotate_axis
    import numpy as np
    from math import cos, sin, pi
    from pytom.agnostic.transform import cut_from_projection
    from pytom.agnostic.io import read

    # Get the size of the original projection
    dim_x = projection.shape[0]
    dim_z = dim_x if dimz is None else dimz

    x, y, z = pick_position
    x = (x + offset[0]) * binning
    y = (y + offset[1]) * binning
    z = (z + offset[2]) * binning

    # Get coordinates of the paricle adjusted for the tilt angle
    yy = y  # assume the rotation axis is around y
    xx = (cos(ang * pi / 180) * (x - dim_x / 2) - sin(ang * pi / 180) * (z - dim_z / 2)) + dim_x / 2

    # Cut the small patch out
    patch = cut_from_projection(projection.squeeze(), [xx, yy], [vol_size, vol_size])
    patch = patch #- patch.mean()

    return patch, xx, yy


def normalize_image(image, mask, p):
    from pytom.agnostic.normalise import meanUnderMask, stdUnderMask
    mp = meanUnderMask(image, mask, p)
    normalized_image = (image - mp) / stdUnderMask(image, mask, mp, p)
    return normalized_image * mask

def run_single_tilt_angle(subtomogram, ang, offset, vol_size, particle_position, particle_rotation,  particle_filename,
                          particle_number, binning, img, create_graphics, fsc_path, dimz, peak_border):
    """
    To run a single tilt angle to allow for parallel computing

    @param ang: the tilt angle
    @type ang: int
    @param subtomogram: the filename of the subtomogram
    @type subtomogram: str
    @param offset: the offset used (x,y,z)
    @type offset: list(int, int, int)
    @param vol_size: the size of the volume to be reconstructed (in pixels)
    @type vol_size: int
    @param particle_position: the position of the particle in vector format,
               as given by particle.pickPosition().toVector()
    @type particle_position: tuple
    @param particle_rotation: the rotation of the particle (Z1/phi, X/the, Z2/psi)
    @type particle_rotation: tuple
    @param particle_filename: the filename of the particle, as given by particle.getfilename()
    @type particle_filename: str
    @param particle_number: the number of the particle, to allow for unique mapping
    @type particle_number: int
    @param binning: the binning factor used
    @type binning: int
    @param img: the filename of the projection to be used
    @type img: str
    @param create_graphics: to flag if images should be created for human inspection of the work done
    @type create_graphics: bool
    @return: the newly found positions of the particle, as a list  in the LOCAL_ALIGNMENT_RESULTS format
    @returntype: list
    """


    from pytom.reconstruction.reconstructionFunctions import alignImageUsingAlignmentResultFile
    from pytom.agnostic.transform import rotate3d, rotate_axis
    from pytom.gui.guiFunctions import datatypeAR, loadstar
    import numpy as np
    from math import cos, sin, pi, sqrt
    from pytom.agnostic.tools import create_circle
    from pytom.agnostic.transform import cut_from_projection
    from pytom.agnostic.filter import applyFourierFilterFull, bandpass_circle
    from pytom.agnostic.io import read, write
    import os
    import time

    t = time.time()
    # print(particle_filename, ang)

    # Filter using FSC
    fsc_mask = None
    k = 0.01
    subtomogram = read(subtomogram) * read(
        '/data/gijsvds/ctem/05_Subtomogram_Analysis/Alignment/GLocal/mask_200_75_5.mrc')

    import os
    from pytom.agnostic.filter import filter_volume_by_profile, profile2FourierVol
    fsc_path = ''
    if os.path.isfile(fsc_path):
        profile = [line.split()[0] for line in open(fsc_path, 'r').readlines()]
        fsc_mask3d = profile2FourierVol(profile, subtomogram.shape)

        subtomogram = applyFourierFilterFull(subtomogram, fsc_mask3d)


    # Cross correlate the templat
    # e and patch, this should give the pixel shift it is after
    from pytom.agnostic.correlation import norm_xcf

    # Get template
    # First rotate the template towards orientation of the particle, then to the tilt angle



    img = read(img)

    rotated1 = rotate3d(subtomogram, phi=particle_rotation[0], the=particle_rotation[1], psi=particle_rotation[2])
    rotated2 = rotate_axis(rotated1, -ang, 'y')  # SWITCHED TO ROTATE AXIS AND ANGLE *-1 THIS IS AN ATTEMPT
    template = rotated2.sum(axis=2)

    # write('pp1_template.mrc', template)


    # img = read(img)
    try:
        patch, xx, yy = cut_patch(img, ang, particle_position, dimz=dimz, vol_size=vol_size, binning=binning )
        mask2d = create_circle(patch.shape, radius=75, sigma=5, num_sigma=2)
        patch *= mask2d
        # write('pp1_patch.mrc', patch)

        if os.path.isfile(fsc_path):
            profile = [line.split()[0] for line in open(fsc_path, 'r').readlines()]
            fsc_mask2d = profile2FourierVol(profile, patch.shape)
            patch = applyFourierFilterFull(patch, fsc_mask2d)


        template = normalize_image(template, mask2d, mask2d.sum())
        patch = normalize_image(patch, mask2d, mask2d.sum())

        if 1:
            ff = xp.ones_like(patch)
        else:
            ff = bandpass_circle(patch.squeeze(), 6, 25, 3)

        ccf = normalised_cross_correlation(template, patch, ff)

        points2d1 = find_sub_pixel_max_value_2d(ccf.copy(), ignore_border=peak_border)
        points2d2 = find_sub_pixel_max_value(ccf.copy(), ignore_border=peak_border)
        points2d = list(points2d2[:2]) + [points2d1[-1]]

        x_diff = points2d[0] - vol_size / 2
        y_diff = points2d[1] - vol_size / 2

        dist = sqrt(x_diff**2+y_diff**2)

        #rx, ry, ii, oox, ooy, x2, y2 = find_sub_pixel_max_value_voltools(ccf.copy())

        #print(rx,ry, x_diff, y_diff, x2, y2, oox, ooy, ii.shape)
        #print(f'{particle_number:3d} {ang:5.1f}, {dist:5.2f} {x_diff} {y_diff} {ccf.max()}')
        if create_graphics:
            rx, ry, ii, oox, ooy, x2, y2 = find_sub_pixel_max_value_voltools(ccf.copy(),k=k)
            print(rx, ry)
            from scipy.ndimage.filters import gaussian_filter
            # Create an image to display and/or test the inner workings of the algorithm

            points = find_sub_pixel_max_value(ccf, ignore_border=peak_border)

            nx, ny, nz = particle_position
            nx += x_diff
            ny += y_diff

            npatch = cut_from_projection(img.squeeze(), [xx + x_diff, yy + y_diff], [vol_size,
                                                                                     vol_size])  # img.sum(axis=2)[int(xx+x_diff-v):int(xx+x_diff+v), int(yy+y_diff-v):int(yy+y_diff+v)]  #
            npatch = normalize_image(npatch, mask2d, mask2d.sum())

            nccf = normalised_cross_correlation(template, npatch.squeeze(), ff)
            npoints = find_sub_pixel_max_value(nccf, ignore_border=peak_border)
            npoints2d = find_sub_pixel_max_value_2d(nccf, ignore_border=peak_border)
            #rx, ry = abs(xp.array(npoints2d[:2]) - 100)

            from scipy.ndimage.filters import gaussian_filter
            # Create an image to display and/or test the inner workings of the algorithm

            points = find_sub_pixel_max_value(ccf, ignore_border=peak_border)

            nx, ny, nz = particle_position
            nx += x_diff
            ny += y_diff

            npatch = cut_from_projection(img.squeeze(), [xx + x_diff, yy + y_diff], [vol_size,
                                                                                     vol_size])  # img.sum(axis=2)[int(xx+x_diff-v):int(xx+x_diff+v), int(yy+y_diff-v):int(yy+y_diff+v)]  #
            npatch = normalize_image(npatch, mask2d, mask2d.sum())

            nccf = normalised_cross_correlation(template, npatch.squeeze(), ff)
            npoints = find_sub_pixel_max_value(nccf, ignore_border=peak_border)
            npoints2d = find_sub_pixel_max_value_2d(nccf, ignore_border=peak_border)
            #rx, ry = abs(xp.array(npoints2d[:2]) - 100)

            m_style = dict(color='tab:blue', linestyle=':', marker='o', markersize=5, markerfacecoloralt='tab:red')
            m_style_alt = dict(color='tab:red', linestyle=':', marker='o', markersize=5, markerfacecoloralt='tab:orange')

            grid = pp.GridSpec(3, 3, wspace=0, hspace=0.35, left=0.05, right=0.95, top=0.90, bottom=0.05)

            ax_0_0 = pp.subplot(grid[0, 0])
            ax_0_1 = pp.subplot(grid[0, 1])
            ax_0_2 = pp.subplot(grid[0, 2])
            ax_1_0 = pp.subplot(grid[1, 0])
            ax_1_1 = pp.subplot(grid[1, 1])
            ax_1_2 = pp.subplot(grid[1, 2])
            ax_2_0 = pp.subplot(grid[2, 0])
            ax_2_1 = pp.subplot(grid[2, 1])
            ax_2_2 = pp.subplot(grid[2, 2])

            ax_0_0.axis('off')
            ax_0_1.axis('off')
            ax_0_2.axis('off')
            ax_1_0.axis('off')
            ax_1_1.axis('off')
            ax_1_2.axis('off')
            ax_2_0.axis('off')
            ax_2_1.axis('off')
            ax_2_2.axis('off')

            axis_title(ax_0_0, "Cutout")
            ax_0_0.imshow(applyFourierFilterFull(patch,xp.fft.fftshift(ff)))
            axis_title(ax_0_1, "Template")
            ax_0_1.imshow(template)
            axis_title(ax_0_2, "Shifted Cutout\n(based on cross correlation)")
            ax_0_2.imshow(npatch.squeeze())

            axis_title(ax_1_0, u"Cross correlation\ncutout × template")
            ax_1_0.imshow(ccf)
            ax_1_0.plot([points[1]], [points[0]], fillstyle='none', **m_style)
            ax_1_0.plot([points2d[1]], [points2d[0]], fillstyle='none', **m_style_alt)
            ax_1_0.plot([vol_size / 2], [vol_size / 2], ",k")

            ax_1_1.text(0.5, 0.8,
                        "Red: 2D spline interpolation\nx: {:f}\ny: {:f}\nBlue: 1D spline interpolation\nx: {:f}\ny: {:f}"
                        "\nBlack: center".format(
                            x_diff, y_diff, points[0] - vol_size / 2, points[1] - vol_size / 2), fontsize=8,
                        horizontalalignment='center', verticalalignment='center', transform=ax_1_1.transAxes)

            axis_title(ax_1_2, u"Cross correlation\nshifted cutout × template")
            ax_1_2.imshow(nccf)
            ax_1_2.plot([npoints[0]], [npoints[1]], fillstyle='none', **m_style)
            ax_1_2.plot([npoints2d[0]], [npoints2d[1]], fillstyle='none', **m_style_alt)
            ax_1_2.plot([vol_size / 2], [vol_size / 2], ",k")

            axis_title(ax_2_0, u"Zoom into red peak\nin CC cutout × template")
            d = 10

            points2d = list(xp.unravel_index(ccf.argmax(), ccf.shape)) + [points2d[-1]]
            peak = ccf[int(points2d[0]) - d:int(points2d[0]) + d, int(points2d[1]) - d:int(points2d[1] + d)]
            ax_2_0.imshow(peak)

            axis_title(ax_2_1, u"Zoom into red peak\nin CC cutout × template\ninterpolated")
            ax_2_1.imshow(ii)

            ax_2_1.plot( [ y2 + (y_diff - ry)/k],[x2 + (x_diff - rx)/k], fillstyle='none', **m_style)
            ax_2_1.plot([y2], [x2], fillstyle='none', **m_style_alt)
            axis_title(ax_2_2, u"Cutout\nGaussian filter σ3")
            import scipy
            ax_2_2.imshow(scipy.ndimage.gaussian_filter(patch, 3))

            pp.savefig("../Images/test/polish_particle_{:04d}_tiltimage_{:05.2f}_shift_{:.1f}.png".format(particle_number, ang, dist))
        del img
    except Exception as e:
        print(e)
        x_diff = y_diff = 0

    print(f'{particle_number:3d} {ang:4d} {x_diff:5.2f} {y_diff:5.2f} {points2d1[0]-patch.shape[0]//2:.2f} {points2d1[1]-patch.shape[0]//2:.2f} {time.time()-t:5.3f}')
    return particle_number, x_diff, y_diff, ang, 0, 0, particle_filename


def axis_title(axis, title):
    """Small helper function to simplify the creation of text on a pylab plot"""
    axis.text(0.5, 1.05, title, fontsize=8, horizontalalignment='center', transform=axis.transAxes)


def run_polished_subtomograms(particle_list_filename, projection_directory, particle_polish_file, binning, offset,
                              vol_size, start_glocal, glocal_jobname, glocal_nodes, glocal_particle_list, dimz):
    """
    Reconstructs subtomograms based on a polished particlelist, writes these to the places as specified in particlelist

    @param particle_list_filename: The name of the file of the particlelist
    @type particle_list_filename: str
    @param projection_directory: The name of the directory containing the projections
    @type projection_directory: str
    @param particle_polish_file: The name of the file containing the polished alignment results
    @type particle_polish_file: str
    @param binning: The binning factor
    @type binning: int
    @param offset: The reconstruction offset
    @type offset: list(int, int, int)
    @param vol_size: The size of the particle
    @type vol_size: int
    @return: void
    @returntype: void
    """
    assert isinstance(particle_list_filename, str)
    assert isinstance(projection_directory, str)
    assert isinstance(particle_polish_file, str)
    assert isinstance(binning, int)
    assert isinstance(offset, list) and len(offset) == 3
    assert isinstance(offset[0], int) and isinstance(offset[1], int) and isinstance(offset[2], int)
    assert isinstance(vol_size, int)

    import os.path
    assert os.path.isfile(particle_list_filename)
    assert os.path.isfile(particle_polish_file)
    assert os.path.isdir(projection_directory)

    cwd = os.getcwd()

    first_batchfile = """#!/usr/bin/bash
#SBATCH --time        1:00:00
#SBATCH -N 1
#SBATCH --partition defq
#SBATCH --ntasks-per-node 20
#SBATCH --job-name    polishedReconstruction                                                                     
#SBATCH --error="../LogFiles/%j-polished_subtomograms.err"
#SBATCH --output="../LogFiles/%j-polished_subtomograms.out"

module load openmpi/2.1.1 python/2.7 lib64/append pytom/dev/dschulte

cd {:s}

reconstructWB.py --particleList {:s} \
--projectionDirectory {:s} \
--coordinateBinning {:d} \
--size {:d} \
--applyWeighting \
--projBinning 1 \
--recOffset {:d},{:d},{:d} \
--particlePolishFile {:s} \
{:s}-n 20"""\
        .format(cwd, particle_list_filename, projection_directory, binning, vol_size, offset[0], offset[1], offset[2],
                particle_polish_file, ("" if dimz is None else ("--dimZ " + str(dimz) + " ")))

    f = open("polished_subtomograms.sh", "w+")
    f.write(first_batchfile)
    f.close()

    import subprocess
    out = subprocess.check_output(['sbatch', 'polished_subtomograms.sh'])
    print("Started reconstruction")

    if start_glocal:
        if out.startswith("Submitted batch job "):
            pid = out.split(" ")[3]

            mask_filename = "/data2/dschulte/BachelorThesis/Data/VPP2/05_Subtomogram_Analysis/Alignment/FRM/test-11-09-19/FRM_mask_200_70_4.mrc"
            iterations = 6
            pixelsize = 2.62
            particleDiameter = 300

            glocal_batchfile = """#!/usr/bin/bash
#SBATCH --time        20:00:00
#SBATCH -N {:d}
#SBATCH --partition defq
#SBATCH --ntasks-per-node 20
#SBATCH --job-name    pGLocalAlign                                                                       
#SBATCH --error="../LogFiles/%j-GLocal_polished_subtomograms.err"
#SBATCH --output="../LogFiles/%j-GLocal_polished_subtomograms.out"
#SBATCH --dependency=afterok:{:s}

module load openmpi/2.1.1 python/2.7 lib64/append pytom/dev/dschulte

cd {:s}

sleep 5 & mpiexec -n {:d} pytom /data2/dschulte/pytom-develop/pytom/bin/GLocalJob.py \
    -p {:s} \
    --mask {:s} \
    --SphericalMask \
    --destination Alignment/GLocal/{:s}/ \
    --numberIterations {:d} \
    --pixelSize {:f} \
    --particleDiameter {:d} \
    --jobName Alignment/GLocal/{:s}/job.xml \
    --noShift""".format(glocal_nodes, pid, cwd, glocal_nodes*20, glocal_particle_list, mask_filename, glocal_jobname, iterations, pixelsize, particleDiameter, glocal_jobname)
            f = open("glocal_align.sh", "w+")
            f.write(glocal_batchfile)
            f.close()

            if not os.path.isdir(cwd + "/Alignment/GLocal/" + glocal_jobname): os.mkdir(cwd + "/Alignment/GLocal/" + glocal_jobname)

            out = subprocess.check_output(['sbatch', 'glocal_align.sh'])

            if out.startswith("Submitted batch job "):
                print("Reconstruction and Glocal alignment scheduled")
            else:
                print("Could not start the Glocal alignment script:\n" + out)
                raise Exception("Could not start the Glocal alignment script: " + out)
        else:
            print("Could not start the reconstruction script:\n" + out)
            raise Exception("Could not start the reconstruction script: " + out)
    else:
        if out.startswith("Submitted batch job "):
            print("Job scheduled succesfully")
        else:
            print("There seems to be a problem with scheduling the job, output: " + out)


def normalised_cross_correlation(first, second, filter_mask=None):
    """
    Do a cross correlation based on numpy

    @param first: The first dataset (numpy 2D)
    @type first: numpy array 2D
    @param second: The second dataset (numpy 2D)
    @type second: numpy array 2D
    @param filter_mask: a filter which is used to filter image 'first'
    @type filter_mask: numpy array 2D
    @return: The cross correlation result
    @returntype: numpy array 2D

    @requires: the shape of first to be equal to the shape of second, and equal t the shape of the filter (if used of course)
    """
    import sys
    from pytom.agnostic.io import write
    assert first.shape == second.shape
    assert len(first.shape) == 2
    if not(filter_mask is None): assert first.shape == filter_mask.shape

    from numpy.fft import fftshift, ifftn, fftn
    import numpy as xp

    prod = 1
    for d in first.shape:
        prod *= d

    if filter_mask is None:
        ffirst = xp.fft.fftn(first)
        fsecond = xp.fft.fftn(second)
    else:
        ffirst = ((xp.fft.fftn(first)) * xp.fft.fftshift(filter_mask))
        fsecond = xp.fft.fftn(second) * xp.fft.fftshift(filter_mask)
    ccmap = xp.real(xp.fft.fftshift(xp.fft.ifftn( fsecond * xp.conj(ffirst)) )) / filter_mask.sum()

    # write('real_ccmap.mrc', ccmap)

    return ccmap

def normalised_cross_correlation_mask(first, second, mask):
    """
    Do cross correlation with a running mask based on numpy

    @param first: The first dataset (numpy 2D)
    @type first: numpy array 2D
    @param second: The second dataset (numpy 2D)
    @type second: numpy array 2D
    @param mask: The mask
    @type mask: numpy array 2D
    @return: The cross correlation result
    @returntype: numpy array 2D

    @requires: the shape of first to be equal to the shape of second and the shape of the mask
    """
    # assert first.shape == second.shape
    # assert first.shape == mask.shape
    # assert len(first.shape) == 2

    import numpy.fft as nf
    import numpy as np

    a = norm_inside_mask(first, mask)
    b = norm_inside_mask(second, mask)

    return np.real(nf.fftshift(nf.ifftn(np.multiply(nf.fftn(b), np.conj(nf.fftn(a)))))) / np.sum(mask)


def norm_inside_mask(inp, mask):
    """
    To normalise a 2D array within a mask

    @param inp: A 2D array to be normalized.
    @type inp: numpy array 2D
    @param mask: A 2D array of the same size as the input to mask of parts of the ixp.
    @type mask: numpy array 2D
    @return: A normalized 2D array.
    @returntype: numpy array 2D

    @requires: the shape of inp to be equal to the shape of the mask
    """
    # assert ixp.shape == mask.shape
    # assert len(ixp.shape) == 2

    import numpy as xp

    mea = xp.divide(xp.sum(xp.multiply(inp, mask)), xp.sum(mask))
    st = xp.sqrt(xp.sum((xp.multiply(mask, mea) + xp.multiply(inp, mask) - mea) ** 2) / xp.sum(mask))
    return xp.multiply((inp - mea) / st, mask)

def find_sub_pixel_max_value_voltools(inp, k=.1, ignore_border=50):
    """
    To find the highest point in a 2D array, with subpixel accuracy based on 1D spline interpolation.
    The algorithm is based on a matlab script "tom_peak.m"

    @param inp: A 2D numpy array containing the data points.
    @type inp: numpy array 2D
    @param k: The smoothing factor used in the spline interpolation, must be 1 <= k <= 5.
    @type k: int
    @return: A list of all points of maximal value in the structure of tuples with the x position, the y position and
        the value.
    @returntype: list
    """
    # assert len(ixp.shape) == 2
    # assert isinstance(k, int) and 1 <= k <= 5

    import cupy
    import numpy as xp
    from scipy.interpolate import InterpolatedUnivariateSpline
    from pytom.voltools import transform

    inp2 = inp[ignore_border:-ignore_border, ignore_border:-ignore_border]

    out = cupy.array(inp,dtype=cupy.float32)
    x,y = xp.array(xp.unravel_index(inp2.argmax(), inp2.shape))+ignore_border

    out = cupy.expand_dims(out,2)
    zoomed = cupy.zeros_like(out,dtype=cupy.float32)
    transform(out, output=zoomed, scale=(k,k,1), translation=[inp.shape[0]//2-x, inp.shape[1]//2-y,0],
              device='gpu:0', interpolation='filt_bspline')

    z = zoomed.get().squeeze()

    x2,y2 = xp.unravel_index(z.argmax(), z.shape)
    return [x + (x2-z.shape[0]//2)*k -z.shape[0]//2, y + (y2-z.shape[1]//2)*k -z.shape[0]//2, z, x, y, x2, y2]



def find_sub_pixel_max_value(inp, k=4, ignore_border=50):
    """
    To find the highest point in a 2D array, with subpixel accuracy based on 1D spline interpolation.
    The algorithm is based on a matlab script "tom_peak.m"

    @param inp: A 2D numpy array containing the data points.
    @type inp: numpy array 2D
    @param k: The smoothing factor used in the spline interpolation, must be 1 <= k <= 5.
    @type k: int
    @return: A list of all points of maximal value in the structure of tuples with the x position, the y position and
        the value.
    @returntype: list
    """
    # assert len(ixp.shape) == 2
    # assert isinstance(k, int) and 1 <= k <= 5

    import numpy as xp
    from scipy.interpolate import InterpolatedUnivariateSpline


    inp = inp[ignore_border:-ignore_border, ignore_border:-ignore_border]

    v = xp.amax(inp)  # the max value
    result = xp.where(inp == v)  # arrays of x and y positions of max values
    output = []

    for xpp, yp in zip(result[0], result[1]):
        # Find the highest point for x (first check if on sides otherwise interpolate)
        if xpp == 1 or xpp == inp.shape[0]:
            x = xpp
            xv = v
        else:
            f = InterpolatedUnivariateSpline(range(0, inp.shape[0]), inp[:, yp], k=k)  # spline interpolation
            cr_pts = f.derivative().roots()
            cr_vals = f(cr_pts)
            val = xp.argmax(cr_vals)
            x = cr_pts[val]
            xv = cr_vals[val]

        # Find the highest point for y (first check if on sides otherwise interpolate)
        if yp == 1 or yp == inp.shape[1]:
            y = yp
            yv = v
        else:
            f = InterpolatedUnivariateSpline(range(0, inp.shape[1]), inp[xpp, :], k=k)  # spline interpolation
            cr_pts = f.derivative().roots()
            cr_vals = f(cr_pts)
            val = xp.argmax(cr_vals)
            y = cr_pts[val]
            yv = cr_vals[val]

        # Calculate the average of the max value to return a single value which is maybe more close to the true value
        output.append([x+ignore_border, y+ignore_border, (xv + yv) / 2])

    return output[0]


def find_sub_pixel_max_value_2d(inp, interpolate_factor=100, smoothing=0, dim=10, border_size=5, ignore_border=37): #ignore_border for 200: 75
    """
    To find the highest point in a given numpy array based on 2d spline interpolation, returns the maximum with subpixel
    precision.

    @param inp: The input data array (numpy 2D)
    @type inp: numpy array 2D
    @param interpolate_factor: The amount of interpolation to be done
    @type interpolate_factor: int
    @param smoothing: The amount of smoothing in the spline interpolation
    @type smoothing: int
    @param dim: The dimensions of the peak cutout, which is interpolated to find the subpixel maximum (initial pixels)
    @type dim: int
    @param border_size: The amount of pixels (initial pixels) to disregard in the peak cutout
    @type border_size: int
    @param ignore_border: The amount of pixels (initial pixels) to disregard in the initial finding of the initial maximum, to force the
       found maximum to be more in the center
    @type ignore_border: int
    @return: The subpixel maximum (x, y, the interpolated peak (excluding the border area))
    @returntype: tuple

     <-------- Vol_Size ------->
    | ignore_border             |
    | |   <------ a ------>     |
    | -> |                 | <- |
    |    | Here the max is |    |
    |    | found           |    |
    |    |    d> <c> <d    |    |
    |    |   |.. max ..|   |    |
    |    |   |... * ...|   |    |
    |    |   |.........|   |    |
    |    |    <-- b -->    |    |
    |    -------------------    |
    |___________________________|

    a: vol_size - 2 * ignore_border     (original pixels)
    b: dim * 2                          (original pixels)
    c: b * interpolate_factor - 2 * d   (interpolated pixels)
    d: border_size * interpolate_factor (interpolated pixels)
    ...: interpolated values
    *: peak found

    """
    # assert len(ixp.shape) == 2
    # assert isinstance(interpolate_factor, int) and interpolate_factor > 0
    # assert isinstance(smoothing, float) and smoothing >= 0
    # assert isinstance(dim, int) and dim > 0
    # assert isinstance(border_size, int) and border_size > 0
    # assert isinstance(ignore_border, int) and ignore_border > 0

    import numpy as xp
    from scipy import interpolate
    import warnings

    border_size = border_size * interpolate_factor

    # Get the position of the initial maximum
    inp_without_border = inp[ignore_border:-ignore_border, ignore_border:-ignore_border]
    initial_max = xp.unravel_index(inp_without_border.argmax(), inp_without_border.shape)
    # Reset the coordinates to be relative to the original inp(ut)
    initial_max = (initial_max[0] + ignore_border, initial_max[1] + ignore_border)
    # imshow(inp[initial_max[0]-10: initial_max[0]+10, initial_max[1]-10: initial_max[1]+10])
    # show()

    # Get the starting points of the peak cutout
    x_dim = inp.shape[0]
    y_dim = inp.shape[1]
    x_start = max([0, initial_max[0] - dim])
    x_end = min([x_dim, initial_max[0] + dim])
    y_start = max([0, initial_max[1] - dim])
    y_end = min([y_dim, initial_max[1] + dim])

    # Create a grid to save the original points and one to save the interpolated points
    x, y = xp.mgrid[0:x_end-x_start,0:y_end-y_start]
    xnew, ynew = xp.mgrid[0:x_end-x_start:1/interpolate_factor, 0:y_end-y_start:1/interpolate_factor]


    #print(x_start, x_end, y_start, y_end)

    # Interpolate the points
    # While catching warnings from the pessimistic algorithm of interpolate,
    # which always thinks it needs too much memory
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tck = interpolate.bisplrep(x, y, inp[x_start:x_end, y_start:y_end], s=smoothing)
        interpolated_grid = interpolate.bisplev(xnew[:, 0], ynew[0,:], tck)
        cropped_inter_grid = interpolated_grid[border_size:-border_size, border_size:-border_size]
        result = xp.unravel_index(cropped_inter_grid.argmax(), cropped_inter_grid.shape)

        rx, ry = result
        ix,iy = interpolated_grid.shape
        # imshow(interpolated_grid)
        # show()
        # Reset the coordinates to point to a place in the original data array
        result = [x_start + (border_size+rx)/interpolate_factor, y_start + (border_size+ry)/interpolate_factor]

        return result[0], result[1], cropped_inter_grid


def create_fsc_mask(fsc, size):
    """Create a 2D mask based on an FSC curve so it can be used to do masked cross correlation"""
    from numpy import meshgrid, arange, sqrt, zeros_like, float32

    X, Y, Z = meshgrid(arange(size), arange(size), arange(size))

    X -= size // 2
    Y -= size // 2
    Z -= size // 2

    R = sqrt(X ** 2 + Y ** 2 + Z ** 2).astype(int)

    out = zeros_like(R).astype(int)

    for n, val in enumerate(fsc):
        out[R == n] = val

    return out[size // 2, :, :]


def gettime():
    """Get the current time, to display in stdout to keep an eye on runtime"""
    from time import gmtime, strftime
    return strftime("%H:%M:%S", gmtime())

