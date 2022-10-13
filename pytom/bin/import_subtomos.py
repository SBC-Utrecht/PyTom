#!/usr/bin/env pytom

import sys, os
import multiprocessing
from pytom.tools.script_helper import ScriptHelper2, ScriptOption2
from pytom.tools.parse_script_options import parse_script_options2
from pytom.agnostic.io import read_star, read, write, read_pixelsize
from pytom.basic.structures import ParticleList, Particle, PickPosition, Rotation, Shift, Wedge
from pytom.basic.score import FLCFScore
from pytom.agnostic.tools import convert_angles
from pytom.agnostic.transform import fourier_reduced2full, fourier_full2reduced
from functools import partial


def convert_subtomos(pl, invert_contrast=False):
    for i, p in enumerate(pl):
        if i > 0 and i % 1000 == 0:
            print(i, ' done')

        if invert_contrast:  # invert contrast if needed
            subtomo_file = p.getFilename()
            folder, filename = os.path.split(subtomo_file)
            basename, ext = os.path.splitext(filename)
            new_subtomo_file = os.path.join(folder, basename + '_pytom' + ext)
            pixel_size = read_pixelsize(subtomo_file)
            subtomo = read(subtomo_file)
            write(new_subtomo_file, -1 * subtomo, pixel_size=pixel_size)
            # p.setFilename(new_subtomo_file)

        # change fourier reduced dimension of 3d ctf
        ctf_file = p.getWedge().getFilename()
        folder, filename = os.path.split(ctf_file)
        basename, ext = os.path.splitext(filename)
        new_ctf_file = os.path.join(folder, basename + '_pytom' + ext)
        pixel_size = read_pixelsize(ctf_file)
        ctf = read(ctf_file)
        ctf[ctf < 0] = 0
        ctf[ctf > 1] = 1
        full_ctf = fourier_reduced2full(ctf, reduced_axis=0, isodd=ctf.shape[1] % 2)
        red_ctf = fourier_full2reduced(full_ctf, reduced_axis=2)
        write(new_ctf_file, red_ctf, pixel_size=pixel_size)
        # p.getWedge.setFilename(new_ctf_file)


def import_subtomos(star_file, xml_file, invert_contrast=False, n_cores=8):
    # read star only handles relion3.0/3.1 currently
    stardata = read_star(star_file)
    particle_list = ParticleList()

    for n in range(len(stardata['CoordinateX'])):
        # set filename
        p = Particle(filename=stardata['ImageName'][n])

        # ===== pick position in tomogram
        p.setPickPosition(PickPosition(stardata['CoordinateX'][n],
                                       stardata['CoordinateY'][n],
                                       stardata['CoordinateZ'][n],
                                       originFilename=stardata['MicrographName'][n]))

        # ======= set shifts
        if 'OriginXAngst' in stardata.dtype.names:
            star_pixel_size = stardata['PixelSize'][n]
            x_shift, y_shift, z_shift = stardata['OriginXAngst'][n], stardata['OriginYAngst'][n], \
                                        stardata['OriginZAngst'][n]

            p.setShift(Shift(x=x_shift / star_pixel_size, y=y_shift / star_pixel_size, z=z_shift / star_pixel_size))

            # factor = binningWarpM / binningPyTom
            # p.getShift().setX(x * factor)
            # p.getShift().setY(y * factor)
            # p.getShift().setZ(z * factor)

        # ======= class info from relion
        if 'GroupNumber' in stardata.dtype.names:
            p.setClass(stardata['GroupNumber'][n])

        # ======= convert and set rotation
        # relion uses internal axis, while pytom external, or the other way around
        # anyway, need to invert the values
        z0 = -stardata['AngleRot'][n]
        y = -stardata['AngleTilt'][n]
        z1 = -stardata['AnglePsi'][n]

        z0, x, z1 = convert_angles((z0, y, z1), rotation_order='zyz', return_order='zxz')

        p.setRotation(Rotation(z1=z0, x=x, z2=z1, paradigm='ZXZ'))

        # ====== set wedge angles
        # need to update the 3d ctf so that pytom can use it
        # this means going reduced2full fourier space in X, and then full2reduced in Z
        ctf_file = stardata['CtfImage'][n]
        p.setWedge(Wedge(wedge_3d_ctf_file=ctf_file))

        # === give particle empty score
        p.setScore(FLCFScore())

        particle_list.append(p)

    print('each process will print some output about its progress')
    print('total number of particles to loop over: ', len(particle_list))
    # map cores to convert the files onto the split particle list
    with multiprocessing.Pool(processes=n_cores) as pool:
        pool.map(partial(convert_subtomos, invert_contrast=invert_contrast), particle_list.splitNSublists(n_cores))

    # update filenames in unsplit particle list
    for p in particle_list:
        if invert_contrast:  # invert contrast if needed
            subtomo_file = p.getFilename()
            folder, filename = os.path.split(subtomo_file)
            basename, ext = os.path.splitext(filename)
            new_subtomo_file = os.path.join(folder, basename + '_pytom' + ext)
            p.setFilename(new_subtomo_file)

        # change fourier reduced dimension of 3d ctf
        ctf_file = p.getWedge().getFilename()
        folder, filename = os.path.split(ctf_file)
        basename, ext = os.path.splitext(filename)
        new_ctf_file = os.path.join(folder, basename + '_pytom' + ext)
        p.getWedge().setFilename(new_ctf_file)

    particle_list.toXMLFile(xml_file)


if __name__ == '__main__':
    helper = ScriptHelper2(
        sys.argv[0].split('/')[-1],  # script name
        description='Import subtomograms from other software. Currently supported: WarpM',
        authors='Marten Chaillet',
        options=[
            ScriptOption2(['-i', '--input-star-file'], 'Input particle list in .star format',
                          'file', 'required'),
            ScriptOption2(['-o', '--output-xml-file'], 'PyTom xml file that will store the particle information',
                          'string', 'optional'),
            ScriptOption2(['--invert'], 'Invert contrast of subtomograms. PyTom wants particles to be negative, '
                                        'while Relion wants positive particles.',
                          'no arguments', 'optional'),
            ScriptOption2(['-c', '--cores'], 'number of cpu cores to split processes on', 'int', 'optional', 1)])

    options = parse_script_options2(sys.argv[1:], helper)
    input_star, output_xml, invert_contrast_flag, n_cores = options

    # make boolean
    invert_contrast_flag = False if invert_contrast_flag is None else True

    # set some paths for the output
    if output_xml is None:
        folder, filename = os.path.split(input_star)
        file_basename = os.path.splitext(filename)[0]
        output_xml = os.path.join(folder, file_basename + '.xml')

    import_subtomos(input_star, output_xml, invert_contrast=invert_contrast_flag, n_cores=n_cores)
