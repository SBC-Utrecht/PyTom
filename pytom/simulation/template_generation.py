#!/usr/bin/env pytom

import numpy as xp
import os
import pytom.simulation.physics as physics


def generate_template(structure_file_path, spacing, binning=1, modify_structure=False, apply_solvent_correction=False,
                      solvent_density=physics.AMORPHOUS_ICE_DENSITY, apply_ctf_correction=False, defocus=3E-6,
                      amplitude_contrast=0.07, voltage=300E3, Cs=2.7E-3, ctf_decay=0.4, display_ctf=False,
                      resolution=30, box_size=None, output_folder=''):
    """

    @param structure_file_path:
    @type structure_file_path: string
    @param spacing:
    @type spacing: float
    @param binning:
    @type binning: int
    @param modify_structure:
    @type modify_structure: bool
    @param solvent_correction:
    @type solvent_correction: bool
    @param solvent_density:
    @type solvent_density: float
    @param apply_ctf_correction:
    @type apply_ctf_correction: bool
    @param defocus:
    @type defocus: float
    @param amplitude_contrast:
    @type amplitude_contrast: float
    @param voltage:
    @type voltage: float
    @param Cs:
    @type Cs: float
    @param ctf_decay:
    @type ctf_decay: float
    @param display_ctf:
    @type display_ctf: bool
    @param apply_lpf:
    @type apply_lpf: bool
    @param resolution:
    @type resolution: float
    @param box_size:
    @type box_size: int
    @return:
    @rtype:
    """
    from pytom.simulation.microscope import create_ctf, display_microscope_function
    from pytom.tompy.transform import resize, fourier_filter
    from pytom.tompy.tools import paste_in_center
    from pytom.simulation.support import create_gaussian_low_pass
    from pytom.simulation.potential import iasa_integration, call_chimera

    assert binning >= 1, 'binning factor smaller than 1 is invalid'
    if solvent_correction:
        assert solvent_density > 0, 'solvent density smaller or equal to 0 is invalid'
        if solvent_density > physics.PROTEIN_DENSITY:
            print(f'WARNING: Solvent density larger than protein density {physics.PROTEIN_DENSITY}')

    if not resolution >= (2 * spacing * binning):
        print(f'Invalid resolution specified, changing to {2*spacing*binning}A')
        resolution = 2 * spacing * binning

    if modify_structure:
        print('Adding symmetry, removing water, and adding hydrogen to pdb for template generation')
        # call chimera to modify the input structure, call_chimera returns the path to the updated structure
        structure_file_path  = call_chimera(structure_file_path, output_folder)

    # generate electrostatic_potential
    # iasa_generation returns a list of the real (and imaginary) part
    template = iasa_integration(structure_file_path, voxel_size=spacing, solvent_exclusion=apply_solvent_correction,
                                 V_sol = physics.V_WATER * (solvent_density/physics.AMORPHOUS_ICE_DENSITY))

    # extend volume to the desired input size before applying convolutions!
    if box_size is not None:
        if box_size > template.shape[0]//binning:
            new_box     = xp.zeros((box_size*binning,)*3)
            template    = paste_in_center(template, new_box)
        elif box_size < template.shape[0]//binning:
            print(f'Box size from electrostatic potential generation is {template.shape[0]//binning} px, which is larger than '
                  f'the requested box size of {box_size} px. We will not reduce the box size with the risk of cutting'
                  f'parts of the molecule.')

    # create ctf function and low pass gaussian if desired
    # for ctf the spacing of pixels/voxels needs to be in meters (not angstrom)
    ctf = create_ctf(template.shape, spacing*1E-10, defocus, amplitude_contrast, voltage, Cs,
                     sigma_decay=ctf_decay) if apply_ctf_correction else 1
    lpf = create_gaussian_low_pass(template.shape, (template.shape[0]*spacing)/resolution)

    # print information back to user
    if apply_ctf_correction: print(f'Applying ctf correction with defocus {defocus*1e6} um')
    print(f'Applying low pass filter to {resolution}A resolution')
    # apply ctf and low pass in fourier space
    filter = lpf * ctf
    if display_ctf:
        print('Displaying combined ctf and lpf frequency modulation')
        display_microscope_function(filter[...,filter.shape[2]//2], form='ctf*lpf')
    template = fourier_filter(template, filter, human=True)

    # binning
    if binning > 1:
        print(f'Binning volume {binning} times')
        template = resize(template, 1/binning, interpolation='Spline')

    return template


if __name__ == '__main__':
    # parameters: file_path, destination, spacing, binning (optional, default is 1), solvent_correction (optional),
    # solvent_density (optional, default 0.93),
    # apply_ctf_correction (optional), defocus (optional, default is 3 um, negative is overfocus),
    # amplitude contrast (optional, default 0.07), voltage (optional, default is 300 keV),
    # Cs (optional, default is 2.7 mm), ctf_decay (optional, default is 0.4), display_ctf (optional),
    # apply_low_pass (optional), resolution_filter (optional, default is 2*spacing*binning), box_size (optional)

    # IN ORDER TO FUNCTION, SCRIPT REQUIRES INSTALLATION OF PYTOM (and dependencies), CHIMERA

    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.tompy.io import write

    # syntax is ScriptOption([short, long], description, requires argument, is optional)
    options = [ScriptOption(['-f', '--file'], 'Protein structure file, either pdb or cif.', True, False),
               ScriptOption(['-d', '--destination'], 'Folder where output should be stored.', True, False),
               ScriptOption(['-s', '--spacing'], 'The pixel spacing of original projections of the dataset in A,'
                                                 ' e.g. 2.62', True, False),
               ScriptOption(['-b', '--binning'], 'Number of times to bin the template. Default is 1 (no binning). If '
                                                 'set to 2 with a spacing of 2.62 the resulting voxel size will '
                                                 'be 5.24', True, True),
               ScriptOption(['-m', '--modify_structure'], 'Activate to call Chimera for adding hydrogen, symmetry and'
                                                          'removing water molecules from the structure.', False, True),
               ScriptOption(['-w', '--solvent_correction'], 'Whether to exclude solvent around each atom as a '
                                                            'correction of the potential.', False, True),
               ScriptOption(['-r', '--solvent_density'], 'Density of solvent, value should not be higher than 1.35 as'
                                                         ' that is the density of proteins. Default is 0.93 g/cm^3.',
                            True, True),
               ScriptOption(['-c', '--ctf_correction'], 'Correct the volume by applying a CTF. Default parameters are '
                                                      'defocus 3 um, amplitude contrast 0.07, voltage 300 keV, '
                                                      'spherical abberation (Cs) 2.7 mm, sigma of gaussian decay 0.4, '
                                                      'optionally plot the CTF to inspect.', False, True),
               ScriptOption(['-z', '--defocus'], 'Defocus in um (negative value is overfocus).', True, True),
               ScriptOption(['-a', '--amplitude_contrast'], 'Amplitude contrast fraction.', True, True),
               ScriptOption(['-v', '--voltage'], 'Acceleration voltage in keV', True, True),
               ScriptOption(['-o', '--Cs'], 'Spherical abberration in mm.', True, True),
               ScriptOption(['-g', '--decay'], 'Sigma of gaussian CTF decay function, 0.4 default.', True, True),
               ScriptOption(['-p', '--plot'], 'Give this option for plotting the CTF for visual inspection.', False,
                            True),
               ScriptOption(['-l', '--lpf_resolution'], 'Specify the resolution of the low pass filter that is applied.'
                                                        'The default value is 2 x spacing x binning (in angstrom), a '
                                                        'smaller resolution than this cannot be selected.', True, True),
               ScriptOption(['-x', '--xyz'], 'Specify a desired size for the output box of the template in number of '
                                          'pixels. By default the molecule is placed in a box with 30A overhang. This '
                                          'usually does not offer enough room to apply a spherical mask.', True, True),
               ScriptOption(['-h', '--help'], 'Help.', False, True)]

    helper = ScriptHelper(sys.argv[0].split('/')[-1], description='Create a template from the specified structure file '
                                                                  'with options for ctf correction and filtering. \n'
                                                                  'Script has dependencies on pytom and chimera.',
                          authors='Marten Chaillet', options=options)
    if len(sys.argv) == 2:
        print(helper)
        sys.exit()
    try:
        filepath, output_folder, spacing, binning, modify_structure, solvent_correction, solvent_density, \
            ctf_correction, defocus, amplitude_contrast, voltage, Cs, sigma_decay, \
            display_ctf, resolution, box_size, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()

    if help:
        print(helper)
        sys.exit()

    if spacing: spacing = float(spacing)

    if binning: binning = int(binning)
    else: binning = 1

    if solvent_density: solvent_density = float(solvent_density)
    else: solvent_density = physics.AMORPHOUS_ICE_DENSITY

    if defocus: defocus = float(defocus)*1E-6
    else: defocus = 3E-6

    if amplitude_contrast: amplitude_contrast = float(amplitude_contrast)
    else: amplitude_contrast = 0.07

    if voltage: voltage = float(voltage)*1E3
    else: voltage = 300E3

    if Cs: Cs = float(Cs)*1E-3
    else: Cs = 2.7E-3

    if sigma_decay: sigma_decay = float(sigma_decay)
    else: sigma_decay = 0.4

    if resolution: resolution = float(resolution)
    else: resolution = 2 * spacing * binning

    if box_size: box_size=int(box_size)

    if not os.path.exists(filepath):
        print('Protein structure file does not exist!')
        sys.exit()
    if not os.path.exists(output_folder):
        print('Output folder does not exist!')
        sys.exit()

    template = generate_template(filepath, spacing,
                                 binning=binning,
                                 modify_structure=modify_structure,
                                 apply_solvent_correction=solvent_correction,
                                 solvent_density=solvent_density,
                                 apply_ctf_correction=ctf_correction,
                                 defocus=defocus,
                                 amplitude_contrast=amplitude_contrast,
                                 voltage=voltage,
                                 Cs=Cs,
                                 ctf_decay=sigma_decay,
                                 display_ctf=display_ctf,
                                 resolution=resolution,
                                 box_size=box_size,
                                 output_folder=output_folder)

    # output structure
    _, file = os.path.split(filepath)
    id, _ = os.path.splitext(file)
    output_filepath = os.path.join(output_folder, f'template_{id}_{spacing*binning:.2f}A_{template.shape[0]}px.mrc')
    print(f'Writing template as {output_filepath}')
    write(output_filepath, template)
