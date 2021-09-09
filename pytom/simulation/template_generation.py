#!/usr/bin/env pytom

import os
import sys
import pytom.simulation.physics as physics
from pytom.tools.script_helper import ScriptHelper, ScriptOption
from pytom.tools.parse_script_options import parse_script_options
from pytom.tompy.io import write
from pytom.simulation.template import generate_template


if __name__ == '__main__':
    # todo template generation could be updated with an option for spherical mask generation

    # parameters: file_path, destination, spacing, binning (optional, default is 1), solvent_correction (optional),
    # solvent_density (optional, default 0.93),
    # apply_ctf_correction (optional), defocus (optional, default is 3 um, negative is overfocus),
    # amplitude contrast (optional, default 0.07), voltage (optional, default is 300 keV),
    # Cs (optional, default is 2.7 mm), ctf_decay (optional, default is 0.4), display_ctf (optional),
    # apply_low_pass (optional), resolution_filter (optional, default is 2*spacing*binning), box_size (optional)

    # IN ORDER TO FUNCTION, SCRIPT REQUIRES INSTALLATION OF PYTOM (and dependencies), CHIMERA

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

