#!/usr/bin/env python

"""
Created on July 12, 2021

This script is a conversion of tom_deconv.m to pytom. tom_deconv was written by Dimitry Tegunov and can be found on
GitHub: https://github.com/dtegunov/tom_deconv.

@author: Marten Chaillet
"""

if __name__ == '__main__':
    # parse command line arguments with ScriptHelper2
    import sys
    from pytom.tools.script_helper import ScriptHelper2, ScriptOption2
    from pytom.tools.parse_script_options import parse_script_options2
    from pytom.agnostic.io import read, write
    from pytom.agnostic.filter import wiener_like_filter
    from pytom.gpu.initialize import xp, device
    from pytom.agnostic.transform import resize

    helper = ScriptHelper2(
        sys.argv[0].split('/')[-1],  # script name
        description='Deconvolute a tomogram or micrograph with wiener-like filter',
        authors='Marten Chaillet (based on tom_deconv, dtegunov)',
        options=[
            ScriptOption2(['-f', '--file'], 'Tomogram or projection file (mrc)', 'file', 'required'),
            ScriptOption2(['-o', '--output'], 'Output file path', 'string', 'required'),
            ScriptOption2(['-s', '--spacing'], 'Spacing of pixels or voxels in angstrom', 'float', 'required'),
            ScriptOption2(['-z', '--defocus'], 'Defocus value of tomogram, defocus is postive, overfocus is '
                                               'negative', 'float', 'required'),
            ScriptOption2(['--snrfalloff'], 'How fast does SNR fall off, i. e. higher values will downweight high '
                                            'frequencies; values like 1.0 or 1.2 seem reasonableFall off snr',
                          'float', 'optional', 1.0),
            ScriptOption2(['--deconvstrength'], 'How much will the signal be deconvoluted overall, i. e. a global '
                                                'scale for SNR; exponential scale: 1.0 is SNR = 1000 at zero '
                                                'frequency, 0.67 is SNR = 100, and so on)', 'float', 'optional', 1.0),
            ScriptOption2(['--highpassnyquist'], 'fraction of Nyquist frequency to be cut off on the lower end ('
                                                 'since it will be boosted the most)', 'float', 'optional', 0.02),
            ScriptOption2(['--phaseflipped'], 'Whether data are already phase flipped', 'no arguments', 'optional'),
            ScriptOption2(['--phaseshift'], 'CTF phase shift in degrees (from a phase plate', 'float', 'optional', .0),
            ScriptOption2(['-v', '--voltage'], 'Acceleration voltage in keV', 'float', 'optional', 300.),
            ScriptOption2(['--Cs'], 'Spherical aberration of microscope in mm', 'float', 'optional', 2.7),
            ScriptOption2(['-a', '--amplitude'], 'Fraction of amplitude contrast', 'float', 'optional', 0.07),
            ScriptOption2(['-b', '--binning'], 'Binning factor for tomogram or projection', 'int', 'optional', 1),
            ScriptOption2(['-g', '--gpuID'], 'GPU ID to run code on', 'int', 'optional')])

    options = parse_script_options2(sys.argv[1:], helper)

    input_file, output_file, spacing_angstrom, defocus_um, snrfalloff, deconvstrength, highpassnyquist, phaseflipped, \
        phaseshift, voltage_kev, Cs_mm, amplitude_contrast, binning, gpu = options

    input = read(input_file)

    if binning > 1:
        input = resize(input, 1.0/binning)
        spacing_angstrom = spacing_angstrom * binning

    filter = wiener_like_filter(input.shape, spacing_angstrom, defocus_um * 1e-6, snrfalloff, deconvstrength,
                                highpassnyquist, voltage=voltage_kev * 1e3, spherical_aberration=Cs_mm * 1e-3,
                                amplitude_contrast=amplitude_contrast, phaseflipped=phaseflipped,
                                phase_shift=phaseshift)

    # TODO use reduced fft
    deconv = xp.fft.ifftn(xp.fft.fftn(input) * xp.fft.ifftshift(filter)).real
    write(output_file, deconv)
