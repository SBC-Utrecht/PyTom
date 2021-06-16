import numpy as xp
import pytom.simulation.physics as physics


def generate_template(structure_file_path, spacing, binning=1, modify_structure=False, apply_solvent_correction=False,
                      solvent_density=physics.AMORPHOUS_ICE_DENSITY, apply_ctf_correction=False, defocus=3E-6,
                      amplitude_contrast=0.07, voltage=300E3, Cs=2.7E-3, ctf_decay=0.4, display_ctf=False,
                      resolution=30, box_size=None, output_folder=''):
    """

    Generating a template for template matching in a set of tomograms. Spacing provided should correspond to original
    spacing of the projections of the dataset. The binning factor can be applied identical to the binning factor used
    for reconstructing the tomograms. It is better to do this because the template will be sampled at higher
    resolution and should therefore be more accurate.

        default values for ctf correction
    defocus                 = 3e-6 m
    amplitude contrast      = 0.07
    voltage                 = 300e3
    spherical aberration    = 2.7e-3
    ctf decay               = 0.4

    @param structure_file_path: full path to file specifying structure in pdb or cif format
    todo some issues with cif format have been reported
    @type  structure_file_path: L{str}
    @param spacing: spacing of voxels in A
    @type  spacing: L{float}
    @param binning: binning factor to apply after generating template first at specified spacing
    @type  binning: L{int}
    @param modify_structure: flag to modify structure by adding hydrogen and adding symmetry using UCSF Chimera.
    Chimera needs to be on PATH in terminal to execute this option.
    @type  modify_structure: L{bool}
    @param apply_solvent_correction: flag to apply solvent correction
    @type  apply_solvent_correction: L{bool}
    @param solvent_density: solvent density to use for the correction, default is amorphous ice density 0.93 g/cm^3
    @type  solvent_density: L{float}
    @param apply_ctf_correction: flag to apply ctf correction to template
    @type  apply_ctf_correction: L{bool}
    @param defocus: defocus value in m, dz > 0 is defocus, dz < 0  is overfocus
    @type  defocus: L{float}
    @param amplitude_contrast: fraction of amplitude contrast
    @type  amplitude_contrast: L{float}
    @param voltage: electron beam acceleration voltage in eV
    @type  voltage: L{float}
    @param Cs: spherical aberration in m
    @type  Cs: L{float}
    @param ctf_decay: sigma of Gaussian ctf decay
    @type  ctf_decay: L{float}
    @param display_ctf: flag to display a plot of the ctf
    @type  display_ctf: L{bool}
    @param resolution: resolution to apply low-pass filter to in A, default is 2 * spacing * binning. values lower
    than the default will be overridden by default value because low-pass filtering the template is always better.
    @type  resolution: L{float}
    @param box_size: force box size of template to be larger than default generation
    @type  box_size: L{int}

    @return: template containg volume, 3d array of floats
    @rtype:  L{np.ndarray}

    @author: Marten Chaillet
    """
    from pytom.simulation.microscope import create_ctf, display_microscope_function
    from pytom.tompy.transform import resize, fourier_filter
    from pytom.tompy.tools import paste_in_center
    from pytom.simulation.support import create_gaussian_low_pass
    from pytom.simulation.potential import iasa_integration, call_chimera

    assert binning >= 1, 'binning factor smaller than 1 is invalid'
    if apply_solvent_correction:
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

