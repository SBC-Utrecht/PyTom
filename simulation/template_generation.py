
import numpy as np


def ctf_array(shape, nyquist):
    """
    @param shape: tuple with the size on each dimension
    @param shape: tuple of ints
    @param spacing: pixel or voxel spacing
    @type spacing: float
    @return: fourier space frequencies
    @rtype: numpy.ndarray
    """
    assert type(shape) == tuple, "shape needs to be a tuple of dimenion sizes, dimenions need to be equal"
    assert len(set(shape)) == 1, "dimensions are not identical"
    assert len(shape) <= 3, "too many dimensions, 3 is max"

    size = shape[0]

    d = np.arange(-nyquist, nyquist, 2. * nyquist / size)

    if len(shape) == 2:
        grids = np.meshgrid(d,d)
    else:
        grids = np.meshgrid(d, d, d)

    # Wave vector and direction
    k = np.sqrt(sum([d**2 for d in grids]))  # xz ** 2 + yz ** 2 + zz ** 2)
    # k2 = k ** 2
    # k4 = k2 ** 2
    return k  # , k2, k4


def create_ctf(shape, spacing, defocus, amplitude_contrast, voltage, spherical_abberration, sigma_decay=0.4):
    """
    Models a CTF.
    INPUT
    @param Dz: defocus value in m
    @type Dz: C{float}
    @param A: Amplitude contrast fraction (e.g., 0.1)
    @type A: C{float}
    @param Voltage: acceleration voltage in eV
    @type Voltage: C{float}
    @param Cs: sphertical abberation in m
    @type Cs: C{float}
    @param k4: frequencies in fourier space to power 4
    @type k4: L{numpy.ndarray}
    @param k2: frequencies in fourier space squared
    @type k2: L{numpy.ndarray}
    @return: CTF in Fourier space
    @rtype: L{numpy.ndarray}
    """
    from pytom.simulation.simulateProjections import wavelength_eV2m
    # k4, k2, k,
    nyquist = 1 / (2 * spacing)
    k = ctf_array(shape, nyquist)

    # lmbd = np.sqrt(150.4 / ((voltage * (1 + voltage / 1022000.)))) * 1E-10
    lmbd = wavelength_eV2m(voltage)

    chi = np.pi * lmbd * defocus * k**2 - 0.5 * np.pi * spherical_abberration * lmbd ** 3 * k ** 4

    ctf = -np.sqrt(1. - amplitude_contrast ** 2) * np.sin(chi) + amplitude_contrast * np.cos(chi)

    if sigma_decay > 0:
        decay = np.exp(-(k / (sigma_decay * nyquist)) ** 2)
    else:
        decay = 1

    return ctf * decay

#     # Electron wavelength
#
#     lmbd = np.sqrt(150.4 / ((voltage * (1 + voltage / 1022000.)))) * 1E-10
#
#     # Phase shift due to lens aberrations and defocus
#     chi = np.pi * lmbd * k2 * Dz - 0.5 * np.pi * lmbd ** 3 * Cs * k4
#
#     # Contrast transfer function
#     ctf = -np.sqrt(1. - A ** 2) * np.sin(chi) + A * np.cos(chi)
#     # TODO what is the point of sqrt(1-A^2)? approximately equal to for 1 for A=0.08
#
#     # CHANGED sign for amplitude contrast to fix phase flip for low frequencies -- SP 7.9.16
#     # originally: ctf = -np.sqrt(1-A**2)*np.sin(chi)-A*np.cos(chi)
#     # Decay function
#     Ny = 1 / (2 * spacing)
#     decay = np.exp(-(k / (sigma_decay_CTF * Ny)) ** 2)
#
#     return ctf * decay
#
#
# def generate_template(structure, pixel_size, solvent_correction=False, ctf_correction=False, lpf=False):
#
#     # electrostatic_potential
#     structure
#     pixel_size
#
#     solvent_correction # with solvent masking (requires oversampling?)
#     solvent_density = 0.93 # default
#
#
#
#     ctf_correction
#     defocus = 3  # default
#     voltage = 300  # default
#
#
#     lpf
#     resolution = 2 * pixel_size * binning # default
#
#     binning
