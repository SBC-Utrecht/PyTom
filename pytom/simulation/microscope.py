import os
from pytom.gpu.initialize import xp, device
import pytom.simulation.physics as physics
from scipy.optimize import curve_fit


def convert_defocusU_defocusV_to_defocus_astigmatism(defocusU, defocusV):
    return 0.5 * (defocusU + defocusV), -0.5 * (defocusU - defocusV)


def convert_defocus_astigmatism_to_defocusU_defocusV(defocus, astigmatism):
    return defocus + astigmatism, defocus - astigmatism


# TODO make this a GPU kernel / numba function as xp.meshgrid has a lot of memory overhead
def fourier_grids(shape, nyquist, indexing='ij', reduced=False):
    """
    Generate a fourier space frequency array where values range from -nyquist to +nyquist, with the center equal
    to zero.

    In combination with angular_grid() the indexing parameter gives:
     - 'ij' produces the angles assuming the x axis is 0
     - 'xy' produces angles assuming the y axis is 0

    @param shape: shape tuple of grid
    @type  shape: L{tuple} -> (L{int},) * n
    @param nyquist: nyquist frequency in fourier space
    @type  nyquist: L{float}
    @param indexing: type of np.meshgrid indexing, either 'ij' or 'xy', if unsure stick with default 'ij'
    @type  indexing: L{str}

    @return: fourier space frequencies, 1d, 2d or 3d array of floats
    @rtype:  L{numpy.ndarray}

    @author: Marten Chaillet, Gijs van der Schot
    """
    assert 1 <= len(shape) <= 3, print('invalid argument for number of dimensions of fourier array')

    # xp.arange(-1, 1, 998) returns an array of length 999, expression xp.arange(size) / (size/2) - 1 solves this
    if reduced:
        d = []
        for i, size in enumerate(shape):
            if i == (len(shape) - 1):
                d.append((xp.arange(size // 2 + 1) / (size // 2) - 1) * nyquist)
            else:
                d.append((xp.arange(size) / (size // 2) - 1) * nyquist)
    else:
        d = [(xp.arange(size) / (size // 2) - 1) * nyquist for size in shape]

    return xp.meshgrid(*d, indexing=indexing)


def normalised_grid(shape, reduced=False):
    grids = fourier_grids(shape, 1, reduced=reduced)
    return xp.sqrt(sum([g**2 for g in grids]))


def ctf_grids(grids):
    """
    Create coordinate grid from fourier_grids (meshgrid) output. Return the grid, grid**2, and grid**4 for ctf calc.

    @param grids: output from fourier grids
    @type  grids: L{tuple} -> (L{np.ndarray},) * n

    @return: three arrays, k, k^2, k^4
    @rtype : L{tuple} -> (L{np.ndarray},) * n

    @author: Marten Chaillet
    """
    k = sum([d ** 2 for d in grids])  # xz ** 2 + yz ** 2 + zz ** 2)
    return xp.sqrt(k), k, k ** 2


def angular_grid(gridx, gridy):
    """
    Create an angular grid, where the x-axis is the 0 angle and the rotation is clockwise.

    From fourier_grids() indexing:
     - 'ij' produces the angles assuming the x axis is 0
     - 'xy' produces angles assuming the y axis is 0

    @param gridx: x grid for fourier_grids (or meshgrid)
    @type  gridx: L{np.ndarray}
    @param gridy: y grid for fourier_grids (or meshgrid)
    @type  gridy: L{np.ndarray}

    @return: angular grid in radians
    @rtype : L{np.ndarray}

    @author: Marten Chaillet
    """
    # correct expression is: x - 1j * y
    # TODO check compared to ctf estimation CTFFIND to see if astigmatism is correct
    return xp.angle(gridx - 1j * gridy)


def precalculate(shape, nyquist):
    grids = fourier_grids(shape, nyquist)
    k, k2, k4 = ctf_grids(grids)
    angles = angular_grid(grids[0], grids[1])
    return k, k2, k4, angles


def astigmatised_defocus_grid(angular_grid, defocusU, defocusV, astigmatism_angle_rad):
    """
    Create an astigmatised defocus grid that can be used in CTF calculations. The astigmatism angle rotation is
    relative to the x-axis.
    NOTE: astigmatised grids are only possible for 2-dimensional CTF

    @param angular_grid: grid with angles, where the x-axis is the 0 angle
    @type  angular_grid: L{np.ndarray}
    @param defocusU: first defocus value
    @type  defocusU: L{float}
    @param defocusV: second defocus value
    @type  defocusV: L{float}
    @param astigmatism_angle_deg: astigmatism angle in degrees
    @type  astigmatism_angle_deg: L{float}

    @return: a defocus grid with astigmatism
    @rtype : L{np.ndarray}

    @author: Marten Chaillet
    """
    return 0.5 * (defocusU + defocusV +
                  (defocusU - defocusV) * xp.cos(2 * (angular_grid - astigmatism_angle_rad)))


def generate_astigmatised_grid(k, Y=None, defocusU=None, defocusV=None, defocus_angle_deg=0., nyquist=None):
    """
    Generate an astigmatised defocus grid with wavevector k and fourier grid in y. If defocusV is not provided a non
    astigmatic grid is created.

    @param k: fourier space wave vector
    @type  k: L{np.ndarray}
    @param Y: fourier grid in y dimension (optional), by default determined from k but providing it is more precise
    @type  Y: L{np.ndarray}
    @param defocusU: defocus in m
    @type  defocusU: L{float}
    @param defocusV: defocus along astigmatism angle in m (optional)
    @type  defocusV: L{float}
    @param defocus_angle_deg: astigmatism angle in degrees
    @type  defocus_angle_deg: L{float}
    @param nyquist: optional nyquist frequencies
    @type  nyquist: L{float}

    @return: return the astigmatised defocus grid
    @rtype:  L{np.ndarray}

    @author: Gijs van der Schot, Marten Chaillet
    """
    size = k.shape[1]
    defocusV = defocusU if defocusV is None else defocusV
    if Y is None:
        if nyquist is None:
            nyquist = k.min(axis=0)
        Y = fourier_grids(k.shape, nyquist)[1]

    theta = xp.arcsin(Y / k)

    z = defocusU * xp.cos(theta - defocus_angle_deg) ** 2 + defocusV * xp.sin(theta - defocus_angle_deg) ** 2
    z[:, -int(size // 2):] = xp.flipud(z[:, -int(size // 2):])

    return z


def sinc_square(x, p1, p2, p3, p4):
    """
    Sinc square function: M{f(x) = p1 * (sin(S{pi} * (x - p2) * p3) / (S{pi} * (x - p2) * p3))^2 + p4 }

    @param x: coordinate
    @type  x: L{float}
    @param p1: parameter 1
    @type  p1: L{float}
    @param p2: parameter 2
    @type  p2: L{float}
    @param p3: parameter 3
    @type  p3: L{float}
    @param p4: parameter 4
    @type  p4: L{float}

    @return: evaluation of sinc square at x
    @rtype:  L{float}

    @author: Marten Chaillet
    """
    return p1 * (xp.sin(xp.pi*(x-p2)*p3) / (xp.pi*(x-p2)*p3))**2 + p4


def fit_sinc_square(xdata, ydata):
    """
    Fit a sinc square function to a set of x,y coordinates provided in xdata and ydata. xdata and ydata should be
    same length and each index in xdata should correspond to the index in ydata.

    @param xdata: x coordinates, 1d array of floats
    @type  xdata: L{np.ndarray}
    @param ydata: y coordinates, 1d array of floats
    @type  ydata: L{np.ndarray}

    @return: the 4 fitted parameters that can be used in sinc_square()
    @rtype:  L{tuple} -> (L{float},) * 4

    @author: Marten Chaillet
    """
    assert len(xdata) == len(ydata), print("length of x and y coordinates lists of data to fit since square to does "
                                           "not match")
    # Here you give the initial parameters for p0 which Python then iterates over
    # to find the best fit
    popt, pcov = curve_fit(sinc_square, xdata, ydata, p0=(1.0, 1.0, 1.0, 1.0))  # THESE PARAMETERS ARE USER DEFINED

    # Performing sum of squares
    p1, p2, p3, p4 = popt[0], popt[1], popt[2], popt[3]
    residuals = ydata - sinc_square(xdata, p1, p2, p3, p4)
    fres = sum(residuals**2)

    print(f'chi-square value for fitting sinc function is {fres}')  # THIS IS YOUR CHI-SQUARE VALUE!

    # Visually inspect fit of function.
    #
    # import matplotlib
    # matplotlib.use('Qt5Agg')
    # import matplotlib.pyplot as plt
    #
    # xaxis = xp.linspace(0,1,100) # we can plot with xdata, but fit will not look good
    # curve_y = sinc_square(xaxis,p1,p2,p3,p4)
    # plt.plot(xdata,ydata,'*')
    # plt.plot(xaxis,curve_y,'-')
    # plt.show()

    return p1, p2, p3, p4


def radial_average(image):
    """
    This calculates the radial average of an image. When used for ctf, dqe, and mtf type display, input should be in
    Fourier space.

    @param image: input to be radially averaged, 2d array of floats
    @type  image: L{np.ndarray}

    @return: coordinates, values
    @rtype:  L{tuple} -> (L{np.ndarray},) * 2

    @author: Marten Chaillet
    """
    assert len(image.shape) == 2, "radial average calculation only works for 2d image arrays"
    assert len(set(image.shape)) == 1, 'differently size dimension, cannot perform radial averaging'

    size = image.shape[0]
    center = (size - 1) / 2
    x, y = xp.meshgrid(xp.arange(size), xp.arange(size))
    R = xp.sqrt((x - center) ** 2 + (y - center) ** 2)

    f = lambda r: image[(R >= r - .5) & (R < r + .5)].mean()
    r = xp.linspace(1, size // 2, num=size // 2)
    mean = xp.vectorize(f)(r)

    return r, mean


def display_microscope_function(image, form='', ylim=(-1, 1), complex=False):
    """
    Display the radial average of a microscope function. If complex flag is set to true the function can also accept
    complex valued inputs.

    todo radial average of non-square images?
    todo complex valued curve should maybe be displayed as amplitude and phase instead of real and imaginary part

    @param image: input to display radial average of, 2d array of floats
    @type  image: L{np.ndarray}
    @param form: name of the type of function that is displayed, will be used as a label for the plot
    @type  form: L{str}
    @param complex: flag for complex valued inputs
    @type  complex: L{bool}

    @return: - (image will be displayed)
    @rtype:  None

    @author: Marten Chaillet
    """
    import matplotlib
    try:
        matplotlib.use('Qt5Agg')
    except:
        pass

    import matplotlib.pyplot as plt

    if complex:
        r1, m1 = radial_average(image.real)
        r2, m2 = radial_average(image.imag)
    else:
        r1, m1 = radial_average(image)

    fig, ax = plt.subplots()
    ax.plot(r1, m1, label=form)
    ax.set_ylim(*ylim)
    if complex: ax.plot(r2, m2, label='imaginary')
    ax.legend()
    # plt.savefig(f'{form}.png')
    # plt.imshow()
    plt.show()
    return


def read_detector_csv(filename):
    """
    Read a csv file containing detector data. The data should be listed as rows with the first column the x
    coordinate and the second the y coordinate. Data should also be comma separated.

    @param filename: path to .csv file to read
    @type  filename: L{str}

    @return: x, y where x and y are both a 1d array with points
    @rtype:  L{tuple} -> (L{np.ndarray},) * 2

    @author: Marten Chaillet
    """
    data = xp.genfromtxt(filename, delimiter=',', dtype=xp.float32)
    x = data[:, 0]
    y = data[:, 1]
    return x, y


def create_detector_response(detector, response_function, image_size, voltage=300E3, oversampling=1, folder='',
                             display=False):
    """
    This function will read a CSV file containing a detector response function at a specific voltage, with format
    {detector}_{response_function}_{int(voltage*1E-3)}kV.csv . This function will be loaded from folder (if provided)
    and otherwise from the current working directory. The CSV file should contain only x and y coordinates,
    where each row has x on first column and y on second, and is comma separated. x values range from 0 to 1 and y
    values range from 0 to 1.

    The function is sampled on the specified image size, assuming a square image. It will be a rotationally symmetrical
    function originating in the center of the image. If oversampling is provided the image will be sampled onto an
    image of size: oversampling * image_size. Which is afterwards cropped (to the center) a number of times equal to
    oversampling. This is required when generating a simulation that is coarse grained, because then the DQE and MTF
    will be larger due consideration of a binned image.

    @param detector: eg. 'K2SUMMIT', 'FALCONII'
    @type  detector: L{str}
    @param response_function: eg. 'DQE' or 'MTF'
    @type  response_function: L{str}
    @param image_size: size of the image to sample the function on, equal x and y dimension
    @type  image_size: L{int}
    @param voltage: voltage of electron beam in eV, default 300e3
    @type  voltage: L{float}
    @param oversampling: number of times function is oversampled, multiple of 1
    @type  oversampling: L{int}
    @param folder: Folder where csv file with detector functions are stored. If not provided assume they are present
    in the current directory. todo add a folder to pytom program where some standard MTF and DQE functions are provided
    @type  folder: L{str}
    @param display: flag for displaying detector function to plot window
    @type  display: L{bool}

    @return: the detector response function, 2d array of floats
    @rtype:  L{np.ndarray}

    @author: Marten Chaillet
    """
    name = f'{detector}_{response_function}_{int(voltage*1E-3)}kV'
    if folder == '':
        filename = f'{name}.csv'
    else:
        filename = os.path.join(folder, f'{name}.csv')
    print(f'Determining {response_function} for {detector}')
    # data is a function of spatial frequency
    qdata, ydata = read_detector_csv(filename)
    params = fit_sinc_square(qdata,ydata)

    sampling_image_size = image_size * oversampling
    # fraction of nyquist maximum
    # Ny = 1
    # R, Y = xp.meshgrid(xp.arange(-Ny, Ny, 2. * Ny / (shape[0])), xp.arange(-Ny, Ny, 2. * Ny / (shape[1])))
    # r = xp.sqrt(R ** 2 + Y ** 2)
    # r = fourier_array(sampling_image_size, 2, 1)

    grids = fourier_grids((sampling_image_size,)*2, 1)  # nyquist is 1, as the fraction of nyquist maximum
    r = xp.sqrt(sum([d**2 for d in grids]))

    detector_response = sinc_square(r, params[0], params[1], params[2], params[3])

    if oversampling > 1:
        # crop the detector function
        cut = sampling_image_size // 2 - image_size // 2
        detector_response = detector_response[cut:-cut, cut:-cut]

    if display:
        display_microscope_function(detector_response, form=response_function, ylim=(0,1), complex=False)

    return detector_response


def transmission_function(sliced_potential, voltage, dz):
    """
    Calculate the transmission function from the sliced potential. The sliced potential is the simulation sample but
    averaged in z dimension per size of the multislice step. Returns:
    M{exp(i * S{sigma}_transfer * sliced_potential * S{delta}f}

    @param sliced_potential: sample averaged in z per step size, 3d array  of floats or complex values
    @type  sliced_potential: L{np.ndarray}
    @param voltage: electron beam voltage in eV, needed for calculating wavelength
    @type  voltage: L{float}
    @param dz: defocus value in m, dz > 0 is defocus, dz < 0  is overfocus
    @type  dz: L{float}

    @return: transmission function, 3d array of complex values
    @rtype:  L{np.ndarray}

    @author: Marten Chaillet
    """
    # wavelength
    Lambda = physics.wavelength_eV2m(voltage)
    # relative mass
    relative_mass = physics.constants["me"] + physics.constants["el"] * voltage / (physics.constants["c"] ** 2)
    # sigma_transfer
    sigma_transfer = 2 * xp.pi * relative_mass * physics.constants["el"] * Lambda / (physics.constants["h"] ** 2)

    return xp.exp(1j * sigma_transfer * sliced_potential * dz)


def fresnel_propagator(image_size, pixel_size, voltage, dz):
    """
    The fresnel propagator describing propagation of the electron wave through each slice of the sample.

    @param image_size: x, y dimension size of square image
    @type  image_size: L{int}
    @param pixel_size: pixel size of image in m
    @type  pixel_size: L{float}
    @param voltage: voltage of electron bream in eV
    @type  voltage: L{float}
    @param dz: defocus value in m, dz > 0 is defocus, dz < 0  is overfocus
    @type  dz: L{float}

    @return: the fresnel propagator function, 2d array of complex values
    @rtype:  L{np.ndarray}

    @author: Marten Chaillet
    """
    Lambda = physics.wavelength_eV2m(voltage)

    nyquist = 1 / (2 * pixel_size)
    grids = fourier_grids((image_size,) * 2, nyquist)
    k = xp.sqrt(sum([d ** 2 for d in grids]))
    # k = fourier_array(image_size, 2, nyquist)

    return xp.exp(-1j * xp.pi * Lambda * (k ** 2) * dz)


def create_ctf_1d(size, spacing, defocus, amplitude_contrast=0.07, voltage=300e3, Cs=2.7e-3, phase_shift_deg=.0,
                  bfactor=.0):
    """
    Create a 1 dimensional ctf curve.
    Based on tom_deconv by Tegunov.

    @param size: number of sampling points
    @type  size: L{int}
    @param spacing: spacing of sampling points in m
    @type  spacing: L{float}
    @param defocus: defocus of CTF in m
    @type  defocus: L{float}
    @param amplitude_contrast: fraction amplitude contrast in CTF
    @type  amplitude_contrast: L{float}
    @param voltage: acceleration voltage in eV
    @type  voltage: L{float}
    @param Cs: spherical aberration
    @type  Cs: L{float}
    @param phase_shift_deg: phaseshift CTF in radians
    @type  phase_shift_deg: L{float}
    @param bfactor: b_factor
    @type  bfactor: L{float}

    @return: 1d numpy array of floats
    @rtype: L{np.ndarray}

    @author: Marten Chaillet
    """
    nyquist = 1 / (2 * spacing)
    k = xp.arange(0, nyquist, nyquist / size)
    k2 = k ** 2
    lmbd = physics.wavelength_eV2m(voltage)

    phase_shift_rad = xp.deg2rad(phase_shift_deg)

    chi = xp.pi * (lmbd * defocus * k2 - 0.5 * Cs * lmbd ** 3 * k2 ** 2) + phase_shift_rad

    acurve = - amplitude_contrast * xp.cos(chi)
    pcurve = - xp.sqrt(1. - amplitude_contrast ** 2) * xp.sin(chi)

    # chi = xp.pi * (0.5 * Cs * lmbd ** 3 * k2 ** 2 + lmbd * defocus * k2) #- phase_shift_rad
    #
    # acurve = - amplitude_contrast * xp.cos(chi)
    # pcurve = xp.sqrt(1. - amplitude_contrast ** 2) * xp.sin(chi)

    bfactor = xp.exp(-bfactor * k2 * 0.25)

    return (pcurve + acurve) * bfactor


def create_ctf(shape, spacing, defocusU, amplitude_contrast, voltage, Cs, sigma_decay=.0,
               display=False, defocusV=None, defocus_angle_deg=.0, phase_shift_deg=.0, zero_cut=-1, precalculated=None):
    """
    This function models a non-complex CTF. It can be used for both 2d or 3d function. It describes a ctf after
    detection (and is therefore not complex).

    This function uses definition given in CTFFIND4, only doing one sin() operation on the grid to reduce
    computational cost.

    @param shape: shape tuple with 2 or 3 elements
    @type  shape: L{tuple} -> (L{int},)
    @param spacing: pixel/voxel spacing in m
    @type  spacing: L{float}
    @param defocusU: defocus value in m, dz > 0 is defocus, dz < 0  is overfocus
    @type  defocusU: L{float}
    @param amplitude_contrast: Amplitude contrast fraction (e.g., 0.1)
    @type  amplitude_contrast: L{float}
    @param voltage: acceleration voltage in eV, e.g. 300e3
    @type  voltage: L{float}
    @param Cs: spherical aberration in m, e.g. 2.7e-3
    @type  Cs: L{float}
    @param sigma_decay: sigma of Gaussian decay function, e.g. 0.4
    @type  sigma_decay: L{float}
    @param display: flag for plotting a radially averaged version of the CTF curve
    @type  display: L{bool}
    @param defocusV: optional astigmatism defocus along defocus angle in m
    @type  defocusV: L{float}
    @param defocus_angle_deg: angle of astigmatism in degrees
    @type  defocus_angle_deg: L{float}
    @param phase_shift_deg: phase shift of ctf function in degrees
    @type  phase_shift_deg: L{float}
    @param precalc: tuple of 3 elements of precalculated fourier grids, k, k**2, k**4, and an angular grid
    @type  precalc: L{tuple} -> (L{np.ndarray},) * 3

    @return: real-valued CTF in Fourier space, 2d or 3d array
    @rtype:  L{np.ndarray}

    @author: Marten Chaillet, Gijs van der Schot
    """
    assert 1 < len(shape) <= 3, 'Invalid shape for create_ctf, only 2d or 3d. For 1d ctf use create_ctf_1d.'

    # get the fourier space coordinate grids and angular grid
    if precalculated is None:
        nyquist = 1 / (2 * spacing)
        k, k2, k4, angles = precalculate(shape, nyquist)
    else:
        k, k2, k4, angles = precalculated

    lmbd = physics.wavelength_eV2m(voltage)

    if defocusV is None:
        # in this case we just need k2 and k4 with one defocus value
        chi = xp.pi * lmbd * defocusU * k2 - 0.5 * xp.pi * Cs * lmbd ** 3 * k4
    else:
        # in this case we need k2, k4 and a defocus grid for the astigmatism
        defocus_grid = astigmatised_defocus_grid(angles, defocusU, defocusV, xp.deg2rad(defocus_angle_deg))
        chi = xp.pi * lmbd * defocus_grid * k2 - 0.5 * xp.pi * Cs * lmbd ** 3 * k4

    # convert phase shift degrees to radians
    phase_shift_rad = xp.deg2rad(phase_shift_deg)
    tan_term = xp.arctan(amplitude_contrast / xp.sqrt(1 - amplitude_contrast ** 2))  # tangent term can also be
    # precomputed
    chi += (phase_shift_rad + tan_term)
    # determine the ctf
    ctf = - xp.sin(chi)
    # ctf = - xp.sqrt(1. - amplitude_contrast ** 2) * xp.sin(chi + phase_shift_rad) - \
    #       amplitude_contrast * xp.cos(chi + phase_shift_rad)

    if sigma_decay > 0:
        decay = xp.exp(-(k / (sigma_decay * nyquist)) ** 2)
        ctf *= decay

    if zero_cut >= 0:  # this part is only needed for template generation
        # one dimensional ctf function will be easier to find the specified zero crossing
        f_prechi = lambda k: xp.pi * lmbd * defocusU * k ** 2 - 0.5 * xp.pi * Cs * lmbd ** 3 * k ** 4
        f_chi = lambda k: - xp.sin(f_prechi(k) + phase_shift_rad + tan_term)
        # generate oversampled range of fourier space frequencies
        k_range = xp.arange(max(k.shape)) / max(k.shape) * nyquist
        # calculate ctf for each k
        values = f_chi(k_range)
        # get the indexes of the zero crossings (the places where the sign changes
        zero_crossings = xp.where(xp.diff(xp.sign(values)))[0]
        # find cutoff frequency
        k_cutoff = k_range[zero_crossings[zero_cut]] if defocusU > 0 else k_range[zero_crossings[zero_cut + 1]]
        # filter the ctf with the cutoff frequency
        flt = (k <= k_cutoff) * 1
        ctf *= flt

    if display:
        if len(shape) == 2:
            display_microscope_function(ctf, form='ctf', complex=False)
        else:
            display_microscope_function(ctf[:, :, shape[2]//2], form='ctf', complex=False)

    return ctf


# todo| In function below image shape should have equal length in each dim, in that case it makes more sense to input
# todo| an image size. But These function could be generated for non square images as well.

def create_simple_complex_ctf(image_shape, pixel_size, defocus, voltage=300E3, Cs=2.7E-3, sigma_decay=0.4,
                              display=False):
    """
    Create a complex valued contrast transfer function of the phase modulation in a 2d array. Dimensions of input image
    shape should be equal.

    @param image_shape: x and y should be same size
    @type  image_shape: L{tuple} -> (L{int},) * 2
    @param pixel_size: pixel size in m
    @type  pixel_size: L{float}
    @param defocus: defocus value in m, dz > 0 is defocus, dz < 0  is overfocus
    @type  defocus: L{float}
    @param voltage: electron beam voltage in eV, default 300e3
    @type  voltage: L{float}
    @param Cs: spherical aberration in m
    @type  Cs: L{float}
    @param sigma_decay: sigma of Gaussian decay function
    @type  sigma_decay: L{float}
    @param display: flag for plotting a radially averaged version of the CTF curve
    @type  display: L{bool}

    @return: phase ctf curve, 2d array of complex values
    @rtype:  L{np.ndarray}

    @author: Marten Chaillet
    """
    assert len(image_shape) == 2, print('image shape should be a tuple of length 2')
    assert len(set(image_shape)) == 1, print('invalid input image/volume for create CTF, dimensions need to be equal.')

    lmbd = physics.wavelength_eV2m(voltage)

    nyquist = 1 / (2 * pixel_size)

    # k = fourier_array(image_shape[0], len(image_shape), nyquist)
    grids = fourier_grids(image_shape, nyquist)
    k = xp.sqrt(sum([d ** 2 for d in grids]))

    complex_ctf = xp.exp(-1j * xp.pi / 2 * (Cs * (lmbd ** 3) * (k ** 4) - 2 * defocus * lmbd * (k ** 2)))

    # Decay function
    decay = 1 if sigma_decay <= 0 else xp.exp(-(k / (sigma_decay * nyquist)) ** 2)
    complex_ctf *= decay

    # visual inspection of CTF
    if display:
        display_microscope_function(complex_ctf, form='ctf', complex=True)

    return complex_ctf


def create_complex_ctf(image_shape, pixel_size, defocus, voltage=300E3, Cs=2.7E-3, Cc=2.7E-3,
                           energy_spread=0.7, illumination_aperture=0.030E-3, objective_diameter=100E-6,
                           focus_length=4.7E-3, astigmatism=0.0, astigmatism_angle=0.0, display=False):
    """
    Create complex valued CTF curve of phase modulation in a 2d array. Adapated from Vulovic et al., 2013.

        # default parameters for extended CTF function (InSilicoTEM, Vulovic, 2013))
    voltage                 = 300E3
    spherical aberration    = 2.7E-3
    chromatic_abberation    = 2.7E-3 # C_c
    energy_spread           = 0.7 # deltaE
    illumination_aperture   = 0.030E-3 # a_i
    objective_diameter      = 100E-6 #
    focus_length            = 4.7E-3 # focal distance
    astigmatism in 0.0E-9
    astigmatism angle in degrees

    @param image_shape: tuple of image shape with equal x and y dimension
    @type  image_shape: L{tuple} -> (L{int},) * 2
    @param pixel_size: pixel size in m
    @type  pixel_size: L{float}
    @param defocus: defocus value in m, dz > 0 is defocus, dz < 0  is overfocus
    @type  defocus: L{float}
    @param voltage: electron beam voltage in eV, default 300e3
    @type  voltage: L{float}
    @param Cs: spherical aberration in m
    @type  Cs: L{float}
    @param Cc: chromatic aberration in m
    @type  Cc: L{float}
    @param energy_spread: spread of electron beam
    @type  energy_spread: L{float}
    @param illumination_aperture: size of aperture in m
    @type  illumination_aperture: L{float}
    @param objective_diameter: diameter of objective lens in m
    @type  objective_diameter: L{float}
    @param focus_length: focal distance of objective lens in m
    @type  focus_length: L{float}
    @param astigmatism: astigmatism in m
    @type  astigmatism: L{float}
    @param astigmatism_angle: angle of astigmatism in degrees
    @type  astigmatism_angle: L{float}
    @param display: flag for displaying function
    @type  display: L{bool}

    @return: fourier space ctf function, 2d array complex valued
    @rtype:  L{np.ndarray}

    @author: Marten Chaillet
    """
    from scipy.ndimage import gaussian_filter

    # print(locals())

    assert len(image_shape) == 2, print('image shape should be a tuple of length 2')
    assert len(set(image_shape)) == 1, print('invalid input image/volume for create CTF, dimensions need to be equal.')

    lmbd = physics.wavelength_eV2m(voltage)

    inratioqqs = xp.sqrt((defocus + astigmatism) / defocus)
    inratioqlq = xp.sqrt(defocus / (defocus - astigmatism))

    # generate fourier space grids with astigmatism
    nyquist = 1 / (2 * pixel_size)
    xx, yy = fourier_grids(image_shape, nyquist)
    qsym = xp.sqrt(xx ** 2 + yy ** 2)
    astigmatism_angle_rad = xp.deg2rad(astigmatism_angle)
    xdot = xx * xp.cos(astigmatism_angle_rad) - yy * xp.sin(astigmatism_angle_rad)
    ydot = xx * xp.sin(astigmatism_angle_rad) + yy * xp.cos(astigmatism_angle_rad)
    q = xp.sqrt((xdot / inratioqlq) ** 2 + (ydot * inratioqqs) ** 2)

    # calculate chi and the CTF
    chi = 0.5 * xp.pi * (Cs * (lmbd ** 3) * (qsym ** 4) - 2 * defocus * lmbd * (q ** 2))
    complex_ctf = xp.cos(chi) - 1j * xp.sin(chi)

    # chromatic envelope
    h = Cc * energy_spread / voltage
    nominator = xp.pi * lmbd * q ** 2 * h
    denominator = 4 * xp.sqrt(xp.log(2))
    chromatic_envelope = xp.exp(- (nominator / denominator) ** 2)
    # spatial envelope
    nums = (xp.pi * Cs * lmbd ** 2 * q ** 3 - xp.pi * defocus * q) ** 2 * illumination_aperture ** 2
    spatial_envelope = xp.exp(- nums / xp.log(2))
    # full envelope
    envelope = chromatic_envelope * spatial_envelope

    # aperture function
    aperture = xp.ones(image_shape)
    qmax = 2 * xp.pi * objective_diameter / (lmbd * focus_length)
    aperture[q > qmax] = 0
    gaussian_filter(aperture, sigma=3, output=aperture)

    # convolute CTF with envelope and aperture
    complex_ctf *= (envelope * aperture)

    if display:
        display_microscope_function(complex_ctf, form='ctf', complex=True)

    return complex_ctf


def test_defocus_grid_time(shape, ntests=100):
    import time

    # Initialize parameters
    df1 = 4e-6
    df2 = 2e-6
    xx, yy = fourier_grids(shape, 1/(2*2.62e-10))
    astigmatism_angle_deg = 30
    astigmatism_angle_rad = xp.deg2rad(astigmatism_angle_deg)
    defocus, astigmatism = convert_defocusU_defocusV_to_defocus_astigmatism(df1, df2)
    inratioqqs = xp.sqrt((defocus + astigmatism) / defocus)
    inratioqlq = xp.sqrt(defocus / (defocus - astigmatism))
    k2 = xx**2 + yy**2
    angles = angular_grid(xx, yy)

    # Test grid coordinates astigmatism without cos over full grid
    t0 = time.time()

    for i in range(ntests):
        xdot = xx * xp.cos(astigmatism_angle_rad) - yy * xp.sin(astigmatism_angle_rad)
        ydot = xx * xp.sin(astigmatism_angle_rad) + yy * xp.cos(astigmatism_angle_rad)
        q = xp.sqrt((xdot / inratioqlq) ** 2 + (ydot * inratioqqs) ** 2) ** 2 * defocus

    t1 = time.time()

    print((t1 - t0) / ntests)

    # Test astigmatism speed with cos over grid but less additional operations
    t0 = time.time()

    for i in range(ntests):
        defocus_grid = astigmatised_defocus_grid(angles, df1, df2, astigmatism_angle_rad)
        q2 = defocus_grid * k2

    t1 = time.time()

    print((t1 - t0) / ntests)

    # output for: test_defocus_grid_time((2000, 3000), ntests=100)
    # >>> 0.2865477156639099
    # >>> 0.2980231332778931
    # speed seems roughly identical for both methods, so the cos method for me is preferred as it is has better
    # correspondence to the mathematical definition of the CTF curve
