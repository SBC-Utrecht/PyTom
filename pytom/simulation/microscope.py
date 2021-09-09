import os
import numpy as xp
import pytom.simulation.physics as physics
from scipy.optimize import curve_fit


def fourier_array(shape, nyquist):
    """
    Generate a fourier space frequency array where values range from -nyquist to +nyquist, with the center equal
    to zero.

    @param shape: tuple with 2 or 3 elements with the size of each dimension
    @type  shape: L{tuple} -> (L{int},)
    @param nyquist: nyquist frequency giving the range of frequencies
    @type  nyquist: L{float}

    @return: fourier space frequencies, 2d or 3d array of floats
    @rtype:  L{numpy.ndarray}

    @author: Marten Chaillet
    """
    assert type(shape) == tuple, "shape needs to be a tuple of dimenion sizes, dimenions need to be equal"
    assert len(set(shape)) == 1, "dimensions are not identical"
    assert len(shape) <= 3, "too many dimensions, 3 is max"

    size = shape[0]

    d = xp.arange(-nyquist, nyquist, 2. * nyquist / size)

    if len(shape) == 2:
        grids = xp.meshgrid(d, d)
    else:
        grids = xp.meshgrid(d, d, d)

    # Wave vector and direction
    k = xp.sqrt(sum([d**2 for d in grids]))  # xz ** 2 + yz ** 2 + zz ** 2)
    # k2 = k ** 2
    # k4 = k2 ** 2
    return k  # , k2, k4


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


def display_microscope_function(image, form='', complex=False):
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
    matplotlib.use('Qt5Agg')
    import matplotlib.pyplot as plt

    if complex:
        r1, m1 = radial_average(image.real)
        r2, m2 = radial_average(image.imag)
    else:
        r1, m1 = radial_average(image)

    fig, ax = plt.subplots()
    ax.plot(r1, m1, label=form)
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
    r = fourier_array((sampling_image_size,)*2, 1)  # nyquist is 1, as the fraction of nyquist maximum

    detector_response = sinc_square(r, params[0], params[1], params[2], params[3])

    if oversampling > 1:
        # crop the detector function
        cut = sampling_image_size // 2 - image_size // 2
        detector_response = detector_response[cut:-cut, cut:-cut]

    if display:
        display_microscope_function(detector_response, form=response_function, complex=False)

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
    k = fourier_array((image_size,)*2, nyquist)

    return xp.exp(-1j * xp.pi * Lambda * (k ** 2) * dz)


def create_ctf(shape, spacing, defocus, amplitude_contrast, voltage, Cs, sigma_decay=0.4,
               display=False):
    """
    This function models a non-complex CTF. It can be used for both 2d or 3d function. It describes a ctf after
    detection (and is therefore not complex).

    @param shape: shape tuple with 2 or 3 elements
    @type  shape: L{tuple} -> (L{int},)
    @param spacing: pixel/voxel spacing in m
    @type  spacing: L{float}
    @param defocus: defocus value in m, dz > 0 is defocus, dz < 0  is overfocus
    @type  defocus: L{float}
    @param amplitude_contrast: Amplitude contrast fraction (e.g., 0.1)
    @type  amplitude_contrast: L{float}
    @param voltage: acceleration voltage in eV
    @type  voltage: L{float}
    @param Cs: spherical aberration in m
    @type  Cs: L{float}
    @param sigma_decay: sigma of Gaussian decay function
    @type  sigma_decay: L{float}
    @param display: flag for plotting a radially averaged version of the CTF curve
    @type  display: L{bool}

    @return: real-valued CTF in Fourier space, 2d or 3d array
    @rtype:  L{np.ndarray}

    @author: Marten Chaillet
    """
    nyquist = 1 / (2 * spacing)
    k = fourier_array(shape, nyquist)
    lmbd = physics.wavelength_eV2m(voltage)

    chi = xp.pi * lmbd * defocus * k**2 - 0.5 * xp.pi * Cs * lmbd ** 3 * k ** 4

    ctf = - xp.sqrt(1. - amplitude_contrast ** 2) * xp.sin(chi) - amplitude_contrast * xp.cos(chi)

    decay = 1 if sigma_decay <= 0 else xp.exp(-(k / (sigma_decay * nyquist)) ** 2)
    ctf *= decay

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

    k = fourier_array(image_shape, nyquist)

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

    assert len(image_shape) == 2, print('image shape should be a tuple of length 2')
    assert len(set(image_shape)) == 1, print('invalid input image/volume for create CTF, dimensions need to be equal.')

    image_size = image_shape[0]

    lmbd = physics.wavelength_eV2m(voltage)

    q_true = 1 / (image_size * pixel_size)

    dfsh = defocus + astigmatism
    dfln = defocus - astigmatism
    inratioqqs = xp.sqrt(dfsh/defocus)
    inratioqlq = xp.sqrt(defocus/dfln)

    x = xp.arange(-(image_size//2), abs(-image_size//2), 1)
    xx, yy = xp.meshgrid(x, x)
    xdot = xx * xp.cos(astigmatism_angle*xp.pi/180) - yy * xp.sin(astigmatism_angle*xp.pi/180)
    ydot = xx * xp.sin(astigmatism_angle*xp.pi/180) - yy * xp.cos(astigmatism_angle*xp.pi/180)
    q = xp.sqrt((xdot / inratioqlq) ** 2 + (ydot * inratioqqs) ** 2) * q_true
    qsym = xp.sqrt(xdot ** 2 + ydot ** 2) * q_true

    # print(r)
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

    complex_ctf *= (envelope * aperture)

    # print(complex_CTF.shape, K.shape, A.shape)
    # visual inspection of CTF
    if display:
        display_microscope_function(complex_ctf, form='ctf', complex=True)

    return complex_ctf
