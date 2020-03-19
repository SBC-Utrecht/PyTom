import numpy as xp
from scipy.optimize import curve_fit
# import matplotlib
# matplotlib.use('Qt5Agg')
# import matplotlib.pyplot as plt
# from pylab import *


def sinc_square(x, p1, p2, p3, p4):
    return p1 * (xp.sin(xp.pi*(x-p2)*p3) / (xp.pi*(x-p2)*p3))**2 + p4


def read_detector_csv(filename):
    data = xp.genfromtxt(filename,delimiter=',',dtype=xp.float32)
    x = data[:,0]
    y = data[:,1]
    return x,y


def fit_sinc_square(xdata,ydata):
    # Here you give the initial parameters for p0 which Python then iterates over
    # to find the best fit
    popt, pcov = curve_fit(sinc_square,xdata,ydata,p0=(1.0, 1.0, 1.0, 1.0)) #THESE PARAMETERS ARE USER DEFINED

    print(popt) # This contains your two best fit parameters

    # Performing sum of squares
    p1 = popt[0]
    p2 = popt[1]
    p3 = popt[2]
    p4 = popt[3]
    residuals = ydata - sinc_square(xdata,p1,p2,p3,p4)
    fres = sum(residuals**2)

    print(f'chi-square value for fitting is {fres}') #THIS IS YOUR CHI-SQUARE VALUE!

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

    # Return the sin parameters
    return p1,p2,p3,p4


def create_detector_response(detector, response_function, image, voltage=300E3, folder=''):
    """
    CSV file containing the detector curve should on the x-axis be expressed as fraction of Nyquist (between 0 and 1)
    and on the y-axis as the MTF or DQE value (also between 0 and 1).

    @param detector: eg. 'FALCONII'
    @type detector: string
    @param response_function: eg. 'DQE'
    @type response_function: string
    @param voltage: Voltage of detector operation (in V not kV)
    @type voltage: float
    @param image: Image is passed for determining the size of the response
    @type image: 2d array (numpy
    @param folder: Folder where csv file with detector functions are stored. If not provided assume they are present
    in the current directory. eg. '/data2/mchaillet/simulation/detectors'
    @type folder: string

    @return: the detector response function
    @rtype 2d array (numpy)

    @author: Marten Chaillet
    """
    name = f'{detector}_{response_function}_{int(voltage*1E-3)}kV'
    if folder == '':
        filename = f'{name}.csv'
    else:
        filename = f'{folder}/{name}.csv'

    # data is a function of spatial frequency
    qdata,ydata = read_detector_csv(filename)
    params = fit_sinc_square(qdata,ydata)

    shape = image.shape
    # fraction of nyquist maximum
    Ny = 1
    R, Y = xp.meshgrid(xp.arange(-Ny, Ny, 2. * Ny / (shape[0])), xp.arange(-Ny, Ny, 2. * Ny / (shape[1])))
    r = xp.sqrt(R ** 2 + Y ** 2)

    return sinc_square(r,params[0],params[1],params[2],params[3])


