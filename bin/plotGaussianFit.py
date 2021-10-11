#!/usr/bin/env pytom

"""
Created on Jun 27, 2011
Updated on Jun 21, 2021

@author: yuxiangchen, Marten Chaillet
"""
# todo add option for prividing a ground truth file for the particle list
# todo increase number of points on roc curve to make more robust?
# todo make sure ValueError are printed with output if data cannot be fit

# plotting
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

# main functionality
import sys
import numpy as np
import scipy.interpolate as interpolate
import traceback
from scipy.optimize import curve_fit
from scipy.special import erf
from pytom.localization.structures import readParticleFile


def check_square_FDR(fdr, recall, epsilon=1e-3):
    """
    @param fdr: list of fdr values
    @type  fdr: L{list}
    @param recall: list of recall values
    @type  recall: L{list}
    @param epsilon: tolerance for closeness to 0 and 1
    @type  epsilon: L{float}

    @return: boolean of whether the FDR is almost square within tolerance
    @rtype:  L{bool}

    @author: Marten Chaillet
    """

    # fdr and recall should contain values very close to 0 and 1, respectively, if function is square
    union = [(f, r) for f, r in zip(fdr, recall) if ((np.abs(0.0 - f) < epsilon) and (np.abs(1.0 - r) < epsilon))]

    # print(union)
    # using list for True False staments would be more pythonic, but I want to force the function to return a boolean
    #  for clarity purposes
    return True if union else False


def distance_to_diag(fdr, recall):
    """

    @param fdr: list of fdr values
    @type  fdr: L{list}
    @param recall: list of recall values
    @type  recall: L{list}

    @return: list of distance of each fdr, recall combination to diagonal line
    @rtype:  L{list}

    @author: Marten Chaillet
    """
    # two point on the diagonal to find the distance to
    lp1, lp2 = (0, 0), (1, 1)
    # list to hold distances
    distance = []
    for f, r in zip(fdr, recall):
        d = np.abs((lp2[0] - lp1[0]) * (lp1[1] - r) - (lp1[0] - f) * (lp2[1] - lp1[1])) / \
            np.sqrt((lp2[0] - lp1[0]) ** 2 + (lp2[1] - lp1[1]) ** 2)
        distance.append(d)
    return distance


if __name__ == '__main__':
    # parse command line arguments with ScriptHelper2

    from pytom.tools.script_helper import ScriptHelper2, ScriptOption2
    from pytom.tools.parse_script_options import parse_script_options2

    helper = ScriptHelper2(
        sys.argv[0].split('/')[-1],  # script name
        description='Fit a Gaussian to particle list peak.',
        authors='Marten Chaillet (based on Yuxiang Chen)',
        options=[
            ScriptOption2(['-f', '--file'], 'Particle list xml file', 'file', 'required'),
            ScriptOption2(['-n', '--numberBins'], 'Numbers of bins to plot the score distribution', 'int',
                          'optional', 10),
            ScriptOption2(['-p', '--gaussianPeak'], 'Index of the peak (p < nbins) to be taken as a starting point '
                                                    'for searching the mean of the particle population', 'int',
                          'optional'),
            ScriptOption2(['--forcePeak'], 'Force the given peak index to be the mean', 'no arguments', 'optional'),
            ScriptOption2(['--cropPlot'], 'Crops distribution on y-axis to not show full height, '
                                          'can make particle peak more clearly visible', 'no arguments', 'optional'),
            ScriptOption2(['-i', '--imageFile'], 'Save plot to an image file; if not provided plot will be shown in '
                                                 'window', 'string', 'optional'),
            ScriptOption2(['-m', '--max'], 'Maximum number to start shifting', 'float', 'optional')])

    options = parse_script_options2(sys.argv[1:], helper)

    pl_filename, num_bins, peak_index, force_peak, crop_plot, imageFile, max_roc_threshold = options

    # parse peak index
    if (peak_index is None) or (peak_index > num_bins):
        print(' - fall back to default peak index estimate, nbins / 2')
        peak_index = num_bins // 2
        if force_peak:
            print(' - WARNING: peak is forced without specifying a peak index, the default peak estimate at bin number '
                  'nbins/2 will be enforced.')

    # read out the scores
    foundParticles = readParticleFile(pl_filename)
    scores = []
    for f in foundParticles:
        scores.append(float(f.score.getValue()))

    # ================================== generate the histogram ========================================================
    matplotlib.rc('font', size=18)
    fig = plt.figure(figsize=(10, 5))
    ax1 = fig.add_subplot(121)
    # ax.hist(x[1:],y,'ro-')
    y, x_hist, _ = ax1.hist(scores, bins=num_bins, histtype='step') #, 'ro-')
    ax1.set_xlabel('Score')
    ax1.set_xlim(x_hist[0], x_hist[-1])
    ax1.set_ylabel('Frequency')

    # ==================================== functions for fitting =======================================================
    # define gaussian function with parameters to fit
    def gauss(x, mu, sigma, A):
        return A * np.exp(-(x - mu) ** 2 / 2 / sigma ** 2)

    # integral of gaussian with certain sigma and A
    def gauss_integral(sigma, A):
        # mu does not influence the integral
        return A * np.abs(sigma) * np.sqrt(2 * np.pi)

    # exponential function
    def exponential(x, a, b, c):
        return a * np.exp(-b * x) + c

    # a cumulative distribution function
    def cumulative_dist(x, Lambda):
        return 1 - np.exp(- Lambda * x)

    # define bimodal function of two gaussians to fit both populations
    def bimodal(x, mu1, sigma1, A1, mu2, sigma2, A2):
        return gauss(x, mu1, sigma1, A1) + gauss(x, mu2, sigma2, A2)

    try:
        # ================================== fit bimodal distribution ==================================================
        # adjust x to center of each bin so len(x)==len(y)
        x = (x_hist[1:] + x_hist[:-1]) / 2

        # expected values
        # left gaussian expectation: mu < x[0] and A > y[0]
        # right gaussian expectation: mu ~ x[half] and A ~ y[half]
        # initial guess for mu_1 should be 0 for better convergence
        expected = (0, .1, y[0], x[peak_index], .1, y[peak_index])
        # force peak of particle population to be at peak index
        if force_peak:
            bounds = ([0, 0, 0, x[peak_index] - 0.01, 0, 0],
                      [.1, 0.2, np.inf, x[peak_index] + 0.01, 0.2, 200])
        else:
            bounds = ([0, 0, 0, 0.1, 0, 0],
                      [.1, 0.2, np.inf, 1.0, 0.2, 200])  # todo 200 is not a good upper bound for
            # second gaussian
        # parameter names for output
        params_names = ['mu_1', 'sigma_1', 'A_1', 'mu_2', 'sigma_2', 'A_2']
        # params_names = ['a', 'b', 'c', 'mu_2', 'sigma_2', 'A_2']
        params, cov = curve_fit(bimodal, x, y, p0=expected, bounds=bounds)  # max iterations argument: maxfev=2000)
        # give sigma of fit for each parameter
        sigma = np.sqrt(np.diag(cov))

        # print information about fit of the model
        print('\nfit of the bimodal model:')
        print('\testimated\t\tsigma')
        for n, p, s in zip(params_names, params, sigma):
            print(f'{n}\t{p:.3f}\t\t{s:.3f}')
        print('\n')

        # plot bimodal model and the gaussian particle population
        ax1.plot(x, bimodal(x, *params), color='blue', lw=2, label='bimodal model')
        population, noise = (params[:3], params[3:6]) if params[0] > params[3] else (params[3:6], params[:3])
        # population = params[2:5]
        ax1.plot(x, gauss(x, *population), color='red', lw=2, label='particle population')
        if crop_plot:
            ax1.set_ylim(0, 3 * population[2])
        ax1.legend()

        # ======================================= Generate a ROC curve =================================================
        roc_steps = 50
        if max_roc_threshold:
            x_roc = np.flip(np.linspace(x[0], max_roc_threshold, roc_steps))
        else:
            x_roc = np.flip(np.linspace(x[0], x[-1], roc_steps))

        # find ratio of hist step vs roc step
        hist_step = (x_hist[-1] - x_hist[0]) / num_bins
        roc_step = (x[-1] - x[0]) / roc_steps
        delta = hist_step / roc_step  # can be used to divide true/false positives by per roc step
        # variable for total number of tp and fp
        n_true_positives = .0
        n_false_positives = .0
        # list for storing probability of true positives and false positives for each cutoff
        recall = []  # recall = TP / (TP + FN); TP + FN is the full area under the Gaussian curve
        fdr = []  # false discovery rate = FP / (TP + FP); == 1 - TP / (TP + FP)

        # find integral of gaussian particle population; NEED TO DIVIDE BY HISTOGRAM BIN STEP
        population_integral = gauss_integral(population[1], population[2]) / hist_step
        print(' - estimation total number of true positives: ', population_integral)

        # should use CDF (cumulative distribution function) of Gaussian, gives probability from -infinity to x
        CDF = lambda x: 0.5 * (1 + erf((x - population[0])/(np.sqrt(2) * population[1])))
        gauss_noise = lambda x: gauss(x, *noise)
        # alternative (less precise) way of determining true positives
        # gauss_pop = lambda x: gauss(x, *population)

        for x_i in x_roc:
            # calculate probability of true positives x_i
            # n_true_positives += gauss_pop(x_i) / delta
            n_true_positives = (1 - CDF(x_i)) * population_integral

            # determine false positives up to this point, could also use CDF
            n_false_positives += gauss_noise(x_i) / delta

            # add probability
            recall.append(n_true_positives / population_integral)
            fdr.append(n_false_positives / (n_true_positives + n_false_positives))

        # ============= attempt to fit ROC curve by cumulative distribution function of exponential distribution =======
        # # ROC curve characterized by f(x) = 1 - e^(-lx)  (l=lambda)
        # params2_names = ['lambda']
        # params2, cov2 = curve_fit(cumulative_dist, [f for f,r in zip(fdr, recall) if f >0.05],
        #                           [r for f, r in zip(fdr, recall) if f > 0.05], p0=[1])
        # # give sigma of fit for each parameter
        # sigma2 = np.sqrt(np.diag(cov2))
        #
        # # print information about fit of the model
        # # print('\nfit of the cumulative distribution:')
        # # print('\testimated\t\tsigma')
        # for n, p, s in zip(params2_names, params2, sigma2):
        #     print(f'{n}\t{p:.3f}\t\t{s:.3f}')
        # print('\n')
        # x_cum = np.linspace(0,1,200)
        # y_cum = cumulative_dist(x_cum, *params2)

        # ============================= use splines to fit ROC curve ===================================================

        # use trigonometry for finding optimal threshold based on given ROC values
        distances = distance_to_diag(fdr, recall)
        id = distances.index(max(distances))
        # plot the threshold on the distribution plot for visual inspection
        ax1.vlines(x_roc[id], 0, max(y), linestyle='dashed', label=f'{x_roc[id]}')
        print(f' - optimal correlation coefficient threshold is {x_roc[id]}.')
        print(f' - this threshold approximately selects {recall[id] * population_integral} particles.')

        # points for plotting
        xs = np.linspace(0, 1, 200)
        # plot the fdr curve
        ax2 = fig.add_subplot(122)
        ax2.scatter(fdr, recall, facecolors='none', edgecolors='r', s=25)
        # add optimal threshold in green
        ax2.scatter(fdr[id], recall[id], s=25, color='green')
        ax2.plot([0, 1], [0, 1], ls="--", c=".3", lw=1)  # transform=ax.transAxes,

        if check_square_FDR(fdr, recall):
            AUC = 1.0
            print('Area Under Curve (AUC): ', AUC)
            # in this case just plot a connecting line between all points
            ax2.plot(fdr, recall, label='spline')
        else:
            recall.append(1.)
            fdr.append(1.)
            recall.insert(0, .0)
            fdr.insert(0, .0)
            # use CubicSpline this fits the (sort of) cumulative dist best)
            spline = interpolate.CubicSpline(fdr, recall)
            # scipy CubicSpline can be easily integrated to find area under curve
            AUC = spline.integrate(0, 1)
            print(' - Area Under Curve (AUC): ', AUC)
            # plot the fitted CubicSpline
            ax2.plot(xs, spline(xs), label='spline')

        # ax2.plot(x_cum, y_cum, label='spline')
        ax2.plot([], [], ' ', label=f'AUC: {AUC:.2f}')
        ax2.legend()
        ax2.set_xlabel('False discovery rate')
        ax2.set_ylabel('Recall')
        ax2.set_xlim(0, 1)
        ax2.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
        ax2.set_ylim(0, 1)
        ax2.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])

    except (RuntimeError, ValueError) as e:
        # runtime error is because the model could not be fit, in that case print error and continue with execution
        traceback.print_exc()

    if imageFile is None:
        plt.tight_layout()
        plt.show()
    else:
        if not ('png' in imageFile or 'PNG' in imageFile):
            imageFile = imageFile + '.png'

        plt.tight_layout()
        plt.savefig(imageFile)
