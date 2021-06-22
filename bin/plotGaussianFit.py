#!/usr/bin/env pytom

"""
Created on Jun 27, 2011
Updated on Jun 21, 2021

@author: yuxiangchen, Marten Chaillet
"""

# plotting
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

# main functionality
import sys
import numpy as np
from scipy.optimize import curve_fit
from pytom.localization.structures import readParticleFile


if __name__ == '__main__':
    # parse command line arguments
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Do the Gaussian fitting on the found particle list.',
                          authors='Yuxiang Chen',
                          options=[ScriptOption(['-f','--file'], 'Particle list after extracting candidates.', True, False),
                                   ScriptOption(['-n','--numberBins'], 'Number of bins of histogram. Default is 10.', True, True),
                                   ScriptOption(['-p','--gaussianPeak'], 'Histogram index of estimated particle '
                                                                         'population peak.', True, True),
                                   ScriptOption(['-c','--numberParticles'], 'Number of particles up to CCC value.', True, True),
                                   ScriptOption(['-i','--imageFile'], 'Save plot to a image file.', True, True),
                                   ScriptOption(['-h', '--help'], 'Help.', False, True)])
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        pl_filename, num_bins, peak_index, ccc_value, imageFile, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()
    if help is True:
        print(helper)
        sys.exit()

    # todo add scripthelper 2 parsing

    # process the arguments
    # parse number of bins
    num_bins = int(num_bins)
    if not num_bins:
        num_bins = 10
    # parse peak index
    if peak_index is not None and not(int(peak_index) > num_bins):
        peak_index = int(peak_index)
    else:
        print('fall back to default peak index, nbins / 2')
        peak_index = num_bins // 2

    
    # read out the scores
    foundParticles = readParticleFile(pl_filename)
    scores = []
    for f in foundParticles:
        scores.append(float(f.score.getValue()))

    # generate the histogram
    matplotlib.rc('font', size=18)
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(121)
    # ax.hist(x[1:],y,'ro-')
    y, x_hist, _ = ax.hist(scores, bins=num_bins, histtype='step') #, 'ro-')
    ax.set_xlabel('Score')
    ax.set_xlim(x_hist[0], x_hist[-1])
    ax.set_ylabel('Frequency')

    # define gaussian function with parameters to fit
    def gauss(x, mu, sigma, A):
        return A * np.exp(-(x - mu) ** 2 / 2 / sigma ** 2)

    def gauss_integral(sigma, A):
        # mu does not influence the integral
        return A * np.sqrt(2 * np.pi * sigma ** 2)

    def exponential(x, a, b, c):
        return a * np.exp(-b * x) + c

    # define bimodal function of two gaussians to fit both populations
    def bimodal(x, mu1, sigma1, A1, mu2, sigma2, A2):
        return gauss(x, mu1, sigma1, A1) + gauss(x, mu2, sigma2, A2)
    # def bimodal(x, a, b, c, mu, sigma, A):
    #     return exponential(x, a, b, c) + gauss(x, mu, sigma, A)

    try:
        # adjust x to center of each bin so len(x)==len(y)
        x = (x_hist[1:] + x_hist[:-1]) / 2

        # expected values
        # left gaussian expectation: mu < x[0] and A > y[0]
        # right gaussian expectation: mu ~ x[half] and A ~ y[half]
        # initial guess for mu_1 should be 0 for better convergence
        expected = (0, .1, y[0], x[peak_index], .1, y[peak_index])
        # expected = (1., 1., 0.0, x[peak_index], .1, y[peak_index])
        bounds = (0, [0.3, 0.2, np.inf, 1.0, 0.2, 200])  # todo 200 is not a good upper bound for second gaussian
        # fit function
        params_names = ['mu_1', 'sigma_1', 'A_1', 'mu_2', 'sigma_2', 'A_2']
        # params_names = ['a', 'b', 'c', 'mu_2', 'sigma_2', 'A_2']
        params, cov = curve_fit(bimodal, x, y, p0=expected, bounds=bounds) #, maxfev=2000)
        # give sigma of fit for each parameter
        sigma = np.sqrt(np.diag(cov))

        # plot bimodal model and the gaussian particle population
        plt.plot(x, bimodal(x, *params), color='blue', lw=2, label='bimodal model')
        population, noise = (params[:3], params[3:6]) if params[0] > params[3] else (params[3:6], params[:3])
        # population = params[2:5]
        plt.plot(x, gauss(x, *population), color='red', lw=2, label='particle population')
        plt.legend()

        # print information about fit of the model
        print('\nfit of the bimodal model:')
        print('\testimated\t\tsigma')
        for n, p, s in zip(params_names, params, sigma):
            print(f'{n}\t{p:.3f}\t\t{s:.3f}')
        print('\n')

        # todo add proposed cutoff based on function overlap?

        # print the estimation of number of true positives
        # x = np.flip(x)
        # gaussian_fnc = lambda x: gauss(x, *population)
        # mu, sigma = population[0], abs(population[1])
        # estimate = 0.
        # for i in x:
        #     if i > mu - sigma:
        #         estimate += gaussian_fnc(i)
        #     else:
        #         break
        # print('One sigma position: %f, particle number estimation: %f' % (mu - sigma, estimate))
        #
        # estimate = 0.
        # for i in x:
        #     if i > mu - 2 * sigma:
        #         estimate += gaussian_fnc(i)
        #     else:
        #         break
        # print( 'Two sigma position: %f, particle number estimation: %f' % (mu - 2 * sigma, estimate))

        # if ccc_value:
        #     ccc_value = float(ccc_value)
        #     estimate = 0.
        #     for i in x:
        #         if i > ccc_value:
        #             estimate += gaussian_fnc(i)
        #         else:
        #             break
        #     print('CCC value position: %f, number of estimation: %f' % (ccc_value, estimate))

        # Generate a ROC curve
        roc_steps = 100
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

        # find integral of gaussian particle population; divide by function step size
        population_integral = gauss_integral(population[1], population[2]) / hist_step
        # print('population integral: ', population_integral)

        for x_i in x_roc:
            # calculate probability of true positives and false positives for x_i
            gauss_pop = lambda x: gauss(x, *population)
            gauss_noise = lambda x: gauss(x, *noise)
            # add the number of true and false positives
            n_true_positives += gauss_pop(x_i) / delta
            n_false_positives += gauss_noise(x_i) / delta
            # add probability
            recall.append(n_true_positives / population_integral)
            fdr.append(n_false_positives / (n_true_positives + n_false_positives))
            # print(n_false_positives / (n_true_positives + n_false_positives), '\t', n_true_positives / population_integral)

        # plot the fdr curve
        ax2 = fig.add_subplot(122)
        ax2.scatter(fdr, recall, facecolors='none', edgecolors='r', s=16)
        ax2.set_xlabel('False discovery rate')
        ax2.set_ylabel('Recall')
        ax2.set_xlim(0, 1)
        ax2.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
        ax2.set_ylim(0, 1)
        ax2.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])

    except RuntimeError as e:
        # runtime error is because the model could not be fit, in that case print error and continue with execution
        print(e)

    if imageFile is None:
        plt.tight_layout()
        plt.show()
    else:
        if not ('png' in imageFile or 'PNG' in imageFile):
            imageFile = imageFile + '.png'

        plt.tight_layout()
        plt.savefig(imageFile)
