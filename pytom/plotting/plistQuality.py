#!/usr/bin/env pytom

"""
Created on Jun 27, 2011

@author: Marten Chaillet, Yuxiang Chen
"""


import os
import numpy as np

# Do not use 'try' around import, this file should only be called if there is a graphical backend
# Or plot could only be written without displaying to screen, I guess this requires another backend.
try:
    import matplotlib
    import matplotlib.pyplot as plt
except:
    import matplotlib
    matplotlib.use('Qt5Agg')
    import matplotlib.pyplot as plt

# TODO evaluation of particle list might also be calculated without plotting


class ScoreHistogramFdr():
    def __init__(self):
        matplotlib.rc('font', size=18)
        self.fig = plt.figure(figsize=(5*2, 5))
        self.hist_ax = self.fig.add_subplot(1, 2, 1)
        self.fdr_ax = self.fig.add_subplot(1, 2, 2)

    def draw_histogram(self, scores, nbins=30, return_bins=False):
        y, x_hist, _ = self.hist_ax.hist(scores, bins=nbins, histtype='step')
        self.hist_ax.set_xlabel('Score')
        self.hist_ax.set_xlim(x_hist[0], x_hist[-1])
        self.hist_ax.set_ylabel('Frequency')
        if return_bins:
            return y, x_hist

    def draw_bimodal(self, x, y1, y2, ymax=None):
        # plot bimodal model and the gaussian particle population
        self.hist_ax.plot(x, y1, color='blue', lw=2, label='bimodal model')
        # population = params[2:5]
        self.hist_ax.plot(x, y2, color='red', lw=2, label='particle population')
        if ymax is not None:
            self.hist_ax.set_ylim(0, ymax)
        self.hist_ax.legend()

    def draw_score_threshold(self, x, ymax):
        self.hist_ax.vlines(x, 0, ymax, linestyle='dashed', label=f'{x:.2f}')
        self.hist_ax.legend()

    def draw_fdr_recall(self, fdr, recall, optimal_id):
        self.fdr_ax.scatter(fdr, recall, facecolors='none', edgecolors='r', s=25)
        # add optimal threshold in green
        self.fdr_ax.scatter(fdr[optimal_id], recall[optimal_id], s=25, color='green')
        self.fdr_ax.plot([0, 1], [0, 1], ls="--", c=".3", lw=1)
        self.fdr_ax.set_xlabel('False discovery rate')
        self.fdr_ax.set_ylabel('Recall')
        self.fdr_ax.set_xlim(0, 1)
        self.fdr_ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
        self.fdr_ax.set_ylim(0, 1)
        self.fdr_ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])

    def draw_fdr_recall_fit(self, fdr, recall, AUC):
        self.fdr_ax.plot(fdr, recall, label='spline')
        self.fdr_ax.plot([], [], ' ', label=f'AUC: {AUC:.2f}')
        self.fdr_ax.legend()

    def write(self, filename, quality=200, transparency=False, bbox='tight'):
        plt.tight_layout()
        plt.savefig(filename, dpi=quality, transparent=transparency, bbox_inches=bbox)

    def display(self):
        plt.tight_layout()
        plt.show()


def check_square_fdr(fdr, recall, epsilon=1e-3):
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


def calculate_histogram(scores, num_steps):
    # construct x and y array according to the given peak index
    # preferably input is already sorted
    scores.sort()  # this is sorted from lowest to highest
    min = scores[0]
    max = scores[-1]

    step = (max - min) / num_steps
    x = []
    for i in range(num_steps):
        x.append(min + i * step)
    x.append(max)

    y = []
    for i in range(num_steps):
        lower = x[i];
        upper = x[i + 1]
        n = len([v for v in scores if lower <= v <= upper])
        y.append(n)

    return x, y


def evaluate_estimates(estimated_positions, ground_truth_positions, tolerance):
    """
    Estimated_positions numpy array, ground truth positions numpy array
    :param estimated_positions:
    :type estimated_positions:
    :param ground_truth_positions:
    :type ground_truth_positions:
    :param tolerance:
    :type tolerance:
    :return:
    :rtype:
    """
    from scipy.spatial.distance import cdist
    n_estimates = estimated_positions.shape[0]
    matrix = cdist(estimated_positions, ground_truth_positions, metric='euclidean')
    correct = [0] * n_estimates
    for i in range(n_estimates):
        if matrix[i].min() < tolerance:
            correct[i] = 1
    return correct


def fdr_recall(correct_particles, scores):
    assert all(i > j for i, j in zip(scores, scores[1:])), print('Scores list should be decreasing.')

    n_true_positives = sum(correct_particles)
    true_positives, false_positives = 0, 0
    fdr, recall = [], []
    for correct, score in zip(correct_particles, scores):
        if correct:
            true_positives += 1
        else:
            false_positives += 1
        if n_true_positives == 0:
            recall.append(0)
        else:
            recall.append(true_positives / n_true_positives)
        fdr.append(false_positives / (true_positives + false_positives))
    return fdr, recall


def get_distance(line, point):

    a1, b1 = line
    x, y = point
    a2 = - (1 / a1)
    b2 = y - a2 * x

    x_int = (b2 - b1) / (a1 - a2)
    y_int = a2 * x_int + b2

    return np.sqrt((x_int - x)**2 + (y_int - y)**2)


def distance_to_random(fdr, recall):
    auc = [0] * len(fdr)
    for i in range(len(fdr)):
        d = get_distance((1, 0), (fdr[i], recall[i]))  # AUC should be 1 at most not 1/2
        if recall[i] > fdr[i]:
            auc[i] = d
        else:
            auc[i] = -d
    return max(auc), np.argmax(auc)


def plist_quality_ground_truth(particle_list, ground_truth, position_tolerance=5, output_figure_name=None, nbins=30):
    # read out the scores, they should be in descending order
    correlation_scores = [-1] * len(particle_list)
    estimated_positions = []
    for i, f in enumerate(particle_list):
        correlation_scores[i] = float(f.getScoreValue())
        estimated_positions.append(f.getPickPosition().toVector())
    correlation_scores.sort(reverse=True)
    correlation_scores = np.array(correlation_scores)
    estimated_positions = np.array(estimated_positions)

    correct = evaluate_estimates(estimated_positions, ground_truth, position_tolerance)
    fdr, recall = fdr_recall(correct, correlation_scores)

    quality, index = distance_to_random(fdr, recall)
    print('Optimal score threshold: ', correlation_scores[index])

    if output_figure_name is not None:
        plot = ScoreHistogramFdr()
        y, _ = plot.draw_histogram(correlation_scores, nbins=nbins, return_bins=True)
        plot.draw_score_threshold(correlation_scores[index], max(y))
        plot.draw_fdr_recall(fdr, recall, index)
        plot.write(output_figure_name, quality=200, transparency=False, bbox='tight')

    return quality / np.cos(np.deg2rad(45))  # divide by max distance to line to normalize between -1 and 1


def plist_quality_gaussian_fit(particle_list, particle_peak_index, force_peak=False,
                               output_figure_name=None, crop_hist=False, num_bins=30):
    # TODO make sure ValueError are printed with output if data cannot be fit
    import scipy.interpolate as interpolate
    import traceback
    from scipy.optimize import curve_fit
    from scipy.special import erf

    # read out the scores
    correlation_scores = [-1] * len(particle_list)
    for i, f in enumerate(particle_list):
        correlation_scores[i] = float(f.getScoreValue())
    correlation_scores.sort(reverse=True)
    correlation_scores = np.array(correlation_scores)

    # draw the histogram
    plot = ScoreHistogramFdr()
    y, x_hist = plot.draw_histogram(correlation_scores, nbins=num_bins, return_bins=True)

    # ==================================== functions for fitting =======================================================
    # define gaussian function with parameters to fit
    def gauss(x, mu, sigma, A):
        return A * np.exp(-(x - mu) ** 2 / 2 / sigma ** 2)

    # integral of gaussian with certain sigma and A
    def gauss_integral(sigma, A):
        # mu does not influence the integral
        return A * np.abs(sigma) * np.sqrt(2 * np.pi)

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
        expected = (0, .1, y[0], x[particle_peak_index], .1, y[particle_peak_index])
        # force peak of particle population to be at peak index
        if force_peak:
            bounds = ([0, 0, 0, x[particle_peak_index] - 0.01, 0, 0],
                      [.1, 0.2, np.inf, x[particle_peak_index] + 0.01, 0.2, 200])
        else:
            bounds = ([0, 0, 0, 0.1, 0, 0],
                      [.1, 0.2, np.inf, 1.0, 0.2, 200])  # todo 200 is not a good upper bound for gaussian 2

        # TODO try fitting, if it does not work dont continue with ROC analysis
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

        y_bimodal = bimodal(x, *params)
        population, noise = (params[:3], params[3:6]) if params[0] > params[3] else (params[3:6], params[:3])
        y_gauss = gauss(x, *population)

        if crop_hist:
            plot.draw_bimodal(x, y_bimodal, y_gauss, ymax=3 * population[2])
        else:
            plot.draw_bimodal(x, y_bimodal, y_gauss)

        # ======================================= Generate a ROC curve =================================================
        roc_steps = 50
        x_roc = np.flip(np.linspace(x[0], x[-1], roc_steps))

        # find ratio of hist step vs roc step
        hist_step = (x_hist[-1] - x_hist[0]) / num_bins
        roc_step = (x[-1] - x[0]) / roc_steps
        delta = hist_step / roc_step  # can be used to divide true/false positives by per roc step
        # variable for total number of tp and fp
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

        for x_i in x_roc:
            # calculate probability of true positives x_i
            # n_true_positives += gauss_pop(x_i) / delta
            n_true_positives = (1 - CDF(x_i)) * population_integral

            # determine false positives up to this point, could also use CDF
            n_false_positives += gauss_noise(x_i) / delta

            # add probability
            recall.append(n_true_positives / population_integral)
            fdr.append(n_false_positives / (n_true_positives + n_false_positives))

        skip = 0
        is_increasing = all([a <= b for a, b in zip(fdr, fdr[1:])])
        while (is_increasing == False) and ((skip + 1) < len(fdr)):
            skip += 1
            is_increasing = all([a <= b for a, b in zip(fdr[skip:], fdr[skip+1:])])
        recall = recall[skip:]
        fdr = fdr[skip:]
        x_roc = x_roc[skip:]

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
        plot.draw_score_threshold(x_roc[id], max(y))
        print(f' - optimal correlation coefficient threshold is {x_roc[id]}.')
        print(f' - this threshold approximately selects {x_roc[id] * population_integral} particles.')

        # points for plotting
        xs = np.linspace(0, 1, 200)
        # plot the fdr curve
        plot.draw_fdr_recall(fdr, recall, id)

        if check_square_fdr(fdr, recall):
            AUC = 1.0
            print('Area Under Curve (AUC): ', AUC)
            # in this case just plot a connecting line between all points
            plot.draw_fdr_recall_fit(fdr, recall, AUC)
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
            plot.draw_fdr_recall_fit(xs, spline(xs), AUC)

    except (RuntimeError, ValueError) as e:
        # runtime error is because the model could not be fit, in that case print error and continue with execution
        traceback.print_exc()

    if output_figure_name is None:
        plot.display()
    else:
        if not ('png' in output_figure_name or 'PNG' in output_figure_name):
            output_figure_name = output_figure_name + '.png'
        plot.write(output_figure_name, )


def read_ground_truth(filename, zoffset):
    from pytom.basic.files import loadtxt
    from pytom.basic.datatypes import SIMULATED_GROUND_TRUTH
    ground_truth = loadtxt(filename, dtype=SIMULATED_GROUND_TRUTH)
    names = ground_truth['ParticleName']
    coordinates = np.vstack((ground_truth['x'], ground_truth['y'], ground_truth['z'])).T
    coordinates[:, 2] += zoffset
    rots = np.vstack((ground_truth['ThetaZ'], ground_truth['PhiX'], ground_truth['PsiZ'])).T
    return names, coordinates, rots


def select_ground_truth(names, coordinates, select):
    return coordinates[names == select]


def evaluate_multiple_lists_gt(particles, particle_list_files, ground_truth_file, ground_truth_zoffset, output='.',
                               nbins=30):
    from pytom.basic.structures import ParticleList
    ground_truth_names, ground_truth_locations, ground_truth_rotations = read_ground_truth(ground_truth_file,
                                                                                           ground_truth_zoffset)
    quality_scores = []
    for name, file in zip(particles, particle_list_files):
        particle_list = ParticleList()
        particle_list.fromXMLFile(file)
        particle_ground_truth = select_ground_truth(ground_truth_names, ground_truth_locations, name)
        quality_scores.append(plist_quality_ground_truth(particle_list, particle_ground_truth,
                                            output_figure_name=os.path.join(output,
                                                                            os.path.splitext(file)[0] + '.png'),
                                            nbins=nbins))
    return quality_scores


def write_quality(results, file_name):
    with open(file_name, 'w') as fstream:
        for r in results:
            fstream.write("{:s} {:.4f}".format(r[0], r[1]))
            fstream.write('\n')


def parse_linker_file_w_gt(input_file, input_name, gt_file, gt_zoffset, output_location, num_bins):
    linker_data = np.genfromtxt(input_file, dtype=[('particle', 'U1000'), ('particle_list_file', 'U1000')],
                                skip_header=1)
    particle_ids = linker_data['particle']
    plist_files = [os.path.join(os.path.split(input_file)[0], f) for f in linker_data['particle_list_file']]
    scores = evaluate_multiple_lists_gt(particle_ids, plist_files, gt_file, gt_zoffset,
                                        output=output_location,
                                        nbins=num_bins)
    results = [(pid, score) for pid, score in zip(particle_ids, scores)]
    write_quality(results, os.path.join(output_location, input_name + '_scores.txt'))


# OLD PLOT GAUSSIAN FIT CODE
# if __name__ == '__main__':
#     # parse command line arguments
#     import sys
#     from pytom.tools.script_helper import ScriptHelper2, ScriptOption2
#     from pytom.tools.parse_script_options import parse_script_options2
#     helper = ScriptHelper2(
#         sys.argv[0].split('/')[-1],  # script name
#         description='Evaluate particle list quality.',
#         authors='Marten Chaillet, Yuxiang Chen',
#         options=[ScriptOption2(['-f', '--file'], 'Particle list with Template Matching candidates or text file that '
#                                                  'links particle ids from ground truth file to particle lists.', 'file',
#                                'optional'),
#                  ScriptOption2(['-o', '--output'], 'Filename of output or directory for output.', 'string',
#                                'optional', '.'),
#                  ScriptOption2(['-g', '--ground_truth'], 'File with simulated ground truth data.', 'file', 'optional'),
#                  ScriptOption2(['--ground_truth_zoffset'], 'Z offset for ground truth locations.', 'int',
#                                'optional', 0),
#                  ScriptOption2(['-n', '--nbins'], 'Number of bins for histogram.', 'int', 'optional', 30)])
#
#     options = parse_script_options2(sys.argv[1:], helper)
#
#     input_file, output_location, gt_file, gt_zoffset, number_of_bins = options
#
#     # read ground truth if available
#     input_name = os.path.splitext(os.path.split(input_file)[1])[0]
#     input_ext = os.path.splitext(input_file)[1]
#     if input_ext == '.txt':
#         linker_data = np.genfromtxt(input_file, dtype=[('particle', 'U1000'), ('particle_list_file', 'U1000')],
#                                     skip_header=1)
#         particle_ids = linker_data['particle']
#         plist_files = [os.path.join(os.path.split(input_file)[0], f) for f in linker_data['particle_list_file']]
#         scores = evaluate_multiple_lists_gt(particle_ids, plist_files, gt_file, gt_zoffset,
#                                             output=output_location,
#                                             nbins=number_of_bins)
#         results = [(pid, score) for pid, score in zip(particle_ids, scores)]
#         write_quality(results, os.path.join(output_location, input_name + '_scores.txt'))
#     elif input_ext == '.xml':
#         # evaluate single list with gaussian fit
#         # read out the scores
#         from pytom.localization.structures import readParticleFile
#         scores = []
#         num_steps = number_of_bins
#         foundParticles = readParticleFile(input_file)
#         for f in foundParticles:
#             scores.append(float(f.score.getValue()))
#
#
#
#         # plot
#         from matplotlib import pyplot
#         import matplotlib
#         import numpy
#
#         matplotlib.rc('lines', linewidth=2)
#         matplotlib.rc('font', size=24)
#         fig = pyplot.figure()
#         plt = fig.add_subplot(111)
#         plt.plot(x[1:], y, 'ro-')
#         plt.set_xlabel('Score')
#         plt.set_ylabel('Frequency')
#
#         # do the fitting
#         from pytom.tools.maths import gaussian_fit
#
#         sigma, mu, a = gaussian_fit(x[peak_index:], y[peak_index - 1:])
#
#         if sigma.__class__ in [complex, numpy.complex128]:
#             sigma = sigma.real
#         if mu.__class__ in [complex, numpy.complex128]:
#             mu = mu.real
#         print('sigma: %f, mu: %f, a: %f' % (sigma, mu, a))
#
#         # plot the Gaussian fitting
#         from math import exp
#
#         gaussian_fnc = lambda x: a * exp(-(x - mu) ** 2 / (2 * sigma ** 2))
#         plt.plot(x[1:], list(map(gaussian_fnc, x[1:])), 'g--')
#
#         # print the estimation of number of true positives
#         x.reverse()
#         while True:
#             if x[-1] > 0:
#                 x.append(x[-1] - step)
#             else:
#                 x = x[:-1]
#                 break
#
#         estimate = 0.
#         for i in x:
#             if i > mu - sigma:
#                 estimate += gaussian_fnc(i)
#             else:
#                 break
#         print('One sigma position: %f, number of estimation: %f' % (mu - sigma, estimate))
#
#         estimate = 0.
#         for i in x:
#             if i > mu - 2 * sigma:
#                 estimate += gaussian_fnc(i)
#             else:
#                 break
#         print('Two sigma position: %f, number of estimation: %f' % (mu - 2 * sigma, estimate))
#
#         if ccc_value:
#             ccc_value = float(ccc_value)
#             estimate = 0.
#             for i in x:
#                 if i > ccc_value:
#                     estimate += gaussian_fnc(i)
#                 else:
#                     break
#             print('CCC value position: %f, number of estimation: %f' % (ccc_value, estimate))
#
#         if imageFile is None:
#             pyplot.show()
#         else:
#             if not ('png' in imageFile or 'PNG' in imageFile):
#                 imageFile = imageFile + '.png'
#
#             pyplot.savefig(imageFile)
#     else:
#         print('Invalid input file format for Particle List Evalulation.')
#         sys.exit(0)
