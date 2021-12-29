#!/usr/bin/env pytom

"""
Created on Jun 27, 2011

@author: yuxiangchen
"""


import os
import numpy as np

# Do not use 'try' around import, this file should only be called if there is a graphical backend
# Or plot could only be written without displaying to screen, I guess this requires another backend.
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt


def init_plot(subplots=2):
    matplotlib.rc('font', size=18)
    fig = plt.figure(figsize=(5*subplots, 5))
    axes = []
    for i in range(subplots):
        axes.append(fig.add_subplot(1, subplots, i + 1))
    return fig, axes


def draw_histogram(ax, scores, nbins=40):
    y, x_hist, _ = ax.hist(scores, bins=nbins, histtype='step')  #, 'ro-')
    ax.set_xlabel('Score')
    ax.set_xlim(x_hist[0], x_hist[-1])
    ax.set_ylabel('Frequency')
    return y, x_hist


def draw_fdr_recall(ax, fdr, recall, optimal_id):
    ax.scatter(fdr, recall, facecolors='none', edgecolors='r', s=25)
    # add optimal threshold in green
    ax.scatter(fdr[optimal_id], recall[optimal_id], s=25, color='green')
    ax.plot([0, 1], [0, 1], ls="--", c=".3", lw=1)
    ax.set_xlabel('False discovery rate')
    ax.set_ylabel('Recall')
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])


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


def plist_quality(particle_list, ground_truth, position_tolerance=5, output_figure_name=None, nbins=30):
    # read out the scores, they should be in descending order
    correlation_scores = [-1] * len(particle_list)
    estimated_positions = []
    for i, f in enumerate(particle_list):
        correlation_scores[i] = float(f.getScoreValue())
        estimated_positions.append(f.getPickPosition().toVector())
    correlation_scores = np.array(correlation_scores)
    estimated_positions = np.array(estimated_positions)

    correct = evaluate_estimates(estimated_positions, ground_truth, position_tolerance)
    fdr, recall = fdr_recall(correct, correlation_scores)

    quality, index = distance_to_random(fdr, recall)

    if output_figure_name is not None:
        fig, axes = init_plot(2)
        _, _ = draw_histogram(axes[0], correlation_scores, nbins=nbins)
        draw_fdr_recall(axes[1], fdr, recall, index)
        plt.savefig(output_figure_name, dpi=200, transparent=True, bbox_inches='tight')

    return quality


def read_ground_truth(filename, zoffset):
    from pytom.basic.files import loadtxt
    from pytom.basic.datatypes import SIMULATED_GROUND_TRUTH
    ground_truth = loadtxt(filename, dtype=SIMULATED_GROUND_TRUTH)
    names = ground_truth['ParticleName']
    coordinates = np.vstack((ground_truth['x'], ground_truth['y'], ground_truth['z'])).T
    coordinates[:, 2] += zoffset
    rots = np.vstack((ground_truth['ThetaZ'], ground_truth['PsiX'], ground_truth['PhiZ'])).T
    return names, coordinates, rots


def select_ground_truth(names, coordinates, select):
    return coordinates[names == select]


def evaluate_template_matching(particles, particle_list_files, ground_truth_file, ground_truth_zoffset,  output='.',
                               nbins=30):
    from pytom.basic.structures import ParticleList
    ground_truth_names, ground_truth_locations, ground_truth_rotations = read_ground_truth(ground_truth_file,
                                                                                           ground_truth_zoffset)
    quality_scores = []
    for name, file in zip(particles, particle_list_files):
        particle_list = ParticleList()
        particle_list.fromXMLFile(file)
        particle_ground_truth = select_ground_truth(ground_truth_names, ground_truth_locations, name)
        quality_scores.append(plist_quality(particle_list, particle_ground_truth,
                                            output_figure_name=os.path.join(output,
                                                                            os.path.splitext(file)[0] + '.png'),
                                            nbins=nbins))
    return quality_scores


if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper2, ScriptOption2
    from pytom.tools.parse_script_options import parse_script_options2
    helper = ScriptHelper2(
        sys.argv[0].split('/')[-1],  # script name
        description='Evaluate particle list quality.',
        authors='Marten Chaillet, Yuxiang Chen',
        options=[ScriptOption2(['-f', '--file'], 'Particle list with Template Matching candidates or text file that '
                                                 'links particle ids from ground truth file to particle lists.', 'file',
                               'optional'),
                 ScriptOption2(['-o', '--output'], 'Filename of output or directory for output.', 'string',
                               'optional', '.'),
                 ScriptOption2(['-g', '--ground_truth'], 'File with simulated ground truth data.', 'file', 'optional'),
                 ScriptOption2(['--ground_truth_zoffset'], 'Z offset for ground truth locations.', 'int',
                               'optional', 0),
                 ScriptOption2(['-n', '--nbins'], 'Number of bins for histogram.', 'int', 'optional', 30)])

    options = parse_script_options2(sys.argv[1:], helper)

    input_file, output_location, gt_file, gt_zoffset, number_of_bins = options

    # read ground truth if available
    input_ext = os.path.splitext(input_file)[1]
    if input_ext == '.txt':
        linker_data = np.genfromtxt(input_file, dtype=[('particle', 'U1000'), ('particle_list_file', 'U1000')],
                                    skip_header=1)
        particle_ids = linker_data['particle']
        plist_files = [os.path.join(os.path.split(input_file)[0], f) for f in linker_data['particle_list_file']]
        scores = evaluate_template_matching(particle_ids, plist_files, gt_file, gt_zoffset,
                                            output=output_location,
                                            nbins=number_of_bins)
        print([(pid, score) for pid, score in zip(particle_ids, scores)])
    elif input_ext == '.xml':
        # evaluate single list with gaussian fit
        # read out the scores
        from pytom.localization.structures import readParticleFile
        scores = []
        num_steps = number_of_bins
        foundParticles = readParticleFile(input_file)
        for f in foundParticles:
            scores.append(float(f.score.getValue()))

        # construct x and y array according to the given peak index
        scores.sort()
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

        # plot
        from matplotlib import pyplot
        import matplotlib
        import numpy

        matplotlib.rc('lines', linewidth=2)
        matplotlib.rc('font', size=24)
        fig = pyplot.figure()
        plt = fig.add_subplot(111)
        plt.plot(x[1:], y, 'ro-')
        plt.set_xlabel('Score')
        plt.set_ylabel('Frequency')

        # do the fitting
        from pytom.tools.maths import gaussian_fit

        sigma, mu, a = gaussian_fit(x[peak_index:], y[peak_index - 1:])

        if sigma.__class__ in [complex, numpy.complex128]:
            sigma = sigma.real
        if mu.__class__ in [complex, numpy.complex128]:
            mu = mu.real
        print('sigma: %f, mu: %f, a: %f' % (sigma, mu, a))

        # plot the Gaussian fitting
        from math import exp

        gaussian_fnc = lambda x: a * exp(-(x - mu) ** 2 / (2 * sigma ** 2))
        plt.plot(x[1:], list(map(gaussian_fnc, x[1:])), 'g--')

        # print the estimation of number of true positives
        x.reverse()
        while True:
            if x[-1] > 0:
                x.append(x[-1] - step)
            else:
                x = x[:-1]
                break

        estimate = 0.
        for i in x:
            if i > mu - sigma:
                estimate += gaussian_fnc(i)
            else:
                break
        print('One sigma position: %f, number of estimation: %f' % (mu - sigma, estimate))

        estimate = 0.
        for i in x:
            if i > mu - 2 * sigma:
                estimate += gaussian_fnc(i)
            else:
                break
        print('Two sigma position: %f, number of estimation: %f' % (mu - 2 * sigma, estimate))

        if ccc_value:
            ccc_value = float(ccc_value)
            estimate = 0.
            for i in x:
                if i > ccc_value:
                    estimate += gaussian_fnc(i)
                else:
                    break
            print('CCC value position: %f, number of estimation: %f' % (ccc_value, estimate))

        if imageFile is None:
            pyplot.show()
        else:
            if not ('png' in imageFile or 'PNG' in imageFile):
                imageFile = imageFile + '.png'

            pyplot.savefig(imageFile)
    else:
        print('Invalid input file format for Particle List Evalulation.')
        sys.exit(0)

