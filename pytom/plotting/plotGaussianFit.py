#!/usr/bin/env pytom

"""
Created on Jun 27, 2011

@author: yuxiangchen
"""
import numpy as np


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


def fdr_recall(correct, scores):
    assert all(i > j for i, j in zip(scores, scores[1:])), print('Scores list should be decreasing.')

    n_true_positives = sum(correct)
    true_positives, false_positives = 0, 0
    fdr, recall = [], []
    for i, score in enumerate(scores):
        if correct[i]:
            true_positives += 1
        else:
            false_positives
        recall.append(true_positives / n_true_positives)
        fdr.append(false_positives / (true_positives + false_positives))
    return recall, fdr


def get_distance(line, point):

    a1, b1 = line
    x, y = point
    a2 = - (1 / a1)
    b2 = y - a2 * x

    x_int = (b2 - b1) / (a1 - a2)
    y_int = a2 * x_int + b2

    return np.sqrt((x_int - x)**2 + (y_int - y)**2)


def distance_to_random(fdr, recall):
    AUC = [0] * len(fdr)
    for i in range(len(fdr)):
        AUC[i] = get_distance((1, 0), (fdr[i], recall[i])) * 2  # AUC should be 1 at most not 1/2
    return max(AUC)


def plist_quality(particle_list, ground_truth, position_tolerance):

    # read out the scores, they should be in descending order
    scores = [-1] * len(particle_list._particleList)
    estimated_positions = []
    for i, f in enumerate(foundParticles):
        scores[i] = float(f.score.getValue())
        estimated_positions.append(f.getPickPosition().toVector())
    scores = np.array(scores)
    estimated_positions = np.array(estimated_positions).T
    ground_truth_positions = np.array(ground_truth).T

    correct = evaluate_estimates(estimated_positions, ground_truth_positions, position_tolerance)
    fdr, recall = fdr_recall(correct, scores)
    quality = distance_to_random(fdr, recall)
    return quality


if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Do the Gaussian fitting on the found particle list.',
                          authors='Yuxiang Chen',
                          options=[ScriptOption(['-f','--file'], 'Particle list after extracting candidates.', True, False),
                                   ScriptOption(['-n','--numberBins'], 'Number of bins of histogram. Default is 10.', True, True),
                                   ScriptOption(['-p','--gaussianPeak'], 'The correspondent index of the gaussian peak.', True, False),
                                   ScriptOption(['-c','--numberParticles'], 'Number of particles up to CCC value.', True, True),
                                   ScriptOption(['-i','--imageFile'], 'Save plot to a image file.', True, True),
                                   ScriptOption(['-h', '--help'], 'Help.', False, True)])
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        pl_filename, num_steps, peak_index, ccc_value, imageFile, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()
    if help is True:
        print(helper)
        sys.exit()
    
    # process the arguments
    num_steps = int(num_steps)
    if not num_steps:
        num_steps = 10
    
    peak_index = int(peak_index)
    
    scores = []
    
    # read out the scores
    from pytom.localization.structures import readParticleFile
    foundParticles = readParticleFile(pl_filename)
    for f in foundParticles:
        scores.append(float(f.score.getValue()))
    
    # construct x and y array according to the given peak index
    scores.sort()
    min = scores[0]
    max = scores[-1]
    
    step = (max-min)/num_steps
    x = []
    for i in range(num_steps):
        x.append(min+i*step)
    x.append(max)
    
    y = []
    for i in range(num_steps):
        lower = x[i]; upper = x[i+1]
        n = len([v for v in scores if lower<=v<=upper])
        y.append(n)
        
    
    # plot
    from matplotlib import pyplot
    import matplotlib
    import numpy
    matplotlib.rc('lines', linewidth=2)
    matplotlib.rc('font', size=24)
    fig = pyplot.figure()
    plt = fig.add_subplot(111)
    plt.plot(x[1:],y,'ro-')
    plt.set_xlabel('Score')
    plt.set_ylabel('Frequency')
    
    # do the fitting
    from pytom.tools.maths import gaussian_fit
    sigma, mu, a = gaussian_fit(x[peak_index:], y[peak_index-1:])
    
    if sigma.__class__ in [complex,numpy.complex128]:
        sigma = sigma.real
    if mu.__class__ in [complex,numpy.complex128]:
        mu = mu.real
    print('sigma: %f, mu: %f, a: %f' % (sigma, mu, a))
    
    # plot the Gaussian fitting
    from math import exp
    gaussian_fnc = lambda x: a*exp(-(x-mu)**2/(2*sigma**2))
    plt.plot(x[1:],list(map(gaussian_fnc, x[1:])),'g--')
    
    # print the estimation of number of true positives
    x.reverse()
    while True:
        if x[-1] > 0:
            x.append(x[-1]-step)
        else:
            x = x[:-1]
            break
    
    estimate = 0.
    for i in x:
        if i > mu-sigma:
            estimate += gaussian_fnc(i)
        else:
            break
    print('One sigma position: %f, number of estimation: %f' % (mu-sigma, estimate))
    
    estimate = 0.
    for i in x:
        if i > mu-2*sigma:
            estimate += gaussian_fnc(i)
        else:
            break
    print('Two sigma position: %f, number of estimation: %f' % (mu-2*sigma, estimate))
    
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