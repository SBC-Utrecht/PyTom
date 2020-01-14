#!/usr/bin/env pytom

'''
Created on Jun 27, 2011

@author: yuxiangchen
'''

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
        print helper
        sys.exit()
    try:
        pl_filename, num_steps, peak_index, ccc_value, imageFile, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print e
        sys.exit()
    if help is True:
        print helper
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
    for i in xrange(num_steps):
        x.append(min+i*step)
    x.append(max)
    
    y = []
    for i in xrange(num_steps):
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
    print 'sigma: %f, mu: %f, a: %f' % (sigma, mu, a)
    
    # plot the Gaussian fitting
    from math import exp
    gaussian_fnc = lambda x: a*exp(-(x-mu)**2/(2*sigma**2))
    plt.plot(x[1:],map(gaussian_fnc, x[1:]),'g--')
    
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
    print 'One sigma position: %f, number of estimation: %f' % (mu-sigma, estimate)
    
    estimate = 0.
    for i in x:
        if i > mu-2*sigma:
            estimate += gaussian_fnc(i)
        else:
            break
    print 'Two sigma position: %f, number of estimation: %f' % (mu-2*sigma, estimate)
    
    if ccc_value:
        ccc_value = float(ccc_value)
        estimate = 0.
        for i in x:
            if i > ccc_value:
                estimate += gaussian_fnc(i)
            else:
                break
        print 'CCC value position: %f, number of estimation: %f' % (ccc_value, estimate)
    
    if imageFile is None:      
        pyplot.show()
    else:
        if not ('png' in imageFile or 'PNG' in imageFile):
            imageFile = imageFile + '.png'
            
        pyplot.savefig(imageFile)