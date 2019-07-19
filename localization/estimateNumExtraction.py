#!/usr/bin/env python

'''
Created on Sep 22, 2010

@author: chen
'''

from pytom.tools.maths import gaussian_fit

def usage():
    print './scriptname -c classified_result_file -t target_particle_name [-s expected_cover_rate_in_sigma (1,2 or 3)] [-m estimate_positive_class_mean_value]'

if __name__ == '__main__':
    import sys, getopt
    
    if len(sys.argv) ==1:
        print "No argument is given!"
        usage()
        sys.exit()        
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hc:t:s:m:", ["help"])
    except getopt.GetoptError:
        print 'Command not right. Exit!'
        sys.exit()
    
    filename = ''
    target_particle = ''
    cover_rate=1
    mean_value=None
    
    for o,a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o in ("-c"):
            filename = a
        if o in ("-t"):
            target_particle = a
        if o in ("-s"):
            cover_rate = int(a)
            if cover_rate not in [1,2,3]:
                print 'Please choose sigma between 1, 2 or 3'
                sys.exit()
        if o in ("-m"):
            mean_value = float(a)
    
    from pytom.localization.structures import readParticleFile
    particles = readParticleFile(filename)
    
    pos_score=[]
    neg_score=[]
    for p in particles:
        if p.classified == target_particle:
            pos_score.append(p.getScore())
        else:
            neg_score.append(p.getScore())
    
    print 'Positive num: %d' % len(pos_score)
    print 'Negative num: %d' % len(neg_score)
    
    if not mean_value:
        mean_value = pos_score[len(pos_score)/2]
    
    num = 10
    total_score = pos_score + neg_score
    total_score.sort()
    min = total_score[0]
    max = total_score[-1]
    
    step = (max-min)/num
    x = []
    for i in xrange(num):
        x.append(min+i*step)
    x.append(max)
    
    y_neg = []
    y_pos = []
    y_total = []
    for i in xrange(num):
        lower = x[i]; upper = x[i+1]
        n1 = len([v for v in neg_score if lower<=v<=upper])
        y_neg.append(n1)
        n2 = len([v for v in pos_score if lower<=v<=upper])
        y_pos.append(n2)
        y_total.append(n1+n2)
    
    # gaussian fit
    l = len([i for i in x if i >= mean_value])
    if 2*l >= len(x):
        xnew = x
        ynew = y_pos[-l:]; ynew.reverse(); ynew = ynew + y_pos[-l:]
        ynew = ynew[-len(x):]
    else:
        xnew = x[-2*l:]
        ynew = y_pos[-l:]; ynew.reverse();
        ynew = ynew + y_pos[-l:]
    
    sigma_pos, mu_pos, a_pos = gaussian_fit(x[1:],y_pos)
    sigma_neg, mu_neg, a_neg = gaussian_fit(x[1:],y_neg)
    sigma_new, mu_new, a_new = gaussian_fit(xnew,ynew)
    
    from math import exp
    y_pos_g = [a_pos*exp(-(v-mu_pos)**2/(2*sigma_pos**2)) for v in x]
    y_neg_g = [a_neg*exp(-(v-mu_neg)**2/(2*sigma_neg**2)) for v in x]
    y_new_g = [a_new*exp(-(v-mu_new)**2/(2*sigma_new**2)) for v in xnew]
    
    lower_score = mu_new-cover_rate*sigma_new
    print 'The estimated lower band of the score is: %.3f' % lower_score
    
    x_est=[]
    s = x[-1]
    while s > lower_score:
        x_est.append(s)
        s = s-step
    x_est.append(lower_score)
    n_est_pos = [a_pos*exp(-(v-mu_pos)**2/(2*sigma_pos**2)) for v in x_est]
    n_est_neg = [a_neg*exp(-(v-mu_neg)**2/(2*sigma_neg**2)) for v in x_est]
    from math import ceil
    n = ceil(sum(n_est_pos))+ceil(sum(n_est_neg))
    print 'The estimated number of extracted peaks is: %d' % n
    
    # plot
    from matplotlib import pyplot
    fig = pyplot.figure()
    plt = fig.add_subplot(111)
    plt.plot(x[1:],y_neg,'ro-', x[1:],y_pos,'go-', x,y_pos_g,'b--', x,y_neg_g,'k--', xnew, y_new_g, 'y--')
    plt.legend(('Classified as negative', 'Classified as positive', 'Positive gaussian fit', 'Negative gaussian fit', 'Estimated correct positive'), 'upper right', shadow=True)
    plt.set_xlabel('Score')
    plt.set_ylabel('Frequency')
    pyplot.show()