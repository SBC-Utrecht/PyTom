#!/usr/bin/env pytom
import matplotlib
matplotlib.use('Qt5Agg')
import numpy
from pylab import *

import os
import lxml.etree as et


def remove_element(el):
    parent = el.getparent()
    if el.tail.strip():
        prev = el.getprevious()
        if prev:
            prev.tail = (prev.tail or '') + el.tail
        else:
            parent.text = (parent.text or '') + el.tail
    parent.remove(el)

def plotTMResults(xmlfiles, show_image=True, outname='', labels=[]):
    
    if type(xmlfiles) == type(''):
        xmlfiles = [xmlfiles]
    print(labels)
        
    fig,ax = subplots(1,1,figsize=(10,6))
    ax.set_xlabel('Particle ID')
    ax.set_ylabel('Correlation Coefficient')
    S = []
    for n, xmlfile in enumerate(xmlfiles):
        if not xmlfile: continue
        tree = et.parse(xmlfile)
        scores = []
        for particle in tree.xpath("Particle"):
            scores.append( float( particle.xpath('Score')[0].get('Value')))

        if len(labels) > n:
            ax.plot(scores,label=labels[n])
        else:
            ax.plot(scores,label=os.path.basename(xmlfile))

        S.append(scores)

    try:
        if numpy.array(S[0][:30]).sum() > numpy.array(S[1][:30]).sum():
            comp = 0
        else:
            comp = 1

        if len(S) == 2:
            for n, sc in enumerate(S[comp]):
                if sc-S[1-comp][n] < -0.0001:
                    ax.vlines(n, 0, sc, label=f'cutoff: {numpy.around(sc,3)}')
                    ax.hlines(sc, 0, n, label=f'num particles: {n+1}')
                    break
    except: pass
    ax.legend()
    fig.tight_layout()
    if outname: 
        savefig(outname)
    if show_image: 
        show()
