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
    ax.legend()
    fig.tight_layout()
    if outname: 
        savefig(outname)
    if show_image: 
        show()
