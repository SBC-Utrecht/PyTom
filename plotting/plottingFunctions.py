#!/usr/bin/env pytom
import matplotlib
matplotlib.use('Qt5Agg')
import numpy
from pylab import *



import os
import lxml.etree as et

def im_show(image):
    try:
        imshow(image.get())
        show()
    except:
        imshow(image)
        show()


def remove_element(el):
    parent = el.getparent()
    if el.tail.strip():
        prev = el.getprevious()
        if prev:
            prev.tail = (prev.tail or '') + el.tail
        else:
            parent.text = (parent.text or '') + el.tail
    parent.remove(el)

def plotTMResults(xmlfiles, show_image=True, outname='', labels=[], plot_cross=False):
    
    if type(xmlfiles) == type(''):
        xmlfiles = [xmlfiles]

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
        if plot_cross:
            if numpy.array(S[0][30:]).sum() > numpy.array(S[1][30:]).sum():
                comp = 0
            else:
                comp = 1

            cid = 0

            if len(S) == 2:
                silence = False
                for n, sc in enumerate(S[comp]):
                    if n < 5: continue
                    if S[1-comp][n] < sc:
                        silence = False

                    if sc-S[1-comp][n] < -0.0001 and not silence:
                        ax.vlines(n, 0, sc, label=f'cutoff {cid}: {numpy.around(sc,3)}')
                        ax.hlines(sc, 0, n, label=f'num particles ({cid}): {n+1}')
                        silence = True
                        cid += 1
    except Exception as e:
        print(e)
        pass
    ax.legend()
    fig.tight_layout()
    if outname: 
        savefig(outname)
    if show_image: 
        show()
