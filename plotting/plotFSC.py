#!/usr/bin/env pytom
import matplotlib
matplotlib.use('Qt5Agg')
import numpy
from pylab import *
import os
from pytom.gui.guiFunctions import loadstar

def plot_FSC(FSCFile, pixelsize, boxsize=0, outname='',show_image=True, c=0.143):
    data = numpy.array(list(map(float, loadstar(FSCFile))))
    dim  = len(data)
    if not boxsize: 
        boxsize = dim
        print('Boxsize set to dimension of FSC data ({})'.format(boxsize))
    x = numpy.arange(dim)*(1/(float(pixelsize)*int(boxsize)))

    fig, ax = subplots(1,1,figsize=(3,5))

    ax.plot( x, data,color='#1989ac',lw=2, label=os.path.basename(FSCFile))

    ax.set_xticks([0,x[dim//4],x[dim//2],x[3*dim//4],x[-1]])
    ax.set_xticklabels(['',1/x[dim//4]//1,1/x[dim//2]//1,1/x[3*dim//4]//1,1/x[-1]//1])
    ax.set_ylabel('Correlation Coefficient' )
    ax.set_xlabel(r'resolution ($\AA$)') 
    if c:
        ax.hlines(c,0,x[-1],lw=1,colors='#54d777', label='Cut off')
    ax.legend()
    ax.set_yticks([ 0, 0.25, 0.5, 0.75, 1])
    ax.set_yticklabels([0,0.25,0.5,0.75,1])
    fig.tight_layout()
    if outname: 
        savefig(outname)
    if show_image: 
        show()


if __name__ == '__main__':
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options

    options = [ScriptOption(['-f','--FCSFile'], 'File with FSC values in single column', True, False),
               ScriptOption(['-b','--boxSize'], 'Box size of average', True, False),
               ScriptOption(['-p','--pixelSize'],'Pixel size of voxels in model.',True, False),
               ScriptOption(['-o','--outputName'],'Save figure under this name.',True, True),
               ScriptOption(['-s','--show'],'Show figure.',False, True),
               ScriptOption(['-c','--FSCCutoff'],'Cut-off used to determine the resolution of your object from the FSC curve. Typical values are 0.5 or 0.143.', True, False),
               ScriptOption(['-h', '--help'], 'Help.', False, True)]


    helper = ScriptHelper(sys.argv[0].split('/')[-1], description='Convert coordinate list to particle list.',
                          authors='Gijs van der Schot', options=options)
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        FSCFile, boxsize, pixelsize, outname, show_image, cutoff, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()

    if help is True:
        print(helper)
        sys.exit()

    if not boxsize: boxsize=0

    if not os.path.exists(FSCFile):
        print('FSC File does not exist')
        sys.exit()

    if not outname:
        show_image = True

    if not cutoff: cutoff = 0
    else: cutoff = float(cutoff)

    if 1:
        plot_FSC(FSCFile, pixelsize, boxsize, show_image=show_image, outname=outname, c=cutoff)
    else:
        print('Error while parsing FSCFile')

