#!/usr/bin/env pytom
import matplotlib

try:
    matplotlib.use('Qt5Agg')
except:
    pass
from pylab import subplots, savefig, show
import os

def plot_FSC(FSCFile, pixelsize, boxsize=0, outname='', show_image=True, c=0.143, directory='', resolution=None, rand=True):
    import numpy, os
    from pytom.gui.guiFunctions import loadstar

    if directory:
        FSCFiles = sorted([os.path.join(directory, fsc) for fsc in os.listdir(directory) if
                           fsc.endswith('.dat') and not 'Filter' in fsc])
        if not FSCFiles: return

    elif type(FSCFile) == str:
        FSCFiles = [FSCFile]
    else:
        FSCFiles = FSCFile

    FSCFile = FSCFiles[0]

    data = numpy.array(list(map(float, loadstar(FSCFile))))
    dim = len(data)
    if not boxsize:
        boxsize = dim
        print('Boxsize set to dimension of FSC data ({})'.format(boxsize))
    x1 = numpy.arange(1,dim+1) * (1 / (float(pixelsize) * int(boxsize)))
    x = numpy.arange(data.shape[0])

    fig, ax = subplots(1, 1, figsize=(10, 5))
    newax = ax.twiny()
    data = numpy.array(list(map(float, loadstar(FSCFile))))
    newax.plot(x+1, data, lw=2)
    for FSCFile in FSCFiles:
        if not rand and ('FSCCorr.dat' in FSCFile or 'FSCRand.dat' in FSCFile):
            continue
        data = numpy.array(list(map(float, loadstar(FSCFile))))
        ax.plot(x, data, lw=2, label=os.path.basename(FSCFile))

    #newax.set_xticks([0, x[dim // 4], x[dim // 2], x[3 * dim // 4], x[-1]])
    #newax.set_xticklabels([0, x[dim // 4] // 1, x[dim // 2] // 1, x[3 * dim // 4] // 1, x[-1] // 1])
    ax.set_ylabel('Correlation Coefficient')
    ax.set_xlabel(r'resolution ($\AA$)')
    if c:
        ax.hlines(c, 0, x[-1], lw=1, colors='#54d777', label='Cut off')
    ax.legend()
    #ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    #ax.set_yticklabels([0, 0.25, 0.5, 0.75, 1])

    ax.set_xticks([0, x[dim // 4], x[dim // 2], x[3 * dim // 4], x[-1]])
    ax.set_xticklabels([f'{(1/x1[0]):.2f}', f'{1 / x1[dim // 4 ]:.2f}', f'{1 / x1[dim // 2]:.2f}', f'{1 / x1[3 * dim // 4]:.2f}', f'{1 / x1[-1]:.2f}'])

    if not resolution is None:
        ax.set_title(f'resolution: {float(resolution):8.3f} A')

    fig.tight_layout()
    if outname:
        savefig(outname)
    if show_image:
        show()


if __name__ == '__main__':
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options

    options = [ScriptOption(['-f', '--FCSFile'], 'File with FSC values in single column', True, True),
               ScriptOption(['-d', '--directory'], 'File with FSC values in single column', True, True),
               ScriptOption(['-b', '--boxSize'], 'Box size of average', True, False),
               ScriptOption(['-p', '--pixelSize'], 'Pixel size of voxels in model.', True, False),
               ScriptOption(['-o', '--outputName'], 'Save figure under this name.', True, True),
               ScriptOption(['-s', '--show'], 'Show figure.', False, True),
               ScriptOption(['-c', '--FSCCutoff'],
                            'Cut-off used to determine the resolution of your object from the FSC curve. Typical values are 0.5 or 0.143.',
                            True, False),
               ScriptOption(['-h', '--help'], 'Help.', False, True)]

    helper = ScriptHelper(sys.argv[0].split('/')[-1], description='Convert coordinate list to particle list.',
                          authors='Gijs van der Schot', options=options)
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        FSCFile, directory, boxsize, pixelsize, outname, show_image, cutoff, help = parse_script_options(sys.argv[1:],
                                                                                                         helper)
    except Exception as e:
        print(e)
        sys.exit()

    if help is True:
        print(helper)
        sys.exit()

    if not boxsize: boxsize = 0

    if not os.path.exists(str(FSCFile)) and not os.path.exists(str(directory)):
        print('FSC File(s) do not exist')
        sys.exit()

    if not outname:
        show_image = True

    if not cutoff:
        cutoff = 0
    else:
        cutoff = float(cutoff)

    directory = '' if directory is None else directory

    if 1:
        plot_FSC(FSCFile, pixelsize, boxsize, show_image=show_image, outname=outname, c=cutoff, directory=directory)
    else:
        print('Error while parsing FSCFile')

