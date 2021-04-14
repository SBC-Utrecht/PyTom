#!/usr/bin/env pytom

import sys
from pytom.gui.mrcOperations import read_mrc, downsample
from pytom_volume import read
from pytom_numpy import vol2npy
import matplotlib
matplotlib.use('Qt5Agg')

from pylab import *
import copy



if __name__=='__main__':


    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Run a localization job. Documentation is available at\n\
                          http://www.pytom.org/doc/pytom/localization.html',
                          authors='Yuxiang Chen, Thomas Hrabe',
                          options=[ScriptOption(['-f','--filename'], 'Specify job.xml filename', arg=True, optional=False),
                                   ScriptOption(['-p','--powerSpectrum'], 'Specify job.xml filename', arg=False, optional=True),
                                   ScriptOption(['-s','--sliceid'], 'Parts you want to split the volume in X dimension', arg=True, optional=True),
                                   ScriptOption(['-b','--factorDownsample'], 'factor for downsampling', arg=True, optional=True),
                                   ScriptOption(['-d','--sliceDirection'], 'slice direction. Options are: x,y,z with z being default.', arg=True, optional=True),
                                   ScriptOption(['-a','--phases'], 'Show phases of the FT of the image', arg=False, optional=True),
                                   ScriptOption(['-h', '--help'], 'Help.', False, True)])
    
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    
    try:
        fname, ps, sliceId, f_ds, sliceDir, phases, b_help = parse_script_options(sys.argv[1:], helper)
    except:
        print(helper)
        sys.exit()

    if b_help:
        print(helper)
        sys.exit()


    if fname.endswith('.em'):
        vol = read(fname)
        data = copy.deepcopy(vol2npy(vol)).T
    else:
        import mrcfile
        data = copy.deepcopy( mrcfile.open(fname,permissive=True).data)

    print( 'max: {} min: {} mean: {} std: {}'.format( data.max(), data.min(), data.mean(), data.std()))
    print(data.shape)
    from scipy.ndimage import gaussian_filter

    print(data.shape)
    
    if sliceDir is None: sliceDir = 'z'
    else: sliceDir = sliceDir.lower()

    if not sliceDir in ['x', 'y', 'z']: 
        sliceDir = 'z'
        print('Slice direction not known. Set to default (z, last axis)')

    if not (sliceId is None):
        if sliceDir == 'x': 
            data = data[int(sliceId),:,:]
        if sliceDir == 'y': 
            data = data[:,int(sliceId),:]
        if sliceDir == 'z': 
            data = data[:,:,int(sliceId)]
        
    if not(f_ds is None):
        data = downsample(data,int(f_ds))
    #dd[dd>dd.mean()+5*dd.std()] = dd.mean()

    if ps or phases:

        ft = fftshift(fftn(data))
        if phases:
            data = angle(ft)
        else:
            data = log10(abs(ft)+0.0001)

    fig,ax = subplots(1,1,figsize=(10,10))
    #ax.imshow( log10(fftshift(abs(fftn(dd)))), cmap='gray')
    ax.imshow(data, cmap='gray')
    ax.set_yticks([])
    ax.set_xticks([])
    fig.tight_layout()
    show()
    #data.close()

