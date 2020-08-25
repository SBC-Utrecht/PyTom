#!/usr/bin/env pytom

'''
Created on Jul 21, 2011

@author: hrabe
'''



if __name__ == '__main__':
    # parse command line argument
    import sys, os
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.basic.correlation import FSC,determineResolution
    from pytom.basic.files import read
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Determine resolution by FSC.',
                          authors='Thomas Hrabe',
                          options=[ScriptOption('--v1', 'First volume path.', arg=True, optional=True),
                                   ScriptOption('--v2', 'Second volume path.', arg=True, optional=True),
                                   ScriptOption('--pl', 'A particleList if v1 and v2 are not available.', arg=True, optional=True),
                                   ScriptOption('--fsc', 'The FSC criterion. Value between 0.0 and 1.0. Standard values are 0.5 or 0.3', arg=True, optional=False),
                                   ScriptOption('--numberBands', 'Number of bands (optional). If not set, numberBands = cubesize/4.', arg=True, optional=True),
                                   ScriptOption(['-m','--mask'], 'Mask (optional, but recomended).', arg=True, optional=True),
                                   ScriptOption('--pixelsize', 'Pixelsize in Angstrom (optional). Will return resolution in Angstrom. ', arg=True, optional=False),
                                   ScriptOption('--xml', 'Output in XML. (optional) ', arg=False, optional=True),
                                   ScriptOption('--randomizePhases', 'Correlation threshold beyond which phases are randomized. Default None. (optional) ', arg=True, optional=True),
                                   ScriptOption('--plot', 'Plot Output. (optional) ', arg=False, optional=True),
                                   ScriptOption('--outputFolder', 'Output Folder. (optional) ', arg=True, optional=True),
                                   ScriptOption(['-v','--verbose'], 'Verbose data. (optional) ', arg=False, optional=True),
                                   ScriptOption(['-h', '--help'], 'Help.', arg=False, optional=True)])


    if len(sys.argv) <=2:
        print(helper)
        sys.exit()
    try:
        v1Filename, v2Filename, particleList, fscCriterion, numberBands, mask, pixelSize, xml, randomize, plot, outdir, verbose, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        sys.exit()
    if help is True:
        print(helper)
        sys.exit()


    try:
        fscCriterion = float(fscCriterion)
    except TypeError:
        fscCriterion = 0.5
    
    try:
        if pixelSize:
            pixelSize = float(pixelSize)
    except ValueError:
        raise ValueError('The value for pixelsize must be a float!')
    
    try:
        mask= read(mask) 
    except:
        mask = None
    
    try:
        if numberBands:
            numberBands = int(numberBands)
    except ValueError:
        raise ValueError('The value for numberBands must be a integer!')

    try:
        if randomize is None:
            randomize= 0
        else:
            randomize = float(randomize)
    except ValueError or randomize > 1 or randomize < 0:
        raise ValueError('The value of randomizePhases should be a float between 0 and 1.')

    if outdir is None:
        outdir = './'

    if v1Filename and v2Filename:    
        v1  = read(v1Filename)
        v2  = read(v2Filename)
        
        if not numberBands:
            numberBands = int(v1.sizeX()//2)
        
        f = FSC(v1,v2,numberBands,mask,verbose)
        #f = [1, 0.9999816330959427, 0.9998727543638058, 0.9986706311763134, 0.9967503089610698, 0.9945473896086455, 0.9953391452631559, 0.9926167503040506, 0.9886499997836082, 0.9846988786130074, 0.9850170987613799, 0.9849222409268831, 0.9820779872082366, 0.981161445747785, 0.978102618561267, 0.9749670311874213, 0.9710488103484851, 0.9686564997495193, 0.9643854526532541, 0.9621314730670771, 0.9568606791204185, 0.9488084842261655, 0.9416015322427522, 0.9316703960630233, 0.9106972097966776, 0.8912878863048055, 0.8795187076235272, 0.8713474813842144, 0.8508750246010772, 0.8195820950483931, 0.8065990902773463, 0.7823660922709397, 0.7521861621134768, 0.7236765089946452, 0.6967229052025852, 0.6815563485825349, 0.6606856188994277, 0.6513326589007892, 0.6340321380485363, 0.6057709663336085, 0.584793117868313, 0.5612341962238455, 0.5623548376402193, 0.5564192235463952, 0.5397585867441497, 0.5089376925831242, 0.4647300273454635, 0.4310334881504598, 0.41741429413879966, 0.4195952490948154, 0.4066609275881065, 0.37019375197265786, 0.31936001816491055, 0.275903152373691, 0.2538399661514517, 0.24197349348826183, 0.20730706794873432, 0.1823204187105925, 0.17041895753522998, 0.16153106667416953, 0.13872290716093824, 0.12428131231796732, 0.09594749780671366, 0.09281895056407187, 0.0950512406930502, 0.07845013400157819, 0.06778918720241832, 0.05699422426625805, 0.04004787096713291, 0.035285330785697615, 0.02761223114687527, 0.029854265039150632, 0.030802679570736933, 0.02855574408867763, 0.04062783248396335, 0.046982827702621556, 0.044667285930674615, 0.03327190294513204, 0.028879433147898908, 0.019113096081122542, 0.018889519864393182, 0.03363102079214279, 0.030416314115916717, 0.015045702588444513, 0.007700419599421394, 0.013662921155622407, 0.02549288977161008, 0.01648898979277964, 0.004577992397744576, 0.003687537468279412, 0.015624522796941348, 0.012150048589636583, 0.013236997547964386, 0.024818980351894827, 0.017881736272355488, 0.008875703095090339, 0.004009836930167128, 0.005169522148403328, 0.013778610598594218, 0.024255111798589142]



        if verbose: print('FSC:\n', f)

        import pytom.tompy.correlation as correlation
        from pytom_numpy import vol2npy
        import numpy as np
        from pytom.tompy.io import write
        if randomize is None:
            for (ii, fscel) in enumerate(f):
                f[ii] = 2.*fscel/(1.+fscel)
            r = determineResolution(f, fscCriterion, verbose)
        else:
            randomizationFrequency    = np.floor(determineResolution(np.array(f), randomize, verbose)[1])
            oddVolumeRandomizedPhase  = correlation.randomizePhaseBeyondFreq(vol2npy(v1), randomizationFrequency)
            evenVolumeRandomizedPhase = correlation.randomizePhaseBeyondFreq(vol2npy(v2), randomizationFrequency)
            write(os.path.join(outdir, 'randOdd.mrc'), oddVolumeRandomizedPhase)
            write(os.path.join(outdir, 'randEven.mrc'), evenVolumeRandomizedPhase)
            oddVolumeRandomizedPhase = read(os.path.join(outdir, 'randOdd.mrc'))
            evenVolumeRandomizedPhase = read(os.path.join(outdir, 'randEven.mrc'))
            fsc_rand = FSC(oddVolumeRandomizedPhase, evenVolumeRandomizedPhase, numberBands, mask, verbose)
            if verbose:
                print('FSC_Random:\n', fsc_rand)
            fsc_corr = list(correlation.calc_FSC_true(np.array(f),np.array(fsc_rand)))
            if verbose:
                print('FSC_true:\n', fsc_corr)

            #for (ii, fscel) in enumerate(fsc_corr):
            #    fsc_corr[ii] = 2.*fscel/(1.+fscel)
            r = determineResolution(fsc_corr,fscCriterion, verbose)
            #os.system('rm randOdd.em randEven.em')

    elif particleList:
        from pytom.basic.structures import ParticleList
        
        pl = ParticleList('.')
        pl.fromXMLFile(particleList)
        
        if len(pl) <= 1:
            raise RuntimeError('There is only 1 or less particles in the particle list. Need at least two! Abort!')
        
        if not numberBands:
            p = pl[0]
            pv = p.getVolume()
            numberBands = int(pv.sizeX()/4)
            
        r = pl.determineResolution(fscCriterion,numberBands,mask,verbose=verbose,plot='',keepHalfsetAverages = True,
                                   halfsetPrefix=os.path.join(outdir, 'plFSC'), randomize=randomize)
        print(f'Even and odd halfsets were written into {outdir} and stored as plFSCeven / odd .em!')
        
    else:
        raise RuntimeError('You must provide either two files or a particlelist as parameters for determining resolution!')
        
    if not xml:
        if v1Filename and v2Filename:
            print('Resolution determined for ', v1Filename, ' and ', v2Filename)
        elif particleList:
            print('Resolution determined for ', particleList, ' ')
            
        print('')
        print('FSC Criterion:   ', fscCriterion)
        print('Number of Bands: ', numberBands)
        print('')
        print('Nyquist: ', r[0])
        print('Band:    ', r[1])
        
        
        if pixelSize:
            from pytom.basic.resolution import bandToAngstrom
            
            print('Resolution determined for pixelsize : ', pixelSize , ' at ', bandToAngstrom(r[1], pixelSize, numberBands), ' Angstrom')
    
    else:
        print('XML')

    if plot and v1Filename and v2Filename:
        import matplotlib
        import numpy
        matplotlib.use('Qt5Agg')
        from pylab import subplots, savefig, show


        fig,ax = subplots(1,1,figsize=(10,5))
        ax.plot(f, label='FSC orig')
        if randomize:
            ax.plot(fsc_rand, label='FSC rand')
            ax.plot(fsc_corr, label='FSC corrected')

        num_pix = len(fsc_rand)
        ax.set_xticks([0,int(round(0.2*num_pix)),int(round(0.4*num_pix)),int(round(0.6*num_pix)),int(round(0.8*num_pix)),int(round(num_pix))])
        ax.set_xticklabels([num_pix*2*pixelSize, num_pix*2*pixelSize/20, num_pix*2*pixelSize/40, numpy.around(num_pix*2*pixelSize/60,2), num_pix*2*pixelSize/80, num_pix*2*pixelSize/100])

        ax.hlines(fscCriterion,0,num_pix,label=f'cutoff = {fscCriterion}')
        ax.set_title(f'Resolution determined for pixelsize : {pixelSize} at {bandToAngstrom(r[1], pixelSize, numberBands)} Angstrom')
        #fig, ax = subplots(1, 1, figsize=(7, 7))
        #ax.plot(f, label='FSC orig')
        #if randomize:
        #    ax.plot(fsc_rand, label='FSC rand')
        #    ax.plot(fsc_corr, label='FSC corrected')
        ax.legend()
        savefig(os.path.join(outdir, 'FSCOrig_FSCRand_FSCTrue.png'))
        show()
    
