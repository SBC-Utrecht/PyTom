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
    from pytom.tompy.correlation import FSC, determineResolution
    from pytom.tompy.io import read
    
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
                                   ScriptOption('--combinedResolution', 'Calculate the resolution for two volumes combined.', arg=False, optional=True),
                                   ScriptOption(['-v','--verbose'], 'Verbose data. (optional) ', arg=False, optional=True),
                                   ScriptOption(['-g', '--gpuID'], 'Which gpu do you want to use?', True, True),
                                   ScriptOption(['-h', '--help'], 'Help.', arg=False, optional=True)])


    if len(sys.argv) <=2:
        print(helper)
        sys.exit()
    try:
        v1Filename, v2Filename, particleList, fscCriterion, numberBands, mask, pixelSize, xml, randomize, plot, outdir, \
        combined_resolution, verbose, gpuIDs, help = parse_script_options(sys.argv[1:], helper)
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
            numberBands = int(v1.shape[0]//2)
        
        f = FSC(v1,v2,numberBands,mask,verbose)

        if combined_resolution:
            for (ii, fscel) in enumerate(f):
                f[ii] = 2. * fscel / (1. + fscel)

        if verbose: print('FSC:\n', f)

        import pytom.tompy.correlation as correlation
        from pytom_numpy import vol2npy
        import numpy as np
        from pytom.tompy.io import write

        if randomize is None or randomize < 1E-3:
            r = determineResolution(f, fscCriterion, verbose)
        else:
            randomizationFrequency    = np.floor(determineResolution(np.array(f), randomize, verbose)[1])
            oddVolumeRandomizedPhase  = correlation.randomizePhaseBeyondFreq(v1, randomizationFrequency)
            evenVolumeRandomizedPhase = correlation.randomizePhaseBeyondFreq(v2, randomizationFrequency)
            # write(os.path.join(outdir, 'randOdd.mrc'), oddVolumeRandomizedPhase)
            # write(os.path.join(outdir, 'randEven.mrc'), evenVolumeRandomizedPhase)
            # oddVolumeRandomizedPhase = read(os.path.join(outdir, 'randOdd.mrc'))
            # evenVolumeRandomizedPhase = read(os.path.join(outdir, 'randEven.mrc'))
            fsc_rand = FSC(oddVolumeRandomizedPhase, evenVolumeRandomizedPhase, numberBands, mask, verbose)

            if combined_resolution:
                for (ii, fscel) in enumerate(fsc_rand):
                    fsc_rand[ii] = 2.*fscel/(1.+fscel)

            if verbose:

                print('FSC_Random:\n', fsc_rand)
            fsc_corr = list(correlation.calc_FSC_true(np.array(f),np.array(fsc_rand)))
            if verbose:
                print('FSC_true:\n', fsc_corr)

            r = determineResolution(fsc_corr, fscCriterion, verbose)
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
            resolution = bandToAngstrom(r[1], pixelSize, numberBands)
            print('Resolution determined for pixelsize : ', pixelSize , ' at ', resolution, ' Angstrom')
    
    else:
        print('XML')

    from pytom.basic.resolution import write_fsc2Ascii

    write_fsc2Ascii(fsc=f, filename=outdir+"/FSCOrig.dat")
    if randomize:
        write_fsc2Ascii(fsc=fsc_rand, filename=outdir + "/FSCRand.dat")
        write_fsc2Ascii(fsc=fsc_corr, filename=outdir + "/FSCCorr.dat")



    if plot and v1Filename and v2Filename:
        from pytom.plotting.plotFSC import plot_FSC
        files = [outdir+"/FSCOrig.dat"]
        if randomize:
            files += [outdir+"/FSCRand.dat",outdir+"/FSCCorr.dat"]
        plot_FSC(files, pixelSize, boxsize=v1.shape[0],
                 outname=os.path.join(outdir, 'FSCOrig_FSCRand_FSCTrue.png'), show_image=True, c=fscCriterion,
                 resolution=resolution, rand=randomize)

        # import matplotlib
        # import numpy
        # matplotlib.use('Qt5Agg')
        # from pylab import subplots, savefig, show
        #
        #
        # fig,ax = subplots(1,1,figsize=(10,5))
        # size = len(f)
        # newax = ax.twiny()
        # x1 = numpy.arange(1, size + 1) * (1 / (float(pixelSize) * int(size*2)))
        # x = numpy.arange(size)
        #
        # newax.plot(x, f, lw=2)
        # ax.plot(x, f, label='FSC orig', lw=2)
        #
        # l = len(f)
        # #size = len(f)*2
        #
        # if randomize:
        #     ax.plot(x, fsc_rand, label='FSC rand')
        #     ax.plot(x, fsc_corr, label='FSC corrected')
        #
        # # xticklabels = [max(2,(r*size))//2 for r in numpy.arange(0,1.1,0.2)]
        # # ax.set_xticks(xticklabels)
        # # ax.set_xticklabels(xticklabels)
        #
        # try: ax.vlines(r[1], 0, fscCriterion, linestyles='dashed', label=f'resolution = {float(resolution):8.3f}', )
        # except Exception as e: print(e)# ax2 = ax.twiny()
        # # ax2.set_xticks(xticklabels)
        # # ax2.set_xticklabels([size*pixelSize/xtick for xtick in xticklabels])
        # ax.set_xticks([0, x[size // 4], x[size // 2], x[3 * size // 4], x[-1]])
        # ax.set_xticklabels([f'{(1/x1[0]):.1f}', f'{1 / x1[size // 4 -1]:.1f}', f'{1 / x1[size // 2-1]:.1f}',
        #                     f'{1 / x1[3 * size // 4-1]:.1f}', f'{1 / x1[-1]:.1f}'])
        #
        # ax.hlines(fscCriterion,0, l,label=f'cutoff = {fscCriterion}')
        # ax.set_title(f'resolution: {float(resolution):8.3f} A')
        # #fig, ax = subplots(1, 1, figsize=(7, 7))
        # #ax.plot(f, label='FSC orig')
        # #if randomize:
        # #    ax.plot(fsc_rand, label='FSC rand')
        # #    ax.plot(fsc_corr, label='FSC corrected')
        # ax.legend()
        # savefig(os.path.join(outdir, 'FSCOrig_FSCRand_FSCTrue.png'))
        # show()
        #
