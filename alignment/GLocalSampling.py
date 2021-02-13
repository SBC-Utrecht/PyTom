'''
Routines for Local Sampling and Reference Filtering using Gold standard FSC.
Created July/Aug 2014

@author: FF
'''
from pytom.basic.structures import PyTomClass
from pytom.gpu.initialize import xp, device
from pytom.angles.localSampling import LocalSampling
from pytom.alignment.alignmentStructures import GLocalSamplingJob, ScoringParameters, FLCFScore, SamplingParameters
from pytom.tompy.mpi import MPI
mpi = MPI()


def mainAlignmentLoop(alignmentJob, verbose=False):
    """
    @param alignmentJob: alignment job
    @type alignmentJob: L{pytom.alignment.GLocalSampling.GLocalSamplingJob}
    @param verbose: verbose mode
    @type verbose: L{bool}

    @author: FF
    """

    if 'gpu' in device:
        from pytom.tompy.structures import Preprocessing, Reference, Rotation
        from pytom.tompy.tools import alignVolumesAndFilterByFSC
        from pytom.basic.resolution import bandToAngstrom, getResolutionBandFromFSC, angleFromResolution, \
            write_fsc2Ascii
        from time import time
        from pytom.angles.angleFnc import differenceAngleOfTwoRotations
        from pytom.tompy.io import write, read_size

    else:
        from pytom.alignment.preprocessing import Preprocessing
        from pytom.alignment.localOptimization import alignVolumesAndFilterByFSC
        from pytom.basic.structures import Reference, Rotation
        from pytom.basic.resolution import bandToAngstrom, getResolutionBandFromFSC, angleFromResolution, write_fsc2Ascii
        from time import time
        from pytom.angles.angleFnc import differenceAngleOfTwoRotations

    assert isinstance(object=alignmentJob, class_or_type_or_tuple=GLocalSamplingJob), \
        "mainAlignmentLoop: alignmentJob must be of type GLocalSamplingJob"
    mpi.begin()
    print("particleList    = "+str(alignmentJob.particleList.getFileName()))
    if alignmentJob.scoringParameters.reference.getFilename():
        print("reference       = "+str(alignmentJob.scoringParameters.reference.getFilename()))
    else:
        print("reference       = using average from particleList")

    print("mask            = "+str(alignmentJob.scoringParameters.mask.getFilename()))
    print("rotations       = "+str(alignmentJob.samplingParameters.rotations))
    print("scoring function= "+str(alignmentJob.scoringParameters.score._type))
    print("symmetries      = "+str(alignmentJob.scoringParameters.symmetries))
    print("destination     = "+str(alignmentJob.destination))
    print("numberIterations= "+str(alignmentJob.max_iter))
    print("binning         = "+str(alignmentJob.samplingParameters.binning))
    print("sampleInfo      = "+str(alignmentJob.samplingParameters.sampleInformation))
    print("weighting       = "+str(alignmentJob.scoringParameters.weighting))
    print("compound Wedge  = "+str(alignmentJob.scoringParameters.compoundWedge))
    print("Calculate on GPU= "+str(alignmentJob.gpu))

    for particle in alignmentJob.particleList:
        particle.setScoreValue(-1000.)
    (odd, even) = alignmentJob.particleList.splitOddEven(verbose=verbose)
    progressBar = True
    setParticleNodesRatio = 2
    neven = len(even)
    nodd = len(odd)
    removeAutocorrelation = alignmentJob.scoringParameters.score.getRemoveAutocorrelation()
    # set CompoundWedgeFilenames to None - may be changed later if flag is used
    evenCompoundWedgeFile = None
    oddCompoundWedgeFile = None

    useExternalRef = False
    if alignmentJob.scoringParameters.reference.getFilename() != '':
        useExternalRef = True

    for ii in range(0, alignmentJob.max_iter):
        print(f'running iteration {ii}/{alignmentJob.max_iter}')
        if 'gpu' in device:
            alignmentJob.scoringParameters.mask = alignmentJob.scoringParameters.mask.convert2numpy()

        tt = time()
        #if verbose:
        print(">>>>>>>>> MPI rank: "+str(mpi.rank)+", Iteration: "+str(ii))
        alignmentJob.scoringParameters.score.setRemoveAutocorrelation(flag=removeAutocorrelation)

        # generate averages - if not external reference provided use average to start with
        t1 = time()
        if useExternalRef == False:
            evenAverage = averageParallel(particleList=even,
                                          averageName=alignmentJob.destination+"/"+str(ii)+'-EvenFiltered.em',
                                          showProgressBar=progressBar, verbose=False, createInfoVolumes=False,
                                          weighting=alignmentJob.scoringParameters.weighting, norm=False,
                                          setParticleNodesRatio=setParticleNodesRatio, gpuIDs=alignmentJob.gpu)
            oddAverage = averageParallel(particleList=odd,
                                          averageName=alignmentJob.destination+"/"+str(ii)+'-OddFiltered.em',
                                          showProgressBar=progressBar, verbose=False, createInfoVolumes=False,
                                          weighting=alignmentJob.scoringParameters.weighting, norm=False,
                                          setParticleNodesRatio=setParticleNodesRatio, gpuIDs=alignmentJob.gpu)
            #write un-filtered averages for info
            if 'gpu' in device:
                print(evenAverage.__class__)
                write(alignmentJob.destination + "/" + str(ii) + '-Even.em', evenAverage.getVolume())
                write(alignmentJob.destination + "/" + str(ii) + '-Odd.em', oddAverage.getVolume())
                nband=read_size(alignmentJob.destination + "/" + str(ii) + '-Even.em','x') //2
                temp_mask = alignmentJob.scoringParameters.mask.getVolume()
            else:
                evenAverage.getVolume().write(alignmentJob.destination+"/"+str(ii)+'-Even.em')
                oddAverage.getVolume().write(alignmentJob.destination+"/"+str(ii)+'-Odd.em')
                nband= evenAverage.getVolume().sizeX()/2
                temp_mask = alignmentJob.scoringParameters.mask.getVolume()

            t2 = time()
            print(">>>>>>>>> averaging done ... took %3.2f seconds" % (t2-t1))
            # filter both volumes by sqrt(FSC)
            (averageEven, averageOdd, fsc, fil, optiRot, optiTrans) = \
                alignVolumesAndFilterByFSC(vol1=evenAverage.getVolume(), vol2=oddAverage.getVolume(),
                                           mask=temp_mask,
                                           nband=nband,
                                           interpolation='linear',
                                           fsc_criterion=alignmentJob.scoringParameters.fsc_criterion,
                                           verbose=True)
            # write average from all particle with correctly rotated odd average
            averageAllVolume = evenAverage.getVolume() + oddAverage.getVolume()
            if 'gpu' in device:
                write(alignmentJob.destination+"/"+str(ii)+'-All.em', averageAllVolume)
            else:
                averageAllVolume.write(alignmentJob.destination+"/"+str(ii)+'-All.em')

            t1 = time()
            print(">>>>>>>>> even and odd averages aligned ... took %3.2f seconds" % (t1-t2))
            try:
                write_fsc2Ascii(fsc=fsc, filename=alignmentJob.destination+"/"+str(ii)+'-FSC.dat')
                # default name of Filter files
                write_fsc2Ascii(fsc=fil, filename=alignmentJob.destination+"/"+str(ii)+'-Filter.dat')
            except:
                pass
            resolutionBand = getResolutionBandFromFSC(fsc, criterion=alignmentJob.scoringParameters.fsc_criterion)
            resolutionAngstrom = bandToAngstrom( band=resolutionBand,
                            pixelSize=alignmentJob.samplingParameters.sampleInformation.getPixelSize(),
                            numberOfBands=len(fsc), upscale=1)
            # read un-corrected averages back in for compoundWedge
            if alignmentJob.scoringParameters.compoundWedge:
                from pytom_volume import read
                averageEven = read(alignmentJob.destination+"/"+str(ii)+'-EvenFiltered-PreWedge.em')
                averageOdd  = read(alignmentJob.destination+"/"+str(ii)+'-OddFiltered-PreWedge.em')
                evenCompoundWedgeFile = alignmentJob.destination+"/"+str(ii)+"-EvenFiltered-WedgeSumUnscaled.em"
                oddCompoundWedgeFile = alignmentJob.destination+"/"+str(ii)+"-OddFiltered-WedgeSumUnscaled.em"

            # adjust rotational sampling
            if (type(alignmentJob.samplingParameters.rotations) == LocalSampling) and \
                    (alignmentJob.samplingParameters.adaptive_res != .0):
                if alignmentJob.samplingParameters.sampleInformation.getParticleDiameter()<0.:
                    raise ValueError("mainAlignmentLoop: set particle diameter or switch off adaptive resolution")
                angularIncrement = angleFromResolution(resolution=resolutionAngstrom,
                            particleDiameter=alignmentJob.samplingParameters.sampleInformation.getParticleDiameter())
                angularIncrement = xp.around(alignmentJob.samplingParameters.adaptive_res * angularIncrement, 1)
                print(">>>>>>>>> Iteration "+str(ii)+": Resolution = %3.2f A; angularIncrement= %2.1f deg." % \
                      (resolutionAngstrom, angularIncrement))
                alignmentJob.samplingParameters.rotations.setIncrement(increment=angularIncrement)
            else:
                print(">>>>>>>>> Iteration "+str(ii)+": Resolution = %3.2f A." % resolutionAngstrom)
        else:
            averageEven = alignmentJob.scoringParameters.reference.getVolume()
            averageOdd  = alignmentJob.scoringParameters.reference.getVolume()
            # override autocorr flag for external ref
            alignmentJob.scoringParameters.score.setRemoveAutocorrelation(flag=False)
            resolutionAngstrom = None

        # external reference only possible for 1st iteration
        useExternalRef = False
        averageAllVolume = averageEven + averageOdd

        resolution = f'_{float(resolutionAngstrom):.3f}' if not (resolutionAngstrom is None) else ''

        if 'gpu' in device:
            write(f'{alignmentJob.destination}/{ii}-AllFiltered{resolution}.em', averageAllVolume)
            write(alignmentJob.destination + "/" + str(ii) + "-EvenFiltered.em", averageEven)
            write(alignmentJob.destination + "/" + str(ii) + "-OddFiltered.em", averageOdd)
        else:
            averageAllVolume.write(f'{alignmentJob.destination}/{ii}-AllFiltered{resolution}.em')
            averageEven.write(alignmentJob.destination + "/" + str(ii) + "-EvenFiltered.em")
            averageOdd.write(alignmentJob.destination + "/" + str(ii) + "-OddFiltered.em")


        currentReferenceEven = Reference( referenceFile=alignmentJob.destination+"/"+str(ii)+"-EvenFiltered.em",
                                          generatedByParticleList=even)
        currentReferenceOdd  = Reference( referenceFile=alignmentJob.destination+"/"+str(ii)+"-OddFiltered.em",
                                          generatedByParticleList=odd)
        # align particleLists
        if verbose:
            print("mainAlignmentLoop: CurrentScore XML: "+str(alignmentJob.scoringParameters.score))
            print("mainAlignmentLoop: CurrentRotations XML: "+str(alignmentJob.samplingParameters.rotations))
            print("mainAlignmentLoop: CurrentMask XML: "+str(alignmentJob.scoringParameters.mask))
        alignmentJob.scoringParameters.score.toXMLFile(filename=alignmentJob.destination+"/"+'CurrentScore.xml')
        alignmentJob.samplingParameters.rotations.toXMLFile(filename=alignmentJob.destination+"/"+'CurrentRotations.xml')
        alignmentJob.scoringParameters.mask.toXMLFile(filename=alignmentJob.destination+"/"+'CurrentMask.xml')
        # split particle lists
        evenSplitList = splitParticleList(particleList=even, setParticleNodesRatio=setParticleNodesRatio)
        oddSplitList  = splitParticleList(particleList=odd, setParticleNodesRatio=setParticleNodesRatio)
        print(">>>>>>>>> Aligning Even ....")

        
        if alignmentJob.gpu is None or alignmentJob.gpu == []:
            print('CPU variant is running')
            bestPeaksEvenSplit = mpi.parfor( alignParticleList,
                                    list(zip(evenSplitList, [currentReferenceEven]*len(evenSplitList),
                                        [evenCompoundWedgeFile]*len(evenSplitList),
                                        [alignmentJob.destination+"/"+'CurrentRotations.xml']*len(evenSplitList),
                                        [alignmentJob.destination+"/"+'CurrentScore.xml']*len(evenSplitList),
                                        [alignmentJob.destination+"/"+'CurrentMask.xml']*len(evenSplitList),
                                        [alignmentJob.scoringParameters.preprocessing]*len(evenSplitList),
                                        [progressBar]*neven, [alignmentJob.samplingParameters.binning]*len(evenSplitList),
                                        [verbose]*len(evenSplitList))))
            print(">>>>>>>>> Aligning Odd  ....")
            bestPeaksOddSplit = mpi.parfor( alignParticleList,
                                    list(zip(oddSplitList, [currentReferenceOdd]*len(oddSplitList),
                                        [oddCompoundWedgeFile]*len(oddSplitList),
                                        [alignmentJob.destination+"/"+'CurrentRotations.xml']*len(oddSplitList),
                                        [alignmentJob.destination+"/"+'CurrentScore.xml']*len(oddSplitList),
                                        [alignmentJob.destination+"/"+'CurrentMask.xml']*len(oddSplitList),
                                        [alignmentJob.scoringParameters.preprocessing]*len(oddSplitList),
                                        [progressBar]*nodd, [alignmentJob.samplingParameters.binning]*len(oddSplitList),
                                        [verbose]*len(oddSplitList))))
            # merge peak lists
            bestPeaksEven = mergeLists(bestPeaksEvenSplit)
            bestPeaksOdd = mergeLists(bestPeaksOddSplit)
            # set orientations, translations, and cc values of particles # better update in alignOneParticle???
            even.updateFromPeaks(peaks=bestPeaksEven)
            odd.updateFromPeaks(peaks=bestPeaksOdd)
            print(">>>>>>>>> Average scores: even %2.3f; odd %2.3f" % (even.averageScore(), odd.averageScore()))
            alignmentJob.particleList.updateFromOddEven(odd, even)

            # reset initial bandpass to zero and pre-process ONLY reference with bandpass in further iterations,
            # not particles
            alignmentJob.scoringParameters.preprocessing = Preprocessing()

            alignmentJob.toXMLFile(filename=alignmentJob.destination+"/"+str(ii)+'-GLocalAlignmentJob.xml')
            alignmentJob.particleList.toXMLFile(filename=alignmentJob.destination+"/"+str(ii)+"-ParticleList.xml")
            even.toXMLFile(filename=alignmentJob.destination+"/"+str(ii)+'-ParticleListEven.xml')
            odd.toXMLFile(filename=alignmentJob.destination+"/"+str(ii)+'-ParticleListOdd.xml')

            ## finally average ALL particles in ONE average and determine FSC
            evenAverage = averageParallel(particleList=even,
                                          averageName=alignmentJob.destination+"/average-Final-Even.em",
                                          showProgressBar=progressBar, verbose=False, createInfoVolumes=False,
                                          weighting=alignmentJob.scoringParameters.weighting, norm=False,
                                          setParticleNodesRatio=setParticleNodesRatio)
            oddAverage = averageParallel(particleList=odd,
                                         averageName=alignmentJob.destination+"/average-Final-Odd.em",
                                         showProgressBar=progressBar, verbose=False, createInfoVolumes=False,
                                         weighting=alignmentJob.scoringParameters.weighting, norm=False,
                                         setParticleNodesRatio=setParticleNodesRatio)
            # filter both volumes by sqrt(FSC)
            (averageEven, averageOdd, fsc, fil, optiRot, optiTrans) = \
            alignVolumesAndFilterByFSC(vol1=evenAverage.getVolume(), vol2=oddAverage.getVolume(),
                                       mask=alignmentJob.scoringParameters.mask.getVolume(),
                                       nband=evenAverage.getVolume().sizeX()/2,
                                       interpolation='linear',
                                       fsc_criterion=alignmentJob.scoringParameters.fsc_criterion,
                                       verbose=False)
            # apply rotations also to odd particle list
            if (differenceAngleOfTwoRotations(rotation1=optiRot, rotation2=Rotation(0, 0, 0)) >
                    alignmentJob.samplingParameters.rotations.getIncrement()):
                odd.addRotation(rot=optiRot.invert())
                odd.addShift(translation=optiTrans.invert())
                alignmentJob.particleList.updateFromOddEven(odd, even)
                print("rotation between averages > increment .... applying rotation and shift to odd particles")
                oddAverage = averageParallel(particleList=odd,
                                             averageName=alignmentJob.destination+"/average-Final-Odd.em",
                                             showProgressBar=progressBar, verbose=False, createInfoVolumes=False,
                                             weighting=alignmentJob.scoringParameters.weighting, norm=False,
                                             setParticleNodesRatio=setParticleNodesRatio)
            final_average = averageParallel(particleList=alignmentJob.particleList,
                                            averageName=alignmentJob.destination+"/average-Final.em",
                                            showProgressBar=progressBar, verbose=False, createInfoVolumes=False,
                                            weighting=alignmentJob.scoringParameters.weighting, norm=False,
                                            setParticleNodesRatio=setParticleNodesRatio)
            from pytom.basic.correlation import FSC
            fsc = FSC(volume1=evenAverage.getVolume(), volume2=oddAverage.getVolume(),
                      numberBands=int(evenAverage.getVolume().sizeX()/2))
            #resolution hokus pokus -> estimate fsc for all particles
            for (ii, fscel) in enumerate(fsc):
                fsc[ii] = 2.*fscel/(1.+fscel)
            try:write_fsc2Ascii(fsc=fsc, filename=alignmentJob.destination+"/FSC-Final.dat")
            except: pass
            resolutionBand = getResolutionBandFromFSC(fsc, criterion=0.143)
            resolutionAngstrom = bandToAngstrom(band=resolutionBand,
                                                pixelSize=alignmentJob.samplingParameters.sampleInformation.getPixelSize(),
                                                numberOfBands=len(fsc), upscale=1)
            print(">>>>>>>>>> Final Resolution = %3.2f A." % resolutionAngstrom)

            # filter final average according to resolution
            from pytom.basic.filter import lowpassFilter
            filtered_final = lowpassFilter(volume=final_average.getVolume(), band=resolutionBand, smooth=resolutionBand/10,
                                           fourierOnly=False)[0]
            filtered_final.write(alignmentJob.destination+"/average-FinalFiltered_"+str(resolutionAngstrom)+".em")
            # clean up temporary files
            from os import remove
            remove(alignmentJob.destination+"/"+'CurrentRotations.xml')
            remove(alignmentJob.destination+"/"+'CurrentScore.xml')
            remove(alignmentJob.destination+"/"+'CurrentMask.xml')

        else:
            from pytom.gpu.gpuFunctions import applyFourierFilter
            from pytom.tompy.io import read

            resultsEven = mpi.parfor(alignParticleListGPU, list(zip(evenSplitList,
                                        [currentReferenceEven] * len(evenSplitList),
                                        [evenCompoundWedgeFile] * len(evenSplitList),
                                        [alignmentJob.destination + "/" + 'CurrentRotations.xml'] * len(evenSplitList),
                                        [alignmentJob.destination + "/" + 'CurrentScore.xml'] * len(evenSplitList),
                                        [alignmentJob.destination + "/" + 'CurrentMask.xml'] * len(evenSplitList),
                                        [alignmentJob.scoringParameters.preprocessing] * len(evenSplitList),
                                        [progressBar] * neven,
                                        [alignmentJob.samplingParameters.binning] * len(evenSplitList),
                                        [verbose] * len(evenSplitList),alignmentJob.gpu)))
            print(">>>>>>>>> Aligning Odd  ....")

            resultsOdd = mpi.parfor(alignParticleListGPU, list(zip(oddSplitList,
                                        [currentReferenceOdd]*len(oddSplitList),
                                        [oddCompoundWedgeFile]*len(oddSplitList),
                                        [alignmentJob.destination+"/"+'CurrentRotations.xml']*len(oddSplitList),
                                        [alignmentJob.destination+"/"+'CurrentScore.xml']*len(oddSplitList),
                                        [alignmentJob.destination+"/"+'CurrentMask.xml']*len(oddSplitList),
                                        [alignmentJob.scoringParameters.preprocessing]*len(oddSplitList),
                                        [progressBar]*nodd, [alignmentJob.samplingParameters.binning]*len(oddSplitList),
                                        [verbose]*len(oddSplitList), alignmentJob.gpu)))

            xp.cuda.Device(alignmentJob.gpu[0]).use()
            bestPeaksEvenSplit, plansEven = zip(*resultsEven)
            bestPeaksOddSplit, plansOdd = zip(*resultsOdd)
            bestPeaksEven = mergeLists(bestPeaksEvenSplit)
            bestPeaksOdd = mergeLists(bestPeaksOddSplit)
            #print(len(bestPeaksEven), len(bestPeaksOdd), len(even), len(odd), resultsEven)
            # set orientations, translations, and cc values of particles # better update in alignOneParticle???
            even.updateFromPeaks(peaks=bestPeaksEven)
            odd.updateFromPeaks(peaks=bestPeaksOdd)
            print(">>>>>>>>> Average scores: even %2.3f; odd %2.3f" % (even.averageScore(), odd.averageScore()))
            alignmentJob.particleList.updateFromOddEven(odd, even)

            # reset initial bandpass to zero and pre-process ONLY reference with bandpass in further iterations,
            # not particles
            alignmentJob.scoringParameters.preprocessing = Preprocessing()

            alignmentJob.toXMLFile(filename=alignmentJob.destination + "/" + str(ii) + '-GLocalAlignmentJob.xml')
            alignmentJob.particleList.toXMLFile(filename=alignmentJob.destination + "/" + str(ii) + "-ParticleList.xml")
            even.toXMLFile(filename=alignmentJob.destination + "/" + str(ii) + '-ParticleListEven.xml')
            odd.toXMLFile(filename=alignmentJob.destination + "/" + str(ii) + '-ParticleListOdd.xml')

            ## finally average ALL particles in ONE average and determine FSC



            average = xp.zeros_like(plansEven[0].sumParticles)
            weight = xp.zeros_like(plansEven[0].sumWeights)


            cvols = {}

            for name in ('Even', 'Odd'):
                if name == 'Even':
                    evenAverage = averageParallel(particleList=even,
                                              averageName=alignmentJob.destination + "/average-Final-Even.em",
                                              showProgressBar=progressBar, verbose=False, createInfoVolumes=False,
                                              weighting=alignmentJob.scoringParameters.weighting, norm=False,
                                              setParticleNodesRatio=setParticleNodesRatio, gpuIDs=alignmentJob.gpu)
                    del plansEven
                else:
                    oddAverage = averageParallel(particleList=odd,
                                             averageName=alignmentJob.destination + "/average-Final-Odd.em",
                                             showProgressBar=progressBar, verbose=False, createInfoVolumes=False,
                                             weighting=alignmentJob.scoringParameters.weighting, norm=False,
                                             setParticleNodesRatio=setParticleNodesRatio, gpuIDs=alignmentJob.gpu)
                    del plansOdd
                # if name == 'Even':
                #     for plan in plansEven:
                #         average += xp.array(plan.sumParticles.get(),dtype=xp.float32)
                #         weight += xp.array(plan.sumWeights.get(),dtype=xp.float32)
                #
                #     del plansEven
                #
                # if name == 'Odd':
                #     for plan in plansOdd:
                #         average += xp.array(plan.sumParticles.get(),dtype=xp.float32)
                #         weight += xp.array(plan.sumWeights.get(),dtype=xp.float32)
                #     del plansOdd
                #
                fname = f"{alignmentJob.destination}/average-Final-{name}.em"
                # # fname2 = f"{alignmentJob.destination}/average-FinalWeight-{name}.em"
                # # fname3 = f"{alignmentJob.destination}/average-FinalSum-{name}.em"
                # # write(fname3, average)
                # # write(fname2, weight)
                #
                #write(fname, applyFourierFilter(average, weight))
                # average *= 0
                # weight *= 0
                cvols[name] = read(fname)

            # filter both volumes by sqrt(FSC)
            (averageEven, averageOdd, fsc, fil, optiRot, optiTrans) = \
                alignVolumesAndFilterByFSC(vol1=cvols['Even'], vol2=cvols['Odd'],
                                           mask=alignmentJob.scoringParameters.mask.getVolume(),
                                           nband=cvols['Even'].shape[0] // 2,
                                           interpolation='linear',
                                           fsc_criterion=alignmentJob.scoringParameters.fsc_criterion,
                                           verbose=False)
            # apply rotations also to odd particle list
            if (differenceAngleOfTwoRotations(rotation1=optiRot, rotation2=Rotation(0, 0, 0)) >
                    alignmentJob.samplingParameters.rotations.getIncrement()):

                odd.addRotation(rot=optiRot.invert())
                odd.addShift(translation=optiTrans.invert().convert2pytomc())
                alignmentJob.particleList.updateFromOddEven(odd, even)
                print("rotation between averages > increment .... applying rotation and shift to odd particles")
                oddAverage = averageParallel(particleList=odd,
                                             averageName=alignmentJob.destination + "/average-Final-Odd.em",
                                             showProgressBar=progressBar, verbose=False, createInfoVolumes=False,
                                             weighting=alignmentJob.scoringParameters.weighting, norm=False,
                                             setParticleNodesRatio=setParticleNodesRatio, gpuIDs=alignmentJob.gpu)
                xp.cuda.Device(alignmentJob.gpu[0]).use()
                oddAverage = oddAverage.getVolume()
            else:
                oddAverage = cvols['Odd']

            final_average = averageParallel(particleList=alignmentJob.particleList,
                                            averageName=alignmentJob.destination + "/average-Final.em",
                                            showProgressBar=progressBar, verbose=False, createInfoVolumes=False,
                                            weighting=alignmentJob.scoringParameters.weighting, norm=False,
                                            setParticleNodesRatio=setParticleNodesRatio,gpuIDs=alignmentJob.gpu)
            xp.cuda.Device(alignmentJob.gpu[0]).use()

            from pytom.tompy.correlation import FSC

            fsc = FSC(volume1=cvols['Even'], volume2=oddAverage,
                      numberBands=int(cvols['Even'].shape[0]// 2))

            # resolution hokus pokus -> estimate fsc for all particles
            for (ii, fscel) in enumerate(fsc):
                fsc[ii] = 2. * fscel / (1. + fscel)

            try:
                write_fsc2Ascii(fsc=fsc, filename=alignmentJob.destination + "/FSC-Final.dat")
            except:
                pass

            resolutionBand = getResolutionBandFromFSC(fsc, criterion=0.143)
            resolutionAngstrom = bandToAngstrom(band=resolutionBand,
                                                pixelSize=alignmentJob.samplingParameters.sampleInformation.getPixelSize(),
                                                numberOfBands=len(fsc), upscale=1)

            print(">>>>>>>>>> Final Resolution = %3.2f A." % resolutionAngstrom)

            # filter final average according to resolution
            from pytom.tompy.filter import bandpass as lowpassFilter
            filtered_final = lowpassFilter(final_average.getVolume(), high=resolutionBand, sigma=resolutionBand / 10)
            write(alignmentJob.destination + f"/average-FinalFiltered_{resolutionAngstrom:.3f}.em", filtered_final)
            # clean up temporary files
            from os import remove
            remove(alignmentJob.destination + "/" + 'CurrentRotations.xml')
            remove(alignmentJob.destination + "/" + 'CurrentScore.xml')
            remove(alignmentJob.destination + "/" + 'CurrentMask.xml')
        print(time()-tt)

    print('destroy mpi!')
    mpi.end()
    print('destroyed mpi')

def alignParticleListWrapper(plFilename, reference, referenceWeighting, rotationsFilename,
                      scoreXMLFilename, maskFilename, preprocessing,
                      progressBar=True, binning=1, verbose=False):
    """
    wrapper to prevent to much parsing - caused MPI freeze :(
    """

def alignParticleList(pl, reference, referenceWeightingFile, rotationsFilename,
                      scoreXMLFilename, maskFilename, preprocessing,
                      progressBar=True, binning=1, verbose=False):
    """
    align a ParticleList

    @param pl: particle list
    @type pl: L{pytom.basic.structures.ParticleList}
    #@param particleFilename:
    @param reference: reference volume
    @type reference: L{pytom_volume.vol}
    @param referenceWeightingFile: File for Fourier weighting of the reference (sum of wedges for instance) = CompoundWedge
    @type referenceWeightingFile: str
    @param rotationsFilename: name of rotations xml file
    @type rotationsFilename: L{str}
    @param maskFilename: name of real-space mask xml file for correlation function
    @type maskFilename: L{str}
    @param scoreXMLFilename: name of XML File of score object
    @type scoreXMLFilename: L{str}  #  L{lxml.etree._Element}
    @param preprocessing: Class storing preprocessing of particle and reference such as bandpass
    @type preprocessing: L{pytom.alignment.preprocessing.Preprocessing}
    @param progressBar: Display progress bar of alignment. False by default.
    @param binning: Is binning applied (currently not properly functioning)
    @param verbose: Print out infos. Writes CC volume to disk!!! Default is False

    @return: Returns the peak list for particle list.
    @author: FF
    """
    from pytom.score.score import fromXMLFile
    from pytom.angles.angle import AngleObject
    from pytom.basic.structures import Mask, ParticleList

    assert type(pl) == ParticleList, "pl is supposed to be a particleList"

    scoreObject = fromXMLFile(filename=scoreXMLFilename)

    rot = AngleObject()
    rotations = rot.fromXMLFile(filename=rotationsFilename)

    mask = Mask()
    mask.fromXMLFile(filename=maskFilename)

    if referenceWeightingFile:
        from pytom_volume import read
        referenceWeighting = read(referenceWeightingFile)
    else:
        referenceWeighting = None

    if verbose:
        print("alignParticleList: rank "+str(mpi.rank))
        print("alignParticleList: angleObject: "+str(rotations))
        print("alignParticleList: scoreObject: "+str(scoreObject))
        print("alignParticleList: mask:        "+str(mask))

    bestPeaks = []
    for particle in pl:
        bestPeak = alignOneParticle(particle=particle, reference=reference, referenceWeighting=referenceWeighting,
                                    rotations=rotations, scoreObject=scoreObject, mask=mask,
                                    preprocessing=preprocessing, progressBar=progressBar, binning=binning,
                                    verbose=verbose)
        bestPeaks.append(bestPeak)

    return bestPeaks

def alignParticleListGPU(pl, reference, referenceWeightingFile, rotationsFilename,
                      scoreXMLFilename, maskFilename, preprocessing,
                      progressBar=True, binning=1, verbose=False, gpuID=0):
    """
    align a ParticleList

    @param pl: particle list
    @type pl: L{pytom.basic.structures.ParticleList}
    #@param particleFilename:
    @param reference: reference volume
    @type reference: L{pytom_volume.vol}
    @param referenceWeightingFile: File for Fourier weighting of the reference (sum of wedges for instance) = CompoundWedge
    @type referenceWeightingFile: str
    @param rotationsFilename: name of rotations xml file
    @type rotationsFilename: L{str}
    @param maskFilename: name of real-space mask xml file for correlation function
    @type maskFilename: L{str}
    @param scoreXMLFilename: name of XML File of score object
    @type scoreXMLFilename: L{str}  #  L{lxml.etree._Element}
    @param preprocessing: Class storing preprocessing of particle and reference such as bandpass
    @type preprocessing: L{pytom.alignment.preprocessing.Preprocessing}
    @param progressBar: Display progress bar of alignment. False by default.
    @param binning: Is binning applied (currently not properly functioning)
    @param verbose: Print out infos. Writes CC volume to disk!!! Default is False

    @return: Returns the peak list for particle list.
    @author: FF
    """
    from pytom.score.score import fromXMLFile
    from pytom.angles.angle import AngleObject
    from pytom.tompy.structures import Wedge, Rotation, Mask, ParticleList
    from pytom.alignment.alignmentFunctions import bestAlignmentGPU
    from pytom.gpu.gpuStructures import GLocalAlignmentPlan
    from time import time
    from pytom.angles.angleFnc import differenceAngleOfTwoRotations
    from pytom.tompy.io import read, read_size
    import numpy as np

    assert type(pl) == ParticleList, "pl must be particleList"

    rot = AngleObject()
    rotations = rot.fromXMLFile(filename=rotationsFilename)

    mask = Mask()
    mask.fromXMLFile(filename=maskFilename)

    wedge = pl[0].getWedge().convert2numpy()

    use_device=f'gpu:{gpuID}'
    print('devices: ', use_device)
    plan = GLocalAlignmentPlan(pl[0], reference, mask, wedge, maskIsSphere=True, cp=xp, device=use_device,
                               interpolation='filt_bspline')

    try:
        preprocessing = preprocessing.convert2numpy()
    except:
        pass

    if verbose:
        print("alignParticleList: rank "+str(mpi.rank))
        print("alignParticleList: angleObject: "+str(rotations))
        print("alignParticleList: scoreObject: "+str(scoreObject))
        print("alignParticleList: mask:        "+str(mask))

    bestPeaks = []
    for particle in pl:
        fname = particle.getFilename()
        t1 = time()
        oldRot = particle.getRotation()
        rotations.setStartRotation(oldRot)
        particleVol = read(particle.getFilename(), deviceID=use_device)
        bestPeak = bestAlignmentGPU(particleVol, rotations, plan, preprocessing=preprocessing, wedgeInfo=wedge)
        bestPeaks += [bestPeak]
        rotations.reset()
        particle.setRotation(bestPeak.getRotation())
        particle.setShift(bestPeak.getShift())
        angDiff = differenceAngleOfTwoRotations(rotation1=bestPeak.getRotation(), rotation2=oldRot)
        t2 = time()
        # print(bestPeak.getRotation())
        # print(bestPeak.getShift())
        print(fname + ": Angular Difference before and after alignment: %2.2f ... took %3.1f seconds ..." % (angDiff, t2 - t1))

    plan.clean()

    return [bestPeaks, plan]

def alignOneParticleWrapper(particle, reference, referenceWeighting=None, rotationsFilename='',
                            scoreXMLFilename='', maskFilename='', preprocessing=None,
                            progressBar=True, binning=1, verbose=False):
    """
    wrapper for alignOneParticle:

    @param particle: particle
    @type particle: L{pytom.basic.structures.Particle}
    #@param particleFilename:
    @param reference: reference volume
    @type reference: L{pytom_volume.vol}
    @param referenceWeighting: Fourier weighting of the reference (sum of wedges for instance)
    @type referenceWeighting: L{pytom.basic.structures.vol}
    @param scoreXMLFilename: name of XML File of score object
    @type scoreXMLFilename: L{str}  #  L{lxml.etree._Element}
    @param rotationsFilename: name of rotations xml file
    @type rotationsFilename: L{str}
    @param maskFilename: name of real-space mask xml file for correlation function
    @type maskFilename: L{str}
    @param preprocessing: Class storing preprocessing of particle and reference such as bandpass
    @type preprocessing: L{pytom.alignment.preprocessing.Preprocessing}
    @param progressBar: Display progress bar of alignment. False by default.
    @param binning: Is binning applied (currently not properly functioning)
    @param verbose: Print out infos. Writes CC volume to disk!!! Default is False

    @return: Returns the best rotation and translation for particle and the corresponding scoring result.
    @rtype: L{pytom.alignment.structures.Peak}
    @author: FF
    """
    from pytom.score.score import fromXMLFile
    from pytom.angles.angle import AngleObject
    from pytom.basic.structures import Mask, Particle

    assert type(particle) == Particle, "particle must be of type Particle"

    scoreObject = fromXMLFile(filename=scoreXMLFilename)

    rot = AngleObject()
    rotations = rot.fromXMLFile(filename=rotationsFilename)

    mask = Mask()
    mask.fromXMLFile(filename=maskFilename)

    if verbose:
        print("alignOneParticleWrapper: rank "+str(mpi.rank))
        print("alignOneParticleWrapper: angleObject: "+str(rotations))
        print("alignOneParticleWrapper: scoreObject: "+str(scoreObject))
        print("alignOneParticleWrapper: mask:        "+str(mask))

    bestPeak = alignOneParticle( particle, reference, referenceWeighting, rotations,
                                 scoreObject, mask, preprocessing, progressBar=progressBar, binning=binning,
                                 verbose=verbose)

    return bestPeak


def alignOneParticle( particle, reference, referenceWeighting, rotations,
                      scoreObject, mask, preprocessing, progressBar=True, binning=1, verbose=False):
    """
    @param particle: particle
    @type particle: L{pytom.basic.Particle}
    @param reference: reference
    @type reference: L{pytom.basic.Reference}
    @param referenceWeighting: Fourier weighting of the reference (sum of wedges, i.e., compound weighting for \
        instance) - if 'False' not performed
    @type referenceWeighting: L{pytom.basic.structures.vol} or str
    @param rotations: All rotations to be scanned
    @type rotations: L{pytom.angles.AngleObject}
    @param scoreObject:
    @type scoreObject: L{pytom.score.score.Score}
    @param mask: real-space mask for correlation function
    @type mask: L{pytom.basic.structures.Particle}
    @param preprocessing: Class storing preprocessing of particle and reference such as bandpass
    @type preprocessing: L{pytom.alignment.preprocessing.Preprocessing}
    @param progressBar: Display progress bar of alignment. False by default.
    @param binning: Is binning applied (currently not properly functioning)
    @param verbose: Print out infos. Writes CC volume to disk!!! Default is False
    @return: Returns the best rotation for particle and the corresponding scoring result.
    @author: FF
    """
    from pytom.basic.structures import Particle, Reference, Mask, Wedge, WedgeInfo
    from pytom.angles.angleList import AngleList
    from pytom.alignment.alignmentFunctions import bestAlignment
    from pytom.score.score import Score
    from pytom.alignment.preprocessing import Preprocessing
    from pytom.angles.angleFnc import differenceAngleOfTwoRotations
    from time import time

    t1 = time()
    fname = particle.getFilename()
    assert type(particle) == Particle, "alignOneParticle: particle not of type Particle"
    partVol = particle.getVolume()

    # prepare reference
    assert type(reference) == Reference, "alignOneParticle: reference not of type Reference"
    if scoreObject.getRemoveAutocorrelation() and reference._generatedByParticleList:
        refVol = reference.subtractParticle( particle=particle, binning=1, verbose=verbose)[0]
    else:
        refVol = reference.getVolume()

    # apply pre-processing to reference, no pre-processing to particles
    if preprocessing is None:
        preprocessing = Preprocessing()
    assert type(preprocessing) == Preprocessing, "alignOneParticle: preprocessing not of type Proprocessing"
    preprocessing.setTaper( taper=refVol.sizeX()/10.)
    refVol = preprocessing.apply(volume=refVol, bypassFlag=True)

    wedge = particle.getWedge()
    assert type(scoreObject) == Score, "alignOneParticle: score not of type Score"
    if mask:
        assert type(mask) == Mask, "alignOneParticle: mask not of type Mask"

    assert type(rotations) == AngleList, "alignOneParticle: rotations not of type AngleList"

    oldRot = particle.getRotation()
    rotations.setStartRotation(oldRot)
    if not referenceWeighting:
        bestPeak = bestAlignment(particle=partVol, reference=refVol,
                                 referenceWeighting='False', wedgeInfo=wedge, rotations=rotations,
                                 scoreObject=scoreObject, mask=mask, preprocessing=None,
                                 progressBar=progressBar, binning=binning, bestPeak=None, verbose=verbose)
    else:
        bestPeak = bestAlignment(particle=partVol, reference=refVol,
                                 referenceWeighting=referenceWeighting, wedgeInfo=wedge, rotations=rotations,
                                 scoreObject=scoreObject, mask=mask, preprocessing=None,
                                 progressBar=progressBar, binning=binning, bestPeak=None, verbose=verbose)
    angDiff = differenceAngleOfTwoRotations(rotation1=bestPeak.getRotation(), rotation2=oldRot)
    t2 = time()
    print(bestPeak.getRotation())
    #print(fname+": Angular Difference before and after alignment: %2.2f ... took %3.1f seconds ..." % (angDiff, t2-t1))
    rotations.reset()

    particle.setRotation(bestPeak.getRotation())
    particle.setShift(bestPeak.getShift())
    return bestPeak




def averageParallel(particleList,averageName, showProgressBar=False, verbose=False,
                    createInfoVolumes=False, weighting=None, norm=False,
                    setParticleNodesRatio=3, gpuIDs=None):
    """
    compute average using parfor
    @param particleList: The particles
    @param averageName: Filename of new average
    @param verbose: Prints particle information. Disabled by default.
    @param createInfoVolumes: Create info data (wedge sum, inverted density) too? False by default.
    @param weighting: weight particles by exp CC in average
    @type weighting: bool
    @param setParticleNodesRatio: minimum number of particles per node
    @type setParticleNodesRatio: L{int}
    @return: A new Reference object
    @rtype: L{pytom.basic.structures.Reference}
    @author: FF

    """
    from pytom_volume import read, complexRealMult
    from pytom.basic.fourier import fft,ifft
    from pytom.basic.filter import lowpassFilter
    from pytom.basic.structures import Reference
    from pytom.alignment.alignmentFunctions import average, invert_WedgeSum

    import os
    splitLists = splitParticleList(particleList, setParticleNodesRatio=setParticleNodesRatio)
    splitFactor = len(splitLists)
    assert splitFactor > 0, "splitFactor == 0, issue with parallelization"

    if 'gpu' in device:
        from pytom.bin.average import averageGPU2 as average
        from pytom.tompy.structures import Reference
        from pytom.tompy.io import read, write
    else:
        gpuIDs = [None,]*splitFactor


    avgNameList = []
    preList = []
    wedgeList = []
    for ii in range(splitFactor):
        avgName = averageName + '_dist' + str(ii) + '.em'
        avgNameList.append(avgName)
        preList.append(averageName + '_dist' + str(ii) + '-PreWedge.em')
        wedgeList.append(averageName + '_dist' + str(ii) + '-WedgeSumUnscaled.em')

    #reference = average(particleList=plist, averageName=xxx, showProgressBar=True, verbose=False,
    # createInfoVolumes=False, weighting=weighting, norm=False)

    averageList = mpi.parfor( average, list(zip(splitLists, avgNameList, [showProgressBar]*splitFactor,
                                           [verbose]*splitFactor, [createInfoVolumes]*splitFactor,
                                           [weighting]*splitFactor, [norm]*splitFactor, gpuIDs)))
    if 'gpu' in device:
        xp.cuda.Device(gpuIDs[0]).use()

    #collect results from files
    unweiAv = read(preList[0])
    wedgeSum = read(wedgeList[0])
    os.system('rm ' + wedgeList[0])
    os.system('rm ' + avgNameList[0])
    os.system('rm ' + preList[0])
    for ii in range(1,splitFactor):
        av = read(preList[ii])
        unweiAv += av
        os.system('rm ' + preList[ii])
        w = read(wedgeList[ii])
        wedgeSum += w
        os.system('rm ' + wedgeList[ii])
        os.system('rm ' + avgNameList[ii])

    if 'gpu' in device:
        from pytom.tompy.tools import invert_WedgeSum
        from pytom.tompy.filter import applyFourierFilter, bandpass
        write(averageName[:len(averageName) - 3] + '-PreWedge.em', unweiAv)
        write(averageName[:len(averageName) - 3] + '-WedgeSumUnscaled.em', wedgeSum)
        wedgeSum = invert_WedgeSum(invol=wedgeSum, r_max=unweiAv.shape[0] / 2 - 2., lowlimit=.05 * len(particleList),
                        lowval=.05 * len(particleList))
        unweiAv = applyFourierFilter(unweiAv, wedgeSum)
        unweiAv = bandpass(unweiAv, high=unweiAv.shape[0]//2-2, sigma=(unweiAv.shape[0]//2-1)/10.)
        write(averageName, unweiAv)
    else:
        unweiAv.write(averageName[:len(averageName)-3]+'-PreWedge.em')
        wedgeSum.write(averageName[:len(averageName)-3] + '-WedgeSumUnscaled.em')

        # convolute unweighted average with inverse of wedge sum
        invert_WedgeSum( invol=wedgeSum, r_max=unweiAv.sizeX()/2-2., lowlimit=.05*len(particleList),
                         lowval=.05*len(particleList))
        fResult = fft(unweiAv)
        r = complexRealMult(fResult,wedgeSum)
        unweiAv = ifft(r)
        unweiAv.shiftscale(0.0,1/float(unweiAv.sizeX()*unweiAv.sizeY()*unweiAv.sizeZ()))
        # low pass filter to remove artifacts at fringes
        unweiAv = lowpassFilter(volume=unweiAv, band=unweiAv.sizeX()/2-2, smooth=(unweiAv.sizeX()/2-1)/10.)[0]

        unweiAv.write(averageName)

    return Reference(averageName,particleList)


def splitParticleList(particleList, setParticleNodesRatio=3):
    """
    @param particleList: The particle list
    @param setParticleNodesRatio: minimum number of particles per node
    @type setParticleNodesRatio: L{int}
    @return: list of particle lists, splitFactor (number of processors or smaller for few particles)
    @rtype: list, L{int}
    @author: FF
    """
    numberOfNodes = mpi.size
    particleNodesRatio = float(len(particleList)) / float(numberOfNodes)
    splitFactor = numberOfNodes
    #make sure each node gets at least setParticleNodesRatio particles.
    if particleNodesRatio < setParticleNodesRatio:
        splitFactor = len(particleList) / int(setParticleNodesRatio)
    assert splitFactor > 0, \
        "splitFactor == 0, too few particles for parallelization - decrease number of processors"
    print(f'splitFactor = {splitFactor}')
    splitLists = particleList.splitNSublists(splitFactor-1)  # somehow better to not include master...
    return splitLists

def mergeLists(lists):
    """
    merge lists that have been split for parallelization
    @param lists: input lists
    @type lists: L{list}
    @return: merged list
    @rtype: L{list}
    @author: FF
    """
    outlist = []
    for tlist in lists:
        for listmember in tlist:
            outlist.append(listmember)
    return outlist

def writeParticleListToUniqueFile(pl, dirName=None):
    """
    write particle list to random filename
    @param pl: particleList
    @type pl: L{pytom.basic.ParticleList}
    @return: filename
    @rtype: str
    """
    from uuid import uuid4
    from pytom.basic.structures import ParticleList
    assert isinstance(object=pl, class_or_type_or_tuple=ParticleList), \
        "writeParticleListToUniqueFile: pl must be of type ParticleList"
    fname = str(uuid4())
    if dirName:
        fname = dirName + '/'+fname
    pl.toXMLFile(filename=fname)
    return fname
