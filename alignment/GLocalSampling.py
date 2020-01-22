'''
Routines for Local Sampling and Reference Filtering using Gold standard FSC.
Created July/Aug 2014

@author: FF
'''

from pytom.basic.structures import PyTomClass

from pytom.angles.localSampling import LocalSampling
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
                                          setParticleNodesRatio=setParticleNodesRatio)
            oddAverage = averageParallel(particleList=odd,
                                          averageName=alignmentJob.destination+"/"+str(ii)+'-OddFiltered.em',
                                          showProgressBar=progressBar, verbose=False, createInfoVolumes=False,
                                          weighting=alignmentJob.scoringParameters.weighting, norm=False,
                                          setParticleNodesRatio=setParticleNodesRatio)
            #write un-filtered averages for info
            evenAverage.getVolume().write(alignmentJob.destination+"/"+str(ii)+'-Even.em')
            oddAverage.getVolume().write(alignmentJob.destination+"/"+str(ii)+'-Odd.em')
            t2 = time()
            print(">>>>>>>>> averaging done ... took %3.2f seconds" % (t2-t1))
            # filter both volumes by sqrt(FSC)
            (averageEven, averageOdd, fsc, fil, optiRot, optiTrans) = \
                alignVolumesAndFilterByFSC(vol1=evenAverage.getVolume(), vol2=oddAverage.getVolume(),
                                           mask=alignmentJob.scoringParameters.mask.getVolume(),
                                           nband=evenAverage.getVolume().sizeX()/2,
                                           interpolation='linear',
                                           fsc_criterion=alignmentJob.scoringParameters.fsc_criterion,
                                           verbose=True)
            # write average from all particle with correctly rotated odd average
            averageAllVolume = evenAverage.getVolume() + oddAverage.getVolume()
            averageAllVolume.write(alignmentJob.destination+"/"+str(ii)+'-All.em')
            t1 = time()
            print(">>>>>>>>> even and odd averages aligned ... took %3.2f seconds" % (t1-t2))
            write_fsc2Ascii(fsc=fsc, filename=alignmentJob.destination+"/"+str(ii)+'-FSC.dat')
            # default name of Filter files
            write_fsc2Ascii(fsc=fil, filename=alignmentJob.destination+"/"+str(ii)+'-Filter.dat')
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
                angularIncrement = round(alignmentJob.samplingParameters.adaptive_res * angularIncrement, 1)
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
        averageEven.write( alignmentJob.destination+"/"+str(ii)+"-EvenFiltered.em")
        averageOdd.write( alignmentJob.destination+"/"+str(ii)+"-OddFiltered.em")
        averageAllVolume = averageEven + averageOdd
        if resolutionAngstrom:
            averageAllVolume.write(alignmentJob.destination+"/"+str(ii)+"-AllFiltered_"+str(resolutionAngstrom)+".em")
        else:
            averageAllVolume.write(alignmentJob.destination+"/"+str(ii)+"-AllFiltered.em")
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
        oddSplitList = splitParticleList(particleList=odd, setParticleNodesRatio=setParticleNodesRatio)
        print(">>>>>>>>> Aligning Even ....")


        if alignmentJob.gpu is None:

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
            write_fsc2Ascii(fsc=fsc, filename=alignmentJob.destination+"/FSC-Final.dat")
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
            bestPeaksOdd, plansOdd = alignParticleListGPU(list(zip(oddSplitList,
                                        [currentReferenceOdd]*len(oddSplitList),
                                        [oddCompoundWedgeFile]*len(oddSplitList),
                                        [alignmentJob.destination+"/"+'CurrentRotations.xml']*len(oddSplitList),
                                        [alignmentJob.destination+"/"+'CurrentScore.xml']*len(oddSplitList),
                                        [alignmentJob.destination+"/"+'CurrentMask.xml']*len(oddSplitList),
                                        [alignmentJob.scoringParameters.preprocessing]*len(oddSplitList),
                                        [progressBar]*nodd, [alignmentJob.samplingParameters.binning]*len(oddSplitList),
                                        [verbose]*len(oddSplitList))))
            bestPeaksEven, plansEven = alignParticleListGPU(list(zip(evenSplitList,
                                        [currentReferenceEven]*len(evenSplitList),
                                        [evenCompoundWedgeFile]*len(evenSplitList),
                                        [alignmentJob.destination+"/"+'CurrentRotations.xml']*len(evenSplitList),
                                        [alignmentJob.destination+"/"+'CurrentScore.xml']*len(evenSplitList),
                                        [alignmentJob.destination+"/"+'CurrentMask.xml']*len(evenSplitList),
                                        [alignmentJob.scoringParameters.preprocessing]*len(evenSplitList),
                                        [progressBar]*neven, [alignmentJob.samplingParameters.binning]*len(evenSplitList),
                                        [verbose]*len(evenSplitList))))


    mpi.end()

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
    from pytom.gpu.plans import prepare_glocal_plan
    
    assert type(pl) == ParticleList, "pl must be particleList"

    rot = AngleObject()
    rotations = rot.fromXMLFile(filename=rotationsFilename)
    #print(np.array(list(rotations)))


    particle = np.zeros(read_size(pl[0].getFilename()))
    reference = vol2npy(reference).astype(no.float32)
    mask = read(maskFilename).astype(np.float32)
    wedge = np.ones_like(particle,dtype=np.float32)

    plan = prepare_glocal_plan(particle, reference, mask, wedge)
    
    if verbose:
        print("alignParticleList: rank "+str(mpi.rank))
        print("alignParticleList: angleObject: "+str(rotations))
        print("alignParticleList: scoreObject: "+str(scoreObject))
        print("alignParticleList: mask:        "+str(mask))

    bestPeaks = []
    for particle in pl:
        plan.particle = ga.to_gpu(read(particle.getFilename()))
        removeCurrent(plan.particle, plan.reference, plan.currentRef)
        oldRot = particle.getRotation()
        setStartRotation(plan.currentRef, oldRot)
        bestPeakList = alignOneParticleGPU([particle], rotations=rotations, binning=binning)
        bestPeaks += bestPeakList
        
    return bestPeaks, plan

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
    print(fname+": Angular Difference before and after alignment: %2.2f ... took %3.1f seconds ..." % (angDiff, t2-t1))
    rotations.reset()

    particle.setRotation(bestPeak.getRotation())
    particle.setShift(bestPeak.getShift())
    return bestPeak


class SamplingParameters(PyTomClass):
    """
    class to store parameters for sampling
    """
    def __init__(self, rotations=LocalSampling(shells=3, increment=3., z1Start=0., z2Start=0., xStart=0.),
                 binning=1, adaptive_res=0.1, sample_info=None):
        """
        @param rotations: rotations for orientation sampling
        @type rotations: L{pytom.angles.angle.AngleObject}
        @param binning: binning FACTOR (1=no binning, 2=2x2x2 voxel-> 1 voxel, etc.)
        @type binning: L{int}
        @param adaptive_res: adaptive resolution offset
        @type adaptive_res: L{float}
        @param sample_info: info about sample
        @type sample_info: L{pytom.basic.structures.SampleInformation}
        @author: FF
        """
        from pytom.angles.angle import AngleObject
        from pytom.basic.structures import SampleInformation

        if rotations == None:
            rotations = LocalSampling(shells=3, increment=3., z1Start=0., z2Start=0., xStart=0.)
        assert type(rotations) == AngleObject, "SamplingParameters: rotations must be of type AngleObject"
        self.rotations = rotations
        assert type(binning) == int, "SamplingParameters: binning must be of type int"
        self.binning = binning
        assert type(adaptive_res) == float, "SamplingParameters: adaptive_res must be of type float"
        self.adaptive_res = adaptive_res
        if sample_info == None:
            sample_info = SampleInformation()
        assert type(sample_info) == SampleInformation, "SamplingParameters: sample_info must be of type SampleInformation"
        self.sampleInformation = sample_info

    def toXML(self):
        """
        generate XML of Sampling Parameters
        @return: xml object of SamplingParameters
        @rtype: L{etree.Element}
        @author: FF
        """
        from lxml import etree

        xmlObj = etree.Element("SamplingParameters")
        xmlObj.append(self.rotations.toXML())
        xmlObj.append(self.sampleInformation.toXML())
        xmlObj.set("Binning", str(self.binning))
        if self.adaptive_res is False:
            xmlObj.set("AdaptiveResolution", '0')
        else:
            xmlObj.set("AdaptiveResolution", str(self.adaptive_res))
        return xmlObj

    def fromXML(self,xmlObj):
        """
        @param xmlObj: xml object
        @type xmlObj: L{lxml.etree._Element}
        @author: FF
        """
        from lxml.etree import _Element
        from pytom.angles.angle import AngleObject
        from pytom.basic.structures import SampleInformation

        if xmlObj.__class__ != _Element:
            raise Exception('You must provide a valid XML object.')
        if xmlObj.tag == "SamplingParameters":
            samplingDescription = xmlObj
        else:
            samplingDescription = xmlObj.xpath('SamplingParameters')
            if len(samplingDescription) == 0:
                raise Exception("This XML is not a SamplingParameters.")
            samplingDescription = samplingDescription[0]

        self.binning = int(samplingDescription.get('Binning'))
        self.adaptive_res = float(samplingDescription.get('AdaptiveResolution'))

        rot = AngleObject()
        self.rotations = rot.fromXML( samplingDescription.xpath('Angles')[0])

        try:
            si = samplingDescription.xpath('SampleInformation')[0]
            self.sampleInformation = SampleInformation()
            self.sampleInformation.fromXML(si)
        except:
            self.sampleInformation = SampleInformation()


from pytom.score.score import FLCFScore


class ScoringParameters(PyTomClass):
    """
    class to store parameters for subtomo alignment score
    @author: FF
    """
    def __init__(self, score=FLCFScore(), ref=None, weighting=False, compoundWedge=False, mask=None,
                 preprocessing=None, fsc_criterion=0.143, symmetries=None):
        """
        @param mask: mask
        @type mask: L{pytom.basic.structures.Mask}
        @param preprocessing: pre-processing parameters (bandpass filter)
        @type preprocessing: L{pytom.alignment.preprocessing.Preprocessing}
        @param weighting: weighting of particles in average by exp(cc)
        @type weighting: bool
        @param compoundWedge: compoundWedge weighting of reference - i.e. used unweighted average for alignment but \
                                rather weight particle accordingly
        @type compoundWedge: bool
        @param symmetries: apply point symmetry for average
        @type symmetries: ???
        @param fsc_criterion: FSC criterion
        @type fsc_criterion: float

        @author: FF
        """
        from pytom.basic.structures import Mask, Reference, Symmetry
        from pytom.alignment.preprocessing import Preprocessing
        from pytom.score.score import Score

        assert type(score) == Score, "ScoringParameters: input score not of pytom type Score"
        self.score = score
        assert (type(ref) == Reference) or (type(ref) == type(None)), \
            "ScoringParameters: input ref not of pytom type Reference or None"
        self.reference = ref
        assert (type(mask) == Mask) or (type(mask) == type(None)), "ScoringParameters: input mask not of pytom type Mask"
        self.mask=mask
        assert (type(preprocessing) == Preprocessing) or (type(preprocessing) == type(None)), \
            "ScoringParameters: input preprocessing not of pytom type Preprocessing"
        self.preprocessing = preprocessing
        assert type(fsc_criterion) == float, "ScoringParameters: input fsc_criterion not of type float"
        self.fsc_criterion = fsc_criterion
        assert type(symmetries) == Symmetry or type(symmetries) == type(None)
        self.symmetries = symmetries
        self.weighting = weighting
        self.compoundWedge = compoundWedge

    def toXML(self):
        """
        generate xml node of Scoring Parameters
        @return: xml object of ScoringParameters
        @rtype: L{etree.Element}
        @author: FF
        """
        from lxml import etree

        xmlObj = etree.Element("ScoringParameters")
        xmlObj.append(self.score.toXML())
        if self.reference:
            xmlObj.append(self.reference.toXML(onlyParticleListFilename=True))
        if self.mask:
            xmlObj.append(self.mask.toXML())
        if self.preprocessing:
            xmlObj.append(self.preprocessing.toXML())
        if self.symmetries != None:
            xmlObj.append(self.symmetries.toXML())
        xmlObj.set("Weighting", str(self.weighting))
        xmlObj.set("CompoundWedge", str(self.compoundWedge))
        return xmlObj

    def fromXML(self,xmlObj):
        """
        read ScoringParameters from xml object
        @param xmlObj: xml object
        @type xmlObj: L{lxml.etree._Element}
        @author: FF
        """
        from lxml.etree import _Element
        from pytom.score.score import fromXML as fromXMLScore
        from pytom.basic.structures import Mask, Reference, MultiSymmetries
        from pytom.alignment.preprocessing import Preprocessing

        if xmlObj.__class__ != _Element:
            raise Exception('You must provide a valid XML object.')
        if xmlObj.tag == "ScoringParameters":
            scoringDescription = xmlObj
        else:
            scoringDescription = xmlObj.xpath('SamplingParameters')
            if len(scoringDescription) == 0:
                raise Exception("This XML is not a SamplingParameters.")
            scoringDescription = scoringDescription[0]

        if scoringDescription.get('Weighting'):
            tmp = str(scoringDescription.get('Weighting'))
            if tmp.lower() == 'true':
                self.weighting = True
            else:
                self.weighting = False
        else:
            self.weighting = False

        if scoringDescription.get('CompoundWedge'):
            tmp = str(scoringDescription.get('CompoundWedge'))
            if tmp.lower() == 'true':
                self.compoundWedge = True
            else:
                self.CompoundWedge = False
        else:
            self.CompoundWedge = False

        if scoringDescription.xpath('Score'):
            self.score = fromXMLScore(scoringDescription.xpath('Score')[0])
        else:
            self.score = FLCFScore()

        self.reference = Reference()
        self.reference.fromXML(xmlObj=scoringDescription.xpath('Reference')[0])

        self.mask = Mask()
        if scoringDescription.xpath('Mask'):
            self.mask.fromXML(xmlObj=scoringDescription.xpath('Mask')[0])

        self.preprocessing = Preprocessing()
        if scoringDescription.xpath('Preprocessing'):
            self.preprocessing.fromXML(xmlObj=scoringDescription.xpath('Preprocessing')[0])

        if scoringDescription.xpath('MultiSymmetries'):
            self.symmetries = MultiSymmetries()
            self.symmetries.fromXML(xmlObj=scoringDescription.xpath('MultiSymmetries')[0])
        else:
            self.symmetries = None


class GLocalSamplingJob(PyTomClass):
    """
    Gold standard Local Sampling Job
    """
    def __init__(self, pl=None, ref=None, mask=None, sample_info=None, rotations=None,
                 preprocessing=None, dest='.', max_iter=10, score=FLCFScore(), weighting=False, compoundWedge=False,
                 binning=1, symmetries=None, adaptive_res=0.1, fsc_criterion=0.143, gpuIDs=None):
        """
        @param pl: particle list
        @type pl: L{pytom.basic.structures.ParticleList}
        @param ref: reference
        @type ref: L{pytom.basic.structures.Reference}
        @param mask: mask
        @type mask: L{pytom.basic.structures.Mask}
        @param sample_info: info about sample ?!
        @param rotations: rotations for orientation sampling
        @type rotations: L{pytom.angles.angle.AngleObject}
        @param preprocessing: pre-processing parameters (mostly bandpass filter)
        @type preprocessing: L{pytom.alignment.preprocessing.Preprocessing}
        @param dest: destination folder
        @type dest: L{str}
        @param max_iter: maximum iteration of quasi ExMax alignment
        @type max_iter: L{int}
        @param score: Scoring object
        @type score: L{pytom.score.score.Score}
        @param weighting: Compound wedge weighting of average
        @type weighting: L{bool}
        @param compoundWedge: compound wedge weighting of reference for alignment
        @type compoundWedge: bool
        @param binning: binning FACTOR (1=no binning, 2=2x2x2 voxel-> 1 voxel, etc.)
        @type binning: L{int}
        @param symmetries: apply point symmetry for average
        @type symmetries: ???
        @param adaptive_res: adaptive resolution offset
        @type adaptive_res: L{float}
        @param fsc_criterion: FSC criterion
        @type fsc_criterion: L{float}

        @author: FF
        """
        from pytom.basic.structures import ParticleList, Reference, Mask, MultiSymmetries
        from pytom.score.score import Score
        from pytom.angles.angle import AngleObject
        from pytom.alignment.preprocessing import Preprocessing

        assert (type(pl) == ParticleList) or (type(pl) == type(None)), \
                "GLocalSamplingJob: pl must be particleList or None"
        self.particleList = pl
        assert type(dest) == str, "GLocalSamplingJob: dest must be a string!"
        self.destination = dest
        assert type(max_iter) == int, "GLocalSamplingJob: max_iter must be an integer!"
        self.max_iter = max_iter

        # set scoring parameters
        assert (type(ref) == Reference) or (type(ref) == type(None)), \
                "GLocalSamplingJob: ref must be Reference or None"
        if score == None:
            score=FLCFScore()
        assert type(score) == Score, "GLocalSamplingJob: score is of type Score"
        assert (type(mask) == Mask) or (type(mask) == type(None)), \
                "GLocalSamplingJob: mask is of type Mask"
        if preprocessing==None:
            self.preprocessing = Preprocessing()
        else:
            self.preprocessing = preprocessing
        assert type(preprocessing) == Preprocessing, \
                "GLocalSamplingJob: preprocessing is of type Preprocessing"
        assert type(fsc_criterion) == float, "GLocalSamplingJob: fsc_criterion is a float"
        assert type(weighting) == bool, "GLocalSamplingJob: weighting is bool"
        assert (type(symmetries) == MultiSymmetries) or (type(symmetries) == type(None))
        self.scoringParameters = ScoringParameters(score=score, ref=ref, weighting=weighting,
                                                   compoundWedge=compoundWedge, mask=mask,
                                                   preprocessing=preprocessing, fsc_criterion=fsc_criterion,
                                                   symmetries=symmetries)
        if rotations == None:
            from pytom.angles.localSampling import LocalSampling
            rotations = LocalSampling(shells=3,increment=3., z1Start=0., z2Start=0., xStart=0.)
        assert type(rotations) == AngleObject, "GLocalSamplingJob: rotations is AngleObject"
        assert type(binning) == int, "GLocalSamplingJob: binning is an int"
        assert type(adaptive_res) == float, "GLocalSamplingJob: adaptive_res is of type float"
        self.samplingParameters = SamplingParameters(rotations=rotations,
                 binning=binning, adaptive_res=adaptive_res, sample_info=sample_info)
        self.gpu=gpuIDs

        self.gpu = gpuIDs

    def fromXML(self, xmlObj):
        """
        @param xmlObj: xml object
        @type xmlObj: L{lxml.etree.Element}
        @author: FF
        """
        from lxml.etree import _Element
        from pytom.basic.structures import ParticleList

        if xmlObj.__class__ != _Element:
            raise Exception('You must provide a valid XML object.')
        if xmlObj.tag == "GLocalSamplingJob":
            jobDescription = xmlObj
        else:  
            jobDescription = xmlObj.xpath('GLocalSamplingJob')
            if len(jobDescription) == 0:
                raise Exception("This XML is not a GLocalSamplingJob.")
            jobDescription = jobDescription[0]

        pl = ParticleList('.')
        particleList_element = jobDescription.xpath('ParticleList')[0]
        filename = particleList_element.get('Filename')
        pl.fromXMLFile(filename=filename)

        self.destination = jobDescription.get('Destination')
        self.max_iter = int(jobDescription.get('MaxIterations'))

        try:
            samplingDescription = jobDescription.xpath('SamplingParameters')[0]
        except IndexError:
            raise Exception("GLocalSamplingJob must contain SamplingParameters")
        self.samplingParameters = SamplingParameters()
        self.samplingParameters.fromXML(xmlObj=samplingDescription)

        try:
            scoringDescription = jobDescription.xpath('ScoringParameters')[0]
        except IndexError:
            raise Exception("GLocalSamplingJob must contain ScoringParameters")
        self.scoringParameters = ScoringParameters()
        self.scoringParameters.fromXML(xmlObj=scoringDescription)


    def toXML(self):
        """
        write Job to XML file

        @author: FF
        """
        from lxml import etree

        jobElement = etree.Element("GLocalSamplingJob")
        if self.particleList:
            particleListElement = etree.Element("ParticleList")
            particleListElement.set("Filename", self.particleList._XMLfilename)
            jobElement.append( particleListElement)
            #if not shortParticleList:
            #    jobElement.append(self.particleList.toXML())
            #else:
            #    jobElement.append(self.particleList.toShortXML())
        jobElement.set("Destination", self.destination)
        jobElement.set("MaxIterations", str(self.max_iter))

        jobElement.append(self.scoringParameters.toXML())
        jobElement.append(self.samplingParameters.toXML())
        return jobElement

    def prettyPrintShortXML(self):
        """
        print job with short XML for particle list
        @author: FF
        """
        from lxml import etree
        return etree.tostring(self.toXML(),pretty_print=True)


    def getParticleList(self):
        """
        get particle list
        @return: particle list
        @rtype: L{pytom.basic.structure.ParticleList}
        @author: FF
        """
        return self.particleList

    @property
    def getParticleListtoXML(self):
        """
        @author: FF
        """
        return self.particleList.toXML()
    
    def check(self):
        """
        @author: FF

        """
        from pytom.tools.files import checkDirExists
        self.particleList.check()
        self.reference.check()
        self.mask.check()
        if not checkDirExists(self.destination):
            raise RuntimeError('Destination path not found! ' + self.destination)


def averageParallel(particleList,averageName, showProgressBar=False, verbose=False,
                    createInfoVolumes=False, weighting=None, norm=False,
                    setParticleNodesRatio=3):
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
                                           [weighting]*splitFactor, [norm]*splitFactor)))
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
