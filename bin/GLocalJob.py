#!/usr/bin/env pytom
'''
Created on Nov 4, 2014

@author: FF
@lastchange: added spherical mask option
'''

if __name__ == '__main__':
    # parse command line arguments

    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.tools.files import checkFileExists,checkDirExists
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1], 
                          description='Create an GLocalSampling job. Documentation is available at\n\
                          http://www.pytom.org/doc/pytom/alignment.html',
                          authors='Friedrich Foerster',
                          options=[ScriptOption(['-p','--particleList'],
                                                'Particle list : xml file storing information to all subvolumes',
                                                arg=True, optional=False),
                                   ScriptOption(['-r','--reference'],
                                                'Reference : the initial reference - if none provided average of particle list',
                                                arg=True, optional=True),
                                   ScriptOption(['-m','--mask'], 'Mask : a mask ', arg=True, optional=False),
                                   ScriptOption(['--SphericalMask'], 'Mask is spherical / speed up!', arg=False,
                                                optional=True),
                                   ScriptOption(['--angleShells'], '# angle shells used for angular refinement. Default= 3',
                                                arg=True, optional=True),
                                   ScriptOption(['--angleIncrement'], 'Angular increment for refinement. Default = 3',
                                                arg=True, optional=True),
                                   ScriptOption(['-s','--scoreObject'], 
                                       'Score Object, i.e., the scoring function. Options: FLCF, NXCF. Default = FLCF',
                                       arg=True, optional=True),
                                   ScriptOption(['--symmetry'], 'PointSymmetry : specify n-fold symmetry (n)',
                                                arg=True, optional=True),
                                   ScriptOption(['--symmetryAngleZ'], 'PointSymmetry axis tilt around Z axis',
                                                arg=True, optional=True),
                                   ScriptOption(['--symmetryAngleX'], 'PointSymmetry axis tilt around X axis',
                                                arg=True, optional=True),
                                   ScriptOption(['-d','--destination'], 'Destination : destination directory',
                                                arg=True, optional=False),
                                   ScriptOption(['-n','--numberIterations'], 'Number of iterations',
                                                arg=True, optional=False),
                                   ScriptOption(['-b','--binning'],
                                                'Perform binning (downscale) of subvolumes by factor. Default=1.',
                                                arg=True, optional=True),
                                   ScriptOption(['--pixelSize'], 'Pixelsize in Angstrom', arg=True, optional=False),
                                   ScriptOption(['--particleDiameter'], 'Particle diameter in Angstrom', arg=True,
                                                optional=False),
                                   ScriptOption(['-w', '--weighting'], 'Weight particles by exp of CC', arg=False,
                                                optional=True),
                                   ScriptOption(['-c', '--compound'], 'Use compound weighting in Fourier space',
                                                arg=False, optional=True),
                                   ScriptOption(['-j','--jobName'], 'Specify job.xml output filename', arg=True,
                                                optional=False),
                                   ScriptOption(['-g','--gpuID'], 'Specify which gpu to use', arg=True,
                                                optional=True),
                                   ScriptOption(['-h', '--help'], 'Help.', arg=False, optional=True)])
    
    if len(sys.argv) <= 2:
        print(helper)
        sys.exit()
    try:
        results = parse_script_options(sys.argv[1:], helper)
        particleList, reference, mask, isSphere, angShells, angleInc, scoreObject, \
        symmetryN, symmetryAxisZ, symmetryAxisX,\
        destination, numberIterations, binning,\
        pixelSize, diameter, weighting, compound, jobName, gpuIDs, help = results
    except Exception as e:
        print(e)
        sys.exit()
        
    if help is True:
        print(helper)
        sys.exit()

    from pytom.alignment.GLocalSampling import GLocalSamplingJob, mainAlignmentLoop
    from pytom.basic.structures import ParticleList, Reference, Mask, SampleInformation, PointSymmetry
    from pytom.score.score import FLCFScore, nxcfScore
    from pytom.angles.localSampling import LocalSampling
    from pytom.alignment.preprocessing import Preprocessing

    #particleList
    if not checkFileExists(particleList):
        raise RuntimeError('ParticleList file ' + particleList + ' does not exist!')
    pl = ParticleList()
    pl.fromXMLFile(particleList)

    if reference:
        if not checkFileExists(reference):
            raise RuntimeError('Reference file ' + reference + ' does not exist!')
        ref = Reference(referenceFile=reference)
    else:
        ref = Reference()
    
    if not checkFileExists(mask):
        raise RuntimeError('Mask file ' + mask + ' does not exist!')
    if isSphere:
        isSphere = True
    else:
        isSphere = False
    m = Mask(filename=mask, isSphere=isSphere)

    if not checkDirExists(destination):
        raise RuntimeError('Destination directory ' + destination + ' does not exist!')

    if not angShells:
        angShells = 3
    if not angleInc:
        angleInc = 3.
    rot = LocalSampling(angShells,angleInc)

    if not pixelSize:
        pixelSize = 1.
    if not diameter:
        diameter = -1
    sampleI = SampleInformation(pixelSize=float(pixelSize), particleDiameter=float(diameter))

    if symmetryN is None or symmetryAxisZ is None or symmetryAxisX is None:
        sym = None
    else:
        sym = PointSymmetry(nfold=int(symmetryN),z2=float(symmetryAxisZ),x=float(symmetryAxisX))

    if not binning:
        binning = 1

    if weighting:
        weighting = True
    else:
        weighting = False

    if compound:
        compound = True
    else:
        compound = False

    gpuIDs = [] if gpuIDs is None else list(map(int, gpuIDs.split(',')))

    ################# fixed parameters #################
    if (scoreObject == None or scoreObject.lower() == 'flcf'):
        score   = FLCFScore()
    elif scoreObject.lower() == 'nxcf':
        score   = nxcfScore()
    else:
        score   = FLCFScore()
    score.setRemoveAutocorrelation(flag=False)
    #pre     = Preprocessing(lowestFrequency = float(lowestFrequency),highestFrequency = float(highestFrequency))
    pre     = Preprocessing()
    adaptive_res = 0.1
    fsc_criterion = 0.143

    locJob = GLocalSamplingJob(pl=pl, ref=ref, mask=m,
                               sample_info=sampleI, rotations=rot,
                               preprocessing=pre,
                               dest=destination, max_iter=int(numberIterations), score=score, binning=int(binning),
                               weighting=weighting, compoundWedge=compound, gpuIDs=gpuIDs,
                               symmetries=None, adaptive_res=adaptive_res, fsc_criterion=fsc_criterion)
    locJob.toXMLFile(jobName)
    # run script
    mainAlignmentLoop( alignmentJob=locJob, verbose=False)
    print('finished')
    
