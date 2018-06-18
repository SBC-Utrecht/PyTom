'''
Created on Sep 30, 2011

@author: hrabe
'''

def createRunscripts(filename,jobFile):
    """
    createRunscripts:
    """
    from pytom.tools.files import dump2TextFile,getPytomPath
    text  = '#!/bin/bash\n'
    text += 'pytom ' + getPytomPath() + '/bin/align.py -j ' + jobFile + '\n'
    
    dump2TextFile(filename, text,append = False)
    
    
def run(parameters,verbose=False):
    """
    run: Create an AlignmentJob
    """
    from pytom.frontend.serverpages.serverMessages import ErrorMessage
    splitParameters = parameters.split('&')
    
    keywords    = []
    arguments   = []
    
    parametersDictionary = {}
    
    if splitParameters.__class__ == list:
        if verbose:
            print splitParameters
        
        for i in xrange(len(splitParameters)):
            
            parameter = splitParameters[i]
            
            split = parameter.split('=')
        
            keyword = split[0]
            argument = split[1]
            
            parametersDictionary[keyword] = argument
    
    response = ''

    try:
        response = interpretRequestParameters(parametersDictionary)
        
    except Exception,err:
        print err
        response = ErrorMessage(err)

    return str(response)

    
def interpretRequestParameters(parameters):
    """
    interpretRequestParameters
    """    
    from pytom.basic.structures import ParticleList,SampleInformation,Reference,Mask,Wedge,PointSymmetry
    from pytom.alignment.preprocessing import Preprocessing 
    from pytom.frontend.serverpages.serverMessages import FileMessage
    
    particleList = ParticleList('.')
    
    if 'plXML' in parameters:
        particleList.fromXMLFile(parameters['plXML'])
    elif 'plDIR' in parameters:
        particleList = ParticleList(parameters['plDIR'])
        particleList.loadDirectory()
    else:
        raise RuntimeError('ParticleList parameter missing in request!')
    
    sampleInfo = SampleInformation()
    if 'pixSize' in parameters:
        sampleInfo.setPixelSize(parameters['pixSize'])
    else:
        raise RuntimeError('Pixelsize parameter missing in request!')
        
    if 'partDia' in parameters:
        sampleInfo.setParticleDiameter(parameters['partDia'])
    else:
        raise RuntimeError('Particle diameter missing in request!')
    
    reference = ''
    
    if 'ref' in parameters:
        reference = Reference(parameters['ref'])
    else:
        raise RuntimeError('Reference parameter missing in request!')
    
    if 'mask' in parameters:
        mask = Mask(parameters['mask'])
    else:
        raise RuntimeError('Mask parameter missing in request!')
    
    angles = None
    
    if 'sampling' in parameters:
        if parameters['sampling'] == 'GLOBAL':
            from pytom.angles.globalSampling import GlobalSampling
            
            if 'angFile' in parameters:
                angles = GlobalSampling(parameters['angFile'])
            else:
                raise RuntimeError('Angle file missing in request!')
        else:
            from pytom.angles.localSampling import LocalSampling
            
            if 'angStart' in parameters:
                startIncrement = int(parameters['angStart'])
            else:
                raise RuntimeError('Angle start missing in request!')
            
            if 'angShells' in parameters:
                shells = int(parameters['angShells'])
            else:
                raise RuntimeError('Angle shells missing in request!')
            
            if 'angInc' in parameters:
                shellIncrement = int(parameters['angInc'])
            else:
                raise RuntimeError('Angle increment missing in request!')
            
            angles = LocalSampling(shells=shells,increment=startIncrement)
    else:
        raise RuntimeError('Sampling completely missing in request!')
    
    
    if not 'lowestF' in parameters:
        raise RuntimeError('Lowest frequency parameter missing in request!')
    
    if not 'highestF' in parameters:
        raise RuntimeError('Highest frequency parameter missing in request!')
    
    if not 'filtSm' in parameters:
        raise RuntimeError('Filter smooth parameter missing in request!')
    
    preprocessing = Preprocessing(float(parameters['lowestF']),float(parameters['highestF']),float(parameters['filtSm']))
    
    adaptive = True
    
    if 'adapt' in parameters:
        fscCriterion=0.5
        resOffset = 0.1
        angleFactor = 0.5
        adaptive = parameters['adapt'] == 'ON'
        
        if parameters['adapt'] == 'ON':
            
            if 'adResC' in parameters:
                fscCriterion = float(parameters['adResC'])
            else:
                raise RuntimeError('Resolution criterion missing in request!')
            
            if 'adResOf' in parameters:
                resOffset = parameters['adResOf']
            else:
                raise RuntimeError('Resolution offset parameter missing in request!')
            
            if 'angFac' in parameters:
                angleFactor = float(parameters['angFac'])
            else:
                raise RuntimeError('Angular factor parameter missing in request!')
        
            
    else:
        raise RuntimeError('Adaptive parameter missing in request!')
    
    score = None
    
    if 'score' in parameters:
        if parameters['score'] == 'xcf':
            from pytom.score.score import xcfScore as scoreClass
        elif parameters['score'] == 'nxcf':
            from pytom.score.score import nxcfScore as scoreClass
        elif parameters['score'] == 'flcf':
            from pytom.score.score import FLCFScore as scoreClass
            
        score = scoreClass()
    else:
        raise RuntimeError('Score parameter missing in request!')
    
    if 'pkPriRad' in parameters:
        radius = float(parameters['pkPriRad'])
    else:
        raise RuntimeError('Peak distance parameter missing in request!')
    
    if 'pkSmooth' in parameters:
        smooth = float(parameters['pkSmooth'])
    else:
        raise RuntimeError('Peak distance smooth missing in request!')
    
    score.setPeakPrior(radius = radius,smooth = smooth)
    
    if 'iter' in parameters:
        iterations = int(parameters['iter'])
    else:
        raise RuntimeError('Number of iterations missing in request!')
    
    if 'binning' in parameters:
        binning = int(parameters['binning'])
    else:
        raise RuntimeError('Scaling parameter missing in request!')
    
    if 'dest' in parameters:
        destination = parameters['dest']
    else:
        raise RuntimeError('Destination parameter missing in request!')
     
    from pytom.alignment.ExMaxAlignment import ExMaxJob
    
    job = ExMaxJob(particleList,destination,reference,score,angles,mask,PointSymmetry(1),1,iterations,preprocessing,-1,binning,sampleInfo,fscCriterion,adaptive,resOffset,angleFactor)
    
    #print job
    jobXMLFile = ''
    
    if 'jobFile' in parameters:
        jobXMLFile = parameters['jobFile']
        job.toXMLFile(jobXMLFile)
        
        if jobXMLFile[-4:] == '.xml':
            jobRunFile = jobXMLFile[0:-4]
        else:
            jobRunFile = jobXMLFile
        createRunscripts(jobRunFile + '.sh',jobXMLFile)
    
    return FileMessage('AlignmentJob',jobXMLFile,'created') 
