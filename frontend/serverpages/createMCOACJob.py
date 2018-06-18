'''
Created on Dec 7, 2011

@author: hrabe
'''
def createRunscripts(filename,jobFile):
    """
    createRunscripts:
    """
    from pytom.tools.files import dump2TextFile,getPytomPath
    text  = '#!/bin/bash\n'
    text += 'pytom ' + getPytomPath() + '/bin/mcoAC.py -j ' + jobFile + '\n'
    
    dump2TextFile(filename, text,append = False)

def run(parameters,verbose=False):
    
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
    
    from pytom.basic.structures import ParticleList,SampleInformation,Reference,Mask,Wedge
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
    
    if 'wa1' in parameters:
        wedgeAngle1 = float(parameters['wa1'])
    else:
        raise RuntimeError('Wedge angle 1 parameter missing in request!')
    
    if 'wa2' in parameters:
        wedgeAngle2 = float(parameters['wa2'])
    else:
        raise RuntimeError('Wedge angle 2 parameter missing in request!')
    
    wedgeInfo = Wedge([wedgeAngle1,wedgeAngle2])
    
    if 'mask' in parameters:
        mask = Mask(parameters['mask'])
    else:
        raise RuntimeError('Mask parameter missing in request!')
    
    if not 'lowestF' in parameters:
        raise RuntimeError('Lowest frequency parameter missing in request!')
    
    if not 'highestF' in parameters:
        raise RuntimeError('Highest frequency parameter missing in request!')
    
    if not 'filtSm' in parameters:
        raise RuntimeError('Filter smooth parameter missing in request!')
    
    preprocessing = Preprocessing(float(parameters['lowestF']),float(parameters['highestF']),float(parameters['filtSm']))
    
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
    
    if 'binning' in parameters:
        binning = int(parameters['binning'])
    else:
        raise RuntimeError('Scaling parameter missing in request!')
    
    if 'classes' in parameters:
        numberClasses = float(parameters['classes'])
    else:
        raise RuntimeError('Number classes parameter missing in request!')
    
    if 'conv' in parameters:
        convergence = float(parameters['conv'])
    else:
        raise RuntimeError('Convergence parameter missing in request!')
    
    if 'dest' in parameters:
        destination = parameters['dest']
    else:
        raise RuntimeError('Destination parameter missing in request!')
    
    sampleInfo = SampleInformation()
    if 'pixSize' in parameters:
        sampleInfo.setPixelSize(float(parameters['pixSize']))
    else:
        raise RuntimeError('Pixelsize parameter missing in request!')
        
    if 'partDia' in parameters:
        sampleInfo.setParticleDiameter(float(parameters['partDia']))
    else:
        raise RuntimeError('Particle diameter missing in request!')
    
    temperature = None
    
    if 'temp' in parameters:
        temperature = parameters['temp']
        
        if 'stemp' in parameters:
            startTemperature = float(parameters['stemp'])
        else:
            raise RuntimeError('Start temperature parameter missing in request!')
        
        if 'astep' in parameters:
            annealingStep = float(parameters['astep'])
        else:
            raise RuntimeError('Annealing step parameter missing in request!')
        
        from pytom.cluster.mcoACStructures import SigmaTemperature
        
        if temperature == 'sigma':
            temperature = SigmaTemperature(startTemperature,annealingStep) 
        
    else:
        raise RuntimeError('Temperature missing in request!')
    
    
    criterion = None
    
    if 'crit' in parameters:
        from pytom.cluster.mcoACStructures import MetropolisCriterion,ThresholdAcceptance 
        criterion = parameters['crit']
        
        if criterion == 'metropolis':
            criterion = MetropolisCriterion()
        elif criterion == 'threshold':
            criterion = ThresholdAcceptance()
    else:
        raise RuntimeError('Criterion missing in request!')
    
    if 'refin' in parameters:
        localSearchIncrement = float(parameters['refin'])
    else:
        raise RuntimeError('Number of refinement rounds missing in request!')
    
    from pytom.cluster.mcoACStructures import MCOACJob
    
    job = MCOACJob(particleList,destination,mask,score,preprocessing,wedgeInfo,binning,sampleInfo,numberClasses,temperature,criterion,convergence,localSearchIncrement,symmetry = None)
        
    jobXMLFile = ''
    
    if 'jobFile' in parameters:
        jobXMLFile = parameters['jobFile']
        job.toXMLFile(jobXMLFile)
        jobRunFile = jobXMLFile[0:-3] 
        createRunscripts(jobRunFile + 'sh',jobXMLFile)
        
    return FileMessage('MCOACJob',jobXMLFile,'created') 
    
    
 