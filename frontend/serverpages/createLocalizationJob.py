'''
Created on Nov 16, 2011

@author: hrabe
'''

def createRunscripts(filename,jobFile, x, y, z):
    """
    createRunscripts:
    """
    from pytom.tools.files import dump2TextFile,getPytomPath
    text  = '#!/bin/bash\n'
    text += 'pytom ' + getPytomPath() + '/bin/localization.py ' + jobFile + ' %d %d %d' % (x, y, z) + '\n'
    
    dump2TextFile(filename, text,append = False)

def run(parameters,verbose=False):
    """
    run: Create an LocalizationJob
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
    from pytom.localization.peak_job import PeakJob
    from pytom.localization.structures import Volume
    from pytom.basic.structures import Mask,Reference,WedgeInfo
    from pytom.angles.globalSampling import GlobalSampling
    from pytom.basic.structures import BandPassFilter
    from pytom.frontend.serverpages.serverMessages import FileMessage
    
    if 'tomo' in parameters:
        vol = Volume(parameters['tomo'])
    else:
        raise RuntimeError('Parameter missing in request!')
    
    if 'ref' in parameters:
        ref = Reference(parameters['ref'])
    else:
        raise RuntimeError('Parameter missing in request!')
    
    if 'mask' in parameters:
        mask = Mask(parameters['mask'])
    else:
        raise RuntimeError('Parameter missing in request!')
    
    if 'angle' in parameters:
        ang = GlobalSampling(parameters['angle'])
    else:
        raise RuntimeError('Parameter missing in request!')
    
    if 'dest' in parameters:
        dest = parameters['dest']
    else:
        dest = './'
    
    if 'low' in parameters:
        low = float(parameters['low'])
    else:
        raise RuntimeError('Parameter missing in request!')
    
    if 'high' in parameters:
        high = float(parameters['high'])
    else:
        raise RuntimeError('Parameter missing in request!')
    
    if 'smooth' in parameters:
        smooth = float(parameters['smooth'])
    else:
        raise RuntimeError('Parameter missing in request!')
    
    if 'w1' in parameters:
        w1 = float(parameters['w1'])
    else:
        raise RuntimeError('Parameter missing in request!')
    
    if 'w2' in parameters:
        w2 = float(parameters['w2'])
    else:
        raise RuntimeError('Parameter missing in request!')
    
    if 'x' in parameters:
        x = float(parameters['x'])
    else:
        x = 0
    
    if 'y' in parameters:
        y = float(parameters['y'])
    else:
        y = 0
    
    if 'z' in parameters:
        z = float(parameters['z'])
    else:
        z = 0
    
    bp = BandPassFilter(low, high, smooth)
    wedge = WedgeInfo([w1, w2])
    
    from pytom.score.score import FLCFScore
    sc = FLCFScore()
    job = PeakJob(vol, ref, mask, wedge, ang, dstDir=dest, score=sc, bandpass=bp)

    jobXMLFile = ''

    if 'jobFile' in parameters:
        jobXMLFile = parameters['jobFile']
        job.toXMLFile(jobXMLFile)
        
        if jobXMLFile[-4:] == '.xml':
            jobRunFile = jobXMLFile[0:-4]
        else:
            jobRunFile = jobXMLFile
            
        createRunscripts(jobRunFile + '.sh', jobXMLFile, x, y, z)
        
    return FileMessage('LocalizationJob',jobXMLFile,'created') 
    
