'''
Created on Apr 25, 2012

@author: hrabe
'''


def createRunscripts(filename,parameters):
    """
    createRunscripts:
    """
    from pytom.tools.files import dump2TextFile,getPytomPath
    text  = '#!/bin/bash\n'
    text += 'pytom ' + getPytomPath() + '/bin/reconstructWB.py ' + parameters
    
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
            print(splitParameters)
        
        for i in range(len(splitParameters)):
            
            parameter = splitParameters[i]
            
            split = parameter.split('=')
        
            keyword = split[0]
            argument = split[1]
            
            parametersDictionary[keyword] = argument
    
    response = ''
    
    try:
        response = interpretRequestParameters(parametersDictionary)
        
    except Exception as err:
        print(err)
        response = ErrorMessage(err)

    return str(response)

    
def interpretRequestParameters(parameters):
    """
    interpretRequestParameters
    """    
    from pytom.reconstruction.reconstructionStructures import ProjectionList
    from pytom.basic.structures import ParticleList
    from pytom.frontend.serverpages.serverMessages import FileMessage
    
    particleList = ParticleList('.')
    
    shellCommand = ''
    
    if 'tomo' in parameters:
        tomogram = parameters['tomo']
        shellCommand += '--tomogram ' + tomogram
    
    if 'plXML' in parameters:
        shellCommand += '--particleList ' + parameters['plXML']  
    
    if shellCommand == '':
        raise RuntimeError('Tomogram or ParticleList parameter missing in request!')
    
    if 'prlDIR' in parameters:
        shellCommand += ' --projectionDirectory ' + parameters['prlDIR']
    elif 'prlXML' in parameters:
        shellCommand += ' --projectionList ' + parameters['prlXML']
    else:
        raise RuntimeError('ProjectionList parameter missing in request!')
    
    if 'ts' in parameters:
        shellCommand += ' --coordinatesScale ' + parameters['ts']
        
    if 'x' in parameters:
        shellCommand += ' --size ' + parameters['x']
    
    if 'xc' in parameters and 'yc' in parameters and 'zc' in parameters:
        shellCommand += ' --recOffset ' + parameters['xc'] + ',' + parameters['yc'] + ','+ parameters['zc']
        
    if 'sb' in parameters:
        shellCommand += ' --projBinning ' + parameters['sb']
        
    # if 'sa' in parameters:
    #     shellCommand += ' --postScale ' + parameters['sa']
        
    if 'jobFile' in parameters:
        jobFile = parameters['jobFile']
    else:
        raise RuntimeError('No job file specified!')
        
    #print shellCommand
    try:
        createRunscripts(jobFile,shellCommand + '\n')
    except Exception as ex:
        print(ex)
    
        
    return FileMessage('ReconstructionJob',jobFile,'created')  

