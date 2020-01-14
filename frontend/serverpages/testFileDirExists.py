'''
Created on Jun 29, 2012

@author: hrabe
'''


def run(parameters,verbose=False):
    """
    run: Answers a http request depending on the provided parameters.
    Parameters can be
    XML=FILENAME
    DIR=DIRNAME
    @param parameters: 
    """
    
    from pytom.tools.files import checkDirExists,checkFileExists 
    from pytom.frontend.serverpages.serverMessages import FileMessage,ErrorMessage
    
    
    if splitParameters.__class__ == list:
        if verbose:
            print(splitParameters)
        
        for i in range(len(splitParameters)):
            
            parameter = splitParameters[i]
            
            split = parameter.split('=')
        
            keyword = split[0]
            argument = split[1]
            
            parametersDictionary[keyword] = argument
    
    
    responseMessage = ErrorMessage('File not found')
    
    if 'File' in parametersDictionary:
        fileName = parametersDictionary['File']
        if checkFileExists(fileName):
            responseMessage('FileExistance',fileName,'YES')
            
    if 'Dir' in parametersDictionary:
        dirName = parametersDictionary['Dir']
        if checkFileExists(dirName):
            responseMessage('DirExistance',dirName,'YES')
            
    return str(responseMessage)    
    
        
        
    
    
    
    
    