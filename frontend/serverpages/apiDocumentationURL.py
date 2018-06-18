'''
Created on Sep 30, 2011

@author: hrabe
'''

def run(parameters,verbose = True):
    """
    run: Determines Path / URL to documentation files
    """
    
    from pytom.tools.files import getPytomPath

    pytomPath = getPytomPath()
    
    apiDocPath = pytomPath + '/doc/index.html'
    
    if verbose:
        print 'API DOC PATH : ', apiDocPath
    
    return apiDocPath