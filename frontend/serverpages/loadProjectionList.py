'''
Created on Apr 13, 2012

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
    
    
    from pytom.reconstruction.reconstructionStructures import ProjectionList
    from pytom.tools.files import checkDirExists,checkFileExists 
    from pytom.frontend.serverpages.serverMessages import ErrorMessage
    
    splitParameters = parameters.split('&')
    
    if splitParameters.__class__ == list:
        if verbose:
            print splitParameters
        
        for i in xrange(len(splitParameters)):
            
            parameter = splitParameters[i]
            
            split = parameter.split('=')
        
            keyword = split[0]
            argument = split[1]
            
            if verbose:
                print 'Keyword : ', keyword
                print 'Arguments : ', argument
                
            if keyword == 'XML':
                from pytom.tools.files import readStringFile,getPytomPath
                import StringIO
                from lxml import etree
                
                if not checkFileExists(argument):
                    print ErrorMessage('File / directory not found!')
                    return str(ErrorMessage('File / directory not found!'))
                    
                
                pl = ProjectionList()
                pl.fromXMLFile(argument)
                
                xsltString = readStringFile(getPytomPath() + '/frontend/html/xslt/ProjectionList.xsl')
                xsltTransform = StringIO.StringIO(xsltString)
                
                transformed = pl.xsltTransform(xsltTransform)
                return etree.tostring(transformed,pretty_print=True)
                
            elif keyword == 'DIR':
                from pytom.tools.files import checkDirExists,getPytomPath,readStringFile
                import StringIO
                from lxml import etree
                
                if not checkDirExists(argument):
                    print ErrorMessage('File / directory not found!')
                    return str(ErrorMessage('File / directory not found!'))
                    
                
                pl = ProjectionList()
                pl.loadDirectory(argument)
                
                xsltString = readStringFile(getPytomPath() + '/frontend/html/xslt/ProjectionList.xsl')
                xsltTransform = StringIO.StringIO(xsltString)
                
                transformed = pl.xsltTransform(xsltTransform)
                return etree.tostring(transformed,pretty_print=True)
                
                
            