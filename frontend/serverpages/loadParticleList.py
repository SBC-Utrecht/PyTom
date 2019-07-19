'''
Created on Sep 30, 2011

@author: hrabe
'''


def run(parameters,verbose=False):
    """
    run: Answers a http request depending on the provided parameters.
    Parameters can be
    XML=FILENAME
    DIR=DIRNAME
    ALIG=NAME
    @param parameters: 
    """
    from pytom.basic.structures import ParticleList
    from pytom.tools.files import checkDirExists,checkFileExists
    
    if verbose:
        print "Parsing particleList request!"
        
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
                    raise IOError('File not found!')
                
                pl = ParticleList('/')
                pl.fromXMLFile(argument)
                
                xsltString = readStringFile(getPytomPath() + '/frontend/html/xslt/ParticleList.xsl')
                xsltTransform = StringIO.StringIO(xsltString)
                
                transformed = pl.xsltTransform(xsltTransform)
                
                return etree.tostring(transformed,pretty_print=True)
                
            elif keyword == 'DIR':
                from pytom.tools.files import checkDirExists,getPytomPath,readStringFile
                import StringIO
                from lxml import etree
                
                if not checkDirExists(argument):
                    raise IOError('File not found!')
                
                pl = ParticleList(argument)
                pl.loadDirectory()
                
                xsltString = readStringFile(getPytomPath() + '/frontend/html/xslt/ParticleList.xsl')
                xsltTransform = StringIO.StringIO(xsltString)
                
                transformed = pl.xsltTransform(xsltTransform)
                
                return etree.tostring(transformed,pretty_print=True)
            
            elif keyword == 'ALIG':
                
                from pytom.tools.files import checkDirExists,getPytomPath,readStringFile
                import StringIO
                from lxml import etree
                
                if not checkDirExists(argument):
                    raise IOError('File not found!')
                
                pl = ParticleList(argument)
                pl.fromAlignmentFile(argument)
                
                xsltString = readStringFile(getPytomPath() + '/frontend/html/xslt/ParticleList.xsl')
                xsltTransform = StringIO.StringIO(xsltString)
                
                transformed = pl.xsltTransform(xsltTransform)
                
                return etree.tostring(transformed,pretty_print=True)
            
            
    elif splitParameters.__class__ == str:
        if verbose:
            print splitParameters

