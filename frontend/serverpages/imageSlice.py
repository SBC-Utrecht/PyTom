'''
Created on Sep 30, 2011

@author: hrabe
'''

def run(parameters,verbose = False):  
    """
    run: Generate an image slice from a density file 
    """
    
    from pytom.tools.files import checkFileExists
    
    if verbose:
        print "Parsing image slice!"
        
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
            
            if keyword == 'File':
                
                from pytom_volume import read
                from pytom.tools.toImage import volumeToPNG
                if not checkFileExists(argument):
                    raise IOError('File not found!')
                
                v = read(argument)
                
                if argument[len(argument)-3:len(argument)] == '.em':
                    pngFilename = argument[0:len(argument)-3]+'.png'
                elif argument[len(argument)-4:len(argument)] in ['.mrc']:
                    pngFilename = argument[0:len(argument)-4]+'.png'
                elif argument[len(argument)-5:len(argument)] in ['.ccp4']:
                    pngFilename = argument[0:len(argument)-5]+'.png'
                    
                volumeToPNG(v,pngFilename)
                
                return pngFilename