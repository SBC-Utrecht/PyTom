'''
Created on 22.06.2012

@author: thomas
'''

'''
Created on Sep 30, 2011

@author: hrabe
'''

def run(parameters,verbose = False):  
    """
    run: Generate an image slice from a density file 
    """
    
    from pytom_volume import read
    from pytom.tools.toImage import volumeToPNG
    from pytom.tools.files import checkFileExists
    from pytom.frontend.serverpages.serverMessage import DataMessage
    
    if verbose:
        print "Parsing image slice!"
        
    splitParameters = parameters.split('&')
    
    if splitParameters.__class__ == list:
        if verbose:
            print splitParameters
        filename = None
        v = None
        
        sliceX = None
        sliceY = None
        sliceZ = None
        
        for i in xrange(len(splitParameters)):
            
            parameter = splitParameters[i]
            
            split = parameter.split('=')
        
            keyword = split[0]
            argument = split[1]
            
            if verbose:
                print 'Keyword : ', keyword
                print 'Arguments : ', argument
            
            if keyword == 'File':
                filename = argument
                
            if keyword == 'SliceX':
                sliceX = argument
            if keyword == 'SliceY':
                sliceY = argument
            if keyword == 'SliceZ':
                sliceZ = argument
                
        if not checkFileExists(filename):
            raise IOError('File not found!')
         
        v = read(filename)
        data = ''
        
        for x in xrange(v.sizeX()):
            for y in xrange(v.sizeX()):
                for z in xrange(v.sizeX()):
                    data += str(v(x,y,z)) + ';'
                    
        return str(DataMessage('volume',data,sizeX = str(v.sizeX),sizeY = str(v.sizeY),sizeZ = str(v.sizeZ)))

