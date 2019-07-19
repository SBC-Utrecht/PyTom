'''
Created on Dec 3, 2009

@author: hrabe
'''

def getPytomPath():
    """
    getPytomPath: Returns absolute path to pytom dir 
    @return: Path
    """

    import pytom
    
    return pytom.__path__[0]

def listFilesOfTypeInDir(directory,extension):
    """
    listFilesOfTypeInDir: Lists all files of one type from a directory
    @param directory: The directory
    @type directory: str
    @param extension: The file extension searched
    @type extension: str
    @return: A list of files found in that folder. Returns empty list if no file is found.
    """
    import os
    
    if not directory[-1] == os.sep:
        directory += os.sep
        
    return [(directory + f) for f in os.listdir(directory) if f.lower().endswith(extension)]

def readStringFile(filename):
    """
    readStringFile: Will read in a string file and return a string
    @param filename: A filename
    @return: String 
    """
    
    if not checkFileExists(filename):
        raise IOError('File ' + filename + ' not found!')
    
    lines = '' 

    f = open(filename)
    try:
        for line in f:
            lines = lines + line
    except :
        print 'Error reading ' + filename + '!'
        assert False
    finally:
        f.close()

    return lines
    

def checkFileExists(filename):
    """
    checksFileExists: Checks whether file is there or not.
    @param filename: The file
    @type filename: str  
    @return: True if file exists, False otherwise
    """
    import os.path
    if not filename or filename == '':
        return False
    return os.path.isfile(filename)
    
def checkDirExists(dirname):
    """
    checkDirExists: Checks whether directory is in filesystem or not.
    @param dirname: The file
    @type dirname: str  
    @return: True if directory exists, False otherwise
    """
    import os.path
    if not dirname or dirname == '':
        return False
    return os.path.isdir(dirname)

def simulationDescriptionToParticleList(directory,prefix = ''):
    """
    simulationDescriptionToParticleList: 
    """
    lenDir = len(directory)
    
    if not directory[lenDir-1] == '/':
        directory = directory + '/'
        
    xmlFile = directory + 'desc.xml'
    
    from lxml import etree
    
    simulationXML = etree.parse(xmlFile)
    
    #print etree.tostring(simulationXML,pretty_print=True)
    particles = simulationXML.xpath('particle')
    parameters = simulationXML.xpath('Simulation_Parameters')
    
    wedge = int(parameters[0].get('wangleEnd'))/2
    
    from pytom.basic.structures import Particle,ParticleList,WedgeInfo

    wi = WedgeInfo(wedge,[0.0,0.0,0.0])
    
    pl = ParticleList(directory)
    
    for particle in particles:
        
        filename = prefix + particle.get('filename')
        rotation = particle.get('rotation')
        rotation = rotation.replace('[','')
        rotation = rotation.replace(']','')
        rotation = rotation.split(',')
        rotation = [float(rotation[0]),float(rotation[1]),float(rotation[2])]
        
        shift = particle.get('shift')
        shift = shift.replace('[','')
        shift = shift.replace(']','')
        shift = shift.split(',')
        shift = [int(round(float(shift[0]))),int(round(float(shift[1]))),int(round(float(shift[2])))]
        
        p = Particle(filename,rotation = rotation,shift = shift,wedge=wi)
        
        pl.append(p)
    
    #al.toXMLFile(directory+'AlignmentList.xml')
    
    return pl



def generateCorrelationMatrixStartScript(startScriptName, xmlFilename):
    
    
    string = '''#!/usr/bin/env pytom
from pytom.cluster.structures import CorrelationMatrixJob
from pytom.cluster.correlationMatrix import distributedCorrelationMatrix 

exj = CorrelationMatrixJob()
exj.fromXMLFile(\'./''' + xmlFilename+ '''\') 
distributedCorrelationMatrix(exj,False) '''

    import os
    
    command = 'echo "' + string + '" > ' + startScriptName 
    
    os.system(command)

def dump2TextFile(filename,text,append = True):
    """
    dump2TextFile:
    @param filename:
    @param text:  
    @author: Thomas Hrabe
    """
    if append:
        f = open(filename,'a')
    else:
        f = open(filename,'w')
        
    f.write(text)
    f.close()

def dumpMsg2Log(filename,text):
    """
    @deprecated: Use dump2TextFile instead!
    """
    dump2TextFile(filename,text,True)

def writeSpider(volume,filename,z1=0,y=0,z2=0,xOff=0,yOff=0,zOff=0):
    """
    writeSpider: Writes volume to disk according to http://www.wadsworth.org/spider_doc/spider/docs/image_doc.html
    @param volume: The volume
    @param filename: Path and name to written file
    @type filename: str  
    @param z1: Z1 Rotation for iangle parameter in header 
    @param y: Y Rotation for iangle parameter in header
    @param z2: Z2 Rotation for iangle parameter in header
    @param xOff: xOff for offset parameter in header
    @param yOff: yOff for offset parameter in header
    @param zOff: zOff for offset parameter in header
    """
    import struct
    from math import ceil
    fileHandle = open(filename,'wb')
    
    fileHandle.write(struct.pack("f", volume.sizeZ())) #nslice
    fileHandle.write(struct.pack("f", volume.sizeX())) #nrows
    fileHandle.write(struct.pack("f", 0)) #irec
    fileHandle.write(struct.pack("f", 0)) #unused
    
    if volume.sizeZ() > 1:
        fileHandle.write(struct.pack("f", 3)) #iform
    else:
        fileHandle.write(struct.pack("f", 1)) #iform
        
    fileHandle.write(struct.pack("f", 0)) #imami
    fileHandle.write(struct.pack("f", 0)) #fmax
    fileHandle.write(struct.pack("f", 0)) #fmin
    fileHandle.write(struct.pack("f", 0)) #av
    fileHandle.write(struct.pack("f", -1)) #sig
    fileHandle.write(struct.pack("f", 0)) #unused
    fileHandle.write(struct.pack("f", volume.sizeY())) #nsam
    fileHandle.write(struct.pack("f", ceil(256./volume.sizeY()))) #labrec
    
    if z1 != 0 or y != 0 or z2 != 0:
        fileHandle.write(struct.pack("f", 1)) #iangle
    else:
        fileHandle.write(struct.pack("f", 0)) #iangle
    
    #write transformation parameters    
    fileHandle.write(struct.pack("f", z1)) #phi
    fileHandle.write(struct.pack("f", y)) #theta
    fileHandle.write(struct.pack("f", z2)) #gamma
    fileHandle.write(struct.pack("f", xOff)) #xoff
    fileHandle.write(struct.pack("f", yOff)) #yoff
    fileHandle.write(struct.pack("f", zOff)) #zoff
    
    fileHandle.write(struct.pack("f", 0)) #scale
    fileHandle.write(struct.pack("f", ceil(256.0/float(volume.sizeY())) * volume.sizeY() * 4)) #labbyt
    fileHandle.write(struct.pack("f", int(volume.sizeY() *4))) #lenbyt

    fillup = int(ceil(256.0/float(volume.sizeY())) * volume.sizeY() * 4 - fileHandle.tell())
    
    
    for i in xrange(fillup):
        fileHandle.write(struct.pack('b',0)) # fill remaining empty bytes till header is full

    for z in xrange(volume.sizeZ()):
        for y in xrange(volume.sizeY()):
            for x in xrange(volume.sizeX()):
                fileHandle.write(struct.pack("f",volume.getV(x,y,z)))
         
    fileHandle.close()
    
def readSpider(filename):
    """
    readSpider: Reads a spider file and returns a volume
    @param filename: The file name 
    @return: L{pytom_volume.vol}
    """
    
    if not checkFileExists(filename):
        raise IOError('File ' + filename + ' not found!')
    
    from pytom_volume import vol
    import struct
    from math import ceil
    
    fileHandle = open(filename,'rb')
    z = int(struct.unpack('f',fileHandle.read(4))[0]) #nslice
    x = int(struct.unpack('f',fileHandle.read(4))[0]) #nrows
    
    fileHandle.seek(44)
    
    y = int(struct.unpack('f',fileHandle.read(4))[0]) #nsam

    volume = vol(x,y,z)
    
    fileHandle.seek(int(ceil(256.0/float(y)) * y * 4))
    
    for z in xrange(volume.sizeZ()):
        for y in xrange(volume.sizeY()):
            for x in xrange(volume.sizeX()):
                value = float(struct.unpack('f',fileHandle.read(4))[0])
                volume(value,x,y,z)
    
    fileHandle.close()
    
    return volume
            
    
    
def writeMatrix2RMatrix(matrix,filename):
    """
    writeMatrix2RMatrix: Will store the matrix provided as a string file so that R understands it.
    @param matrix: The source matrix
    @param filename: Name of string file matrix written    
    """
    from pytom_volume import vol
    if matrix.__class__ != vol:
        raise TypeError('Parameter provided must be a image of pytom_volume.vol type')
    
    if matrix.sizeZ() > 1:
        raise RuntimeError('Parameter provided must be a N x M x 1 matrix, not a volume!')
    
    fileHandle = open(filename,'w')
    
    for x in xrange(matrix.sizeX()):
        for y in xrange(matrix.sizeY()):
            fileHandle.write(str(matrix(x,y,0)))
            fileHandle.write(' ')
        fileHandle.write('\n')
        
    fileHandle.close()
            
    



