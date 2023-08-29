'''
Created on Dec 7, 2010

@author: hrabe
'''


def readProxy(fileName, subregion1=0, subregion2=0, subregion3=0, 
              subregion4=0, subregion5=0, subregion6=0, sampling1=0, 
               sampling2=0, sampling3=0, binning1=0, binning2=0,  binning3=0):
    """
    readProxy: Proxy function to easily replace pytom_volume calls with read \
    function below  
    """

    return read(fileName,subregion=[subregion1,subregion2,subregion3, 
                                    subregion4,subregion5,subregion6], 
                         sampling=[sampling1,sampling2,sampling3],  
                         binning=[binning1,binning2,binning3])


def read(file,subregion=[0,0,0,0,0,0],sampling=[0,0,0],binning=[0,0,0]):
    """
    read: Reads a file
    @param file: Path to file. Supports EM, MRC and CCP4 files 
    @type file: str
    @param subregion: Specifies a subregion to be read. The first three 
    values specify the upper left coordinates within the large volume, 
    the last three the length of the subregion along each dimension. 
    @type subregion: List of 6 integers  
    @param sampling: Change read sampling. Read every second (2), 
    third (3) pixel along each dimension.
    @type sampling:  List of 3 integers
    @param binning: Bin volume along each dimension. Note, 1 will do nothing,
    2 will bin with a kernelsize of 2 pixels along each dimension, 3 will bin
    with a kernelsize of 3 pixels along each dimension and so forth. 
    @type binning:  List of 3 integers
    @return: A volume object. L{pytom_volume.vol}
    @author: Thomas Hrabe
    """
    from pytom.tools.files import checkFileExists
    from pytom_volume import read
    
    if not file.__class__ == str:
        raise TypeError('File parameter must be a string!') 
    
    if not checkFileExists(file):
        raise IOError('File not found or path is wrong: ' + file)

    print(subregion)

    try:                                                                                                     
        f = read(file,subregion[0],subregion[1],subregion[2],subregion[3], 
                 subregion[4],subregion[5],sampling[0],sampling[1],sampling[2], 
                 binning[0],binning[1],binning[2])
        return f
    except RuntimeError as e:
        #redundant to code above, but just in case it goes through
        if "Wrong file format or file doesn't exist!" in e.message:
            raise IOError('File not found or path is wrong: ' + file)
        else:
            raise


def readSubvolumeFromFourierspaceFile(filename,size_x,size_y,size_z):
    """
    readSubvolumeFromFourierspaceFile: This function is required when data \
    (in real space) is read in binned mode and a related fourier space file 
    like a wedge needs to be read alongside. 
    Works only if fourier file is reduced complex without any shift applied.      
    @param filename: The fourier space file name
    @param size_x: X final size of subvolume if it was complete 
    (what L{pytom.basic.structures.Wedge.returnWedgeVolume} with 
    humanUnderstandable == True returns)
    @param size_y: Y final size of subvolume if it was complete 
    (what L{pytom.basic.structures.Wedge.returnWedgeVolume} 
    with humanUnderstandable == True returns)
    @param size_z: Z final size of subvolume if it was complete 
    (what L{pytom.basic.structures.Wedge.returnWedgeVolume} 
    with humanUnderstandable == True returns)
    @return: A subvolume 
    @author: Thomas Hrabe
    """
    from pytom_volume import vol,subvolume,paste
    from pytom.basic.fourier import fourierSizeOperation
    [newX,newY,newZ] = fourierSizeOperation(size_x,size_y,size_z, 
                                            reducedToFull = False)
    newVolume = vol(newX,newY,newZ)
    newVolume.setAll(0)
    newX = newX / 2
    newY = newY / 2

    if filename.__class__ == str:
        originalVolume  = read(filename)
    elif filename.__class__ == vol:
        #open a backdoor for this function to take volumes, but 
        #this should be rather an exception -> not fully documented
        originalVolume = filename
    else:
        raise TypeError('Filename must be a string')
        
    originalSizeX   = int(originalVolume.size_x())
    originalSizeY   = int(originalVolume.size_y())
    
    #the original volume is reduced complex without shift -> 
    #zero frequency is in outer corner (0,0,0)
    #read a subvolume around every corner with a subvolume 
    #of half x,y of the final volume with constant z
    
    
    if filename.__class__ == str:
        firstSubvolume  = read(filename,subregion = [0,0,0,newX,newY,newZ])
    else:
        firstSubvolume  = subvolume(filename,0,0,0,newX,newY,newZ)
    
    if filename.__class__ == str:
        secondSubvolume = read(filename,subregion = [originalSizeX - newX, 
                                                     0,0,newX,newY,newZ])
    else:
        secondSubvolume  = subvolume(filename,originalSizeX - newX,0,0,
                                     newX,newY,newZ)

    if filename.__class__ == str:
        thirdSubvolume  = read(filename,subregion = [0,originalSizeY - newY,0, 
                                                     newX,newY,newZ])
    else:
        thirdSubvolume  = subvolume(filename,0,originalSizeY - newY,0, 
                                    newX,newY,newZ)

    if filename.__class__ == str:
        fourthSubvolume = read(filename,subregion = [originalSizeX - newX, 
                                                     originalSizeY - newY,0, 
                                                     newX,newY,newZ])
    else:
        fourthSubvolume  = subvolume(filename,originalSizeX - newX, 
                                     originalSizeY - newY,0,newX,newY,newZ)


    #merge the volumes to the final volume
    paste(firstSubvolume,newVolume,0,0,0)
    paste(secondSubvolume,newVolume,newX,0,0)
    paste(thirdSubvolume,newVolume,0,newY,0)
    paste(fourthSubvolume,newVolume,newX,newY,0)

    return newVolume


class NaiveAtom:

    def __init__(self, atomSeq, atomType, x, y, z, resSeq, resType):
        
        self._atomSeq = atomSeq
        self._atomType = atomType
        self._x = x
        self._y = y
        self._z = z
        self._resSeq = resSeq
        self._resType = resType
        
    def getAtomType(self):
        return self._atomType

    def getAtomSeq(self):
        return self._atomSeq

    def getX(self):
        return self._x

    def getY(self):
        return self._y

    def getZ(self):
        return self._z

    def setX(self,value):
        if value.__class__ == str:
            value = float(value)

        if value.__class__ != int and value.__class__ != float:
            raise TypeError('You must provide an int or float to NaiveAtom.setX()')

        self._x = value

    def setY(self,value):
        if value.__class__ == str:
            value = float(value)

        if value.__class__ != int and value.__class__ != float:
            raise TypeError('You must provide an int or float to NaiveAtom.setY()')
        
        self._y = value

    def setZ(self,value):
        if value.__class__ == str:
            value = float(value)

        if value.__class__ != int and value.__class__ != float:
            raise TypeError('You must provide an int or float to NaiveAtom.setZ()')

        self._z = value


def naivePDBParser(pdbPath,chainName=None):

    atomList = []

    pdbFile = open(pdbPath, 'r')
    try:
        for line in pdbFile:
        
            '''
PDB example
ATOM   4366  OXT SER I 456      10.602  32.380  -1.590  1.00 53.05           O
            '''
        
            name = line[:6]
            
            if name == "ATOM  ":
                chain   = line[21]
                
                if chainName is not None and chain != chainName:
                    continue
                
                atomSeq = line[5:11]
                atomType = line[11:17].strip()
                resType = line[17:20]
                resSeq = line[22:26]
                x = float(line[29:38])
                y = float(line[38:46])
                z = float(line[46:54])
        
                atomList.append(NaiveAtom(atomSeq, atomType, x, y, z, resSeq, resType))    
    finally:
        pdbFile.close()
    
    return atomList


def mmCIFParser(filePath,chainName = None):
    """
    mmCIFParser: Parses mmCIF files from PDB and returns a list of Atom coordinates
    @param filePath: Path to the file
    @param chainName: Focus on one specific chain. Optional, if not specified, the whole file will be used. (This is the oposite to L{pytom.files.naivePDBParser}).
    @return: List of L{pytom.files.NaiveAtom}s. 
    """
    import re
    
    mmCIFFile = open(filePath, 'r')
    lines = mmCIFFile.readlines()
    mmCIFFile.close()
    
    atoms = []
    
    for line in lines:
        try:
            if line[:4] == 'ATOM':
                parts = re.sub( '\s+', ' ', line ).split(' ')
                chain = parts[6]
                
                if chainName is not None and chain != chainName:
                    continue
                
                atomType = parts[2].strip()
                x = float(parts[10])
                y = float(parts[11])
                z = float(parts[12])
                atom = NaiveAtom('', atomType, x, y, z, '', '')
                atoms.append(atom)
                
        except:
            continue
        finally:
            pass

    return atoms


def initSphere(cubeSize,radius,smoothing=0,centerX=None,centerY=None,centerZ=None):
    """
    initSphere: Initilizes a volume with a sphere
    @param cubeSize: The size of the whole volume
    @param radius: Radius of the sphere
    @param smoothing: Smoothing at the edges of the sphere
    @param centerX: Center of shpere along X axis
    @param centerY: Center of shpere along Y axis
    @param centerZ: Center of shpere along Z axis  
    """
    from pytom_volume import vol,initSphere
    
    sphere = vol(cubeSize,cubeSize,cubeSize)
    sphere.setAll(0)
    
    if centerX is None:
        centerX = cubeSize / 2 - 1
    
    if centerY is None:
        centerY = cubeSize / 2 - 1
        
    if centerZ is None:
        centerZ = cubeSize / 2 - 1
        
    initSphere(sphere,radius,smoothing,0,centerX,centerY,centerZ)

    return sphere


def cutParticlesFromTomogram(particleList,cubeSize,sourceTomogram,coordinateBinning = 0,binningFactorOut = 0):
    """
    cutParticlesFromTomogram: Cuts out particles from source tomogram. Source tomograms set for each particle will be used. If they dont exist, sourceTomogram will be set otherwise from Particle.getSource() .
    @param particleList: 
    @type particleList: L{pytom.basic.structures.ParticleList}
    @param cubeSize: Size of cut out cubes
    @type cubeSize: int 
    @param sourceTomogram: The source tomogram (either file name or volume). 
    @type sourceTomogram: L{str} or L{pytom_volume.vol}  
    @param coordinateBinning: binning factor affecting coordinates. was template matching processed on binned data? use fractions (1/2 , 1/4) if tomogram is binned and coordinates are from unbinned...
    @param binningFactorOut: binning factor for final cubes
    @author: Thomas Hrabe
    """
    
    from pytom.tools.files import checkDirExists
    
    destinationDirectory = particleList.getDirectory()
    
    if not checkDirExists(destinationDirectory):
        raise RuntimeError('The destination directory ' + destinationDirectory + ' does not exist. Create directory first.')
    
    from pytom_volume import read 
    
    cubeRadius = cubeSize / 2
    
    for particle in particleList:
        #
        pickPosition = particle.getPickPosition()
        x = int(pickPosition.getX()) 
        y = int(pickPosition.getY())
        z = int(pickPosition.getZ())
        originTomogram = pickPosition.getOriginFilename()
        
        if originTomogram == '':
            originTomogram = sourceTomogram
            
        
        #
        if coordinateBinning > 0:
            x = x * coordinateBinning
            y = y * coordinateBinning
            z = z * coordinateBinning
        
        #
        if not originTomogram or originTomogram == '':
            originTomogram = sourceTomogram
        
        try:
            subcube = read(originTomogram,x - cubeRadius,y - cubeRadius, 
                           z - cubeRadius,cubeSize,cubeSize,cubeSize,0,0,0, 
                           binningFactorOut,binningFactorOut,binningFactorOut)
        except Exception:
            print(particle)
            assert False
        
        subcube.write(particle.getFilename())
        
    
def int32toint8(n):
    """
    @param n: integer
    @type n: int32
    @return: list of int8s
    @rtype: 4-dim list
    """
    bytearr = []
    while n:
        n, d = divmod(n, 256)
        bytearr.append(d)
    return bytearr


import numpy as np
def int8toint32(bytearr):
    """
    @param bytearr: list of bytes
    @rtype: int32
    """
    n = np.int32(0)
    n = bytearr[3]*(256**3) + bytearr[2]*(256**2) + bytearr[1]*256 + bytearr[0]
    return n


class EMHeader():
    def __init__(self):
        self.raw_data = np.zeros(128, dtype='int32')
        self.raw_data[0] = 83886086 # '0x5000006', TODO: hard-coded, to be changed!
    
    def set_dim(self, x, y, z):
        """
        @param x: x-dimension
        @type x: int
        @param y: y-dimension
        @type y: int
        @param z: z-dimension
        @type z: int
        """
        self.raw_data[1] = x
        self.raw_data[2] = y
        self.raw_data[3] = z
    
    def get_dim(self):
        return self.raw_data[1:4]

    def get_1st4bytes(self):
        """
        get first byte of header which specifies endianness and datatype
        """
        n = self.raw_data[0]
        inibytes = int32toint8(n)
        return inibytes

    def set_1st4bytes(self, datatype=None, machinetype=None, verbose=False):
        """
	set 1st 4 bytes of header and alter values for datatype and 
	machinetype if required

	@param datatype: 
	@type datatype: numpy type
	@param machinetype: 
	@type machinetype: string
        """
        inibytes = self.get_1st4bytes()
        if datatype:
            if datatype == float:
                inibytes[3] = np.int8(5)
                if verbose:
                    print("datatype set to " + str(datatype))
            elif datatype == np.int8:
                inibytes[3] = np.int8(1)
                if verbose:
                    print("datatype set to " + str(datatype))
            elif datatype == np.int16:
                inibytes[3] = np.int8(2)
                if verbose:
                    print("datatype set to " + str(datatype))
            elif datatype == np.int32:
                inibytes[3] = np.int8(4)
                if verbose:
                    print("datatype set to " + str(datatype))
            elif datatype == np.complex:
                inibytes[3] = np.int8(8)
                if verbose:
                    print("datatype set to " + str(datatype))
            elif datatype == np.double:
                inibytes[3] = np.int8(9)
                if verbose:
                    print("datatype set to " + str(datatype))
            elif datatype == np.complex64:
                inibytes[3] = np.int8(10)
                if verbose:
                    print("datatype set to " + str(datatype))
            else:
                print("datatype not implemented for EM file")
                print(datatype)
                raise ValueError

        if machinetype:
            if machinetype == 'OS-9':
                inibytes[0] = 0
                if verbose:
                    print("machinetype set to " + machinetype)
            elif machinetype == 'VAX':
                inibytes[0] = 1
                if verbose:
                    print("machinetype set to " + machinetype)
            elif machinetype == 'Convex':
                inibytes[0] = 2
                if verbose:
                    print("machinetype set to " + machinetype)
            elif machinetype == 'SGI':
                inibytes[0] = 3
                if verbose:
                    print("machinetype set to " + machinetype)
            elif machinetype == 'SUN':
                inibytes[0] = 4
                if verbose:
                    print("machinetype set to " + machinetype)
            elif machinetype == 'Mac':
                inibytes[0] = 5
                if verbose:
                    print("machinetype set to " + machinetype)
            elif machinetype == 'PC':
                inibytes[0] = 6
                if verbose:
                    print("machinetype set to " + machinetype)
            else:
                print("machinetype not implemented for EM file")
                print(machinetype)
                raise ValueError

        self.raw_data[0] = int8toint32(inibytes)

    def get_datatype(self):
        """
        get type of data from header
        @return: numpy type
            """
        inibytes = self.get_1st4bytes()
        dtype = inibytes[3]
        if dtype == 1:
            return np.int8
        if dtype == 2:
            return np.int16
        if dtype == 4:
            return np.int32
        if dtype == 5:
            return np.float
        if dtype == 8:
            return np.complex
        if dtype == 9:
            return np.double
        if dtype == 10:
            return np.complex64
        else:
            raise ValueError
    
    def get_machinetype(self):
        """
        get type of data from header
        @return: numpy type
        """
        inibytes = self.get_1st4bytes()
        mtype = inibytes[0]
        if mtype == 0:
            return 'OS-9'
        if mtype == 1:
            return 'VAX'
        if mtype == 2:
            return 'Convex'
        if mtype == 3:
            return 'SGI'
        if mtype == 4:
            return 'SUN'
        if mtype == 5:
            return 'Mac'
        if mtype == 6:
            return 'PC'
        else:
            raise ValueError

    def set_tiltangle(self, angle):
        """
        @param angle: tilt angle in deg
        """
        self.raw_data[24+18] = np.int32(angle*1000) # 19th
    
    def get_tiltangle(self):
        """
        @return: tilt angle in deg
        """
        return self.raw_data[24+18]/1000.
    
    def to_binary(self):
        return self.raw_data.tostring()
    
    def from_binary(self, data):
        self.raw_data = data


def read_em_header(filename):
    """
    read_em_header: Reads the EM header only.
    @param filename: The em file
    @type filename: str
    @returns: L{pytom.basic.files.EMHeader}
    """

    from pytom.tools.files import checkFileExists

    if not checkFileExists(filename):
        raise IOError('readEMHeader: File not found! ' + filename)

    f = open(filename, 'r')
    try:
        header_data = np.fromfile(f, np.dtype('int32'), 128)
        header = EMHeader()
        header.from_binary(header_data)
    finally:
        f.close()

    return header


def read_em(filename, binning=None):
    """Read the whole EM file: the header as well as the data.
    @param filename: filename
    @type filename: str
    @param binning: binning factor in x, y, z (default: 1)
    @type binning: 3-dim array
    @return: [data, header]
    """
    from pytom_volume import read
    
    # read the header
    header = read_em_header(filename)

    # read the data
    if binning:
        v = read(filename, 0,0,0,0,0,0,0,0,0, binning[0], binning[1], binning[2])
    else:
        v = read(filename)
    
    if binning:
        dim = header.get_dim()
        header.set_dim( dim[0]/binning[0], dim[1]/binning[1], dim[2]/binning[2] )

    # files are always casted as float irrespective of original datatype
    # => ensure header is correct
    header.set_1st4bytes( datatype=type(1.), machinetype=None)

    return v, header


def write_em(filename, data, header=None):
    """Write the EM header as well as the data into one file.
    @param filename: filename
    @type filename: str
    @param data: volume data
    @type data: pytom_volume.vol
    @param header: em header
    @type header: EMHeader
    """
    data.write(filename) # write the data first
    
    if header: # write the header
        header.set_dim(data.size_x(), data.size_y(), data.size_z()) # set the dimension
        try:
            f = open(filename, 'rb+')
            f.write(header.to_binary())
        finally:
            f.close()


def write_em_header(filename, header):
    try:
        f = open(filename, 'rb+')
        f.write(header.to_binary())
    finally:
        f.close()


def read_size(filename):
    emfile = filename.endswith('.em')*1
    f = open(filename, 'r')
    try:
        dt_header = np.dtype('int32')
        header = np.fromfile(f, dt_header, 4)
        x = header[0+emfile]
        y = header[1+emfile]
        z = header[2+emfile]
    except:
        raise Exception("reading of MRC file failed")

    f.close()

    return [x,y,z]


# Conversion between fileformats


def atomList2em(atomList, pixelSize, cubeSize, densityNegative=False):
    """
    atomList2em:
    @param atomList:
    @param pixelSize:
    @param cubeSize:
    @param densityNegative:
    @return:
    """
    from math import floor
    from pytom_volume import vol

    if len(atomList) == 0:
        raise RuntimeError('atomList2em : Your atom list is empty!')

    # get map
    volume = vol(cubeSize, cubeSize, cubeSize)
    volume.setAll(0.0)

    centroidX = 0
    centroidY = 0
    centroidZ = 0

    for i in range(len(atomList)):
        centroidX += atomList[i].getX()
        centroidY += atomList[i].getY()
        centroidZ += atomList[i].getZ()

    centroidX = centroidX / len(atomList)
    centroidY = centroidY / len(atomList)
    centroidZ = centroidZ / len(atomList)

    centerX = floor(float(cubeSize) / 2.0)
    centerY = floor(float(cubeSize) / 2.0)
    centerZ = floor(float(cubeSize) / 2.0)

    shiftX = centroidX - centerX
    shiftY = centroidY - centerY
    shiftZ = centroidZ - centerZ

    for i in range(len(atomList)):
        atomList[i].setX(round(atomList[i].getX() / pixelSize) + centerX)
        atomList[i].setY(round(atomList[i].getY() / pixelSize) + centerY)
        atomList[i].setZ(round(atomList[i].getZ() / pixelSize) + centerZ)

    periodicTableAvailable = True
    try:
        # searching for periodic table library http://pypi.python.org/pypi/periodictable
        from periodictable import elements
    except ImportError:
        periodicTableAvailable = False

    for i in range(len(atomList)):
        x = int(atomList[i].getX())
        y = int(atomList[i].getY())
        z = int(atomList[i].getZ())

        if x >= cubeSize or y >= cubeSize or z >= cubeSize:
            raise RuntimeError('Cube size is too small. Please specify a larger cube for PDB structure!')

        currentValue = volume(x, y, z)

        if periodicTableAvailable:
            atomName = atomList[i].getAtomType()[0]
            element = elements.symbol(atomName)
            mass = element.mass
            volume.setV(currentValue + mass, x, y, z)
        else:
            if atomList[i].getAtomType()[0] == 'H':  ##maybe take this out
                volume.setV(currentValue + 1.0, x, y, z)
            elif atomList[i].getAtomType()[0] == 'C':
                volume.setV(currentValue + 6.0, x, y, z)
            elif atomList[i].getAtomType()[0] == 'N':
                volume.setV(currentValue + 7.0, x, y, z)
            elif atomList[i].getAtomType()[0] == 'O':
                volume.setV(currentValue + 8.0, x, y, z)
            elif atomList[i].getAtomType()[0] == 'P':
                volume.setV(currentValue + 15.0, x, y, z)
            elif atomList[i].getAtomType()[0] == 'S':
                volume.setV(currentValue + 16.0, x, y, z)

    if densityNegative:
        volume = volume * -1

    return volume


def pdb2em(pdbPath, pixelSize, cubeSize, chain=None, densityNegative=False):
    """
    pdb2em: Creates an volume out of a PDB file
    @param pdbPath: Path to PDB file or PDB id for online download
    @param pixelSize: The pixel size to convert to
    @param cubeSize: Resulting cube size
    @return: A volume
    @author: Thomas Hrabe & Luis Kuhn
    """
    from math import floor
    from pytom_volume import vol

    atomList = naivePDBParser(pdbPath, chain)

    return atomList2em(atomList, pixelSize, cubeSize, densityNegative)


def mmCIF2em(mmCIFPath, pixelSize, cubeSize, chain=None, densityNegative=False):
    """
    pdb2em: Creates an volume out of a mmCIF file
    @param mmCIFPath: Path to mmCIF file
    @param pixelSize: The pixel size to convert to
    @param cubeSize: Resulting cube size
    @return: A volume
    @author: Thomas Hrabe
    """
    from math import floor
    from pytom_volume import vol

    atomList = mmCIFParser(mmCIFPath, chain)

    return atomList2em(atomList, pixelSize, cubeSize, densityNegative)


def ccp42em(filename, target):
    from pytom_volume import read
    from pytom.tools.files import checkFileExists, checkDirExists
    import os

    if not checkFileExists(filename):
        raise RuntimeError('CCP4 file not found! ', filename)

    if not checkDirExists(target):
        raise RuntimeError('Destination directory not found! ', target)

    emfile = read(filename)

    newFilename = name_to_format(filename, target, "em")

    emfile.write(newFilename, 'em')


def ccp42mrc(filename, target):
    from pytom_volume import read
    from pytom.tools.files import checkFileExists, checkDirExists
    import os

    if not checkFileExists(filename):
        raise RuntimeError('CCP4 file not found! ', filename)

    if not checkDirExists(target):
        raise RuntimeError('Destination directory not found! ', target)

    emfile = read(filename)

    newFilename = name_to_format(filename, target, "mrc")

    emfile.write(newFilename, 'mrc')


def em2mrc(filename, target):
    from pytom_volume import read
    from pytom.tools.files import checkFileExists, checkDirExists
    import os

    if not checkFileExists(filename):
        raise RuntimeError('EM file not found! ', filename)

    if not checkDirExists(target):
        raise RuntimeError('Destination directory not found! ', target)

    emfile = read(filename)

    newFilename = name_to_format(filename, target, "mrc")

    emfile.write(newFilename, 'mrc')


def em2ccp4(filename, target):
    from pytom_volume import read
    from pytom.tools.files import checkFileExists, checkDirExists
    import os

    if not checkFileExists(filename):
        raise RuntimeError('EM file not found! ', filename)

    if not checkDirExists(target):
        raise RuntimeError('Destination directory not found! ', target)

    emfile = read(filename)

    newFilename = name_to_format(filename, target, "ccp4")

    emfile.write(newFilename, 'ccp4')


def mrc2ccp4(filename, target):
    from pytom_volume import read
    from pytom.tools.files import checkFileExists, checkDirExists
    import os

    if not checkFileExists(filename):
        raise RuntimeError('MRC file not found! ', filename)

    if not checkDirExists(target):
        raise RuntimeError('Destination directory not found! ', target)

    emfile = read(filename)

    newFilename = name_to_format(filename, target, "ccp4")

    emfile.write(newFilename, 'ccp4')


def mrc2em(filename, target):
    from pytom_volume import read
    from pytom.tools.files import checkFileExists, checkDirExists
    import os
    if not checkFileExists(filename):
        raise RuntimeError('MRC file not found! ', filename)

    if not checkDirExists(target):
        raise RuntimeError('Destination directory not found! ', target)

    emfile = read(filename)

    newFilename = name_to_format(filename, target, "em")

    emfile.write(newFilename, 'em')


def convertCoords2PL(coordinate_file, particleList_file, subtomoPrefix=None,
        wedge_angle=None):
    pl = ParticleList()
    pl.loadCoordinateFile( filename=coordinate_file, name_prefix=subtomoPrefix,
        wedge_angle=wedge_angle)
    pl.toXMLFile(particleList_file)


def name_to_format(filename, target, extension):
    import os
    basename = os.path.basename(filename)
    return target + ("" if target.endswith(os.sep) else os.sep) + ".".join(basename.split(".")[:-1]) + '.' + extension
