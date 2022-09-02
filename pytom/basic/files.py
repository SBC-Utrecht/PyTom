"""
Created on Dec 7, 2010

@author: hrabe
"""
import os

def get_install_folder():
    """
    Get path to pytom install folder.

    @return: full path
    @rtype:  L{str}
    """
    import os
    return os.path.dirname(os.path.dirname(os.popen('which pytom').read()[:-1]))


def readProxy(fileName, subregion1=0, subregion2=0, subregion3=0,
              subregion4=0, subregion5=0, subregion6=0, sampling1=0,
              sampling2=0, sampling3=0, binning1=0, binning2=0, binning3=0):
    """
    readProxy: Proxy function to easily replace pytom_volume calls with read \
    function below
    """

    return read(fileName, subregion=[subregion1, subregion2, subregion3,
                                     subregion4, subregion5, subregion6],
                sampling=[sampling1, sampling2, sampling3],
                binning=[binning1, binning2, binning3])


def read(file, subregion=[0, 0, 0, 0, 0, 0], sampling=[0, 0, 0], binning=[0, 0, 0]):
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

    if not isinstance(file, str):
        raise TypeError('File parameter must be a string!')

    if not checkFileExists(file):
        raise IOError('File not found or path is wrong: ' + file)

    try:
        f = read(file, subregion[0], subregion[1], subregion[2], subregion[3],
                 subregion[4], subregion[5], sampling[0], sampling[1], sampling[2],
                 binning[0], binning[1], binning[2])
        return f
    except Exception as e:
        # redundant to code above, but just in case it goes through
        if "Wrong file format or file doesn't exist!" in e:
            raise IOError('File not found or path is wrong: ' + file)
        else:
            raise


def readSubvolumeFromFourierspaceFile(filename, sizeX, sizeY, sizeZ):
    """
    readSubvolumeFromFourierspaceFile: This function is required when data \
    (in real space) is read in binned mode and a related fourier space file
    like a wedge needs to be read alongside.
    Works only if fourier file is reduced complex without any shift applied.
    @param filename: The fourier space file name
    @param sizeX: X final size of subvolume if it was complete
    (what L{pytom.basic.structures.Wedge.returnWedgeVolume} with
    humanUnderstandable == True returns)
    @param sizeY: Y final size of subvolume if it was complete
    (what L{pytom.basic.structures.Wedge.returnWedgeVolume}
    with humanUnderstandable == True returns)
    @param sizeZ: Z final size of subvolume if it was complete
    (what L{pytom.basic.structures.Wedge.returnWedgeVolume}
    with humanUnderstandable == True returns)
    @return: A subvolume
    @author: Thomas Hrabe
    """
    from pytom_volume import vol, subvolume, paste
    from pytom.basic.fourier import fourierSizeOperation
    [newX, newY, newZ] = fourierSizeOperation(sizeX, sizeY, sizeZ,
                                              reducedToFull=False)
    newVolume = vol(newX, newY, newZ)
    newVolume.setAll(0)
    newX = newX / 2
    newY = newY / 2

    if filename.__class__ == str:
        originalVolume = read(filename)
    elif filename.__class__ == vol:
        # open a backdoor for this function to take volumes, but
        # this should be rather an exception -> not fully documented
        originalVolume = filename
    else:
        raise TypeError('Filename must be a string')

    originalSizeX = int(originalVolume.sizeX())
    originalSizeY = int(originalVolume.sizeY())

    # the original volume is reduced complex without shift ->
    # zero frequency is in outer corner (0,0,0)
    # read a subvolume around every corner with a subvolume
    # of half x,y of the final volume with constant z

    if isinstance(filename, str):
        firstSubvolume = read(filename, subregion=[0, 0, 0, newX, newY, newZ])
    else:
        firstSubvolume = subvolume(filename, 0, 0, 0, newX, newY, newZ)

    if isinstance(filename, str):
        secondSubvolume = read(filename, subregion=[originalSizeX - newX,
                                                    0, 0, newX, newY, newZ])
    else:
        secondSubvolume = subvolume(filename, originalSizeX - newX, 0, 0,
                                    newX, newY, newZ)

    if isinstance(filename, str):
        thirdSubvolume = read(filename, subregion=[0, originalSizeY - newY, 0,
                                                   newX, newY, newZ])
    else:
        thirdSubvolume = subvolume(filename, 0, originalSizeY - newY, 0,
                                   newX, newY, newZ)

    if isinstance(filename, str):
        fourthSubvolume = read(filename, subregion=[originalSizeX - newX,
                                                    originalSizeY - newY, 0,
                                                    newX, newY, newZ])
    else:
        fourthSubvolume = subvolume(filename, originalSizeX - newX,
                                    originalSizeY - newY, 0, newX, newY, newZ)

    # merge the volumes to the final volume
    paste(firstSubvolume, newVolume, 0, 0, 0)
    paste(secondSubvolume, newVolume, newX, 0, 0)
    paste(thirdSubvolume, newVolume, 0, newY, 0)
    paste(fourthSubvolume, newVolume, newX, newY, 0)

    return newVolume


class NaiveAtom:

    def __init__(self, atomSeq, atomType, x, y, z, resSeq, resType, occupancy, tempFact, element):

        self._atomSeq = atomSeq
        self._atomType = atomType
        self._x = x
        self._y = y
        self._z = z
        self._resSeq = resSeq
        self._resType = resType
        self._occupancy = occupancy
        self._tempFact = tempFact
        self._element = element

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

    def getOccupancy(self):
        return self._occupancy

    def getTempFact(self):
        return self._tempFact

    def getElement(self):
        return self._element

    def setX(self, value):
        if value.__class__ == str:
            value = float(value)

        if value.__class__ != int and value.__class__ != float:
            raise TypeError('You must provide an int or float to NaiveAtom.setX()')

        self._x = value

    def setY(self, value):
        if value.__class__ == str:
            value = float(value)

        if value.__class__ != int and value.__class__ != float:
            raise TypeError('You must provide an int or float to NaiveAtom.setY()')

        self._y = value

    def setZ(self, value):
        if value.__class__ == str:
            value = float(value)

        if value.__class__ != int and value.__class__ != float:
            raise TypeError('You must provide an int or float to NaiveAtom.setZ()')

        self._z = value


def naivePDBParser(pdbPath, chainName=None):
    atomList = []

    pdbFile = open(pdbPath, 'r')
    try:
        for line in pdbFile:

            '''
PDB example
ATOM   4366  OXT SER I 456      10.602  32.380  -1.590  1.00 53.05           O
            '''

            '''
CIF example
ATOM   3458 O OXT . GLU B 1 220 ? 28.062  59.037 64.587 1.00 43.20  ? 220 GLU B OXT 1 
            '''

            name = line[:6]

            if name == "ATOM  ":
                chain = line[21]
                if chainName != 'all' and chainName != 'All':
                    if chainName is not None and chain != chainName:
                        continue
                atomdata = line.split()
                if len(atomdata) > 17: # For PDB this is never true...
                    line = atomdata
                    atomSeq  = line[1]
                    atomType = line[3]
                    resType  = line[5]
                    resSeq   = line[8]
                    x        = float(line[10])
                    y        = float(line[11])
                    z        = float(line[12])
                    occupancy = float(line[13])
                    tempFact = float(line[14])
                    element = line[15]

                else:
                    atomSeq  = line[5:11]
                    atomType = line[12:17].strip()
                    resType  = line[17:20]
                    resSeq   = line[22:26]
                    x        = float(line[30:38])
                    y        = float(line[38:46])
                    z        = float(line[46:54])
                    occupancy = float(line[54:60])
                    tempFact = float(line[60:66])
                    element  = line[76:78].strip()

                atomList.append(NaiveAtom(atomSeq, atomType, x, y, z, resSeq, resType, occupancy, tempFact, element))
    finally:
        pdbFile.close()

    return atomList


def mmCIFParser(filePath, chainName=None):
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
                parts = re.sub('\s+', ' ', line).split(' ')
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


def recenterVolume(volume, densityNegative=False):
    from scipy.ndimage import center_of_mass
    from pytom.agnostic.io import read, write
    from pytom.agnostic.tools import paste_in_center
    from pytom.gpu.initialize import xp
    from pytom_numpy import vol2npy
    import os

    try:
        a = vol2npy(volume).copy()
        vol =True
    except:
        a = volume
        vol = False

    if densityNegative:
        a *= -1

    x, y, z = list(map(int, center_of_mass(a)))
    cx, cy, cz = a.shape[0] // 2, a.shape[1] // 2, a.shape[2] // 2

    sx = min(x, a.shape[0] - x)
    sy = min(y, a.shape[0] - y)
    sz = min(z, a.shape[0] - z)

    ac = a[x - sx:x + sx, y - sy:y + sy, z - sz:z + sz]
    b = xp.zeros_like(a)

    b = paste_in_center(ac, b)

    if densityNegative: b *= -1

    if vol:
        write('recenteredDBV21.em', b)
        from pytom.basic.files import read
        vol = read('recenteredDBV21.em')
        os.system('rm recenteredDBV21.em')
        return vol
    else:
        return b


def initSphere(cubeSize, radius, smoothing=0, centerX=None, centerY=None, centerZ=None):
    """
    initSphere: Initilizes a volume with a sphere
    @param cubeSize: The size of the whole volume
    @param radius: Radius of the sphere
    @param smoothing: Smoothing at the edges of the sphere
    @param centerX: Center of shpere along X axis
    @param centerY: Center of shpere along Y axis
    @param centerZ: Center of shpere along Z axis
    """
    from pytom_volume import vol, initSphere

    sphere = vol(cubeSize, cubeSize, cubeSize)
    sphere.setAll(0)

    if centerX is None:
        centerX = cubeSize / 2 - 1

    if centerY is None:
        centerY = cubeSize / 2 - 1

    if centerZ is None:
        centerZ = cubeSize / 2 - 1

    initSphere(sphere, radius, smoothing, 0, centerX, centerY, centerZ)

    return sphere


def cutParticlesFromTomogram(particleList, cubeSize, sourceTomogram, coordinateBinning=0, binningFactorOut=0):
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
        raise RuntimeError(
            'The destination directory ' + destinationDirectory + ' does not exist. Create directory first.')

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
            subcube = read(originTomogram, x - cubeRadius, y - cubeRadius,
                           z - cubeRadius, cubeSize, cubeSize, cubeSize, 0, 0, 0,
                           binningFactorOut, binningFactorOut, binningFactorOut)
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
    n = bytearr[3] * (256 ** 3) + bytearr[2] * (256 ** 2) + bytearr[1] * 256 + bytearr[0]
    return n


class EMHeader():
    def __init__(self):
        self.raw_data = np.zeros(128, dtype='int32')
        self.raw_data[0] = 83886086  # '0x5000006', TODO: hard-coded, to be changed!

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
        self.raw_data[24 + 18] = np.int32(angle * 1000)  # 19th

    def get_tiltangle(self):
        """
        @return: tilt angle in deg
        """
        return self.raw_data[24 + 18] / 1000.

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
        v = read(filename, 0, 0, 0, 0, 0, 0, 0, 0, 0, binning[0], binning[1], binning[2])
    else:
        v = read(filename)

    if binning:
        dim = header.get_dim()
        header.set_dim(dim[0] / binning[0], dim[1] / binning[1], dim[2] / binning[2])

    # files are always casted as float irrespective of original datatype
    # => ensure header is correct
    header.set_1st4bytes(datatype=type(1.), machinetype=None)

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
    data.write(filename)  # write the data first

    if header:  # write the header
        header.set_dim(data.sizeX(), data.sizeY(), data.sizeZ())  # set the dimension
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


# Conversion scripts and helper functions

# These scripts all have same parameters

def mdoc2meta(filename, target, prefix=None, outname=''):
    import numpy
    from pytom.gui.guiFunctions import datatype, headerText as header, fmt
    from pytom.agnostic.io import read, write, read_tilt_angle
    from pytom.tools.files import checkFileExists, checkDirExists
    import os
    if not checkFileExists(filename):
        raise RuntimeError('MDOC file not found! ', filename)

    if not checkDirExists(target):
        raise RuntimeError('Destination directory not found! ', target)


    newFilename = name_to_format(filename, target, "meta") if (not outname is None and not outname) else outname

    mdocEntries = open(filename, 'r').read().split('ZValue')[1:]
    metadata = numpy.zeros((len(mdocEntries)), dtype=datatype)
    metadata['Magnification'] = 1

    datadict = {'TiltAngle': 0, 'Magnification': 79000, 'Intensity': 1.0, 'PixelSpacing': 1.75,
                'DefocusV': 3, 'DefocusU': 3, 'InPlaneRotation': 0, 'FileName': '', 'Voltage': 300, 'MarkerDiameter': 100,
                'SphericalAberration': 2.7, 'AmplitudeContrast': 0.08 }

    for k in datadict.keys():
        metadata[k] = datadict[k]

    queries = ['TiltAngle', 'PixelSpacing', 'Magnification', 'Defocus', 'SubFramePath', 'Voltage']
    metaname = {'TiltAngle':['TiltAngle'], 'PixelSpacing': ['PixelSpacing'], 'Magnification':['Magnification'],
                'Defocus':['DefocusU', 'DefocusV'], 'SubFramePath':['FileName'], 'Voltage': ['Voltage']}

    for n, entry in enumerate(mdocEntries):
        metadata['AcquisitionOrder'][n] = n
        lines = entry.split('\n')
        for line in lines:
            point = [[q, line.split(' ')[-1]] for q in queries if line.startswith(q)]
            for q,p in point:
                for mname in metaname[q]:
                    try:
                        metadata[mname][n] = float(p)
                    except:
                        metadata[mname][n] = p

    metadata['DefocusU'] *= -1
    metadata['DefocusV'] *= -1

    metadata = numpy.sort(metadata, order='TiltAngle')

    for n, fname in enumerate(metadata['FileName']):
        metadata['FileName'][n] =os.path.basename(metadata['FileName'][n].replace('\\', '/'))

    if not prefix is None:
        for n, fname in enumerate(metadata['FileName']):
            metadata['FileName'][n] = f'{prefix}{n:02d}.mrc'

    numpy.savetxt(newFilename, metadata, fmt=fmt, header=header)

def ccp42em(filename, target, prefix='sorted_', outname=None):
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

def ccp42mrc(filename, target, prefix='sorted_', outname=None):
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

def em2mrc(filename, target, prefix='sorted_', outname=None):
    from pytom.agnostic.io import read, write, read_tilt_angle
    from pytom.tools.files import checkFileExists, checkDirExists
    import os

    if not checkFileExists(filename):
        raise RuntimeError('EM file not found! ', filename)

    if not checkDirExists(target):
        raise RuntimeError('Destination directory not found! ', target)

    data = read(filename)
    tilt_angle = read_tilt_angle(filename)

    newFilename = name_to_format(filename, target, "mrc")

    write(newFilename, data, tilt_angle=tilt_angle)

def em2ccp4(filename, target, prefix='sorted_', outname=None):
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

def mrc2ccp4(filename, target, prefix='sorted_', outname=None):
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

def mrc2em(filename, target, prefix='sorted_', outname=None):
    from pytom.agnostic.io import read, write, read_tilt_angle
    from pytom.tools.files import checkFileExists, checkDirExists
    import os
    if not checkFileExists(filename):
        raise RuntimeError('MRC file not found! ', filename)

    if not checkDirExists(target):
        raise RuntimeError('Destination directory not found! ', target)

    data = read(filename)
    tilt_angle = read_tilt_angle(filename)

    newFilename = name_to_format(filename, target, "em")

    write(newFilename, data, tilt_angle=tilt_angle)

def mrcs2mrc(filename, target, prefix='sorted_', outname=None):
    import mrcfile

    data = mrcfile.open(filename, permissive=True).data.copy()

    for sliceID in range(data.shape[0]):
        outname = os.path.join(target, f'{prefix}{sliceID:02d}.mrc')
        mrcfile.new(outname, data[sliceID,:,:].astype('float32'), overwrite=True)

def files2mrcs(folder, target, prefix='sorted_', outname=None):
    import mrcfile
    files = [os.path.join(folder, fname) for fname in os.listdir(folder) if fname.startswith(prefix)]
    data = mrcfile.open(filename, permissive=True).data.copy().squeeze()
    if len(files):
        stack = numpy.zeros((len(files), *data.shape), dtype=numpy.float32)
        for n, filename in enumerate(files):
            stack[n, :, :] = mrcfile.open(filename, permissive=True).data.copy().squeeze()
        mrcfile.new(target, stack, overwrite=True)


def pdb2em(filename, target, prefix='', pixelSize=1, cubeSize=200, chain=None, invertDensity=False, fname='', recenter=True, outname=None):
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

    pdbPath = filename

    atomList = naivePDBParser(pdbPath, chain)

    vol = atomList2em(atomList, pixelSize, cubeSize, invertDensity)

    if recenter:
        vol = recenterVolume(vol, invertDensity)

    if fname:
        vol.write(fname)

    else:
        return vol

def pdb2mrc(pdbPath, pixelSize=1, cubeSize=200, chain=None, invertDensity=False, fname='', recenter=True, outname=None):
    """
    pdb2em: Creates an volume out of a PDB file
    @param pdbPath: Path to PDB file or PDB id for online download
    @param pixelSize: The pixel size to convert to
    @param cubeSize: Resulting cube size
    @return: A volume
    @author: Thomas Hrabe & Luis Kuhn
    """

    vol = pdb2em(pdbPath, pixelSize, cubeSize, chain=chain, invertDensity=invertDensity, recenter=recenter)

    if fname:
        from pytom.agnostic.io import write as writeNPY
        from pytom_numpy import vol2npy

        writeNPY(fname, vol2npy(vol))

    else:
        return vol

def mmCIF2em(mmCIFPath, pixelSize=1, cubeSize=200, chain=None, densityNegative=False, fname='', recenter=True, outname=None):
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

    vol = atomList2em(atomList, pixelSize, cubeSize, densityNegative)

    if recenter:
        vol = recenterVolume(vol, densityNegative)

    if fname:
        vol.write(fname)

    else:
        return vol

def mmCIF2mrc(mmCIFPath, pixelSize=1, cubeSize=200, chain=None, densityNegative=False, fname='', recenter=True, outname=None):
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

    vol = atomList2em(atomList, pixelSize, cubeSize, densityNegative)

    if recenter:
        vol = recenterVolume(vol, densityNegative)

    if fname:
        from pytom.agnostic.io import write as writeNPY
        from pytom_numpy import vol2npy

        writeNPY(fname, vol2npy(vol))

    else:
        return vol


def txt2pl(filename, target, prefix='', outname='', subtomoPrefix=None, wedgeAngle=None):
    from pytom.basic.structures import ParticleList
    particleList_file = outname
    coordinate_file = outname if outname else name_to_format(filename, target, "pl")
    pl = ParticleList()
    pl.loadCoordinateFile( filename=coordinate_file, name_prefix=subtomoPrefix,
        wedgeAngle=wedgeAngle)
    pl.toXMLFile(particleList_file)


def  pl2star(filename, target, prefix='', pixelsize=1., binningPyTom=1., binningWarpM=1., outname='', wedgeAngles=None, prexf='', sorted_folder=''):
    from pytom.agnostic.tools import zxz2zyz
    from pytom.basic.structures import ParticleList
    from pytom.basic.datatypes import RELION31_PICKPOS_STAR, fmtR31S, headerRelion31Subtomo
    import os
    import numpy as np

    # If filename is a directory combine all xmls, otherwise read filename
    if os.path.isdir(filename):
        pl = ParticleList()
        xmls = [os.path.join(filename, fname) for fname in os.listdir(filename) if fname.endswith('.xml')]
        for fname in xmls:
            try:
                tempxml = ParticleList()
                tempxml.fromXMLFile(fname)
                pl = pl + tempxml
            except:
                pass
    else:
        pl = ParticleList()
        pl.fromXMLFile(filename)

    stardata = np.zeros((len(pl)), dtype=RELION31_PICKPOS_STAR)

    for n, p in enumerate(pl):
        x, y, z = p.getPickPosition().toVector()

        stardata['CoordinateX'][n] = x * binningPyTom / binningWarpM
        stardata['CoordinateY'][n] = y * binningPyTom / binningWarpM
        stardata['CoordinateZ'][n] = z * binningPyTom / binningWarpM

        stardata['MicrographName'][n] = p.getPickPosition().getOriginFilename()

        stardata['Magnification'][n] = 1

        stardata['DetectorPixelSize'][n] = pixelsize

        stardata['GroupNumber'][n] = p.getClass()

        z0, z1, x = p.getRotation().toVector()

        z0, y, z1 = zxz2zyz(z0, x, z1)

        stardata['AngleRot'][n] = -z0
        stardata['AngleTilt'][n] = y
        stardata['AnglePsi'][n] = -z1

    newFilename = name_to_format(filename if outname == '' else outname, target, "star")


    np.savetxt(newFilename, stardata, fmt=fmtR31S, header=headerRelion31Subtomo, comments='')

def star2xml(filename, target, prefix='', pixelsize=1., binningPyTom=1., binningWarpM=1., outname='', wedgeAngles=None, prexf='', sorted_folder=''):
    star2pl(filename, target, prefix, pixelsize, binningPyTom, binningWarpM, outname, wedgeAngles)

def  star2pl(filename, target, prefix='', pixelsize=1., binningPyTom=1., binningWarpM=1., outname='', wedgeAngles=None, prexf='', sorted_folder=''):
    from pytom.agnostic.tools import zxz2zyz, zyz2zxz
    from pytom.basic.structures import ParticleList, Particle
    from pytom.basic.datatypes import RELION31_PICKPOS_STAR, fmtR31S, headerRelion31Subtomo
    import os
    import numpy as np

    stardata= loadtxt(filename, dtype=RELION31_PICKPOS_STAR)

    pl = ParticleList()

    for n in range(len(stardata['CoordinateX'])):

        x, y, z = stardata['CoordinateX'][n], stardata['CoordinateX'][n], stardata['CoordinateX'][n]

        p = Particle(stardata['MicrographName'][n])

        factor = binningWarpM / binningPyTom

        p.getShift().setX(x * factor)
        p.getShift().setY(y * factor)
        p.getShift().setZ(z * factor)


        p.setClass(stardata['GroupNumber'][n])

        z0, z1, x = p.getRotation().toVector()

        z0 = -stardata['AngleRot'][n]
        y  = stardata['AngleTilt'][n]
        z1 = -stardata['AnglePsi'][n]


        z0, x, z1 = zyz2zxz(z0, y, z1)

        p.getRotation().setZ1(z0)
        p.getRotation().setZ2(z1)
        p.getRotation().setX(x)
        if not wedgeAngles is None:
            p.getWedge().setWedgeAngles(wedgeAngles)

        pl.append(p)

    newFilename = name_to_format(filename if outname == '' else outname, target, "xml")

    pl.toXMLFile(newFilename)

def log2txt(filename, target, prefix='', pixelsize=1., binningPyTom=1., binningWarpM=1., outname='', wedgeAngles=None, prexf='', sorted_folder=''):
    import numpy as np
    from pytom.basic.datatypes import DATATYPE_TASOLUTION as dtype_ta, DATATYPE_ALIGNMENT_RESULTS_RO, FMT_ALIGNMENT_RESULTS_RO, HEADER_ALIGNMENT_RESULTS_RO
    from pytom.voltools.utils import transform_matrix
    import os

    ta_fname = filename
    shift_fname = prexf
    prefix, filetype, folder = 'sorted_', 'mrc', sorted_folder

    if folder == '': folder =target

    shift = np.loadtxt(shift_fname)
    ta = loadtxt(ta_fname, dtype=dtype_ta, skip_header=3)

    NUM = len(ta)
    ar = np.zeros((NUM), dtype=DATATYPE_ALIGNMENT_RESULTS_RO)
    ar['TiltAngle'] = ta['Tilt']
    ar['InPlaneRotation'] = - ta['Rotation']
    ar['Magnification'] = 1 / ta['Mag']
    ar['OperationOrder'] = "TRS"

    for n, ((shx, shy), mat) in enumerate(zip(shift[:,-2:], shift[:, :4])):
        # fill the rotation matrix with rotation around z-axis
        m = transform_matrix(rotation=(ta['Rotation'][n], 0, 0),
                             rotation_order='rzxz',
                             scale=(ta['Mag'][n], ) * 3)[:2, :2]
        # m[0, 0] = numpy.cos(ar['InPlaneRotation'][n]*numpy.pi/180.)  # why multiply with pp (==0) before??
        # m[0, 1] = numpy.sin(ar['InPlaneRotation'][n]*numpy.pi/180.)
        # m[1, 0] = -numpy.sin(ar['InPlaneRotation'][n]*numpy.pi/180.)
        # m[1, 1] = numpy.cos(ar['InPlaneRotation'][n]*numpy.pi/180.)
        # m[0, 0], m[0, 1], m[1, 0], m[1, 1] = mat[0], mat[2], mat[1], mat[3]

        # dot multiply the rotation with the translation
        # puts the translation in the reference frame??
        # shift[n, -2:] = numpy.dot(m, numpy.array((shx, shy)))
        # Warp does (-shx, -shy) and then dot(m.T, shifts)
        shift_t = np.dot(m.T, np.array((shx, shy)))

        # denom = mat[0] ** 2 + mat[1] ** 2
        # scale = np.sqrt(denom)
        # scaleY = (mat[0] * mat[3] - mat[2] * mat[1]) / scaleX
        # skew = np.rad2deg(np.arctan2(mat[0] * mat[2] + mat[1] * mat[3], denom))
        # skewY = 0
        # ar['InPlaneRotation'][n] = -zxz[2]

        # set output
        # ar['InPlaneRotation'][n] = - np.rad2deg(np.arctan2(mat[1], mat[0]))
        # ar['Magnification'][n] = scale
        ar['AlignmentTransX'][n] = shift_t[0]
        ar['AlignmentTransY'][n] = shift_t[1]

    # set sorted filenames in pytomproject folder
    ar['FileName'] = sorted([os.path.join(folder, fname) for fname in os.listdir(folder) if fname.startswith(prefix) and fname.endswith('.'+filetype)])

    # save alignmentresults text files in the target location
    savetxt(os.path.join(target, 'alignmentResults.txt'), ar, fmt=FMT_ALIGNMENT_RESULTS_RO, header=HEADER_ALIGNMENT_RESULTS_RO)

def txt2wimp(fname, target, prefix, outname=''):
    '''This functions creates and saves a file named outname which has the wimp format.
    @parm fname: path to markerfile.txt file in PyTom DATATYPE_MARKERFILE format
    @param outname: path to output file in wimp format, if empty the name of fname will be used
    '''
    from pytom.basic.datatypes import DATATYPE_MARKERFILE
    import numpy, os

    data = loadtxt(fname, dtype=DATATYPE_MARKERFILE)

    outname = outname if outname else os.path.join(target, fname.split('/')[-1][:-4] + '.wimp')

    ids = numpy.unique(data['MarkerIndex'])

    max_num_obj = num_obj = len(ids)
    num_mode = len(data['MarkerIndex'])

    wimpfile = f'''
     Model file name........................{outname}
     max # of object.......................   {max_num_obj}
     # of node.............................  {num_mode+max_num_obj}
     # of object...........................   {num_obj}
      Object sequence : 
    '''

    headermarker = '''  Object #:           {}
     # of point:           {}
     Display switch:1  247
         #    X       Y       Z      Mark    Label 
    '''

    ii = data['MarkerIndex']
    x = data['PositionX']
    y = data['PositionY']
    a = data['TiltAngle']

    dict_angs = {}
    angs = numpy.sort(numpy.unique(a))

    for l, aaa in enumerate(angs):
        dict_angs[f'{aaa:.2f}'] = l

    line = '{:7d} {:7.2f} {:7.2f} {:7.2f} {:3d}\n'

    cntr = 0

    for n, id in enumerate(ids):
        num_points = len(ii[ii == id])
        wimpfile += headermarker.format(n, num_points)

        xx = x[ii == id]
        yy = y[ii == id]
        aa = a[ii == id]

        for m in range(num_points):
            wimpfile += line.format(cntr, xx[m], yy[m], dict_angs[f'{aa[m]:.2f}'], 0)
            cntr += 1

        cntr += 1

    wimpfile += '''
      END
    '''

    out = open(outname, 'w')
    out.write(wimpfile)
    out.close()

def txt2fid(fname, target, prefix, outname=''):
    import os
    outname = outname if outname else os.path.join(target, fname.split('/')[-1][:-4] + '.wimp')
    txt2wimp(fname,outname=outname)
    os.system('wimp')

def name_to_format(filename, target, extension):
    import os
    basename = os.path.basename(filename)
    return target + ("" if target.endswith(os.sep) else os.sep) + ".".join(basename.split(".")[:-1]) + '.' + extension


def headerline(line):
    if line.startswith('data_') or line.startswith('loop_') or line.startswith('_') or line.startswith('#'):
        return False
    else:
        return True


def loadtxt(filename, dtype='float32', usecols=None, skip_header=0, max_rows=None):
    import numpy
    with open(filename, 'r') as f:
        stop = 1E9 if max_rows is None else max_rows

        lines = [line for n, line in enumerate(f) if n >= skip_header and headerline(line) and n < stop]
        arr = numpy.genfromtxt(lines, dtype=dtype, usecols=usecols, max_rows=max_rows)
    return arr


def loadstar(filename, dtype='float32', usecols=None, skip_header=0, max_rows=None):
    return loadtxt(filename, dtype=dtype, usecols=usecols, skip_header=skip_header, max_rows=max_rows)


def savetxt(filename, arr, header='', fmt='', comments='#'):
    import numpy
    numpy.savetxt(filename, arr, comments=comments, header=header, fmt=fmt)

