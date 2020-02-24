'''
started on Mar 10, 2011

@author: luiskuhn, foerster
'''

from pytom.basic.structures import PyTomClass
import numpy as np
import mrcfile
from pytom_numpy import vol2npy

class Projection(PyTomClass):
    """
    Projection: 
    """
    def __init__(self, filename=None, tiltAngle=None, offsetX=None,offsetY=None,
            alignmentTransX=None,alignmentTransY=None,
            alignmentRotation=None,alignmentMagnification=1.,
            alignedFilename=None, index=1):
        """
        Initialize projection class

        @param filename: The filename
        @type filename: string
        @param tiltAngle: Tilt angle of a projection
        @type tiltAngle: float
        @param offsetX: Offset X is used for subtomogram reconstruction from cut out projections. Takes care of offset due to positions determined in binned tomograms.
        @type offsetX: float
        @param offsetY:
        @type offsetY: float
        @param alignmentTransX: X Shift determined from alignment method
        @type alignmentTransX: float
        @param alignmentTransY: Y Shift determined from alignment method
        @type alignmentTransY: float
        @param alignmentRotation: Rotation determined from fine alignment method
        @type alignmentRotation: float
        @param alignmentMagnification: Magnification determined from fine alignment method      
        @type alignmentMagnification: float
        @param alignedFilename: filename of aligned projection
        @type alignedFilename: string
        @param index: Index of projection in filename (default: 1)
        @type index: int
    
        """
        self._filename = filename
        self._xSize = None
        self._ySize = None
        
        if filename and tiltAngle == None:
            self._loadTiltAngle()
#            print "TiltAngle = "+str(self._tiltAngle)+" deg"
        else:
            self._tiltAngle = tiltAngle

        self.setIndex( index)
        
        if offsetX == None:    
            self._offsetX = 0
        else:
            self._offsetX = offsetX    
        
        if offsetY == None:
            self._offsetY = 0
        else:
            self._offsetY = offsetY
        
        if alignmentTransX:
            self._alignmentTransX = alignmentTransX 
        else:
            self._alignmentTransX = 0.
            
        if alignmentTransY:
            self._alignmentTransY = alignmentTransY
        else:
            self._alignmentTransY = 0.
        
        if alignmentRotation:
            self._alignmentRotation = alignmentRotation
        else:
            self._alignmentRotation = 0.
        
        self._alignmentMagnification = alignmentMagnification

    def info(self):
        tline = (" %3d: " %self.getIndex())
        tline = tline + ("Filename: %20s" %self._filename)
        tline = tline + (", TiltAngle= %6.2f" %self._tiltAngle)
        return tline
        
    def _loadTiltAngle(self):
        """
        _loadTiltAngle: Loads projection EM file and parses header for tilt angle. Assigns it to the object. Works for EM files only!
        """
        import struct

        if self._filename.split('.')[-1] == 'em':
            fh = open(self._filename,'rb')
            fh.seek(168)
            self._tiltAngle = float(struct.unpack('i',fh.read(4))[0])/1000
            fh.close()
        else:
            from pytom.gui.mrcOperations import read_mrc_header
            fh = read_mrc_header(self._filename)
            self._tiltAngle = fh['user'][0]/1000. #TODO find correct location of the tiltangle and autowrite it.

    def getDimensions(self):
        """
        get dimensions of projection

        @rtype: array (3-dim)
        """
        if not self._filename:
            "Projection Filename not set :("
            return None
        from pytom.basic.files import readProxy as read
        
        proj = read(self._filename)
        return proj.sizeX(), proj.sizeY(), proj.sizeZ()

    def getFilename(self):
        return self._filename
        
    def getTiltAngle(self):
        return self._tiltAngle

    def getOffsetX(self):
        return self._offsetX
    
    def getOffsetY(self):
        return self._offsetY
    
    def getAlignmentTransX(self):
        return self._alignmentTransX
    
    def getAlignmentTransY(self):
        return self._alignmentTransY
    
    def getAlignmentRotation(self):
        return self._alignmentRotation
    
    def getAlignmentMagnification(self):
        return self._alignmentMagnification

    def getAlignedFilename(self):
        """
        get filename of aligned projection

        @rtype: string
        """
        return self._alignedFilename

    def setIndex(self, index):
        """
        set index of projection

        @param index: index of projection (as indicated in filename)
        @type index: int
        """
        self._index = index

    def getIndex(self):
        """
        get index of projection

        @rtype: int
        """
        return self._index
    
    def returnProjectionImage(self,rescale=1.):
        """
        @param rescale: rescale image by factor
	    @type rescale: int
        """
        from pytom.basic.files import readProxy as read
        return read(self.getFilename(),0,0,0,0,0,0,0,0,0,rescale,rescale,1)
    
    def setFilename(self, filename):
        """
        set filename of projection

        @param filename: Filename
        @type filename: float
        """
        self._filename = filename
    
    def setTiltAngle(self, tiltAngle):
        """
        set tilt angle in projection

        @param tiltAngle: tilt angle (deg)
        @type tiltAngle: float
        """
        self._tiltAngle = tiltAngle
    
    def setAlignmentTransX(self,value):
        """
        set x-translation in projection

        @param value: translation X
        @type value: float
        """
        self._alignmentTransX = float(value)
    
    def setAlignmentTransY(self,value):
        """
        set y-translation in projection

        @param value: translation Y
        @type value: float
        """
        self._alignmentTransY = float(value)
    
    def setAlignmentRotation(self,value):
        """
        set Rotation value in projection

        @param value: rotation angle (deg)
        @type value: float
        """
        self._alignmentRotation = float(value)
    
    def setAlignmentMagnification(self,value):
        """
        set Magnification value in projection

        @param value: magnification
        @type value: float
        """
        self._alignmentMagnification = float(value)

    def setAlignedFilename(self, alignedFilename):
        """
        set filename of aligned projection

        @type alignedFilename: string
        """
        self._alignedFilename = alignedFilename

    def _setSize(self):
        """
        set X- and Y-size of proj according to file
        """
        from pytom.basic.files import readProxy as read
        p = read(self._filename)
        self._xSize = p.sizeX()
        self._ySize = p.sizeY()
        
    def getXSize(self):
        """
        get X-dimension of proj

        @return: x-dim
        @rtype: float
        """
        if not self._xSize:
            self._setSize()
            
        return self._xSize
    
    def getYSize(self):
        """
        get Y-dimension of proj

        @return: y-dim
        @rtype: float
        """
        if not self._ySize:
            self._setSize()
            
        return self._ySize
    
    def toXML(self):
        """
        print projection metadata into xml
        """
        from lxml import etree
        projection_element = etree.Element('Projection', Filename=str(self._filename),
                                           TiltAngle=str(float(self._tiltAngle)),OffsetX=str(float(self._offsetX)),
                                           OffsetY=str(float(self._offsetY)),
                                           AlignmentTransX=str(float(self._alignmentTransX)),
                                           AlignmentTransY=str(float(self._alignmentTransY)),
                                           AlignmentAngle=str(float(self._alignmentRotation)),
                                           AlignmentMagnification=str(float(self._alignmentMagnification)))
        return projection_element        
        
    def fromXML(self,xmlObj):
        """
        @param xmlObj: A xml object  
        @type xmlObj: L{lxml.etree._Element}
        @author: Thomas Hrabe 
        """
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element :
            raise Exception('Is not a lxml.etree._Element! You must provide a valid XMLobject.')
        
        if xmlObj.tag == 'Projection':
            projection_element = xmlObj
        else:
            Exception('Is not a Projection! You must provide a valid Projection object.')
        
        self._filename = projection_element.get('Filename')
        self._tiltAngle = float(projection_element.get('TiltAngle'))
        self._offsetX = float(projection_element.get('OffsetX'))
        self._offsetY = float(projection_element.get('OffsetY'))
        self._alignmentMagnification = float(projection_element.get('AlignmentMagnification'))
        
        
        if projection_element.get('AlignmentTransX') == [] or projection_element.get('AlignmentTransX') == None:
            self._alignmentTransX = float(projection_element.get('AlignmentOffsetX'))
        else: 
            self._alignmentTransX = float(projection_element.get('AlignmentTransX'))
            
        if projection_element.get('AlignmentTransY') == [] or projection_element.get('AlignmentTransY') == None:
            self._alignmentTransY = float(projection_element.get('AlignmentOffsetY'))
        else: 
            self._alignmentTransY = float(projection_element.get('AlignmentTransY'))
            
        self._alignmentRotation = float(projection_element.get('AlignmentAngle'))
        
        
class ProjectionList(PyTomClass):
    """
    ProjectionList: List of projections
    """
    def __init__(self, list = None):
        self._list = list or []

    def info(self):
        tline = ""
        for projection in self._list:
            tline = tline + str(projection) + "\n"
        return tline

    def loadDirectory(self, directory, metafile='', prefix=''):
        """
        loadDirectory: Will take all projections in a directory, determine their respective tilt angle from the header and populate this list object.
        @param directory: directory where files are located
        @type directory: L{str}
        """
        from pytom.tools.files import checkDirExists
        from pytom.gui.guiFunctions import datatype, loadstar
        import numpy

        if metafile:
            metadata = loadstar(metafile, dtype=datatype)
            tiltAngles = metadata['TiltAngle']
        if not checkDirExists(directory):
            raise RuntimeError('Directory ' + directory + ' does not exist!')

        import os

        if directory[-1] != os.sep:
            # append / to dir
            directory += os.sep

        files = [line for line in sorted(os.listdir(directory)) if prefix in line]
        self.tilt_angles = []

        for n, file in enumerate(files):
            if file[len(file) - 3:len(file)] == '.em' or file[-4:] == '.mrc':
                projection = Projection(directory + file)
                if metafile:
                    N = int(file.split('_')[-1].split('.')[0])
                    projection.setTiltAngle(int(round(tiltAngles[N])))
                self.tilt_angles.append(projection._tiltAngle)
                self.append(projection)

    def sort(self):
        self._list.sort(key=lambda p: p.getTiltAngle())
        
    def append(self, projection):
        """
        append projection to projection list

        @param projection: a projection
        @type projection: L{pytom.reconstructionStructures.Projection}
        """
        self._list.append(projection)
        
    def __getitem__(self,key):
        """
        __getitem__: Retrieve projection at position definded by key
        @param key:
        @type key: int
        @rtype: L{pytom.reconstructionStructures.Projection}. 
        """
        if not isinstance(key, int):
            raise TypeError('Provide an int')
        
        if len(self) <= key:
            raise IndexError('Index out of range.')
        
        return self._list[key]
        
    def __len__(self):
        return len(self._list)
    
    def generateVolumes(self, particles, cubeSize, binning=1, 
            applyWeighting=False, showProgressBar=False, verbose=False,
            preScale=1, postScale=1):
        """
        generateVolumes (old function name, now replaced by reconstructVolumes)
        @deprecated: 
        """
        print('Generate Volumes is deprecated, use reconstructVolumes instead!')
        self.reconstructVolumes(particles, cubeSize, binning, applyWeighting,
	        showProgressBar ,verbose,preScale,postScale)

    def reconstructVolume(self, dims=[512, 512, 128], reconstructionPosition=[0, 0, 0],
                          binning=1, applyWeighting=False, alignResultFile='', gpu=-1):
        """
        reconstruct a single 3D volume from weighted and aligned projections

        @param dims: 3D Size of reconstructed volume
        @type dims: 3d list
        @param binning: down-sizing of projections prior to 3d rec
        @type binning: int
        @param reconstructionPosition: offset center of reconstruction (BEFORE possible binning)
        @type reconstructionPosition: 3-dim list
        @param applyWeighting: apply weighting
        @type applyWeighting: bool
        @param gpu: run reconstruction on gpu
        @type gpu: int
        @author: FF
        last change: include binning in reconstructionPosition
        """
        from pytom_volume import vol, backProject

        vol_bp = vol(dims[0], dims[1], dims[2])
        vol_bp.setAll(0.0)

        # stacks for images, projections angles etc.
        if not alignResultFile:
            [vol_img, vol_phi, vol_the, vol_offsetProjections] = self.toProjectionStack(
                binning=binning, applyWeighting=applyWeighting, showProgressBar=False,
                verbose=False)
        else:
            [vol_img, vol_phi, vol_the, vol_offsetProjections] = self.toProjectionStackFromAlignmentResultFile(
                alignResultFile, binning=binning, applyWeighting=applyWeighting, showProgressBar=False, verbose=False)

        # volume storing reconstruction offset from center (x,y,z)
        recPosVol = vol(3, vol_img.sizeZ(), 1)
        recPosVol.setAll(0.0)
        for iproj in range(0, vol_img.sizeZ()):
            for ii in range(0, 3):
                recPosVol.setV(float(reconstructionPosition[ii] / binning), ii, iproj, 0)

        # finally backproject into volume
        import time
        s = time.time()

        if gpu > -1:
            from pytom.tompy.reconstruction_functions import backProjectGPU as backProject
            import cupy as cp
            import voltools

            interpolation = voltools.Interpolations.FILT_BSPLINE
            projections = vol2npy(vol_img)
            projections = cp.array(projections)
            proj_angles = vol2npy(vol_the)[0, 0, :]
            vol_bp = cp.zeros((vol_bp.sizeX(), vol_bp.sizeY(), vol_bp.sizeZ()), dtype=cp.float32)
            vol_bp = backProject(projections, vol_bp, vol_phi, proj_angles, recPosVol, vol_offsetProjections, interpolation).get()
        else: 
            from pytom_volume import backProject
            backProject(vol_img, vol_bp, vol_phi, vol_the, recPosVol, vol_offsetProjections)
        print(f'backproject time: {time.time()-s}')

        return vol_bp

    def reconstructVolumes(self, particles, cubeSize, binning=1, applyWeighting=False,
                           showProgressBar=False, verbose=False, preScale=1, postScale=1, num_procs=5,
                           alignResultFile='', polishResultFile=''):
        """
        reconstructVolumes: reconstruct a subtomogram given a particle object.

        @param particles: A list of particles
        @type particles: L{pytom.basic.structures.ParticleList}
        @param cubeSize: 3D Size of reconstructed particle (here, x == y == z)
        @param binning: binning of projections
        @param applyWeighting: Will apply weighting for weighted backprojection to projection. Default is off (False)!
        @param showProgressBar: False by default
        @param verbose: False by default
        @param preScale: Scale (bin) the projections BEFORE reconstruction
        @param postScale: Scale the tomogram AFTER reconstruction
        """

        from pytom_volume import vol, paste, backProject, complexRealMult
        from pytom.basic.files import readProxy as read
        from pytom.tools.ProgressBar import FixedProgBar
        from pytom.basic.fourier import fft, ifft
        from pytom.basic.filter import circleFilter, rampFilter, fourierFilterShift
        import time, os
        from pytom_numpy import vol2npy
        from pytom.basic.files import write_em
        from multiprocessing import Process
        print('start')

        if len(particles) == 0:
            raise RuntimeError('ParticleList is empty!')

        if cubeSize.__class__ != int:
            raise TypeError('reconstructVolumes: Parameter cubeSize must be of type int!')

        if showProgressBar:
            progressBar = FixedProgBar(0, len(particles), 'Particle volumes generated ')
            progressBar.update(0)
            numberParticleVolumes = 0

        # stacks for images, projections angles etc.
        if not alignResultFile:
            resultProjstack = self.toProjectionStack(binning=binning, applyWeighting=applyWeighting,
                                                     showProgressBar=False,
                                                     verbose=False, num_procs=num_procs)
        else:
            from pytom.reconstruction.reconstructionFunctions import toProjectionStackFromAlignmentResultsFile
            resultProjstack = toProjectionStackFromAlignmentResultsFile(alignResultFile, weighting=applyWeighting,
                                                                             num_procs=num_procs, binning=binning,
                                                                             circleFilter=True)
        # if verbose: return

        procs = []
        num_finished, num_started = 0, 0
        polishIndices = range(len(self))
        if polishResultFile:
            try:
                from pytom.gui.guiFunctions import LOCAL_ALIGNMENT_RESULTS, loadstar
                import pytom.basic.combine_transformations as ct
                ppf = loadstar(polishResultFile, dtype=LOCAL_ALIGNMENT_RESULTS)
                num = ppf['ParticleIndex'][-1] + 1
                rppf = ppf.reshape(num, ppf.shape[0] // num)
                polishangles = rppf[0]['TiltAngle']
                polishIndices = []
                for angle in self.tilt_angles:
                    for n, pangle in enumerate(polishangles):
                        if abs(pangle-angle) < 0.01:
                            polishIndices.append(n)
                            break

                if len(polishIndices) != len(self):
                    print(self.tilt_angles)
                    print(polishIndices)
                    raise Exception('data from polishfile does not contain the same angles as the angles from tiltimages')

                polishResults = rppf
                print(polishIndices)

            except Exception as e:
                print(e)
                raise Exception('data from polishfile does not contain the same angles as the angles from tiltimages')

        submitted = 0
        if num_procs:
            for n, result in enumerate(resultProjstack):
                if not os.path.exists(os.path.dirname(particles[0].getFilename())):
                    os.mkdir(os.path.dirname(particles[0].getFilename()))
                write_em(os.path.dirname(particles[0].getFilename()) + '/.temp_{}.em'.format(n), result)

            for particleIndex in range(len(particles)):
                p = particles[particleIndex]

                while len(procs) >= num_procs:
                    time.sleep(0.1)
                    procs = [proc for proc in procs if proc.is_alive()]
                    num_finished += num_started - num_finished - len(procs)
                if polishResultFile:
                    particlePolishResults = polishResults[particleIndex][polishIndices]
                else:
                    particlePolishResults = None
                proc = Process(target=self.extract_single_particle, args=(p, particleIndex, verbose, binning, postScale,
                                                                          cubeSize, particlePolishResults))
                procs.append(proc)
                proc.start()
                submitted += 1
                num_started += 1
                time.sleep(0.3)
            # print('Subtomogram reconstruction for particle {} has started (batch mode).'.format(particleIndex))
        else:
            vol_bp = vol(cubeSize, cubeSize, cubeSize)
            for particleIndex in range(len(particles)):
                if polishResultFile:
                    particlePolishResults = polishResults[particleIndex][polishIndices]
                else:
                    particlePolishResults = None

                p = particles[particleIndex]
                # self.extract_single_particle(p, particleIndex, verbose, binning, postScale, cubeSize)
                # print('Subtomogram reconstruction for particle {} has started.'.format(particleIndex))
                if showProgressBar:
                    progressBar.update(particleIndex)

                if verbose:
                    print(p)

                vol_bp.setAll(0.0)

                reconstructionPosition = vol(3, len(self), 1)
                reconstructionPosition.setAll(0.0)

                # adjust coordinates of subvolumes to binned reconstruction
                if particlePolishResults is None:
                    for i in range(num_projections):
                        reconstructionPosition(float(p.getPickPosition().getX() / binning), 0, i, 0)
                        reconstructionPosition(float(p.getPickPosition().getY() / binning), 1, i, 0)
                        reconstructionPosition(float(p.getPickPosition().getZ() / binning), 2, i, 0)
                else:
                    import pytom.basic.combine_transformations as ct
                    x, y, z = float(p.getPickPosition().getX() / binning), \
                              float(p.getPickPosition().getY() / binning), \
                              float(p.getPickPosition().getZ() / binning)

                    for i in range(len(self)):
                        offsetX = polishResults['AlignmentTransX'][i]
                        offsetY = polishResults['AlignmentTransY'][i]

                        pick_position = np.matrix([x, y, z]).T

                        rotation = ct.matrix_rotate_3d_y(self._list[i].getTiltAngle())

                        shift = np.matrix([offsetX, offsetY, 0]).T

                        new_position = rotation * pick_position

                        new_position_shifted = new_position + shift
                        reposition = np.linalg.inv(rotation) * new_position_shifted
                        reposition = np.array(reposition.T)[0]

                        reconstructionPosition(float(reposition[0] / binning), 0, i, 0)
                        reconstructionPosition(float(reposition[1] / binning), 1, i, 0)
                        reconstructionPosition(float(reposition[2] / binning), 2, i, 0)

                if verbose:
                    print((p.getPickPosition().getX() / binning,
                           p.getPickPosition().getY() / binning,
                           p.getPickPosition().getZ() / binning))

                [vol_img, vol_phi, vol_the, vol_offsetProjections] = resultProjstack

                backProject(vol_img, vol_bp, vol_phi, vol_the, reconstructionPosition, vol_offsetProjections)

                if postScale > 1:
                    volumeRescaled = vol(cubeSize / postScale, cubeSize / postScale, cubeSize / postScale)
                    rescaleSpline(vol_bp, volumeRescaled)
                    volumeRescaled.write(p.getFilename())
                else:
                    vol_bp.write(p.getFilename())

        while len(procs):
            time.sleep(0.01)
            procs = [proc for proc in procs if proc.is_alive()]
        # for n, result in enumerate(resultProjstack):

        os.system('rm -f {}'.format(os.path.dirname(particles[0].getFilename()) + '/.temp_*.em'))
        progressBar.update(len(particles))
        print('\n Subtomogram reconstructions have finished.\n\n')

    def extract_single_particle(self, p, pid, verbose, binning, postScale, cubeSize, polishResults=None):
        import os
        from pytom.basic.files import read
        from pytom_volume import vol, backProject, rescaleSpline
        from pytom_numpy import vol2npy
        import numpy as np

        try:
            folder = os.path.dirname(p.getFilename())

            results = []
            for index in range(4):
                results.append(read('{}/.temp_{}.em'.format(folder, index)))

            [vol_img, vol_phi, vol_the, vol_offsetProjections] = results
            num_projections = vol_img.sizeZ()

            if verbose:
                print(p)

            vol_bp = vol(cubeSize, cubeSize, cubeSize)
            vol_bp.setAll(0.0)

            reconstructionPosition = vol(3, num_projections, 1)
            reconstructionPosition.setAll(0.0)

            # adjust coordinates of subvolumes to binned reconstruction
            if polishResults is None:
                for i in range(num_projections):
                    reconstructionPosition(float(p.getPickPosition().getX() / binning), 0, i, 0)
                    reconstructionPosition(float(p.getPickPosition().getY() / binning), 1, i, 0)
                    reconstructionPosition(float(p.getPickPosition().getZ() / binning), 2, i, 0)
            else:
                import pytom.basic.combine_transformations as ct
                x, y, z = float(p.getPickPosition().getX() / binning), \
                          float(p.getPickPosition().getY() / binning), \
                          float(p.getPickPosition().getZ() / binning)

                for i in range(len(self)):
                    offsetX = polishResults['AlignmentTransX'][i]
                    offsetY = polishResults['AlignmentTransY'][i]

                    pick_position = np.matrix([x, y, z]).T

                    rotation = ct.matrix_rotate_3d_y(self._list[i].getTiltAngle())

                    shift = np.matrix([offsetX, offsetY, 0]).T

                    new_position = rotation * pick_position

                    new_position_shifted = new_position + shift
                    reposition = np.linalg.inv(rotation) * new_position_shifted
                    reposition = np.array(reposition.T)[0]

                    reconstructionPosition(float(reposition[0] / binning), 0, i, 0)
                    reconstructionPosition(float(reposition[1] / binning), 1, i, 0)
                    reconstructionPosition(float(reposition[2] / binning), 2, i, 0)

            if verbose:
                print((p.getPickPosition().getX() / binning, p.getPickPosition().getY() / binning,
                       p.getPickPosition().getZ() / binning))

            backProject(vol_img, vol_bp, vol_phi, vol_the, reconstructionPosition, vol_offsetProjections)

            if postScale > 1:
                volumeRescaled = vol(cubeSize / postScale, cubeSize / postScale, cubeSize / postScale)
                rescaleSpline(vol_bp, volumeRescaled)
                volumeRescaled.write(p.getFilename())
            else:
                vol_bp.write(p.getFilename())

            for a in [vol_img, vol_phi, vol_the, vol_offsetProjections]: del a
            del results


        except Exception as e:
            print('Caught exception in worker thread (x = %d):' % pid)
            print()
            raise e

    def saveParticleProjections(self, particles, projectionSize, binning=1,
                                applyWeighting=False, showProgressBar=False, verbose=False, outputScale=1):
        """
        saveParticleProjections: Saves all small projections of particle.

        @param particles:
        @type particles: L{pytom.basic.structures.ParticleList}
        @param projectionSize: x,y size of the resulting small projections
        @param binning: Necessary for positions in projections. Were L{pytom.basic.structures.PickLocation} determined in a binned volume and do you want to reconstruct a unbinned subtomogram? Use libtomc binning notation!
        @param applyWeighting: Will apply weighting for weighted backprojection to projection. Default is off!
        @param showProgressBar:
        @param verbose:
        @param outputScale: Scale the output by a factor. 1 will leave as it is, 2 will scale by 2, 3 by 3, 4 by 4...
        """
        from pytom_volume import vol, complexRealMult, rescaleSpline
        from pytom.basic.files import readProxy as read
        from pytom.reconstruction.reconstructionFunctions import positionInProjection
        from pytom.tools.files import writeSpider, checkDirExists
        from pytom.tools.ProgressBar import FixedProgBar
        from os import mkdir
        from pytom.basic.fourier import fft, ifft
        from pytom.basic.filter import circleFilter, rampFilter, fourierFilterShift
        from pytom.tools.maths import ZRotationMatrix
        from pytom.basic.normalise import mean0std1

        if showProgressBar:
            progressBar = FixedProgBar(0, len(particles) * len(self._list), 'Particle projections generated ')
            progressBar.update(0)
            numberParticleProjections = 0

        tmpImage = read(self._list[0].getFilename())
        imgSizeX = tmpImage.sizeX()
        imgSizeY = tmpImage.sizeY()

        if applyWeighting:
            weightSlice = fourierFilterShift(rampFilter(projectionSize, projectionSize))
            circleFilterRadius = tmpImage.sizeX() / 2
            circleSlice = fourierFilterShift(circleFilter(projectionSize, projectionSize, circleFilterRadius))

        for p in particles:

            particleName = p.getFilename()
            directoryName = particleName[0:len(particleName) - 3] + '/'

            if not checkDirExists(directoryName):
                mkdir(directoryName)

            particleProjectionList = ProjectionList([])
            particleProjectionListSpider = ProjectionList([])

            if verbose:
                print(p)

            vector = [p.getPickPosition().getX() * binning, p.getPickPosition().getY() * binning,
                      p.getPickPosition().getZ() * binning]

            if verbose:
                print('Original position: ', vector)

            for i in range(len(self._list)):

                if verbose:
                    print(self._list[i])

                # position in aligned projection
                # the angle MUST be negative because we want the inverse rotation! The matrix within the function is the forward rotation!
                position2D = positionInProjection(vector, -self._list[i].getTiltAngle(), 'Y')
                if verbose:
                    print('Position in projection: ', position2D)

                # position in unaligned projection

                vector2 = [position2D[0] * 1 / self._list[i].getAlignmentMagnification(),
                           position2D[1] * 1 / self._list[i].getAlignmentMagnification(), 0.0]
                if verbose:
                    print('Position after magnification: ', vector2)

                # position2D was extended to [x,y,0] -> 3D, can now be multiplied to a 3D rotation matrix. only the 2D part is used
                rotationMatrix = ZRotationMatrix(-self._list[i].getAlignmentRotation())

                newPosition = rotationMatrix * vector2

                if verbose:
                    print('Position after inplane rotation: ', newPosition)

                position2D[0] = newPosition[0] - self._list[i].getAlignmentTransX()
                position2D[1] = newPosition[1] - self._list[i].getAlignmentTransY()

                if verbose:
                    print(position2D)

                # catch io error (RuntimeError) exception when cut out piece is outside of big file.
                # in case of error, omit this projection, print warning
                try:
                    # changing origin, after that in upper left corner
                    # x and y are integers!
                    x = (int(position2D[0]) - (projectionSize / 2) + imgSizeX / 2)
                    y = (int(position2D[1]) - (projectionSize / 2) + imgSizeY / 2)

                    if verbose:
                        print('Final position in projection: ', x, y)

                    # determining error from pick center -> offset
                    offsetX = position2D[0] - float(int(position2D[0]))
                    offsetY = position2D[1] - float(int(position2D[1]))

                    if verbose:
                        print('Position offset: ', offsetX, offsetY)

                    projection = read(self._list[i].getFilename(), x, y, 0, projectionSize, projectionSize, 1, 0, 0, 0,
                                      0, 0, 0)

                    if applyWeighting:
                        projection = ifft(complexRealMult(complexRealMult(fft(projection), weightSlice), circleSlice))

                    try:
                        mean0std1(projection)
                    except ValueError as e:
                        print(str(e.message))
                        continue

                    if outputScale > 1:
                        projectionRescaled = vol(projectionSize / outputScale, projectionSize / outputScale, 1)
                        rescaleSpline(projection, projectionRescaled)
                        projection = projectionRescaled

                    projection.write(directoryName + 'projection_' + str(i) + '.em')
                    writeSpider(projection, directoryName + 'projection_' + str(i) + '.spi')

                    particleProjection = Projection(directoryName + 'projection_' + str(i) + '.em',
                                                    self._list[i].getTiltAngle(), offsetX, offsetY,
                                                    self._list[i].getAlignmentTransX(),
                                                    self._list[i].getAlignmentTransY(),
                                                    self._list[i].getAlignmentRotation(),
                                                    self._list[i].getAlignmentMagnification())
                    particleProjectionList.append(particleProjection)

                    particleProjectionSpider = Projection(directoryName + 'projection_' + str(i) + '.spi',
                                                          self._list[i].getTiltAngle(), offsetX, offsetY,
                                                          self._list[i].getAlignmentTransX(),
                                                          self._list[i].getAlignmentTransY(),
                                                          self._list[i].getAlignmentRotation(),
                                                          self._list[i].getAlignmentMagnification())
                    particleProjectionListSpider.append(particleProjectionSpider)

                except RuntimeError:
                    print('Warning : Particle out of bounds (' + str(x) + ',' + str(y) + ') for projection!')
                    print(p)
                    print(self._list[i])
                # assert False

                if showProgressBar:
                    numberParticleProjections = numberParticleProjections + 1
                    progressBar.update(numberParticleProjections)

            particleProjectionList.toXMLFile(directoryName + 'projectionList.xml')
            particleProjectionListSpider.toXMLFile(directoryName + 'projectionListSpider.xml')

    def toXML(self):
        """
        """
        from lxml import etree

        projectionList_element = etree.Element('ProjectionList')

        ##self._list = sorted(self._list, key=lambda Projection: Projection._tiltAngle)

        for projection in self._list:
            projectionList_element.append(projection.toXML())

        return projectionList_element

    def fromXML(self, xmlObj):
        """
        """
        from lxml.etree import _Element

        if xmlObj.__class__ != _Element:
            raise Exception('Is not a lxml.etree._Element! You must provide a valid XMLobject.')

        if xmlObj.tag == 'ProjectionList':
            projectionList_element = xmlObj
        else:
            Exception('Is not a ProjectionList! You must provide a valid ProjectionList object.')

        for projection in projectionList_element:
            newProjection = Projection()
            newProjection.fromXML(projection)
            self.append(newProjection)

        self._list = sorted(self._list, key=lambda Projection: Projection._tiltAngle)

    def toProjectionStack(self, binning=1, applyWeighting=False, tiltAngle=None, showProgressBar=False, verbose=False,
                          num_procs=1):
        """
        toProjectionStack:

        @param binning: binning factor
        @type binning: int
        @param applyWeighting: applyWeighting
        @type applyWeighting: bool
        @param showProgressBar: show pretty bar
        @type showProgressBar: bool
        @param verbose: talkative
        @type verbose: bool

        @return: Will return a stack of projections - [Imagestack,phiStack,thetaStack,offsetStack]
        """
        from pytom_numpy import vol2npy
        from pytom_volume import vol, paste, complexRealMult
        from pytom.basic.files import readProxy as read
        from pytom.tools.ProgressBar import FixedProgBar
        from pytom.basic.fourier import fft, ifft
        from pytom.basic.filter import circleFilter, rampFilter, exactFilter, rotateFilter, fourierFilterShift, \
            fourierFilterShift_ReducedComplex
        from multiprocessing import Process
        import time

        # determine image dimensions according to first image in projection list
        imgDim = read(self._list[0].getFilename(), 0, 0, 0, 0, 0, 0, 0, 0, 0, binning, binning, 1).sizeX()

        stack = vol(imgDim, imgDim, len(self._list))
        stack.setAll(0.0)

        phiStack = vol(1, 1, len(self._list))
        phiStack.setAll(0.0)

        thetaStack = vol(1, 1, len(self._list))
        thetaStack.setAll(0.0)

        offsetStack = vol(1, 2, len(self._list))
        offsetStack.setAll(0.0)

        if int(applyWeighting):
            weightSlice = fourierFilterShift(rampFilter(imgDim, imgDim))

            circleFilterRadius = imgDim // 2
            circleSlice = fourierFilterShift_ReducedComplex(circleFilter(imgDim, imgDim, circleFilterRadius))

        if showProgressBar:
            progressBar = FixedProgBar(0, len(self._particleList), 'Particle volumes generated ')
            progressBar.update(0)

        self.tilt_angles = []

        for (i, projection) in enumerate(self._list):
            self.tilt_angles.append(projection._tiltAngle)

        self.tilt_angles = sorted(self.tilt_angles)

        if num_procs:
            for (i, projection) in enumerate(self._list):

                if verbose:
                    print(projection)

                if int((applyWeighting)) >= 1:

                    if int(applyWeighting) != 2:
                        weightSlice = fourierFilterShift(exactFilter(self.tilt_angles, projection._tiltAngle,
                                                                     imgDim, imgDim, imgDim))

                    else:
                        weightSlice = fourierFilterShift(rotateFilter(self.tilt_angles, projection._tiltAngle,
                                                                      imgDim, imgDim, imgDim))

                image = read(projection.getFilename(), 0, 0, 0, 0, 0, 0, 0, 0, 0, binning, binning, 1)

                if int(applyWeighting):
                    image = ifft(complexRealMult(complexRealMult(fft(image), weightSlice), circleSlice), scaling=True)

                if verbose:
                    print(projection.getTiltAngle(), projection.getOffsetX(), projection.getOffsetY())
                
                thetaStack(int(round(projection.getTiltAngle())), 0, 0, i)
                offsetStack(projection.getOffsetX(), 0, 0, i)
                offsetStack(projection.getOffsetY(), 0, 1, i)
                paste(image, stack, 0, 0, i)

                if showProgressBar:
                    progressBar.update(i)
        else:
            procs = []
            todo = range(len(self._list))
            self.temp_stack = stack
            if applyWeighting: self.temp_weightSlice = weightSlice
            self.temp_circleSlice = circleSlice
            self.temp_thetaStack = thetaStack
            self.temp_offsetStack = offsetStack

            for pid in range(num_procs):
                args = (pid, num_procs, imgDim, applyWeighting, verbose, binning)
                proc = Process(target=self.read_projections, args=args)
                procs.append(proc)
                proc.start()

            while len(procs):
                time.sleep(0.4)
                procs = [proc for proc in procs if proc.is_alive()]
            stack, thetaStack, offsetStack = self.temp_stack, self.temp_thetaStack, self.temp_offsetStack
        return [stack, phiStack, thetaStack, offsetStack]

    def toProjectionStackFromAlignmentResultsFile(self, alignmentResultsFile, weighting=None, lowpassFilter=0.9,
                                                  binning=1, circleFilter=False, num_procs=1):
        """read image and create aligned projection stack, based on the results described in the alignmentResultFile.

           @param alignmentResultsFile: result file generate by the alignment script.
           @type datatypeAR: gui.guiFunction.datatypeAR
           @param weighting: weighting (<0: analytical weighting, >1: exact weighting, 0/None: no weighting )
           @type weighting: float
           @param lowpassFilter: lowpass filter (in Nyquist)
           @type lowpassFilter: float
           @param binning: binning (default: 1 = no binning). binning=2: 2x2 pixels -> 1 pixel, binning=3: 3x3 pixels -> 1 pixel, etc.

           @author: GS
        """
        print('weighting: ', weighting)
        import numpy
        from pytom_numpy import vol2npy
        from pytom.basic.files import read_em, write_em
        from pytom.basic.functions import taper_edges
        from pytom.basic.transformations import general_transform2d
        from pytom.basic.fourier import ifft, fft
        from pytom.basic.filter import filter as filterFunction, bandpassFilter
        from pytom.basic.filter import circleFilter, rampFilter, exactFilter, fourierFilterShift, fourierFilterShift_ReducedComplex
        from pytom_volume import complexRealMult, vol, paste
        import pytom_freqweight
        from pytom.basic.transformations import resize, rotate
        from pytom.gui.guiFunctions import fmtAR, headerAlignmentResults, datatype, datatypeAR, loadstar
        from pytom.reconstruction.reconstructionStructures import Projection, ProjectionList
        from pytom_numpy import vol2npy

        print("Create aligned images from alignResults.txt")

        alignmentResults = loadstar(alignmentResultsFile, dtype=datatypeAR)
        imageList = alignmentResults['FileName']
        tilt_angles = alignmentResults['TiltAngle']

        a = mrcfile.open(imageList[0], permissive=True)
        imdim = a.data.T.shape[0]

        if binning > 1:
            imdim = int(float(imdim) / float(binning) + .5)
        else:
            imdim = imdim

        sliceWidth = imdim

        # pre-determine analytical weighting function and lowpass for speedup
        if (weighting != None) and (float(weighting) < -0.001):
            weightSlice = fourierFilterShift(rampFilter(imdim, imdim))

        if circleFilter:
            circleFilterRadius = imdim // 2
            circleSlice = fourierFilterShift_ReducedComplex(circleFilter(imdim, imdim, circleFilterRadius))
        else:
            circleSlice = vol(imdim, imdim // 2 + 1, 1)
            circleSlice.setAll(1.0)

        # design lowpass filter
        if lowpassFilter:
            if lowpassFilter > 1.:
                lowpassFilter = 1.
                print("Warning: lowpassFilter > 1 - set to 1 (=Nyquist)")
            # weighting filter: arguments: (angle, cutoff radius, dimx, dimy,
            lpf = pytom_freqweight.weight(0.0, lowpassFilter * imdim // 2, imdim, imdim // 2 + 1, 1,
                                          lowpassFilter / 5. * imdim)
            # lpf = bandpassFilter(volume=vol(imdim, imdim,1),lowestFrequency=0,highestFrequency=int(lowpassFilter*imdim/2),
            #                     bpf=None,smooth=lowpassFilter/5.*imdim,fourierOnly=False)[1]

        projectionList = ProjectionList()
        for n, image in enumerate(imageList):
            atx = alignmentResults['AlignmentTransX'][n]
            aty = alignmentResults['AlignmentTransY'][n]
            rot = alignmentResults['InPlaneRotation'][n]
            mag = alignmentResults['Magnification'][n]
            print(imageList[n], tilt_angles[n], atx, aty, rot, mag)
            projection = Projection(imageList[n], tiltAngle=tilt_angles[n], alignmentTransX=atx, alignmentTransY=aty,
                                    alignmentRotation=rot, alignmentMagnification=mag)
            projectionList.append(projection)

        stack = vol(imdim, imdim, len(imageList))
        stack.setAll(0.0)

        phiStack = vol(1, 1, len(imageList))
        phiStack.setAll(0.0)

        thetaStack = vol(1, 1, len(imageList))
        thetaStack.setAll(0.0)

        offsetStack = vol(1, 2, len(imageList))
        offsetStack.setAll(0.0)

        for (ii, projection) in enumerate(projectionList):
            if projection._filename.split('.')[-1] == 'st':
                from pytom.basic.files import EMHeader, read
                idx = projection._index
                image = read(file=projection._filename,
                             subregion=[0, 0, idx - 1, imdim, imdim, 1],
                             sampling=[0, 0, 0], binning=[0, 0, 0])
                if not (binning == 1) or (binning == None):
                    image = resize(volume=image, factor=1 / float(binning))[0]
            else:
                # read projection files
                from pytom.basic.files import EMHeader, read, read_em_header
                image = read(str(projection._filename))
                # image = rotate(image,180.,0.,0.)
                image = resize(volume=image, factor=1 / float(binning))[0]

            if lowpassFilter:
                filtered = filterFunction(volume=image, filterObject=lpf, fourierOnly=False)
                image = filtered[0]

            tiltAngle = projection._tiltAngle

            # normalize to contrast - subtract mean and norm to mean
            immean = vol2npy(image).mean()
            image = (image - immean) / immean

            # smoothen borders to prevent high contrast oscillations
            image = taper_edges(image, imdim // 30)[0]

            # transform projection according to tilt alignment
            transX = projection._alignmentTransX / binning
            transY = projection._alignmentTransY / binning
            rot = float(projection._alignmentRotation)
            mag = float(projection._alignmentMagnification)

            image = general_transform2d(v=image, rot=rot, shift=[transX, transY], scale=mag, order=[2, 1, 0], crop=True)

            # smoothen once more to avoid edges
            image = taper_edges(image, imdim // 30)[0]

            # analytical weighting
            if (weighting != None) and (weighting < 0):
                #image = (ifft(complexRealMult(fft(image), w_func)) / (image.sizeX() * image.sizeY() * image.sizeZ()))
                image = ifft(complexRealMult(complexRealMult(fft(image), weightSlice), circleSlice), scaling=True)


            elif (weighting != None) and (weighting > 0):
                weightSlice = fourierFilterShift(exactFilter(tilt_angles, tiltAngle, imdim, imdim, sliceWidth))
                #image = (ifft(complexRealMult(fft(image), w_func)) / (image.sizeX() * image.sizeY() * image.sizeZ()))
                image = ifft(complexRealMult(complexRealMult(fft(image), weightSlice), circleSlice), scaling=True)

            thetaStack(int(round(projection.getTiltAngle())), 0, 0, ii)
            offsetStack(int(round(projection.getOffsetX())), 0, 0, ii)
            offsetStack(int(round(projection.getOffsetY())), 0, 1, ii)
            paste(image, stack, 0, 0, ii)
            fname = 'sorted_aligned_novel_{:02d}.mrc'.format(ii)
            #write_em(fname.replace('mrc', 'em'), image)
            #mrcfile.new(fname, vol2npy(image).copy().astype('float32').T, overwrite=True)

        return [stack, phiStack, thetaStack, offsetStack]

    def read_projections(self, pid, num_procs, imgDim, applyWeighting, verbose, binning):
        from pytom_volume import vol, paste, complexRealMult
        from pytom.basic.files import readProxy as read
        from pytom.tools.ProgressBar import FixedProgBar
        from pytom.basic.fourier import fft, ifft
        from pytom.basic.filter import circleFilter, rampFilter, exactFilter, fourierFilterShift, \
            fourierFilterShift_ReducedComplex

        reduced_list = self._list[pid::num_procs]

        for (i, projection) in enumerate(reduced_list):
            print(i * num_procs + pid)

            if int(applyWeighting) >= 1:
                print(self.tilt_angles, projection._tiltAngle, imgDim)

                self.temp_weightSlice = fourierFilterShift(exactFilter(self.tilt_angles, projection._tiltAngle,
                                                                       imgDim, imgDim, imgDim))

            image = read(projection.getFilename(), 0, 0, 0, 0, 0, 0, 0, 0, 0, binning, binning, 1)

            if int(applyWeighting):
                image = ifft(complexRealMult(complexRealMult(fft(image), self.temp_weightSlice), self.temp_circleSlice))

            self.temp_thetaStack(int(round(projection.getTiltAngle())), 0, 0, pid + num_procs * i)
            self.temp_offsetStack(projection.getOffsetX(), 0, 0, i * num_procs + pid)
            self.temp_offsetStack(projection.getOffsetY(), 0, 1, i * num_procs + pid)
            paste(image, self.temp_stack, 0, 0, i * num_procs + pid)

    def saveAsProjectionStack(self,filename,scale=1,applyWeighting=False,showProgressBar=False,verbose=False):
        """
        saveAsProjectionStack: Saves this projection list as a stack of projections
        """
        stack = self.toProjectionStack(scale,applyWeighting,showProgressBar,verbose)
        stack.write(filename)

    def setAlignmentParameters(self,alignmentXMLFile,verbose=False):
        """
        setAlignmentParameters:

        @param alignmentXMLFile:
        @param verbose: verbose level
	    @type verbose: bool
        """
        from lxml import etree
        from pytom.tools.files import readStringFile
        
        string = readStringFile(alignmentXMLFile)
        alignXML = etree.fromstring(string)
        
        if verbose:
            print(etree.tostring(alignXML,pretty_print = True))
        
        for projection in self:
            
            projectionName = projection.getFilename()
            
            names = alignXML.xpath('/root/projection/name') 
            
            for name in names:
                id = -1
                
                if str(name.text).find(projectionName) > -1 or projectionName.find(str(name.text)) > -1:
                    id = int(name.get('idx'))
                else:
                    continue
                
                if verbose:
                    print(etree.tostring(name,pretty_print = True))
                
                query = '/root/projection/alpha[@idx=' +str(id)+ ']'
                if verbose:
                    print(query)
                angle = alignXML.xpath(query)
                projection.setAlignmentRotation(float(angle[0].text))
                
                query = '/root/projection/tx[@idx="'+str(id)+'"]'
                if verbose:
                    print(query)
                x = alignXML.xpath(query)
                projection.setAlignmentTransX(float(x[0].text))
                
                query = '/root/projection/ty[@idx="'+str(id)+'"]'
                if verbose:
                    print(query)
                y = alignXML.xpath(query)
                projection.setAlignmentTransY(float(y[0].text))
                
                query = '/root/projection/s[@idx="'+str(id)+'"]'
                if verbose:
                    print(query)
                magnification = alignXML.xpath(query)
                projection.setAlignmentMagnification(float(magnification[0].text))
    
            
class Reconstruction(PyTomClass):            

                
    def __init__(self,projectionList,particleList):
        """
        __init__
        """

        self._projectionList = None
        self._particleList = None

        if projectionList:
            self._projectionList = projectionList
                
        if particleList:
            self._particleList = particleList    
        
    
    def apply(self,cubeSize,preScale=1,postScale=1,target="",coordinatesScale=1):
        """
        apply: Apply reconstruction method
        @param cubeSize: Size of resulting cube
        @param preScale: Downscale projections BEFORE reconstruction. Default = 1. Can be 1x, 2x
        @param postScale:   Downscale projections AFTER reconstruction. Default = 1. Can be 1x, 2x
        @param target: Target volume if no particle list is given.
        @param coordinatesScale: Do the determined coordinates stem from a "binned" volume?  
        """
        
        raise RuntimeError('You must run this function on a child of Reconstruction')


class WeightedBackprojection(PyTomClass):
    """
    WeightedBackprojection: A class that can do everything related to WB
    """          
                
    def apply(self,cubeSize,preScale=1,postScale=1,target="",coordinatesScale=1,applyWeighting=False,showProgressBar=False,verbose=False):
        """
        @param applyWeighting: Are the projections weighted. Apply weighting if not. Is False by default 
        """

        from pytom.tools.ProgressBar import FixedProgBar
        from pytom_volume import vol, backProject
        
        if not self._projectionList or len(self._projectionList) == 0:
            raise RuntimeError('This reconstruction object does not have a valid ProjectionList')
        
        if not self._particleList:
            from pytom.basic.structures import ParticleList,Particle
            
            p = Particle(target)
            pl = ParticleList('.')
            
            self._particleList = pl
            
        targetVolume = vol(cubeSize, cubeSize, cubeSize)
        targetVolume.setAll(0.0)
        
        [projectionStack,phiStack,thetaStack,projectionOffsetStack] = self._projectionStack.toProjectionStack() 

        if showProgressBar:
            progressBar = FixedProgBar(0,len(self._particleList),'Particle volumes generated ')
            progressBar.update(0)
            numberParticleVolumes = 0
        
        
        particleOffsetStack = vol(3, len(self) ,1)
        particleOffsetStack.setAll(0.0)
        
        
        for particleIndex in range(len(self._particleList)):    
            p = self._particleList[particleIndex]
            
            if verbose:
                print(p)
            
            for i in range(len(self)):
                particleOffsetStack.setV(p.getPickPosition().getX()*coordinatesScale, 0, i, 0) 
                particleOffsetStack.setV(p.getPickPosition().getY()*coordinatesScale, 1, i, 0)
                particleOffsetStack.setV(p.getPickPosition().getZ()*coordinatesScale, 2, i, 0)   
            
            backProject(projectionStack, targetVolume, phiStack, thetaStack, projectionOffsetStack,projectionOffsetStack)
            targetVolume.write(p.getFilename())

            if showProgressBar:
                numberParticleVolumes = numberParticleVolumes + 1
                progressBar.update(numberParticleVolumes)
            
        
    def reconstructTomogram(self,cubeSize,preScale=1,postScale=1,target="",coordinatesScale=1,applyWeighting=False,showProgressBar=False,verbose=False):
        """
        reconstructTomogram: Will reconstruct a WB tomogram
        """
        
                
                
        
    
    
    
    
    
