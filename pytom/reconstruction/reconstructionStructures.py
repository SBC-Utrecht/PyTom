'''
started on Mar 10, 2011

@author: luiskuhn, foerster, Marten Chaillet
'''
import sys
from pytom.basic.structures import PyTomClass, PyTomClassError
from pytom.agnostic.tools import convert_operation_order_str2list, convert_operation_order_list2str

# TODO  - what should be the default operation order? => default order seems to be [0, 1, 2]
# TODO this was in previous versions already available as parameter to basic.transformations -> general_transform()
# TODO  - projectionlist should also be able to read projections from a stacked mrc file

# ============================= helper function parallel cpu alignment ====================

def init_pool_processes(shared_stack_):
    """
    Make a shared array inheritable for all the workers in multiprocessing Pool.
    @param shared_stack_:
    @type shared_stack_:
    @return:
    @rtype:
    """
    global shared_stack
    shared_stack = shared_stack_


def put_in_shared_stack(stack_dtype, stack_shape, image, index):
    """
    For accessing shared arrays in python multiprocessing.

    @param stack_dtype:
    @type stack_dtype:
    @param stack_shape:
    @type stack_shape:
    @param image:
    @type image:
    @param index:
    @type index:
    @return:
    @rtype:
    """
    import numpy as np
    from pytom.lib.pytom_numpy import vol2npy
    # create a np view to global shared_stack
    _shared_stack_np = np.frombuffer(shared_stack.get_obj(),
                                                       dtype=stack_dtype).reshape(stack_shape, order='F')
    _shared_stack_np[:, :, index] = vol2npy(image).copy(order='F')


def align_single_projection(index, projection_list, tilt_angles, shared_dtype, shared_stack_shape, weighting, binning,
                            low_pass_freq, apply_circle_filter, scale_factor_particle, angle_specimen, verbose,
                            show_progress, imdim, imdimX, imdimY, slice_width, offset):
    """
    For parallel alignment of tilt series.

    @param index: index of a projection in the list
    @type index: L{int}
    @param weighting: ramp weighting (== -1), none (== 0), and exact weighting (== 1)
    @type weighting: L{int}
    @param binning: binning factor to apply to projections (default: 1 = no binning). binning=2: 2x2 pixels -> 1
    pixel, binning=3: 3x3 pixels -> 1 pixel, etc.)
    @type binning: L{int}
    @param low_pass_freq: lowpass filter (in Nyquist)
    @type low_pass_freq: L{float}
    @param apply_circle_filter: optional circular filter in fourier space
    @type circleFilter: L{bool}
    @param scale_factor_particle: scaling factor for particles
    @type scale_factor_particle: L{float}
    @param angle_specimen: additional rotation of the specimen
    @type angle_specimen: L{float}
    @param verbose: print output default=False
    @type verbose: L{bool}
    @param imdim: x and y dimension of final image
    @type imdim: L{int}
    @param imdimX: input projection x dimension size
    @type imdimX: L{int}
    @param imdimY: input projection y dimension size
    @type imdimY: L{int}

    @return: Will return a list of stacks: projections, phi angles (around y axis), theta angles (around x axis)
    and an offsetStack (user specified x, y offset)
    @rtype: L{list}

    @author: Marten Chaillet
    """
    from pytom.lib.pytom_volume import complexRealMult, vol, pasteCenter, mean
    from pytom.basic.functions import taper_edges
    from pytom.basic.transformations import general_transform2d, resize
    from pytom.basic.fourier import ifft, fft
    from pytom.basic.filter import filter, circleFilter, rampFilter, exactFilter, fourierFilterShift, \
        fourierFilterShift_ReducedComplex
    from pytom.basic.files import read
    import pytom.lib.pytom_freqweight as pytom_freqweight

    projection = projection_list[index + offset]
    tilt_angle = tilt_angles[index + offset]

    # pre-determine analytical weighting function and lowpass for speedup
    if weighting == -1:
        weightSlice = fourierFilterShift(rampFilter(imdim, imdim))
    elif weighting == 1:
        weightSlice = fourierFilterShift(exactFilter(tilt_angles, tilt_angle,
                                                     imdim, imdim, slice_width))

    if apply_circle_filter:
        circleFilterRadius = imdim // 2
        circleSlice = fourierFilterShift_ReducedComplex(circleFilter(imdim, imdim, circleFilterRadius))
    else:
        circleSlice = vol(imdim, imdim // 2 + 1, 1)
        circleSlice.setAll(1.0)

    # design lowpass filter
    if low_pass_freq > 0:
        sigma = 0.015
        lpf = pytom_freqweight.weight(0.0, low_pass_freq * imdim // 2, imdim, imdim // 2 + 1, 1,
                                      sigma * imdim)
        # lpf = pytom_freqweight.weight(0.0, low_pass_freq * imdim // 2, imdim, imdim // 2 + 1, 1,
        #                               low_pass_freq / 5. * imdim)

    if projection._filename.split('.')[-1] == 'st':
        image = read(file=projection.getFilename(),
                     subregion=[0, 0, index - 1 + offset, imdim, imdim, 1],
                     sampling=[0, 0, 0], binning=[0, 0, 0])
    else:
        image = read(projection.getFilename())

    if binning > 1:
        image = resize(volume=image, factor=1 / float(binning))[0]

    # normalize to contrast - subtract mean and norm to mean
    immean = mean(image)
    image = (image - immean) / immean

    # smoothen borders to prevent high contrast oscillations
    image = taper_edges(image, imdim // 30)[0]

    # transform projection according to tilt alignment
    transX = projection.getAlignmentTransX() / binning
    transY = projection.getAlignmentTransY() / binning
    rot = projection.getAlignmentRotation()
    mag = projection.getAlignmentMagnification() * scale_factor_particle
    order = projection.getOperationOrder()

    # 3 -- square if needed
    if imdimY != imdimX:
        newImage = vol(imdim, imdim, 1)
        newImage.setAll(0)
        pasteCenter(image, newImage)
        image = newImage

    # IMOD Rotates around size/2 -0.5, wheras pytom defaults to rotating around size/2
    # so IMOD default order is scaling > rotation > translation, that seems a weird order...
    center = ((image.size_x() - 1) / 2., (image.size_y() - 1) / 2., 0) if order == [1, 0, 2] else None

    if not (rot == 0 and transX == 0 and transY == 0 and mag == 1):
        image = general_transform2d(v=image, rot=rot, shift=[transX, transY], scale=mag, order=order,
                                    crop=True, center=center)

    if low_pass_freq > 0:
        filtered = filter(volume=image, filterObject=lpf, fourierOnly=False)
        image = filtered[0]

    # smoothen once more to avoid edges
    image = taper_edges(image, imdim // 30)[0]

    # apply weighting
    if weighting in [-1, 1]:
        image = ifft(complexRealMult(complexRealMult(fft(image), weightSlice), circleSlice), scaling=True)

    # place back in the stack that is accesible for the main process
    put_in_shared_stack(shared_dtype, shared_stack_shape, image, index)

    if show_progress:
        print(f'parallel-pool: finished aligning projection {index}')

    if verbose:
        print(tilt_angle, projection.getOffsetX(), projection.getOffsetY())

    return index, (projection.getAlignmentAxisAngle(), tilt_angle - angle_specimen,
                   (projection.getOffsetX(), projection.getOffsetY()))


# ============================== Projection Class =========================================

class Projection(PyTomClass):
    """
    Projection: 
    """
    def __init__(self, filename=None, tiltAngle=None, offsetX=0, offsetY=0, alignmentTransX=0.,
                 alignmentTransY=0., alignmentRotation=0., alignmentMagnification=1.,
                 alignmentAxisAngle=0., operationOrder=[2, 0, 1], index=0):
        """
        Initialize projection class

        @param filename: The filename
        @type filename: string
        @param tiltAngle: Tilt angle of a projection, rotation around microscope y-axis
        @type tiltAngle: float
        @param offsetX: offset X is used for subtomogram reconstruction from cut out projections
        @type offsetX: float
        @param offsetY: offset Y is used for subtomogram reconstruction from cut out projections
        @type offsetY: float
        @param alignmentTransX: X Shift determined from alignment method
        @type alignmentTransX: float
        @param alignmentTransY: Y Shift determined from alignment method
        @type alignmentTransY: float
        @param alignmentRotation: in plane rotation from alignment, around microscope z-axis
        @type alignmentRotation: float
        @param alignmentMagnification: Magnification determined from fine alignment method
        @type alignmentMagnification: float
        @param alignmentAxisAngle: rotation angle around microscope system x-axis
        @type alignmentAxisAngle: float
        @param operationOrder: tuple of alignment operation order (r, t, s) => (2, 0, 1) means translation first,
        then scaling, then rotation
        @type operationOrder: list
        @param index: Index of projection in filename (default: 1)
        @type index: int
        """

        self._filename = str(filename) if filename is not None else filename
        self._checkValidFile()

        self._xSize = None
        self._ySize = None

        if filename and tiltAngle is None:
            self._loadTiltAngle()
        else:
            self._tiltAngle = float(tiltAngle)

        self._index = int(index)
        self._offsetX = float(offsetX)
        self._offsetY = float(offsetY)
        self._alignmentTransX = float(alignmentTransX)
        self._alignmentTransY = float(alignmentTransY)
        self._alignmentRotation = float(alignmentRotation)
        self._alignmentMagnification = float(alignmentMagnification)
        self._alignmentAxisAngle = float(alignmentAxisAngle)
        if isinstance(operationOrder, str):
            self._operationOrder = convert_operation_order_str2list(operationOrder)
        else:
            self._operationOrder = operationOrder  # default = (2, 0, 1) translation, scaling, rotation
            # imod default is also (2, 0, 1) but set default to (1, 0, 2) is identical but then we can recognize
            # imod alignments in results files.
            # common in literature is translation, scaling, rotation => (2, 0, 1)

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
        from pytom.gui.mrcOperations import read_mrc_header

        try:
            if self._filename.split('.')[-1] == 'em':
                fh = open(self._filename,'rb')
                fh.seek(168)
                self._tiltAngle = float(struct.unpack('i',fh.read(4))[0])/1000
                fh.close()
            else:
                fh = read_mrc_header(self._filename)
                self._tiltAngle = fh['user'][0]/1000.  # TODO find correct location of the tiltangle and autowrite it.
        except:
            print(f'could not read tilt angle for projection {self._filename}')
            self._tiltAngle = None

    def _checkValidFile(self):
        """
        Class private method that is called whenever a Projection() is initialized or the filename is set,
        it will warn the user if the current filename is valid
        """
        from pytom.tools.files import checkFileExists
        if not checkFileExists(self._filename):
            print('Warning: Projection() initialized without valid filename')
            sys.exit(0)

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
        return proj.size_x(), proj.size_y(), proj.size_z()

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

    def getAlignmentAxisAngle(self):
        """
        tilt around microscope X axis

        @rtype: float
        """
        return self._alignmentAxisAngle

    def getOperationOrder(self, as_string=False):
        """
        get order of applying operations

        @rtype: tuple -> (r, t, s)
        """
        if as_string:
            return convert_operation_order_list2str(self._operationOrder)
        else:
            return self._operationOrder

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
        self._checkValidFile()

    def setTiltAngle(self, tiltAngle):
        """
        set tilt angle in projection

        @param tiltAngle: tilt angle (deg)
        @type tiltAngle: float
        """
        self._tiltAngle = float(tiltAngle)

    def setOffsetX(self, offsetX):
        """
        set offset in X direction

        @param offsetX: offset (pixels)
        @type offsetX: float
        """
        self._offsetX = float(offsetX)

    def setOffsetY(self, offsetY):
        """
        set offset in Y direction

        @param offsetY: offset (pixels)
        @type offsetY float
        """
        self._offsetY = float(offsetY)
    
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

    def setAlignmentAxisAngle(self, alignmentAxisAngle):
        """
        set tilt around microscope X axis

        @type axisAngle: float
        """
        self._alignmentAxisAngle = float(alignmentAxisAngle)

    def setOperationOrder(self, operationOrder):
        """
        set operation order for aligning tilts as a list

        @type operationOrder: list -> [r, t, s]
        """
        if isinstance(operationOrder, str):
            self._operationOrder = convert_operation_order_str2list(operationOrder)
        else:
            self._operationOrder = operationOrder

    def _setSize(self):
        """
        set X- and Y-size of proj according to file
        """
        from pytom.basic.files import readProxy as read
        p = read(self._filename)
        self._xSize = p.size_x()
        self._ySize = p.size_y()
        
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
        self._checkValidFile()
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


# ============================== ProjectionList Class =====================================

class ProjectionList(PyTomClass):
    """
    ProjectionList: List of projections
    """
    def __init__(self, projection_list=None):
        if isinstance(projection_list, ProjectionList):
            self._list = projection_list._list
            self._init_tilt_angles() if projection_list._tilt_angles is None else projection_list._tilt_angles
        elif isinstance(projection_list, list):
            self._list = projection_list
            self._init_tilt_angles()
        else:
            self._list = []
            self._tilt_angles = []

    def _init_tilt_angles(self):
        self._tilt_angles = [p.getTiltAngle() for p in self._list]

    def info(self):
        tline = ""
        for projection in self._list:
            tline = tline + str(projection) + "\n"
        return tline

    def loadDirectory(self, directory, metafile=None, prefix=''):
        """
        loadDirectory: Will take all projections in a directory, determine their respective tilt angle from the
        header and populate this list object.

        @param directory: directory where files are located
        @type directory: L{str}
        @param metafile: metafile with information about tilt series
        @type metafile: L{str}
        @param prefix: prefix of projection files in dir, will be used to filter files
        @type prefix: L{str}
        """
        from pytom.tools.files import checkDirExists
        from pytom.basic.files import loadstar
        from pytom.basic.datatypes import DATATYPE_METAFILE
        import os

        if not checkDirExists(directory):
            raise RuntimeError('Directory ' + directory + ' does not exist!')

        if metafile is None:
            possible_files = [line for line in sorted(os.listdir(directory))
                              if '.meta' == os.path.splitext(line)[1]]
            if len(possible_files) == 1:  # if we find a metafile, assign it
                metafile = os.path.join(directory, possible_files[0])
        if metafile is not None:  # if we have a metafile, load its tilt angles
            metadata = loadstar(metafile, dtype=DATATYPE_METAFILE)
            tiltAngles = metadata['TiltAngle']

        files = [line for line in sorted(os.listdir(directory)) if prefix in line]

        # if metafile is provided we could just append directory with filename in metafile ...
        # metafile could also contain the full path
        for file in files:
            if os.path.splitext(file)[1] in ['.em', '.mrc']:
                i = int(file.split('_')[-1].split('.')[0])  # TODO this is hackish => but leave it for now
                if metafile is not None:
                    projection = Projection(filename=os.path.join(directory, file), tiltAngle=tiltAngles[i], index=i)
                else:
                    projection = Projection(filename=os.path.join(directory, file), index=i)
                self.append(projection)

        if len(self._list) != 0:
            if metafile is not None or not any([p.getTiltAngle() is None for p in self._list]):  # with metafile we
                # have the tilts, otherwise tilts might have been readable from the em/mrc files but we will need to
                # them all
                self.sort()
            else:
                self.sort(key=lambda p: p.getIndex())
        else:
            raise PyTomClassError("No projections in directory that could be initialised for ProjectionList.")

    def load_alignment(self, alignment_file):
        """
        Load the alignment parameters for the projections from an alignment results file.
        If the Projections have not yet been initialized the filename from the alignment results file will be used
        to initialize the projections.

        @param alignment_file: alignment file txt
        @type alignment_file: L{str}
        """
        from pytom.basic.datatypes import DATATYPE_ALIGNMENT_RESULTS_RO, DATATYPE_ALIGNMENT_RESULTS
        from pytom.basic.files import loadstar

        try:
            alignment_results = loadstar(alignment_file, dtype=DATATYPE_ALIGNMENT_RESULTS_RO)
        except:
            alignment_results = loadstar(alignment_file, dtype=DATATYPE_ALIGNMENT_RESULTS)

        # WARNING: the alignment results can be loaded on a specified directory with tilts, if so...
        #           ensure the alignments are appropriate for the projections
        if len(self) == 0:
            for i in range(len(alignment_results)):
                # set offsetx and offsety ??
                projection = Projection(filename=alignment_results[i]['FileName'],
                                        tiltAngle=alignment_results[i]['TiltAngle'],
                                        alignmentTransX=alignment_results[i]['AlignmentTransX'],
                                        alignmentTransY=alignment_results[i]['AlignmentTransY'],
                                        alignmentRotation=alignment_results[i]['InPlaneRotation'],
                                        alignmentMagnification=alignment_results[i]['Magnification'],
                                        index=i)
                if 'OperationOrder' in alignment_results.dtype.names:
                    projection.setOperationOrder(alignment_results[i]['OperationOrder'])
                if 'AxisAngle' in alignment_results.dtype.names:
                    projection.setAlignmentAxisAngle(float(alignment_results[i]['AxisAngle']))
                self.append(projection)
        else:
            assert len(self) == len(alignment_results), print('Current number of projections in ProjectionList does '
                                                              'not correspond with the amount of projections in '
                                                              'the alignment parameters file. Exiting.')
            for i in range(len(self)):
                self[i].setTiltAngle(alignment_results[i]['TiltAngle'])
                self[i].setAlignmentTransX(alignment_results[i]['AlignmentTransX'])
                self[i].setAlignmentTransY(alignment_results[i]['AlignmentTransY'])
                self[i].setAlignmentRotation(alignment_results[i]['InPlaneRotation'])
                self[i].setAlignmentMagnification(alignment_results[i]['Magnification'])
                if 'OperationOrder' in alignment_results.dtype.names:
                    self[i].setOperationOrder(alignment_results[i]['OperationOrder'])
                if 'AxisAngle' in alignment_results.dtype.names:
                    self[i].setAlignmentAxisAngle(alignment_results[i]['AxisAngle'])

        # sort by tilt angle
        self.sort()

    def ready_for_reconstruction(self):
        """
        Check if the projection list has sufficient information to start reconstruction.
        @return: True or False
        @rtype: Bool
        """
        from pytom.tools.files import checkFileExists

        # check if files all exist
        files_ready = all([checkFileExists(p._filename) for p in self._list]) if len(self._list) > 0 else False

        # if files exist but not tilt angles are known try to initialize them
        if files_ready and len(self._tilt_angles) == 0:
            self._init_tilt_angles()

        # if we still dont have them this will be False and we wont be ready
        angles_ready = all([a is not None for a in self._tilt_angles]) if len(self._tilt_angles) > 0 else False

        return files_ready and angles_ready

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

        # sort the projections by tilt angle
        if not any([p.getTiltAngle() is None for p in self._list]):
            self.sort()
        else:
            self.sort(key=lambda p: p.getIndex())

    def sort(self, key=lambda p: p.getTiltAngle()):
        """
        Sort the ProjectionList based on key.
        @param key: lambda expression to select which item of the Projection to sort on
        @type key: L{function}
        """
        sorted_list = sorted(self._list, key=key)
        if sorted_list != self._list:
            self._list = sorted_list
            [p.setIndex(i) for i, p in enumerate(self._list) if i != p.getIndex()]  # update indices of the
            # if they have been sorted with key .getIndex() this will not change the indices of the projections
        if not any([p.getTiltAngle() is None for p in self._list]):
            self._init_tilt_angles()
        else:
            print('Warning: ProjectionList could not initialize any or all of the tilt angles')
        
    def append(self, projection):
        """
        append projection to projection list

        @param projection: a projection
        @type projection: L{pytom.reconstructionStructures.Projection}
        """
        self._list.append(projection)
        self._tilt_angles.append(projection.getTiltAngle())

    def __getitem__(self, key):
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

    # Backward compatability
    def generateVolumes(self, particles, cubeSize, binning=1, 
            applyWeighting=False, showProgressBar=False, verbose=False,
            preScale=1, postScale=1):
        """
        generateVolumes (old function name, now replaced by reconstructVolumes)
        @deprecated: 
        """
        print('Generate Volumes is deprecated, use reconstructVolumes instead!')
        w = -1 if applyWeighting else 0
        self.reconstructVolumes(particles, cubeSize, binning=binning, weighting=w,
                                show_progress_bar=showProgressBar,
                                verbose=verbose, pre_scale_factor=preScale, post_scale=postScale)

    # ================== Tomogram reconstruction ============================
    def reconstructVolume(self, dims=[512, 512, 128], reconstructionPosition=[0, 0, 0],
                          recon_interpolation='filt_bspline', binning=1, weighting=0, specimen_angle=0,
                          low_pass_ny_fraction=0.9, apply_circle_filter=True, scale_factor=1., tilt_range=None,
                          show_progress_bar=True, verbose=False, num_procs=1):
        """ reconstruct a a single 3D volume from weighted and aligned projections on either a gpu or a cpu

        @param dims: 3D Size of reconstructed volume
        @type dims: list
        @param reconstructionPosition: offset center of reconstruction (BEFORE possible binning)
        @type reconstructionPosition: list
        @param read_only: only reading the file, not processing
        @type reconstructionPosition: bool
        @param binning: down-sizing of projections prior to 3d rec
        @type binning: int
        @param weighting: apply weighting
        @type weighting: int
        @param alignmentResultsFile: filename of alignment results file (see pytom.basic.datatypes).
        @type alignmentResultsFile: txt
        @param gpu: run reconstruction on gpu
        @type gpu: int
        @param specimen_angle: angle of the specimen, the tilt angle will be corrected by this angle
        @type specimen_angle: float

        @return: reconstructed volume
        @rtype: L{pytom.lib.pytom_volume.vol} or L{numpy.ndarray}
        """
        from pytom.gpu.initialize import device

        assert len(dims) == 3, "invalid dimensions for tomogram, needs to be 3 dimensional"

        if 'gpu' in device:
            vol = self.reconstructVolumeGPU(dims=dims, reconstructionPosition=reconstructionPosition,
                                            recon_interpolation=recon_interpolation, binning=binning,
                                            weighting=weighting, specimen_angle=specimen_angle,
                                            low_pass_ny_fraction=low_pass_ny_fraction,
                                            apply_circle_filter=apply_circle_filter, scale_factor=scale_factor,
                                            tilt_range=tilt_range, show_progress_bar=show_progress_bar, verbose=verbose)
        else:
            vol = self.reconstructVolumeCPU(dims=dims, reconstructionPosition=reconstructionPosition,
                                            recon_interpolation=recon_interpolation, binning=binning,
                                            weighting=weighting, specimen_angle=specimen_angle,
                                            low_pass_ny_fraction=low_pass_ny_fraction,
                                            apply_circle_filter=apply_circle_filter, scale_factor=scale_factor,
                                            tilt_range=tilt_range, show_progress_bar=show_progress_bar,
                                            verbose=verbose, num_procs=num_procs)
        return vol

    def reconstructVolumeGPU(self, dims=[512, 512, 128], reconstructionPosition=[0, 0, 0],
                             recon_interpolation='filt_bspline', binning=1, weighting=0,
                             specimen_angle=0, low_pass_ny_fraction=0.9, apply_circle_filter=True, scale_factor=1.,
                             tilt_range=None, show_progress_bar=True, verbose=False):
        """
        reconstruct a single 3D volume from weighted and aligned projections on a gpu

        @param dims: 3D Size of reconstructed volume
        @type dims: list
        @param reconstructionPosition: offset center of reconstruction (BEFORE possible binning)
        @type reconstructionPosition: list
        @param interpolation: interpolation method used (default='filt_bspline')
        @type reconstructionPosition: str
        @param binning: down-sizing of projections prior to 3d rec
        @type binning: int
        @param weighting: apply weighting
        @type weighting: int
        @param alignmentResultsFile: filename of alignment results file (see pytom.basic.datatypes).
        @type alignmentResultsFile: txt
        @param gpu: run reconstruction on gpu
        @type gpu: int
        @param specimen_angle: angle of the specimen, the tilt angle will be corrected by this angle
        @type specimen_angle: float
        @param read_only: only reading the file, not processing
        @type reconstructionPosition: bool

        @return: reconstructed volume
        @rtype: L{numpy.ndarray}

        @author: FF
        last change: include binning in reconstructionPosition
        """
        from pytom.gpu.initialize import device, xp
        from pytom.agnostic.reconstruction_functions import backProjectGPU
        import time

        ready = self.ready_for_reconstruction()
        if not ready:
            print('ProjectionList is not initialized with sufficient parameters to start reconstruction.')
            sys.exit(0)

        s = time.time()

        # TODO what to do with scale factor?? should also be applied to positions for reconstructions
        projections, vol_phi, proj_angles, vol_offsetProjections = self.to_projection_stack_gpu(
            weighting=weighting, binning=binning, low_pass_freq=low_pass_ny_fraction,
            apply_circle_filter=apply_circle_filter, scale_factor_particle=scale_factor, angle_specimen=specimen_angle,
            tilt_range=tilt_range, show_progress_bar=show_progress_bar, verbose=verbose)

        # volume storing reconstruction offset from center (x,y,z)
        recPosVol = xp.zeros((projections.shape[2], 3), dtype=xp.float32)

        for ii in range(0,3):
            recPosVol[:, ii] = float(reconstructionPosition[ii] / binning)

        print(f'alignment time: {time.time()-s}')
        s = time.time()

        # finally backproject into volume
        # TODO interpolation methods is also not used in GPU back project
        vol_bp = xp.zeros(dims, dtype=xp.float32)
        backProjectGPU(projections, vol_bp, vol_phi, proj_angles, recPosVol, vol_offsetProjections,
                             recon_interpolation)

        print(f'backproject time: {time.time()-s}')

        return vol_bp

    def reconstructVolumeCPU(self, dims=[512, 512, 128], reconstructionPosition=[0, 0, 0],
                             recon_interpolation='filt_bspline', binning=1, weighting=0,
                             specimen_angle=0, low_pass_ny_fraction=0.9, apply_circle_filter=True, scale_factor=1.,
                             tilt_range=None, show_progress_bar=True, verbose=False, num_procs=1):
        """
        reconstruct a single 3D volume from weighted and aligned projections on a CPU

        @param dims: 3D Size of reconstructed volume
        @type dims: list
        @param reconstructionPosition: offset center of reconstruction (BEFORE possible binning)
        @type reconstructionPosition: list
        @param read_only: only reading the file, not processing
        @type reconstructionPosition: bool
        @param binning: down-sizing of projections prior to 3d rec
        @type binning: int
        @param weighting: apply weighting
        @type weighting: int
        @param alignmentResultsFile: filename of alignment results file (see pytom.basic.datatypes).
        @type alignmentResultsFile: txt
        @param gpu: run reconstruction on gpu
        @type gpu: int
        @param specimen_angle: angle of the specimen, the tilt angle will be corrected by this angle
        @type specimen_angle: float

        @return: reconstructed volume
        @rtype: L{pytom.lib.pytom_volume.vol}

        @author: FF
        last change: include binning in reconstructionPosition
        """
        from pytom.lib.pytom_volume import vol, backProject
        import time
        import sys

        ready = self.ready_for_reconstruction()
        if not ready:
            print('ProjectionList is not initialized with sufficient parameters to start reconstruction.')
            sys.exit(0)

        if recon_interpolation != 'filt_bspline':
            print('WARNING: other interpolation methods for reconstruction on CPU currently not available')

        s = time.time()

        vol_bp = vol(dims[0], dims[1], dims[2])
        vol_bp.setAll(0.0)

        if num_procs == 1:  # start alignment function based on available cpu cores
            vol_img, vol_phi, vol_the, vol_offsetProjections = self.to_projection_stack(
                weighting=weighting, binning=binning, low_pass_freq=low_pass_ny_fraction,
                apply_circle_filter=apply_circle_filter, scale_factor_particle=scale_factor,
                angle_specimen=specimen_angle, tilt_range=tilt_range, show_progress_bar=show_progress_bar, \
                                                                                      verbose=verbose)

        elif isinstance(num_procs, int) and num_procs > 1:
            vol_img, vol_phi, vol_the, vol_offsetProjections = self.to_projection_stack_parallel(
                weighting=weighting, binning=binning, low_pass_freq=low_pass_ny_fraction,
                apply_circle_filter=apply_circle_filter, scale_factor_particle=scale_factor, tilt_range=tilt_range,
                angle_specimen=specimen_angle, num_procs=num_procs, show_progress_bar=show_progress_bar, verbose=verbose)

        else:
            print('Invalid number of processes for aligning projections. Exiting...')
            sys.exit(0)

        # volume storing reconstruction offset from center (x,y,z)
        recPosVol = vol(3, vol_img.size_z(), 1)
        recPosVol.setAll(0.0)
        for iproj in range(0, vol_img.size_z()):
            for ii in range(0, 3):
                recPosVol.setV(float(reconstructionPosition[ii] / binning), ii, iproj, 0)

        print(f'alignment time: {time.time() - s}')
        s = time.time()

        # backproject into volume
        backProject(vol_img, vol_bp, vol_phi, vol_the, recPosVol, vol_offsetProjections)

        print(f'backproject time: {time.time() - s}')

        return vol_bp

    # ====================== Subtomogram reconstruction =============
    def reconstructVolumes(self, particles, cube_size, binning=1, weighting=0, low_pass_ny_fraction=0.9,
                           apply_circle_filter=True,  pre_scale_factor=1., post_scale=1., num_procs=1,
                           ctfcenter=None, polishResultFile='', specimen_angle=0.,
                           tilt_range=None, show_progress_bar=False,
                           verbose=False, gpuIDs=None):
        """
        reconstructVolumes: reconstruct a subtomogram given a particle object.

        @param particles: A list of particles
        @type particles: L{pytom.basic.structures.ParticleList}
        @param cube_size: 3D Size of reconstructed particle (here, x == y == z)
        @param binning: binning of projections
        @param applyWeighting: Will apply weighting for weighted backprojection to projection. Default is off (False)!
        @param showProgressBar: False by default
        @param verbose: False by default
        @param preScale: Scale (bin) the projections BEFORE reconstruction
        @param post_scale: Scale the tomogram AFTER reconstruction
        @param num_procs: Number of parallel processes
        @param alignmentResultsFile: filename of alignment results file (see pytom.basic.datatypes).
        @type alignmentResultsFile: txt
        @param polishResultFile: result file generate by the particle polishing (see pytom.basic.datatypes).
        @type polishResultFile: txt
        @param scaleFactorParticle: scaling factor for the reconstruction default=1
        @type scaleFactorParticle: float
        @param gpuIDs: List of gpu-ids on which reconstruction will take place.
        @type gpuIDs: list
        @param specimen_angle: angle of the specimen, the tilt angle will be corrected by this angle
        """
        import time, os, sys
        from pytom.agnostic.io import write
        # from pytom.basic.files import write_em
        from pytom.lib.pytom_volume import vol, backProject, rescaleSpline
        from pytom.tools.ProgressBar import FixedProgBar
        from multiprocessing import Process, set_start_method

        ready = self.ready_for_reconstruction()
        if not ready:
            print('ProjectionList is not initialized with sufficient parameters to start reconstruction.')
            sys.exit(0)

        # important to set to allow child processes on multiple GPU's (otherwise cupy init is not working)
        set_start_method('spawn')

        if len(particles) == 0:
            raise RuntimeError('ParticleList is empty!')

        if cube_size.__class__ != int:
            raise TypeError('reconstructVolumes: Parameter cube_size must be of type int!')

        if gpuIDs is not None:
            if isinstance(gpuIDs, list):
                num_procs = len(gpuIDs)
            elif isinstance(gpuIDs, int):
                num_procs = 1
            else:
                print('invalid argument for gpuIDs, exiting...')
                sys.exit(0)

        # stacks for images, projections angles etc.
        if gpuIDs is None:
            # Running on the CPU
            if num_procs > 1:
                stacks = self.to_projection_stack_parallel(weighting=weighting, binning=binning,
                                                            low_pass_freq=low_pass_ny_fraction,
                                                            apply_circle_filter=apply_circle_filter,
                                                            scale_factor_particle=pre_scale_factor,
                                                            angle_specimen=specimen_angle, tilt_range=tilt_range,
                                                            show_progress_bar=show_progress_bar, verbose=verbose,
                                                            num_procs=num_procs)
            else:
                stacks = self.to_projection_stack(weighting=weighting, binning=binning,
                                                   low_pass_freq=low_pass_ny_fraction,
                                                   apply_circle_filter=apply_circle_filter,
                                                   scale_factor_particle=pre_scale_factor,
                                                   angle_specimen=specimen_angle, tilt_range=tilt_range,
                                                   show_progress_bar=show_progress_bar, verbose=verbose)
        else:
            # Running on the GPU
            stacks = self.to_projection_stack_gpu(weighting=weighting, binning=binning,
                                                   low_pass_freq=low_pass_ny_fraction,
                                                   apply_circle_filter=apply_circle_filter,
                                                   scale_factor_particle=pre_scale_factor,
                                                   angle_specimen=specimen_angle, tilt_range=tilt_range,
                                                   show_progress_bar=show_progress_bar, verbose=verbose)
            # stacks = [x.get() for x in stacks]

        # create list for processes
        procs = []

        # set the temp dir from particle storage path
        temp_dir = os.path.dirname(particles[0].getFilename())

        if num_procs > 1 or gpuIDs is not None:

            # temp storage of projections
            for n in range(4):
                if not os.path.exists(temp_dir):
                    os.mkdir(temp_dir)
                write(os.path.join(temp_dir, '.temp_' + str(n) + '.em'), stacks[n])
            del stacks

            # set gpu or cpu extraction
            # prepare gpuIDs as a list to assign to processes
            if gpuIDs is None:
                extract = self.extract_single_particle
                # set gpus to none for each process
                gpuIDs = [None, ] * num_procs
            else:
                extract = self.extract_particles_on_gpu
                if isinstance(gpuIDs, int):
                    gpuIDs = [gpuIDs, ]

            # create each process
            for i in range(num_procs):
                # select the particles for the process
                ps = particles[i::num_procs]

                # if empty skip
                if not len(ps) > 0:
                    continue

                # start all processes
                proc = Process(target=extract, args=(ps, i, verbose, binning, post_scale, cube_size, polishResultFile,
                                                     gpuIDs[i], temp_dir))
                procs.append(proc)
                proc.start()

                # print('Subtomogram reconstruction for particle {} has started (batch mode).'.format(particleIndex))
        else:
            if show_progress_bar:
                progressBar = FixedProgBar(0, len(particles), 'Particle volumes generated ')
                progressBar.update(0)

            t = time.time()

            vol_bp = vol(cube_size, cube_size, cube_size)
            [vol_img, vol_phi, vol_the, vol_offsetProjections] = stacks

            for particleIndex in range(len(particles)):
                t = time.time()
                p = particles[particleIndex]
                # self.extract_single_particle(p, particleIndex, verbose, binning, post_scale, cube_size)
                # print('Subtomogram reconstruction for particle {} has started.'.format(particleIndex))
                if show_progress_bar:
                    progressBar.update(particleIndex)

                if verbose:
                    print(p)

                vol_bp.setAll(0.0)

                reconstructionPosition = vol(3, vol_img.size_z(), 1)
                reconstructionPosition.setAll(0.0)

                # adjust coordinates of subvolumes to binned reconstruction
                for i in range(vol_img.size_z()):
                    reconstructionPosition(float(p.getPickPosition().getX() / binning), 0, i, 0)
                    reconstructionPosition(float(p.getPickPosition().getY() / binning), 1, i, 0)
                    reconstructionPosition(float(p.getPickPosition().getZ() / binning), 2, i, 0)

                if verbose:
                    print((p.getPickPosition().getX() / binning, p.getPickPosition().getY() / binning,
                           p.getPickPosition().getZ() / binning))

                backProject(vol_img, vol_bp, vol_phi, vol_the, reconstructionPosition, vol_offsetProjections)

                if post_scale > 1:
                    volumeRescaled = vol(cube_size / post_scale, cube_size / post_scale, cube_size / post_scale)
                    rescaleSpline(vol_bp, volumeRescaled)
                    volumeRescaled.write(p.getFilename())
                else:
                    vol_bp.write(p.getFilename())

            print(f'recon time: {time.time()-t:.3f} sec')

        while len(procs):
            time.sleep(0.1)
            procs = [proc for proc in procs if proc.is_alive()]

        if num_procs > 1 or gpuIDs is not None:
            for n in range(4):
                os.system('rm -f {}'.format(os.path.join(temp_dir, '.temp_' + str(n) + '.em')))

        print('\n Subtomogram reconstructions have finished.\n\n')

    def extract_particles_on_gpu(self, particles, pid, verbose, binning, post_scale, cube_size, filename_ppr,
                                 gpuID, temp_dir):
        import os
        from pytom.agnostic.io import read, write
        from pytom.gpu.initialize import xp, device
        from pytom.agnostic.transform import resize
        from pytom.agnostic.reconstruction_functions import backProjectGPU
        import time

        xp.cuda.Device(gpuID).use()

        print(f'start recon in process: {pid}')
        t = time.time()

        try:
            # read the temp stacks
            stacks = []
            for n in range(4):
                stacks.append(read(os.path.join(temp_dir, '.temp_' + str(n) + '.em'), keepnumpy=True))

            # put in gpu mem
            [projections, vol_phi, vol_the, vol_offsetProjections] = [xp.asarray(x) for x in stacks]
            num_projections = projections.shape[2]

            vol_bp = xp.zeros((cube_size, cube_size, cube_size), dtype=xp.float32)
            reconstructionPosition = xp.zeros((num_projections, 3), dtype=xp.float32)

            interpolation = 'filt_bspline'

            proj_angles = vol_the.squeeze()

            for p in particles:

                # set back project volume and reconstruction position back to zero
                reconstructionPosition.fill(.0)
                vol_bp.fill(.0)

                # adjust coordinates of subvolumes to binned reconstruction
                reconstructionPosition[:, 0] = float(p.getPickPosition().getX() / binning)
                reconstructionPosition[:, 1] = float(p.getPickPosition().getY() / binning)
                reconstructionPosition[:, 2] = float(p.getPickPosition().getZ() / binning)

                # run the back projection for the particle
                backProjectGPU(projections, vol_bp, vol_phi.squeeze(), proj_angles,
                               reconstructionPosition, vol_offsetProjections, interpolation)

                # do post scaling
                if post_scale > 1:
                    # scale and write
                    scaled = resize(vol_bp, 1. / post_scale, interpolation='Spline')
                    write(p.getFileName(), scaled)
                else:
                    write(p.getFilename(), vol_bp)

        except Exception as e:
            print('Caught exception in worker thread (x = %d):' % pid)
            print()
            raise e

        print(f'recon time in process {pid}: {time.time()-t:.3f} sec')

    def extract_single_particle(self, particles, pid, verbose, binning, post_scale, cube_size, filename_ppr,
                                gpuID, temp_dir):
        from pytom.lib.pytom_volume import vol, backProject, rescaleSpline
        from pytom.basic.files import read
        import numpy as np
        import time
        import os

        print(f'start recon in process: {pid}')
        t = time.time()

        try:
            # read the temp stacks
            stacks = []
            for n in range(4):
                stacks.append(read(os.path.join(temp_dir, '.temp_' + str(n) + '.em')))

            # put in variables
            [vol_img, vol_phi, vol_the, vol_offsetProjections] = stacks
            num_projections = vol_img.size_z()

            vol_bp = vol(cube_size, cube_size, cube_size)
            vol_bp.setAll(0.0)

            reconstructionPosition = vol(3, num_projections, 1)
            reconstructionPosition.setAll(0.0)
            start_index = pid * len(self)

            for p in particles:
                # reset
                reconstructionPosition.setAll(0.0)
                vol_bp.setAll(0.0)

                # adjust coordinates of subvolumes to binned reconstruction
                if not filename_ppr:
                    for i in range(num_projections):
                        reconstructionPosition(float(p.getPickPosition().getX() / binning), 0, i, 0)
                        reconstructionPosition(float(p.getPickPosition().getY() / binning), 1, i, 0)
                        reconstructionPosition(float(p.getPickPosition().getZ() / binning), 2, i, 0)

                else:
                    # possible particle polish results loading
                    x, y, z = float(p.getPickPosition().getX() / binning), float(
                        p.getPickPosition().getY() / binning), float(p.getPickPosition().getZ() / binning)

                    for i in range(len(self)):
                        from pytom.gui.guiFunctions import LOCAL_ALIGNMENT_RESULTS, loadstar
                        import pytom.basic.combine_transformations as ct

                        particle_polish_file = loadstar(filename_ppr, dtype=LOCAL_ALIGNMENT_RESULTS)

                        offsetX = float(particle_polish_file['AlignmentTransX'][start_index + i])
                        offsetY = float(particle_polish_file['AlignmentTransY'][start_index + i])

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

                backProject(vol_img, vol_bp, vol_phi, vol_the, reconstructionPosition, vol_offsetProjections)

                if post_scale > 1:
                    volumeRescaled = vol(cube_size / post_scale, cube_size / post_scale, cube_size / post_scale)
                    rescaleSpline(vol_bp, volumeRescaled)
                    volumeRescaled.write(p.getFilename())
                else:
                    vol_bp.write(p.getFilename())

        except Exception as e:
            print('Caught exception in worker thread (x = %d):' % pid)
            print()
            raise e

        print(f'recon time in process {pid}: {time.time()-t:.3f} sec')

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
        from pytom.lib.pytom_volume import vol, complexRealMult, rescaleSpline
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
        imgSizeX = tmpImage.size_x()
        imgSizeY = tmpImage.size_y()

        if applyWeighting:
            weightSlice = fourierFilterShift(rampFilter(projectionSize, projectionSize))
            circleFilterRadius = tmpImage.size_x() / 2
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

    # ====== new general to_projection_stack functions ===================
    def to_projection_stack(self, weighting=0, binning=1, low_pass_freq=0.9, apply_circle_filter=True,
                            scale_factor_particle=1., angle_specimen=0., particle_diameter=None, pixel_size=2.62,
                            tilt_range=None, show_progress_bar=False, verbose=False):
        """
        Create a projection stack from current Projections in the ProjectionList by weighting the projections.
        Alignment is applied if aligned parameters have been initialized on the Projections.

        @param weighting: ramp weighting (== -1), none (== 0), and exact weighting (== 1)
        @type weighting: L{int}
        @param binning: binning factor to apply to projections (default: 1 = no binning). binning=2: 2x2 pixels -> 1
        pixel, binning=3: 3x3 pixels -> 1 pixel, etc.)
        @type binning: L{int}
        @param low_pass_freq: lowpass filter (in Nyquist)
        @type low_pass_freq: L{float}
        @param apply_circle_filter: optional circular filter in fourier space
        @type circleFilter: L{bool}
        @param scale_factor_particle: scaling factor for particles
        @type scale_factor_particle: L{float}
        @param angle_specimen: additional rotation of the specimen
        @type angle_specimen: L{float}
        @param show_progress_bar: write progress to terminal
        @type show_progress_bar: L{bool}
        @param verbose: print output default=False
        @type verbose: L{bool}

        @return: Will return a list of stacks: projections, phi angles (around y axis), theta angles (around x axis)
        and an offsetStack (user specified x, y offset)
        @rtype: L{list}

        @author: Marten Chaillet
        """
        from pytom.lib.pytom_volume import complexRealMult, vol, paste, pasteCenter, mean
        from pytom.basic.functions import taper_edges
        from pytom.basic.transformations import general_transform2d, resize
        from pytom.basic.fourier import ifft, fft
        from pytom.basic.filter import filter, circleFilter, rampFilter, exactFilter, fourierFilterShift, \
            fourierFilterShift_ReducedComplex
        from pytom.basic.files import read
        from pytom.tools.ProgressBar import FixedProgBar
        import pytom.lib.pytom_freqweight as pytom_freqweight

        if low_pass_freq > 1.:
            low_pass_freq = 1.
            print("Warning: lowpassFilter > 1 - set to 1 (=Nyquist)")

        if tilt_range is None:
            tilt_range = (int(round(min(self._tilt_angles))), int(round(max(self._tilt_angles))))

        imdimX = self[0].getXSize()
        imdimY = self[0].getYSize()

        if binning > 1:
            imdimX = int(float(imdimX) / float(binning) + .5)
            imdimY = int(float(imdimY) / float(binning) + .5)

        imdim = max(imdimX, imdimY)

        # get a list of tilt angles that are in the range
        tilts_in_range = [t for t in self._tilt_angles if tilt_range[0] <= int(round(t)) <= tilt_range[1]]

        # find the offset of the lower limit of provided tilt-range
        offset = 0
        while int(round(self._tilt_angles[offset])) < tilt_range[0]:
            offset += 1

        # prep stacks for output
        stack = vol(imdim, imdim, len(tilts_in_range))
        stack.setAll(0.0)

        phiStack = vol(1, 1, len(tilts_in_range))
        phiStack.setAll(0.0)

        thetaStack = vol(1, 1, len(tilts_in_range))
        thetaStack.setAll(0.0)

        offsetStack = vol(1, 2, len(tilts_in_range))
        offsetStack.setAll(0.0)

        if particle_diameter is None:  # if no particle diameter provided we set slice width relative to full image size
            slice_width = 1
        elif particle_diameter <= 0: # Particle diameter should be bigger than 0
            raise ValueError("particle diameter must be greater than 0")
        else:  # ideal slice width for exact filter is relative to particle diameter
            slice_width = (pixel_size * binning * imdim) / particle_diameter

        # pre-determine analytical weighting function and lowpass for speedup
        if weighting == -1:
            weightSlice = fourierFilterShift(rampFilter(imdim, imdim))

        if apply_circle_filter:
            circleFilterRadius = imdim // 2
            circleSlice = fourierFilterShift_ReducedComplex(circleFilter(imdim, imdim, circleFilterRadius))
        else:
            circleSlice = vol(imdim, imdim // 2 + 1, 1)
            circleSlice.setAll(1.0)

        # design lowpass filter
        if low_pass_freq > 0:
            sigma = 0.015
            lpf = pytom_freqweight.weight(0.0, low_pass_freq * imdim // 2, imdim, imdim // 2 + 1, 1,
                                          sigma * imdim)
            # lpf = pytom_freqweight.weight(0.0, low_pass_freq * imdim // 2, imdim, imdim // 2 + 1, 1,
            #                               low_pass_freq / 5. * imdim)

        if show_progress_bar:
            progressBar = FixedProgBar(0, len(tilts_in_range), 'Creating aligned and weighted projections')
            progressBar.update(0)

        if verbose:
            print('projections to go over: ', [i for i, _ in enumerate(tilts_in_range)])

        for i, (tilt_angle, projection) in enumerate(zip(self._tilt_angles, self._list)):

            if not (tilt_range[0] <= int(round(tilt_angle)) <= tilt_range[1]):
                continue

            if projection._filename.split('.')[-1] == 'st':
                image = read(file=projection.getFilename(),
                             subregion=[0, 0, i - 1, imdim, imdim, 1],
                             sampling=[0, 0, 0], binning=[0, 0, 0])
            else:
                image = read(projection.getFilename())

            if binning > 1:
                image = resize(volume=image, factor=1 / float(binning))[0]

            # normalize to contrast - subtract mean and norm to mean
            immean = mean(image)
            image = (image - immean) / immean

            # smoothen borders to prevent high contrast oscillations
            image = taper_edges(image, imdim // 30)[0]

            # transform projection according to tilt alignment
            transX = projection.getAlignmentTransX() / binning
            transY = projection.getAlignmentTransY() / binning
            rot = projection.getAlignmentRotation()
            mag = projection.getAlignmentMagnification() * scale_factor_particle
            order = projection.getOperationOrder()

            # 3 -- square if needed
            if imdimY != imdimX:
                newImage = vol(imdim, imdim, 1)
                newImage.setAll(0)
                pasteCenter(image, newImage)
                image = newImage

            # IMOD Rotates around size/2 -0.5, wheras pytom defaults to rotating around size/2
            # so IMOD default order is scaling > rotation > translation, default order usually
            center = ((image.size_x() - 1) / 2., (image.size_y() - 1) / 2., 0) if order == [1, 0, 2] else None

            if not (rot == 0 and transX == 0 and transY == 0 and mag == 1):
                image = general_transform2d(v=image, rot=rot, shift=[transX, transY], scale=mag, order=order,
                                            crop=True, center=center)

            if low_pass_freq > 0:
                filtered = filter(volume=image, filterObject=lpf, fourierOnly=False)
                image = filtered[0]

            # smoothen once more to avoid edges
            image = taper_edges(image, imdim // 30)[0]

            # weighting
            if weighting == -1:  # ramp filter
                image = ifft(complexRealMult(complexRealMult(fft(image), weightSlice), circleSlice), scaling=True)
            elif weighting == 1:  # exact filter
                weightSlice = fourierFilterShift(exactFilter(self._tilt_angles, tilt_angle, imdim, imdim, slice_width))
                image = ifft(complexRealMult(complexRealMult(fft(image), weightSlice), circleSlice), scaling=True)

            phiStack(projection.getAlignmentAxisAngle(), 0, 0, i - offset)
            thetaStack(tilt_angle - angle_specimen, 0, 0, i - offset)  # TODO why round + int?
            offsetStack(int(round(projection.getOffsetX())), 0, 0, i - offset)
            offsetStack(int(round(projection.getOffsetY())), 0, 1, i - offset)
            paste(image, stack, 0, 0, i - offset)

            if show_progress_bar:
                progressBar.update(i - offset)

            if verbose:
                print(tilt_angle, projection.getOffsetX(), projection.getOffsetY())

        return [stack, phiStack, thetaStack, offsetStack]

    def _create_shared_stack(self, imdim, n_tilts):
        import numpy as np
        import ctypes
        import multiprocessing as mp
        # create shared arrays
        self._shared_dtype = np.float32
        self._shared_stack_shape = (imdim, imdim, n_tilts)
        self._shared_stack = mp.Array(ctypes.c_float, imdim * imdim * n_tilts)
        self._shared_stack_np = np.frombuffer(self._shared_stack.get_obj(),
                                                dtype=self._shared_dtype).reshape(self._shared_stack_shape,
                                                                                  order='F')
        self._shared_stack_np[:] = np.zeros(self._shared_stack_shape, dtype=self._shared_dtype, order='F')

    def _retrieve_shared_stack(self, stack_vol):
        from pytom.lib.pytom_numpy import npy2vol
        copy = self._shared_stack_np.copy(order='F')
        del self._shared_stack_np, self._shared_stack
        stack_vol.copyVolume(npy2vol(copy, len(copy.shape)))

    def to_projection_stack_parallel(self, weighting=0, binning=1, low_pass_freq=0.9, apply_circle_filter=True,
                                     scale_factor_particle=1., angle_specimen=0., particle_diameter=None,
                                     pixel_size=2.62, tilt_range=None, num_procs=1, show_progress_bar=False,
                                     verbose=False):
        """
        Parallel version.
        Create a projection stack from current Projections in the ProjectionList by weighting the projections.
        Alignment is applied if aligned parameters have been initialized on the Projections.

        @param weighting: ramp weighting (== -1), none (== 0), and exact weighting (== 1)
        @type weighting: L{int}
        @param binning: binning factor to apply to projections (default: 1 = no binning). binning=2: 2x2 pixels -> 1
        pixel, binning=3: 3x3 pixels -> 1 pixel, etc.)
        @type binning: L{int}
        @param low_pass_freq: lowpass filter (in Nyquist)
        @type low_pass_freq: L{float}
        @param apply_circle_filter: optional circular filter in fourier space
        @type circleFilter: L{bool}
        @param scale_factor_particle: scaling factor for particles
        @type scale_factor_particle: L{float}
        @param angle_specimen: additional rotation of the specimen
        @type angle_specimen: L{float}
        @param show_progress_bar: write progress to terminal
        @type show_progress_bar: L{bool}
        @param verbose: print output default=False
        @type verbose: L{bool}
        @param num_procs: number of available processes for parallelization
        @type num_procs: L{int}

        @return: Will return a list of stacks: projections, phi angles (around y axis), theta angles (around x axis)
        and an offsetStack (user specified x, y offset)
        @rtype: L{list}

        @author: Marten Chaillet
        """
        from pytom.lib.pytom_volume import vol
        from functools import partial
        import multiprocessing as mp
        from contextlib import closing

        if low_pass_freq > 1.:
            low_pass_freq = 1.
            print("Warning: lowpassFilter > 1 - set to 1 (=Nyquist)")

        if tilt_range is None:
            tilt_range = (int(round(min(self._tilt_angles))), int(round(max(self._tilt_angles))))

        imdimX = self[0].getXSize()
        imdimY = self[0].getYSize()

        if binning > 1:
            imdimX = int(float(imdimX) / float(binning) + .5)
            imdimY = int(float(imdimY) / float(binning) + .5)

        imdim = max(imdimX, imdimY)

        if particle_diameter is None:  # if no particle diameter provided we set slice width relative to full image size
            slice_width = 1
        elif particle_diameter <= 0: # Particle diameter should be bigger than 0
            raise ValueError("particle diameter must be greater than 0")
        else:  # ideal slice width for exact filter is relative to particle diameter
            slice_width = (pixel_size * binning * imdim) / particle_diameter

        # get a list of tilt angles that are in the range
        tilts_in_range = [t for t in self._tilt_angles if tilt_range[0] <= int(round(t)) <= tilt_range[1]]

        # find the offset of the lower limit of provided tilt-range
        offset = 0
        while int(round(self._tilt_angles[offset])) < tilt_range[0]:
            offset += 1

        # prepare the image stack shared between the processes
        self._create_shared_stack(imdim, len(tilts_in_range))

        if verbose:
            print('projections to go over: ', [i for i, _ in enumerate(tilts_in_range)])

        with closing(mp.Pool(num_procs, initializer=init_pool_processes, initargs=(self._shared_stack, ))) as p:
            results = p.map(partial(align_single_projection, projection_list=self._list,
                                    tilt_angles=self._tilt_angles, shared_dtype=self._shared_dtype,
                                    shared_stack_shape=self._shared_stack_shape, weighting=weighting, binning=binning,
                                    low_pass_freq=low_pass_freq, apply_circle_filter=apply_circle_filter,
                                    scale_factor_particle=scale_factor_particle, angle_specimen=angle_specimen,
                                    verbose=verbose, show_progress=show_progress_bar, imdim=imdim, imdimX=imdimX,
                                    imdimY=imdimY, slice_width=slice_width, offset=offset),
                            [i for i, _ in enumerate(tilts_in_range)])
        p.join()
        # map could also iter over range(len(self._list)), however i explicitly zip _tilt_angles and the _list to
        # make sure that they are both present

        if show_progress_bar:
            print('========== All processes finished alignment ============')

        # get the stack as a pytom volume
        stack = vol(imdim, imdim, len(tilts_in_range))
        self._retrieve_shared_stack(stack)

        phiStack = vol(1, 1, len(tilts_in_range))
        phiStack.setAll(0.0)

        thetaStack = vol(1, 1, len(tilts_in_range))
        thetaStack.setAll(0.0)

        offsetStack = vol(1, 2, len(tilts_in_range))
        offsetStack.setAll(0.0)

        # fill the stacks
        for i, result in results:  # why round(int())?
            phiStack(result[0], 0, 0, i)
            thetaStack(result[1], 0, 0, i)
            offsetStack(int(round(result[2][0])), 0, 0, i)
            offsetStack(int(round(result[2][1])), 0, 1, i)

        return [stack, phiStack, thetaStack, offsetStack]

    def to_projection_stack_gpu(self, weighting=0, binning=1, low_pass_freq=0.9, apply_circle_filter=True,
                                scale_factor_particle=1., angle_specimen=0., particle_diameter=None, pixel_size=2.62,
                                tilt_range=None, show_progress_bar=False, verbose=False):
        from pytom.agnostic.io import read
        from pytom.agnostic.tools import taper_edges, paste_in_center
        from pytom.agnostic.filter import circle_filter, ramp_filter, exact_filter, ellipse_filter, bandpass_circle
        import pytom.voltools as vt
        from pytom.gpu.initialize import xp, device
        from pytom.agnostic.transform import resize
        from pytom.tools.ProgressBar import FixedProgBar
        from pytom.basic.transformations import general_transform_matrix
        from pytom.lib.pytom_numpy import vol2npy
        import numpy as np

        if low_pass_freq > 1.:
            low_pass_freq = 1.
            print("Warning: lowpassFilter > 1 - set to 1 (=Nyquist)")

        if tilt_range is None:
            tilt_range = (int(round(min(self._tilt_angles))), int(round(max(self._tilt_angles))))

        imdimX = self[0].getXSize()
        imdimY = self[0].getYSize()

        if binning > 1:
            imdimX = int(float(imdimX) / float(binning) + .5)
            imdimY = int(float(imdimY) / float(binning) + .5)

        imdim = max(imdimY, imdimX)

        # get a list of tilt angles that are in the range
        tilts_in_range = [t for t in self._tilt_angles if tilt_range[0] <= int(round(t)) <= tilt_range[1]]

        # find the offset of the lower limit of provided tilt-range
        offset = 0
        while int(round(self._tilt_angles[offset])) < tilt_range[0]:
            offset += 1

        # prep stacks for output
        stack = xp.zeros((imdim, imdim, len(tilts_in_range)), dtype=xp.float32)
        phiStack = xp.zeros((len(tilts_in_range)), dtype=xp.float32)
        thetaStack = xp.zeros((len(tilts_in_range)), dtype=xp.float32)
        offsetStack = xp.zeros((len(tilts_in_range), 2), dtype=xp.int32)

        if particle_diameter is None: # if no particle diameter provided we set slice width relative to full image size
            slice_width = 1
        elif particle_diameter <= 0: # Particle diameter should be bigger than 0
            raise ValueError("particle diameter must be greater than 0")
        else:  # ideal slice width for exact filter is relative to particle diameter
            slice_width = (pixel_size * binning * imdim) / particle_diameter

        # pre-determine analytical weighting function and lowpass for speedup
        if weighting == -1:
            weightSlice = xp.fft.fftshift(ramp_filter(imdim, imdim))

        if apply_circle_filter:
            # use ellipse filter if images are not squared?
            # circleSlice = xp.fft.fftshift(ellipse_filter(imdim, imdim, circleFilterRadiusX, circleFilterRadiusY))
            circle_filter_radius = imdim // 2
            circleSlice = xp.fft.fftshift(circle_filter(imdim, imdim, circle_filter_radius))
        else:
            circleSlice = xp.ones((imdim, imdim))

        # initialize taper_mask variable
        taper_mask = None
        outputImage = None

        if verbose:
            print('projections to go over: ', [i for i, _ in enumerate(tilts_in_range)])

        if show_progress_bar:
            progressBar = FixedProgBar(0, len(tilts_in_range), 'Creating aligned and weighted projections')
            progressBar.update(0)

        for i, (tilt_angle, projection) in enumerate(zip(self._tilt_angles, self._list)):

            if not (tilt_range[0] <= int(round(tilt_angle)) <= tilt_range[1]):
                continue

            if projection._filename.split('.')[-1] == 'st':
                # image = read(projection.getFilename(), deviceID=device)[:, :, i].squeeze()
                image = read(projection.getFilename())[:, :, i].squeeze()
            else:
                # image = read(projection.getFilename(), deviceID=device).squeeze()
                image = read(projection.getFilename()).squeeze()

            # x,y needs to be flipped for correct orientations after reading
            # reading with pytom.basic.files not same a reading with pytom.agnostic.io
            image = image.T
            # TODO this also changes image dimension information, needs to be adjusted

            if binning > 1:
                image = resize(volume=image, factor=1 / float(binning))

            # 1 -- normalize to contrast - subtract mean and norm to mean
            immean = image.mean()
            image = (image - immean) / immean  # maybe only division by mean and not subtraction

            # 2 -- smoothen borders to prevent high contrast oscillations
            if taper_mask is None:
                image, taper_mask = taper_edges(image, imdim // 30)
            else:
                image *= taper_mask

            # transform projection according to tilt alignment
            transX = projection.getAlignmentTransX() / binning
            transY = projection.getAlignmentTransY() / binning
            rot = projection.getAlignmentRotation()
            mag = projection.getAlignmentMagnification() * scale_factor_particle
            order = projection.getOperationOrder()

            # 3 -- square if needed
            if imdimY != imdimX:
                newImage = xp.zeros((imdim, imdim), dtype=xp.float32)
                image = paste_in_center(image, newImage)

            # TODO lots of arrays created, can this be more effiecient?
            inputImage = xp.expand_dims(image, 2)

            if outputImage is None:
                outputImage = xp.zeros_like(inputImage, dtype=xp.float32)
            else:
                outputImage *= 0

            center = ((inputImage.shape[0] - 1) / 2., (inputImage.shape[1] - 1) / 2., 0) if \
                order == [1, 0, 2] else None

            # TODO provide operation order directly to vt and let vt calculate transform matrix
            mtx = general_transform_matrix(inputImage.shape, rot=[rot, 0, 0], shift=[transX, transY, 0],
                                           scale=[mag, mag, 1], order=order, center=center)
            # apply matrix inversion because: (1) dot(a, mtx) == (2) dot(inv(mtx), a)
            # pytom uses (2) and voltools uses (1)
            mtx = np.linalg.inv(vol2npy(mtx._matrix).copy(order='F'))

            # pass the matrix to voltools
            vt.transform(inputImage.astype(xp.float32), matrix=mtx, output=outputImage,
                         device=device, interpolation='filt_bspline')
            # vt.transform(inputImage.astype(xp.float32), rotation=[0, 0, rot], rotation_order='rxyz', output=outputImage,
            #              device=device, translation=[transX, transY, 0], scale=[mag, mag, 1],
            #              interpolation='filt_bspline')
            del inputImage
            image = outputImage.squeeze()

            # 5 -- lowpass filter (optional)
            if low_pass_freq > 0:  # no mutliply in fourier space??
                sigma = 0.015
                image = bandpass_circle(image, high=(low_pass_freq * imdim // 2), sigma=(sigma * imdim))
                # image = bandpass_circle(image, high=low_pass_freq * imdim // 2, sigma=low_pass_freq / 5. *
                #                                                                                 imdim)

            # 6 -- smoothen once more to avoid edges
            image *= taper_mask

            # 7 -- weighting
            if weighting == -1:
                image = xp.fft.ifftn(xp.fft.fftn(image) * weightSlice * circleSlice).real
            elif weighting == 1:
                weightSlice = xp.fft.fftshift(exact_filter(self._tilt_angles, tilt_angle, imdim, imdim, slice_width))
                image = xp.fft.ifftn(xp.fft.fftn(image) * weightSlice * circleSlice).real

            # thetaStack[ii] = float(round(tilt_angle - angle_specimen))  # TODO why round + int?
            thetaStack[i - offset] = tilt_angle - angle_specimen
            offsetStack[i - offset, :] = xp.array([int(round(projection.getOffsetX())),
                                          int(round(projection.getOffsetY()))])
            stack[:, :, i - offset] = image

            if show_progress_bar:
                progressBar.update(i - offset)

            if verbose:
                print(tilt_angle, projection.getOffsetX(), projection.getOffsetY())

        return [stack, phiStack, thetaStack, offsetStack]
