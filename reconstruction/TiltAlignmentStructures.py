'''
Structures for alignment of tilt series using fiducial markers
'''
from pytom_volume import vol, read, transform
from pytom_numpy import vol2npy
import numpy
import os
from pytom.reconstruction.reconstructionStructures import Projection, ProjectionList
from pytom.reconstruction.tiltAlignmentFunctions import markerResidual, alignmentFixMagRot
from pytom.basic.structures import PyTomClass
from pytom_volume import vol
from pytom_numpy import vol2npy
from pytom.gui.guiFunctions import loadstar

class TiltSeries(PyTomClass):
    """
    tilt series

    @author: FF
    """

    def __init__(self, tiltSeriesName='', TiltAlignmentParas='', alignedTiltSeriesName='',
                 markerFileName='', firstProj=1, lastProj=None, projIndices=None,
                 tiltSeriesFormat='em', verbose=False, alignmentResultFile=''):
        """
        initialize TiltSeries 

        @param tiltSeriesName: Name of Tilt series files
        @type tiltSeriesName: string
        @param TiltAlignmentParas: Tilt series alignment parameters
        @type TiltAlignmentParas: TiltAlignmentParameters
        @param firstProj: Index of first projection
        @type firstProj: int
        @param lastProj: Index of last projection
        @type lastProj: int
        @param projIndices: Projection indices - can be used to ommit specific projections from tilt series
        @type projIndices: array
        @param markerFileName: Name of file containing marker coordinates in projections
        @type markerFileName: string
        @param alignedTiltSeriesName: Name of Aligned Tilt series files
        @type alignedTiltSeriesName: string
        @param tiltSeriesFormat: format of tilt series (em/mrc)
        @type tiltSeriesFormat: string
        @param verbose: Verbose mode on / off
        @type verbose: bool
        @author: FF
        """
        self._tiltSeriesName = tiltSeriesName
        self.verbose = verbose
        print(TiltAlignmentParas)
        if TiltAlignmentParas:
            self.init_vars_markerbased_alignment(tiltSeriesName, TiltAlignmentParas, alignedTiltSeriesName,
                                                 markerFileName, firstProj=firstProj, lastProj=lastProj,
                                                 tiltSeriesFormat=tiltSeriesFormat, verbose=False)

    def init_vars_markerbased_alignment(self, tiltSeriesName, TiltAlignmentParas, alignedTiltSeriesName, markerFileName,
                                        firstProj=1, lastProj=None, tiltSeriesFormat='em', verbose=False):
        if tiltSeriesName:
            self._firstProj = firstProj
            self._lastProj = lastProj

            folder = os.path.dirname(tiltSeriesName)
            prefix = os.path.basename(tiltSeriesName)
            files = [line for line in os.listdir(folder) if line.endswith(tiltSeriesFormat) and line.startswith(prefix)]
            self._projIndices = [line.split('_')[-1].split('.')[0] for line in files ]
            self._projIndices.sort(key=int)
            self._lenPI = len(self._projIndices)

            self._projIndices = numpy.array(self._projIndices)
            self._tiltSeriesFormat = tiltSeriesFormat
            self._TiltAlignmentParas = TiltAlignmentParas
            self._alignedTiltSeriesName = alignedTiltSeriesName

            print(markerFileName, len(files))
            if markerFileName.endswith('.em'): self.mf = vol2npy(read(markerFileName))

            elif markerFileName.endswith('.txt'): self.mf = self.txt2markerfile(markerFileName, len(files))
            else: raise Exception('Unknown file type for markerfile.\n Please submit either a .txt file or a .em file')
            #print(self.mf[0,:,0])
            #print self.mf[0,:,0]
            
            # set Projection List
            self._firstIndex, self._lastIndex = -1, len(self._projIndices)
            projs = []


            if tiltSeriesFormat != 'st':
                for cnt, ii in enumerate(self._projIndices):

                    if int(ii) < self._firstProj or int(ii) > self._lastProj: continue
                    if self._firstIndex < 0: self._firstIndex = cnt
                    self._lastIndex = cnt+1
                    fname = tiltSeriesName + "_" + str( ii ) + "." + tiltSeriesFormat

                    if alignedTiltSeriesName:
                        proj = Projection(filename=fname,
                                          alignedFilename=alignedTiltSeriesName + "_" + str(ii) + "." +tiltSeriesFormat,
                                          index=ii, tiltAngle=self.mf[0,cnt,0],
                                          offsetX=0., offsetY=0.,
                                          alignmentTransX=0., alignmentTransY=0.,
                                          alignmentRotation=0., alignmentMagnification=1.)
                        

                    else:
        
                        proj = Projection(filename=fname,
                                          alignedFilename=None,
                                          index=ii, tiltAngle=self.mf[0,cnt,0],
                                          offsetX=0., offsetY=0.,
                                          alignmentTransX=0., alignmentTransY=0.,
                                          alignmentRotation=0., alignmentMagnification=1.)
                    projs.append(proj)
                    if self.verbose:
                        print(("Projection " + fname + " appended ..."))
            else:
                tmp = tiltSeriesName.split('.')[-1]
                if tmp != 'st':
                    fname = tiltSeriesName + ".st"
                else:
                    fname = tiltSeriesName

                for cnt, ii in enumerate(self._projIndices):
                    if int(ii) < self._firstProj or int(ii) > self._lastProj: continue
                    if self._firstIndex < 0: self.firstIndex = cnt
                    self._lastIndex = cnt+1

                    if alignedTiltSeriesName:
                        proj = Projection(filename=fname,
                                          alignedFilename=alignedTiltSeriesName + "_" + str(ii) + tiltSeriesFormat,
                                          index=ii, tiltAngle=self.mf[0,ii,0],
                                          offsetX=0., offsetY=0.,
                                          alignmentTransX=0., alignmentTransY=0.,
                                          alignmentRotation=0., alignmentMagnification=1.)
                        
                    else:
                        
                        proj = Projection(filename=fname,
                                          alignedFilename=None,
                                          index=ii, tiltAngle=self.mf[0,ii,0],
                                          offsetX=0., offsetY=0.,
                                          alignmentTransX=0., alignmentTransY=0.,
                                          alignmentRotation=0., alignmentMagnification=1.)
                    projs.append(proj)
                    if self.verbose:
                        print(("Projection " + fname + " appended ..."))
            self._ProjectionList = ProjectionList(projs)
            self._imdim = projs[0].getDimensions()[0]

            # read markerFile if set
            self._markerFileName = markerFileName
            self._projIndices = self._projIndices[self._firstIndex:self._lastIndex]

        self._Markers = []

        if markerFileName:
            if self.verbose:
                print(("reading marker file: " + str(markerFileName)))
            self._Markers = self.readMarkerFile(markerFileName)
        else:
            if self.verbose:
                print("No marker file specified")

    def info(self):
        """
        print current projections and their alignment

        @author: FF
        """
        tline = ""
        for (ii, projection) in enumerate(self._ProjectionList):
            tiltAngle = projection._tiltAngle
            transX = -projection._alignmentTransX
            transY = -projection._alignmentTransY
            rot = -(projection._alignmentRotation + 90.)
            mag = projection._alignmentMagnification
            tline = tline + ("%3d: " % ii)
            tline = tline + ("%15s; " % projection._filename)
            tline = tline + ("tiltAngle=%9.3f; " % tiltAngle)
            tline = tline + ("transX=%9.3f; " % transX)
            tline = tline + ("transY=%9.3f; " % transY)
            tline = tline + ("rot=%9.3f; " % rot)
            tline = tline + ("mag=%9.3f\n" % mag)
        print(tline)

    def toXML(self):
        """
        create xml of tilt series and its alignment

        @author: FF
        """
        from lxml import etree

        tiltseries_element = etree.Element("TiltSeries")
        for proj in self._ProjectionList:
            tiltseries_element.append(proj.toXML())

        return tiltseries_element

    def createEmptyProjections(self, imdim=2048):
        """
        @param imdim: image dimensions in pixels (default: 2048)
        @type imdim: int

        @author: FF
        """
        projs = []
        for ii in self._projIndices:
            proj = Projection(filename=None,
                              alignedFilename=None,
                              index=ii, tiltAngle=None,
                              offsetX=0., offsetY=0.,
                              alignmentTransX=0., alignmentTransY=0.,
                              alignmentRotation=0., alignmentMagnification=1.)
            projs.append(proj)
        self._ProjectionList = ProjectionList(projs)
        self._imdim = imdim

    def removeProjection(self, iremove):
        """
        remove projection with specified index from TiltSeries and corresponding marker

        @param iremove: index of removed projection (in FileName)
        @type iremove: int

        @author: FF
        """
        # check that dims of Projections and Markers are the same
        nprojProjs = len(self._ProjectionList._list)
        nprojMarker = len(self.Markers[0].xProj)
        if (nprojProjs != nprojMarker):
            "Houston: we have a problem!"
            "Numbers of projections in Markers and Projections do not match"
        kk = -1
        for proj in self._ProjectionList._list:
            kk = kk + 1
            ii = proj._index
            if (ii == iremove):
                break
        self._ProjectionList._list.pop(kk)
        self._projIndices.remove(iremove)
        for Marker in self.Markers:
            Marker.xProj.pop(kk)
            Marker.yProj.pop(kk)
            Marker._projIndices.remove(iremove)
        if self.verbose:
            print(("Projection " + str(iremove) + " removed from TiltSeries"))

    def txt2markerfile(self, filename, num_tilt_images):
        data = loadstar(filename)
        datalen = data.shape[0]
        id = data[-1][0]
        if not id:
            id = data[-2][0]
        id = int(round(id)) +1
        start = 0
        for i in range(len(data)):
            if data[i][0] == 0:
                start += 1
        print(start)

        x, y = int(round(id)), datalen // id#, num_tilt_images
        markerfile = data.reshape(x, y, 4)[:, :, 1:].transpose(2, 1, 0)
        #markerfile[1:,:,:] = markerfile[1:,:,:][::-1,:]
        return  markerfile

    def readMarkerFile(self, markerFileName):
        """
        read MarkerFile and update self._Markers

        @param markerFileName: File name of marker file
        @type markerFileName: string

        @author: FF
        """
        from math import pi

        if markerFileName.endswith('.em') or markerFileName.endswith('.mrc'):
            markerFileVol = read(markerFileName)
            nproj = markerFileVol.sizeY()
            nproj -= self._firstIndex
            nproj -= self._lenPI - self._lastIndex
            nmark = markerFileVol.sizeZ()
            markerFile = vol2npy(markerFileVol)
            markerFile = markerFile[:, self._firstIndex:self._lastIndex, :]

        else:
            markerFile = self.txt2markerfile(markerFileName, len(self._ProjectionList._list) )
            nproj = markerFile.shape[1]
            nproj -= self._firstIndex
            nproj -= self._lenPI - self._lastIndex
            nmark = markerFile.shape[2]
            markerFile = markerFile[:, self._firstIndex:self._lastIndex, :]

        # make sure that nproj matches number of Projections in self._ProjectionList
        if (nproj != len(self._ProjectionList._list)):
            print("Number of projections specified in TiltSeries and MarkerFileName do not match!")
            print(("  Markerfile: " + str(nproj)))
            print(("  TiltSeries: " + str(len(self._ProjectionList._list))))
            print("Please fix!")

        # check that tilt angles in marker file and projections are the same
        for (iproj, proj) in enumerate(self._ProjectionList._list):
            tiltAngle = markerFile[0, iproj, 0]
            proj._tiltAngle = tiltAngle
            if abs(tiltAngle - proj._tiltAngle) > 0.01:
                print(("Warning: tilt angles in Projection " + str(iproj + 1)
                      + " differ in markerFile"))
                print(("MarkerFile: %4.1f" % tiltAngle) + (" vs. image Header: %4.1f" % proj._tiltAngle))

        # delete old markerFiles
        if (len(self._Markers) > 0):
            print("overwriting pre-loaded Markers")
        self._Markers = []
        for imark in range(0, nmark):
            self._Markers.append(Marker(self._projIndices))
            x = markerFile[1, 0:nproj, imark]
            y = markerFile[2, 0:nproj, imark]
            print(imark, x, y)
            self._Markers[imark].set_xProjs(x)
            self._Markers[imark].set_yProjs(y)
        return self._Markers

    def readIMODwimp(self, markerFileName, prexgfile=None, tltfile=None, preBin=1, verbose=False):
        """
        read Marker coordinates from IMOD wimp file

        @param markerFileName: name of winp file
        @type markerFileName: str
        @param prexgfile: file containing pre-shifts (IMOD way of doing things)
        @type prexgfile: str
        @param tltfile: file containing tilt angles - typically end on .tlt or .rawtlt
        @type tltfile: str
        @param preBin: pre-binning in IMOD
        @type preBin: int
        @param verbose: verbose?
        @type verbose: bool

        @author: FF
        """
        from pytom.reconstruction.tiltAlignmentFunctions import readIMODmarkerfile
        from pytom.reconstruction.tiltAlignmentFunctions import getIMODpreshifts
        from pytom.reconstruction.tiltAlignmentFunctions import applyPreshiftsToMarkers
        self._Markers = readIMODmarkerfile(markerFileName, binning=preBin)
        if verbose:
            print("########################")
            print("Markers before prexgfile")
            print("########################")
            for marker in self._Markers:
                marker.info()
        if tltfile:
            self.getTiltAnglesFromIMODfile(tltfile=tltfile)

        if prexgfile:
            (shiftX, shiftY) = getIMODpreshifts(prexgfile)
            applyPreshiftsToMarkers(markers=self._Markers, shiftX=shiftX, shiftY=shiftY)
            if verbose:
                print("########################")
                print("Markers after prexgfile ")
                print("########################")
                for marker in self._Markers:
                    marker.info()

    def getTiltAnglesFromIMODfile(self, tltfile):
        """
        get tiltangles from IMOD file and update them in TiltSeries

        @param tltfile: file containing tilt angles - typically end on .tlt or .rawtlt
        @type tltfile: str
        @return: tiltAngles
        @rtype: list

        @author: FF
        """
        from pytom.reconstruction.tiltAlignmentFunctions import readIMODtiltAngles
        from math import sin, cos, pi

        tiltAngles = readIMODtiltAngles(tltfile)
        if len(tiltAngles) != len(self._ProjectionList):
            print("Number of tilt angles in tiltfile:   " + len(tiltAngles))
            print("Number of tilt angles in tiltseries: " + len(self._ProjectionList))
            raise IndexError('Number of tilt angles in tltfile does not match TiltSeries')
        for (iproj, proj) in enumerate(self._ProjectionList):
            proj._tiltAngle = tiltAngles[iproj]
        return tiltAngles

    def writeMarkerFile(self, markerFileName):
        """
        write MarkerFile using data in self._Markers

        @param markerFileName: File name of marker file
        @type markerFileName: string

        @author: FF
        """
        markerFileVol = vol(12, len(self._ProjectionList._list), len(self._Markers))
        markerFileVol.setAll(0.)
        for (imark, Marker) in enumerate(self._Markers):
            for (itilt, proj) in enumerate(self._ProjectionList._list):
                markerFileVol.setV(Marker.get_xProj(itilt), 1, itilt, imark)
                markerFileVol.setV(Marker.get_yProj(itilt), 2, itilt, imark)
                if imark == 0:
                    markerFileVol.setV(int(round(proj._tiltAngle)), 0, int(itilt), int(imark))
        markerFileVol.write(markerFileName)


    def write_aligned_projs(self, weighting=None, lowpassFilter=None, binning=1, verbose=False, write_images=True):
        """
        write aligned (and weighted) projections

        @param weighting: weighting (<0: analytical weighting, >1 exact weighting (value corresponds to object diameter\
                    in pixel AFTER binning)
        @type weighting: float
        @param lowpassFilter: lowpass filter (in Nyquist)
        @type lowpassFilter: float
        @param binning: binning (default: 1 = no binning). binning=2: 2x2 pixels -> 1 pixel, \
        binning=3: 3x3 pixels -> 1 pixel, etc.
        @type binning: int or float
        @param verbose: verbose mode?
        @type verbose: bool

        @author: FF
        """
        from pytom.reconstruction.writeAlignedProjections import writeAlignedProjections

        writeAlignedProjections(self, weighting=weighting, lowpassFilter=lowpassFilter, binning=binning,
                                verbose=verbose, write_images=write_images)

    def write_projections(self, filename, binning=1, lowpassFilter=None, verbose=False):
        """
        write projection files to disk
        @param filename: projections will be stored as <filename>_index.em
        @type filename: str
        @param binning: binning (default: 1 = no binning). binning=2: 2x2 pixels -> 1 pixel, \
        binning=3: 3x3 pixels -> 1 pixel, etc.
        @type binning: int
        @param lowpassFilter: lowpass filter (in Nyquist)
        @type lowpassFilter: float
        @param verbose: verbose mode?
        @type verbose: bool

        @author: FF
        """
        from pytom.basic.files import read_em, write_em
        from pytom.basic.filter import filter as filterFunction
        import pytom_freqweight
        from pytom.basic.transformations import resize

        if binning:
            imdim = int(float(self._imdim) / float(binning) + .5)
        else:
            imdim = self._imdim

        # design lowpass filter
        if lowpassFilter:
            if lowpassFilter > 1.:
                lowpassFilter = 1.
                print("Warning: lowpassFilter > 1 - set to 1 (=Nyquist)")
            # weighting filter: arguments: (angle, cutoff radius, dimx, dimy,
            lpf = pytom_freqweight.weight(0.0, lowpassFilter*imdim/2, imdim, imdim/2+1, 1, lowpassFilter/5.*imdim)

        for (ii, projection) in enumerate(self._ProjectionList):
            if projection._filename.split('.')[-1] == 'st':
                from pytom.basic.files import EMHeader, read
                header = EMHeader()
                header.set_dim(x=imdim, y=imdim, z=1)
                idx = projection._index
                if verbose:
                    print("reading in projection %d" % idx)
                image = read(file=projection._filename, subregion=[0, 0, idx - 1, self._imdim, self._imdim, 1],
                             sampling=[0, 0, 0], binning=[0, 0, 0])
            else:
                # read projection files
                (image, header) = read_em(projection._filename)
            if not (binning == 1) or (binning == None):
                image = resize(volume=image, factor=1 / float(binning))[0]
            if lowpassFilter:
                filtered = filterFunction(volume=image, filterObject=lpf, fourierOnly=False)
                image = filtered[0]

            tiltAngle = projection._tiltAngle
            if verbose:
                print("tiltAngle=%2.2f" % tiltAngle)
            header.set_tiltangle(tiltAngle)
            newFilename = (filename + "_" + str(projection.getIndex()) + '.em')
            write_em(filename=newFilename, data=image, header=header)

    def reconstructVolume(self, dims=[512, 512, 128], reconstructionPosition=[0, 0, 0], binning=1, alignResultFile=''):
        """
        reconstruct a single 3D volume from weighted and aligned projections

        @param dims: 3D Size of reconstructed volume
        @type dims: 3d list
        @param binning: down-sizing of projections prior to 3d rec
        @type binning: int
        @param reconstructionPosition: offset center of reconstruction
        @type reconstructionPosition: 3-dim list

        @author: FF
        """
        # set alignedProjectionList
        projs = []
        if not alignResultFile:
            for (kk, ii) in enumerate(self._projIndices):
                tiltAngle = self._ProjectionList[kk]._tiltAngle
                proj = Projection(filename=self._alignedTiltSeriesName + "_" + str(ii) + "." + self._tiltSeriesFormat,
                                  alignedFilename=self._alignedTiltSeriesName + "_" + str(
                                      ii) + "." + self._tiltSeriesFormat,
                                  index=ii, tiltAngle=tiltAngle,
                                  offsetX=0., offsetY=0.,
                                  alignmentTransX=0., alignmentTransY=0.,
                                  alignmentRotation=0., alignmentMagnification=1.)
                projs.append(proj)

        self._alignedProjectionList = ProjectionList(projs)
        # reconstruct tomogram
        vol_bp = self._alignedProjectionList.reconstructVolume(dims=dims, alignResultFile=alignResultFile,
                                                               reconstructionPosition=reconstructionPosition,
                                                               binning=binning)
        return vol_bp


class TiltAlignment:
    """
    class containing alignment of a tilt series (e.g., shifts, rotations)

    @author: FF
    """

    def __init__(self, TiltSeries_, verbose=False):
        """
        @param TiltSeries_: Tilt series 
        @type TiltSeries_: TiltSeries
        @param verbose: verbose mode?
        @type verbose: L{bool}

        @author: FF
        """
        self.TiltSeries_ = TiltSeries_
        self._ntilt = len(TiltSeries_._projIndices)
        self._projIndices = TiltSeries_._projIndices
        self._alignmentTransX = numpy.array(self._ntilt * [0.])
        self._alignmentTransY = numpy.array(self._ntilt * [0.])
        self._alignmentRotations = numpy.array(self._ntilt * [0.])
        self._alignmentMagnifications = numpy.array(self._ntilt * [1.])
        self._alignmentBeamTilt = 0.
        # get tilt angles from projections, stored in self._tiltAngles, self._cTilt, self._sTilt
        self._tiltAngles = numpy.array(self._ntilt * [0.])
        self._cTilt = numpy.array(self._ntilt * [0.])
        self._sTilt = numpy.array(self._ntilt * [0.])
        self.getTiltAnglesFromTiltSeries(TiltSeries_)
        self.Psi = 0.
        self._dmagnfoc = 0.
        self._drotfoc = 0.
        self._dbeam = None
        self._Markers = []
        self.getMarkersFromTiltSeries(TiltSeries_)
        self.verbose = verbose

    def info(self):
        """
        print current alignment

        @author: FF
        """
        tline = ""
        for ii in range(0, self._ntilt):
            tiltAngle = self._tiltAngles[ii]
            transX = -self._alignmentTransX[ii]
            transY = -self._alignmentTransY[ii]
            rot = -(self._alignmentRotations[ii] + 90.)
            mag = self._alignmentMagnifications[ii]

            tline = tline + ("%3d: " % ii)
            tline = tline + ("tiltAngle=%8.3f; " % tiltAngle)
            tline = tline + ("transX=%9.3f; " % transX)
            tline = tline + ("transY=%9.3f; " % transY)
            tline = tline + ("rot=%9.3f; " % rot)
            tline = tline + ("mag=%9.3f\n" % mag)
        tline = tline + ("================================================================================\n")
        for (imark, marker) in enumerate(self._Markers):
            tline = tline + ("Marker %2d =" % imark)
            tline = tline + ("%6.1f, " % marker.get_r()[0])
            tline = tline + ("%6.1f, " % marker.get_r()[1])
            tline = tline + ("%6.1f\n" % marker.get_r()[2])

        print(tline)

    def resetAlignmentCenter(self):
        """
        set alignment center according to dimensions of images
        """
        cent = self.TiltSeries_._TiltAlignmentParas.cent
        imdim = self.TiltSeries_._imdim
        print(imdim)
        if cent[0] != imdim//2+1:
            #rint "Centers do not match: cent="+str(cent)+", imdim="+str(imdim)
            self.TiltSeries_._TiltAlignmentParas.cent = [imdim//2+1, imdim//2+1]

    def getMarkersFromTiltSeries(self, TiltSeries_):
        """
        get marker coordinates from TiltSeries and update in TiltAlignment

        @param TiltSeries_: tilt series
        @type TiltSeries_: TiltSeries

        @author: FF
        """
        self._Markers = TiltSeries_._Markers
        return self._Markers

    def setMarkersInTiltSeries(self, TiltSeries_):
        """
        set marker coordinates in TiltSeries

        @param TiltSeries_: tilt series
        @type TiltSeries_: TiltSeries

        @author: FF
        """
        TiltSeries_._Markers = self._Markers

    def printTransRotMag(self):
        """
        print TransRotMag
        @author: FF
        """

    def getTranslationsFromTiltSeries(self, TiltSeries_):
        """
        get translations from TiltSeries and update

        @param TiltSeries_: tilt series
        @type TiltSeries_: L{pytom.reconstruction.TiltAlignmentStructures.TiltSeries}
        @return: x- and y-translation
        @rtype: list

        @author: FF
        """
        # initialize alignment in separate array - easier for optimization
        nprojs = len(TiltSeries_._ProjectionList._list)
        self._alignmentTransX = nprojs * [0.]
        self._alignmentTransY = nprojs * [0.]
        for (kk, proj) in enumerate(TiltSeries_._ProjectionList._list):
            self._alignmentTransX[kk] = proj.getAlignmentTransX()
            self._alignmentTransY[kk] = proj.getAlignmentTransY()
        return self._alignmentTransX, self._alignmentTransY

    def setTranslationsInTiltSeries(self, TiltSeries_):
        """
        set translations in TiltSeries according to stored values

        @param TiltSeries_: tilt series
        @type TiltSeries_: L{pytom.reconstruction.TiltAlignmentStructures.TiltSeries}

        @author: FF
        """
        for (kk, proj) in enumerate(TiltSeries_._ProjectionList._list):
            proj.setAlignmentTransX(self._alignmentTransX[kk])
            proj.setAlignmentTransY(self._alignmentTransY[kk])

    def getRotationsFromTiltSeries(self, TiltSeries_):
        """
        get rotations from TiltSeries and update

        @param TiltSeries_: tilt series
        @type TiltSeries_: L{pytom.reconstruction.TiltAlignmentStructures.TiltSeries}
        @return: rotation angles of tilt axes
        @rtype: array

        @author: FF
        """
        # initialize alignment in seperate array - easier for optimization
        self._alignmentRotations = numpy.array( len(TiltSeries_._ProjectionList._list) * [0.] )
        for (kk, proj) in enumerate(TiltSeries_._ProjectionList._list):
            self._alignmentRotations[kk] = proj.getAlignmentRotation()
        return self._alignmentRotations

    def setRotationsInTiltSeries(self, TiltSeries_):
        """
        set rotations in TiltSeries 

        @param TiltSeries_: tilt series
        @type TiltSeries_: L{pytom.reconstruction.TiltAlignmentStructures.TiltSeries}

        @author: FF
        """
        for (kk, proj) in enumerate(TiltSeries_._ProjectionList._list):
            proj.setAlignmentRotation(self._alignmentRotations[kk])

    def getMagnificationsFromTiltSeries(self, TiltSeries_):
        """
        get magnifications from TiltSeries and update

        @param TiltSeries_: tilt series
        @type TiltSeries_: L{pytom.reconstruction.TiltAlignmentStructures.TiltSeries}

        @author: FF
        """
        # initialize alignment in seperate array - easier for optimization
        self._alignmentMagnifications = len(TiltSeries_._ProjectionList._list) * [1.]
        for (kk, proj) in enumerate(TiltSeries_._ProjectionList._list):
            self._alignmentMagnifications[kk] = proj.getAlignmentMagnification()
        return self._alignmentMagnifications

    def setMagnificationsInTiltSeries(self, TiltSeries_):
        """
        set magnifications in TiltSeries

        @param TiltSeries_: tilt series
        @type TiltSeries_: L{pytom.reconstruction.TiltAlignmentStructures.TiltSeries}

        @author: FF
        """
        kk = 0
        for proj in TiltSeries_._ProjectionList._list:
            proj.setAlignmentMagnification(self._alignmentMagnifications[kk])
            kk = kk + 1

    def getTiltAnglesFromTiltSeries(self, TiltSeries_):
        """
        get tilt angles from TiltSeries -
        also computes sin and cos and stores both

        @param TiltSeries_: tilt series
        @type TiltSeries_: L{pytom.reconstruction.TiltAlignmentStructures.TiltSeries}

        @rtype: array

        @author: FF
        """
        from math import sin, cos, pi

        if len(TiltSeries_._projIndices) != self._ntilt:
            print("Houston: We have a problem!")
        # initialize alignment in seperate array - easier for optimization
        # sin and cos
        self._tiltAngles = numpy.array(self._ntilt * [0.])
        self._cTilt = numpy.array(self._ntilt * [0.])
        self._sTilt = numpy.array(self._ntilt * [0.])
        for (kk, proj) in enumerate(TiltSeries_._ProjectionList._list):
            the = proj.getTiltAngle()
            self._tiltAngles[kk] = the
            self._sTilt[kk] = sin(the / 180. * pi)
            self._cTilt[kk] = cos(the / 180. * pi)
        return self._tiltAngles

    def getTiltAnglesFromIMODfile(self, tltfile):
        """
        get tiltangles from IMOD file

        @param tltfile: file containing tilt angles - typically end on .tlt or .rawtlt
        @type tltfile: str
        @return: tiltAngles
        @rtype: list

        @author: FF
        """
        from pytom.reconstruction.tiltAlignmentFunctions import readIMODtiltAngles
        from math import sin, cos, pi

        tiltAngles = readIMODtiltAngles(tltfile)
        self._ntilt = len(tiltAngles)
        self._tiltAngles = numpy.array(self._ntilt * [0.])
        self._cTilt = numpy.array(self._ntilt * [0.])
        self._sTilt = numpy.array(self._ntilt * [0.])
        for ii in range(0, len(tiltAngles)):
            self._tiltAngles[ii] = tiltAngles[ii]
            self._sTilt[ii] = sin(tiltAngles[ii] / 180. * pi)
            self._cTilt[ii] = cos(tiltAngles[ii] / 180. * pi)
        return self._tiltAngles

    def randomize(self, amplitude):
        """
        add random distortion to current alignment

        @param amplitude: amplitude of distortion
        @type amplitude: L{float}
        """
        from random import gauss
        # update alignment from projections
        self.getMarkersFromTiltSeries(self.TiltSeries_)
        self.getTranslationsFromTiltSeries(self.TiltSeries_)
        self.getRotationsFromTiltSeries(self.TiltSeries_)
        self.getMagnificationsFromTiltSeries(self.TiltSeries_)

        # add random amplitude
        # markers
        irefmark = self.TiltSeries_._TiltAlignmentParas.irefmark
        for (imark, mark) in enumerate(self._Markers):
            if (imark + 1) != irefmark:
                r = mark.get_r()
                for ii in range(0, 3):
                    r[ii] = gauss(0, amplitude) + r[ii]
                self._Markers[imark].set_r(r)

        ireftilt = self.TiltSeries_._TiltAlignmentParas.ireftilt
        # translations, rotations, magnifications
        for itilt in range(0, self._ntilt):
            self._alignmentTransX[ii] = self._alignmentTransX[ii] + gauss(0, amplitude)
            self._alignmentTransY[ii] = self._alignmentTransY[ii] + gauss(0, amplitude)
            self._alignmentRotations[ii] = self._alignmentRotations[ii] + gauss(0, amplitude)
            if (self._projIndices[itilt] != ireftilt):
                self._alignmentMagnifications[ii] = self._alignmentMagnifications[ii] + gauss(0, amplitude)

        # update alignment from projections
        self.setMarkersInTiltSeries(self.TiltSeries_)
        self.getTranslationsFromTiltSeries(self.TiltSeries_)
        self.getRotationsFromTiltSeries(self.TiltSeries_)
        self.getMagnificationsFromTiltSeries(self.TiltSeries_)

    def alignmentResidual(self,cut=-1):
        """
        calculate residual of a marker model given the marker coords
        """
        from pytom.reconstruction.tiltAlignmentFunctions import markerResidual
        if cut == -1 or cut +1 > len(self._Markers):
            start = 0
            end = len(self._Markers)
        else:
            start = cut
            end = start+1

        residual = markerResidual(cent=self.TiltSeries_._TiltAlignmentParas.cent,
                                  Markers_=self._Markers[start:end],
                                  cTilt=self._cTilt, sTilt=self._sTilt,
                                  transX=self._alignmentTransX, transY=self._alignmentTransY,
                                  rotInPlane=self._alignmentRotations,
                                  isoMag=self._alignmentMagnifications, dBeam=self._alignmentBeamTilt,
                                  dMagnFocus=None, dRotFocus=None, equationSet=False)
        return residual

    def alignmentScore(self, optimizableVariables):
        """
        compute alignment score for given parameters

        @param optimizableVariables: array of variables that are subject to optimization
        @type optimizableVariables: numpy array
        @return: alignment score
        """
        from pytom.reconstruction.tiltAlignmentFunctions import markerResidual
        self.setOptimizableVariables(self.TiltSeries_._TiltAlignmentParas, optimizableVariables)
        if self.TiltSeries_._TiltAlignmentParas.leastsq == True:
            score = markerResidual(self.TiltSeries_._TiltAlignmentParas.cent,
                                   Markers_=self._Markers,
                                   cTilt=self._cTilt, sTilt=self._sTilt,
                                   transX=self._alignmentTransX, transY=self._alignmentTransY,
                                   rotInPlane=self._alignmentRotations,
                                   isoMag=self._alignmentMagnifications, dBeam=self._alignmentBeamTilt,
                                   dMagnFocus=None, dRotFocus=None, equationSet=True)
        else:
            score = markerResidual(self.TiltSeries_._TiltAlignmentParas.cent,
                                   Markers_=self._Markers,
                                   cTilt=self._cTilt, sTilt=self._sTilt,
                                   transX=self._alignmentTransX, transY=self._alignmentTransY,
                                   rotInPlane=self._alignmentRotations,
                                   isoMag=self._alignmentMagnifications, dBeam=self._alignmentBeamTilt,
                                   dMagnFocus=None, dRotFocus=None, equationSet=False)
        return score

    def getOptimizableVariables(self, TiltAlignmentParameters_):
        """
        generate numpy array of optimizable variables from tilt series and given alignment parameters

        @param TiltAlignmentParameters_: parameters for tilt alignment
        @type TiltAlignmentParameters_: TiltAlignmentParameters
        @return: optimizableVariables (marker coords, translations, rotations, magnifications, beam tilt)
        @rtype: numpy array
        """
        ntilt = self._ntilt
        nmark = len(self._Markers)

        # coordinates of reference (imark) mark and reference projection (iproj) are fixed
        #nopti = (nmark-1)*3 + (ntilt-1)*2
        nopti = (nmark-1) * 3 + (ntilt) * 2

        #check that irefmark and ireftilt are set properly
        if not (TiltAlignmentParameters_.irefmark in range(nmark)):
            TiltAlignmentParameters_.irefmark = 1
            print("Warning: irefmark must be 1<= irefmark <=nmark")
            print("New irefmark: " + str(TiltAlignmentParameters_.irefmark))

        if not (TiltAlignmentParameters_.ireftilt in self._projIndices.astype(int)):
            TiltAlignmentParameters_.ireftilt = abs(self._tiltAngles).argmin() + 1
            print("Warning: ireftilt must be in range of projection indices")
            print("New ireftilt: " + str(TiltAlignmentParameters_.ireftilt))

        #variable rotation for projections 
        if TiltAlignmentParameters_.drot:
            nopti = nopti + ntilt
        else:
            nopti = nopti + 1

        #variable magnifications for projections 
        if TiltAlignmentParameters_.dmag:
            nopti = nopti + ntilt - 1

        # beam tilt
        if TiltAlignmentParameters_.dbeam:
            nopti = nopti + 1

        ## gradient on image rotation and magnification in projections
        #if TiltAlignmentParameters_.dGradRotMag:
        #    nopti = nopti + 2

        print(nopti)
        optimizableVariables = numpy.array(nopti * [0.])

        # marker 3D coords

        ivar = 0
        for (imark, Marker) in enumerate(self._Markers):
            # reference marker irefmark is fixed to standard value
            if ((imark ) != TiltAlignmentParameters_.irefmark):
                r = Marker.get_r()
                optimizableVariables[ivar] = r[0]
                optimizableVariables[ivar + 1] = r[1]
                optimizableVariables[ivar + 2] = r[2]
                ivar = ivar + 3

        # translations
        for itilt in range(0, ntilt):
            # translation in reference projection is zero
            #if self._projIndices[itilt] != TiltAlignmentParameters_.ireftilt:
            optimizableVariables[ivar] = self._alignmentTransX[itilt]
            optimizableVariables[ivar + 1] = self._alignmentTransY[itilt]
            ivar = ivar + 2

        # image rotations
        if TiltAlignmentParameters_.drot:
            for itilt in range(0, ntilt):
                optimizableVariables[ivar] = self._alignmentRotations[itilt]
                ivar = ivar + 1
        # all rotations are the same - take the first one
        else:
            optimizableVariables[ivar] = self._alignmentRotations[0]
            ivar = ivar + 1

        # magnification changes
        if TiltAlignmentParameters_.dmag:
            for itilt in range(0, ntilt):
                # magnification of reference projection is 1.
                if int(self._projIndices[itilt]) != TiltAlignmentParameters_.ireftilt:
                    optimizableVariables[ivar] = self._alignmentMagnifications[itilt]
                    ivar = ivar + 1
        # beam inclination
        if TiltAlignmentParameters_.dbeam:
            optimizableVariables[ivar] = self._alignmentBeamTilt
            ivar = ivar + 1

        # focus gradient (TODO)
        #if TiltAlignmentParameters_.dGradRotMag:
        #    optimizableVariables[ivar]   = self._alignmentMagnFoc
        #    optimizableVariables[ivar+1] = self._alignmentRotFoc

        return optimizableVariables

    def setOptimizableVariables(self, TiltAlignmentParameters_, optimizableVariables):
        """
        set values in tilt alignment according to specified optimizable variables and alignment parameters

        @param TiltAlignmentParameters_: parameters for tilt alignment
        @type TiltAlignmentParameters_: TiltAlignmentParameters
        @param optimizableVariables: optimizable alignment variables
        @type optimizableVariables: numpy array
        @return: alignment variables
        @rtype: 
        """
        ntilt = self._ntilt
        nmark = len(self._Markers)

        # coordinates of reference (imark) mark and reference projection (iproj) are fixed
        #FFnopti = (nmark-1)*3 + (ntilt-1)*2
        nopti = (nmark - 1) * 3 + (ntilt) * 2

        #variable rotation for projections 
        if TiltAlignmentParameters_.drot:
            nopti = nopti + ntilt
        else:
            nopti = nopti + 1

        #variable magnifications for projections 
        if TiltAlignmentParameters_.dmag:
            nopti = nopti + ntilt - 1

        # beam tilt
        if TiltAlignmentParameters_.dbeam:
            nopti = nopti + 1

        ## gradient on image rotation and magnification in projections
        #if TiltAlignmentParameters_.dGradRotMag:
        #    nopti = nopti + 2

        # check that number of variables is ok
        if len(optimizableVariables) != nopti:
            print("Length optimizableVariables: " + str(len(optimizableVariables)))
            print("N optmization: " + str(nopti))
            raise IndexError('length of optimizableVariables does not match TiltAlignmentParameters')

            # marker 3D coords
        ivar = 0
        for (imark, Marker) in enumerate(self._Markers):
            # reference marker irefmark is fixed to standard value
            if ((imark ) != TiltAlignmentParameters_.irefmark):
                r = numpy.array([optimizableVariables[ivar],
                                 optimizableVariables[ivar + 1], optimizableVariables[ivar + 2]])
                self._Markers[imark].set_r(r)
                ivar = ivar + 3

        # translations
        for itilt in range(0, ntilt):
            # translation in reference projection is zero
            #FFif (self._projIndices[itilt] != TiltAlignmentParameters_.ireftilt):
            self._alignmentTransX[itilt] = optimizableVariables[ivar]
            self._alignmentTransY[itilt] = optimizableVariables[ivar + 1]
            ivar = ivar + 2

        # image rotations
        if TiltAlignmentParameters_.drot:
            for itilt in range(0, ntilt):
                self._alignmentRotations[itilt] = optimizableVariables[ivar]
                ivar = ivar + 1
        # all rotations are the same - take the first one
        else:
            self._alignmentRotations[0] = optimizableVariables[ivar]
            ivar = ivar + 1

        # magnification changes
        if TiltAlignmentParameters_.dmag:
            for itilt in range(0, ntilt):
                # magnification of reference projection is 1.
                if (int(self._projIndices[itilt]) != TiltAlignmentParameters_.ireftilt):
                    self._alignmentMagnifications[itilt] = optimizableVariables[ivar]
                    ivar = ivar + 1

        # beam inclination
        if TiltAlignmentParameters_.dbeam:
            self._alignmentBeamTilt = optimizableVariables[ivar]
            ivar = ivar + 1

            # focus gradient (TODO)
            #if TiltAlignmentParameters_.dGradRotMag:
            #    optimizableVariables[ivar]   = self._alignmentMagnFoc
            #    optimizableVariables[ivar+1] = self._alignmentRotFoc

    def alignFromFiducials(self, mute=True):
        """
        align tilt series

        all necessary information is automatically obtained.
        Alignment seeded by:
          - marker 3D-coordinates in Markers
          - translations in projections of tilt series
          - rotations in projections of tilt series
          - magnifications in projections of tilt series
        @param mute: L{bool}
        @return: alignment score (=residual of markers)
        @rtype: float

        @author: FF
        """
        from math import sqrt
        import scipy.optimize
        from pytom.reconstruction.tiltAlignmentFunctions import markerResidual

        if self.TiltSeries_._TiltAlignmentParas.optimizer == 'fmin':
            optimizer = scipy.optimize.fmin
            if not mute:
                print("using scipy fmin optimizer")
        elif self.TiltSeries_._TiltAlignmentParas.optimizer == 'fmin_slsqp':
            optimizer = scipy.optimize.fmin_slsqp
            if not mute:
                print("using scipy fmin_slsqp (Sequential Least SQuares Programming) optimizer")
        elif self.TiltSeries_._TiltAlignmentParas.optimizer == 'fmin_cg':
            optimizer = scipy.optimize.fmin_cg
            if not mute:
                print("using scipy fmin_cg (conjugate gradients) optimizer")
        elif self.TiltSeries_._TiltAlignmentParas.optimizer == 'leastsq':
            optimizer = scipy.optimize.leastsq
            if not mute:
                print("using scipy leastsq optimizer - optimize matrix instead of scalar function")
            self.TiltSeries_._TiltAlignmentParas.leastsq = True
        elif self.TiltSeries_._TiltAlignmentParas.optimizer == 'fmin_powell':
            optimizer = scipy.optimize.fmin_powell
            if not mute:
                print("using scipy fmin_powell optimizer")
        else:
            if not mute:
                print(("optimizer " + str(self.TiltSeries_._TiltAlignmentParas.optimizer) +
                      " not known"))
        # first update alignment from projections
        self.getMarkersFromTiltSeries(self.TiltSeries_)
        self.getTranslationsFromTiltSeries(self.TiltSeries_)
        self.getRotationsFromTiltSeries(self.TiltSeries_)
        self.getMagnificationsFromTiltSeries(self.TiltSeries_)
        optimizableVariables0 = self.getOptimizableVariables(self.TiltSeries_._TiltAlignmentParas)

        # alignment score before optimization
        score = markerResidual(self.TiltSeries_._TiltAlignmentParas.cent,
            Markers_=self._Markers,
            cTilt=self._cTilt, sTilt=self._sTilt,
            transX=self._alignmentTransX, transY=self._alignmentTransY,
            rotInPlane=self._alignmentRotations,
            isoMag=self._alignmentMagnifications, dBeam=self._alignmentBeamTilt,
            dMagnFocus=None, dRotFocus=None, equationSet=False)

        if not mute:
            print(( "Alignment score before optimization (square root of residual): "
                   + str(sqrt(score)) ))

        # optimize scoring function
        if ((self.TiltSeries_._TiltAlignmentParas.optimizer == 'fmin') or
                (self.TiltSeries_._TiltAlignmentParas.optimizer == 'fmin_powell')):
            optimizableVariables = optimizer(self.alignmentScore, optimizableVariables0,
                                             xtol=0.000001, ftol=0.000001,
                                             maxiter=self.TiltSeries_._TiltAlignmentParas.maxIter, maxfun=None)
        elif self.TiltSeries_._TiltAlignmentParas.optimizer == 'fmin_cg':
            optimizableVariables = optimizer(self.alignmentScore, optimizableVariables0,
                                             gtol=0.0000001,
                                             maxiter=self.TiltSeries_._TiltAlignmentParas.maxIter)
        elif self.TiltSeries_._TiltAlignmentParas.optimizer == 'fmin_slsqp':
            optimizableVariables = optimizer(self.alignmentScore, optimizableVariables0,
                                             iter=self.TiltSeries_._TiltAlignmentParas.maxIter, acc=1e-08)
        elif self.TiltSeries_._TiltAlignmentParas.optimizer == 'leastsq':
            optimizableVariables, success = optimizer(self.alignmentScore, optimizableVariables0,
                                                      maxfev=self.TiltSeries_._TiltAlignmentParas.maxIter, epsfcn=0.0,
                                                      factor=10)

        score = markerResidual(self.TiltSeries_._TiltAlignmentParas.cent,
            Markers_=self._Markers,
            cTilt=self._cTilt, sTilt=self._sTilt,
            transX=self._alignmentTransX, transY=self._alignmentTransY,
            rotInPlane=self._alignmentRotations,
            isoMag=self._alignmentMagnifications, dBeam=self._alignmentBeamTilt,
            dMagnFocus=None, dRotFocus=None, equationSet=False)
        if not mute:
            print("Alignment Score after optimization: " + str(sqrt(score)))

        self.setOptimizableVariables(self.TiltSeries_._TiltAlignmentParas,
                                     optimizableVariables)
        # finally set values in tilt series
        self.setMarkersInTiltSeries(self.TiltSeries_)
        self.setTranslationsInTiltSeries(self.TiltSeries_)
        self.setRotationsInTiltSeries(self.TiltSeries_)
        self.setMagnificationsInTiltSeries(self.TiltSeries_)
        return sqrt(score)

    def alignmentResidualGradient(self, TiltSeries_):
        """
        """

    def optimizeAlignment(self, TiltSeries_):
        """
        """

    def computeCoarseAlignment(self, TiltSeries_, mute=True, outfile=''):
        """
        compute alignment analytically (constant mag. and tilt axis)

        @param TiltSeries_: Tiltseries
        @type TiltSeries_: L{Tiltseries}
        @param mute: turn output silent
        @type mute: L{bool}

        @author: FF
        """
        #print('ref index: ', numpy.argwhere( self._projIndices.astype(int) == TiltSeries_._TiltAlignmentParas.ireftilt)[0][0], TiltSeries_._TiltAlignmentParas.ireftilt )
        (psiindeg, shiftX, shiftY, x, y, z, distLine, diffX, diffY,
         shiftVarX, shiftVarY) = alignmentFixMagRot(
            Markers_=self._Markers, cTilt=self._cTilt, sTilt=self._sTilt,
            ireftilt=numpy.argwhere( self._projIndices.astype(int) == TiltSeries_._TiltAlignmentParas.ireftilt)[0][0],
            irefmark=TiltSeries_._TiltAlignmentParas.irefmark,
            r=TiltSeries_._TiltAlignmentParas.r, imdim=TiltSeries_._imdim,
            handflip=TiltSeries_._TiltAlignmentParas.handflip, mute=mute, writeResults=outfile)
        if not mute:
            print(("Tilt Axis: %.2f" % psiindeg))
        # copy parameters to TiltSeries
        self._alignmentRotations = numpy.array(self._ntilt * [psiindeg])
        self.setRotationsInTiltSeries(TiltSeries_)
        self._alignmentTransX = shiftX
        self._alignmentTransY = shiftY
        self.set_TranslationsInTiltSeries(TiltSeries_)
        self.Psi = psiindeg

        for (imark, Marker) in enumerate(self._Markers):
            Marker.set_r(numpy.array([x[imark], y[imark], z[imark]]))
            
    def set_TranslationsInTiltSeries(self, TiltSeries_):
        """
        set translations in TiltSeries

        @author: FF
        """
        for (kk, Proj) in enumerate(TiltSeries_._ProjectionList):
            Proj._alignmentTransX = self._alignmentTransX[kk]
            Proj._alignmentTransY = self._alignmentTransY[kk]

    def get_RotationsFromTiltSeries(self, TiltSeries_):
        """
        get rotations from TiltSeries

        @author: FF
        """
        # initialize alignment 
        self.rotInPlane = len(TiltSeries_.Projections) * [0.]
        kk = 0
        for Proj in TiltSeries_.Projections:
            self.rotInPlane[kk] = Proj.rotInPlane
            kk = kk + 1
        return self.rotInPlane

    def set_RotationsInTiltSeries(self, TiltSeries_):
        """
        set Rotations in TiltSeries

        @author: FF
        """
        kk = 0
        for Proj in TiltSeries_.Projections:
            Proj.rotInPlane = self.rotInPlane[kk]
            kk = kk + 1

    def refineMarkerPositions(self, TiltSeries_, dimBox=32,
                              finealigfile=None, verbose=False):
        """
        refine coordinates of markers

        @param TiltSeries_: Tilt series
        @type TiltSeries_: L{TiltSeries}
        @param dimBox: dimension of box for markers
        @type dimBox: L{int}
        @param finealigfile: output file with refined coordinates
        @type finealigfile: L{str}
        @param verbose: verbose for debugging
        @type verbose: L{bool}
    
        @author: FF
        """
        from pytom.image2D.imageStructures import ImageStack
        from pytom.basic.functions import initSphere, taper_edges

        # prepare mask
        #mask = initSphere(sizeX=dimBox, sizeY=dimBox, sizeZ=1, radius=dimBox/5.,
        #          smooth=dimBox/5., maxradius=0, cent=None)
        mask = None

        # coarse alignment to get errors
        print(">>> Alignment error before marker position fitting:")
        self.computeCoarseAlignment(TiltSeries_, mute=True)

        for (imark, marker) in enumerate(TiltSeries_._Markers):
            markerImageStack = ImageStack(verbose=verbose)
            for itilt in range(0, len(marker._projIndices)):
                x = round(marker.get_xProj(itilt))
                y = round(marker.get_yProj(itilt))
                if (x > -1) and ( y > -1):
                    proj = TiltSeries_._ProjectionLists[itilt]
                    filename = proj.getFilename()
                    #copy data to ImageStack
                    markerImageStack.addImageFromFile(filename=filename,
                                                      boxCoords=[x - dimBox / 2 - 2, y - dimBox / 2 - 2, 0],
                                                      dims=[dimBox, dimBox, 1], shiftX=0, shiftY=0, rotation=0,
                                                      appliedShiftX=0, appliedShiftY=0, index=itilt)
                else:
                    proj = TiltSeries_._ProjectionList[itilt]
                    print("marker not clicked in " + proj.getFilename())
            #markerImageStack.normalize( normtype="StdMeanInMask", mask=mask)
            markerImageStack.bandpass(lowfreq=dimBox / 20., hifreq=7. * dimBox / 32.,
                                      smooth=3. * dimBox / 32., bpf=None)
            markerImageStack.normalize(normtype="StdMean")
            # smoothen edges
            markerImageStack.taper_edges(width=dimBox / 8.)
            markerImageStack.exMaxAlign(niter=10, mask=mask)
            if verbose:
                markerImageStack.writeWorkingCopies(filename="Marker" + str(imark) + '.em')
                markerImageStack.writeAverage(filename='av_' + str(imark) + ".em")

            # set new coordinates
            irun = 0
            for itilt in range(0, len(marker._projIndices)):
                x = round(marker.get_xProj(itilt))
                y = round(marker.get_yProj(itilt))
                if ( x > -1 ) and ( y > -1 ):
                    marker.set_xProj(itilt, x + markerImageStack.images[irun].shiftX)
                    marker.set_yProj(itilt, y + markerImageStack.images[irun].shiftY)
                    irun = irun + 1
        print(">>> Alignment error after marker position fitting:")
        self.computeCoarseAlignment(TiltSeries_, mute=True)
        #write new coordinates
        if finealigfile:
            TiltSeries_.writeMarkerFile(markerFileName=finealigfile)


class Marker:
    """
    Marker in tilt series
    """

    def __init__(self, projIndices):
        """
        initialize Marker class

        @param projIndices: projection Indices
        @type projIndices: list of integers

        @author: FF
        """
        self._nproj = len(projIndices)
        self._projIndices = projIndices
        # coords in projections
        self.xProj = numpy.array(self._nproj * [-1.])
        self.yProj = numpy.array(self._nproj * [-1.])
        # 3D coordinates of marker
        self._r = numpy.array([0., 0., 0.])

    def get_numberOfProjections(self):
        """
        get number of projections for marker (also update nproj in structure)

        @return: number of projections
        @rtype: L{int}
        """
        self._nproj = len(self._projIndices)
        return self._nproj

    def info(self):
        """
        """
        r = self.get_r()
        tline = ("3D model: x=%7.1f, " % r[0])
        tline = tline + ("y=%7.1f, " % r[1])
        tline = tline + ("z=%7.1f\n" % r[2])
        tline = tline + "============================================\n"
        tline = tline + "Coordinates in projections: \n"
        for iproj in range(0, self._nproj):
            tline = tline + ("%3d: " % (iproj + 1))
            tline = tline + ("%7.1f, " % self.get_xProj(iproj))
            tline = tline + ("%7.1f\n" % self.get_yProj(iproj))
        print(tline)

    def set_xProj(self, iproj, xProj):
        """
        set x-value of marker in projection iproj

        @param iproj: index of projection
        @type iproj: L{int}
        @param xProj: x-coordinate
        @type xProj: L{float}

        @author: FF
        """
        self.xProj[iproj] = xProj

    def get_xProj(self, iproj):
        """
        get x-value of marker in projection iproj

        @author: FF
        """
        return self.xProj[iproj]

    def set_yProj(self, iproj, yProj):
        """
        set y-value of marker in projection iproj

        @param iproj: index of projection
        @type iproj: L{int}
        @param yProj: y-coordinate
        @type yProj: L{float}

        @author: FF
        """
        self.yProj[iproj] = yProj

    def get_yProj(self, iproj):
        """
        get y-value of marker in projection iproj

        @author: FF
        """
        return self.yProj[iproj]

    def set_xProjs(self, xProjs):
        """
        set x-values of marker in all projections of tilt series

        @param xProjs: x-coordinates of markers in projection
        @type xProjs: numpy array

        @author: FF
        """
        for (ii, xProj) in enumerate(xProjs):
            self.xProj[ii] = xProj

    def set_yProjs(self, yProjs):
        """
        set y-values of marker in all projections of tilt series

        @param yProjs: y-coordinates of markers in projection
        @type yProjs: numpy array

        @author: FF
        """
        for (ii, yProj) in enumerate(yProjs):
            self.yProj[ii] = yProj

    def set_r(self, r):
        """
        set 3d coordinates r of marker

        @param r: marker
        @type r: 3d numpy array

        @author: FF
        """
        self._r = r

    def get_r(self):
        """
        get 3d coordinates r of marker

        @return: coordinates of marker
        @rtype: 3d numpy array

        @author: FF
        """
        return self._r


class TiltAlignmentParameters:
    """
    class to store alignment of projections in tilt series

    @author: FF
    """

    def __init__(self, dmag=True, drot=True, dbeam=True,
                 finealig=True, finealigfile='xxx.dat',
                 grad=False,
                 irefmark=1, ireftilt=0, r=[0., 0., 0.],
                 cent=[2049, 2049],
                 handflip=False,
                 optimizer='fmin', maxIter=2000, leastsq=False):
        """initialize TiltAlignmentParameters class

        @param dmag: different magnification of tilts (default: true)
        @type dmag: bool
        @param finealig: fine-center gold beads, true/false(default)
        @type finealig: bool
        @param finealigfile: filename of re-centered gold beads (only used when 'finealig' is true)
        @type finealigfile: string
        @param drot: different in-plane rotations of projections (default: True)
        @type drot: bool
        @param dbeam: beam inclination (default: False)
        @type dbeam: bool
        @param grad: compute gradients of residual for optimization (default: False)
        @type grad: bool
        @param irefmark: index of reference marker (default:1)
        @type irefmark: int
        @param ireftilt: index of reference tilt angle (default: 0)
        @type ireftilt: int
        @param cent: lateral center (default: [2049,2049])
        @type cent: 2-dim array
        @param r: coordinates assigned to reference marker (default: [x(ireftilt),y(ireftilt),1025] ); if [0,0,0] default values are taken 
        @type r: array
        @param handflip: add pi to intially determined tilt axis
        @type handflip: bool
        @param maxIter:  Maximum number of iteration for optimization (default: 2000)
        @type maxIter: int
        @param optimizer: optimzation algorithm 
        @type optimizer: str
        @param leastsq: least-squared optimization
        @type leastsq: bool

        @author: Friedrich Foerster
        """
        self.dmag = dmag
        self.drot = drot
        self.dbeam = dbeam
        self.finealig = finealig
        self.finealigfile = finealigfile
        self.grad = grad
        self.irefmark = irefmark
        self.ireftilt = ireftilt
        self.r = r
        self.cent = cent
        self.handflip = handflip
        self.optimizer = optimizer
        self.maxIter = maxIter
        self.leastsq = leastsq

