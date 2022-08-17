'''
functions underlying 3D reconstruction
Created on Dec 7, 2010

@lastchange: Nov 2015, FF
@author: hrabe, ff
'''

def adjustPickPositionInParticleList(particleList,offsetZ,binning,centerX,centerY,centerZ):
    """
    adjustPickPositionInParticleList:
    @param offsetZ: Z direction distance where the tomogram used for template matching was cut out
    @param binning: Binning factor of tomogram and thus of determined pick positions. Was the tomo 1x binned -> binning = 2^1 = 2 and so forth
    @param centerX: The center coordinate in the original sized tomogram
    @param centerY: The center coordinate in the original sized tomogram
    @param centerZ: The center coordinate in the original sized tomogram
    @return: Particle list with adjusted coordinates, ready for reconstruction   
    """
    
    for particle in particleList:
        
        pickPosition = particle.getPickPosition()
        x = pickPosition.getX()
        y = pickPosition.getY()
        z = pickPosition.getZ()
        
        z = z + offsetZ
        
        x = x * binning - centerX
        y = y * binning - centerY
        z = z * binning - centerZ
        
        pickPosition.setX(x)
        pickPosition.setY(y)
        pickPosition.setZ(z)
        
        particle.setPickPosition(pickPosition)
    
    return particleList
    
        
def positionInProjection(location3D, tiltAngle, tiltAxis='Y'):
    """
    positionInProjection:
    @param location3D: [x,y,z] vector of particle positions in tomogram
    @param tiltAngle: tilt angle
    @type tiltAngle: float
    @param tiltAxis: Set tilt axis of projections. 'Y' is default.
    @type tiltAxis: str
    @return: 2D Position vector [x,y]
    @author: Thomas Hrabe  
    """
    if tiltAxis == 'X':
        from pytom.tools.maths import XRotationMatrix
        rotationMatrix = XRotationMatrix(tiltAngle) 
    elif tiltAxis == 'Y':
        from pytom.tools.maths import YRotationMatrix
        rotationMatrix = YRotationMatrix(tiltAngle)
    elif tiltAxis == 'Z':
        from pytom.tools.maths import ZRotationMatrix
        rotationMatrix = ZRotationMatrix(tiltAngle)

    position2D = rotationMatrix * location3D

    return [position2D[0],position2D[1]]


def positionsInProjections(location3D,tiltStart,tiltIncrement,tiltEnd,tiltAxis='Y'):
    """
    positionsInProjections: Returns 2D positions in projections of a particle according to its 3D location in tomogram  
    @param location3D: 
    @param tiltStart:
    @param tiltIncrement:
    @param tiltEnd:
    @param tiltAxis: Set tilt axis of projections. Y is default.
    @return: List of positions starting with position for projection at tiltStart
    @author: Thomas Hrabe
    """
    positionList = []
    
    for angle in xrange(tiltStart,tiltEnd + tiltIncrement,tiltIncrement): 
        positionList.append(positionInProjection(location3D,angle,tiltAxis))
        
    return positionList


def alignWeightReconstruct(tiltSeriesName, markerFileName, lastProj, tltfile=None, prexgfile=None, preBin=None,
                           volumeName=None, volumeFileType='em', alignResultFile='',
                           voldims=None, recCent=[0,0,0], tiltSeriesFormat='st', firstProj=0, irefmark=1, ireftilt=1,
                           handflip=False, alignedTiltSeriesName='align/myTilt', weightingType=-1,
                           lowpassFilter=1., projBinning=1, outMarkerFileName=None, verbose=False, outfile='',
                           write_images=True, shift_markers=True, logfile_residual='', refMarkTomo='', profile=False,
                           gpuID=-1,specimen_angle=0):
    """
    @param tiltSeriesName: Name of tilt series (set of image files in .em or .mrc format) or stack file (ending '.st').\
    Note: the actual file ending should NOT be provided.
    @type tiltSeriesName: str
    @param markerFileName: name of EM markerfile or IMOD wimp File containing marker coordinates
    @type markerFileName: str
    @param lastProj: index of last projection
    @param tltfile: ascii file containing the tilt angles
    @type tltfile: str
    @type lastProj: int
    @param prexgfile: file containing pre-shifts (IMOD way of doing things)
    @type prexgfile: str
    @param preBin: pre-binning in IMOD used for marker determination
    @type preBin: int
    @param volumeName: Filename of volume
    @type volumeName: str
    @param volumeFileType: type of output volume - em or mrc
    @type volumeFileType: str
    @param voldims: dimensions of reconstruction volume, e.g. [512,512,128] - if None chosen only aligned projs are \
    written
    @type voldims: list
    @param recCent: offset from center of volume - for example choose z!=0 to shift in z (post-binning coordinates) - \
    default (0,0,0)
    @type recCent: list
    @param tiltSeriesFormat: file format of tilt series: 'em', 'mrc' or 'st'
    @type tiltSeriesFormat: str
    @param firstProj: index of first projection - default 1
    @type firstProj: int
    @param irefmark: index of reference marker for alignment
    @type irefmark: int
    @param ireftilt: index of reference tilt
    @type ireftilt: int
    @param handflip: add 180 deg to default tilt axis resulting in change of handedness?
    @type handflip: bool
    @param alignedTiltSeriesName: name of ailgned tilt series
    @type alignedTiltSeriesName: str
    @param weightingType: type of weighting - -1 is analytical 'ramp' weighting and 0 is no weighting
    @type weightingType: int
    @param lowpassFilter: lowpass filter in Nyquist
    @type lowpassFilter: float
    @param projBinning: binning of projection (factor rather than power of two!)
    @type projBinning: int or float
    @param outMarkerFileName: filename of output marker file
    @type outMarkerFileName: str
    @param verbose: verbose?
    @type verbose: bool

    @author: FF
    """
    from pytom.reconstruction.TiltAlignmentStructures import TiltSeries, TiltAlignment, TiltAlignmentParameters
    from pytom.gui.guiFunctions import savestar, headerMarkerResults, fmtMR
    import numpy

    if profile:
        import time
        s = time.time()

    if not alignResultFile:  # we need to find the alignment, optionally write out aligned images and reconstruction
        if verbose:
            print("Function alignWeightReconstruct started")
            mute = False
        else:
            mute = True

        if not shift_markers:
            irefmarkInit = irefmark
            irefmark = refMarkTomo if refMarkTomo else 0

        tiltParas = TiltAlignmentParameters(dmag=True, drot=True, dbeam=False, finealig=True,
                                            finealigfile='xxx.txt', grad=False,
                                            irefmark=irefmark, ireftilt=ireftilt, r=None, cent=[2049, 2049],
                                            handflip=handflip, optimizer='leastsq', maxIter=1000)

        markerFileType = markerFileName.split('.')[-1]
        # align with wimpfile
        if not preBin:
            preBin=1
        if markerFileType == 'wimp':
            if verbose:
                print(" WIMP file used for alignment")
            tiltSeries = TiltSeries(tiltSeriesName=tiltSeriesName, TiltAlignmentParas=tiltParas,
                                    alignedTiltSeriesName=alignedTiltSeriesName,
                                    markerFileName=None, firstProj=firstProj, lastProj=lastProj,
                                    tiltSeriesFormat=tiltSeriesFormat)
            tiltSeries.readIMODwimp(markerFileName=markerFileName, prexgfile=prexgfile, tltfile=tltfile, preBin=preBin,
                                    verbose=False)
        else:
            if verbose:
                format = markerFileName.split('.')[-1].upper()
                print(f"{format} markerfile file used for alignment")
            tiltSeries = TiltSeries(tiltSeriesName=tiltSeriesName, TiltAlignmentParas=tiltParas,
                                    alignedTiltSeriesName=alignedTiltSeriesName,
                                    markerFileName=markerFileName, firstProj=firstProj, lastProj=lastProj,
                                    tiltSeriesFormat=tiltSeriesFormat)
            if tltfile:
                tiltSeries.getTiltAnglesFromIMODfile(tltfile=tltfile)
        tiltAlignment = TiltAlignment(TiltSeries_=tiltSeries)
        if outMarkerFileName:
            tiltSeries.writeMarkerFile(markerFileName=outMarkerFileName)
            tiltSeries._markerFileName = outMarkerFileName
        tiltAlignment.resetAlignmentCenter()  # overrule cent in Paras



        tiltAlignment.computeCoarseAlignment(tiltSeries, mute=mute, optimizeShift=shift_markers)
        ireftilt = int(
            numpy.argwhere(tiltAlignment._projIndices.astype(int) == tiltSeries._TiltAlignmentParas.ireftilt)[0][0])


        tiltAlignment.alignFromFiducials(mute=mute, shift_markers=True, logfile_residual=logfile_residual)




        if not shift_markers:
            tiltSeries._TiltAlignmentParas.irefmark = irefmarkInit
            tiltParas.irefmark = irefmarkInit

            for (imark, Marker) in enumerate(tiltAlignment._Markers):
                r = Marker.get_r()
                # Marker.set_r(numpy.array([r[0] - tiltSeries._ProjectionList[ireftilt]._alignmentTransX,
                #                           r[1] - tiltSeries._ProjectionList[ireftilt]._alignmentTransY, r[2]]))

            tiltAlignment.alignFromFiducials(mute=mute, shift_markers=False, logfile_residual=logfile_residual)

        values = []

        if outfile:  # this is false with the default name ''
            # Retrieve the index of reference image
            if verbose:
                print(f'Writing reference marker file: {outfile}')

            ireftilt = int(numpy.argwhere( tiltAlignment._projIndices.astype(int) == tiltSeries._TiltAlignmentParas.ireftilt)[0][0])

            # for n, i in enumerate(tiltAlignment._projIndices.astype(int)):
            #
            #
            #     print("{:3.0f} {:6.0f} {:6.0f} {:6.0f} {:6.0f}".format(tiltAlignment._tiltAngles[n],
            #           tiltAlignment._Markers[tiltSeries._TiltAlignmentParas.irefmark].xProj[n],
            #           tiltAlignment._Markers[tiltSeries._TiltAlignmentParas.irefmark].yProj[n],
            #           tiltSeries._ProjectionList[int(n)]._alignmentTransX,
            #           tiltSeries._ProjectionList[int(n)]._alignmentTransY))

            # Create center point (3D)
            cent = tiltAlignment.TiltSeries_._TiltAlignmentParas.cent
            cent.append(float(tiltSeries._imdim // 2 + 1))

            # Create shift
            shift = [tiltSeries._ProjectionList[ireftilt]._alignmentTransX,
                     tiltSeries._ProjectionList[ireftilt]._alignmentTransY, 0]

            # Retrieve the two relevant angles
            inPlaneAng = tiltAlignment._alignmentRotations[ireftilt]
            tiltAng = tiltAlignment._tiltAngles[ireftilt]

            # Retrieve the original pick position
            xx = tiltAlignment._Markers[tiltSeries._TiltAlignmentParas.irefmark].xProj[ireftilt]
            yy = tiltAlignment._Markers[tiltSeries._TiltAlignmentParas.irefmark].yProj[ireftilt]
            dimx, dimy = float(tiltSeries._imdimX // 2 + 1), float(tiltSeries._imdimY // 2 + 1)

            if abs(90 - (inPlaneAng % 180)) < 45:
                dx, dy = dimx, dimy
            else:
                dx, dy = dimy, dimx

            dz = min(dx,dy)

            # Create ref marker point
            ref = tiltAlignment._Markers[tiltSeries._TiltAlignmentParas.irefmark].get_r()
            ref = rotate_vector3D(ref, [0,0,0], shift, inPlaneAng, tiltAng)

            for (imark, Marker) in enumerate(tiltAlignment._Markers):
                # reference marker irefmark is fixed to standard value
                r = rotate_vector3D(Marker.get_r(), [0,0,0], shift, inPlaneAng, tiltAng)
                values.append([imark, r[0]-ref[0], r[1]-ref[1], r[2]-ref[2], r[0]+dx, r[1]+dy, r[2]+dz])

            savestar(outfile, numpy.array(values), header=headerMarkerResults, fmt=fmtMR)

        # creating dir for aligned tilt series if default filename
        if alignedTiltSeriesName == 'align/myTilt':
            from os import mkdir
            try:
                mkdir('align')
            except OSError:
                print(" dir 'align' already exists - writing aligned files into existing dir")

        tiltSeries.write_aligned_projs(weighting=weightingType, lowpassFilter=lowpassFilter, binning=projBinning,
                                       verbose=verbose, write_images=write_images)
        if voldims:
            # overrule tiltSeriesFormat - aligned tiltseries is always a series of em files
            #tiltSeries._tiltSeriesFormat = 'em'
            vol_bp = tiltSeries.reconstructVolume(dims=voldims, reconstructionPosition=recCent, binning=1,
                                                  specimen_angle=specimen_angle, gpu=gpuID)

            try:
                vol_bp.write(volumeName, volumeFileType)
            except:
                from pytom.agnostic.io import write
                write(volumeName, vol_bp) #vol_bp.write(volumeName, volumeFileType)

    else:  # an alignment file is supplied: align the images and reconstruct the volume

        tiltSeries = TiltSeries(tiltSeriesName=tiltSeriesName, tiltSeriesFormat=tiltSeriesFormat)
        vol_bp = tiltSeries.reconstructVolume(dims=voldims, reconstructionPosition=recCent, binning=projBinning,
                                              alignResultFile=alignResultFile, applyWeighting=weightingType, gpu=gpuID,
                                              specimen_angle=specimen_angle, read_only=False)

        try:
            vol_bp.write(volumeName, volumeFileType)
        except:
            from pytom.agnostic.io import write
            write(volumeName, vol_bp)

    if profile:
        print(f'alignWeightReconstruct took {(time.time()-s):.1f} sec')


# TODO why is this function here?????
def rotate_vector3D(point3D, centerPoint3D, shift3D, inPlaneRotAngle, tiltAngle):
    from pytom.tools.maths import rotate_vector2d
    from numpy import sin, cos, pi

    newPoint = []

    # Recenter and shift
    for p,c,s in zip(point3D, centerPoint3D, shift3D):
        newPoint.append(p-c-s)

    # In-plane rotation of x and y coordinate
    inPlaneRotAngle = -inPlaneRotAngle-90
    cpsi = cos((inPlaneRotAngle%360) / 180. * pi)
    spsi = sin((inPlaneRotAngle%360) / 180. * pi)
    xmod, ymod = rotate_vector2d(newPoint[:2], cpsi, spsi)
    zmod = newPoint[2]

    # In-plane rotation of x and z coordinate over tiltAngle
    sTilt = sin(tiltAngle / 180. * pi)
    cTilt = cos(tiltAngle / 180. * pi)
    xmod,zmod = rotate_vector2d([xmod,zmod], cTilt, sTilt)

    rotatedNewPoint = [xmod, ymod, zmod]

    # Add center
    transformedPoint = []
    for p,c in zip(rotatedNewPoint, centerPoint3D):
        transformedPoint.append(p+c)

    return transformedPoint


# TODO this one could directly be placed in the projection list structure as
# TODO toProjectionStackFromAlignmentResultsFileGPU

def alignImagesUsingAlignmentResultFile(alignmentResultsFile, weighting=None, lowpassFilter=0.9, binning=1,
                                        circleFilter=True, angle_specimen=0):
    from pytom.gui.guiFunctions import fmtAR, headerAlignmentResults, datatype, datatypeAR, loadstar
    from pytom.reconstruction.reconstructionStructures import Projection, ProjectionList
    from pytom.agnostic.io import read, write, read_size
    from pytom.agnostic.tools import taper_edges, create_circle, paste_in_center as pasteCenter
    from pytom.agnostic.filter import circle_filter, ramp_filter, exact_filter, ellipse_filter
    import pytom.voltools as vt
    from pytom.gpu.initialize import xp, device
    from pytom.agnostic.transform import resize

    print("Create aligned images from alignResults.txt")

    alignmentResults = loadstar(alignmentResultsFile, dtype=datatypeAR)
    imageList = alignmentResults['FileName']
    tilt_angles = alignmentResults['TiltAngle']

    imdimX = read_size(imageList[0], 'x')
    imdimY = read_size(imageList[0], 'y')

    if binning > 1:
        imdimX = int(float(imdimX) / float(binning) + .5)
        imdimY = int(float(imdimY) / float(binning) + .5)

    imdim = max(imdimY, imdimX)
    sliceWidth = imdim

    if (weighting != None) and (float(weighting) < -0.001):
        cfreq = abs(tilt_angles[1:]-tilt_angles[:-1])
        cfreq = float(1/xp.sin(cfreq.min()*xp.pi/180))//1
        weightSlice = xp.fft.fftshift(ramp_filter(imdim, imdim, cfreq, len(tilt_angles)))

    if circleFilter:
        circleFilterRadiusX = imdim // 2
        circleFilterRadiusY = imdim // 2
        circleSlice = xp.fft.fftshift(ellipse_filter(imdim, imdim, circleFilterRadiusX, circleFilterRadiusY))
    else:
        circleSlice = xp.ones((imdim, imdim))

    # design lowpass filter
    if lowpassFilter:
        if lowpassFilter > 1.:
            lowpassFilter = 1.
            print("Warning: lowpassFilter > 1 - set to 1 (=Nyquist)")

    projectionList = ProjectionList()
    for n, image in enumerate(imageList):
        atx = alignmentResults['AlignmentTransX'][n]
        aty = alignmentResults['AlignmentTransY'][n]
        rot = alignmentResults['InPlaneRotation'][n]
        mag = 1 / (alignmentResults['Magnification'][n])
        projection = Projection(imageList[n], tiltAngle=tilt_angles[n], alignmentTransX=atx, alignmentTransY=aty,
                                alignmentRotation=rot, alignmentMagnification=mag)
        projectionList.append(projection)

    stack = xp.zeros((imdim, imdim, len(imageList)), dtype=xp.float32)
    phiStack = xp.zeros((len(imageList)), dtype=xp.float32)
    thetaStack = xp.zeros((len(imageList)), dtype=xp.float32)
    offsetStack = xp.zeros((len(imageList), 2), dtype=xp.float32)

    for (ii, projection) in enumerate(projectionList):
        #print(f'Align {projection._filename}')
        image = read(str(projection._filename)).squeeze()

        if binning > 1:
            image = resize(image, 1 / binning)

        tiltAngle = projection._tiltAngle

        # 1 -- normalize to contrast - subtract mean and norm to mean
        immean = image.mean()
        image = (image - immean) / immean

        # 2 -- smoothen borders to prevent high contrast oscillations
        if ii == 0:
            image, taper_mask = taper_edges(image, imdim // 30)

        else:
            image *= taper_mask #taper_edges(image,imdim//30,taper_mask)[0]

        # 3 -- square if needed
        if imdimY != imdimX:
            newImage = xp.zeros((imdim, imdim), dtype=xp.float32)
            image = pasteCenter(image, newImage)

        # 4 -- transform projection according to tilt alignment
        transX = projection._alignmentTransX / binning
        transY = projection._alignmentTransY / binning
        rot = float(projection._alignmentRotation)
        mag = float(projection._alignmentMagnification)

        inputImage = xp.expand_dims(image, 2)

        if ii ==0: outputImage = xp.zeros_like(inputImage, dtype=xp.float32)
        else: outputImage *= 0

        # TODO provide operation order parameter to voltools
        vt.transform(inputImage.astype(xp.float32), rotation=[0, 0, rot], rotation_order='rxyz', output=outputImage,
                     device=device, translation=[transX, transY, 0], scale=[mag, mag, 1], interpolation='filt_bspline')
        del inputImage
        image = outputImage.squeeze()

        # 5 -- lowpass filter (optional)
        if lowpassFilter:
            from pytom.agnostic.filter import bandpass_circle
            image = bandpass_circle(image, high=lowpassFilter * imdim // 2, sigma=lowpassFilter /5. * imdim)

        # 6 -- smoothen once more to avoid edges
        if imdimY == imdimX: image *= taper_mask
        else:
           image = taper_edges(image, imdim // 30)[0]

        # 7 -- weighting
        if (weighting != None) and (weighting < 0):
            # image = (ifft(complexRealMult(fft(image), w_func)) / (image.sizeX() * image.sizeY() * image.sizeZ()))
            image = xp.fft.ifftn(xp.fft.fftn(image) * weightSlice * circleSlice).real
        elif (weighting != None) and (weighting > 0):
            weightSlice = xp.fft.fftshift(exact_filter(tilt_angles, tiltAngle, imdim, imdim, sliceWidth))
            image = xp.fft.ifftn(xp.fft.fftn(image) * weightSlice * circleSlice).real

        thetaStack[ii] = float(round(projection.getTiltAngle() - angle_specimen))
        offsetStack[ii,:] = xp.array([int(round(projection.getOffsetX())), int(round(projection.getOffsetY()))])
        stack[:, :, ii] = image

    return [stack, phiStack, thetaStack, offsetStack]
