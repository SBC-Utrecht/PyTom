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

    if not alignResultFile:
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




        if outfile:
            # Retrieve the index of reference image
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

    else:

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


def alignImageUsingAlignmentResultFile(alignmentResultsFile, indexImage, weighting=None, lowpassFilter=0.9, binning=1,
                                       circleFilter=False):

    from pytom.gui.guiFunctions import fmtAR, headerAlignmentResults, datatype, datatypeAR, loadstar
    from pytom.reconstruction.reconstructionStructures import Projection, ProjectionList
    from pytom.agnostic.io import read, write, read_size
    from pytom.agnostic.tools import taper_edges, create_circle, paste_in_center as pasteCenter
    from pytom.agnostic.filter import circle_filter, ramp_filter, exact_filter, ellipse_filter
    import pytom.voltools as vt
    from pytom.gpu.initialize import xp, device

    # print("Create aligned images from alignResults.txt")

    alignmentResults = loadstar(alignmentResultsFile, dtype=datatypeAR)
    imageList = alignmentResults['FileName']
    tilt_angles = alignmentResults['TiltAngle']

    imdimX = read_size(imageList[0], 'x')
    imdimY = read_size(imageList[0], 'y')

    imdim = max(imdimX, imdimY)

    if binning > 1:
        imdimX = int(float(imdimX) / float(binning) + .5)
        imdimY = int(float(imdimY) / float(binning) + .5)

    imdim = max(imdimY, imdimX)
    sliceWidth = imdimX

    if (weighting != None) and (float(weighting) < -0.001):
        weightSlice = xp.fft.fftshift(ramp_filter(imdim, imdim))

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

        # weighting filter: arguments: (()dimx, dimy), cutoff radius, sigma
        # lpf = xp.fft.fftshift(create_circle((imdimX,imdimY),lowpassFilter*(imdim//2), sigma=0.4*lowpassFilter*(imdim//2)))

    projectionList = ProjectionList()
    for n, image in enumerate(imageList):
        atx = alignmentResults['AlignmentTransX'][n]
        aty = alignmentResults['AlignmentTransY'][n]
        rot = alignmentResults['InPlaneRotation'][n]
        mag = 1 / (alignmentResults['Magnification'][n] )
        projection = Projection(imageList[n], tiltAngle=tilt_angles[n], alignmentTransX=atx, alignmentTransY=aty,
                                alignmentRotation=rot, alignmentMagnification=mag)
        projectionList.append(projection)


    for (ii, projection) in enumerate(projectionList):
        if not ii == indexImage:
            continue
        from pytom.agnostic.transform import resize

        # print(f'read {projection._filename}')
        image = read(str(projection._filename)).squeeze()

        if binning > 1:
            image = resize(image, 1 / binning)

        #write(f'test/image_{ii}.mrc', image, tilt_angle=tilt_angles[ii])


        tiltAngle = projection._tiltAngle

        # 1 -- normalize to contrast - subtract mean and norm to mean
        immean = image.mean()
        image = (image - immean) / immean

        # 2 -- smoothen borders to prevent high contrast oscillations
        image = taper_edges(image, imdim // 30)[0]

        # 3 -- square if needed
        if imdimY != imdimX:
            newImage = xp.zeros((imdim, imdim), dtype=xp.float32)
            image = pasteCenter(image, newImage)



        # 4 -- transform projection according to tilt alignment
        transX = projection._alignmentTransX / binning
        transY = projection._alignmentTransY / binning
        rot = float(projection._alignmentRotation)
        mag = float(projection._alignmentMagnification)

        inputImage = xp.expand_dims(image, 2).copy()
        outputImage = xp.zeros_like(inputImage, dtype=xp.float32)


        vt.transform(inputImage.astype(xp.float32), rotation=[ 0, 0, rot], rotation_order='rxyz', output=outputImage,
                     center=[inputImage.shape[0]//2,inputImage.shape[1]//2,0],
                     device=device, translation=[transX, transY, 0], scale=[mag, mag, 1], interpolation='filt_bspline')

        del image
        image = outputImage.squeeze()

        # 5 -- Optional Low Pass Filter
        if lowpassFilter:
            from pytom.agnostic.filter import bandpass_circle

            image = bandpass_circle(image, high=lowpassFilter * (min(imdimX, imdimY) // 2),
                                    sigma=0.4 * lowpassFilter * (min(imdimX, imdimY) // 2))
            # image = xp.abs((xp.fft.ifftn(xp.fft.fftn(image) * lpf)))

        # 6 -- smoothen once more to avoid edges
        image = taper_edges(image, imdim // 30)[0]

        # 7 -- analytical weighting
        if (weighting != None) and (weighting < 0):
            print(weightSlice.shape, circleSlice.shape)
            # image = (ifft(complexRealMult(fft(image), w_func)) / (image.sizeX() * image.sizeY() * image.sizeZ()))
            image = xp.fft.ifftn(xp.fft.fftn(image) * weightSlice.T * circleSlice).real


        elif (weighting != None) and (weighting > 0):
            weightSlice = xp.fft.fftshift(exact_filter(tilt_angles, tiltAngle, imdim, imdim, sliceWidth))
            image = xp.fft.ifftn(xp.fft.fftn(image) * weightSlice * circleSlice).real

        del inputImage, outputImage, circleSlice

        return image.astype(xp.float32)


def alignImagesUsingAlignmentResultFile(alignmentResultsFile, weighting=None, lowpassFilter=0.9, binning=1,
                                        circleFilter=True, angle_specimen=0):
    import os
    from pytom.basic.files import read as readCVol
    from pytom_numpy import vol2npy, npy2vol
    from pytom.gui.guiFunctions import fmtAR, headerAlignmentResults, datatype, datatypeAR, loadstar
    from pytom.reconstruction.reconstructionStructures import Projection, ProjectionList
    from pytom.agnostic.io import read, write, read_size
    from pytom.agnostic.tools import taper_edges, create_circle, paste_in_center as pasteCenter
    from pytom.agnostic.filter import circle_filter, ramp_filter, exact_filter, ellipse_filter
    import pytom.voltools as vt
    from pytom.gpu.initialize import xp, device
    from pytom.agnostic.transform import resize
    from pytom.basic.files import read as RR
    from pytom.basic.transformations import general_transform2d
    import pytom_freqweight
    from pytom.basic.filter import filter as filterFunction
    from pytom_volume import complexRealMult
    from pytom.basic.fourier import ifft, fft

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


        # import matplotlib
        # matplotlib.use('Qt5Agg')
        # from pylab import imshow, show
        # imshow(image.get().real, cmap='gray')
        # show()

    arrays = []

    for fname, arr in (
    ('stack.mrc', stack), ('offsetStack.mrc', phiStack), ('thetaStack.mrc', thetaStack), ('phiStack.mrc', offsetStack)):
    #     if 'gpu' in device:
    #         arr = arr.get()
    #     import numpy as np
    #     if 'cpu' in device:
    #         res = npy2vol(np.array(arr, dtype='float32', order='F'), 3)
    #         arrays.append(res)
    #     else:
        arrays.append(arr)
    #
    #     write('stack.mrc', stack)
    #     stack = readCVol('stack.mrc')
    # write('offsetstack.mrc', offsetStack)
    # offsetStack = readCVol('offsetstack.mrc')
    # write('thetastack.mrc', thetaStack)
    # thetaStack = readCVol('thetastack.mrc')
    # write('phistack.mrc', phiStack)
    # phiStack = readCVol('phistack.mrc')
    #
    # os.remove('stack.mrc')
    # os.remove('offsetstack.mrc')
    # os.remove('thetastack.mrc')
    # os.remove('psistack.mrc')




    return arrays


def toProjectionStackFromAlignmentResultsFile(alignmentResultsFile, weighting=None, lowpassFilter=0.9,
                                              binning=1, circleFilter=False, num_procs=1, outdir='',
                                              prefix='sorted_aligned'):
    """read image and create aligned projection stack, based on the results described in the alignmentResultFile.

       @param alignmentResultsFile: result file generate by the alignment script.
       @type datatypeAR: gui.guiFunction.datatypeAR
       @param weighting: weighting (<0: analytical weighting, >1: exact weighting, 0/None: no weighting )
       @type weighting: float
       @param lowpassFilter: lowpass filter (in Nyquist)
       @type lowpassFilter: float
       @param binning: binning (default: 1 = no binning). binning=2: 2x2 pixels -> 1 pixel, binning=3: 3x3 pixels -> 1 pixel, etc.

       @author: GvdS
    """
    print('weighting: ', weighting)
    import numpy
    from pytom_numpy import vol2npy
    from pytom.basic.files import read_em, write_em
    from pytom.basic.functions import taper_edges
    from pytom.basic.transformations import general_transform2d
    from pytom.basic.fourier import ifft, fft
    from pytom.basic.filter import filter as filterFunction, bandpassFilter
    from pytom.basic.filter import circleFilter, rampFilter, exactFilter, fourierFilterShift, \
        fourierFilterShift_ReducedComplex
    from pytom_volume import complexRealMult, vol, paste
    import pytom_freqweight
    from pytom.basic.transformations import resize, rotate
    from pytom.gui.guiFunctions import fmtAR, headerAlignmentResults, datatype, datatypeAR, loadstar
    from pytom.reconstruction.reconstructionStructures import Projection, ProjectionList
    from pytom_numpy import vol2npy
    import mrcfile
    from pytom.agnostic.io import write, read_size
    import os

    print("Create aligned images from alignResults.txt")

    alignmentResults = loadstar(alignmentResultsFile, dtype=datatypeAR)
    imageList = alignmentResults['FileName']
    tilt_angles = alignmentResults['TiltAngle']

    imdim = int(read_size(imageList[0], 'x'))

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
    imageList = []
    tilt_angles = []
    for n, image in enumerate(alignmentResults['FileName']):
        atx = alignmentResults['AlignmentTransX'][n]
        aty = alignmentResults['AlignmentTransY'][n]
        rot = alignmentResults['InPlaneRotation'][n]
        mag = alignmentResults['Magnification'][n]
        # print(image, alignmentResults['TiltAngle'][n])
        # if abs(alignmentResults['TiltAngle'][n]) > 20:
        #     continue
        tilt_angles.append(alignmentResults['TiltAngle'][n])
        imageList.append(image)
        projection = Projection(imageList[-1], tiltAngle=tilt_angles[-1], alignmentTransX=atx, alignmentTransY=aty,
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

        #print(ii, immean, projection._filename)

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
            # image = (ifft(complexRealMult(fft(image), w_func)) / (image.sizeX() * image.sizeY() * image.sizeZ()))
            image = ifft(complexRealMult(complexRealMult(fft(image), weightSlice), circleSlice), scaling=True)


        elif (weighting != None) and (weighting > 0):
            weightSlice = fourierFilterShift(exactFilter(tilt_angles, tiltAngle, imdim, imdim, sliceWidth))
            # image = (ifft(complexRealMult(fft(image), w_func)) / (image.sizeX() * image.sizeY() * image.sizeZ()))
            image = ifft(complexRealMult(complexRealMult(fft(image), weightSlice), circleSlice), scaling=True)

        thetaStack(int(round(projection.getTiltAngle())), 0, 0, ii)
        offsetStack(int(round(projection.getOffsetX())), 0, 0, ii)
        offsetStack(int(round(projection.getOffsetY())), 0, 1, ii)
        paste(image, stack, 0, 0, ii)
        fname = '{}_{:02d}.mrc'.format(prefix, int(imageList[ii].split('_')[-1].split('.')[0]))
        if outdir:
            import mrcfile
            # write_em(os.path.join(outdir, fname.replace('mrc', 'em')), image)
            write(os.path.join(outdir, fname), vol2npy(image).copy().astype('float32'))
            print('written file: ', fname)

    return [stack, phiStack, thetaStack, offsetStack]
