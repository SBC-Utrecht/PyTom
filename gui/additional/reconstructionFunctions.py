'''
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
                           voldims=None, recCent=[0,0,0], tiltSeriesFormat='st', firstProj=1, irefmark=1, ireftilt=1,
                           handflip=False, alignedTiltSeriesName='align/myTilt', weightingType=-1,
                           lowpassFilter=1., projBinning=1, outMarkerFileName=None, verbose=False, outfile='',
                           write_images=True):
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
    from pytom.gui.additional.TiltAlignmentStructures import TiltSeries, TiltAlignment, TiltAlignmentParameters

    if not alignResultFile:
        if verbose:
            print("Function alignWeightReconstruct started")
            mute = False
        else:
            mute = True
        from pytom.gui.additional.TiltAlignmentStructures import TiltSeries, TiltAlignment, TiltAlignmentParameters
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
                print(" EM markerfile file used for alignment")
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
        print(outfile)
        tiltAlignment.computeCoarseAlignment(tiltSeries, mute=mute, outfile=outfile)
        tiltAlignment.alignFromFiducials(mute=mute)
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
            vol_bp = tiltSeries.reconstructVolume(dims=voldims, reconstructionPosition=recCent, binning=1)

            vol_bp.write(volumeName, volumeFileType)

    else:
        print('new code')
        tiltSeries = TiltSeries(tiltSeriesName=tiltSeriesName)
        vol_bp = tiltSeries.reconstructVolume(dims=voldims, reconstructionPosition=recCent, binning=projBinning,
                                              alignResultFile=alignResultFile)

        vol_bp.write(volumeName, volumeFileType)


