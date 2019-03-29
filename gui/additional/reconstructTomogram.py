#!/usr/bin/env pytom
"""
Created on Jul 20, 2013

@author: FF
"""

if __name__ == '__main__':
    import sys
    #from pytom.reconstruction.TiltAlignmentStructures import TiltAlignmentParameters, TiltSeries, TiltAlignment
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.gui.additional.reconstructionFunctions import alignWeightReconstruct

    options=[ScriptOption(['--tiltSeriesName'], 'Name tilt series - either prefix of sequential tilt series files \
             expected as "tiltSeriesName_index.em/mrc" or full name of stack "tiltSeriesName.st"',
                          arg=True, optional=False),
             ScriptOption(['--tiltSeriesFormat'], 'Format of tilt series (series of "em" or "mrc" images or "st" stack).',
                          arg=True, optional=True),
             ScriptOption(['--firstIndex'], 'Index of first projection.', arg=True, optional=True),
             ScriptOption(['--lastIndex'], 'Index of last projection.', arg=True, optional=True),
             ScriptOption(['--projIndices'], 'Use numbering in filename as index', arg=False, optional=True),
             ScriptOption(['--tltFile'], 'tltFile containing tilt angles.', arg=True, optional=True),
             ScriptOption(['--prexgFile'], 'prexgFile containing pre-shifts from IMOD.', arg=True, optional=True),
             ScriptOption(['--preBin'], 'pre-Binning in IMOD prior to marker determination.', arg=True, optional=True),
             ScriptOption(['--referenceIndex'], 'Index of reference projection used for alignment.', arg=True,
                          optional=False),
             ScriptOption(['--markerFile'], 'Name of EM markerfile or IMOD wimp File containing marker coordinates.',
                          arg=True, optional=False),
             ScriptOption(['--referenceMarkerIndex'], 'Index of reference marker to set up coordinate system.',
                          arg=True, optional=False),
             ScriptOption(['--handflip'], 'Is your tilt series outside of 0-180deg (Specify if yes).', arg=False,
                          optional=True),
             ScriptOption(['--projectionTargets'],
                          'Relative or absolute path to the aligned projections that will be generated + file prefix.\
                          default: "align/myTilt"', arg=True, optional=True),
             ScriptOption(['--fineAlignFile'],
                          'Relative or absolute path to the file with fineAlign parameters (type should be *.dat).',
                          arg=True, optional=True),
             ScriptOption(['--projectionBinning'], 'Binning of projections during read - default: 1.', arg=True,
                          optional=True),
             ScriptOption(['--lowpassFilter'], 'Lowpass filter in Nyquist after binning.', arg=True, optional=False),
             ScriptOption(['--tomogramFile'],
                          'Relative or absolute path to final tomogram (no tomogram written if not specified).',
                          arg=True, optional=True),
             ScriptOption(['--fileType'], 'File type (can be em or mrc - no tomogram written if not specified).',
                          arg=True, optional=True),
             ScriptOption(['--tomogramSizeX'], 'Size of tomogram in x (no tomogram written if not specified).',
                          arg=True, optional=True),
             ScriptOption(['--tomogramSizeY'], 'Size of tomogram in y (no tomogram written if not specified).',
                          arg=True, optional=True),
             ScriptOption(['--tomogramSizeZ'], 'Size of tomogram in z (no tomogram written if not specified).',
                          arg=True, optional=True),
             ScriptOption(['--reconstructionCenterX'],
                          'Center where tomogram will be reconstructed (no tomogram written if not specified).',
                          arg=True, optional=True),
             ScriptOption(['--reconstructionCenterY'],
                          'Center where tomogram will be reconstructed (no tomogram written if not specified).',
                          arg=True, optional=True),
             ScriptOption(['--reconstructionCenterZ'],
                          'Center where tomogram will be reconstructed (no tomogram written if not specified).',
                          arg=True, optional=True),
             ScriptOption(['--weightingType'], 'Type of weighting (-1 default r-weighting, 0 no weighting)', arg=True,
                          optional=True),
             ScriptOption(['--verbose'], 'Enable verbose mode', arg=False, optional=True),
             ScriptOption(['-h', '--help'], 'Help.', False, True)]
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1],
                          description='Align and weight projections, save them and reconstruct tomogram (optional). \n\
                                      See http://pytom.org/doc/pytom/reconstructTomograms.html for documentation.',
                          authors='Friedrich Foerster',
                          options = options)

    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        tiltSeriesName, tiltSeriesFormat, firstProj, lastProj, projIndices,\
        tltFile, prexgFile, preBin, referenceIndex, markerFileName, referenceMarkerIndex, handflip, \
        projectionTargets, fineAlignFile, projBinning, lowpassFilter, \
        volumeName, filetype, \
        tomogramSizeX, tomogramSizeY, tomogramSizeZ, \
        reconstructionCenterX, reconstructionCenterY, reconstructionCenterZ, \
        weightingType, verbose, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(sys.version_info)
        print(e)
        sys.exit()
        
    if help is True:
        print(helper)
        sys.exit()
    # input parameters
    #tiltSeriesName = tiltSeriesPath + tiltSeriesPrefix  # "../projections/tomo01_sorted" # ending is supposed to be tiltSeriesName_index.em (or mrc)
    if not tiltSeriesFormat:
        tiltSeriesFormat = 'em'
    if firstProj:
        firstProj = int(firstProj)  # 1 # index of first projection
    else:
        firstProj = 1
    lastProj = int(lastProj)  # 41 # index of last projection
    ireftilt = int(referenceIndex)  # 21 # reference projection (used for alignment)
    if referenceMarkerIndex:
        irefmark = int(referenceMarkerIndex)  # reference marker (defines 3D coordinate system)
    else:
        irefmark = 1

    handflip = handflip is not None  # False # is your tilt axis outside 0-180 deg?
    # output parameters
    if projectionTargets:
        # weighted and aligned projections are stored as alignedTiltSeriesName_index.em
        alignedTiltSeriesName = projectionTargets
    else:
        alignedTiltSeriesName = 'align/myTilt'
    if projBinning:
        projBinning = int(projBinning)  # binning factor
    else:
        projBinning = 1
    if lowpassFilter:
        lowpassFilter = float(lowpassFilter)  # lowpass filter in Nyquist (post-binning)
    else:
        lowpassFilter = 1.

    if weightingType is None:
        weightingType = -1
    else:
        weightingType = int(weightingType)
    # only write projections and do NOT reconstruct tomogram (following parameters would be obsolete)
    onlyWeightedProjections = False
    if not volumeName:
        onlyWeightedProjections = True
         
    if filetype is None:
        if volumeName:
            filetype = volumeName.split('.')[-1]
        else:
            filetype = tiltSeriesFormat

    if tomogramSizeX is not None or tomogramSizeY is not None or tomogramSizeZ is not None:
        # dimensions of reconstructed tomogram
        voldims = [int(tomogramSizeX), int(tomogramSizeY), int(tomogramSizeZ)]
    else:
        onlyWeightedProjections = True
        voldims = [0, 0, 0]       # offset from center of volume - for example choose z!=0 to shift in z (post-binning coordinates)
    if reconstructionCenterX:
        reconstructionCenterX = int(reconstructionCenterX)
    else:
        reconstructionCenterX = 0
    if reconstructionCenterY:
        reconstructionCenterY = int(reconstructionCenterY)
    else:
        reconstructionCenterY = 0
    if reconstructionCenterZ:
        reconstructionCenterZ = int(reconstructionCenterZ)
    else:
        reconstructionCenterZ = 0
    reconstructionPosition = [reconstructionCenterX, reconstructionCenterY, reconstructionCenterZ]
    if preBin:
        preBin=int(preBin)

    outMarkerFileName = 'MyMarkerFile.em'
    if verbose:
        print("Tilt Series: "+str(tiltSeriesName)+", "+str(firstProj)+"-"+str(lastProj))
        print("Index of Reference Projection: "+str(referenceIndex))
        print("Marker Filename: "+str(markerFileName))
        print("TltFile: "+str(tltFile))
        print("prexgFile: "+str(prexgFile))
        print("Index of Reference Marker: "+str(referenceMarkerIndex))
        print("Handflip: "+str(handflip))
        print("Projection Targets: "+str(projectionTargets))
        print("FineAlignmentFile: "+str(fineAlignFile))
        print("Binning Factor of Projections: "+str(projBinning)+", lowpass filter (in Ny): "+str(lowpassFilter))
        print("Name of Reconstruction Volume: "+str(volumeName)+" of Filetype: "+str(filetype))
        print("Reconstruction size: "+str(voldims))
        print("Reconstruction center: "+str(reconstructionPosition))
        print("write only aligned projections out: "+str(onlyWeightedProjections))

    alignWeightReconstruct(tiltSeriesName=tiltSeriesName, markerFileName=markerFileName, lastProj=lastProj,
                           tltfile=tltFile, prexgfile=prexgFile, preBin=preBin,
                           volumeName=volumeName, volumeFileType=filetype,
                           voldims=voldims, recCent=reconstructionPosition,
                           tiltSeriesFormat=tiltSeriesFormat, firstProj=firstProj, irefmark=irefmark, ireftilt=ireftilt,
                           handflip=handflip,
                           alignedTiltSeriesName=alignedTiltSeriesName,
                           weightingType=weightingType,
                           lowpassFilter=lowpassFilter, projBinning=projBinning,
                           outMarkerFileName=outMarkerFileName, verbose=True)


