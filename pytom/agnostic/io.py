#! /usr/bin/env python
# -*- coding: utf-8 -*-

from pytom.gpu.initialize import xp, device
import numpy as np
#typing imports
import pytom.lib as pytom_lib

# Reading functions
from typing import Union
from pytom.gpu.initialize import xpt

def read(filename, ndarray=True, order='F', keepnumpy=False, deviceID=None, dtype=None, read_optics_group=False) -> Union[xpt.NDArray[float], pytom_lib.pytom_volume.vol]:
    """General reading function. Can read em, mrc, st, rec, txt, log and star file. For EM and MRC files: only support read the type float32 on little-endian machines.

    @param filename: file name to read.
    @param ndarray: if False data is converted to PyTom vol. Default = True
    @param order: order of the ndarray. Default = Fortran
    @param keepnumpy: return numpy.ndarray if True, else class depends on the init off pytom. Could be cupy.ndarray as well. Default = False
    @param deviceID: optional kw: GPU_ID of device to which the data is copied if requested. Default = None
    @param dtype: optional kw, used for reading log, star and txt files. Default = None

    @return The data from the EM file in ndarray
    """

    import os

    allowed_formats = {'em': read_em,
                       'mrc': read_mrc, 'rec': read_mrc, 'st': read_mrc,
                       'txt': read_txt, 'log': read_txt, 'meta': read_txt,
                       'star': read_star if not read_optics_group else read_star_optics_group}

    name, ext = os.path.splitext(filename)
    ext = ext.lower()[1:]

    assert filename
    assert os.path.exists(filename)
    assert ext in allowed_formats.keys()

    data = allowed_formats[ext](filename, order, keepnumpy, deviceID, dtype)

    if ndarray:
        return data
    else:
        return n2v(data)


def read_star(filename, order='F', keepnumpy=False, deviceID=None, dtype=None):
    '''Read star file for relion 3.1 or relion3.0 or user defined dtype.'''
    from pytom.basic.files import loadtxt

    size_optics_group = 0

    if dtype is None:

        relion_31_options = {'AccumMotionEarly': np.double,
                             'AccumMotionLate': np.double,
                             'AccumMotionTotal': np.double,
                             'AccuracyRotations': np.double,
                             'AccuracyTranslations': np.double,
                             'AccuracyTranslationsAngst': np.double,
                             'AdaptiveOversampleFraction': np.double,
                             'AdaptiveOversampleOrder': 'i4',
                             'AmplitudeContrast': np.double,
                             'AmplitudeCorrelationMaskedMaps': np.double,
                             'AmplitudeCorrelationUnmaskedMaps': np.double,
                             'AnglePsi': np.double,
                             'AnglePsiFlip': bool,
                             'AnglePsiFlipRatio': np.double,
                             'AnglePsiPrior': np.double,
                             'AngleRot': np.double,
                             'AngleRotFlipRatio': np.double,
                             'AngleRotPrior': np.double,
                             'AngleTilt': np.double,
                             'AngleTiltPrior': np.double,
                             'AngstromResolution': np.double,
                             'AreaId': 'i4',
                             'AreaName': 'U200',
                             'AutoLocalSearchesHealpixOrder': 'i4',
                             'AutopickFigureOfMerit': np.double,
                             'AvailableMemory': np.double,
                             'AverageNrOfFrames': 'i4',
                             'AveragePmax': np.double,
                             'AverageValue': np.double,
                             'BeamTiltClass': 'i4',
                             'BeamTiltX': np.double,
                             'BeamTiltY': np.double,
                             'BestResolutionThusFar': np.double,
                             'BfactorUsedForSharpening': np.double,
                             'BodyKeepFixed': 'i4',
                             'BodyMaskName': 'U200',
                             'BodyReferenceName': 'U200',
                             'BodyRotateDirectionX': np.double,
                             'BodyRotateDirectionY': np.double,
                             'BodyRotateDirectionZ': np.double,
                             'BodyRotateRelativeTo': 'i4',
                             'BodySigmaAngles': np.double,
                             'BodySigmaOffset': np.double,
                             'BodySigmaOffsetAngst': np.double,
                             'BodySigmaPsi': np.double,
                             'BodySigmaRot': np.double,
                             'BodySigmaTilt': np.double,
                             'BodyStarFile': 'U200',
                             'ChangesOptimalClasses': np.double,
                             'ChangesOptimalOffsets': np.double,
                             'ChangesOptimalOrientations': np.double,
                             'ChromaticAberration': np.double,
                             'ClassDistribution': np.double,
                             'ClassNumber': 'i4',
                             'ClassPriorOffsetX': np.double,
                             'ClassPriorOffsetY': np.double,
                             'CoarseImageSize': 'i4',
                             'Comment': 'U200',
                             'ConvergenceCone': np.double,
                             'CoordinateX': np.double,
                             'CoordinateY': np.double,
                             'CoordinateZ': np.double,
                             'CorrectedFourierShellCorrelationPhaseRandomizedMaskedMaps': np.double,
                             'CorrelationFitGuinierPlot': np.double,
                             'CtfAstigmatism': np.double,
                             'CtfBfactor': np.double,
                             'CtfDataAreCtfPremultiplied': bool,
                             'CtfDataArePhaseFlipped': bool,
                             'CtfFigureOfMerit': np.double,
                             'CtfImage': 'U200',
                             'CtfMaxResolution': np.double,
                             'CtfPowerSpectrum': 'U200',
                             'CtfScalefactor': np.double,
                             'CtfValidationScore': np.double,
                             'CtfValue': np.double,
                             'CurrentImageSize': 'i4',
                             'CurrentIteration': 'i4',
                             'CurrentResolution': np.double,
                             'DataDimensionality': 'i4',
                             'DataType': 'i4',
                             'DefocusAngle': np.double,
                             'DefocusU': np.double,
                             'DefocusV': np.double,
                             'DetectorPixelSize': np.double,
                             'Diff2RandomHalves': np.double,
                             'DifferentialPhaseResidualMaskedMaps': np.double,
                             'DifferentialPhaseResidualUnmaskedMaps': np.double,
                             'DoAutoRefine': bool,
                             'DoCorrectCtf': bool,
                             'DoCorrectMagnification': bool,
                             'DoCorrectNorm': bool,
                             'DoCorrectScale': bool,
                             'DoExternalReconstruct': bool,
                             'DoFastSubsetOptimisation': bool,
                             'DoHelicalRefine': bool,
                             'DoIgnoreCtfUntilFirstPeak': bool,
                             'DoMapEstimation': bool,
                             'DoOnlyFlipCtfPhases': bool,
                             'DoRealignMovies': bool,
                             'DoSkipAlign': bool,
                             'DoSkipRotate': bool,
                             'DoSolventFlattening': bool,
                             'DoSolventFscCorrection': bool,
                             'DoSplitRandomHalves': bool,
                             'DoStochasticEM': bool,
                             'DoStochasticGradientDescent': bool,
                             'DoZeroMask': bool,
                             'EERGrouping': 'i4',
                             'EERUpsampling': 'i4',
                             'Enabled': bool,
                             'EnergyLoss': np.double,
                             'EstimatedResolution': np.double,
                             'EvenZernike': np.double,
                             'ExperimentalDataStarFile': 'U200',
                             'ExtReconsDataImag': 'U200',
                             'ExtReconsDataReal': 'U200',
                             'ExtReconsResult': 'U200',
                             'ExtReconsResultStarfile': 'U200',
                             'ExtReconsWeight': 'U200',
                             'FinalResolution': np.double,
                             'FittedInterceptGuinierPlot': np.double,
                             'FittedSlopeGuinierPlot': np.double,
                             'FixSigmaNoiseEstimates': bool,
                             'FixSigmaOffsetEstimates': bool,
                             'FixTauEstimates': bool,
                             'FourierCompleteness': np.double,
                             'FourierMask': 'U200',
                             'FourierShellCorrelation': np.double,
                             'FourierShellCorrelationCorrected': np.double,
                             'FourierShellCorrelationMaskedMaps': np.double,
                             'FourierShellCorrelationParticleMaskFraction': np.double,
                             'FourierShellCorrelationParticleMolWeight': np.double,
                             'FourierShellCorrelationUnmaskedMaps': np.double,
                             'FourierSpaceInterpolator': 'i4',
                             'GoldStandardFsc': np.double,
                             'GroupName': 'U200',
                             'GroupNrParticles': 'i4',
                             'GroupNumber': 'i4',
                             'GroupScaleCorrection': np.double,
                             'HasConverged': bool,
                             'HasHighFscAtResolLimit': bool,
                             'HasLargeSizeIncreaseIterationsAgo': 'i4',
                             'HealpixOrder': 'i4',
                             'HealpixOrderOriginal': 'i4',
                             'HelicalCentralProportion': np.double,
                             'HelicalKeepTiltPriorFixed': bool,
                             'HelicalMaskTubeInnerDiameter': np.double,
                             'HelicalMaskTubeOuterDiameter': np.double,
                             'HelicalOffsetStep': np.double,
                             'HelicalRise': np.double,
                             'HelicalRiseInitial': np.double,
                             'HelicalRiseInitialStep': np.double,
                             'HelicalRiseMax': np.double,
                             'HelicalRiseMin': np.double,
                             'HelicalSigmaDistance': np.double,
                             'HelicalSymmetryLocalRefinement': bool,
                             'HelicalTrackLength': np.double,
                             'HelicalTrackLengthAngst': np.double,
                             'HelicalTubeID': 'i4',
                             'HelicalTubePitch': np.double,
                             'HelicalTwist': np.double,
                             'HelicalTwistInitial': np.double,
                             'HelicalTwistInitialStep': np.double,
                             'HelicalTwistMax': np.double,
                             'HelicalTwistMin': np.double,
                             'HighresLimitExpectation': np.double,
                             'HighresLimitSGD': np.double,
                             'IgnoreHelicalSymmetry': bool,
                             'ImageDimensionality': 'i4',
                             'ImageId': 'i4',
                             'ImageName': 'U200',
                             'ImageOriginalName': 'U200',
                             'ImagePixelSize': np.double,
                             'ImageSize': 'i4',
                             'ImageSizeX': 'i4',
                             'ImageSizeY': 'i4',
                             'ImageSizeZ': 'i4',
                             'ImageWeight': np.double,
                             'IncrementImageSize': 'i4',
                             'Is3DSampling': bool,
                             'Is3DTranslationalSampling': bool,
                             'IsFlip': bool,
                             'IsHelix': bool,
                             'JobIsContinue': bool,
                             'JobOptionDefaultValue': 'U200',
                             'JobOptionDirectoryDefault': 'U200',
                             'JobOptionFilePattern': 'U200',
                             'JobOptionGUILabel': 'U200',
                             'JobOptionHelpText': 'U200',
                             'JobOptionMenuOptions': 'U200',
                             'JobOptionSliderMax': np.double,
                             'JobOptionSliderMin': np.double,
                             'JobOptionSliderStep': np.double,
                             'JobOptionValue': 'U200',
                             'JobOptionVariable': 'U200',
                             'JobType': 'i4',
                             'JobTypeName': 'U200',
                             'JoboptionType': 'i4',
                             'JoinHalvesUntilThisResolution': np.double,
                             'KullbackLeiblerDivergence': np.double,
                             'KurtosisExcessValue': np.double,
                             'LensStability': np.double,
                             'LocalSymmetryFile': 'U200',
                             'LogAmplitudesIntercept': np.double,
                             'LogAmplitudesMTFCorrected': np.double,
                             'LogAmplitudesOriginal': np.double,
                             'LogAmplitudesSharpened': np.double,
                             'LogAmplitudesWeighted': np.double,
                             'LogLikeliContribution': np.double,
                             'LogLikelihood': np.double,
                             'LongitudinalDisplacement': np.double,
                             'LowresLimitExpectation': np.double,
                             'MagMat00': np.double,
                             'MagMat01': np.double,
                             'MagMat10': np.double,
                             'MagMat11': np.double,
                             'Magnification': np.double,
                             'MagnificationCorrection': np.double,
                             'MagnificationSearchRange': np.double,
                             'MagnificationSearchStep': np.double,
                             'MaskName': 'U200',
                             'Matrix_1_1': np.double,
                             'Matrix_1_2': np.double,
                             'Matrix_1_3': np.double,
                             'Matrix_2_1': np.double,
                             'Matrix_2_2': np.double,
                             'Matrix_2_3': np.double,
                             'Matrix_3_1': np.double,
                             'Matrix_3_2': np.double,
                             'Matrix_3_3': np.double,
                             'MaxNumberOfPooledParticles': 'i4',
                             'MaxValueProbDistribution': np.double,
                             'MaximumCoarseImageSize': 'i4',
                             'MaximumValue': np.double,
                             'MicrographBinning': np.double,
                             'MicrographDefectFile': 'U200',
                             'MicrographDoseRate': np.double,
                             'MicrographEndFrame': 'i4',
                             'MicrographFrameNumber': 'i4',
                             'MicrographGainName': 'U200',
                             'MicrographId': 'i4',
                             'MicrographMetadata': 'U200',
                             'MicrographMovieName': 'U200',
                             'MicrographName': 'U200',
                             'MicrographNameNoDW': 'U200',
                             'MicrographOriginalPixelSize': np.double,
                             'MicrographPixelSize': np.double,
                             'MicrographPreExposure': np.double,
                             'MicrographShiftX': np.double,
                             'MicrographShiftY': np.double,
                             'MicrographStartFrame': 'i4',
                             'MicrographTiltAngle': np.double,
                             'MicrographTiltAxisDirection': np.double,
                             'MicrographTiltAxisOutOfPlane': np.double,
                             'MinRadiusNnInterpolation': 'i4',
                             'MinimumValue': np.double,
                             'ModelStarFile': 'U200',
                             'ModelStarFile2': 'U200',
                             'MolecularWeight': np.double,
                             'MotionModelCoeff': np.double,
                             'MotionModelCoeffsIdx': 'i4',
                             'MotionModelVersion': 'i4',
                             'MovieFrameNumber': 'i4',
                             'MovieFramesRunningAverage': 'i4',
                             'MtfFileName': 'U200',
                             'MtfValue': np.double,
                             'NormCorrection': np.double,
                             'NormCorrectionAverage': np.double,
                             'NrBodies': 'i4',
                             'NrClasses': 'i4',
                             'NrGroups': 'i4',
                             'NrHelicalAsymUnits': 'i4',
                             'NrHelicalNStart': 'i4',
                             'NrOfFrames': 'i4',
                             'NrOfSignificantSamples': 'i4',
                             'NumberOfIterWithoutChangingAssignments': 'i4',
                             'NumberOfIterWithoutResolutionGain': 'i4',
                             'NumberOfIterations': 'i4',
                             'OddZernike': np.double,
                             'OffsetRange': np.double,
                             'OffsetRangeOriginal': np.double,
                             'OffsetStep': np.double,
                             'OffsetStepOriginal': np.double,
                             'OpticsGroup': 'i4',
                             'OpticsGroupName': 'U200',
                             'OpticsStarFile': 'U200',
                             'OrientSamplingStarFile': 'U200',
                             'OrientationDistribution': np.double,
                             'OrientationalPriorMode': 'i4',
                             'OrientationsID': 'i4',
                             'OriginX': np.double,
                             'OriginXAngst': np.double,
                             'OriginXPrior': np.double,
                             'OriginXPriorAngst': np.double,
                             'OriginY': np.double,
                             'OriginYAngst': np.double,
                             'OriginYPrior': np.double,
                             'OriginYPriorAngst': np.double,
                             'OriginZ': np.double,
                             'OriginZAngst': np.double,
                             'OriginZPrior': np.double,
                             'OriginZPriorAngst': np.double,
                             'OriginalImageSize': 'i4',
                             'OriginalParticleName': 'U200',
                             'OutputRootName': 'U200',
                             'OverallAccuracyRotations': np.double,
                             'OverallAccuracyTranslations': np.double,
                             'OverallAccuracyTranslationsAngst': np.double,
                             'OverallFourierCompleteness': np.double,
                             'PaddingFactor': np.double,
                             'ParticleBoxFractionMolecularWeight': np.double,
                             'ParticleBoxFractionSolventMask': np.double,
                             'ParticleDiameter': np.double,
                             'ParticleFigureOfMerit': np.double,
                             'ParticleId': 'i4',
                             'ParticleName': 'U200',
                             'ParticleNumber': 'i4',
                             'ParticleSelectZScore': np.double,
                             'PerFrameCumulativeWeight': np.double,
                             'PerFrameRelativeWeight': np.double,
                             'PhaseShift': np.double,
                             'PipeLineEdgeFromNode': 'U200',
                             'PipeLineEdgeProcess': 'U200',
                             'PipeLineEdgeToNode': 'U200',
                             'PipeLineJobCounter': 'i4',
                             'PipeLineNodeName': 'U200',
                             'PipeLineNodeType': 'i4',
                             'PipeLineProcessAlias': 'U200',
                             'PipeLineProcessName': 'U200',
                             'PipeLineProcessStatus': 'i4',
                             'PipeLineProcessType': 'i4',
                             'PixelSize': np.double,
                             'PsiStep': np.double,
                             'PsiStepOriginal': np.double,
                             'RadiusMaskExpImages': 'i4',
                             'RadiusMaskMap': 'i4',
                             'RandomSeed': 'i4',
                             'RandomSubset': 'i4',
                             'RandomiseFrom': np.double,
                             'ReconstructImageName': 'U200',
                             'ReferenceDimensionality': 'i4',
                             'ReferenceImage': 'U200',
                             'ReferenceSigma2': np.double,
                             'ReferenceSpectralPower': np.double,
                             'ReferenceTau2': np.double,
                             'RefsAreCtfCorrected': bool,
                             'Resolution': np.double,
                             'ResolutionInversePixel': np.double,
                             'ResolutionSquared': np.double,
                             'SGDGradientImage': 'U200',
                             'SamplingPerturbFactor': np.double,
                             'SamplingPerturbInstance': np.double,
                             'SamplingRate': np.double,
                             'SamplingRateX': np.double,
                             'SamplingRateY': np.double,
                             'SamplingRateZ': np.double,
                             'ScheduleBooleanVariableName': 'U200',
                             'ScheduleBooleanVariableResetValue': bool,
                             'ScheduleBooleanVariableValue': bool,
                             'ScheduleCurrentNodeName': 'U200',
                             'ScheduleEdgeBooleanVariable': 'U200',
                             'ScheduleEdgeInputNodeName': 'U200',
                             'ScheduleEdgeIsFork': bool,
                             'ScheduleEdgeNumber': 'i4',
                             'ScheduleEdgeOutputNodeName': 'U200',
                             'ScheduleEdgeOutputNodeNameIfTrue': 'U200',
                             'ScheduleEmailAddress': 'U200',
                             'ScheduleFloatVariableName': 'U200',
                             'ScheduleFloatVariableResetValue': np.double,
                             'ScheduleFloatVariableValue': np.double,
                             'ScheduleJobHasStarted': bool,
                             'ScheduleJobMode': 'U200',
                             'ScheduleJobName': 'U200',
                             'ScheduleJobNameOriginal': 'U200',
                             'ScheduleName': 'U200',
                             'ScheduleOperatorInput1': 'U200',
                             'ScheduleOperatorInput2': 'U200',
                             'ScheduleOperatorName': 'U200',
                             'ScheduleOperatorOutput': 'U200',
                             'ScheduleOperatorType': 'U200',
                             'ScheduleOriginalStartNodeName': 'U200',
                             'ScheduleStringVariableName': 'U200',
                             'ScheduleStringVariableResetValue': 'U200',
                             'ScheduleStringVariableValue': 'U200',
                             'Selected': 'i4',
                             'SgdFinalIterations': 'i4',
                             'SgdFinalResolution': np.double,
                             'SgdFinalSubsetSize': 'i4',
                             'SgdInBetweenIterations': 'i4',
                             'SgdInitialIterations': 'i4',
                             'SgdInitialResolution': np.double,
                             'SgdInitialSubsetSize': 'i4',
                             'SgdMaxSubsets': 'i4',
                             'SgdMuFactor': np.double,
                             'SgdSigma2FudgeHalflife': 'i4',
                             'SgdSigma2FudgeInitial': np.double,
                             'SgdSkipAnneal': bool,
                             'SgdStepsize': np.double,
                             'SgdSubsetSize': 'i4',
                             'SgdWriteEverySubset': 'i4',
                             'Sigma2Noise': np.double,
                             'SigmaOffsets': np.double,
                             'SigmaOffsetsAngst': np.double,
                             'SigmaPriorPsiAngle': np.double,
                             'SigmaPriorRotAngle': np.double,
                             'SigmaPriorTiltAngle': np.double,
                             'SignalToNoiseRatio': np.double,
                             'SkewnessValue': np.double,
                             'SmallestChangesClasses': 'i4',
                             'SmallestChangesOffsets': np.double,
                             'SmallestChangesOrientations': np.double,
                             'SolventMask2Name': 'U200',
                             'SolventMaskName': 'U200',
                             'SortedIndex': 'i4',
                             'SpectralIndex': 'i4',
                             'SpectralOrientabilityContribution': np.double,
                             'SphericalAberration': np.double,
                             'SsnrMap': np.double,
                             'StandardDeviationValue': np.double,
                             'StarFileMovieParticles': 'U200',
                             'SymmetryGroup': 'U200',
                             'Tau2FudgeFactor': np.double,
                             'TauSpectrumName': 'U200',
                             'TiltAngleLimit': np.double,
                             'TransversalDisplacement': np.double,
                             'UnfilteredMapHalf1': 'U200',
                             'UnfilteredMapHalf2': 'U200',
                             'UnknownLabel': 'U200',
                             'UseTooCoarseSampling': bool,
                             'Voltage': np.double,
                             'WidthMaskEdge': 'i4'}

        relion_30_options = {'AccumMotionEarly': np.double,
                             'AccumMotionLate': np.double,
                             'AccumMotionTotal': np.double,
                             'AccuracyRotations': np.double,
                             'AccuracyTranslations': np.double,
                             'AdaptiveOversampleFraction': np.double,
                             'AdaptiveOversampleOrder': 'i4',
                             'AmplitudeContrast': np.double,
                             'AmplitudeCorrelationMaskedMaps': np.double,
                             'AmplitudeCorrelationUnmaskedMaps': np.double,
                             'AnglePsi': np.double,
                             'AnglePsiFlipRatio': np.double,
                             'AnglePsiPrior': np.double,
                             'AngleRot': np.double,
                             'AngleRotPrior': np.double,
                             'AngleTilt': np.double,
                             'AngleTiltPrior': np.double,
                             'AngstromResolution': np.double,
                             'AreaId': 'i4',
                             'AreaName': 'U200',
                             'AutoLocalSearchesHealpixOrder': 'i4',
                             'AutopickFigureOfMerit': np.double,
                             'AvailableMemory': np.double,
                             'AverageNrOfFrames': 'i4',
                             'AveragePmax': np.double,
                             'AverageValue': np.double,
                             'BeamTiltClass': 'i4',
                             'BeamTiltX': np.double,
                             'BeamTiltY': np.double,
                             'BestResolutionThusFar': np.double,
                             'BfactorUsedForSharpening': np.double,
                             'BodyKeepFixed': 'i4',
                             'BodyMaskName': 'U200',
                             'BodyReferenceName': 'U200',
                             'BodyRotateDirectionX': np.double,
                             'BodyRotateDirectionY': np.double,
                             'BodyRotateDirectionZ': np.double,
                             'BodyRotateRelativeTo': 'i4',
                             'BodySigmaAngles': np.double,
                             'BodySigmaOffset': np.double,
                             'BodySigmaPsi': np.double,
                             'BodySigmaRot': np.double,
                             'BodySigmaTilt': np.double,
                             'BodyStarFile': 'U200',
                             'ChangesOptimalClasses': np.double,
                             'ChangesOptimalOffsets': np.double,
                             'ChangesOptimalOrientations': np.double,
                             'ChromaticAberration': np.double,
                             'ClassDistribution': np.double,
                             'ClassNumber': 'i4',
                             'ClassPriorOffsetX': np.double,
                             'ClassPriorOffsetY': np.double,
                             'CoarseImageSize': 'i4',
                             'Comment': 'U200',
                             'ConvergenceCone': np.double,
                             'CoordinateX': np.double,
                             'CoordinateY': np.double,
                             'CoordinateZ': np.double,
                             'CorrectedFourierShellCorrelationPhaseRandomizedMaskedMaps': np.double,
                             'CorrelationFitGuinierPlot': np.double,
                             'CtfAstigmatism': np.double,
                             'CtfBfactor': np.double,
                             'CtfDataAreCtfPremultiplied': bool,
                             'CtfDataArePhaseFlipped': bool,
                             'CtfFigureOfMerit': np.double,
                             'CtfImage': 'U200',
                             'CtfMaxResolution': np.double,
                             'CtfScalefactor': np.double,
                             'CtfValidationScore': np.double,
                             'CtfValue': np.double,
                             'CurrentImageSize': 'i4',
                             'CurrentIteration': 'i4',
                             'CurrentResolution': np.double,
                             'DataDimensionality': 'i4',
                             'DataType': 'i4',
                             'DefocusAngle': np.double,
                             'DefocusU': np.double,
                             'DefocusV': np.double,
                             'DetectorPixelSize': np.double,
                             'Diff2RandomHalves': np.double,
                             'DifferentialPhaseResidualMaskedMaps': np.double,
                             'DifferentialPhaseResidualUnmaskedMaps': np.double,
                             'DoAutoRefine': bool,
                             'DoCorrectCtf': bool,
                             'DoCorrectMagnification': bool,
                             'DoCorrectNorm': bool,
                             'DoCorrectScale': bool,
                             'DoFastSubsetOptimisation': bool,
                             'DoHelicalRefine': bool,
                             'DoIgnoreCtfUntilFirstPeak': bool,
                             'DoMapEstimation': bool,
                             'DoOnlyFlipCtfPhases': bool,
                             'DoRealignMovies': bool,
                             'DoSkipAlign': bool,
                             'DoSkipRotate': bool,
                             'DoSolventFlattening': bool,
                             'DoSolventFscCorrection': bool,
                             'DoSplitRandomHalves': bool,
                             'DoStochasticGradientDescent': bool,
                             'DoZeroMask': bool,
                             'Enabled': bool,
                             'EnergyLoss': np.double,
                             'EstimatedResolution': np.double,
                             'ExperimentalDataStarFile': 'U200',
                             'FinalResolution': np.double,
                             'FittedInterceptGuinierPlot': np.double,
                             'FittedSlopeGuinierPlot': np.double,
                             'FixSigmaNoiseEstimates': bool,
                             'FixSigmaOffsetEstimates': bool,
                             'FixTauEstimates': bool,
                             'FourierCompleteness': np.double,
                             'FourierShellCorrelation': np.double,
                             'FourierShellCorrelationCorrected': np.double,
                             'FourierShellCorrelationMaskedMaps': np.double,
                             'FourierShellCorrelationUnmaskedMaps': np.double,
                             'FourierSpaceInterpolator': 'i4',
                             'GoldStandardFsc': np.double,
                             'GroupName': 'U200',
                             'GroupNrParticles': 'i4',
                             'GroupNumber': 'i4',
                             'GroupScaleCorrection': np.double,
                             'HasConverged': bool,
                             'HasHighFscAtResolLimit': bool,
                             'HasLargeSizeIncreaseIterationsAgo': 'i4',
                             'HealpixOrder': 'i4',
                             'HelicalCentralProportion': np.double,
                             'HelicalKeepTiltPriorFixed': bool,
                             'HelicalMaskTubeInnerDiameter': np.double,
                             'HelicalMaskTubeOuterDiameter': np.double,
                             'HelicalOffsetStep': np.double,
                             'HelicalRise': np.double,
                             'HelicalRiseInitial': np.double,
                             'HelicalRiseInitialStep': np.double,
                             'HelicalRiseMax': np.double,
                             'HelicalRiseMin': np.double,
                             'HelicalSigmaDistance': np.double,
                             'HelicalSymmetryLocalRefinement': bool,
                             'HelicalTrackLength': np.double,
                             'HelicalTubeID': 'i4',
                             'HelicalTubePitch': np.double,
                             'HelicalTwist': np.double,
                             'HelicalTwistInitial': np.double,
                             'HelicalTwistInitialStep': np.double,
                             'HelicalTwistMax': np.double,
                             'HelicalTwistMin': np.double,
                             'HighresLimitExpectation': np.double,
                             'HighresLimitSGD': np.double,
                             'IgnoreHelicalSymmetry': bool,
                             'ImageDimensionality': 'i4',
                             'ImageId': 'i4',
                             'ImageName': 'U200',
                             'ImageOriginalName': 'U200',
                             'ImageSize': 'i4',
                             'ImageSizeX': 'i4',
                             'ImageSizeY': 'i4',
                             'ImageSizeZ': 'i4',
                             'ImageWeight': np.double,
                             'IncrementImageSize': 'i4',
                             'Is3DSampling': bool,
                             'Is3DTranslationalSampling': bool,
                             'IsFlip': bool,
                             'IsHelix': bool,
                             'JoinHalvesUntilThisResolution': np.double,
                             'KullbackLeiblerDivergence': np.double,
                             'KurtosisExcessValue': np.double,
                             'LensStability': np.double,
                             'LocalSymmetryFile': 'U200',
                             'LogAmplitudesIntercept': np.double,
                             'LogAmplitudesMTFCorrected': np.double,
                             'LogAmplitudesOriginal': np.double,
                             'LogAmplitudesSharpened': np.double,
                             'LogAmplitudesWeighted': np.double,
                             'LogLikeliContribution': np.double,
                             'LogLikelihood': np.double,
                             'LongitudinalDisplacement': np.double,
                             'Magnification': np.double,
                             'MagnificationCorrection': np.double,
                             'MagnificationSearchRange': np.double,
                             'MagnificationSearchStep': np.double,
                             'MaskName': 'U200',
                             'Matrix_1_1': np.double,
                             'Matrix_1_2': np.double,
                             'Matrix_1_3': np.double,
                             'Matrix_2_1': np.double,
                             'Matrix_2_2': np.double,
                             'Matrix_2_3': np.double,
                             'Matrix_3_1': np.double,
                             'Matrix_3_2': np.double,
                             'Matrix_3_3': np.double,
                             'MaxNumberOfPooledParticles': 'i4',
                             'MaxValueProbDistribution': np.double,
                             'MaximumCoarseImageSize': 'i4',
                             'MaximumValue': np.double,
                             'MicrographBinning': np.double,
                             'MicrographDefectFile': 'U200',
                             'MicrographDoseRate': np.double,
                             'MicrographEndFrame': 'i4',
                             'MicrographFrameNumber': 'i4',
                             'MicrographGainName': 'U200',
                             'MicrographId': 'i4',
                             'MicrographMetadata': 'U200',
                             'MicrographMovieName': 'U200',
                             'MicrographName': 'U200',
                             'MicrographNameNoDW': 'U200',
                             'MicrographOriginalPixelSize': np.double,
                             'MicrographPreExposure': np.double,
                             'MicrographShiftX': np.double,
                             'MicrographShiftY': np.double,
                             'MicrographStartFrame': 'i4',
                             'MicrographTiltAngle': np.double,
                             'MicrographTiltAxisDirection': np.double,
                             'MicrographTiltAxisOutOfPlane': np.double,
                             'MinRadiusNnInterpolation': 'i4',
                             'MinimumValue': np.double,
                             'ModelStarFile': 'U200',
                             'ModelStarFile2': 'U200',
                             'MotionModelCoeff': np.double,
                             'MotionModelCoeffsIdx': 'i4',
                             'MotionModelVersion': 'i4',
                             'MovieFrameNumber': 'i4',
                             'MovieFramesRunningAverage': 'i4',
                             'MtfValue': np.double,
                             'NormCorrection': np.double,
                             'NormCorrectionAverage': np.double,
                             'NrBodies': 'i4',
                             'NrClasses': 'i4',
                             'NrGroups': 'i4',
                             'NrHelicalAsymUnits': 'i4',
                             'NrOfFrames': 'i4',
                             'NrOfSignificantSamples': 'i4',
                             'NumberOfIterWithoutChangingAssignments': 'i4',
                             'NumberOfIterWithoutResolutionGain': 'i4',
                             'NumberOfIterations': 'i4',
                             'OffsetRange': np.double,
                             'OffsetStep': np.double,
                             'OrientSamplingStarFile': 'U200',
                             'OrientationDistribution': np.double,
                             'OrientationalPriorMode': 'i4',
                             'OrientationsID': 'i4',
                             'OriginX': np.double,
                             'OriginXPrior': np.double,
                             'OriginY': np.double,
                             'OriginYPrior': np.double,
                             'OriginZ': np.double,
                             'OriginZPrior': np.double,
                             'OriginalImageSize': 'i4',
                             'OriginalParticleName': 'U200',
                             'OutputRootName': 'U200',
                             'OverallAccuracyRotations': np.double,
                             'OverallAccuracyTranslations': np.double,
                             'OverallFourierCompleteness': np.double,
                             'PaddingFactor': np.double,
                             'ParticleDiameter': np.double,
                             'ParticleFigureOfMerit': np.double,
                             'ParticleId': 'i4',
                             'ParticleName': 'U200',
                             'ParticleNumber': 'i4',
                             'ParticleSelectZScore': np.double,
                             'PerFrameCumulativeWeight': np.double,
                             'PerFrameRelativeWeight': np.double,
                             'PhaseShift': np.double,
                             'PipeLineEdgeFromNode': 'U200',
                             'PipeLineEdgeProcess': 'U200',
                             'PipeLineEdgeToNode': 'U200',
                             'PipeLineJobCounter': 'i4',
                             'PipeLineNodeName': 'U200',
                             'PipeLineNodeType': 'i4',
                             'PipeLineProcessAlias': 'U200',
                             'PipeLineProcessName': 'U200',
                             'PipeLineProcessStatus': 'i4',
                             'PipeLineProcessType': 'i4',
                             'PixelSize': np.double,
                             'PsiStep': np.double,
                             'RadiusMaskExpImages': 'i4',
                             'RadiusMaskMap': 'i4',
                             'RandomSeed': 'i4',
                             'RandomSubset': 'i4',
                             'RandomiseFrom': np.double,
                             'ReconstructImageName': 'U200',
                             'ReferenceDimensionality': 'i4',
                             'ReferenceImage': 'U200',
                             'ReferenceSigma2': np.double,
                             'ReferenceSpectralPower': np.double,
                             'ReferenceTau2': np.double,
                             'RefsAreCtfCorrected': bool,
                             'Resolution': np.double,
                             'ResolutionInversePixel': np.double,
                             'ResolutionSquared': np.double,
                             'SGDGradientImage': 'U200',
                             'SamplingPerturbFactor': np.double,
                             'SamplingPerturbInstance': np.double,
                             'SamplingRate': np.double,
                             'SamplingRateX': np.double,
                             'SamplingRateY': np.double,
                             'SamplingRateZ': np.double,
                             'Selected': 'i4',
                             'SgdFinalIterations': 'i4',
                             'SgdFinalResolution': np.double,
                             'SgdFinalSubsetSize': 'i4',
                             'SgdInBetweenIterations': 'i4',
                             'SgdInitialIterations': 'i4',
                             'SgdInitialResolution': np.double,
                             'SgdInitialSubsetSize': 'i4',
                             'SgdMaxSubsets': 'i4',
                             'SgdMuFactor': np.double,
                             'SgdSigma2FudgeHalflife': 'i4',
                             'SgdSigma2FudgeInitial': np.double,
                             'SgdSkipAnneal': bool,
                             'SgdStepsize': np.double,
                             'SgdSubsetSize': 'i4',
                             'SgdWriteEverySubset': 'i4',
                             'Sigma2Noise': np.double,
                             'SigmaOffsets': np.double,
                             'SigmaPriorPsiAngle': np.double,
                             'SigmaPriorRotAngle': np.double,
                             'SigmaPriorTiltAngle': np.double,
                             'SignalToNoiseRatio': np.double,
                             'SkewnessValue': np.double,
                             'SmallestChangesClasses': 'i4',
                             'SmallestChangesOffsets': np.double,
                             'SmallestChangesOrientations': np.double,
                             'SolventMask2Name': 'U200',
                             'SolventMaskName': 'U200',
                             'SortedIndex': 'i4',
                             'SpectralIndex': 'i4',
                             'SpectralOrientabilityContribution': np.double,
                             'SphericalAberration': np.double,
                             'SsnrMap': np.double,
                             'StandardDeviationValue': np.double,
                             'StarFileMovieParticles': 'U200',
                             'SymmetryGroup': 'U200',
                             'Tau2FudgeFactor': np.double,
                             'TauSpectrumName': 'U200',
                             'TiltAngleLimit': np.double,
                             'TransversalDisplacement': np.double,
                             'UnfilteredMapHalf1': 'U200',
                             'UnfilteredMapHalf2': 'U200',
                             'UseTooCoarseSampling': bool,
                             'Voltage': np.double,
                             'WidthMaskEdge': 'i4'}

        for relion_options in (relion_31_options, relion_30_options):
            try:
                dtype, dtype_optics_group, size_optics_group = construct_dtype_relion(filename, relion_options)
                break
            except Exception as e:
                print(e)
                pass

    if dtype is None:
        raise Exception('Parsing Relion File Failed: header is not of relion3.0 or relion3.1 format.')

    return loadtxt(filename, dtype=dtype, skip_header=size_optics_group+1)

def read_star_optics_group(filename, order='F', keepnumpy=False, deviceID=None, dtype=None):
    '''Read star file for relion 3.1 or relion3.0 or user defined dtype.'''
    from pytom.basic.files import loadtxt

    size_optics_group = 0

    if dtype is None:

        relion_31_options = {'AccumMotionEarly': np.double,
                             'AccumMotionLate': np.double,
                             'AccumMotionTotal': np.double,
                             'AccuracyRotations': np.double,
                             'AccuracyTranslations': np.double,
                             'AccuracyTranslationsAngst': np.double,
                             'AdaptiveOversampleFraction': np.double,
                             'AdaptiveOversampleOrder': 'i4',
                             'AmplitudeContrast': np.double,
                             'AmplitudeCorrelationMaskedMaps': np.double,
                             'AmplitudeCorrelationUnmaskedMaps': np.double,
                             'AnglePsi': np.double,
                             'AnglePsiFlip': bool,
                             'AnglePsiFlipRatio': np.double,
                             'AnglePsiPrior': np.double,
                             'AngleRot': np.double,
                             'AngleRotFlipRatio': np.double,
                             'AngleRotPrior': np.double,
                             'AngleTilt': np.double,
                             'AngleTiltPrior': np.double,
                             'AngstromResolution': np.double,
                             'AreaId': 'i4',
                             'AreaName': 'U200',
                             'AutoLocalSearchesHealpixOrder': 'i4',
                             'AutopickFigureOfMerit': np.double,
                             'AvailableMemory': np.double,
                             'AverageNrOfFrames': 'i4',
                             'AveragePmax': np.double,
                             'AverageValue': np.double,
                             'BeamTiltClass': 'i4',
                             'BeamTiltX': np.double,
                             'BeamTiltY': np.double,
                             'BestResolutionThusFar': np.double,
                             'BfactorUsedForSharpening': np.double,
                             'BodyKeepFixed': 'i4',
                             'BodyMaskName': 'U200',
                             'BodyReferenceName': 'U200',
                             'BodyRotateDirectionX': np.double,
                             'BodyRotateDirectionY': np.double,
                             'BodyRotateDirectionZ': np.double,
                             'BodyRotateRelativeTo': 'i4',
                             'BodySigmaAngles': np.double,
                             'BodySigmaOffset': np.double,
                             'BodySigmaOffsetAngst': np.double,
                             'BodySigmaPsi': np.double,
                             'BodySigmaRot': np.double,
                             'BodySigmaTilt': np.double,
                             'BodyStarFile': 'U200',
                             'ChangesOptimalClasses': np.double,
                             'ChangesOptimalOffsets': np.double,
                             'ChangesOptimalOrientations': np.double,
                             'ChromaticAberration': np.double,
                             'ClassDistribution': np.double,
                             'ClassNumber': 'i4',
                             'ClassPriorOffsetX': np.double,
                             'ClassPriorOffsetY': np.double,
                             'CoarseImageSize': 'i4',
                             'Comment': 'U200',
                             'ConvergenceCone': np.double,
                             'CoordinateX': np.double,
                             'CoordinateY': np.double,
                             'CoordinateZ': np.double,
                             'CorrectedFourierShellCorrelationPhaseRandomizedMaskedMaps': np.double,
                             'CorrelationFitGuinierPlot': np.double,
                             'CtfAstigmatism': np.double,
                             'CtfBfactor': np.double,
                             'CtfDataAreCtfPremultiplied': bool,
                             'CtfDataArePhaseFlipped': bool,
                             'CtfFigureOfMerit': np.double,
                             'CtfImage': 'U200',
                             'CtfMaxResolution': np.double,
                             'CtfPowerSpectrum': 'U200',
                             'CtfScalefactor': np.double,
                             'CtfValidationScore': np.double,
                             'CtfValue': np.double,
                             'CurrentImageSize': 'i4',
                             'CurrentIteration': 'i4',
                             'CurrentResolution': np.double,
                             'DataDimensionality': 'i4',
                             'DataType': 'i4',
                             'DefocusAngle': np.double,
                             'DefocusU': np.double,
                             'DefocusV': np.double,
                             'DetectorPixelSize': np.double,
                             'Diff2RandomHalves': np.double,
                             'DifferentialPhaseResidualMaskedMaps': np.double,
                             'DifferentialPhaseResidualUnmaskedMaps': np.double,
                             'DoAutoRefine': bool,
                             'DoCorrectCtf': bool,
                             'DoCorrectMagnification': bool,
                             'DoCorrectNorm': bool,
                             'DoCorrectScale': bool,
                             'DoExternalReconstruct': bool,
                             'DoFastSubsetOptimisation': bool,
                             'DoHelicalRefine': bool,
                             'DoIgnoreCtfUntilFirstPeak': bool,
                             'DoMapEstimation': bool,
                             'DoOnlyFlipCtfPhases': bool,
                             'DoRealignMovies': bool,
                             'DoSkipAlign': bool,
                             'DoSkipRotate': bool,
                             'DoSolventFlattening': bool,
                             'DoSolventFscCorrection': bool,
                             'DoSplitRandomHalves': bool,
                             'DoStochasticEM': bool,
                             'DoStochasticGradientDescent': bool,
                             'DoZeroMask': bool,
                             'EERGrouping': 'i4',
                             'EERUpsampling': 'i4',
                             'Enabled': bool,
                             'EnergyLoss': np.double,
                             'EstimatedResolution': np.double,
                             'EvenZernike': np.double,
                             'ExperimentalDataStarFile': 'U200',
                             'ExtReconsDataImag': 'U200',
                             'ExtReconsDataReal': 'U200',
                             'ExtReconsResult': 'U200',
                             'ExtReconsResultStarfile': 'U200',
                             'ExtReconsWeight': 'U200',
                             'FinalResolution': np.double,
                             'FittedInterceptGuinierPlot': np.double,
                             'FittedSlopeGuinierPlot': np.double,
                             'FixSigmaNoiseEstimates': bool,
                             'FixSigmaOffsetEstimates': bool,
                             'FixTauEstimates': bool,
                             'FourierCompleteness': np.double,
                             'FourierMask': 'U200',
                             'FourierShellCorrelation': np.double,
                             'FourierShellCorrelationCorrected': np.double,
                             'FourierShellCorrelationMaskedMaps': np.double,
                             'FourierShellCorrelationParticleMaskFraction': np.double,
                             'FourierShellCorrelationParticleMolWeight': np.double,
                             'FourierShellCorrelationUnmaskedMaps': np.double,
                             'FourierSpaceInterpolator': 'i4',
                             'GoldStandardFsc': np.double,
                             'GroupName': 'U200',
                             'GroupNrParticles': 'i4',
                             'GroupNumber': 'i4',
                             'GroupScaleCorrection': np.double,
                             'HasConverged': bool,
                             'HasHighFscAtResolLimit': bool,
                             'HasLargeSizeIncreaseIterationsAgo': 'i4',
                             'HealpixOrder': 'i4',
                             'HealpixOrderOriginal': 'i4',
                             'HelicalCentralProportion': np.double,
                             'HelicalKeepTiltPriorFixed': bool,
                             'HelicalMaskTubeInnerDiameter': np.double,
                             'HelicalMaskTubeOuterDiameter': np.double,
                             'HelicalOffsetStep': np.double,
                             'HelicalRise': np.double,
                             'HelicalRiseInitial': np.double,
                             'HelicalRiseInitialStep': np.double,
                             'HelicalRiseMax': np.double,
                             'HelicalRiseMin': np.double,
                             'HelicalSigmaDistance': np.double,
                             'HelicalSymmetryLocalRefinement': bool,
                             'HelicalTrackLength': np.double,
                             'HelicalTrackLengthAngst': np.double,
                             'HelicalTubeID': 'i4',
                             'HelicalTubePitch': np.double,
                             'HelicalTwist': np.double,
                             'HelicalTwistInitial': np.double,
                             'HelicalTwistInitialStep': np.double,
                             'HelicalTwistMax': np.double,
                             'HelicalTwistMin': np.double,
                             'HighresLimitExpectation': np.double,
                             'HighresLimitSGD': np.double,
                             'IgnoreHelicalSymmetry': bool,
                             'ImageDimensionality': 'i4',
                             'ImageId': 'i4',
                             'ImageName': 'U200',
                             'ImageOriginalName': 'U200',
                             'ImagePixelSize': np.double,
                             'ImageSize': 'i4',
                             'ImageSizeX': 'i4',
                             'ImageSizeY': 'i4',
                             'ImageSizeZ': 'i4',
                             'ImageWeight': np.double,
                             'IncrementImageSize': 'i4',
                             'Is3DSampling': bool,
                             'Is3DTranslationalSampling': bool,
                             'IsFlip': bool,
                             'IsHelix': bool,
                             'JobIsContinue': bool,
                             'JobOptionDefaultValue': 'U200',
                             'JobOptionDirectoryDefault': 'U200',
                             'JobOptionFilePattern': 'U200',
                             'JobOptionGUILabel': 'U200',
                             'JobOptionHelpText': 'U200',
                             'JobOptionMenuOptions': 'U200',
                             'JobOptionSliderMax': np.double,
                             'JobOptionSliderMin': np.double,
                             'JobOptionSliderStep': np.double,
                             'JobOptionValue': 'U200',
                             'JobOptionVariable': 'U200',
                             'JobType': 'i4',
                             'JobTypeName': 'U200',
                             'JoboptionType': 'i4',
                             'JoinHalvesUntilThisResolution': np.double,
                             'KullbackLeiblerDivergence': np.double,
                             'KurtosisExcessValue': np.double,
                             'LensStability': np.double,
                             'LocalSymmetryFile': 'U200',
                             'LogAmplitudesIntercept': np.double,
                             'LogAmplitudesMTFCorrected': np.double,
                             'LogAmplitudesOriginal': np.double,
                             'LogAmplitudesSharpened': np.double,
                             'LogAmplitudesWeighted': np.double,
                             'LogLikeliContribution': np.double,
                             'LogLikelihood': np.double,
                             'LongitudinalDisplacement': np.double,
                             'LowresLimitExpectation': np.double,
                             'MagMat00': np.double,
                             'MagMat01': np.double,
                             'MagMat10': np.double,
                             'MagMat11': np.double,
                             'Magnification': np.double,
                             'MagnificationCorrection': np.double,
                             'MagnificationSearchRange': np.double,
                             'MagnificationSearchStep': np.double,
                             'MaskName': 'U200',
                             'Matrix_1_1': np.double,
                             'Matrix_1_2': np.double,
                             'Matrix_1_3': np.double,
                             'Matrix_2_1': np.double,
                             'Matrix_2_2': np.double,
                             'Matrix_2_3': np.double,
                             'Matrix_3_1': np.double,
                             'Matrix_3_2': np.double,
                             'Matrix_3_3': np.double,
                             'MaxNumberOfPooledParticles': 'i4',
                             'MaxValueProbDistribution': np.double,
                             'MaximumCoarseImageSize': 'i4',
                             'MaximumValue': np.double,
                             'MicrographBinning': np.double,
                             'MicrographDefectFile': 'U200',
                             'MicrographDoseRate': np.double,
                             'MicrographEndFrame': 'i4',
                             'MicrographFrameNumber': 'i4',
                             'MicrographGainName': 'U200',
                             'MicrographId': 'i4',
                             'MicrographMetadata': 'U200',
                             'MicrographMovieName': 'U200',
                             'MicrographName': 'U200',
                             'MicrographNameNoDW': 'U200',
                             'MicrographOriginalPixelSize': np.double,
                             'MicrographPixelSize': np.double,
                             'MicrographPreExposure': np.double,
                             'MicrographShiftX': np.double,
                             'MicrographShiftY': np.double,
                             'MicrographStartFrame': 'i4',
                             'MicrographTiltAngle': np.double,
                             'MicrographTiltAxisDirection': np.double,
                             'MicrographTiltAxisOutOfPlane': np.double,
                             'MinRadiusNnInterpolation': 'i4',
                             'MinimumValue': np.double,
                             'ModelStarFile': 'U200',
                             'ModelStarFile2': 'U200',
                             'MolecularWeight': np.double,
                             'MotionModelCoeff': np.double,
                             'MotionModelCoeffsIdx': 'i4',
                             'MotionModelVersion': 'i4',
                             'MovieFrameNumber': 'i4',
                             'MovieFramesRunningAverage': 'i4',
                             'MtfFileName': 'U200',
                             'MtfValue': np.double,
                             'NormCorrection': np.double,
                             'NormCorrectionAverage': np.double,
                             'NrBodies': 'i4',
                             'NrClasses': 'i4',
                             'NrGroups': 'i4',
                             'NrHelicalAsymUnits': 'i4',
                             'NrHelicalNStart': 'i4',
                             'NrOfFrames': 'i4',
                             'NrOfSignificantSamples': 'i4',
                             'NumberOfIterWithoutChangingAssignments': 'i4',
                             'NumberOfIterWithoutResolutionGain': 'i4',
                             'NumberOfIterations': 'i4',
                             'OddZernike': np.double,
                             'OffsetRange': np.double,
                             'OffsetRangeOriginal': np.double,
                             'OffsetStep': np.double,
                             'OffsetStepOriginal': np.double,
                             'OpticsGroup': 'i4',
                             'OpticsGroupName': 'U200',
                             'OpticsStarFile': 'U200',
                             'OrientSamplingStarFile': 'U200',
                             'OrientationDistribution': np.double,
                             'OrientationalPriorMode': 'i4',
                             'OrientationsID': 'i4',
                             'OriginX': np.double,
                             'OriginXAngst': np.double,
                             'OriginXPrior': np.double,
                             'OriginXPriorAngst': np.double,
                             'OriginY': np.double,
                             'OriginYAngst': np.double,
                             'OriginYPrior': np.double,
                             'OriginYPriorAngst': np.double,
                             'OriginZ': np.double,
                             'OriginZAngst': np.double,
                             'OriginZPrior': np.double,
                             'OriginZPriorAngst': np.double,
                             'OriginalImageSize': 'i4',
                             'OriginalParticleName': 'U200',
                             'OutputRootName': 'U200',
                             'OverallAccuracyRotations': np.double,
                             'OverallAccuracyTranslations': np.double,
                             'OverallAccuracyTranslationsAngst': np.double,
                             'OverallFourierCompleteness': np.double,
                             'PaddingFactor': np.double,
                             'ParticleBoxFractionMolecularWeight': np.double,
                             'ParticleBoxFractionSolventMask': np.double,
                             'ParticleDiameter': np.double,
                             'ParticleFigureOfMerit': np.double,
                             'ParticleId': 'i4',
                             'ParticleName': 'U200',
                             'ParticleNumber': 'i4',
                             'ParticleSelectZScore': np.double,
                             'PerFrameCumulativeWeight': np.double,
                             'PerFrameRelativeWeight': np.double,
                             'PhaseShift': np.double,
                             'PipeLineEdgeFromNode': 'U200',
                             'PipeLineEdgeProcess': 'U200',
                             'PipeLineEdgeToNode': 'U200',
                             'PipeLineJobCounter': 'i4',
                             'PipeLineNodeName': 'U200',
                             'PipeLineNodeType': 'i4',
                             'PipeLineProcessAlias': 'U200',
                             'PipeLineProcessName': 'U200',
                             'PipeLineProcessStatus': 'i4',
                             'PipeLineProcessType': 'i4',
                             'PixelSize': np.double,
                             'PsiStep': np.double,
                             'PsiStepOriginal': np.double,
                             'RadiusMaskExpImages': 'i4',
                             'RadiusMaskMap': 'i4',
                             'RandomSeed': 'i4',
                             'RandomSubset': 'i4',
                             'RandomiseFrom': np.double,
                             'ReconstructImageName': 'U200',
                             'ReferenceDimensionality': 'i4',
                             'ReferenceImage': 'U200',
                             'ReferenceSigma2': np.double,
                             'ReferenceSpectralPower': np.double,
                             'ReferenceTau2': np.double,
                             'RefsAreCtfCorrected': bool,
                             'Resolution': np.double,
                             'ResolutionInversePixel': np.double,
                             'ResolutionSquared': np.double,
                             'SGDGradientImage': 'U200',
                             'SamplingPerturbFactor': np.double,
                             'SamplingPerturbInstance': np.double,
                             'SamplingRate': np.double,
                             'SamplingRateX': np.double,
                             'SamplingRateY': np.double,
                             'SamplingRateZ': np.double,
                             'ScheduleBooleanVariableName': 'U200',
                             'ScheduleBooleanVariableResetValue': bool,
                             'ScheduleBooleanVariableValue': bool,
                             'ScheduleCurrentNodeName': 'U200',
                             'ScheduleEdgeBooleanVariable': 'U200',
                             'ScheduleEdgeInputNodeName': 'U200',
                             'ScheduleEdgeIsFork': bool,
                             'ScheduleEdgeNumber': 'i4',
                             'ScheduleEdgeOutputNodeName': 'U200',
                             'ScheduleEdgeOutputNodeNameIfTrue': 'U200',
                             'ScheduleEmailAddress': 'U200',
                             'ScheduleFloatVariableName': 'U200',
                             'ScheduleFloatVariableResetValue': np.double,
                             'ScheduleFloatVariableValue': np.double,
                             'ScheduleJobHasStarted': bool,
                             'ScheduleJobMode': 'U200',
                             'ScheduleJobName': 'U200',
                             'ScheduleJobNameOriginal': 'U200',
                             'ScheduleName': 'U200',
                             'ScheduleOperatorInput1': 'U200',
                             'ScheduleOperatorInput2': 'U200',
                             'ScheduleOperatorName': 'U200',
                             'ScheduleOperatorOutput': 'U200',
                             'ScheduleOperatorType': 'U200',
                             'ScheduleOriginalStartNodeName': 'U200',
                             'ScheduleStringVariableName': 'U200',
                             'ScheduleStringVariableResetValue': 'U200',
                             'ScheduleStringVariableValue': 'U200',
                             'Selected': 'i4',
                             'SgdFinalIterations': 'i4',
                             'SgdFinalResolution': np.double,
                             'SgdFinalSubsetSize': 'i4',
                             'SgdInBetweenIterations': 'i4',
                             'SgdInitialIterations': 'i4',
                             'SgdInitialResolution': np.double,
                             'SgdInitialSubsetSize': 'i4',
                             'SgdMaxSubsets': 'i4',
                             'SgdMuFactor': np.double,
                             'SgdSigma2FudgeHalflife': 'i4',
                             'SgdSigma2FudgeInitial': np.double,
                             'SgdSkipAnneal': bool,
                             'SgdStepsize': np.double,
                             'SgdSubsetSize': 'i4',
                             'SgdWriteEverySubset': 'i4',
                             'Sigma2Noise': np.double,
                             'SigmaOffsets': np.double,
                             'SigmaOffsetsAngst': np.double,
                             'SigmaPriorPsiAngle': np.double,
                             'SigmaPriorRotAngle': np.double,
                             'SigmaPriorTiltAngle': np.double,
                             'SignalToNoiseRatio': np.double,
                             'SkewnessValue': np.double,
                             'SmallestChangesClasses': 'i4',
                             'SmallestChangesOffsets': np.double,
                             'SmallestChangesOrientations': np.double,
                             'SolventMask2Name': 'U200',
                             'SolventMaskName': 'U200',
                             'SortedIndex': 'i4',
                             'SpectralIndex': 'i4',
                             'SpectralOrientabilityContribution': np.double,
                             'SphericalAberration': np.double,
                             'SsnrMap': np.double,
                             'StandardDeviationValue': np.double,
                             'StarFileMovieParticles': 'U200',
                             'SymmetryGroup': 'U200',
                             'Tau2FudgeFactor': np.double,
                             'TauSpectrumName': 'U200',
                             'TiltAngleLimit': np.double,
                             'TransversalDisplacement': np.double,
                             'UnfilteredMapHalf1': 'U200',
                             'UnfilteredMapHalf2': 'U200',
                             'UnknownLabel': 'U200',
                             'UseTooCoarseSampling': bool,
                             'Voltage': np.double,
                             'WidthMaskEdge': 'i4'}

        relion_30_options = {'AccumMotionEarly': np.double,
                             'AccumMotionLate': np.double,
                             'AccumMotionTotal': np.double,
                             'AccuracyRotations': np.double,
                             'AccuracyTranslations': np.double,
                             'AdaptiveOversampleFraction': np.double,
                             'AdaptiveOversampleOrder': 'i4',
                             'AmplitudeContrast': np.double,
                             'AmplitudeCorrelationMaskedMaps': np.double,
                             'AmplitudeCorrelationUnmaskedMaps': np.double,
                             'AnglePsi': np.double,
                             'AnglePsiFlipRatio': np.double,
                             'AnglePsiPrior': np.double,
                             'AngleRot': np.double,
                             'AngleRotPrior': np.double,
                             'AngleTilt': np.double,
                             'AngleTiltPrior': np.double,
                             'AngstromResolution': np.double,
                             'AreaId': 'i4',
                             'AreaName': 'U200',
                             'AutoLocalSearchesHealpixOrder': 'i4',
                             'AutopickFigureOfMerit': np.double,
                             'AvailableMemory': np.double,
                             'AverageNrOfFrames': 'i4',
                             'AveragePmax': np.double,
                             'AverageValue': np.double,
                             'BeamTiltClass': 'i4',
                             'BeamTiltX': np.double,
                             'BeamTiltY': np.double,
                             'BestResolutionThusFar': np.double,
                             'BfactorUsedForSharpening': np.double,
                             'BodyKeepFixed': 'i4',
                             'BodyMaskName': 'U200',
                             'BodyReferenceName': 'U200',
                             'BodyRotateDirectionX': np.double,
                             'BodyRotateDirectionY': np.double,
                             'BodyRotateDirectionZ': np.double,
                             'BodyRotateRelativeTo': 'i4',
                             'BodySigmaAngles': np.double,
                             'BodySigmaOffset': np.double,
                             'BodySigmaPsi': np.double,
                             'BodySigmaRot': np.double,
                             'BodySigmaTilt': np.double,
                             'BodyStarFile': 'U200',
                             'ChangesOptimalClasses': np.double,
                             'ChangesOptimalOffsets': np.double,
                             'ChangesOptimalOrientations': np.double,
                             'ChromaticAberration': np.double,
                             'ClassDistribution': np.double,
                             'ClassNumber': 'i4',
                             'ClassPriorOffsetX': np.double,
                             'ClassPriorOffsetY': np.double,
                             'CoarseImageSize': 'i4',
                             'Comment': 'U200',
                             'ConvergenceCone': np.double,
                             'CoordinateX': np.double,
                             'CoordinateY': np.double,
                             'CoordinateZ': np.double,
                             'CorrectedFourierShellCorrelationPhaseRandomizedMaskedMaps': np.double,
                             'CorrelationFitGuinierPlot': np.double,
                             'CtfAstigmatism': np.double,
                             'CtfBfactor': np.double,
                             'CtfDataAreCtfPremultiplied': bool,
                             'CtfDataArePhaseFlipped': bool,
                             'CtfFigureOfMerit': np.double,
                             'CtfImage': 'U200',
                             'CtfMaxResolution': np.double,
                             'CtfScalefactor': np.double,
                             'CtfValidationScore': np.double,
                             'CtfValue': np.double,
                             'CurrentImageSize': 'i4',
                             'CurrentIteration': 'i4',
                             'CurrentResolution': np.double,
                             'DataDimensionality': 'i4',
                             'DataType': 'i4',
                             'DefocusAngle': np.double,
                             'DefocusU': np.double,
                             'DefocusV': np.double,
                             'DetectorPixelSize': np.double,
                             'Diff2RandomHalves': np.double,
                             'DifferentialPhaseResidualMaskedMaps': np.double,
                             'DifferentialPhaseResidualUnmaskedMaps': np.double,
                             'DoAutoRefine': bool,
                             'DoCorrectCtf': bool,
                             'DoCorrectMagnification': bool,
                             'DoCorrectNorm': bool,
                             'DoCorrectScale': bool,
                             'DoFastSubsetOptimisation': bool,
                             'DoHelicalRefine': bool,
                             'DoIgnoreCtfUntilFirstPeak': bool,
                             'DoMapEstimation': bool,
                             'DoOnlyFlipCtfPhases': bool,
                             'DoRealignMovies': bool,
                             'DoSkipAlign': bool,
                             'DoSkipRotate': bool,
                             'DoSolventFlattening': bool,
                             'DoSolventFscCorrection': bool,
                             'DoSplitRandomHalves': bool,
                             'DoStochasticGradientDescent': bool,
                             'DoZeroMask': bool,
                             'Enabled': bool,
                             'EnergyLoss': np.double,
                             'EstimatedResolution': np.double,
                             'ExperimentalDataStarFile': 'U200',
                             'FinalResolution': np.double,
                             'FittedInterceptGuinierPlot': np.double,
                             'FittedSlopeGuinierPlot': np.double,
                             'FixSigmaNoiseEstimates': bool,
                             'FixSigmaOffsetEstimates': bool,
                             'FixTauEstimates': bool,
                             'FourierCompleteness': np.double,
                             'FourierShellCorrelation': np.double,
                             'FourierShellCorrelationCorrected': np.double,
                             'FourierShellCorrelationMaskedMaps': np.double,
                             'FourierShellCorrelationUnmaskedMaps': np.double,
                             'FourierSpaceInterpolator': 'i4',
                             'GoldStandardFsc': np.double,
                             'GroupName': 'U200',
                             'GroupNrParticles': 'i4',
                             'GroupNumber': 'i4',
                             'GroupScaleCorrection': np.double,
                             'HasConverged': bool,
                             'HasHighFscAtResolLimit': bool,
                             'HasLargeSizeIncreaseIterationsAgo': 'i4',
                             'HealpixOrder': 'i4',
                             'HelicalCentralProportion': np.double,
                             'HelicalKeepTiltPriorFixed': bool,
                             'HelicalMaskTubeInnerDiameter': np.double,
                             'HelicalMaskTubeOuterDiameter': np.double,
                             'HelicalOffsetStep': np.double,
                             'HelicalRise': np.double,
                             'HelicalRiseInitial': np.double,
                             'HelicalRiseInitialStep': np.double,
                             'HelicalRiseMax': np.double,
                             'HelicalRiseMin': np.double,
                             'HelicalSigmaDistance': np.double,
                             'HelicalSymmetryLocalRefinement': bool,
                             'HelicalTrackLength': np.double,
                             'HelicalTubeID': 'i4',
                             'HelicalTubePitch': np.double,
                             'HelicalTwist': np.double,
                             'HelicalTwistInitial': np.double,
                             'HelicalTwistInitialStep': np.double,
                             'HelicalTwistMax': np.double,
                             'HelicalTwistMin': np.double,
                             'HighresLimitExpectation': np.double,
                             'HighresLimitSGD': np.double,
                             'IgnoreHelicalSymmetry': bool,
                             'ImageDimensionality': 'i4',
                             'ImageId': 'i4',
                             'ImageName': 'U200',
                             'ImageOriginalName': 'U200',
                             'ImageSize': 'i4',
                             'ImageSizeX': 'i4',
                             'ImageSizeY': 'i4',
                             'ImageSizeZ': 'i4',
                             'ImageWeight': np.double,
                             'IncrementImageSize': 'i4',
                             'Is3DSampling': bool,
                             'Is3DTranslationalSampling': bool,
                             'IsFlip': bool,
                             'IsHelix': bool,
                             'JoinHalvesUntilThisResolution': np.double,
                             'KullbackLeiblerDivergence': np.double,
                             'KurtosisExcessValue': np.double,
                             'LensStability': np.double,
                             'LocalSymmetryFile': 'U200',
                             'LogAmplitudesIntercept': np.double,
                             'LogAmplitudesMTFCorrected': np.double,
                             'LogAmplitudesOriginal': np.double,
                             'LogAmplitudesSharpened': np.double,
                             'LogAmplitudesWeighted': np.double,
                             'LogLikeliContribution': np.double,
                             'LogLikelihood': np.double,
                             'LongitudinalDisplacement': np.double,
                             'Magnification': np.double,
                             'MagnificationCorrection': np.double,
                             'MagnificationSearchRange': np.double,
                             'MagnificationSearchStep': np.double,
                             'MaskName': 'U200',
                             'Matrix_1_1': np.double,
                             'Matrix_1_2': np.double,
                             'Matrix_1_3': np.double,
                             'Matrix_2_1': np.double,
                             'Matrix_2_2': np.double,
                             'Matrix_2_3': np.double,
                             'Matrix_3_1': np.double,
                             'Matrix_3_2': np.double,
                             'Matrix_3_3': np.double,
                             'MaxNumberOfPooledParticles': 'i4',
                             'MaxValueProbDistribution': np.double,
                             'MaximumCoarseImageSize': 'i4',
                             'MaximumValue': np.double,
                             'MicrographBinning': np.double,
                             'MicrographDefectFile': 'U200',
                             'MicrographDoseRate': np.double,
                             'MicrographEndFrame': 'i4',
                             'MicrographFrameNumber': 'i4',
                             'MicrographGainName': 'U200',
                             'MicrographId': 'i4',
                             'MicrographMetadata': 'U200',
                             'MicrographMovieName': 'U200',
                             'MicrographName': 'U200',
                             'MicrographNameNoDW': 'U200',
                             'MicrographOriginalPixelSize': np.double,
                             'MicrographPreExposure': np.double,
                             'MicrographShiftX': np.double,
                             'MicrographShiftY': np.double,
                             'MicrographStartFrame': 'i4',
                             'MicrographTiltAngle': np.double,
                             'MicrographTiltAxisDirection': np.double,
                             'MicrographTiltAxisOutOfPlane': np.double,
                             'MinRadiusNnInterpolation': 'i4',
                             'MinimumValue': np.double,
                             'ModelStarFile': 'U200',
                             'ModelStarFile2': 'U200',
                             'MotionModelCoeff': np.double,
                             'MotionModelCoeffsIdx': 'i4',
                             'MotionModelVersion': 'i4',
                             'MovieFrameNumber': 'i4',
                             'MovieFramesRunningAverage': 'i4',
                             'MtfValue': np.double,
                             'NormCorrection': np.double,
                             'NormCorrectionAverage': np.double,
                             'NrBodies': 'i4',
                             'NrClasses': 'i4',
                             'NrGroups': 'i4',
                             'NrHelicalAsymUnits': 'i4',
                             'NrOfFrames': 'i4',
                             'NrOfSignificantSamples': 'i4',
                             'NumberOfIterWithoutChangingAssignments': 'i4',
                             'NumberOfIterWithoutResolutionGain': 'i4',
                             'NumberOfIterations': 'i4',
                             'OffsetRange': np.double,
                             'OffsetStep': np.double,
                             'OrientSamplingStarFile': 'U200',
                             'OrientationDistribution': np.double,
                             'OrientationalPriorMode': 'i4',
                             'OrientationsID': 'i4',
                             'OriginX': np.double,
                             'OriginXPrior': np.double,
                             'OriginY': np.double,
                             'OriginYPrior': np.double,
                             'OriginZ': np.double,
                             'OriginZPrior': np.double,
                             'OriginalImageSize': 'i4',
                             'OriginalParticleName': 'U200',
                             'OutputRootName': 'U200',
                             'OverallAccuracyRotations': np.double,
                             'OverallAccuracyTranslations': np.double,
                             'OverallFourierCompleteness': np.double,
                             'PaddingFactor': np.double,
                             'ParticleDiameter': np.double,
                             'ParticleFigureOfMerit': np.double,
                             'ParticleId': 'i4',
                             'ParticleName': 'U200',
                             'ParticleNumber': 'i4',
                             'ParticleSelectZScore': np.double,
                             'PerFrameCumulativeWeight': np.double,
                             'PerFrameRelativeWeight': np.double,
                             'PhaseShift': np.double,
                             'PipeLineEdgeFromNode': 'U200',
                             'PipeLineEdgeProcess': 'U200',
                             'PipeLineEdgeToNode': 'U200',
                             'PipeLineJobCounter': 'i4',
                             'PipeLineNodeName': 'U200',
                             'PipeLineNodeType': 'i4',
                             'PipeLineProcessAlias': 'U200',
                             'PipeLineProcessName': 'U200',
                             'PipeLineProcessStatus': 'i4',
                             'PipeLineProcessType': 'i4',
                             'PixelSize': np.double,
                             'PsiStep': np.double,
                             'RadiusMaskExpImages': 'i4',
                             'RadiusMaskMap': 'i4',
                             'RandomSeed': 'i4',
                             'RandomSubset': 'i4',
                             'RandomiseFrom': np.double,
                             'ReconstructImageName': 'U200',
                             'ReferenceDimensionality': 'i4',
                             'ReferenceImage': 'U200',
                             'ReferenceSigma2': np.double,
                             'ReferenceSpectralPower': np.double,
                             'ReferenceTau2': np.double,
                             'RefsAreCtfCorrected': bool,
                             'Resolution': np.double,
                             'ResolutionInversePixel': np.double,
                             'ResolutionSquared': np.double,
                             'SGDGradientImage': 'U200',
                             'SamplingPerturbFactor': np.double,
                             'SamplingPerturbInstance': np.double,
                             'SamplingRate': np.double,
                             'SamplingRateX': np.double,
                             'SamplingRateY': np.double,
                             'SamplingRateZ': np.double,
                             'Selected': 'i4',
                             'SgdFinalIterations': 'i4',
                             'SgdFinalResolution': np.double,
                             'SgdFinalSubsetSize': 'i4',
                             'SgdInBetweenIterations': 'i4',
                             'SgdInitialIterations': 'i4',
                             'SgdInitialResolution': np.double,
                             'SgdInitialSubsetSize': 'i4',
                             'SgdMaxSubsets': 'i4',
                             'SgdMuFactor': np.double,
                             'SgdSigma2FudgeHalflife': 'i4',
                             'SgdSigma2FudgeInitial': np.double,
                             'SgdSkipAnneal': bool,
                             'SgdStepsize': np.double,
                             'SgdSubsetSize': 'i4',
                             'SgdWriteEverySubset': 'i4',
                             'Sigma2Noise': np.double,
                             'SigmaOffsets': np.double,
                             'SigmaPriorPsiAngle': np.double,
                             'SigmaPriorRotAngle': np.double,
                             'SigmaPriorTiltAngle': np.double,
                             'SignalToNoiseRatio': np.double,
                             'SkewnessValue': np.double,
                             'SmallestChangesClasses': 'i4',
                             'SmallestChangesOffsets': np.double,
                             'SmallestChangesOrientations': np.double,
                             'SolventMask2Name': 'U200',
                             'SolventMaskName': 'U200',
                             'SortedIndex': 'i4',
                             'SpectralIndex': 'i4',
                             'SpectralOrientabilityContribution': np.double,
                             'SphericalAberration': np.double,
                             'SsnrMap': np.double,
                             'StandardDeviationValue': np.double,
                             'StarFileMovieParticles': 'U200',
                             'SymmetryGroup': 'U200',
                             'Tau2FudgeFactor': np.double,
                             'TauSpectrumName': 'U200',
                             'TiltAngleLimit': np.double,
                             'TransversalDisplacement': np.double,
                             'UnfilteredMapHalf1': 'U200',
                             'UnfilteredMapHalf2': 'U200',
                             'UseTooCoarseSampling': bool,
                             'Voltage': np.double,
                             'WidthMaskEdge': 'i4'}

        for relion_options in (relion_31_options, relion_30_options):
            try:
                dtype, dtype_optics_group, size_optics_group = construct_dtype_relion(filename, relion_options)
                break
            except Exception as e:
                pass

    if dtype is None:
        raise Exception('Parsing Relion File Failed: header is not of relion3.0 or relion3.1 format.')

    return loadtxt(filename, dtype=dtype_optics_group, max_rows=size_optics_group)

def read_txt(filename, order='F', keepnumpy=False, deviceID=None, dtype=None):
    '''Read a text based file, could be a star file, a pytom support file, or general txt file that can be parsed by loadtxt'''
    from pytom.basic.datatypes import datatype_list
    from pytom.gui.guiFunctions import datatypeAR, datatype, datatypeMR, datatype0, ALIGNMENT_ERRORS
    from pytom.basic.files import loadtxt

    datatype_list += [datatypeAR, datatype, datatypeMR, datatype0, ALIGNMENT_ERRORS, dtype]
    success = False
    data = None
    columnIDS = []
    skip_header=0

    info = open(filename)

    try:
        for line in info.readlines():
            if line.startswith('#') and len(line) > 4 and not 'loop' in line and not 'data' in line:
                if '_' in line:
                    columnIDS.append(line.split('_')[1].split()[0])
                elif '#' in line:
                    columnIDS.append(line.split('#')[1].split()[0])
                else:
                    columnIDS.append(line.split(' ')[1].split()[0])
            if line.startswith(' view '):
                for i in line.split():
                    columnIDS.append(i.capitalize())

                columnIDS = columnIDS[:-2] + ['MeanResid']

                skip_header = 3

    except Exception as e:
        print(e)
        raise(f'PyTom cannot interpret header of {filename}.')

    # We assume it is either a txt file with pytom_defined format or easy to read with loadtxt (dtype=None)
    # Lets try all pytom-defined datatypes.

    info.close()


    for DTYPE in datatype_list:

        try:
            data = loadtxt(filename, dtype=DTYPE, skip_header=skip_header)
            for name in data.dtype.names:
                if not name in columnIDS:
                    raise Exception(f'not right {name}')

            success = True
        except Exception as e:
            pass
        if success:
            return data

    if not success: raise Exception(f'No valid interpreter found for {filename}')


def read_mrc(filename, order='F', keepnumpy=False, deviceID=None, dtype=None):
    import numpy as np
    f = open(filename, 'r')
    try:
        dt_header = np.dtype('int32')
        header = np.fromfile(f, dt_header, 256)
        x = header[0]
        y = header[1]
        z = header[2]

        if header[23]:
            extended_header = np.fromfile(f, np.float32, header[23] // 4)

        # read the data type
        dt = header[3]
        default_type = False
        if dt == 0:  # byte
            dt_data = np.dtype('int8')
        elif dt == 1:  # short
            dt_data = np.dtype('int16')
        elif dt == 2:  # long
            dt_data = np.dtype('float32')  # little-endian, float32
            default_type = True
        elif dt == 4:  # float32
            dt_data = np.dtype('complex32')
        elif dt == 6:  # float complex
            dt_data = np.dtype('uint16')
        elif dt == 12:
            dt_data = np.dtype('float16')
        else:
            raise Exception("Data type not supported yet!")

        v = np.fromfile(f, dt_data, x * y * z)
    finally:
        f.close()

    if keepnumpy:
        volume = np.array(v.reshape((x, y, z), order=order), dtype='float32').copy()  # fortran-order array
    else:
        if not deviceID == None:
            id = int(deviceID.split(":")[1])
            xp.cuda.Device(id).use()
        volume = xp.array(v.reshape((x, y, z), order=order), dtype=xp.float32).copy()  # fortran-order array

    return volume


def read_em(filename, order='F', keepnumpy=False, deviceID=None, dtype=None):
    """Read EM file. Now only support read the type float32 on little-endian machines.

    @param filename: file name to read.

    @return The data from the EM file in ndarray
    """
    f = open(filename, 'r')

    try:
        dt_header = np.dtype('int32')
        header = np.fromfile(f, dt_header, 128)
        x = header[1]
        y = header[2]
        z = header[3]

        # read the data type
        dt = int(hex(header[0])[2])
        default_type = False

        if dt == 1:  # byte
            raise Exception("Data type not supported yet!")
        elif dt == 2:  # short
            dt_data = np.dtype('<i2')
        elif dt == 4:  # long
            dt_data = np.dtype('<i4')
        elif dt == 5:  # float32
            dt_data = np.dtype('<f4')  # little-endian, float32
            default_type = True
        elif dt == 8:  # float complex
            raise Exception("Data type not supported yet!")
        elif dt == 9:  # double
            dt_data = np.dtype('<f8')
        elif dt == 10:  # double complex
            raise Exception("Data type not supported yet!")
        else:
            raise Exception("Data type not supported yet!")

        v = np.fromfile(f, dt_data, x * y * z)
    finally:
        f.close()

    if keepnumpy:
        volume = np.array(v.reshape((x, y, z), order=order), dtype='float32').copy()  # fortran-order array
    else:  # if the input data is not the default type, convert
        if not deviceID == None:
            id = int(deviceID.split(":")[1])
            xp.cuda.Device(id).use()
        volume = xp.array(v.reshape((x, y, z), order=order), dtype='float32').copy()  # fortran-order array

    return volume


def read_size(filename, dim=''):
    emfile = filename.endswith('.em') * 1
    f = open(filename, 'r')
    try:
        dt_header = np.dtype('int32')
        header = np.fromfile(f, dt_header, 4)
        x = header[0 + emfile]
        y = header[1 + emfile]
        z = header[2 + emfile]
    except:
        raise Exception("reading of MRC file failed")

    f.close()

    if dim == 'x':
        return x
    elif dim == 'y':
        return y
    elif dim == 'z':
        return z
    else:
        return [x, y, z]


def read_header(filename):
    emfile = filename.endswith('.em') * 1
    f = open(filename, 'rb')
    try:
        if emfile:
            header_data = np.fromfile(f, np.dtype('int32'), 128)
        else:
            """
            header_data = []

            header_data += list(np.fromfile(f, np.dtype('int32'), 10))
            header_data += list(np.fromfile(f, np.dtype('float32'), 3))
            header_data += list(np.fromfile(f, np.dtype('int32'), 6))
            header_data += list(np.fromfile(f, np.dtype('float32'), 3))
            header_data += list(np.fromfile(f, np.dtype('int32'), 2))
            header_data += list(np.fromfile(f, np.dtype('float32'), 25))
            header_data += list(np.fromfile(f, np.dtype('int32'), 5))
            header_data += list(np.fromfile(f, np.dtype('float32'), 1))
            header_data += list(np.fromfile(f, np.dtype('int32'), 201))


            print(header_data)
            """
            header_data = np.fromfile(f, np.dtype('float32'), 256)
            # """
    finally:
        f.close()

    return header_data


def read_tilt_angle(filename):
    """Reads the EM header only.
    @param filename: The em file
    """
    emfile = filename.endswith('.em') * 1
    f = open(filename, 'rb')
    try:
        if emfile:
            header_data = np.fromfile(f, np.dtype('int32'), 128)
            tilt_angle = header_data[42] / 1000.
        else:
            header_data = np.fromfile(f, np.dtype('float32'), 256)
            tilt_angle = header_data[43]
    finally:
        f.close()

    return tilt_angle


def read_pixelsize(filename, dim=''):
    import sys
    import numpy as np

    f = open(filename, 'rb')
    try:
        header_data = []

        header_data += list(np.fromfile(f, np.dtype('int32'), 10))
        header_data += list(np.fromfile(f, np.dtype('float32'), 3))
        header_data += list(np.fromfile(f, np.dtype('int32'), 6))
        header_data += list(np.fromfile(f, np.dtype('float32'), 3))
        header_data += list(np.fromfile(f, np.dtype('int32'), 2))
        header_data += list(np.fromfile(f, np.dtype('float32'), 25))
        header_data += list(np.fromfile(f, np.dtype('int32'), 5))
        header_data += list(np.fromfile(f, np.dtype('float32'), 1))
        header_data += list(np.fromfile(f, np.dtype('int32'), 201))

    finally:
        f.close()

    if filename.endswith('.em'):
        sx, sy, sz = header_data[1:4]
    else:
        sx, sy, sz = header_data[:3]

    x, y, z = np.array(header_data[10:13]) / np.array((sx, sy, sz), np.float32)
    if dim == 'x':
        return x
    elif dim == 'y':
        return y
    elif dim == 'z':
        return z
    else:
        return [x, y, z]


def readSubvolumeFromFourierspaceFile(filename, size_x, size_y, size_z):
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
    from pytom.agnostic.io import read
    from pytom.basic.fourier import fourierSizeOperation
    [newX, newY, newZ] = fourierSizeOperation(size_x, size_y, size_z,
                                              reducedToFull=False)

    if filename.__class__ == str:
        originalVolume = read(filename)
    elif filename.__class__ == xp.array:
        # open a backdoor for this function to take volumes, but
        # this should be rather an exception -> not fully documented
        originalVolume = filename
    else:
        raise TypeError('Filename must be a string')

    originalVolume = xp.fft.fftshift(originalVolume, axes=(0, 1))
    newVolume = originalVolume[size_x // 2 - newX // 2:size_x // 2 + newX // 2,
                size_y // 2 - newY // 2:size_y // 2 + newY // 2, :newZ]
    newVolume = xp.fft.fftshift(newVolume, axes=(0, 1))

    return newVolume


# Support Functions Reading
def construct_dtype_relion(filename, relion_options):
    '''Constructs a structured array dtype from header of star file
    @param filename: absolute or relative path to star file
    @param relion_options: dictionary with all header names as defined in relion version
    @return: returns a dtype definition to create numpy structured array'''

    options = relion_options.keys()

    file = open(filename, 'r')
    data = file.readlines()
    file.close()
    DTYPE_OPTICS_GROUP, DTYPE = [], []

    add_optics, add_data = False, True

    size_optics_group = 0

    # Create dtype for parsing relion file
    for n, line in enumerate(data):
        if 'data_optics' in line:
            add_optics, add_data = True, False

        if add_optics:
            size_optics_group = n

        if 'data_particles' in line or 'data_micro' in line:
            add_optics, add_data = False, True

        parts = line.split()
        if not len(parts) > 1:
            continue
        if parts[0].startswith('_rln'):
            name = parts[0][4:]
            if name in options:
                if add_data:
                    DTYPE.append((name, relion_options[name]))

                    if len(DTYPE) != int(parts[-1].replace('#', '')):
                        raise Exception('Parsing Relion File Failed: some header params are not recognised')
                if add_optics:
                    DTYPE_OPTICS_GROUP.append((name, relion_options[name]))


    return DTYPE, DTYPE_OPTICS_GROUP, size_optics_group


# Writing functions

def write(filename, data, tilt_angle=0, pixel_size=1, order='F', fmt=None, header=None, rotation_angles=None):
    """Write data to file. Now only support written in type float32 on little-endian machines.

    @param filename: file name.
    @param data: data which should be written to file. Can be ndarray or pytom volume
    @param tilt_angle: if data is a projection, this is the tilt angle of the projection.
    @param pixel_size: size of pixels/voxels in Angstrom).
    @param order: what is the order of the data (if ndarray).
    @param fmt: if writing txt or star file, fmt defines the output layout. If None, no formatting is used.
    @param header: if writing txt, log or star file, header is header of file. If None, no header is written.
    @param data: data to write.

    """
    import os
    from pytom.lib.pytom_volume import vol

    # Define the allowed file formats and related write function.
    write_functions = {'em': write_em,
                       'mrc': write_mrc, 'rec': write_mrc, 'st': write_mrc, 'mrcs': write_mrc,
                       'txt': write_txt, 'log': write_txt, 'star': write_star, 'meta': write_txt}

    # Extension determines which write function is called
    d, ext = os.path.splitext(filename)
    ext = ext.lower()[1:]

    assert filename
    assert ext in write_functions.keys()

    # If data is instance of vol datatype, convert data to numpy array. Needed because header is defined differently.
    if isinstance(data, vol):
        from pytom.lib.pytom_numpy import vol2npy
        try:
            data_npy = vol2npy(data).copy()
            # Write data to file, using respective write function
            write_functions[ext](filename, data_npy, tilt_angle=tilt_angle, pixel_size=pixel_size, order=order, fmt=fmt,
                                 header=header, rotation_angles=rotation_angles)
        except Exception as e:
            print(data.__class__, e)
            raise Exception('Invalid data type of data. Please provide an ndarray or a pytom volume object.')
    else:
        # Write data to file, using respective write function
        write_functions[ext](filename, data, tilt_angle=tilt_angle, pixel_size=pixel_size, order=order, fmt=fmt,
                             header=header, rotation_angles=rotation_angles)


def write_star(filename, data, tilt_angle=0, pixel_size=1, inplanerot=0, magnification=1., dx=0., dy=0.,
               current_tilt_angle=999, rotation_angles=None,
               order='F', fmt=None, header=None):
    from pytom.basic.files import savetxt

    if header is None:
        header = '\ndata_\n\nloop_\n\n'
        for n, h in enumerate(data.dtype.names):
            header += '_rln{} #{}\n'.format(h, n + 1)

        header = header[:-1]

    if fmt is None:
        fmt = ''

        dtypes = [str(data.dtype[n]) for n in data.dtype.names]

        for dtype in dtypes:
            if dtype in (np.float64, np.float32, 'float32', 'float64'):
                fmt += '%-9.4f '
            elif dtype in (np.int32, 'int32'):
                fmt += '%6d '
            elif ('<U' in str(dtype)):
                fmt += '%60s '

        fmt = fmt[:-1] if len(fmt) else None

    savetxt(filename, data, fmt=fmt, header=header,comments='')


def write_txt(filename, data, tilt_angle=0, pixel_size=1, inplanerot=0, magnification=1., dx=0., dy=0.,
              current_tilt_angle=999, rotation_angles=None,
              order='F', fmt=None, header=None):
    from pytom.basic.files import savetxt

    savetxt(filename, data, fmt=fmt, header=header)


def write_mrc(filename, data, tilt_angle=0, pixel_size=1, inplanerot=0, magnification=1., dx=0., dy=0.,
              current_tilt_angle=999, rotation_angles=None,
              order='F', fmt=None, header=None):
    import numpy as np
    try:
        data = data.get()
    except:
        pass
    if data.dtype != np.dtype('float32'):  # if the origin data type is not float32, convert
        data = data.astype(np.float32)

    assert len(data.shape) < 4 and len(data.shape) > 0

    size = np.ones((3), dtype=np.int32)
    size[:len(data.shape)] = data.shape

    mode = 2
    cell_angles = (0, 0, 0)
    cell_size = pixel_size * size
    dmin, dmax, dmean = data.min(), data.max(), data.mean()

    if current_tilt_angle == 999:
        current_tilt_angle = tilt_angle

    extra = [0, 0, 0, 20140] + [0] * 13 + [dx, dy, tilt_angle, inplanerot, magnification, current_tilt_angle, 0,
                                           0]  # set tiltangle at 43th byte

    if not rotation_angles is None:
        extra[0:3] = rotation_angles

    if np.little_endian:
        machst = 0x00004144
    else:
        machst = 0x11110000

    rms = data.std()

    strings = [
        binary_string(size, np.int32),  # nx, ny, nz
        binary_string(mode, np.int32),  # mode
        binary_string((0, 0, 0), np.int32),  # nxstart, nystart, nzstart
        binary_string(size, np.int32),  # mx, my, mz
        binary_string(cell_size, np.float32),  # cella
        binary_string(cell_angles, np.int32),  # cellb
        binary_string((1, 2, 3), np.int32),  # mapc, mapr, maps
        binary_string((dmin, dmax, dmean), np.float32),  # dmin, dmax, dmean
        binary_string(0, np.int32),  # ispg
        binary_string(0, np.int32),  # nsymbt
        binary_string(extra, np.float32),  # extra
        binary_string((0., 0., 0.), np.int32),  # origin
        binary_string(0, np.int32),  # map
        binary_string(machst, np.int32),  # machst
        binary_string(rms, np.float32),  # rms
        binary_string(0, np.int32),  # nlabl
        binary_string([0] * 200, np.int32),
    ]

    header = b"".join(strings)

    f = open(filename, 'wb')
    try:
        f.write(header)
        f.write(data.tobytes(order=order))  # fortran-order array
    finally:
        f.close()


def write_em(filename, data, tilt_angle=0, pixel_size=1, inplanerot=0, magnification=1., dx=0., dy=0.,
             current_tilt_angle=999, rotation_angles=None,
             order='F', fmt=None, header=None):
    """Write EM file. Now only support written in type float32 on little-endian machines.

    @param filename: file name.
    @param data: data to write.
    """

    try:
        data = data.get()
    except:
        pass

    if data.dtype != np.dtype('float32'):  # if the origin data type is not float32, convert
        data = np.array(data, dtype='float32')

    header = np.zeros(128, dtype='int32')
    header[0] = 83886086  # '0x5000006', TODO: hard-coded, to be changed!
    header[24 + 18] = int(tilt_angle * 1000)  # set tilt angle

    if not rotation_angles is None:
        header[43] = float(rotation_angles[0])
        header[44] = float(rotation_angles[1])
        header[45] = float(rotation_angles[2])

    if len(data.shape) == 3:
        header[1:4] = data.shape
    elif len(data.shape) == 2:
        header[1:3] = data.shape
        header[3] = 1
    elif len(data.shape) == 1:
        header[1] = data.shape[0]
        header[2:4] = 1
    else:
        raise Exception("Input data shape invalid!")

    f = open(filename, 'wb')
    try:
        f.write(header.tobytes())
        f.write(data.tobytes(order=order))  # fortran-order array
    finally:
        f.close()


# Supporting functions write

def binary_string(values, type):
    import numpy as np
    return np.array(values, type).tobytes()


def n2v(data):
    """Transfer Numpy array into a Pytom volume.
    Note the data will be shared between Numpy and Pytom!

    @param data: data to convert.
    """
    try:
        from pytom.lib.pytom_volume import vol
        from pytom.lib.pytom_numpy import npy2vol
    except:
        raise ImportError("Pytom library is not installed or set properly!")

    if data.dtype != np.dtype("float32"):
        data = np.array(data, dtype="float32")

    if len(data.shape) == 3:
        vv = data.copy(order='F')
        v = npy2vol(vv, 3)
    elif len(data.shape) == 2:
        vv = data.reshape(data.shape[0], data.shape[1], 1).copy(order='F')
        v = npy2vol(vv, 2)
    else:
        raise Exception("Data shape invalid!")

    return v
