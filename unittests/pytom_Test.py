def pytom_TestSuite():
    """
    
    """
    from unittest import TestSuite
    suite = TestSuite()
    
    import pytom_LibraryTest
    suite.addTest(pytom_LibraryTest.pytom_LibraryTest("library_Test"))
    #---------------------
    import pytom_IOTest
    suite.addTest(pytom_IOTest.pytom_IOTest("read_Test"))
    
    #---------------------
    import pytom_AngleTest
    suite.addTest(pytom_AngleTest.pytom_AngleTest("EulerAngleList_oldRotations_Test"))
    suite.addTest(pytom_AngleTest.pytom_AngleTest("EulerAngleList_numberRotations_Test"))
    suite.addTest(pytom_AngleTest.pytom_AngleTest("LocalSampling_numberRotations_Test"))
    suite.addTest(pytom_AngleTest.pytom_AngleTest("LocalSampling_reset_Test"))
    suite.addTest(pytom_AngleTest.pytom_AngleTest("LocalSampling_focusRotation_Test"))
    suite.addTest(pytom_AngleTest.pytom_AngleTest("LocalSampling_init_Test"))
    suite.addTest(pytom_AngleTest.pytom_AngleTest("AV3Sampling_Test"))
    suite.addTest(pytom_AngleTest.pytom_AngleTest("GlobalLocalCombined_Test"))
    suite.addTest(pytom_AngleTest.pytom_AngleTest("AngleEMList_focusRotation_Test"))
    suite.addTest(pytom_AngleTest.pytom_AngleTest("AngleEMList_focusRotation_Test"))
    suite.addTest(pytom_AngleTest.pytom_AngleTest("angleFromResolution_Test"))
    suite.addTest(pytom_AngleTest.pytom_AngleTest("pointRotateZXZ_Tests"))
    suite.addTest(pytom_AngleTest.pytom_AngleTest("EquidistantList_Test"))
    suite.addTest(pytom_AngleTest.pytom_AngleTest("zxzToMatToZXZ_Test"))
    suite.addTest(pytom_AngleTest.pytom_AngleTest("axisAngleToMatToAxisAngle_Test"))
    suite.addTest(pytom_AngleTest.pytom_AngleTest("indexAccess_Test"))
    suite.addTest(pytom_AngleTest.pytom_AngleTest("AngleList_sliceAccess_Test"))
    suite.addTest(pytom_AngleTest.pytom_AngleTest("LocalSamplingNumberAngles_Test"))
    suite.addTest(pytom_AngleTest.pytom_AngleTest("matrixMult_Test"))
    #---------------------    
    import pytom_MathTest
    suite.addTest(pytom_MathTest.pytom_MathTest("Matrix_Test"))
    suite.addTest(pytom_MathTest.pytom_MathTest("Identity_Test"))
    suite.addTest(pytom_MathTest.pytom_MathTest("RotationMatrix_Test"))
    suite.addTest(pytom_MathTest.pytom_MathTest("sumRotations_Test"))
    suite.addTest(pytom_MathTest.pytom_MathTest("Pcacov_Test"))
    #---------------------
    import pytom_XMLTest
    suite.addTest(pytom_XMLTest.pytom_XMLTest("ParticleXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("ReferenceXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("ParticleListXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("ReferenceListXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("ExpectationMaximisationJobXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("MaximisationJobXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("MaximisationResultXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("ExpectationJobXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("ExpectationJobMsgXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("ExpectationResultXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("ExpectationResultMsgXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("EulerAngleListXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("AngleEquidistantXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("AngleListFromEMXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("PreprocessingXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("WedgeInfoXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("AsymmetricWedgeInfoXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("StatusMessageXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("AlignmentListXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("ScoresFromXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("GrowingAverageJobXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("SymmetryXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("PeakPriorXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("EquidistantListXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("MaskXML_Test"))
    ##suite.addTest(pytom_XMLTest.pytom_XMLTest("CorrelationMatrixJobXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("ReferenceListXML_Test"))
    ##suite.addTest(pytom_XMLTest.pytom_XMLTest("CorrelationVectorJobXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("ShiftXML_Test"))
    ##suite.addTest(pytom_XMLTest.pytom_XMLTest("MultiRefAlignJobXML_Test"))
    ##suite.addTest(pytom_XMLTest.pytom_XMLTest("SimulatedAnnealingAlignJobXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("MCOEXMXJobXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("AnnealingTemperature_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("AnnealingJob_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("ClusterSwap_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("SwapList_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("ProjectionXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("ProjectionListXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("AnnealingCriterionXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("SingleTiltWedgeXML_Test"))
    suite.addTest(pytom_XMLTest.pytom_XMLTest("DoubleTiltWedgeXML_Test"))
    #---------------------
    
    import pytom_CorrelationTest 
    suite.addTest(pytom_CorrelationTest.pytom_CorrelationTest("XCF_Test"))
    suite.addTest(pytom_CorrelationTest.pytom_CorrelationTest("NXCC_Test"))
    suite.addTest(pytom_CorrelationTest.pytom_CorrelationTest("NXCF_Test"))
    suite.addTest(pytom_CorrelationTest.pytom_CorrelationTest("FLCF_Test"))
    suite.addTest(pytom_CorrelationTest.pytom_CorrelationTest("SOC_Test"))
    suite.addTest(pytom_CorrelationTest.pytom_CorrelationTest("FSC_Test"))
    suite.addTest(pytom_CorrelationTest.pytom_CorrelationTest("SubPixelPeak_Test"))
    #---------------------
    
    import pytom_NormTest
    suite.addTest(pytom_NormTest.pytom_NormTest("mean0std1_Test"))

    #---------------------
    import pytom_StructuresTest
    suite.addTest(pytom_StructuresTest.pytom_StructuresTest("Wedge_getWedgeVolume_Test"))
    suite.addTest(pytom_StructuresTest.pytom_StructuresTest("Wedge_apply_Test"))
    suite.addTest(pytom_StructuresTest.pytom_StructuresTest("DoubleTiltWedge_Test"))
    suite.addTest(pytom_StructuresTest.pytom_StructuresTest("ParticleList_getParticleByFilename_Test"))
    suite.addTest(pytom_StructuresTest.pytom_StructuresTest("ParticleList_sliceTest"))
    suite.addTest(pytom_StructuresTest.pytom_StructuresTest("ParticleList_assignTest"))
    suite.addTest(pytom_StructuresTest.pytom_StructuresTest("ParticleList_addTest"))
    suite.addTest(pytom_StructuresTest.pytom_StructuresTest("ParticleList_checkTest"))
    suite.addTest(pytom_StructuresTest.pytom_StructuresTest("Mask_checkTest"))
    suite.addTest(pytom_StructuresTest.pytom_StructuresTest("Reference_checkTest"))
    suite.addTest(pytom_StructuresTest.pytom_StructuresTest("ReferenceList_accessTest"))
    suite.addTest(pytom_StructuresTest.pytom_StructuresTest("AlignmentList_sortByParticleListTest"))
    suite.addTest(pytom_StructuresTest.pytom_StructuresTest("Shift_OperatorTest"))
    
    #---------------------
    
    import pytom_FilterTest
    suite.addTest(pytom_FilterTest.pytom_FilterTest("ones_Test"))
    suite.addTest(pytom_FilterTest.pytom_FilterTest("complexDiv_Test"))
    suite.addTest(pytom_FilterTest.pytom_FilterTest("profile_Test"))
    suite.addTest(pytom_FilterTest.pytom_FilterTest("wedgeFilter_Test"))
    suite.addTest(pytom_FilterTest.pytom_FilterTest("wedgeRotation_Test"))
    suite.addTest(pytom_FilterTest.pytom_FilterTest("ramp_Test"))

    #---------------------

    import pytom_ScoreTest
    suite.addTest(pytom_ScoreTest.pytom_ScoreTest("xcfScore_Test"))
    suite.addTest(pytom_ScoreTest.pytom_ScoreTest("nxcfScore_Test"))
    suite.addTest(pytom_ScoreTest.pytom_ScoreTest("flcfScore_Test"))
    suite.addTest(pytom_ScoreTest.pytom_ScoreTest("socScore_Test"))
    #---------------------
    
    import pytom_SimulationTest
    suite.addTest(pytom_SimulationTest.pytom_SimulationTest("gaussianNoise_Test"))
    #---------------------
    
    import pytom_LocalizationTest
    suite.addTest(pytom_LocalizationTest.pytom_LocalTest("Volume_Test"))
    suite.addTest(pytom_LocalizationTest.pytom_LocalTest("Orientation_Test"))
    suite.addTest(pytom_LocalizationTest.pytom_LocalTest("FoundParticle_Test"))
    suite.addTest(pytom_LocalizationTest.pytom_LocalTest("PeakJob_Test"))
    suite.addTest(pytom_LocalizationTest.pytom_LocalTest("PeakResult_Test"))
    suite.addTest(pytom_LocalizationTest.pytom_LocalTest("PeakJobMsg_Test"))
    suite.addTest(pytom_LocalizationTest.pytom_LocalTest("PeakResultMsg_Test"))
    suite.addTest(pytom_LocalizationTest.pytom_LocalTest("BackwardCompatibility_Test"))

    #---------------------
    import pytom_SymmetryTest
    suite.addTest(pytom_SymmetryTest.pytom_SymmetryTest("HelicalSymmetry_ApplyToParticleTest"))
    #suite.addTest(pytom_SymmetryTest.pytom_SymmetryTest("HelicalSymmetry_ApplyToParticleList"))
    #---------------------
    import pytom_AlignTest
    
    suite.addTest(pytom_AlignTest.pytom_AlignTest("ExpectationMaximisationJobCheck_Test"))
    suite.addTest(pytom_AlignTest.pytom_AlignTest("Peak_Test"))
    suite.addTest(pytom_AlignTest.pytom_AlignTest("bestAlignment_Test"))
    suite.addTest(pytom_AlignTest.pytom_AlignTest("BackwardCompatibility_Test"))
    ## suite.addTest(pytom_AlignTest.pytom_AlignTest("Aligner_Test"))
    ##suite.addTest(pytom_AlignTest.pytom_AlignTest("Expectation_Test"))
    ##suite.addTest(pytom_AlignTest.pytom_AlignTest("Symmetry_Test"))
    #---------------------
    
    ##   import pytom_FittingTest
    ##    suite.addTest(pytom_FittingTest.pytom_FittingTest("RefinementFitting_Test"))
    #---------------------
       
    ##import pytom_ResultsTest
    ##suite.addTest(pytom_ResultsTest.pytom_ResultsTest("assessClassification_Test"))
    ##suite.addTest(pytom_ResultsTest.pytom_ResultsTest("determineClassSwaps_Test"))
    #---------------------
    #import pytom_NumpyTest
    ## suite.addTest(pytom_NumpyTest.pytom_NumpyTest("numpy_Test"))
    #---------------------
    
    import pytom_ToolsTest
    suite.addTest(pytom_ToolsTest.pytom_ToolsTest("VolumeStorage_Test"))
    
    #---------------------
    
    import pytom_ReconstructionTest
    #suite.addTest(pytom_ReconstructionTest.pytom_ReconstructionTest("goldReconstruction_Test"))
    
    #---------------------
    import pytom_ServerpagesTest
    suite.addTest(pytom_ServerpagesTest.pytom_ServerpagesTest("ReconstructionCall_Test"))
    suite.addTest(pytom_ServerpagesTest.pytom_ServerpagesTest("LocalizationCall_Test"))
    #suite.addTest(pytom_ServerpagesTest.pytom_ServerpagesTest("AlignmentCall_Test"))
    #suite.addTest(pytom_ServerpagesTest.pytom_ServerpagesTest("MCOEXMXCall_Test"))
    #suite.addTest(pytom_ServerpagesTest.pytom_ServerpagesTest("MCOACCall_Test"))

    #---------------------
    import pytom_BinTest
    #suite.addTest(pytom_BinTest.pytom_BinTest("localization_Test"))
    #suite.addTest(pytom_BinTest.pytom_BinTest("alignment_Test"))
    #suite.addTest(pytom_BinTest.pytom_BinTest("mcoEXMX_Test"))
    #suite.addTest(pytom_BinTest.pytom_BinTest("mcoAC_Test"))
    #---------------------
    import pytom_TiltAlignTest
    #suite.addTest(pytom_TiltAlignTest.pytom_TiltAlignTest("testConsistency"))

    #---------------------
    import pytom_FRMTest
    suite.addTest(pytom_FRMTest.pytom_FRMTest("FRM_Test"))

    #---------------------
    import pytom_TransformationTest
    suite.addTest(pytom_TransformationTest.pytom_TransformationTest("test_resize2D"))
    suite.addTest(pytom_TransformationTest.pytom_TransformationTest("test_resize3D"))
    suite.addTest(pytom_TransformationTest.pytom_TransformationTest("testOrigin"))

    #---------------------
    import pytom_LocalOptimizationTest
    suite.addTest(pytom_LocalOptimizationTest.pytom_LocalOptimizationTest("localOpti_Test"))

    #---------------------
    #---------------------
    #---------------------
    #---------------------
    #---------------------
    #---------------------
    
    return suite

def run():


    pytom_Suite = pytom_TestSuite()
    from unittest import TextTestRunner
    
    runner = TextTestRunner()
    runner.run(pytom_Suite)
