import unittest
import pytom
from pytom.reconstruction.TiltAlignmentStructures import TiltAlignmentParameters, TiltSeries, TiltAlignment
from pytom.reconstruction.tiltAlignmentFunctions import simulate_markers
from pytom.tools.files import checkDirExists

class pytom_TiltAlignTest(unittest.TestCase):
    def setUp(self):
        """set up"""
        self.MyTiltAlignmentParas=TiltAlignmentParameters(dmag=True, drot=True, dbeam=False,
            finealig=True, finealigfile='xxx.dat',
            grad=False,
            irefmark=1, ireftilt=21, r=None, cent= [1025,1025],
            handflip=False,
            optimizer='leastsq', maxIter=0)
        self.MyTiltSeries= TiltSeries(tiltSeriesName="./testData/tiltSeries/RM_03", 
	        TiltAlignmentParas=self.MyTiltAlignmentParas, 
                alignedTiltSeriesName="./output/AliTest",
                markerFileName="./testData/tiltSeries/markfile_temp.em.finealig", firstProj=1, lastProj=41)
        self.MyTiltAlignment = TiltAlignment(self.MyTiltSeries)
        self.MyTiltAlignment.computeCoarseAlignment( self.MyTiltSeries)

        if not checkDirExists('output'):
            import os
            os.mkdir('output')

    #def optimizableParametersTest(self):
    #    """
#	test copy to optimizable parameters
#        """
#        print self.MyTiltAlignment
#        opti_pars = self.MyTiltAlignment.getOptimizableVariables( self.MyTiltAlignmentParas)
#	print opti_pars
#	opti_pars = range(0,len(opti_pars))
#        self.MyTiltAlignment.setOptimizableVariables( self.MyTiltAlignmentParas, opti_pars)
#        print self.MyTiltAlignment

    def testConsistency(self):

        from numpy import array
        markCoords = []
        tiltAngles = array(len(self.MyTiltSeries._ProjectionList)*[0.])
        for (imark, Marker) in enumerate(self.MyTiltSeries._Markers):
            markCoords.append(Marker.get_r())
        for (itilt, proj) in enumerate(self.MyTiltSeries._ProjectionList):
            tiltAngles[itilt]=proj.getTiltAngle()
            (MyMarkers, MyTiltSeries, MyTiltAlignmentParas) = simulate_markers(markCoords, 
                tiltAngles, tiltAxis=-76.71, ireftilt=None,
                ampTrans=100., ampRot=.5, ampMag=.01, dBeam=None, dMagnFocus=None)

        MyTiltSeries.writeMarkerFile('output/markfile_test.em')
        # print MyTiltSeries
        MyTiltSeries= TiltSeries(tiltSeriesName="./testData/tiltSeries/RM_03", 
    	        TiltAlignmentParas=self.MyTiltAlignmentParas, 
                    alignedTiltSeriesName="./output/AliTest",
                    markerFileName="./output/markfile_test.em", firstProj=1, lastProj=41)
        MyTiltAlignment = TiltAlignment(MyTiltSeries)
        MyTiltAlignment.computeCoarseAlignment( MyTiltSeries)
        score = MyTiltAlignment.alignFromFiducials()
        # print MyTiltSeries
        self.assertTrue( score < 0.01, 'Tilt Alignment Consistency damaged!')

    def testIMODread(self):
        """
	test reading IMOD alignment file
        """
        self.MyTiltSeries.readIMODwimp( 
	        markerFileName='./testData/ctt_coll_1_13_markers.wimp', 
	        prexgfile='./testData/ctt_coll_1_13.prexg', verbose=True)


    #def runTest(self):
    #    self.optimizableParametersTest()

    def tearDown(self):
        if checkDirExists('output'):
            import os
            os.system('rm -rf output')

if __name__ == '__main__':
    unittest.main()



