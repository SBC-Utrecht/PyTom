import pytom
from pytom.reconstruction.TiltAlignmentStructures import TiltAlignmentParameters, TiltSeries, TiltAlignment

MyTiltAlignmentParas=TiltAlignmentParameters(dmag=True, drot=True, dbeam=False,
        finealig=True, finealigfile='xxx.dat',
        grad=False,
        irefmark=1, ireftilt=21, r=None, cent= [1025,1025],
        handflip=False,
        optimizer='leastsq', maxIter=0)
        #optimizer='fmin_slsqp', maxIter=20000)
        #optimizer='fmin_cg', maxIter=20000)
        #optimizer='fmin_powell', maxIter=20000)
        #optimizer='fmin', maxIter=20000)
MyTiltSeries= TiltSeries(tiltSeriesName="./testData/tiltSeries/RM_03", TiltAlignmentParas=MyTiltAlignmentParas, 
                alignedTiltSeriesName="./output/AliTest",
                markerFileName="./testData/tiltSeries/markfile_temp.em.finealig", firstProj=1, lastProj=41)
MyTiltAlignment = TiltAlignment(MyTiltSeries)
MyTiltAlignment.computeCoarseAlignment( MyTiltSeries)


#print MyTiltAlignment
MyTiltAlignment.alignFromFiducials()
#print MyTiltSeries
#MyTiltSeries.write_aligned_projs( weighting=None, lowpassFilter=.5, binning=4)
MyTiltSeries.write_aligned_projs( weighting=-1, lowpassFilter=.5, binning=4)
vol_bp = MyTiltSeries.reconstructVolume( dims=[512,512,128], reconstructionPosition=[0,0,0], binning=1)
vol_bp.write("output/testvol.em")

