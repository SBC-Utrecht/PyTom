#!/usr/bin/env python
'''
script to align projections based on fiducials, generate weighted and 
aligned projections, and reconstruct a tomogram using weighted 
backprojection

@author: FF
11.10.2012
'''
import pytom
from pytom.reconstruction.TiltAlignmentStructures import TiltAlignmentParameters, TiltSeries, TiltAlignment

# input parameters
tiltSeriesName="../../../../pytom-testrepo/master/testData/tiltSeries/RM_03" # ending is supposed to be tiltSeriesName_index.em (or mrc)
firstProj=1  # index of first projection
lastProj=41  # index of last projection
ireftilt=21  # reference projection (used for alignment)
markerFileName="../../../../pytom-testrepo/master/testData/tiltSeries/markfile_temp.em.finealig" # marker file (EM/TOM format
irefmark=1   # reference marker (defines 3D coordinate system)
handflip=False  # is your tilt axis outside 0-180 deg?

# output parameters
alignedTiltSeriesName="./aliTest"  # weighted and aligned projections are stored as alignedTiltSeriesName_index.em
projBinning = 4                # binning factor
lowpassFilter=.5               # lowpass filter in Nyquist (post-binning)
onlyWeightedProjections=False  # only write projections and do NOT reconstruct tomogram (following parameters would be obsolete)
voldims=[512,512,128]          # dimensions of reconstructed tomogram
reconstructionPosition=[0,0,0] # offset from center of volume - for example choose z!=0 to shift in z (post-binning coordinates)
volumeName="./testvol.em"      # final tomogram
filetype='em'                  # filetype of final tomogram - can be 'em' or 'mrc'


####################################################################
######## do not edit unless you know what your are doing ###########
####################################################################
MyTiltAlignmentParas=TiltAlignmentParameters(
        dmag=True, drot=True, dbeam=False,
        finealig=True, finealigfile='xxx.dat',
        grad=False,
        irefmark=irefmark, ireftilt=ireftilt, r=None, cent= [1025,1025],
        handflip=handflip,
        optimizer='leastsq', maxIter=0)
MyTiltSeries= TiltSeries(tiltSeriesName=tiltSeriesName, 
                TiltAlignmentParas=MyTiltAlignmentParas, 
                alignedTiltSeriesName=alignedTiltSeriesName,
                markerFileName=markerFileName, 
		firstProj=firstProj, lastProj=lastProj)
MyTiltAlignment = TiltAlignment(MyTiltSeries)
MyTiltAlignment.computeCoarseAlignment( MyTiltSeries)
MyTiltAlignment.alignFromFiducials()
MyTiltSeries.write_aligned_projs( weighting=-1, lowpassFilter=lowpassFilter, 
        binning=projBinning)
if not onlyWeightedProjections:
    vol_bp = MyTiltSeries.reconstructVolume( dims=voldims, 
        reconstructionPosition=reconstructionPosition, binning=1)
    vol_bp.write(volumeName, filetype)

