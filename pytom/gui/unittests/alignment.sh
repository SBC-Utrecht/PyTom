cd UnitTestProject/03_Tomographic_Reconstruction/tomogram_001

/Users/gijs/Documents/PyTomPrivate/bin/pytom /Users/gijs/Documents/PyTomPrivate/gui/additional/reconstructTomogram.py \
    --tiltSeriesName sorted/sorted  \
    --firstIndex 0 \
    --lastIndex 36 \
    --referenceIndex 18 \
    --markerFile alignment/markerfile.em \
    --referenceMarkerIndex 1 \
    --expectedRotationAngle 0 \
    --projectionTargets alignment/unweighted_unbinned_marker_1/sorted_aligned \
    --projectionBinning 1 \
    --lowpassFilter 0.9 \
    --tiltSeriesFormat mrc \
    --weightingType 0 