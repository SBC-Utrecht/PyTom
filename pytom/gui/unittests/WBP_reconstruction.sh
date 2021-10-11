cd UnitTestProject/03_Tomographic_Reconstruction/tomogram_001

/Users/gijs/Documents/PyTomPrivate/bin/pytom /Users/gijs/Documents/PyTomPrivate/gui/additional/reconstructTomogram.py \
    --tiltSeriesName sorted/sorted \
    --firstIndex 0 \
    --lastIndex 36 \
    --referenceIndex 18 \
    --markerFile reconstruction/WBP/markerfile.em \
    --referenceMarkerIndex 1 \
    --expectedRotationAngle 0 \
    --projectionTargets reconstruction/WBP/temp_files_unweighted/sorted_aligned \
    --projectionBinning 8 \
    --lowpassFilter 0.9  \
    --weightingType 1  \
    --tomogramFile reconstruction/WBP/tomogram_001_WBP.mrc \
    --tiltSeriesFormat mrc \
    --fileType mrc  \
    --tomogramSizeX 232  \
    --tomogramSizeY 232 \
    --tomogramSizeZ 232 \
    --reconstructionCenterX 0 \
    --reconstructionCenterY 0 \
    --reconstructionCenterZ 0


unlink ../../04_Particle_Picking/Tomograms/tomogram_001_WBP.mrc
ln -s UnitTestProject/03_Tomographic_Reconstruction/tomogram_001/reconstruction/WBP/tomogram_001_WBP.mrc ../../04_Particle_Picking/Tomograms/tomogram_001_WBP.mrc