cd UnitTestProject/03_Tomographic_Reconstruction/tomogram_001

/Users/gijs/Documents/PyTomPrivate/bin/pytom /Users/gijs/Documents/PyTomPrivate/gui/additional/reconstructTomogram.py \
    --tiltSeriesName reconstruction/INFR/sorted  \
    --firstIndex 0 \
    --lastIndex 37 \
    --referenceIndex 19 \
    --markerFile reconstruction/INFR/markerfile.em \
    --referenceMarkerIndex 1 \
    --expectedRotationAngle 0 \
    --projectionTargets reconstruction/INFR/temp_files_unweighted/sorted_aligned \
    --projectionBinning 8 \
    --lowpassFilter 0.9 \
    --weightingType 0 \
    --tiltSeriesFormat em \
    --fileType em

unlink ../../04_Particle_Picking/Tomograms/tomogram_001_INFR.em
/Users/gijs/Documents/PyTomPrivate/bin/pytom /Users/gijs/Documents/PyTomPrivate/reconstruction/reconstruct_INFR.py -d reconstruction/INFR/temp_files_unweighted/ -o reconstruction/INFR/tomogram_001_INFR.em 
ln -s UnitTestProject/03_Tomographic_Reconstruction/tomogram_001/reconstruction/INFR/tomogram_001_INFR.em ../../04_Particle_Picking/Tomograms/tomogram_001_INFR.em