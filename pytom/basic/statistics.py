'''
a couple of functions for statistical analysis of volumes

'''
def mean(volume,mask,fVolume=0,fMask=0,numberOfMaskVoxels=-1):
    """
    suggestion frido: rename to mean_moving_mask
    
    mean: Determines the mean of volume under mask. 
    Each assigned value is determined through the value of the voxels surrounding the current voxel that are covered by the mask.
    @param volume: The volume of interest
    @param mask: A mask under which the mean is determined
    @param fVolume: Optional - the fouriertransformed volume (fft must not be repeated)
    @param fMask: Optional - the fouriertransformed mask (fft must not be repeated)     
    @param numberOfMaskVoxels: Optional - the number of voxels != 0 
    """
    import pytom_volume;
    if not fVolume.__class__ == pytom_volume.vol_comp:
        from pytom.basic.fourier import fft;
        fVolume = fft(volume);
    
    if not fMask.__class__ == pytom_volume.vol_comp:
        from pytom.basic.fourier import fft;
        fMask = fft(mask);

    if numberOfMaskVoxels <0:
        numberOfMaskVoxels = pytom_volume.numberSetVoxels(mask);

    fMeanVol = fVolume * fMask;
    from pytom.basic.fourier import ifft;
    meanVolume = ifft(fMeanVol);
    
    from pytom.basic.fourier import iftshift;
    iftshift(meanVolume);
    
    meanVolume.shiftscale(0,1.0/(numberOfMaskVoxels*meanVolume.numelem()));

    return meanVolume;


def std(volume,mask,meanVolume=0,fMask=0,numberOfMaskVoxels=-1):
    """
    suggestion frido: rename to std_moving_mask
    
    std: Determines the std of volume under moving mask. 
    Each assigned value is determined through the value of the voxels surrounding the current voxel that are covered by the mask.
    @param volume: The volume of interest
    @param mask: A mask under which the std is determined
    @param meanVolume: Optional - the meanVolume determined by mean. (mean must not be recalculated)
    @param fMask: Optional - the fouriertransformed mask (fft must not be repeated)
    @param numberOfMaskVoxels: Optional - the number of voxels != 0     
    """
    from pytom.basic.fourier import fft,ifft,iftshift;
    from pytom_volume import power;
    import pytom_volume;
    
    if not fMask.__class__ == pytom_volume.vol_comp:
        fMask = fft(mask);
    
    if not meanVolume.__class__ == pytom_volume.vol:
        meanVolume = mean(volume,mask);
    
    if numberOfMaskVoxels <0:
        numberOfMaskVoxels = pytom_volume.numberSetVoxels(mask);

    volumeSqrd = volume * volume;
    fVolumeSqrd = fft(volumeSqrd);
    
    fVolumeSqrd = fVolumeSqrd * fMask;
    
    varianceVolume = iftshift(ifft(fVolumeSqrd));
    varianceVolume.shiftscale(0,1.0/(numberOfMaskVoxels*varianceVolume.numelem()));
    
    power(meanVolume,2);
    
    varianceVolume = varianceVolume - meanVolume;
    
    power(varianceVolume,0.5);
    
    return varianceVolume;

#def mean(volume, mask=None):
#    """calculate mean of volume - if mask specified only under mask
#    
#       @param volume: 
#       @param mask:
#    """
#    
#    resV = mean_moving_mask(volume, mask)
#    res = resV.getV(resV.size_x()/2,resV.size_y()/2,resV.size_z()/2)
#    
#    return res
#    
#def std(volume, mask=None):
#    """calculate standard deviation of volume - if mask specified only under mask
#    
#       @param volume: 
#       @param mask:
#    """
    
def statisticOfDistanceMatrix(distanceMatrix):
    """
    statisticOfDistanceMatrix: Determines max, mean and deviation of a distance matrix
    @param distanceMatrix: Must be a square matrix!
    @return: [max,mean,std] of distanceMatrix
    @author: Thomas Hrabe
    """
    from pytom_volume import vol,max
    from pytom.tools.maths import listMean,listStd

    if not distanceMatrix.__class__ == vol:
        raise TypeError('Parameter must be a pytom_volume.vol!')

    if distanceMatrix.size_z() > 1:
        raise RuntimeError('Parameter must be a 2D pytom_volume.vol!')

    if distanceMatrix.size_x() != distanceMatrix.size_y():
        raise RuntimeError('Matrix must be a square! Size x != y')

    values = []

    for i in range(distanceMatrix.size_x()):
        if i < distanceMatrix.size_x():
            for j in range(i+1,distanceMatrix.size_y()):
                values.append(distanceMatrix(i,j,0))

    m = listMean(values)
    s = listStd(values, m)

    return [max(distanceMatrix),m,s]
    

def averagePlanes(volume,sliceStart,sliceEnd,sliceStep=1,axis='Z'):
    """
    averagePlanes:
    @param volume:
    @param sliceStart:
    @param sliceEnd:
    @param axis:
    @author:     
    """
    
    from pytom_volume import subvolume
    
    planes = []
    
    
    for i in range(sliceStart,sliceEnd,sliceStep):
        
        if axis == 'X':
            v = subvolume(volume,0,0,i,volume.size_x(),volume.size_y(),1)
            planes.append(v)
            
        if axis == 'Y':
            v = subvolume(volume,0,0,i,volume.size_x(),volume.size_y(),1)
            planes.append(v)
        
        if axis == 'Z':
            v = subvolume(volume,0,0,i,volume.size_x(),volume.size_y(),1)
            planes.append(v)
    
    if len(planes)>0:
        p = planes[0]
        
        for i in range(len(planes)-1):
            
            p += planes[i+1]
            
        return p
    
    else:
        
        raise RuntimeError('Could not select on single plane from your parameters.')
    
            
