def add(volume,SNR=1):
    """
    add Adds white noise to volume
    @param volume: A volume
    @param SNR: Signal to Noise ratio of result
    @type SNR: int or float > 0 
    @return: Volume containing noise with SNR == SNR
    @author: Thomas Hrabe  
    """
    
    if(SNR < 0):
        return volume
    
    from math import sqrt
    from pytom_volume import vol,mean,variance,gaussianNoise
    
    m = mean(volume)
    s = sqrt(variance(volume,False)/SNR) # SNR = Var(signal) / Var(noise)
    
#    s = sqrt(variance(volume,False)/SNR)-variance(volume,False)
#    
#    if s<0:
#        return volume
#    elif s==0:
#        s =1

    noise = vol(volume.sizeX(),volume.sizeY(),volume.sizeZ())
    
    gaussianNoise(noise,m,s) # s is actually the std

    result = volume + noise

    return result
    
    
    
    