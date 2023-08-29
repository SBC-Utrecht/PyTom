def power(volume,exponent,inplace=False):
    """
    power: Pixelwise power 
    @param volume: The volume
    @type volume: L{pytom.lib.pytom_volume.vol}
    @param exponent: The exponent
    @type exponent: L{float}
    @param inplace: Perform power inplace? Default is False   
    @type inplace: L{bool}
    @return: volume
    @rtype: L{pytom.lib.pytom_volume.vol}
    """

    if inplace:
        from pytom.lib.pytom_volume import power
        
        power(volume,exponent)

    else:
        #return new volume object
        from pytom.lib.pytom_volume import vol,power
        
        volume2 = vol(volume.size_x(),volume.size_y(),volume.size_z())
        volume2.copyVolume(volume)
        
        power(volume2,exponent)
        
        return volume2
    
    
def determineRotationCenter(particle, binning):
    """
    determineRotationCenter:
    @param particle: The particle 
    @type particle: Either L{pytom.lib.pytom_volume.vol} or string specifying the particle file name
    @param binning: Binning factor
    @return: [centerX,centerY,centerZ]  

    @author: Thomas Hrabe
    """
    if particle.__class__ == str:
        from pytom.lib.pytom_volume import read
        particle = read(particle)
    
    centerX = particle.size_x() / 2.0 * (1.0/float(binning)) 
    centerY = particle.size_y() / 2.0 * (1.0/float(binning))
    centerZ = particle.size_z() / 2.0 * (1.0/float(binning))
    # 
    #if binning > 1:
    #   centerX = centerX - 0.25*(binning-1)
    #    centerY = centerY - 0.25*(binning-1)
    #    centerZ = centerZ - 0.25*(binning-1)
     
    return [centerX,centerY,centerZ]
    
    
    
    
