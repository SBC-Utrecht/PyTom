'''
Created on Mar 2, 2011

@author: luiskuhn
'''

def theta_vol_projection(vol_src, theta):
        
    from pytom.lib.pytom_volume import vol, rotate, subvolume, sum
    
    vol_src_dim_x = vol_src.size_x()
    vol_src_dim_y = vol_src.size_y()
    vol_src_dim_z = vol_src.size_z()
    
    vol_dst = vol(vol_src_dim_x, vol_src_dim_y, vol_src_dim_z)
    vol_dst.setAll(0.0)
    
    rotate(vol_src, vol_dst, 270, 90, theta)

    vol_img = vol(vol_src_dim_x, vol_src_dim_y, 1)
    vol_img.setAll(0.0)
    
    for i in range(vol_src_dim_x):
        for j in range(vol_src_dim_y):
            
            vol_img.setV(sum(subvolume(vol_dst, i, j, 0, 1, 1, vol_src_dim_z)), i, j, 0)
    
    return vol_img
    
    
  


def getWeightedProjectionCube(sourceVolume, thetaAngles):
    
    from pytom.lib.pytom_volume import vol, paste, complexRealMult
    
    from pytom.basic.fourier import fft
    from pytom.basic.fourier import ifft
    
    from pytom.basic.filter import fourierFilterShift
    from pytom.basic.filter import circleFilter
    from pytom.basic.filter import rampFilter
    
    
    dimX = sourceVolume.size_x()
    dimY = sourceVolume.size_y()
        
    imageCube = vol(dimX, dimY, len(thetaAngles))
    imageCube.setAll(0.0)
        
    weightSlice = fourierFilterShift(rampFilter(sourceVolume))
    
    circleFilterRadius = dimX/2
    circleSlice = fourierFilterShift(circleFilter(sourceVolume, circleFilterRadius))

    for k in range(len(thetaAngles)):

        angle = thetaAngles[k]    
    
        image = theta_vol_projection(sourceVolume, angle)
    
        ##filtered_img = ifft( complexRealMult(complexRealMult(complexRealMult(fft(image), filter_slice), circle_filter_slice), ctf) ) 
        filteredImage = ifft( complexRealMult( complexRealMult( fft(image), weightSlice), circleSlice) ) 
    
        paste(filteredImage, imageCube, 0, 0, k)
        
        #for i in range(dimX):
        #    for j in range(dimY):        
        #        imageCube.setV(filteredImage.getV(i, j, 0), i, j, k)
                
    
    return imageCube
                
     
