'''
Created on May 12, 2011

@author: luiskuhn
'''

def particleVolume(particleList, templateVolume, dimX, dimY, dimZ, volume=None):
    
    from pytom_volume import vol, paste, rotate, subvolume
    
    if volume == None:
        volume = vol(dimX, dimY, dimZ)
        volume.setAll(0.0)
    
    for p in particleList:
        
        x = p.getPickPosition().getX()
        y = p.getPickPosition().getY()
        z = p.getPickPosition().getZ()
    
        z1 = p.getRotation().getZ1()
        z2 = p.getRotation().getZ2()
        x1 = p.getRotation().getX()        
              
        tempTemplate = vol(templateVolume.size_x(), templateVolume.size_y(), templateVolume.size_z())
        tempTemplate.setAll(0.0)
        rotate(templateVolume, tempTemplate, z1, z2, x1)
        tempTemplate = tempTemplate + subvolume(volume, int(x -templateVolume.size_x()/2), int(y -templateVolume.size_y()/2), int(z -templateVolume.size_z()/2), templateVolume.size_x(), templateVolume.size_y(), templateVolume.size_z())        
        paste(tempTemplate, volume, int(x -templateVolume.size_x()/2), int(y -templateVolume.size_y()/2), int(z -templateVolume.size_z()/2))
  
    return volume


def multipleParticlesVolume(particleList, templateMap, dimX, dimY, dimZ):
    
    from pytom_volume import read, vol
    
    volumeMap = {}
    keys = list(templateMap.keys())
    for key in keys:        
        volume = vol(dimX, dimY, dimZ)
        volume.setAll(0.0)        
        volumeMap[key] = volume
        
        templateVolume = read(templateMap[key])
        templateMap[key] = templateVolume        
            
        
    for i in range(len(particleList)):
        
        p = particleList[i]
        templateVolume = templateMap[p.getClassName()]
        volume = volumeMap[p.getClassName()]
        
        volume = particleVolume(particleList[i:i+1], templateVolume, 0, 0, 0, volume)
          
    return volumeMap

