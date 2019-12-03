# -*- coding: utf-8 -*-
'''
Created on Oct 29, 2010

@author:Thomas Hrabe
'''

def _display(volume,sliceNumber = 0,projectionAxis='z'):
    """
    _display: Generate image for volume display
    @author: Thomas Hrabe
    @param volume: The image / volume
    @param sliceNumber: Slice number at which the volume will be cut
    @param projectionAxis: Defines the axis. Default is z, means you look on the xy plane. (x equals to 'yz' plane, y to 'xz') 
    """
    import Image
    from pytom_volume import min,max,mean,variance,limit,subvolume
    from math import sqrt
    
    if projectionAxis == 'x':
        size1 = volume.sizeY()
        size2 = volume.sizeZ()
    if projectionAxis == 'y':
        size1 = volume.sizeX()
        size2 = volume.sizeZ()    
    if projectionAxis == 'z':
        size1 = volume.sizeX()
        size2 = volume.sizeY()

    img = Image.new('L',(size1,size2))
    
    if sliceNumber > 0:
        
        if projectionAxis == 'x':
            volume = subvolume(volume,sliceNumber,0,0,1,size1,size2)
        elif projectionAxis == 'y':
            volume = subvolume(volume,0,sliceNumber,0,size1,1,size2)
        elif projectionAxis == 'z':
            volume = subvolume(volume,0,0,sliceNumber,size1,size2,1)
         
    elif sliceNumber == 0 and volume.sizeZ() > 1:
        sliceNumber = int(volume.sizeZ()/2)
        
        if projectionAxis == 'x':
            volume = subvolume(volume,sliceNumber,0,0,1,size1,size2)
        elif projectionAxis == 'y':
            volume = subvolume(volume,0,sliceNumber,0,size1,1,size2)
        elif projectionAxis == 'z':
            volume = subvolume(volume,0,0,sliceNumber,size1,size2,1)
    
    minValue = min(volume)
    volume.shiftscale(-minValue,1)
    maxValue = max(volume)
    if maxValue == 0:
        maxValue = 1
    volume.shiftscale(0,255.0/maxValue)
        
    for i in range(size1):
        for j in range(size2):
            if projectionAxis == 'x':
                value = volume(0,i,j)
            elif projectionAxis == 'y':
                value = volume(i,0,j)
            elif projectionAxis == 'z':
                value = volume(i,j,0)
            
            img.putpixel((i,j),int(value))
    
    return img

def display(volume,sliceNumber = 0,projectionAxis='z'):
    """
    display: Will quickly display image / volume on screen
    @param volume: The image / volume
    @param sliceNumber: Slice number at which the volume will be cut
    @param projectionAxis: Defines the axis. Default is z, means you look \
    on the xy plane. (x equals to 'yz' plane, y to 'xz') 
    @author: Thomas Hrabe
    """
    img = _display(volume,sliceNumber,projectionAxis)
    img.show()
    
def volumeToPNG(volume,pngPath,sliceNumber = 0,projectionAxis='z'):
    """
    volumeToPNG: Takes a volume and writes a PNG file 
    @param volume: 
    @param pngPath:
    @param sliceNumber: Slice number at which the volume will be cut
    @param projectionAxis: Defines the axis. Default is z, means you look on the xy plane. (x equals to 'yz' plane, y to 'xz') 
    @author: Thomas Hrabe
    """
    png = _display(volume,sliceNumber,projectionAxis)
    png.save(pngPath,format = 'PNG')
    

def volumeSequenceToPNGs(volume,sequenceDirectory,sliceStep,projectionAxis='z'):
    """
    volumeSequenceToPNGs: Creates a slice sequence along an specified axis 
    @param volume:
    @param sequenceDirectory:
    @param sliceStep:
    @param projectionAxis: Defines the axis. Default is z, means you look on the xy plane. (x equals to 'yz' plane, y to 'xz') 
    @author: Thomas Hrabe
    """
    if projectionAxis == 'x':
        size = volume.sizeX()
    elif projectionAxis == 'y':
        size = volume.sizeY()
    elif projectionAxis == 'z':
        size = volume.sizeZ()
    else:
        raise ValueError('projectionAxis must either be x,y, or z!')
    
    for i in range(0,size,sliceStep):
        volumeToPNG(volume,sequenceDirectory + str(i) + '.png',i,projectionAxis)


def _dspcub(volume,sigma=None,projectionAxis='z'):
    """
    _dspcub: Creates a dspcub image
    dspcub2PNG: display z-slices in 2D image
    @param volume: The input volume
    @type volume: L{pytom_volume.vol}
    @parameter sigma: thresholding as multiples of standard deviation sigma
    @type sigma: float
    @param projectionAxis: Defines the axis. Default is z, means you look on the xy plane. (x equals to 'yz' plane, y to 'xz') 
    @return: Image
    @rtype: Image
    @@author: Pia Unverdorben 
    """
    
    if volume.sizeZ() == 1:
        raise Exception('You can use ')
    
    import Image
    from pytom_volume import min,max,mean,variance,limit,subvolume
    from math import sqrt, ceil, floor
    
    sizeX=volume.sizeX()
    sizeY=volume.sizeY()
    sizeZ=volume.sizeZ()
    
    if projectionAxis == 'x':
        imagesPerRow=int(floor(sqrt(sizeX)))
        size1 = float(sizeX)
        sizeI = sizeY
        sizeJ = sizeZ
    elif projectionAxis == 'y':
        imagesPerRow=int(floor(sqrt(sizeY)))
        size1 = float(sizeY)
        sizeI = sizeX
        sizeJ = sizeZ
    elif projectionAxis == 'z':
        imagesPerRow=int(floor(sqrt(sizeZ)))
        size1=float(sizeZ)
        sizeI = sizeX
        sizeJ = sizeY
        
    numberIterations=imagesPerRow*imagesPerRow
    
    if size1<numberIterations:    
        iterationSteps=float(numberIterations/size1)
    else:
        iterationSteps=float(size1/numberIterations) 
    
    iterationSteps = int(iterationSteps)
    
    if projectionAxis == 'x':
        img=Image.new('L', (sizeY*imagesPerRow, sizeZ*imagesPerRow))
    elif projectionAxis == 'y':
        img=Image.new('L', (sizeX*imagesPerRow, sizeZ*imagesPerRow))
    elif projectionAxis == 'z':
        img=Image.new('L', (sizeX*imagesPerRow, sizeY*imagesPerRow))    

    # normalize according to standard deviation if sigma is specified
    if sigma:
        meanv = mean(volume)
        stdv = sqrt(variance(volume,False))
        minValue = meanv - float(sigma) * stdv
        maxValue = meanv + float(sigma) * stdv
    else:
        minValue = min(volume)
        maxValue = max(volume)
        
    for sliceNumber in range(0,numberIterations, iterationSteps):
        
        if projectionAxis == 'x':
            png=Image.new('L', (sizeY, sizeZ))
        elif projectionAxis == 'y':
            png=Image.new('L', (sizeX, sizeZ))
        elif projectionAxis == 'z':
            png=Image.new('L', (sizeX, sizeY))
            
        (yvalue,xvalue)=divmod(sliceNumber, imagesPerRow)
        
        if projectionAxis == 'x':
            vol = subvolume(volume,sliceNumber,0,0,1,sizeY,sizeZ)
        elif projectionAxis == 'y':
            vol = subvolume(volume,0,sliceNumber,0,sizeX,1,sizeZ)
        elif projectionAxis == 'z':
            vol = subvolume(volume,0,0,sliceNumber,sizeX,sizeY,1)
            
        vol.shiftscale(-minValue,1)
        vol.shiftscale(0,255/maxValue)
        
        for i in range(sizeI):
            for j in range(sizeJ):
            
                if projectionAxis == 'x':
                    value = vol(0,i,j)
                elif projectionAxis == 'y':
                    value = vol(i,0,j)
                elif projectionAxis == 'z':
                    value = vol(i,j,0)
                             
                png.putpixel((i,j),int(value))
                      
        img.paste(png, (xvalue*sizeI, yvalue*sizeJ)) 
    
    return img

def dspcub(volume,sigma=None,projectionAxis='z'):
    """
    _dspcub: Displays a volume in dspcub
    dspcub2PNG: display z-slices in 2D image based on Pia's code
    @param volume: The input volume
    @type volume: L{pytom_volume.vol}
    @parameter sigma: thresholding as multiples of standard deviation sigma
    @param projectionAxis: Defines the axis. Default is z, means you look on the xy plane. (x equals to 'yz' plane, y to 'xz')
    @type sigma: float
    """
    
    img = _dspcub(volume,sigma,projectionAxis)
    img.show()
    
def dspcub2PNG(volume,pngPath,sigma=None,projectionAxis='z'):
    """
    dspcub2PNG: display z-slices in 2D image and save as PNG
    @param volume: The input volume
    @type volume: L{pytom_volume.vol}
    @param pngPath: Path to file written to disk
    @type pngPath: str    
    @parameter sigma: thresholding as multiples of standard deviation sigma
    @param projectionAxis: Defines the axis. Default is z, means you look on the xy plane. (x equals to 'yz' plane, y to 'xz')
    @type sigma: float
    @author: Pia Unverdorben 
    @date: 18/11/2010
    """
    
    img = _dspcub(volume,sigma,projectionAxis)
                     
    img.save(pngPath, format = 'PNG')

    
    
    
