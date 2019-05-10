'''
Created on Sep 6, 2011
last change: aug 31, 2012 FF

@author: hrabe, foerster
'''
def initSphere(sizeX, sizeY, sizeZ, radius, smooth=0, maxradius=0, cent=None, filename=''):
    """
    initSphere: Initializes a sphere
    @param sizeX: X size of cube
    @type sizeX: int
    @param sizeY: Y size of cube
    @type sizeY: int
    @param sizeZ: Z size of cube
    @type sizeZ: int
    @param radius: Radius of sphere mask
    @type radius: float
    @param smooth: Smooth of sphere mask
    @type smooth: float
    @param maxradius: maximum radius - data will be set to zero beyond this radius - default 0 does not do anything
    @param cent: centre vector
    @type cent: array (3-dim)
    @param filename: If specified by the user, the spherical mask will be written to disk.    
    """
    from pytom_volume import initSphere, vol
    
        
    v = vol(sizeX,sizeY,sizeZ)
    
    if cent:
        initSphere(v, radius, smooth, maxradius, cent[0], cent[1], cent[2])
    else:
        initSphere(v, radius, smooth, maxradius, sizeX/2, sizeY/2, sizeZ/2)
    
    if filename != '':
        v.write(filename)
    
    return v

def taper_edges(image, width, taper_mask=None):
    """
    taper edges of image (or volume) with cos function

    @param image: input image (or volume)
    @type image: pytom volume
    @param width: width of edge
    @type width: int
    @param taper_mask: mask for tapering - if None it will be generated
    @type taper_mask: L{pytom_volume.vol}

    @return: image with smoothened edges, taper_mask
    @rtype: array-like

    @author: FF
    """
    import numpy
    from pytom_volume import vol
    from math import cos, pi

    assert type(image) == vol, "taper_edges: image must be type pytom_volume.vol"
    dims = [image.sizeX(), image.sizeY(), image.sizeZ()]
    assert dims[0]>1, "taper_edges: image must be 2D or 3D"
    assert dims[1]>1, "taper_edges: image must be 2D or 3D"

    width = int(width)

    if taper_mask==None:
        # pre-compute scale factor at borders
        fact = numpy.array(list(range(1,width+1)),dtype='float32')
        for (ii,ifact) in enumerate(fact):
            fact[ii] = cos(ifact*pi/(2.*(width)))
        #2-dim
        if dims[2] == 1:
            taper_mask = vol(dims[0],dims[1],1)
            taper_mask.setAll(1.)
            for ix in range(0,dims[0]):
                ixleft  = ix-width+1
                ixright = dims[0]-ix-width
                for iy in range(0,dims[1]):
                    iyleft  = iy-width+1
                    iyright = dims[1]-iy-width
                    if ( ixleft<=0 or ixright<=0 or iyleft<=0 or iyright<=0):
                        ii = abs(min( [ ixleft, ixright, iyleft, iyright] ))
                        taper_mask.setV(float(fact[ii]),ix,iy,0)
        #3-dim
        else:
            taper_mask = vol(dims[0],dims[1],dims[2])
            taper_mask.setAll(1.)
            for ix in range(0,dims[0]):
                ixleft  = ix-width+1
                ixright = dims[0]-ix-width
                for iy in range(0,dims[1]):
                    iyleft  = iy-width+1
                    iyright = dims[1]-iy-width
                    for iz in range(0,dims[2]):
                        izleft  = iz-width+1
                        izright = dims[2]-iz-width
                        if ( ixleft<=0 or ixright<=0 or iyleft<=0 or iyright<=0 
                                    or izleft<=0 or izright<=0):
                            ii = abs(min( [ixleft,ixright,iyleft,iyright,
                                    izleft,izright] ))
                            taper_mask.setV(float(fact[ii]),ix,iy,iz)
    image = image * taper_mask
    return image, taper_mask

def limit_in_sphere( invol, r_max=None, lowlimit=None, lowval=0., hilimit=None, hival=0., reduced=True):
    """
    limit grey values of volume within center radius

    @param invol: input volume
    @type invol: L{pytom_volume.vol} or L{pytom_volume.vol_comp}
    @param r_max: radius
    @type r_max: L{int}
    @param lowlimit: lower limit - all values below this value that lie in the specified radius will be replaced \
                by lowval
    @type lowlimit: L{float}
    @param lowval: replacement value
    @type lowval: L{float}
    @param hilimit: upper limit - all values above this value that lie in the specified radius will be replaced \
                by hival
    @type hilimit: L{float}
    @param hival: replacement value
    @type hival: L{float}

    @param reduced: reduced complex?
    @type reduced: L{bool}

    @author: FF
    """
    from math import sqrt
    if not r_max:
        r_max=invol.sizeY()/2-1
    if reduced:
        centX1 = 0
        centX2 = invol.sizeX()-1
        centY1 = 0
        centY2 = invol.sizeY()-1
        centZ  = 0
    else:
        centX1 = int(invol.sizeX()/2+.5)
        centY1 = int(invol.sizeY()/2+.5)
        centZ  = int(invol.sizeZ()/2+.5)
    for ix in range(0,invol.sizeX()):
        for iy in range(0,invol.sizeY()):
            for iz in range(0,invol.sizeZ()):
                if reduced:
                    d1 = (ix-centX1)**2
                    d2 = (ix-centX2)**2
                    dx = min(d1,d2)
                    d1 = (iy-centY1)**2
                    d2 = (iy-centY2)**2
                    dy = min(d1,d2)
                else:
                    dx = (ix-centX1)**2
                    dy = (iy-centY1)**2
                dz = (iz-centZ)**2
                r = sqrt(dx+dy+dz)
                if r < r_max:
                    v = invol.getV( ix, iy, iz)
                    if lowlimit:
                        if v < lowlimit:
                            invol.setV( lowval, ix, iy, iz)
                    if hilimit:
                        if v > hilimit:
                            invol.setV( hival, ix, iy, iz)

def scale(volume,factor,interpolation='Spline'):
    """
    scale: Scale (enlarge/shrink) a volume by a factor (=change size) . 
    @param volume:
    @param factor: a factor > 0. Factors < 1 will downscale the volume, factors > 1 will upscale.
    @param interpolation: Can be Spline (default), Lagrange or Linear
    @return: The scaled volume
    @author: Thomas Hrabe  
    @deprecated: Use L{pytom.basic.transformations.scale} instead!
    """
    print('pytom.basic.functions.scale: This function is deprecated. Use pytom.basic.transformations.scale instead!')
    if factor <=0:
        raise RuntimeError('Scaling factor must be > 0!')
    
    from pytom_volume import vol
    from math import ceil
    
    if interpolation == 'Spline':
        from pytom_volume import rescaleSpline as rescale
    elif interpolation == 'Lagrange':
        from pytom_volume import rescaleCubic as rescale
    elif interpolation == 'Linear':
        from pytom_volume import rescale
         
    
    sizeX = volume.sizeX()
    sizeY = volume.sizeY()
    sizeZ = volume.sizeZ()
    newSizeX = int(ceil(sizeX * float(factor)))
    newSizeY = int(ceil(sizeY * float(factor)))
    newSizeZ = int(ceil(sizeZ * float(factor)))

    newVolume = vol(newSizeX,newSizeY,newSizeZ)
    rescale(volume,newVolume)
    
    return newVolume

def mirror(volume,axis = 'x',copyFlag = True):
    """
    mirror: Mirrors a volume at defined axis (x,y,z)
    @param volume: 
    @param axis:
    @param copyFlag: If true, will mirror volume into a new volume (default). If false, will do mirror inplace.
    @return: The mirrored volume
    @author: Thomas Hrabe  
    @deprecated: Use L{pytom.basic.transformations.mirror} instead!
    """
    print('pytom.basic.functions.mirror: This function is deprecated. Use pytom.basic.transformations.mirror instead!')
    if axis == 'x':
        transformation = [-1,1,1]
    elif axis == 'y':
        transformation = [1,-1,1]
    elif axis == 'z':
        transformation = [1,1,-1]
    
    centerX = volume.sizeX()
    centerY = volume.sizeY()
    centerZ = volume.sizeZ()
    
    if copyFlag:
        from pytom_volume import vol
        returnVolume = vol(volume.sizeX(),volume.sizeY(),volume.sizeZ())
        
        for x in range(volume.sizeX()):
            for y in range(volume.sizeY()):
                for z in range(volume.sizeZ()):
                    
                    xMirrored = (x-centerX) * transformation[0] + centerX
                    yMirrored = (y-centerY) * transformation[1] + centerY
                    zMirrored = (z-centerZ) * transformation[2] + centerZ
                    
                    returnVolume.setV(volume.getV(x,y,z),xMirrored,yMirrored,zMirrored)
        return returnVolume
        
    else:
        for x in range(volume.sizeX()):
            for y in range(volume.sizeY()):
                for z in range(volume.sizeZ()):
                    
                    tmp = volume.getV(x,y,z)
                    xMirrored = (x-centerX) * transformation[0] + centerX
                    yMirrored = (y-centerY) * transformation[1] + centerY
                    zMirrored = (z-centerZ) * transformation[2] + centerZ
                    
                    volume.setV(volume.getV(xMirrored,yMirrored,zMirrored),x,y,z)
                    
                    volume.setV(tmp,xMirrored,yMirrored,zMirrored)
