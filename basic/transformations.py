'''
Created on Apr 29, 2011

@author: hrabe
'''


def shift(volume,shiftX,shiftY,shiftZ,imethod='fourier',twice=False):
    """
    shift: Performs a shift on a volume
    @param volume: the volume
    @param shiftX: shift in x direction
    @param shiftY: shift in y direction
    @param shiftZ: shift in z direction
    @param imethod: Select interpolation method. Real space : linear, cubic, spline . Fourier space: fourier
    @param twice: Zero pad volume into a twice sized volume and perform calculation there.
    @return: The shifted volume.   
    @author: Yuxiang Chen and Thomas Hrabe 
    """
    if imethod == 'fourier':
        from pytom_volume import vol_comp,reducedToFull,fullToReduced,shiftFourier
        from pytom.basic.fourier import fft,ifft,ftshift, iftshift
        
        fvolume = fft(volume)
        fullFVolume = reducedToFull(fvolume)

        destFourier = vol_comp(fullFVolume.sizeX(),fullFVolume.sizeY(),fullFVolume.sizeZ())
        
        shiftFourier(fullFVolume,destFourier,shiftX,shiftY,shiftZ)
        
        resFourier = fullToReduced(destFourier)
        
        return ifft(resFourier)/volume.numelem()
        
    else:
        from pytom_volume import vol
        if imethod == 'linear':
            from pytom_volume import transform
        elif imethod == 'cubic':
            from pytom_volume import transformCubic as transform
        elif imethod == 'spline':
            from pytom_volume import transformSpline as transform
            
        centerX = volume.sizeX()/2
        centerY = volume.sizeY()/2
        centerZ = volume.sizeZ()/2
        
        res = vol(volume.sizeX(),volume.sizeY(),volume.sizeZ())
        transform(volume,res,0,0,0,centerX,centerY,centerZ,shiftX,shiftY,shiftZ,0,0,0)
        
        return res
    
def rotate(volume,rotation,x=None,z2=None,imethod='spline',twice=False):
    """
    rotate: Rotates a volume by a specified rotation
    @param volume: The volume
    @param rotation: The rotation (L{pytom.basic.structures.Rotation}) OR the z1 rotation. 
    @type rotation: Should be of type L{pytom.basic.structures.Rotation}. If not and a scalar is provided, you must provide x, z2 as the remaining two ZXZ rotation angles.
    @param imethod: Interpolation method. Can be: linear, cubic, spline, fourierSpline    
    @param twice: Zero pad volume into a twice sized volume and perform calculation there.
    @return: The rotated volume.   

    @author: Yuxiang Chen and Thomas Hrabe
    """
    
    from pytom.basic.structures import Rotation
    if rotation.__class__ == Rotation and x==None and z2==None:
        z1 = rotation.getZ1()
        x  = rotation.getX()
        z2 = rotation.getZ2()
    elif (rotation.__class__ != float and rotation.__class__ != int) or (x.__class__ != float and x.__class__ != int) or (z2.__class__ != float and z2.__class__ != int):
        
        raise TypeError('Rotation parameter must be a Rotation object or provide z1,x,z2 as floats!')
    else:
        z1 = rotation
        
    if imethod == 'fourier':
        
        return transformFourierSpline(volume,z1,z2,x,0,0,0,twice)
        
    else:
        from pytom_volume import vol
        if imethod == 'linear':
            from pytom_volume import transform
        elif imethod == 'cubic':
            from pytom_volume import transformCubic as transform
        elif imethod == 'spline':
            from pytom_volume import transformSpline as transform
            
        centerX = volume.sizeX()/2
        centerY = volume.sizeY()/2
        centerZ = volume.sizeZ()/2
        
        res = vol(volume.sizeX(),volume.sizeY(),volume.sizeZ())
        transform(volume,res,z1,z2,x,centerX,centerY,centerZ,0,0,0,0,0,0)
        
        return res
    
def transformFourierSpline(volume,z1,z2,x,shiftX,shiftY,shiftZ,twice=False):
    """
    transformFourierSpline: Rotate and shift a volume in fourierspace
    @param volume:
    @param z1:
    @param z2:
    @param x:
    @param shiftX: Shift after rotation
    @param shiftY: Shift after rotation 
    @param shiftZ: Shift after rotation
    @param twice: Zero pad volume into a twice sized volume and perform calculation there.
    @return: The transformed volume.   
    @author: Yuxiang Chen and Thomas Hrabe
    """
    from pytom.basic.fourier import fft, ifft, ftshift, iftshift
    from pytom_volume import vol, pasteCenter, subvolume, transformFourierSpline
    
    if z1 == 0 and z2 == 0 and x == 0:
        return vol(volume)
    
    if twice:
        # paste into a twice sized volume
        v = vol(volume.sizeX()*2, volume.sizeY()*2, volume.sizeZ()*2)
        pasteCenter(volume, v)
    else:
        v = volume
    
    fvol = fft(iftshift(v, inplace=False)) # these steps have to be done in python level because of the fft

    resF = transformFourierSpline(fvol,z1,z2,x,shiftX,shiftY,shiftZ)
    
    res = ftshift(ifft(resF),inplace=False) / v.numelem() # reverse it back
    
    if twice:
        # cut the center part back
        res = subvolume(res, (volume.sizeX()+1)/2, (volume.sizeY()+1)/2, (volume.sizeZ()+1)/2, volume.sizeX(), volume.sizeY(), volume.sizeZ())
    
    return res
    
def rotateFourierSpline(volume,z1,z2,x,twice=False):
    """
    rotateFourierSpline:
    @param volume:
    @param z1:
    @param z2:
    @param x:     
    @param twice: Zero pad volume into a twice sized volume and perform calculation there (=oversampling).
    @type twice: L{bool}
    """
    
#    from pytom.basic.fourier import fft, ifft, ftshift, iftshift
#    from pytom_volume import vol, pasteCenter, subvolume, rotateSplineInFourier
#    
#    if z1 == 0 and z2 == 0 and x == 0:
#        return vol(volume)
#    
#    if twice:
#        # paste into a twice sized volume
#        v = vol(volume.sizeX()*2, volume.sizeY()*2, volume.sizeZ()*2)
#        pasteCenter(volume, v)
#    else:
#        v = volume
#    
#    fvol = fft(iftshift(v, inplace=False)) # these steps have to be done in python level because of the fft
#    
#    resF = rotateSplineInFourier(fvol,z1,z2,x)
#    
#    res = ftshift(ifft(resF),inplace=False) / v.numelem() # reverse it back
#    
#    if twice:
#        # cut the center part back
#        res = subvolume(res, (volume.sizeX()+1)/2, (volume.sizeY()+1)/2, (volume.sizeZ()+1)/2, volume.sizeX(), volume.sizeY(), volume.sizeZ())
#    
#    return res
    
    return transformFourierSpline(volume,z1,z2,x,0,0,0,twice)

def scale(volume,factor,interpolation='Spline'):
    """
    scale: Scale a volume by a factor in REAL SPACE - see also resize function for more accurate operation in Fourier \
    space
    @param volume: input volume
    @type volume: L{pytom_volume.vol}
    @param factor: a factor > 0. Factors < 1 will de-magnify the volume, factors > 1 will magnify.
    @type factor: L{float}
    @param interpolation: Can be Spline (default), Cubic or Linear
    @type interpolation: L{str}

    @return: The scaled volume
    @author: Thomas Hrabe  
    """
    
    if factor <=0:
        raise RuntimeError('Scaling factor must be > 0!')
    
    from pytom_volume import vol
    from math import ceil
    
    if interpolation == 'Spline':
        from pytom_volume import rescaleSpline as rescale
    elif interpolation == 'Cubic':
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

def resize(volume, factor, interpolation='Fourier'):
    """
    resize volume in real or Fourier space
    @param volume: input volume
    @type volume: L{pytom_volume.vol}
    @param factor: a factor > 0. Factors < 1 will de-magnify the volume, factors > 1 will magnify.
    @type factor: L{float}
    @param interpolation: Can be 'Fourier' (default), 'Spline', 'Cubic' or 'Linear'
    @type interpolation: L{str}

    @return: The re-sized volume
    @rtype: L{pytom_volume.vol}
    @author: FF
    """

    if (interpolation == 'Spline') or (interpolation == 'Cubic') or (interpolation == 'Linear'):
        return scale(volume, factor, interpolation='Spline')
    else:
        from pytom.basic.fourier import fft, ifft

        fvol = fft(data=volume)
        newfvol = resizeFourier(fvol=fvol, factor=factor)
        newvol = ifft(newfvol)

        return newvol, newfvol

def resizeFourier(fvol, factor):
    """
    resize Fourier transformed by factor

    @param fvol: Fourier transformed of a volume  - reduced complex
    @type fvol: L{pytom_volume.vol_comp}
    @param factor:  a factor > 0. Factors < 1 will de-magnify the volume, factors > 1 will magnify.
    @type factor: float
    @return: resized Fourier volume (deruced complex)
    @rtype: L{pytom_volume.vol_comp}
    @author: FF
    """
    from pytom_volume import vol_comp

    assert isinstance(fvol, vol_comp), "fvol must be reduced complex"

    oldFNx = fvol.sizeX()
    oldFNy = fvol.sizeY()
    oldFNz = fvol.sizeZ()
    oldNx = oldFNx
    oldNy = fvol.getFtSizeY()
    oldNz = fvol.getFtSizeZ()
    # new dims in real and Fourier space
    newFNx = int(float(oldFNx*factor)+0.5)
    newNx = newFNx
    newNy = int(float(oldNy*factor)+0.5)
    # check 2D images
    if oldNz > 1:
        newNz = int(float(oldNz*factor)+0.5)
        newFNz = newNz // 2 +1
        newFNy = newNy
    else:
        newNz = 1
        newFNz = 1
        newFNy = newNy //2 + 1

    newfvol = vol_comp(newFNx, newFNy, newFNz)
    newfvol.setFtSizeX(newNx)
    newfvol.setFtSizeY(newNy)
    newfvol.setFtSizeZ(newNz)
    scf = 1./(newNx*newNy*newNz)
    # magnify image
    if factor >= 1.:
        newfvol.setAll(0.)
        oldFNx_2 = oldFNx//2+1
        if oldNz == 1:
            iz = 0
            for iy in range(oldFNy):
                for ix in range(oldFNx_2):
                    newfvol.setV( fvol.getV(ix, iy, iz)*scf, ix, iy, iz)
                for ix in range(oldFNx_2-1, oldFNx):
                    ixNew = ix + (newFNx - oldFNx)
                    newfvol.setV( fvol.getV(ix, iy, iz)*scf, ixNew, iy, iz)
        else:
            oldFNy_2 = oldFNy//2+1
            for iz in range(oldFNz):
                for iy in range(oldFNy_2):
                    for ix in range(oldFNx_2):
                        newfvol.setV( fvol.getV(ix, iy, iz)*scf, ix, iy, iz)
                    for ix in range(oldFNx_2-1, oldFNx):
                        ixNew = ix + (newFNx - oldFNx)
                        newfvol.setV( fvol.getV(ix, iy, iz)*scf, ixNew, iy, iz)
                for iy in range(oldFNy_2-1, oldFNy):
                    iyNew = iy + (newFNy - oldFNy)
                    for ix in range(oldFNx_2):
                        newfvol.setV( fvol.getV(ix, iy, iz)*scf, ix, iyNew, iz)
                    for ix in range(oldFNx_2-1, oldFNx):
                        ixNew = ix + (newFNx - oldFNx)
                        newfvol.setV( fvol.getV(ix, iy, iz)*scf, ixNew, iyNew, iz)
    # de-magnify image
    else:
        newFNx_2 = newFNx//2+1
        if oldFNz == 1:
            iz = 0
            for iy in range(newFNy):
                for ix in range(newFNx_2):
                    newfvol.setV( fvol.getV(ix, iy, iz)*scf, ix, iy, iz)
                for ix in range(newFNx_2-1, newFNx):
                    ixOld = ix + (oldFNx - newFNx)
                    newfvol.setV( fvol.getV(ixOld, iy, iz)*scf, ix, iy, iz)
        else:
            newFNy_2 = newFNy//2+1
            for iz in range(newFNz):
                for iy in range(newFNy_2):
                    for ix in range(newFNx_2):
                        newfvol.setV( fvol.getV(ix, iy, iz)*scf, ix, iy, iz)
                    for ix in range(newFNx_2-1, newFNx):
                        ixOld = ix + (oldFNx - newFNx)
                        newfvol.setV( fvol.getV(ixOld, iy, iz)*scf, ix, iy, iz)
                for iy in range(newFNy_2-1, newFNy):
                    iyOld = iy + (oldFNy - newFNy)
                    for ix in range(newFNx_2):
                        newfvol.setV( fvol.getV(ix, iyOld, iz)*scf, ix, iy, iz)
                    for ix in range(newFNx_2-1, newFNx):
                        ixOld = ix + (oldFNx - newFNx)
                        newfvol.setV( fvol.getV(ixOld, iyOld, iz)*scf, ix, iy, iz)


    return newfvol

def mirror(volume,axis = 'x',copyFlag = True):
    """
    mirror: Mirrors a volume at defined axis (x,y,z)
    @param volume: 
    @param axis:
    @param copyFlag: If true, will mirror volume into a new volume (default). If false, will do mirror inplace.
    @return: The mirrored volume
    @author: Thomas Hrabe  
    """
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

def general_transform_crop(v, rot=None, shift=None, scale=None, order=[0, 1, 2]):
    """
    Perform general transformation using 3rd order spline interpolation 
    using volume of identical size. The origin stays invariant upon rotation
    AND magnification!

    @param v: volume
    @type v: L{pytom_volume.vol}
    @param rot: rotate
    @type rot: pytom.basic.structures.Rotate or list
    @param shift: shift
    @type shift: pytom.basic.structures.Shift or list
    @param scale: scale / magnification along each dimension
    @type scale: list
    @param order: the order in which the three operations are performed (smaller means first) \
    e.g.: [2,1,0]: rotation: 2, scale: 1, translation:0 => translation 1st
    @type order: list
    @return: pytom_volume

    @author: FF
    """
    from pytom.basic.structures import Rotation, Shift
    from pytom_volume import vol
    if not isinstance(v,vol):
        raise TypeError('general_transform_crop: v must be of type pytom_volume.vol! Got ' + str(v.__class__) + ' instead!')

    if rot is None:
        rot = Rotation(0.,0.,0.)
    if rot.__class__ == list and len(rot) == 3:
        rot = Rotation(rot[0], rot[1], rot[2])
    if shift is None:
        shift = Shift(0.,0.,0.)
    if shift.__class__ == list and len(shift) == 3:
        shift = Shift(shift[0], shift[1], shift[2])
    if scale is None:
        scale = [1.0, 1.0, 1.0]
    if scale.__class__ == list and len(scale) == 3:
        #check
        if int(v.sizeX()*scale[0]) < 1 or int(v.sizeY()*scale[1]) < 1 or int(v.sizeZ()*scale[2]) < 1:
            raise RuntimeError("Scale not possible! Please check all the dimension after scaling is bigger than 1!")
    else:
        raise TypeError("Scale parameter invalid! Should be a list of 3 values!")
    
    from pytom_volume import vol, general_transform
    from pytom.tools.maths import Matrix
    
    # invert matrix
    rotM = rot.toMatrix(True)
    shiftM = shift.toMatrix()
    scaleM = Matrix(4,4)
    scaleM[0,0] = scale[0]
    scaleM[1,1] = scale[1]
    scaleM[2,2] = scale[2]
    scaleM[3,3] = 1
    
    # multiply them according to the order
    rotCenter1 = Shift(-int(v.sizeX()/2), -int(v.sizeY()/2), -int(v.sizeZ()/2)).toMatrix()
    rotCenter2 = Shift(int(v.sizeX()/2), int(v.sizeY()/2), int(v.sizeZ()/2)).toMatrix()

    # multiply them according to the order
    all_mtx = [None, None, None]
    all_mtx[order[0]] = rotCenter2 * (rotM * rotCenter1) # for the rotation center!
    all_mtx[order[1]] = shiftM
    all_mtx[order[2]] = rotCenter2 * (scaleM* rotCenter1) # for the magnification center!
    mtx = all_mtx[2] * (all_mtx[1] * all_mtx[0])

    res = vol(v.sizeX(), v.sizeY(), v.sizeZ())
    general_transform(v, res, mtx._matrix)
    return res

def general_transform(v, rot=None, shift=None, scale=None, order=[0, 1, 2]):
    """Perform general transformation using 3rd order spline interpolation.
    @param v: volume
    @type v: L{pytom_volume.vol}
    @param rot: rotate
    @type rot: pytom.basic.structures.Rotate or list
    @param shift: shift
    @type shift: pytom.basic.structures.Shift or list
    @param scale: scale / magnification along each dimension
    @type scale: list
    @param order: the order in which the three operations are performed (smaller means first)
    @type order: list
    @return: pytom_volume
    """
    from pytom.basic.structures import Rotation, Shift

    from pytom_volume import vol
    if not isinstance(v,vol):
        raise TypeError('general_transform: v must be of type pytom_volume.vol! Got ' + str(v.__class__) + ' instead!')

    if rot is None:
        rot = Rotation(0.,0.,0.)
    if rot.__class__ == list and len(rot) == 3:
        rot = Rotation(rot[0], rot[1], rot[2])
    if shift is None:
        shift = Shift(0.,0.,0.)
    if shift.__class__ == list and len(shift) == 3:
        shift = Shift(shift[0], shift[1], shift[2])
    if scale is None:
        scale = [1.0, 1.0, 1.0]
    if scale.__class__ == list and len(scale) == 3:
        #check
        if int(v.sizeX()*scale[0]) < 1 or int(v.sizeY()*scale[1]) < 1 or int(v.sizeZ()*scale[2]) < 1:
            raise Exception("Scale not possible! Please check all the dimension after scaling is bigger than 1!")
    else:
        raise Exception("Scale parameter invalid! Should be a list of 3 values!")
    
    from pytom_volume import vol, general_transform
    from pytom.tools.maths import Matrix
    
    # invert matrix
    rotM = rot.toMatrix(True)
    shiftM = shift.toMatrix()
    scaleM = Matrix(4,4)
    scaleM[0,0] = scale[0]
    scaleM[1,1] = scale[1]
    scaleM[2,2] = scale[2]
    scaleM[3,3] = 1
    
    # multiply them according to the order
    all_mtx = [None, None, None]
    try:
        if order[2] > order[0]: # rotation first
            rotCenter1 = Shift(-int(v.sizeX()/2), -int(v.sizeY()/2), -int(v.sizeZ()/2)).toMatrix()
            rotCenter2 = Shift(int(v.sizeX()/2), int(v.sizeY()/2), int(v.sizeZ()/2)).toMatrix()
        else: # scale first, so the center is different
            rotCenter1 = Shift(-int(v.sizeX()*scale)/2, -int(v.sizeY()*scale)/2, -int(v.sizeZ()*scale)/2).toMatrix()
            rotCenter2 = Shift(int(v.sizeX()*scale)/2, int(v.sizeY()*scale)/2, int(v.sizeZ()*scale)/2).toMatrix()
        all_mtx[order[0]] = rotCenter2 * (rotM * rotCenter1) # for the rotation center!
        all_mtx[order[1]] = shiftM
        all_mtx[order[2]] = scaleM
    except:
        raise Exception("The given order is wrong! Should be a list of 0,1,2!")
    
    mtx = all_mtx[2] * (all_mtx[1] * all_mtx[0])
    
    res = vol(int(v.sizeX()*scale[0]), int(v.sizeY()*scale[1]), int(v.sizeZ()*scale[2]))
    general_transform(v, res, mtx._matrix)
    return res

def general_transform2d(v, rot=None, shift=None, scale=None, order=[0, 1, 2], crop=True):
    """Perform general transformation of 2D data using 3rd order spline interpolation.
    @param v: volume in 2d (it's got to be an image)
    @type v: L{pytom_volume.vol}
    @param rot: rotate
    @type rot: float
    @param shift: shift
    @type shift: L{pytom.basic.structures.Shift} or list [x,y]
    @param scale: uniform scale
    @type scale: float
    @param order: the order in which the three operations are performed (smaller means first) \
    e.g.: [2,1,0]: rotation: 2, scale: 1, translation:0 => translation 1st
    @type order: list
    @param crop: crop the resulting volume to have the same size as original volume (default: True)
    @type crop: boolean
    @return: pytom_volume
    """
    from pytom.basic.structures import Shift

    from pytom_volume import vol
    if not isinstance(v,vol):
        if isinstance(v,tuple) and isinstance(v[0],vol):
            v = v[0]
        else:
            raise TypeError('general_transform2d: v must be of type pytom_volume.vol! Got ' + str(v.__class__) + ' instead!')
    if v.sizeZ() != 1:
        raise RuntimeError('general_transform2d: v must be 2D!')
    if rot is None:
        rot = [0., 0., 0.]
    elif rot.__class__ != int and rot.__class__ != float:
        raise RuntimeError("Please provide only 1 angle because you are dealing with 2D data!")
    else:
        rot = [rot, 0., 0.]

    if shift.__class__ == list and len(shift) == 2:
        shift = Shift(shift[0], shift[1], 0.)
    if shift.__class__ == list and len(shift) == 3:
        shift = Shift(shift[0], shift[1], 0.) # you cannot shift 2D image along z axis, can you?

    if crop:
        vv= general_transform_crop(v, rot=rot, shift=shift, scale=[scale, scale, 1], order=order)
    else:
        vv = general_transform(v, rot, shift, [scale, scale, 1], order)
    return vv

def project(v, rot, verbose=False):
    """
    rotate and subsequently project volume along z
    @param v: volume
    @type v: L{pytom_volume.vol}
    @param rot: rotation - either Rotation object or single angle interpreted as rotation along y-axis
    @type rot: L{pytom.basic.structures.Rotation} or float
    @return: projection (2D image)
    @rtype: L{pytom_volume.vol}
    @author: FF
    """
    from pytom_volume import vol
    from pytom_numpy import vol2npy, npy2vol
    from numpy import ndarray, float32
    from pytom.basic.transformations import general_transform_crop
    from pytom.basic.structures import Rotation
    from numpy import sum as proj

    if not isinstance(v, vol):
        raise TypeError('project: v must be of type pytom_volume.vol! Got ' + str(v.__class__) + ' instead!')
    if isinstance(rot, float):
        rot = Rotation(z1=90., z2=270., x=rot, paradigm='ZXZ')
    if not isinstance(rot, Rotation):
        raise TypeError('project: rot must be of type Rotation or float! Got ' + str(v.__class__) + ' instead!')
    if verbose:
        print("project: Rotation for projection: "+str(rot))
    rotvol = general_transform_crop(v=v, rot=rot, shift=None, scale=None, order=[0, 1, 2])
    # slightly awkward: projection in numpy ...
    npvol = vol2npy(rotvol)
    npprojection = ndarray([v.sizeX(), v.sizeY(), 1], dtype=float32, buffer=None, offset=0, strides=npvol.strides,
                           order='F')
    proj(npvol, axis=2, dtype=float32, out=npprojection)
    projection = npy2vol(npprojection,3)
    return projection
    
