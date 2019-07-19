
def matToZYZ(rotMatrix):
    """
    matToZYZ: Converts a rotation matrix to ZYZ rotation angles
    @param rotMatrix:
    @type rotMatrix: L{pytom.tools.maths.Matrix} 
    @return: [z1,y,z2] Rotation
    """
    from pytom.angles.angle import rad2deg
    from pytom.tools.maths import checkEpsilon
    from math import atan2,acos
    from pytom.tools.maths import epsilon
    
    zyz = [0, 0, 0]

    if abs(rotMatrix[2,2] -1) <= epsilon:
        zyz[0] = 0
        zyz[1] = 0
        zyz[2] = atan2(rotMatrix[0, 1], rotMatrix[0, 0]) * rad2deg
    else:
        zyz[0] = atan2(rotMatrix[2, 1],rotMatrix[2, 0]) * rad2deg
        zyz[1] = acos(rotMatrix[2, 2]) * rad2deg
        zyz[2] = atan2(rotMatrix[1, 2],-rotMatrix[0, 2]) * rad2deg
    
    
    zyz[0] = checkEpsilon(zyz[0])
    zyz[1] = checkEpsilon(zyz[1])
    zyz[2] = checkEpsilon(zyz[2])
        
    return zyz

def zxzToZYZ(z1,z2,x,isRad=False):
    """
    zxzToZYZ
    """
    
    rotMatrix = zxzToMat(z1, z2, x, isRad)
                    
    return matToZYZ(rotMatrix)   
    
def zyzToMat(z1,z2,y,isRad = False):
    """
    zyzToMat
    """
    from pytom.tools.maths import YRotationMatrix,ZRotationMatrix
    
    zm1 = ZRotationMatrix(z1,isRad)
    xm = YRotationMatrix(y,isRad)
    zm2= ZRotationMatrix(z2,isRad)
    
    return (zm2 * (xm * zm1))
    

def angleFromResolutionDeprecated(band,cubeSize,pixelSize=-1,particleDiameter=-1,shrinkFactor=1):
    """
    angleFromResolution: determines the smallest possible rotation angle for a current resolution (see parameters) 
    @param band: The band (aka ring) up to where resolution is good. (in Pixel)
    @param cubeSize: Size of the particle cube. (in Pixel)
    @param pixelSize: Pixelsize in volume. (in Angstrom)
    @param particleDiameter: Diameter of object. (in Angstrom)
    @param shrinkFactor: How many voxels are combined for processing? Default is 1  
    @return: Angle in degrees, or cubeSize / 4 as angle if any of the parameters are invalid
    @deprecated: Yes
    @author: Thomas Hrabe    
    """
    from pytom.angles.angle import rad2deg
    
    if pixelSize == -1:
        print 'Error pixelSize in angleFromResolution, returning default value'
        return cubeSize / 4
    
    if particleDiameter == -1:
        print 'Error particleDiameter in angleFromResolution, returning default value'
        return cubeSize / 4
    
    if not band or band == 0: 
        print 'Error band in angleFromResolution, returning default value'
        return cubeSize / 4
    
    if not cubeSize or cubeSize == 0:
        print 'Error cubeSize in angleFromResolution, returning default value'
        return cubeSize / 4
    
    if shrinkFactor <= 0:
        raise RuntimeError('There is something wrong with your binning parameter! It is <= 0! Please check!')
    
    angle = rad2deg * (cubeSize*pixelSize*shrinkFactor) / (band * particleDiameter)
    
    return angle
    

def zxzToMat(z1,z2=None,x=None,isRad=False):
    """
    zxzToMat : Converts ZXZ angles phi,psi,theta to a rotation matrix
    @param z1: scalar or Rotation 
    @param z2: scalar or None
    @param x: scalar or None 
    @param isRad: Are z1,z2,x in radians?
    @return: The corresponding rotation matrix  
    @author: Thomas Hrabe
    """
    from pytom.tools.maths import XRotationMatrix,ZRotationMatrix
    from pytom.basic.structures import Rotation
    
    if z1.__class__ == Rotation:
        x   = z1.getX()
        z2  = z1.getZ2()
        z1  = z1.getZ1()
    elif z1 == None or z2 == None or x == None:
        raise TypeError('Please provide either a Rotation as z1 or values for z1,z2,x')
    else:
        pass
    zm1 = ZRotationMatrix(z1,isRad)
    xm = XRotationMatrix(x,isRad)
    zm2= ZRotationMatrix(z2,isRad)
    
    return (zm2 * (xm * zm1))

def matToZXZ(rotMatrix,inRad=False):
    """
    matToZXZ : Converts a rotation matrix to zxz angles z1,z2,x.
    @param rotMatrix: The rotation matrix
    @param inRad: Do you want the returned angles be in radians?
    @return: [z1,z2,x] 
    @author: Friedrich Forster
    """
    if not (rotMatrix.getSizeX() == 3 and rotMatrix.getSizeY() == 3):
        raise RuntimeError, 'Input matrix must have a shape of (3x3).'

    import math
    from pytom.basic.structures import Rotation
    from pytom.angles.angle import rad2deg
    from numpy import sign
    #from pytom.tools.maths import epsilon
    epsilon = .001
    if rotMatrix.isIdentity():
        return Rotation(0,0,0)

    # determine X-rotation angle
    cosX = rotMatrix[2,2]
    if cosX >= 1.:
        x = 0.
    elif cosX <= -1.:
        x = math.pi
    else:
        x = math.acos(rotMatrix[2,2])
    sinX = math.sin(x)
    if abs(sinX) >= epsilon:
        cosZ1 = rotMatrix[2,1]/sinX
        if cosZ1 > 1:
            cosZ1 = 1.
        elif cosZ1 < -1.:
            cosZ1 = -1.
        sinZ1 = rotMatrix[2,0]/sinX
        #z1    = math.acos(cosZ1)
        #if sinZ1 < 0.:
        #    z1 = 2*math.pi - z1
        #now z2
        z1    =  arctan(sign(sinX)*rotMatrix[2,0], sign(sinX)*rotMatrix[2,1])
        cosZ2 = -rotMatrix[1,2]/sinX
        sinZ2 = rotMatrix[0,2]/sinX
        if cosZ2 > 1:
            cosZ2 = 1.
        elif cosZ2 < -1.:
            cosZ2 = -1.
        #z2    = math.acos(cosZ2)
        #if sinZ2 < 0.:
        #    z2 = 2*math.pi - z2
        z2    =  arctan(sign(sinX)*rotMatrix[0,2], -sign(sinX)*rotMatrix[1,2])
    else:
        # set z1=0
        z1 = 0.
        cosZ2 = rotMatrix[0,0]
        if cosZ2 > 1.:
            z2 = 0
        elif cosZ2 < -1.:
            z2 = math.pi
        else:
            z2 = math.acos(cosZ2)
        # x=0 deg
        if cosX > 0:
            if rotMatrix[0,1] < 0.:
                z2 = 2*math.pi - z2
        # x=180 deg
        else:
            if rotMatrix[0,1] > 0.:
                z2 = 2*math.pi - z2
        
    if not inRad:
        from pytom.angles.angle import rad2deg
        z1 = z1 *rad2deg
        z2 = z2 *rad2deg
        x = x *rad2deg
    
    return Rotation(z1,z2,x)


def arctan(sinX, cosX):
    """
    compute arctan as function of sin(x) and cos(x)
    @param sinX: sine
    @type sinX: float
    @param cosX: cosine
    @type cosX: float
    @return: atan - automatically map to correct segment
    @rtype: float
    """
    from math import atan, pi
    if cosX == 0:
        if sinX >= 0:
           x = pi/2.
        else:
           x = 3.*pi/2.
    else:
        x = atan(sinX/cosX)
    if cosX < 0:
        x = x + pi
    if x < 0:
        x = x + 2.*pi
    return x


def matToAxisAngle(matrix):
    """
    matToAxisAngle: Converts rotation matrix to axis/angle notation. Formula from (http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToAngle/index.htm)
    @param matrix: 
    @type matrix: L{pytom.tools.maths.Matrix}
    @return: Returns list of angle and [x,y,z]. angle is in degrees!!! [x,y,z] represent a vector of length == 1 
    @rtype: [angle,[x,y,z]]  
    """
    from math import acos,sqrt
    from pytom.angles.angle import rad2deg
    
    if matrix.isSymmetric():
        #means all non diagonal entries are equal
        
        if matrix.isIdentity():
            #return "no rotation"
            return [0,[1,0,0]]
        
        from math import sqrt
        from pytom.tools.maths import epsilon
        #return rotation about 180
        angle = 180
        xx = (matrix[0,0]+1)/2
        yy = (matrix[1,1]+1)/2
        zz = (matrix[2,2]+1)/2
        xy = (matrix[0,1]+matrix[1,0])/4
        xz = (matrix[0,2]+matrix[2,0])/4
        yz = (matrix[1,2]+matrix[2,1])/4
        
        if (xx > yy) and (xx > zz):
            if xx < epsilon: 
                x = 0
                y = 0.7071
                z = 0.7071
            else:
                x = sqrt(xx)
                y = xy/x
                z = xz/x
        elif (yy > zz):
            if yy< epsilon:
                x = 0.7071
                y = 0
                z = 0.7071
            else:
                y = sqrt(yy)
                x = xy/y
                z = yz/y
        else:
            if zz< epsilon:
                x = 0.7071
                y = 0.7071
                z = 0
            else:
                z = sqrt(zz)
                x = xz/z
                y = yz/z
        
        return [angle,[x,y,z]]
    
    value = (matrix[0,0]+matrix[1,1]+matrix[2,2]-1)/2
    if abs(value) -1 < 0.000001: 
        value = float("%.5f" % value)
    else:
        raise RuntimeError()
                 
    angle = acos(value)
    
    div = sqrt(pow(matrix[2,1]-matrix[1,2],2) + pow(matrix[0,2]-matrix[2,0],2) + pow(matrix[1,0]-matrix[0,1],2))
    
    x = (matrix[2,1]-matrix[1,2])/div
    y = (matrix[0,2]-matrix[2,0])/div
    z = (matrix[1,0]-matrix[0,1])/div
    
    return [angle * rad2deg,[x,y,z]]
    
def zxzToAxisAngle(z1,z2,x,isRad=False):
    """
    zxzToAxisAngle: Computes axis/angle notation from ZXZ eulerian angles. 
    @param z1: Phi 
    @param z2: Psi
    @param x: Theta
    @param isRad: Are angles in radians? Default = False. 
    @return: Returns [angle (in rad) , axis[x,y,z]]
    """
    rotMatrix = zxzToMat(z1,z2,x,isRad)
    return matToAxisAngle(rotMatrix)
    
def axisAngleToMat(axis,angle,isRad=False):    
    """
    axisAngleToMat: Converts axis angle representation to a rotation matrix
    @param axis:
    @param angle: 
    @param isRad: Angle in radians? False by default. 
    """
    from pytom.tools.maths import Matrix
    from math import sin,cos
    from pytom.angles.angle import deg2rad
    
    m = Matrix(3,3)
    
    if not isRad:
        angle = angle * deg2rad
    
    cosAngle = cos(angle)
    sinAngle = sin(angle)
    
    x = axis[0]/(axis[0]**2+axis[1]**2+axis[2]**2)**0.5
    y = axis[1]/(axis[0]**2+axis[1]**2+axis[2]**2)**0.5
    z = axis[2]/(axis[0]**2+axis[1]**2+axis[2]**2)**0.5
    
    m[0,0] = (1-cosAngle) *x*x + cosAngle
    m[0,1] = (1-cosAngle) *x*y - z*sinAngle 
    m[0,2] = (1-cosAngle) *x*z + y*sinAngle
    
    m[1,0] = (1-cosAngle) *x*y + z*sinAngle
    m[1,1] = (1-cosAngle) *y*y + cosAngle
    m[1,2] = (1-cosAngle) *y*z - x*sinAngle
    
    m[2,0] = (1-cosAngle) *x*z - y*sinAngle
    m[2,1] = (1-cosAngle) *y*z + x*sinAngle
    m[2,2] = (1-cosAngle) *z*z + cosAngle
    
    return m

def axisAngleToZXZ(axis,angle,isRad=False):
    """
    axisAngleToZXZ:
    @param axis:
    @param angle: 
    @param isRad: Angle in radians? False by default.
    @rtype: L{pytom.basic.structure.Rotation}
    @todo: UnitTest
    """
    m = axisAngleToMat(axis,angle,isRad)
    return matToZXZ(m)

def angleToUnitSphereVector(longitude,latitude,isRadians = False):
    """
    angleToUnitSphereVector:
    @param longitude: ranges from 0 to 360 (set to phi)
    @param latitude: ranges from 0 to 180 (set to theta)
    @param isRadians: default = false
    """
    from math import sin,cos
    from pytom.tools.maths import epsilon
    
    if not isRadians:
        from pytom.angles.angle import deg2rad
        longitude = float(longitude) * deg2rad
        latitude = float(latitude) * deg2rad
        
    v = [0,0,0]
    v[0] = -sin(longitude)
    v[1] = cos(latitude)*cos(longitude)
    v[2] = cos(longitude)*sin(latitude)
    
    if v[0] < epsilon:
        v[0] == 0
    
    if v[1] < epsilon:
        v[1] == 0
        
    if v[2] < epsilon:
        v[2] == 0
    
    return v

def rotationDistances(rotation1,rotation2,isRad=False):
    """
    rotationDistances: Calculates rotation distance on unit sphere and for psi
    @param rotation1: First rotation
    @param rotation2: Second rotation
    @param isRad: optional (default = False)
    @return: List of distances [Distance on unit sphere , psi distance]  
    """
    from math import sqrt,cos,sin,pow,acos,pi
    from pytom.angles.angleFnc import angleToUnitSphereVector
    from pytom.angles.angle import deg2rad,rad2deg
    
    if not isRad:
        
        rotation1[0] = deg2rad* rotation1[0]
        rotation1[1] = deg2rad* rotation1[1]
        rotation1[2] = deg2rad* rotation1[2]
        
        rotation2[0] = deg2rad* rotation2[0]
        rotation2[1] = deg2rad* rotation2[1]
        rotation2[2] = deg2rad* rotation2[2]
       
    m1 = zxzToMat(rotation1[0],rotation1[1],rotation1[2],True)
    m2 = zxzToMat(rotation2[0],rotation2[1],rotation2[2],True)
    
    rot1 = matToZXZ(m1)
    rot2 = matToZXZ(m2)
    
    v1 = angleToUnitSphereVector(rot1[0],rot1[2],False)
    v2 = angleToUnitSphereVector(rot2[0],rot2[2],False)
    
    a = sqrt(pow(v1[0]-v2[0],2)+pow(v1[1]-v2[1],2)+pow(v1[2]-v2[2],2))
    b = sqrt(pow(v1[0],2)+pow(v1[1],2)+pow(v1[2],2))
    c = sqrt(pow(v2[0],2)+pow(v2[1],2)+pow(v2[2],2))
    
    rotationDistance = rad2deg * acos((b*b+c*c-a*a)/(2*b*c))
    
    psiDistance =  abs((rad2deg*rotation1[1])%(360)-(rad2deg*rotation2[1])%(360))
    
    return rotationDistance,psiDistance

def differenceAngleOfTwoRotations(rotation1, rotation2):
    """
    compute difference angles between two rotations
    @param rotation1: 1st rotation
    @type rotation1: L{pytom.basic.structures.Rotation}
    @param rotation2: 2nd rotation
    @type rotation2: L{pytom.basic.structures.Rotation}
    @return: difference angle (in deg)
    @rtype: L{float}
    @author: FF
    """
    from pytom.basic.structures import Rotation
    from math import acos, pi

    assert type(rotation1) == Rotation, "rotation1 must be Rotation!"
    assert type(rotation2) == Rotation, "rotation2 must be Rotation!"
    rot = rotation1*rotation2.invert()
    mat = rot.toMatrix()
    trace = mat[0,0] + mat[1,1] + mat[2,2]
    cosAng =  .5*(trace - 1)
    if cosAng > 1:
        cosAng = 1.
    if cosAng < -1.:
        cosAng = -1.
    the = acos(cosAng) * 180./pi

    return the
    

def angleObjectCheckDistance(angleObject):
    """
    angleObjectCheckDistance : Determines Mean and Standard Deviation of the succesive angular distances in Quaternion space for a given angle object
    @param angleObject: Angle Object
    @type angleObject: L{pytom.angles.angle.AngleObject} 
    @return: [mean,std] - meanDistance and stdDistance
    @author: Thomas Hrabe
    """
    from pytom.angles.quaternions import  Quaternion

    oldRotation = angleObject.nextRotation()
    rotation = angleObject.nextRotation()
    
    distances = []

    while not rotation == [-1,-1,-1]:
        q1 = Quaternion(oldRotation[0],oldRotation[1],oldRotation[2])
        q2 = Quaternion(rotation[0],rotation[1],rotation[2])

        distances.append(q1.distance(q2))
        
        oldRotation = rotation
        rotation = angleObject.nextRotation()
    
    n = len(distances)
    
    mean = 0
    
    for i in xrange(n):
        
        mean = mean + float(distances[i] / n)
    
    std = 0
    from math import sqrt
    for i in xrange(n):
        
        std = std + float((distances[i] - mean) * (distances[i] - mean) / n)
    
    std = sqrt(std)
    
    return [mean,std]
    
        
def z2XFromCenterVector(vector,center=[0,0,0],r = None):       
    """
    z2XFromCenterVector: Determines Z2 and X rotation from a center / vector pair. Helps to determine orientations of object prior to alignment. Formula from Springer Mathematische Formeln - p. 248
    @param vector: Vector pointing to a object. 
    @param center: Center of object 
    @param r: Distance of vector from center (optional)
    @return: [Z2,X] 
    """
    from math import atan2,pi,acos,sqrt
    from pytom.angles.angle import rad2deg
    #shift volume center to center   
    vector[0] = vector[0] - center[0]
    vector[1] = vector[1] - center[1]
    vector[2] = vector[2] - center[2]
    
    radius = r or sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2])
    
    x = acos(vector[2]/radius)
    
    z2 = atan2(vector[1],vector[0]) + pi/2
    #Shift by pi/2 because formula determines left rotated angle
     
    return [rad2deg * z2,rad2deg * x]
    
    
        
def pointRotateZXZ(vector,z1,z2=None,x=None,isDeg = True):
    """
    pointRotateZXZ: Rotates a point / vector around [0,0,0] by z1,z2,x
    @param vector: The point / vector
    @param z1: Either a angle float or a Rotation object
    @type z1: float or L{pytom.basic.structures.Rotation}
    @param z2: float or None
    @param x: float or None
    @param isDeg: True by default      
    """
    from pytom.basic.structures import Rotation
    from pytom.tools.maths import XRotationMatrix,ZRotationMatrix
 
    assert len(vector) == 3
    
    if z1.__class__ == Rotation:
        z2 = z1.getZ2()
        x  = z1.getX()
        z1 = z1.getZ1()

    
    z1Matrix = ZRotationMatrix(z1, not isDeg)
    xMatrix  = XRotationMatrix(x,  not isDeg)
    z2Matrix = ZRotationMatrix(z2, not isDeg)
    
    rotationMatrix = xMatrix * z1Matrix
    rotationMatrix = z2Matrix * rotationMatrix
       
    res = rotationMatrix * vector
    
    return res
    
    
def rotationDistancesFromAngleList(angleList):
    """
    innerClassRotationDistance: Computes distance of one rotation to any other, for all rotations. 
    @param angleList: The current class number
    @type angleList: L{pytom.angles.angleList.AngleList}
    @return: distances[i][j] - distance of rotation i to rotation j
    """
    
    numberRotations = angleList.numberRotations()
    
    listDistances = [[None] * numberRotations for _ in xrange(numberRotations)]
    
    for i in xrange(numberRotations):
        
        quat1 = angleList[i].toQuaternion()
        
        for j in xrange(i+1,numberRotations):
            
            quat2 = angleList[j].toQuaternion()
            
            listDistances[i][j] = quat1.distance(quat2) 
            
    return listDistances   

def meanAndStdRotationDistance(distances):
    """
    meanAndStdRotationDistance: Mean and standart deviation of rotations. Returns mean=0, std=0 if class is empty
    @param distances: as computed by rotationDistancesFromAngleList
    """
    from math import sqrt
    dists = []
    
    counter = 0
    
    for i in xrange(len(distances)):
        for j in xrange(len(distances[i])):
            
            if distances[i][j] and (not i == j):
                dists.append(distances[i][j])
                counter = counter +1
    
    if counter > 0:
        mean = sum(dists) / counter
    else:
        mean = sum(dists)
        #print 'Warning - no class member here!'
        
    std =0    
    for i in xrange(len(dists)):
        std = std + sqrt((dists[i]-mean)*(dists[i]-mean))
    
    if counter > 0:
        std = std / counter
    else:
        std = 0
        #print 'Warning - no class member here!'
    
    return [mean,std]
    

def vector2euler(vec, reference_vec=[0,0,1]):
    """Transform a vector to an Euler angle representation.
    Or: find the Euler angle to rotate the reference vector to the given vector.
    Note there are infinite possible ways to do this, since the inplane rotation is not specified.
    """
    from pytom.tools.maths import euclidianDistance
    if(euclidianDistance(vec, [0,0,0]) == 0):
        raise RuntimeError("Vector length should be bigger than 0!")
    
    if(euclidianDistance(vec, reference_vec) == 0):
        from pytom.basic.structures import Rotation
        return Rotation(0,0,0)
    
    import numpy as np
    vec = vec/np.linalg.norm(vec)
    axis = np.cross(reference_vec, vec)
    angle = np.math.acos(np.dot(vec, reference_vec))
    mat = axisAngleToMat(axis, angle, True)
    rotation = matToZXZ(mat)
    
    return rotation
    

