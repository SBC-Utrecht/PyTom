'''
Created on Sep 24, 2009

@author: hrabe
'''

"""
@var epsilon: value taken from numpy for float64 using this command 
numpy.finfo(float).eps
epsilon = 2.2204460492503131e-16
"""

"""
@var epsilon: value is a float epsilon now
"""
epsilon = 0.0000005

def rotate_vector(r, rotation):
    """
    rotate 2d or 3d vector by specified rotation

    @param r: vector
    @type r: 2d or 3d array
    @param rotation: rotation (float for 2d or 3d vector containing Euler angles Z1,X1,Z2; all in deg)
    @type rotation: float or 3d array

    """
    from math import cos, sin, pi
    if len(r) == 2:
        if type(rotation) != float:
	    raise TypeError('rotation has to be float for 2D input vector!')
        cpsi = cos(rotation/180.*pi)
        spsi = sin(rotation/180.*pi)
        return rotate_vector2d(r, cpsi, spsi)
    if len(r) == 3:
        if len(rotation) != 3:
	    raise TypeError('rotation has to be 3d array for 3D input vector!')
	else:
	    print("Function not implemented, yet")

def rotate_vector2d(r, cpsi, spsi):
    """
    rotate 2d or 3d vector by specified rotation

    @param r: vector
    @type r: 2d array
    @param cpsi: cos of rotation angle
    @type cpsi: float
    @param spsi: cos of rotation angle
    @type spsi: float
    @return: 2d vector
    """
    from numpy import array
    x = cpsi*r[0] - spsi*r[1]
    y = spsi*r[0] + cpsi*r[1]
    return array([x,y])

def checkEpsilon(value):    
    """
    check if value is below epsilon
    @param value: variable to be checked
    @type value: float
    """
    if abs(value) < epsilon:
        return 0.0
    else:
        return float(value)
    
def euclidianDistance(v1,v2):
    """
    euclidianDistance: Calculates euclidian distance of v1,v2
    @param v1:
    @param v2:
    @return: Euclidian distance of v1,v2  
    """

    assert len(v1) == len(v2), "euclideanDistance: lengths of v1 and v2 different"
    
    from math import sqrt
    
    sum = 0.0
    
    for i in range(len(v1)):
        sum += pow(float(v1[i]) - float(v2[i]),2)
    
    return checkEpsilon(sum**0.5)
    

def listMean(theList):    
    """
    listMean: calculates the mean of a python list
    @param theList: a list of numbers
    @return: The mean of the list
    """
    if not theList.__class__ == list:
        theList = [theList]
    
    mean = 0
    for i in range(len(theList)):
        l = theList[i]
        mean = mean + l
    
    return checkEpsilon(mean/len(theList))

def listStd(theList,mean=None):
    """
    listStd: calculates the std of a python list
    @param theList: a list of numbers
    @param mean: The mean value if available. Optional.
    @return: Standart deviation of list 
    """
    from math import sqrt
    
    if not theList.__class__ == list:
        theList = [theList]
        
    mean = mean or listMean(theList)
    
    std = 0
     
    for i in range(len(theList)):
        std = std + (theList[i]-mean)**2
        
    return checkEpsilon(sqrt(std/len(theList)))
    
def scale(v,value): 
    """
    scale: Scale a vector by a value
    @param v: vector to be scaled
    @param value: scaling factor

    last change: 10.12.12 FF: fixed bug!
    """
    assert v.__class__ == list
    assert value.__class__ in [int,float]
   
    newV = [0] * len(v)
   
    for i in range(len(v)):
        if not (v[i].__class__ in [int,float]):
            raise RuntimeError('Please provide a list of ints or floats to pytom.tools.maths.scale!' + str(v))
       
        newV[i] = checkEpsilon(v[i] * value)
    
    return newV
   
def plus(v1,v2):
    """
    plus: Add two lists
    """
    
    assert len(v1) == len(v2)
    newV = [0] * len(v1)

    for i in range(len(v1)):
        newV[i] = checkEpsilon(v1[i] + v2[i])
    
    return newV

def mult(v1,v2):
    """
    mult: Multiply two lists
    """
    
    assert len(v1) == len(v2)
    newV = [0] * len(v1)

    for i in range(len(v1)):
        newV[i] = checkEpsilon(v1[i] * v2[i])
    
    return newV

def scalar(v1,v2):
    """
    scalar: Inner product of two vectors
    """
    assert len(v1) == len(v2)
    
    scal = 0
    for i in range(len(v1)):
        scal = scal + v1[i] * v2[i]
    
    return checkEpsilon(scal)


class Matrix(object):
    """
    Matrix : Index convention according to http://mathworld.wolfram.com/Matrix.html. First index is Row index , Second index is Column index.
    @todo: Add Unit Test
    """
    
    def __init__(self,sizeX,sizeY=None):
        """
        @param sizeX: size in X
        @type sizeX: int
        @param sizeY: size in Z
        @type sizeY: int
        """
        from pytom_volume import vol
        
        if sizeX.__class__ == vol:
            self._matrix = sizeX
        elif sizeY == None:
            raise TypeError('Size x and y of the matrix must be set!')
        else:
            self._matrix = vol(sizeX,sizeY,1)
            self._matrix.setAll(0)
        
    def getMatrix(self):
        return self._matrix
    
    def getSizeX(self):
        return self._matrix.sizeX()
    
    def getSizeY(self):
        return self._matrix.sizeY()
    
    def __getitem__(self,key):
        assert key[0] < self.getSizeX()
        assert key[1] < self.getSizeY()
        
        return self._matrix(key[0],key[1],0)
        
    def __setitem__(self,key,value):
        assert key[0] < self.getSizeX()
        assert key[1] < self.getSizeY()
    
    
        self._matrix(float(value),key[0],key[1],0)
    
    def __str__(self):
        strng = ''
        for x in range(self.getSizeX()):
            strng = strng + str(self.getRow(x)) + '\n'
            
        return strng
    
    def getRow(self,index):
        
        vector = []
        
        for y in range(self.getSizeY()):
            vector.append(self[index,y])
            
        return vector
    
    def getColumn(self,index):
        
        vector = []
        
        for x in range(self.getSizeX()):
            vector.append(self[x,index])
        
        return vector
    
    def __eq__(self,otherMatrix):
        equal = True
        
        for x in range(self.getSizeX()):
            for y in range(self.getSizeY()):
                equal = equal and abs(self[x,y] - otherMatrix[x,y]) <= epsilon
                
        return equal
    
    def isSymmetric(self):
        
        symmetric = True
        
        for x in range(self.getSizeX()):
            for y in range(self.getSizeY()):
        
                if not x == y:
                    symmetric = symmetric and (abs(self[x,y] - self[y,x]) < epsilon)
        
        return symmetric
        
    def isIdentity(self):    
        
        identity = True
        
        for x in range(self.getSizeX()):
            for y in range(self.getSizeY()):
        
                if not x == y:
                    identity = identity and ( abs(self[x,y]) < epsilon)
                else:
                    identity = identity and ( abs(self[x,y] - 1) < epsilon)
                    
        return identity
        
    
    def trace(self):
        tmp = 0.0
        for i in range(self.getSizeX()):
            for j in range(self.getSizeY()):
                tmp += self[i, j]
        
        return tmp
    
    def transpose(self):
        tmp = Matrix(self.getSizeY(), self.getSizeX())
        for i in range(self.getSizeX()):
            for j in range(self.getSizeY()):
                tmp[i,j] = self[j,i]
        
        return tmp
    
    def elementwise_mult(self, m):
        if self._matrix.sizeX() != m.getSizeX() or self._matrix.sizeY() != m.getSizeY():
            raise RuntimeError('Matrices must be of same size for elementwise multiplication!')
        
        result = self._matrix * m.getMatrix()
        
        return Matrix(result)
        
    def __mul__(self,value):
        if isinstance(value, int) or value.__class__ == float:
            #multiply self by scalar
            newMatrix = Matrix(self.getSizeX(),self.getSizeY())
            
            for x in range(self.getSizeX()):
                for y in range(self.getSizeY()):
                    newMatrix[x,y] = checkEpsilon(self[x,y] * value)
            return newMatrix
        
        elif value.__class__ == list:
            #multiply by vector
            vector = value
            
            assert self.getSizeY() == len(vector)

            if self.isIdentity():
                return vector
            
            newVector = [0 for _ in range(len(vector))]
            
            for x in range(self.getSizeX()):
                row = self.getRow(x)
            
                newVector[x] = scalar(row,vector)
                 
            return newVector
        
        else:
            #multiply by matrix
            
            otherMatrix = value
            assert self._matrix.sizeX() == otherMatrix.getSizeY()
            assert self._matrix.sizeY() == otherMatrix.getSizeX()
            
            if self.isIdentity():
                return otherMatrix
            
            newMatrix = Matrix(self.getSizeX(),self.getSizeY())
            
            for x in range(self.getSizeX()):
                for y in range(self.getSizeY()):
                    
                    row = self.getRow(x)
                    column = otherMatrix.getColumn(y)
                    newMatrix[x,y] = scalar(row,column)

                    
            return newMatrix
    
class Identity(Matrix):
    """
    Identity: Identity matrix of arbitraty size
    """
    
    def __init__(self,sizeX,sizeY):
        super(Identity,self).__init__(sizeX,sizeY)

        for x in range(sizeX):
            self._matrix(1,x,x,0)
    
    def __setitem__(self,key,value):
        """
        __setitem__: This matrix is read only 
        """
        pass

    def isIdentity(self):
        """
        isIdentity: Overrides parent method
        """
        return True
    
    
class XRotationMatrix(Matrix):
    """
    XRotationMatrix:
    """
    
    def __init__(self,angle,angleInRadians=False):
        super(XRotationMatrix,self).__init__(3,3)
        
        if not angleInRadians:
            from pytom.angles.angle import deg2rad
            angle = deg2rad * angle
        
        from math import sin,cos
        self._matrix(checkEpsilon(cos(angle)),1,1,0)
        self._matrix(checkEpsilon(sin(angle)),2,1,0)
        self._matrix(checkEpsilon(cos(angle)),2,2,0)
        self._matrix(checkEpsilon(-sin(angle)),1,2,0)
        self._matrix(1,0,0,0)    
    
class YRotationMatrix(Matrix):
    """
    YRotationMatrix:
    """
    
    def __init__(self,angle,angleInRadians=False):
    
        super(YRotationMatrix,self).__init__(3,3)
        
        if not angleInRadians:
            from pytom.angles.angle import deg2rad
            angle = deg2rad * angle
        
        from math import sin,cos
        self._matrix(checkEpsilon(cos(angle)),0,0,0)
        self._matrix(checkEpsilon(-sin(angle)),2,0,0)
        self._matrix(checkEpsilon(cos(angle)),2,2,0)
        self._matrix(checkEpsilon(sin(angle)),0,2,0)
        self._matrix(1,1,1,0)
        
class ZRotationMatrix(Matrix):
    """
    ZRotationMatrix:
    """
    
    def __init__(self,angle,angleInRadians=False):
        super(ZRotationMatrix,self).__init__(3,3)
        
        if not angleInRadians:
            from pytom.angles.angle import deg2rad
            angle = deg2rad * angle
        
        from math import sin,cos
        self._matrix(checkEpsilon(cos(angle)),0,0,0)
        self._matrix(checkEpsilon(sin(angle)),1,0,0)
        self._matrix(checkEpsilon(cos(angle)),1,1,0)
        self._matrix(checkEpsilon(-sin(angle)),0,1,0)
        self._matrix(1,2,2,0)
        

def pcacov(matrix):
    """
    pcacov: Will determine eigenVector, eigenValue, eigenValueNormed for a matrix consistent to the matlab pcacov function.
    @param matrix: The matrix
    @return: [eigenVector, eigenValue, eigenVectorNormed] - eigenvectors are [nth vector, value] 
    @author: Thomas Hrabe 
    """
    
    from numpy import sum,max,diag
    from scipy.linalg import svd
    
    from pytom_volume import vol
    
    if matrix.__class__ == vol:
        from pytom_numpy import vol2npy
        matrix = vol2npy(matrix)
    
    [x,latent,coeff] = svd(matrix)
    
    totalvar = sum(diag(latent))
    
    explained = 100*latent / totalvar
    
    [sizeX,sizeY] = coeff.shape
    
    abscoeff = abs(coeff)
       
    maxIndex = abscoeff.argmax(1)

    for x in range(sizeX):
        if coeff[x,maxIndex[x]] < 0:
            coeff[x,:] = coeff[x,:] * -1;

    return [coeff,latent,explained]

    
def normalize2D(data, start=None, end=None, byRow=True):
    """
    normalize2D: normalize 2D data
    """
    from scipy import mean, std, shape
    row, col = shape(data)
    
    if byRow: # normalize by row
        if not start:
            start = 0
        if not end:
            end = row
        
        ndata = data[:start]
        for d in data[start:end]:
            m = mean(d)
            scale = std(d)
            ndata.append([(i-m)/scale for i in d])
        ndata.extend(data[end:])
    else: # normalize by column
        if not start:
            start = 0
        if not end:
            end = col
        
        ndata = []
        m = [0 for i in range(col)]
        scale = [1 for i in range(col)]
        
        for i in range(start, end):
            tmp = [d[i] for d in data]
            m[i] = mean(tmp)
            scale[i] = std(tmp)
        
        for d in data:
            nd = []
            for i in range(col):
                nd.append((d[i]-m[i])/scale[i])
            ndata.append(nd)
    
    return ndata

def rotation_distance(ang1, ang2):
    """
    rotation_distance: given two angles (lists), return the euler distance (degree).
    @param ang1: 
    @type ang1: 3-dim list or L{pytom.basic.structures.Rotation}
    @param ang2:
    @type ang2: 3-dim list or L{pytom.basic.structures.Rotation}
    @return: distance in deg
    @rtype: float
    """
    from pytom.basic.structures import Rotation
    matrix1 = Rotation(ang1).toMatrix()
    matrix2 = Rotation(ang2).toMatrix()
    res = matrix1.elementwise_mult(matrix2)
    trace = res.trace()
    
    from math import pi, acos
    temp=0.5*(trace-1.0)
    if temp >= 1.0:
        return 0.0
    if temp <= -1.0:
        return 180
    return acos(temp)*180/pi


class TransformationMatrix(Matrix):
    """
    TransformationMatrix: Class that combines any rotation and translation to one matrix (first rotation and then translation).
    """
    
    def __init__(self, rotation, translation):
        from pytom.basic.structures import Rotation, Shift
        if rotation.__class__ == Rotation:
            pass
        elif rotation.__class__ == list:
            rotation = Rotation(rotation)
        else:
            raise RuntimeError("Rotation should be of type list or pytom.basic.structures.Rotation!")
        if not isinstance(translation, (Shift, list)):
            raise RuntimeError("Translation should be of type list or pytom.basic.structures.Shift!")
        
        super(TransformationMatrix, self).__init__(4,4)
        self.setRotationCoef(rotation.toMatrix())
        self.setTranslationCoef(translation)
    
    def setRotationCoef(self, mat):
        self._matrix(mat[0,0], 0,0,0)
        self._matrix(mat[0,1], 0,1,0)
        self._matrix(mat[0,2], 0,2,0)
        self._matrix(mat[1,0], 1,0,0)
        self._matrix(mat[1,1], 1,1,0)
        self._matrix(mat[1,2], 1,2,0)
        self._matrix(mat[2,0], 2,0,0)
        self._matrix(mat[2,1], 2,1,0)
        self._matrix(mat[2,2], 2,2,0)
    
    def setTranslationCoef(self, tran):
        self._matrix(tran[0], 0,3,0)
        self._matrix(tran[1], 1,3,0)
        self._matrix(tran[2], 2,3,0)
        self._matrix(1, 3,3,0)
    
    def getRotationMatrix(self):
        mat = Matrix(3, 3)
        mat[0,0] = self._matrix(0,0,0)
        mat[0,1] = self._matrix(0,1,0)
        mat[0,2] = self._matrix(0,2,0)
        mat[1,0] = self._matrix(1,0,0)
        mat[1,1] = self._matrix(1,1,0)
        mat[1,2] = self._matrix(1,2,0)
        mat[2,0] = self._matrix(2,0,0)
        mat[2,1] = self._matrix(2,1,0)
        mat[2,2] = self._matrix(2,2,0)
        
        return mat
    
    def getTranslation(self):
        return [self._matrix(0,3,0), self._matrix(1,3,0), self._matrix(2,3,0)]
    
    def __mul__(self, otherMat):
        mat = super(TransformationMatrix, self).__mul__(otherMat)
        
        res = TransformationMatrix([0,0,0], [0,0,0])
        res._matrix = mat._matrix
        return res


def gaussian_fit(x, y, h=0.2):
    """Apply Gaussian data fitting: y=A*exp(-(x-mu)^2/(2*sigma^2)), using polyfit
    @param x: x
    @type x: list
    @param y: y
    @type y: list
    @param h: the fraction threshold of the maximum y height that the data is taken into account of fitting
    @type h: float between 0-1
    
    @return: [sigma,mu,a]
    @rtype: list
    """
    if len(x)!=len(y):
        raise RuntimeError("X and Y have different sizes!")
    if len(x)<=3:
        raise RuntimeError("Input data too short!")
    
    # cutting
    import numpy as np
    ymax = np.max(y)
    ynew = []; xnew = []
    for i in range(len(y)):
        if y[i] > ymax*h and y[i]>0.:
            ynew.append(y[i])
            xnew.append(x[i])
    
    # fitting
    from math import log, sqrt, exp
    from pylab import polyfit
    ylog = [log(a) for a in ynew]
    xlog = xnew
    
    p = polyfit(xlog,ylog,2)
    a2 = p[0]; a1=p[1]; a0=p[2]
    
    from cmath import sqrt
    sigma=sqrt(-1/(2*a2)) # can be complex numbers
    mu=a1*sigma**2
    a=exp(a0+mu**2/(2*sigma**2))
    
    return [sigma, mu, a]

