'''
Created on Mar 16, 2011

@author: luiskuhn
'''

def getTuples(xmippAngleList, transformList, lineOffset):
    """
    getTuples: Constrain parameters into correct xmipp ranges  
    @return: xmipp doc lines needed for xmipp reconstruction
    @author: luiskuhn
    """
    from pytom.tools.maths import epsilon
    
    tmp = []
    
    for i in range(len(xmippAngleList)):
        
        angles = xmippAngleList[i]
        
        flip = 0
        rot = angles[0]
        psi = 0
        tilt = 0
        
        
        if angles[2] < epsilon:
            psi = 360 + angles[2]
        else:
            psi = angles[2]
        
        if angles[1] > 90.0+epsilon:
            flip = 1
            rot = -(180 - angles[0])
            tilt = 180 - angles[1]
            psi = 180 - psi
        else:
            tilt = angles[1]
        
        if rot < -180.0-epsilon:
            rot =  360 + rot
        if rot > 180.0+epsilon:
            rot = rot - 360
        
        if psi < epsilon:
            psi = 360 + psi
            
        shiftx = transformList[i][1][0]
        shifty = transformList[i][1][1]
        
        filename = transformList[i][2]
        
        tmp.append([filename, [lineOffset + i, 8, rot, tilt, psi, shiftx, shifty, -1, flip, 0]])
        
    return tmp
        
        
      
        

def checkAngles(angleList, flag,verbose = False):
    """
    checkAngles:
    @author: luiskuhn
    """ 
    if verbose:
        print 'not really checking angles, yet'
    
        

def getParticleTransformLines(particle, projectionList, lineOffset):
    """
    getParticleTransformLines: Weird function name indeed, converts particle transformation parameters (rotation&shift) into xmipp convention 
    @param particle: The particle
    @param projectionList: A projection list (small projection per particle)
    @param lineOffset: Some kind of index offset used
    @author: luiskuhn
    """
    from pytom.basic.structures import Rotation
    from pytom.reconstruction.reconstructionStructures import ProjectionList
    from pytom.angles.angleFnc import zxzToZYZ, pointRotateZXZ
    
      
    alignmentRotation = particle.getRotation()
    
    invertedAlignmentRotation = alignmentRotation.invert()
    
    #print invertedAlignmentRotation
    alignmentShift = particle.getShift()
    
    transformList = []
    
    for projection in projectionList:
        
        theta = projection.getTiltAngle()      
        offsetX = projection.getOffsetX()
        offsetY = projection.getOffsetY()
        
        projectionAlignmentRotation = projection.getAlignmentRotation() 
        
        #projectionRotation = Rotation(270 - projectionAlignmentRotation , 90, theta)
        projectionRotation = Rotation(270 , 90, theta) 
        rotationSum = invertedAlignmentRotation + projectionRotation        
        
        invertedRotationSum = rotationSum.invert()       
        #rotatedShift = pointRotateZXZ([-alignmentShift.getX(), -alignmentShift.getY(), -alignmentShift.getZ()], -90, -(270 - projectionAlignmentRotation), -theta)
        rotatedShift = pointRotateZXZ([-alignmentShift.getX(), -alignmentShift.getY(), -alignmentShift.getZ()], -90, -270 , -theta)
        #projectedShift = [round(rotatedShift[0]), round(rotatedShift[1])]
        projectedShift = [rotatedShift[0]-offsetX, rotatedShift[1]-offsetY]
        filename = projection.getFilename()
                
        transformList.append([invertedRotationSum, projectedShift, filename])
        
    xmippAngleList = []
    
    for i in range(len(transformList)):
        
        invertedRotationSum = transformList[i][0]               
        
        xmippAngleList.append( zxzToZYZ(invertedRotationSum.getZ1(), invertedRotationSum.getZ2(), invertedRotationSum.getX()) )
        
    
    ##print xmippAngleList
    
    checkAngles(xmippAngleList, True)
    
    tuples = getTuples(xmippAngleList, transformList, lineOffset)
        
    return tuples


def setReference(tupleList):
    """
    setReference: Do something (???)
    @author: luiskuhn
    """
    refMap = {}
    
    refCount = 1
    
    for i in range(len(tupleList)):
        
        tuple = tupleList[i]
        
        rot = round(tuple[1][2], 15)
        tilt = round(tuple[1][3], 15)
        strRot = "%10.15f" % rot
        strTilt = "%10.15f" % tilt
        key = strRot + '_' + strTilt
                
        if refMap.has_key(key):
            
            tuple[1][7] = refMap[key]
        else:
            
            tuple[1][7] = refCount
            refMap[key] = refCount
            
            refCount = refCount + 1            
            
            
            

def setupXMIPP(particleList, projectionLists, docFilename, libFilename):
    """
    setupXMIPP: Setup xmipp doc and lib files for xmipp reconstruction
    @param particleList: A list of particles
    @param projectionLists: A list of projectionsLists
    @param docFilename: The resulting xmipp doc
    @param libFilename: The resulting xmipp lib
    @author: luiskuhn
    """
    
    tupleList = []
    lineOffset = 1
    
    for i in range(len(particleList)):
        
        particle = particleList[i]
        projectionList = projectionLists[i]
        
        particleTupleList = getParticleTransformLines(particle, projectionList, lineOffset)
        
        lineOffset = lineOffset + len(particleTupleList)
        
        tupleList.extend(particleTupleList)
        
    setReference(tupleList)
       
    writeDoc(docFilename, tupleList)
    writeLib(libFilename, tupleList)
    
        
def copyProjectionsForXMIPP(particleList, projectionLists, xmippProjectionDirectory,target,showProgressBar = False,verbose=False):
    """
    
    """
    from pytom.tools.ProgressBar import FixedProgBar
    from pytom.tools.files import readSpider,writeSpider
    import os 
        
    
    projectionPropertiesList = []
    lineOffset = 1
    
    for i in range(len(particleList)):
        
        particle = particleList[i]

        projectionList = projectionLists[i]
        
        particleTupleList = getParticleTransformLines(particle, projectionList, lineOffset)
        
        lineOffset = lineOffset + len(particleTupleList)
        
        projectionPropertiesList.extend(particleTupleList)
    
    if showProgressBar:
        progressBar = FixedProgBar(0,len(projectionPropertiesList),'Projections modified ')
        progressBar.update(0)
    
    for i in xrange(len(projectionPropertiesList)):
        
        projectionName = projectionPropertiesList[i][0]
        
        filename = projectionName.rsplit('/')
        filename = filename[len(filename)-1]
        
        newFilename = xmippProjectionDirectory + os.sep + filename[0:len(filename)-4]+ '_P' +str(i) + '.spi'
        
        projection = readSpider(projectionName)
        
        #print projectionPropertiesList[i][1][2], projectionPropertiesList[i][1][3], projectionPropertiesList[i][1][4], projectionPropertiesList[i][1][5], projectionPropertiesList[i][1][6], 0
        writeSpider(projection,newFilename, projectionPropertiesList[i][1][2], projectionPropertiesList[i][1][3], projectionPropertiesList[i][1][4], projectionPropertiesList[i][1][5], projectionPropertiesList[i][1][6], 0)
        
        if showProgressBar:
            progressBar.update(i)
        

    cmd = 'ls ' + xmippProjectionDirectory + '/*.spi | awk \'{print $1 \" 1\"}\' > ' + xmippProjectionDirectory + '/reconstruction.sel'
    os.system(cmd)

    cmd = 'xmipp_reconstruct_fourier -i ' + xmippProjectionDirectory + '/reconstruction.sel -o ' + xmippProjectionDirectory + '/result.out -sym C1 -thr 1  -pad_vol 3'
    os.system(cmd)
    
    v = readSpider(xmippProjectionDirectory + '/result.out')
    v.write(target)
    
def writeDoc(filename, lines):
    """
    writeDoc: Writes xmipp doc file required for xmipp fourier reconstruction
    @author: luiskuhn
    """
    nameFormat = ' ; %s\n'
    valuesFormat = '%5d %d  %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n'
    
    header = ' ; Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Ref (6), Flip (7), maxCC (8)\n'
    
    file = open(filename, 'w')
    
    file.write(header)
    for i in range(len(lines)):
        
        name = lines[i][0]
        values = tuple(lines[i][1])
        
        nameString = nameFormat % name
        valuesString = valuesFormat % values
        
        file.write(nameString)
        file.write(valuesString)
        
    file.close()


def writeLib(filename, lines):
    """
    writeLib: Writes xmipp lib file required for xmipp fourier reconstruction
    @author: luiskuhn
    """     
    valuesFormat = '%5d %d  %10.5f %10.5f %10.5f\n'
        
    file = open(filename, 'w')
        
    for i in range(len(lines)):
                
        values = lines[i][1]
        valuesTuple = (i+1, 3, values[2], values[3], 0)
        
        valuesString = valuesFormat % valuesTuple
       
        file.write(valuesString)
        
    file.close()
    
        
def xmippFourierReconstruction(resultEMFile, docFile = 'xmipp.doc', libFile = 'xmipp.lib', tempDir = 'xmippTempDir'):
    """
    xmippFourierReconstruction: Calls xmipp fourier reconstruction
    @param resultEMFile: The filename of the reconstructed tomogram
    @param docFile: The required xmipp doc file 
    @param libFile: The required xmipp lib file
    @param tempDir: Temporary directory where xmipp can work in (will be deleted automatically after xmipp is done)
    @author: luiskuhn
    """
    import os
    from pytom.tools.files import readSpider
    
    dir = tempDir
    cmd = 'mkdir ' + dir
    os.system(cmd)
    cmd = 'chmod ugo+rwx ' + dir
    
    dir = tempDir + '/logs/'
    cmd = 'mkdir ' + dir
    os.system(cmd)
    cmd = 'chmod ugo+rwx ' + dir
    os.system(cmd)    
    
    dir = tempDir + '/classes/'
    cmd = 'mkdir ' + dir
    os.system(cmd)
    cmd = 'chmod ugo+rwx ' + dir
    os.system(cmd)
    
    dir = tempDir + '/vols/'
    cmd = 'mkdir ' + dir
    os.system(cmd)
    cmd = 'chmod ugo+rwx ' + dir
    os.system(cmd)
    
    dir = tempDir + '/classes/ProjMatchClasses/'
    cmd = 'mkdir ' + dir
    os.system(cmd)
    cmd = 'chmod ugo+rwx ' + dir
    os.system(cmd)
    
    cmd = 'xmipp_angular_class_average -i ' + docFile + ' -lib ' + libFile + ' -dont_write_selfiles -o ' + tempDir + '/classes/ProjMatchClasses/proj_match'
    print cmd
    os.system(cmd)
    
    cmd = 'ls ' + tempDir + '/classes/ProjMatchClasses/proj_match_class*.xmp | awk \'{print $1 \" 1\"}\' > ' + tempDir + '/classes/ProjMatchClasses/reconstruction.sel'
    print cmd
    os.system(cmd)

    cmd = 'xmipp_reconstruct_fourier -i ' + tempDir + '/classes/ProjMatchClasses/reconstruction.sel -o ' + tempDir + '/vols/out.spi -thr 1 -weight'
#    cmd = 'xmipp_reconstruct_wbp -i ' + tempDir + '/classes/ProjMatchClasses/reconstruction.sel -o ' + tempDir + '/vols/out.spi  -weight'
    print cmd
    os.system(cmd)
    
    out = readSpider(tempDir + '/vols/out.spi')
    out.write(resultEMFile)

    cmd = 'rm -rf ' + tempDir
    os.system(cmd)

    cmd = 'rm ' + docFile
    os.system(cmd)

    cmd = 'rm ' + libFile
    os.system(cmd)

def callXMIPPReconstruction(resultEMFile, docFile = 'xmipp.doc', libFile = 'xmipp.lib', tempDir = 'xmippTempDir', reconstructionType='fourier'):
    """
    callXMIPPReconstruction:
    @param resultEMFile: The filename of the reconstructed tomogram
    @param docFile: The required xmipp doc file 
    @param libFile: The required xmipp lib file
    @param tempDir: Temporary directory where xmipp can work in (will be deleted automatically after xmipp is done)
    @param reconstructionType:    
    @author: Thomas Hrabe
    """
    
    if reconstructionType == 'fourier':
        xmippFourierReconstruction(resultEMFile, docFile, libFile, tempDir)
    else:
        raise RuntimeError('The reconstruction method specified is not supported yet!')
    