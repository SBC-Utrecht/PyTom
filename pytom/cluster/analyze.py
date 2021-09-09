'''
Created on Jan 11, 2010

@author: hrabe
'''

def determineRotationsByClassAssignment(classifiedParticleList,rotationsParticleList):
    """
    determineRotationsByClassAssignment: Will determine rotations (stored in rotationsParticleList) \
    to clusters determined in classifiedParticleList

    @param classifiedParticleList:
    @param rotationsParticleList:
    @return: [rotationsParticleList clustered according to classifiedParticleList, \
    list of L{pytom.angles.angleList.AngleList} storing the only cluster rotations]

    @author: Thomas Hrabe     
    """
    
    from pytom.angles.angleList import AngleList
    
    rotationsCopy = rotationsParticleList.copy()
    rotationsCopy.setClassesFromList(classifiedParticleList)
    
    classLists = rotationsCopy.splitByClass()
    
    clusteredRotations = []
    
    for classMembers in classLists:
        
        angleList = AngleList()
        
        for particle in classMembers:
            angleList.append(particle.getRotation())
        
        clusteredRotations.append(angleList)
        
    return [classLists,clusteredRotations]
    

def extractClusterRotationsFromParticleList(particleList):
    """
    extractClusterRotationsFromParticleList: Splits a particleList according to the clusters stored
    @param particleList: The particleList 
    @type particleList: either str or L{pytom.alignment.structures.ParticleList}. 
    @return: Returns a list of  L{pytom.angles.angleList.AngleList} storing cluster rotations 
    """
    from pytom.angles.angleList import AngleList
    
    particleClasses = particleList.splitByClass()
    
    classes = []
    
    for i in range(len(particleClasses)):
        className = particleClasses[i][0].getClassName()
        classes.append(className)
    
    classRotations = []
    
    for className in classes:
        
        particles = particleList.particleFromClass(className)
        angleList = AngleList()
        
        for particle in particles:    
            angleList.append(particle.getRotation())
        
        classRotations.append(angleList)    
        
    return classRotations
 


    
    
    
    
    
