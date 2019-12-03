'''
Created on Dec 22, 2010

@author: hrabe
'''



def randomiseParticleListClasses(particleList, numberClasses):
    """
    randomiseParticleListClasses: Randomises class membership of particles in a ParticleList
    @param particleList:
    @param numberClasses:  
    @author: Thomas Hrabe
    """
    
    import random
    
    random.seed()
    
    for i in range(len(particleList)):
        
        particle = particleList[i]
    
        particle.setClass(int(random.random()*numberClasses))
        
        particleList[i] = particle
        
    return particleList


def determineClassSwaps(particleList1,particleList2):
    """
    determineClassSwaps
    @param particleList1: The previous particleList
    @param particleList2:  The current particleList
    @return: List of particles that swapped class and a corresponding list of [previous class, current class].
    @rtype: [L{pytom.basic.structures.ParticleList},[previousClass,currentClass]]
    """
    from pytom.basic.structures import ParticleList
    
    returnList = ParticleList(particleList1.getDirectory())
    classSwaps = []
    
    for particle in particleList1:
        
        particleName = particle.getFilename()
        particleClass = particle.getClassName()
        
        otherParticle = particleList2.getParticleByFilename(particleName)
        otherParticleClass = otherParticle.getClassName()
        
        if not particleClass == otherParticleClass:
            returnList.append(otherParticle)
            
        classSwaps.append([particleClass,otherParticleClass])
            
    return [returnList,classSwaps]
        

    
    