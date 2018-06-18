'''
Created on Feb 1, 2011

@author: hrabe
'''

def createClassificationResultDictionaries(classifiedParticleList,groundTruthParticleList,verbose=False):
    """
    assessClassification: Comment in LB 31.1.2011
    @param classifiedParticleList: list of classified particles
    @param groundTruthParticleList: ground truth
    @param verbose: If True, will print the class dictionaries generated. Default is false.
    @return:  A dictionary that maps new classnames to groundTruthClasses
    """
    from pytom.basic.structures import Particle,ParticleList
    
    #from file if filename
    if classifiedParticleList.__class__ == str:
        classifiedParticleListFile = classifiedParticleList
        
        classifiedParticleList = ParticleList('/')
        classifiedParticleList.fromXMLFile(classifiedParticleListFile)
    
    #from file if filename
    if groundTruthParticleList.__class__ == str:
        
        groundTruthParticleList = ParticleList('/')
        groundTruthParticleList.fromXMLFile(groundTruthParticleList)
    
    newClassesToGroundTruthMap = {} #maps new class name to a ground truth class name 
    gtClassNames = {} #maps if one class has been determined 
    
    for gtClass in groundTruthParticleList.splitByClass():
        particle = gtClass[0]
        gtClassNames[particle.getClassName()] = False
        
    if verbose:
        print 'gtClassNames : ',gtClassNames
    
    groundTruthParticleListXML = groundTruthParticleList.toXML()
    
    for newClass in classifiedParticleList.splitByClass():
        
        if verbose:
            print ''
            print 'newClass : ', newClass
            
        particlesFromGroundTruth = ParticleList('/')
        
        
        for particle in newClass:
            #collect all particles that were assigned to newClass
            particleName = particle.getFilename()
            #gtParticle = groundTruthParticleList.getParticleByFilename(particleName)
            particleXML = groundTruthParticleListXML.xpath('/ParticleList/Particle[@Filename="'+str(particleName)+'"]')
            gtParticle = Particle('a')
            gtParticle.fromXML(particleXML[0])
            
            particlesFromGroundTruth.append(gtParticle)
    
        if verbose:
            print 'len(particlesFromGroundTruth) : ',len(particlesFromGroundTruth) 
        
        #sort classes according to size to descending order
        sortedClasses = sorted(particlesFromGroundTruth.splitByClass(),key = lambda x: len(x), reverse=True)
        
        classWasAssigned = False
        classIndex = 0
        
        if verbose:
            print 'len(sortedClasses) : ', len(sortedClasses)
            
        while not classWasAssigned and classIndex < len(sortedClasses):
            sortedClass = sortedClasses[classIndex]            
            className = sortedClass[0].getClassName()
            
            if verbose:
                print 'className : ' + className
                print 'len(sortedClass) : ' , len(sortedClass)
                
            classWasAssigned = not gtClassNames[className]
            
            if verbose:
                print 'classWasAssigned : ', classWasAssigned

            if not classWasAssigned:
                classIndex = classIndex + 1
            else:
                gtClassNames[className] = True
                newClassesToGroundTruthMap[newClass[0].getClassName()] = className
                
                if verbose:
                    print 'gtClassNames : ', gtClassNames
                    print 'newClassesToGroundTruthMap : ' , newClassesToGroundTruthMap

    return newClassesToGroundTruthMap

def assessClassification(classifiedParticleList,groundTruthParticleList,verbose=False):
    """
    assessClassification: Comment in LB 31.1.2011
    @param classifiedParticleList: list of classified particles
    @param groundTruthParticleList: ground truth
    @param verbose: If True, will print the class dictionaries generated. Default is false.
    @return:  [trueHits,falseHits,trueHitsPercent,falseHitsPercent,numberClusters,[clusterSizes]]
    """
    from pytom.basic.structures import Particle,ParticleList
    
    #from file if filename
    if classifiedParticleList.__class__ == str:
        classifiedParticleListFile = classifiedParticleList
        
        classifiedParticleList = ParticleList('/')
        classifiedParticleList.fromXMLFile(classifiedParticleListFile)
    
    #from file if filename
    if groundTruthParticleList.__class__ == str:
        groundTruthParticleListFile = groundTruthParticleList
        
        groundTruthParticleList = ParticleList('/')
        groundTruthParticleList.fromXMLFile(groundTruthParticleListFile)
     
    gtClassNamesAssigned = {} #maps if one class has been determined 
    
    for gtClass in groundTruthParticleList.splitByClass():
        particle = gtClass[0]
        gtClassNamesAssigned[particle.getClassName()] = False
        
    
    if verbose:
        print 'GT Classes ',gtClassNamesAssigned    
        
    
    gtClassesPerClass = {}
    
    newClasses = classifiedParticleList.splitByClass()
    newClassNamesAssigned = {}
    numberClasses = len(newClasses)
    classSizes = []
    
    for i in xrange(len(newClasses)):
        classSizes.append(len(newClasses[i]))
    
    for newClass in newClasses:
    
        newClassParticleList = ParticleList(newClass.getDirectory())
        newClassNamesAssigned[newClass[0].getClassName()] = False
        
        for particle in newClass:
            pp = groundTruthParticleList[particle.getFilename()]
            newClassParticleList.append(pp)
            
        gtClassSizeDictionary = {}
        
        for gtClass in newClassParticleList.splitByClass():
            particle = gtClass[0]
            gtParticle = groundTruthParticleList[particle.getFilename()]
            gtClassSizeDictionary[gtParticle.getClassName()] = len(gtClass)
            
        gtClassesPerClass[newClass[0].getClassName()] = [newClassParticleList,gtClassSizeDictionary]
    
    if verbose:
        print 'Class distribution dictionary'
        for k in gtClassesPerClass.keys():
            print k,gtClassesPerClass[k]
     
    
    gtToClassDictionary = {}
    
    for gtName in gtClassNamesAssigned.keys():
    
        newClassIndex = 0
        classSizeList = []
        
        largestClass  = -1
        maxClassName = 'unknown'
        assigned = False
        for newClassName in gtClassesPerClass.keys():
        
            l = gtClassesPerClass[newClassName]
            gtClassSizeDictionary = l[1]
            if verbose:
                print 'GT Name',gtName,' New Class Name',newClassName
                print 'GT Name Size',gtClassSizeDictionary
            
            try:
                if verbose:
                    print gtClassSizeDictionary[gtName]
                if largestClass < gtClassSizeDictionary[gtName] and not newClassNamesAssigned[newClassName]:
                    largestClass = gtClassSizeDictionary[gtName]
                    maxClassName = newClassName
                    if verbose:
                        print 'SWAP'
                
            except KeyError:
                pass
            
        gtToClassDictionary[gtName] = maxClassName
        newClassNamesAssigned[maxClassName] = True
        gtClassNamesAssigned[gtName] = True            
        
        for newClassName in gtClassesPerClass.keys():
            try:
                l = gtClassesPerClass[newClassName]
                gtClassSizeDictionary = l[1]
            
                del gtClassSizeDictionary[gtName]
                l[1] = gtClassSizeDictionary
            
                gtClassesPerClass[newClassName] = l
            except KeyError:
                pass
    if verbose:            
        print 'GT to New Dictionary'            
        print gtToClassDictionary
    
    
    trueHits = 0
    falseHits = 0
    
    classifiedParticleListXML = classifiedParticleList.toXML()
    
    for gtParticle in groundTruthParticleList:
        
        particleXML = classifiedParticleListXML.xpath('/ParticleList/Particle[@Filename="'+str(gtParticle.getFilename())+'"]')
        particle = Particle('a')
        particle.fromXML(particleXML[0])

        
        if particle.getClassName() == gtToClassDictionary[gtParticle.getClassName()]:
            trueHits += 1
        else:
            falseHits += 1
        
    return  [trueHits,falseHits,float(trueHits)/len(classifiedParticleList),float(falseHits)/len(classifiedParticleList),numberClasses,classSizes]


