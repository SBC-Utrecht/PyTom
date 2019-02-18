"""
functions for Constrained PCA
"""
import numpy as np

def subTomoClust(particleListFilename, classifiedParticleListFilename, 
        cccName, neig, nclass, verbose=False):
    """
    subtomogram clustering using CCC and kmeans
    @param particleListFilename: particle list filename
    @type particleListFilename: str
    @param classifiedParticleListFilename: output particle list filename
    @type classifiedParticleListFilename: str
    @param cccName: Name of (Constrained) Correlation Coefficient (CCC) matrix
    @type cccName: str
    @param neig: number of eigenvectors
    @type neig: int
    @param nclass: number of classes
    @type nclass: int

    @author: FF Jan 2013
    """
    from pytom.basic.structures import ParticleList
    pl = ParticleList('.')
    pl.fromXMLFile(particleListFilename)
    if verbose:
        print("Particle List read in")
    ccc = readCCC(cccName)
    if verbose:
        print("CCC read in")
    coeffs, eigvalues = SVD_analysis(ccc)
    if verbose:
        print("Eigen analysis done")
    labels = kmeansCluster(coeff=coeffs, neig=neig, nclass=nclass)
    if verbose:
        print("kmeans clustering done")
    for (ipart,part) in enumerate(pl):
        part.setClass(className=str(labels[ipart]))
    if verbose:
        print("Class labels assigned")
    pl.toXMLFile(classifiedParticleListFilename)
    if verbose:
        print("File written")

def readCCC(CCCName):
    """
    read (Constrained) Correlation Coefficient (CCC) matrix
    @param CCCName: Name of file
    @type CCCName: str

    @return: CCC
    @rtype: L{numpy.array}

    @author: FF Jan 2013
    """
    ccc = np.genfromtxt(CCCName, delimiter=',', skip_header=0)
    return ccc

def SVD_analysis(ccc):
    """
    compute eigenvectors of CCC matrix
    @param ccc: (Constrained) Correlation Coefficient (CCC) matrix
    @type ccc: L{numpy.array}

    @return: coefficients, eigenvalues
    @author: FF Jan 2013
    @rtype: L{numpy.array}, L{numpy.array}
    """
    from math import sqrt
    coeff, eigvalues, v = np.linalg.svd(ccc)
    #rescale vectors so they include variance
    for ii in range(0,len(eigvalues)):
        coeff[:,ii] = sqrt(abs(eigvalues[ii]))*coeff[:,ii]
    return coeff, eigvalues

def kmeansCluster(coeff, neig, nclass):
    """
    cluster eigen coefficients using k-means clustering
    @param coeff: coefficients of eigen analysis
    @type coeff: L{numpy.array}
    @param neig: number of eigenvectors
    @type neig: int
    @param nclass: number of classes
    @type nclass: int

    @return: labels
    @rtype: L{numpy.array}

    @author: FF Jan 2013
    """
    import scipy.cluster
    centroids,labels = scipy.cluster.vq.kmeans2(data=coeff[:,0:neig],k=nclass,iter=50)
    return labels

def averageClasses(particleListFilename, avName):
    """
    write class averages of classified particles
    @param particleListFilename: particle list filename
    @type particleListFilename: str
    @param avName: Name for class averages (<avName>_iclass.em)
    @type avName: str

    @author: FF
    @date: Jan 2013
    """
    from pytom.basic.structures import ParticleList
    pl = ParticleList()
    pl.fromXMLFile(particleListFilename)
    pl.sortByClassLabel()
    pls = pl.splitByClass()

    for cl in pls:
        className = cl[0].getClassName()
	cl.average(avName + "_" + str(className) + '.em')
	print(className, ' contains ' , len(cl) , ' particles') 


