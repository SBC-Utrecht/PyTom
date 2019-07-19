'''
Created on Apr 15, 2010

@author: hrabe
'''

import os

_proc_status = '/proc/%d/status' % os.getpid()

_scale = {'kB': 1024.0, 'mB': 1024.0*1024.0,
          'KB': 1024.0, 'MB': 1024.0*1024.0}

def _VmB(VmKey):
    '''Private.
    '''
    global _proc_status, _scale
     # get pseudo file  /proc/<pid>/status
    try:
        t = open(_proc_status)
        v = t.read()
        t.close()
    except:
        return 0.0  # non-Linux?
    #print v# get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
    i = v.index(VmKey)
    v = v[i:].split(None, 3)  # whitespace
    if len(v) < 3:
        return 0.0  # invalid format?
     # convert Vm value to bytes
    return float(v[1]) * _scale[v[2]]


def memory(since=0.0):
    '''Return memory usage in bytes.
    '''
    return _VmB('VmSize:') - since


def resident(since=0.0):
    '''Return resident memory usage in bytes.
    '''
    return _VmB('VmRSS:') - since


def stacksize(since=0.0):
    '''Return stack size in bytes.
    '''
    return _VmB('VmStk:') - since

def printProcessStats():
    """
    printProcessStats: Prints current memory usage of this pytom process to screen.
    If in mpi mode, will display mpi_id, too
    """
    import os, pytom_mpi
    
    if pytom_mpi.isInitialised():
        print 'Usage of MPI-ID (',pytom_mpi.rank(),')'
        
    pid = os.getpid()
    
    psCommand = 'ps -v -p ' + str(pid)  
    os.system(psCommand)
    #print memory()

def cleanUp():
    from pytom.tools.memory import deleteWedgeStorage
    deleteWedgeStorage()

class Singleton(object):
    """
    Singleton: Abstract class used for singleton classes. Most children will be used for efficient memory management to speed up code.
    @deprecated: Use python global class attribute instead such as for volume storage! 
    """
    
    def __new__(cls, *a, **k):
        if not hasattr(cls, '_inst'):
            cls._inst = super(Singleton, cls).__new__(cls, *a, **k)
            
        return cls._inst

class VolumeStorage(object):
    """
    VolumeStorage: volumes will store volume once read from disk in memory. 
    This class is quite useful when reading and binning is repeated many times. Use for speed up! 
    """
    
    _volumeList = {}
    
    def append(self,volume,key):
        self._volumeList[key] = volume
    
    def get(self,key):
        """
        get: return a volume volume 
        """
        
        volume = self._volumeList[key]
        
        from pytom_volume import vol
        
        returnVolume = vol(volume.sizeX(),volume.sizeY(),volume.sizeZ())
        
        returnVolume.copyVolume(volume)
        
        return returnVolume
    
        
        

def read(filename,subregionX=0,subregionY=0,subregionZ=0,subregionXL=0,subregionYL=0,subregionZL=0,samplingX=0,samplingY=0,samplingZ=0,binningX=0,binningY=0,binningZ=0):
    """
    read: Will read a file either from disk, or copy the identical volume if file has already been read from disk. -> no disk access for data. Saves a lot of time if binning is set! Consumes a lot of memory on the other hand!!!
    @param filename: Filename of volume file. Will be used as key for later storing. 
    @param subregionX: Subregion start in X (default is 0)
    @param subregionY: Subregion start in Y (default is 0)
    @param subregionZ: Subregion start in Z (default is 0)
    @param subregionXL: Subregion X length in pixels (default is 0)
    @param subregionYL: Subregion Y length in pixels (default is 0)
    @param subregionZL: Subregion Z length in pixels (default is 0)
    @param samplingX: Sampling in X (default is 0)
    @param samplingY: Sampling in Y (default is 0)
    @param samplingZ: Sampling in Z (default is 0)
    @param binningX: Binning in X (default is 0)
    @param binningY: Binning in Y (default is 0)
    @param binningZ: Binning in Z (default is 0)
    """
    #import time
    #print time.time()
    volumeStorage = VolumeStorage()
    key = filename + str(subregionX) + str(subregionY) + str(subregionZ) + str(subregionXL) + str(subregionYL) + str(subregionZL) + str(samplingX) + str(samplingY) + str(samplingZ) + str(binningX) + str(binningY) + str(binningZ)
    
    volume = None
    
    try:
        volume = volumeStorage.get(key)
    except:
    
        from pytom_volume import read as readFromDisk
        
        volume = readFromDisk(filename,int(subregionX),int(subregionY),int(subregionZ),
                                       int(subregionXL),int(subregionYL),int(subregionZL),
                                       int(samplingX),int(samplingY),int(samplingZ),
                                       int(binningX),int(binningY),int(binningZ))
        
        volumeStorage.append(volume,key)
    
    return volume
        
        
class WedgeStorage(VolumeStorage):
    """
    WedgeStorage: A storage used for wedge volumes only. Subclass of volume storage.
    """
    
    
class FilterStorage(VolumeStorage):
    """
    FilterStorage: A storage used for filters. 
    """
    
    