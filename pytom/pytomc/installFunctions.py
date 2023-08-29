import os
import sys

def help():
    print('PyTom compile script')
    print('')
    print('OPTIONS')
    print('--libDir             Specify multiple directories for dynamic library search')
    print('--includeDir         Specify multiple directories for header file search')
    print('--exeDir             Specify multiple directories for excecuteable search')
    print('--pythonVersion      Specify a python version you want to link against if you have many (2.5, 2.6 ...)')
    print('--target             Specify which build target you like')
    print('--minicondaEnvDir    Specify the directory in which your miniconda environment is installed. i.e. /path/to/miniconda/envs/pytom_env')
    print('                     all         : build everything')
    print('                     swig        : build swig modules only (part of all)')
    print('                     libtomc     : build C++ libraries only (part of all)')
    print('                     np          : build numpy interface (part of all)')
    print('                     smpi        : build mpi module (part of all)')
    print('                     clean       : clean all files')
    print('                     cleanswig   : clean all swig files')
    print('                     cleanlibtomc: clean C++ libraries')
    print('                     check       : check if compile was correct')
    print('')
    print('Good luck!')
    
def parseArguments(args):
    libParameter        = '--libDir'
    incParameter        = '--includeDir'
    exeParameter        = '--exeDir'
    pythonVersion       = '--pythonVersion'
    target              = '--target'
    minconEnvDir        = '--minicondaEnvDir' 
    
    libIndex    = None
    incIndex    = None
    exeIndex    = None
    pythonIndex = None
    targetIndex = None
    minconIndex = None
    
    for i in range(len(args)):
        
        if libParameter in args[i]:    
            libIndex = i
        elif incParameter in args[i]:    
            incIndex = i
        elif exeParameter in args[i]:    
            exeIndex = i
        elif pythonVersion in args[i]:    
            pythonIndex = i
        elif target in args[i]:    
            targetIndex = i
        elif minconEnvDir in args[i]:
            minconIndex = i
            
    libPaths = []
    incPaths = []
    exePaths = []

    if libIndex != None:
        for libIterator in range(libIndex+1,len(args)):
            otherKeyword = incParameter in args[libIterator] or exeParameter in args[libIterator] or pythonVersion in args[libIterator] or target in args[libIterator]
            
            if not otherKeyword:
                if os.path.exists(args[libIterator]):
                    libPaths.append(args[libIterator])
            else:
                break
    if incIndex  != None:  
        for incIterator in range(incIndex+1,len(args)):
            otherKeyword = libParameter in args[incIterator] or exeParameter in args[incIterator] or pythonVersion in args[incIterator] or target in args[incIterator]
            
            if not otherKeyword:
                if os.path.exists(args[incIterator]):
                    incPaths.append(args[incIterator])
            else:
                break    
    if exeIndex  != None:            
        for exeIterator in range(exeIndex+1,len(args)):
            otherKeyword = incParameter in args[exeIterator] or libParameter in args[exeIterator] or pythonVersion in args[exeIterator] or target in args[exeIterator]
            
            if not otherKeyword:
                if os.path.exists(args[exeIterator]):
                    exePaths.append(args[exeIterator])
            else:
                break
                
    pythonVersion = ''
    if pythonIndex  != None:     
        pythonVersion = args[pythonIndex + 1]
    else:
        import sys
        pythonVersion = f'{sys.version_info[0]}.{sys.version_info[1]}'
    
    target = ''
    if targetIndex != None:
        target = args[targetIndex + 1]

    envdir=''
    if minconIndex != None:
        libPaths.append(os.path.join(args[minconIndex+1], 'lib'))
        exePaths.append(os.path.join(args[minconIndex+1], 'bin'))
        incPaths.append(os.path.join(args[minconIndex+1], 'include'))
        print('Updating the exe, inc and lib paths using minicondaEnvDir')
        
        for f in ('c++', 'cc', 'gcc', 'g++' ):
            fname = os.path.join(args[minconIndex+1],'bin/x86_64-conda_cos6-linux-gnu-'+f)
            fnm = os.path.join(args[minconIndex+1],'bin/' + f)
            if not os.path.exists(fnm) and os.path.exists(fname):
                os.system(f'ln -s {fname} {fnm}')
            
            
        envdir = args[minconIndex+1]
                
        
    return [libPaths,incPaths,exePaths,pythonVersion,envdir, target]
    
# find the file recursively, return the dir
def find_file(name, dir, required=''):
    """
    find a file within a directory
    """
    for root, dirnames, filenames in os.walk(dir):
        if name in filenames and required in root:
            return root

    return None

def find_dir(dirname, dir):
    """
    find a file within a directory
    """
    for root, dirnames, filenames in os.walk(dir):
        
        for directory in dirnames:
            if dirname in directory:
                return root + os.sep + directory

    return None

def findDir(dirname,dirs):
    for directory in dirs:
        return_dir = find_dir(dirname,directory)
        if return_dir:
            return return_dir

def find(obj, searchDir, required=''):
    """
    findObj: find a file within a list of directories
    @param obj: One or list of objects 
    @type obj: Each obj entry is a str
    @param searchDir: List of potential directories 
    """
        
    if isinstance(obj,str):
        obj = [obj]
        
    returnValue = None,None
    for o in obj:
        dir = findObj(o,searchDir, required=required)
        
        if dir is not None:
            returnValue = [o,dir]

    return returnValue

def findObj(obj, searchDir, required=''):
    """
    findObj: find a file within a list of directories
    @param obj: One object 
    @type obj: str
    @param searchDir: List of potential directories 
    """
    
    for dir in searchDir:
        parentDir = find_file(obj, dir, required=required)
        if parentDir:
            print('Searching : ' , obj, '\t\t Found : ', True)
            return parentDir
    else:
        print('Searching : ' , obj, '\t\t Found : ', False)
        return None

def adjustLibraryVersions(library,versionList,flag,extension, search_dir):
    """
    find a library and determine its proper version. return according library path and linker flag
    """
    lib = None
    libraryFlag = flag
    versionIncrement = 0
    while not lib:
      
        version = versionList[versionIncrement]
        
        lib = find(library + version + extension,search_dir)
        libraryFlag = flag + version
        versionIncrement += 1
    
    return [lib,libraryFlag]
    
def generateExecuteables(libPaths=None,binPaths=None,pyPaths=None,python_version=''):
    
    os.chdir('..')
    pytomDirectory = os.getcwd()
    os.chdir(pytomDirectory + os.sep + 'bin')
    os.chdir(pytomDirectory + os.sep + 'pytomc')
    generatePyTomScript(pytomDirectory, python_version)
    generateIPyTomScript(pytomDirectory)
    generatePyTomGuiScript(pytomDirectory, python_version)

def generatePyTomScript(pytomDirectory,python_version):
    
    pytomCommand = f"""#!/usr/bin/env bash
cat {sys.prefix}/pytom_data/LICENSE.txt

FID=0
export PYTOM_GPU=-1

for a in $*
do
    if [ $FID -eq 1 ]; then
        FID=0
        export PYTOM_GPU="$a"
    fi

    if [ "$a" = '-g' ] || [ "$a" = '--gpuID' ]; then
        export PYTOM_GPU=1
        FID=1
    fi
done

{sys.executable} -O $*
"""

    f = open(pytomDirectory + os.sep + 'bin' + os.sep + 'pytom','w')
    f.write(pytomCommand)
    f.close()
    
    os.system('chmod 755 ' + pytomDirectory + os.sep + 'bin' + os.sep + 'pytom')

def generatePyTomGuiScript(pytomDirectory, python_version):
    pytomguiCommand = '# !/bin/bash\n'
    pytomguiCommand += f'{sys.executable} {pytomDirectory}/gui/pytomGUI.py $1\n'

    f = open(pytomDirectory + os.sep + 'bin' + os.sep + 'pytomGUI', 'w')
    f.write(pytomguiCommand)
    f.close()

    os.system('chmod 755 ' + pytomDirectory + os.sep + 'bin' + os.sep + 'pytomGUI')

def generateIPyTomScript(pytomDirectory):

    ipytomCommand = '#!/usr/bin/env bash\n'
    ipytomCommand += f'cat {sys.prefix}/pytom_data/LICENSE.txt\n'
    ipytomCommand += 'ipython $* -i\n'
    
    f = open(pytomDirectory + os.sep + 'bin' + os.sep + 'ipytom','w')
    f.write(ipytomCommand)
    f.close()
    
    os.system('chmod 755 ' + pytomDirectory + os.sep + 'bin' + os.sep + 'ipytom')

    
def checkCompileDirsExist():
    
    import os.path
    
    paths =  ['./libs/libtomc/libs','./swigCpp','./swigLibs']
    
    for path in paths:
        if not os.path.isdir(path):
            os.mkdir(path)

def _readStringFile(filename):
    """
    readStringFile: Helper function. Will read in a string file and return a string
    @param filename: A filename
    @return: String 
    """

    lines = '' 

    f = open(filename)
    try:
        for line in f:
            lines = lines + line
    except :
        print('Error reading ' + filename + '!')
        assert False
    finally:
        f.close()

    return lines

def checkSwigVersion():
    """
    checkSwigVersion: Returns true if swig version is below 3.0.12
    """
    from os import system
    
    system('swig -version | grep Version > swigVersion.txt')
    txt = _readStringFile('swigVersion.txt')
    system('rm swigVersion.txt')
    
    txts = txt.split(' ')
    versionString = txts[-1].split('\n')[0]
    
    return versionString < '3.0.12'


def checkPythonVersion():
    """
    checkSwigVersion: Returns true if swig version is below 3.0.12
    """
    import sys
    versionString = f'{sys.version_info[0]}.{sys.version_info[1]}'

    return versionString < '3.0.12'
