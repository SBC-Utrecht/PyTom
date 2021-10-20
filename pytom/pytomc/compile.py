import os
import sys
import platform

from installFunctions import *

if checkSwigVersion():
    print('Your Swig Version is less than 3.0.12.')
    print('Please install 3.0.12 or newer.')
    print('Compilation of PyTom will fail otherwise!')
    sys.exit(1)

if checkPythonVersion():
    print('Your Python Version is less than 3.6.')
    print('Please install 3.6 or newer.')
    print('Compilation of PyTom will fail otherwise!')
    sys.exit(1)

args = sys.argv[1:]

if len(args) == 0:
    help()
    sys.exit(1)

[libPaths , includePaths , exePaths, pythonVersion , minicondaDir,  phony_target] = parseArguments(args)


def check4specialchars(path):
    special_chars = [' ', ')', '(' ]

    outpath = path[:1]
    
    for n, char in enumerate(path[1:]):
        if char in special_chars and path[n] != '\\':
            outpath += "\\" + char
        else:
            outpath += char
    print(path)
    print(outpath)
    return outpath


if 'LD_LIBRARY_PATH' in os.environ:
    os.environ['LD_LIBRARY_PATH'] = check4specialchars(os.environ['LD_LIBRARY_PATH'])
    for path in os.environ['LD_LIBRARY_PATH'].split(':'):
        libPaths.append(path)
else:
    print('Warning: Your system does not have LD_LIBRARY_PATH set.')
    print('PyTom assumes this to be the standard environment variable to store system paths!')
    print('')

if 'INCLUDE_PATH' in os.environ: # most OS does not set this, ignore the warning
    os.environ['INCLUDE_PATH'] = check4specialchars(os.environ['INCLUDE_PATH'])
    for path in os.environ['INCLUDE_PATH'].split(':'):
        includePaths.append(path)
    
if 'PATH' in os.environ:
    os.environ['PATH'] = check4specialchars(os.environ['PATH'])
    for path in os.environ['PATH'].split(':'):
        exePaths.append(path)
else:
    print('Warning: Your system does not have PATH set.')
    print('PyTom assumes this to be the standard environment variable to store system paths!')
    print('')

for path in libPaths:
    newPath = ''
    if '/lib' in path:
        newPath = path.replace('lib','include')
    elif '/lib64' in path:
        newPath = path.replace('lib64','include')   
        
    if os.path.exists(newPath):
        includePaths.append(newPath)

print(includePaths)

# prerequisites
need_exes = ["mpic++"]
need_headers = ["mpi.h", "Python.h", "fftw3.h","ndarrayobject.h"]
need_libs = ["libmpi", "libpython", "libfftw3"]

if platform.system() == "Linux":
    dynamicExtension = '.so'
elif platform.system() == 'Darwin':
    dynamicExtension = '.dylib'
    #dynamicExtension = '.so'
else:
    #win is not supported anyway..
    dynamicExtension = '.dll'

need_libs = [lib+dynamicExtension for lib in need_libs]

#-------------------------
# variables required for compile
mpi_exe = None
lib_mpi = None
include_mpi = None

lib_python = None
include_python = None

lib_fftw = None
include_fftw = None

include_boost = None

include_numpy = None
#-------------------------

# mpi_exe,mpi_exePath = find(["mpicc","mpic++","openmpic++","mpicxx-openmpi-gcc45"],exePaths)
mpi_exe,mpi_exePath = find(["mpicc","mpic++"],exePaths)
if mpi_exePath is None:
    print('MPI compiler not found!')
    exit(1)
        
libraryFile,lib_mpi = find("libmpi" + dynamicExtension,libPaths)
if lib_mpi is None:
    print('MPI library path not found!')
    exit(1)
    
includeFile,include_mpi = find("mpi.h",includePaths)
if include_mpi is None:
    print('MPI include path not found!')
    exit(1)

if 0 and include_python is None:
    try:
        import numpy
        p = numpy.__path__[0].split('lib')[0]
        for a in range(3):

            if [os.path.join(p,d) for d in os.listdir(p) if d == 'include' and os.path.isdir(d)]:
                import glob
                inc = os.path.join(p,'include')
                a = glob.glob(f'{inc}/Python.h') + glob.glob(f'{inc}/*/Python.h')

                if a:
                    inc = os.path.dirname(a[0])
                    includePaths = [inc] + includePaths
                    print(includePaths)
                    break

            p = os.path.dirname(p)

        includeFile, include_python = find("Python.h", includePaths, required=f'python{pythonVersion}')
    except:
        pass

includeFile,include_python = find("Python.h",includePaths, required=f'python{pythonVersion}')
if include_python is None:
    print('Python include path not found!')
    exit(1)

if pythonVersion == '':
    [lib_python, lib_pythonFlag] = adjustLibraryVersions("libpython",['2.6','2.5','2.7', '3.7'],'lpython',dynamicExtension,libPaths)
else:
    libraryFile,lib_python      = find("libpython" + pythonVersion + dynamicExtension,libPaths)
    lib_pythonFlag  = 'lpython' + pythonVersion
    
#if include_fftw is None:
#    
#    includeNew = os.path.dirname(os.popen('locate fftw3.h').read().split()[0])
#    if includeNew:
#        includePaths = [includeNew] + includePaths

includeFile,include_fftw = find("fftw3.h",includePaths)
    
if include_fftw is None:
    print('FFTW include path not found!')
    exit(1)

libraryFile,lib_fftw = find("libfftw3" + dynamicExtension,libPaths)

#if lib_fftw is None:    
#    libPathsNew = os.path.dirname(os.popen('locate libfftw3'+dynamicExtension).read().split()[0])
#    if libPathsNew:
#        libPaths = [libPathsNew] + libPaths

libraryFile, lib_fftw = find("libfftw3"+ dynamicExtension, libPaths)

if lib_fftw is None:
    print('FFTW library path not found!')
    exit(1)

includeFile,include_boost = find("type_traits.hpp",includePaths) 
if include_boost is None:
    print('Boost include path not found!')
    exit(1)
    
if include_boost:
    include_boost = os.path.abspath(os.path.join(include_boost, os.pardir)) + os.sep

includeFile,include_numpy = find("ndarrayobject.h",includePaths)    
if include_numpy is None:
    #includePathsNew = [p for p in os.popen('locate ndarrayobject.h').read().split() if pythonVersion in p]
    #if includePathsNew:
    #    includePaths = [includePathsNew] + includePaths
    #    includeFile,include_numpy = find("ndarrayobject.h",includePaths)    
    if include_numpy is None:
        try:
            import numpy
            path = os.path.join(numpy.__path__[0], 'core/include/numpy')
            includePaths = [path] + includePaths
            includeFile,include_numpy = find("ndarrayobject.h",includePaths)    
        except Exception as e:
            print(e)
            pass
if include_numpy is None:
    print('Numpy include path not found!')
    exit(1)

setenv_line     = ""
setflags_line   = "export PYTOMC_DIR='" + os.getcwd() +"'" 
setflags_line   += " && export LDFLAGS_ADDITIONAL='-D_GLIBCXX_DEBUG'"

if lib_mpi.__class__ == str:
    setenv_line += "export LIBPATH='" + lib_mpi + "'"
    
if mpi_exePath.__class__ == str:
    if len(setenv_line) > 0:
        setenv_line += " && "
        
    setenv_line += "export PATH='" + mpi_exePath + "':" + os.environ['PATH']
    setenv_line += " && export MPICXX='" + mpi_exe + "'"
    
if include_fftw.__class__ == str:
    setflags_line += " && export INCLUDE_FFTW='-I" + include_fftw + "'"

if lib_fftw.__class__ == str:
    setflags_line += " && export LDFLAGS_FFTW='-L" + lib_fftw + " -lfftw3 -lfftw3f'"

if include_python.__class__ == str:
    setflags_line += " && export INCLUDE_PYTHON='-I" + include_python + "'"    

if include_numpy.__class__ == str:
    setflags_line += " && export INCLUDE_NUMPY='-I" + include_numpy + "'"

if lib_python.__class__ == str:
    setflags_line += " && export LDFLAGS_PYTHON='-L" + lib_python + " -" + lib_pythonFlag + "'"

if include_mpi.__class__ == str:
    setflags_line += " && export INCLUDE_MPI='-I" + include_mpi +"'"

if lib_mpi.__class__ == str:
    setflags_line += " && export LDFLAGS_MPI='-L" + lib_mpi + " -lmpi'"

if include_boost.__class__ == str:
    setflags_line += " && export INCLUDE_BOOST=-I'" + include_boost + "'"


if platform.system() == 'Darwin':
    setflags_line += ' && export OSX=" "'

#generate local compile directories
checkCompileDirsExist()

#print found dependencies
print("Check compiling prerequisites:\n")
print("")
print("Python     - lib_python: "+str(lib_python))
print("Python     - include_python: "+str(include_python))
print("FFTW3      - lib_fftw: "+str(lib_fftw))
print("FFTW3      - include_fftw: "+str(include_fftw))
print("OpenMPI    - mpic++: "+str(mpi_exe))
print("OpenMPI    - libmpi: "+str(lib_mpi))
print("OpenMPI    - include_mpi: "+str(include_mpi))
print("Numpy      - include_numpy: "+str(include_numpy))
print('Boost     - include_boost: ' + str(include_boost))

print('')
print('Flags determined: ')
print(setflags_line)
print('')

nosh       = False # you can choose not to compile the SH Alignment library
nompi4py   = False  # you can choose not to compile the mpi4py library
nonfft     = False  # you can choose not to compile the NFFT library
novoltools = True


if phony_target == "swig":
    os.system(setenv_line+" && "+setflags_line+" && "+"make swig")
elif phony_target == "np":
    os.system(setenv_line+" && "+setflags_line+" && "+"make np")
elif phony_target == "smpi":
    os.system(setenv_line+" && "+setflags_line+" && "+"make smpi")
elif phony_target == "libtomc":
    os.system(setenv_line+" && "+setflags_line+" && "+"make libtomc")
elif phony_target == "check":
    os.system(setenv_line+" && "+setflags_line+" && "+"make check")
elif phony_target == "cleanlibtomc":
    os.system("make cleanlibtomc")
elif phony_target == "cleanswig":
    os.system("make cleanswig")
elif phony_target == "clean":
    os.system("make clean")
    exit(1)
elif phony_target == "all":
    os.system(setenv_line+" && "+setflags_line+" && "+"make all")
elif phony_target == "nosh":
    os.system(setenv_line+" && "+setflags_line+" && "+"make all")
    nosh = True
elif phony_target == "nompi4py":
    os.system(setenv_line+" && "+setflags_line+" && "+"make all")
    nompi4py = True
elif phony_target == "nonfft":
    os.system(setenv_line+" && "+setflags_line+" && "+"make all")
    nonfft = True
else:
    raise RuntimeError("The phony target is empty!")


### compile the SH Alignment
sh_ld_library_paths = []
sh_python_paths = None
if nosh is False:
    try:
        print()
        print("############ Start to compile the SH Alignment Library ############")
        print()
        from sh_alignment.compile import compile
        sh_ld_library_paths, sh_python_paths = compile(include_path=includePaths, library_path=libPaths, python_version='python'+pythonVersion)
        sh_ld_library_paths = sh_ld_library_paths.split(':')
        sh_python_paths = sh_python_paths.split(':')
    except Exception as e:
        print(e)
        print("Compilation of SH Alignment failed! Disable this functionality.")


### compile external libraries
# mp4py
if nompi4py is False:
    try:
        print()
        print("############ Start to compile mpi4py Library ############")
        print()
        mpicc_exe, mpicc_exePath = find("mpicc", exePaths)
        mpicc_exe = mpicc_exePath + '/' + mpicc_exe
        os.system("cd ../external/src/mpi4py/ && python"+str(pythonVersion)+" setup.py build --mpicc="+mpicc_exe)
        os.system("cd ../external/src/mpi4py/ && python"+str(pythonVersion)+" setup.py install --prefix=../../")
        print()
        # add into the path
        path1 = os.path.abspath(os.path.join(os.getcwd(), os.pardir)) + '/external/lib/'
        path2 = os.path.abspath(os.path.join(os.getcwd(), os.pardir)) + '/external/lib/python'+str(pythonVersion)+'/site-packages/'
        if sh_python_paths is None:
            sh_python_paths = [path1, path2]
        elif sh_python_paths.__class__ == list:
            sh_python_paths.append(path1)
            sh_python_paths.append(path2)
        else:
            raise Exception()
    except Exception as e:
        print(e)
        print("Compilation of mpi4py failed! Disable this functionality.")

if novoltools is False:
    try:
        print()
        print("############ Start to compile voltools Library ############")
        print()

        os.system("cd ../external/src/voltools/ && python"+str(pythonVersion)+" setup.py install --prefix=../../")
        print()
        
        # add into the path                                                                                                                                                                 
        path1 = os.path.abspath(os.path.join(os.getcwd(), os.pardir)) + '/external/lib/'
        path2 = os.path.abspath(os.path.join(os.getcwd(), os.pardir)) + '/external/lib/python'+str(pythonVersion)+'/site-packages/'
        if sh_python_paths is None:
            sh_python_paths = [path1, path2]
        elif sh_python_paths.__class__ == list:
            if not path1 in sh_python_paths:
                sh_python_paths.append(path1)
            if not path2 in sh_python_paths:
                sh_python_paths.append(path2)
        else:
            raise Exception()

    except Exception as e:
        print(e)
        print("Compilation of voltools failed! Disable this functionality.")


# NFFT
if nonfft is False:
    try:
        print()
        print("############ Start to compile NFFT Library ############")
        print()

        extra = '' 
        if minicondaDir:
            import os
            tt= '../external/src/nfft-3.1.3/include/fftw3.h'
            if os.path.exists(tt):
                os.system(f'unlink {tt}')
                
            os.system(f'ln -s {minicondaDir}/include/fftw3.h {tt}')
            extra = f" export C_INCLUDE_PATH='{minicondaDir}/include'  &&" 
        
        install_path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))+"/external/"
        os.system(f"cd ../external/src/nfft-3.1.3/ && {extra} chmod +x ./configure && ./configure --prefix="+install_path)

    
        os.system("cd ../external/src/nfft-3.1.3/ && make && make install")
        print()
        # compile the swig part
        set_flags = 'export PYTHON_INCLUDE_PATH="%s"' % include_python
        # get the parent directory for numpy header file
        if include_numpy[-5:] == "numpy":
            parent_include_numpy = include_numpy[:-5]
        elif include_numpy[-6:] == "numpy/":
            parent_include_numpy = include_numpy[:-6]
        else:
            raise RuntimeError("Numpy header file arrayobject.h is in the wrong place!")
        set_flags += ' && export NUMPY_INCLUDE_PATH="%s"' % parent_include_numpy
        set_flags += ' && export FFTW_INCLUDE_PATH="%s"' % include_fftw
        set_flags += ' && export FFTW_LIB_PATH="%s"' % lib_fftw
        #set_flags += ' && export NFFT_LIB_PATH="%s"' % ("/usr/local/lib/")
        #set_flags += ' && export NFFT_INCLUDE_PATH="%s"' % ("/Users/gijs/Documents/PostDocUtrecht/nfft/include/")
        set_flags += ' && export NFFT_LIB_PATH="%s"' % (install_path+"lib/")
        set_flags += ' && export NFFT_INCLUDE_PATH="%s"' % (install_path+"src/nfft-3.1.3/include/")
        set_flags += ' && export PYTHON_LIB_PATH="%s"' % lib_python
        set_flags += ' && export PYTHON_VERSION="python%s"' % pythonVersion
        
        command = set_flags+" && cd nufft/ && chmod +x mkswig.sh && ./mkswig.sh"
        print(command)
        os.system(command)
        # add into the path
        py_path = os.path.abspath(os.path.join(os.getcwd(), os.pardir)) + '/pytomc/nufft/'
        if sh_python_paths is None:
            sh_python_paths = [py_path]
        elif sh_python_paths.__class__ == list:
            sh_python_paths.append(py_path)
        else:
            raise Exception()
        ld_path = os.path.abspath(os.path.join(os.getcwd(), os.pardir)) + '/external/lib/'
        sh_ld_library_paths.append(ld_path)
        sh_ld_library_paths.append(py_path)
    except Exception as e:
        print(e)
        print("Compilation of NFFT failed! Disable this functionality.")


if os.path.isfile("../lib/_pytom_fftplan.so") \
    and os.path.isfile("../lib/_pytom_freqweight.so") \
    and os.path.isfile("../lib/_pytom_mpi.so") \
    and os.path.isfile("../lib/_pytom_numpy.so") \
    and os.path.isfile("../lib/_pytom_volume.so") \
    and phony_target and not ('clean' in phony_target):
    print('Generating executables:')
    print('../bin/pytom')
    print('../bin/ipytom')
    print('../bin/pytomGUI')


    genexelibs = list(set([lib_mpi, lib_fftw, lib_python] + sh_ld_library_paths[:1]))
    genexeincl = sh_python_paths

    
    
    generateExecuteables(genexelibs, exePaths, genexeincl, python_version=pythonVersion)

if minicondaDir:
    print('link c-libs', os.getcwd(), minicondaDir)
    os.system(f'cp lib/lib*.so {minicondaDir}/lib')
    os.system(f'cp lib/_pytom*.so {minicondaDir}/lib/python3.8/site-packages/')
    os.system(f'cp lib/*.py {minicondaDir}/lib/python3.8/site-packages/')
