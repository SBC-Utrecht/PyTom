#!/usr/bin/env csh
if ($?LD_LIBRARY_PATH>0) then
setenv LD_LIBRARY_PATH '/usr/local/Cellar/open-mpi/3.1.0/lib/:/usr/local/Cellar/fftw/3.3.7_1/lib/:/usr/local/Cellar/python/3.7.2_1/Frameworks/Python.framework/Versions/3.7/lib/:/usr/local/Cellar/fftw/3.3.7_1/lib/:/Users/gijs/Documents/PyTomPrivate/pytomc/sh_alignment/SpharmonicKit27/:/Users/gijs/Documents/PyTomPrivate/pytomc/sh_alignment/frm/swig/:/Users/gijs/Documents/PyTomPrivate/external/lib/:/Users/gijs/Documents/PyTomPrivate/pytomc/nufft/:/Users/gijs/Documents/PyTomPrivate/pytomc/libs/libtomc/libs':$LD_LIBRARY_PATH
else
setenv LD_LIBRARY_PATH '/usr/local/Cellar/open-mpi/3.1.0/lib/:/usr/local/Cellar/fftw/3.3.7_1/lib/:/usr/local/Cellar/python/3.7.2_1/Frameworks/Python.framework/Versions/3.7/lib/:/usr/local/Cellar/fftw/3.3.7_1/lib/:/Users/gijs/Documents/PyTomPrivate/pytomc/sh_alignment/SpharmonicKit27/:/Users/gijs/Documents/PyTomPrivate/pytomc/sh_alignment/frm/swig/:/Users/gijs/Documents/PyTomPrivate/external/lib/:/Users/gijs/Documents/PyTomPrivate/pytomc/nufft/:/Users/gijs/Documents/PyTomPrivate/pytomc/libs/libtomc/libs'
endif

if ($?PATH>0) then
setenv PATH '/usr/local/Cellar/gcc/7.3.0_1/bin/:/Users/gijs/Documents/PyTomPrivate/bin:/usr/local/opt/llvm/bin:/usr/local/bin:/usr/bin/:/usr/local/opt/qt/bin:/Users/gijs/Documents/PyTomPrivate/bin:/usr/local/opt/llvm/bin:/usr/local/bin:/usr/bin/:/usr/local/opt/qt/bin:/Applications/IMOD/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/Library/TeX/texbin:/opt/X11/bin':$PATH
else
setenv PATH '/usr/local/Cellar/gcc/7.3.0_1/bin/:/Users/gijs/Documents/PyTomPrivate/bin:/usr/local/opt/llvm/bin:/usr/local/bin:/usr/bin/:/usr/local/opt/qt/bin:/Users/gijs/Documents/PyTomPrivate/bin:/usr/local/opt/llvm/bin:/usr/local/bin:/usr/bin/:/usr/local/opt/qt/bin:/Applications/IMOD/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/Library/TeX/texbin:/opt/X11/bin'
endif

if ($?PYTHONPATH>0) then
setenv PYTHONPATH '/Users/gijs/Documents:/Users/gijs/Documents/PyTomPrivate/pytomc:/Users/gijs/Documents/PyTomPrivate/pytomc/sh_alignment/frm/swig/:/Users/gijs/Documents/PyTomPrivate/external/lib/:/Users/gijs/Documents/PyTomPrivate/external/lib/python3.7/site-packages/:/Users/gijs/Documents/PyTomPrivate/pytomc/nufft/:/Users/gijs/Documents:/Users/gijs/Documents/PyTomPrivate/pytomc/swigModules':$PYTHONPATH
else
setenv PYTHONPATH '/Users/gijs/Documents:/Users/gijs/Documents/PyTomPrivate/pytomc:/Users/gijs/Documents/PyTomPrivate/pytomc/sh_alignment/frm/swig/:/Users/gijs/Documents/PyTomPrivate/external/lib/:/Users/gijs/Documents/PyTomPrivate/external/lib/python3.7/site-packages/:/Users/gijs/Documents/PyTomPrivate/pytomc/nufft/:/Users/gijs/Documents:/Users/gijs/Documents/PyTomPrivate/pytomc/swigModules'
endif

