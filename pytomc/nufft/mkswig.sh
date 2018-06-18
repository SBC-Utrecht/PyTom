#!/bin/sh

swig -python nufft.i
gcc -fPIC -O3 -c nufft.c nufft_wrap.c -I$PYTHON_INCLUDE_PATH -I$NUMPY_INCLUDE_PATH -I$FFTW_INCLUDE_PATH -I$NFFT_INCLUDE_PATH
cc -shared nufft.o nufft_wrap.o -L$FFTW_LIB_PATH -lfftw3 -lm -L$NFFT_LIB_PATH -lnfft3 -o _swig_nufft.so
