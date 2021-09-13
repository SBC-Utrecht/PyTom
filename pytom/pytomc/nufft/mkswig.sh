#!/bin/sh

echo $PYTHON_VERSION
echo $NUMPY_INCLUDE_PATH
echo $FFTW_INCLUDE_PATH
echo $NFFT_INCLUDE_PATH
echo $FFTW_LIB_PATH
echo $NFFT_LIB_PATH
echo $PYTHON_LIB_PATH


echo `which swig`

swig -python nufft.i
gcc -fPIC -O3 -c nufft.c nufft_wrap.c -I$PYTHON_INCLUDE_PATH -I$NUMPY_INCLUDE_PATH -I$FFTW_INCLUDE_PATH -I$NFFT_INCLUDE_PATH
cc -shared nufft.o nufft_wrap.o -L$FFTW_LIB_PATH -lfftw3 -lm -L$NFFT_LIB_PATH -lnfft3 -L$PYTHON_LIB_PATH -l$PYTHON_VERSION -o ../../lib/_swig_nufft.so
