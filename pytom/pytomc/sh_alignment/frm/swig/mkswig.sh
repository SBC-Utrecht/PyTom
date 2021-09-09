#!/bin/sh

swig -python frm.i
gcc -O3 -fPIC -g -c frm.c frm_wrap.c -I../src -I$PYTHON_INCLUDE_PATH -I$NUMPY_INCLUDE_PATH -I$FFTW_INCLUDE_PATH -I../../SpharmonicKit27
cc -shared -L$PYTHON_LIB_PATH -l$PYTHON_VERSION frm.o frm_wrap.o -L$FFTW_LIB_PATH -lfftw3 -lm ../src/lib_vio.o ../src/lib_pio.o ../src/lib_std.o ../src/lib_eul.o ../src/lib_pwk.o ../src/lib_vec.o ../src/lib_vwk.o ../src/lib_tim.o -L../../SpharmonicKit27 -lsphkit -o _swig_frm.so
