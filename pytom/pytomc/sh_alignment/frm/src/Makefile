CC = gcc

CFLAGS = -O3 -fPIC

# CFLAGS = -O3 -fPIC -g -I$(FFTW2_INCLUDE_PATH) 
# LDFLAGS =  -L$(FFTW2_LIB_PATH) -lfftw -lm
# LDSRFLAGS = -lsrfftw -lsfftw -lm
# SPHKITFLAG= -L../../SpharmonicKit27 -lsphkit # sphkit lib

SITUSLIB = lib_vio.o lib_pio.o lib_std.o lib_eul.o \
           lib_pwk.o lib_vec.o lib_vwk.o lib_tim.o

PROGS = frmr frmv frm6 

default: frm_rad frm_3d frm_6d 

lib: $(SITUSLIB)

# frm_rad:$(SITUSLIB) frm_rad.o
# 	$(CC) $(CFLAGS) $(SITUSLIB) frm_rad.o $(SPHKITFLAG) \
#         $(LDFLAGS) -o frmr

# frm_3d: $(SITUSLIB) frm_3d.o
# 	$(CC) $(CFLAGS) $(SITUSLIB) frm_3d.o $(SPHKITFLAG) \
#         $(LDFLAGS) -o frmv

# frm_6d: $(SITUSLIB) frm_6d.o
# 	$(CC) $(CFLAGS) $(SITUSLIB) frm_6d.o $(SPHKITFLAG) \
#         $(LDSRFLAGS) -o frm6

install: 
	mv -f $(PROGS) ../bin

clean: 
	rm -f *.o $(PROGS) *~ core

