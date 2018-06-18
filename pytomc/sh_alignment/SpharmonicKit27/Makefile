############################################################################
############################################################################
##
##		Spherical Harmonic Transform Kit 2.5
##
## Sean Moore, Dennis Healy, Dan Rockmore, Peter Kostelec
## smoore@bbn.com, {healy,rockmore,geelong}@cs.dartmouth.edu
##
## Contact: Peter Kostelec
##          geelong@cs.dartmouth.edu
##
##
## Copyright 1997, 1998 Sean Moore, Dennis Healy,
##                      Dan Rockmore, Peter Kostelec
##
##
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
##
##
## Commercial use is absolutely prohibited.
##
## See the accompanying LICENSE file for details.
##
############################################################################
############################################################################

CC = gcc
F77 = f77

# define WALLCLOCK on the CFLAGS line if want to time
# walltime and not cpu time (cpu time is default); also
# define optimization flags suitable for your platform

# define FFTPACK if want to use fftpack or your own fft
# and dct routines; *not* defining FFTPACK means the
# code will use slower fft and dct routines provided
# in this distribution
# CFLAGS = -O3 -fp_reorder -DWALLCLOCK
# CFLAGS = -O3 -fp_reorder -DFFTPACK

#CFLAGS = -O2
CFLAGS = -O3 -g -fPIC
LDFLAGS = -lm

# If FFTPACK *is* defined, uncomment the FFTFLAGS line below
# and set the path to where the library lives, so the compiler
# knows where to link to

# FFTFLAGS = -L/net/puddleby/geelong/prcstmp/modified_fftpack -lfftpack
FFTFLAGS =

# for naive
NAIVESRC = naive_synthesis.c weights.c csecond.c

NAIVEOBJ = naive_synthesis.o weights.o csecond.o


# for semi-naive
SEMISRC = OURperms.c OURmods.c newFCT.c weights.c \
	cospmls.c seminaive.c oddweights.c \
	csecond.c primitive.c

SEMIOBJ = OURperms.o OURmods.o newFCT.o weights.o \
	cospmls.o seminaive.o oddweights.o \
	csecond.o primitive.o

# for FLT_classic
FLTCLASSICSRC = OURperms.c OURmods.c weights.c \
	  primitive.c precomp_flt_classic.c \
	  newFCT.c flt_classic.c \
	  csecond.c

FLTCLASSICOBJ = OURperms.o OURmods.o weights.o \
	  primitive.o precomp_flt_classic.o \
	  newFCT.o flt_classic.o \
	  csecond.o

# for FLT_dhmid
FLTDHMIDSRC = OURperms.c OURmods.c weights.c \
	  primitive.c precomp_flt_dhmid.c \
	  newFCT.c flt_dhmid.c \
	  csecond.c

FLTDHMIDOBJ = OURperms.o OURmods.o weights.o \
	  primitive.o precomp_flt_dhmid.o \
	  newFCT.o flt_dhmid.o \
	  csecond.o

# for use with FLT_hybrid
FLTHYBRIDSRC = $(SEMISRC) \
	precomp_flt_hybrid.c flt_hybrid.c

FLTHYBRIDOBJ = $(SEMIOBJ) \
	precomp_flt_hybrid.o flt_hybrid.o 


# seminaive spherical transform and convolution
FSTSEMISRC = $(SEMISRC) naive_synthesis.c \
	primitive_FST.c permroots.c indextables.c FFTcode.c \
	fft_grids.c

FSTSEMIOBJ = $(SEMIOBJ) naive_synthesis.o \
	primitive_FST.o permroots.o indextables.o FFTcode.o \
	fft_grids.o

CONVSEMISRC = $(FSTSEMISRC) MathFace.c

CONVSEMIOBJ = $(FSTSEMIOBJ) MathFace.o


# hybrid/seminaive spherical transform and convolution
FSTHYBRIDSRC =  $(FLTHYBRIDSRC) \
	  primitive_FST.c permroots.c indextables.c FFTcode.c \
	  fft_grids.c

FSTHYBRIDOBJ = $(FLTHYBRIDOBJ) \
	  primitive_FST.o permroots.o indextables.o FFTcode.o \
	  fft_grids.o

FSTHYBRIDXSRC =  $(FLTHYBRIDSRC) \
	  primitive_FST.c permroots.c indextables.c FFTcode.c \
	  fft_grids.c

FSTHYBRIDXOBJ = $(FLTHYBRIDOBJ) \
	  primitive_FST.o permroots.o indextables.o FFTcode.o \
	  fft_grids.o

CONVHYBRIDSRC = $(FSTHYBRIDSRC) \
	  MathFace.c

CONVHYBRIDOBJ = $(FSTHYBRIDOBJ) \
	  MathFace.o


###################################################################
##################################################################
######
######              things that can be made
######
##################################################################
##################################################################

# first some shortcuts

all:
	make \
	legendre \
	sphere \
	shared


legendre:
	make \
	test_naive \
	test_stability_naive \
	test_semi \
	test_flt_classic \
	test_flt_dhmid \
	test_flt_hybrid

sphere:
	make \
	FST_precomp2disk \
	test_FST_semi_memo \
	test_FST_semi_disk \
	test_FST_semi_fly \
	CONV_SEMI_DISK \
	CONV_SEMI_MEMO \
	CONV_SEMI_FLY \
	test_FST_hybrid_disk \
	test_FST_hybrid_memo \
	test_FST_hybrid_memoX \
	CONV_HYBRID_DISK \
	CONV_HYBRID_MEMO

clean: 
	rm *.o

shared:
	cc -shared OURperms.o OURmods.o newFCT.o weights.o cospmls.o seminaive.o oddweights.o csecond.o primitive.o FST_semi_memo.o naive_synthesis.o primitive_FST.o permroots.o indextables.o FFTcode.o fft_grids.o -o libsphkit.so

# now the make definitions for the individual executables

test_naive: $(NAIVEOBJ) test_naive.o
	$(CC) $(CFLAGS) $(NAIVEOBJ) test_naive.o \
	$(LDFLAGS) -o test_naive

test_stability_naive: $(NAIVEOBJ) test_stability_naive.o
	$(CC) $(CFLAGS) $(NAIVEOBJ) test_stability_naive.o \
	$(LDFLAGS) -o test_stability_naive

test_semi: $(SEMIOBJ) test_semi.o
	$(CC) $(CFLAGS) $(SEMIOBJ) test_semi.o \
	$(FFTFLAGS) $(LDFLAGS) -o test_semi

test_semi_roundtrip: $(SEMIOBJ) naive_synthesis.o test_semi_roundtrip.o
	$(CC) $(CFLAGS) $(SEMIOBJ) naive_synthesis.o test_semi_roundtrip.o \
	$(FFTFLAGS) $(LDFLAGS) -o test_semi_roundtrip

test_flt_classic: $(FLTCLASSICOBJ) test_flt_classic.o
	$(CC) $(CFLAGS) $(FLTCLASSICOBJ) test_flt_classic.o \
	$(FFTFLAGS) $(LDFLAGS) -o test_flt_classic

test_flt_dhmid: $(FLTDHMIDOBJ) test_flt_dhmid.o
	$(CC) $(CFLAGS) $(FLTDHMIDOBJ) test_flt_dhmid.o \
	$(FFTFLAGS) $(LDFLAGS) -o test_flt_dhmid

test_flt_hybrid: $(FLTHYBRIDOBJ) test_flt_hybrid.o
	$(CC) $(CFLAGS) $(FLTHYBRIDOBJ) test_flt_hybrid.o \
	$(FFTFLAGS) $(LDFLAGS) -o test_flt_hybrid

FST_precomp2disk: $(FSTHYBRIDOBJ) FST_precomp2disk.o
	$(CC) $(CFLAGS)  $(FSTHYBRIDOBJ) \
	FST_precomp2disk.o $(FFTFLAGS) $(LDFLAGS) -o FST_precomp2disk

test_FST_semi_memo: $(FSTSEMIOBJ) FST_semi_memo.o test_FST_semi_memo.o
	$(CC) $(CFLAGS) $(FSTSEMIOBJ) FST_semi_memo.o test_FST_semi_memo.o \
	$(FFTFLAGS) $(LDFLAGS) -o test_FST_semi_memo

test_FST_semi_disk: $(FSTSEMIOBJ) FST_semi_disk.o test_FST_semi_disk.o
	$(CC) $(CFLAGS) $(FSTSEMIOBJ) FST_semi_disk.o test_FST_semi_disk.o \
	$(FFTFLAGS) $(LDFLAGS) -o test_FST_semi_disk

test_FST_semi_fly: $(FSTSEMIOBJ) FST_semi_fly.o test_FST_semi_fly.o
	$(CC) $(CFLAGS) $(FSTSEMIOBJ) FST_semi_fly.o test_FST_semi_fly.o \
	$(FFTFLAGS) $(LDFLAGS) -o test_FST_semi_fly

CONV_SEMI_DISK: $(CONVSEMIOBJ) FST_semi_disk.o CONV_SEMI_DISK.o
	$(CC) $(CFLAGS) $(CONVSEMIOBJ) FST_semi_disk.o CONV_SEMI_DISK.o \
	$(FFTFLAGS) $(LDFLAGS) -o CONV_SEMI_DISK

CONV_SEMI_MEMO: $(CONVSEMIOBJ) FST_semi_memo.o CONV_SEMI_MEMO.o
	$(CC) $(CFLAGS) $(CONVSEMIOBJ) FST_semi_memo.o CONV_SEMI_MEMO.o \
	$(FFTFLAGS) $(LDFLAGS) -o CONV_SEMI_MEMO

CONV_SEMI_FLY: $(CONVSEMIOBJ) FST_semi_fly.o CONV_SEMI_FLY.o
	$(CC) $(CFLAGS) $(CONVSEMIOBJ) FST_semi_fly.o CONV_SEMI_FLY.o \
	$(FFTFLAGS) $(LDFLAGS) -o CONV_SEMI_FLY

test_FST_hybrid_disk: $(FSTHYBRIDOBJ) FST_hybrid_disk.o test_FST_hybrid_disk.o
	$(CC) $(CFLAGS) \
	$(FSTHYBRIDOBJ) FST_hybrid_disk.o test_FST_hybrid_disk.o \
	$(FFTFLAGS) $(LDFLAGS) -o test_FST_hybrid_disk

test_FST_hybrid_memo: $(FSTHYBRIDOBJ) FST_hybrid_memo.o test_FST_hybrid_memo.o
	$(CC) $(CFLAGS) \
	$(FSTHYBRIDOBJ) FST_hybrid_memo.o test_FST_hybrid_memo.o \
	$(FFTFLAGS) $(LDFLAGS) -o test_FST_hybrid_memo

test_FST_hybrid_memoX: $(FSTHYBRIDXOBJ) FST_hybrid_memoX.o test_FST_hybrid_memoX.o
	$(CC) $(CFLAGS) \
	$(FSTHYBRIDXOBJ) FST_hybrid_memoX.o test_FST_hybrid_memoX.o \
	$(FFTFLAGS) $(LDFLAGS) -o test_FST_hybrid_memoX

CONV_HYBRID_DISK: $(CONVHYBRIDOBJ) FST_hybrid_disk.o CONV_HYBRID_DISK.o
	$(CC) $(CFLAGS) \
	$(CONVHYBRIDOBJ) FST_hybrid_disk.o CONV_HYBRID_DISK.o \
	$(FFTFLAGS) $(LDFLAGS) -o CONV_HYBRID_DISK

CONV_HYBRID_MEMO: $(CONVHYBRIDOBJ) FST_hybrid_memo.o CONV_HYBRID_MEMO.o
	$(CC) $(CFLAGS) \
	$(CONVHYBRIDOBJ) FST_hybrid_memo.o CONV_HYBRID_MEMO.o \
	$(FFTFLAGS) $(LDFLAGS) -o CONV_HYBRID_MEMO



# and now for LOTS OF dependencies ...

# DO NOT DELETE THIS LINE -- make depend depends on it.

CONV_HYBRID_DISK.o: CONV_HYBRID_DISK.c FST_hybrid_disk.h MathFace.h \
 cospmls.h precomp_flt_hybrid.h

CONV_HYBRID_MEMO.o: CONV_HYBRID_MEMO.c FST_hybrid_memo.h MathFace.h \
 cospmls.h precomp_flt_hybrid.h

CONV_SEMI_DISK.o: CONV_SEMI_DISK.c FST_semi_disk.h MathFace.h \
 cospmls.h

CONV_SEMI_FLY.o: CONV_SEMI_FLY.c FST_semi_fly.h MathFace.h cospmls.h

CONV_SEMI_MEMO.o: CONV_SEMI_MEMO.c FST_semi_memo.h MathFace.h \
 cospmls.h

FFTcode.o: FFTcode.c indextables.h permroots.h

FST_hybrid_disk.o: FST_hybrid_disk.c FFTcode.h config.h cospmls.h \
 flt_hybrid.h primitive.h primitive_FST.h seminaive.h

FST_hybrid_memo.o: FST_hybrid_memo.c FFTcode.h cospmls.h flt_hybrid.h \
 primitive.h primitive_FST.h seminaive.h

FST_hybrid_memoX.o: FST_hybrid_memoX.c FFTcode.h cospmls.h \
 flt_hybrid.h primitive.h primitive_FST.h seminaive.h

FST_precomp2disk.o: FST_precomp2disk.c config.h cospmls.h \
 precomp_flt_hybrid.h

FST_semi_disk.o: FST_semi_disk.c FFTcode.h config.h cospmls.h \
 naive_synthesis.h primitive.h primitive_FST.h seminaive.h

FST_semi_fly.o: FST_semi_fly.c FFTcode.h cospmls.h naive_synthesis.h \
 primitive.h primitive_FST.h seminaive.h

FST_semi_memo.o: FST_semi_memo.c FFTcode.h cospmls.h naive_synthesis.h \
 primitive.h primitive_FST.h seminaive.h

MathFace.o: MathFace.c

OURmods.o: OURmods.c

OURperms.o: OURperms.c

cospmls.o: cospmls.c newFCT.h primitive.h

csecond.o: csecond.c

flt_classic.o: flt_classic.c csecond.h newFCT.h precomp_flt_classic.h \
 primitive.h weights.h

flt_dhmid.o: flt_dhmid.c csecond.h newFCT.h precomp_flt_dhmid.h \
 primitive.h weights.h

flt_hybrid.o: flt_hybrid.c csecond.h cospmls.h newFCT.h oddweights.h \
 seminaive.h weights.h

indextables.o: indextables.c

naive_synthesis.o: naive_synthesis.c csecond.h weights.h

newFCT.o: newFCT.c OURperms.h OURmods.h newFCT.h

oddweights.o: oddweights.c

permroots.o: permroots.c

precomp_flt_classic.o: precomp_flt_classic.c newFCT.h primitive.h \
 weights.h

precomp_flt_dhmid.o: precomp_flt_dhmid.c newFCT.h primitive.h \
 weights.h

precomp_flt_hybrid.o: precomp_flt_hybrid.c cospmls.h newFCT.h \
 primitive.h

primitive.o: primitive.c

primitive_FST.o: primitive_FST.c

semi.o: semi.c cospmls.h seminaive.h

seminaive.o: seminaive.c cospmls.h csecond.h newFCT.h oddweights.h \
 weights.h

test_FST_hybrid_disk.o: test_FST_hybrid_disk.c FST_hybrid_disk.h \
 cospmls.h csecond.h precomp_flt_hybrid.h primitive_FST.h

test_FST_hybrid_memo.o: test_FST_hybrid_memo.c FST_hybrid_memo.h \
 cospmls.h csecond.h precomp_flt_hybrid.h primitive_FST.h

test_FST_hybrid_memoX.o: test_FST_hybrid_memoX.c FST_hybrid_memo.h \
 cospmls.h csecond.h precomp_flt_hybrid.h primitive_FST.h

test_FST_semi_disk.o: test_FST_semi_disk.c FST_semi_disk.h cospmls.h \
 csecond.h primitive_FST.h

test_FST_semi_fly.o: test_FST_semi_fly.c FST_semi_fly.h cospmls.h \
 csecond.h primitive_FST.h seminaive.h

test_FST_semi_memo.o: test_FST_semi_memo.c FST_semi_memo.h csecond.h \
 cospmls.h primitive_FST.h seminaive.h

test_flt_classic.o: test_flt_classic.c flt_classic.h primitive.h

test_flt_dhmid.o: test_flt_dhmid.c flt_dhmid.h primitive.h

test_flt_hybrid.o: test_flt_hybrid.c cospmls.h flt_hybrid.h \
 precomp_flt_hybrid.h primitive.h

test_naive.o: test_naive.c naive_synthesis.h

test_semi.o: test_semi.c cospmls.h primitive.h seminaive.h

test_semi_roundtrip.o: test_semi_roundtrip.c naive_synthesis.h \
 cospmls.h seminaive.h

test_stability_naive.o: test_stability_naive.c naive_synthesis.h

weights.o: weights.c
