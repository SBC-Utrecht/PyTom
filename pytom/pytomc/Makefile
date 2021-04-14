CC = g++

LIBTOMC_BASEDIR=./libs/libtomc

LDFLAGS_LIBTOMC=-L$(PYTOMC_DIR)/libs/libtomc/libs/ -ltomc
LIBTOMC_STATIC=shared
LIBS = -L$(PYTOMC_DIR)/libs/libtomc/libs/
LDFLAGS_TOMC = -lstdc++ -shared -ltomc -fpermissive

     

libtomc:
	make -C $(LIBTOMC_BASEDIR) shared


ifdef LDFLAGS_FFTW
LDFLAGS := $(LDFLAGS_TOMC) $(LDFLAGS_FFTW)
else
LDFLAGS := $(LDFLAGS_TOMC)
endif

ifdef LDFLAGS_PYTHON
LDFLAGS += $(LDFLAGS_PYTHON)
endif

ifdef LDFLAGS_MPI
LDFLAGS += $(LDFLAGS_MPI)
endif


ifndef GCCSWIGFLAGS
GCCSWIGFLAGS :=-fPIC -c -g3 -Wall -Wno-long-long -pedantic -ansi --std=c++98
endif 

PYMODULESPATH=./swigModules/


SWIGCODE=./swigCode/
SWIGCPP=./swigCpp/
SWIGLIBS=./swigLibs/
SWIGTEST=./test/
SWIGFLAGS=-Wall -python -c++ -outdir $(PYMODULESPATH)


ifdef LDFLAGS_ADDITIONAL
GCCSWIGFLAGS += $(LDFLAGS_ADDITIONAL)
endif

ifdef INCLUDE_FFTW
INCLUDE := -I$(LIBTOMC_BASEDIR)/include -I./include $(INCLUDE_FFTW) 
else
INCLUDE := -I$(LIBTOMC_BASEDIR)/include -I./include
endif

ifdef INCLUDE_PYTHON
INCLUDE += $(INCLUDE_PYTHON)
endif

ifdef INCLUDE_BOOST
INCLUDE += $(INCLUDE_BOOST)
endif

ifdef INCLUDE_PYTHON
INCLUDE += $(INCLUDE_PYTHON)
endif

ifdef INCLUDE_NUMPY
INCLUDE += $(INCLUDE_NUMPY)
endif


HEADERS=$(LIBTOMC_BASEDIR)/include/tom/*.h* $(LIBTOMC_BASEDIR)/include/tom/*/*.h* ./include/*h*

SRCS=./src/*c*


SWIGCLASSES=$(SWIGLIBS)swigVolume.o

$(SWIGLIBS)swigVolume.o: src/swigVolume.cpp
	$(CC) $(GCCSWIGFLAGS) $(INCLUDE) $< -o $@

SWIGOBJECTS= 	libtomc \
				$(SWIGLIBS)swigVolume.o \
				$(PYMODULESPATH)_pytom_volume.so \
				$(PYMODULESPATH)_pytom_fftplan.so \
				$(PYMODULESPATH)_pytom_freqweight.so
					
SWIGINCLUDE= $(INCLUDE) -I./src

NUMPYOBJECTS= libtomc \
			  $(SWIGLIBS)swigVolume.o \
			  $(PYMODULESPATH)_pytom_volume.so \
			  $(PYMODULESPATH)_pytom_numpy.so	


MPIOBJECTS= $(SWIGLIBS)swigMPI.o \
			$(PYMODULESPATH)_pytom_mpi.so


$(SWIGLIBS)swigMPI.o: src/swigMPI.cpp
	$(MPICXX) $(GCCSWIGFLAGS) $(INCLUDE) $(INCLUDE_MPI) $< -o $@	
$(SWIGCPP)pytom_mpi_wrap.cpp: $(SWIGCODE)pytom_mpi.i $(HEADERS) ./src/swigMPI.cpp
	swig $(SWIGFLAGS) $(SWIGINCLUDE) $(INCLUDE_MPI) -o $@ $<
$(SWIGLIBS)pytom_mpi_wrap.o: $(SWIGCPP)pytom_mpi_wrap.cpp
	$(MPICXX) $(GCCSWIGFLAGS) $(INCLUDE) $< -o $@
$(PYMODULESPATH)_pytom_mpi.so: $(SWIGLIBS)pytom_mpi_wrap.o
	$(MPICXX)  $< -o $@  $(LDFLAGS) $(LDFLAGS_MPI) $(LIBS) $(INCLUDE_MPI)

$(SWIGCPP)pytom_volume_wrap.cpp: $(SWIGCODE)pytom_volume.i $(HEADERS) ./src/swigVol*c*
	swig $(SWIGFLAGS) $(SWIGINCLUDE) -o $@ $< 	
$(SWIGLIBS)pytom_volume_wrap.o: $(SWIGCPP)pytom_volume_wrap.cpp
	$(CC) $(GCCSWIGFLAGS) $(INCLUDE) $< -o $@
$(PYMODULESPATH)_pytom_volume.so: $(SWIGLIBS)pytom_volume_wrap.o
	$(CC)  $< -o $@  $(LDFLAGS) $(LIBS)

$(SWIGCPP)pytom_fftplan_wrap.cpp: $(SWIGCODE)pytom_fftplan.i $(HEADERS) ./src/swigFft*c*
	swig $(SWIGFLAGS) $(SWIGINCLUDE) -o $@ $< 	
$(SWIGLIBS)pytom_fftplan_wrap.o: $(SWIGCPP)pytom_fftplan_wrap.cpp
	$(CC) $(GCCSWIGFLAGS) $(INCLUDE) $< -o $@
$(PYMODULESPATH)_pytom_fftplan.so: $(SWIGLIBS)pytom_fftplan_wrap.o
	$(CC)  $< -o $@  $(LDFLAGS) $(LIBS)
	
$(SWIGCPP)pytom_freqweight_wrap.cpp: $(SWIGCODE)pytom_freqweight.i $(HEADERS) ./src/swigFreqWeight*c*
	swig $(SWIGFLAGS) $(SWIGINCLUDE) -o $@ $< 	
$(SWIGLIBS)pytom_freqweight_wrap.o: $(SWIGCPP)pytom_freqweight_wrap.cpp
	$(CC) $(GCCSWIGFLAGS) $(INCLUDE) $< -o $@
$(PYMODULESPATH)_pytom_freqweight.so: $(SWIGLIBS)pytom_freqweight_wrap.o
	$(CC)  $< -o $@  $(LDFLAGS) $(LIBS)
	
$(SWIGCPP)pytom_numpy_wrap.cpp: $(SWIGCODE)pytom_numpy.i $(HEADERS) ./src/swigTomNumpy*c* ./src/swigVol*c*
	swig $(SWIGFLAGS) $(SWIGINCLUDE) -o $@ $< 	
$(SWIGLIBS)pytom_numpy_wrap.o: $(SWIGCPP)pytom_numpy_wrap.cpp
	$(CC) -fPIC -c -g -Wall -Wno-long-long  -ansi --std=c++98 $(INCLUDE) $< -o $@ #-I$(NUMPYINCLUDE) 
$(PYMODULESPATH)_pytom_numpy.so: $(SWIGLIBS)pytom_numpy_wrap.o
	$(CC)  $< -o $@  $(LDFLAGS) $(LIBS)
	
ct: $(HEADERS) $(SRCS)
	$(CC) $(LDFLAGS_LIBTOMC) $(INCLUDE) src/swigVolume.cpp src/swigVolumeFnc.cpp ct.cpp -o c.x -g
	
s : swig
swig: $(SWIGOBJECTS)
np : $(NUMPYOBJECTS)
smpi: $(MPIOBJECTS)
all : swig smpi np

default: libtomc all
		
sc: swigclean
	
swigclean:
	rm -rf $(SWIGCPP)*
	rm -rf $(SWIGLIBS)*
	rm -rf $(PYMODULESPATH)*.py*
	rm -rf $(PYMODULESPATH)*so
	rm -rf $(SWIGTEST)*pyc
	
clean:
	rm -rf ./libs/libtomc/libs/*
	rm -rf $(SWIGCPP)*
	rm -rf $(SWIGLIBS)*
	rm -rf $(PYMODULESPATH)*.py*
	rm -rf $(PYMODULESPATH)*so
	rm -rf $(SWIGTEST)*pyc
		
check:
	ldd -r -d $(PYMODULESPATH)*so
	
svnup:
	svn update
	make all
	echo 'import PyTom_Test;PyTom_Test.run()' | python
	
	

unexport LDFLAGS
unexport LDFLAGS_LIBTOMC
