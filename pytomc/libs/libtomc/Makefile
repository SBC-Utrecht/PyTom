

include ./Makefile.params



OBJECTS=    $(LIBDIR)/tom__FreqWeight.o \
            $(LIBDIR)/tom__volume.o \
            $(LIBDIR)/tom__volume_fcn.o \
            $(LIBDIR)/tom__fftw__fftw_plan.o \
            $(LIBDIR)/tom__io__io.o \
            $(LIBDIR)/tom__io__io_cpp.o \
            $(LIBDIR)/tom__tools__Log.o \
            $(LIBDIR)/tom__tools__ParameterSet.o \
            $(LIBDIR)/tom__tools__snippets.o \
            $(LIBDIR)/tom__transf__transform_cpp.o \

LIBTOMC_SHARED=$(LIBDIR)/libtomc.so

.PHONY: default

default: shared 


.PHONY: all
all: shared_all
.PHONY: shared 
shared: shared_tomc
.PHONY: shared_tomc
shared_tomc: $(LIBTOMC_SHARED)



$(LIBTOMC_SHARED): $(OBJECTS)
	$(CXX) $(SHAREFLAGS) $(WLFLAGS) -o $@ $^ $(LDFLAGS_)




# OBJECTS
$(LIBDIR)/tom__volume.o: src/tom/volume.cpp $(HEADERS_)
	$(CXX) $(CFLAGS_) -c -o $@ $< $(INCLUDE_)
$(LIBDIR)/tom__volume_fcn.o: src/tom/volume_fcn.cpp	$(HEADERS_)
	$(CXX) $(CFLAGS_) -c -o $@ $< $(INCLUDE_)
$(LIBDIR)/tom__FreqWeight.o: src/tom/FreqWeight.cpp $(HEADERS_)
	$(CXX) $(CFLAGS_) -c -o $@ $< $(INCLUDE_)
$(LIBDIR)/tom__fftw__fftw_plan.o: src/tom/fftw/fftw_plan.cpp $(HEADERS_)
	$(CXX) $(CFLAGS_) -c -o $@ $< $(INCLUDE_)
$(LIBDIR)/tom__io__io_cpp.o: src/tom/io/io.cpp $(HEADERS_)
	$(CXX) $(CFLAGS_) -c -o $@ $< $(INCLUDE_)
$(LIBDIR)/tom__io__io.o: src/tom/io/io.c $(HEADERS_)
	$(CXX) $(CFLAGS_) -c -o $@ $< $(INCLUDE_)
$(LIBDIR)/tom__transf__transform_cpp.o: src/tom/transf/transform.cpp $(HEADERS_)
	$(CXX) $(CFLAGS_) -c -o $@ $< $(INCLUDE_)
$(LIBDIR)/tom__tools__Log.o: src/tom/tools/Log.cpp $(HEADERS_)
	$(CXX) $(CFLAGS_) -c -o $@ $< $(INCLUDE_)
$(LIBDIR)/tom__tools__ParameterSet.o: src/tom/tools/ParameterSet.cpp $(HEADERS_)
	$(CXX) $(CFLAGS_) -c -o $@ $< $(INCLUDE_)
$(LIBDIR)/tom__tools__snippets.o: src/tom/tools/snippets.cpp $(HEADERS_)
	$(CXX) $(CFLAGS_) -c -o $@ $< $(INCLUDE_)





.PHONY: doxygen_rebuild
doxygen_rebuild: doxygen_clean doxygen
.PHONY: doxygen_clean
doxygen_clean:
	rm -rf ./doc/html ./doc/latex
.PHONY: 
doxygen: doxygen
	doxygen src/Doxyfile

.PHONY: rebuild
rebuild: clean default


.PHONY: clean
clean:
	


