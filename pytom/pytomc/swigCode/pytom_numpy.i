%module pytom_numpy

%include "exception.i"

%exception {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

%{
#include "swigTomNumpy.hpp"
#include "../src/swigTomNumpy.cpp"
%}

%include "swigTomNumpy.hpp"
%include "swigTomNumpy.cpp"




%template(vol2npy) swigTom::volume_to_numpy<float,float>;
%template(vol2npy) swigTom::volume_to_numpy<std::complex<float>,float>;
%template(npy2vol) swigTom::numpy_to_volume<float,float>;
%template(npy2vol) swigTom::numpy_to_volume<std::complex<float>,float>;

//%template(vol2npy) swigTom::volume_to_numpy<double,double>;
//%template(vol2npy) swigTom::volume_to_numpy<std::complex<double>,double>;
//%template(npy2vol) swigTom::numpy_to_volume<double,double>;
//%template(npy2vol) swigTom::numpy_to_volume<std::complex<double>,double>;