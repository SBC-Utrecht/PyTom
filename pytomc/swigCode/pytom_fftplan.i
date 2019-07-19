%module pytom_fftplan
%include "std_string.i"
%include "exception.i"

%{
#include "swigFftPlan.hpp"
#include "swigVolume.hpp"
#include "../src/swigFftFnc.cpp"
%}

%include "swigFftPlan.hpp"
%include "swigVolume.hpp"
%include "swigFftFnc.cpp"

%exception {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}


%template(plan) swigTom::swigFftPlan<float,float>;
%template(fftShift) swigTom::fftShift<float, float>;
%template(fftShift) swigTom::fftShift<std::complex<float>, float>;
//%template(reducedToFull) swigTom::reducedToFull<float, float>;
//%template(reducedToFull) swigTom::reducedToFull<std::complex<float>, float>;
//%template(fullToReduced) swigTom::fullToReduced<float, float>;


//%template(plan) swigTom::swigFftPlan<double,double>;
//%template(fftShift) swigTom::fftShift<double, double>;
//%template(fftShift) swigTom::fftShift<std::complex<double>, double>;
//%template(reducedToFull) swigTom::reducedToFull<double, double>;
//%template(reducedToFull) swigTom::reducedToFull<std::complex<double>, double>;
//%template(fullToReduced) swigTom::fullToReduced<double, double>;
