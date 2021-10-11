%module pytom_freqweight
%include "std_string.i"
%include "exception.i"

%{
#include "swigFreqWeight.hpp"
#include "../src/swigFreqWeight.cpp"
%}

%include "swigFreqWeight.hpp"
%include "swigVolume.hpp"
%include "../src/swigFreqWeight.cpp"



%exception {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}


%template(weight) swigTom::swigFreqWeight<float>;

//%template(weight) swigTom::swigFreqWeight<double>;
