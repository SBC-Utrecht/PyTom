%module pytom_mpi
%include "std_string.i"
%include "exception.i"

%exception {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

%{
#include "swigMPI.hpp"
#include "../src/swigMPI.cpp"
%}


%include "swigMPI.hpp"


