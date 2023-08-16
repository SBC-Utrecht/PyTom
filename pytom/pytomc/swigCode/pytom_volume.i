%module pytom_volume
%include "std_string.i"
%include "std_complex.i"
%include "exception.i"

%{
#include "swigVolume.hpp"
#include "../src/swigVolume.cpp"
#include "../src/swigVolumeFnc.cpp"
#include <exception>
#include <iostream>
#include <complex>


%}

%ignore swigTom::swigVolume<T>::swigVolume(tom::Volume<T>& v);
%include "swigVolume.hpp"



%include "swigVolume.cpp"
%include "swigVolumeFnc.cpp"




%exception {
  try {
    $action
  } catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}

%typemap(out) tom::st_idx{
	PyObject * list = PyList_New(0);
	
	PyList_Insert(list, 0, PyInt_FromLong($1.x));
	PyList_Insert(list, 1, PyInt_FromLong($1.y));
	PyList_Insert(list, 2, PyInt_FromLong($1.z));
	
	$result = list;
}
%typemap(out) std::tuple<std::size_t, std::size_t, std::size_t>{
       PyObject * tuple = PyTuple_New(3);
       
       PyTuple_SetItem(tuple, 0, PyInt_FromLong(std::get<0>($1)));
       PyTuple_SetItem(tuple, 1, PyInt_FromLong(std::get<1>($1)));
       PyTuple_SetItem(tuple, 2, PyInt_FromLong(std::get<2>($1)));
       
       $result = tuple;
}
%nodefaultdtor swigTom::swigVolume;


%feature("autodoc", "1");
%template(vol) swigTom::swigVolume<float,float>;
%template(vol_comp) swigTom::swigVolume<std::complex<float>,float>;
%template(read) swigTom::read<float,float>;
%template(conj_mult) swigTom::conj_mult<std::complex<float>,float>;
%template(initSphere) swigTom::initSphere<float,float>;
%template(maskNorm) swigTom::maskNorm<float,float>;
%template(rotate) swigTom::rotate<float,float>;
%template(rotateCubic) swigTom::rotateCubic<float,float>;
%template(rotateSpline) swigTom::rotateSpline<float,float>;
%template(rotateSplineInFourier) swigTom::rotateSplineInFourier<std::complex<float >,float>;
%template(shift) swigTom::shift<float,float>;
%template(shiftFourier) swigTom::shiftFourier<std::complex<float>,float>;
%template(transform) swigTom::transform<float,float>;
%template(general_transform) swigTom::general_transform<float,float>;
%template(transformCubic) swigTom::transformCubic<float,float>;
%template(transformSpline) swigTom::transformSpline<float,float>;
%template(transformFourierSpline) swigTom::transformFourierSpline<std::complex<float >,float>;
%template(peak) swigTom::peak<float>;
%template(conjugate) swigTom::conjugate<float>;
%template(paste) swigTom::pasteIn<float,float>;
%template(pasteCenter) swigTom::pasteCenter<float,float>;
%template(power) swigTom::power<float,float>;
%template(power) swigTom::power<std::complex<float> ,float>;
%template(numberSetVoxels) swigTom::numberSetVoxels<float>;
%template(sum) swigTom::sum<float,float>;
%template(sum) swigTom::sum<std::complex<float>,float>;
%template(projectSum) swigTom::projectSum<float,float>;
%template(mean) swigTom::mean<float,float>;
%template(variance) swigTom::variance<float,float>;
%template(min) swigTom::min<float,float>;
%template(max) swigTom::max<float,float>;
%template(limit) swigTom::limit<float,float>;
%template(abs) swigTom::abs<float,float>;
%template(abs) swigTom::abs<std::complex<float>,float>;
%template(reducedToFull) swigTom::reducedToFull<float,float>;
%template(reducedToFull) swigTom::reducedToFull<std::complex<float>,float>;
%template(gaussianNoise) swigTom::gaussianNoise<float,float>;
%template(complexDiv) swigTom::complexDiv<float>;
%template(vectorize) swigTom::vectorize<float>;
%template(subvolume) swigTom::getSubregion<float,float>;
%template(subvolume) swigTom::getSubregion<std::complex<float>,float >;
%template(putSubVolume) swigTom::putSubregion<float,float>;
%template(putSubVolume) swigTom::putSubregion<std::complex<float>,float>;
%template(writeSubregion) swigTom::writeSubregion<float>;
%template(updateResFromIdx) swigTom::updateResFromIdx<float, float>;
%template(updateResFromVol) swigTom::updateResFromVol<float, float>;
%template(mirrorVolume) swigTom::mirrorVolume<float, float>;
%template(rescale) swigTom::rescale<float, float>;
%template(rescaleCubic) swigTom::rescaleCubic<float, float>;
%template(rescaleSpline) swigTom::rescaleSpline<float, float>;
%template(interpolate) swigTom::interpolate<float, float>;
%template(interpolateCubic) swigTom::interpolateCubic<float, float>;
%template(interpolateSpline) swigTom::interpolateSpline<float, float>;
%template(interpolateFourierSpline) swigTom::interpolateFourierSpline<std::complex<float>, float>;
%template(backProject) swigTom::backProject<float,float>;
%template(complexRealMult) swigTom::complexRealMult<float,float>;
%template(real) swigTom::real<float,float>;
%template(imag) swigTom::imag<float,float>;
%template(mergeRealImag) swigTom::mergeRealImag<float,float>;
%template(volTOsf) swigTom::volTOsf<float,float>;
%template(fvolTOsf) swigTom::fvolTOsf<float,float>;
%template(sfFourierShift) swigTom::sfFourierShift<float,float>;
%template(fullToReduced) swigTom::fullToReduced<float,float>;
%template(fullToReduced) swigTom::fullToReduced<std::complex<float>,float>;
