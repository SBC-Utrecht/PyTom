/*
 * swigNumpy.hpp
 *
 *  Created on: Jan 20, 2009
 *      Author: hrabe
 */

#ifndef SWIGTOMNUMPY_HPP_
#define SWIGTOMNUMPY_HPP_

#include <Python.h>
#include <swigVolume.hpp>
#include <ndarrayobject.h>


namespace swigTom{
	template<typename T, typename TSHIFT_SCALE>
	PyObject* volume_to_numpy(const swigTom::swigVolume<T,TSHIFT_SCALE>& vol);

	template<typename T,typename TSHIFT_SCALE>
	swigTom::swigVolume<T,TSHIFT_SCALE>* numpy_to_volume(const PyObject* array,const std::size_t dim);

}
#endif /* SWIGTOMNUMPY_HPP_ */

