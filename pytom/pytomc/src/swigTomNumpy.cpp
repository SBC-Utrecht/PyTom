/*
 * swigTomNumpy.cpp
 *
 *  Created on: Jan 20, 2009
 *      Author: hrabe
 */

#include <swigTomNumpy.hpp>
#include <swigVolume.hpp>
#include <tom/volume_fcn.hpp>
#include <Python.h>
#include <ndarrayobject.h>
#include <exception>
#include <iostream>

#define FORTRAN_ORDER 1  //corresponds to F order ( according to numpy naming)
#define C_ORDER 0 //corresponds to C order ( according to numpy naming)

namespace swigTom{
template<typename T,typename TSHIFT_SCALE>
PyObject* volume_to_numpy(const swigTom::swigVolume<T,TSHIFT_SCALE>& vol){

	//PyArray_API is in __multiarray_api.h of numpy
	if(PyArray_API == NULL){
		int impArrReturn =_import_array();

		if(impArrReturn < 0){
			PyErr_Print();
			throw std::runtime_error("Fatal error: could not initialise numpy. Array exchange will not be possible!");
		}
	}

	std::size_t dimensions = 0;
	if(vol.size_x() > 1)
		dimensions = 1;
	if(vol.size_y() > 1)
		dimensions = 2;
	if(vol.size_z() > 1)
		dimensions = 3;
	if(dimensions == 0)
		throw std::runtime_error("Volume has size 0,0,0. Aborting!");

	npy_intp* dims = (npy_intp*)PyMem_Malloc(sizeof(npy_intp)*dimensions);
	npy_intp* strides = (npy_intp*)PyMem_Malloc(sizeof(npy_intp)*dimensions);


	if(dimensions >= 1){
		dims[0] = vol.size_x();
		strides[0] = vol.strideX();
	}
	if(dimensions >= 2){
		dims[1] = vol.size_y();
		strides[1] = vol.strideY();
	}
	if(dimensions >= 3){
		dims[2] = vol.size_z();
		strides[2] = vol.strideZ();
	}

	char type = 0;

	if(tom::is_float<T>())
		type = NPY_FLOATLTR;
	else if(tom::is_double<T>())
		type = NPY_DOUBLELTR;
	else if(tom::is_float_complex<T>())
		type = NPY_CFLOATLTR;
	else if(tom::is_double_complex<T>())
		type = NPY_CDOUBLELTR;


	PyArray_Descr *descr = PyArray_DescrFromType(type);
	PyObject* obj = PyArray_NewFromDescr(&PyArray_Type,descr,
			dimensions, dims,
			strides, (char*) &vol.get(),
			NPY_CARRAY, NULL);

	return obj;
}

template<typename T,typename TSHIFT_SCALE>
swigTom::swigVolume<T,TSHIFT_SCALE>* numpy_to_volume(const PyObject* array,const std::size_t dim){
	npy_intp* dims = PyArray_DIMS(array);
	npy_intp* strides = PyArray_STRIDES(array);
    
    short order = FORTRAN_ORDER; 
    if(dim == 3){
        
        if(strides[0] > strides[2]){
            order = C_ORDER;
        }
        
    }
    
	if(dim == 2){
		dims[2] = 0;
		strides[2] = 0;
        
        if(strides[0] > strides[1]){
            order = C_ORDER;
        }
	}
	if(dim == 1){
		dims[1] = 0;
		strides[1] = 0;
	}
	if(dim <1)
		throw std::runtime_error("Volume has 0 dimension. Aborting!");

    
    if(order == FORTRAN_ORDER){
    	// create the volume with existing data. will not free the data.
        // tom::Volume<T> tomVol((T*) PyArray_DATA(array), dims[0], dims[1], dims[2],strides[0],strides[1],strides[2], false, NULL);
        // return new swigTom::swigVolume<T,TSHIFT_SCALE>(tomVol); // copy the memory, create a new one and return

        // share the memory and let the pytom volume free the memory
        return new swigTom::swigVolume<T,TSHIFT_SCALE>((T*) PyArray_DATA(array), dims[0], dims[1], dims[2],strides[0],strides[1],strides[2]);
    }else{
        throw std::runtime_error("Your numpy array is in C order. Use Fortran order instead!");
    }
}

}


