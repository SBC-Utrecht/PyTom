/*
 * swigFftFnc.cpp
 *
 *  Created on: Dec 16, 2008
 *      Author: hrabe
 */
#include <tom/volume.hpp>
#include <tom/volume_fcn.hpp>
#include <math.h>

namespace swigTom{

template<typename T,typename TSCALE_SHIFT>
void fftShift(swigVolume<T,TSCALE_SHIFT>& vol, bool isIfftShift){
	tom::fftshift(vol,isIfftShift);
}
/*
template<typename T,typename TSCALE_SHIFT>
swigVolume<T,TSCALE_SHIFT> reducedToFull(const swigVolume<T,TSCALE_SHIFT>& source){
	swigVolume<T,TSCALE_SHIFT> destination(source.getSizeX(),source.getSizeY(),(source.getSizeZ()-1)*2);
	tom::hermitian_symmetry_to_full(source,destination);
	return destination;
}

template<typename T,typename TSCALE_SHIFT>
swigVolume<T,TSCALE_SHIFT> fullToReduced(const swigVolume<T,TSCALE_SHIFT>& source){

	std::size_t x,y,z;

	x = source.getSizeX();
	y = source.getSizeY();
	z = source.getSizeZ()/2+1;

	tom::Volume<T> destination(x,y,z,NULL,NULL);

	destination.setValues(tom::Volume<T>(const_cast<T *>(&source.get()),x,y,z,source.getStrideX(),source.getStrideY(),source.getStrideZ(),false,NULL));

	return swigTom::swigVolume<T,TSCALE_SHIFT>(destination);

}
*/
}
