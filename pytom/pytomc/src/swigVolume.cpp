/*
 * swigVolume.cpp
 *
 *  Created on: Nov 27, 2008
 *      Author: hrabe
 */

#include <swigVolume.hpp>

namespace swigTom{

	/**
	 *	\brief Element wise + operator for swigVolume
	 */
	template<typename T, typename TSCALE_SHIFT>
	swigVolume<T,TSCALE_SHIFT> swigVolume<T,TSCALE_SHIFT>::operator+(const swigVolume<T,TSCALE_SHIFT> &otherVol) const{
		swigVolume<T,TSCALE_SHIFT> newVol(*this);
		tom::element_wise_add(newVol, otherVol);
		return newVol;
	}
	/**
	 *	\brief + operator for swigVolume
	 */
	template<typename T, typename TSCALE_SHIFT>
	swigVolume<T,TSCALE_SHIFT> swigVolume<T,TSCALE_SHIFT>::operator+(const TSCALE_SHIFT &value) const{
		swigVolume<T,TSCALE_SHIFT> newVol(*this);
		newVol.shift_scale(value,(TSCALE_SHIFT) 1);
		return newVol;
	}

	/**
	 *	\brief Element wise * operator for swigVolume
	 */
	template<typename T, typename TSCALE_SHIFT>
	swigVolume<T,TSCALE_SHIFT> swigVolume<T,TSCALE_SHIFT>::operator*(const swigVolume<T,TSCALE_SHIFT> &otherVol) const{
		swigVolume<T,TSCALE_SHIFT> newVol(*this);
		tom::element_wise_multiply(newVol, otherVol);
		return newVol;
	}


	/**
	 *	\brief * operator for swigVolume
	 */
	template<typename T, typename TSCALE_SHIFT>
	swigVolume<T,TSCALE_SHIFT> swigVolume<T,TSCALE_SHIFT>::operator*(const TSCALE_SHIFT &value) const{
		swigVolume<T,TSCALE_SHIFT> newVol(*this);
		newVol.shift_scale((TSCALE_SHIFT) 0,value);
		return newVol;
	}
/*
    template<typename T, typename TSCALE_SHIFT>
	swigVolume<std::complex<T>,TSCALE_SHIFT> swigVolume<T,TSCALE_SHIFT>::operator*(const std::complex<T> &value) const{
        
		swigVolume<std::complex<T>,TSCALE_SHIFT> complexVol = swigVolume<std::complex<T>,TSCALE_SHIFT>(this->size_x(),this->size_y(),this->size_z());
        
        T val=0;
        
        for(std::size_t x = 0; x<this->size_x();x++){
            for(std::size_t y = 0; y<this->size_y();y++){
                for(std::size_t z = 0; z<this->size_z();z++){
                    val = this->get(x,y,z);
                    complexVol(val*value,x,y,z);
                }
            }
        }
        
        return complexVol;
	}
*/    
	/**
	 *	\brief Element wise - operator for swigVolume
	 */
	template<typename T, typename TSCALE_SHIFT>
	swigVolume<T,TSCALE_SHIFT> swigVolume<T,TSCALE_SHIFT>::operator-(const swigVolume<T,TSCALE_SHIFT> &otherVol) const{
		swigVolume<T,TSCALE_SHIFT> newVol(*this);
		tom::element_wise_sub(newVol, otherVol);
		return newVol;
	}

	/**
	 *	\brief * operator for swigVolume
	 */
	template<typename T, typename TSCALE_SHIFT>
	swigVolume<T,TSCALE_SHIFT> swigVolume<T,TSCALE_SHIFT>::operator-( const TSCALE_SHIFT &value) const{
		swigVolume<T,TSCALE_SHIFT> newVol(*this);
		newVol.shift_scale(-value,(TSCALE_SHIFT) 1);
		return newVol;
	}
	/**
	 *	\brief Element wise / operator for swigVolume
	 */
	template<typename T, typename TSCALE_SHIFT>
	swigTom::swigVolume<T,TSCALE_SHIFT> swigTom::swigVolume<T,TSCALE_SHIFT>::operator/(const swigVolume<T,TSCALE_SHIFT> &otherVol) const{
		swigVolume<T,TSCALE_SHIFT> newVol(*this);
		tom::element_wise_div(newVol, otherVol,(T)0);
		return newVol;
	}
	/**
	 *	\brief Element wise / operator for swigVolume

	template<typename T, typename TSCALE_SHIFT>
	swigTom::swigVolume<std::complex<T>,TSCALE_SHIFT> swigTom::swigVolume<std::complex<T>,TSCALE_SHIFT>::operator/(const swigVolume<TSCALE_SHIFT,TSCALE_SHIFT> &otherVol) const{
		swigVolume<T,TSCALE_SHIFT> newVol(*this);
		tom::element_wise_div(newVol, otherVol,(std::complex<T>)0);
		return newVol;
	}
	*/

	/**
	 *	\brief * operator for swigVolume
	 */
	template<typename T, typename TSCALE_SHIFT>
	swigVolume<T,TSCALE_SHIFT> swigVolume<T,TSCALE_SHIFT>::operator/( const TSCALE_SHIFT &value) const{
		swigVolume<T,TSCALE_SHIFT> newVol(*this);
		if(value ==(TSCALE_SHIFT) 0.0)
			throw std::runtime_error("Division by zero! Abort!");
		
		newVol.shift_scale((TSCALE_SHIFT) 0.0,(TSCALE_SHIFT)1/value);
		return newVol;
	}

	template<typename T, typename TSCALE_SHIFT>
	const T  swigVolume<T,TSCALE_SHIFT>::operator()(const std::size_t &x,const std::size_t &y,const std::size_t &z){
		//not tested yet
		return this->get(x,y,z);
	}
	template<typename T, typename TSCALE_SHIFT>
	void swigVolume<T,TSCALE_SHIFT>::operator()(const T &value,const std::size_t &x,const std::size_t &y,const std::size_t &z){
		this->get(x,y,z) = value;
	}

}



