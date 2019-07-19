/****************************************************************************//**
 * \file swigVolume.hpp
 * \brief Contains fft related functions.
 * \author  Thomas Hrabe
 * \version 0.2
 * \date    16.12.2008
 *******************************************************************************/

#ifndef SWIGFFTFNC_HPP_
#define SWIGFFTFNC_HPP_


#include <swigVolume.hpp>

namespace swigTom{

/****************************************************************************//**
 * 	\brief	Performs inplace fftshift on volume
 * 	\param[in] vol The volume
 *******************************************************************************/
template<typename T,typename TSCALE_SHIFT>
void fftShift(swigVolume<T,TSCALE_SHIFT>& vol, bool isIfftShift);


/****************************************************************************//**
 * 	\brief	Copies a reduced "complex" volume to full size
 * 	\param[in] source The source volume
 * 	\param[in,out] destination The destination volume
 *******************************************************************************/
//template<typename T,typename TSCALE_SHIFT>
//swigVolume<T,TSCALE_SHIFT> reducedToFull(const swigVolume<T,TSCALE_SHIFT>& source);

/****************************************************************************//**
 * 	\brief	Copies a full "complex" volume to reduced size
 * 	\param[in] source The source volume
 * 	\param[out] destination The destination volume
 *******************************************************************************/
//template<typename T,typename TSCALE_SHIFT>
//swigVolume<T,TSCALE_SHIFT> fullToReduced(const swigVolume<T,TSCALE_SHIFT>& source);


#endif /* SWIGFFTFNC_HPP_ */
