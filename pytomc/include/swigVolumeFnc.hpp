/****************************************************************************//**
 * \file swigVolumeFnc.cpp
 * \brief Header file for swigVolume functions
 * \author  Thomas Hrabe
 * \version 0.2
 * \date    1.12.2008
 *******************************************************************************/

#ifndef SWIGVOLUMEFNC_HPP_
#define SWIGVOLUMEFNC_HPP_



#include <swigVolume.hpp>
#include <tom/transf/transform.hpp>

namespace swigTom{
	template<typename T,typename TSCALE_SHIFT>
	swigTom::swigVolume<T,TSCALE_SHIFT> read(std::string fileName);

	template<typename T,typename TSCALE_SHIFT>
	swigTom::swigVolume<T,TSCALE_SHIFT> read(std::string fileName,std::size_t subregion1,std::size_t subregion2,std::size_t subregion3,
																  std::size_t subregion4,std::size_t subregion5,std::size_t subregion6,
																  std::size_t sampling1,std::size_t sampling2,std::size_t sampling3,
																  std::size_t binning1,std::size_t binning2,std::size_t binning3);

	template<typename T,typename TSCALE_SHIFT>
	void conj_mult(swigTom::swigVolume<T,TSCALE_SHIFT>& v1,const swigTom::swigVolume<T,TSCALE_SHIFT>& v2);

	template<typename T,typename TSCALE_SHIFT>
	void initSphere(swigVolume<T,TSCALE_SHIFT>& vol ,float 	radius,float 	sigma, float 	max_radius, float centerx, float centery, float centerz);

	template<typename T,typename TSCALE_SHIFT>
	void maksNorm(swigVolume<T,TSCALE_SHIFT>& vol, swigVolume<T,TSCALE_SHIFT>& mask,bool is_boolean_mask);

	template<typename T,typename TSCALE_SHIFT>
	void transform(swigVolume<T,TSCALE_SHIFT>& src,swigVolume<T,TSCALE_SHIFT>& dst,double phi,double psi,double theta,double centerX,double centerY,double centerZ,double preX,double preY,double preZ,double postX,double postY,double postZ);

	template<typename T,typename TSCALE_SHIFT>
	void rotate(swigVolume<T,TSCALE_SHIFT>& src,swigVolume<T,TSCALE_SHIFT>& dst,double phi,double psi,double theta);

	template<typename T,typename TSCALE_SHIFT>
	void transformCubic(swigVolume<T,TSCALE_SHIFT>& src,swigVolume<T,TSCALE_SHIFT>& dst,double phi,double psi,double theta,double centerX,double centerY,double centerZ,double preX,double preY,double preZ,double postX,double postY,double postZ);

	template<typename T,typename TSCALE_SHIFT>
	void rotateCubic(swigVolume<T,TSCALE_SHIFT>& src,swigVolume<T,TSCALE_SHIFT>& dst,double phi,double psi,double theta);

	template<typename T,typename TSCALE_SHIFT>
	void transformSpline(swigVolume<T,TSCALE_SHIFT>& src,swigVolume<T,TSCALE_SHIFT>& dst,double phi,double psi,double theta,double centerX,double centerY,double centerZ,double preX,double preY,double preZ,double postX,double postY,double postZ);

	template<typename T,typename TSCALE_SHIFT>
	swigVolume<T,TSCALE_SHIFT> transformFourierSpline(swigVolume<T,TSCALE_SHIFT>& src,double phi,double psi,double theta,double shiftX,double shiftY,double shiftZ);

	template<typename T,typename TSCALE_SHIFT>
	void rotateSpline(swigVolume<T,TSCALE_SHIFT>& src,swigVolume<T,TSCALE_SHIFT>& dst,double phi,double psi,double theta);

    template<typename T,typename TSCALE_SHIFT>
	swigVolume<T,TSCALE_SHIFT> rotateSplineInFourier(swigVolume<T,TSCALE_SHIFT>& src,double phi,double psi,double theta);
    
	template<typename T,typename TSCALE_SHIFT>
	void shift(swigVolume<T,TSCALE_SHIFT>& src,swigVolume<T,TSCALE_SHIFT>& dst,double shiftX,double shiftY,double shiftZ);

	template<typename T,typename TSCALE_SHIFT>
	void shiftFourier(swigVolume<T,TSCALE_SHIFT>& src,swigVolume<T,TSCALE_SHIFT>& dst,double preX,double preY,double preZ);

	template<typename T>
	tom::st_idx peak(const swigVolume<T,T>& volume);

	template<typename T>
	tom::st_idx peak(const swigVolume<T,T>& volume,const swigVolume<T,T>& mask);

	template<typename T>
	void conjugate(swigVolume<std::complex<T>,T>& volume);

	template<typename T,typename TSHIFT_SCALE>
	void pasteIn(const swigVolume<T,TSHIFT_SCALE>& volume,swigVolume<T,TSHIFT_SCALE>& destination,signed long x,signed long y,signed long z);

	template<typename T,typename TSHIFT_SCALE>
	void pasteCenter(const swigVolume<T,TSHIFT_SCALE>& volume,swigVolume<T,TSHIFT_SCALE>& destination);

	template<typename T, typename T2>
	void power(swigVolume<T,T2>& volume,const T2 & exponent);

	template<typename T>
	std::size_t numberSetVoxels(const swigVolume<T,T>& volume);

	template<typename T, typename TSCALE_SHIFT>
	swigVolume<T,TSCALE_SHIFT> abs(const swigVolume<T,TSCALE_SHIFT>& volume);

	template<typename T,typename TSCALE_SHIFT>
	T sum(const swigVolume<T,TSCALE_SHIFT>& volume);

	template<typename T,typename TSCALE_SHIFT>
	swigVolume<T,TSCALE_SHIFT> projectSum(const swigVolume<T,TSCALE_SHIFT>& volume, int axis);

	template<typename T,typename TSCALE_SHIFT>
	double mean(const swigVolume<T,TSCALE_SHIFT>& volume);

	template<typename T,typename TSCALE_SHIFT>
	double variance(const swigVolume<T,TSCALE_SHIFT>& volume,bool use_sample_standard_deviation);

	template<typename T,typename TSCALE_SHIFT>
	T min(const swigVolume<T,TSCALE_SHIFT>& volume);

	template<typename T,typename TSCALE_SHIFT>
	T max(const swigVolume<T,TSCALE_SHIFT>& volume);

	template<typename T,typename TSCALE_SHIFT>
	void limit(swigVolume<T,TSCALE_SHIFT>& volume, T lowerBound, T lowerReplacement, T upperBound, T upperReplacement,bool doLower, bool doUpper);

	template<typename T,typename TSCALE_SHIFT>
	swigTom::swigVolume<T,TSCALE_SHIFT> reducedToFull(swigVolume<T,TSCALE_SHIFT>& volume);

	template<typename T, typename TSCALE_SHIFT>
    swigVolume<T,TSCALE_SHIFT> real(const swigVolume<std::complex<T>,TSCALE_SHIFT> &vol);

	template<typename T, typename TSCALE_SHIFT>
    swigVolume<T,TSCALE_SHIFT> imag(const swigVolume<std::complex<T>,TSCALE_SHIFT> &vol);

	template<typename T,typename TSCALE_SHIFT>
	void gaussianNoise(swigVolume<std::complex<T>,TSCALE_SHIFT>& volume,double mean, double sigma);

	template<typename T>
	swigVolume<std::complex<T>,T > complexDiv(const swigVolume<std::complex<T>,T> &volume, const swigVolume<T,T> &div);

	template<typename T>
	swigVolume<T,T> vectorize(swigVolume<T,T> &source);

	template<typename T,typename TSHIFT_SCALE>
	swigVolume<T,TSHIFT_SCALE> getSubregion(swigVolume<T,TSHIFT_SCALE> &source,std::size_t startX,std::size_t startY,std::size_t startZ,std::size_t endX,std::size_t endY,std::size_t endZ);

	template<typename T,typename TSHIFT_SCALE>
	void putSubregion(swigVolume<T,TSHIFT_SCALE> &source, swigVolume<T,TSHIFT_SCALE> &destination,std::size_t positionX,std::size_t positionY,std::size_t positionZ);

	template<typename T>
	void writeSubregion(swigVolume<T,T> &source, std::string filename, std::size_t positionX, std::size_t positionY, std::size_t positionZ);

	template<typename T,typename TSHIFT_SCALE>
	void updateResFromIdx(swigVolume<T,TSHIFT_SCALE>& resV,const swigVolume<T,TSHIFT_SCALE>& newV, swigVolume<T,TSHIFT_SCALE>& orientV, const std::size_t orientIdx);

	template<typename T,typename TSHIFT_SCALE>
	void updateResFromVol(swigVolume<T,TSHIFT_SCALE>& resV,const swigVolume<T,TSHIFT_SCALE>& newV, swigVolume<T,TSHIFT_SCALE>& orientV, const swigVolume<T,TSHIFT_SCALE>& neworientV);

	template<typename T,typename TSHIFT_SCALE>
	void mirrorVolume(swigVolume<T,TSHIFT_SCALE>& src, swigVolume<T,TSHIFT_SCALE>& des);

	template<typename T,typename TSHIFT_SCALE>
	void rescale(const swigVolume<T,TSHIFT_SCALE>& source,swigVolume<T,TSHIFT_SCALE>& destination);

	template<typename T,typename TSHIFT_SCALE>
	void rescaleCubic(const swigVolume<T,TSHIFT_SCALE>& source,swigVolume<T,TSHIFT_SCALE>& destination);

	template<typename T,typename TSHIFT_SCALE>
	void rescaleSpline(const swigVolume<T,TSHIFT_SCALE>& source,swigVolume<T,TSHIFT_SCALE>& destination);

	template<typename T,typename TSHIFT_SCALE>
	T interpolate(const swigVolume<T,TSHIFT_SCALE>& source,const float x,const float y,const float z);

	template<typename T,typename TSHIFT_SCALE>
	T interpolateCubic(const swigVolume<T,TSHIFT_SCALE>& source,const float x,const float y,const float z);

	template<typename T,typename TSHIFT_SCALE>
	T interpolateSpline(const swigVolume<T,TSHIFT_SCALE>& source,const float x,const float y,const float z);

	template<typename T,typename TSHIFT_SCALE>
	T interpolateFourierSpline(const swigVolume<T,TSHIFT_SCALE>& source,const float x,const float y,const float z);

	template<typename T,typename TSCALE_SHIFT>
	void backProject(swigVolume<T,TSCALE_SHIFT>& src,
					 swigVolume<T,TSCALE_SHIFT>& dst,
					 swigVolume<T,TSCALE_SHIFT>& phi,
					 swigVolume<T,TSCALE_SHIFT>& theta,
					 swigVolume<T,TSCALE_SHIFT>& offset,
					 swigVolume<T,TSCALE_SHIFT>& offsetProjections);

	template<typename T,typename TSCALE_SHIFT>
	swigVolume<std::complex<T>, TSCALE_SHIFT> complexRealMult(const swigVolume<std::complex<T>,TSCALE_SHIFT> &vol,
															  const swigVolume<T,TSCALE_SHIFT> &otherVol);

	template<typename T,typename TSCALE_SHIFT>
	void backProjectExtended(swigVolume<T,TSCALE_SHIFT>& src,
					 swigVolume<T,TSCALE_SHIFT>& dst,
					 swigVolume<T,TSCALE_SHIFT>& phi,
					 swigVolume<T,TSCALE_SHIFT>& theta,
					 swigVolume<T,TSCALE_SHIFT>& psi,
					 swigVolume<T,TSCALE_SHIFT>& offset,
					 swigVolume<T,TSCALE_SHIFT>& offsetProjections);
    
    template<typename T, typename TSCALE_SHIFT>
    swigVolume<std::complex<T>,TSCALE_SHIFT> mergeRealImag(const swigVolume<T,TSCALE_SHIFT> &real,const swigVolume<T,TSCALE_SHIFT> &imag);

    template<typename T,typename TSCALE_SHIFT>
	swigTom::swigVolume<T,TSCALE_SHIFT> fullToReduced(swigVolume<T,TSCALE_SHIFT>& volume);
}

#endif /* SWIGVOLUMEFNC_HPP_ */
