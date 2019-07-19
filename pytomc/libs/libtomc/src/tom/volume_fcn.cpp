#include <tom/volume_fcn.hpp>

#include <limits>

#include <cmath>
#include <complex>
#include <cassert>
#include <ctime>



#include <boost/lambda/lambda.hpp>
#include <boost/lambda/if.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/casts.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include <boost/type_traits/is_floating_point.hpp>

#include <tom/fftw/fftw_plan.hpp>

#include <assert.h>

#include <iostream>

#define PI 3.141592653589793238




namespace {
template<typename T, typename TMASK, typename TPRECISION>
struct for_each__tom__norm_mask__get_mean {
    for_each__tom__norm_mask__get_mean(): n(0), sum(0) { }
    std::size_t n; TPRECISION sum;
    inline void operator()(const T &v, const TMASK &mask, std::size_t, std::size_t, std::size_t) {
        if (mask) {
            this->n ++;
            this->sum += v;
        }
    }
};
template<typename T, typename TMASK, typename TPRECISION>
struct for_each__tom__norm_mask__apply_mask_no_variance {
    for_each__tom__norm_mask__apply_mask_no_variance(TPRECISION mean): sum(0), mean(mean) { }
    TPRECISION sum;
    TPRECISION const mean;
    inline void operator()(T &v, const TMASK &mask, std::size_t, std::size_t, std::size_t) {
        if (mask) {
            const TPRECISION v2 = (v - this->mean) * mask;
            v = v2;
            this->sum += v2;
        } else {
            v = 0;
        }
    }
};
template<typename T, typename TMASK, typename TPRECISION>
struct for_each__tom__norm_mask__set_mean0 {
    for_each__tom__norm_mask__set_mean0(TPRECISION mean): mean(mean) { }
    const TPRECISION mean;
    inline void operator()(T &v, const TMASK &mask, std::size_t, std::size_t, std::size_t) {
        if (mask) {
            v = v - this->mean;
        }
    }
};
template<typename T, typename TMASK, typename TPRECISION>
struct for_each__tom__norm_mask__set_mean0_all {
    for_each__tom__norm_mask__set_mean0_all(TPRECISION mean): mean(mean) { }
    const TPRECISION mean;
    inline void operator()(T &v, const TMASK &mask, std::size_t, std::size_t, std::size_t) {
        if (mask) {
            v = v - this->mean;
        } else {
            v = 0;
        }
    }
};
template<typename T, typename TMASK, typename TPRECISION>
struct for_each__tom__norm_mask__sub_meanv {
    for_each__tom__norm_mask__sub_meanv(TPRECISION mean): sum(0), sum2(0), mean(mean) { }
    TPRECISION sum, sum2;
    const TPRECISION mean;
    inline void operator()(T &v, const TMASK &mask, std::size_t, std::size_t, std::size_t) {
        if (mask) {
            const TPRECISION v2 = (v - this->mean) * mask;
            v = v2;
            this->sum += v2;
            this->sum2 += v2*v2;
        } else {
            v = 0;
        }
    }
};
template<typename T, typename TMASK, typename TPRECISION>
struct for_each__tom__norm_mask__sub_meanv_bool {
    for_each__tom__norm_mask__sub_meanv_bool(TPRECISION mean): sum(0), sum2(0), mean(mean) { }
    TPRECISION sum, sum2;
    const TPRECISION mean;
    inline void operator()(T &v, const TMASK &mask, std::size_t, std::size_t, std::size_t) {
        if (mask) {
            const TPRECISION v2 = v - this->mean;
            v = v2;
            this->sum += v2;
            this->sum2 += v2*v2;
        } else {
            v = 0;
        }
    }
};
template<typename T, typename TMASK, typename TPRECISION>
struct for_each__tom__norm_under_mask_5 {
    for_each__tom__norm_under_mask_5(TPRECISION mean, TPRECISION stddev): mean(mean), stddev(stddev) { assert(stddev != 0); }
    const TPRECISION mean, stddev;
    inline void operator()(T &v, const TMASK &mask, std::size_t, std::size_t, std::size_t) {
        if (mask) {
            v = (v - this->mean) / this->stddev;
        }
    }
};
} // namespace


template<typename T>
void tom::back_project(tom::Volume<T> &src,
					   tom::Volume<T> &dst,
					   int volumeSizeX, int volumeSizeY, int volumeSizeZ,
					   float Z1, float Y,float Z2,
					   int projectionSizeX, int projectionSizeY,
					   float x_offset, float y_offset, float z_offset,
					   float projectionOffsetX,float projectionOffsetY){
	
	float Tr11, Tr12, Tr13;
	float Tr21, Tr22, Tr23;
	float Tr31, Tr32, Tr33;
	
	float interpolationPositionX,interpolationPositionY,interpolationOffsetX,interpolationOffsetY;
	
	int reconstructionCentreX,reconstructionCentreY,reconstructionCentreZ;
	int projectionCentreX,projectionCentreY;
	int projectionPositionX,projectionPositionY;
	
	/*  Calculate Angles to Radians*/
	//std::cout << Z1 << " " << Y << " " << Z2 << std::endl;
	Z1=(float)tom::math::deg2rad((double)Z1);
	Y=(float)tom::math::deg2rad((double)Y);
	Z2=(float)tom::math::deg2rad((double)Z2);
	//std::cout << Z1 << " " << Y << " " << Z2 << std::endl;
	

	//first row of ZYZ Matrix
	Tr11 = (float)(cos(Y)*cos(Z1)*cos(Z2)-sin(Z1)*sin(Z2));
	Tr21 = (float)(cos(Y)*sin(Z1)*cos(Z2)+cos(Z1)*sin(Z2));
	Tr31 = (float)(-sin(Y)*cos(Z2));
	//std::cout << Tr11 << " " << Tr21 << " " << Tr31 << std::endl;
	//second row of ZYZ Matrix
	Tr12 = (float)(-cos(Y)*cos(Z1)*sin(Z2)-sin(Z1)*cos(Z2));
	Tr22 = (float)(-cos(Y)*sin(Z1)*sin(Z2)+cos(Z1)*cos(Z2));
	Tr32 = (float)(sin(Y)*sin(Z2));
	//std::cout << Tr12 << " " << Tr22 << " " << Tr32 << std::endl;	
	//third row of ZYZ Matrix
	Tr13 = (float)(sin(Y)*cos(Z1));
	Tr23 = (float)(sin(Y)*sin(Z1));
	Tr33 = (float)(cos(Y));
	
	/*
	//Reconstruction works with these settings!
	Tr11 = (float)cos(Y)*(float)cos(Z1);
	Tr21 = (float)cos(Y)*(float)sin(Z1);
	Tr31 = -(float)sin(Y);
	std::cout << Tr11 << " " << Tr21 << " " << Tr31 << std::endl;
	Tr12 = -(float)sin(Z1);
	Tr22 = (float)cos(Z1);
	Tr32 = 0;
	std::cout << Tr12 << " " << Tr22 << " " << Tr32 << std::endl;
	Tr13 = 0;
	Tr23 = 0;
	Tr33 = 0;
	
	assert(false);
	*/
	/*C  Define centre of the reconstruction body*/
	reconstructionCentreX=volumeSizeX/2+1 - x_offset;
	reconstructionCentreY=volumeSizeY/2+1 - y_offset;
	reconstructionCentreZ=volumeSizeZ/2+1 - z_offset;
	
	projectionCentreX=projectionSizeX/2+1 + projectionOffsetX;
	projectionCentreY=projectionSizeY/2+1 + projectionOffsetY;
	
	/*C
	 C  Main loop over all pixels of the object*/
	
	float valueX1=0.0;
	float valueX2=0.0;
	float valueY=0.0;
	
	for (int indexZ=1;indexZ<=volumeSizeZ;indexZ++){
		
		//transformedZPosition = (indexZ-reconstructionCentreZ)*Tr31;
		
		for (int indexY=1;indexY<=volumeSizeY;indexY++){
			
			//transformedXPosition = (indexY-reconstructionCentreY)*Tr21 + transformedZPosition;
			//transformedYPosition = (indexY-reconstructionCentreY)*Tr22;
			
			for (int indexX=1;indexX<=volumeSizeX;indexX++){
				
				/*C  Calculate point (interpolationPositionX,interpolationPositionY) in the projection plane corresponding to
				 C  the object point (indexX,indexY,indexZ)*/
				//interpolationPositionX = (indexX-reconstructionCentreX)*Tr11 + transformedXPosition + projectionCentreX;
				//interpolationPositionY = (indexX-reconstructionCentreX)*Tr12 + transformedYPosition + projectionCentreY;
				
				interpolationPositionX = Tr11 * (indexX-reconstructionCentreX) + Tr21 * (indexY-reconstructionCentreY) + Tr31 * (indexZ-reconstructionCentreZ) + projectionCentreX;
				interpolationPositionY = Tr12 * (indexX-reconstructionCentreX) + Tr22 * (indexY-reconstructionCentreY) + Tr32 * (indexZ-reconstructionCentreZ) + projectionCentreY;
				
				projectionPositionX = (int)(interpolationPositionX);
				projectionPositionY = (int)(interpolationPositionY);
		
				
				interpolationOffsetX = (interpolationPositionX-projectionPositionX);
				interpolationOffsetY = (interpolationPositionY-projectionPositionY);
				/*C  In case of projection point outside frame of projection do nothing.
				 IF ((projectionPositionX.LT.1).OR.(projectionPositionY.LT.1).OR.(projectionPositionX.GT.projectionSizeX).OR.(projectionPositionY.GT.projectionSizeY)) */
				
				if ( !(projectionPositionX < 1 || projectionPositionY < 1 || projectionPositionX > projectionSizeX || projectionPositionY > projectionSizeY) ){
					/*  Perform bilinear interpolation */
					
					/*  projectionPositionX last column of projection array?*/
					
					if (projectionPositionX==projectionSizeX){
						valueX1 = src.get(projectionPositionX -1, projectionPositionY -1, 0);
					}else{
						valueX1 = src.get(projectionPositionX -1, projectionPositionY -1, 0) + ( src.get(projectionPositionX, projectionPositionY -1, 0) - src.get(projectionPositionX -1, projectionPositionY -1, 0) )*interpolationOffsetX;
					}
					
					/*  projectionPositionY last row of projection array?*/
					
					if (projectionPositionY==projectionSizeY){
						interpolationOffsetY=0.0;
						valueX2=0.0;
					}else if(projectionPositionX==projectionSizeX){
						valueX2 = src.get(projectionPositionX -1, projectionPositionY, 0);
					}else{
						valueX2 = src.get(projectionPositionX -1, projectionPositionY, 0) + ( src.get(projectionPositionX, projectionPositionY, 0) - src.get(projectionPositionX -1, projectionPositionY, 0) )*interpolationOffsetX;
					}

					valueY = valueX1 + (valueX2-valueX1) *interpolationOffsetY;
					
					dst.get(indexX -1, indexY -1, indexZ -1) = dst.get(indexX -1, indexY -1, indexZ -1) + valueY;
					
				}
			}
		}
	}
}



/****************************************************************************//**
 * \brief Normalises the input volume to mean 0 and standard deviation 1 with mask
 *
 * \param[in,out] v The volume to be normalised.
 * \param[in] mask A mask under which \a v is to be normalized. \a mask must have the same
 *   size as the volume \a v. The mask should contain numbers between 0 and 1.
 * \param[in] stddev_type How to finally scale the volume. It can be one of the integer
 *   constants NORM_NO_NORM, NORM_STDDEV_1 and NORM_STDDEV_SAMPLE_1.
 * \param[out] variance If not NULL, the variance of the volume is returned.
 *   If \a stddev_type == NORM_NO_NORM then it is the variance
 *   of the volume after executing the function (with divisor N, not (N-1)). If \a stddev_type is
 *   NORM_STDDEV_1 or NORM_STDDEV_SAMPLE_1 it always 1.
 *
 * All computations are done with the precision of the template type \a TPRECISION. \n
 * What this function does is the following: first the volume \a v is scaled to have
 * a mean of 0 under the \a mask (i.e. the voxels of \a v where the corresponding
 * voxel of the \a mask is != 0). Then each element of the \a v is multiplied with
 * \a mask and the volume is again shifted to have a mean of 0 under the mask.
 * If \a stddev_type == NORM_NO_NORM
 * then the variance is computed and the function exits. Otherwise a final scaling is done
 * to set the standard deviation. If NORM_STDDEV_1 it takes the standard
 * deviation with divisor N (N the number of voxels). If NORM_STDDEV_SAMPLE_1
 * it takes the sample standard deviation with divisor (N-1).
 *******************************************************************************/
template<typename T, typename TMASK, typename TPRECISION>
void tom::norm_mask(tom::Volume<T> &v, const tom::Volume<TMASK> &mask, tom::norm::ntype stddev_type, TPRECISION *variance, bool is_boolean_mask) {

    if (stddev_type != tom::norm::NORM_NO_NORM &&
        stddev_type != tom::norm::NORM_STDDEV_1 &&
        stddev_type != tom::norm::NORM_STDDEV_SAMPLE_1) {
        throw std::invalid_argument("The normalisation type of the deviation is not defined.");
    }


    std::size_t val_n;
    const std::size_t numel = v.getSizeX() * v.getSizeY() * v.getSizeZ();
    TPRECISION val_sum, val_sum2, val_mean, val_stddev, val_variance;

    val_n = 0;
    val_sum = 0.;
    val_sum2 = 0.;
    {
        // First compute the number of elements with a mask-value ~= 0 and its mean.
        for_each__tom__norm_mask__get_mean<T, TMASK, TPRECISION> s;
        tom::loop::for_each<const tom::Volume<T>, const tom::Volume<TMASK>, ::for_each__tom__norm_mask__get_mean<T, TMASK, TPRECISION> &>(v, mask, s);
        val_n = s.n;
        val_sum = s.sum;
    }
    val_mean = val_sum / static_cast<TPRECISION>(val_n);

    if (stddev_type == tom::norm::NORM_NO_NORM && !variance) {
        // First compute the number of elements with a mask-value ~= 0 and its mean.
        if (is_boolean_mask) {
            if (val_mean != 0.) {
                // Finally normalise again under the mask.
                tom::loop::for_each<tom::Volume<T>, const tom::Volume<TMASK>, ::for_each__tom__norm_mask__set_mean0_all<T, TMASK, TPRECISION> >(v, mask, ::for_each__tom__norm_mask__set_mean0_all<T, TMASK, TPRECISION>(val_mean));
            }
        } else {
            for_each__tom__norm_mask__apply_mask_no_variance<T, TMASK, TPRECISION> s(val_mean);
            tom::loop::for_each<tom::Volume<T>, const tom::Volume<TMASK>, ::for_each__tom__norm_mask__apply_mask_no_variance<T, TMASK, TPRECISION> &>(v, mask, s);
            val_sum = s.sum;
            val_mean = val_sum / static_cast<TPRECISION>(val_n);
            if (val_mean != 0.) {
                // Finally normalise again under the mask.
                tom::loop::for_each<tom::Volume<T>, const tom::Volume<TMASK>, ::for_each__tom__norm_mask__set_mean0<T, TMASK, TPRECISION> >(v, mask, ::for_each__tom__norm_mask__set_mean0<T, TMASK, TPRECISION>(val_mean));
            }
        }
    } else {
        // Shift inside the mask to give mean==0 and multiply with the mask.
        if (is_boolean_mask) {
            for_each__tom__norm_mask__sub_meanv_bool<T, TMASK, TPRECISION> s(val_mean);
            tom::loop::for_each<tom::Volume<T>, const tom::Volume<TMASK>, ::for_each__tom__norm_mask__sub_meanv_bool<T, TMASK, TPRECISION> &>(v, mask, s);
            val_sum = s.sum;
            val_sum2 = s.sum2;
        } else {
            for_each__tom__norm_mask__sub_meanv<T, TMASK, TPRECISION> s(val_mean);
            tom::loop::for_each<tom::Volume<T>, const tom::Volume<TMASK>, ::for_each__tom__norm_mask__sub_meanv<T, TMASK, TPRECISION> &>(v, mask, s);
            val_sum = s.sum;
            val_sum2 = s.sum2;
        }

        val_mean = val_sum / static_cast<TPRECISION>(val_n);
        val_variance = (val_sum2 - static_cast<TPRECISION>(val_n)*val_mean*val_mean) / static_cast<TPRECISION>(numel - (stddev_type==tom::norm::NORM_STDDEV_SAMPLE_1 ? 1 : 0));

        if (stddev_type != tom::norm::NORM_NO_NORM) {
            val_stddev = sqrt(val_variance);
            if (val_mean!=0. || (val_stddev!=1.&&val_stddev!=0.)) {
                // Finally normalise again under the mask.
                if (val_stddev!=1. && val_stddev!=0.) {
                    tom::loop::for_each<tom::Volume<T>, const tom::Volume<TMASK>, ::for_each__tom__norm_under_mask_5<T, TMASK, TPRECISION> >(v, mask, ::for_each__tom__norm_under_mask_5<T, TMASK, TPRECISION>(val_mean, val_stddev));
                } else {
                    tom::loop::for_each<tom::Volume<T>, const tom::Volume<TMASK>, ::for_each__tom__norm_mask__set_mean0<T, TMASK, TPRECISION> >(v, mask, ::for_each__tom__norm_mask__set_mean0<T, TMASK, TPRECISION>(val_mean));
                }
            }
            if (val_stddev) {
                val_variance = 1.;
            } else {
                val_variance = 0.;
            }
        } else {
            if (val_mean!=0.) {
                // Finally normalise again under the mask.
                tom::loop::for_each<tom::Volume<T>, const tom::Volume<TMASK>, ::for_each__tom__norm_mask__set_mean0<T, TMASK, TPRECISION> >(v, mask, ::for_each__tom__norm_mask__set_mean0<T, TMASK, TPRECISION>(val_mean));
            }
        }
        if (variance) {
            *variance = val_variance;
        }
    }
}




/****************************************************************************//**
 * \brief Adds two volumes elementwise.
 *
 * \param[in,out] v1 The left side operand of the elementwise addition and the
 *   destination volume at the same time.
 * \param[in] v2 The right side operand of the addition.
 *
 * Self assignment is not a problem. In that case the volume is added to itself.
 *******************************************************************************/
 
template<typename T1, typename T2>
void tom::element_wise_add(tom::Volume<T1> &v1, const tom::Volume<T2> &v2) {
    using boost::lambda::_1;
    using boost::lambda::_2;
    tom::loop::for_each_no_idx(v1, v2, _1 += _2);
}
/****************************************************************************//**
 * \brief Subtracts two volumes elementwise.
 *
 * \param[in,out] v1 The left side operand of the elementwise subtraction and the
 *   destination volume at the same time.
 * \param[in] v2 The right side operand of the subtraction.
 *
 * Self assignment is not a problem. In that case the volume is subtracted by itself. Yields 0 volume.
 *******************************************************************************/
template<typename T1, typename T2>
void tom::element_wise_sub(tom::Volume<T1> &v1, const tom::Volume<T2> &v2) {
    using boost::lambda::_1;
    using boost::lambda::_2;
    tom::loop::for_each_no_idx(v1, v2, _1 -= _2);
}




namespace {
template<typename T1, typename T2>
struct for_each__tom__element_wise_div {
    for_each__tom__element_wise_div(const T1 &inf_value): _inf_value(inf_value) { }
    const T1 &_inf_value;
    inline void operator()(T1 &v1, const T2 &v2, std::size_t, std::size_t, std::size_t) const {
        if (v2) {
            v1 = v1 / v2;
        } else {
            v1 = _inf_value;
        }
    }
};
template<>
struct for_each__tom__element_wise_div<std::complex<float>, std::complex<float> >{
    for_each__tom__element_wise_div(const std::complex<float> &inf_value): _inf_value(inf_value) { }
    const std::complex<float> &_inf_value;
    inline void operator()(std::complex<float> &v1, const std::complex<float> &v2, std::size_t, std::size_t, std::size_t) const {
        if (v2.real() != 0 || v2.imag() != 0) {
            v1 = v1 / v2;
        } else {
            v1 = _inf_value;
        }
    }
};
template<>
struct for_each__tom__element_wise_div<std::complex<double>, std::complex<double> >{
    for_each__tom__element_wise_div(const std::complex<double> &inf_value): _inf_value(inf_value) { }
    const std::complex<double> &_inf_value;
    inline void operator()(std::complex<double> &v1, const std::complex<double> &v2, std::size_t, std::size_t, std::size_t) const {
        if (v2.real() != 0 || v2.imag() != 0) {
            v1 = v1 / v2;
        } else {
            v1 = _inf_value;
        }
    }
};
}
/****************************************************************************//**
 * \brief Elementwise division
 *
 * \param[in,out] b The left side operand of the elementwise division (numerator)
 *   and the destination volume at the same time.
 * \param[in] a The right side operand of the addition (denominator)
 * \param[in] inf_value Value to assign if the corresponding voxel of @a a is 0.
 *******************************************************************************/
template<typename T1, typename T2>
void tom::element_wise_div(tom::Volume<T1> &v1, const tom::Volume<T2> &v2, T1 inf_value) {
    tom::loop::for_each<tom::Volume<T1>, const tom::Volume<T2>, ::for_each__tom__element_wise_div<T1, T2> >(v1, v2, ::for_each__tom__element_wise_div<T1, T2>(inf_value));
}




/****************************************************************************//**
 * \brief Set lower limit of the volume.
 *
 * \param[in,out] a The volume to be limited.
 * \param[in] threshold The threshold
 * \param[in] value The value.
 *
 * Every value of the volume \a a which is smaller! then threshold is set
 * to value.
 *******************************************************************************/
template<typename T>
void tom::element_wise_set_below_threshold(tom::Volume<T> &a, T threshold, T value) {
    tom::loop::for_each_no_idx(a, boost::lambda::if_then(boost::lambda::_1<threshold,
                                                         boost::lambda::_1 = value));
}






namespace local {
namespace functor {
namespace {
template<typename T>
struct minmax {
    T &min_, &max_;
    minmax(const T &first_val, T &min, T &max): min_(min), max_(max) {
        min_ = max_ = first_val;
    }
    void operator()(const T &val) {
        if (val < min_) {
            min_ = val;
        } else if (val > max_) {
            max_ = val;
        }
    }
}; // struct minmax
template<typename T>
struct min {
    T &min_;
    min(const T &first_val, T &min): min_(min) {
        min_ = first_val;
    }
    void operator()(const T &val) {
        if (val < min_) {
            min_ = val;
        }
    }
}; // struct min
template<typename T>
struct max {
    T &max_;
    max(const T &first_val, T &max): max_(max) {
        max_ = first_val;
    }
    void operator()(const T &val) {
        if (val > max_) {
            max_ = val;
        }
    }
}; // struct max
template<typename T, typename TPRECISION>
struct stat_1 {
    TPRECISION sum_, sum2_;
    stat_1(): sum_(0), sum2_(0) { }
    void operator()(const T &a) {
        sum_  += a;
        sum2_ += a*a;
    }
}; // struct stat_1
template<typename T, typename TPRECISION>
struct stat_2 {
    TPRECISION sum_, sum2_;
    T &min_, &max_;
    stat_2(const T &first_val, T &min, T &max): sum_(0), sum2_(0), min_(min), max_(max) { }
    void operator()(const T &a) {
        if (a < min_) {
            min_ = a;
        } else if (a > max_) {
            max_ = a;
        }
        sum_  += a;
        sum2_ += a*a;
    }
}; // struct stat_2
template<typename T>
struct limit{
	T lowerBound_;
	T lowerReplacement_;
	T upperBound_;
	T upperReplacement_;
	bool doLower_;
	bool doUpper_;

	limit(T &first_val, T lb, T lr, T ub, T ur,bool l, bool u):lowerBound_(lb), lowerReplacement_(lr), upperBound_(ub), upperReplacement_(ur) , doLower_(l) , doUpper_(u) {}

	void operator()(T &a){

		if(doLower_)
			if(a <= lowerBound_)
				a = lowerReplacement_;

		if(doUpper_)
			if(a >= upperBound_)
				a = upperReplacement_;

	}
}; // struct limit

} // namespace
} // namespace functor
} // namespace local
/****************************************************************************//**
 * \brief Find the minimum and the maximum of the volume.
 *******************************************************************************/
template<typename T>
void tom::minmax(const Volume<T> &v, T &minimum, T &maximum) {
    tom::loop::for_each_no_idx(v, local::functor::minmax<T>(v.get(), minimum, maximum));
}
/****************************************************************************//**
 * \brief Find the minimum of the volume.
 *******************************************************************************/
template<typename T>
T tom::min(const Volume<T> &v) {
    T minimum;
    tom::loop::for_each_no_idx(v, local::functor::min<T>(v.get(), minimum) );
    return minimum;
}
/****************************************************************************//**
 * \brief Find the maximum of the volume.
 *******************************************************************************/
template<typename T>
T tom::max(const Volume<T> &v) {
    T maximum;
    tom::loop::for_each_no_idx(v, local::functor::max<T>(v.get(), maximum));
    return maximum;
}
/****************************************************************************//**
 * \brief Replace values within a volume.
 *******************************************************************************/
template<typename T>
void tom::limit(Volume<T> &v,T lowerBound, T lowerReplacement, T upperBound, T upperReplacement, bool doLower, bool doUpper){
	tom::loop::for_each_no_idx(v,local::functor::limit<T>(v.get(),lowerBound,lowerReplacement,upperBound,upperReplacement, doLower, doUpper));
}

/****************************************************************************//**
 * \brief Compute the mean and the variance of a volume.
 *******************************************************************************/
template<typename T, typename TPRECISION>
void tom::stat(const Volume<T> &v, TPRECISION &mean, TPRECISION &variance, bool use_sample_standard_deviation) {
    BOOST_STATIC_ASSERT((boost::is_floating_point<TPRECISION>::value));
    local::functor::stat_1<T, TPRECISION> f;
    tom::loop::for_each_no_idx<const tom::Volume<T>, local::functor::stat_1<T, TPRECISION> &>(v, f);

    mean = f.sum_ / static_cast<TPRECISION>(v.numel());
    variance = (f.sum2_ - static_cast<TPRECISION>(v.numel())*mean*mean) / static_cast<TPRECISION>(v.numel() - (use_sample_standard_deviation ? 1 : 0));
}
/****************************************************************************//**
 * \brief Return the variance, the mean and the min, max of the volume.
 *******************************************************************************/
template<typename T, typename TPRECISION>
void tom::stat(const Volume<T> &v, TPRECISION &mean, TPRECISION &variance, bool use_sample_standard_deviation, T &min, T &max) {
    BOOST_STATIC_ASSERT((boost::is_floating_point<TPRECISION>::value));
    local::functor::stat_2<T, TPRECISION> f(v.get(), min, max);
    tom::loop::for_each_no_idx<const tom::Volume<T>, local::functor::stat_2<T, TPRECISION> &>(v, f);
    mean = f.sum_ / static_cast<TPRECISION>(v.numel());
    variance = (f.sum2_ - static_cast<TPRECISION>(v.numel())*mean*mean) /
                (static_cast<TPRECISION>(v.numel()) - static_cast<TPRECISION>(use_sample_standard_deviation ? 1 : 0));
}
/****************************************************************************//**
 * \brief Calculate the squared sum of a volume (2-norm).
 *******************************************************************************/
template<typename T, typename TPRECISION>
TPRECISION tom::sum_squared(const tom::Volume<T> &v) {
    using boost::lambda::_1;
    using boost::lambda::ll_static_cast;
    TPRECISION sum2 = 0;
    tom::loop::for_each_no_idx(v, sum2 += ll_static_cast<TPRECISION>(_1)*ll_static_cast<TPRECISION>(_1));
    return sum2;
}
/****************************************************************************//**
 * \brief Calculate the sum over a volume
 * \param[in] v Volume.
 *******************************************************************************/
template<typename T, typename TPRECISION>
TPRECISION tom::sum(const Volume<T> &v) {
    using boost::lambda::_1;
    TPRECISION sum = 0;
    tom::loop::for_each_no_idx(v, sum += _1);
    return sum;
}





namespace local {
namespace functor {
namespace {
template<typename T1, typename T2> inline void element_wise_multiply_f(T1 &v1, const T2 &v2) { v1 = v1 * v2; }
#if __xlC__ == 0x0800
// in the complex header from the currently used xlC compiler the methods real and imag to not return a reference...
// strange :)
template<> inline void element_wise_multiply_f<std::complex<float >, double>(std::complex<float > &v1, const double &v2) { v1.real(v1.real()*v2); v1.imag(v1.imag()*v2); }
template<> inline void element_wise_multiply_f<std::complex<double>, float >(std::complex<double> &v1, const float  &v2) { v1.real(v1.real()*v2); v1.imag(v1.imag()*v2); }
#else
template<> inline void element_wise_multiply_f<std::complex<float >, double>(std::complex<float > &v1, const double &v2) { v1.real() = v1.real()*v2; v1.imag() = v1.imag()*v2; }
template<> inline void element_wise_multiply_f<std::complex<double>, float >(std::complex<double> &v1, const float  &v2) { v1.real() = v1.real()*v2; v1.imag() = v1.imag()*v2; }
#endif
template<typename T1, typename T2>
struct element_wise_multiply {
    void operator()(T1 &v1, const T2 &v2) { element_wise_multiply_f<T1, T2>(v1, v2); }
};
} } } // namespace local::functor::<unnamed>
/****************************************************************************//**
 * \brief Elementwise multiplication
 *
 * \param[in,out] v1 The left side operand of the elementwise multiplication
 *   and the destination volume at the same time.
 * \param[in] v2 The right side operand of the multiplication.
 *******************************************************************************/
template<typename T1, typename T2>
void tom::element_wise_multiply(tom::Volume<T1> &v1, const tom::Volume<T2> &v2) {
    tom::loop::for_each_no_idx(v1, v2, local::functor::element_wise_multiply<T1, T2>());
}







namespace {
template<typename T>
struct for_each__tom__conjugate {
    inline void operator()(T &v1,std::size_t, std::size_t, std::size_t) { v1 = conj(v1); }
};
}
/****************************************************************************//**
 * \brief Conjugate complex input (inplace)
 *
 * \param[in,out] v The complex volume which will be conjugated.
 *******************************************************************************/
template<typename T>
void tom::conjugate(tom::Volume<T> &v) {
    tom::loop::for_each<tom::Volume<T>, ::for_each__tom__conjugate<T> >(v,::for_each__tom__conjugate<T>());
}





namespace tom {
namespace functor {
namespace {
template<typename T1, typename T2, typename T3>
struct element_wise_conj_multiply_copy_ {
    void operator()(T1 &v1, const T2 &v2, const T3 &v3, std::size_t, std::size_t, std::size_t) const {
        v1.DO_NOT_INSTANTIATE_THIS_TYPE;
        v1 = v2 * v3;
    }
};

template<typename T1, typename T2, typename T3>
struct element_wise_conj_multiply_copy_<T1, std::complex<T2>, T3> {
    void operator()(T1 &v1, const std::complex<T2> &v2, const T3 &v3, std::size_t, std::size_t, std::size_t) const {
        v1 = std::conj(v2) * v3;
    }
};
} // namespace
} // namespace functor
}
/***********************************************************************//**
 * \file volume_fcn.cpp
 * \brief Implementation of functions manipulating the volume class.
 * \author  Thomas Haller
 * \version 0.2
 * \date    26.11.2008
 **************************************************************************/

namespace tom {
namespace functor {
namespace {
template<typename T1, typename T2>
struct element_wise_conj_multiply_inplace_ {
    void operator()(T1 &v1, const T2 &v2, std::size_t, std::size_t, std::size_t) const {
        v1.DO_NOT_INSTANTIATE_THIS_TYPE;
        v1 = v1 * v2;
    }
};
template<typename T1, typename T2>
struct element_wise_conj_multiply_inplace_<std::complex<T1>, T2> {
    void operator()(std::complex<T1> &v1, const T2 &v2, std::size_t, std::size_t, std::size_t) const {
        v1 = std::conj(v1) * v2;
    }
};
} // namespace
} // namespace functor
} // namespace tom
/****************************************************************************//**
 * \brief Elementwise multiplication
 *
 * \param[in,out] v1 The left side operand of the elementwise multiplication
 *   and the destination volume at the same time. The conjugate complex
 *   value of b is taken.
 * \param[in] v2 The right side operand of the multiplication.
 *******************************************************************************/
template<typename T1, typename T2>
void tom::element_wise_conj_multiply(tom::Volume<T1> &v1, const tom::Volume<T2> &v2) {
    tom::loop::for_each(v1, v2, tom::functor::element_wise_conj_multiply_inplace_<T1, T2>());
}






// namespace tom
/****************************************************************************//**
 * \brief Elementwise multiplication of two volumes.
 *
 * \param[out] v1 Resulting volume which is filled with conj(v2).*v3
 * \param[in] v2 Operand of the multiplication. Its conjugate complex is taken.
 * \param[in] v3 Second operand of the multiplication.
 *
 * All volumes must have the same size.
 *******************************************************************************/
template<typename T1, typename T2, typename T3>
void tom::element_wise_conj_multiply(tom::Volume<T1> &v1, const tom::Volume<T2> &v2, const tom::Volume<T3> &v3) {
    tom::loop::for_each(v1, v2, v3, tom::functor::element_wise_conj_multiply_copy_<T1, T2, T3>());
}




/****************************************************************************//**
 * \brief Shift zero-frequency component to center of spectrum.
 *
 * \param[in] vsrc The source volume.
 * \param[out] v The shifted volume.
 * \param[in] is_ifftshift If true it behaves as \c ifftshift from MATLAB.
 *   Otherwise as \c fftshift. In case of even volume sizes, there is not
 *   difference. Otherwise ifftshift is the inverse operation of the fftshift.
 *
 * See the documentation of fftshift/ifftshift from MATLAB :)
 * \TODO The current implementation could be done better.
 *******************************************************************************/
template<typename T>
void tom::fftshift(const Volume<T> &vsrc, Volume<T> &v, bool is_ifftshift) {

    if (vsrc.getSizeX() != v.getSizeX() || vsrc.getSizeY() != v.getSizeY() || vsrc.getSizeZ() != v.getSizeZ()) {
        throw std::runtime_error("The volumes must have the same size.");
    }

    const std::size_t sizex = v.getSizeX();
    const std::size_t sizey = v.getSizeY();
    const std::size_t sizez = v.getSizeZ();
    const std::size_t sizex2 = sizex/2 + (is_ifftshift ? 0 : sizex%2);
    const std::size_t sizey2 = sizey/2 + (is_ifftshift ? 0 : sizey%2);
    const std::size_t sizez2 = sizez/2 + (is_ifftshift ? 0 : sizez%2);

    std::size_t x, y, z;
    std::size_t z_tmp, y_tmp, x_tmp;
    if (v.isContiguous() && vsrc.isContiguous()) {
        const T *pvsrc = &vsrc.get();
        T *pv = &v.get();
        std::size_t i;
        i = 0;
        for (z=0; z<sizez; z++) {
            z_tmp = ((z+sizez2)%sizez) * sizey;
            for (y=0; y<sizey; y++) {
                y_tmp = (z_tmp + (y+sizey2)%sizey) * sizex;
                for (x=0; x<sizex; x++) {
                    x_tmp = y_tmp + (x+sizex2)%sizex;
                    pv[i++] = pvsrc[x_tmp];
                }
            }
        }
    } else {
        for (z=0; z<sizez; z++) {
            z_tmp = (z+sizez2)%sizez;
            for (y=0; y<sizey; y++) {
                y_tmp = (y+sizey2)%sizey;
                for (x=0; x<sizex; x++) {
                    v.get(x,y,z) = vsrc.get((x+sizex2)%sizex, y_tmp, z_tmp);
                }
            }
        }
    }
}




namespace {
template<typename T>
struct for_each__tom__number_set_voxels {
	for_each__tom__number_set_voxels() {numel=0;}
	std::size_t numel;
    inline void operator()(const T &a, std::size_t, std::size_t, std::size_t) {
        if (a != 0) {
        	numel++;
        }
    }
};
}
/****************************************************************************//**
 * \brief Determines the number of voxels != 0 in volume
 *******************************************************************************/
template<typename T>
std::size_t tom::number_set_voxels(const tom::Volume<T> &v) {
    using boost::lambda::_1;
    using boost::lambda::if_then;
    using boost::lambda::var;
    std::size_t n = 0;
    tom::loop::for_each_no_idx(v, if_then(_1!=0, var(n)++));

    #ifndef NDEBUG
    ::for_each__tom__number_set_voxels<T> f;
    tom::loop::for_each<const tom::Volume<T>, ::for_each__tom__number_set_voxels<T> &>(v, f);
    assert(f.numel == n);
    #endif

    return n;
}




namespace lc {
namespace functor {
namespace {
template<typename T1, typename T2>
struct peakv {
    peakv(T1 max): max_(max) { }
    T1 max_;
    std::vector<tom::st_idx> res_;
    void operator()(const T1 &v1, const T2 &v2, const tom::loop::index_type &x, const tom::loop::index_type &y, const tom::loop::index_type &z) {
        if (v2 && !(v1<max_)) {
            tom::st_idx idx;
            idx.x = x;
            idx.y = y;
            idx.z = z;
            if (v1 > max_) {
                max_ = v1;
                res_.clear();
            }
            res_.push_back(idx);
        }
    }
};
} } } // namespace lc::functor::<unnamed>
/****************************************************************************//**
 * Gets the index of the maximum value.
 *******************************************************************************/
template<typename T1, typename T2>
std::vector<tom::st_idx> tom::peak(const tom::Volume<T1> &v, const tom::Volume<T2> &mask) {
    BOOST_STATIC_ASSERT((boost::is_floating_point<T1>::value || boost::is_integral<T1>::value));
    T1 s0;
    if (boost::is_floating_point<T1>::value) {
        s0 = - std::numeric_limits<T1>::max();
    } else if (boost::is_integral<T1>::value) {
        s0 = std::numeric_limits<T1>::min();
    } else {
    }
    lc::functor::peakv<T1, T2> s(s0);
    tom::loop::for_each<const tom::Volume<T1>, const tom::Volume<T2>, lc::functor::peakv<T1, T2> &>(v, mask, s);
    return s.res_;
}

/****************************************************************************//**
 * \brief Expands the reduced complex volume to its full dimension.
 * \param[in] vsrc The complex volume. It can be either reduced
 *   (with dimension z only of size z/2+1) or full. However only
 *   the upper half is considered.
 * \param[out] v The complex volume. It has the dimension of the full
 *   (output) volume. The upper half (0:sizex-1, 0:sizey-1, 0:sizez/2+1)
 *   contains the complex values as returned by a R2C fftw.
 *   The second half (0:sizex-1, 0:sizey-1, sizez/2+2:sizez) is filled with
 *   the conjugate complex values to fullfill the hermitian symmetry.
 *******************************************************************************/
template<typename T>
void tom::hermitian_symmetry_to_full(const tom::Volume<T> &vsrc, tom::Volume<T> &v) {
    if (vsrc.getSizeZ() != v.getSizeZ() && vsrc.getSizeZ() != v.getSizeZ()/2+1) {
        throw std::invalid_argument("The input volume must be either reduced or full complex.");
    }
    std::size_t sizez = v.getSizeZ()/2+1;
    if (&vsrc != &v) {
    	//pointer to the first element in data source
        T *data = const_cast<T *>(&vsrc.get());

        tom::Volume<T>(v, NULL, v.getSizeX(), v.getSizeY(), sizez, v.getStrideX(), v.getStrideY(), v.getStrideZ()).setValues(
            tom::Volume<T>(data, vsrc.getSizeX(), vsrc.getSizeY(), sizez, vsrc.getStrideX(), vsrc.getStrideY(), vsrc.getStrideZ(), false, NULL));
    }
    tom::hermitian_symmetry_to_full(v);
}









namespace {
template<typename T>
inline T local_fcn_conjugate_complex(const T &a) {
    return conj(a);
}
template<> inline float  local_fcn_conjugate_complex(const float  &a) { return a; }
template<> inline double local_fcn_conjugate_complex(const double &a) { return a; }
}
/****************************************************************************//**
 * \brief Expands the reduced volume to its full dimension.
 *
 * \param[in,out] v The volume. It has the dimension of the full
 *   (output) volume. The upper half (0:sizex-1, 0:sizey-1, 0:sizez/2+1)
 *   contains the already initialized values.
 *   This function copies that half into the second half accoding to
 *   the properties of the hermitian symetrie.\n
 *   For example the fftw real to complex transformations return a the half
 *   volume of complex numbers.
 *
 * Can also be a pure real volume (T==float or double)
 *******************************************************************************/
template<typename T>
void tom::hermitian_symmetry_to_full(tom::Volume<T> &v) {
    typedef T local_TYPE;
    local_TYPE *local_A = &v.get();
    local_TYPE *v_data = &v.get();
    std::size_t local_sizex = v.getSizeX();
    std::size_t local_sizey = v.getSizeY();
    std::size_t local_sizez = v.getSizeZ();
    std::size_t local_stridex = v.getStrideX();
    std::size_t local_stridey = v.getStrideY();
    std::size_t local_stridez = v.getStrideZ();

    {
        std::size_t x, y, z;
        std::size_t x2, y2, z2;
        std::size_t y2_cnt, x2_cnt;

        tom::loop::ptr_add_byte_offset(local_A, (local_sizez/2+1)*local_stridez);
        if (!(local_stridex%sizeof(local_TYPE)) &&
            !(local_stridey%sizeof(local_TYPE)) &&
            !(local_stridez%sizeof(local_TYPE))) {

            local_stridex /= sizeof(local_TYPE);
            local_stridey /= sizeof(local_TYPE);
            local_stridez /= sizeof(local_TYPE);

            if (local_stridex == 1) {
                if (local_stridey == local_sizex &&
                    local_stridez == local_sizex*local_sizey) {
                    std::size_t i = 0;
                    for (z=local_sizez/2+1, z2=(local_sizez+1)/2-1; z<local_sizez; z++,z2--) {
                        for (y=0, y2_cnt=local_sizey; y<local_sizey; y++, y2_cnt--) {
                            y2 = y2_cnt % local_sizey;
                            for (x=0, x2_cnt=local_sizex; x<local_sizex; x++, i++, x2_cnt--) {
                                x2 = x2_cnt % local_sizex;
                                local_A[i] = local_fcn_conjugate_complex<local_TYPE>(v_data[z2*local_stridez + y2*local_stridey + x2]);
                            }
                        }
                    }
                } else {
                    local_stridez -= local_sizey*local_stridey;
                    for (z=local_sizez/2+1, z2=(local_sizez+1)/2-1; z<local_sizez; z++,z2--) {
                        for (y=0, y2_cnt=local_sizey; y<local_sizey; y++, y2_cnt--) {
                            y2 = y2_cnt % local_sizey;
                            for (x=0, x2_cnt=local_sizex; x<local_sizex; x++, x2_cnt--) {
                                x2 = x2_cnt % local_sizex;
                                local_A[x] = local_fcn_conjugate_complex<local_TYPE>(v_data[z2*local_stridez + y2*local_stridey + x2]);
                            }
                            local_A += local_stridey;
                        }
                        local_A += local_stridez;
                    }
                }
            } else {
                local_stridez -= local_sizey*local_stridey;
                local_stridey -= local_sizex*local_stridex;
                for (z=local_sizez/2+1, z2=(local_sizez+1)/2-1; z<local_sizez; z++,z2--) {
                    for (y=0, y2_cnt=local_sizey; y<local_sizey; y++, y2_cnt--) {
                        y2 = y2_cnt % local_sizey;
                        for (x=0, x2_cnt=local_sizex; x<local_sizex; x++, x2_cnt--) {
                            x2 = x2_cnt % local_sizex;
                            local_A[0] = local_fcn_conjugate_complex<local_TYPE>(v_data[z2*local_stridez + y2*local_stridey + x2*local_stridex]);
                            local_A += local_stridex;
                        }
                        local_A += local_stridey;
                    }
                    local_A += local_stridez;
                }
            }
        } else {
            local_stridez -= local_sizey*local_stridey;
            local_stridey -= local_sizex*local_stridex;
            for (z=local_sizez/2+1, z2=(local_sizez+1)/2-1; z<local_sizez; z++,z2--) {
                for (y=0, y2_cnt=local_sizey; y<local_sizey; y++, y2_cnt--) {
                    y2 = y2_cnt % local_sizey;
                    for (x=0, x2_cnt=local_sizex; x<local_sizex; x++, x2_cnt--) {
                        x2 = x2_cnt % local_sizex;
                        *local_A = local_fcn_conjugate_complex<local_TYPE>(* reinterpret_cast<local_TYPE *>((reinterpret_cast<char *>(v_data)) + z2*local_stridez + y2*local_stridey + x2*local_stridex));
                        tom::loop::ptr_add_byte_offset(local_A, local_stridex);
                    }
                    tom::loop::ptr_add_byte_offset(local_A, local_stridey);
                }
                tom::loop::ptr_add_byte_offset(local_A, local_stridez);
            }
        }
    }
}







namespace {
template<typename T>
struct for_each__tom__make_binary {
    void operator()(T &a, std::size_t, std::size_t, std::size_t) const {
        a = (a==0) ? 0 : 1;
    }
};
}
/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
void tom::make_binary(tom::Volume<T> &v) {
    tom::loop::for_each(v, ::for_each__tom__make_binary<T>());
}










namespace tom {
namespace functor {
namespace {
template<typename T>
struct element_wise_power_pow_ {
    element_wise_power_pow_(T exponent): exponent_(exponent) { }
    T exponent_;
    void operator()(T &a, std::size_t, std::size_t, std::size_t) const {
		if(a < 0 && exponent_ > -1 && exponent_ < 1)
			throw std::runtime_error("Encountered negative value in volume with exponent between -1 and 1! ABORT!");
		
        a = tom::math::pow(a, exponent_);
    }
};
} // namespace
} // namespace functor
} // namespace tom
/****************************************************************************//**
 * \brief Calculates the power of a volume to the exponent
 *******************************************************************************/
template<typename T>
void tom::element_wise_power(tom::Volume<T> &v, T exponent) {
	/*
    if (exponent == 0.5) {
        tom::element_wise_operation(v, &tom::math::sqrt<T>);
    } else {
        tom::loop::for_each(v, tom::functor::element_wise_power_pow_<T>(exponent));
    }
	*/
	
	tom::loop::for_each(v, tom::functor::element_wise_power_pow_<T>(exponent));
}


namespace tom {
namespace functor {
namespace {
template<typename T,typename T2>
struct element_wise_power_pow_2_ {
    element_wise_power_pow_2_(T2 exponent): exponent_(exponent) { }
    T2 exponent_;
    void operator()(T &a, std::size_t, std::size_t, std::size_t) const {
        a = pow(a, exponent_);
    }
};
} // namespace
} // namespace functor
} // namespace tom
/****************************************************************************//**
 * \brief Calculates the power of a volume to the exponent
 *******************************************************************************/
template<typename T,typename T2>
void tom::element_wise_power(tom::Volume<T> &v, T2 exponent) {
	tom::loop::for_each(v, tom::functor::element_wise_power_pow_2_<T,T2>(exponent));
}




namespace {
template<typename T>
struct for_each__tom__element_wise_max {
	for_each__tom__element_wise_max(T val): val(val) { }
	const T val;
	inline void operator()(T &a, std::size_t, std::size_t, std::size_t) { if (a < this->val) { a = this->val; } }
};
}
/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
void tom::element_wise_max(tom::Volume<T> &v, T val) {
	tom::loop::for_each(v, ::for_each__tom__element_wise_max<T>(val));
}








namespace {
template<typename T>
struct for_each__tom__peak {
    for_each__tom__peak(T val0): res(), max(val0) { }
    for_each__tom__peak(const for_each__tom__peak& v): res(), max(0) {
        assert(!"I SHOULD NOT BE CALLED!!");
    }
    std::vector<tom::st_idx> res;
    T max;
    inline void operator()(const T &a, std::size_t x, std::size_t y, std::size_t z) {
        if (!(a < this->max)) {
            tom::st_idx idx;
            idx.x = x;
            idx.y = y;
            idx.z = z;
            if (a > this->max) {
                this->max = a;
                res.clear();
            }
            this->res.push_back(idx);
        }
    }
};
}
/****************************************************************************//**
 * \brief Find the peak of a volume.
 *******************************************************************************/
template<typename T>
std::vector<tom::st_idx> tom::peak(const tom::Volume<T> &v) {
    ::for_each__tom__peak<T> p(v.get());
    tom::loop::for_each<const tom::Volume<T>, ::for_each__tom__peak<T> &>(v, p);
    return p.res;
}



namespace tom {
namespace functor {
namespace {
template<typename T>
struct isfinite {
    bool isfinite_;
    isfinite(): isfinite_(true) { }
    bool operator()(const T &a, std::size_t, std::size_t, std::size_t);
};
template<> inline bool isfinite<float                >::operator()(const float                &a, std::size_t, std::size_t, std::size_t) { return (isfinite_ = std::isfinite(a)                                  ); }
template<> inline bool isfinite<double               >::operator()(const double               &a, std::size_t, std::size_t, std::size_t) { return (isfinite_ = std::isfinite(a)                                  ); }
template<> inline bool isfinite<std::complex<float > >::operator()(const std::complex<float > &a, std::size_t, std::size_t, std::size_t) { return (isfinite_ = std::isfinite(a.real()) && std::isfinite(a.imag())); }
template<> inline bool isfinite<std::complex<double> >::operator()(const std::complex<double> &a, std::size_t, std::size_t, std::size_t) { return (isfinite_ = std::isfinite(a.real()) && std::isfinite(a.imag())); }
} // namespace
} // namespace functor
} // namespace tom
/****************************************************************************//**
 * \brief Determines whether there are non-finite floating point values in the volume.
 *
 * Operating with NaN or Inf floating point volumes will lead to failure.
 * Most (all?) function from libtomc don't check for these exception for performance
 * and assume that all values are finite. Use this function to check for non-finite
 * values. It calles \c isfinite from \c math.h.
 *******************************************************************************/
template<typename T>
bool tom::isfinite(const tom::Volume<T> &v) {
    tom::functor::isfinite<T> op;
    tom::loop::for_each_while<const tom::Volume<T>, tom::functor::isfinite<T> &>(v, op);
    return op.isfinite_;
}




/****************************************************************************//**
 * \brief Computes the variance and the mean of a volume under a given mask.
 *
 * \param[in] v Volume of which the mean and the variance is computed.
 * \param[out] mean Here the mean under the mask will be returned.
 * \param[out] variance Here the variance under the mask will be returned.
 * \param[in]  use_sample_standard_deviation If true, the variance is computed
 *    using the standard deviation of the sample (with nominator N-1).
 * \param[in] mask Only those voxels in the source volume are considered where
 *    the mask has values != 0.
 *    However the mean and the variance are normalised relative to the whole image!!
 *******************************************************************************/
template<typename T, typename TMASK, typename TPRECISION>
void tom::stat(const Volume<T> &v, TPRECISION &mean, TPRECISION &variance, bool use_sample_standard_deviation, const Volume<TMASK> &mask) {
    using boost::lambda::if_then;
    using boost::lambda::_1;
    using boost::lambda::_2;
    using boost::lambda::ll_static_cast;
    TPRECISION sum  = 0;
    TPRECISION sum2 = 0;
    tom::loop::for_each_no_idx(v, mask, if_then(_2!=0, (sum+=_1, sum2+=ll_static_cast<TPRECISION>(_1)*ll_static_cast<TPRECISION>(_1)) ) );

    mean = sum / static_cast<TPRECISION>(v.numel());
    variance = (sum2 - static_cast<TPRECISION>(v.numel())*mean*mean) / static_cast<TPRECISION>(v.numel()-(use_sample_standard_deviation ? 1 : 0));
}




namespace {
template<typename TPRECISION, typename T>
struct for_each_step__tom__init_sphere__sigma {
    for_each_step__tom__init_sphere__sigma(TPRECISION centerx, TPRECISION centery, TPRECISION centerz,
                                         TPRECISION radius, TPRECISION sigma, TPRECISION max_radius)
                                         : centerx(centerx), centery(centery), centerz(centerz),
                                           radius(radius), sigma(sigma), max_radius(max_radius),
                                           current_z(0), current_yz(0) { }
    TPRECISION centerx, centery, centerz, radius, sigma, max_radius;
    TPRECISION current_z, current_yz;
    inline void stepz(TPRECISION z) { current_z  = z - centerz; current_z  *= current_z; }
    inline void stepy(TPRECISION y) { current_yz = y - centery; current_yz  = current_yz*current_yz + current_z; }
    inline void operator()(T &a, TPRECISION x) {
        x -= centerx;
        x = tom::math::sqrt<TPRECISION>(x*x + current_yz);
        if (x <= this->radius) {
            a = 1;
        } else if (x > this->max_radius) {
            a = 0;
        } else {
            x = ((x - radius)/sigma);
            a = tom::math::exp<TPRECISION>(- x*x);
        }
    }
};
template<typename TPRECISION>
struct for_each_step__tom__init_sphere__sigma<TPRECISION, char> {
    for_each_step__tom__init_sphere__sigma(TPRECISION centerx, TPRECISION centery, TPRECISION centerz,
                                         TPRECISION radius, TPRECISION sigma, TPRECISION max_radius_input)
                                         : centerx(centerx), centery(centery), centerz(centerz),
                                           radius(radius), sigma(sigma), max_radius(max_radius_input*max_radius_input),
                                           current_z(0), current_yz(0) { }
    TPRECISION centerx, centery, centerz, radius, sigma, max_radius;
    TPRECISION current_z, current_yz;
    inline void stepz(TPRECISION z) {
        current_z = z - centerz;
        current_z *= current_z;
    }
    inline void stepy(TPRECISION y) {
        current_yz = y - centery;
        current_yz = current_yz*current_yz + current_z;
    }
    inline void operator()(char &a, TPRECISION x) {
        x -= centerx;
        // Sigma with an integer type is only the same as using max_radius.
        a = x*x + current_yz <= max_radius;
    }
};
template<typename TPRECISION, typename T>
struct for_each_step__tom__init_sphere {
    for_each_step__tom__init_sphere(TPRECISION centerx, TPRECISION centery, TPRECISION centerz, TPRECISION radius)
                                         : centerx_(centerx), centery_(centery), centerz_(centerz), radius_(radius*radius), z_(0), y_(0) { }
    TPRECISION centerx_, centery_, centerz_, radius_;
    TPRECISION z_, y_;
    inline void stepz(TPRECISION z) {
        z_ = z - centerz_;
        z_ = z_*z_; }
    inline void stepy(TPRECISION y) {
        y_ = y - centery_;
        y_ = y_*y_ + z_; }
    inline void operator()(T &a, TPRECISION x) {
        x -= centerx_;
        const TPRECISION dist = x*x + y_;
        a = (dist <= radius_);
        //a = (tom::math::sqrt<TPRECISION>(x*x + this->current_y) <= this->radius);
    }
};
}
/****************************************************************************//**
 * \brief initialises volume with sphere mask.
 * \param[out] mask The mask which will be initialized.
 * \param[in] radius The radius of the sphere. Inside this radius the mask has
 *   value 0.
 * \param[in] sigma If > 0 there is a gaussian descent of the values with deviation
 *   sigma.
 * \param[in] max_radius Ignored if sigma <= 0. Otherwise it specifies, that outside
 *   of max_radius all elements are set to 0. Setting max_radius smaller than
 *   radius defaults to max_radius = radius+sqrt(2*sigma).
 * \param[in] centerx
 * \param[in] centery
 * \param[in] centerz
 *
 * Inside of the radius, the values are set to 1(!!)
 *******************************************************************************/
template<typename T>
void tom::init_sphere(tom::Volume<T> &mask, float radius, float sigma, float max_radius, float centerx, float centery, float centerz) {
    typedef float TPRECISION;
    if (radius < 0.) {
        radius = 0.;
    }
    if (sigma > 0) {
        if (max_radius < radius) {
            max_radius = radius + sqrt(2.*sigma);
        }
        tom::loop::for_each_step(mask, ::for_each_step__tom__init_sphere__sigma<TPRECISION, T>(centerx, centery, centerz, radius, sigma, max_radius));
    } else if (radius > 0.) {
        tom::loop::for_each_step(mask, ::for_each_step__tom__init_sphere       <TPRECISION, T>(centerx, centery, centerz, radius                   ));
    } else {
        mask.setValues(0);
    }
}






namespace tom {
namespace functor {
namespace {
template<typename T, typename TPRECISION>
struct init_gaussian_sphere {
    init_gaussian_sphere(TPRECISION _centerx, TPRECISION _centery, TPRECISION _centerz, TPRECISION _sigma)
        :   centerx(_centerx), centery(_centery), centerz(_centerz),
            sigma_square_twice_minus(-2. * _sigma*_sigma), current_z(0), current_yz(0), sum_(0) {
    }
    TPRECISION centerx, centery, centerz, sigma_square_twice_minus, factor;
    TPRECISION current_z, current_yz;
    double sum_;
    void stepz(TPRECISION z) { current_z  = z - centerz; current_z  *= current_z; }
    void stepy(TPRECISION y) { current_yz = y - centery; current_yz  = current_yz*current_yz + current_z; }
    void operator()(T &a, TPRECISION x) {
        x -= centerx;
        x = x*x + current_yz;
        x = tom::math::exp<TPRECISION>(x / sigma_square_twice_minus);
        sum_ += x;
        a     = x;
    }
};
} // namespace
} // namespace functor
} // namespace tom
/****************************************************************************//**
 * \brief initialises volume with a Gausian sphere mask.
 * \param[out] mask The mask which will be initialized.
 * \param[in] radius The radius of the sphere. Inside this radius the mask has
 *   value 0.
 * \param[in] sigma If > 0 there is a gaussian descent of the values with deviation
 *   sigma.
 * \param[in] max_radius Ignored if sigma <= 0. Otherwise it specifies, that outside
 *   of max_radius all elements are set to 0. Setting max_radius smaller than
 *   radius defaults to max_radius = radius+sqrt(2*sigma).
 * \param[in] centerx
 * \param[in] centery
 * \param[in] centerz
 *
 * The sum under the sphere is set to 1.
 *******************************************************************************/
template<typename T>
void tom::init_gaussian_sphere(tom::Volume<T> &v, double sigma, double centerx, double centery, double centerz) {
    typedef double TPRECISION;
    if (sigma < 0) {
        std::invalid_argument("tom::init_gaussian_sphere: sigma can not be negative.");
    } else if (sigma == 0) {
        v.setValues(0.);
        const std::size_t centerx_i = static_cast<std::size_t>(centerx);
        const std::size_t centery_i = static_cast<std::size_t>(centery);
        const std::size_t centerz_i = static_cast<std::size_t>(centerz);
        if (centerx_i==centerx && centery_i==centery && centerz_i==centerz &&
            centerx_i<v.getSizeX() && centery_i<v.getSizeY() && centerz_i<v.getSizeZ()) {
            v.get(centerx_i, centery_i, centerz_i) = 1.;
        } else {
            std::runtime_error("tom::init_gaussian_sphere: Cannot set Dirac pulse outside grid point of the discrete volume.");
        }
    } else {
        tom::functor::init_gaussian_sphere<T, TPRECISION> s(centerx, centery, centerz, sigma);
        tom::loop::for_each_step<tom::Volume<T>, tom::functor::init_gaussian_sphere<T, TPRECISION> &>(v, s);
        v.template shift_scale<TPRECISION>(0, 1. / s.sum_);
    }
}




namespace tom {
namespace functor {
namespace {
template<typename T>
struct fourier_shell_correlation_ {
    typedef double TPRECISION;
    typedef float  TPRECISION_RING;

    fourier_shell_correlation_(std::size_t n_shells, std::size_t size):
        n_shells(n_shells),
        ampl1(NULL),
        ampl2(NULL),
        ampl_diff(NULL),
        phares1(NULL),
        phares2(NULL),
        n_shell(NULL),
        center((size+1)/2),
        sy(0),
        sz(0) {
        assert(size>0 && n_shells>0);
        try {
            ampl1       = new double[n_shells];
            ampl2       = new double[n_shells];
            ampl_diff   = new double[n_shells];
            phares1     = new double[n_shells];
            phares2     = new double[n_shells];
            n_shell     = new std::size_t[n_shells];
        } catch (...) {
            delete[] ampl1;
            delete[] ampl2;
            delete[] ampl_diff;
            delete[] phares1;
            delete[] phares2;
            delete[] n_shell;
            throw;
        }
        for (std::size_t i=0; i<n_shells; i++) {
            ampl1[i]        = 0;
            ampl2[i]        = 0;
            ampl_diff[i]    = 0;
            phares1[i]      = 0;
            phares2[i]      = 0;
            n_shell[i]      = 0;
        }
        shell_thickness = static_cast<TPRECISION_RING>(size-center) / static_cast<TPRECISION_RING>(n_shells);
        size3 = size*size*size;
    }
    ~fourier_shell_correlation_() {
        delete[] ampl1;
        delete[] ampl2;
        delete[] ampl_diff;
        delete[] phares1;
        delete[] phares2;
        delete[] n_shell;
    }
    inline void stepz(signed long z) {
        sz = z - center;
        sz = sz * sz;
    }
    inline void stepy(signed long y) {
        sy = y - center;
        sy = sz + sy*sy;
    }
    inline void operator()(const std::complex<T> &f1, const std::complex<T> &f2, signed long x) {
        x -= center;
        if ((x = sy + x*x)) {
            const TPRECISION_RING shelld = tom::math::sqrt<TPRECISION_RING>(x) / shell_thickness;
            if (shelld < n_shells) {
                const std::size_t shell = static_cast<std::size_t>(shelld);
                assert(shelld>=0 && shell<n_shells);

                n_shell[shell]++;

                const TPRECISION f1_real = f1.real();
                const TPRECISION f1_imag = f1.imag();
                const TPRECISION f2_real = f2.real();
                const TPRECISION f2_imag = f2.imag();
                const TPRECISION amplitude1 = (f1_real*f1_real + f1_imag*f1_imag) / size3;
                const TPRECISION amplitude2 = (f2_real*f2_real + f2_imag*f2_imag) / size3;
                const TPRECISION f12_real = f1_real - f2_real;
                const TPRECISION f12_imag = f1_imag - f2_imag;
                const TPRECISION amplitude_diff = (f12_real*f12_real + f12_imag*f12_imag) / size3;
                const TPRECISION phares_v = tom::math::sqrt<TPRECISION>(amplitude1) + tom::math::sqrt<TPRECISION>(amplitude2); // Here there is a difference from tom_comparec.c

                ampl1[shell] += amplitude1;
                ampl2[shell] += amplitude2;
                ampl_diff[shell] += amplitude_diff;
                phares1[shell] += phares_v;

                TPRECISION arg = 2*tom::math::sqrt<TPRECISION>(amplitude1*amplitude2); /* calculation of cos(theta) for Phase Residuum */
                if (arg > 0) {
                    arg = (amplitude1 + amplitude2 - amplitude_diff)/arg;
                    if (arg>1) {
                        arg = 1;
                    } else if (arg<-1) {
                        arg = -1;
                    }
                }
                const TPRECISION delta = tom::math::acos<TPRECISION>(arg);    /* Phaseshift */
                phares2[shell] += phares_v*delta*delta;      /* Phase Residuum */
            }
        }
    }


    std::size_t n_shells;
    double *ampl1, *ampl2, *ampl_diff, *phares1, *phares2;
    std::size_t *n_shell;
    TPRECISION_RING shell_thickness;
    signed long center;
    signed long sy, sz;
    TPRECISION size3;
private:
    fourier_shell_correlation_(const fourier_shell_correlation_ &v)
        :   n_shells(0),
            ampl1(0),
            ampl2(0),
            ampl_diff(0),
            phares1(0),
            phares2(0),
            n_shell(0),
            center(0),
            sy(0),
            sz(0) {
        assert(!"HIDDEN");
    }
};
} // namespace
} // namespace functor
} // namespace tom
/****************************************************************************//**
 * The function work fastest, it the volumes are already fft-shifted and
 * reduced complex. Otherwise temporary copies have to be made.
 *
 * The F1_ and F2_ must not be normalized (e.g. to mean0+1std).
 *******************************************************************************/
template<typename T>
void tom::fourier_shell_correlation(const tom::Volume<std::complex<T> > &F1_, bool is_fft_shifted1, const tom::Volume<std::complex<T> > &F2_, bool is_fft_shifted2, std::size_t size, std::size_t n_shells, std::vector<double> &out) {

    if (F1_.getSizeX()!=size || F2_.getSizeX()!=size || F1_.getSizeY()!=size || F2_.getSizeY()!=size || (F1_.getSizeZ()!=size && F1_.getSizeZ()!=size/2+1) || (F2_.getSizeZ()!=size && F2_.getSizeZ()!=size/2+1)) {
        std::invalid_argument("The input volumes must have the same size (the z-dimension can opionally be reduced).");
    }
    const tom::Volume<std::complex<T> > *F1 = &F1_;
    const tom::Volume<std::complex<T> > *F2 = &F2_;

    // Prepare the input volumes so that they are reduced complex and fftshifted
    std::auto_ptr<tom::Volume<std::complex<T> > > F1p, F1p_half, F2p, F2p_half, vtmp;
    if (F1_.getSizeZ() == size/2+1) {
        if (!is_fft_shifted1) {
            vtmp.reset(    new tom::Volume<std::complex<T> >(size, size, size, NULL,NULL));
            F1p.reset(     new tom::Volume<std::complex<T> >(size, size, size, NULL,NULL));
            F1p_half.reset(new tom::Volume<std::complex<T> >(*F1p, NULL, size, size, size/2+1, F1p->getStrideX(), F1p->getStrideY(), F1p->getStrideZ()));
            tom::hermitian_symmetry_to_full(F1_, *vtmp);
            tom::fftshift(*vtmp, *F1p, false);
            F1 = F1p_half.get();
        }
    } else {
        if (is_fft_shifted1) {
            F1p_half.reset(new tom::Volume<std::complex<T> >(const_cast<std::complex<T> *>(&F1_.get()), size, size, size/2+1, F1_.getStrideX(), F1_.getStrideY(), F1_.getStrideZ(), false, NULL));
        } else {
            F1p.reset(     new tom::Volume<std::complex<T> >(size, size, size, NULL,NULL));
            F1p_half.reset(new tom::Volume<std::complex<T> >(*F1p, NULL, size, size, size/2+1, F1p->getStrideX(), F1p->getStrideY(), F1p->getStrideZ()));
            tom::fftshift(F1_, *F1p, false);
        }
        F1 = F1p_half.get();
    }
    if (F2_.getSizeZ() == size/2+1) {
        if (!is_fft_shifted2) {
            if (!vtmp.get()) { vtmp.reset(    new tom::Volume<std::complex<T> >(size, size, size, NULL,NULL)); }
            F2p.reset(     new tom::Volume<std::complex<T> >(size, size, size, NULL,NULL));
            F2p_half.reset(new tom::Volume<std::complex<T> >(*F2p, NULL, size, size, size/2+1, F2p->getStrideX(), F2p->getStrideY(), F2p->getStrideZ()));
            tom::hermitian_symmetry_to_full(F2_, *vtmp);
            tom::fftshift(*vtmp, *F2p, false);
            F2 = F2p_half.get();
        }
    } else {
        if (is_fft_shifted2) {
            F2p_half.reset(new tom::Volume<std::complex<T> >(const_cast<std::complex<T> *>(&F2_.get()), size, size, size/2+1, F2_.getStrideX(), F2_.getStrideY(), F2_.getStrideZ(), false, NULL));
        } else {
            F2p.reset(     new tom::Volume<std::complex<T> >(size, size, size, NULL,NULL));
            F2p_half.reset(new tom::Volume<std::complex<T> >(*F2p, NULL, size, size, size/2+1, F2p->getStrideX(), F2p->getStrideY(), F2p->getStrideZ()));
            tom::fftshift(F2_, *F2p, false);
        }
        F2 = F2p_half.get();
    }
    vtmp.reset(NULL);

    tom::functor::fourier_shell_correlation_<T> s(n_shells, size);
    tom::loop::for_each_step<   const tom::Volume<std::complex<T> >,
                                const tom::Volume<std::complex<T> >,
                                tom::functor::fourier_shell_correlation_<T> &
                            >(*F1, *F2, s);

    std::size_t ps = 0;
    out.resize(10 * n_shells);
    for (std::size_t i=0; i<n_shells; i++) {
        double ccc, mean;

        ccc = s.ampl1[i] * s.ampl2[i];
        assert(ccc >= 0);
        ccc  = (ccc>0)          ? (s.ampl1[i] + s.ampl2[i] - s.ampl_diff[i]) / sqrt(ccc) / 2.
                                : 1.;
        mean = (s.n_shell[i])   ? (s.ampl1[i] + s.ampl2[i] - s.ampl_diff[i]) / (s.n_shell[i]*2)
                                : 0.;
        if(mean < 0) {
            mean = -mean;
        }
        mean = sqrt(mean);
        const double rmsd = tom::math::sqrt<double>(s.ampl_diff[i] / (2*(s.ampl1[i] + s.ampl2[i]) - s.ampl_diff[i]));
        s.ampl1[i] = tom::math::sqrt<double>(s.ampl1[i] / s.n_shell[i]);
        s.ampl2[i] = tom::math::sqrt<double>(s.ampl2[i] / s.n_shell[i]);
        s.ampl_diff[i] = tom::math::sqrt<double>(s.ampl_diff[i]/s.n_shell[i]);
        const double sqrn = 2/tom::math::sqrt<double>(s.n_shell[i]);                 /* 2 / sqrt(N) */
        const double phares = tom::math::sqrt<double>(s.phares2[i]/s.phares1[i])*57.29577951308232087665461840231273527024; /* 180/pi */      /* Phase Residuum in degrees */

        out[i+n_shells*0] = s.n_shell[i];
        out[i+n_shells*1] = 2.*static_cast<double>(n_shells)/static_cast<double>(i+1);               /* pixel resolution */
        out[i+n_shells*2] = s.ampl_diff[i];
        out[i+n_shells*3] = rmsd;
        out[i+n_shells*4] = s.ampl1[i];
        out[i+n_shells*5] = s.ampl2[i];
        out[i+n_shells*6] = mean;
        out[i+n_shells*7] = ccc;
        out[i+n_shells*8] = sqrn;
        out[i+n_shells*9] = phares;

        ps += s.n_shell[i];
    }
}




/***************************************************************************//**
 * \param[in,out] v Volume
 * \param[in] sigma Sigma (sqrt(variance) of the gaussian noise. The noise has
 *  a mean of zero.
 *
 * Adds to each voxel of the volume \c v Gaussian noise with a standard deviation
 * of \c sigma.
 ******************************************************************************/
template<typename T>
void tom::gaussian_noise_add(tom::Volume<T> &v, double sigma) {

    // Create a Mersenne twister random number generator
    // that is seeded once with #seconds since 1970
    static boost::mt19937 rng(std::time(0));

    // select Gaussian probability distribution
    boost::normal_distribution<T> norm_dist(0, sigma);

    // bind random number generator to distribution, forming a function
    typedef boost::variate_generator<boost::mt19937 &, boost::normal_distribution<T> > t_gen;
    t_gen normal_sampler(rng, norm_dist);

    using boost::lambda::_1;
    using boost::lambda::bind;
    using boost::ref;

    // - The static_cast is necessary to choose the proper overloaded version of the operator().
    // - To call the operator() of normal_sampler, it must be bind with boost::lambda::bind.
    // - Moreover, the object parameter must be either a pointer or wrapped into boost::ref
    //   http://www.boost.org/doc/libs/1_38_0/doc/html/lambda/le_in_details.html#member_functions_as_targets
    tom::loop::for_each_no_idx( v,
        _1 += bind(static_cast<typename t_gen::result_type (t_gen::*)()>(&t_gen::operator()), ref(normal_sampler)) );

}

/***************************************************************************//**
 * \param[in,out] v Volume
 * \param[in] mean Mean of the gaussian noise.
 * \param[in] sigma Sigma (sqrt(variance) of the gaussian noise.
 *
 * Assigns noise to each voxel of the volume \c v Gaussian noise with a mean of \c mean and standard deviation
 * of \c sigma.
 ******************************************************************************/
template<typename T>
void tom::gaussian_noise(tom::Volume<T> &v, double mean, double sigma) {

	// Create a Mersenne twister random number generator
	// that is seeded once with #seconds since 1970
	static boost::mt19937 rng(std::time(0));

	// select Gaussian probability distribution
	boost::normal_distribution<T> norm_dist(mean, sigma);

    // bind random number generator to distribution, forming a function
    typedef boost::variate_generator<boost::mt19937 &, boost::normal_distribution<T> > t_gen;
    t_gen normal_sampler(rng, norm_dist);

    using boost::lambda::_1;
    using boost::lambda::bind;
    using boost::ref;


    // - The static_cast is necessary to choose the proper overloaded version of the operator().
    // - To call the operator() of normal_sampler, it must be bind with boost::lambda::bind.
    // - Moreover, the object parameter must be either a pointer or wrapped into boost::ref
    //   http://www.boost.org/doc/libs/1_38_0/doc/html/lambda/le_in_details.html#member_functions_as_targets
    tom::loop::for_each_no_idx( v,
            _1 = bind(static_cast<typename t_gen::result_type (t_gen::*)()>(&t_gen::operator()), ref(normal_sampler)) );


}




/***************************************************************************//**
 *
 ******************************************************************************/
template<typename T>
void tom::element_wise_add(tom::Volume<T> &v1, T factor_for_v1, const tom::Volume<T> &v2, T factor_for_v2) {

    using boost::lambda::_1;
    using boost::lambda::_2;
    if (factor_for_v1==0 && factor_for_v2==0) {
        v1.setValues(0);
    } else if (factor_for_v1==0) {
        tom::loop::for_each_no_idx(v1, v2, _1 = _2 * factor_for_v2);
    } else if (factor_for_v2==0) {
        v1.template scale<T>(factor_for_v1);
    } else {
        if (factor_for_v1==1 && factor_for_v2==1) {
            tom::element_wise_add(v1, v2);
        } else if (factor_for_v1==1) {
            tom::loop::for_each_no_idx(v1, v2, _1 = _1               + _2*factor_for_v2);
        } else if (factor_for_v2==1) {
            tom::loop::for_each_no_idx(v1, v2, _1 = _1*factor_for_v1 + _2              );
        } else {
            tom::loop::for_each_no_idx(v1, v2, _1 = _1*factor_for_v1 + _2*factor_for_v2);
        }
    }

}




namespace local {
namespace functor {
namespace {
template<typename T>
struct fill_noise_from_annulus_first_run {
    typedef boost::int_fast32_t tint;
    fill_noise_from_annulus_first_run(std::size_t sizex, std::size_t sizey, std::size_t sizez,
        const double &inner_radius, const double &outer_radius)
        :   inner_radius_(inner_radius), outer_radius_(outer_radius),
            centerx_(sizex/2), centery_(sizey/2), centerz_(sizez/2),
            current_z_(0), current_yz_(0),
            n_(0), sum_(0), sum2_(0) { }
    const double &inner_radius_, &outer_radius_;
    tint centerx_, centery_, centerz_;
    tint current_z_, current_yz_;
    std::size_t n_;
    double sum_, sum2_;
    inline void stepz(tint z) { current_z_  = z - centerz_; current_z_ *= current_z_; }
    inline void stepy(tint y) { current_yz_ = y - centery_; current_yz_ = current_yz_*current_yz_ + current_z_; }
    inline void operator()(const T &v, char &mask, tint x) {
        x -= centerx_;
        const double radius = tom::math::sqrt<double>(x*x + current_yz_);
        if (!(mask=(radius > outer_radius_)) && radius>=inner_radius_) {
            n_++;
            const double v2 = v;
            sum_  += v2;
            sum2_ += v2*v2;
        }
    }
};
} // namespace
} // namespace functor
} // namespace local
/***************************************************************************//**
 *
 * A radius equal to inner_radius or outer_radius is considered as inside the
 * annulus. Hence, the voxel is used for calculating the variance.
 ******************************************************************************/
template<typename T>
void tom::fill_noise_from_annulus(tom::Volume<T> &v, double inner_radius, double outer_radius, double sigma_factor) {
    assert(!(inner_radius<0 || outer_radius<0));
    assert(inner_radius==0 || !(inner_radius>=outer_radius));

    std::size_t sx = v.getSizeX()/2;
    std::size_t sy = v.getSizeY()/2;
    std::size_t sz = v.getSizeZ()/2;

    if (outer_radius > 0 && outer_radius < tom::math::sqrt<double>(sx*sx+sy*sy+sz*sz)) {

        tom::Volume<char> mask(v.getSizeX(), v.getSizeY(), v.getSizeZ(), 0,0);
        local::functor::fill_noise_from_annulus_first_run<T> f1(v.getSizeX(), v.getSizeY(), v.getSizeZ(), inner_radius, outer_radius);
        tom::loop::for_each_step<   const tom::Volume<T>,
                                    tom::Volume<char>,
                                    local::functor::fill_noise_from_annulus_first_run<T> &
                                >(v, mask, f1);
        if (f1.n_ > 0) {
            const double mean     = f1.sum_ / f1.n_;
            using boost::lambda::_1;
            using boost::lambda::_2;
            using boost::lambda::bind;
            using boost::lambda::if_then;
            using boost::ref;
            if (sigma_factor <= 0) {
                tom::loop::for_each_no_idx( v, mask, if_then(_2!=0, _1 = static_cast<const T>(mean)) );
            } else {
                const bool use_sample_standard_deviation = false;
                const double variance = (f1.sum2_ - f1.n_*mean*mean) / static_cast<double>(f1.n_ - (use_sample_standard_deviation ? 1 : 0));

                // Create a Mersenne twister random number generator
                // that is seeded once with #seconds since 1970
                static boost::mt19937 rng(std::time(0));

                // select Gaussian probability distribution
                boost::normal_distribution<T> norm_dist(mean, sigma_factor * tom::math::sqrt(variance));

                // bind random number generator to distribution, forming a function
                typedef boost::variate_generator<boost::mt19937 &, boost::normal_distribution<T> > t_gen;
                t_gen normal_sampler(rng, norm_dist);

                // - The static_cast is necessary to choose the proper overloaded version of the operator().
                // - To call the operator() of normal_sampler, it must be bind with boost::lambda::bind.
                // - Moreover, the object parameter must be either a pointer or wrapped into boost::ref
                //   http://www.boost.org/doc/libs/1_38_0/doc/html/lambda/le_in_details.html#member_functions_as_targets
                tom::loop::for_each_no_idx( v, mask,
                    if_then(_2!=0, _1 = bind(static_cast<typename t_gen::result_type (t_gen::*)()>(&t_gen::operator()), ref(normal_sampler)) ) );
            }
        }
    }
}

/****************************************************************************//**
 * \brief Returns a copy !!! of input volume as a vector with dimensions x*y*z,1,1
 * \param[in] v Volume
 * This function is useful for later eigenvalue / singular value decomposition ....
 *******************************************************************************/
template<typename T>
tom::Volume<T> tom::vectorize(tom::Volume<T> &source){

	const std::size_t x = source.getSizeX();
	const std::size_t y = source.getSizeY();
	const std::size_t z = source.getSizeZ();

	tom::Volume<T> returnVolume(x*y*z,1,1,NULL,NULL);

	tom::Volume<T> viewWindow(source,(void*)&source.get(),x*y*z,1,1,0,0,0);

	returnVolume.setValues(viewWindow);

	return returnVolume;
}


/****************************************************************************//**
 * \brief Returns subregion of the source
 * \param[in] v Volume
 * \param[in] subregion coordinates = selfexplaining
 * Access easyly a subregion of source
 *******************************************************************************/
template<typename T>
tom::Volume<T> tom::getSubregion(tom::Volume<T> &source,std::size_t startX,std::size_t startY,std::size_t startZ,std::size_t sizeX,std::size_t sizeY,std::size_t sizeZ){

	if(startX < 0) throw std::runtime_error(std::string("tom::getSubregion - StartX is <0 ?!"));
	if(startY < 0) throw std::runtime_error(std::string("tom::getSubregion - StartY is <0 ?!"));
	if(startZ < 0) throw std::runtime_error(std::string("tom::getSubregion - StartZ is <0 ?!"));

	if(startX + sizeX > source.getSizeX()) throw std::runtime_error(std::string("tom::getSubregion - startX +sizeX is larger than source.sizeX!"));
	if(startY + sizeY > source.getSizeY()) throw std::runtime_error(std::string("tom::getSubregion - startY +sizeY is larger than source.sizeY!"));
	if(startZ + sizeZ > source.getSizeZ()) throw std::runtime_error(std::string("tom::getSubregion - startZ +sizeZ is larger than source.sizeZ!"));

	tom::Volume<T> returnVolume(sizeX,sizeY,sizeZ,NULL,NULL);

	tom::Volume<T> viewWindow(source,(void*)&source.get(startX,startY,startZ),sizeX,sizeY,sizeZ,source.getStrideX(),source.getStrideY(),source.getStrideZ());

	returnVolume.setValues(viewWindow);
	return returnVolume;
}

/****************************************************************************//**
 * \brief Writes one volume into the region of the other
 * \param[in] v Volume
 * \param[in] destination coordinates = selfexplaining
 * Write smaller volume into a bigger one!
 *******************************************************************************/
template<typename T>
void tom::putSubregion(tom::Volume<T> &source, tom::Volume<T> &destination,std::size_t positionX,std::size_t positionY,std::size_t positionZ){

	if(positionX < 0) throw std::runtime_error(std::string("tom::putSubregion - positionX is <0 ?!"));
	if(positionY < 0) throw std::runtime_error(std::string("tom::putSubregion - positionY is <0 ?!"));
	if(positionZ < 0) throw std::runtime_error(std::string("tom::putSubregion - positionZ is <0 ?!"));

	if(positionX > destination.getSizeX()) throw std::runtime_error(std::string("tom::putSubregion - positionX is larger than destination.sizeX!"));
	if(positionY > destination.getSizeY()) throw std::runtime_error(std::string("tom::putSubregion - positionY is larger than destination.sizeY!"));
	if(positionZ > destination.getSizeZ()) throw std::runtime_error(std::string("tom::putSubregion - positionZ is larger than destination.sizeZ!"));

	if(positionX + source.getSizeX() > destination.getSizeX()) throw std::runtime_error(std::string("tom::putSubregion - positionX + source.getSizeX is larger than destination.sizeX!"));
	if(positionY + source.getSizeY() > destination.getSizeY()) throw std::runtime_error(std::string("tom::putSubregion - positionY + source.getSizeY is larger than destination.sizeY!"));
	if(positionZ + source.getSizeZ() > destination.getSizeZ()) throw std::runtime_error(std::string("tom::putSubregion - positionZ + source.getSizeZ is larger than destination.sizeZ!"));

	tom::Volume<T> viewWindow(destination,(void*)&destination.get(positionX,positionY,positionZ),source.getSizeX(),source.getSizeY(),source.getSizeZ(),destination.getStrideX(),destination.getStrideY(),destination.getStrideZ());

	viewWindow.setValues(source);
}

/****************************************************************************//**
 * \brief Calculates the sum of the elementwise product of two volumes.
 *
 * This function is useful as it returns the CC between two volumes (times numel)
 *******************************************************************************/
template<typename T, typename TPRECISION>
TPRECISION tom::sum_of_element_wise_product(const tom::Volume<T> &v1, const tom::Volume<T> &v2) {
    using boost::lambda::_1;
    using boost::lambda::_2;
    using boost::lambda::ll_static_cast;
    TPRECISION sum = 0.;
    tom::loop::for_each_no_idx(v1, v2, sum += ll_static_cast<TPRECISION>(_1)*ll_static_cast<TPRECISION>(_2));
    return sum;
}




// template instantiation
template int tom::get_tom_io_type<int                  >();
template int tom::get_tom_io_type<double               >();
template int tom::get_tom_io_type<fftw_complex         >();
template int tom::get_tom_io_type<std::complex<double> >();

template void tom::back_project<float>(tom::Volume<float> &src, tom::Volume<float> &dst,
									   int volumeSizeX, int volumeSizeY, int volumeSizeZ,
									   float Z1, float Y, float Z2,
									   int projectionSizeX, int projectionSizeY,
									   float x_offset, float y_offset, float z_offset,
									   float projectionOffsetX,float projectionOffsetY);


template void tom::norm_mask<float , float , float>(tom::Volume<float > &v, const tom::Volume<float > &mask, tom::norm::ntype stddev_type, float *variance, bool is_boolean_mask);
template void tom::norm_mask<float , float , double>(tom::Volume<float > &v, const tom::Volume<float > &mask, tom::norm::ntype stddev_type, double *variance, bool is_boolean_mask);
template void tom::norm_mask<double, double, double>(tom::Volume<double> &v, const tom::Volume<double> &mask, tom::norm::ntype stddev_type, double *variance, bool is_boolean_mask);
template void tom::norm_mask<double, char  , double>(tom::Volume<double> &v, const tom::Volume<char  > &mask, tom::norm::ntype stddev_type, double *variance, bool is_boolean_mask);

template void tom::fftshift<char                 >(const Volume<char                 > &vsrc, Volume<char                 > &v, bool is_ifftshift);
template void tom::fftshift<float                >(const Volume<float                > &vsrc, Volume<float                > &v, bool is_ifftshift);
template void tom::fftshift<double               >(const Volume<double               > &vsrc, Volume<double               > &v, bool is_ifftshift);
template void tom::fftshift<std::complex<float > >(const Volume<std::complex<float > > &vsrc, Volume<std::complex<float > > &v, bool is_ifftshift);
template void tom::fftshift<std::complex<double> >(const Volume<std::complex<double> > &vsrc, Volume<std::complex<double> > &v, bool is_ifftshift);

template std::vector<tom::st_idx> tom::peak<float >(const tom::Volume<float > &v);
template std::vector<tom::st_idx> tom::peak<double>(const tom::Volume<double> &v);
template std::vector<tom::st_idx> tom::peak<float , char>(const tom::Volume<float > &v, const tom::Volume<char> &mask);
template std::vector<tom::st_idx> tom::peak<double, char>(const tom::Volume<double> &v, const tom::Volume<char> &mask);
template std::vector<tom::st_idx> tom::peak<float , float >(const tom::Volume<float > &v, const tom::Volume<float > &mask);
template std::vector<tom::st_idx> tom::peak<double, double>(const tom::Volume<double> &v, const tom::Volume<double> &mask);


template void tom::hermitian_symmetry_to_full<float                >(const tom::Volume<float                > &vsrc, tom::Volume<float                > &v);
template void tom::hermitian_symmetry_to_full<double               >(const tom::Volume<double                > &vsrc, tom::Volume<double              > &v);
template void tom::hermitian_symmetry_to_full<std::complex<float > >(const tom::Volume<std::complex<float > > &vsrc, tom::Volume<std::complex<float > > &v);
template void tom::hermitian_symmetry_to_full<std::complex<double> >(const tom::Volume<std::complex<double> > &vsrc, tom::Volume<std::complex<double> > &v);

template void tom::init_sphere<char  >(tom::Volume<char  > &mask, float radius, float sigma, float max_radius, float centerx, float centery, float centerz);
template void tom::init_sphere<float >(tom::Volume<float > &mask, float radius, float sigma, float max_radius, float centerx, float centery, float centerz);
template void tom::init_sphere<double>(tom::Volume<double> &mask, float radius, float sigma, float max_radius, float centerx, float centery, float centerz);

template void tom::init_gaussian_sphere<float >(tom::Volume<float > &v, double sigma, double centerx, double centery, double centerz);
template void tom::init_gaussian_sphere<double>(tom::Volume<double> &v, double sigma, double centerx, double centery, double centerz);




template void tom::element_wise_add<float               , float                >(tom::Volume<float                > &b, const tom::Volume<float                > &a);
template void tom::element_wise_add<double              , double               >(tom::Volume<double               > &b, const tom::Volume<double               > &a);
template void tom::element_wise_add<std::complex<float >, std::complex<float > >(tom::Volume<std::complex<float > > &b, const tom::Volume<std::complex<float > > &a);
template void tom::element_wise_add<std::complex<double>, std::complex<double> >(tom::Volume<std::complex<double> > &b, const tom::Volume<std::complex<double> > &a);

template void tom::element_wise_add<float >(tom::Volume<float > &v1, float  factor_for_v1, const tom::Volume<float > &v2, float  factor_for_v2);
template void tom::element_wise_add<double>(tom::Volume<double> &v1, double factor_for_v1, const tom::Volume<double> &v2, double factor_for_v2);


template void tom::element_wise_sub<float , float >(tom::Volume<float > &b, const tom::Volume<float > &a);
template void tom::element_wise_sub<double, double>(tom::Volume<double> &b, const tom::Volume<double> &a);
template void tom::element_wise_sub<std::complex<float >, std::complex<float > >(tom::Volume<std::complex<float > > &b, const tom::Volume<std::complex<float > > &a);
template void tom::element_wise_sub<std::complex<double>, std::complex<double> >(tom::Volume<std::complex<double> > &b, const tom::Volume<std::complex<double> > &a);

template void tom::element_wise_multiply<std::complex<float >, float >(tom::Volume<std::complex<float > > &b, const tom::Volume<float > &a);
template void tom::element_wise_multiply<std::complex<float >, double>(tom::Volume<std::complex<float > > &b, const tom::Volume<double> &a);
template void tom::element_wise_multiply<std::complex<double>, float >(tom::Volume<std::complex<double> > &b, const tom::Volume<float > &a);
template void tom::element_wise_multiply<std::complex<double>, double>(tom::Volume<std::complex<double> > &b, const tom::Volume<double> &a);

template void tom::element_wise_multiply<float,float>(tom::Volume<float > &b, const tom::Volume<float > &a);
template void tom::element_wise_multiply<double,double>(tom::Volume<double > &b, const tom::Volume<double > &a);

template void tom::element_wise_multiply<std::complex<float >, std::complex<float > >(tom::Volume<std::complex<float > > &b, const tom::Volume<std::complex<float > > &a);
template void tom::element_wise_multiply<std::complex<double>, std::complex<double> >(tom::Volume<std::complex<double> > &b, const tom::Volume<std::complex<double> > &a);


template void tom::element_wise_div<std::complex<float >, float >(tom::Volume<std::complex<float > > &b, const tom::Volume<float > &a, std::complex<float > inf_value);
template void tom::element_wise_div<std::complex<double>, double>(tom::Volume<std::complex<double> > &b, const tom::Volume<double> &a, std::complex<double> inf_value);
template void tom::element_wise_div<float , float>(tom::Volume<float>  &b, const tom::Volume<float> &a, float inf_value);
template void tom::element_wise_div<double , double>(tom::Volume<double> &b, const tom::Volume<double> &a, double inf_value);
template void tom::element_wise_div<std::complex<float >, std::complex<float> >(tom::Volume<std::complex<float > > &b, const tom::Volume<std::complex<float> > &a, std::complex<float > inf_value);
template void tom::element_wise_div<std::complex<double>, std::complex<double> >(tom::Volume<std::complex<double> > &b, const tom::Volume<std::complex<double> > &a, std::complex<double> inf_value);


//template void tom::element_wise_conj_multiply<float,float>(tom::Volume<float > &b, const tom::Volume<float > &a);
//template void tom::element_wise_conj_multiply<double,double>(tom::Volume<double>&b, const tom::Volume<double> &a);
template void tom::element_wise_conj_multiply<std::complex<float >, std::complex<float > >(tom::Volume<std::complex<float > > &b, const tom::Volume<std::complex<float > > &a);
template void tom::element_wise_conj_multiply<std::complex<double>, std::complex<double> >(tom::Volume<std::complex<double> > &b, const tom::Volume<std::complex<double> > &a);


//template void tom::element_wise_conj_multiply<float , float , float >(tom::Volume<float> &v1, const tom::Volume<float> &v2, const tom::Volume<float>&v3);
//template void tom::element_wise_conj_multiply<double, double, double>(tom::Volume<double> &v1, const tom::Volume<double> &v2, const tom::Volume<double>&v3);
template void tom::element_wise_conj_multiply<std::complex<float >, std::complex<float >, std::complex<float > >(tom::Volume<std::complex<float > > &v1, const tom::Volume<std::complex<float > > &v2, const tom::Volume<std::complex<float > > &v3);
template void tom::element_wise_conj_multiply<std::complex<double>, std::complex<double>, std::complex<double> >(tom::Volume<std::complex<double> > &v1, const tom::Volume<std::complex<double> > &v2, const tom::Volume<std::complex<double> > &v3);


template void tom::conjugate<std::complex<float> >(tom::Volume<std::complex<float> > &v);
template void tom::conjugate<std::complex<double> >(tom::Volume<std::complex<double> > &v);


template void tom::element_wise_power<float >(tom::Volume<float > &v, float  exponent);
template void tom::element_wise_power<double>(tom::Volume<double> &v, double exponent);

template void tom::element_wise_power<std::complex<float>, float >(tom::Volume<std::complex<float> > &v, float  exponent);
//template void tom::element_wise_power<std::complex<double>, double>(tom::Volume<std::complex<double> > &v, double exponent);

template void tom::element_wise_max<float >(tom::Volume<float > &v, float  val);
template void tom::element_wise_max<double>(tom::Volume<double> &v, double val);

template std::size_t tom::number_set_voxels(const tom::Volume<float> &v);
template std::size_t tom::number_set_voxels(const tom::Volume<double> &v);
/*
template std::size_t tom::number_set_voxels(const tom::Volume<std::complex<float> > &v);
template std::size_t tom::number_set_voxels(const tom::Volume<std::complex<double> > &v);
*/

template void tom::make_binary<float >(tom::Volume<float > &v);
template void tom::make_binary<double>(tom::Volume<double> &v);


template void tom::fourier_shell_correlation<float >(const tom::Volume<std::complex<float > > &F1, bool is_fft_shifted1, const tom::Volume<std::complex<float > > &F2, bool is_fft_shifted2, std::size_t size, std::size_t n_shells, std::vector<double> &r);
template void tom::fourier_shell_correlation<double>(const tom::Volume<std::complex<double> > &F1, bool is_fft_shifted1, const tom::Volume<std::complex<double> > &F2, bool is_fft_shifted2, std::size_t size, std::size_t n_shells, std::vector<double> &r);

template void tom::element_wise_set_below_threshold<float >(tom::Volume<float > &a, float  threshold, float  value);
template void tom::element_wise_set_below_threshold<double>(tom::Volume<double> &a, double threshold, double value);


template bool tom::isfinite<float                >(const tom::Volume<float                > &v);
template bool tom::isfinite<double               >(const tom::Volume<double               > &v);
template bool tom::isfinite<std::complex<float > >(const tom::Volume<std::complex<float > > &v);
template bool tom::isfinite<std::complex<double> >(const tom::Volume<std::complex<double> > &v);

template float  tom::sum_squared<float , float >(const tom::Volume<float > &v);
template double tom::sum_squared<float , double>(const tom::Volume<float > &v);
template double tom::sum_squared<double, double>(const tom::Volume<double> &v);


template void tom::gaussian_noise_add<float >(tom::Volume<float > &v, double sigma);
template void tom::gaussian_noise_add<double>(tom::Volume<double> &v, double sigma);

template void tom::gaussian_noise<float >(tom::Volume<float > &v, double mean, double sigma);
template void tom::gaussian_noise<double>(tom::Volume<double> &v, double mean, double sigma);

template tom::Volume<float> tom::vectorize(tom::Volume<float> &v);
template tom::Volume<double> tom::vectorize(tom::Volume<double> &v);

template tom::Volume<float> tom::getSubregion(tom::Volume<float> &source,std::size_t startX,std::size_t startY,std::size_t startZ,std::size_t sizeX,std::size_t sizeY,std::size_t sizeZ);
template tom::Volume<double> tom::getSubregion(tom::Volume<double> &source,std::size_t startX,std::size_t startY,std::size_t startZ,std::size_t sizeX,std::size_t sizeY,std::size_t sizeZ);

template tom::Volume<std::complex<float> > tom::getSubregion(tom::Volume<std::complex<float> > &source,std::size_t startX,std::size_t startY,std::size_t startZ,std::size_t sizeX,std::size_t sizeY,std::size_t sizeZ);
template tom::Volume<std::complex<double> > tom::getSubregion(tom::Volume<std::complex<double> > &source,std::size_t startX,std::size_t startY,std::size_t startZ,std::size_t sizeX,std::size_t sizeY,std::size_t sizeZ);

template void tom::putSubregion(tom::Volume<float > &source, tom::Volume<float>& destination,std::size_t positionX,std::size_t positionY,std::size_t positionZ);
template void tom::putSubregion(tom::Volume<double> &source, tom::Volume<double> &destination,std::size_t positionX,std::size_t positionY,std::size_t positionZ);

template void tom::putSubregion(tom::Volume<std::complex<float  > > &source, tom::Volume<std::complex<float> >& destination,std::size_t positionX,std::size_t positionY,std::size_t positionZ);
template void tom::putSubregion(tom::Volume<std::complex<double > > &source, tom::Volume<std::complex<double> >& destination,std::size_t positionX,std::size_t positionY,std::size_t positionZ);


template std::complex<float > tom::sum<std::complex<float > , std::complex<float > >(const Volume<std::complex<float > > &v);
template std::complex<double> tom::sum<std::complex<double> , std::complex<double> >(const Volume<std::complex<double> > &v);


template float          tom::sum<float , float        >(const Volume<float > &v);
template double         tom::sum<float , double       >(const Volume<float > &v);
template float          tom::sum<double, float        >(const Volume<double> &v);
template double         tom::sum<double, double       >(const Volume<double> &v);
template unsigned int   tom::sum<char  , unsigned int >(const Volume<char  > &v);
template unsigned long  tom::sum<char  , unsigned long>(const Volume<char  > &v);
template char           tom::min(const Volume<char  > &v);
template double         tom::min(const Volume<double> &v);
template float          tom::min(const Volume<float > &v);
template char           tom::max(const Volume<char  > &v);
template int            tom::max(const Volume<int   > &v);
template float          tom::max(const Volume<float > &v);
template double         tom::max(const Volume<double> &v);
template void tom::minmax(const Volume<float > &v, float  &min, float  &max);
template void tom::minmax(const Volume<double> &v, double &min, double &max);
template void tom::stat  (const Volume<float > &v, double &mean, double &variance, bool use_sample_standard_deviation);
template void tom::stat  (const Volume<double> &v, double &mean, double &variance, bool use_sample_standard_deviation);
template void tom::stat  (const Volume<double> &v, double &mean, double &variance, bool use_sample_standard_deviation, double &min, double &max);
template void tom::stat  (const Volume<double> &v, double &mean, double &variance, bool use_sample_standard_deviation, const Volume<char  > &mask);
template void tom::stat  (const Volume<float > &v, double &mean, double &variance, bool use_sample_standard_deviation, const Volume<float > &mask);
template void tom::stat  (const Volume<double> &v, double &mean, double &variance, bool use_sample_standard_deviation, const Volume<double> &mask);

template void tom::limit (Volume<float > &v,float  lowerBound,float  lowerReplacement,float  upperBound,float  upperReplacement, bool doLower, bool doUpper);
template void tom::limit (Volume<double> &v,double lowerBound,double lowerReplacement,double upperBound,double upperReplacement, bool doLower, bool doUpper);


template void tom::fill_noise_from_annulus(Volume<float > &v, double inner_radius, double outer_radius, double sigma_factor);
template void tom::fill_noise_from_annulus(Volume<double> &v, double inner_radius, double outer_radius, double sigma_factor);


template float  tom::sum_of_element_wise_product(const tom::Volume<float > &v1, const tom::Volume<float > &v2);
template double tom::sum_of_element_wise_product(const tom::Volume<double> &v1, const tom::Volume<double> &v2);
