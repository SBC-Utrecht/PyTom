/****************************************************************************//**
 * \file volume_fcn.hpp
 * \brief The header file for single functions manipulating tom::Volume.
 * \author  Thomas Haller
 * \version 0.1
 * \date    10.12.2007
 *******************************************************************************/
#ifndef ___INCLUDE_CORE__VOLUME_FCN_HPP__
#define ___INCLUDE_CORE__VOLUME_FCN_HPP__


#include <vector>
#include <complex>
#include <cmath>
#include <cstdlib>

#include <boost/type_traits/remove_reference.hpp>
#include <boost/type_traits/is_floating_point.hpp>

#include <fftw3.h>

#include <tom/volume.hpp>
#include <tom/tools/triple.hpp>
#include <tom/volume_loop.hpp>


/****************************************************************************//**
 * \brief base namespace of libtomc.
 *******************************************************************************/
namespace tom {


typedef tom::tools::triple<std::size_t, std::size_t, std::size_t> st_idx;


template<typename T> bool is_double();
template<typename T> bool is_float();
template<typename T> bool is_int();
template<typename T> bool is_float_complex();
template<typename T> bool is_double_complex();


template<typename T, typename TPRECISION> TPRECISION sum               (const Volume<T> &v);
template<typename T>                      double     sum               (const Volume<T> &v);
template<typename T, typename TPRECISION> TPRECISION sum_squared       (const Volume<T> &v);
template<typename T                     > double     sum_squared       (const Volume<T> &v);
template<typename T                     > double     mean              (const Volume<T> &v                                    );
template<typename T, typename TPRECISION> TPRECISION mean              (const Volume<T> &v                                    );
template<typename T                     > double     variance          (const Volume<T> &v, bool use_sample_standard_deviation);
template<typename T, typename TPRECISION> TPRECISION variance          (const Volume<T> &v, bool use_sample_standard_deviation);
template<typename T                     > double     variance_mean_free(const Volume<T> &v, bool use_sample_standard_deviation);
template<typename T, typename TPRECISION> TPRECISION variance_mean_free(const Volume<T> &v, bool use_sample_standard_deviation);
template<typename T                     > T          min               (const Volume<T> &v);
template<typename T                     > T          max               (const Volume<T> &v);
template<typename T                     > void       minmax            (const Volume<T> &v, T &min, T &max);
template<typename T, typename TPRECISION> void       stat              (const Volume<T> &v, TPRECISION &mean, TPRECISION &variance, bool use_sample_standard_deviation);
template<typename T, typename TPRECISION> void       stat              (const Volume<T> &v, TPRECISION &mean, TPRECISION &variance, bool use_sample_standard_deviation, T &min, T &max);

template<typename T, typename TMASK, typename TPRECISION>
    void stat(const Volume<T> &v, TPRECISION &mean, TPRECISION &variance, bool use_sample_standard_deviation, const Volume<TMASK> &mask);




template<typename T> void make_binary(tom::Volume<T> &v);


template<typename T> void init_sphere(tom::Volume<T> &mask, float radius, float sigma, float max_radius);
template<typename T> void init_sphere(tom::Volume<T> &mask, float radius, float sigma, float max_radius, float centerx, float centery, float centerz);

template<typename T> void limit(tom::Volume<T> &volume,T lowerBound, T lowerReplacement, T upperBound, T upperReplacement,bool doLower, bool doUpper);

template<typename T> void init_gaussian_sphere(tom::Volume<T> &mask, double sigma);
template<typename T> void init_gaussian_sphere(tom::Volume<T> &mask, double sigma, double centerx, double centery, double centerz);


template<typename T> void hermitian_symmetry_to_full(      tom::Volume<T> &v                      );
template<typename T> void hermitian_symmetry_to_full(const tom::Volume<T> &vsrc, tom::Volume<T> &v);



template<typename T> void fftshift(      tom::Volume<T> &v   ,                       bool is_ifftshift);
template<typename T> void fftshift(const tom::Volume<T> &vsrc, tom::Volume<T> &v   , bool is_ifftshift);

template<typename T> void  fftshift(      tom::Volume<T> &v                      );
template<typename T> void  fftshift(const tom::Volume<T> &vsrc, tom::Volume<T> &v);
template<typename T> void ifftshift(      tom::Volume<T> &v                      );
template<typename T> void ifftshift(const tom::Volume<T> &vsrc, tom::Volume<T> &v);


//template<typename T> void shift(tom::Volume<std::complex<T> > &v, std::size_t sizez, double shiftx, double shifty, double shiftz);

namespace norm {
    enum ntype {
        NORM_NO_NORM,
        NORM_STDDEV_1,
        NORM_STDDEV_SAMPLE_1
    };
}
template<typename T, typename T2, typename T3>
void norm_mask(tom::Volume<T> &v, const tom::Volume<T2> &mask, tom::norm::ntype stddev_type, T3 *variance, bool is_boolean_mask);

template<typename T, typename TIDX>
void update_correlation_volume(const Volume<T> &v, Volume<T> &cc, Volume<TIDX> &vindex, TIDX index);


template<typename T>              std::vector<tom::st_idx> peak(const tom::Volume<T> &v);
template<typename T, typename T2> std::vector<tom::st_idx> peak(const tom::Volume<T> &v, const tom::Volume<T2> &mask);


template<typename T> bool isfinite(const tom::Volume<T> &v);


template<typename T> void conjugate(tom::Volume<T> &v);


template<typename T> void element_wise_abs(tom::Volume<T> &v);
template<typename T1,typename T2> void element_wise_abs(tom::Volume<T1> &v1,const tom::Volume<T2> &v2);


template<typename T, typename T2> void element_wise_add          (tom::Volume<T> &b, const tom::Volume<T2> &a);
template<typename T             > void element_wise_add          (tom::Volume<T> &v1, T factor_for_v1, const tom::Volume<T> &v2, T factor_for_v2);
template<typename T, typename T2> void element_wise_sub          (tom::Volume<T> &b, const tom::Volume<T2> &a);
template<typename T, typename T2> void element_wise_multiply     (tom::Volume<T> &b, const tom::Volume<T2> &v);
template<typename T, typename T2> void element_wise_div          (tom::Volume<T> &b, const tom::Volume<T2> &a, T inf_value);
template<typename T, typename T2> void element_wise_conj_multiply(tom::Volume<T> &b, const tom::Volume<T2> &a);
template<typename T1, typename T2, typename T3> void element_wise_conj_multiply(tom::Volume<T1> &v1, const tom::Volume<T2> &v2, const tom::Volume<T3> &v3);

template<typename T, typename TPRECISION> TPRECISION sum_of_element_wise_product(const tom::Volume<T> &v1, const tom::Volume<T> &v2);


template<typename T> void element_wise_set_below_threshold(tom::Volume<T> &a, T threshold, T value);

template<typename T> std::size_t number_set_voxels(const tom::Volume<T> &v);

template<typename T> void element_wise_power(tom::Volume<T> &v, T exponent);
template<typename T,typename T2> void element_wise_power(tom::Volume<T> &v, T2 exponent);

template<typename T> void element_wise_max(tom::Volume<T> &v, T val);


template<typename T> void fourier_shell_correlation(const tom::Volume<std::complex<T> > &F1, bool is_fft_shifted1, const tom::Volume<std::complex<T> > &F2, bool is_fft_shifted2, std::size_t size, std::size_t n_shells, std::vector<double> &r);


template<typename T, typename TOP> void element_wise_operation(tom::Volume<T> &v, TOP o);


template<typename T> void gaussian_noise_add(tom::Volume<T> &v, double sigma);
template<typename T> void gaussian_noise(tom::Volume<T> &v, double mean, double sigma);

template<typename T>
void fill_noise_from_annulus(tom::Volume<T> &v, double inner_radius, double outer_radius, double sigma_factor);

template<typename T>
Volume<T> vectorize(tom::Volume<T> &source);

template<typename T>
Volume<T> getSubregion(tom::Volume<T> &source,std::size_t startX,std::size_t startY,std::size_t startZ,std::size_t endX,std::size_t endY,std::size_t endZ);

template<typename T>
void putSubregion(tom::Volume<T> &source, tom::Volume<T>& destination,std::size_t positionX,std::size_t positionY,std::size_t positionZ);

namespace math {
double deg2rad(double x);
double rad2deg(double x);
template<typename T> T log (T x);
template<typename T> T exp (T x);
template<typename T> T sqrt(T x);
template<typename T> T acos(T x);
template<typename T> T abs (T x);
template<typename T> T abs (std::complex<T> x);
template<typename T> T ceil(T x);
template<typename T> T pow (T x,T y);
}


template<typename T1, typename T2> struct equal_;
template<typename T1, typename T2>
inline bool equal(const T1 &v1, const T2 &v2) { return equal_<T1, T2>::eval(v1, v2); }


template<typename T> void back_project(tom::Volume<T> &src,
									   tom::Volume<T> &dst,
									   int Ox_max, int Oy_max, int Oz_max,
									   float Z1, float Y,float Z2,
									   int Ix_max, int Iy_max,
									   float x_offset, float y_offset, float z_offset,
									   float projectionOffsetX,float projectionOffsetY);

}








namespace tom {
/****************************************************************************//**
 * \brief contains inline template functions to coose the proper underlying mathematical operation.
 *
 * Contains template functions to call the appropriate floating point operation
 * depending on the data type (for example \c sin vs. \c sinf).
 *******************************************************************************/
namespace math {
template<          > inline float  log (float  x) { return ::logf(x); }
template<          > inline double log (double x) { return std::log(x); }
template<          > inline float  exp (float  x) { return ::expf(x); }
template<          > inline double exp (double x) { return std::exp(x); }
template<          > inline float  sqrt(float  x) { return ::sqrtf(x); }
template<          > inline double sqrt(double x) { return std::sqrt(x);  }
template<          > inline float  acos(float  x) { return ::acosf(x); }
template<          > inline double acos(double x) { return std::acos(x); }
template<          > inline float  ceil(float  x) { return ::ceilf(x); }
template<          > inline double ceil(double x) { return std::ceil(x); }
template< 		   > inline float  pow(float x,float y) { return powf(x,y); }
template< 		   > inline double pow(double x,double y) { return pow(x,y); }

template<          > inline long int abs (long int x) { return std::labs(x); }
template<          > inline int      abs (int      x) { return std::abs(x); }
template<          > inline float    abs (float    x) { return ::fabsf(x); }
template<          > inline double   abs (double   x) { return std::fabs(x); }

template<          > inline float    abs (std::complex<float >   x) { return std::abs(x); }
template<          > inline double   abs (std::complex<double>   x) { return std::abs(x); }




}


template<typename T1, typename T2> struct equal_                        { static bool eval(const T1               &v1, const T2  &v2) { return v1 == v2; } };
template<typename T1             > struct equal_<std::complex<T1>, int> { static bool eval(const std::complex<T1> &v1, const int &v2) { return v1.imag()==0 && v1.real()==v2; } };


}


inline double tom::math::deg2rad(double x) { return x *  0.017453292519943295769236907684886127134; }
inline double tom::math::rad2deg(double x) { return x * 57.295779513082320876798154814105170332   ; }


namespace tom {
/** \brief Returns true if the template function  is for int32 */
template <typename T> inline bool is_int     ()        { return false; }
template <          > inline bool is_int<int>()        { return true ; }
/** \brief Returns true if the template function  is for float */
template <typename T> inline bool is_float       ()         { return false; }
template <          > inline bool is_float<float>()         { return true ; }
/** \brief Returns true if the template function  is for double */
template <typename T> inline bool is_double        ()        { return false; }
template <          > inline bool is_double<double>()        { return true ; }
/** \brief Returns true if the template function is for float complex numbers */
template <typename T> inline bool is_float_complex                      () { return false; }
template <          > inline bool is_float_complex<std::complex<float> >() { return true ; }
/** \brief Returns true if the template function  is for double complex numbers */
template <typename T> inline bool is_double_complex                       ()  { return false; }
template <          > inline bool is_double_complex<std::complex<double> >()  { return true ; }
}


/****************************************************************************//**
 * \brief Makes an fftshift of the volume inplace.
 * \param[in,out] v
 * \param[in] is_ifftshift If true, an inverse fft-shift is performed (as in
 *      the \c ifftshift function in MATLAB). This differs in the choise (+/-0) of the
 *      centre for volumes with an even number of elements per dimension.
 *
 * \todo The current implementation copies the volume. Make a better implementation!!
 *******************************************************************************/
template<typename T>
inline void tom::fftshift(tom::Volume<T> &v, bool is_ifftshift) {
    tom::Volume<T> v2(v);
    tom::fftshift(v2, v, is_ifftshift);
}



/****************************************************************************//**
 * \brief initialises volume with sphere mask.
 * \param[out] mask
 * \param[in] radius
 * \param[in] sigma
 * \param[in] max_radius
 *
 * Calls init_sphere with the center parameters set to the center of the
 * volume. Beware, that the center is always integer and is rounded towards
 * infinity, in case of odd mask size.\n
 *******************************************************************************/
template<typename T>
inline void tom::init_sphere(tom::Volume<T> &mask, float radius, float sigma, float max_radius) {
    ::tom::init_sphere<T>(mask, radius, sigma, max_radius, mask.getSizeX() / 2, mask.getSizeY() / 2, mask.getSizeZ() / 2);
}



/****************************************************************************//**
 * \brief initialises volume with a Gausian sphere mask.
 * \param[out] mask
 * \param[in] radius
 * \param[in] sigma
 * \param[in] max_radius
 *
 * Calls init_gaussian_sphere with the center parameters set to the center of the
 * volume. Beware, that the center is always integer and is rounded towards
 * infinity, in case of odd mask size.\n
 *******************************************************************************/
template<typename T>
inline void tom::init_gaussian_sphere(tom::Volume<T> &mask, double sigma) {
    ::tom::init_gaussian_sphere<T>(mask, sigma, mask.getSizeX() / 2, mask.getSizeY() / 2, mask.getSizeZ() / 2);
}




namespace tom {
namespace functor {
namespace {
template<typename T, typename TOP>
struct element_wise_operation_ {
    element_wise_operation_(TOP &o): o_(o) { }
    void operator()(T &v, std::size_t, std::size_t, std::size_t) {
        v = o_(v);
    }
    TOP &o_;
};
} // namespace
} // namespace functor
} // namespace tom
/****************************************************************************//**
 * \brief Template function to perform an operation on each element.
 *
 *  \code
 *      tom::Volume<float> v(...);
 *      tom::element_wise_operation(v, &tom::math::sqrt<float>);
 *  \endcode
 *******************************************************************************/
template<typename T, typename TOP>
inline void tom::element_wise_operation(tom::Volume<T> &v, TOP o) {
    typedef typename boost::remove_reference<TOP>::type TOP_NOREF;
    tom::loop::for_each(v, tom::functor::element_wise_operation_<T, TOP_NOREF>(o));
}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
inline double tom::sum(const tom::Volume<T> &v) {
    return tom::sum<T, double>(v);
}
/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
inline double tom::sum_squared(const tom::Volume<T> &v) {
    return tom::sum_squared<T, double>(v);
}




/****************************************************************************//**
 * \brief Assigns absolute values to v1, inplace.
 *
 * \param[in,out] v1 The left side operand of the elementwise abs and the
 *   destination volume.
 *
 *******************************************************************************/
template<typename T>
inline void tom::element_wise_abs(tom::Volume<T> &v) {
    tom::element_wise_operation<T, T(*)(T)>(v, static_cast<T(*)(T)>(&tom::math::abs<T>));
}






namespace tom {
namespace functor {
namespace {
template<typename T1, typename T2>
struct element_wise_abs_ {
    void operator()(T1 &v1, const T2 &v2, std::size_t, std::size_t, std::size_t) { v1 = tom::math::abs(v2); }
};
} // namespace
} // namespace functor
} // namespace tom
/****************************************************************************//**
 * \brief Assigns the absolute value of v2 to v1.
 *
 * \param[in,out] v1 The destination volume.
 * \param[in] v2 The source for the abs function.
 *
 * v2 may be complex. The assigned values will be the magnitude of each complex voxel.
 *******************************************************************************/
template<typename T1, typename T2>
void tom::element_wise_abs(tom::Volume<T1> &v1, const tom::Volume<T2> &v2) {
    tom::loop::for_each(v1, v2, tom::functor::element_wise_abs_<T1, T2>());
}






namespace tom {
namespace functor {
namespace {
template<typename T, typename TIDX>
struct update_correlation_volume_ {
    update_correlation_volume_(TIDX index): index_(index) { }
    const TIDX index_;
    inline void operator()(const T &v, T &vmax, TIDX & vindex, std::size_t, std::size_t, std::size_t) const {
        if (vmax < v) {
            vmax = v;
            vindex = index_;
        }
    }
};
} // namespace
} // namespace functor
} // namespace tom
/***************************************************************************//**
 * \brief Take the elementwise maximum of two volumes.
 *
 * The result is saved in \c cc and if a new maximum at such a position is found,
 * \c index is saved at that particular position in \c vindex.
 ******************************************************************************/
template<typename T, typename TIDX>
inline void tom::update_correlation_volume(const tom::Volume<T> &v, tom::Volume<T> &cc, tom::Volume<TIDX> &vindex, TIDX index) {
    tom::loop::for_each(v, cc, vindex, tom::functor::update_correlation_volume_<T, TIDX>(index));
}




/***************************************************************************//**
 *
 ******************************************************************************/
template<typename T>
inline void tom::fftshift(tom::Volume<T> &v) {
    fftshift(v, false);
}


/***************************************************************************//**
 *
 ******************************************************************************/
template<typename T>
inline void tom::fftshift(const tom::Volume<T> &vsrc, tom::Volume<T> &v) {
    fftshift(vsrc, v, false);
}


/***************************************************************************//**
 *
 ******************************************************************************/
template<typename T>
inline void tom::ifftshift(tom::Volume<T> &v) {
    fftshift(v, true);
}


/***************************************************************************//**
 *
 ******************************************************************************/
template<typename T>
inline void tom::ifftshift(const tom::Volume<T> &vsrc, tom::Volume<T> &v) {
    fftshift(vsrc, v, true);
}



/****************************************************************************//**
 * \brief Compute the variance of the volume
 *
 * \param[in] v Volume
 * \param[in] use_sample_standard_deviation If true the denominator is N-1 instead
 *    of N (N number of voxels).
 *
 * \return The variance of the volume.
 *******************************************************************************/
template<typename T, typename TPRECISION>
inline TPRECISION tom::variance(const Volume<T> &v, bool use_sample_standard_deviation) {
    TPRECISION variance, mean;
    stat(v, mean, variance, use_sample_standard_deviation);
    return variance;
}


/****************************************************************************//**
 * \brief Compute the variance of the volume
 *
 * \param[in] v Volume
 * \param[in] use_sample_standard_deviation If true the denominator is N-1 instead
 *    of N (N number of voxels).
 *
 * \return The variance of the volume.
 *******************************************************************************/
template<typename T>
inline double tom::variance(const Volume<T> &v, bool use_sample_standard_deviation) {
    return variance<T, double>(v, use_sample_standard_deviation);
}



/****************************************************************************//**
 * \brief Compute the mean of the volume with precision double.
 *******************************************************************************/
template<typename T>
inline double tom::mean(const Volume<T> &v) {
    return mean<T, double>(v);
}


/****************************************************************************//**
 * \brief Compute the mean of the volume with precision double.
 *******************************************************************************/
template<typename T, typename TPRECISION>
inline TPRECISION tom::mean(const Volume<T> &v) {
    BOOST_STATIC_ASSERT((boost::is_floating_point<TPRECISION>::value));
    return sum<T, double>(v) / static_cast<TPRECISION>(v.numel());
}

/****************************************************************************//**
 * Compute the variance of a volume where it is known that the mean is zero.
 * Performs the computation with double precision.
 *******************************************************************************/
template<typename T>
inline double tom::variance_mean_free(const Volume<T> &v, bool use_sample_standard_deviation) {
    return variance_mean_free<T, double>(v, use_sample_standard_deviation);
}



/****************************************************************************//**
 * \brief Compute the variance of a volume where it is known that the mean is zero.
 *
 * \warning if the mean is not really zero, the result will be wrong.
 *******************************************************************************/
template<typename T, typename TPRECISION>
inline TPRECISION tom::variance_mean_free(const Volume<T> &v, bool use_sample_standard_deviation) {
    BOOST_STATIC_ASSERT((boost::is_floating_point<TPRECISION>::value));
    return sum_squared<T, TPRECISION>(v) / static_cast<TPRECISION>(v.numel() - (use_sample_standard_deviation ? 1 : 0));
}





#endif

