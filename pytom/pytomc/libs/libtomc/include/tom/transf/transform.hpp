/***********************************************************************//**
 * \file transform.hpp
 * \brief
 * \author  Thomas Haller
 * \version 0.2
 * \date    12.12.2007
 **************************************************************************/
#ifndef ___INCLUDE__TOM__TRANSFORM__TRANSFORM_HPP__
#define ___INCLUDE__TOM__TRANSFORM__TRANSFORM_HPP__


#include <fftw3.h>

#include <tom/volume.hpp>


namespace tom {


/***********************************************************************//**
 * \brief namespace for transformation, interpolation and rotation.
 *
 * namespace for transformation, interpolation and rotation.
 **************************************************************************/
namespace transf {



int invert_4x4(const double *P, double *Pinv, int N);

void sum_rotation(double *P, bool homogenous, bool continue_P, double *angles_out, double *shifts_out, std::size_t n, const double *angles, const int *axes, const double *shifts);





/***********************************************************************//**
 * \brief exception class.
 * The interpolations functions should throw this type of exception
 * during \c setvolume, if they don't support the memory alignmet (stride
 * parameters) of the volume.
 **************************************************************************/
class memory_alignment_not_supported: public std::invalid_argument {
public:
    memory_alignment_not_supported(const std::string &s): std::invalid_argument(s) {
    }
};



/***********************************************************************//**
 * Class for tom::transf::transform which performs nearest neighbour interpolation.
 **************************************************************************/
template<typename T>
class InterpolNearestNeighbour {

public:
    typedef float idx_type;
    InterpolNearestNeighbour(T defaultval);
    void setvolume(const tom::Volume<T> &src);
    void setvolume(const T *data,const std::size_t sizex,const std::size_t sizey,const std::size_t sizez,const std::size_t stridey,const std::size_t stridez);
    T interpolate(typename tom::transf::InterpolNearestNeighbour<T>::idx_type x, typename tom::transf::InterpolNearestNeighbour<T>::idx_type y, typename tom::transf::InterpolNearestNeighbour<T>::idx_type z);

private:
    const T *data;
    std::size_t sizex, sizey, sizez;
    std::size_t stridey, stridez;

    T defaultval;
};

/***********************************************************************//**
 * Class for tom::transf::transform which performs trilinear interpolation.
 **************************************************************************/
template<typename T>
class InterpolTriLinear {

public:
    typedef float idx_type;
    InterpolTriLinear(T defaultval);
    void setvolume(const tom::Volume<T> &src);
    void setvolume(const T *data,const std::size_t sizex,const std::size_t sizey,const std::size_t sizez,const std::size_t stridey,const std::size_t stridez);
    T interpolate(const typename tom::transf::InterpolTriLinear<T>::idx_type &x, const typename tom::transf::InterpolTriLinear<T>::idx_type &y, const typename tom::transf::InterpolTriLinear<T>::idx_type &z);

private:
    const T *data;
    std::size_t sizex, sizey, sizez;
    typename tom::transf::InterpolTriLinear<T>::idx_type sizex_m1, sizey_m1, sizez_m1;
    std::size_t stridey, stridez;

    T defaultval;
};



/***********************************************************************//**
 * Class for tom::transf::transform which performs tricubic interpolation.
 **************************************************************************/
template<typename T>
class InterpolTriCubic {

public:
    typedef float idx_type;
    InterpolTriCubic(T defaultval);
    void setvolume(const tom::Volume<T> &src);
    void setvolume(const T *data,const std::size_t sizex,const std::size_t sizey,const std::size_t sizez,const std::size_t stridey,const std::size_t stridez);
    T interpolate(const typename tom::transf::InterpolTriCubic<T>::idx_type &x, const typename tom::transf::InterpolTriCubic<T>::idx_type &y, const typename tom::transf::InterpolTriCubic<T>::idx_type &z);

private:
    const T *data;
    std::size_t sizex, sizey, sizez;
    std::size_t stridey, stridez;

    T defaultval;
};


/***********************************************************************//**
 * \brief Initialise the interpolation object with the volume.
 *
 * If the interpolation does not support the stride parameters it should
 * throw tom::transf::memory_alignment_not_supported.
 **************************************************************************/
template<typename T>
void tom::transf::InterpolTriLinear<T>::setvolume(const tom::Volume<T> &src) {

    const std::size_t stridex = src.getStrideX();
    const std::size_t stridey = src.getStrideY();
    const std::size_t stridez = src.getStrideZ();

    if (stridex!=sizeof(T) || stridey%sizeof(T) || stridez%sizeof(T)) {
        throw tom::transf::memory_alignment_not_supported("The volume of InterpolTriLinear must be aligned with sizeof(T) and the fastest dimension must be contigous.");
    }

    this->sizex = src.getSizeX();
    this->sizey = src.getSizeY();
    this->sizez = src.getSizeZ();
    this->sizex_m1 = this->sizex - 1;
    this->sizey_m1 = this->sizey - 1;
    this->sizez_m1 = this->sizez - 1;
    this->stridey = stridey/sizeof(T);
    this->stridez = stridez/sizeof(T);
    this->data = &src.get();

}

template<typename T>
void tom::transf::InterpolTriLinear<T>::setvolume(const T *data,const std::size_t sizex,const std::size_t sizey,const std::size_t sizez,const std::size_t stridey,const std::size_t stridez) {

    this->sizex = sizex;
    this->sizey = sizey;
    this->sizez = sizez;
    this->stridey = stridey;
    this->stridez = stridez;
    this->data = data;
}


/***********************************************************************//**
* \brief Initialise the interpolation object with the volume.
*
* If the interpolation does not support the stride parameters it should
* throw tom::transf::memory_alignment_not_supported.
**************************************************************************/
template<typename T>
void tom::transf::InterpolTriCubic<T>::setvolume(const tom::Volume<T> &src) {

    const std::size_t stridex = src.getStrideX();
    const std::size_t stridey = src.getStrideY();
    const std::size_t stridez = src.getStrideZ();

    if (stridex!=sizeof(T) || stridey%sizeof(T) || stridez%sizeof(T)) {
        throw tom::transf::memory_alignment_not_supported("The volume of InterpolTriCubic must be aligned with sizeof(T) and the fastest dimension must be contigous.");
    }

    this->sizex = src.getSizeX();
    this->sizey = src.getSizeY();
    this->sizez = src.getSizeZ();
    this->stridey = stridey/sizeof(T);
    this->stridez = stridez/sizeof(T);
    this->data = &src.get();

}

template<typename T>
void tom::transf::InterpolTriCubic<T>::setvolume(const T *data,const std::size_t sizex,const std::size_t sizey,const std::size_t sizez,const std::size_t stridey,const std::size_t stridez) {

    this->sizex = sizex;
    this->sizey = sizey;
    this->sizez = sizez;
    this->stridey = stridey;
    this->stridez = stridez;
    this->data = data;
}

/***********************************************************************//**
* Class for tom::transf::transform which performs cubic spline interpolation.
**************************************************************************/
template<typename T>
class InterpolCubicSpline {

public:
	typedef float idx_type;
	InterpolCubicSpline(T defaultval);
	void setvolume(const tom::Volume<T> &src);
	void setvolume(const T *data,const std::size_t sizex,const std::size_t sizey,const std::size_t sizez,const std::size_t stridey,const std::size_t stridez);
	T interpolate(const typename tom::transf::InterpolCubicSpline<T>::idx_type &x, const typename tom::transf::InterpolCubicSpline<T>::idx_type &y, const typename tom::transf::InterpolCubicSpline<T>::idx_type &z);

private:
	const T *data;
	std::size_t sizex, sizey, sizez;
	std::size_t stridey, stridez;

	T defaultval;
};


/// Default type for interpolation.
template<typename T>
struct InterpolDefault {
public:
    typedef InterpolTriLinear<T> type;
};


/***********************************************************************//**
* Class for tom::transf::transform which performs fourier spline interpolation.
**************************************************************************/
template<typename T>
class InterpolFourierSpline {
    
public:
    typedef float idx_type;
    InterpolFourierSpline(T defaultval);
    void setvolume(const tom::Volume<T> &src);
    void setvolume(const T *data,const std::size_t sizex,const std::size_t sizey,const std::size_t sizez,const std::size_t stridey,const std::size_t stridez);
    T interpolate(typename tom::transf::InterpolFourierSpline<T>::idx_type x, typename tom::transf::InterpolFourierSpline<T>::idx_type y, typename tom::transf::InterpolFourierSpline<T>::idx_type z);
    
private:
    const T *data;
    std::size_t sizex, sizey, sizez;
    std::size_t stridey, stridez;
    
    T defaultval;
};
    



template<typename T, typename TINTERP>
void transform(const tom::Volume<T> &src, tom::Volume<T> &dst, const double *P, bool is_affinity, T valinf, TINTERP interp);

template<typename T>
void transformSpline(const tom::Volume<T> &src, tom::Volume<T> &v_rot, const tom::Volume<T> &mtx, T defaultval=0);

template<typename T, typename TINTERP>
void transformFourier(const tom::Volume<T> &src, tom::Volume<T> &dst, const double *P, bool is_affinity, T valinf, TINTERP interp,double shiftX,double shiftY,double shiftZ);

template<typename T>
void shiftFourier(const tom::Volume<T > &src, tom::Volume<T > &dst, double shift_x, double shift_y, double shift_z);

template<typename T> void rotate(const tom::Volume<T> &src, tom::Volume<T> &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, T defaultval=0);
template<typename T> void rotate(const tom::Volume<T> &src, tom::Volume<T> &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, double shiftx_pre, double shifty_pre, double shiftz_pre, double shiftx_post, double shifty_post, double shiftz_post, T defaultval=0);
template<typename T> void rotate(const tom::Volume<T> &src, tom::Volume<T> &v_rot, double phi, double psi, double theta, T defaultval=0);

template<typename T> void rotateCubic(const tom::Volume<T> &src, tom::Volume<T> &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, T defaultval=0);
template<typename T> void rotateCubic(const tom::Volume<T> &src, tom::Volume<T> &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, double shiftx_pre, double shifty_pre, double shiftz_pre, double shiftx_post, double shifty_post, double shiftz_post, T defaultval=0);
template<typename T> void rotateCubic(const tom::Volume<T> &src, tom::Volume<T> &v_rot, double phi, double psi, double theta, T defaultval=0);

template<typename T> void rotateSpline(const tom::Volume<T> &src, tom::Volume<T> &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, T defaultval=0);
template<typename T> void rotateSpline(const tom::Volume<T> &src, tom::Volume<T> &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, double shiftx_pre, double shifty_pre, double shiftz_pre, double shiftx_post, double shifty_post, double shiftz_post, T defaultval=0);
template<typename T> void rotateSpline(const tom::Volume<T> &src, tom::Volume<T> &v_rot, double phi, double psi, double theta, T defaultval=0);
    
template<typename T> void rotateFourierSpline(const tom::Volume<T > &src, tom::Volume<T > &v_rot, double phi, double psi, double theta, T defaultval=0);
template<typename T> void rotateFourierSpline(const tom::Volume<T > &src, tom::Volume<T > &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, T defaultval = 0);
template<typename T> void transformFourierSpline(const tom::Volume<T > &src, tom::Volume<T > &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, double shiftX, double shiftY, double shiftZ,T defaultval);
    
template<typename T> void rescale_fs(const tom::Volume<T               > &vsrc,                        tom::Volume<T               > &vdst,                        unsigned fftw_flag=FFTW_MEASURE);
template<typename T> void rescale_fs(const tom::Volume<T               > &vsrc,                        tom::Volume<std::complex<T> > &vdst, std::size_t sizez_dst, unsigned fftw_flag=FFTW_MEASURE);
template<typename T> void rescale_fs(const tom::Volume<std::complex<T> > &vsrc, std::size_t sizez_src, tom::Volume<T               > &vdst,                        unsigned fftw_flag=FFTW_MEASURE);
template<typename T> void rescale_fs(const tom::Volume<std::complex<T> > &vsrc, std::size_t sizez_src, tom::Volume<std::complex<T> > &vdst, std::size_t sizez_dst                                 );
template<typename T> void rescale_fs(const tom::Volume<T               > &vsrc,                        tom::Volume<T               > &vdst,                        double shiftx, double shifty, double shiftz, unsigned fftw_flag=FFTW_MEASURE);
template<typename T> void rescale_fs(const tom::Volume<T               > &vsrc,                        tom::Volume<std::complex<T> > &vdst, std::size_t sizez_dst, double shiftx, double shifty, double shiftz, unsigned fftw_flag=FFTW_MEASURE);
template<typename T> void rescale_fs(const tom::Volume<std::complex<T> > &vsrc, std::size_t sizez_src, tom::Volume<T               > &vdst,                        double shiftx, double shifty, double shiftz, unsigned fftw_flag=FFTW_MEASURE);
template<typename T> void rescale_fs(const tom::Volume<std::complex<T> > &vsrc, std::size_t sizez_src, tom::Volume<std::complex<T> > &vdst, std::size_t sizez_dst, double shiftx, double shifty, double shiftz                                 );

template<typename T> void rescale(const tom::Volume<T> &vsrc, tom::Volume<T> &vdst);
template<typename T, typename TINTERP> void rescale(const tom::Volume<T> &vsrc, tom::Volume<T> &vdst, TINTERP interp);

template<typename T> void rescale_rebin(const tom::Volume<T> &vsrc, tom::Volume<T> &vdst, std::size_t sizex, std::size_t sizey, std::size_t sizez, const T *defaultval, unsigned fftw_flag);
template<typename T> void rescale_rebin(const tom::Volume<T> &vsrc, std::size_t binsrcx, std::size_t binsrcy, std::size_t binsrcz,
                                              tom::Volume<T> &vdst, std::size_t bindstx, std::size_t bindsty, std::size_t bindstz,
                                              std::size_t sizex, std::size_t sizey, std::size_t sizez, const T *defaultval, unsigned fftw_flag);


template<typename T, typename TBIN> void bin(const tom::Volume<T> &vsrc, tom::Volume<TBIN> &vdst);
template<typename T, typename TBIN> void bin(const tom::Volume<T> &vsrc, tom::Volume<TBIN> &vdst, std::size_t binx);
template<typename T, typename TBIN> void bin(const tom::Volume<T> &vsrc, tom::Volume<TBIN> &vdst, std::size_t binx, std::size_t biny, std::size_t binz);


template<typename T1, typename T2> void paste(const tom::Volume<T1> &v1, tom::Volume<T2> &v2, const T2 *outside_value);
template<typename T1, typename T2> void paste(const tom::Volume<T1> &v1, tom::Volume<T2> &v2, signed long cornerx, signed long cornery, signed long cornerz, const T2 *outside_value);

} // namespace transf
} // namespace tom







// INLINE FUNCTIONS

/***********************************************************************//**
 * \brief Constructor of the nearest neighbour interpolaion.
 *
 * Constructor of the nearest neighbour interpolaion.
 **************************************************************************/
template<typename T>
inline tom::transf::InterpolNearestNeighbour<T>::InterpolNearestNeighbour(T defaultval_p):
    data(NULL),
    sizex(0), sizey(0), sizez(0), stridey(0), stridez(0),
    defaultval(defaultval_p) {
}


/***********************************************************************//**
 * \brief Interpolate a point in the volume
 *
 * Interpolate a point in the volume
 **************************************************************************/
template<typename T>
inline T tom::transf::InterpolNearestNeighbour<T>::interpolate(typename tom::transf::InterpolNearestNeighbour<T>::idx_type x, typename tom::transf::InterpolNearestNeighbour<T>::idx_type y, typename tom::transf::InterpolNearestNeighbour<T>::idx_type z) {

    x += 0.5;
    y += 0.5;
    z += 0.5;
    if (x<0 || y<0 || z<0 || x>=this->sizex || y>=this->sizey || z>=this->sizez) {
        return this->defaultval;
    }
    return this->data[static_cast<std::size_t>(z)*this->stridez +
                       static_cast<std::size_t>(y)*this->stridey +
                       static_cast<std::size_t>(x)];

}


/***********************************************************************//**
 * \brief Constructor of the trilinear interpolaion.
 *
 * Constructor of the trilinear interpolaion.
 **************************************************************************/
template<typename T>
inline tom::transf::InterpolTriLinear<T>::InterpolTriLinear(T defaultval_p):
    data(NULL),
    sizex(0), sizey(0), sizez(0), sizex_m1(0), sizey_m1(0), sizez_m1(0), stridey(0), stridez(0),
    defaultval(defaultval_p) {
}




/***********************************************************************//**
 * \brief Interpolate a point in the volume
 *
 * Interpolate a point in the volume
 **************************************************************************/
template<typename T>
inline T tom::transf::InterpolTriLinear<T>::interpolate(const typename tom::transf::InterpolTriLinear<T>::idx_type &x, const typename tom::transf::InterpolTriLinear<T>::idx_type &y, const typename tom::transf::InterpolTriLinear<T>::idx_type &z) {

    typedef typename tom::transf::InterpolTriLinear<T>::idx_type FLOATTYPE;

    if (x<0. || y<0. || z<0. ||
        x > (this->sizex_m1) ||
        y > (this->sizey_m1) ||
        z > (this->sizez_m1)) {
        return this->defaultval;
    }
    const size_t floorx = static_cast<std::size_t>(x);
    const size_t floory = static_cast<std::size_t>(y);
    const size_t floorz = static_cast<std::size_t>(z);
    const FLOATTYPE xoffseth = x - floorx;
    const FLOATTYPE yoffseth = y - floory;
    const FLOATTYPE zoffseth = z - floorz;

    const T *src = &this->data[floorz * this->stridez + floory*this->stridey + floorx];
	//printf("Coordinates %f %f %f \n",x,y,z);
    switch (static_cast<int>(xoffseth > 0.) | (static_cast<int>(yoffseth > 0.)<<1) | (static_cast<int>(zoffseth > 0.)<<2)) {
        case 0x07: {
                FLOATTYPE x00x, x01x, x10x, x11x, x0yx;
                x00x = src[0] + (src[1]-src[0])*xoffseth;
                src += this->stridey;
                x01x = src[0] + (src[1]-src[0])*xoffseth;
                src += this->stridez;
                x11x = src[0] + (src[1]-src[0])*xoffseth;
                src -= this->stridey;
                x10x = src[0] + (src[1]-src[0])*xoffseth;
                x0yx = x00x + (x01x - x00x) * yoffseth;
                return x0yx + ((x10x + (x11x - x10x) * yoffseth) - x0yx) * zoffseth;
            }
        case 0x06: {
                const FLOATTYPE x0y0 = src[0] + (src[this->stridey]-src[0])*yoffseth;
                src += this->stridez;
                return x0y0 + (src[0] + (src[this->stridey]-src[0])*yoffseth - x0y0)*zoffseth;
            }
        case 0x05: {
                const FLOATTYPE x00x = src[0] + (src[1]-src[0])*xoffseth;
                src += this->stridez;
                return x00x + (src[0] + (src[1]-src[0])*xoffseth - x00x)*zoffseth;
            }
        case 0x03: {
                const FLOATTYPE x00x = src[0] + (src[1]-src[0])*xoffseth;
                src += this->stridey;
                return x00x + (src[0] + (src[1]-src[0])*xoffseth - x00x)* yoffseth;
            }
        case 0x04:
            return src[0] + (src[this->stridez]-src[0])*zoffseth;
        case 0x02:
            return src[0] + (src[this->stridey]-src[0])*yoffseth;
        case 0x01:
            return src[0] + (src[1]            -src[0])*xoffseth;
    }
    return src[0];
}


/***********************************************************************//**
 * \brief Constructor of the tricubic interpolaion.
 *
 * Constructor of the tricubic interpolaion.
 **************************************************************************/
template<typename T>
inline tom::transf::InterpolTriCubic<T>::InterpolTriCubic(T defaultval_p):
    data(NULL),
    sizex(0), sizey(0), sizez(0), stridey(0), stridez(0),
    defaultval(defaultval_p) {
}


/*macro for cubic interpolation*/
/*calculates the polynomial value P_i(x)*y_i according to */
/*http://mathworld.wolfram.com/LagrangeInterpolatingPolynomial.html*/
//warning : arguments with operations must be in (), undefined behaviour otherwise
//x : float position at interpolation point
//y : the value
//xj : current voxel
//x2 : the other coordinate
//x3 : the other coordinate
#define CUB_INT(x,y,xj,x2,x3,val){\
val = y*(x-x2)/(xj-x2)*(x-x3)/(xj-x3);\
}

/***********************************************************************//**
 * \brief Interpolate a point in the volume
 *
 * Interpolate a point in the volume
 **************************************************************************/
template<typename T>
T tom::transf::InterpolTriCubic<T>::interpolate(const typename tom::transf::InterpolTriCubic<T>::idx_type &x, const typename tom::transf::InterpolTriCubic<T>::idx_type &y, const typename tom::transf::InterpolTriCubic<T>::idx_type &z) {

    typedef typename tom::transf::InterpolTriCubic<T>::idx_type FLOATTYPE;

    bool is3D = this->sizez > 1;

	// is the current position in the data or outside. return default value if outside
    if (is3D &&
    	(x<1. || y<1. || z<1. ||
        x > (this->sizex)-2 ||
        y > (this->sizey)-2 ||
        z > (this->sizez)-2) ) {
    		InterpolTriLinear<T> linearInterpolation = InterpolTriLinear<T>(this->defaultval);
    		linearInterpolation.setvolume(this->data,this->sizex,this->sizey,this->sizez,this->stridey,this->stridez);
    		return linearInterpolation.interpolate(x,y,z);
    }

	if(!is3D &&
    	(x<1. || y<1. ||
        x > (this->sizex)-2 ||
        y > (this->sizey)-2) ){
			InterpolTriLinear<T> linearInterpolation = InterpolTriLinear<T>(this->defaultval);
			linearInterpolation.setvolume(this->data,this->sizex,this->sizey,this->sizez,this->stridey,this->stridez);
			return linearInterpolation.interpolate(x,y,z);
	}

    const long floorx = static_cast<long>(x);//integer position along x
    const long floory = static_cast<long>(y);//integer position along y
    const long floorz = static_cast<long>(z);//integer position along z
    const FLOATTYPE xoffseth = x - floorx;//floating point offset along x
    const FLOATTYPE yoffseth = y - floory;//floating point offset along y
    const FLOATTYPE zoffseth = z - floorz;//floating point offset along z

    const T *src = &this->data[floorz * this->stridez + floory*this->stridey + floorx];//point into the closest position to interpolation
	//printf("Coordinates %f %f %f \n",x,y,z);

    switch (static_cast<int>(xoffseth > 0.) | (static_cast<int>(yoffseth > 0.)<<1) | (static_cast<int>(zoffseth > 0.)<<2)) {
		case 0x00:{
			//all match grid, no interpolation
			//printf("0x00 \n");
			return src[0];
		}
		default: {

			/*
			interpolation square
				p1  p2  p3
				p4  p5  p6
					  P
				p7  p8  p9
			 */

			const long px_pl_1 = floorx + 1; //all voxels plus 1 from x -> p3,p6,p9
			const long px_mi_1 = floorx - 1; //all voxels minus 1 from x -> p1,p4,p7

			const long py_pl_1 = floory + 1; //all voxels plus 1 from y -> p1,p2,p3
			const long py_mi_1 = floory - 1; //all voxels minus 1 from x -> p7,p8,p9

			//those will hold the current voxel values
			T v1 = this->defaultval;
			T v2 = this->defaultval;
			T v3 = this->defaultval;
			T v4 = this->defaultval;
			T v5 = this->defaultval;
			T v6 = this->defaultval;
			T v7 = this->defaultval;
			T v8 = this->defaultval;
			T v9 = this->defaultval;

			T line1 = this->defaultval; //interpolation values for each line (y)
			T line2 = this->defaultval;
			T line3 = this->defaultval;

			int lowerLayerBound = -1;
			int upperLayerBound = 1;
			int layerOffset = 1;
			T layerValues[3];

			if(!is3D){
				lowerLayerBound = 0;
				upperLayerBound = 0;
				layerOffset = 0;
			}

			//interpolation values for each layer (z)
			for (int zIteration=lowerLayerBound; zIteration <= upperLayerBound; zIteration++) {
				//current position in memory plus current z layer offset in voxels (of type T)
				//first will be negative (-1), second 0 (same layer), third is 1, next layer

				//load the pixel values
				v1 = *(src-1-this->stridey + zIteration*this->stridez); //one line up in y direction, one position back in x
				v2 = *(src  -this->stridey + zIteration*this->stridez); //one line up in y direction
				v3 = *(src+1-this->stridey + zIteration*this->stridez); //one line up in y direction, one position forward in x
				v4 = *(src-1 + zIteration*this->stridez); //same line in y
				v5 = *(src + zIteration*this->stridez); //...
				v6 = *(src+1 + zIteration*this->stridez);
				v7 = *(src-1+this->stridey + zIteration*this->stridez);
				v8 = *(src  +this->stridey + zIteration*this->stridez);
				v9 = *(src+1+this->stridey + zIteration*this->stridez);


				//printf("Value %f %f %f %f %f %f %f %f %f \n",v1,v2,v3,v4,v5,v6,v7,v8,v9);

				//interpolate first row 1 2 3
				CUB_INT(x,v1,px_mi_1,floorx,px_pl_1,line1); //px1,px2,px3
				CUB_INT(x,v2,floorx,px_pl_1,px_mi_1,line2); //px2,px3,px1
				CUB_INT(x,v3,px_pl_1,px_mi_1,floorx,line3); //px3,px1,px2
				//store values into v1
				//printf("Line 1 %f %f %f\n",line1,line2,line3);
				v1 = line1+line2+line3;


				//same for the next rows
				CUB_INT(x,v4,px_mi_1,floorx,px_pl_1,line1); //px4,px5,px6
				CUB_INT(x,v5,floorx,px_pl_1,px_mi_1,line2); //px5,px6,px4
				CUB_INT(x,v6,px_pl_1,px_mi_1,floorx,line3); //px6,px5,px4
				//printf("Line 2 %f %f %f\n",line1,line2,line3);
				v2 = line1+line2+line3;

				CUB_INT(x,v7,px_mi_1,floorx,px_pl_1,line1); //px7,px8,px9
				CUB_INT(x,v8,floorx,px_pl_1,px_mi_1,line2); //px8,px9,px7
				CUB_INT(x,v9,px_pl_1,px_mi_1,floorx,line3); //px9,px7,px8
				//printf("Line 3 %f %f %f\n",line1,line2,line3);
				v3 = line1+line2+line3;

				//interpolate col 2 5 8 in y direction
				CUB_INT(y,v1,py_mi_1,floory,py_pl_1,line1); //py2,py5,py8
				CUB_INT(y,v2,floory,py_pl_1,py_mi_1,line2); //py5,py8,py2
				CUB_INT(y,v3,py_pl_1,py_mi_1,floory,line3); //py8,py2,py5
				//printf("Row 1%f %f %f\n",line1,line2,line3);

				layerValues[zIteration + layerOffset] = line1+line2+line3;
			}
            //printf("Layer Values %f %f %f \n",layerValues[0],layerValues[1],layerValues[2]);
			//printf("FloorZ %d %d %d \n",floorz-1,floorz,floorz +1);

			if(is3D){
				CUB_INT(z,layerValues[0],(floorz-1),floorz,(floorz+1),line1);
				CUB_INT(z,layerValues[1],floorz,(floorz+1),(floorz-1),line2);
				CUB_INT(z,layerValues[2],(floorz+1),(floorz-1),floorz,line3);
			}
		    //printf("Layer 1 %f %f %f\n",line1,line2,line3);

			return (line1 + line2 + line3);

		}


		return this->defaultval;
    }
    return this->defaultval;
}


/***********************************************************************//**
* \brief Constructor of the trilinear interpolaion.
*
* Constructor of the trilinear interpolaion.
**************************************************************************/
template<typename T>
inline tom::transf::InterpolCubicSpline<T>::InterpolCubicSpline(T defaultval_p):
data(NULL),
sizex(0), sizey(0), sizez(0), stridey(0), stridez(0),
defaultval(defaultval_p) {
}


/*prestep for cspline interpolation for the calculation of D1 and D2 according to source*/
/*http://mathworld.wolfram.com/CubicSpline.html*/
/*macro for spline interpolation in a row*/
/*3*(f1*(v2-v1)+f2*(v3-v1)+f3*(v4-v2)+f4*(v4-v3))*/
#define CSPL_INT(f1,f2,f3,f4,a,b,c,d,val){\
val = ((T)3)*(f1*(b-a)+f2*(c-a)+f3*(d-b)+f4*(d-c));\
}

/*macro to calculate the result of cspline interpolation using 4 features*/
/*c2 and c3 are calculated by CSPL_INT*/
/*value = c2+D1*(r_y-py6)+(3*(c3-c2)-2*D1-D2)*(r_y-py6)*(r_y-py6)+(2*(c2-c3)+D1+D2)*(r_y-py6)*(r_y-py6)*(r_y-py6);*/
/*c2~c2, c3~c3,D1~D1,D2~D2*,off~(r_y-py6) - distance from nearest point/*/
#define CSPL_CALC(c2,c3,D1,D2,off,val){\
val = c2+D1*off+(((T)3)*(c3-c2)-((T)2)*D1-D2)*off*off+(((T)2)*(c2-c3)+D1+D2)*off*off*off;\
}

/***********************************************************************//**
* \brief Interpolate a point in the volume
*
* Interpolate a point in the volume
**************************************************************************/
template<typename T>
T tom::transf::InterpolCubicSpline<T>::interpolate(const typename tom::transf::InterpolCubicSpline<T>::idx_type &x, const typename tom::transf::InterpolCubicSpline<T>::idx_type &y, const typename tom::transf::InterpolCubicSpline<T>::idx_type &z) {

    typedef typename tom::transf::InterpolCubicSpline<T>::idx_type FLOATTYPE;

    bool is3D = this->sizez > 1;

	// is the current position in the data or outside. return default value if outside
    if (is3D &&
    	(x<2. || y<2. || z<2. ||
        x > (this->sizex)-3 ||
        y > (this->sizey)-3 ||
        z > (this->sizez)-3) ) {
    		InterpolTriCubic<T> cubicInterpolation = InterpolTriCubic<T>(this->defaultval);
    		cubicInterpolation.setvolume(this->data,this->sizex,this->sizey,this->sizez,this->stridey,this->stridez);
    	    return cubicInterpolation.interpolate(x,y,z);
    }

	if(!is3D &&
    	(x<2. || y<2. ||
        x > (this->sizex)-3 ||
        y > (this->sizey)-3) ){
			InterpolTriCubic<T> cubicInterpolation = InterpolTriCubic<T>(this->defaultval);
			cubicInterpolation.setvolume(this->data,this->sizex,this->sizey,this->sizez,this->stridey,this->stridez);
		    return cubicInterpolation.interpolate(x,y,z);
	}

    const long floorx = static_cast<long>(x);//integer position along x
    const long floory = static_cast<long>(y);//integer position along y
    const long floorz = static_cast<long>(z);//integer position along z
    const FLOATTYPE xoffseth = x - floorx;//floating point offset along x
    const FLOATTYPE yoffseth = y - floory;//floating point offset along y
    const FLOATTYPE zoffseth = z - floorz;//floating point offset along z

    const T *src = &this->data[floorz * this->stridez + floory*this->stridey + floorx];//point into the closest position to interpolation
	//printf("Coordinates %f %f %f \n",x,y,z);

    switch (static_cast<int>(xoffseth > 0.) | (static_cast<int>(yoffseth > 0.)<<1) | (static_cast<int>(zoffseth > 0.)<<2)) {
		case 0x00:{
			//all match grid, no interpolation
			return src[0];
		}
		default:{
			/*interpolation square

			 p1   p2   p3 p4
			 p5   p6   p7 p8
			         P
			 p9   p10  p11 p12
			 p13  p14  p15 p16

			 */
			T f1= -0.1556;
			T f2=0.3111;
			T f3=-0.0889;
			T f4=0.0444;

			//those will hold the current voxel values
			T v1 = this->defaultval;
			T v2 = this->defaultval;
			T v3 = this->defaultval;
			T v4 = this->defaultval;
			T v5 = this->defaultval;
			T v6 = this->defaultval;
			T v7 = this->defaultval;
			T v8 = this->defaultval;
			T v9 = this->defaultval;
			T v10 = this->defaultval;
			T v11 = this->defaultval;
			T v12 = this->defaultval;
			T v13 = this->defaultval;
			T v14 = this->defaultval;
			T v15 = this->defaultval;
			T v16 = this->defaultval;

			T line1 = this->defaultval; //interpolation values for each line (y)
			T line2 = this->defaultval;
			T line3 = this->defaultval;
			T line4 = this->defaultval;

			T D1 = this->defaultval;
			T D2 = this->defaultval;

			int lowerLayerBound = -1;
			int upperLayerBound = 2;
			int layerOffset = 1;
			T layerValues[4];

			if(!is3D){
				lowerLayerBound = 0;
				upperLayerBound = 0;
				layerOffset = 0;
			}
			//interpolation values for each layer (z)
			for (int zIteration=lowerLayerBound; zIteration <= upperLayerBound; zIteration++) {

				//load the pixel values
				v1 = *(src-1-this->stridey + zIteration*this->stridez); //one line up in y direction, one position back in x
				v2 = *(src  -this->stridey + zIteration*this->stridez); //one line up in y direction
				v3 = *(src+1-this->stridey + zIteration*this->stridez); //one line up in y direction, one position forward in x
				v4 = *(src+2-this->stridey + zIteration*this->stridez);

				v5 = *(src-1 + zIteration*this->stridez); //same line in y
				v6 = *(src + zIteration*this->stridez); //...
				v7 = *(src+1 + zIteration*this->stridez);
				v8 = *(src+2 + zIteration*this->stridez);

				v9 = *(src-1+this->stridey + zIteration*this->stridez);
				v10 = *(src  +this->stridey + zIteration*this->stridez);
				v11 = *(src+1+this->stridey + zIteration*this->stridez);
				v12 = *(src+2+this->stridey + zIteration*this->stridez);

				v13 = *(src-1+2*this->stridey + zIteration*this->stridez);
				v14 = *(src  +2*this->stridey + zIteration*this->stridez);
				v15 = *(src+1+2*this->stridey + zIteration*this->stridez);
				v16 = *(src+2+2*this->stridey + zIteration*this->stridez);

				/*calculate spline value for line 1 2 3 4 above pixel of interest */
                CSPL_INT(f1,f2,f3,f4,v1,v2,v3,v4,D1);
                CSPL_INT(f4,f3,f2,f1,v1,v2,v3,v4,D2);
                CSPL_CALC(v2,v3,D1,D2,xoffseth,line1);

				/*calculate spline value for line 5 6 7 8 above pixel of interest */
                CSPL_INT(f1,f2,f3,f4,v5,v6,v7,v8,D1);
                CSPL_INT(f4,f3,f2,f1,v5,v6,v7,v8,D2);
                CSPL_CALC(v6,v7,D1,D2,xoffseth,line2);

				/*calculate spline value for line 9 10 11 12 above pixel of interest */
                CSPL_INT(f1,f2,f3,f4,v9,v10,v11,v12,D1);
                CSPL_INT(f4,f3,f2,f1,v9,v10,v11,v12,D2);
                CSPL_CALC(v10,v11,D1,D2,xoffseth,line3);

				/*calculate spline value for line 13 14 15 16 above pixel of interest */
                CSPL_INT(f1,f2,f3,f4,v13,v14,v15,v16,D1);
                CSPL_INT(f4,f3,f2,f1,v13,v14,v15,v16,D2);
                CSPL_CALC(v14,v15,D1,D2,xoffseth,line4);

				/*finaly, calculate spline into y direction and save into value[z]*/
                CSPL_INT(f1,f2,f3,f4,line1,line2,line3,line4,D1);
                CSPL_INT(f4,f3,f2,f1,line1,line2,line3,line4,D2);
                CSPL_CALC(line2,line3,D1,D2,yoffseth,layerValues[zIteration + layerOffset]);

			}

			if (is3D) {
				/*calculate spline value for z direction*/
				CSPL_INT(f1,f2,f3,f4,layerValues[0],layerValues[1],layerValues[2],layerValues[3],D1);
				CSPL_INT(f4,f3,f2,f1,layerValues[0],layerValues[1],layerValues[2],layerValues[3],D2);

				CSPL_CALC(layerValues[1],layerValues[2],D1,D2,zoffseth,D1);

				return D1;
			}else {
				return layerValues[0];
			}
		}
	}
    return this->defaultval;
}


/***********************************************************************
* \brief Constructor of the nearest neighbour interpolaion.
*
* Constructor of the nearest neighbour interpolaion.
**************************************************************************/
template<typename T>
inline tom::transf::InterpolFourierSpline<T>::InterpolFourierSpline(T defaultval_p):
data(NULL),
sizex(0), sizey(0), sizez(0), stridey(0), stridez(0),
defaultval(defaultval_p) {
}


/***********************************************************************
* \brief Interpolate a point in the volume
*
* Interpolate a point in the volume
**************************************************************************/
template<typename T>
inline T tom::transf::InterpolFourierSpline<T>::interpolate(typename tom::transf::InterpolFourierSpline<T>::idx_type x, typename tom::transf::InterpolFourierSpline<T>::idx_type y, typename tom::transf::InterpolFourierSpline<T>::idx_type z) {
    
    typedef typename tom::transf::InterpolFourierSpline<T>::idx_type FLOATTYPE;
    
    bool is3D = this->sizez > 1;

	// is the current position in the data or outside. return default value if outside
    if (is3D &&
    	(x<2. || y<2. || z<2. ||
          x > (this->sizex)-3 ||
          y > (this->sizey)-3 ||
          z > (this->sizez)-3) ) {
    	    return (T)0;
        }
    const long floorx = static_cast<long>(x);//integer position along x
    const long floory = static_cast<long>(y);//integer position along y
    const long floorz = static_cast<long>(z);//integer position along z
    const FLOATTYPE xoffseth = x - floorx;//floating point offset along x
    const FLOATTYPE yoffseth = y - floory;//floating point offset along y
    const FLOATTYPE zoffseth = z - floorz;//floating point offset along z
    
    const T *src = &this->data[floorz * this->stridez + floory*this->stridey + floorx];//point into the closest position to interpolation
	//printf("Coordinates %f %f %f \n",x,y,z);
    
    switch (static_cast<int>(xoffseth > 0.) | (static_cast<int>(yoffseth > 0.)<<1) | (static_cast<int>(zoffseth > 0.)<<2)) {
		case 0x00:{
			//all match grid, no interpolation
			//printf("0x00 \n");
			return src[0];
		}
		default:{
            T f1= -0.1556;
            T f2=0.3111;
            T f3=-0.0889;
            T f4=0.0444;
            
            //those will hold the current voxel values
            T v1 = this->defaultval;
            T v2 = this->defaultval;
            T v3 = this->defaultval;
            T v4 = this->defaultval;
            T v5 = this->defaultval;
            T v6 = this->defaultval;
            T v7 = this->defaultval;
            T v8 = this->defaultval;
            T v9 = this->defaultval;
            T v10 = this->defaultval;
            T v11 = this->defaultval;
            T v12 = this->defaultval;
            T v13 = this->defaultval;
            T v14 = this->defaultval;
            T v15 = this->defaultval;
            T v16 = this->defaultval;
            
            T line1 = this->defaultval; //interpolation values for each line (y)
            T line2 = this->defaultval;
            T line3 = this->defaultval;
            T line4 = this->defaultval;
            
            T D1 = this->defaultval;
            T D2 = this->defaultval;
            
            int lowerLayerBound = -1;
            int upperLayerBound = 2;
            int layerOffset = 1;
            T layerValues[4];
            
            if(!is3D){
                lowerLayerBound = 0;
                upperLayerBound = 0;
                layerOffset = 0;
            }
            //interpolation values for each layer (z)
            for (int zIteration=lowerLayerBound; zIteration <= upperLayerBound; zIteration++) {
                
                //load the pixel values
                v1 = *(src-1-this->stridey + zIteration*this->stridez); //one line up in y direction, one position back in x
                v2 = *(src  -this->stridey + zIteration*this->stridez); //one line up in y direction
                v3 = *(src+1-this->stridey + zIteration*this->stridez); //one line up in y direction, one position forward in x
                v4 = *(src+2-this->stridey + zIteration*this->stridez);
                
                v5 = *(src-1 + zIteration*this->stridez); //same line in y
                v6 = *(src + zIteration*this->stridez); //...
                v7 = *(src+1 + zIteration*this->stridez);
                v8 = *(src+2 + zIteration*this->stridez);
                
                v9 = *(src-1+this->stridey + zIteration*this->stridez);
                v10 = *(src  +this->stridey + zIteration*this->stridez);
                v11 = *(src+1+this->stridey + zIteration*this->stridez);
                v12 = *(src+2+this->stridey + zIteration*this->stridez);
                
                v13 = *(src-1+2*this->stridey + zIteration*this->stridez);
                v14 = *(src  +2*this->stridey + zIteration*this->stridez);
                v15 = *(src+1+2*this->stridey + zIteration*this->stridez);
                v16 = *(src+2+2*this->stridey + zIteration*this->stridez);
                
                /*calculate spline value for line 1 2 3 4 above pixel of interest */
                CSPL_INT(f1,f2,f3,f4,v1,v2,v3,v4,D1);
                CSPL_INT(f4,f3,f2,f1,v1,v2,v3,v4,D2);
                CSPL_CALC(v2,v3,D1,D2,xoffseth,line1);
                
                /*calculate spline value for line 5 6 7 8 above pixel of interest */
                CSPL_INT(f1,f2,f3,f4,v5,v6,v7,v8,D1);
                CSPL_INT(f4,f3,f2,f1,v5,v6,v7,v8,D2);
                CSPL_CALC(v6,v7,D1,D2,xoffseth,line2);
                
                /*calculate spline value for line 9 10 11 12 above pixel of interest */
                CSPL_INT(f1,f2,f3,f4,v9,v10,v11,v12,D1);
                CSPL_INT(f4,f3,f2,f1,v9,v10,v11,v12,D2);
                CSPL_CALC(v10,v11,D1,D2,xoffseth,line3);
                
                /*calculate spline value for line 13 14 15 16 above pixel of interest */
                CSPL_INT(f1,f2,f3,f4,v13,v14,v15,v16,D1);
                CSPL_INT(f4,f3,f2,f1,v13,v14,v15,v16,D2);
                CSPL_CALC(v14,v15,D1,D2,xoffseth,line4);
                
                /*finaly, calculate spline into y direction and save into value[z]*/
                CSPL_INT(f1,f2,f3,f4,line1,line2,line3,line4,D1);
                CSPL_INT(f4,f3,f2,f1,line1,line2,line3,line4,D2);
                CSPL_CALC(line2,line3,D1,D2,yoffseth,layerValues[zIteration + layerOffset]);
                
            }
            
            if (is3D) {
                /*calculate spline value for z direction*/
                CSPL_INT(f1,f2,f3,f4,layerValues[0],layerValues[1],layerValues[2],layerValues[3],D1);
                CSPL_INT(f4,f3,f2,f1,layerValues[0],layerValues[1],layerValues[2],layerValues[3],D2);
                
                CSPL_CALC(layerValues[1],layerValues[2],D1,D2,zoffseth,D1);
                
                return D1;
            }else {
                return layerValues[0];
            }
        }
    }
    
    

    
    
    
}


/****************************************************************************//**
 * \brief Rotates the volume around the its center.
 *
 * Beware in case of volume with even size, the center is rounded towards to
 * next integer (size/2)
 *******************************************************************************/
template <typename T>
inline void tom::transf::rotate(const tom::Volume<T> &src, tom::Volume<T> &v_rot, double phi, double psi, double theta, T defaultval) {
    tom::transf::rotate(src, v_rot, phi, psi, theta, src.getSizeX()/2, src.getSizeY()/2, src.getSizeZ()/2, defaultval);
}


/****************************************************************************//**
 * \brief Rotates the volume.
 *******************************************************************************/
template<typename T>
inline void tom::transf::rotate(const tom::Volume<T> &src, tom::Volume<T> &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, T defaultval) {
    tom::transf::rotate(src, v_rot, phi, psi, theta, centerx, centery, centerz, 0,0,0, 0,0,0, defaultval);
}

/****************************************************************************//**
* \brief Rotates the volume around the its center. Uses TriCubic interpolation!
*
* Beware in case of volume with even size, the center is rounded towards to
* next integer (size/2)
*******************************************************************************/
template <typename T>
inline void tom::transf::rotateCubic(const tom::Volume<T> &src, tom::Volume<T> &v_rot, double phi, double psi, double theta, T defaultval) {
    tom::transf::rotateCubic(src, v_rot, phi, psi, theta, src.getSizeX()/2, src.getSizeY()/2, src.getSizeZ()/2, defaultval);
}


/****************************************************************************//**
* \brief Rotates the volume. Uses TriCubic interpolation!
*******************************************************************************/
template<typename T>
inline void tom::transf::rotateCubic(const tom::Volume<T> &src, tom::Volume<T> &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, T defaultval) {
	tom::transf::rotateCubic(src, v_rot, phi, psi, theta, centerx, centery, centerz, 0,0,0, 0,0,0, defaultval);
}


/****************************************************************************
* \brief Rotates the volume around the its center. Uses CubicSpline interpolation!
*
* Beware in case of volume with even size, the center is rounded towards to
* next integer (size/2)
*******************************************************************************/
template <typename T>
inline void tom::transf::rotateSpline(const tom::Volume<T> &src, tom::Volume<T> &v_rot, double phi, double psi, double theta, T defaultval) {
    tom::transf::rotateSpline(src, v_rot, phi, psi, theta, src.getSizeX()/2, src.getSizeY()/2, src.getSizeZ()/2, defaultval);
}


/****************************************************************************//**
* \brief Rotates the volume. Uses CubicSpline interpolation!
*******************************************************************************/
template<typename T>
inline void tom::transf::rotateSpline(const tom::Volume<T> &src, tom::Volume<T> &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, T defaultval) {
	tom::transf::rotateSpline(src, v_rot, phi, psi, theta, centerx, centery, centerz, 0,0,0, 0,0,0, defaultval);
}

/****************************************************************************
 * \brief Rotates the volume around the its center. Uses Fourier Spline interpolation!
 *
 * Beware in case of volume with even size, the center is rounded towards to
 * next integer (size/2)
 *******************************************************************************/
template <typename T>
inline void tom::transf::rotateFourierSpline(const tom::Volume<T > &src, tom::Volume<T > &v_rot, double phi, double psi, double theta, T defaultval) {
    tom::transf::rotateFourierSpline(src, v_rot, phi, psi, theta, src.getSizeX()/2, src.getSizeY()/2, src.getSizeZ()/2, defaultval);
}



/****************************************************************************//**
 * \brief Rescale a volume in time domain.
 *
 * \param[in] vsrc The source volume in Fourier space.
 * \param[in] vdst The resulting volume in time domain.
 *******************************************************************************/
template<typename T>
inline void tom::transf::rescale(const tom::Volume<T> &vsrc, tom::Volume<T> &vdst) {
    tom::transf::rescale<T, typename InterpolDefault<T>::type>(vsrc, vdst, typename InterpolDefault<T>::type(0));
}

/****************************************************************************//**
* \brief Rescale a volume in time domain.Interpolates tricubic.
*
* \param[in] vsrc The source volume in Fourier space.
* \param[in] vdst The resulting volume in time domain.
******************************************************************************
template<typename T>
inline void tom::transf::rescaleCubic(const tom::Volume<T> &vsrc, tom::Volume<T> &vdst) {
    tom::transf::rescale<T, typename InterpolTriCubic<T> >(vsrc, vdst, typename InterpolTriCubic<T>(0));
}*/
/****************************************************************************//**
 * \brief Rescale a volume in time domain.
 *
 * \param[in] vsrc The source volume in Fourier space.
 * \param[in] vdst The resulting volume in time domain.
 * \param[in] interp Interpolation functor.
 *******************************************************************************/
template<typename T, typename TINTERP>
inline void tom::transf::rescale(const tom::Volume<T> &vsrc, tom::Volume<T> &vdst, TINTERP interp) {
    const double P[16] = {  static_cast<double>(vdst.getSizeX())/(vsrc.getSizeX()), 0., 0., 0. ,
                            0., static_cast<double>(vdst.getSizeY())/(vsrc.getSizeY()), 0., 0. ,
                            0., 0., static_cast<double>(vdst.getSizeZ())/(vsrc.getSizeZ()), 0. ,
                            0., 0., 0., 1. };
    tom::transf::transform<T, TINTERP>(vsrc, vdst, P, true, 0, interp);
}




/****************************************************************************//**
 * \brief Binning of a volume.
 *
 * The binning factor is deduced by the size of the ouput volume.
 *******************************************************************************/
template<typename T, typename TBIN>
inline void tom::transf::bin(const tom::Volume<T> &vsrc, tom::Volume<TBIN> &vdst) {
    bin<T, TBIN>(vsrc, vdst, vsrc.getSizeX()/vdst.getSizeX(), vsrc.getSizeY()/vdst.getSizeY(), vsrc.getSizeZ()/vdst.getSizeZ());
}

/****************************************************************************//**
 * \brief Binning of a volume.
 *
 * The binning factor is deduced by the size of the ouput volume.
 *******************************************************************************/
template<typename T, typename TBIN>
inline void tom::transf::bin(const tom::Volume<T> &vsrc, tom::Volume<TBIN> &vdst, std::size_t binf) {
    bin<T, TBIN>(vsrc, vdst, binf, binf, binf);
}



template<typename T>
inline void tom::transf::rescale_fs(const tom::Volume<T> &vsrc, tom::Volume<T> &vdst, unsigned fftw_flag) {
    rescale_fs(vsrc, vdst, 0.,0.,0., fftw_flag);
}


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
inline void tom::transf::rescale_fs(const tom::Volume<T> &vsrc, tom::Volume<std::complex<T> > &vdst, std::size_t sizez_dst, unsigned fftw_flag) {
    rescale_fs(vsrc, vdst, sizez_dst, 0.,0.,0., fftw_flag);
}


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
inline void tom::transf::rescale_fs(const tom::Volume<std::complex<T> > &vsrc, std::size_t sizez_src, tom::Volume<T> &vdst, unsigned fftw_flag) {
    rescale_fs(vsrc, sizez_src, vdst, 0.,0.,0., fftw_flag);
}


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
inline void tom::transf::rescale_fs(const tom::Volume<std::complex<T> > &vsrc, std::size_t sizez_src, tom::Volume<std::complex<T> > &vdst, std::size_t sizez_dst) {
    rescale_fs(vsrc, sizez_src, vdst, sizez_dst, 0.,0.,0.);
}



/****************************************************************************//**
 * \brief Wrappes rescale_rebin where the binning factor is deduced.
 *******************************************************************************/
template<typename T>
inline void tom::transf::rescale_rebin(const tom::Volume<T> &vsrc, tom::Volume<T> &vdst, std::size_t sizex, std::size_t sizey, std::size_t sizez, const T *defaultval, unsigned fftw_flag) {
    rescale_rebin(  vsrc, sizex/vsrc.getSizeX(), sizey/vsrc.getSizeY(), sizez/vsrc.getSizeZ(),
                    vdst, sizex/vdst.getSizeX(), sizey/vdst.getSizeY(), sizez/vdst.getSizeZ(), sizex, sizey, sizez, defaultval, fftw_flag);
}



/****************************************************************************//**
 * \brief Wrapper for tom::transf::paste which maps the centers to each other.
 *******************************************************************************/
template<typename T1, typename T2>
inline void tom::transf::paste(const tom::Volume<T1> &v1, tom::Volume<T2> &v2, const T2 *outside_value) {
    paste(v1, v2,   static_cast<signed long>(v2.getSizeX()/2) - static_cast<signed long>(v1.getSizeX()/2),
                    static_cast<signed long>(v2.getSizeY()/2) - static_cast<signed long>(v1.getSizeY()/2),
                    static_cast<signed long>(v2.getSizeZ()/2) - static_cast<signed long>(v1.getSizeZ()/2), outside_value);
}



#endif



