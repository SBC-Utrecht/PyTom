/***********************************************************************//**
 * \file transform.cpp
 * \brief Implementations of transform.hpp
 * \author  Thomas Haller
 * \version 0.2
 * \date    18.12.2007
 *
 * Contains the implementations of transform.hpp
 **************************************************************************/
#include <tom/transf/transform.hpp>
#include <tom/tools/Precision.hpp>


#include <algorithm>

#include <tom/volume_fcn.hpp>
#include <tom/fftw/fftw_plan.hpp>

#define PI 3.141592653589793238

namespace {
template<typename T> inline bool is_simple_interpolationfcn                                 () { return false; }
template<>           inline bool is_simple_interpolationfcn<tom::transf::InterpolTriLinear<float > >() { return true;  }
template<>           inline bool is_simple_interpolationfcn<tom::transf::InterpolTriLinear<double> >() { return true;  }
template<>           inline bool is_simple_interpolationfcn<tom::transf::InterpolNearestNeighbour<float > >() { return true;  }
template<>           inline bool is_simple_interpolationfcn<tom::transf::InterpolNearestNeighbour<double> >() { return true;  }
template<>           inline bool is_simple_interpolationfcn<tom::transf::InterpolTriCubic<float > >() { return true;  }
template<>           inline bool is_simple_interpolationfcn<tom::transf::InterpolTriCubic<double> >() { return true;  }
template<>           inline bool is_simple_interpolationfcn<tom::transf::InterpolCubicSpline<float > >() { return true;  }
template<>           inline bool is_simple_interpolationfcn<tom::transf::InterpolCubicSpline<double> >() { return true;  }
template<>           inline bool is_simple_interpolationfcn<tom::transf::InterpolFourierSpline<float > >() { return true;  }
template<>           inline bool is_simple_interpolationfcn<tom::transf::InterpolFourierSpline<double> >() { return true;  }
}









/***********************************************************************//**
 * \brief Initialise the interpolation object with the volume.
 *
 * If the interpolation does not support the stride parameters if should
 * throw tom::transf::memory_alignment_not_supported.
 **************************************************************************/
template<typename T>
void tom::transf::InterpolNearestNeighbour<T>::setvolume(const tom::Volume<T> &src) {

    const std::size_t stridex = src.getStrideX();
    const std::size_t stridey = src.getStrideY();
    const std::size_t stridez = src.getStrideZ();

    if (stridex!=sizeof(T) || stridey%sizeof(T) || stridez%sizeof(T)) {
        throw tom::transf::memory_alignment_not_supported("The volume of InterpolNearestNeighbour must be aligned with sizeof(T) and the fastest dimension must be contigous.");
    }

    this->sizex = src.getSizeX();
    this->sizey = src.getSizeY();
    this->sizez = src.getSizeZ();
    this->stridey = stridey/sizeof(T);
    this->stridez = stridez/sizeof(T);
    this->data = &src.get();
}

template<typename T>
void tom::transf::InterpolNearestNeighbour<T>::setvolume(const T *data,const std::size_t sizex,const std::size_t sizey,const std::size_t sizez,const std::size_t stridey,const std::size_t stridez) {

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
void tom::transf::InterpolCubicSpline<T>::setvolume(const tom::Volume<T> &src) {

    const std::size_t stridex = src.getStrideX();
    const std::size_t stridey = src.getStrideY();
    const std::size_t stridez = src.getStrideZ();

    if (stridex!=sizeof(T) || stridey%sizeof(T) || stridez%sizeof(T)) {
        throw tom::transf::memory_alignment_not_supported("The volume of InterpolCubicSpline must be aligned with sizeof(T) and the fastest dimension must be contigous.");
    }

    this->sizex = src.getSizeX();
    this->sizey = src.getSizeY();
    this->sizez = src.getSizeZ();
    this->stridey = stridey/sizeof(T);
    this->stridez = stridez/sizeof(T);
    this->data = &src.get();

}


template<typename T>
void tom::transf::InterpolCubicSpline<T>::setvolume(const T *data,const std::size_t sizex,const std::size_t sizey,const std::size_t sizez,const std::size_t stridey,const std::size_t stridez) {

    this->sizex = sizex;
    this->sizey = sizey;
    this->sizez = sizez;
    this->stridey = stridey;
    this->stridez = stridez;
    this->data = data;
}


/***********************************************************************
* \brief Initialise the interpolation object with the volume.
*
* If the interpolation does not support the stride parameters it should
* throw tom::transf::memory_alignment_not_supported.
**************************************************************************/
template<typename T>
void tom::transf::InterpolFourierSpline<T>::setvolume(const tom::Volume<T> &src) {
    
    const std::size_t stridex = src.getStrideX();
    const std::size_t stridey = src.getStrideY();
    const std::size_t stridez = src.getStrideZ();
    
    if (stridex!=sizeof(T) || stridey%sizeof(T) || stridez%sizeof(T)) {
        throw tom::transf::memory_alignment_not_supported("The volume of InterpolFourierSpline must be aligned with sizeof(T) and the fastest dimension must be contigous.");
    }
    
    this->sizex = src.getSizeX();
    this->sizey = src.getSizeY();
    this->sizez = src.getSizeZ();
    this->stridey = stridey/sizeof(T);
    this->stridez = stridez/sizeof(T);
    this->data = &src.get();
    
}


template<typename T>
void tom::transf::InterpolFourierSpline<T>::setvolume(const T *data,const std::size_t sizex,const std::size_t sizey,const std::size_t sizez,const std::size_t stridey,const std::size_t stridez) {
    
    this->sizex = sizex;
    this->sizey = sizey;
    this->sizez = sizez;
    this->stridey = stridey;
    this->stridez = stridez;
    this->data = data;
}

/****************************************************************************
 * \brief Performs a linear (projective) transform of a volume.
 *
 * \param[in] src Source volume.
 * \param[out] dst Destination volume. Can be of different size than \c src
 * \param[in] P Transformation matrix. It is a pointer to a 4x4 matrix,
 *      i.e. 16 double values. If the parameter \c is_affinity is true,
 *      the last row is assumed to be [0, 0, 0, 1] and \c P can legally point
 *      to only 12 values. \c P must not be singular (det != 0).
 * \param[in] is_affinity Specifies that the last row of \c P is [0 0 0 1].
 *      Thereby is is faster to perform the rotation. Moreover no point
 *      in the destination volume was originally at infinity.
 *      Particularly isometries (rotation+shift) are affine transformation.
 * \param[in] valinf The value to set the voxel in the destination volume
 *      if the original point was at infinity. This can only happen if the
 *      last (4th) row evaluates to 0.
 * \param[in] interp The interpolation class which implements the
 *      actual interpolation. It must support the methods \c setvolume and \c interpolate.
 *
 * Both volumes have points at coordinates 0, ..., size-1 along each dimension.
 * The original coordinate X_s of point X_d = [x_d, y_d, z_d] is computed by
 * [x_s', y_s', z_s', w_s']^T = P*[x_d, y_d, z_d, 1]^T;
 * X_s = [x_s'/w_s', y_s'/w_s', z_s'/w_s'].\\ X_s is than the coordinate in \c src
 * from where the value in the destination volume at point X_d is interpolated.
 * If the value lies outside of the source volume, \c interp deciedes how to set
 * the volume.
 *
 * In the current implementation, the volume is copied if the fastest dimension (X)
 * is not aligned (getStrideX() != sizeof(T).
 *
 * The function is self assignment save (and creates a copy of src in case of self assignment)
 * However it does only detect self assignment if the volumes have the same data pointer.
 **************************************************************************/
template<typename T, typename TINTERP>
void tom::transf::transform(const tom::Volume<T> &src, tom::Volume<T> &dst, const double *P, bool is_affinity, T valinf, TINTERP interp) {

    double Pinv[16];
    { // Compute the inverse of the transformation matrix.
        Pinv[ 0] = P[ 0];   Pinv[ 1] = P[ 1];   Pinv[ 2] = P[ 2];   Pinv[ 3] = P[ 3];
        Pinv[ 4] = P[ 4];   Pinv[ 5] = P[ 5];   Pinv[ 6] = P[ 6];   Pinv[ 7] = P[ 7];
        Pinv[ 8] = P[ 8];   Pinv[ 9] = P[ 9];   Pinv[10] = P[10];   Pinv[11] = P[11];
        if (is_affinity) {
            Pinv[12] = 0.;      Pinv[13] = 0.;      Pinv[14] = 0.;      Pinv[15] = 1.;
        } else {
            Pinv[12] = P[12];   Pinv[13] = P[13];   Pinv[14] = P[14];   Pinv[15] = P[15];
            is_affinity = are_same(Pinv[12],0.) && are_same(Pinv[13],0.) && are_same(Pinv[14],0.) && are_same(Pinv[15],0.);
            if (is_affinity && Pinv[15]!=1.) {
                // Rescale the matrix to have an affinity (last row [0 0 0 1]
                for (int i=0; i<12; i++) {
                    Pinv[i] /= Pinv[15];
                }
                Pinv[15] = 1.;
            }
        }
        if (!invert_4x4(Pinv, Pinv, 4)) {
            throw std::invalid_argument("The transformation matrix is non invertible");
        }
        if (is_affinity) {
            Pinv[12] = 0.;      Pinv[13] = 0.;      Pinv[14] = 0.;      Pinv[15] = 1.;
        }
    }

    // Pointers to make a copy of the src or dst if the alignment does is not supported.
    std::auto_ptr<tom::Volume<T> > src_copy;
    const tom::Volume<T> *src_p = &src;

    if (&src.get() == &dst.get()) {
        // Avoid self assignment.
        // ATTENTION: self assignment is not detected if the volumes only share
        // some common memory (subvolume).
        src_copy.reset(new tom::Volume<T>(src));
        src_p = src_copy.get();
    }


    // See whether interp supports the stride parameters of the source volume.
    try {
        interp.setvolume(*src_p);
    } catch (tom::transf::memory_alignment_not_supported &e) {
        // In case of an exception, make a copy...
        if (src_copy.get()) {
            // Made already a copy...
            throw;
        }
        src_copy.reset(new tom::Volume<T>(src));
        src_p = src_copy.get();

        // If there is still an exception, don't catch it anymore.
        interp.setvolume(*src_p);
    }

    std::auto_ptr<tom::Volume<T> > dst_copy;
    tom::Volume<T> *dst_p = &dst;

    std::size_t dst_stridex = dst.getStrideX();
    std::size_t dst_stridey = dst.getStrideY();
    std::size_t dst_stridez = dst.getStrideZ();
    if (dst_stridex!=sizeof(T) || dst_stridey%sizeof(T) || dst_stridez%sizeof(T)) {
        // Not yet implemented... therefore simply make a copy.
        dst_copy.reset(new tom::Volume<T>(dst));
        dst_p = dst_copy.get();

        dst_stridex = dst_p->getStrideX();
        dst_stridey = dst_p->getStrideY();
        dst_stridez = dst_p->getStrideZ();
    }
    dst_stridex /= sizeof(T); dst_stridey /= sizeof(T); dst_stridez /= sizeof(T);

    const std::size_t dst_sizex = dst_p->getSizeX();
    const std::size_t dst_sizey = dst_p->getSizeY();
    const std::size_t dst_sizez = dst_p->getSizeZ();

    if (is_affinity  &&
        are_same(Pinv[ 0], 1.) && are_same(Pinv[ 1], 0.) && are_same(Pinv[ 2], 0.) && are_same(Pinv[ 3], 0.) &&
        are_same(Pinv[ 4], 0.) && are_same(Pinv[ 5], 1.) && are_same(Pinv[ 6], 0.) && are_same(Pinv[ 7], 0.) &&
        are_same(Pinv[ 8], 0.) && are_same(Pinv[ 9], 0.) && are_same(Pinv[10], 1.) && are_same(Pinv[11], 0.) &&
        ::is_simple_interpolationfcn<TINTERP>() && dst_p->is_equal_size(*src_p)) {
        /* Projection is identity. Do only copy
        Only in case of certain interpolation classes (is_simple_interpolationfcn).
        Otherwise it may be desired to do the "interpolation" nonetheless,
        for example, if the interpolation-function is a filter in time domain over
        several neighbours in the volume space. */
        dst_p->setValues(*src_p);
    } else {

        typedef typename TINTERP::idx_type FLOATTYPE;


        T *pdst = &dst_p->get();

        /* Copy the projective transformation into single variables. */
        const FLOATTYPE p00 = Pinv[ 0];
        const FLOATTYPE p01 = Pinv[ 1];
        const FLOATTYPE p02 = Pinv[ 2];
        const FLOATTYPE p03 = Pinv[ 3];
        const FLOATTYPE p10 = Pinv[ 4];
        const FLOATTYPE p11 = Pinv[ 5];
        const FLOATTYPE p12 = Pinv[ 6];
        const FLOATTYPE p13 = Pinv[ 7];
        const FLOATTYPE p20 = Pinv[ 8];
        const FLOATTYPE p21 = Pinv[ 9];
        const FLOATTYPE p22 = Pinv[10];
        const FLOATTYPE p23 = Pinv[11];
        const FLOATTYPE p30 = Pinv[12];
        const FLOATTYPE p31 = Pinv[13];
        const FLOATTYPE p32 = Pinv[14];
        const FLOATTYPE p33 = Pinv[15];
        FLOATTYPE tmpz_00, tmpz_01, tmpz_02;
        FLOATTYPE tmpy_00, tmpy_01, tmpy_02;
        std::size_t dstx, dsty, dstz;

        interp.setvolume(*src_p);

        dst_stridez -= dst_sizey*dst_stridey;
        dst_stridey -= dst_sizex;
        if (is_affinity) {
            for (dstz=0; dstz<dst_sizez; dstz++) {
                tmpz_00 = p02*dstz + p03;
                tmpz_01 = p12*dstz + p13;
                tmpz_02 = p22*dstz + p23;
                for (dsty=0; dsty<dst_sizey; dsty++) {
                    tmpy_00 = p01*dsty + tmpz_00;
                    tmpy_01 = p11*dsty + tmpz_01;
                    tmpy_02 = p21*dsty + tmpz_02;
                    for (dstx=0; dstx<dst_sizex; dstx++) {
                        *pdst++ = interp.interpolate(p00*dstx + tmpy_00, p10*dstx + tmpy_01, p20*dstx + tmpy_02);
                    }
                    pdst += dst_stridey;
                }
                pdst += dst_stridez;
            }
        } else {
            FLOATTYPE pw, tmpz_03, tmpy_03;
            for (dstz=0; dstz<dst_sizez; dstz++) {
                tmpz_00 = p02*dstz + p03;
                tmpz_01 = p12*dstz + p13;
                tmpz_02 = p22*dstz + p23;
                tmpz_03 = p32*dstz + p33;
                for (dsty=0; dsty<dst_sizey; dsty++) {
                    tmpy_00 = p01*dsty + tmpz_00;
                    tmpy_01 = p11*dsty + tmpz_01;
                    tmpy_02 = p21*dsty + tmpz_02;
                    tmpy_03 = p31*dsty + tmpz_03;
                    for (dstx=0; dstx<dst_sizex; dstx++) {
                        if((pw = p30*dstx + tmpy_03) == 0.) {
                            *pdst++ = valinf;
                        } else {
                            *pdst++ = interp.interpolate((p00*dstx + tmpy_00) / pw, (p10*dstx + tmpy_01) / pw, (p20*dstx + tmpy_02) / pw);
                        }
                    }
                    pdst += dst_stridey;
                }
                pdst += dst_stridez;
            }
        }
    }
}

/****************************************************************************
* \brief Transform according to the given matrix.
*******************************************************************************/
template<typename T>
void tom::transf::transformSpline(const tom::Volume<T> &src, tom::Volume<T> &v_rot, const tom::Volume<T> &mtx, T defaultval) {

    double P[16] = { mtx.get(0,0,0), mtx.get(0,1,0), mtx.get(0,2,0), mtx.get(0,3,0),
		             mtx.get(1,0,0), mtx.get(1,1,0), mtx.get(1,2,0), mtx.get(1,3,0),
		             mtx.get(2,0,0), mtx.get(2,1,0), mtx.get(2,2,0), mtx.get(2,3,0),
		             mtx.get(3,0,0), mtx.get(3,1,0), mtx.get(3,2,0), mtx.get(3,3,0)};

    tom::transf::transform<T, InterpolCubicSpline<T> >(src, v_rot, P, true, 0, InterpolCubicSpline<T>(defaultval));

}

/***************************************************************************
 * \brief Performs a linear (projective) transform of a volume.
 *
 * \param[in] src Source volume.
 * \param[out] dst Destination volume. Can be of different size than \c src
 * \param[in] P Transformation matrix. It is a pointer to a 4x4 matrix,
 *      i.e. 16 double values. If the parameter \c is_affinity is true,
 *      the last row is assumed to be [0, 0, 0, 1] and \c P can legally point
 *      to only 12 values. \c P must not be singular (det != 0).
 * \param[in] is_affinity Specifies that the last row of \c P is [0 0 0 1].
 *      Thereby is is faster to perform the rotation. Moreover no point
 *      in the destination volume was originally at infinity.
 *      Particularly isometries (rotation+shift) are affine transformation.
 * \param[in] valinf The value to set the voxel in the destination volume
 *      if the original point was at infinity. This can only happen if the
 *      last (4th) row evaluates to 0.
 * \param[in] interp The interpolation class which implements the
 *      actual interpolation. It must support the methods \c setvolume and \c interpolate.
 *
 * Both volumes have points at coordinates 0, ..., size-1 along each dimension.
 * The original coordinate X_s of point X_d = [x_d, y_d, z_d] is computed by
 * [x_s', y_s', z_s', w_s']^T = P*[x_d, y_d, z_d, 1]^T;
 * X_s = [x_s'/w_s', y_s'/w_s', z_s'/w_s'].\\ X_s is than the coordinate in \c src
 * from where the value in the destination volume at point X_d is interpolated.
 * If the value lies outside of the source volume, \c interp deciedes how to set
 * the volume.
 *
 * In the current implementation, the volume is copied if the fastest dimension (X)
 * is not aligned (getStrideX() != sizeof(T).
 *
 * The function is self assignment save (and creates a copy of src in case of self assignment)
 * However it does only detect self assignment if the volumes have the same data pointer.
 **************************************************************************/
template<typename T, typename TINTERP>
void tom::transf::transformFourier(const tom::Volume<T> &src, tom::Volume<T> &dst, const double *P, bool is_affinity, T valinf, TINTERP interp,double shiftX,double shiftY,double shiftZ) {

    double Pinv[16];
    { // Compute the inverse of the transformation matrix.
        Pinv[ 0] = P[ 0];   Pinv[ 1] = P[ 1];   Pinv[ 2] = P[ 2];   Pinv[ 3] = P[ 3];
        Pinv[ 4] = P[ 4];   Pinv[ 5] = P[ 5];   Pinv[ 6] = P[ 6];   Pinv[ 7] = P[ 7];
        Pinv[ 8] = P[ 8];   Pinv[ 9] = P[ 9];   Pinv[10] = P[10];   Pinv[11] = P[11];
        if (is_affinity) {
            Pinv[12] = 0.;      Pinv[13] = 0.;      Pinv[14] = 0.;      Pinv[15] = 1.;
        } else {
            Pinv[12] = P[12];   Pinv[13] = P[13];   Pinv[14] = P[14];   Pinv[15] = P[15];
            is_affinity = are_same(Pinv[12],0.) && are_same(Pinv[13],0.) && are_same(Pinv[14],0.) && are_same(Pinv[15],0.);
            if (is_affinity && Pinv[15]!=1.) {
                // Rescale the matrix to have an affinity (last row [0 0 0 1]
                for (int i=0; i<12; i++) {
                    Pinv[i] /= Pinv[15];
                }
                Pinv[15] = 1.;
            }
        }
        if (!invert_4x4(Pinv, Pinv, 4)) {
            throw std::invalid_argument("The transformation matrix is non invertible");
        }
        if (is_affinity) {
            Pinv[12] = 0.;      Pinv[13] = 0.;      Pinv[14] = 0.;      Pinv[15] = 1.;
        }
    }

    // Pointers to make a copy of the src or dst if the alignment does is not supported.
    std::auto_ptr<tom::Volume<T> > src_copy;
    const tom::Volume<T> *src_p = &src;

    if (&src.get() == &dst.get()) {
        // Avoid self assignment.
        // ATTENTION: self assignment is not detected if the volumes only share
        // some common memory (subvolume).
        src_copy.reset(new tom::Volume<T>(src));
        src_p = src_copy.get();
    }


    // See whether interp supports the stride parameters of the source volume.
    try {
        interp.setvolume(*src_p);
    } catch (tom::transf::memory_alignment_not_supported &e) {
        // In case of an exception, make a copy...
        if (src_copy.get()) {
            // Made already a copy...
            throw;
        }
        src_copy.reset(new tom::Volume<T>(src));
        src_p = src_copy.get();

        // If there is still an exception, don't catch it anymore.
        interp.setvolume(*src_p);
    }

    std::auto_ptr<tom::Volume<T> > dst_copy;
    tom::Volume<T> *dst_p = &dst;

    std::size_t dst_stridex = dst.getStrideX();
    std::size_t dst_stridey = dst.getStrideY();
    std::size_t dst_stridez = dst.getStrideZ();
    if (dst_stridex!=sizeof(T) || dst_stridey%sizeof(T) || dst_stridez%sizeof(T)) {
        // Not yet implemented... therefore simply make a copy.
        dst_copy.reset(new tom::Volume<T>(dst));
        dst_p = dst_copy.get();

        dst_stridex = dst_p->getStrideX();
        dst_stridey = dst_p->getStrideY();
        dst_stridez = dst_p->getStrideZ();
    }
    dst_stridex /= sizeof(T); dst_stridey /= sizeof(T); dst_stridez /= sizeof(T);

    const std::size_t dst_sizex = dst_p->getSizeX();
    const std::size_t dst_sizey = dst_p->getSizeY();
    const std::size_t dst_sizez = dst_p->getSizeZ();

	typedef typename TINTERP::idx_type FLOATTYPE;


	T *pdst = &dst_p->get();
	T value;
	/* Copy the projective transformation into single variables. */
	const FLOATTYPE p00 = Pinv[ 0];
	const FLOATTYPE p01 = Pinv[ 1];
	const FLOATTYPE p02 = Pinv[ 2];
	const FLOATTYPE p03 = Pinv[ 3];
	const FLOATTYPE p10 = Pinv[ 4];
	const FLOATTYPE p11 = Pinv[ 5];
	const FLOATTYPE p12 = Pinv[ 6];
	const FLOATTYPE p13 = Pinv[ 7];
	const FLOATTYPE p20 = Pinv[ 8];
	const FLOATTYPE p21 = Pinv[ 9];
	const FLOATTYPE p22 = Pinv[10];
	const FLOATTYPE p23 = Pinv[11];
	const FLOATTYPE p30 = Pinv[12];
	const FLOATTYPE p31 = Pinv[13];
	const FLOATTYPE p32 = Pinv[14];
	const FLOATTYPE p33 = Pinv[15];
	FLOATTYPE tmpz_00, tmpz_01, tmpz_02;
	FLOATTYPE tmpy_00, tmpy_01, tmpy_02;
	std::size_t dstx, dsty, dstz;

	interp.setvolume(*src_p);

	dst_stridez -= dst_sizey*dst_stridey;
	dst_stridey -= dst_sizex;

	// calculate the center position (given the full volume)
	std::size_t center_x = dst_sizex/2;
	std::size_t center_y = dst_sizey/2;
	std::size_t center_z = dst_sizez/2;

	int xx,yy,zz;

	if (is_affinity) {
		// loop only over the half size of the fourier volume and asign the complex conjugate to the other index
		for (dstz=0; dstz<dst_sizez/2+1; dstz++) {
			tmpz_00 = p02*dstz + p03;
			tmpz_01 = p12*dstz + p13;
			tmpz_02 = p22*dstz + p23;
			for (dsty=0; dsty<dst_sizey; dsty++) {
				tmpy_00 = p01*dsty + tmpz_00;
				tmpy_01 = p11*dsty + tmpz_01;
				tmpy_02 = p21*dsty + tmpz_02;
				for (dstx=0; dstx<dst_sizex; dstx++) {

                    xx = dstx - center_x;
                    yy = dsty - center_y;
                    zz = dstz - center_z;

					value = interp.interpolate(p00*dstx + tmpy_00, p10*dstx + tmpy_01, p20*dstx + tmpy_02);

					T shift = T(cos(2*PI/dst_sizex*xx*shiftX), -sin(2*PI/dst_sizex*xx*shiftX)) *
							  T(cos(2*PI/dst_sizey*yy*shiftY), -sin(2*PI/dst_sizey*yy*shiftY)) *
							  T(cos(2*PI/dst_sizez*zz*shiftZ), -sin(2*PI/dst_sizez*zz*shiftZ));

					*pdst++ = value * shift;
					try {
						dst_p->get(2*center_x-dstx, 2*center_y-dsty, 2*center_z-dstz) = std::conj(value * shift );
					}
					catch(...) // do not catch the index out of bound exception
					{
					}
				}
				pdst += dst_stridey;
			}
			pdst += dst_stridez;
		}
	} else {
		FLOATTYPE pw, tmpz_03, tmpy_03;
		// for (dstz=0; dstz<dst_sizez; dstz++) {
		for (dstz=0; dstz<dst_sizez/2+1; dstz++) {
			tmpz_00 = p02*dstz + p03;
			tmpz_01 = p12*dstz + p13;
			tmpz_02 = p22*dstz + p23;
			tmpz_03 = p32*dstz + p33;
			for (dsty=0; dsty<dst_sizey; dsty++) {
				tmpy_00 = p01*dsty + tmpz_00;
				tmpy_01 = p11*dsty + tmpz_01;
				tmpy_02 = p21*dsty + tmpz_02;
				tmpy_03 = p31*dsty + tmpz_03;
				for (dstx=0; dstx<dst_sizex; dstx++) {
					if((pw = p30*dstx + tmpy_03) == 0.) {
						*pdst++ = valinf;
					} else {
						value = interp.interpolate((p00*dstx + tmpy_00) / pw, (p10*dstx + tmpy_01) / pw, (p20*dstx + tmpy_02) / pw);
						*pdst++ = value;
						try {
							dst_p->get(2*center_x-dstx, 2*center_y-dsty, 2*center_z-dstz) = std::conj(value);
						}
						catch(...) // do not catch the index out of bound exception
						{
						}
					}
				}
				pdst += dst_stridey;
			}
			pdst += dst_stridez;
		}
	}
}



#ifdef __USE_LAPACK
    #if defined _WIN32
        #define dgetri_ dgetri
        #define dgetrf_ dgetrf
    #endif
    #ifdef __cplusplus
        #define __EXTERN_LAPACK extern "C"
    #else
        #define __EXTERN_LAPACK extern
    #endif

/****************************************************************************//**
 * \brief Prototype for LAPACK functions needed for matrix inversion.
 *
 * See http://www.netlib.org/lapack/double/dgetri.f
 *******************************************************************************/
__EXTERN_LAPACK int dgetri_(int *N, double *A, int *LDA, int *IPIV, double *WORK, int *LWORK, int *INFO );
/** \copydoc dgetri_ */
__EXTERN_LAPACK int dgetrf_(int *N, int *M, double *A, int *LDA, int *IPIV, int *INFO);
#endif





/****************************************************************************//**
 * \brief Compute inverse of a matrix.
 *
 * \param[in]  P the input matrix.
 * \param[out] Pinv the output matrix inverted.
 * \param[in]  N The size of the matrix.
 * \returns 0 if the input arguments are non valid or the matrix is
 *   singular.
 *
 * P and Pinv are matrices of size N and lie continous in memory
 * (for example the third element in the second row is
 * P[1,2] = P[1*N + 2]).
 * Is self assignment safe.
 * If the function returns 0 the content of Pinv is undefined.
 *******************************************************************************/
int tom::transf::invert_4x4(const double *P, double *Pinv, int N) {

#ifdef __USE_LAPACK
    int INFO;
    int LWORK = N*N;
    int *IPIV = 0;
    double *WORK = 0;

    if (!(IPIV = (int *)malloc(sizeof(int)*N)) ||
        !(WORK = (double *)malloc(sizeof(double)*LWORK))) {
        if (IPIV) { free(IPIV); }
        return 0;
    }

    if (P != Pinv) {
        size_t i;
        size_t NN = N*N;
        for (i=0; i<NN; i++) {
            Pinv[i] = P[i];
        }
    }

    dgetrf_(&N, &N, Pinv, &N, IPIV, &INFO);

    if (INFO) {
        free(IPIV);
        free(WORK);
        return 0;
    }

    dgetri_(&N, Pinv, &N, IPIV, WORK, &LWORK, &INFO);

    free(IPIV);
    free(WORK);
#else
    const double *Plocal = P;
    int i;
    double div;
    double *Pl = 0;
    if (N==1) {
        Pinv[0] = 1./P[0];
    } else if (N==2) {
        if (P == Pinv) {
            Pl = (double *)malloc(sizeof(double) * N*N);
            if (!Pl) {
                return 0;
            }
            for (i=0; i<N*N; i++) {
                Pl[i] = P[i];
            }
            Plocal = Pl;
        }
        div = Plocal[0]*Plocal[3] - Plocal[1]*Plocal[2];
        if (!div) {
            return 0;
        }
        Pinv[0] =   Plocal[3] / div;
        Pinv[1] = - Plocal[1] / div;
        Pinv[2] = - Plocal[2] / div;
        Pinv[3] =   Plocal[0] / div;
        if (Pl) { free(Pl); }
    } else if (N==3) {
        if (P == Pinv) {
            Pl = (double *)malloc(sizeof(double) * N*N);
            if (!Pl) {
                return 0;
            }
            for (i=0; i<N*N; i++) {
                Pl[i] = P[i];
            }
            Plocal = Pl;
        }
        #if 0
        | a11 a12 a13 |-1             |   a33a22-a32a23  -(a33a12-a32a13)   a23a12-a22a13  |
        | a21 a22 a23 |    =  1/DET * | -(a33a21-a31a23)   a33a11-a31a13  -(a23a11-a21a13) |
        | a31 a32 a33 |               |   a32a21-a31a22  -(a32a11-a31a12)   a22a11-a21a12  |

        with DET  =  a11(a33a22-a32a23)-a21(a33a12-a32a13)+a31(a23a12-a22a13)
        #endif
        div =        + Plocal[0]*(Plocal[8]*Plocal[4] - Plocal[7]*Plocal[5])
                     - Plocal[3]*(Plocal[8]*Plocal[1] - Plocal[7]*Plocal[2])
                     + Plocal[6]*(Plocal[5]*Plocal[1] - Plocal[4]*Plocal[2]);
        if (!div) {
            return 0;
        }
        Pinv[0] =   (Plocal[4]*Plocal[8] - Plocal[5]*Plocal[7]) / div;
        Pinv[1] =   (Plocal[7]*Plocal[2] - Plocal[8]*Plocal[1]) / div;
        Pinv[2] =   (Plocal[1]*Plocal[5] - Plocal[2]*Plocal[4]) / div;
        Pinv[3] =   (Plocal[5]*Plocal[6] - Plocal[3]*Plocal[8]) / div;
        Pinv[4] =   (Plocal[0]*Plocal[8] - Plocal[2]*Plocal[6]) / div;
        Pinv[5] =   (Plocal[2]*Plocal[3] - Plocal[0]*Plocal[5]) / div;
        Pinv[6] =   (Plocal[3]*Plocal[7] - Plocal[4]*Plocal[6]) / div;
        Pinv[7] =   (Plocal[1]*Plocal[6] - Plocal[0]*Plocal[7]) / div;
        Pinv[8] =   (Plocal[0]*Plocal[4] - Plocal[1]*Plocal[3]) / div;
        if (Pl) { free(Pl); }
    } else if (N==4) {
        double Ainv[9];
        double D_mCAinvB__inv;
        double CAinv[3];
        double AinvB[3];

        Ainv[0] = P[ 0];
        Ainv[1] = P[ 1];
        Ainv[2] = P[ 2];
        Ainv[3] = P[ 4];
        Ainv[4] = P[ 5];
        Ainv[5] = P[ 6];
        Ainv[6] = P[ 8];
        Ainv[7] = P[ 9];
        Ainv[8] = P[10];
        if (!invert_4x4(Ainv, Ainv, 3)) {
            return 0;
        }
        CAinv[0] = P[12]*Ainv[0] + P[13]*Ainv[3] + P[14]*Ainv[6];
        CAinv[1] = P[12]*Ainv[1] + P[13]*Ainv[4] + P[14]*Ainv[7];
        CAinv[2] = P[12]*Ainv[2] + P[13]*Ainv[5] + P[14]*Ainv[8];
        AinvB[0] = Ainv[0]*P[ 3] + Ainv[1]*P[ 7] + Ainv[2]*P[11];
        AinvB[1] = Ainv[3]*P[ 3] + Ainv[4]*P[ 7] + Ainv[5]*P[11];
        AinvB[2] = Ainv[6]*P[ 3] + Ainv[7]*P[ 7] + Ainv[8]*P[11];
        D_mCAinvB__inv = P[15] - ( CAinv[0]*P[ 3] + CAinv[1]*P[ 7] + CAinv[2]*P[11] );
        if (!D_mCAinvB__inv) {
            return 0;
        }
        D_mCAinvB__inv = 1. / D_mCAinvB__inv;

        Pinv[ 0] = Ainv[0] + AinvB[0]*D_mCAinvB__inv*CAinv[0];
        Pinv[ 1] = Ainv[1] + AinvB[0]*D_mCAinvB__inv*CAinv[1];
        Pinv[ 2] = Ainv[2] + AinvB[0]*D_mCAinvB__inv*CAinv[2];
        Pinv[ 3] = - AinvB[0] * D_mCAinvB__inv;
        Pinv[ 4] = Ainv[3] + AinvB[1]*D_mCAinvB__inv*CAinv[0];
        Pinv[ 5] = Ainv[4] + AinvB[1]*D_mCAinvB__inv*CAinv[1];
        Pinv[ 6] = Ainv[5] + AinvB[1]*D_mCAinvB__inv*CAinv[2];
        Pinv[ 7] = - AinvB[1] * D_mCAinvB__inv;
        Pinv[ 8] = Ainv[6] + AinvB[2]*D_mCAinvB__inv*CAinv[0];
        Pinv[ 9] = Ainv[7] + AinvB[2]*D_mCAinvB__inv*CAinv[1];
        Pinv[10] = Ainv[8] + AinvB[2]*D_mCAinvB__inv*CAinv[2];
        Pinv[11] = - AinvB[2] * D_mCAinvB__inv;
        Pinv[12] = - D_mCAinvB__inv * CAinv[0];
        Pinv[13] = - D_mCAinvB__inv * CAinv[1];
        Pinv[14] = - D_mCAinvB__inv * CAinv[2];
        Pinv[15] = D_mCAinvB__inv;
    } else {
        //printf("ERROR: NOT IMPLEMENTED MATRIX INVERSION OF ORDER>4 (%s %d)", __FILE__, __LINE__);
        assert(0);
        return 0;
    }
#endif
    return 1;

}





/****************************************************************************//**
 * \brief Concatenates a series of shifts and rotations into one transformation.
 *
 * \param[out] P If not 0, it must point to 3x4 doubles or 3x3 doubles. This output will be
 *  filled with the resulting transformation matrix. The matrix
 *  elements are at position P[num_row*4 + num_column] or P[num_row*3 + num_column],
 *  depending on the value of \c homogenous.
 * \param[in] homogenous If true, \c P points to a homogenous 3x4 matrix. Then,
 *  the 4th column contains the resulting shifts of the transformation. Otherwise
 *  it points to a 3x3 matrix, containing only the rotation matrix.
 * \param[in] continue_P If true, the content in \c P is taken as starting
 *  transformation. Otherwise the matrix is initialised with the identity matrix
 *  before starting to apply the transformations.
 *  If the input in \c P is not a isometrie (orthogonal matrix), the result is
 *  as applying a isometrie to an affinity transformation. However, the result
 *  will not be a projective transformation, as the 4th row of \c P is assumed
 *  to be [0,0,0,1] and is not touched.
 * \param[out] angles_out If not 0, it must point to 3 double values which will
 *  be filled with the angles [phi,psi,theta] in radians, according to the TOM
 *  convention.
 * \param[out] shifts_out If not 0, it must point to 3 double values which will
 *  hold the shifts of the transformation in X, Y and Z direction, respectively.
 * \param[in] n Number of rotations/shifts (determines the length of the
 *  following parameters).
 * \param[in] angles If not 0, it points to \c n rotation angles in radians.
 * \param[in] axes if not 0, it determine the axes about which the rotation is
 *  performed (0==x-axes, 1==y-axes, 2==z-axes). If one of \c axes or \c angles
 *  is 0, no rotation is performed.
 * \param[in] shifts If not 0, it points to \c 3*n displacements, where
 *  shifts[i*3+0], shifts[i*3+1], shifts[i*3+2] are the x,y and z displacements.
 *
 * First the transformatio rotates as given by angles/axes, and then(!!) the shift
 * is performed.
 *
 * In R^3, coordinate system rotations of the x-, y-, and z-axes in a
 * counterclockwise direction when looking towards the origin give the matrices.
 * http://mathworld.wolfram.com/RotationMatrix.html.
 * Note that this is different from the TOM-convention, which rotates in this
 * notation about -Z,-X,-Z (\see tom::transf::rotate).
 *******************************************************************************/
void tom::transf::sum_rotation(double *P, bool homogenous, bool continue_P, double *angles_out, double *shifts_out, std::size_t n, const double *angles, const int *axes, const double *shifts) {

    if (!P && !angles_out && !shifts_out) {
        return;
    }

    double sina, cosa, p0, p1;
    double P_[12];
    double shifts_out_[3];
    if (!P) {
        P = P_;
        homogenous = false;
    }
    if (!shifts_out) {
        shifts_out = shifts_out_;
    }
    if (!axes) {
        angles = 0;
    }
    //determine size of matrix (set by homogenous) and assign step the y-size
    const int step = homogenous ? 4 : 3;
    if (continue_P && P!=P_) {
        if (homogenous) {
            shifts_out[0] = P[0*step+3];
            shifts_out[1] = P[1*step+3];
            shifts_out[2] = P[2*step+3];
        } else {
            shifts_out[0] = 0.;
            shifts_out[1] = 0.;
            shifts_out[2] = 0.;
        }
    } else {
        P[0*step+0] = 1.;    P[0*step+1] = 0.;    P[0*step+2] = 0.;
        P[1*step+0] = 0.;    P[1*step+1] = 1.;    P[1*step+2] = 0.;
        P[2*step+0] = 0.;    P[2*step+1] = 0.;    P[2*step+2] = 1.;
        shifts_out[0] = 0.;
        shifts_out[1] = 0.;
        shifts_out[2] = 0.;
    }

    std::size_t ni;
    int i;
    //loop over all matrix lines
    //remember, this happens line by line => this will be processed n times!
    for (ni=0; ni<n; ni++) {
        //if current entry in the angles array is >0 then...
        if (angles && angles[ni]!=0.) {
            //precalculate
            sina = sin(angles[ni]);
            cosa = cos(angles[ni]);
//            std::cout << std::endl;
//            std::cout << ni << std::endl;
//            std::cout << sina << " " << cosa << std::endl;

            //apply rotation matrix combination based on axis index
            switch (axes[ni] % 3) {
                case 0:
                    //comments identical for the lines below:
                    //apply rotation to vector element in the current line
                    for (i=0; i<3; i++) {
                        p0 = P[1*step+i];
                        p1 = P[2*step+i];
                        P[1*step+i] = + p0*cosa + p1*sina;
                        P[2*step+i] = - p0*sina + p1*cosa;
                    }
                    //if shifts have been specified, apply a transformation to them as well
                    if (shifts) {
                        p0 = shifts_out[1];
                        p1 = shifts_out[2];
//                        std::cout << "1 "<< p0 << " " << p1 << std::endl;
                        shifts_out[1] = + p0*cosa + p1*sina;
                        shifts_out[2] = - p0*sina + p1*cosa;
//                        std::cout << shifts_out[1] << " " << shifts_out[2]<< std::endl;
                    }
                    break;
                case 1:
                    for (i=0; i<3; i++) {
                        p0 = P[0*step+i];
                        p1 = P[2*step+i];
                        P[0*step+i] = + p0*cosa - p1*sina;
                        P[2*step+i] = + p0*sina + p1*cosa;
                    }
                    if (shifts) {
                        p0 = shifts_out[0];
                        p1 = shifts_out[2];
//                        std::cout << "2 "<< p0 << " " << p1 << std::endl;
                        shifts_out[0] = + p0*cosa - p1*sina;
                        shifts_out[2] = + p0*sina + p1*cosa;
//                        std::cout << shifts_out[1] << " " << shifts_out[2]<< std::endl;
                    }
                    break;
                case 2:
                    for (i=0; i<3; i++) {
                        p0 = P[0*step+i];
                        p1 = P[1*step+i];
                        P[0*step+i] = + p0*cosa + p1*sina;
                        P[1*step+i] = - p0*sina + p1*cosa;
                    }
                    if (shifts) {
                        p0 = shifts_out[0];
                        p1 = shifts_out[1];
//                        std::cout << "3 "<< p0 << " " << p1 << std::endl;
                        shifts_out[0] = + p0*cosa + p1*sina;
                        shifts_out[1] = - p0*sina + p1*cosa;
//                        std::cout << shifts_out[1] << " " << shifts_out[2]<< std::endl;
                    }
                    break;
            }
        }

        if (shifts) {
            shifts_out[0] += shifts[3*ni+0];
            shifts_out[1] += shifts[3*ni+1];
            shifts_out[2] += shifts[3*ni+2];
        }
    }


    if (homogenous && P!=P_) {
        P[3] = shifts_out[0];
        P[7] = shifts_out[1];
        P[11] = shifts_out[2];
    }

    if (angles_out) {
        // Recompute the rotation angles...
        angles_out[0] = atan2(P[2*step+0],  P[2*step+1]);
        angles_out[1] = atan2(P[0*step+2], -P[1*step+2]);
        angles_out[2] = acos (P[2*step+2]);
    }
}





/****************************************************************************//**
 * \brief First rotate the volume around a point and then shifts it.
 *
 * Uses the tom-convection for rotation. It does not rotate the object,
 * but its coordinate system with angles phi, theta, psi (in that order).
 *******************************************************************************/
template<typename T>
void tom::transf::rotate(const tom::Volume<T> &src, tom::Volume<T> &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, double shiftx_pre, double shifty_pre, double shiftz_pre, double shiftx_post, double shifty_post, double shiftz_post, T defaultval) {

    double P[16] = { 0., 0., 0., 0.,
                     0., 0., 0., 0.,
                     0., 0., 0., 0.,
                     0., 0., 0., 1. };

    const int axes[4] = {   0, 2, 0, 2 };
    const double angles[4] = { 0., -phi, -theta, -psi };
    const double shifts[3*4] = { shiftx_pre - centerx, shifty_pre - centery, shiftz_pre - centerz, 0,0,0, 0,0,0, centerx+shiftx_post, centery+shifty_post, centerz+shiftz_post };

    tom::transf::sum_rotation(P, true, false, 0, 0, 4, angles, axes, shifts);

    tom::transf::transform<T, typename InterpolDefault<T>::type>(src, v_rot, P, true, 0, typename InterpolDefault<T>::type(defaultval));

}



/****************************************************************************//**
* \brief First rotate the volume around a point and then shifts it.
*
* Uses the tom-convection for rotation. It does not rotate the object,
* but its coordinate system with angles phi, theta, psi (in that order).
*******************************************************************************/
template<typename T>
void tom::transf::rotateCubic(const tom::Volume<T> &src, tom::Volume<T> &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, double shiftx_pre, double shifty_pre, double shiftz_pre, double shiftx_post, double shifty_post, double shiftz_post, T defaultval) {

    double P[16] = { 0., 0., 0., 0.,
		0., 0., 0., 0.,
		0., 0., 0., 0.,
		0., 0., 0., 1. };

    const int axes[4] = {   0, 2, 0, 2 };
    const double angles[4] = { 0., -phi, -theta, -psi };
    const double shifts[3*4] = { shiftx_pre - centerx, shifty_pre - centery, shiftz_pre - centerz, 0,0,0, 0,0,0, centerx+shiftx_post, centery+shifty_post, centerz+shiftz_post };

    tom::transf::sum_rotation(P, true, false, 0, 0, 4, angles, axes, shifts);

    tom::transf::transform<T, InterpolTriCubic<T> >(src, v_rot, P, true, 0, InterpolTriCubic<T>(defaultval));

}

/****************************************************************************
* \brief First rotate the volume around a point and then shifts it.
*
* Uses the tom-convection for rotation. It does not rotate the object,
* but its coordinate system with angles phi, theta, psi (in that order).
*******************************************************************************/
template<typename T>
void tom::transf::rotateSpline(const tom::Volume<T> &src, tom::Volume<T> &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, double shiftx_pre, double shifty_pre, double shiftz_pre, double shiftx_post, double shifty_post, double shiftz_post, T defaultval) {

    double P[16] = { 0., 0., 0., 0.,
		0., 0., 0., 0.,
		0., 0., 0., 0.,
		0., 0., 0., 1. };

    const int axes[4] = {   0, 2, 0, 2 };
    const double angles[4] = { 0., -phi, -theta, -psi };
    const double shifts[3*4] = { shiftx_pre - centerx, shifty_pre - centery, shiftz_pre - centerz, 0,0,0, 0,0,0, centerx+shiftx_post, centery+shifty_post, centerz+shiftz_post };

    tom::transf::sum_rotation(P, true, false, 0, 0, 4, angles, axes, shifts);

    tom::transf::transform<T, InterpolCubicSpline<T> >(src, v_rot, P, true, 0, InterpolCubicSpline<T>(defaultval));

}

/****************************************************************************
 * \brief First rotate the volume around a point and then shifts it.
 *
 * Uses the tom-convection for rotation. It does not rotate the object,
 * but its coordinate system with angles phi, theta, psi (in that order).
 *******************************************************************************/
template<typename T>
//void tom::transf::rotateFourierSpline(const tom::Volume<T > &src, tom::Volume<T > &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, double shiftx_pre, double shifty_pre, double shiftz_pre, double shiftx_post, double shifty_post, double shiftz_post, T defaultval) {
void tom::transf::rotateFourierSpline(const tom::Volume<T > &src, tom::Volume<T > &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, T defaultval) {
    double P[16] = { 0., 0., 0., 0.,
		0., 0., 0., 0.,
		0., 0., 0., 0.,
		0., 0., 0., 1. };
    
    const int axes[4] = {   0, 2, 0, 2 };
    const double angles[4] = { 0., -phi, -theta, -psi };
    const double shifts[3*4] = { 0 - centerx, 0 - centery, 0 - centerz, 0,0,0, 0,0,0, centerx+0 , centery+0, centerz+0};
 

    tom::transf::sum_rotation(P, true, false, 0, 0, 4, angles, axes, shifts);
    
    //remove shift influence on transformation by setting P values to 0
    P[12]=P[13]=P[14]=0.0;
    P[15]=1;

    tom::transf::transformFourier<T, InterpolFourierSpline<T > >(src, v_rot, P, true, 0, InterpolFourierSpline<T >(defaultval), 0,0,0);
    
}

/****************************************************************************
 * \brief First rotate the volume around a point and then shifts it.
 *
 * Uses the tom-convection for rotation. It does not rotate the object,
 * but its coordinate system with angles phi, theta, psi (in that order).
 *******************************************************************************/
template<typename T>
//void tom::transf::rotateFourierSpline(const tom::Volume<T > &src, tom::Volume<T > &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, double shiftx_pre, double shifty_pre, double shiftz_pre, double shiftx_post, double shifty_post, double shiftz_post, T defaultval) {
void tom::transf::transformFourierSpline(const tom::Volume<T > &src, tom::Volume<T > &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, double shiftX, double shiftY, double shiftZ,T defaultval) {
    double P[16] = { 0., 0., 0., 0.,
		0., 0., 0., 0.,
		0., 0., 0., 0.,
		0., 0., 0., 1. };

    const int axes[4] = {   0, 2, 0, 2 };
    const double angles[4] = { 0., -phi, -theta, -psi };
    const double shifts[3*4] = { 0 - centerx, 0 - centery, 0 - centerz, 0,0,0, 0,0,0, centerx+0 , centery+0, centerz+0};


    tom::transf::sum_rotation(P, true, false, 0, 0, 4, angles, axes, shifts);

    //remove shift influence on transformation by setting P values to 0
    P[12]=P[13]=P[14]=0.0;
    P[15]=1;

    tom::transf::transformFourier<T, InterpolFourierSpline<T > >(src, v_rot, P, true, 0, InterpolFourierSpline<T >(defaultval), shiftX, shiftY, shiftZ);

}

template<typename T>
void tom::transf::shiftFourier(const tom::Volume<T > &src, tom::Volume<T > &dst, double shift_x, double shift_y, double shift_z) {
    const std::size_t sizex = src.getSizeX();
    const std::size_t sizey = src.getSizeY();
    const std::size_t sizez = src.getSizeZ();

    int ii,jj,kk;

    for (std::size_t i=0; i<sizex; i++)
        for (std::size_t j=0; j<sizey; j++)
            for (std::size_t k=0; k<sizez; k++)
            {
                T value = src.get(i,j,k);

                ii = (i<sizex/2)?i:i-sizex;
                jj = (j<sizey/2)?j:j-sizey;
                kk = (k<sizez/2)?k:k-sizez;
                
                T shift = T(cos(2*PI/sizex*ii*shift_x), -sin(2*PI/sizex*ii*shift_x)) *
                          T(cos(2*PI/sizey*jj*shift_y), -sin(2*PI/sizey*jj*shift_y)) *
                          T(cos(2*PI/sizez*kk*shift_z), -sin(2*PI/sizez*kk*shift_z));
                
                dst.get(i,j,k) = value*shift;
            } 
}


namespace {
template<typename T>
struct for_each__tom__rescale_fs {
    for_each__tom__rescale_fs(const tom::Volume<std::complex<T> > &_vsrc, std::size_t _src_sizez, const tom::Volume<std::complex<T> > &_vdst, std::size_t _dst_sizez, double shiftx, double shifty, double shiftz, bool divide_by_numel)
        :   vsrc(_vsrc),
            vdst(_vdst),
            src_sizez_full(_src_sizez),
            dst_sizez_full(_dst_sizez) {
        src_sizex = vsrc.getSizeX();
        src_sizey = vsrc.getSizeY();
        src_sizez = vsrc.getSizeZ();
        dst_sizex = vdst.getSizeX();
        dst_sizey = vdst.getSizeY();
        dst_sizez = vdst.getSizeZ();
        assert(src_sizez_full/2+1 == src_sizez);
        assert(dst_sizez_full/2+1 == dst_sizez);
        src_sizex_offset = src_sizex/2      + src_sizex%2;
        src_sizey_offset = src_sizey/2      + src_sizey%2;
        src_sizez_offset = src_sizez_full/2 + src_sizez_full%2;
        dst_sizex_offset = dst_sizex/2;
        dst_sizey_offset = dst_sizey/2;
        dst_sizez_offset = dst_sizez_full/2;
        double f = (static_cast<double>(src_sizex)/dst_sizex) * (static_cast<double>(src_sizey)/dst_sizey) * (static_cast<double>(src_sizez_full)/dst_sizez_full);
        if (divide_by_numel) {
            f *= dst_sizex * dst_sizey * dst_sizez_full;
        }
        factor = 1./f;

        do_shift_ = shiftx!=0. || shifty!=0. || shiftz!=0.;
        if (do_shift_) {
            shiftx_ = shiftx / vdst.getSizeX();
            shifty_ = shifty / vdst.getSizeY();
            shiftz_ = shiftz / _dst_sizez;
        }
    }
    void do_src_ifftshift(std::size_t &x, std::size_t &y, std::size_t &z) const {
        x = (x + src_sizex_offset) % src_sizex;
        y = (y + src_sizey_offset) % src_sizey;
        z = (z + src_sizez_offset) % src_sizez_full;
    }
    void do_dst_fftshift(std::size_t &x, std::size_t &y, std::size_t &z) const {
        x = (x + dst_sizex_offset) % dst_sizex;
        y = (y + dst_sizey_offset) % dst_sizey;
        z = (z + dst_sizez_offset) % dst_sizez_full;
    }
    #define _W(cmd)
    bool tranf_coord(std::size_t &x, std::size_t &y, std::size_t &z) const {
        // Compute the fft-shifted destination coordinate...
        _W(std::cout << "[" << std::setw(5) << x << ',' << std::setw(5) << y << ',' << std::setw(5) << z << ']';)
        _W(std::size_t x_ = x;)
        _W(std::size_t y_ = y;)
        _W(std::size_t z_ = z;)

        do_dst_fftshift(x, y, z);
        _W(std::cout << " -> [" << std::setw(5) << x << ',' << std::setw(5) << y << ',' << std::setw(5) << z << ']';)
        x = x + src_sizex/2;
        y = y + src_sizey/2;
        z = z + src_sizez_full/2;
        if (x < dst_sizex/2 || y < dst_sizey/2 || z < dst_sizez_full/2) {
        } else {
            x = x - dst_sizex/2;
            y = y - dst_sizey/2;
            z = z - dst_sizez_full/2 ;
            _W(std::cout << " -> [" << std::setw(5) << x << ',' << std::setw(5) << y << ',' << std::setw(5) << z << ']';)
            if (x >= src_sizex || y>=src_sizey || z>=src_sizez_full) {
            } else {
                do_src_ifftshift(x, y, z);
                if (z < src_sizez) {
                    _W(std::cout << " -> [" << std::setw(5) << x << ',' << std::setw(5) << y << ',' << std::setw(5) << z << ']' << ((x_==x&&y_==y&&z_==z)?'*':' ') << std::endl;)
                    return true;
                }
            }
        }
        _W(std::cout << std::endl;)
        return false;
    }
    #undef _W
    void operator()(std::complex<T> &v, std::size_t x, std::size_t y, std::size_t z) const {
        if (do_shift_) {
            shift_tmpx_ = x;
            shift_tmpy_ = y;
            shift_tmpz_ = z;
        }
        if (tranf_coord(x, y, z)) {
            v = vsrc.get(x, y, z);
            if (do_shift_) {
                if (shift_tmpx_ > dst_sizex/2) {
                    shift_tmp_ = shiftx_ * (static_cast<double>(shift_tmpx_) - dst_sizex);
                } else {
                    shift_tmp_ = shiftx_ * shift_tmpx_;
                }
                if (shift_tmpy_ > dst_sizey/2) {
                    shift_tmp_ += shifty_ * (static_cast<double>(shift_tmpy_) - dst_sizey);
                } else {
                    shift_tmp_ += shifty_ * shift_tmpy_;
                }
                if (shift_tmpz_ > dst_sizez_full/2) {
                    shift_tmp_ += shiftz_ * (static_cast<double>(shift_tmpz_) - dst_sizez_full);
                } else {
                    shift_tmp_ += shiftz_ * shift_tmpz_;
                }
                shift_tmp_ *= -6.2831853071795864769252867665590057684;
                v *= std::complex<T>(std::cos(shift_tmp_), std::sin(shift_tmp_));
            }
            v *= factor;
        } else {
            v = 0;
        }
    }
private:
    const tom::Volume<std::complex<T> > &vsrc;
    const tom::Volume<std::complex<T> > &vdst;
    std::size_t src_sizez_full;
    std::size_t dst_sizez_full;
    std::size_t src_sizex, src_sizey, src_sizez;
    std::size_t dst_sizex, dst_sizey, dst_sizez;
    std::size_t src_sizex_offset, src_sizey_offset, src_sizez_offset;
    std::size_t dst_sizex_offset, dst_sizey_offset, dst_sizez_offset;
    double shiftx_, shifty_, shiftz_;
    bool do_shift_;
    mutable double shift_tmp_;
    mutable std::size_t shift_tmpx_, shift_tmpy_, shift_tmpz_;
    T factor;
};
// Hilfsfuntion fuer rescale_fs.
template<typename T>
void tom__rescale_fs(const tom::Volume<std::complex<T> > &vsrc, std::size_t sizez_src, tom::Volume<std::complex<T> > &vdst, std::size_t sizez_dst, double shiftx, double shifty, double shiftz, bool divide_by_numel) {

    if (vsrc.getSizeZ() != sizez_src/2+1) {
        throw std::invalid_argument("rescale_fs: source volume has wrong size.");
    }
    if (vdst.getSizeZ() != sizez_dst/2+1) {
        throw std::invalid_argument("rescale_fs: destination volume has wrong size.");
    }

    if (&vsrc == &vdst) {
        // Do nothing.
    } else if (vsrc.is_equal_size(vdst) && sizez_src==sizez_dst) {
        vdst.setValues(vsrc);
        if (divide_by_numel) {
            vdst.shift_scale(0., 1./(vdst.getSizeX()*vdst.getSizeY()*sizez_dst));
        }
    } else {
        ::for_each__tom__rescale_fs<T> op(vsrc, sizez_src, vdst, sizez_dst, shiftx, shifty, shiftz, divide_by_numel);
        tom::loop::for_each<tom::Volume<std::complex<T> >, ::for_each__tom__rescale_fs<T> &>(vdst, op);
    }
}
} // namespace
/****************************************************************************//**
 * \brief Rescale a volume in Fourier space.
 *
 * \param[in] vsrc The source volume in time domain.
 * \param[in] sizez_src The size of the full volume.
 * \param[in] vdst The resulting volume in Fourier space.
 * \param[in] sizez_dst The size of the full volume.
 * \param[in] shiftx Shift the resulting volume in Fourier space.
 * \param[in] shifty Shift the resulting volume in Fourier space.
 * \param[in] shiftz Shift the resulting volume in Fourier space.
 *
 * analog to tom_rescale from the TOM-Toolbox.
 *
 * \todo The current implemenation is sub-optimal :)
 *******************************************************************************/
template<typename T>
void tom::transf::rescale_fs(const tom::Volume<std::complex<T> > &vsrc, std::size_t sizez_src, tom::Volume<std::complex<T> > &vdst, std::size_t sizez_dst, double shiftx, double shifty, double shiftz) {
    ::tom__rescale_fs(vsrc, sizez_src, vdst, sizez_dst, shiftx, shifty, shiftz, false);
}



/****************************************************************************//**
 * \brief Rescale a volume in Fourier space.
 *
 * \param[in] vsrc The source volume in time domain.
 * \param[in] vdst The resulting volume.
 * \param[in] shiftx Shift the resulting volume in Fourier space.
 * \param[in] shifty Shift the resulting volume in Fourier space.
 * \param[in] shiftz Shift the resulting volume in Fourier space.
 * \param[in] fftw_flag The planner flag for the fftw.
 *
 * analog to tom_rescale from the TOM-Toolbox.
 *******************************************************************************/
template<typename T>
void tom::transf::rescale_fs(const tom::Volume<T> &vsrc, tom::Volume<T> &vdst, double shiftx, double shifty, double shiftz, unsigned fftw_flag) {

    tom::Volume<T>                 vsrc2_rs(vsrc.getSizeX(), vsrc.getSizeY(), vsrc.getSizeZ()    , tom::fftw::fftw_malloc<T>(), tom::fftw::fftw_free<T>());
    tom::Volume<std::complex<T> >  vsrc2_fs(vsrc.getSizeX(), vsrc.getSizeY(), vsrc.getSizeZ()/2+1, tom::fftw::fftw_malloc<T>(), tom::fftw::fftw_free<T>());
    tom::fftw::Plan<T> plan_src_t2f(vsrc2_rs, vsrc2_fs, fftw_flag | FFTW_DESTROY_INPUT);

    vsrc2_rs.setValues(vsrc);
    plan_src_t2f.execute(vsrc2_rs, vsrc2_fs);

    tom::transf::rescale_fs(vsrc2_fs, vsrc.getSizeZ(), vdst, shiftx, shifty, shiftz, fftw_flag);
}


/****************************************************************************//**
 * \brief Rescale a volume in Fourier space.
 *
 * \param[in] vsrc The source volume in time domain.
 * \param[in] vdst The resulting volume in Fourier space.
 * \param[in] sizez_dst The size of the full volume.
 * \param[in] shiftx Shift the resulting volume in Fourier space.
 * \param[in] shifty Shift the resulting volume in Fourier space.
 * \param[in] shiftz Shift the resulting volume in Fourier space.
 * \param[in] fftw_flag The planner flag for the fftw.
 *
 * analog to tom_rescale from the TOM-Toolbox.
 *******************************************************************************/
template<typename T>
void tom::transf::rescale_fs(const tom::Volume<T> &vsrc, tom::Volume<std::complex<T> > &vdst, std::size_t sizez_dst, double shiftx, double shifty, double shiftz, unsigned fftw_flag) {

    tom::Volume<T>                 vsrc_ts(vsrc.getSizeX(), vsrc.getSizeY(), vsrc.getSizeZ()    , tom::fftw::fftw_malloc<T>(), tom::fftw::fftw_free<T>());
    tom::Volume<std::complex<T> >  vsrc_fs(vsrc.getSizeX(), vsrc.getSizeY(), vsrc.getSizeZ()/2+1, tom::fftw::fftw_malloc<T>(), tom::fftw::fftw_free<T>());
    tom::fftw::Plan<T> plan_t2f(vsrc_ts, vsrc_fs, fftw_flag | FFTW_DESTROY_INPUT);

    vsrc_ts.setValues(vsrc);
    plan_t2f.execute(vsrc_ts, vsrc_fs);

    ::tom__rescale_fs(vsrc_fs, vsrc.getSizeZ(), vdst, sizez_dst, shiftx, shifty, shiftz, false);
}


/****************************************************************************//**
 * \brief Rescale a volume in Fourier space.
 *
 * \param[in] vsrc The source volume in Fourier space.
 * \param[in] sizez_src The size of the full volume.
 * \param[in] vdst The resulting volume in time domain.
 * \param[in] shiftx Shift the resulting volume in Fourier space.
 * \param[in] shifty Shift the resulting volume in Fourier space.
 * \param[in] shiftz Shift the resulting volume in Fourier space.
 * \param[in] fftw_flag The planner flag for the fftw.
 *
 * analog to tom_rescale from the TOM-Toolbox.
 *******************************************************************************/
template<typename T>
void tom::transf::rescale_fs(const tom::Volume<std::complex<T> > &vsrc, std::size_t sizez_src, tom::Volume<T> &vdst, double shiftx, double shifty, double shiftz, unsigned fftw_flag) {
    tom::Volume<std::complex<T> > vdst_fs(vdst.getSizeX(), vdst.getSizeY(), vdst.getSizeZ()/2+1, tom::fftw::fftw_malloc<T>(), tom::fftw::fftw_free<T>());
    tom::fftw::Plan<T> plan_f2t(vdst_fs, vdst, fftw_flag | FFTW_DESTROY_INPUT);
    ::tom__rescale_fs(vsrc, sizez_src, vdst_fs, vdst.getSizeZ(), shiftx, shifty, shiftz, true);

    plan_f2t.execute(vdst_fs, vdst);
}






/****************************************************************************//**
 * \brief Rescale a binned volume to a new binning size.
 *
 * \param[in] vsrc The source volume.
 * \param[in] binsrcx Binning factor with which the source volume was binned.
 * \param[in] binsrcy Binning factor with which the source volume was binned.
 * \param[in] binsrcz Binning factor with which the source volume was binned.
 * \param[out] vdst The resulting volume.
 * \param[in] bindstx Binning factor with which the destination volume should be binned.
 * \param[in] bindsty Binning factor with which the destination volume should be binned.
 * \param[in] bindstz Binning factor with which the destination volume should be binned.
 * \param[in] sizex X-Size of the source volume as it was originally. The
 *  size of \c vsrc and \c vdst must correspond to a possible binned version of
 *  this size.
 * \param[in] sizey analog to \c sizex
 * \param[in] sizez analog to \c sizey
 * \param[in] default_val This is the value pasted into the fully enlarged volume
 *  where values are thrown away while binning. Setting to 0, makes rescale_rebin use
 *  the mean value of vsrc.
 *
 * The current implementation enlarges first the volume \c vsrc using
 * \c rescale_fs to the original size of \c sizex x \c sizey x \c sizez (if
 * the original sized volume has voxels at the border which were dropped
 * during binning, these border voxels are filled with \c default_val.
 * If the destination volume is also a binned version, the volume is scaled down
 * again using \c tom::transf::bin. This upscaling to the full size could be
 * problematic if the orignal volume is very large.
 *
 * \TODO Implement a better version which ommits the middle step of up
 * scaling to full size.
 *******************************************************************************/
template<typename T>
void tom::transf::rescale_rebin(const tom::Volume<T> &vsrc, std::size_t binsrcx, std::size_t binsrcy, std::size_t binsrcz, tom::Volume<T> &vdst, std::size_t bindstx, std::size_t bindsty, std::size_t bindstz, std::size_t sizex, std::size_t sizey, std::size_t sizez, const T *defaultval, unsigned fftw_flag) {

    const std::size_t binsrc[3] = { binsrcx?binsrcx:1, binsrcy?binsrcy:1, binsrcz?binsrcz:1 };
    const std::size_t bindst[3] = { bindstx?bindstx:1, bindsty?bindsty:1, bindstz?bindstz:1 };


    std::ostringstream ss;

    if (vsrc.getSizeX()>sizex || sizex%vsrc.getSizeX()>=binsrc[0] || sizex/binsrc[0]!=vsrc.getSizeX() ||
        vsrc.getSizeY()>sizey || sizey%vsrc.getSizeY()>=binsrc[1] || sizey/binsrc[1]!=vsrc.getSizeY() ||
        vsrc.getSizeZ()>sizez || sizez%vsrc.getSizeZ()>=binsrc[2] || sizez/binsrc[2]!=vsrc.getSizeZ()) {
        ss << "tom::transf::rescale_rebin: source volume of size " << vsrc.getSizeX() << 'x' << vsrc.getSizeY() << 'x' << vsrc.getSizeZ() << " does not result from binning a " << sizex << 'x' << sizey << 'x' << sizez << " volume.";
        throw std::invalid_argument(ss.str());
    }
    if (vdst.getSizeX()>sizex || sizex%vdst.getSizeX()>=bindst[0] || sizex/bindst[0]!=vdst.getSizeX() ||
        vdst.getSizeY()>sizey || sizey%vdst.getSizeY()>=bindst[1] || sizey/bindst[1]!=vdst.getSizeY() ||
        vdst.getSizeZ()>sizez || sizez%vdst.getSizeZ()>=bindst[2] || sizez/bindst[2]!=vdst.getSizeZ()) {
        ss << "tom::transf::rescale_rebin: destination volume of size " << vdst.getSizeX() << 'x' << vdst.getSizeY() << 'x' << vdst.getSizeZ() << " does not result from binning a " << sizex << 'x' << sizey << 'x' << sizez << " volume.";
        throw std::invalid_argument(ss.str());
    }

    if (vsrc.is_equal_size(vdst)) {
        vdst.setValues(vsrc);
    } else {
        std::auto_ptr<tom::Volume<T> > v2_;
        const tom::Volume<T> *vlarge;
        if (vsrc.is_equal_size(sizex, sizey, sizez)) {
            vlarge = &vsrc;
        } else {
            tom::Volume<T> *vlarge2, *vlarge3;
            if (vdst.is_equal_size(sizex, sizey, sizez)) {
                vlarge2 = &vdst;
            } else {
                v2_.reset(new tom::Volume<T>(sizex, sizey, sizez, tom::fftw::fftw_malloc<T>(), tom::fftw::fftw_free<T>()));
                vlarge2 = v2_.get();
            }
            std::auto_ptr<tom::Volume<T> > v3_;
            vlarge = vlarge3 = vlarge2;
            if (sizex%binsrc[0] || sizey%binsrc[1] || sizez%binsrc[2]) {
                // The src volume is smaller because the size got truncated
                // while binning. Cut out sub-volumes to set the region
                // outside of the original volume to zero.
                T defaultval_ = defaultval ? *defaultval : tom::mean<T, double>(vsrc);
                tom::Volume<T> v_x(*vlarge2, &vlarge2->get(sizex-sizex%binsrc[0], 0, 0), sizex%binsrc[0], sizey, sizez, vlarge2->getStrideX(),vlarge2->getStrideY(),vlarge2->getStrideZ());
                v_x.setValues(defaultval_);
                tom::Volume<T> v_y(*vlarge2, &vlarge2->get(0, sizey-sizey%binsrc[1], 0), sizex, sizey%binsrc[1], sizez, vlarge2->getStrideX(),vlarge2->getStrideY(),vlarge2->getStrideZ());
                v_y.setValues(defaultval_);
                tom::Volume<T> v_z(*vlarge2, &vlarge2->get(0, 0, sizez-sizez%binsrc[2]), sizex, sizey, sizez%binsrc[2], vlarge2->getStrideX(),vlarge2->getStrideY(),vlarge2->getStrideZ());
                v_z.setValues(defaultval_);
                v3_.reset(new tom::Volume<T>(*vlarge2, 0, sizex-sizex%binsrc[0], sizey-sizey%binsrc[1], sizez-sizez%binsrc[2], vlarge2->getStrideX(),vlarge2->getStrideY(),vlarge2->getStrideZ()));
                vlarge3 = v3_.get();
            }
            tom::transf::rescale_fs(vsrc, *vlarge3, ((binsrc[0]-1)/2.), ((binsrc[1]-1)/2.), ((binsrc[2]-1)/2.), fftw_flag);
        }
        if (!vlarge->is_equal_size(vdst)) {
            assert(vlarge != &vdst);
            tom::transf::bin(*vlarge, vdst);
        }
    }

}




namespace {
template<typename T, typename TBIN>
struct for_each__tom__transf__bin {
    for_each__tom__transf__bin(tom::Volume<TBIN> &vdst, const std::size_t bin[3])
        :   v_(vdst) {
        bin_[0] = bin[0];
        bin_[1] = bin[1];
        bin_[2] = bin[2];
    }
    void operator()(const T val, std::size_t x, std::size_t y, std::size_t z) const {
        assert(x/bin_[0]<v_.getSizeX() && y/bin_[1]<v_.getSizeY() && z/bin_[2]<v_.getSizeZ());
        v_.get_unsafe(x/bin_[0], y/bin_[1], z/bin_[2]) += val;
    }
private:
    tom::Volume<TBIN> &v_;
    std::size_t bin_[3];
};
}
/****************************************************************************//**
 * \brief Binning of a volume.
 *
 * The size of the output volume must correspond to the binned volume size.
 *******************************************************************************/
template<typename T, typename TBIN>
void tom::transf::bin(const tom::Volume<T> &vsrc, tom::Volume<TBIN> &vdst, std::size_t binx, std::size_t biny, std::size_t binz) {

    const std::size_t bin[3] = { binx?binx:1, biny?biny:1, binz?binz:1 };
    if (vdst.getSizeX()>vsrc.getSizeX() || vsrc.getSizeX()%vdst.getSizeX()>=bin[0] ||
        vdst.getSizeY()>vsrc.getSizeY() || vsrc.getSizeY()%vdst.getSizeY()>=bin[1] ||
        vdst.getSizeZ()>vsrc.getSizeZ() || vsrc.getSizeZ()%vdst.getSizeZ()>=bin[2]) {
        std::ostringstream ss;
        ss << "tom::transf::bin: source volume of size " << vsrc.getSizeX() << 'x' << vsrc.getSizeY() << 'x' << vsrc.getSizeZ() << " can not be binned to a volume of size " << vdst.getSizeX() << 'x' << vdst.getSizeY() << 'x' << vdst.getSizeZ() << ".";
        throw std::invalid_argument(ss.str());
    }
    if (bin[0]==1 && bin[1]==1 && bin[2]==1) {
        if (!vsrc.is_equal_memory(vdst)) {
            vdst.setValues(vsrc);
        }
        return;
    }

    const tom::Volume<T> *pvsrc = &vsrc;
    std::auto_ptr<tom::Volume<T> > vsrc_copy;
    if (vsrc.getSizeX()%bin[0] || vsrc.getSizeY()%bin[1] || vsrc.getSizeZ()%bin[2]) {
        vsrc_copy.reset(new tom::Volume<T>(const_cast<T *>(&vsrc.get()), vsrc.getSizeX()-vsrc.getSizeX()%bin[0], vsrc.getSizeY()-vsrc.getSizeY()%bin[1], vsrc.getSizeZ()-vsrc.getSizeZ()%bin[2], vsrc.getStrideX(), vsrc.getStrideY(), vsrc.getStrideZ(), false, 0));
        pvsrc = vsrc_copy.get();
    }
    vdst.setValues(0);
    tom::loop::for_each<const tom::Volume<T>, const for_each__tom__transf__bin<T, TBIN> &>(*pvsrc, for_each__tom__transf__bin<T, TBIN>(vdst, bin));
    vdst.template shift_scale<TBIN>(0, 1./(bin[0]*bin[1]*bin[2]));
}




/****************************************************************************//**
 * \brief Paste one volume into an other.
 *
 * \param[in] v1 The volume which is pasted into \c v2.
 * \param[out] v2 The destination volume.
 * \param[in] cornerx The first corner [0,0,0] of \c v1 is pasted to this coordinate.
 * \param[in] cornery The first corner [0,0,0] of \c v1 is pasted to this coordinate.
 * \param[in] cornerz
 * \param[in] outside_value If 0, the values in \c v2 which are not overwritten
 *  by \c v1 are left unchanged. Otherwise \c *outside_value is inserted.
 * The sizes of the volumes can be arbitrary (i.e. paste a smaller into a larger
 * volume is ok, just as the other way around).
 *******************************************************************************/
template<typename T1, typename T2>
void tom::transf::paste(const tom::Volume<T1> &v1, tom::Volume<T2> &v2, signed long cornerx, signed long cornery, signed long cornerz, const T2 *outside_value) {

    typedef signed long TIDX;

    TIDX v1_startx, v1_starty, v1_startz;
    TIDX v2_startx, v2_starty, v2_startz;
    TIDX sizex, sizey, sizez;

    v1_startx = std::min<signed long>(std::max<signed long>(0, -cornerx), v1.getSizeX());
    v1_starty = std::min<signed long>(std::max<signed long>(0, -cornery), v1.getSizeY());
    v1_startz = std::min<signed long>(std::max<signed long>(0, -cornerz), v1.getSizeZ());

    v2_startx = std::min<signed long>(std::max<signed long>(0, cornerx), v2.getSizeX());
    v2_starty = std::min<signed long>(std::max<signed long>(0, cornery), v2.getSizeY());
    v2_startz = std::min<signed long>(std::max<signed long>(0, cornerz), v2.getSizeZ());

    sizex = std::min<signed long>(v1.getSizeX()-v1_startx, v2.getSizeX()-v2_startx);
    sizey = std::min<signed long>(v1.getSizeY()-v1_starty, v2.getSizeY()-v2_starty);
    sizez = std::min<signed long>(v1.getSizeZ()-v1_startz, v2.getSizeZ()-v2_startz);

    if (sizex<=0 || sizey<=0 || sizez<=0) {
        if (outside_value) {
            v2.setValues(*outside_value);
        }
        return;
    }
    {
        tom::Volume<T1> v1_(const_cast<T1 *>(&v1.get(v1_startx, v1_starty, v1_startz)), sizex, sizey, sizez, v1.getStrideX(), v1.getStrideY(), v1.getStrideZ(), false, 0);
        tom::Volume<T2> v2_(v2, &v2.get(v2_startx, v2_starty, v2_startz), sizex, sizey, sizez, v2.getStrideX(), v2.getStrideY(), v2.getStrideZ());
        v2_.setValues(v1_);
    }
    if (outside_value) {
        T2 outside_value_ = *outside_value;
        if (v2_startx > 0) {
            tom::Volume<T2> v2_(v2, &v2.get(0,         0,         0), v2_startx              , v2.getSizeY()          , v2.getSizeZ(), v2.getStrideX(), v2.getStrideY(), v2.getStrideZ());
            v2_.setValues(outside_value_);
        }
        if (v2_starty > 0) {
            tom::Volume<T2> v2_(v2, &v2.get(v2_startx, 0,         0), v2.getSizeX()-v2_startx, v2_starty              , v2.getSizeZ(), v2.getStrideX(), v2.getStrideY(), v2.getStrideZ());
            v2_.setValues(outside_value_);
        }
        if (v2_startz > 0) {
            tom::Volume<T2> v2_(v2, &v2.get(v2_startx, v2_starty, 0), v2.getSizeX()-v2_startx, v2.getSizeY()-v2_starty, v2_startz    , v2.getStrideX(), v2.getStrideY(), v2.getStrideZ());
            v2_.setValues(outside_value_);
        }
        if (static_cast<std::size_t>(v2_startx+sizex)<v2.getSizeX()) {
            tom::Volume<T2> v2_(v2, &v2.get(v2_startx+sizex, v2_starty, v2_startz), v2.getSizeX()-v2_startx-sizex, v2.getSizeY()-v2_starty, v2.getSizeZ()-v2_startz, v2.getStrideX(), v2.getStrideY(), v2.getStrideZ());
            v2_.setValues(outside_value_);
        }
        if (static_cast<std::size_t>(v2_starty+sizey)<v2.getSizeY()) {
            tom::Volume<T2> v2_(v2, &v2.get(v2_startx, v2_starty+sizey, v2_startz), sizex, v2.getSizeY()-v2_starty-sizey, v2.getSizeZ()-v2_startz, v2.getStrideX(), v2.getStrideY(), v2.getStrideZ());
            v2_.setValues(outside_value_);
        }
        if (static_cast<std::size_t>(v2_startz+sizez)<v2.getSizeZ()) {
            tom::Volume<T2> v2_(v2, &v2.get(v2_startx, v2_starty, v2_startz+sizez), sizex, sizey, v2.getSizeZ()-v2_startz-sizez, v2.getStrideX(), v2.getStrideY(), v2.getStrideZ());
            v2_.setValues(outside_value_);
        }
    }
}






// Template instantiations.


template void tom::transf::transform<float , tom::transf::InterpolNearestNeighbour<float > >(const tom::Volume<float > &src, tom::Volume<float > &dst, const double *P, bool is_affinity, float  valinf, tom::transf::InterpolNearestNeighbour<float > interp);
template void tom::transf::transform<double, tom::transf::InterpolNearestNeighbour<double> >(const tom::Volume<double> &src, tom::Volume<double> &dst, const double *P, bool is_affinity, double valinf, tom::transf::InterpolNearestNeighbour<double> interp);
template void tom::transf::transform<float , tom::transf::InterpolTriLinear<float > >(const tom::Volume<float > &src, tom::Volume<float > &dst, const double *P, bool is_affinity, float  valinf, tom::transf::InterpolTriLinear<float > interp);
template void tom::transf::transform<double, tom::transf::InterpolTriLinear<double> >(const tom::Volume<double> &src, tom::Volume<double> &dst, const double *P, bool is_affinity, double valinf, tom::transf::InterpolTriLinear<double> interp);
template void tom::transf::transform<float , tom::transf::InterpolTriCubic<float > >(const tom::Volume<float > &src, tom::Volume<float > &dst, const double *P, bool is_affinity, float valinf, tom::transf::InterpolTriCubic<float > interp);
template void tom::transf::transform<double, tom::transf::InterpolTriCubic<double> >(const tom::Volume<double> &src, tom::Volume<double> &dst, const double *P, bool is_affinity, double valinf, tom::transf::InterpolTriCubic<double> interp);

template void tom::transf::transformSpline(const tom::Volume<float> &src, tom::Volume<float> &v_rot, const tom::Volume<float> &mtx, float defaultval);

template void tom::transf::transform<std::complex<float> , tom::transf::InterpolFourierSpline<std::complex<float> > >(const tom::Volume<std::complex<float> > &src, tom::Volume<std::complex<float> > &dst, const double *P, bool is_affinity, std::complex<float> valinf, tom::transf::InterpolFourierSpline<std::complex<float> > interp);


template void tom::transf::rotate(const tom::Volume<float > &src, tom::Volume<float > &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, double shiftx_pre, double shifty_pre, double shiftz_pre, double shiftx_post, double shifty_post, double shiftz_post, float  defaultval);
template void tom::transf::rotate(const tom::Volume<double> &src, tom::Volume<double> &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, double shiftx_pre, double shifty_pre, double shiftz_pre, double shiftx_post, double shifty_post, double shiftz_post, double defaultval);

template void tom::transf::rotateCubic(const tom::Volume<float > &src, tom::Volume<float > &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, double shiftx_pre, double shifty_pre, double shiftz_pre, double shiftx_post, double shifty_post, double shiftz_post, float  defaultval);
template void tom::transf::rotateCubic(const tom::Volume<double> &src, tom::Volume<double> &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, double shiftx_pre, double shifty_pre, double shiftz_pre, double shiftx_post, double shifty_post, double shiftz_post, double defaultval);

template void tom::transf::rotateSpline(const tom::Volume<float > &src, tom::Volume<float > &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, double shiftx_pre, double shifty_pre, double shiftz_pre, double shiftx_post, double shifty_post, double shiftz_post, float  defaultval);
template void tom::transf::rotateSpline(const tom::Volume<double> &src, tom::Volume<double> &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, double shiftx_pre, double shifty_pre, double shiftz_pre, double shiftx_post, double shifty_post, double shiftz_post, double defaultval);

//template void tom::transf::rotateFourierSpline(const tom::Volume<std::complex<float> > &src, tom::Volume<std::complex<float> > &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, double shiftx_pre, double shifty_pre, double shiftz_pre, double shiftx_post, double shifty_post, double shiftz_post, std::complex<float>  defaultval);
//template void tom::transf::rotateFourierSpline(const tom::Volume<std::complex<float> > &src, tom::Volume<std::complex<float> > &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, double shiftx_pre, double shifty_pre, double shiftz_pre, double shiftx_post, double shifty_post, double shiftz_post, std::complex<float> defaultval);

template void tom::transf::rotateFourierSpline(const tom::Volume<std::complex<float> > &src, tom::Volume<std::complex<float> > &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, std::complex<float> defaultval);

template void tom::transf::transformFourierSpline(const tom::Volume<std::complex<float>  > &src, tom::Volume<std::complex<float>  > &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, double shiftX, double shiftY, double shiftZ,std::complex<float>  defaultval);

//template void tom::transf::rotateFourierSpline(const tom::Volume<float > &src, tom::Volume<float > &v_rot, double phi, double psi, double theta, double centerx, double centery, double centerz, double shiftx_pre, double shifty_pre, double shiftz_pre, double shiftx_post, double shifty_post, double shiftz_post, float defaultval);
template void tom::transf::shiftFourier(const tom::Volume<std::complex<float> > &src, tom::Volume<std::complex<float> > &dst, double shift_x, double shift_y, double shift_z);


template void tom::transf::rescale_fs(const tom::Volume<float                > &vsrc,                        tom::Volume<float                > &vdst,                        double shiftx, double shifty, double shiftz, unsigned fftw_flag);
template void tom::transf::rescale_fs(const tom::Volume<double               > &vsrc,                        tom::Volume<double               > &vdst,                        double shiftx, double shifty, double shiftz, unsigned fftw_flag);
template void tom::transf::rescale_fs(const tom::Volume<float                > &vsrc,                        tom::Volume<std::complex<float > > &vdst, std::size_t sizez_dst, double shiftx, double shifty, double shiftz, unsigned fftw_flag);
template void tom::transf::rescale_fs(const tom::Volume<double               > &vsrc,                        tom::Volume<std::complex<double> > &vdst, std::size_t sizez_dst, double shiftx, double shifty, double shiftz, unsigned fftw_flag);
template void tom::transf::rescale_fs(const tom::Volume<std::complex<float > > &vsrc, std::size_t sizez_src, tom::Volume<float                > &vdst,                        double shiftx, double shifty, double shiftz, unsigned fftw_flag);
template void tom::transf::rescale_fs(const tom::Volume<std::complex<double> > &vsrc, std::size_t sizez_src, tom::Volume<double               > &vdst,                        double shiftx, double shifty, double shiftz, unsigned fftw_flag);
template void tom::transf::rescale_fs(const tom::Volume<std::complex<float > > &vsrc, std::size_t sizez_src, tom::Volume<std::complex<float > > &vdst, std::size_t sizez_dst, double shiftx, double shifty, double shiftz                    );
template void tom::transf::rescale_fs(const tom::Volume<std::complex<double> > &vsrc, std::size_t sizez_src, tom::Volume<std::complex<double> > &vdst, std::size_t sizez_dst, double shiftx, double shifty, double shiftz                    );

template void tom::transf::rescale_rebin<float >(const tom::Volume<float > &, std::size_t binsrcx, std::size_t binsrcy, std::size_t binsrcz, tom::Volume<float > &, std::size_t bindstx, std::size_t bindsty, std::size_t bindstz, std::size_t, std::size_t, std::size_t, const float  *default_val, unsigned fftw_flag);
template void tom::transf::rescale_rebin<double>(const tom::Volume<double> &, std::size_t binsrcx, std::size_t binsrcy, std::size_t binsrcz, tom::Volume<double> &, std::size_t bindstx, std::size_t bindsty, std::size_t bindstz, std::size_t, std::size_t, std::size_t, const double *default_val, unsigned fftw_flag);

template void tom::transf::bin<float , double>(const tom::Volume<float > &vsrc, tom::Volume<double> &vdst);
template void tom::transf::bin<double, double>(const tom::Volume<double> &vsrc, tom::Volume<double> &vdst);

template void tom::transf::paste<float , float >(const tom::Volume<float > &v1, tom::Volume<float > &v2, signed long cornerx, signed long cornery, signed long cornerz, const float  *outside_value);
template void tom::transf::paste<double, double>(const tom::Volume<double> &v1, tom::Volume<double> &v2, signed long cornerx, signed long cornery, signed long cornerz, const double *outside_value);

 
