/****************************************************************************//**
 * \file fftw_plan.cpp
 * \brief Contains implementations of the functions related to the fourier-transformation, defined in fftw_plan.hpp
 * \author  Thomas Haller
 * \version 0.2
 * \date    12.11.2007
 *******************************************************************************/
#include <tom/fftw/fftw_plan.hpp>



#include <cstring>
#include <cstdlib>
#include <cassert>
#include <exception>
#include <iostream>
#include <cctype>
#include <algorithm>
#include <typeinfo>
#include <unistd.h>

#include <boost/lexical_cast.hpp>
#include <boost/static_assert.hpp>

#include <tom/volume_fcn.hpp>



/****************************************************************************//**
 * \brief Variable for saving the wisdom-filename.
 *******************************************************************************/
char *wisdom_filenamef = 0;
char *wisdom_filename = 0;



namespace {
template<typename T>
class plan_destroyer {
public:
	void operator()(void *p);
};
template<>
class plan_destroyer<float > {
public:
	void operator()(void *p) {
		if (p) {
			fftwf_destroy_plan(reinterpret_cast<fftwf_plan>(p));
		}
	}
};
template<>
class plan_destroyer<double> {
public:
	void operator()(void *p) {
		if (p) {
			fftw_destroy_plan (reinterpret_cast<fftw_plan >(p));
		}
	}
};

}




/****************************************************************************//**
 *
 *******************************************************************************/
std::string tom::fftw::flag2str(unsigned flag) {
    flag = flag & (FFTW_ESTIMATE | FFTW_MEASURE | FFTW_PATIENT | FFTW_EXHAUSTIVE);
    if (flag == FFTW_ESTIMATE) {
        return "FFTW_ESTIMATE";
    } else if (flag == FFTW_MEASURE) {
        return "FFTW_MEASURE";
    } else if (flag == FFTW_PATIENT) {
        return "FFTW_PATIENT";
    } else if (flag == FFTW_EXHAUSTIVE) {
        return "FFTW_EXHAUSTIVE";
    } else {
        throw std::runtime_error("unrecognised fftw-flag \"" + boost::lexical_cast<std::string>(flag) + "\"");
    }
}


/****************************************************************************//**
 *
 *******************************************************************************/
unsigned tom::fftw::str2flag(const std::string &s) {
    std::string S(s);
    std::transform(S.begin(), S.end(), S.begin(), &toupper);
    if (S == "FFTW_PATIENT") {
        return FFTW_PATIENT;
    } else if (S == "FFTW_MEASURE") {
        return FFTW_MEASURE;
    } else if (S == "FFTW_ESTIMATE") {
        return FFTW_ESTIMATE;
    } else if (S == "FFTW_EXHAUSTIVE") {
        return FFTW_EXHAUSTIVE;
    } else {
        throw std::runtime_error("unrecognised string for fftw-flag \"" + s + "\"");
    }
}


/****************************************************************************//**
 * \brief Returns the currently set wisdom filename.
 *
 * Do not modify the resulting pointer.\n
 * These wisdom-functions are not thread save!
 *******************************************************************************/
template<typename T>
const char *tom::fftw::get_wisdom_name() {
    assert(tom::is_float<T>() || tom::is_double<T>());
    if (tom::is_float<T>()) {
        return wisdom_filenamef;
    }
    return wisdom_filename;
}


/****************************************************************************//**
 * \brief Load the wisdom from file and memorize its name in
 *
 * \param[in] fname The name of the wisdom file. Although this name is
 *   saved, the handled argument can be savely freed/modified afterwards.
 *   It is save to pass get_fftw_wisdom_name() as input (self-assignment).
 * \param[in] load If true, the wisdom is reloaded from the file. Otherwise
 *   its name is just remembered.
 * \return 1 in case of success, 0 otherwise. In case of failure,
 *   probably load was true, but the file contains invalid data.
 *   Or maybe (less likely), allocating memory for the copy of fname failed.
 *
 * For example you can set the name, and load it later be calling
 * \code
 * set_fftw_wisdom_name(name, 0);
 * ...
 * set_fftw_wisdom_name(get_fftw_wisdom_name(), 1);
 * \endcode
 *
 * These wisdom-functions are not thread save!
 *******************************************************************************/
template<typename T>
int tom::fftw::set_wisdom_name(const char *fname, int load) {
    assert(tom::is_float<T>() || tom::is_double<T>());

    int res = 0;
    char *wisdom_filename_local = NULL;


    if (fname) {
        wisdom_filename_local = (char *)malloc(sizeof(char)*strlen(fname)+1);
        if (!wisdom_filename_local) {
            return 0;
        }
        strcpy(wisdom_filename_local, fname);
    }
    tom::fftw::clear_wisdom_name<T>();
    if (tom::is_float<T>()) {
        wisdom_filenamef = wisdom_filename_local;
    } else {
        wisdom_filename = wisdom_filename_local;
    }
    if (wisdom_filename_local && load) {
        FILE *f = fopen(fname, "rb");
        if (f) {
            res = tom::is_float<T>() ? fftwf_import_wisdom_from_file(f) : fftw_import_wisdom_from_file(f);
            fclose(f);
        }
    }
    return res;
}


/****************************************************************************//**
 * \brief Saves the wisdom to the file with name get_fftw_wisdom_name()
 *
 * \return 1 in case of success, 0 otherwise. Failure can indicate, that
 *   the filename was not set with set_fftw_wisdom_name or the file could not
 *   be written.
 *
 * These wisdom-functions are not thread save!
 *******************************************************************************/
template<typename T>
int tom::fftw::save_wisdom() {
    assert(tom::is_float<T>() || tom::is_double<T>());
    int res = 0;
    char *filename = tom::is_float<T>() ? wisdom_filenamef : wisdom_filename;
    if (filename) {
        FILE *f = fopen(filename, "wb");
        if (f) {
            if (tom::is_float<T>()) {
                fftwf_export_wisdom_to_file(f);
            } else {
                fftw_export_wisdom_to_file(f);
            }
            res = !ferror(f);
            fclose(f);
        }
    }
    return res;
}



/****************************************************************************//**
 * \brief Saves the wisdom to the file with name get_fftw_wisdom_name()
 *
 * \return 1 in case of success, 0 otherwise. Failure can indicate, that
 *   the filename was not set with set_fftw_wisdom_name or the file could not
 *   be written.
 *
 * These wisdom-functions are not thread save!
 *******************************************************************************/
template<typename T>
void tom::fftw::clear_wisdom_name() {
    assert(tom::is_float<T>() || tom::is_double<T>());
    if (tom::is_float<T>()) {
        if (wisdom_filenamef) {
            free(wisdom_filenamef);
            wisdom_filenamef = NULL;
        }
    } else {
        if (wisdom_filename) {
            free(wisdom_filename);
            wisdom_filename = NULL;
        }
    }
}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
std::string tom::fftw::create_wisdom_name(const std::string &dir, const std::string &hostname) {
    assert(tom::is_float<T>() || tom::is_double<T>());
    if (dir.empty()) {
        return "";
    }
    std::stringstream ss;
    ss << dir;
    if (dir[dir.size()-1] != '/') {
        ss << '/';
    }
    ss << fftw_version << "_";
    if (hostname.empty()) {
        #ifndef HOST_NAME_MAX
        #define HOST_NAME_MAX 256
        #endif
        #define LEN___ (HOST_NAME_MAX>255 ? HOST_NAME_MAX : 256)
        char name[LEN___];
        gethostname(name, LEN___);
        name[LEN___-1] = 0;
        #undef LEN___
        ss << name;
    } else {
        ss << hostname;
    }
    ss << ".wisdom" << (tom::is_float<T>() ? "f" : "");

    return ss.str();
}






/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
int tom::fftw::Plan<T>::is_valid_dft_r2c(const tom::Volume<T> &vsrc, const tom::Volume<std::complex<T> > &vdst) {
    if (vsrc.getStrideX()%sizeof(T) == 0 &&
        vsrc.getStrideY()%sizeof(T) == 0 &&
        vsrc.getStrideZ()%sizeof(T) == 0 &&
        vdst.getStrideX()%sizeof(std::complex<T>) == 0 &&
        vdst.getStrideY()%sizeof(std::complex<T>) == 0 &&
        vdst.getStrideZ()%sizeof(std::complex<T>) == 0) {
        if (vsrc.getSizeX()      == vdst.getSizeX() &&
             vsrc.getSizeY()      == vdst.getSizeY() &&
            (vsrc.getSizeZ()/2+1) == vdst.getSizeZ()) {
            return TOM_FFTW_PLAN_3D;
        } else if (vsrc.getSizeX()      == vdst.getSizeX() &&
                   (vsrc.getSizeY()/2+1) == vdst.getSizeY() &&
                   vsrc.getSizeZ()==1 && vdst.getSizeZ()==1) {
            return TOM_FFTW_PLAN_2D;
        }
    }
    return TOM_FFTW_PLAN_INVALID;
}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
int tom::fftw::Plan<T>::is_valid_dft_c2r(const tom::Volume<std::complex<T> > &vsrc, const tom::Volume<T> &vdst) {
    if (vsrc.getStrideX()%sizeof(std::complex<T>) == 0 &&
        vsrc.getStrideY()%sizeof(std::complex<T>) == 0 &&
        vsrc.getStrideZ()%sizeof(std::complex<T>) == 0 &&
        vdst.getStrideX()%sizeof(T) == 0 &&
        vdst.getStrideY()%sizeof(T) == 0 &&
        vdst.getStrideZ()%sizeof(T) == 0) {
        if (vsrc.getSizeX() == vdst.getSizeX() &&
            vsrc.getSizeY() == vdst.getSizeY() &&
            vsrc.getSizeZ() == (vdst.getSizeZ()/2+1)) {
            return TOM_FFTW_PLAN_3D;
        } else if (vsrc.getSizeX() == vdst.getSizeX() &&
                   vsrc.getSizeY() == (vdst.getSizeY()/2+1) &&
                   vsrc.getSizeZ()==1 && vdst.getSizeZ()==1) {
            return TOM_FFTW_PLAN_2D;
        }
    }
    return TOM_FFTW_PLAN_INVALID;
}


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
tom::fftw::Plan<T>::Plan(const tom::fftw::Plan<T> &originalPlan){

	if(originalPlan.plan.get() == NULL)
		throw std::invalid_argument("tom::fftw::Plan(tom::fftw::Plan<T> &originalPlan) - The original plan is empty!");

	std::auto_ptr<tom::Volume<T               > > p_time0(new tom::Volume<T               >(*originalPlan.p_time0,NULL,originalPlan.p_time0->getSizeX(),originalPlan.p_time0->getSizeY(),originalPlan.p_time0->getSizeZ(),originalPlan.p_time0->getStrideX(),originalPlan.p_time0->getStrideY(),originalPlan.p_time0->getStrideZ()));

    std::auto_ptr<tom::Volume<std::complex<T> > > p_freq0(new tom::Volume<std::complex<T> >(*originalPlan.p_freq0,NULL,originalPlan.p_freq0->getSizeX(),originalPlan.p_freq0->getSizeY(),originalPlan.p_freq0->getSizeZ(),originalPlan.p_freq0->getStrideX(), originalPlan.p_freq0->getStrideY(), originalPlan.p_freq0->getStrideZ()));

	this->p_time0 = p_time0.release();
	this->p_freq0 = p_freq0.release();
	this->p_time1 = NULL;
	this->p_freq1 = NULL;

	this->type = originalPlan.type;
	this->plan = originalPlan.plan;
}


/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
tom::fftw::Plan<T>::Plan(tom::Volume<T> &vsrc, tom::Volume<std::complex<T> > &vdst, unsigned flags){

    const int type = this->is_valid_dft_r2c(vsrc, vdst);
    if (type == TOM_FFTW_PLAN_INVALID) {
        throw std::invalid_argument("tom::fftw::Plan() - The input volumes dont have not the right size (complex must be z/2+1 of real) or alignment for dft_r2c");
    }

    int rank = (type==TOM_FFTW_PLAN_3D ? 3 : 2);
    int howmany_rank = 0;
    fftw_iodim iodim[3];
    iodim[0].n  = vsrc.getSizeX();
    iodim[0].is = vsrc.getStrideX()/sizeof(T);
    iodim[0].os = vdst.getStrideX()/sizeof(std::complex<T>);
    iodim[1].n  = vsrc.getSizeY();
    iodim[1].is = vsrc.getStrideY()/sizeof(T);
    iodim[1].os = vdst.getStrideY()/sizeof(std::complex<T>);
    iodim[2].n  = vsrc.getSizeZ();
    iodim[2].is = vsrc.getStrideZ()/sizeof(T);
    iodim[2].os = vdst.getStrideZ()/sizeof(std::complex<T>);
    fftw_iodim *howmany_dims = NULL;


    boost::shared_ptr<void> p(static_cast<void *>(0), ::plan_destroyer<T>());
    void *pp;
    if (tom::is_float<T>()) {
        pp = static_cast<void *>(fftwf_plan_guru_dft_r2c(rank, iodim, howmany_rank, howmany_dims, (float  *)&vsrc.get(), (fftwf_complex *)&vdst.get(), flags));
    } else {
        pp = static_cast<void *>(fftw_plan_guru_dft_r2c (rank, iodim, howmany_rank, howmany_dims, (double *)&vsrc.get(), (fftw_complex  *)&vdst.get(), flags));
    }
    if (!pp) {
        throw std::runtime_error("tom::fftw::Plan() - could not create fftw plan dft_r2c. Check for fftw flags?");
    }
    p.reset(pp, ::plan_destroyer<T>());

    std::auto_ptr<tom::Volume<std::complex<T> > > p_freq0(new tom::Volume<std::complex<T> >(vdst, NULL, vdst.getSizeX(), vdst.getSizeY(), vdst.getSizeZ(), vdst.getStrideX(), vdst.getStrideY(), vdst.getStrideZ()));
    std::auto_ptr<tom::Volume<T               > > p_time0(new tom::Volume<T               >(vsrc, NULL, vsrc.getSizeX(), vsrc.getSizeY(), vsrc.getSizeZ(), vsrc.getStrideX(), vsrc.getStrideY(), vsrc.getStrideZ()));

    this->plan.swap(p);
    this->type = tom::fftw::Plan<T>::FFTW_DFT_R2C;
    this->p_freq0 = p_freq0.release();
    this->p_time0 = p_time0.release();
    this->p_freq1 = NULL;
    this->p_time1 = NULL;
}




/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
tom::fftw::Plan<T>::~Plan() {
    delete this->p_freq0;
    delete this->p_time0;
    delete this->p_freq1;
    delete this->p_time1;

}

/****************************************************************************//**
 * \brief Constructor of a C2R plan.
 *
 * \param[in,out] vsrc The source volume.
 * \param[in,out] vdst The destination volume.
 * \param[in] flags Flags passed directly to fftw_plan_guru_dft_c2r.
 *
 * Calles fftw_plan_guru_dft_c2r to create an fftw_plan. Depending on the flags
 * the input volumes are destroyed. See the documentation of fftw_plan_guru_dft_c2r.\n
 * Unfortunately there seems to be a bug in fftw3.1 concerning fftw_plan_guru_dft_c2r with
 * FFTW_PATIENT and FFTW_EXHAUSTIVE. The function writes outside the boundaries of
 * vdst. One sollution could be to initialize vdst with last dimension vsrc.getSizeZ()
 * large, instead of vsrc.getSizeZ()/2+1. If there is already wisdom about this type
 * of transformation, the error also don't happens.\n
 * It is not possible to use fftw_plan_dft_c2r_3d instead of the guru version because
 * fftw_plan_dft_c2r_3d expects the first dimension (vsrc.getSizeX()) to be reduced
 * instead of the last. Also the dimensions x,y,z are swapped. I think thats
 * an error in the FFTW-api, because its not consistent with the guru interface nor
 * the documentation.
 * With fftw3.2 this bug seams to be fixed.
 *******************************************************************************/
template<typename T>
tom::fftw::Plan<T>::Plan(tom::Volume<std::complex<T> > &vsrc, tom::Volume<T> &vdst, unsigned flags) {

    const int type = this->is_valid_dft_c2r(vsrc, vdst);
    if (type == TOM_FFTW_PLAN_INVALID) {
        throw std::invalid_argument("tom::fftw::Plan() - The input volumes have not the right size (complex must be z/2+1 of real) or alignment for dft_c2r");
    }

    const int rank = type==TOM_FFTW_PLAN_3D ? 3 : 2;
    int howmany_rank = 0;
    fftw_iodim iodim[3];
    iodim[0].n  = vdst.getSizeX();
    iodim[0].is = vsrc.getStrideX()/sizeof(std::complex<T>);
    iodim[0].os = vdst.getStrideX()/sizeof(T);
    iodim[1].n  = vdst.getSizeY();
    iodim[1].is = vsrc.getStrideY()/sizeof(std::complex<T>);
    iodim[1].os = vdst.getStrideY()/sizeof(T);
    iodim[2].n  = vdst.getSizeZ();
    iodim[2].is = vsrc.getStrideZ()/sizeof(std::complex<T>);
    iodim[2].os = vdst.getStrideZ()/sizeof(T);
    fftw_iodim *howmany_dims = NULL;
    std::auto_ptr<tom::Volume<std::complex<T> > > auto_vol_fourier;
    std::auto_ptr<tom::Volume<T               > > auto_vol_real;


    boost::shared_ptr<void> p(static_cast<void *>(0), ::plan_destroyer<T>());
    void *pp = 0;
    if (tom::is_float<T>()) {
        pp = fftwf_plan_guru_dft_c2r(rank, iodim, howmany_rank, howmany_dims, reinterpret_cast<fftwf_complex *>(&vsrc.get()), reinterpret_cast<float  *>(&vdst.get()), flags);
    } else {
        pp = fftw_plan_guru_dft_c2r (rank, iodim, howmany_rank, howmany_dims, reinterpret_cast<fftw_complex  *>(&vsrc.get()), reinterpret_cast<double *>(&vdst.get()), flags);
    }
    if (!pp) {
        throw std::runtime_error("tom::fftw::Plan() - could not create fftw plan dft_c2r. Check for fftw flags?");
    }
    p.reset(pp, ::plan_destroyer<T>());

    auto_vol_fourier.reset(new tom::Volume<std::complex<T> >(vsrc, NULL, vsrc.getSizeX(), vsrc.getSizeY(), vsrc.getSizeZ(), vsrc.getStrideX(), vsrc.getStrideY(), vsrc.getStrideZ()));
    auto_vol_real   .reset(new tom::Volume<T               >(vdst, NULL, vdst.getSizeX(), vdst.getSizeY(), vdst.getSizeZ(), vdst.getStrideX(), vdst.getStrideY(), vdst.getStrideZ()));

    this->plan.swap(p);
    this->type = tom::fftw::Plan<T>::FFTW_DFT_C2R;
    this->p_freq0 = auto_vol_fourier.release();
    this->p_time0 = auto_vol_real   .release();
    this->p_freq1 = NULL;
    this->p_time1 = NULL;
}





/****************************************************************************//**
 * \brief Execute the plan
 *
 * \param[in,out] vsrc The source volume.
 * \param[in,out] vdst The destination volume.
 *
 * It calles simply fftw_execute with the plan created upon instantiation.
 * The parameters are not acctually needed, its more for making clear what
 * is going to happen. i.e. you can not pass different volumes than at
 * construction time!
 *******************************************************************************/
template<typename T>
void tom::fftw::Plan<T>::execute(tom::Volume<T> &vsrc, tom::Volume<std::complex<T> > &vdst) const {

	assert(this->type!=tom::fftw::Plan<T>::FFTW_DFT_R2C || (this->p_time0&&this->p_freq0));
    if (this->type != tom::fftw::Plan<T>::FFTW_DFT_R2C ||
        !this->p_freq0->is_equal_memory(vdst) ||
        !this->p_time0->is_equal_memory(vsrc)) {
        throw std::runtime_error("tom::Volume::execute (R2C) - Call execute with the volumes used for construction of the plan.");
    }
    if (tom::is_float<T>()) {
        fftwf_execute((fftwf_plan)this->plan.get());
    } else {
        fftw_execute ((fftw_plan )this->plan.get());
    }
}

/****************************************************************************//**
 * \brief Execute the plan
 *
 * \param[in,out] vsrc The source volume.
 * \param[in,out] vdst The destination volume.
 *
 * It calles simply fftw_execute with the plan created upon instantiation.
 * The parameters are not acctually needed, its more for making clear what
 * is going to happen. i.e. you can not pass different volumes than at
 * construction time!
 *******************************************************************************/
template<typename T>
void tom::fftw::Plan<T>::execute(tom::Volume<std::complex<T> > &vsrc, tom::Volume<T> &vdst) const {

    if (this->type != tom::fftw::Plan<T>::FFTW_DFT_C2R ||
        !this->p_freq0->is_equal_memory(vsrc) ||
        !this->p_time0->is_equal_memory(vdst)) {
        throw std::runtime_error("tom::Volume::execute (C2R) - Call execute with the volumes used for construction of the plan.");
    }
    if (tom::is_float<T>()) {
        fftwf_execute(reinterpret_cast<fftwf_plan>(this->plan.get()));
    } else {
        fftw_execute (reinterpret_cast<fftw_plan >(this->plan.get()));
    }
}


/****************************************************************************//**
 * \brief Get access to the source volume.
 *
 * Template method because depending on the Type, either a complex or a real
 * valued volume is returned. If the wrong templated function is called, and
 * std::runtime_error is thrown.
 *
 * This volume is not the same object as when created the plan.
 * However, they use the same memory.
 *******************************************************************************/
template<typename T> template<typename T2>
tom::Volume<T2> &tom::fftw::Plan<T>::getSrcVol() {
    switch (type) {
        case Plan<T>::FFTW_DFT_R2C:
            if (!tom::is_double<T2>() && !tom::is_float<T2>()) {
                throw std::runtime_error("the fftw-plan contains a real valued source volume (R2C).");
            }
            assert(p_time0 && typeid(*p_time0) == typeid(tom::Volume<T2>));
            return *reinterpret_cast<tom::Volume<T2> *>(p_time0);
        case Plan<T>::FFTW_DFT_C2R:
            break;
        case Plan<T>::FFTW_DFT_FORWARD:
        case Plan<T>::FFTW_DFT_BACKWARD:
        default:
            assert(0);
    }
    assert(type == Plan<T>::FFTW_DFT_C2R);
    if (!tom::is_double_complex<T2>() && !tom::is_float_complex<T2>()) {
        throw std::runtime_error("the fftw-plan contains a complex valued source volume (C2R).");
    }
    assert(p_freq0 && typeid(*p_freq0) == typeid(tom::Volume<T2>));
    return *reinterpret_cast<tom::Volume<T2> *>(p_freq0);
 }

/****************************************************************************//**
 * \brief Get access to the destination volume.
 *
 * Template method because depending on the Type, either a complex or a real
 * valued volume is returned. If the wrong templated function is called, and
 * std::runtime_error is thrown.
 *
 * This volume is not the same object as when created the plan.
 * However, they use the same memory.
 *******************************************************************************/
template<typename T> template<typename T2>
tom::Volume<T2> &tom::fftw::Plan<T>::getDstVol() {
    switch (type) {
        case Plan<T>::FFTW_DFT_R2C:
            if (!tom::is_double_complex<T2>() && !tom::is_float_complex<T2>()) {
                throw std::runtime_error("the fftw-plan contains a complex valued destination volume (R2C).");
            }
            assert(p_freq0 && typeid(*p_freq0) == typeid(tom::Volume<T2>));
            return *reinterpret_cast<tom::Volume<T2> *>(p_freq0);
        case Plan<T>::FFTW_DFT_C2R:
            break;
        case Plan<T>::FFTW_DFT_FORWARD:
        case Plan<T>::FFTW_DFT_BACKWARD:
        default:
            assert(0);
    }
    assert(type == Plan<T>::FFTW_DFT_C2R);
    if (!tom::is_double<T2>() && !tom::is_float<T2>()) {
        throw std::runtime_error("the fftw-plan contains a real valued destination volume (C2R).");
    }
    assert(p_time0 && typeid(*p_time0) == typeid(tom::Volume<T2>));
    return *reinterpret_cast<tom::Volume<T2> *>(p_time0);
 }




template class tom::fftw::Plan<float >;
template class tom::fftw::Plan<double>;

template tom::Volume<float                > &tom::fftw::Plan<float >::getSrcVol();
template tom::Volume<std::complex<float > > &tom::fftw::Plan<float >::getSrcVol();
template tom::Volume<double               > &tom::fftw::Plan<double>::getSrcVol();
template tom::Volume<std::complex<double> > &tom::fftw::Plan<double>::getSrcVol();
template tom::Volume<float                > &tom::fftw::Plan<float >::getDstVol();
template tom::Volume<std::complex<float > > &tom::fftw::Plan<float >::getDstVol();
template tom::Volume<double               > &tom::fftw::Plan<double>::getDstVol();
template tom::Volume<std::complex<double> > &tom::fftw::Plan<double>::getDstVol();


template const char *tom::fftw::get_wisdom_name<float >();
template const char *tom::fftw::get_wisdom_name<double>();
template int         tom::fftw::set_wisdom_name<float >(const char *fname, int load);
template int         tom::fftw::set_wisdom_name<double>(const char *fname, int load);
template void        tom::fftw::clear_wisdom_name<float >();
template void        tom::fftw::clear_wisdom_name<double>();
template int         tom::fftw::save_wisdom<float >();
template int         tom::fftw::save_wisdom<double>();
template std::string tom::fftw::create_wisdom_name<float >(const std::string &dir, const std::string &hostname);
template std::string tom::fftw::create_wisdom_name<double>(const std::string &dir, const std::string &hostname);










