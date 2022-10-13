/***********************************************************************//**
 * \file FreqWeight.cpp
 * \brief Contains the implementation of FreqWeight.hpp
 * \author  Thomas Haller
 * \version 0.3
 * \date    30.07.2010
 *
 * Contains the implementation of FreqWeight.hpp
 **************************************************************************/
#include <tom/FreqWeight.hpp>

#include <string>
#include <exception>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <math.h>

#include <tom/volume_fcn.hpp>
#include <tom/transf/transform.hpp>

/****************************************************************************//**
 * \brief Create a FreqWeight object from the description.
 *
 * Template method to construct the FreqWeight from its description.
 * The data type is determined by the template type.
 * The difference to getWeight is, that the returned pointer must be deleted
 * by the caller.
 *******************************************************************************/
template<typename T> tom::FreqWeight<T> *tom::FreqWeightDesc::createWeight(std::size_t sizex, std::size_t sizey, std::size_t sizez) const {
    FreqWeightDesc::data_type t;
    if (tom::is_double<T>()) {
        t = FreqWeightDesc::type_double;
    } else {
        // Should only be decleared for float and double.
        assert(tom::is_float<T>());
        t = FreqWeightDesc::type_float;
    }

    // There is not way to assure that the pointer returned is really of
    // the assumed type. The class implementor of createWeightV must guarantee that.
    return static_cast<tom::FreqWeight<T> *>(createWeightV(sizex, sizey, sizez, t));
}



/****************************************************************************//**
 * \brief Create a FreqWeight object from the description.
 *
 * Template method to construct the FreqWeight from its description.
 * The data type is determined by the template type.
 * The difference to createWeight is, that the returned pointer is still owned
 * by the FreqWeightDesc and must not be deleted. Successive calls of getWeight
 * don't create the object again, but return instead the same pointer as before.
 * This is usefull if several volumes share the same FreqWeightDesc. Using
 * getWeight to access the actuall weighting, allows to share also the weighting
 * object.
 *******************************************************************************/
template<typename T> tom::FreqWeight<T> &tom::FreqWeightDesc::getWeight(std::size_t sizex, std::size_t sizey, std::size_t sizez) {

    // There is not way to assure that the pointer returned is really of
    // the assumed type. The class implementor of createWeightV must guarantee that.
    if (tom::is_double<T>()) {
        if (getWeight_double_ && !getWeight_double_->is_equal_size(sizex, sizey, sizez)) {
            delete getWeight_double_;
            getWeight_double_ = 0;
        }
        if (!getWeight_double_) {
            getWeight_double_ = static_cast<tom::FreqWeight<double> *>(createWeightV(sizex, sizey, sizez, FreqWeightDesc::type_double));
        }
        return *reinterpret_cast<tom::FreqWeight<T> *>(getWeight_double_);
    }
    // Should only be decleared for float and double.
    assert(tom::is_float<T>());

    if (getWeight_float_ && !getWeight_float_->is_equal_size(sizex, sizey, sizez)) {
        delete getWeight_float_;
        getWeight_float_ = 0;
    }
    if (!getWeight_float_) {
        getWeight_float_ = static_cast<tom::FreqWeight<float> *>(createWeightV(sizex, sizey, sizez, FreqWeightDesc::type_float));
    }
    return *reinterpret_cast<tom::FreqWeight<T> *>(getWeight_float_);
}



/****************************************************************************//**
 * \brief Static constructor of descriptions from a text.
 *
 * Create a weighting description from a string (this format can be used in
 * configuration (text) files.
 *******************************************************************************/
tom::FreqWeightDesc *tom::FreqWeightDesc::fromString(const std::string s) {

    std::string s1, s2, s3, s4;
    {
        // Check that there is only one line.
        std::istringstream iss(s);
        if (!std::getline(iss, s1).eof()) {
            throw std::runtime_error("More than one line found (FreqWeightDesc::fromString).");
        }
    }
    std::istringstream iss(s);
    iss >> s1 >> std::ws;
    std::getline(iss, s2);
    assert(iss.eof());
    iss.clear();

    std::auto_ptr<tom::FreqWeightDesc> res;

    if (s1 == "" || s1=="allpass" || s1=="nowedge") {
        //return "allpass";
        if (!s2.empty()) {
            throw std::runtime_error((s1==""?std::string():s1) + " can have no parameters.");
        }
        res.reset(new tom::FreqWeightDesc_AllPass());
    } else if (s1 == "emfile") {
        //ss << "emfile  \"" << filename_ << "\"";
        bool stretch = false;
        bool resize = true;
        {
            iss.str(s2);
            iss >> std::ws >> s3;
            if (s3=="stretch" || s3=="nostretch") {
                stretch = s3=="stretch";
                iss >> std::ws;
                std::getline(iss, s2);
                assert(iss.eof());
                iss.clear();
                iss.str(s2);
                iss >> std::ws >> s3;
            }
            if (s3=="resize" || s3=="noresize") {
                resize = s3=="resize";
                iss >> std::ws;
                std::getline(iss, s2);
                assert(iss.eof());
                iss.clear();
            }
        }
        if (s2.size() > 2 && s2[0]=='"' && s2[s2.size()-1]=='"') {
            s2 = s2.substr(1, s2.size()-2);
        }
        if (s2.empty()) {
            throw std::runtime_error("emfile wedge expects the filename as parameter.");
        }
        res.reset(new tom::FreqWeightDesc_EmFile(s2, stretch,  resize));
    } else if (s1=="singleaxiswedge" || s1=="simple") {
        //ss << "singleaxiswedge" << (getAngle()*(180./3.141592653589793238512808959406186204433)) << " " << cutoff_radius_nobinning_;
        iss.str(s2);
        iss >> s2 >> std::ws >> s3 >> std::ws >> s4 >> std::ws;
        if (!iss.eof()) {
            throw std::runtime_error(s1 + " expects two to three parameters. The (half) opening angle in degree, the cut-off radius and optionally the corresponding Nyquist frequency of the x-dimension.");
        }
        double angle, cutoff_radius=0.;
        std::size_t nyquist_x = 0;
        try {
            angle = boost::lexical_cast<double>(s2);
            if (!s3.empty()) {
                cutoff_radius = boost::lexical_cast<double>(s3);
            }
            if (!s4.empty()) {
                nyquist_x = boost::lexical_cast<std::size_t>(s4);
            }
        } catch (boost::bad_lexical_cast &e) {
            throw std::runtime_error(s1 + " expects two or three arguments (i.e. the angle in DEG, the cutoff-radius and the corresponding Nyquist frequency along x).");
        }
        res.reset(new tom::FreqWeightDesc_SingleAxisWedge(tom::math::deg2rad(angle), cutoff_radius, nyquist_x));
    } else {
        throw std::runtime_error("Unknown description for weighting function (\"" + s + "\").");
    }
    return res.release();
}





/****************************************************************************//**
 * \brief Constructor of FreqWeight_SingleAxisWedge
 *
 * The angles must be in radians.
 * The cutoff_radius is absolute and relative to the volume size.
 * (e.g for sizex==32, cutoff_radius=15 would make a sphere with radius 15 in the
 * Fourier space of size [sizex,sizey,sizez].
 *******************************************************************************/
template<typename T>
tom::FreqWeight_SingleAxisWedge<T>::FreqWeight_SingleAxisWedge(std::size_t sizex, std::size_t sizey, std::size_t sizez, double angle, double cutoff_radius)
    :   FreqWeight<T>(sizex, sizey, sizez),
        w_redu_(NULL),
        w_full_(NULL),
        cphi_redu_(0.),
        cpsi_redu_(0.),
        cthe_redu_(0.),
        cphi_full_(0.),
        cpsi_full_(0.),
        cthe_full_(0.),
        w_full_valid_(false),
        angle_(angle),
        phi_(0.),
        psi_(0.),
        the_(0.),
        cutoff_radius_(cutoff_radius) {
    if (angle_ <= 0. || cutoff_radius_ < 0.) {
        throw std::invalid_argument("The SingleAxisWedge must have an opening angle > 0 and a cutoff_radius >= 0 (0 means no cut off).");
    }

}



/****************************************************************************//**
 * \brief Returns the volume for the weighting.
 *******************************************************************************/
template<typename T>
const tom::Volume<T> *tom::FreqWeight_SingleAxisWedge<T>::get_weight(bool reduced_complex) const {
    init_wedge_volume(reduced_complex);
    return reduced_complex ? w_redu_ : w_full_;
}


/****************************************************************************//**
 * Checks whether v is a complex volume of the proper size to be weighted.
 * v can be reduced complex. If it is, true is returned.
 *******************************************************************************/
template<typename T>
bool tom::FreqWeight<T>::is_reduced(const tom::Volume<typename FreqWeight<T>::Tcomplex> &v) const {
    const bool reduced_complex = v.getSizeZ() == this->sizez_/2+1;
    if (this->sizex_!=v.getSizeX() || this->sizey_!=v.getSizeY() || (!reduced_complex && this->sizez_!=v.getSizeZ())) {
        throw std::invalid_argument("The size of the volume does not match the weighting size");
    }
    return reduced_complex;
}




/****************************************************************************//**
 * \brief Applies the weighting to a complex volume.
 *******************************************************************************/
template<typename T>
bool tom::FreqWeight_SingleAxisWedge<T>::weight(tom::Volume<typename FreqWeight<T>::Tcomplex> &vsrc) const {
    const bool reduced_complex = this->is_reduced(vsrc);
    init_wedge_volume(reduced_complex);
    tom::element_wise_multiply<typename FreqWeight<T>::Tcomplex, T>(vsrc, *(reduced_complex ? w_redu_ : w_full_));
    return true;
}



/****************************************************************************//**
 * \brief Initializes the (rotated) weighting with the either zero or one (depending on the frequency).
 *******************************************************************************/
template<typename T>
void tom::FreqWeight_SingleAxisWedge<T>::init_wedge_volume(bool reduced) const {
    #define EQUAL_CROT(type) (cphi_##type##_==phi_ && cpsi_##type##_==psi_ && cthe_##type##_==the_)
    std::auto_ptr<tom::Volume<T> > w, w2;

    assert(!w_redu_ || w_redu_->is_equal_size(this->sizex_, this->sizey_, this->sizez_/2+1));
    assert(!w_full_ || w_full_->is_equal_size(this->sizex_, this->sizey_, this->sizez_    ));

    if (reduced) {
        if (!w_redu_) {
            if (w_full_) {
                w_redu_ = new tom::Volume<T>(*w_full_, NULL, this->sizex_, this->sizey_, this->sizez_/2+1, w_full_->getStrideX(), w_full_->getStrideY(), w_full_->getStrideZ());
                if (!EQUAL_CROT(full)) {
                    init_wedge_volume(*w_redu_);
                    cphi_full_ = phi_;
                    cpsi_full_ = psi_;
                    cthe_full_ = the_;
                    w_full_valid_ = false;
                }
            } else {
                w_redu_ = new tom::Volume<T>(this->sizex_, this->sizey_, this->sizez_/2+1, 0,0);
                init_wedge_volume(*w_redu_);
            }
        } else if (EQUAL_CROT(redu)) {
            // Reduced volume is already in the right orientation.
            return;
        } else {
            init_wedge_volume(*w_redu_);
            if (w_full_ && w_redu_->is_shared_memory(*w_full_)) {
                // Part of the full volume was overwritten.
                cphi_full_ = phi_;
                cpsi_full_ = psi_;
                cthe_full_ = the_;
                w_full_valid_ = false;
            }
        }
        cphi_redu_ = phi_;
        cpsi_redu_ = psi_;
        cthe_redu_ = the_;
    } else {
        if (!w_full_) {
            w_full_ = new tom::Volume<T>(this->sizex_, this->sizey_, this->sizez_, 0,0);
            std::auto_ptr<tom::Volume<T> > v2(new tom::Volume<T>(*w_full_, NULL, this->sizex_, this->sizey_, this->sizez_/2+1, w_full_->getStrideX(), w_full_->getStrideY(), w_full_->getStrideZ()));
            if (w_redu_ && EQUAL_CROT(redu)) {
                v2->setValues(*w_redu_);
                delete w_redu_;
                w_redu_ = v2.release();
            } else {
                if (w_redu_) {
                    delete w_redu_;
                }
                w_redu_ = v2.release();
                init_wedge_volume(*w_redu_);
            }
        } else {
            assert(w_redu_ && w_redu_->is_shared_memory(*w_full_));
            init_wedge_volume(*w_redu_);
        }
        ::tom::hermitian_symmetry_to_full<T>(*w_full_);
        cphi_redu_ = phi_;
        cpsi_redu_ = psi_;
        cthe_redu_ = the_;
        cphi_full_ = phi_;
        cpsi_full_ = psi_;
        cthe_full_ = the_;
        w_full_valid_ = true;
    }
    #undef EQUAL_CROT
}
#if defined(___TEST_PERFORMANCE_1_) || defined(___TEST_PERFORMANCE_2_)
#include <iostream>
#include <time.h>
#include <boost/lexical_cast.hpp>
/*******************************************************************************
 * A main routine to make a performance test whether it is faster to call
 * init_wedge_volume on the full volume of to call it on the reduced on and
 * use hermitian_symmetry_to_full to "copy" the values.
 * This is relevant in FreqWeight_SingleAxisWedge<T>::init_wedge_volume
 *
 * execute with
 *  alias X='rm .libs/libtomc.a .libs/tom__FreqWeight.o -f; CFLAGS_="-D___TEST_PERFORMANCE_${V}_" make && g++ .libs/libtomc.a  -lfftw3 -lfftw3f -lm && echo -e "\n\n-------------------" && ./a.out'
 *  export V=2; X
 ******************************************************************************/
typedef float TFLOAT;
std::complex<TFLOAT> rand1(std::complex<TFLOAT>) { return std::complex<TFLOAT>(std::rand(), std::rand()); }
int main(int argc, char **argv) {
    std::cout << "Make performance test for init_wedge_volume... (" << __FILE__ << ":" << __LINE__ << ")." << std::endl;
    const std::size_t size[3] = { 200, 200, 200};
    const std::size_t nrun = 15;

    std::cout.precision(20);

    tom::Volume<std::complex<TFLOAT> > v(size[0], size[1], size[2], 0,0);
    v.printInfo("vol");
    srand(5);
    tom::element_wise_operation(v, &::rand1);

    clock_t t1, t2;

    t1 = clock();
    for (std::size_t i=0; i<nrun; i++) {
        tom::FreqWeight_SingleAxisWedge<TFLOAT> w(tom::math::deg2rad(30.), size[0]/2.0);
        w.rotate(rand()*i*1.,rand()*3.*i,rand()*5.*i);
        w.weight(v, size[2]);
        #if 0
        w.get_weight(size[0], size[1], size[2], size[2])->write_to_em("./vol"
        #ifdef ___TEST_PERFORMANCE_1_
        "1"
        #else
        "2"
        #endif
        "_" + boost::lexical_cast<std::string>(i) + ".em", 0);
        #endif
    }
    t2 = clock();
    std::cout << "...Took " << (double(t2-t1)/CLOCKS_PER_SEC) << " s" << std::endl;
    return 0;
}
#endif




/****************************************************************************//**
 * \brief Initialise a volume with the right mask for a SingleAxisWedge.
 *
 * This method is static, to allow a usage from outside :)
 *******************************************************************************/
template<typename T>
void tom::FreqWeight_SingleAxisWedge<T>::init_wedge_volume(tom::Volume<T> &v) const {


    assert(this->sizex_ == v.getSizeX());
    assert(this->sizey_ == v.getSizeY());
    assert(this->sizez_ == v.getSizeZ() || this->sizez_/2+1 == v.getSizeZ());

    T *vdata = &v.get();
    const std::size_t v_sizez = v.getSizeZ();
    const std::size_t offsetx = this->sizex_/2 + this->sizex_%2;
    const std::size_t offsety = this->sizey_/2 + this->sizey_%2;
    const std::size_t offsetz = this->sizez_/2 + this->sizez_%2; // For the fftshift.

    assert(v.getStrideX()==sizeof(T) && v.getStrideY()==this->sizex_*v.getStrideX() && v.getStrideZ()==this->sizey_*v.getStrideY());

    typedef float TCOMP;

    TCOMP centerx = ceil((this->sizex_-1) / 2.);
    TCOMP centery = ceil((this->sizey_-1) / 2.);
    TCOMP centerz = ceil((this->sizez_-1) / 2.);

    TCOMP P[16] = {   0., 0., 0., 0.,
                      0., 0., 0., 0.,
                      0., 0., 0., 0.,
                      0., 0., 0., 1. };
    {
        double Plocal[16] = {   0., 0., 0., 0.,
                                0., 0., 0., 0.,
                                0., 0., 0., 0.,
                                0., 0., 0., 1. };
        const int axes[4] = {   0, 2, 0, 2 };
        const double angles[4] = { 0., psi_, the_, phi_ };
        // Only shift to the center, not back...
        const double shifts[3*4] = { -centerx, -centery, -centerz, 0,0,0, 0,0,0, 0,0,0 };
        tom::transf::sum_rotation(Plocal, true, false, 0, 0, 4, angles, axes, shifts);
        for (int i=0; i<16; i++) {
            P[i] = Plocal[i];
        }
    }
    const TCOMP maxdist = sqrt(centerx*centerx + centery*centery + centerz*centerz);

    const TCOMP tan_angle = tan(getAngleRad());

    const TCOMP cutoff_radius = cutoff_radius_>0 && cutoff_radius_<maxdist ? cutoff_radius_ : 0.;

    const TCOMP z0_threshold = 1e-4;

    std::size_t x, y, z;
    std::size_t xs, ys, zs;
    TCOMP x_, y_, z_;
    TCOMP x_tmpz, y_tmpz, z_tmpz;
    TCOMP x_tmpy, y_tmpy, z_tmpy;

    if (cutoff_radius > 0.) {
        const TCOMP cutoff_radius_squared = cutoff_radius*cutoff_radius;
        for (z=0; z<v_sizez; z++) {
            zs = (z + offsetz)%this->sizez_;
            x_tmpz = P[ 2]*zs + P[ 3];
            y_tmpz = P[ 6]*zs + P[ 7];
            z_tmpz = P[10]*zs + P[11];
            for (y=0; y<this->sizey_; y++) {
                ys = (y + offsety)%this->sizey_;
                x_tmpy = P[ 1]*ys + x_tmpz;
                y_tmpy = P[ 5]*ys + y_tmpz;
                z_tmpy = P[ 9]*ys + z_tmpz;
                for (x=0; x<this->sizex_; x++) {
                    xs = (x + offsetx)%this->sizex_;
                    x_ =                P[ 0]*xs + x_tmpy;
                    y_ =                P[ 4]*xs + y_tmpy;
                    z_ = tom::math::abs(P[ 8]*xs + z_tmpy);

                    *vdata++ = (x_*x_ + y_*y_ + z_*z_ <= cutoff_radius_squared) &&
                               ((z_<z0_threshold) || (tan_angle <= (tom::math::abs<TCOMP>(x_) / z_)));
                }
            }
        }
    } else {
        for (z=0; z<v_sizez; z++) {
            zs = (z + offsetz)%this->sizez_;
            x_tmpz = P[ 2]*zs + P[ 3];
            //y_tmpz = P[ 6]*zs + P[ 7];
            z_tmpz = P[10]*zs + P[11];
            for (y=0; y<this->sizey_; y++) {
                ys = (y + offsety)%this->sizey_;
                x_tmpy = P[ 1]*ys + x_tmpz;
                //y_tmpy = P[ 5]*ys + y_tmpz;
                z_tmpy = P[ 9]*ys + z_tmpz;
                for (x=0; x<this->sizex_; x++) {
                    xs = (x + offsetx)%this->sizex_;
                    x_ = P[ 0]*xs + x_tmpy;
                    //y_ = P[ 4]*xs + y_tmpy;
                    z_ = P[ 8]*xs + z_tmpy;
                    z_ = tom::math::abs<TCOMP>(z_);

                    *vdata++ = (z_<z0_threshold) || (tan_angle <= (tom::math::abs<TCOMP>(x_) / z_));
                }
            }
        }
    }
}


/****************************************************************************//**
 * \brief Constructor of FreqWeight_AsymetricSingleAxisWedge
 *
 * The angles must be in radians.
 * The cutoff_radius is absolute and relative to the volume size.
 * (e.g for sizex==32, cutoff_radius=15 would make a sphere with radius 15 in the
 * Fourier space of size [sizex,sizey,sizez].
 *******************************************************************************/
template<typename T>
tom::FreqWeight_AsymmetricSingleAxisWedge<T>::FreqWeight_AsymmetricSingleAxisWedge(std::size_t sizex, std::size_t sizey, std::size_t sizez, double angle1, double angle2, double cutoff_radius)
    :   FreqWeight<T>(sizex, sizey, sizez),
        w_redu_(NULL),
        w_full_(NULL),
        cphi_redu_(0.),
        cpsi_redu_(0.),
        cthe_redu_(0.),
        cphi_full_(0.),
        cpsi_full_(0.),
        cthe_full_(0.),
        w_full_valid_(false),
        angle1_(angle1),
        angle2_(angle2),
        phi_(0.),
        psi_(0.),
        the_(0.),
        cutoff_radius_(cutoff_radius) {
    if (angle1_ <= 0. || angle2_ <= 0.) {
        throw std::invalid_argument("The AssymetricSingleAxisWedge must have two opening angles > 0.");
    }

    if(cutoff_radius_ < 0. ){
    	throw std::invalid_argument("The AssymetricSingleAxisWedge must have a cutoff_radius >= 0 (0 means no cut off).");
    }

}



/****************************************************************************//**
 * \brief Returns the volume for the weighting.
 *******************************************************************************/
template<typename T>
const tom::Volume<T> *tom::FreqWeight_AsymmetricSingleAxisWedge<T>::get_weight(bool reduced_complex) const {
    init_wedge_volume(reduced_complex);
    return reduced_complex ? w_redu_ : w_full_;
}


/****************************************************************************//**
 * \brief Applies the weighting to a complex volume.
 *******************************************************************************/
template<typename T>
bool tom::FreqWeight_AsymmetricSingleAxisWedge<T>::weight(tom::Volume<typename FreqWeight<T>::Tcomplex> &vsrc) const {
    const bool reduced_complex = this->is_reduced(vsrc);
    init_wedge_volume(reduced_complex);
    tom::element_wise_multiply<typename FreqWeight<T>::Tcomplex, T>(vsrc, *(reduced_complex ? w_redu_ : w_full_));
    return true;
}



/****************************************************************************//**
 * \brief Initializes the (rotated) weighting with the either zero or one (depending on the frequency).
 *******************************************************************************/
template<typename T>
void tom::FreqWeight_AsymmetricSingleAxisWedge<T>::init_wedge_volume(bool reduced) const {
    #define EQUAL_CROT(type) (cphi_##type##_==phi_ && cpsi_##type##_==psi_ && cthe_##type##_==the_)
    std::auto_ptr<tom::Volume<T> > w, w2;

    assert(!w_redu_ || w_redu_->is_equal_size(this->sizex_, this->sizey_, this->sizez_/2+1));
    assert(!w_full_ || w_full_->is_equal_size(this->sizex_, this->sizey_, this->sizez_    ));

    if (reduced) {
        if (!w_redu_) {
            if (w_full_) {
                w_redu_ = new tom::Volume<T>(*w_full_, NULL, this->sizex_, this->sizey_, this->sizez_/2+1, w_full_->getStrideX(), w_full_->getStrideY(), w_full_->getStrideZ());
                if (!EQUAL_CROT(full)) {
                    init_wedge_volume(*w_redu_);
                    cphi_full_ = phi_;
                    cpsi_full_ = psi_;
                    cthe_full_ = the_;
                    w_full_valid_ = false;
                }
            } else {
                w_redu_ = new tom::Volume<T>(this->sizex_, this->sizey_, this->sizez_/2+1, 0,0);
                init_wedge_volume(*w_redu_);
            }
        } else if (EQUAL_CROT(redu)) {
            // Reduced volume is already in the right orientation.
            return;
        } else {
            init_wedge_volume(*w_redu_);
            if (w_full_ && w_redu_->is_shared_memory(*w_full_)) {
                // Part of the full volume was overwritten.
                cphi_full_ = phi_;
                cpsi_full_ = psi_;
                cthe_full_ = the_;
                w_full_valid_ = false;
            }
        }
        cphi_redu_ = phi_;
        cpsi_redu_ = psi_;
        cthe_redu_ = the_;
    } else {
        if (!w_full_) {
            w_full_ = new tom::Volume<T>(this->sizex_, this->sizey_, this->sizez_, 0,0);
            std::auto_ptr<tom::Volume<T> > v2(new tom::Volume<T>(*w_full_, NULL, this->sizex_, this->sizey_, this->sizez_/2+1, w_full_->getStrideX(), w_full_->getStrideY(), w_full_->getStrideZ()));
            if (w_redu_ && EQUAL_CROT(redu)) {
                v2->setValues(*w_redu_);
                delete w_redu_;
                w_redu_ = v2.release();
            } else {
                if (w_redu_) {
                    delete w_redu_;
                }
                w_redu_ = v2.release();
                init_wedge_volume(*w_redu_);
            }
        } else {
            assert(w_redu_ && w_redu_->is_shared_memory(*w_full_));
            init_wedge_volume(*w_redu_);
        }
        ::tom::hermitian_symmetry_to_full<T>(*w_full_);
        cphi_redu_ = phi_;
        cpsi_redu_ = psi_;
        cthe_redu_ = the_;
        cphi_full_ = phi_;
        cpsi_full_ = psi_;
        cthe_full_ = the_;
        w_full_valid_ = true;
    }
    #undef EQUAL_CROT
}



/****************************************************************************//**
 * \brief Initialise a volume with the right mask for a SingleAxisWedge.
 *
 * This method is static, to allow a usage from outside :)
 *******************************************************************************/
template<typename T>
void tom::FreqWeight_AsymmetricSingleAxisWedge<T>::init_wedge_volume(tom::Volume<T> &v) const {


    assert(this->sizex_ == v.getSizeX());
    assert(this->sizey_ == v.getSizeY());
    assert(this->sizez_ == v.getSizeZ() || this->sizez_/2+1 == v.getSizeZ());

    T *vdata = &v.get();
    const std::size_t v_sizez = v.getSizeZ();
    const std::size_t offsetx = this->sizex_/2 + this->sizex_%2;
    const std::size_t offsety = this->sizey_/2 + this->sizey_%2;
    const std::size_t offsetz = this->sizez_/2 + this->sizez_%2; // For the fftshift.

    assert(v.getStrideX()==sizeof(T) && v.getStrideY()==this->sizex_*v.getStrideX() && v.getStrideZ()==this->sizey_*v.getStrideY());

    typedef float TCOMP;

    TCOMP centerx = ceil((this->sizex_-1) / 2.);
    TCOMP centery = ceil((this->sizey_-1) / 2.);
    TCOMP centerz = ceil((this->sizez_-1) / 2.);

    TCOMP P[16] = {   0., 0., 0., 0.,
                      0., 0., 0., 0.,
                      0., 0., 0., 0.,
                      0., 0., 0., 1. };
    {
        double Plocal[16] = {   0., 0., 0., 0.,
                                0., 0., 0., 0.,
                                0., 0., 0., 0.,
                                0., 0., 0., 1. };
        const int axes[4] = {   0, 2, 0, 2 };
        const double angles[4] = { 0., psi_, the_, phi_ };
        // Only shift to the center, not back...
        const double shifts[3*4] = { -centerx, -centery, -centerz, 0,0,0, 0,0,0, 0,0,0 };
        tom::transf::sum_rotation(Plocal, true, false, 0, 0, 4, angles, axes, shifts);
        for (int i=0; i<16; i++) {
            P[i] = Plocal[i];
        }
    }
    const TCOMP maxdist = sqrt(centerx*centerx + centery*centery + centerz*centerz);

    const TCOMP tan_angle1 = tan(angle1_);
	const TCOMP tan_angle2 = tan(angle2_);

    const TCOMP cutoff_radius = cutoff_radius_>0 && cutoff_radius_<maxdist ? cutoff_radius_ : 0.;

    const TCOMP z0_threshold = 1e-4;

    std::size_t x, y, z;
    std::size_t xs, ys, zs;
    TCOMP x_, y_, z_;
    TCOMP x_tmpz, y_tmpz, z_tmpz;
    TCOMP x_tmpy, y_tmpy, z_tmpy;
    if (cutoff_radius > 0.) {
        const TCOMP cutoff_radius_squared = cutoff_radius*cutoff_radius;

        bool quad1,quad3;

        for (z=0; z<v_sizez; z++) {
            zs = (z + offsetz)%this->sizez_;
            x_tmpz = P[ 2]*zs + P[ 3];
            y_tmpz = P[ 6]*zs + P[ 7];
            z_tmpz = P[10]*zs + P[11];
            for (y=0; y<this->sizey_; y++) {
                ys = (y + offsety)%this->sizey_;
                x_tmpy = P[ 1]*ys + x_tmpz;
                y_tmpy = P[ 5]*ys + y_tmpz;
                z_tmpy = P[ 9]*ys + z_tmpz;
                for (x=0; x<this->sizex_; x++) {
                    xs = (x + offsetx)%this->sizex_;
                    x_ =                P[ 0]*xs + x_tmpy;
                    y_ =                P[ 4]*xs + y_tmpy;
                    z_ = tom::math::abs(P[ 8]*xs + z_tmpy);

                    //determine whether x is on the left or on the right part. (and lower quadrants)

                    quad1 = x_ >=0. && z_ >=0.;
                    quad3 = x_ <=0. && z_ <=0.;

                    //then use the according angle
                    if(quad1 || quad3){
                    	*vdata++ = (x_*x_ + y_*y_ + z_*z_ <= cutoff_radius_squared) &&
								   ((z_<z0_threshold) || (tan_angle1 <= (tom::math::abs<TCOMP>(x_) / z_)));
                    }else{
                    	*vdata++ = (x_*x_ + y_*y_ + z_*z_ <= cutoff_radius_squared) &&
                    	           ((z_<z0_threshold) || (tan_angle2 <= (tom::math::abs<TCOMP>(x_) / z_)));
                    }
                }
            }
        }
    } else {
    	bool quad1,quad3;

        for (z=0; z<v_sizez; z++) {
            zs = (z + offsetz)%this->sizez_;
            x_tmpz = P[ 2]*zs + P[ 3];
            //y_tmpz = P[ 6]*zs + P[ 7];
            z_tmpz = P[10]*zs + P[11];
            for (y=0; y<this->sizey_; y++) {
                ys = (y + offsety)%this->sizey_;
                x_tmpy = P[ 1]*ys + x_tmpz;
                //y_tmpy = P[ 5]*ys + y_tmpz;
                z_tmpy = P[ 9]*ys + z_tmpz;
                for (x=0; x<this->sizex_; x++) {
                    xs = (x + offsetx)%this->sizex_;
                    x_ = P[ 0]*xs + x_tmpy;
                    //y_ = P[ 4]*xs + y_tmpy;
                    z_ = P[ 8]*xs + z_tmpy;
                    z_ = tom::math::abs<TCOMP>(z_);

                    //determine whether x is on the left or on the right part. (and lower quadrants)

                    quad1 = x_ >=0. && z_ >=0.;
					quad3 = x_ <=0. && z_ <=0.;

					//then use the according angle
                    if(quad1 || quad3){
						*vdata++ = tan_angle1 <= (tom::math::abs<TCOMP>(x_) / z_);
					}else{
						*vdata++ = tan_angle2 <= (tom::math::abs<TCOMP>(x_) / z_);
					}

                }
            }
        }
    }
}



//----------------------------------------------------------------------------------
/****************************************************************************//**
 * \brief Constructor of FreqWeight_SmoothedAsymmetricSingleAxisWedge
 *
 * The angles must be in radians.
 * The cutoff_radius is absolute and relative to the volume size.
 *
 * (e.g for sizex==32, cutoff_radius=15 would make a sphere with radius 15 in the
 * Fourier space of size [sizex,sizey,sizez].
 *******************************************************************************/
template<typename T>
tom::FreqWeight_SmoothedAsymmetricSingleAxisWedge<T>::FreqWeight_SmoothedAsymmetricSingleAxisWedge(std::size_t sizex, std::size_t sizey, std::size_t sizez, double angle1, double angle2, double cutoff_radius, double smooth)
    :   FreqWeight<T>(sizex, sizey, sizez),
        w_redu_(NULL),
        w_full_(NULL),
        cphi_redu_(0.),
        cpsi_redu_(0.),
        cthe_redu_(0.),
        cphi_full_(0.),
        cpsi_full_(0.),
        cthe_full_(0.),
        w_full_valid_(false),
        angle1_(angle1),
        angle2_(angle2),
        phi_(0.),
        psi_(0.),
        the_(0.),
        cutoff_radius_(cutoff_radius),
        smooth_(smooth){
    if (angle1_ <= 0. || angle2_ <= 0.) {
        throw std::invalid_argument("The AssymetricSingleAxisWedge must have two opening angles > 0.");
    }

    if(cutoff_radius_ < 0. ){
    	throw std::invalid_argument("The AssymetricSingleAxisWedge must have a cutoff_radius >= 0 (0 means no cut off).");
    }

    if(smooth_ < 0. ){
    	throw std::invalid_argument("The AssymetricSingleAxisWedge must have a smooth >= 0 (0 means no smoothing).");
    }

}



/****************************************************************************//**
 * \brief Returns the volume for the weighting.
 *******************************************************************************/
template<typename T>
const tom::Volume<T> *tom::FreqWeight_SmoothedAsymmetricSingleAxisWedge<T>::get_weight(bool reduced_complex) const {
    init_wedge_volume(reduced_complex);
    return reduced_complex ? w_redu_ : w_full_;
}


/****************************************************************************//**
 * \brief Applies the weighting to a complex volume.
 *******************************************************************************/
template<typename T>
bool tom::FreqWeight_SmoothedAsymmetricSingleAxisWedge<T>::weight(tom::Volume<typename FreqWeight<T>::Tcomplex> &vsrc) const {
    const bool reduced_complex = this->is_reduced(vsrc);
    init_wedge_volume(reduced_complex);
    tom::element_wise_multiply<typename FreqWeight<T>::Tcomplex, T>(vsrc, *(reduced_complex ? w_redu_ : w_full_));
    return true;
}



/****************************************************************************//**
 * \brief Initializes the (rotated) weighting with the either zero or one (depending on the frequency).
 *******************************************************************************/
template<typename T>
void tom::FreqWeight_SmoothedAsymmetricSingleAxisWedge<T>::init_wedge_volume(bool reduced) const {
    #define EQUAL_CROT(type) (cphi_##type##_==phi_ && cpsi_##type##_==psi_ && cthe_##type##_==the_)
    std::auto_ptr<tom::Volume<T> > w, w2;

    assert(!w_redu_ || w_redu_->is_equal_size(this->sizex_, this->sizey_, this->sizez_/2+1));
    assert(!w_full_ || w_full_->is_equal_size(this->sizex_, this->sizey_, this->sizez_    ));

    if (reduced) {
        if (!w_redu_) {
            if (w_full_) {
                w_redu_ = new tom::Volume<T>(*w_full_, NULL, this->sizex_, this->sizey_, this->sizez_/2+1, w_full_->getStrideX(), w_full_->getStrideY(), w_full_->getStrideZ());
                if (!EQUAL_CROT(full)) {
                    init_wedge_volume(*w_redu_);
                    cphi_full_ = phi_;
                    cpsi_full_ = psi_;
                    cthe_full_ = the_;
                    w_full_valid_ = false;
                }
            } else {
                w_redu_ = new tom::Volume<T>(this->sizex_, this->sizey_, this->sizez_/2+1, 0,0);
                init_wedge_volume(*w_redu_);
            }
        } else if (EQUAL_CROT(redu)) {
            // Reduced volume is already in the right orientation.
            return;
        } else {
            init_wedge_volume(*w_redu_);
            if (w_full_ && w_redu_->is_shared_memory(*w_full_)) {
                // Part of the full volume was overwritten.
                cphi_full_ = phi_;
                cpsi_full_ = psi_;
                cthe_full_ = the_;
                w_full_valid_ = false;
            }
        }
        cphi_redu_ = phi_;
        cpsi_redu_ = psi_;
        cthe_redu_ = the_;
    } else {
        if (!w_full_) {
            w_full_ = new tom::Volume<T>(this->sizex_, this->sizey_, this->sizez_, 0,0);
            std::auto_ptr<tom::Volume<T> > v2(new tom::Volume<T>(*w_full_, NULL, this->sizex_, this->sizey_, this->sizez_/2+1, w_full_->getStrideX(), w_full_->getStrideY(), w_full_->getStrideZ()));
            if (w_redu_ && EQUAL_CROT(redu)) {
                v2->setValues(*w_redu_);
                delete w_redu_;
                w_redu_ = v2.release();
            } else {
                if (w_redu_) {
                    delete w_redu_;
                }
                w_redu_ = v2.release();
                init_wedge_volume(*w_redu_);
            }
        } else {
            assert(w_redu_ && w_redu_->is_shared_memory(*w_full_));
            init_wedge_volume(*w_redu_);
        }
        ::tom::hermitian_symmetry_to_full<T>(*w_full_);
        cphi_redu_ = phi_;
        cpsi_redu_ = psi_;
        cthe_redu_ = the_;
        cphi_full_ = phi_;
        cpsi_full_ = psi_;
        cthe_full_ = the_;
        w_full_valid_ = true;
    }
    #undef EQUAL_CROT
}



/****************************************************************************//**
 * \brief Initialise a volume with the right mask for a SingleAxisWedge.
 *
 * This method is static, to allow a usage from outside :)
 *******************************************************************************/
template<typename T>
void tom::FreqWeight_SmoothedAsymmetricSingleAxisWedge<T>::init_wedge_volume(tom::Volume<T> &v) const {

    assert(this->sizex_ == v.getSizeX());
    assert(this->sizey_ == v.getSizeY());
    assert(this->sizez_ == v.getSizeZ() || this->sizez_/2+1 == v.getSizeZ());

    T *vdata = &v.get();
    const std::size_t v_sizez = v.getSizeZ();
    const std::size_t offsetx = this->sizex_/2 + this->sizex_%2;
    const std::size_t offsety = this->sizey_/2 + this->sizey_%2;
    const std::size_t offsetz = this->sizez_/2 + this->sizez_%2; // For the fftshift.
//    std::cout << "Offset" << offsety << " " << offsety << " " << offsetz << std::endl;
    assert(v.getStrideX() == sizeof(T) && v.getStrideY() == this->sizex_*v.getStrideX() && v.getStrideZ() == this->sizey_*v.getStrideY());

    typedef float TCOMP;

    TCOMP centerx = ceil((this->sizex_-1) / 2.);
    TCOMP centery = ceil((this->sizey_-1) / 2.);
    TCOMP centerz = ceil((this->sizez_-1) / 2.);

    TCOMP P[16] = {   0., 0., 0., 0.,
                      0., 0., 0., 0.,
                      0., 0., 0., 0.,
                      0., 0., 0., 1. };
    {
        double Plocal[16] = {   0., 0., 0., 0.,
                                0., 0., 0., 0.,
                                0., 0., 0., 0.,
                                0., 0., 0., 1. };
        const int axes[4] = {   0, 2, 0, 2 };
        const double angles[4] = { 0., psi_, the_, phi_ };
        // Only shift to the center, not back...
        const double shifts[3*4] = { -centerx, -centery, -centerz, 0,0,0, 0,0,0, 0,0,0 };
        tom::transf::sum_rotation(Plocal, true, false, 0, 0, 4, angles, axes, shifts);
//        std::cout << "Matrix" << std::endl;
        for (int i=0; i<16; i++) {
            P[i] = Plocal[i];
//            std::cout << P[i] << " ";
//            if((i+1)%4 == 0 && (i+1)/4 >0)
//            	std::cout << std::endl;
        }
    }

    const TCOMP maxdist = sqrt(centerx*centerx + centery*centery + centerz*centerz);

    const TCOMP tan_angle1 = tan(angle1_);
	const TCOMP tan_angle2 = tan(angle2_);

//	const TCOMP tan_angle1Smooth = tan(angle1_ + smooth_);
//	const TCOMP tan_angle2Smooth = tan(angle2_ + smooth_);
	const TCOMP range_angle1Smooth = smooth_/sin(angle1_);
	const TCOMP range_angle2Smooth = smooth_/sin(angle2_);

	const TCOMP cutoff_radius = cutoff_radius_>0 && cutoff_radius_<maxdist ? cutoff_radius_ : 0.;
    // what is this parameter?
    const TCOMP z0_threshold = 1e-4;

    std::size_t x, y, z;
    std::size_t xs, ys, zs;
    TCOMP x_, y_, z_;
    TCOMP x_tmpz, y_tmpz, z_tmpz;
    TCOMP x_tmpy, y_tmpy, z_tmpy;
	
//	std::cout << "Cutoff: " << cutoff_radius << std::endl;
    if (cutoff_radius > 0.) {
        const TCOMP cutoff_radius_squared = cutoff_radius*cutoff_radius;
        bool quad1,quad3;

        for (z=0; z<v_sizez; z++) {
            zs = (z + offsetz)%this->sizez_;
            x_tmpz = P[ 2]*zs + P[ 3];
            y_tmpz = P[ 6]*zs + P[ 7];
            z_tmpz = P[10]*zs + P[11];
            for (y=0; y<this->sizey_; y++) {
                ys = (y + offsety)%this->sizey_;
                x_tmpy = P[ 1]*ys + x_tmpz;
                y_tmpy = P[ 5]*ys + y_tmpz;
                z_tmpy = P[ 9]*ys + z_tmpz;
                for (x=0; x<this->sizex_; x++) {
                    xs = (x + offsetx)%this->sizex_;
                    x_ =                P[ 0]*xs + x_tmpy;
                    y_ =                P[ 4]*xs + y_tmpy;
                    z_ = tom::math::abs(P[ 8]*xs + z_tmpy);

                    //determine whether x is on the left or on the right part. (and lower quadrants)
                    quad1 = x_ >=0. && z_ >=0.;
                    quad3 = x_ <=0. && z_ <=0.;
                    // changed FF - not clear why asymmetrical
                    //TH: I think it just came from the logical approach to avoid ambiguities, but x < 0 and z == 0 would not have been cought
                    //TH: That may have been it. Did you try your tests for the bug with this code?
                    
                    if(quad1 || quad3){
                    	if ((x_*x_ + y_*y_ + z_*z_ <= cutoff_radius_squared) && ((z_<z0_threshold) || (tan_angle1 <= (tom::math::abs<TCOMP>(x_) / z_)))){
							*vdata++ = 1;
                    	}
						else if( (x_*x_ + y_*y_ + z_*z_ <= cutoff_radius_squared) && ((z_<z0_threshold) || abs(z_-(x_/tan(angle1_)))<range_angle1Smooth) ){
							double dist = abs(z_-(x_/tan(angle1_)))*sin(angle1_)/smooth_;
							*vdata++ = 1 - dist;
						}
                    	else{
                    		*vdata++ = 0;
                    	}
                    }
                    else{
                    	if ((x_*x_ + y_*y_ + z_*z_ <= cutoff_radius_squared) && ((z_<z0_threshold) || (tan_angle2 <= (tom::math::abs<TCOMP>(x_) / z_)))){
							*vdata++ = 1;
                    	}
						else if( (x_*x_ + y_*y_ + z_*z_ <= cutoff_radius_squared) && ((z_<z0_threshold) || abs(z_+(x_/tan(angle2_)))<range_angle2Smooth) ){
							double dist = abs(z_+(x_/tan(angle2_)))*sin(angle2_)/smooth_;
							*vdata++ = 1 - dist;
						}
                    	else{
                    		*vdata++ = 0;
                    	}
					}
                }
            }
        }
    }
    else {
    	bool quad1,quad3;

        for (z=0; z<v_sizez; z++) {
            zs = (z + offsetz)%this->sizez_;
            x_tmpz = P[ 2]*zs + P[ 3];
            y_tmpz = P[ 6]*zs + P[ 7];
            z_tmpz = P[10]*zs + P[11];
            for (y=0; y<this->sizey_; y++) {
                ys = (y + offsety)%this->sizey_;
                x_tmpy = P[ 1]*ys + x_tmpz;
                y_tmpy = P[ 5]*ys + y_tmpz;
                z_tmpy = P[ 9]*ys + z_tmpz;
                for (x=0; x<this->sizex_; x++) {
                    xs = (x + offsetx)%this->sizex_;
                    x_ = P[ 0]*xs + x_tmpy;
                    y_ = P[ 4]*xs + y_tmpy;
                    z_ = P[ 8]*xs + z_tmpy;
                    z_ = tom::math::abs<TCOMP>(z_);

                    //determine whether x is on the left or on the right part. (and lower quadrants)
                    quad1 = x_ >=0. && z_ >=0.;
					quad3 = x_ <=0. && z_ <=0.;

					//then use the according angle
                    if(quad1 || quad3){
                    	if(tan_angle1 <= (tom::math::abs<TCOMP>(x_) / z_))
                    		*vdata++ = 1;
                    	else if(abs(z_-(x_/tan(angle1_)))<range_angle1Smooth){
                    		//smoothing applied as offset for wedge angle
                    		double dist = abs(z_-(x_/tan(angle1_)))*sin(angle1_)/smooth_;
							*vdata++ = 1 - dist;
                    	}
                    	else{
                    		*vdata++ = 0;
                    	}
					}else{
                    	if(tan_angle2 <= (tom::math::abs<TCOMP>(x_) / z_))
                    		*vdata++ = 1;
                    	else if(abs(z_+(x_/tan(angle2_)))<range_angle2Smooth){
                    		//smoothing applied as offset for wedge angle
                    		double dist = abs(z_+(x_/tan(angle2_)))*sin(angle2_)/smooth_;
							*vdata++ = 1 - dist;
                    	}
                    	else{
                    		*vdata++ = 0;
                    	}
					}
                }
            }
        }
    }
}
//----------------------------------------------------------------------------------

/****************************************************************************//**
 * \brief Constructor
 *
 * \param[in,out] w Volume to initialise the weighting. May be modified by
 *  the caller if \c copy is false.
 * \param[in] copy If true, the input volume is copied. Otherwise the weighting
 *  shares the volume with \c w. In this case, \c w is altered by the constructor
 *  and should not be altered anymore by the caller.
 *  However, the caller is allowed to delete \c w, because the constructor
 *  uses \c share_memory to access the same memory, it does not keep the
 *  valus of the object \c w itself.
 *******************************************************************************/
template<typename T>
tom::FreqWeight_Volume<T>::FreqWeight_Volume(tom::Volume<T> &w, bool copy):
    FreqWeight<T>(1, 1, 1),
    w_(0),
    w_full_(0),
    w_redu_(0),
    cphi_(0),
    cpsi_(0),
    cthe_(0),
    phi_(0),
    psi_(0),
    the_(0) {

    std::auto_ptr<tom::Volume<T> > ww;
    if (copy) {
        ww.reset(new tom::Volume<T>(w));
    } else {
        ww.reset(new tom::Volume<T>(w, 0, w.getSizeX(), w.getSizeY(), w.getSizeZ(), w.getStrideX(), w.getStrideY(), w.getStrideZ()));
    }
    this->resetVolume(ww);
}


/****************************************************************************//**
 * \brief Constructor
 *
 * \param[in] w Volume to initialise the weighting. The volume is never modified
 *   and a copy of it is made instead.
 *******************************************************************************/
template<typename T>
tom::FreqWeight_Volume<T>::FreqWeight_Volume(const tom::Volume<T> &w)
    :   FreqWeight<T>(1,1,1),
        w_(0),
        all_pass_(false),
        w_full_(0),
        w_redu_(0),
        cphi_(0),
        cpsi_(0),
        cthe_(0),
        phi_(0),
        psi_(0),
        the_(0) {

    std::auto_ptr<tom::Volume<T> > vv(new tom::Volume<T>(w));
    resetVolume(vv);
    //autov->template shift_scale<double>(0, 1./static_cast<double>(w_max));
}


/****************************************************************************//**
 * \brief Constructor
 *
 * \param[in] w Volume to initialise the weighting.
 *
 * The pointer in \c w is passed to the weighting and deleted by its destructor.
 *******************************************************************************/
template<typename T>
void tom::FreqWeight_Volume<T>::resetVolume(std::auto_ptr<tom::Volume<T> > &w) {

    delete w_;
    w_ = 0;
    delete w_full_;
    w_full_ = 0;
    delete w_redu_;
    w_redu_ = 0;
    all_pass_ = false;
    cphi_ = 0;
    cpsi_ = 0;
    cthe_ = 0;
    phi_ = 0;
    psi_ = 0;
    the_ = 0;

    this->sizex_ = w->getSizeX();
    this->sizey_ = w->getSizeY();
    this->sizez_ = w->getSizeZ();

    assert(w.get());

    if (!tom::isfinite(*w)) {
        throw std::invalid_argument("volume for frequency weighting contains non-finite values.");
    }
    T min, max;
    tom::minmax(*w, min, max);
    if (min < 0) {
        throw std::invalid_argument("volume for frequency weighting contains negative values (" + boost::lexical_cast<std::string>(min) + ").");
    }
    if (!(max > 0)) {
        throw std::invalid_argument("volume for frequency weighting contains no positive values.");
    }

    // Maybe fuzzy comparison would be better... :)
    all_pass_ = /*0 &&*/ (min==max) && (min==1);

    w_ = w.release();
    //w_->template shift_scale<double>(0, 1./static_cast<double>(max));
}



/****************************************************************************//**
 * \brief Applies the rotated weighting to a volume.
 *
 * Applies the rotated weighting to a volume.
 *******************************************************************************/
template<typename T>
bool tom::FreqWeight_Volume<T>::weight(tom::Volume<typename FreqWeight<T>::Tcomplex> &vsrc) const {

    assert(w_);

    const bool reduced_complex = this->is_reduced(vsrc);

    if (all_pass_) { return false; }

    init_wedge_volume();
    assert((reduced_complex && w_redu_) || (!reduced_complex && w_full_));

    if (reduced_complex) {
        tom::element_wise_multiply<typename FreqWeight<T>::Tcomplex, T>(vsrc, *w_redu_);
    } else {
        tom::element_wise_multiply<typename FreqWeight<T>::Tcomplex, T>(vsrc, *w_full_);
    }
    return true;
}



/****************************************************************************//**
 * \brief Returns the (rotated) weighting.
 *
 * \see tom::FreqWeight\<T\>::get_weight
 *******************************************************************************/
template<typename T>
const tom::Volume<T> *tom::FreqWeight_Volume<T>::get_weight(bool reduced_complex) const {
    assert(w_);
    init_wedge_volume();
    return reduced_complex ? w_redu_ : w_full_;
}



/****************************************************************************//**
 * \brief Protected default constructor.
 *
 * For derived classes to initialise the weighting half way in their
 * constructors.
 *******************************************************************************/
template<typename T>
tom::FreqWeight_Volume<T>::FreqWeight_Volume():
    FreqWeight<T>(1,1,1),
    w_(0),
    all_pass_(true),
    w_full_(),
    w_redu_(),
    cphi_(0),
    cpsi_(0),
    cthe_(0),
    phi_(0),
    psi_(0),
    the_(0) {
}



/****************************************************************************//**
 *
 *******************************************************************************/
template<typename T>
void tom::FreqWeight_Volume<T>::init_wedge_volume() const {

    assert(w_);

    if (w_full_ && cphi_==phi_ && cpsi_==psi_ && cthe_==the_) {

    } else {
        if (!w_full_) {
            assert(!w_redu_);
            w_full_ = new tom::Volume<T>(w_->getSizeX(), w_->getSizeY(), w_->getSizeZ(), 0, 0);
            w_redu_ = new tom::Volume<T>(*w_full_, NULL, w_full_->getSizeX(), w_full_->getSizeY(), w_full_->getSizeZ()/2+1, w_full_->getStrideX(), w_full_->getStrideY(), w_full_->getStrideZ());
        }

        if (phi_==0. && psi_==0. && the_==0.) {
            tom::fftshift(*w_, *w_full_, true);
        } else {
            tom::Volume<T> v(w_->getSizeX(), w_->getSizeY(), w_->getSizeZ(), NULL, NULL);
            tom::transf::rotate(*w_, v, phi_, psi_, the_);
            tom::fftshift(v, *w_full_, true);
        }

        cphi_ = phi_;
        cpsi_ = psi_;
        cthe_ = the_;

    }
}



/****************************************************************************//**
 * \brief Constructor
 *
 * Constructor. The emfile is fftshifted, i.e. it has the zero frequency at
 * ceil(center)
 *******************************************************************************/
template<typename T>
tom::FreqWeight_EmFile<T>::FreqWeight_EmFile(const std::string &filename, size_t sizex, size_t sizey, size_t sizez, bool stretch, bool  resize)
    :   FreqWeight_Volume<T>(),
        filename_(filename),
        stretch_(stretch),
        resize_(resize) {

    std::ostringstream ss;

    tom::io::VolumeEM v_em(filename);

    const bool sizes_differ = !v_em.is_equal_size(sizex, sizey, sizez);

    if (!resize_ && sizes_differ) {
        ss << "Error reading weighting from emfile \"" << filename << "\": size does not match (" << v_em->dims[0] << 'x' << v_em->dims[1] << 'x' << v_em->dims[2] << " vs. " << sizex << 'x' << sizey << 'x' << sizez << ')';
        throw std::runtime_error(ss.str());
    }
    this->sizex_ = sizex;
    this->sizey_ = sizey;
    this->sizez_ = sizez;
    std::auto_ptr<tom::Volume<T> > v;

    if (sizes_differ) {
        v.reset(new tom::Volume<T>(sizex, sizey, sizez, 0,0));
        std::auto_ptr<tom::Volume<T> > v0 = v_em.read_volume<T>(0, 0, 0);

        if (stretch) {
            assert(!"TODO: IMPLEMENT RESCALING OF WEIGHTING");
            throw "TODO: IMPLEMENT RESCALING OF WEIGHTING";
        } else {
            T outside_value = 0;
            tom::transf::paste(*v0, *v, &outside_value);
        }
    } else {
        v = v_em.read_volume<T>(0, 0, 0);
    }
    assert(v.get() && tom::isfinite(*v));

    try {
        this->resetVolume(v);
    } catch (std::invalid_argument &e) {
        throw std::invalid_argument("While reading weighting \"" + filename + "\": " + e.what());
    }
}




/****************************************************************************//**
 * \brief Default constructor.
 *
 * All values must be <=1 (positive only).
 *******************************************************************************/
template<typename T>
tom::FreqWeight_Bandpass<T>::FreqWeight_Bandpass(const float lowestFrequency, const float highestFrequency, const std::size_t x,const std::size_t y,const std::size_t z,const float smooth):
FreqWeight<T>(1,1,1),
b_(0),
all_pass_(false)
{

    this->lowestFrequency = lowestFrequency;
    this->highestFrequency = highestFrequency;
    this->size_x = x;
    this->size_y = y;
    this->size_z = z;
    this->smooth = smooth;
}


namespace {
template<typename T>
struct for_each_step__tom__init_bandpass {

    std::size_t posX,posY,posZ;

    std::size_t widthX,widthY,widthZ;

    std::size_t centerX,centerY,centerZ;

    std::size_t cubeSizeX,cubeSizeY,cubeSizeZ;

    float rx,ry,rz;

    T smooth;

    float lowestFrequency,highestFrequency;
    
    for_each_step__tom__init_bandpass(const float lowestFrequency, const float highestFrequency,const std::size_t wx,const std::size_t wy,const std::size_t wz, const T smooth){
        this->lowestFrequency = lowestFrequency;
        this->highestFrequency = highestFrequency;
        
        this->widthX =wx;
        this->centerX = wx/2;
        this->cubeSizeX = wx/2;


        if(wz == 1){

        	this->widthY =(wy-1)*2;
        	this->centerY = wy-1;
        	this->cubeSizeY = wy-1;

        	this->widthZ =1;
        	this->centerZ = 0;
        	this->cubeSizeZ = 0;

        }else{

        	this->widthY =wy;
        	this->centerY = wy/2;
        	this->cubeSizeY = wy/2;

        	this->widthZ =(wz-1)*2;
        	this->centerZ = wz-1;
        	this->cubeSizeZ = wz-1;

        }


        this->posX =0;
        this->posY =0;
        this->posZ =0;

        this->smooth = smooth;

    }

    void stepz(std::size_t z) {
        posZ = z;
        rz = (posZ+cubeSizeZ)%widthZ;
    }

    void stepy(std::size_t y) {
        posY = y;
        ry = (posY+cubeSizeY)%widthY;
    }

    void operator()(T &a, std::size_t x) {
        posX = x;

        rx = (posX+cubeSizeX)%widthX;
        //where am I relative to the small cubes (remember, we're in Fourier space, thats why we have to consider the shifted cubes)


        //radius = distance from the center
        float radius = sqrt((rx-centerX)*(rx-centerX)+(ry-centerY)*(ry-centerY)+(rz-centerZ)*(rz-centerZ));

        T val=0;

        //determine if (x,y,z) in any band
        if(this->is_in_band(radius)){
            a = 1;
        }
        //is the radius in the smoothing area?
        else if(this->smooth > 0 && this->is_in_smooth(radius,val)){
        	a = val;
        	val=0;
        }
        else{
            a = 0;
        }
    }

    bool is_in_smooth(const float radius,T & val){

    	//radius is lager than lowerBandEnd - smooth but smaller than lowerBandEnd; same below for upperEnd but inverse operators
        if(this->lowestFrequency - 3*this->smooth <= radius && radius < this->lowestFrequency){
            float dist = this->lowestFrequency - radius;
            T pos = ((T)dist)/this->smooth;
            //set value assessed by gaussian distribution
            val = exp(-pow(pos,2)/2);

            return true;
    	}
        else if(this->highestFrequency+3*this->smooth >= radius && radius > this->highestFrequency){
    		float dist = radius - this->highestFrequency;
    		T pos = ((T)dist)/this->smooth;
    		val = exp(-pow(pos,2)/2);

            return true;
        }
    	return false;
    }

    bool is_in_band(const float radius){
        
        if(this->lowestFrequency <= radius && radius < this->highestFrequency){
            return true;
        }
        
        return false;
    }

};
}

/****************************************************************************//**
 * \brief Initialises bandpass volume.
 *
 * If not set, the bandpass volume will be initialised.
 *******************************************************************************/
template<typename T>
void tom::FreqWeight_Bandpass<T>::init_bandpass_volume() const{

    if(this->b_ != NULL)
        return;

    this->b_ = new tom::Volume<T>(this->size_x,this->size_y,this->size_z,NULL,NULL);

    tom::loop::for_each_step(*this->b_, ::for_each_step__tom__init_bandpass<T>(this->lowestFrequency,this->highestFrequency,this->size_x,this->size_y,this->size_z,this->smooth));

}

/****************************************************************************//**
 * \brief Weigting function.
 *
 * * \see tom::FreqWeight\<T\>::weight
 *******************************************************************************/
template<typename T>
bool tom::FreqWeight_Bandpass<T>::weight(tom::Volume<typename FreqWeight<T>::Tcomplex> &vsrc) const{

    init_bandpass_volume();

    tom::element_wise_multiply<typename FreqWeight<T>::Tcomplex, T>(vsrc, *this->b_);

    return true;
}

/****************************************************************************//**
 * \brief Returns the weighting volume.
 *
 * * \see tom::FreqWeight\<T\>::get_weight
 *******************************************************************************/
template<typename T>
const tom::Volume<T> *tom::FreqWeight_Bandpass<T>::get_weight(bool b) const{
    if(this->b_ == NULL)
        init_bandpass_volume();
    return b_;
}

// template instantiations.
template class tom::FreqWeight<float >;
template class tom::FreqWeight<double>;
template class tom::FreqWeight_SingleAxisWedge<float >;
template class tom::FreqWeight_SingleAxisWedge<double>;
template class tom::FreqWeight_AsymmetricSingleAxisWedge<float >;
template class tom::FreqWeight_AsymmetricSingleAxisWedge<double>;
template class tom::FreqWeight_SmoothedAsymmetricSingleAxisWedge<float >;
template class tom::FreqWeight_SmoothedAsymmetricSingleAxisWedge<double>;
template class tom::FreqWeight_Volume<float >;
template class tom::FreqWeight_Volume<double>;
template class tom::FreqWeight_Bandpass<float >;
template class tom::FreqWeight_Bandpass<double>;
template class tom::FreqWeight_EmFile<float >;
template class tom::FreqWeight_EmFile<double>;

template tom::FreqWeight<float > *tom::FreqWeightDesc::createWeight(std::size_t sizex, std::size_t sizey, std::size_t sizez) const;
template tom::FreqWeight<double> *tom::FreqWeightDesc::createWeight(std::size_t sizex, std::size_t sizey, std::size_t sizez) const;
template tom::FreqWeight<float > &tom::FreqWeightDesc::getWeight(std::size_t sizex, std::size_t sizey, std::size_t sizez);
template tom::FreqWeight<double> &tom::FreqWeightDesc::getWeight(std::size_t sizex, std::size_t sizey, std::size_t sizez);


