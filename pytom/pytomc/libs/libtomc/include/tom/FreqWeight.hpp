/****************************************************************************//**
 * \file FreqWeight.hpp
 * \brief The header file for wedge FreqWeight.
 * \author  Thomas Haller
 * \version 0.3
 * \date    30.07.2010
 *
 * FreqWeight is an abstract base class for all forms of weighting in Fourier space.
 * The missing wedge is a specialized frequency weighting and there are several types
 * implemented.
 *******************************************************************************/
#ifndef ___INCLUDE__TOM__FREQWEIGTH_HPP__
#define ___INCLUDE__TOM__FREQWEIGTH_HPP__




#include <tom/volume_fcn.hpp>
#include <vector>

#ifdef ___TEST_PERFORMANCE_1_
// forward declaration for a test (main) routine decleared in FreqWeight.cpp
int main(int argc, char **argv);
#endif



namespace tom {


typedef enum { identifier_NULL, identifier_AllPass, identifier_SingleAxisWedge, identifier_EmFile } Identifier_FreqWeight;



/****************************************************************************//**
 * \brief Abstract base class for all types of frequency weighting.
 *******************************************************************************/
template <typename T>
class FreqWeight {

public:

    /***********************************************************************//**
     * \brief Typedef of the element type (will be either float or double)
     *
     * Typedef of the real element type (will be either float or double).
     **************************************************************************/
    typedef T Treal;

    /***********************************************************************//**
     * \brief Typedef of the complex element type
     *
     * Typedef of the complex element type.
     **************************************************************************/
    typedef std::complex<T> Tcomplex;

    /// Destructor. No operation in the base class.
    virtual ~FreqWeight() { };

    /***********************************************************************//**
     * \brief Rotates the weighting.
     * \param[in] phi Rotation angle in radians!!
     * \param[in] psi Rotation angle in radians!!
     * \param[in] theta Rotation angle in radians!!
     *
     * Abstract method that rotates the weighting.
     * The class tom::FreqWeight and their descended classes support the weighting
     * of volumes rotated about their centre <CODE>int(size/2)</CODE>. \n
     * When a volume is rotated in real space about its centre (using
     * tom::transf::rotate), its frequencies in Fourier space are also rotated
     * about the zero-frequency. Moreover there occures a phase shift, which however
     * is not relevant for the weighting, as it only weights the amplitudes, ignoring
     * the phase.\n
     * Calling \c rotate on the FreqWeight object brings the weighting in the same
     * orientation as the volume, which can afterwards be weigted calling \c weight. \n
     * Some descended class might implement a rotation invariant filter and therefore
     * perform no operations when calling rotate (e.g. a low-pass filter).
     **************************************************************************/
    virtual void rotate(double phi, double psi, double theta) = 0;

    /***********************************************************************//**
     * \brief If true, the weighting does not actually weight any frequency.
     *
     * Calling all_pass allows to determine whether a call af \c weight actually
     * changes any frequencies in the volume. This might depend on the
     * current orientation (as set by calling \c rotate).
     **************************************************************************/
    virtual bool all_pass() const = 0;

    /****************************************************************************//**
     * \brief Get a volume with the coefficients of the weighting.
     *
     * \param[in] reduced_complex Return only the upper half of the weighting.
     * \returns A volume of size sizex*sizey*sizez_reduced_or_full with the weightings.
     *      Do not delete this volume and don't use it after calling an other method
     *      of the object. If 0 is returned, this corresponds to an all-pass filter
     *      (a volume filled with ones everywhere).
     *
     * The output volume should not be modified (or deleted) and is only valid
     * until \c this gets destructed or an other method is called.
     *
     * The volume is not Fourier shifted and has the zero frequency at position [0,0,0]
     *******************************************************************************/
    virtual const tom::Volume<T> *get_weight(bool reduced_complex) const = 0;


    /****************************************************************************//**
     * \brief Applay the weighting to the volume in fourier space.
     *
     * \param[in,out] vsrc The Fourier transformed volume of complex elements.
     *   The volume must have the same size as the weighting, except SizeZ, which
     *   can be SizeZ/2+1 for reduced complex volumes.
     * \return A boolean volume, specifying whether vsrc was changed by the wedge.
     *   This can be usefull, to decide whether a normalisation done before is still
     *   valid.
     *******************************************************************************/
    virtual bool weight(tom::Volume<typename FreqWeight<T>::Tcomplex> &vsrc) const = 0;

    std::size_t getSizeX() const { return sizex_; }
    std::size_t getSizeY() const { return sizey_; }
    std::size_t getSizeZ() const { return sizez_; }

    bool is_equal_size(std::size_t sizex, std::size_t sizey, std::size_t sizez) const {
        return sizex == sizex_ && sizey == sizey_ && sizez == sizez_;
    }

    bool is_reduced(const tom::Volume<typename FreqWeight<T>::Tcomplex> &v) const;


protected:
    FreqWeight(std::size_t sizex, std::size_t sizey, std::size_t sizez) {
        if (sizex <= 0 || sizey <= 0 || sizez <= 0) {
            throw std::invalid_argument("weighing must have a positive volume size.");
        }
        sizex_ = sizex;
        sizey_ = sizey;
        sizez_ = sizez;
    };
    std::size_t sizex_, sizey_, sizez_;

private:
    FreqWeight(const FreqWeight<T> &f)               { assert(!"HIDDEN");               }
    FreqWeight<T> &operator=(const FreqWeight<T> &f) { assert(!"HIDDEN"); return *this; }

};





/***************************************************************************//**
 * \brief A class to handle the description of a FreqWeight object.
 *
 * From this class an actual FreqWeight object can be created.
 * The description can be generated from a configuration string and represented
 * as a string.
 ******************************************************************************/
class FreqWeightDesc {
public:
    /// Constructor.
    FreqWeightDesc()
        :   getWeight_float_(0),
            getWeight_double_(0) {
    }

    template<typename T> tom::FreqWeight<T> *createWeight(std::size_t sizex, std::size_t sizey, std::size_t sizez) const;
    template<typename T> tom::FreqWeight<T> &getWeight(std::size_t sizex, std::size_t sizey, std::size_t sizez);

    /// Virtual desctructor.
    virtual ~FreqWeightDesc() {
        clearWeight();
    };

    /// Releases internal weightings allocated by getWeight.
    void clearWeight() {
        delete getWeight_double_;
        getWeight_double_ = 0;
        delete getWeight_float_;
        getWeight_float_ = 0;
    };

    /// Returns the line from the configuration file which describes the weighting.
    virtual std::string toString() const = 0;

    /// Static function to initialise the weighting from a string.
    static FreqWeightDesc *fromString(const std::string s);

    /// Returns a unique Id depending on the type, starting at 0.
    virtual Identifier_FreqWeight getIdentifier() const = 0;
protected:

    tom::FreqWeight<float > *getWeight_float_;
    tom::FreqWeight<double> *getWeight_double_;


    /// constant to describe the data type (single or double).
    typedef enum { type_float, type_double } data_type;

    /***********************************************************************//**
    * \brief Called by createWeight to have a (abstract) virtual function to construct
    *   the weighting with the proper data type.
    * \param[in] sizex Size of the created weighting.
    * \param[in] sizey Size of the created weighting.
    * \param[in] sizez Size of the created weighting.
    * \param[in] t Determines the data type for which the weighting is constructed.
    * \returns A pointer to the weigting. The implementor of the descended class must
    *   take care to returns either a pointer to \c tom::FreqWeight\<double\> or
    *   to \c tom::FreqWeight\<float\>.
    ***************************************************************************/
    virtual void *createWeightV(std::size_t sizex, std::size_t sizey, std::size_t sizez, FreqWeightDesc::data_type t) const = 0;

};






/****************************************************************************//**
 * \brief Dummy weighting which is an all-pass filter.
 *******************************************************************************/
template <typename T>
class FreqWeight_AllPass: public FreqWeight<T> {

public:

    /// Constructor.
    FreqWeight_AllPass(std::size_t sizex, std::size_t sizey, std::size_t sizez)
        : FreqWeight<T>(sizex, sizey, sizez) {
        /* NOP */
    }

    /// Virtual destructor (no operation).
    virtual ~FreqWeight_AllPass() { /* NOP */ }

    /// Rotate the volume (no operation).
    virtual void rotate(double phi, double psi, double theta) { /* NOP */ }


    /// This class is allways an all pass.
    virtual bool all_pass() const {
        return true;
    }

    /// Returns always 0.
    virtual const tom::Volume<T> *get_weight(bool reduced_complex) const {
        return 0;
    }

    /// Applies the weighting (no operation).
    virtual bool weight(tom::Volume<typename FreqWeight<T>::Tcomplex> &vsrc) const {
        // throw an exception if the volume sizes do not match.
        this->is_reduced(vsrc);
        return false;
    }
};





/****************************************************************************//**
 * \brief Descriptor for FreqWeightDesc.
 *******************************************************************************/
class FreqWeightDesc_AllPass: public FreqWeightDesc {
public:
    /// Constructor (no operation)
    FreqWeightDesc_AllPass(): FreqWeightDesc() { }

    /// Destructor (no operation)
    virtual ~FreqWeightDesc_AllPass() { };

    /// String to describe the weighting.
    virtual std::string toString() const {
        return "allpass";
    }

    /// Returns a unique Id depending on the type, starting at 0.
    virtual Identifier_FreqWeight getIdentifier() const { return identifier_AllPass; };

protected:
    /// Creates the weighting as void*, depending on the type.
    virtual void *createWeightV(std::size_t sizex, std::size_t sizey, std::size_t sizez, FreqWeightDesc::data_type t) const {
        if (t == FreqWeightDesc::type_float) {
            return new tom::FreqWeight_AllPass<float>(sizex, sizey, sizez);
        }
        assert(t == FreqWeightDesc::type_double);
        return new tom::FreqWeight_AllPass<double>(sizex, sizey, sizez);
    }

};



/****************************************************************************//**
 * \brief A frequency weighting for a single axis wedge.
 *
 * The wedge has two parameters. One is the (half) opening angle and the other
 * an \c cutoff_radius for a low-pass filter.
 *******************************************************************************/
template<typename T>
class FreqWeight_SingleAxisWedge: public FreqWeight<T> {

public:

    FreqWeight_SingleAxisWedge(std::size_t sizex, std::size_t sizey, std::size_t sizez, double angle, double cutoff_radius);

    /// Destructor
    virtual ~FreqWeight_SingleAxisWedge() {
        delete w_redu_;
        delete w_full_;
    }

    /// Rotate the weighting into a certain orientation (angles in radians).
    virtual void rotate(double phi, double psi, double the) {
        phi_ = phi;
        psi_ = psi;
        the_ = the;
    }

    virtual const tom::Volume<T> *get_weight(bool reduced_complex) const;

    /// Is the weighting an all-pass for a volume of this size?
    /// Returns always false, because the opening angle is > 0. and therefore always some frequencies
    /// are filtered out.
    virtual bool all_pass() const {
        return  false;
    }

    virtual bool weight(tom::Volume<typename FreqWeight<T>::Tcomplex> &vsrc) const;

    /// Returns the opening angle in degree.
    double getAngleDeg() const     { return tom::math::rad2deg(angle_); }
    /// Returns the opening angle in radians.
    double getAngleRad() const     { return angle_; }
    /// Returns the cut off radius (0. means no cut off).
    double getCutoffRadius() const { return cutoff_radius_; }


private:
    FreqWeight_SingleAxisWedge(const FreqWeight_SingleAxisWedge<T> &c): FreqWeight<T>(1,1,1) { assert(!"HIDDEN"); } //hide copy constructor.
    FreqWeight_SingleAxisWedge &operator=(const FreqWeight_SingleAxisWedge &c)               { assert(!"HIDDEN"); return *this; } // hide assignment constructor.

    // Caches (are mutable)
    mutable tom::Volume<T> *w_redu_;
    mutable tom::Volume<T> *w_full_;

    // Current oritentation of the cache.
    mutable double cphi_redu_;
    mutable double cpsi_redu_;
    mutable double cthe_redu_;
    mutable double cphi_full_;
    mutable double cpsi_full_;
    mutable double cthe_full_;

    // Means that at least the upper half of the full volume is non valid.
    // It can only become invalid if the reduced weighting shares the volume
    // and overwrites its half.
    mutable bool w_full_valid_;

    // Opening angle. In radians.
    double angle_;

    // Current orientation.
    double phi_;
    double psi_;
    double the_;

    double cutoff_radius_;

    void init_wedge_volume(bool reduced) const;
    void init_wedge_volume(tom::Volume<T> &vol) const;

    #ifdef ___TEST_PERFORMANCE_1_
    // friend declaration for a test (main) routine decleared in FreqWeight.cpp
    //friend int ::main(int argc, char **argv);
    #endif
};



/****************************************************************************//**
 * \brief Descriptor class for SingleAxisWedge.
 *******************************************************************************/
class FreqWeightDesc_SingleAxisWedge: public FreqWeightDesc {
public:
    /***********************************************************************//**
    * \brief Constructor
    *
    * \param[in] angle The opening angle is in radians.
    * \param[in] cutoff_radius Frequencies larger than this cutoff_radius are set
    *   to 0 and therefore this is the cutoff frequency of a low-pass filter.
    *   If \c nyquist_x is 0, the cutoff_radius is meant as an absolute value,
    *   regardless of what the actuall volume size may be. If \c nyquist_x
    *   is a positive number, the cutoff_radius of the created weigthing is scaled
    *   and taken relatively.
    *   Setting \c cutoff_radius to 0, means no cutoff.
    * \param[in] nyquist_x The intended Nyquist frequency of the weighting in X-direction.
    *   If 0, the \c cutoff_radius is taken as an absolute value of frequencies.
    *   This would be for example useful, if you know that a volume has information
    *   up to a certain frequency. Rescaling/binning the volume, does not actually
    *   change the cutoff_radius, i.e. the region in Fourier space where there
    *   is any information.
    *   By setting \c nyquist_x to a volume size, \c cutoff_radius is taken relatively
    *   to the actual volume size of the weighting (see \c getResultingCutoffRadius).
    ***************************************************************************/
    FreqWeightDesc_SingleAxisWedge(double angle, double cutoff_radius, std::size_t nyquist_x)
        :   FreqWeightDesc(),
            cutoff_radius_(cutoff_radius),
            nyquist_x_(nyquist_x),
            angle_(angle) {
        if (angle_ <= 0. || cutoff_radius_<0.) {
            throw std::invalid_argument("The SingleAxisWedge must have an opening angle > 0 and a cutoff_radius >= 0 (0 means no cut off).");
        }
    }

    /// Destructor (no operation).
    virtual ~FreqWeightDesc_SingleAxisWedge() { };

    /// Returns the opening angle in degrees.
    double getAngleDeg() const     { return tom::math::rad2deg(angle_); }
    /// Returns the opening angle in radians.
    double getAngleRad() const     { return angle_; }
    /// Returns the cut off radius (0 means no cut off).
    double getCutoffRadius() const { return cutoff_radius_; }
    /// Returns the cut off radius which will be applied to a weight of a certain size (0 means no cut off).
    double getResultingCutoffRadius(std::size_t sizex, std::size_t sizey, std::size_t sizez) const {
        if (nyquist_x_) {
            return cutoff_radius_ / static_cast<double>(nyquist_x_) * static_cast<double>(sizex/2);
        }
        return cutoff_radius_;
    }
    /// Returns the Nyquist frequency for dimension X.
    std::size_t getNyquistX() const { return nyquist_x_; }

    /// Return a description.
    virtual std::string toString() const {
        std::stringstream ss;
        ss << "singleaxiswedge " << getAngleDeg() << " " << cutoff_radius_;
        if (nyquist_x_) {
            ss << " " << nyquist_x_;
        }
        return ss.str();
    }

    /// Returns a unique Id depending on the type, starting at 0.
    virtual Identifier_FreqWeight getIdentifier() const { return identifier_SingleAxisWedge; };

protected:
    /***********************************************************************//**
     * Init a weighting from the description.
     * Because template functions can not be virtual, this function returns a
     * pointer to void.
     **************************************************************************/
    virtual void *createWeightV(std::size_t sizex, std::size_t sizey, std::size_t sizez, FreqWeightDesc::data_type t) const {
        if (t == FreqWeightDesc::type_float) {
            return new tom::FreqWeight_SingleAxisWedge<float>(sizex, sizey, sizez, getAngleRad(), getResultingCutoffRadius(sizex, sizey, sizez));
        }
        assert(t == FreqWeightDesc::type_double);
        return new tom::FreqWeight_SingleAxisWedge<double>(sizex, sizey, sizez, getAngleRad(), getResultingCutoffRadius(sizex, sizey, sizez));
    }

private:
    double cutoff_radius_;
    std::size_t nyquist_x_;

    // half opening angle in radians
    double angle_;
};


/****************************************************************************//**
 * \brief A frequency weighting for a asymmetric, single axis wedge.
 *
 * The wedge has three parameters. The first two define the opening angles. and the other
 * an \c cutoff_radius for a low-pass filter.
 *******************************************************************************/
template<typename T>
class FreqWeight_AsymmetricSingleAxisWedge: public FreqWeight<T> {

public:

	FreqWeight_AsymmetricSingleAxisWedge(std::size_t sizex, std::size_t sizey, std::size_t sizez, double angle1, double angle2, double cutoff_radius);

    /// Destructor
    virtual ~FreqWeight_AsymmetricSingleAxisWedge() {
        delete w_redu_;
        delete w_full_;
    }

    /// Rotate the weighting into a certain orientation (angles in radians).
    virtual void rotate(double phi, double psi, double the) {
        phi_ = phi;
        psi_ = psi;
        the_ = the;
    }

    virtual const tom::Volume<T> *get_weight(bool reduced_complex) const;

    /// Is the weighting an all-pass for a volume of this size?
    /// Returns always false, because the opening angle is > 0. and therefore always some frequencies
    /// are filtered out.
    virtual bool all_pass() const {
        return false;
    }

    virtual bool weight(tom::Volume<typename FreqWeight<T>::Tcomplex> &vsrc) const;

    /// Returns the opening angle in degree.
    double getFirstAngleDeg() const     { return tom::math::rad2deg(angle1_); }
    /// Returns the opening angle in radians.
    double getFirstAngleRad() const     { return angle1_; }

    /// Returns the opening angle in degree.
    double getSecondAngleDeg() const     { return tom::math::rad2deg(angle2_); }
	/// Returns the opening angle in radians.
	double getSecondAngleRad() const     { return angle2_; }

    /// Returns the cut off radius (0. means no cut off).
    double getCutoffRadius() const { return cutoff_radius_; }


private:
	//FreqWeight_SingleAxisWedge(const FreqWeight_SingleAxisWedge<T> &c): FreqWeight<T>(1,1,1) { assert(!"HIDDEN"); } // hide copy constructor.
	FreqWeight_AsymmetricSingleAxisWedge(const FreqWeight_AsymmetricSingleAxisWedge<T> &c): FreqWeight<T>(1,1,1) { assert(!"HIDDEN"); } // hide copy constructor.
	FreqWeight_AsymmetricSingleAxisWedge &operator=(const FreqWeight_AsymmetricSingleAxisWedge<T> &c)               { assert(!"HIDDEN"); return *this; } // hide assignment constructor.

    // Caches (are mutable)
    mutable tom::Volume<T> *w_redu_;
    mutable tom::Volume<T> *w_full_;

    // Current oritentation of the cache.
    mutable double cphi_redu_;
    mutable double cpsi_redu_;
    mutable double cthe_redu_;
    mutable double cphi_full_;
    mutable double cpsi_full_;
    mutable double cthe_full_;

    // Means that at least the upper half of the full volume is non valid.
    // It can only become invalid if the reduced weighting shares the volume
    // and overwrites its half.
    mutable bool w_full_valid_;

    // Opening angles. In radians.
    double angle1_;
    double angle2_;

    // Current orientation.
    double phi_;
    double psi_;
    double the_;

    double cutoff_radius_;

    void init_wedge_volume(bool reduced) const;
    void init_wedge_volume(tom::Volume<T> &vol) const;

};


/****************************************************************************//**
 * \brief A frequency weighting for a asymmetric, single axis wedge.
 *
 * The wedge has three parameters. The first two define the opening angles. and the other
 * an \c cutoff_radius for a low-pass filter.
 *******************************************************************************/
template<typename T>
class FreqWeight_SmoothedAsymmetricSingleAxisWedge: public FreqWeight<T> {

public:

	FreqWeight_SmoothedAsymmetricSingleAxisWedge(std::size_t sizex, std::size_t sizey, std::size_t sizez, double angle1, double angle2, double cutoff_radius,double smooth);

    /// Destructor
    virtual ~FreqWeight_SmoothedAsymmetricSingleAxisWedge() {
        delete w_redu_;
        delete w_full_;
    }

    /// Rotate the weighting into a certain orientation (angles in radians).
    virtual void rotate(double phi, double psi, double the) {
        phi_ = phi;
        psi_ = psi;
        the_ = the;
    }

    virtual const tom::Volume<T> *get_weight(bool reduced_complex) const;

    /// Is the weighting an all-pass for a volume of this size?
    /// Returns always false, because the opening angle is > 0. and therefore always some frequencies
    /// are filtered out.
    virtual bool all_pass() const {
        return false;
    }

    virtual bool weight(tom::Volume<typename FreqWeight<T>::Tcomplex> &vsrc) const;

    /// Returns the opening angle in degree.
    double getFirstAngleDeg() const     { return tom::math::rad2deg(angle1_); }
    /// Returns the opening angle in radians.
    double getFirstAngleRad() const     { return angle1_; }

    /// Returns the opening angle in degree.
    double getSecondAngleDeg() const     { return tom::math::rad2deg(angle2_); }
	/// Returns the opening angle in radians.
	double getSecondAngleRad() const     { return angle2_; }

    /// Returns the cut off radius (0. means no cut off).
    double getCutoffRadius() const { return cutoff_radius_; }

    /// Returns smooth value
    double getSmoothRange() const { return smooth_;}

private:
	FreqWeight_SmoothedAsymmetricSingleAxisWedge(const FreqWeight_SmoothedAsymmetricSingleAxisWedge<T> &c): FreqWeight<T>(1,1,1) { assert(!"HIDDEN"); } // hide copy constructor.
	FreqWeight_SmoothedAsymmetricSingleAxisWedge &operator=(const FreqWeight_SmoothedAsymmetricSingleAxisWedge<T> &c)            { assert(!"HIDDEN"); return *this; } // hide assignment constructor.

    // Caches (are mutable)
    mutable tom::Volume<T> *w_redu_;
    mutable tom::Volume<T> *w_full_;

    // Current oritentation of the cache.
    mutable double cphi_redu_;
    mutable double cpsi_redu_;
    mutable double cthe_redu_;
    mutable double cphi_full_;
    mutable double cpsi_full_;
    mutable double cthe_full_;

    // Means that at least the upper half of the full volume is non valid.
    // It can only become invalid if the reduced weighting shares the volume
    // and overwrites its half.
    mutable bool w_full_valid_;

    // Opening angles. In radians.
    double angle1_;
    double angle2_;

    // Current orientation.
    double phi_;
    double psi_;
    double the_;

    double cutoff_radius_;

	double smooth_;

    void init_wedge_volume(bool reduced) const;
    void init_wedge_volume(tom::Volume<T> &vol) const;

};

/****************************************************************************//**
 * \brief Weighting initialised from a tom::Volume.
 *
 * Applying the weighting to a different sized volume, an exception is thrown.
 *******************************************************************************/
template<typename T>
class FreqWeight_Volume: public FreqWeight<T> {

public:

    FreqWeight_Volume(tom::Volume<T> &w, bool copy);
    FreqWeight_Volume(const tom::Volume<T> &w);

    /// Destructor.
    virtual ~FreqWeight_Volume() {
        delete w_full_;
        delete w_redu_;

        // Even if w_ is not a copy (according to the constructor), this
        // tom::Volume-object is not shared, only the underlying data.
        delete w_;
    }

    /// Gives access to the not rotated volume. Do not modify or delete the pointer!
    const tom::Volume<T> *get_volume() const {
        assert(w_);
        return w_;
    }

    /***********************************************************************//**
     * \brief Rotates the weighting.
     *
     * The volume is not actually rotated. This is postponed until
     * it is needed the first time.
     * Rotating a volume in real space around the centre, corresponds to a Rotation in Fourier space
     * of the weighting (where the phase is ignored).
     **************************************************************************/
    virtual void rotate(double phi, double psi, double the) {
        phi_ = phi;
        psi_ = psi;
        the_ = the;
    }


    virtual bool weight(tom::Volume<typename FreqWeight<T>::Tcomplex> &vsrc) const;

    /// Return whether the weighting is an all pass.
    virtual bool all_pass() const {
        return all_pass_;
    }

    virtual const tom::Volume<T> *get_weight(bool reduced_complex) const;

protected:
    FreqWeight_Volume(); // Needed for specializations such as FreqWeight_EmFile.

    tom::Volume<T> *w_;

    void resetVolume(std::auto_ptr<tom::Volume<T> > &w);

    bool all_pass_;

private:
    FreqWeight_Volume(const FreqWeight_Volume<T> &c): FreqWeight<T>(1,1,1) { assert(!"HIDDEN");               } // hide copy constructor.
    FreqWeight_Volume &operator=(const FreqWeight_Volume &c)               { assert(!"HIDDEN"); return *this; } // hide assignment constructor.

    void init_wedge_volume() const;

    // Cache for the rotated volumes.
    mutable tom::Volume<T> *w_full_;
    mutable tom::Volume<T> *w_redu_;

    // For the current orientation as saved in wedge_full_ and wedge_reduced_
    mutable double cphi_;
    mutable double cpsi_;
    mutable double cthe_;

    // The orientation last set by rotate.
    double phi_;
    double psi_;
    double the_;
};




/****************************************************************************//**
 * \brief A FreqWeight_Volume which is read from an EM-file.
 *******************************************************************************/
template <typename T>
class FreqWeight_EmFile: public FreqWeight_Volume<T> {

public:
    FreqWeight_EmFile(const std::string &filename, size_t sizex, size_t sizey, size_t sizez, bool stretch, bool resize);

    /// Destructor (no operation).
    virtual ~FreqWeight_EmFile() { }

private:
    FreqWeight_EmFile(const FreqWeight_EmFile<T> &c) { assert(0); } // hide copy constructor.
    FreqWeight_EmFile &operator=(const FreqWeight_EmFile &c) { assert(0); return *this; } // hide assignment constructor.

    std::string filename_;
    bool stretch_;
    bool resize_;
};



/****************************************************************************//**
 * \brief A descriptor for the FreqWeight_EmFile.
 *******************************************************************************/
class FreqWeightDesc_EmFile: public FreqWeightDesc {
public:
    /// Constructor
    FreqWeightDesc_EmFile(const std::string &filename, bool stretch, bool resize)
        :   FreqWeightDesc(),
            stretch_(stretch),
            filename_(filename),
            resize_(resize) {
    };

    /// Destructor (no operation).
    virtual ~FreqWeightDesc_EmFile() { };

    /// Get the config entry text
    virtual std::string toString() const {
        std::ostringstream ss;
        ss << "emfile " << (stretch_ ? "" : "no") << "stretch " << (resize_?"":"no") << "resize   \"" << filename_ << "\"";
        return ss.str();
    }

    /// Get the name of the EM-File
    std::string getFilename() const {
        return filename_;
    }

    /// Does the weighting stretch/resize the volume or cut out the relevant part?
    bool getStretch() const {
        return stretch_;
    }

    /// The weight throws an exception if the weighting is constructed for a
    /// different size than the volume in the emfile.
    bool getResize() const {
        return resize_;
    }

    void setResize(bool resize) {
        resize_ = resize;
    }

    /// Returns a unique Id depending on the type, starting at 0.
    virtual Identifier_FreqWeight getIdentifier() const { return identifier_EmFile; };

    typedef std::pair<float,float> band;

protected:
    virtual void *createWeightV(std::size_t sizex, std::size_t sizey, std::size_t sizez, FreqWeightDesc::data_type t) const {
        if (t == FreqWeightDesc::type_float) {
            return new tom::FreqWeight_EmFile<float>(getFilename(), sizex, sizey, sizez, getStretch(), getResize());
        }
        assert(t == FreqWeightDesc::type_double);
        return new tom::FreqWeight_EmFile<double>(getFilename(), sizex, sizey, sizez, getStretch(), getResize());
    }
private:

	std::vector<band> band_vector_;

    bool stretch_;
    std::string filename_;
    bool resize_;


};



/****************************************************************************//**
 * \brief Bandpass filter class.
 *
 * Bandpass filter class.
 *******************************************************************************/
template<typename T>
class FreqWeight_Bandpass: public FreqWeight<T> {

public:

	FreqWeight_Bandpass(const float lowestFrequency, const float highestFrequency,const std::size_t x,const std::size_t y,const std::size_t z, const float smooth);

    /// Destructor.
    virtual ~FreqWeight_Bandpass() {
        // Even if w_ is not a copy (according to the constructor), this
        // tom::Volume-object is not shared, only the underlying data.

        delete b_;
    }

    /// Gives access to the not rotated volume. Do not modify or delete the pointer!
    const tom::Volume<T> *get_volume() const {
        assert(b_);
        return b_;
    }

    virtual void rotate(double phi, double psi, double the) {};


    virtual bool weight(tom::Volume<typename FreqWeight<T>::Tcomplex> &vsrc) const;

    /// Return whether the weighting is an all pass filter.
    virtual bool all_pass() const {
        return all_pass_;
    }

    virtual const tom::Volume<T> *get_weight(bool b) const;

    void setSmooth(float smooth){this->smooth=smooth;};
    float getSmooth(){return smooth;};

    void setLowestFrequency(float f){this->lowestFrequency=f;};
    float getLowestFrequency(){return this->lowestFrequency;};

    
    void setHighestFrequency(float f){this->highestFrequency=f;};
    float getHighestFrequency(){return this->highestFrequency;};

    

protected:

	mutable tom::Volume<T> *b_;

    bool all_pass_;

private:
	FreqWeight_Bandpass &operator=(const FreqWeight_Bandpass &c)               { assert(!"HIDDEN"); return *this; } // hide assignment constructor.

	std::size_t size_x;
	std::size_t size_y;
	std::size_t size_z;

    void init_bandpass_volume() const;
    float lowestFrequency;
    float highestFrequency;
    float smooth;

};

} // namespace tom


// Inline functions.


/// Operator<< to output the FreqWeight.
/// Calls the virtual method toString and acts therefore polymorph.
template<typename T>
inline std::ostream &operator<<(std::ostream &os, const tom::FreqWeight<T> &w) {
    os << w.toString();
    return os;
}





#endif









