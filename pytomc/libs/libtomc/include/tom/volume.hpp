/****************************************************************************//**
 * \file volume.hpp
 * \brief The header file for the class tom::Volume.
 * \author  Thomas Haller
 * \version 0.2
 * \date    22.11.2007
 *******************************************************************************/
#ifndef ___INCLUDE__TOM__VOLUME_HPP__
#define ___INCLUDE__TOM__VOLUME_HPP__


#include <cstddef>
#include <stdexcept>
#include <complex>
#include <vector>
#include <assert.h>




#ifndef NDEBUG
#   include <sys/types.h>
#   include <unistd.h>
#   include <iostream>
#   define DBGMSGE(cmd) { std::cerr << getpid() << ":" << __FILE__ << ":" <<__LINE__ << ": " << cmd << std::endl; }
#   define DBGMSG(cmd)  { std::cout << getpid() << ":" << __FILE__ << ":" <<__LINE__ << ": " << cmd << std::endl; }
#else
#   define DBGMSGE(cmd)
#   define DBGMSG(cmd)
#endif


namespace tom {
template<typename T> class Volume;
}

#include <tom/io/io.hpp>



namespace tom {


void setFcnMem( void *(*fcn_malloc)(size_t n), void (*fcn_free)(void *ptr));
void *(*(getFcnMalloc()))(size_t n);
void (*(getFcnFree()))(void *ptr);
template <typename T> void *fcn_malloc_new(size_t size);
template <typename T> void fcn_free_delete(void *ptr);



/** \brief Forward declaration of st_volume_memory from tom_volume.cpp. */
struct st_volume_memory;


#ifndef NO_DEFINE_DEPRECATED
#   define DEFINE_DEPRECATED            __attribute__ ((deprecated))
#else
#   define DEFINE_DEPRECATED
#endif



/****************************************************************************//**
 * \brief Class to save a volume.
 *
 * A generic class which holds a pointer to a volume. The volume must be
 * created/initialised when the object is instanciated and can not be changed.\n
 * The class counts its references to the same memory and frees it automatically.\n
 * The function to allocate memory can be set using tom::setFcnMem().\n
 * Volume can share the same memory using the right constructor. However doing
 * so the original size of the volume can not be exceeded (out_of_range).
 *******************************************************************************/
template<typename T>
class Volume {

public:

    typedef T element_type;

    Volume(std::size_t sizex, std::size_t sizey, std::size_t sizez, void *(*fcn_malloc)(size_t size), void (*fcn_free)(void *ptr));
    Volume(T *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez, bool free, void (*fcn_free)(void *ptr));
    Volume(const Volume<T> &v);
    //Volume(Volume<T> &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez);
    template<typename T2> Volume(const Volume<T2> &v);
    template<typename T2> Volume(Volume<T2> &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez);
    template<typename T2> Volume(Volume<T2> &v, bool get_real_part);

    ~Volume();

    std::size_t getSizeX() const;
    std::size_t getSizeY() const;
    std::size_t getSizeZ() const;
    std::size_t getStrideX() const;
    std::size_t getStrideY() const;
    std::size_t getStrideZ() const;

    bool to_free_memory() const;
    int number_memory_used() const;
    void (*get_free_function())(void *ptr);
    std::size_t get_memsize() const;

    bool isContiguous() const;

    const T &get() const;
    const T &get(std::size_t x, std::size_t y, std::size_t z) const;
    T &get();
    T &get(std::size_t x, std::size_t y, std::size_t z);
    const T &get_unsafe(std::size_t x, std::size_t y, std::size_t z) const;
    T &get_unsafe(std::size_t x, std::size_t y, std::size_t z);

    /// Don't use these functions but instead the one from tom::io (io.hpp)
    void write_to_em(const std::string &filename, const tom_io_em_header *header                             ) const DEFINE_DEPRECATED;
    void write_to_em(const std::string &filename, const tom_io_em_header *header, const uint32_t* first_voxel) const DEFINE_DEPRECATED;

    void setValues(T val);
    template <typename T2> void setValues(const Volume<T2> &v);

    template<typename TPRECISION> void shift(TPRECISION factor_shift);
    template<typename TPRECISION> void scale(TPRECISION factor_scale);
    template<typename TPRECISION> void shift_scale(TPRECISION factor_shift, TPRECISION factor_scale);
    template<typename TPRECISION> void scale_shift(TPRECISION factor_scale, TPRECISION factor_shift);

    bool is_equal_size(std::size_t s) const;
    bool is_equal_size(std::size_t x, std::size_t y, std::size_t z) const;
    template<typename T2> bool is_equal_size(const tom::Volume<T2> &v2) const;


    // Avoid the following functions. They are not meaningful for every template type
    // and should be removed from the interface of tom::Volume.
    // To compute the variance and the mean, use the functions from volume_fcn instead.
    // For example for std::complex<T> these functions are not meaningful and should therefore
    // not be part of the interface.
    T min() const                                                                                                                       DEFINE_DEPRECATED;
    T max() const                                                                                                                       DEFINE_DEPRECATED;
    void minmax(T &min, T &max) const                                                                                                   DEFINE_DEPRECATED;
    double mean() const                                                                                                                 DEFINE_DEPRECATED;
    double variance(bool use_sample_standard_deviation) const                                                                           DEFINE_DEPRECATED;
    double variance_mean_free(bool use_sample_standard_deviation) const                                                                 DEFINE_DEPRECATED;
    void stat(double &m, double &variance, bool use_sample_standard_deviation) const                                                    DEFINE_DEPRECATED;
    template<typename T2> void stat(double &m, double &variance, bool use_sample_standard_deviation, const tom::Volume<T2> &mask) const DEFINE_DEPRECATED;
    void stat(double &m, double &variance, T &min, T &max, bool use_sample_standard_deviation) const                                    DEFINE_DEPRECATED;

    bool operator==(const T &val) const;
    bool operator==(const tom::Volume<T> &v) const;
    void printInfo() const;
    void printInfo(const std::string &name) const;

    std::size_t getByteOffset(std::size_t x, std::size_t y, std::size_t z) const;

    bool is_equal_memory(const tom::Volume<T> &vol) const;
    template<typename T2> bool is_shared_memory(const tom::Volume<T2> &vol) const;

    std::size_t numel() const;



private:

    T *data;
    std::size_t sizex;
    std::size_t sizey;
    std::size_t sizez;
    std::size_t stridex;
    std::size_t stridey;
    std::size_t stridez;
    ::tom::st_volume_memory *volume_memory;

    Volume(); // Hide empty constructor.
    tom::Volume<T> &operator=(const tom::Volume<T> &v) { throw std::runtime_error("Assignment operator not allowed"); } // Hide assignment operator.

    void stat(double *m, double *variance, T *min, T *max, bool use_sample_standard_deviation) const /*DEFINE_DEPRECATED*/;

    template<typename T2> void share_memory(Volume<T2> &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez);

    void initialize(std::size_t sizex, std::size_t sizey, std::size_t sizez, void *(*fcn_malloc)(size_t size), void (*fcn_free)(void *ptr));


    T &get_(std::size_t x, std::size_t y, std::size_t z);
    const T &get_(std::size_t x, std::size_t y, std::size_t z) const;


    friend class tom::Volume<float           >;
    friend class tom::Volume<double          >;
    friend class tom::Volume<std::complex<float>  >;
    friend class tom::Volume<std::complex<double> >;

}; // class tom::Volume<T>


template<typename T> bool less_ptr_Volume(const tom::Volume<T> * v1, const tom::Volume<T> *v2);



/****************************************************************************//**
 * A class providing as operator a comparison of two pointers to tom::Volume
 * by calling tom::less_ptr_Volume.
 *******************************************************************************/
template<class T>
struct volume_less: public std::binary_function<const tom::Volume<T> *, const tom::Volume<T> *, bool> {
    bool operator()(const tom::Volume<T> *__x, const tom::Volume<T> *__y) const {
        return tom::less_ptr_Volume<T>(__x, __y);
    }
};

} // namespace tom









// Inline functions...



/** \brief Returns the size of the volume along X. */
template<typename T>
inline std::size_t tom::Volume<T>::getSizeX() const {
    return this->sizex;
}

/** \brief Returns the size of the volume along Y. */
template<typename T>
inline std::size_t tom::Volume<T>::getSizeY() const {
    return this->sizey;
}

/****************************************************************************//**
 * \brief Returns the size of the volume along Z.
 *
 * It is the fastest dimension. First all the voxels along X apear in memory
 * with all the same Y and Z coordinate.
 *******************************************************************************/
template<typename T>
inline std::size_t tom::Volume<T>::getSizeZ() const {
    return this->sizez;
}



/****************************************************************************//**
 * \brief Returns the size of elements in the volume.
 *******************************************************************************/
template<typename T>
inline std::size_t tom::Volume<T>::numel() const {
    return this->getSizeX()*this->getSizeY()*this->getSizeZ();
}



/****************************************************************************//**
 * \brief Copies a volume (templated "copy constructor")
 *
 * \param[in] v A reference to the source volume.
 *
 * Allocates new contiguous memory holding all elements of \c v .
 * The elements are copied using setValues. Beware, that is uses the functions
 * saved in tom::setFcnMem() to allocate the memory (where NULL defaults to
 * the new[]/delete[] operator).
 *******************************************************************************/
template<typename T> template<typename T2>
tom::Volume<T>::Volume(const Volume<T2> &v) {
    this->initialize(v.getSizeX(), v.getSizeY(), v.getSizeZ(), NULL, NULL);
    this->setValues<T2>(v);
}



/****************************************************************************//**
 * \brief Returns the distance between two neighbouring elements along X.
 *
 * The distance is measured in bytes, thus it must be at least sizeof(T).
 *******************************************************************************/
template<typename T>
inline std::size_t tom::Volume<T>::getStrideX() const {
    return this->stridex;
}

/****************************************************************************//**
 * \brief Returns the distance between two neighbouring elements along Y.
 *
 * The distance is measured in bytes. It must be at least getSizeX()*getStrideX().
 *******************************************************************************/
template<typename T>
inline std::size_t tom::Volume<T>::getStrideY() const {
    return this->stridey;
}

/****************************************************************************//**
 * \brief Returns the distance between two neighbouring elements along Z.
 *
 * The distance is measured in bytes. It must be at least getSizeY()*getStrideY().
 *******************************************************************************/
template<typename T>
inline std::size_t tom::Volume<T>::getStrideZ() const {
    return this->stridez;
}



/** \brief Returns true if there is not gap between two voxels along each direction. */
template<typename T>
inline bool tom::Volume<T>::isContiguous() const {
    return  this->getStrideX()==sizeof(T) &&
            this->getStrideY()==this->getSizeX()*this->getStrideX() &&
            this->getStrideZ()==this->getSizeY()*this->getStrideY();
}





/****************************************************************************//**
 * \brief Constructor of the class. Allocates memory.
 *
 * \param[in] sizex The size of the volume. Must be positive.
 * \param[in] sizey The size of the volume. Must be positive.
 * \param[in] sizez The size of the volume. Must be positive.
 * \param[in] fcn_malloc The allocation function used to allocate the memory.
 *   Setting this to \c NULL defaults to tom::getFcnMalloc(). Setting the
 *   default function to \c NULL results in using the new[] operator.
 * \param[in] fcn_free The function to de-allocate the memory afterwards.
 *   \c NULL defaults to tom::getFcnFree().
 *
 * It is guaranteed that the allocated memory is contiguous.\n
 * fcn_malloc and fcn_free must correspond to each other. They both must be
 * set to NULL or both must be set to a valid function.\n
 * The memory is owned by the volume, that means it's destructor will
 * call fcn_free. However if there are more then one Volumes which point
 * to this volume only the last one will free it.
 *******************************************************************************/
template <typename T>
inline tom::Volume<T>::Volume(std::size_t sizex, std::size_t sizey, std::size_t sizez, void *(*fcn_malloc)(size_t size), void (*fcn_free)(void *ptr))  {
    this->initialize(sizex, sizey, sizez, fcn_malloc, fcn_free);
}





/****************************************************************************//**
 * \brief Copies a volume (copy constructor)
 *
 * \param[in] v A reference to the source volume.
 *
 * Allocates new contiguous memory holding all elements of \c v .
 * The elements are copied using setValues. Beware, that is uses the functions
 * saved in tom::setFcnMem() to allocate the memory (where NULL defaults to
 * the new[]/delete[] operator).
 *******************************************************************************/
template<typename T>
inline tom::Volume<T>::Volume(const Volume<T> &v) {

    this->initialize(v.getSizeX(), v.getSizeY(), v.getSizeZ(), NULL, NULL);
    this->setValues<T>(v);
}


/** Const access to the first voxel. */
template <typename T>
inline const T &tom::Volume<T>::get() const {
    return this->data[0];
}

/** Access to the first voxel. */
template <typename T>
inline T &tom::Volume<T>::get() {
    return this->data[0];
}

/****************************************************************************//**
 * \brief Const access to the data.
 *
 * This function return a reference to a voxel of the allocated memory.
 * An exception is thrown when accessing out of range.\n
 * Notice that the coordinate is specified as (X,Y,Z) instead of (Z,Y,X).
 * This way it conforms to all other prototypes where X comes first too.\n
 *******************************************************************************/
template <typename T>
inline const T &tom::Volume<T>::get(std::size_t x, std::size_t y, std::size_t z) const {
    if (x >= this->getSizeX() || y >= this->getSizeY() || z >= this->getSizeZ()) {
        throw std::out_of_range("get() const tries to access out of range.");
    }
    return (reinterpret_cast<const T *>(reinterpret_cast<char *>(this->data) + this->getByteOffset(x,y,z)))[0];
}

/****************************************************************************//**
 * \brief Access to the data.
 *
 * This function return a reference to a voxel of the allocated memory.
 * An exception is thrown when accessing out of range.\n
 * Notice that the coordinate is specified as (X,Y,Z) instead of (Z,Y,X).
 * This way it conforms to all other prototypes where X comes first too.\n
 *******************************************************************************/
template <typename T>
inline T &tom::Volume<T>::get(std::size_t x, std::size_t y, std::size_t z) {
    if (x >= this->getSizeX() || y >= this->getSizeY() || z >= this->getSizeZ()) {
        throw std::out_of_range("get() tries to access out of range.");
    }
    return (reinterpret_cast<T *>(reinterpret_cast<char *>(this->data) + this->getByteOffset(x,y,z)))[0];
}




/****************************************************************************//**
 * \brief Access to the data without out of range check.
 *******************************************************************************/
template <typename T>
inline const T &tom::Volume<T>::get_unsafe(std::size_t x, std::size_t y, std::size_t z) const {
    return (reinterpret_cast<T *>(reinterpret_cast<char *>(this->data) + this->getByteOffset(x,y,z)))[0];
}

/****************************************************************************//**
 * \brief Access to the data without out of range check.
 *******************************************************************************/
template <typename T>
inline T &tom::Volume<T>::get_unsafe(std::size_t x, std::size_t y, std::size_t z) {
    return (reinterpret_cast<T *>(reinterpret_cast<char *>(this->data) + this->getByteOffset(x,y,z)))[0];
}


/****************************************************************************//**
 * \brief Returns byte offset of a voxel from the base pointer.
 *
 * Returns the distance of a voxel-coordinate (zero based) from the beginning
 * of the data (in bytes). Does not perform out of range check.
 *******************************************************************************/
template <typename T>
inline std::size_t tom::Volume<T>::getByteOffset(std::size_t x, std::size_t y, std::size_t z) const {
    return x*this->getStrideX() + y*this->getStrideY() + z*this->getStrideZ();
}



/****************************************************************************//**
 * \brief Constructor which cuts out a subvolume of an already existing one.
 *
 * \param[in] v A reference to the volume which memory will be taken.
 * \param[in] data A pointer to the first element of the new volume. Setting
 *   this parameter to NULL, means taking the first element of v. (= &v.get()).
 * \param[in] sizex Size of the new volume.
 * \param[in] sizey Size of the new volume.
 * \param[in] sizez Size of the new volume.
 * \param[in] stridex Stride along x (in bytes).
 * \param[in] stridey Stride along y (in bytes).
 * \param[in] stridez Stride along z (in bytes).
 *
 * This constructor allows to share the same memory between different volumes.
 * For example if you have a volume of complex numbers fftwf_complex and want
 * a subregion of the imaginary part, you can create a new volume, using the
 * same memory. \n
 * Again, ommitting the stride parameters by setting them to 0, defaults to
 * contiguous memory (i.e. stridex=sizeof(T), stridex=sizex*stridex, ...).
 * Be carefull, in case of different type: that may not be what you want.
 *******************************************************************************/
template<typename T> template<typename T2>
inline tom::Volume<T>::Volume(Volume<T2> &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez) {
    this->share_memory(v, data, sizex, sizey, sizez, stridex, stridey, stridez);
}

/****************************************************************************//**
 * \brief Constructor which cuts out a subvolume of an already existing one.
 *
 * \param[in] v A reference to the volume which memory will be taken.
 * \param[in] data A pointer to the first element of the new volume. Setting
 *   this parameter to NULL, means taking the first element of v. (= &v.get()).
 * \param[in] sizex Size of the new volume.
 * \param[in] sizey Size of the new volume.
 * \param[in] sizez Size of the new volume.
 * \param[in] stridex Stride along x (in bytes).
 * \param[in] stridey Stride along y (in bytes).
 * \param[in] stridez Stride along z (in bytes).
 *
 * This constructor allows to share the same memory between different volumes.
 * For example if you have a volume of complex numbers fftwf_complex and want
 * a subregion of the imaginary part, you can create a new volume, using the
 * same memory. \n
 * Again, ommitting the stride parameters by setting them to 0, defaults to
 * contiguous memory (i.e. stridex=sizeof(T), stridex=sizex*stridex, ...).
 * Be carefull, in case of different type: that may not be what you want.
template<typename T>
inline tom::Volume<T>::Volume(Volume<T> &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez) {
    this->share_memory(v, data, sizex, sizey, sizez, stridex, stridey, stridez);
}
*******************************************************************************/

/****************************************************************************//**
 * \brief Constructor which takes the real (or imaginary) part of a complex volume
 * \brief Does not work!
 *******************************************************************************/
template<typename T> template<typename T2>
inline tom::Volume<T>::Volume(Volume<T2> &v, bool get_real_part) {
    // T2 is not strictly set to std::complex<T> to allow other type of volumes
    // as for example int, where no std::complex<int> exists.
    // Thus T2 is a template-type.
    void *data = get_real_part ? &v.get().real() : &v.get().imag();
    this->share_memory(v, data, v.getSizeX(), v.getSizeY(), v.getSizeZ(), v.getStrideX(), v.getStrideY(), v.getStrideZ());
}












/****************************************************************************//**
 * \brief Checks whether two volumes have the same dimension.
 *******************************************************************************/
template<typename T> template<typename T2>
inline bool tom::Volume<T>::is_equal_size(const tom::Volume<T2> &v) const {
    return  this->getSizeX() == v.getSizeX() &&
            this->getSizeY() == v.getSizeY() &&
            this->getSizeZ() == v.getSizeZ();
}


/****************************************************************************//**
 * \brief Checks whether the volume has a cubic size.
 *******************************************************************************/
template<typename T>
inline bool tom::Volume<T>::is_equal_size(std::size_t s) const {
    return is_equal_size(s,s,s);
}


/****************************************************************************//**
 * \brief Checks whether the volume has a certain volume size.
 *******************************************************************************/
template<typename T>
inline bool tom::Volume<T>::is_equal_size(std::size_t x, std::size_t y, std::size_t z) const {
    return this->getSizeX()==x && this->getSizeY()==y && this->getSizeZ()==z;
}





/****************************************************************************//**
 * \brief Check whether this volume points to the same volume vol and size and stride have the same values.
 *******************************************************************************/
template<typename T>
inline bool tom::Volume<T>::is_equal_memory(const tom::Volume<T> &vol) const{
    return this->getSizeX() == vol.getSizeX() &&
            this->getSizeY() == vol.getSizeY() &&
            this->getSizeZ() == vol.getSizeZ() &&
            this->getStrideX() == vol.getStrideX() &&
            this->getStrideY() == vol.getStrideY() &&
            this->getStrideZ() == vol.getStrideZ() &&
            &this->get() == &vol.get();
}


/****************************************************************************//**
 * \brief Check whether this volume shares the memory with an other one.
 *
 * ATTENTION: It does only check whether both volumes use the same memory-structure
 * as it is shared by the \c share_memory. Otherwise a complete test whether two volumes
 * use any common elements is quite difficult due to the possibly different
 * memory layouts.
 *******************************************************************************/
template<typename T> template<typename T2>
bool tom::Volume<T>::is_shared_memory(const tom::Volume<T2> &vol) const {
    return volume_memory == vol.volume_memory;
}





/****************************************************************************//**
 * \brief A strong total order relation for volumes.
 *
 * \param[in] v1
 * \param[in] v2
 * \return True if \c v1 < \a v2. Defining a strong total order relation.
 *
 * This operator can be used for compare-operator as needed in associative
 * containers such as std::map. Even if it is not possible to use the tom::Volume
 * as key elements itself (because of the missing copy constructor), it can be
 * used inside a compare operator which uses pointers to volumes.\n
 * \see http://uw713doc.sco.com/en/SDK_c++/_Key_and_Value_Types.html.
 * Problems happen, if floating point volume contains NaNs!
 *******************************************************************************/
template<typename T>
inline bool tom::less_ptr_Volume(const tom::Volume<T> *v1, const tom::Volume<T> *v2) {
    if (v1 == v2) {
        return false;
    }
    if (v1) {
        if (v2) {
            const T *dv1 = &v1->get();
            const T *dv2 = &v2->get();
            if (v1->isContiguous() && v2->isContiguous()) {
                const std::size_t nv1 = v1->getSizeX()*v1->getSizeY()*v1->getSizeZ();
                const std::size_t nv2 = v2->getSizeX()*v2->getSizeY()*v2->getSizeZ();
                const std::size_t nvx = nv1<nv2 ? nv1 : nv2;
                register std::size_t i = 0;
                while (i<nvx && dv1[i]==dv2[i]) {
                    i++;
                }
                if (i < nvx) {
                    return dv1[i] < dv2[i];
                } else {
                    return nv1 < nv2;
                }
            } else {
                const std::size_t v1sx = v1->getSizeX();
                const std::size_t v1sy = v1->getSizeY();
                const std::size_t v1sz = v1->getSizeZ();
                const std::size_t v2sx = v2->getSizeX();
                const std::size_t v2sy = v2->getSizeY();
                const std::size_t v2sz = v2->getSizeZ();
                const std::size_t nv1 = v1sx*v1sy*v1sz;
                const std::size_t nv2 = v2sx*v2sy*v2sz;
                const std::size_t nvx = nv1<nv2 ? nv1 : nv2;
                std::size_t i=0;
                std::size_t i1x=0, i1y=0, i1z=0, i2x=0, i2y=0, i2z=0;
                while (i<nvx) {
                    assert(i1x<v1sx && i1y<v1sy && i1z<v1sz && i2x<v2sx && i2y<v2sy && i2z<v2sz);
                    const T &div1 = *reinterpret_cast<const T *>(reinterpret_cast<const char *>(dv1) + v1->getByteOffset(i1x,i1y,i1z));
                    const T &div2 = *reinterpret_cast<const T *>(reinterpret_cast<const char *>(dv2) + v2->getByteOffset(i2x,i2y,i2z));

                    if (div1 != div2) {
                        return div1 < div2;
                    }
                    i++;
                    i1x++;
                    i2x++;
                    if (i1x >= v1sx) { i1x = 0; i1y++; }
                    if (i1y >= v1sy) { i1y = 0; i1z++; }
                    if (i2x >= v2sx) { i2x = 0; i2y++; }
                    if (i2y >= v2sy) { i2y = 0; i2z++; }
                }
                return nv1 < nv2;
            }
        }
        // v1!=NULL && v2==NULL
        return false;
    }
    // v1==NULL && v2!=NULL
    return true;
}



#endif


