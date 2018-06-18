/***********************************************************************//**
 * \file volume.cpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1
 * \date    22.11.2007
 **************************************************************************/
#include <tom/volume.hpp>




#include <typeinfo>
#include <limits>
#include <iostream>
#include <cstring>
#include <cassert>
#include <cerrno>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/casts.hpp>
#include <boost/lambda/if.hpp>


#include <fftw3.h>



#include <tom/volume_fcn.hpp>


/*
#ifndef NDEBUG
#if FPU_SENDS_SIGFPE
#include <fpu_control.h>
namespace lc {
namespace debug {
namespace {
*/
/***************************************************************************//**
 * Dump core in case of an floating point exception (for debugging).
 * \see http://www.yosefk.com/blog/the-c-sucks-series-the-quest-for-the-entry-point.html
 *
 * \attention So you call this function somewhere during your program's initialization
 * sequence, and sure enough, computations producing NaN after the call to fpu_setup
 * result in core dumps. Then one day someone computes a NaN before the call to fpu_setup,
 * and you get a core dump the first time you try to use the FPU after that point.
 * Because that's how x86 maintains its "illegal operation" flags and that's how it uses
 * them to signal exceptions.
 ******************************************************************************/
 /*
struct fpu_setup {
    fpu_setup() {
        unsigned short cw;
        _FPU_GETCW(cw);
        cw &= ~_FPU_MASK_ZM;//Divide by zero
        cw &= ~_FPU_MASK_IM;//Invalid operation
        _FPU_SETCW(cw);
    }
} fpu_setup_init;
} } } // namespace lc::debug
#endif
#endif
*/


namespace tom {
/****************************************************************************//**
 * \brief Structure where every tom::Volume remembers whether to free its memory.
 *
 * Volume can share the same memory using the constructor. All these objects point
 * to the same st_volume_memory to know how many other objects shares the memory
 * with them, and whether and how to free it.
 *******************************************************************************/
struct st_volume_memory {
    void *ptr;                      /**< Pointer to the memory which can be passed to free. Must not be the same as the base pointer of its volume. */
    std::size_t size;               /**< The maximum size allocated in bytes, beginning from ptr. */
    bool free;                      /**< If true, the destructor of the last object frees the memory calling fcn_free. */
    int cnt;                        /**< The number of objects sharing the same memory. */
    void (*fcn_free)(void *ptr);    /**< A pointer to the deallocation function. if !free it is NULL, otherwise it must be set. */
};
}


/****************************************************************************//**
 * \brief Returns a structure containing the memory infos of the volume.
 *******************************************************************************/
template<typename T>
std::size_t tom::Volume<T>::get_memsize() const {
    return this->volume_memory->size - ((const char *)this->data - (const char *)this->volume_memory->ptr);
}





/* saves default malloc-function.
 * Returned by tom::getFcnMalloc() */
void *(*fcn_malloc_default)(size_t n ) = NULL;

/* saves default free-function.
 * Returned by tom::getFcnMalloc() */
void  (*fcn_free_default  )(void *ptr) = NULL;






/****************************************************************************//**
 * \brief Sets the default memory allocation functions.
 *
 * \param[in] fcn_malloc Function to allocate memory. Synopsis must be the same
 *   as for malloc() from ISO C89.
 * \param[in] fcn_free Function to de-allocate memory. Synopsis must be the same
 *   as for free() from ISO C89.
 *
 * Sets the allocation and de-allocation functions as used in some cases to
 * allocate memory on the heap. The user must take care, that they correspond,
 * that means, memory allocated by \c fcn_malloc can be freed calling \c fcn_free and
 * vice versa.
 * Both parameters can be set to \c NULL. In that case the methods from tom_volume.cpp
 * use the new[] and delete[] operaters from C++. Not that these call
 * the default constructor of T, while malloc only allocates the memory.\n
 * Upon program start these functions are initialized to \c NULL.
 * In case of failure malloc should return NULL or throw std::bad_alloc. free should
 * fail or throw exceptions.
 *******************************************************************************/
void tom::setFcnMem( void *(*fcn_malloc)(size_t n),
                void (*fcn_free)(void *ptr)) {
    if ((!fcn_malloc && fcn_free) || (fcn_malloc && !fcn_free)) {
        throw std::invalid_argument("You must specify both fcn_malloc() and fcn_free() or leave both unspecified.");
    }
    fcn_malloc_default = fcn_malloc;
    fcn_free_default = fcn_free;
}



/****************************************************************************//**
 * \brief Returns the default memory allocation function.
 *
 * Can be set by tom::setMemFcn() (even to NULL).
 *******************************************************************************/
void *(*(tom::getFcnMalloc()))(size_t n) {
    return fcn_malloc_default;
}

/****************************************************************************//**
 * \brief Returns the default memory de-allocation function.
 *
 * Can be set by tom::setMemFcn() (even to NULL).
 *******************************************************************************/
void (*(tom::getFcnFree()))(void *) {
    return fcn_free_default;
}




/****************************************************************************//**
 * \brief Returns true, if the objects destructor will free the memory.
 *
 * Of course the memory is only freed it there are not other volumes left which
 * use the same data. See also number_memory_used().
 *******************************************************************************/
template<typename T>
bool tom::Volume<T>::to_free_memory() const {
    return this->volume_memory->free;
}


/****************************************************************************//**
 * \brief Returns how many objects in total share this data.
 *
 * This value will always be larger or equal to 1. The destructor of the object
 * decrements this value. If to_free_memory() is true and number_memory_used()
 * is 0, the destructor calls the de-allocation function.
 *******************************************************************************/
template<typename T>
int tom::Volume<T>::number_memory_used() const {
    return this->volume_memory->cnt;
}


/****************************************************************************//**
 * \brief Destructor of the volume.
 *
 * The destructor of the object the number of objects sharing the same data.
 * If to_free_memory() is true and number_memory_used()
 * is 0, the destructor calls the de-allocation function.
 *******************************************************************************/
template<typename T>
tom::Volume<T>::~Volume() {
    //std::cout << "call destructor: " << std::endl;
    //this->printInfo("v(destr)");

    if (--this->volume_memory->cnt < 1) {
        if (this->volume_memory->free) {
            //std::cout << "  +++++ do free memory" << std::endl;
            this->volume_memory->fcn_free(this->volume_memory->ptr);
        }
        //std::cout << "  +++++ free volume_memory" << std::endl;
        delete this->volume_memory;
    }
}




/****************************************************************************//**
 *
 *******************************************************************************/
template <typename T>
void (*tom::Volume<T>::get_free_function())(void *ptr) {
    return this->volume_memory->fcn_free;
}




/****************************************************************************//**
 * \brief Initializes the volume by taking the data from an existing one.
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
 *******************************************************************************/
template<typename T> template<typename T2>
void tom::Volume<T>::share_memory(Volume<T2> &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez) {

    if (sizex<=0 || sizey<=0 || sizez<=0) {
        throw std::invalid_argument("tom::Volume::share_memory - Empty volume is not allowed.");
    }
    if (!tom_io_check_stride(sizeof(T), sizex, sizey, sizez, &stridex, &stridey, &stridez)) {
        throw std::invalid_argument("tom::Volume::share_memory - Stride parameter too small according to data-type and volume size. PyTom rather supports Fortran ordered memory. Maybe you have a C order memory?");
    }

    std::size_t i;

    const char *database = reinterpret_cast<char *>(v.volume_memory->ptr);
    char *datac = reinterpret_cast<char *>(data);
    if (!datac) {
        datac = reinterpret_cast<char *>(&v.get());
    }

    if ((  datac < database) ||
        (i=datac - database) > v.volume_memory->size ||
        (v.volume_memory->size - i) < ((sizez-1) *stridez + (sizey-1)*stridey + (sizex-1) * stridex + sizeof(T))) {
        throw std::out_of_range("tom::Volume::share_memory - New data lies outside the (known) allocated memory.");
    }


    this->volume_memory = v.volume_memory;
    this->volume_memory->cnt++;

    this->data = const_cast<T *>(reinterpret_cast<const T *>(datac));
    this->sizex = sizex;
    this->sizey = sizey;
    this->sizez = sizez;
    this->stridex = stridex;
    this->stridey = stridey;
    this->stridez = stridez;

}


/****************************************************************************//**
 * \brief Create an empty volume.
 * This constructor is private, since no empty volumes are allowed.
 * This constructor can be used, to create an non initialized volume within
 * a member function or friend.
 *******************************************************************************/
template<typename T>
tom::Volume<T>::Volume()
    :   data(NULL),
        sizex(0),
        sizey(0),
        sizez(0),
        stridex(0),
        stridey(0),
        stridez(0),
        volume_memory(NULL) {
}


/****************************************************************************//**
 * \brief Creates a volume from already allocated memory.
 *
 * \param[in] data A pointer to the volume.
 * \param[in] sizex The size of the volume. Must be positive.
 * \param[in] sizey The size of the volume. Must be positive.
 * \param[in] sizez The size of the volume. Must be positive.
 * \param[in] free Whether the last volume should free the allocated memory in its
 *   destructor.
 * \param[in] stridex The distance in memory between two x-elements. This is
 *   measured in bytes, thus for contiguous memory, this equals to sizeof(T).
 *   Setting to zero defaults to contiguous
 * \param[in] stridey The distance in memory between two neighoubing y-elements.
 *   in bytes. For contiguous memory, this equals to sizex*stridex.
 * \param[in] stridez Same for z.
 * \param[in] fcn_free The function to de-allocate the memory afterwards.
 *   \c NULL defaults to tom::getFcnFree(), which by itself defaults to the
 *   C++ delete[] T operator, if set to \c NULL. In that case the function
 *   which is returned from tom::getFcnFree() at construction time is used,
 *   not the one from deconstruction time.
 *
 * You are responsible, that the memory pointed by data is large enough to
 * hold all values. Further the pointer must be valid as long as the volume
 * exists, and must be freed by the volume itself (in case of \c free ) or
 * by the caller.\n
 * \c fcn_free must be the right function to deallocate the memory. If the
 * memory was allocated using the new[] operator, set \c fcn_free to \c NULL
 * and tom::setFcnMem() also.\n
 * If you want to share memory with an other volume, use the appropriate
 * constructor instead.
 *
 * If the constructor throws an exception, the volume is not freed!!
 * In this case the caller must ensure exception safety.
 *******************************************************************************/
template<typename T>
tom::Volume<T>::Volume(T *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez, bool free, void (*fcn_free)(void *ptr)) {

    if (!data || sizex<=0 || sizey<=0 || sizez<=0) {
        throw std::invalid_argument("Volume copy constructor. Source volume has size 0,0,0!");
    }

    if (!tom_io_check_stride(sizeof(T), sizex, sizey, sizez, &stridex, &stridey, &stridez)) {
        throw std::invalid_argument("Stride parameter too small according to data-type and volume size. PyTom rather supports Fortran ordered memory. Maybe you have a C order memory?");
    }

    if (free) {
        if (!fcn_free) {
            if (!(fcn_free = tom::getFcnFree())) {
                fcn_free = &tom::fcn_free_delete<T>;
            }
        }
    } else {
        fcn_free = NULL;
    }

    this->volume_memory = new st_volume_memory();
    this->volume_memory->ptr = data;
    this->volume_memory->size = (sizez-1) *stridez + (sizey-1)*stridey + sizex * stridex;
    this->volume_memory->free = free;

    this->volume_memory->cnt = 1;
    this->volume_memory->fcn_free = fcn_free;

    this->data = data;
    this->sizex = sizex;
    this->sizey = sizey;
    this->sizez = sizez;
    this->stridex = stridex;
    this->stridey = stridey;
    this->stridez = stridez;
}




/****************************************************************************//**
 * \brief Initialize volume by allocating new (continous) memory.
 *******************************************************************************/
template<typename T>
void tom::Volume<T>::initialize(std::size_t sizex, std::size_t sizey, std::size_t sizez, void *(*fcn_malloc)(size_t size), void (*fcn_free)(void *ptr)) {
	//std::cout << "INIT" << std::endl;
    if (sizex<=0 || sizey<=0 || sizez<=0) {
        throw std::invalid_argument("Initialize of volume failed. Specified size is 0,0,0!");
    }

    if (!fcn_malloc && !fcn_free) {
        fcn_malloc = tom::getFcnMalloc();
        fcn_free = tom::getFcnFree();
        if (!fcn_malloc || !fcn_free) {
            fcn_malloc = &tom::fcn_malloc_new<T>;
            fcn_free = &tom::fcn_free_delete<T>;
        }
    } else if (!(fcn_malloc && fcn_free)) {
        throw std::invalid_argument("You must specify both fcn_malloc() and fcn_free() or leave both unspecified.");
    }

    size_t numel_bytes = sizex*sizey*sizez*sizeof(T);

    try {
        data = (T *)fcn_malloc(numel_bytes);
    } catch (...) {
        throw std::bad_alloc();
    }
    if (!data) {
        throw std::bad_alloc();
    }

    try {
        this->volume_memory = new st_volume_memory();
    } catch (...) {
        fcn_free(data);
        throw; /* Rethrow exception. */
    }

    this->volume_memory->ptr = data;
    this->volume_memory->size = numel_bytes;
    this->volume_memory->free = true;
    this->volume_memory->cnt = 1;
    this->volume_memory->fcn_free = fcn_free;

    this->data = data;
    this->sizex = sizex;
    this->sizey = sizey;
    this->sizez = sizez;
    this->stridex = sizeof(T);
    this->stridey = this->sizex * this->stridex;
    this->stridez = this->sizey * this->stridey;

}


/****************************************************************************//**
 * \brief Writes the volume as em-file.
 *
 * \param[in] filename The name of the em-file.
 * \param[in] header A pointer to the used em-header. If set to NULL, a default
 *   header is used with only the needed fields set according to the calling
 *   object. Even if a header is given, the fields \c dims and \c type are set
 *   to match the data from the volume.
 *
 * Calles tom_io_em_write to write the data. If an error occures
 * there, its integer error status is thrown as exception.
 *
 * This method is misplaced in this class. Use the functions from tom::io instead.
 *
 * \deprecated Don't use these methods, use the tom::io::write_to_em instead.
 *******************************************************************************/
template<typename T>
void tom::Volume<T>::write_to_em(const std::string &filename, const tom_io_em_header *header) const {
    // DEPRECATED use tom::io::write_to_em instead.
    #if 1
    tom_io_em_header header_local;

    if (this->getSizeX() > std::numeric_limits<uint32_t>::max() ||
        this->getSizeY() > std::numeric_limits<uint32_t>::max() ||
        this->getSizeZ() > std::numeric_limits<uint32_t>::max()) {
        throw std::runtime_error("The volume is too large to save it as em-format.");
    }


    if (header) {
        header_local = *header;
    } else {
        memset(&header_local, 0, sizeof(header_local));
        header_local.machine = 6;
    }

    if (!tom_io_em_set_iotype_header(&header_local, tom::get_tom_io_type<T>())) {
        throw std::runtime_error("This type can not be saved to em-file.");
    }
    header_local.dims[0] = this->getSizeX();
    header_local.dims[1] = this->getSizeY();
    header_local.dims[2] = this->getSizeZ();

    int res  = tom_io_em_write(filename.c_str(), &header_local, this->data, tom_io_em_get_iotype_header(&header_local), this->getStrideX(), this->getStrideY(), this->getStrideZ());
    if (res != TOM_ERR_OK) {
        const int errno_ = errno;
        std::stringstream ss;
        std::vector<char> buf(1024);
        ss << "Error writing em-file \"" << filename << "\": " << tom_io_strerror(res, errno_, &buf[0], buf.size());
        throw std::runtime_error(ss.str());
    }
    #else
    tom::io::write_to_em<T>(*this, filename, header);
    #endif
}

/****************************************************************************//**
 * \brief Writes the volume into an existing em-file.
 *
 * \param[in] filename The name of the em-file.
 * \param[in] header A pointer to the used em-header. If set to NULL, a default
 *   header is used with only the needed fields set according to the calling
 *   object. Even if a header is given, the fields \c dims and \c type are set
 *   to match the data from the volume.
 * \param[in] first_voxel position where the volume is written to.
 *
 * Calles tom_io_em_write to write the data. If an error occures
 * there, its integer error status is thrown as exception.\n
 * Uses tom_io_em_write_paste to write the volume into an existing em-file.
 *
 * This method is misplaced in this class. Use the functions from tom::io instead.
 *
 * \deprecated Don't use these methods, use tom::io::write_to_em_paste instead.
 *******************************************************************************/
template<typename T>
void tom::Volume<T>::write_to_em(const std::string &filename, const tom_io_em_header *header,const uint32_t* first_voxel) const {

    tom_io_em_header header_local;

    if (this->getSizeX() > std::numeric_limits<uint32_t>::max() ||
        this->getSizeY() > std::numeric_limits<uint32_t>::max() ||
        this->getSizeZ() > std::numeric_limits<uint32_t>::max()) {
        throw std::runtime_error("tom::Volume<T>::write_to_em-paste : The volume is too large to save it as em-format.");
    }


    if (header) {
        header_local = *header;
    } else {
        memset(&header_local, 0, sizeof(header_local));
        header_local.machine = 6;
    }

    if (!tom_io_em_set_iotype_header(&header_local, tom::get_tom_io_type<T>())) {
        throw std::runtime_error("tom::Volume<T>::write_to_em-paste : This type can not be saved to em-file.");
    }
    header_local.dims[0] = this->getSizeX();
    header_local.dims[1] = this->getSizeY();
    header_local.dims[2] = this->getSizeZ();

    int header_read;

    int res  = tom_io_em_write_paste(filename.c_str(),this->data,tom_io_em_get_iotype_header(&header_local),
                                    this->getSizeX(), this->getSizeY(), this->getSizeZ(),
                                    this->getStrideX(), this->getStrideY(),this->getStrideZ(),
                                    &header_local, &header_read, true, first_voxel, NULL);

    if (res != TOM_ERR_OK) {
        throw res;
    }
}


/****************************************************************************//**
 * \brief A malloc function using internally the C++ new[] operator.
 *
 * \param[in] size The number of bytes to be allocated.
 *
 * Creates <tt> ceil(size/sizeof(T)) </tt> objects of \c T with its default
 * constructor. Which should make not a big difference for elementary datatypes.
 *******************************************************************************/
template<typename T> void *tom::fcn_malloc_new(size_t size) {
    size_t remainder = size % sizeof(T);
    size /= sizeof(T);
    if (remainder > 0) {
        size++;
    }
    return new T[size];
}

/****************************************************************************//**
 * \brief A free function using internally the C++ delete[] operator.
 *******************************************************************************/
template<typename T> void tom::fcn_free_delete(void *ptr) {
    delete[] (T *)ptr;
}




/****************************************************************************//**
 * \brief Writes debugging info to stdout.
 *******************************************************************************/
template<typename T>
void tom::Volume<T>::printInfo() const {
    std::cout << "Volume unnamed:" << std::endl <<
                 "  size[zyx] = [" << this->getSizeZ()   << "," << this->getSizeY()   << "," << this->getSizeX()   << "]; of tom_io_type " << tom::get_tom_io_type<T>() << std::endl <<
                 "  stride[zyx] = [" << this->getStrideZ() << "," << this->getStrideY() << "," << this->getStrideX() << "]; gap[zyx] = [" << (this->getStrideZ()-this->getStrideY()*this->getSizeY()) << "," << (this->getStrideY()-this->getStrideX()*this->getSizeX()) << "," << (this->getStrideX()-sizeof(T)) << "]; " << std::endl;
    if (tom::is_double<T>() || tom::is_float<T>()) {
        double mean, variance, variance_false;
        T min, max;
        //tom::stat(*this, mean, variance, true, min, max);
        this->stat(&mean, &variance, &min, &max, true);
        this->stat(0, &variance_false, 0, 0, false);
        std::cout << "  range = [" << min << " .. " << mean << " .. " << max << "]" << std::endl
                  << "  stddev = " << sqrt(variance_false) << " (stddev_sample = " << sqrt(variance) << ")" << std::endl;
    }
    std::cout << "  ptr = " << &this->get() << ", baseptr = " << this->volume_memory->ptr << " (offset=" << (((size_t)&this->get()) - ((size_t)(this->volume_memory->ptr))) << ")" << std::endl <<
                 "  pointed to and ";
    if (!this->to_free_memory()) { std::cout << "not "; }
    std::cout << "to free by " << this->number_memory_used() << " volumes using " << ((void *)((size_t)this->volume_memory->fcn_free)) << " (";
    void (*fcn_free_delete_double)(void*) = &tom::fcn_free_delete<double>;
    void (*fcn_free_delete_float)(void*) = &tom::fcn_free_delete<float>;
    if (this->volume_memory->fcn_free == &fftw_free) {
        std::cout << "fftw_free";
    } else if (this->volume_memory->fcn_free == &fftwf_free) {
        std::cout << "fftwf_free";
    } else if (this->volume_memory->fcn_free == fcn_free_delete_double) {
        std::cout << "fcn_free_delete<double>";
    } else if (this->volume_memory->fcn_free == fcn_free_delete_float) {
        std::cout << "fcn_free_delete<float>";
    } else if (reinterpret_cast<void (*)(void)>(this->volume_memory->fcn_free) == reinterpret_cast<void (*)(void)>(&free)) {
        std::cout << "free";
    } else if (!this->volume_memory->fcn_free) {
        std::cout << "NULL";
    } else {
        std::cout << "unknown";
    }
    std::cout << ")" << std::endl;
}


/****************************************************************************//**
 * \brief Writes debugging info to stdout.
 *******************************************************************************/
template<typename T>
void tom::Volume<T>::printInfo(const std::string &name) const {
    std::cout << "Volume " << name << ":" << std::endl <<
                 "  size[zyx] = [" << this->getSizeZ()   << "," << this->getSizeY()   << "," << this->getSizeX()   << "]; of tom_io_type " << tom::get_tom_io_type<T>() << std::endl <<
                 "  stride[zyx] = [" << this->getStrideZ() << "," << this->getStrideY() << "," << this->getStrideX() << "]; gap[zyx] = [" << (this->getStrideZ()-this->getStrideY()*this->getSizeY()) << "," << (this->getStrideY()-this->getStrideX()*this->getSizeX()) << "," << (this->getStrideX()-sizeof(T)) << "]; " << std::endl;
    if (tom::is_double<T>() || tom::is_float<T>()) {
        double mean, variance, variance_false;
        T min, max;
        //tom::stat(*this, mean, variance, true, min, max);
        this->stat(&mean, &variance, &min, &max, true);
        this->stat(0, &variance_false, 0, 0, false);
        std::cout << "  range = [" << min << " .. " << mean << " .. " << max << "]" << std::endl
                  << "  stddev = " << sqrt(variance_false) << " (stddev_sample = " << sqrt(variance) << ")" << std::endl;
    }
    std::cout << "  ptr = " << &this->get() << ", baseptr = " << this->volume_memory->ptr << " (offset=" << (((size_t)&this->get()) - ((size_t)(this->volume_memory->ptr))) << ")" << std::endl <<
                 "  pointed to and ";
    if (!this->to_free_memory()) { std::cout << "not "; }
    std::cout << "to free by " << this->number_memory_used() << " volumes using " << ((void *)((size_t)this->volume_memory->fcn_free)) << " (";
    void (*fcn_free_delete_double)(void*) = &tom::fcn_free_delete<double>;
    void (*fcn_free_delete_float)(void*) = &tom::fcn_free_delete<float>;
    if (this->volume_memory->fcn_free == &fftw_free) {
        std::cout << "fftw_free";
    } else if (this->volume_memory->fcn_free == &fftwf_free) {
        std::cout << "fftwf_free";
    } else if (this->volume_memory->fcn_free == fcn_free_delete_double) {
        std::cout << "fcn_free_delete<double>";
    } else if (this->volume_memory->fcn_free == fcn_free_delete_float) {
        std::cout << "fcn_free_delete<float>";
    } else if (reinterpret_cast<void (*)(void)>(this->volume_memory->fcn_free) == reinterpret_cast<void (*)(void)>(&free)) {
        std::cout << "free";
    } else if (!this->volume_memory->fcn_free) {
        std::cout << "NULL";
    } else {
        std::cout << "unknown";
    }
    std::cout << ")" << std::endl;
}






namespace {
template<typename T> inline void setValues__contiguous_copy(T                    *dest, const T                    *src, std::size_t n) { for (std::size_t i=0; i<n; i++) { dest[i] = src[i]; } } // Dont use memcpy for calling the assignment operator.
template<          > inline void setValues__contiguous_copy(int                  *dest, const int                  *src, std::size_t n) { memcpy(dest, src, n*sizeof(int                 )); }
template<          > inline void setValues__contiguous_copy(float                *dest, const float                *src, std::size_t n) { memcpy(dest, src, n*sizeof(float               )); }
template<          > inline void setValues__contiguous_copy(double               *dest, const double               *src, std::size_t n) { memcpy(dest, src, n*sizeof(double              )); }
template<          > inline void setValues__contiguous_copy(std::complex<float > *dest, const std::complex<float > *src, std::size_t n) { memcpy(dest, src, n*sizeof(std::complex<float >)); }
template<          > inline void setValues__contiguous_copy(std::complex<double> *dest, const std::complex<double> *src, std::size_t n) { memcpy(dest, src, n*sizeof(std::complex<double>)); }
template<typename T1, typename T2> inline void for_each__tom__volume__setValues2_(T1 &v1, const T2 &v2) { v1 = v2; }
template<typename T1, typename T2>
struct for_each__tom__volume__setValues2 {
    inline void operator()(T1 &v1, const T2 &v2, std::size_t, std::size_t, std::size_t) { for_each__tom__volume__setValues2_<T1,T2>(v1, v2); }
};
}
/****************************************************************************//**
 * \brief Copy the values from one Volume to an other with type conversion.
 *
 * \param[in] v The source volume, of the same size as *this. The datatypes can
 *   differ.
 *******************************************************************************/
template<typename T>
template<typename T2>
void tom::Volume<T>::setValues(const Volume<T2> &v) {
    if (!this->is_equal_size(v)) {
        throw std::invalid_argument("Both volumes must have the same size.");
    }
    if (static_cast<const void *>(&v.get()) == static_cast<void *>(&this->get())) {
        if (typeid(T) != typeid(T2)) {
            throw std::invalid_argument("Self assigment of different datatypes.");
        }
        return;
    }

    if (typeid(T)==typeid(T2) && this->isContiguous() && v.isContiguous()) {
        ::setValues__contiguous_copy<T>(&this->get(), reinterpret_cast<const T *>(&v.get()), this->getSizeX()*this->getSizeY()*this->getSizeZ());
    } else {
        tom::loop::for_each(*this, v, ::for_each__tom__volume__setValues2<T, T2>());
    }
}







namespace {
template<typename T>
struct for_each__tom__volume__setValues {
    const T &_val;
    inline for_each__tom__volume__setValues(const T &val): _val(val) { }
    inline void operator()(T &a, std::size_t, std::size_t, std::size_t) const { a = _val; }
};
}
/****************************************************************************//**
 * \brief Set all voxels to a specific value.
 *
 * \param[in] val The value which is assigned to each voxel of this.
 *******************************************************************************/
template<typename T>
void tom::Volume<T>::setValues(T val) {
    tom::loop::for_each(*this, ::for_each__tom__volume__setValues<T>(val));
}




namespace lc {
namespace functor {
namespace {
template<typename T>
struct operator_equal_value {
    operator_equal_value(const T &val): val_(val), equal_(true) { }
    const T val_;
    bool equal_;
    bool operator()(const T &a, const tom::loop::index_type &, const tom::loop::index_type &, const tom::loop::index_type &) {
        return (equal_ = (a == val_));
    }
};
} } } // namespace lc::functor::<unnamed>
/****************************************************************************//**
 * \brief Checks wether all elements of the volume equal one value.
 *
 * \param[in] val The value to be compared with the volume.
 * \returns true if all elements of the volume are equal the scalar.
 *   Otherwise false.
 *
 * If the first element not equal val is found, false is returned. otherwise
 * every value has to be checked.
 *******************************************************************************/
template<typename T>
bool tom::Volume<T>::operator==(const T &val) const {
    lc::functor::operator_equal_value<T> s(val);
    tom::loop::for_each_while<const tom::Volume<T>, lc::functor::operator_equal_value<T> &>(*this, s);
    return s.equal_;
}




namespace tom {
namespace functor {
namespace {
template<typename T>
struct operator_equal_ {
    operator_equal_(): equal_(true) { }
    bool equal_;
    bool operator()(const T &v1, const T &v2, std::size_t, std::size_t, std::size_t) {
        if (v1 != v2) {
            equal_ = false;
        }
        return equal_;
    }
};
} // namespace
} // namespace functor
} // namespace tom
/****************************************************************************//**
 * \brief Checks wether all elements are equal to each other.
 *
 * \param[in] vol Volume to be compared.
 * \returns true if all elements of the volume equal.
 *   Otherwise false. If the volumes have different size, false is returned.
 *
 * If the first element not equal val is found, false is returned. otherwise
 * every value has to be checked.
 *******************************************************************************/
template<typename T>
bool tom::Volume<T>::operator==(const tom::Volume<T> &v) const {
    if (!this->is_equal_size<T>(v)) {
        return false;
    }
    if (this==&v || (&this->get()==&v.get() && this->getStrideX()==v.getStrideX() && this->getStrideY()==v.getStrideY() && this->getStrideZ()==v.getStrideZ())) {
        return true;
    }
    tom::functor::operator_equal_<T> f;
    tom::loop::for_each_while<const Volume<T>, const Volume<T>, tom::functor::operator_equal_<T> &>(*this, v, f);
    return f.equal_;
}






namespace {
template<typename T, typename TPRECISION>
struct for_each__tom__volume__stat__minmaxsum2 {
    for_each__tom__volume__stat__minmaxsum2(T init): min(init), max(init), sum(0), sum2(0) { }
    T min, max;
    TPRECISION sum, sum2;
    inline void operator()(const T &a, std::size_t, std::size_t, std::size_t) {
        if (a < this->min) {
            this->min = a;
        } else if (a > this->max) {
            this->max = a;
        }
        this->sum += a;
        this->sum2 += a*a;
    }
};
template<typename T, typename TPRECISION>
struct for_each__tom__volume__stat__minmaxsum {
    for_each__tom__volume__stat__minmaxsum(T init): min(init), max(init), sum(0) { }
    T min, max;
    TPRECISION sum;
    inline void operator()(const T &a, std::size_t, std::size_t, std::size_t) {
        if (a < this->min) {
            this->min = a;
        } else if (a > this->max) {
            this->max = a;
        }
        this->sum += a;
    }
};
template<typename T>
struct for_each__tom__volume__stat__minmax {
    for_each__tom__volume__stat__minmax(T init): min(init), max(init) { }
    T min, max;
    inline void operator()(const T &a, std::size_t, std::size_t, std::size_t) {
        if (a < this->min) {
            this->min = a;
        } else if (a > this->max) {
            this->max = a;
        }
    }
};
template<typename T>
struct for_each__tom__volume__stat__min {
    for_each__tom__volume__stat__min(T init): min(init) { }
    T min;
    inline void operator()(const T &a, std::size_t, std::size_t, std::size_t) { if (a < this->min) { this->min = a; } }
};
template<typename T>
struct for_each__tom__volume__stat__max {
    for_each__tom__volume__stat__max(T init): max(init) { }
    T max;
    inline void operator()(const T &a, std::size_t, std::size_t, std::size_t) { if (a > this->max) { this->max = a; } }
};
template<typename T, typename TPRECISION>
struct for_each__tom__volume__stat__sum {
    for_each__tom__volume__stat__sum(): sum(0) { }
    TPRECISION sum;
    inline void operator()(const T &a, std::size_t, std::size_t, std::size_t) { this->sum += a; }
};
template<typename T, typename TPRECISION>
struct for_each__tom__volume__stat__sum2 {
    for_each__tom__volume__stat__sum2(): sum(0), sum2(0) { }
    TPRECISION sum, sum2;
    inline void operator()(const T &a, std::size_t, std::size_t, std::size_t) { this->sum += a; this->sum2 += a*a; }
};
}
/****************************************************************************//**
 * \brief Private method which computes the mean, variance and min/max.
 *
 * Depending on which parameter are needed, different things are calculated.
 *******************************************************************************/
template<typename T>
void tom::Volume<T>::stat(double *m, double *variance, T *min, T *max, bool use_sample_standard_deviation) const {

    typedef double TPRECISION;

    const std::size_t numel = this->getSizeX() * this->getSizeY() * this->getSizeZ();

    TPRECISION lmean = 0.;
    TPRECISION lvariance = 0.;
    T lmin, lmax;
    lmin = lmax = this->get();

    if (m && !variance && !min && !max) {
        ::for_each__tom__volume__stat__sum<T, TPRECISION> s;
        tom::loop::for_each<const tom::Volume<T>, ::for_each__tom__volume__stat__sum<T, TPRECISION> &>(*this, s);
        lmean = s.sum / static_cast<TPRECISION>(numel);
    } else if (variance && !min && !max) {
        ::for_each__tom__volume__stat__sum2<T, TPRECISION> s;
        tom::loop::for_each<const tom::Volume<T>, ::for_each__tom__volume__stat__sum2<T, TPRECISION> &>(*this, s);
        lmean = s.sum / static_cast<TPRECISION>(numel);
        lvariance = (s.sum2 - static_cast<TPRECISION>(numel)*lmean*lmean) / (static_cast<TPRECISION>(numel)- (use_sample_standard_deviation ? 1 : 0));
    } else if (!m && !variance && min && !max) {
        ::for_each__tom__volume__stat__min<T> s(this->get());
        tom::loop::for_each<const tom::Volume<T>, ::for_each__tom__volume__stat__min<T> &>(*this, s);
        lmin = s.min;
    } else if (!m && !variance && !min && max) {
        ::for_each__tom__volume__stat__max<T> s(this->get());
        tom::loop::for_each<const tom::Volume<T>, ::for_each__tom__volume__stat__max<T> &>(*this, s);
        lmax = s.max;
    } else if (!m && !variance && min && max) {
        ::for_each__tom__volume__stat__minmax<T> s(this->get());
        tom::loop::for_each<const tom::Volume<T>, ::for_each__tom__volume__stat__minmax<T> &>(*this, s);
        lmin = s.min;
        lmax = s.max;
    } else if (!variance) {
        ::for_each__tom__volume__stat__minmaxsum<T, TPRECISION> s(this->get());
        tom::loop::for_each<const tom::Volume<T>, ::for_each__tom__volume__stat__minmaxsum<T, TPRECISION> &>(*this, s);
        lmin = s.min;
        lmax = s.max;
        lmean = s.sum / static_cast<TPRECISION>(numel);
    } else if (m || variance || min || max) {
        ::for_each__tom__volume__stat__minmaxsum2<T, TPRECISION> s(this->get());
        tom::loop::for_each<const tom::Volume<T>, ::for_each__tom__volume__stat__minmaxsum2<T, TPRECISION> &>(*this, s);
        lmin = s.min;
        lmax = s.max;
        lmean = s.sum / static_cast<TPRECISION>(numel);
        lvariance = (s.sum2 - static_cast<TPRECISION>(numel)*lmean*lmean) / (static_cast<TPRECISION>(numel)- (use_sample_standard_deviation ? 1 : 0));
    }
    if (max) { *max = lmax; }
    if (min) { *min = lmin; }
    if (m) { *m = lmean; }
    if (variance) { *variance = lvariance > 0. ? lvariance : 0.; }
}
/** \cond __HIDE_FUNCTIONS */
namespace tom {
template<> void Volume<std::complex<float > >::stat(double *m, double *variance, std::complex<float > *min, std::complex<float > *max, bool use_sample_standard_deviation) const { throw std::runtime_error("function stat not defined for complex volume"); }
template<> void Volume<std::complex<double> >::stat(double *m, double *variance, std::complex<double> *min, std::complex<double> *max, bool use_sample_standard_deviation) const { throw std::runtime_error("function stat not defined for complex volume"); }
}
/** \endcond __HIDE_FUNCTIONS */





/****************************************************************************//**
 * \brief Computes the variance and the mean of a volume under a given mask.
 *
 * \param[out] mean Here the mean under the mask will be returned.
 * \param[out] variance Here the variance under the mask will be returned.
 * \param[in]  use_sample_standard_deviation If true, the variance is computed
 *    using the standard deviation of the sample (with nominator N-1).
 * \param[in] mask Only those voxels in the source volume are considered where
 *    the mask has values != 0.
 *    However the mean and the variance are normalised relative to the whole image!!
 *******************************************************************************/
template<typename T>
template<typename T2>
void tom::Volume<T>::stat(double &m, double &variance, bool use_sample_standard_deviation, const tom::Volume<T2> &mask) const {
    tom::stat(*this, m, variance, use_sample_standard_deviation, mask);
}





namespace {
template<typename T, typename TPRECISION>
struct for_each__tom__volume__variance_mean_free {
    for_each__tom__volume__variance_mean_free(): sum2(0) { }
    TPRECISION sum2;
    inline void operator()(const T &a, std::size_t, std::size_t, std::size_t) { this->sum2 += a*a; }
};
}
/****************************************************************************//**
 * \brief Computes the variance of a volume which has a mean of zero.
 *
 * \param[in] use_sample_standard_deviation If true, the variance is computed
 *    using the standard deviation of the sample (with nominator N-1).
 *
 * \return The variance of the volume.
 *
 * If the mean of the volume is not 0, the result is wrong. This case is not
 * checked.
 *******************************************************************************/
template<typename T>
double tom::Volume<T>::variance_mean_free(bool use_sample_standard_deviation) const {

    typedef double TPRECISION;

    ::for_each__tom__volume__variance_mean_free<T, TPRECISION> s;
    tom::loop::for_each<const tom::Volume<T>, ::for_each__tom__volume__variance_mean_free<T, TPRECISION> &>(*this, s);

    return s.sum2 / static_cast<TPRECISION>(this->getSizeX() * this->getSizeY() * this->getSizeZ() - (use_sample_standard_deviation ? 1 : 0));
}
/** \cond __HIDE_FUNCTIONS */
namespace tom {
template<> double Volume<std::complex<float > >::variance_mean_free(bool use_sample_standard_deviation) const { throw std::runtime_error("function stat not defined for complex volume"); }
template<> double Volume<std::complex<double> >::variance_mean_free(bool use_sample_standard_deviation) const { throw std::runtime_error("function stat not defined for complex volume"); }
}
/** \endcond __HIDE_FUNCTIONS */



/****************************************************************************//**
 * \brief Compute the minimum and the maximum of the volume.
 *******************************************************************************/
template<typename T>
void tom::Volume<T>::minmax(T &minimum, T &maximum) const {
    this->stat(NULL, NULL, &minimum, &maximum, false);
    //tom::minmax(*this, minimum, maximum);
}
/****************************************************************************//**
 * \brief Return the mean and the variance of the volume
 *******************************************************************************/
template <typename T>
void tom::Volume<T>::stat(double &m, double &variance, bool use_sample_standard_deviation) const {
    this->stat(&m, &variance, NULL, NULL, use_sample_standard_deviation);
    //tom::stat(*this, m, variance, use_sample_standard_deviation);
}
/****************************************************************************//**
 * \brief Return the variance, the mean and the min, max of the volume.
 *******************************************************************************/
template <typename T>
void tom::Volume<T>::stat(double &m, double &variance, T &min, T &max, bool use_sample_standard_deviation) const {
    this->stat(&m, &variance, &min, &max, use_sample_standard_deviation);
    //tom::stat(*this, m, variance, use_sample_standard_deviation, min, max);
}
/****************************************************************************//**
 * \brief Returns the minumum value.
 *******************************************************************************/
template <typename T>
T tom::Volume<T>::min() const {
    T min; this->stat(NULL, NULL, &min, NULL, false); return min;
    //return tom::min(*this);
}
/****************************************************************************//**
 * \brief Return the maximum value.
 *******************************************************************************/
template <typename T>
T tom::Volume<T>::max() const {
    T max; this->stat(NULL, NULL, NULL, &max, false); return max;
    //return tom::max(*this);
}
/****************************************************************************//**
 * \brief Compute the variance of the volume
 *
 * \param[in] use_sample_standard_deviation If true the denominator is N-1 instead
 *    of N (N number of voxels).
 *
 * \return The variance of the volume.
 *******************************************************************************/
template <typename T>
double tom::Volume<T>::variance(bool use_sample_standard_deviation) const {
    double variance; this->stat(NULL, &variance, NULL, NULL, use_sample_standard_deviation); return variance;
    //return tom::variance(*this, use_sample_standard_deviation);
}
/****************************************************************************//**
 * \brief Return the mean of the volume
 *******************************************************************************/
template <typename T>
double tom::Volume<T>::mean() const {
    double mean; this->stat(&mean, NULL, NULL, NULL, false); return mean;
    //return tom::sum<T, double>(*this) / static_cast<double>(numel());
}











namespace tom { namespace functor { namespace {
template<typename T, typename TPR> struct shift {
    shift(const TPR &factor): factor_(factor) { }
    const TPR &factor_;
    void operator()(T &val, const tom::loop::packed_idx &) {
        val = static_cast<TPR>(val) + factor_;
    }
};
template<typename T, typename TPR> struct shift<std::complex<T>, TPR> {
    shift(const TPR &factor): factor_(factor) { }
    const TPR &factor_;
    void operator()(std::complex<T> &val, const tom::loop::packed_idx &) {
        val = static_cast<std::complex<TPR> >(val) + factor_;
    }
};
template<typename T, typename TPR> struct shift<std::complex<T>, std::complex<TPR> > {
    shift(const std::complex<TPR> &factor): factor_(factor) { }
    const std::complex<TPR> &factor_;
    void operator()(std::complex<T> &val, const tom::loop::packed_idx &) {
        val = static_cast<std::complex<TPR> >(val) + factor_;
    }
};
} } } // namespace tom::functor::<unnamed>
/****************************************************************************//**
 * \brief Adds a scalar to the volume.
 *******************************************************************************/
template<typename T>
template<typename TPRECISION>
void tom::Volume<T>::shift(TPRECISION factor_shift) {
    if (!tom::equal(factor_shift, 0)) {
        tom::loop::for_each_packed_idx(*this, tom::functor::shift<T, TPRECISION>(factor_shift));
    }
}





namespace tom { namespace functor { namespace {
template<typename T, typename TPR> struct scale {
    scale(const TPR &factor): factor_(factor) { }
    const TPR &factor_;
    void operator()(T &val) {
        val = static_cast<TPR>(val) * factor_;
    };
};
template<typename T, typename TPR> struct scale<std::complex<T>, TPR> {
    scale(const TPR &factor): factor_(factor) { }
    const TPR &factor_;
    void operator()(std::complex<T> &val) {
        val = static_cast<std::complex<TPR> >(val) * factor_;
    };
};
template<typename T, typename TPR> struct scale<std::complex<T>, std::complex<TPR> > {
    scale(const std::complex<TPR> &factor): factor_(factor) { }
    const std::complex<TPR> &factor_;
    void operator()(std::complex<T> &val) {
        val = static_cast<std::complex<TPR> >(val) * factor_;
    };
};
} } } // namespace tom::functor::<unnamed>
/****************************************************************************//**
 * \brief Multiplies the volume with a factor.
 *******************************************************************************/
template<typename T>
template<typename TPRECISION>
void tom::Volume<T>::scale(TPRECISION factor_scale) {
    if (tom::equal(factor_scale, 0)) {
        setValues(0);
    } else if (!tom::equal(factor_scale, 1)) {
        tom::loop::for_each_no_idx(*this, tom::functor::scale<T, TPRECISION>(factor_scale));
    }
}


namespace tom { namespace functor { namespace {
template<typename T, typename TPR> struct shift_scale {
    shift_scale(const TPR &factor_shift, const TPR &factor_scale): factor_shift_(factor_shift), factor_scale_(factor_scale) { }
    const TPR factor_shift_, factor_scale_;
    void operator()(T &val, const tom::loop::packed_idx &) {
        val = (static_cast<TPR>(val) + factor_shift_) * factor_scale_;
    }
};
template<typename T, typename TPR> struct shift_scale<std::complex<T>, TPR> {
    shift_scale(const TPR &factor_shift, const TPR &factor_scale): factor_shift_(factor_shift), factor_scale_(factor_scale) { }
    const TPR &factor_shift_, &factor_scale_;
    void operator()(std::complex<T> &val, const tom::loop::packed_idx &) {
        val = (static_cast<std::complex<TPR> >(val) + factor_shift_) * factor_scale_;
    }
};
template<typename T, typename TPR> struct shift_scale<std::complex<T>, std::complex<TPR> > {
    shift_scale(const std::complex<TPR> &factor_shift, const std::complex<TPR> &factor_scale): factor_shift_(factor_shift), factor_scale_(factor_scale) { }
    const std::complex<TPR> &factor_shift_, &factor_scale_;
    void operator()(std::complex<T> &val, const tom::loop::packed_idx &) {
        val = (static_cast<std::complex<TPR> >(val) + factor_shift_) * factor_scale_;
    }
};
} } } // namespace tom::functor::<unnamed>
/****************************************************************************//**
 * \brief Adds and scales the volume by a scalar.
 *******************************************************************************/
template<typename T>
template<typename TPRECISION>
void tom::Volume<T>::shift_scale(TPRECISION factor_shift, TPRECISION factor_scale) {
    if (tom::equal(factor_scale, 0)) {
        setValues(0);
    } else if (!tom::equal(factor_shift, 0) && !tom::equal(factor_scale, 1)) {
        tom::loop::for_each_packed_idx(*this, tom::functor::shift_scale<T, TPRECISION>(factor_shift, factor_scale));
    } else if (!tom::equal(factor_shift, 0)) {
        this->shift(factor_shift);
    } else if (!tom::equal(factor_scale, 1)) {
        this->scale(factor_scale);
    }
}





namespace tom { namespace functor { namespace {
template<typename T, typename TPR> struct scale_shift {
    scale_shift(const TPR &factor_scale, const TPR &factor_shift): factor_shift_(factor_shift), factor_scale_(factor_scale) { }
    const TPR &factor_shift_, &factor_scale_;
    void operator()(T &val, const tom::loop::packed_idx &) {
        val = (static_cast<TPR>(val) * factor_scale_) + factor_shift_;
    }
};
template<typename T, typename TPR> struct scale_shift<std::complex<T>, TPR> {
    scale_shift(const TPR &factor_scale, const TPR &factor_shift): factor_shift_(factor_shift), factor_scale_(factor_scale) { }
    const TPR &factor_shift_, &factor_scale_;
    void operator()(std::complex<T> &val, const tom::loop::packed_idx &) {
        val = (static_cast<std::complex<TPR> >(val) * factor_scale_) + factor_shift_;
    }
};
template<typename T, typename TPR> struct scale_shift<std::complex<T>, std::complex<TPR> > {
    scale_shift(const std::complex<TPR> &factor_scale, const std::complex<TPR> &factor_shift): factor_shift_(factor_shift), factor_scale_(factor_scale) { }
    const std::complex<TPR> &factor_shift_, &factor_scale_;
    void operator()(std::complex<T> &val, const tom::loop::packed_idx &) {
        val = (static_cast<std::complex<TPR> >(val) * factor_scale_) + factor_shift_;
    }
};
} } } // namespace tom::functor::<unnamed>

/****************************************************************************//**
 * \brief First, multiplies the volume with a scalar and then adds a constant.
 *******************************************************************************/
template<typename T>
template<typename TPRECISION>
void tom::Volume<T>::scale_shift(TPRECISION factor_scale, TPRECISION factor_shift) {
    if (tom::equal(factor_scale, 0)) {
        setValues(factor_shift);
    } else if (!tom::equal(factor_shift, 0) && !tom::equal(factor_scale, 1)) {
        tom::loop::for_each_packed_idx(*this, tom::functor::scale_shift<T, TPRECISION>(factor_scale, factor_shift));
    } else if (!tom::equal(factor_shift, 0)) {
        this->shift(factor_shift);
    } else if (!tom::equal(factor_scale, 1)) {
        this->scale(factor_scale);
    }
}







// template instantiation.
/** \cond __HIDE_FUNCTIONS */

template class tom::Volume<char         >;
template class tom::Volume<int16_t      >;
template class tom::Volume<uint16_t     >;
template class tom::Volume<int          >;
template class tom::Volume<float        >;
template class tom::Volume<double       >;
template class tom::Volume<std::complex<float > >;
template class tom::Volume<std::complex<double> >;


//Copy Constructor
template tom::Volume<int16_t              >::Volume(Volume<int16_t              > &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez);
template tom::Volume<int                  >::Volume(Volume<int                  > &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez);
template tom::Volume<float                >::Volume(Volume<float                > &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez);
template tom::Volume<float                >::Volume(Volume<std::complex<float > > &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez);
template tom::Volume<double               >::Volume(Volume<double               > &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez);
template tom::Volume<double               >::Volume(Volume<std::complex<double> > &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez);
template tom::Volume<std::complex<float > >::Volume(Volume<std::complex<float > > &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez);
template tom::Volume<std::complex<double> >::Volume(Volume<std::complex<double> > &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez);

//Complex to Real Constructor
//template tom::Volume<float>::Volume(Volume<std::complex<float> > &v, bool get_real_part);
                                                    
                                                    
                                                    
template void tom::Volume<float >::stat<float >(double &m, double &variance, bool use_sample_standard_deviation, const tom::Volume<float > &mask) const;
template void tom::Volume<double>::stat<double>(double &m, double &variance, bool use_sample_standard_deviation, const tom::Volume<double> &mask) const;


template void tom::Volume<float                >::setValues<char                 >(const Volume<char                 > &v);
template void tom::Volume<float                >::setValues<int                  >(const Volume<int                  > &v);
template void tom::Volume<float                >::setValues<uint16_t             >(const Volume<uint16_t             > &v);
template void tom::Volume<int                  >::setValues<int                  >(const Volume<int                  > &v);
template void tom::Volume<int16_t              >::setValues<float                >(const Volume<float                > &v);
template void tom::Volume<float                >::setValues<int16_t              >(const Volume<int16_t              > &v);
template void tom::Volume<int                  >::setValues<float                >(const Volume<float                > &v);
template void tom::Volume<float                >::setValues<float                >(const Volume<float                > &v);
template void tom::Volume<float                >::setValues<double               >(const Volume<double               > &v);
template void tom::Volume<double               >::setValues<char                 >(const Volume<char                 > &v);
template void tom::Volume<double               >::setValues<int                  >(const Volume<int                  > &v);
template void tom::Volume<double               >::setValues<float                >(const Volume<float                > &v);
template void tom::Volume<double               >::setValues<double               >(const Volume<double               > &v);
template void tom::Volume<std::complex<float > >::setValues<std::complex<float > >(const Volume<std::complex<float > > &v);
template void tom::Volume<std::complex<float > >::setValues<std::complex<double> >(const Volume<std::complex<double> > &v);
template void tom::Volume<std::complex<double> >::setValues<std::complex<float > >(const Volume<std::complex<float > > &v);
template void tom::Volume<std::complex<double> >::setValues<std::complex<double> >(const Volume<std::complex<double> > &v);


template void tom::Volume<float                >::scale      (float  factor_scale);
template void tom::Volume<float                >::scale      (double factor_scale);
template void tom::Volume<double               >::scale      (double factor_scale);
template void tom::Volume<std::complex<float > >::scale      (float factor_scale);
template void tom::Volume<std::complex<float > >::scale      (double factor_scale);
template void tom::Volume<std::complex<double> >::scale      (double factor_scale);
template void tom::Volume<float                >::shift      (float  factor_shift);
template void tom::Volume<float                >::shift      (double factor_shift);
template void tom::Volume<std::complex<float>                >::shift       (std::complex<float> factor_shift);
template void tom::Volume<double               >::shift      (double factor_shift);
template void tom::Volume<std::complex<double>               >::shift      (std::complex<double> factor_shift);
template void tom::Volume<float                >::shift_scale(float  factor_shift, float  factor_scale);
template void tom::Volume<float                >::shift_scale(double factor_shift, double factor_scale);
template void tom::Volume<double               >::shift_scale(double factor_shift, double factor_scale);
template void tom::Volume<std::complex<float > >::shift_scale(float  factor_shift, float  factor_scale);
template void tom::Volume<std::complex<float > >::shift_scale(std::complex<float >  factor_shift, std::complex<float >  factor_scale);
template void tom::Volume<std::complex<float > >::shift_scale(double factor_shift, double factor_scale);
template void tom::Volume<std::complex<double> >::shift_scale(double factor_shift, double factor_scale);
template void tom::Volume<std::complex<double> >::shift_scale(std::complex<double> factor_shift, std::complex<double> factor_scale);
template void tom::Volume<float                >::scale_shift(float  factor_scale, float  factor_shift);
template void tom::Volume<float                >::scale_shift(double factor_scale, double factor_shift);
template void tom::Volume<double               >::scale_shift(double factor_scale, double factor_shift);


template void tom::Volume<char  >::share_memory<char  >(Volume<char  > &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez);
template void tom::Volume<int   >::share_memory<int   >(Volume<int   > &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez);
template void tom::Volume<double>::share_memory<float >(Volume<float > &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez);
template void tom::Volume<float >::share_memory<double>(Volume<double> &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez);



/** \endcond __HIDE_FUNCTIONS */

