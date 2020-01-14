/***********************************************************************//**
 * \file io.cpp
 * \brief IO-Functions. Contains implementations of io.hpp
 * \author  Thomas Haller
 * \version 0.2
 * \date    12.01.2009
 **************************************************************************/
#include <tom/io/io.hpp>


#include <sstream>
#include <cerrno>
#include <stdexcept>
#include <cstdio>
#include <cstring>
#include <limits>
#include <iostream>

#include <tom/volume_fcn.hpp>



namespace tom {
namespace io {
/****************************************************************************//**
 * \brief Returns an integer specifying the data-type.
 *
 * \returns Number as specified in \link tom_defines.h \endlink for defines
 *   such as TOM_IO_TYPE_DOUBLE. Returns TOM_IO_TYPE_UNKNOWN
 *   if the type is not defined there.
 *******************************************************************************/
template<typename T> int get_tom_io_type                       () { return TOM_IO_TYPE_UNKNOWN;   }
template<          > int get_tom_io_type<char                 >() { return TOM_IO_TYPE_INT8;      }
template<          > int get_tom_io_type<int8_t               >() { return TOM_IO_TYPE_INT8;      }
template<          > int get_tom_io_type<int16_t              >() { return TOM_IO_TYPE_INT16;     }
template<          > int get_tom_io_type<int32_t              >() { return TOM_IO_TYPE_INT32;     }
template<          > int get_tom_io_type<float                >() { return TOM_IO_TYPE_FLOAT;     }
template<          > int get_tom_io_type<fftwf_complex        >() { return TOM_IO_TYPE_COMPLEX32; }
template<          > int get_tom_io_type<std::complex<float>  >() { return TOM_IO_TYPE_COMPLEX32; }
template<          > int get_tom_io_type<double               >() { return TOM_IO_TYPE_DOUBLE;    }
template<          > int get_tom_io_type<fftw_complex         >() { return TOM_IO_TYPE_COMPLEX64; }
template<          > int get_tom_io_type<std::complex<double> >() { return TOM_IO_TYPE_COMPLEX64; }
template<          > int get_tom_io_type<int64_t              >() { return TOM_IO_TYPE_INT64;     }

template<          > int get_tom_io_type<std::complex<int16_t>  >() { return TOM_IO_TYPE_COMPLEX16; }
template<          > int get_tom_io_type<uint16_t             >() { return TOM_IO_TYPE_UINT16;    }
} // namespace io
} // namespace tom


/***************************************************************************//**
 * Read the header.
 ******************************************************************************/
void tom::io::em_read_header(const std::string &filename, tom_io_em_header &header) {

    const int tom_errnum = tom_io_em_read_header(filename.c_str(), 0, &header, 0);

    if (tom_errnum != TOM_ERR_OK) {
        const int errno_ = errno;
        std::ostringstream ss;
        char cbuf[2048];
        ss << "Error reading header of emfile \"" << filename << "\": " << tom_io_strerror(tom_errnum, errno_, cbuf, sizeof(cbuf));
        throw std::runtime_error(ss.str());
    }
}



/***************************************************************************//**
 * Constructor for the tom::io::Volume class
 *
 * \param[in] filename Name of the em-file.
 * \param[in] fcn_malloc Pointer to an malloc function used allocate the memory
 *      when reading the volume.
 * \param[in] fcn_free Pointer to the free function for freeing the memory.
 *      This must correspond to the \c fcn_malloc function.
 * \param[in] close_after_first_read The emfile is always opened uppon
 *      construction of this (with throwing an exception in case of an error).
 *      If this parameter is true, the file will be closed after reading the
 *      volume the first time. In that case, several calls of read_volume
 *      must each time reopen the file. If false, the file stays open until
 *      the object is destructed.
 *
 * Leaving both functions 0, defaults to tom::getFcnMalloc and tom::getFcnFree,
 * which, if not set, defaults to new T[]/ delete[] T
 ******************************************************************************/
tom::io::VolumeEM::VolumeEM(const std::string &filename, void *(*fcn_malloc)(size_t size), void (*fcn_free)(void *ptr), bool close_after_first_read)
    :   fcn_malloc_(fcn_malloc),
        fcn_free_(fcn_free),
        headerp_(new tom_io_em_header),
        header_(*headerp_),
        filename_(filename),
        f_(0),
        f_is_closed_(true),
        close_after_first_read_(close_after_first_read),
        seeked_at_data_start_(false) {
    if ((!fcn_malloc && fcn_free) || (fcn_malloc && !fcn_free)) {
        throw std::invalid_argument("Both fcn_malloc() and fcn_free() must be specified or leave unspecified (0).");
    }
    open_and_seek_file(true);
}


/***************************************************************************//**
 * Constructor for the tom::io::Volume class
 *
 * \param[in] filename Name of the em-file.
 * \param[in,out] header tom_io_em_header which is used for reading the emfile.
 *      This structure must exist as long as the VolumeEM exists and shall not
 *      be modified. Modifying or deleting the header is erroneous.
 * \param[in] fcn_malloc Pointer to an malloc function used allocate the memory
 *      when reading the volume.
 * \param[in] fcn_free Pointer to the free function for freeing the memory.
 *      This must correspond to the \c fcn_malloc function.
 * \param[in] close_after_first_read The emfile is always opened uppon
 *      construction of this (with throwing an exception in case of an error).
 *      If this parameter is true, the file will be closed after reading the
 *      volume the first time. In that case, several calls of read_volume
 *      must each time reopen the file. If false, the file stays open until
 *      the object is destructed.
 *
 * Leaving both functions 0, defaults to tom::getFcnMalloc and tom::getFcnFree,
 * which, if not set, defaults to new T[]/ delete[] T
 ******************************************************************************/
tom::io::VolumeEM::VolumeEM(const std::string &filename, tom_io_em_header &header, void *(*fcn_malloc)(size_t size), void (*fcn_free)(void *ptr), bool close_after_first_read)
    :   fcn_malloc_(fcn_malloc),
        fcn_free_(fcn_free),
        headerp_(0),
        header_(header),
        filename_(filename),
        f_(0),
        f_is_closed_(true),
        close_after_first_read_(close_after_first_read),
        seeked_at_data_start_(false) {
    if ((!fcn_malloc && fcn_free) || (fcn_malloc && !fcn_free)) {
        throw std::invalid_argument("Both fcn_malloc() and fcn_free() must be specified or leave unspecified (0).");
    }
    open_and_seek_file(true);
}


/***************************************************************************//**
 * Destructor for the tom::io::Volume class
 ******************************************************************************/
tom::io::VolumeEM::~VolumeEM() {
    if (!f_is_closed_) {
        std::fclose(f_);
    }
    delete headerp_;
}


/***************************************************************************//**
 * Opens the file.
 ******************************************************************************/
void tom::io::VolumeEM::open_and_seek_file(bool use_volume_header) {
    if (f_is_closed_) {
        std::auto_ptr<tom_io_em_header> header2;
        tom_io_em_header *header;
        if (use_volume_header) {
            header = &header_;
        } else {
            header2.reset(header = new tom_io_em_header());
        }
        const int tom_errnum = tom_io_em_read_header(filename_.c_str(), 0, header, &f_);
        if (tom_errnum != TOM_ERR_OK) {
            const int errno_ = errno;
            std::ostringstream ss;
            char cbuf[2048];
            ss << "Error opening emfile \"" << filename_ << "\": " << tom_io_strerror(tom_errnum, errno_, cbuf, sizeof(cbuf));
            throw std::runtime_error(ss.str());
        }
        f_is_closed_ = false;
        seeked_at_data_start_ = true;
        if (!use_volume_header && std::memcmp(header, &header_, sizeof(header_))!=0) {
            throw std::runtime_error("Re-reading the emheader yields different results.");
        }
    } else {
        if (!seeked_at_data_start_) {
            fseek(f_, sizeof(tom_io_em_header), SEEK_SET);
            seeked_at_data_start_ = true;
        }
    }
}


/***************************************************************************//**
 * Reads the volume from file and returns a newly created volume.
 ******************************************************************************/
template<typename T>
std::auto_ptr<tom::Volume<T> > tom::io::VolumeEM::read_volume(const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning) {

    std::auto_ptr<tom::Volume<T> > v;

    try {
        // Open/reopen the file and seek to the start of the data (after the header)
        open_and_seek_file(false);

        // Prepare the memory function.
        void *(*fcn_malloc)(size_t size) = fcn_malloc_;
        void (*fcn_free)(void *ptr)      = fcn_free_;
        assert((fcn_malloc && fcn_free) || (!fcn_malloc && !fcn_free));
        if (!fcn_malloc) {
            fcn_malloc = tom::getFcnMalloc();
            fcn_free = tom::getFcnFree();
            if (!fcn_malloc || !fcn_free) {
                fcn_malloc = &tom::fcn_malloc_new<T>;
                fcn_free = &tom::fcn_free_delete<T>;
            }
        }
        int i;
        uint32_t dims[3];
        uint32_t binning_[3]  = {1, 1, 1};
        uint32_t sampling_[3] = {1, 1, 1};
        uint32_t subregion_[6];

        /* Get the subregion and check the dimension. */
        if ((i=tom_io_calculate_sizes(header_.dims, subregion, subregion_, sampling, sampling_, binning, binning_, 0, dims)) != TOM_ERR_OK) {
            const int errno_ = errno;
            std::ostringstream ss;
            char cbuf[1024];
            ss << "Error reading emfile \"" << filename_ << "\": " << tom_io_strerror(i, errno_, cbuf, sizeof(cbuf));
            throw std::runtime_error(ss.str());
        }

        v.reset(new tom::Volume<T>(dims[0], dims[1], dims[2], fcn_malloc, fcn_free));

        if ((i=tom_io_read_vol(f_, header_.dims, tom_io_em_get_iotype_header(&header_), tom_io_em_is_swaped(&header_), subregion_, sampling_, binning_, &v->get(), tom::io::get_tom_io_type<T>())) != TOM_ERR_OK) {
            const int errno_ = errno;
            std::ostringstream ss;
            char cbuf[1024];
            ss << "Error reading emfile \"" << filename_ << "\": " << tom_io_strerror(i, errno_, cbuf, sizeof(cbuf));
            throw std::runtime_error(ss.str());
        }

    } catch (...) {
        std::fclose(f_);
        f_is_closed_ = true;
        seeked_at_data_start_ = false;
        throw;
    }

    if (close_after_first_read_) {
        // Close the file.
        std::fclose(f_);
        f_is_closed_ = true;
    }
    seeked_at_data_start_ = false;

    return v;
}





/***************************************************************************//**
 * Initializes a volume while reading it from file.
 * The volume \c v must have the proper size, otherwise an exception is thrown.
 ******************************************************************************/
template<typename T>
void tom::io::VolumeEM::read_volume(tom::Volume<T> &v, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning) {

    try {
        // Open/reopen the file and seek to the start of the data (after the header)
        open_and_seek_file(false);

        int i;
        uint32_t dims[3];
        uint32_t binning_[3]  = {1, 1, 1};
        uint32_t sampling_[3] = {1, 1, 1};
        uint32_t subregion_[6];

        std::ostringstream ss;

        /* Get the subregion and check the dimension. */
        if ((i=tom_io_calculate_sizes(header_.dims, subregion, subregion_, sampling, sampling_, binning, binning_, 0, dims)) != TOM_ERR_OK) {
            const int errno_ = errno;
            char cbuf[1024];
            ss << "Error reading emfile \"" << filename_ << "\": " << tom_io_strerror(i, errno_, cbuf, sizeof(cbuf));
            throw std::runtime_error(ss.str());
        }

        if (!v.is_equal_size(dims[0], dims[1], dims[2])) {
            ss << "Error: reading emfile \"" << filename_ << "\" would result in a volume size of " << dims[0] << 'x' << dims[1] << 'x' << dims[2] << " but a volume of size " << v.getSizeX() << 'x' << v.getSizeY() << 'x' << v.getSizeZ() << " was given.";
            throw std::runtime_error(ss.str());
        }

        if ((i=tom_io_read_vol(f_, header_.dims, tom_io_em_get_iotype_header(&header_), tom_io_em_is_swaped(&header_), subregion_, sampling_, binning_, &v.get(), tom::io::get_tom_io_type<T>())) != TOM_ERR_OK) {
            const int errno_ = errno;
            char cbuf[1024];
            ss << "Error reading emfile \"" << filename_ << "\": " << tom_io_strerror(i, errno_, cbuf, sizeof(cbuf));
            throw std::runtime_error(ss.str());
        }

    } catch (...) {
        std::fclose(f_);
        f_is_closed_ = true;
        seeked_at_data_start_ = false;
        throw;
    }

    if (close_after_first_read_) {
        // Close the file.
        std::fclose(f_);
        f_is_closed_ = true;
    }
    seeked_at_data_start_ = false;
}




/****************************************************************************//**
 * \brief Reads a volume from em-file.
 *
 * \param[out] v A reference to a pointer where the volume is returned.
 * \param[in] filename The filename of the emfile.
 * \param[in] subregion The subregion to read. See tom_io_em_read.
 * \param[in] sampling The sampling factor. See tom_io_em_read.
 * \param[in] binning The binning factor. See tom_io_em_read.
 * \param[out] header A pointer to the em-header. Set to \c NULL if you are not
 *   interested in the header.
 * \param[in] fcn_malloc A pointer to a function with which the memory will be
 *   allocated. Setting to \c NULL defaults to tom::getFcnMalloc(). If
 *   getFcnMalloc returns \c NULL too, the new[] for type T is used.
 * \param[in] fcn_free A pointer to a function with which the memory allocated with
 *   fcn_malloc can be freed again. Setting to \c NULL defaults to
 *   tom::getFcnFree(), and here too, to delete[].
 *
 * Calls tom_io_em_read to read the file, and throws the integer error status
 * returned by it, in case of an error.
 *
 * \TODO there is still a time consuming assert, to check that the new implementation
 *  which ought to replace tom::read_from_em yields the same results.
 *******************************************************************************/
template<typename T>
void tom::io::read_from_em(tom::Volume<T> *&v, const std::string &filename, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning, tom_io_em_header *header, void *(*fcn_malloc)(size_t size), void (*fcn_free)(void *ptr)) {

    std::auto_ptr<tom::Volume<T> > v_ = tom::io::VolumeEM(filename, fcn_malloc, fcn_free).read_volume<T>(subregion, sampling, binning);

    #ifndef NDEBUG

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


    void *vdata = NULL;
    uint32_t dims[3];
    int restype;
    std::auto_ptr<tom::Volume<T> > v_local;

    try {
        int res = ::tom_io_em_read(filename.c_str(), subregion, sampling, binning, header, &vdata, dims, &restype, fcn_malloc, fcn_free);

        if (res != TOM_ERR_OK) {
            const int errno_ = errno;
            std::stringstream ss;
            std::vector<char> buf(1024);
            ss << "Error reading em-file \"" << filename << "\": " << tom_io_strerror(res, errno_, &buf[0], buf.size());
            throw std::runtime_error(ss.str());
        }

        if (restype == tom::get_tom_io_type<T>()) {
            v_local.reset(new tom::Volume<T>((T *)vdata, dims[0], dims[1], dims[2], 0, 0, 0, true, fcn_free));
            vdata = NULL;
        } else {
            v_local.reset(new tom::Volume<T>(dims[0], dims[1], dims[2], fcn_malloc, fcn_free));
            switch (restype) {
                case TOM_IO_TYPE_FLOAT: {
                        tom::Volume<float> v2((float *)vdata, dims[0], dims[1], dims[2], 0, 0, 0, true, fcn_free);
                        vdata = NULL; /* The data from vdata is now owned by v2. Leave clean up to local object. */
                        if (tom::is_double<T>()) {
                            tom::Volume<double> *v_local_typed = reinterpret_cast<tom::Volume<double> *>(v_local.get());
                            v_local_typed->setValues(v2);
                        } 
						else {
                            throw std::runtime_error("Reading EM file: problem with type conversion from float your volume type!");
                        }
					
                    }
                    break;
                case TOM_IO_TYPE_COMPLEX32: {
                        tom::Volume<std::complex<float> > v2((std::complex<float> *)vdata, dims[0], dims[1], dims[2], 0, 0, 0, true, fcn_free);
                        vdata = NULL; /* The data from vdata is now owned by v2. Leave it to him to clean up. */
                        if (tom::is_double_complex<T>()) {
                            tom::Volume<std::complex<double> > *v_local_typed = reinterpret_cast<tom::Volume<std::complex<double> > *>(v_local.get());
                            v_local_typed->setValues(v2);
                        } else {
                            throw std::runtime_error("Reading EM file: problem with type conversion from complex float your volume type!");
                        }
                    }
                    break;
                case TOM_IO_TYPE_DOUBLE: {
                        tom::Volume<double> v2((double *)vdata, dims[0], dims[1], dims[2], 0, 0, 0, true, fcn_free);
                        vdata = NULL; /* The data from vdata is now owned by v2. Leave it to him to clean up. */
                        if (tom::is_float<T>()) {
                            tom::Volume<float> *v_local_typed = reinterpret_cast<tom::Volume<float> *>(v_local.get());
                            v_local_typed->setValues(v2);
                        } else {
                            throw std::runtime_error("Reading EM file: problem with type conversion from double your volume type!");
                        }
                    }
                    break;
                case TOM_IO_TYPE_COMPLEX64: {
                        tom::Volume<std::complex<double> > v2((std::complex<double> *)vdata, dims[0], dims[1], dims[2], 0, 0, 0, true, fcn_free);
                        vdata = NULL; /* The data from vdata is now owned by v2. Leave it to him to clean up. */
                        if (tom::is_float_complex<T>()) {
                            tom::Volume<std::complex<float> > *v_local_typed = reinterpret_cast<tom::Volume<std::complex<float> > *>(v_local.get());
                            v_local_typed->setValues(v2);
                        } else {
                            throw std::runtime_error("Reading EM file: problem with type conversion from complex double your volume type!");
                        }
                    }
                    break;
                case TOM_IO_TYPE_INT8: {
						tom::Volume<char> v2((char *)vdata, dims[0], dims[1], dims[2], 0, 0, 0, true, fcn_free);
						vdata = NULL; /* The data from vdata is now owned by v2. Leave clean up to local object. */
						if (tom::is_float<T>()) {
							tom::Volume<float> *v_local_typed = reinterpret_cast<tom::Volume<float> *>(v_local.get());
							v_local_typed->setValues(v2);
						} else {
							throw std::runtime_error("Reading EM file: problem with type conversion from int8 your volume type!");
						}
					}
					break;
                case TOM_IO_TYPE_INT16: {
						tom::Volume<int16_t> v2((int16_t *)vdata, dims[0], dims[1], dims[2], 0, 0, 0, true, fcn_free);
						vdata = NULL; /* The data from vdata is now owned by v2. Leave clean up to local object. */
						if (tom::is_float<T>()) {
								tom::Volume<float> *v_local_typed = reinterpret_cast<tom::Volume<float> *>(v_local.get());
								v_local_typed->setValues(v2);
						} else {
							throw std::runtime_error("Reading EM file: problem with type conversion from int16 your volume type!");
						}
					}
					break;
                case TOM_IO_TYPE_INT32:
					{
						//type of file is int32! generate a int32 copy of the data
						tom::Volume<int32_t> v2((int32_t *)vdata, dims[0], dims[1], dims[2], 0, 0, 0, true, fcn_free);
						vdata = NULL; /* The data from vdata is now owned by v2. Leave clean up to local object. */
						if (tom::is_float<T>()) {
							tom::Volume<float> *v_local_typed = reinterpret_cast<tom::Volume<float> *>(v_local.get());
							//copy from v2 
							v_local_typed->setValues(v2);
						} else {
							throw std::runtime_error("Reading EM file: problem with type conversion from int32 your volume type!");
						}
					}
					break;
				default:
                    throw std::logic_error("Reading EM file: the datatype in the file is not supported!");
            }
        }
    } catch (...) {
        /* Clean up... */
        if  (vdata) { fcn_free(vdata); }
        throw;
    }
    assert(*v_ == *v_local);
    #endif

    v = v_.release();
}



/***************************************************************************//**
 * Gets the resulting size of the volume for a certain \c subregion, \c sampling
 * and \c binning. If the parameters are out of range (e.g. the subregion), an
 * error number is returned and the size parameters are left unchanged.
 * Otherwise TOM_ERR_OK is returned.
 ******************************************************************************/
int tom::io::VolumeEM::get_size(const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning, std::size_t &sizex, std::size_t &sizey, std::size_t &sizez) const {
    int i;
    uint32_t dims[3];

    /* Get the subregion and check the dimension. */
    if ((i=tom_io_calculate_sizes(header_.dims, subregion, 0, sampling, 0, binning, 0, 0, dims)) != TOM_ERR_OK) {
        return i;
    }
    sizex = dims[0];
    sizey = dims[1];
    sizez = dims[2];
    return TOM_ERR_OK;
}





/****************************************************************************//**
 * \brief Writes the volume as em-file.
 *
 * \param[in] v Input volume.
 * \param[in] filename The name of the em-file.
 * \param[in] header A pointer to the used em-header. If set to NULL, a default
 *   header is used with only the needed fields set according to the calling
 *   object. Even if a header is given, the relevant fields (\c dims and \c type)
 *   are set to match the data from the volume.
 *
 * Existing files are overwritten. The function can perform a type conversion
 * from the datatype of the volume (\c TVOL) to the datatype in em-file
 * (\c TFILE)
 *
 * Calles tom_io_em_write to write the data. If an error occures
 * there, its integer error status is thrown as exception.
 *******************************************************************************/
template<typename TVOL, typename TFILE>
void tom::io::write_to_em(const tom::Volume<TVOL> &v, const std::string &filename, const tom_io_em_header *header) {

    if (v.getSizeX() > std::numeric_limits<uint32_t>::max() ||
        v.getSizeY() > std::numeric_limits<uint32_t>::max() ||
        v.getSizeZ() > std::numeric_limits<uint32_t>::max()) {
        throw std::runtime_error("tom::io::write_to_em: The volume is too large to save it as em-format.");
    }

    tom_io_em_header header_local;
    if (header) {
        header_local = *header;
    } else {
        memset(&header_local, 0, sizeof(header_local));
        header_local.machine = 6;
    }

    assert(!"UNSOUPPORTED SRC-TYPE" || tom::get_tom_io_type<TVOL >() != TOM_IO_TYPE_UNKNOWN );
    assert(!"UNSOUPPORTED DST-TYPE" || tom::get_tom_io_type<TFILE>() != TOM_IO_TYPE_UNKNOWN );

    if (!tom_io_em_set_iotype_header(&header_local, tom::get_tom_io_type<TFILE>())) {
        assert(!"Can not save the volume as this type. NOT SUPPORTED (either by io.c or by EM-format.");
    }
    header_local.dims[0] = v.getSizeX();
    header_local.dims[1] = v.getSizeY();
    header_local.dims[2] = v.getSizeZ();

    const int res  = tom_io_em_write(filename.c_str(), &header_local, &v.get(), tom::get_tom_io_type<TVOL>(),
                                     v.getStrideX(), v.getStrideY(), v.getStrideZ());

    if (res != TOM_ERR_OK) {
        const int errno_ = errno;
        std::stringstream ss;
        std::vector<char> buf(1024);
        ss << "Error writing em-file \"" << filename << "\". " << tom_io_strerror(res, errno_, &buf[0], buf.size());
        throw std::runtime_error(ss.str());
    }
}



/*
template<typename TVOL>
void tom::io::write_to_em_paste(const tom::Volume<TVOL> &v, const std::string &filename, const tom_io_em_header *header, const TVOL *first_voxel) {

    assert(!" tom::io::write_to_em_paste THIS FUNCTION WAS NOT YET NEEDED, NOR TESTED. CHECK IT NOW AS YOU GET THIS EXCEPTION");

    if (v.getSizeX() > std::numeric_limits<uint32_t>::max() ||
        v.getSizeY() > std::numeric_limits<uint32_t>::max() ||
        v.getSizeZ() > std::numeric_limits<uint32_t>::max()) {
        throw std::runtime_error("tom::io::write_to_em_paste: The volume is too large to save it as em-format.");
    }

    uint32_t first_voxelv[3];
    if (first_voxel) {
        if (first_voxel < &v.get()) {
            throw std::invalid_argument("first voxel is not inside the volume");
        }
        std::size_t idxx, idxy, idxz;
        std::ptrdiff_t offset = first_voxel - &v.get();

        std::size_t size_1 = v.getSizeX()*v.getStrideX();
        std::size_t size_2 = v.getSizeY()*v.getStrideY() * size_1;

        idxz = offset / size_2;
        idxy = (offset % size_2) / size_1;
        idxx = ((offset % size_2) % size_1) / v.getStrideX();
        if (((offset % size_2) % size_1) % v.getStrideX() != 0 ||
            idxx >= v.getSizeX() ||
            idxy >= v.getSizeY() ||
            idxz >= v.getSizeZ()) {
            throw std::invalid_argument("first voxel is not inside the volume");
        }
        first_voxelv[0] = idxx;
        first_voxelv[1] = idxy;
        first_voxelv[2] = idxz;
    } else {
        first_voxelv[0] =
        first_voxelv[1] =
        first_voxelv[2] = 0;
    }



    tom_io_em_header header_local;
    if (header) {
        header_local = *header;
    } else {
        memset(&header_local, 0, sizeof(header_local));
        header_local.machine = 6;
    }

    assert(!"UNSOUPPORTED TYPE" || tom::get_tom_io_type<TVOL >() != TOM_IO_TYPE_UNKNOWN );
    if (!tom_io_em_set_iotype_header(&header_local, tom::get_tom_io_type<TVOL>())) {
        assert(!"Can not save the volume as this type. NOT SUPPORTED (either by io.c or by EM-format.");
    }
    header_local.dims[0] = v.getSizeX();
    header_local.dims[1] = v.getSizeY();
    header_local.dims[2] = v.getSizeZ();

    int header_read;


    int res  = tom_io_em_write_paste(filename.c_str(), &v.get(), tom::get_tom_io_type<TVOL>(),
                                    v.getSizeX(), v.getSizeY(), v.getSizeZ(),
                                    v.getStrideX(), v.getStrideY(),v.getStrideZ(),
                                    &header_local, &header_read, true, first_voxelv, 0);
    if (res != TOM_ERR_OK) {
        const int errno_ = errno;
        std::stringstream ss;
        std::vector<char> buf(1024);
        ss << "Error pasting to em-file \"" << filename << "\". " << tom_io_strerror(res, errno_, &buf[0], buf.size());
        throw std::runtime_error(ss.str());
    }
}
*/

/**
 * @brief Write a volume as pasting into an existing file.
 *
 * @param[in] v Input volume.
 * @param[in] filename The name of the file.
 * @param[in] first_voxel The index of the first voxel which the volume should start from.
 */
template<typename TVOL>
void tom::io::write_to_em_paste(const tom::Volume<TVOL> &v, const std::string &filename, const uint32_t *first_voxel) {

    if (v.getSizeX() > std::numeric_limits<uint32_t>::max() ||
        v.getSizeY() > std::numeric_limits<uint32_t>::max() ||
        v.getSizeZ() > std::numeric_limits<uint32_t>::max()) {
        throw std::runtime_error("tom::io::write_to_em_paste: The volume is too large to save it as em-format.");
    }

    int res  = tom_io_em_write_paste(filename.c_str(), &v.get(), tom::get_tom_io_type<TVOL>(),
                                      v.getSizeX(), v.getSizeY(), v.getSizeZ(),
                                      v.getStrideX(), v.getStrideY(),v.getStrideZ(),
                                      NULL, NULL, true, first_voxel, 0);
    if (res != TOM_ERR_OK) {
        const int errno_ = errno;
        std::stringstream ss;
        std::vector<char> buf(1024);
        ss << "Error pasting to em-file \"" << filename << "\". " << tom_io_strerror(res, errno_, &buf[0], buf.size());
        throw std::runtime_error(ss.str());
    }
}

/**
 * @brief Read the MRC header from the file
 *
 * @param[in] filename Name of the file
 * @param[out] header MRC header
 */
void tom::io::mrc_read_header(const std::string &filename, tom_io_mrc_header &header) {

    const int tom_errnum = tom_io_mrc_read_header(filename.c_str(), 0, &header, 0);

    if (tom_errnum != TOM_ERR_OK) {
        const int errno_ = errno;
        std::ostringstream ss;
        char cbuf[2048];
        ss << "Error reading header of mrc file \"" << filename << "\": " << tom_io_strerror(tom_errnum, errno_, cbuf, sizeof(cbuf));
        throw std::runtime_error(ss.str());
    }
}

/**
 * @brief Print some information of the MRC header
 */
void tom::io::mrc_print_header_info(const tom_io_mrc_header &header)
{
	std::cout << "number of columns: " << header.nx << std::endl;
	std::cout << "number of rows: " << header.ny << std::endl;
	std::cout << "number of sections: " << header.nz << std::endl;
	std::cout << "data type: " << header.mode << std::endl;
	std::cout << "number of first column in map: " << header.nxstart << std::endl;
	std::cout << "number of first row in map: " << header.nystart << std::endl;
	std::cout << "number of first section in map: " << header.nzstart << std::endl;
	std::cout << "number of intervals along X: " << header.mx << std::endl;
	std::cout << "number of intervals along Y: " << header.my << std::endl;
	std::cout << "number of intervals along Z: " << header.mz << std::endl;
//	std::cout << "cell dimensions in angstroms: " << header.cella[0] << header.cella[1] << header.cella[2] << std::endl;
//	std::cout << "cell angles in degrees: " << header.cellb[0] << header.cellb[1] << header.cellb[2] << std::endl;
	std::cout << "axis corresp to cols: " << header.mapc << std::endl;
	std::cout << "axis corresp to rows: " << header.mapr << std::endl;
	std::cout << "axis corresp to sections: " << header.maps << std::endl;
	std::cout << "minimum density value: " << header.dmin << std::endl;
	std::cout << "maximum density value: " << header.dmax << std::endl;
	std::cout << "mean density value: " << header.dmean << std::endl;
	std::cout << "space group number: " << header.ispg << std::endl;
	std::cout << "number of bytes used for symmetry data: " << header.nsymbt << std::endl;
//	std::cout << "extra space used for anything: " << header.nx << std::endl;
//	std::cout << "origin in X,Y,Z used for transforms: " << header.origin[0] << header.origin[1] << header.origin[2] << std::endl;
//	std::cout << "character string 'MAP' to identify file type: " << header.map << std::endl;
	std::cout << "machine stamp: " << header.machst << std::endl;
	std::cout << "rms deviation of map from mean density: " << header.rms << std::endl;
	std::cout << "number of labels being used: " << header.nlabl << std::endl;
}

/**
 * @brief Reads a volume as well as the header from mrc-file.
 *
 * @param[out] v A reference to a pointer where the volume is returned.
 * @param[in] filename The filename of the file.
 * @param[in] subregion The subregion to read.
 * @param[in] sampling The sampling factor.
 * @param[in] binning The binning factor.
 * @param[out] header A pointer to the mrc-header. Set to \c NULL if you are not
 *   interested in the header.
 * @param[in] fcn_malloc A pointer to a function with which the memory will be
 *   allocated. Setting to \c NULL defaults to tom::getFcnMalloc(). If
 *   getFcnMalloc returns \c NULL too, the new[] for type T is used.
 * @param[in] fcn_free A pointer to a function with which the memory allocated with
 *   fcn_malloc can be freed again. Setting to \c NULL defaults to
 *   tom::getFcnFree(), and here too, to delete[].
 *
 * Calls tom_io_mrc_read to read the file, and throws the integer error status
 * returned by it, in case of an error.
 */
template<typename T>
void tom::io::read_from_mrc(tom::Volume<T> *&v, const std::string &filename, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning, tom_io_mrc_header *header, void *(*fcn_malloc)(size_t size), void (*fcn_free)(void *ptr))
{

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


    void *vdata = NULL;
    uint32_t dims[3];
    int restype;
    std::auto_ptr<tom::Volume<T> > v_local;

    try {
        int res = ::tom_io_mrc_read(filename.c_str(), subregion, sampling, binning, header, &vdata, dims, &restype, fcn_malloc, fcn_free);

        if (res != TOM_ERR_OK) {
            const int errno_ = errno;
            std::stringstream ss;
            std::vector<char> buf(1024);
            ss << "Error reading mrc-file \"" << filename << "\": " << tom_io_strerror(res, errno_, &buf[0], buf.size());
            throw std::runtime_error(ss.str());
        }

        if (restype == tom::get_tom_io_type<T>()) {
            v_local.reset(new tom::Volume<T>((T *)vdata, dims[0], dims[1], dims[2], 0, 0, 0, true, fcn_free));
            vdata = NULL;
        } else { // Type casting is only allowed for float and double
            v_local.reset(new tom::Volume<T>(dims[0], dims[1], dims[2], fcn_malloc, fcn_free));
            switch (restype) {
                case TOM_IO_TYPE_FLOAT: {
                        tom::Volume<float> v2((float *)vdata, dims[0], dims[1], dims[2], 0, 0, 0, true, fcn_free);
                        vdata = NULL; /* The data from vdata is now owned by v2. Leave clean up to local object. */
                        if (tom::is_double<T>()) {
                            tom::Volume<double> *v_local_typed = reinterpret_cast<tom::Volume<double> *>(v_local.get());
                            v_local_typed->setValues(v2);
                        } else {
                            throw std::runtime_error("Reading MRC file: problem with type conversion from float your volume type!");
                        }
                    }
                    break;
                case TOM_IO_TYPE_DOUBLE: {
                        tom::Volume<double> v2((double *)vdata, dims[0], dims[1], dims[2], 0, 0, 0, true, fcn_free);
                        vdata = NULL; /* The data from vdata is now owned by v2. Leave it to him to clean up. */
                        if (tom::is_float<T>()) {
                            tom::Volume<float> *v_local_typed = reinterpret_cast<tom::Volume<float> *>(v_local.get());
                            v_local_typed->setValues(v2);
                        } else {
                            throw std::runtime_error("Reading MRC file: problem with type conversion from double your volume type!");
                        }
                    }
                    break;
                case TOM_IO_TYPE_INT8: {
						tom::Volume<char> v2((char *)vdata, dims[0], dims[1], dims[2], 0, 0, 0, true, fcn_free);
						vdata = NULL;
                        if (tom::is_float<T>()) {
                            tom::Volume<float> *v_local_typed = reinterpret_cast<tom::Volume<float> *>(v_local.get());
                            v_local_typed->setValues(v2);
                        } else {
                            throw std::runtime_error("Reading MRC file: problem with type conversion from int8 your volume type!");
                        }
					}
					break;
                case TOM_IO_TYPE_INT16: {
						tom::Volume<int16_t> v2((int16_t *)vdata, dims[0], dims[1], dims[2], 0, 0, 0, true, fcn_free);
						vdata = NULL; /* The data from vdata is now owned by v2. Leave clean up to local object. */
						if (tom::is_float<T>()) {
								tom::Volume<float> *v_local_typed = reinterpret_cast<tom::Volume<float> *>(v_local.get());
								v_local_typed->setValues(v2);
						} else {
							throw std::runtime_error("Reading MRC file: problem with type conversion from int16 your volume type!");
						}
					}
					break;
				case TOM_IO_TYPE_UINT16: {
						tom::Volume<uint16_t> v2((uint16_t *)vdata, dims[0], dims[1], dims[2], 0, 0, 0, true, fcn_free);
						vdata = NULL; /* The data from vdata is now owned by v2. Leave clean up to local object. */
						if (tom::is_float<T>()) {
								tom::Volume<float> *v_local_typed = reinterpret_cast<tom::Volume<float> *>(v_local.get());
								v_local_typed->setValues(v2);
						} else {
							throw std::runtime_error("Reading MRC file: problem with type conversion from uint16 your volume type!");
						}
					}
					break;
				case TOM_IO_TYPE_INT32:
					{
						//type of file is int32! generate a int32 copy of the data
						tom::Volume<int32_t> v2((int32_t *)vdata, dims[0], dims[1], dims[2], 0, 0, 0, true, fcn_free);
						vdata = NULL; /* The data from vdata is now owned by v2. Leave clean up to local object. */
						if (tom::is_float<T>()) {
							tom::Volume<float> *v_local_typed = reinterpret_cast<tom::Volume<float> *>(v_local.get());
							//copy from v2
							v_local_typed->setValues(v2);
						} else {
							throw std::runtime_error("Reading MRC file: problem with type conversion from int32 your volume type!");
						}
					}
					break;
                default:
                    throw std::logic_error("Reading MRC file: the datatype in the file is not supported!");
            }
        }
    } catch (...) {
        /* Clean up... */
        if  (vdata) { fcn_free(vdata); }
        throw;
    }

    v = v_local.release();
}
/*
template<>
void tom::io::write_to_mrc(const tom::Volume<std::complex<float> > &v, const std::string &filename, const tom_io_mrc_header *header)
{
	throw std::runtime_error("tom::io::write_to_mrc does not support complex numbers!");
}
*/

/**
 * @brief Writes the volume as mrc-file.
 *
 * @param[in] v Input volume.
 * @param[in] filename The name of the file.
 * @param[in] header A pointer to the used mrc-header. If set to NULL, a default
 *   header is used with only the needed fields set according to the calling
 *   object. Even if a header is given, the relevant fields (\c dims and \c type)
 *   are set to match the data from the volume.
 *
 * Existing files are overwritten. The function can perform a type conversion
 * from the datatype of the volume (\c TVOL) to the datatype in mrc-file
 * (\c TFILE)
 *
 * Calles tom_io_mrc_write to write the data. If an error occures
 * there, its integer error status is thrown as exception.
 *
 * Note: Swapping is not considered in this case.
 * And UINT_16 and complex 16-bit integers are not supported yet!
 */

template<typename TVOL, typename TFILE>
void tom::io::write_to_mrc(const tom::Volume<TVOL> &v, const std::string &filename, const tom_io_mrc_header *header) {

    if (v.getSizeX() > std::numeric_limits<uint32_t>::max() ||
        v.getSizeY() > std::numeric_limits<uint32_t>::max() ||
        v.getSizeZ() > std::numeric_limits<uint32_t>::max()) {
        throw std::runtime_error("tom::io::write_to_mrc: The volume is too large to save it as mrc-format.");
    }

    tom_io_mrc_header header_local;
    if (header) {
        header_local = *header;
    } else {
        memset(&header_local, 0, sizeof(header_local));
    }

    assert(!"UNSOUPPORTED SRC-TYPE" || tom::get_tom_io_type<TVOL >() != TOM_IO_TYPE_UNKNOWN );
    assert(!"UNSOUPPORTED DST-TYPE" || tom::get_tom_io_type<TFILE>() != TOM_IO_TYPE_UNKNOWN );

    // set the data format and other stuff in the header file
    switch(tom::get_tom_io_type<TFILE>())
    {
		case TOM_IO_TYPE_INT8:
			header_local.mode = 0;
			break;
		case TOM_IO_TYPE_INT16:
			header_local.mode = 1;
			break;
		case TOM_IO_TYPE_FLOAT:
			header_local.mode = 2;
			break;
		case TOM_IO_TYPE_COMPLEX32:
			header_local.mode = 4;
			break;
		default:
			std::string errorString("This data type is not supported for MRC! Supported ones are TOM_IO_TYPE_INT8 TOM_IO_TYPE_INT16 TOM_IO_TYPE_FLOAT TOM_IO_TYPE_COMPLEX32.");
			throw std::runtime_error(errorString);
    }
    header_local.nx = v.getSizeX();
    header_local.ny = v.getSizeY();
    header_local.nz = v.getSizeZ();

    const int res  = tom_io_mrc_write(filename.c_str(), &header_local, &v.get(), tom::get_tom_io_type<TVOL>(),
                                     v.getStrideX(), v.getStrideY(), v.getStrideZ());

    if (res != TOM_ERR_OK) {
        const int errno_ = errno;
        std::stringstream ss;
        std::vector<char> buf(1024);
        ss << "Error writing mrc-file \"" << filename << "\". " << tom_io_strerror(res, errno_, &buf[0], buf.size());
        throw std::runtime_error(ss.str());
    }
}

/**
 * @brief Write a volume as pasting into an existing file.
 *
 * @param[in] v Input volume.
 * @param[in] filename The name of the file.
 * @param[in] first_voxel The index of the first voxel which the volume should start from.
 */
template<typename TVOL>
void tom::io::write_to_mrc_paste(const tom::Volume<TVOL> &v, const std::string &filename, const uint32_t *first_voxel) {

    if (v.getSizeX() > std::numeric_limits<uint32_t>::max() ||
        v.getSizeY() > std::numeric_limits<uint32_t>::max() ||
        v.getSizeZ() > std::numeric_limits<uint32_t>::max()) {
        throw std::runtime_error("tom::io::write_to_mrc_paste: The volume is too large to save it as mrc-format.");
    }

    int res  = tom_io_mrc_write_paste(filename.c_str(), &v.get(), tom::get_tom_io_type<TVOL>(),
                                      v.getSizeX(), v.getSizeY(), v.getSizeZ(),
                                      v.getStrideX(), v.getStrideY(),v.getStrideZ(),
                                      NULL, NULL, true, first_voxel, 0);
    if (res != TOM_ERR_OK) {
        const int errno_ = errno;
        std::stringstream ss;
        std::vector<char> buf(1024);
        ss << "Error pasting to mrc-file \"" << filename << "\". " << tom_io_strerror(res, errno_, &buf[0], buf.size());
        throw std::runtime_error(ss.str());
    }
}

/**
 * @brief Read the CCP4 header from the file
 *
 * @param[in] filename Name of the file
 * @param[out] header CCP4 header
 */
void tom::io::ccp4_read_header(const std::string &filename, tom_io_ccp4_header &header) {

    const int tom_errnum = tom_io_ccp4_read_header(filename.c_str(), 0, &header, 0);

    if (tom_errnum != TOM_ERR_OK) {
        const int errno_ = errno;
        std::ostringstream ss;
        char cbuf[2048];
        ss << "Error reading header of ccp4 file \"" << filename << "\": " << tom_io_strerror(tom_errnum, errno_, cbuf, sizeof(cbuf));
        throw std::runtime_error(ss.str());
    }
}

/**
 * @brief Print some information of the CCP4 header
 */
void tom::io::ccp4_print_header_info(const tom_io_ccp4_header &header)
{
	std::cout << "number of columns: " << header.nc << std::endl;
	std::cout << "number of rows: " << header.nr << std::endl;
	std::cout << "number of sections: " << header.ns << std::endl;
	std::cout << "data type: " << header.mode << std::endl;
	std::cout << "number of first column in map: " << header.ncstart << std::endl;
	std::cout << "number of first row in map: " << header.nrstart << std::endl;
	std::cout << "number of first section in map: " << header.nsstart << std::endl;
	std::cout << "number of intervals along X: " << header.nx << std::endl;
	std::cout << "number of intervals along Y: " << header.ny << std::endl;
	std::cout << "number of intervals along Z: " << header.nz << std::endl;
	std::cout << "cell dimensions in angstroms: " << header.x_length << " "<< header.y_length << " " << header.z_length << std::endl;
	std::cout << "cell angles in degrees: " << header.alpha << " " << header.beta << " " << header.gamma << std::endl;
	std::cout << "axis corresp to cols: " << header.mapc << std::endl;
	std::cout << "axis corresp to rows: " << header.mapr << std::endl;
	std::cout << "axis corresp to sections: " << header.maps << std::endl;
	std::cout << "minimum density value: " << header.amin << std::endl;
	std::cout << "maximum density value: " << header.amax << std::endl;
	std::cout << "mean density value: " << header.amean << std::endl;
	std::cout << "space group number: " << header.ispg << std::endl;
	std::cout << "number of bytes used for storing symmetry operators: " << header.nsympg << std::endl;
	std::cout << "flag for skew transformation: " << header.lskflg << std::endl;
	std::cout << "machine stamp: " << header.machst << std::endl;
	std::cout << "rms deviation of map from mean density: " << header.arms << std::endl;
	std::cout << "number of labels being used: " << header.nlabl << std::endl;
}

/**
 * @brief Reads a volume as well as the header from ccp4-file.
 *
 * @param[out] v A reference to a pointer where the volume is returned.
 * @param[in] filename The filename of the file.
 * @param[in] subregion The subregion to read.
 * @param[in] sampling The sampling factor.
 * @param[in] binning The binning factor.
 * @param[out] header A pointer to the ccp4-header. Set to \c NULL if you are not
 *   interested in the header.
 * @param[in] fcn_malloc A pointer to a function with which the memory will be
 *   allocated. Setting to \c NULL defaults to tom::getFcnMalloc(). If
 *   getFcnMalloc returns \c NULL too, the new[] for type T is used.
 * @param[in] fcn_free A pointer to a function with which the memory allocated with
 *   fcn_malloc can be freed again. Setting to \c NULL defaults to
 *   tom::getFcnFree(), and here too, to delete[].
 *
 * Calls tom_io_ccp4_read to read the file, and throws the integer error status
 * returned by it, in case of an error.
 */
template<typename T>
void tom::io::read_from_ccp4(tom::Volume<T> *&v, const std::string &filename, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning, tom_io_ccp4_header *header, void *(*fcn_malloc)(size_t size), void (*fcn_free)(void *ptr))
{

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


    void *vdata = NULL;
    uint32_t dims[3];
    int restype;
    std::auto_ptr<tom::Volume<T> > v_local;

    try {
        int res = ::tom_io_ccp4_read(filename.c_str(), subregion, sampling, binning, header, &vdata, dims, &restype, fcn_malloc, fcn_free);

        if (res != TOM_ERR_OK) {
            const int errno_ = errno;
            std::stringstream ss;
            std::vector<char> buf(1024);
            ss << "Error reading ccp4-file \"" << filename << "\": " << tom_io_strerror(res, errno_, &buf[0], buf.size());
            throw std::runtime_error(ss.str());
        }

        if (restype == tom::get_tom_io_type<T>()) {
            v_local.reset(new tom::Volume<T>((T *)vdata, dims[0], dims[1], dims[2], 0, 0, 0, true, fcn_free));
            vdata = NULL;
        } else { // Type casting is only allowed for float and double
            v_local.reset(new tom::Volume<T>(dims[0], dims[1], dims[2], fcn_malloc, fcn_free));
            switch (restype) {
                case TOM_IO_TYPE_FLOAT: {
                        tom::Volume<float> v2((float *)vdata, dims[0], dims[1], dims[2], 0, 0, 0, true, fcn_free);
                        vdata = NULL; /* The data from vdata is now owned by v2. Leave clean up to local object. */
                        if (tom::is_double<T>()) {
                            tom::Volume<double> *v_local_typed = reinterpret_cast<tom::Volume<double> *>(v_local.get());
                            v_local_typed->setValues(v2);
                        } else {
                            throw std::runtime_error("Reading CCP4 file: problem with type conversion from float your volume type!");
                        }
                    }
                    break;
                case TOM_IO_TYPE_DOUBLE: {
                        tom::Volume<double> v2((double *)vdata, dims[0], dims[1], dims[2], 0, 0, 0, true, fcn_free);
                        vdata = NULL; /* The data from vdata is now owned by v2. Leave it to him to clean up. */
                        if (tom::is_float<T>()) {
                            tom::Volume<float> *v_local_typed = reinterpret_cast<tom::Volume<float> *>(v_local.get());
                            v_local_typed->setValues(v2);
                        } else {
                            throw std::runtime_error("Reading CCP4 file: problem with type conversion from double your volume type!");
                        }
                    }
                    break;
                case TOM_IO_TYPE_INT8:
                case TOM_IO_TYPE_INT16:
                case TOM_IO_TYPE_INT32:
                default:
                    throw std::logic_error("Reading CCP4 file: the datatype in the file is not supported!");
            }
        }
    } catch (...) {
        /* Clean up... */
        if  (vdata) { fcn_free(vdata); }
        throw;
    }

    v = v_local.release();
}

/**
 * @brief Writes the volume as ccp4-file.
 *
 * @param[in] v Input volume.
 * @param[in] filename The name of the file.
 * @param[in] header A pointer to the used ccp4-header. If set to NULL, a default
 *   header is used with only the needed fields set according to the calling
 *   object. Even if a header is given, the relevant fields (\c dims and \c type)
 *   are set to match the data from the volume.
 *
 * Existing files are overwritten. The function can perform a type conversion
 * from the datatype of the volume (\c TVOL) to the datatype in mrc-file
 * (\c TFILE)
 *
 * Calles tom_io_ccp4_write to write the data. If an error occures
 * there, its integer error status is thrown as exception.
 *
 * Note: Swapping is not considered in this case.
 */
template<typename TVOL, typename TFILE>
void tom::io::write_to_ccp4(const tom::Volume<TVOL> &v, const std::string &filename, const tom_io_ccp4_header *header) {

    if (v.getSizeX() > std::numeric_limits<uint32_t>::max() ||
        v.getSizeY() > std::numeric_limits<uint32_t>::max() ||
        v.getSizeZ() > std::numeric_limits<uint32_t>::max()) {
        throw std::runtime_error("tom::io::write_to_ccp4: The volume is too large to save it as ccp4-format.");
    }

    tom_io_ccp4_header header_local;
    if (header) {
        header_local = *header;
    } else {
        memset(&header_local, 0, sizeof(header_local));
    }

    assert(!"UNSOUPPORTED SRC-TYPE" || tom::get_tom_io_type<TVOL >() != TOM_IO_TYPE_UNKNOWN );
    assert(!"UNSOUPPORTED DST-TYPE" || tom::get_tom_io_type<TFILE>() != TOM_IO_TYPE_UNKNOWN );

    // set the data format and other stuff in the header file
    switch(tom::get_tom_io_type<TFILE>())
    {
		case TOM_IO_TYPE_INT8:
			header_local.mode = 0;
			break;
		case TOM_IO_TYPE_INT16:
			header_local.mode = 1;
			break;
		case TOM_IO_TYPE_FLOAT:
			header_local.mode = 2;
			break;
		case TOM_IO_TYPE_COMPLEX32:
			header_local.mode = 4;
			break;
		default:
			std::string errorString("This data type is not supported for CCP4! Supported ones are TOM_IO_TYPE_INT8 TOM_IO_TYPE_INT16 TOM_IO_TYPE_FLOAT TOM_IO_TYPE_COMPLEX32.");
			throw std::runtime_error(errorString);
    }
    header_local.nc = v.getSizeX();
    header_local.nr = v.getSizeY();
    header_local.ns = v.getSizeZ();

    const int res  = tom_io_ccp4_write(filename.c_str(), &header_local, &v.get(), tom::get_tom_io_type<TVOL>(),
                                     v.getStrideX(), v.getStrideY(), v.getStrideZ());

    if (res != TOM_ERR_OK) {
        const int errno_ = errno;
        std::stringstream ss;
        std::vector<char> buf(1024);
        ss << "Error writing ccp4-file \"" << filename << "\". " << tom_io_strerror(res, errno_, &buf[0], buf.size());
        throw std::runtime_error(ss.str());
    }
}

/**
 * @brief Write a volume as pasting into an existing file.
 *
 * @param[in] v Input volume.
 * @param[in] filename The name of the file.
 * @param[in] first_voxel The index of the first voxel which the volume should start from.
 */
template<typename TVOL>
void tom::io::write_to_ccp4_paste(const tom::Volume<TVOL> &v, const std::string &filename, const uint32_t *first_voxel) {

    if (v.getSizeX() > std::numeric_limits<uint32_t>::max() ||
        v.getSizeY() > std::numeric_limits<uint32_t>::max() ||
        v.getSizeZ() > std::numeric_limits<uint32_t>::max()) {
        throw std::runtime_error("tom::io::write_to_ccp4_paste: The volume is too large to save it as ccp4-format.");
    }

    int res  = tom_io_ccp4_write_paste(filename.c_str(), &v.get(), tom::get_tom_io_type<TVOL>(),
                                       v.getSizeX(), v.getSizeY(), v.getSizeZ(),
                                       v.getStrideX(), v.getStrideY(),v.getStrideZ(),
                                       NULL, NULL, true, first_voxel, 0);
    if (res != TOM_ERR_OK) {
        const int errno_ = errno;
        std::stringstream ss;
        std::vector<char> buf(1024);
        ss << "Error pasting to ccp4-file \"" << filename << "\". " << tom_io_strerror(res, errno_, &buf[0], buf.size());
        throw std::runtime_error(ss.str());
    }
}

/**
 * @brief Check the name of the given file and retrieve the type of the file according to the suffix (case insensitive)
 *
 * @param[in] filename The name of the file
 *
 * @return Enum type specified in the header file
 */
int tom::io::check_name_suffix(const std::string &filename)
{
	if(filename.find(std::string(".em"))!=std::string::npos || filename.find(std::string(".EM"))!=std::string::npos)
		return EM;
	else if(filename.find(std::string(".mrc"))!=std::string::npos || filename.find(std::string(".MRC"))!=std::string::npos)
		return MRC;
	else if(filename.find(std::string(".cpp"))!=std::string::npos || filename.find(std::string(".CPP"))!=std::string::npos)
		return CCP4;
	else
		return UNKNOWN;
}

/**
 * @brief Get the file's type by checking the validity and the surffix of the file
 *
 * @param[in] filename The name of the file
 *
 * @return Enum type specified in the header file
 */
int tom::io::get_file_type(const std::string &filename)
{
	int file_suffix = check_name_suffix(filename);

	bool matchEM = false;
	bool matchMRC = false;
	bool matchCCP4 = false;

    tom_io_em_header e_header;
    if(tom_io_em_read_header(filename.c_str(), "rb", &e_header, NULL) == TOM_ERR_OK)
    	matchEM = true;

    tom_io_ccp4_header c_header;
    if(tom_io_ccp4_read_header(filename.c_str(), "rb", &c_header, NULL) == TOM_ERR_OK)
    	matchCCP4 = true;

    tom_io_mrc_header m_header;
    if(tom_io_mrc_read_header(filename.c_str(), "rb", &m_header, NULL) == TOM_ERR_OK)
    	matchMRC = true;

    switch(file_suffix)
    {
		case EM:
			if(matchEM == true)
				return EM;
			else
				break;
		case MRC:
			if(matchMRC == true)
				return MRC;
			else
				break;
		default:
			break;
    }

    if(matchCCP4 == true)
    	return CCP4;
    else if(matchMRC == true)
    	return MRC;
    else if(matchEM == true)
    	return EM;
    else
    	return UNKNOWN;
}

/**
 * @brief Wrap function for reading all kinds of files (EM, MRC, CCP4).
 * Cannot retrieve the header information!
 *
 * @param[out] v A reference to a pointer where the volume is returned.
 * @param[in] filename The filename of the file.
 * @param[in] file_type Read the file as type specified here, default is unknown.
 * You can also use the surffix of the file to specify the file type you want to achieve.
 * @param[in] subregion The subregion to read.
 * @param[in] sampling The sampling factor.
 * @param[in] binning The binning factor.
 */
template<typename T>
void tom::io::read_from_file(tom::Volume<T> *&v, const std::string &filename, const int file_type, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning)
{
	int f_local = file_type;

	if(f_local == UNKNOWN)
		f_local = get_file_type(filename);

	switch(f_local)
	{
		case UNKNOWN:
			char errorMsg[256];
			sprintf(errorMsg, "%s: Wrong file format or file doesn't exist!", filename.c_str());
			throw std::logic_error(std::string(errorMsg));
			return;
		case CCP4:
			read_from_ccp4(v, filename, subregion, sampling, binning, NULL, NULL, NULL);
			return;
		case MRC:
			read_from_mrc(v, filename, subregion, sampling, binning, NULL, NULL, NULL);
			return;
		case EM:
			tom::io::read_from_em(v, filename, subregion, sampling, binning, NULL, NULL, NULL);
			return;
		default:
			return;
	}
}

/**
 * @brief Wrap function for writing all kinds of files (EM, MRC, CCP4).
 * Cannot give the header to be written!
 *
 * @param[in] v Input volume.
 * @param[in] filename The filename of the file.
 * @param[in] file_type Write the file as type specified here, default is unknown.
 * You can also use the surffix of the file to specify the file type you want to achieve.
 */
template<typename TVOL, typename TFILE>
void tom::io::write_to_file(const tom::Volume<TVOL> &v, const std::string &filename, const int file_type)
{
	int f_local = file_type;

	if(f_local == UNKNOWN)
		f_local = check_name_suffix(filename);

	switch(f_local)
	{
		case CCP4:
			write_to_ccp4<TVOL, TFILE>(v, filename, NULL);
			return;
		case MRC:
			write_to_mrc<TVOL, TFILE>(v, filename, NULL);
			return;
		case EM:
			write_to_em<TVOL, TFILE>(v, filename, NULL);
			return;
		default:
			throw std::logic_error("Must specify a valid type!");
			return;
	}
}

/**
 * @brief Write a volume as pasting into an existing file.
 *
 * @param[in] v Input volume.
 * @param[in] filename The name of the file.
 * @param[in] first_voxel The index of the first voxel which the volume should start from.
 */
template<typename TVOL>
void tom::io::write_to_file_paste(const tom::Volume<TVOL> &v, const std::string &filename, const uint32_t *first_voxel, const int file_type)
{
	int f_local = file_type;

	if(f_local == UNKNOWN)
		f_local = check_name_suffix(filename);

	switch(f_local)
	{
		case CCP4:
			write_to_ccp4_paste<TVOL>(v, filename, first_voxel);
			return;
		case MRC:
			write_to_mrc_paste<TVOL>(v, filename, first_voxel);
			return;
		case EM:
			write_to_em_paste<TVOL>(v, filename, first_voxel);
			return;
		default:
			throw std::logic_error("Must specify a valid type!");
			return;
	}
}

// TEMPLATE INSTANTIATIONS

template void tom::io::read_from_em<int32_t>(tom::Volume<int32_t> *&v, const std::string &filename, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning, tom_io_em_header *header, void *(*fcn_malloc)(size_t size), void (*fcn_free)(void *ptr));
template void tom::io::read_from_em<float  >(tom::Volume<float  > *&v, const std::string &filename, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning, tom_io_em_header *header, void *(*fcn_malloc)(size_t size), void (*fcn_free)(void *ptr));
template void tom::io::read_from_em<double >(tom::Volume<double > *&v, const std::string &filename, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning, tom_io_em_header *header, void *(*fcn_malloc)(size_t size), void (*fcn_free)(void *ptr));


template std::auto_ptr<tom::Volume<float                > > tom::io::VolumeEM::read_volume<float                >(const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning);
template std::auto_ptr<tom::Volume<double               > > tom::io::VolumeEM::read_volume<double               >(const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning);
template std::auto_ptr<tom::Volume<std::complex<float > > > tom::io::VolumeEM::read_volume<std::complex<float > >(const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning);
template std::auto_ptr<tom::Volume<std::complex<double> > > tom::io::VolumeEM::read_volume<std::complex<double> >(const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning);
template void tom::io::VolumeEM::read_volume<float                >(tom::Volume<float                > &v, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning);
template void tom::io::VolumeEM::read_volume<double               >(tom::Volume<double               > &v, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning);
template void tom::io::VolumeEM::read_volume<std::complex<float > >(tom::Volume<std::complex<float > > &v, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning);
template void tom::io::VolumeEM::read_volume<std::complex<double> >(tom::Volume<std::complex<double> > &v, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning);


template void tom::io::write_to_em<char                , char                 >(const tom::Volume<char                 > &v, const std::string &filename, const tom_io_em_header *header);
template void tom::io::write_to_em<float               , float                >(const tom::Volume<float                > &v, const std::string &filename, const tom_io_em_header *header);
template void tom::io::write_to_em<float               , double               >(const tom::Volume<float                > &v, const std::string &filename, const tom_io_em_header *header);
template void tom::io::write_to_em<double              , float                >(const tom::Volume<double               > &v, const std::string &filename, const tom_io_em_header *header);
template void tom::io::write_to_em<double              , double               >(const tom::Volume<double               > &v, const std::string &filename, const tom_io_em_header *header);
template void tom::io::write_to_em<std::complex<float >, std::complex<float > >(const tom::Volume<std::complex<float > > &v, const std::string &filename, const tom_io_em_header *header);
template void tom::io::write_to_em<std::complex<double>, std::complex<float > >(const tom::Volume<std::complex<double> > &v, const std::string &filename, const tom_io_em_header *header);
template void tom::io::write_to_em<std::complex<double>, std::complex<double> >(const tom::Volume<std::complex<double> > &v, const std::string &filename, const tom_io_em_header *header);

template void tom::io::write_to_em_paste(const tom::Volume<float> &v, const std::string &filename, const uint32_t *first_voxel);
template void tom::io::write_to_em_paste(const tom::Volume<double> &v, const std::string &filename, const uint32_t *first_voxel);

template void tom::io::read_from_mrc<char   >(tom::Volume<char   > *&v, const std::string &filename, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning, tom_io_mrc_header *header, void *(*fcn_malloc)(size_t size), void (*fcn_free)(void *ptr));
template void tom::io::read_from_mrc<int16_t>(tom::Volume<int16_t> *&v, const std::string &filename, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning, tom_io_mrc_header *header, void *(*fcn_malloc)(size_t size), void (*fcn_free)(void *ptr));
template void tom::io::read_from_mrc<double >(tom::Volume<double > *&v, const std::string &filename, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning, tom_io_mrc_header *header, void *(*fcn_malloc)(size_t size), void (*fcn_free)(void *ptr));
template void tom::io::read_from_mrc<float  >(tom::Volume<float  > *&v, const std::string &filename, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning, tom_io_mrc_header *header, void *(*fcn_malloc)(size_t size), void (*fcn_free)(void *ptr));

template void tom::io::write_to_mrc<char                , char                 >(const tom::Volume<char                 > &v, const std::string &filename, const tom_io_mrc_header *header);
template void tom::io::write_to_mrc<float               , float                >(const tom::Volume<float                > &v, const std::string &filename, const tom_io_mrc_header *header);
template void tom::io::write_to_mrc<double              , float                >(const tom::Volume<double               > &v, const std::string &filename, const tom_io_mrc_header *header);
template void tom::io::write_to_mrc<float               , double               >(const tom::Volume<float                > &v, const std::string &filename, const tom_io_mrc_header *header);
template void tom::io::write_to_mrc<double              , double               >(const tom::Volume<double               > &v, const std::string &filename, const tom_io_mrc_header *header);
template void tom::io::write_to_mrc<std::complex<float >, std::complex<float > >(const tom::Volume<std::complex<float > > &v, const std::string &filename, const tom_io_mrc_header *header);
template void tom::io::write_to_mrc<std::complex<double>, std::complex<float > >(const tom::Volume<std::complex<double> > &v, const std::string &filename, const tom_io_mrc_header *header);
template void tom::io::write_to_mrc<std::complex<double>, std::complex<double> >(const tom::Volume<std::complex<double> > &v, const std::string &filename, const tom_io_mrc_header *header);
template void tom::io::write_to_mrc_paste(const tom::Volume<float> &v, const std::string &filename, const uint32_t *first_voxel);
template void tom::io::write_to_mrc_paste(const tom::Volume<double> &v, const std::string &filename, const uint32_t *first_voxel);

template void tom::io::read_from_ccp4<char   >(tom::Volume<char   > *&v, const std::string &filename, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning, tom_io_ccp4_header *header, void *(*fcn_malloc)(size_t size), void (*fcn_free)(void *ptr));
template void tom::io::read_from_ccp4<double >(tom::Volume<double > *&v, const std::string &filename, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning, tom_io_ccp4_header *header, void *(*fcn_malloc)(size_t size), void (*fcn_free)(void *ptr));
template void tom::io::read_from_ccp4<float  >(tom::Volume<float  > *&v, const std::string &filename, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning, tom_io_ccp4_header *header, void *(*fcn_malloc)(size_t size), void (*fcn_free)(void *ptr));

template void tom::io::write_to_ccp4<char                , char                 >(const tom::Volume<char                 > &v, const std::string &filename, const tom_io_ccp4_header *header);
template void tom::io::write_to_ccp4<float               , float                >(const tom::Volume<float                > &v, const std::string &filename, const tom_io_ccp4_header *header);
template void tom::io::write_to_ccp4<double              , float                >(const tom::Volume<double               > &v, const std::string &filename, const tom_io_ccp4_header *header);
template void tom::io::write_to_ccp4<float               , double               >(const tom::Volume<float                > &v, const std::string &filename, const tom_io_ccp4_header *header);
template void tom::io::write_to_ccp4<double              , double               >(const tom::Volume<double               > &v, const std::string &filename, const tom_io_ccp4_header *header);
template void tom::io::write_to_ccp4<std::complex<float >, std::complex<float > >(const tom::Volume<std::complex<float > > &v, const std::string &filename, const tom_io_ccp4_header *header);
template void tom::io::write_to_ccp4<std::complex<double>, std::complex<float > >(const tom::Volume<std::complex<double> > &v, const std::string &filename, const tom_io_ccp4_header *header);
template void tom::io::write_to_ccp4<std::complex<double>, std::complex<double> >(const tom::Volume<std::complex<double> > &v, const std::string &filename, const tom_io_ccp4_header *header);
template void tom::io::write_to_ccp4_paste(const tom::Volume<float> &v, const std::string &filename, const uint32_t *first_voxel);
template void tom::io::write_to_ccp4_paste(const tom::Volume<double> &v, const std::string &filename, const uint32_t *first_voxel);

template void tom::io::read_from_file<char   >(tom::Volume<char   > *&v, const std::string &filename, const int file_type, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning);
template void tom::io::read_from_file<float  >(tom::Volume<float  > *&v, const std::string &filename, const int file_type, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning);
template void tom::io::read_from_file<double >(tom::Volume<double > *&v, const std::string &filename, const int file_type, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning);
template void tom::io::read_from_file<int16_t>(tom::Volume<int16_t> *&v, const std::string &filename, const int file_type, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning);
template void tom::io::read_from_file<int32_t>(tom::Volume<int32_t> *&v, const std::string &filename, const int file_type, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning);

template void tom::io::write_to_file<char                , char                 >(const tom::Volume<char                 > &v, const std::string &filename, const int file_type);
template void tom::io::write_to_file<float               , float                >(const tom::Volume<float                > &v, const std::string &filename, const int file_type);
template void tom::io::write_to_file<double              , float                >(const tom::Volume<double               > &v, const std::string &filename, const int file_type);
template void tom::io::write_to_file<float               , double               >(const tom::Volume<float                > &v, const std::string &filename, const int file_type);
template void tom::io::write_to_file<double              , double               >(const tom::Volume<double               > &v, const std::string &filename, const int file_type);

template void tom::io::write_to_file_paste<float >(const tom::Volume<float > &v, const std::string &filename, const uint32_t *first_voxel, const int file_type);
template void tom::io::write_to_file_paste<double>(const tom::Volume<double> &v, const std::string &filename, const uint32_t *first_voxel, const int file_type);


