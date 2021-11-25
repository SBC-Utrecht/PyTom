/****************************************************************************//**
 * \file io.hpp
 * \brief The header file of the io-functions in tom_io.cpp
 * \author  Thomas Haller
 * \version 0.2
 * \date    12.01.2009
 *
 * Contains C++ functions to enhance tom/io/io.c.
 *******************************************************************************/
#ifndef ___INCLUDE__TOM__IP__IO_HPP__
#define ___INCLUDE__TOM__IP__IO_HPP__

#include <string>
#include <memory>

#include <tom/io/io.h>
#include <tom/volume.hpp>


namespace tom {

namespace io {




/***************************************************************************//**
 * \brief Class for reading the image.
 *
 * This class opens an em-file and allows to read the header and the volume.
 * This can be used to read the header in a first step and the data in a second,
 * without closing and reopening the file.
 * The volume can be read only once.
 ******************************************************************************/
class VolumeEM {

public:
    VolumeEM(const std::string &filename, void *(*fcn_malloc)(size_t size)=0, void (*fcn_free)(void *ptr)=0, bool close_after_first_read=true);
    VolumeEM(const std::string &filename, tom_io_em_header &header, void *(*fcn_malloc)(size_t size)=0, void (*fcn_free)(void *ptr)=0, bool close_after_first_read=true);
    ~VolumeEM();

    const tom_io_em_header &get_header() const;
    const tom_io_em_header *operator->() const;
    template<typename T> std::auto_ptr<tom::Volume<T> > read_volume(const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning);
    template<typename T> void read_volume(tom::Volume<T> &v, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning);

    template<typename T> bool is_equal_size(const tom::Volume<T> &v) const;
    bool is_equal_size(std::size_t size) const;
    bool is_equal_size(std::size_t sizex, std::size_t sizey, std::size_t sizez) const;

    template<typename T> bool is_equal_size(const tom::Volume<T> &v, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning) const;
    bool is_equal_size(std::size_t sizex, std::size_t sizey, std::size_t sizez, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning) const;
    int get_size(const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning, std::size_t &sizex, std::size_t &sizey, std::size_t &sizez) const;

    const std::string &get_filename() const;

private:
    VolumeEM(const VolumeEM &v): header_(v.header_) { assert(!"HIDDEN");               }
    VolumeEM &operator=(const VolumeEM &v)          { assert(!"HIDDEN"); return *this; }

    void open_and_seek_file(bool use_volume_header);

    void *(*fcn_malloc_)(size_t size);
    void (*fcn_free_)(void *ptr);

    tom_io_em_header *headerp_;
    tom_io_em_header &header_;
    std::string filename_;
    FILE *f_;
    bool f_is_closed_;
    bool close_after_first_read_;
    bool seeked_at_data_start_;

}; // class tom::io::Volume


void em_read_header(const std::string &filename, tom_io_em_header &header);


template<typename T> void read_from_em(tom::Volume<T> *&v, const std::string &filename, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning, tom_io_em_header *header, void *(*fcn_malloc)(size_t size)=0, void (*fcn_free)(void *ptr)=0);
template<typename T> int get_tom_io_type();


template<typename TVOL, typename TFILE> void write_to_em(const tom::Volume<TVOL> &v, const std::string &filename, const tom_io_em_header *header);
template<typename TVOL>                 void write_to_em(const tom::Volume<TVOL> &v, const std::string &filename, const tom_io_em_header *header);

//template<typename TVOL> void write_to_em_paste(const tom::Volume<TVOL> &v, const std::string &filename, const tom_io_em_header *header, const TVOL *first_voxel);
template<typename TVOL> void write_to_em_paste(const tom::Volume<TVOL> &v, const std::string &filename, const uint32_t *first_voxel);

// Here begins the mrc file part
void mrc_read_header(const std::string &filename, tom_io_mrc_header &header);
void mrc_print_header_info(const tom_io_mrc_header &header);
template<typename T> void read_from_mrc(tom::Volume<T> *&v, const std::string &filename, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning, tom_io_mrc_header *header, void *(*fcn_malloc)(size_t size)=0, void (*fcn_free)(void *ptr)=0);
template<typename TVOL, typename TFILE> void write_to_mrc(const tom::Volume<TVOL> &v, const std::string &filename, const tom_io_mrc_header *header);
template<typename TVOL>                 void write_to_mrc(const tom::Volume<TVOL> &v, const std::string &filename, const tom_io_mrc_header *header);
template<typename TVOL> void write_to_mrc_paste(const tom::Volume<TVOL> &v, const std::string &filename, const uint32_t *first_voxel);

// Here begins the ccp4 file part
void ccp4_read_header(const std::string &filename, tom_io_ccp4_header &header);
void ccp4_print_header_info(const tom_io_ccp4_header &header);
template<typename T> void read_from_ccp4(tom::Volume<T> *&v, const std::string &filename, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning, tom_io_ccp4_header *header, void *(*fcn_malloc)(size_t size)=0, void (*fcn_free)(void *ptr)=0);
template<typename TVOL, typename TFILE> void write_to_ccp4(const tom::Volume<TVOL> &v, const std::string &filename, const tom_io_ccp4_header *header);
template<typename TVOL>                 void write_to_ccp4(const tom::Volume<TVOL> &v, const std::string &filename, const tom_io_ccp4_header *header);
template<typename TVOL> void write_to_ccp4_paste(const tom::Volume<TVOL> &v, const std::string &filename, const uint32_t *first_voxel);

enum FILE_TYPE
{
	UNKNOWN,
	EM,
	MRC,
	CCP4
};

int check_name_suffix(const std::string &filename);
int get_file_type(const std::string &filename);
template<typename T> void read_from_file(tom::Volume<T> *&v, const std::string &filename, const int file_type=UNKNOWN, const uint32_t *subregion=NULL, const uint32_t *sampling=NULL, const uint32_t *binning=NULL);
template<typename TVOL, typename TFILE> void write_to_file(const tom::Volume<TVOL> &v, const std::string &filename, const int file_type=UNKNOWN);
template<typename TVOL                > void write_to_file(const tom::Volume<TVOL> &v, const std::string &filename, const int file_type=UNKNOWN);

template<typename TVOL> void write_to_file_paste(const tom::Volume<TVOL> &v, const std::string &filename, const uint32_t *first_voxel, const int file_type=UNKNOWN);
} // namespace io


// Deprecated versions in the wrong namespace :)
template<typename T> void read_from_em(tom::Volume<T> *&v, const std::string &filename, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning, tom_io_em_header *header, void *(*fcn_malloc)(size_t size)=0, void (*fcn_free)(void *ptr)=0);
template<typename T> int get_tom_io_type();


} // namespace tom






// INLINE DECLARATIONS





/***************************************************************************//**
 * Access the header of the volume.
 ******************************************************************************/
inline const tom_io_em_header &tom::io::VolumeEM::get_header() const {
    return header_;
}


/***************************************************************************//**
 * Bring tom::io::read_from_em from the namespace tom::io to tom, for supporting
 * older definitions.
 ******************************************************************************/
template<typename T>
inline void tom::read_from_em(tom::Volume<T> *&v, const std::string &filename, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning, tom_io_em_header *header, void *(*fcn_malloc)(size_t size), void (*fcn_free)(void *ptr)) {
    tom::io::read_from_em(v, filename, subregion, sampling, binning, header, fcn_malloc, fcn_free);
}



/***************************************************************************//**
 * Bring tom::io::get_tom_io_type from the namespace tom::io to tom, for supporting
 * older definitions.
 ******************************************************************************/
template<typename T>
inline int tom::get_tom_io_type() {
    return tom::io::get_tom_io_type<T>();
}


/***************************************************************************//**
 * Operator to access the header of the volume.
 ******************************************************************************/
inline const tom_io_em_header *tom::io::VolumeEM::operator->() const {
    return &header_;
}


/***************************************************************************//**
 * Check that the VolumeEM has an equal size as a given tom::Volume
 ******************************************************************************/
template<typename T>
inline bool tom::io::VolumeEM::is_equal_size(const tom::Volume<T> &v) const {
    return v.is_equal_size(header_.dims[0], header_.dims[1], header_.dims[2]);
}


/***************************************************************************//**
 * Check the size of the volume
 ******************************************************************************/
inline bool tom::io::VolumeEM::is_equal_size(std::size_t size) const {
    return is_equal_size(size, size, size);
}


/***************************************************************************//**
 * Check the size of the volume
 ******************************************************************************/
inline bool tom::io::VolumeEM::is_equal_size(std::size_t sizex, std::size_t sizey, std::size_t sizez) const {
    return sizex==header_.dims[0] && sizey==header_.dims[1] && sizez==header_.dims[2];
}


/***************************************************************************//**
 * Check whether the em-file has a corresponging size to the volume \c v using a
 * certain \c subregion, \c sampling and \c binning.
 * Uses \c get_size to determine the resulting size and returns also false if the
 * subregion/binning/sampling parameters are out of range.
 ******************************************************************************/
template<typename T>
inline bool tom::io::VolumeEM::is_equal_size(const tom::Volume<T> &v, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning) const {
    std::size_t sizex, sizey, sizez;
    return get_size(subregion, sampling, binning, sizex, sizey, sizez)==TOM_ERR_OK && v.is_equal_size(sizex, sizey, sizez);
}


/***************************************************************************//**
 * Check whether the em-file has a certain size for a given
 * \c subregion, \c sampling and \c binning.
 * Uses \c get_size to determine the resulting size and returns also false if the
 * subregion/binning/sampling parameters are out of range.
 ******************************************************************************/
inline bool tom::io::VolumeEM::is_equal_size(std::size_t sizex, std::size_t sizey, std::size_t sizez, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning) const {
    std::size_t sizex_, sizey_, sizez_;
    return get_size(subregion, sampling, binning, sizex_, sizey_, sizez_)==TOM_ERR_OK && sizex==sizex_ && sizey==sizey_ && sizez==sizez_;
}


/***************************************************************************//**
 * Check whether the em-file has a certain size for a given
 * \c subregion, \c sampling and \c binning.
 * Uses \c get_size to determine the resulting size and returns also false if the
 * subregion/binning/sampling parameters are out of range.
 ******************************************************************************/
inline const std::string &tom::io::VolumeEM::get_filename() const {
    return filename_;
}



/***************************************************************************//**
 * Save volume to em-file.
 *
 * Wrapper for the type converting version of write_to_em. The data-type on disk
 * is equal to the one of the volume.
 *
 * No type conversion is done!
 ******************************************************************************/
template<typename TVOL>
inline void tom::io::write_to_em(const tom::Volume<TVOL> &v, const std::string &filename, const tom_io_em_header *header) {
    write_to_em<TVOL, TVOL>(v, filename, header);
}


/***************************************************************************//**
 * Save volume to mrc-file.
 *
 * Wrapper for the type converting version of write_to_mrc. The data-type on disk
 * is equal to the one of the volume.
 *
 * No type conversion is done!
 ******************************************************************************/

template<typename TVOL>
inline void tom::io::write_to_mrc(const tom::Volume<TVOL> &v, const std::string &filename, const tom_io_mrc_header *header) {
    write_to_mrc<TVOL, TVOL>(v, filename, header);
}

/***************************************************************************//**
 * Save volume to ccp4-file.
 *
 * Wrapper for the type converting version of write_to_ccp4. The data-type on disk
 * is equal to the one of the volume.
 *
 * No type conversion is done!
 ******************************************************************************/
template<typename TVOL>
inline void tom::io::write_to_ccp4(const tom::Volume<TVOL> &v, const std::string &filename, const tom_io_ccp4_header *header) {
    write_to_ccp4<TVOL, TVOL>(v, filename, header);
}

template<typename TVOL>
inline void tom::io::write_to_file(const tom::Volume<TVOL> &v, const std::string &filename, const int file_type)
{
	write_to_file<TVOL, TVOL>(v, filename, file_type);
}

#endif
