/****************************************************************************//**
 * \file io.h
 * \brief The header file of the io-functions in tom_io.c
 * \author  Thomas Haller
 * \version 0.1
 * \date    12.11.2007
 *
 * Contains functions to read and write data from/to a file.
 * The used format is EM
 *******************************************************************************/
#ifndef ___INCLUDE_CORE__IO_H__
#define ___INCLUDE_CORE__IO_H__

#define _FILE_OFFSET_BITS 64
#ifndef __LARGEFILE_SOURCE
    #define __LARGEFILE_SOURCE
#endif
#define _FILE_OFFSET_BITS 64


#ifdef __cplusplus
extern "C" {
#endif


#include <stddef.h>
#include <stdio.h>
#if defined __GNUC__ || defined __xlC__
    #include <stdint.h>
#else
    /* for systems that do not provide stdint.h */
    #include <system/pstdint.h>
#endif

#if defined _WIN32 && defined __LCC__
    #define int64_t int64_T
    #define uint64_t uint64_T
#endif




/* Error/returning codes by (some) functions. Make
 * these values nagtive!!!! */
#define TOM_ERR_OK                        ((int)( 0))
#define TOM_ERR_WRONG_INPUT               ((int)(-1))
#define TOM_ERR_WRONG_HEADER              ((int)(-2))
#define TOM_ERR_SUBREGION_OUT             ((int)(-3))
#define TOM_ERR_NOT_YET_IMPLEMENTED       ((int)(-4))
#define TOM_ERR_MALLOC                    ((int)(-5))
#define TOM_ERR_OPEN_FILE                 ((int)(-6))
#define TOM_ERR_FILESIZE_MISMATCH         ((int)(-7))
#define TOM_ERR_READ_FILE                 ((int)(-8))
#define TOM_ERR_BINNING_TOO_HIGH          ((int)(-9))
#define TOM_ERR_IOTYPE_NOT_SUPPORTED      ((int)(-10))
#define TOM_ERR_WRONG_IOTYPE_CONVERSION   ((int)(-11))
#define TOM_ERR_WRITE_FILE                ((int)(-12))
#define TOM_ERR_NO_COMPLEX_BINNING        ((int)(-13))
#define TOM_ERR_WRONG_DATA_SIZE           ((int)(-14))
#define TOM_ERR_VOLUME_TOO_LARGE          ((int)(-15))
#define TOM_ERR_VOLUME_TOO_LARGE_FOR_EM   ((int)(-16))




/* These are integer numbers defining the raw datatype.
 * The numerical value corresponds to the parameter in the em-header.
 * This is arbitrary. Don't rely on the numerical values of the defines.
 * Only make them >= 1. See also iotype_datasize */
#define TOM_IO_TYPE_UNKNOWN                    ((int)0)
#define TOM_IO_TYPE_INT8                       ((int)1)
#define TOM_IO_TYPE_INT16                      ((int)2)
#define TOM_IO_TYPE_INT32                      ((int)4)
#define TOM_IO_TYPE_FLOAT                      ((int)5)
#define TOM_IO_TYPE_COMPLEX32                  ((int)8)
#define TOM_IO_TYPE_DOUBLE                     ((int)9)
#define TOM_IO_TYPE_COMPLEX64                  ((int)10)
#define TOM_IO_TYPE_INT64                      ((int)11)

#define TOM_IO_TYPE_UINT16                     ((int)6)
#define TOM_IO_TYPE_COMPLEX16                  ((int)3)




/****************************************************************************//**
 * \brief Structure for the em-header.
 *
 * A binary compatible structure to the header of the em-file.
 * It contains the exactly same bits as stored in the file, exept the case
 * when the integers are swaped.
 *******************************************************************************/
typedef struct tom_io_em_header {
    int8_t machine;             /**< Byte 1: Machine Coding
                                            (OS-9;      0),
                                            (VAX;       1),
                                            (Convex;    2),
                                            (SGI;       3),
                                            (Mac;       5),
                                            (PC;        6). */
    int8_t byte2;               /**< General purpose. On OS-9 system: 0 old version 1 is new version. */
    int8_t byte3;               /**< Not used in standard EM-format, if this byte is 1 the header is abandoned. */
    int8_t type;                /**< Data Type Coding. */
    uint32_t dims[3];           /**< Three long integers (3x4 bytes) are image size in x, y, z Dimension. */
    int8_t comment[80];         /**< 80 Characters as comment. */
    int32_t emdata[40];         /**< 40 long integers (4 x 40 bytes) are user defined parameters. */
    int8_t  userdata[256];      /**< 256 Byte with userdata, i.e. the username. */
} tom_io_em_header;



/****************************************************************************//**
 * \brief Complex element corresponding to the complex em-type.
 *
 * Format to save complex data in em-format.
 *******************************************************************************/
typedef struct {
    float re; /**< real part */
    float im; /**< imaginary part */
} tom_io_em_complex32;

/**************************************************************
**************//**
 * \brief Complex element corresponding to the complex em-type.
 *
 * Format to save complex data in em-format.
 *******************************************************************************/
typedef struct {
    double re; /**< real part */
    double im; /**< imaginary part */
} tom_io_em_complex64;



int tom_io_iotype_datasize(int type);


#define TOM_IO_MACHINE_BYTEORDER_BIG_ENDIAN      0
#define TOM_IO_MACHINE_BYTEORDER_LITTLE_ENDIAN   1
int tom_io_machine_byteorder();
void tom_io_swap(void *x, size_t size);

int tom_io_em_read_header(const char *fname, const char *mode, tom_io_em_header *header, FILE **fhandle);
int tom_io_em_is_valid_header(const void *header, size_t size);
int tom_io_em_is_swaped(const tom_io_em_header *header);
int tom_io_em_swap_header(tom_io_em_header *header);
int tom_io_em_iotype2emtype(int);
int tom_io_em_emtype2iotype(int);
int tom_io_em_get_iotype_header(const tom_io_em_header *header);
int tom_io_em_set_iotype_header(tom_io_em_header *header, int iotype);
void tom_io_em_get_microscope_name(const tom_io_em_header *header, char *name, size_t bytes_reserved);

int tom_io_em_read( const char *fname, const uint32_t *subregion_, const uint32_t *sampling, const uint32_t *binning, tom_io_em_header *header,
                    void **data, uint32_t *dims, int *restype,
                    void *(*fcn_malloc)(size_t size),
                    void  (*fcn_free)(void *ptr));
int tom_io_read_vol_sac(FILE *f, const uint32_t voldims[], int iotype, int swaped,
                    const uint32_t *subregion_, const uint32_t *sampling_, const uint32_t *binning_,
                    void **data, uint32_t *dims_, int *restype_,
                    void *(*fcn_malloc)(size_t size),
                    void  (*fcn_free)(void *ptr));
int tom_io_read_vol(FILE *f, const uint32_t voldims[], int iotype, int swaped,
                    const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning,
                    void *pdata, int restype);
int tom_io_calculate_sizes(const uint32_t *voldims,
                    const uint32_t *subregion_, uint32_t *subregion_out,
                    const uint32_t *sampling_, uint32_t *sampling_out,
                    const uint32_t *binning_, uint32_t *binning_out,
                    uint32_t *dims_sampled_out,
                    uint32_t *dims_out);

/*int tom_io_write_vol(FILE *f, int iotype_write, const void *data, int iotype, size_t sizex, size_t sizey, size_t sizez, size_t stridex, size_t stridey, size_t stridez, int swap_bytes);*/
int tom_io_write_vol(       FILE *f, int iotype_write, const void *data, int iotype,
                            size_t sizex, size_t sizey, size_t sizez,
                            size_t stridex, size_t stridey, size_t stridez, int swap_bytes,
                            size_t fstridex, size_t fstridey, size_t fstridez);

int tom_io_em_write(const char *filename, const tom_io_em_header *header, const void *data, int iotype, size_t stridex, size_t stridey, size_t stridez);
int tom_io_em_write_append_stack(const char *filename, const void *data, int iotype, uint32_t sizex, uint32_t sizey, uint32_t sizez, size_t stridex, size_t stridey, size_t stridez, tom_io_em_header *header_out, int *header_read, int allow_conversion);
int tom_io_em_write_paste(const char *filename, const void *data, int iotype, uint32_t sizex, uint32_t sizey, uint32_t sizez, size_t stridex, size_t stridey, size_t stridez, tom_io_em_header *header_out, int *header_read, int allow_conversion, const uint32_t *first_voxel, const uint32_t *sampling);

int tom_io_check_stride(size_t size_of_type, size_t sizex, size_t sizey, size_t sizez, size_t *stridex, size_t *stridey, size_t *stridez);

char *tom_io_strerror(int tom_errnum, int errno_, char *buf, size_t buflen);


/****************************************************************************//**
 * @brief Structure of MRC header
 * The format can be available at
 * http://ami.scripps.edu/software/mrctools/mrc_specification.php
 *******************************************************************************/
typedef struct tom_io_mrc_header {

	uint32_t nx; ///< number of columns
	uint32_t ny; ///< number of rows
	uint32_t nz; ///< number of sections

	/**
	 * MODE     data type :
	 *  0        image : signed 8-bit bytes range -128 to 127
	 *  1        image : 16-bit halfwords
	 *  2        image : 32-bit reals
	 *  3        transform : complex 16-bit integers
	 *  4        transform : complex 32-bit reals
	 *  6        image : unsigned 16-bit range 0 to 65535
	 */
	uint32_t mode;
	int32_t nxstart; ///< number of first column in map
	int32_t nystart; ///< number of first row in map
	int32_t nzstart; ///< number of first section in map

	int32_t mx; ///< number of intervals along X
	int32_t my; ///< number of intervals along Y
	int32_t mz; ///< number of intervals along Z

	int32_t cella[3]; ///< cell dimensions in angstroms
	int32_t cellb[3]; ///< cell angles in degrees

	int32_t mapc; ///< axis corresp to cols (1,2,3 for X,Y,Z)
	int32_t mapr; ///< axis corresp to rows (1,2,3 for X,Y,Z)
	int32_t maps; ///< axis corresp to sections (1,2,3 for X,Y,Z)

	int32_t dmin; ///< minimum density value
	int32_t dmax; ///< maximum density value
	int32_t dmean; ///< mean density value

	uint32_t ispg; ///< space group number 0 or 1
	uint32_t nsymbt; ///< number of bytes used for symmetry data
	int32_t extra[25]; ///< extra space used for anything
	int32_t origin[3]; ///< origin in X,Y,Z used for transforms

	int8_t map[4]; ///< character string 'MAP' to identify file type
	int32_t machst; ///< machine stamp
	int32_t rms; ///< rms deviation of map from mean density
	uint32_t nlabl; ///< number of labels being used

	int8_t label[800]; ///< LABEL(20,10) 10 80-character text labels

} tom_io_mrc_header;

int tom_io_mrc_read_header(const char *filename, const char *mode, tom_io_mrc_header *header, FILE **fhandle);
int tom_io_mrc_get_iotype_header(const tom_io_mrc_header *header);
int tom_io_mrc_read(const char *filename, const uint32_t *subregion_, const uint32_t *sampling, const uint32_t *binning, tom_io_mrc_header *header,
                    void **data, uint32_t *dims, int *restype,
                    void *(*fcn_malloc)(size_t size),
                    void  (*fcn_free)(void *ptr));

int tom_io_mrc_write(const char *filename, const tom_io_mrc_header *header, const void *data, int iotype, size_t stridex, size_t stridey, size_t stridez);
int tom_io_mrc_write_paste(const char *filename, const void *data, int iotype, uint32_t sizex, uint32_t sizey, uint32_t sizez, size_t stridex, size_t stridey, size_t stridez, tom_io_mrc_header *header_out, int *header_read, int allow_conversion, const uint32_t *first_voxel, const uint32_t *sampling);


/**
 * @brief Sturcture of CCP4 header
 * The format can be available at
 * http://www.ccp4.ac.uk/html/maplib.html#description
 */
typedef struct tom_io_ccp4_header {

	uint32_t nc; ///< Number of Columns (fastest changing in map)
	uint32_t nr; ///< Number of rows
	uint32_t ns; ///< Number of sections (slowest changing in map)

	/**
	 * MODE     Data type :
	 *  0        envelope stored as signed bytes (from -128 lowest to 127 highest)
	 *  1        image : stored as Integer*2 (16 bit integer)
	 *  2        image : stored as Reals (32 bit real)
	 *  3        transform : stored as Complex Integer*2
	 *  4        transform : stored as Complex Reals
	 *  5        NONE
	 *  Note: Mode 2 is the normal mode used in the CCP4 programs. Other modes than 2 and 0 may NOT WORK!
	 */
	int32_t mode;
	int32_t ncstart; ///< Number of first COLUMN  in map
	int32_t nrstart; ///< Number of first ROW     in map
	int32_t nsstart; ///< Number of first SECTION in map

	int32_t nx; ///< Number of intervals along X
	int32_t ny; ///< Number of intervals along Y
	int32_t nz; ///< Number of intervals along Z

	int32_t x_length; ///< Cell Dimensions (Angstroms)
	int32_t y_length; ///< Cell Dimensions (Angstroms)
	int32_t z_length; ///< Cell Dimensions (Angstroms)

	int32_t alpha; ///< Cell Angles     (Degrees)
	int32_t beta;  ///< Cell Angles     (Degrees)
	int32_t gamma; ///< Cell Angles     (Degrees)

	int32_t mapc; ///< Which axis corresponds to Cols.  (1,2,3 for X,Y,Z)
	int32_t mapr; ///< Which axis corresponds to Rows   (1,2,3 for X,Y,Z)
	int32_t maps; ///< Which axis corresponds to Sects. (1,2,3 for X,Y,Z)

	int32_t amin;  ///< Minimum density value
	int32_t amax;  ///< Maximum density value
	int32_t amean; ///< Mean    density value

	int32_t ispg; ///< Space group number
	int32_t nsympg; ///< Number of bytes used for storing symmetry operators

	int32_t lskflg; ///< Flag for skew transformation, =0 none, =1 if foll
	int32_t skwmat[9]; ///< Skew matrix S (in order S11, S12, S13, S21 etc) if LSKFLG .ne. 0.
	/**
	 * Skew translation t if LSKFLG .ne. 0.
	 * Skew transformation is from standard orthogonal
	 * coordinate frame (as used for atoms) to orthogonal map frame, as
	 * Xo(map) = S * (Xo(atoms) - t)
	 */
	int32_t skwtrn[3];

	int32_t future_use[15]; ///< All set to zero by default
	int8_t map[4]; ///< Character string 'MAP ' to identify file type
	int32_t machst; ///< Machine stamp indicating the machine type which wrote file

	int32_t arms; ///< Rms deviation of map from mean density
	int32_t nlabl; ///< Number of labels being used
	int8_t label[800]; ///< LABEL(20,10) 10 80-character text labels

} tom_io_ccp4_header;

int tom_io_ccp4_read_header(const char *filename, const char *mode, tom_io_ccp4_header *header, FILE **fhandle);
int tom_io_ccp4_get_iotype_header(const tom_io_ccp4_header *header);
int tom_io_ccp4_read(const char *filename, const uint32_t *subregion_, const uint32_t *sampling, const uint32_t *binning, tom_io_ccp4_header *header,
                    void **data, uint32_t *dims, int *restype,
                    void *(*fcn_malloc)(size_t size),
                    void  (*fcn_free)(void *ptr));

int tom_io_ccp4_write(const char *filename, const tom_io_ccp4_header *header, const void *data, int iotype, size_t stridex, size_t stridey, size_t stridez);
int tom_io_ccp4_write_paste(const char *filename, const void *data, int iotype, uint32_t sizex, uint32_t sizey, uint32_t sizez, size_t stridex, size_t stridey, size_t stridez, tom_io_ccp4_header *header_out, int *header_read, int allow_conversion, const uint32_t *first_voxel, const uint32_t *sampling);

// Header transformation
void tom_io_em2mrc_header(const tom_io_em_header *em_header, tom_io_mrc_header *mrc_header);
void tom_io_em2ccp4_header(const tom_io_em_header *em_header, tom_io_ccp4_header *ccp4_header);

void tom_io_mrc2em_header(const tom_io_mrc_header *mrc_header, tom_io_em_header *em_header);
void tom_io_mrc2ccp4_header(const tom_io_mrc_header *mrc_header, tom_io_ccp4_header *ccp4_header);

void tom_io_ccp42em_header(const tom_io_ccp4_header *ccp4_header, tom_io_em_header *em_header);
void tom_io_ccp42mrc_header(const tom_io_ccp4_header *ccp4_header, tom_io_mrc_header *mrc_header);

#ifdef __cplusplus
} /* extern "C" */
#endif




#endif
