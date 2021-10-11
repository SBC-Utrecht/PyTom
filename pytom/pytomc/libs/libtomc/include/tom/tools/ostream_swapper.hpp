/****************************************************************************//**
 * \file ostream_swapper.hpp
 * \brief Contains the class tom::tools::ostream_swapper.
 * \author  Thomas Haller
 * \version 0.1.0
 * \date    24.01.2008
 *
 * Taken from \url http://groups.google.com/group/borland.public.cppbuilder.language/browse_thread/thread/52e6bab7756110b0/dfc46e2bfaaa6cbc?lnk=st&#dfc46e2bfaaa6cbc
 *******************************************************************************/
#ifndef ___INCLUDE_OSTREAM_SWAPPER_HPP__
#define ___INCLUDE_OSTREAM_SWAPPER_HPP__



#include <ostream>


namespace tom {
namespace tools {


/****************************************************************************//**
 * \brief Sets the rdbuf of an std::ostreams to redirect the output of a steam.
 *
 * The destructor resets the original buffers.
 *******************************************************************************/
class ostream_swapper {

public:
    /****************************************************************************//**
    * \param[in,out] stream The rdbuf of this stream is set to the one of \c sreplace.
    *   Thereby writing to stream redirects the output to the one of \c sreplace.
    * \param[in,out] sreplace
    * stream must be have a longer live time than the ostream_swapper because the
    * destructor tries to reset the old buffer.
    *******************************************************************************/
    ostream_swapper(std::ostream &stream, std::ostream &sreplace)
        : swapped_back_(true), rdbuf_(stream.rdbuf()), stream_(stream) {
        std::streambuf *srdbuf = sreplace.rdbuf();
        stream.rdbuf(srdbuf);
        swapped_back_ = false;
    }
    ~ostream_swapper() {
        this->swap_back();
    }
    /****************************************************************************//**
    * Reset the original content of the stream before destruction of the object.
    * Afterwords it is save that the stream is destroyed even before the
    * ostream_swapper object.
    *******************************************************************************/
    void swap_back() {
        if (!swapped_back_) {
            // set the status variable first, in case of an exception...
            swapped_back_ = true;
            stream_.rdbuf(rdbuf_);
        }
    }
private:
    bool swapped_back_;
    std::streambuf *rdbuf_;
    std::ostream &stream_;
}; // class tom::tools::ostream_swapper


} // namespace tools
} // namespace tom


#endif

