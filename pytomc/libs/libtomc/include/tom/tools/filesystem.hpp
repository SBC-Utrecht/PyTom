/****************************************************************************//**
 * \file filesystem.hpp
 * \author  Thomas Haller
 * \version 0.2
 * \date    05.03.2007
 * \brief Wrapper for boost::filesystem.
 *
 * Contains the namespace tom::tools::fs.
 * This only uses the namespace of boost::filesystem.
 * It is intended for systems where the boost library is not available
 * for one or the other reasen.\\
 * Defining NO_BOOST_FILESYSTEM would allow an own implementation of the filesystem
 * functions. However this is not yet implemented and probably never will :). \\
 * This file was introduced because at first there was a system where boost::filesystem
 * did not compile.
 *******************************************************************************/
#ifndef ___INCLUDE__TOM__TOOLS__FILESYSTEM_HPP__
#define ___INCLUDE__TOM__TOOLS__FILESYSTEM_HPP__


#ifndef NO_BOOST_FILESYSTEM

#include <boost/filesystem.hpp>


namespace tom {
namespace tools {

    /// Wrappes the boost::filesystem
    namespace fs = boost::filesystem;
} // namespace tools
} // namespace tom


#else
#error not_yet_implemented

// create a basic implementation of the boost/filesystem api.
namespace tom {
namespace tools {
namespace fs {

// TODO

} // namespace fs
} // namespace tools
} // namespace tom

#endif



#endif


