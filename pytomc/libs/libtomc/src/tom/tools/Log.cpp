/***************************************************************************//**
 * \file Log.cpp
 * \brief Implementations of Log.hpp
 * \author  Michael Stoelken
 * \version 0.1
 * \date    29.01.2009
 *
 * Implementations of Log.hpp
 ******************************************************************************/
#include <tom/tools/Log.hpp>

#include <iostream>


#define COUT std::cout << __FILE__ << ":" << __LINE__ << ": "



#include <tom/tools/snippets.hpp>


extern const tom::tools::s_Endl tom::tools::endl = tom::tools::s_Endl();


/***************************************************************************//**
 * Constructor
 ******************************************************************************/
tom::tools::Log::Log(std::ostream *s_default, char separator)
    :   s_default_(s_default),
        s_(0),
        swp_(0),
        entry_closed_(true),
        is_redirected_(false),
        old_level_(INFO),
        old_msgid_(""),
        setw_tag_(0),
        separator_(separator) {

    std::auto_ptr<std::ostringstream> s(new std::ostringstream());
    s_ = s.release();
}




/***************************************************************************//**
 * Destructor
 ******************************************************************************/
tom::tools::Log::~Log() {

    if (!entry_closed_) {
        *s_ << '\n';
    }
    if (!is_redirected_ && s_default_) {
        *s_default_ << static_cast<std::ostringstream *>(s_)->str();
    }
    delete swp_;
    delete s_;
}


/***************************************************************************//**
 *
 ******************************************************************************/
void tom::tools::Log::redirect(const std::string &filename) {
    if (is_redirected_) {
        throw std::runtime_error("Logging class already redirected.");
    }
    std::auto_ptr<std::ostream> f(new std::ofstream(filename.c_str(), std::ios::trunc));
    const int errno_ = errno;
    if (!f->good()) {
        std::ostringstream ss;
        ss << "could not open log file \"" << filename << "\" for writing: " << tom::tools::strerror(errno_);
        throw std::runtime_error(ss.str());
    }
    assert(dynamic_cast<std::ostringstream *>(s_));
    *f << static_cast<std::ostringstream *>(s_)->str();
    delete s_;
    s_ = f.release();
    is_redirected_ = true;
}


/***************************************************************************//**
 *
 ******************************************************************************/
void tom::tools::Log::redirect(std::ostream &stream) {
    if (is_redirected_) {
        throw std::runtime_error("Logging class already redirected.");
    }
    assert(dynamic_cast<std::ostringstream *>(s_));
    stream << static_cast<std::ostringstream *>(s_)->str();
    swp_ = new tom::tools::ostream_swapper(*s_, stream);
    is_redirected_ = true;
}

