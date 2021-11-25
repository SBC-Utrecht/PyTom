/****************************************************************************//**
 * \file ParameterSet.cpp
 * \brief Contains the implementations of ParameterSet.hpp
 * \author  Thomas Haller
 * \version 0.2
 * \date    21.11.2008
 *
 *
 * See the documentation of ParameterSet.hpp
 *******************************************************************************/
#include <tom/tools/ParameterSet.hpp>




#include <cassert>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <errno.h>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <tom/tools/snippets.hpp>




/****************************************************************************//**
 * \brief Constructor to initialize the parameterset from config file.
 *
 * \param[in] filename Name of the configuration (text) file.
 *******************************************************************************/
tom::tools::ParameterSet::ParameterSet(const std::string &filename)
    :   filename_(filename),
        mtag_(),
        mtagit_(),
        midx_() {

    std::vector<std::pair<std::string, std::string> > m = parse_file(filename);

    init(m);
}




/****************************************************************************//**
 * \brief Constructor to initialize the ParameterSet.
 *
 * \param[in] m Vector containing the tuples (tag, value). The order of the
 *  Vector might be relevant for the one that uses the ParameterSet. Therefore
 *  the class does not mix it up.
 *******************************************************************************/
tom::tools::ParameterSet::ParameterSet(const std::vector<std::pair<std::string, std::string> > &m)
    :   filename_(),
        mtag_(),
        mtagit_(),
        midx_() {
    init(m);
}


namespace {
template<typename T> inline T           _parameterset_cast(const std::string &s) { return boost::lexical_cast<T>(boost::trim_copy(s));  }
template<          > inline std::string _parameterset_cast(const std::string &s) { std::string s_(s); boost::trim(s_); return s_;       }
}
/****************************************************************************//**
 * \brief Return the values associated with the tag.
 *
 * \param[in] tag Tag (identifier) of the entry that is looked up.
 * \param[in] first_idx Only parameters with a higher (or equal) index are considered.
 * \param[in] dflt Pointer to the default value if the tag is missing. If 0, an exception
 *      is thrown if the parameter is missing. *
 * \return
 *******************************************************************************/
template<typename T>
T tom::tools::ParameterSet::pop(const std::string &tag, std::size_t first_idx, const T *dflt) {


    std::auto_ptr<tom::tools::ParameterSet::btuple> entry = pop(tag, first_idx);
    std::ostringstream ss;

    T val;
    if (entry.get()) {
        try {
            val = ::_parameterset_cast<T>(VAL(*entry));
        } catch (boost::bad_lexical_cast &e) {
            ss << "Could not cast '" << tag << "'='" << VAL(*entry) << "' to type.";
            throw std::runtime_error(ss.str());
        }
    } else {
        if (dflt) {
            val = *dflt;
        } else {
            entry = get(tag, false, 0);
            if (entry.get()) {
                throw std::runtime_error("No more key '" + tag + "' in parameter set.");
            } else {
                throw std::runtime_error("No key '" + tag + "' in parameter set.");
            }
        }
    }
    return val;
}



/****************************************************************************//**
 * \brief Return the values associated with the tag.
 *
 * \param[in] tag Tag (identifier) of the entry that is looked up.
 * \param[in] first_idx Only parameters with a higher (or equal) index are considered.
 * \return A auto_ptr to the btuple. It contains 0, if there is no entry with the
 *  given tag, or if all entries are already touched/retrieved.
 *
 * Pops the first found matching parameter from the set (with higher index than first_idx).
 * A once retrieved element can not be poped again (however it can still be accessed
 * with \c get).
 * After calling \c pop several times, the function will return 0 as all values
 * are marked as retrieved.
 *******************************************************************************/
std::auto_ptr<tom::tools::ParameterSet::btuple> tom::tools::ParameterSet::pop(const std::string &tag, std::size_t first_idx) {

    std::auto_ptr<btuple> result;

    std::pair<std::multimap<std::string, btuple>::iterator, std::multimap<std::string, btuple>::iterator> range = mtag_.equal_range(tag);

    std::multimap<std::string, btuple>::iterator it = range.first;
    std::multimap<std::string, btuple>::iterator best_hit;
    for (; it!=range.second; ++it) {
        const btuple &value = it->second;
        if (UNT(value) &&                                                       // The value was not yet retrieved or search any.
            (first_idx <= IND(value)) &&                                        // The index is not smaller than allowed.
            (!result.get() || IND(value) < IND(*result)) ) {                    // Take the smallest index found.
            result.reset(new btuple(value));
            best_hit = it;
        }
    }
    if (result.get()) {
        UNT(best_hit->second) = false;
    }
    return result;
}



/****************************************************************************//**
 * \brief Return the values associated with the tag.
 *
 * \param[in] tag String describing the tag/name of the entry (case sensitive).
 * \param[in] search_only_untouched If true, get only returnes those who are not
 *  yet retrieved using pop. Otherwise it also those which were already poped from
 *  the set are regarded.
 * \param[in] first_idx Only parameters with a higher (or equal) index are considered.
 *  By repeatedly calling get with increasing first_idx all parameters with the
 *  same name can be returned.
 *******************************************************************************/
std::auto_ptr<tom::tools::ParameterSet::btuple> tom::tools::ParameterSet::get(const std::string &tag, bool search_only_untouched, std::size_t first_idx) const {

    std::auto_ptr<btuple> result;

    std::pair<std::multimap<std::string, btuple>::const_iterator, std::multimap<std::string, btuple>::const_iterator> range = mtag_.equal_range(tag);
    std::multimap<std::string, btuple>::const_iterator it = range.first;
    for (; it!=range.second; ++it) {
        const btuple &value = it->second;
        if ((!search_only_untouched || UNT(value)) &&                           // The value was not yet retrieved or search any.
            (first_idx <= IND(value)) &&                                        // The index is not smaller than allowed.
            (!result.get() || IND(value) < IND(*result)) ) {                    // Take the smallest index found.
            // Remember the found value.
            result.reset(new btuple(value));
        }
    }
    return result;
}





/****************************************************************************//**
 * \brief Return the values associated with the tag.
 *******************************************************************************/
std::auto_ptr<std::pair<std::string, tom::tools::ParameterSet::btuple> > tom::tools::ParameterSet::pop(std::size_t first_idx) {
    std::multimap<std::string, boost::tuple<std::size_t, bool, std::string> >::iterator it_best = mtag_.end();
    std::multimap<std::string, boost::tuple<std::size_t, bool, std::string> >::iterator it = mtagit_;
    bool found_no_untouched = true;
    bool found_best = false;
    while (it != mtag_.end()) {
        if (UNT(mtagit_->second)) {
            if (IND(mtagit_->second)>=first_idx && (!found_best || IND(it_best->second)>IND(mtagit_->second))) {
                found_best = true;
                it_best = it;
            }
            found_no_untouched = false;
        } else if (found_no_untouched) {
            mtagit_++;
        }
        it++;
    }
    if (it_best != mtag_.end()) {
        UNT(it_best->second) = false;
        return std::auto_ptr<std::pair<std::string, tom::tools::ParameterSet::btuple> >(new std::pair<std::string, tom::tools::ParameterSet::btuple>(it_best->first, it_best->second));
    }
    return std::auto_ptr<std::pair<std::string, tom::tools::ParameterSet::btuple> >(0);
}







/****************************************************************************//**
 * \brief Parse the configuration file.
 *
 * \returns Vector containing the tuples (tag, value). The order of the
 *  Vector might be relevant for the one that uses the ParameterSet. Therefore
 *  the class does not mix it up.
 *******************************************************************************/
std::vector<std::pair<std::string, std::string> > tom::tools::ParameterSet::parse_file(const std::string &filename) {

    std::vector<std::pair<std::string, std::string> > result;

    std::ifstream istr(filename.c_str());
    {
        const int errno_ = errno;
        if (!istr.good()) {
            throw std::runtime_error("Error opening config file \"" + filename + "\": " + tom::tools::strerror(errno_));
        }
    }
    std::string s, s2;
    while (std::getline(istr, s)) {
        boost::trim_left(s);
        if (s.empty() || s[0]=='#') {
            continue;
        }
        std::string::size_type idx = s.find_first_of('=');
        if (idx == std::string::npos) {
            boost::trim_right(s);
            result.push_back(std::pair<std::string, std::string>(s, ""));
        } else {
            std::string stag = s.substr(0, idx);
            boost::trim_right(stag);
            s.erase(0, idx+1);
            while (!s.empty() && s[s.length()-1]=='\\') {
                s[s.length()-1] = '\n';
                if (std::getline(istr, s2)) {
                    s += s2;
                }
            }
            result.push_back(std::pair<std::string, std::string>(stag, s));
        }
    }

    return result;
}




/****************************************************************************//**
 * \brief initialize the ParameterSet.
 *
 * \param[in] m Vector containing the tuples (tag, value). The order of the
 *  Vector might be relevant for the one that uses the ParameterSet. Therefore
 *  the class does not mix it up.
 *******************************************************************************/
std::vector<std::pair<std::string, tom::tools::ParameterSet::btuple> > tom::tools::ParameterSet::get_all(bool search_only_untouched) const {

    std::vector<std::pair<std::string, tom::tools::ParameterSet::btuple> > res;
    std::multimap<std::string, boost::tuple<std::size_t, bool, std::string> >::const_iterator it = mtag_.begin();

    for (; it!=mtag_.end(); it++) {
        if (!search_only_untouched || UNT(it->second)) {
            res.push_back(*it);
        }
    }
    return res;
}




/****************************************************************************//**
 * \brief initialize the ParameterSet.
 *
 * \param[in] m Vector containing the tuples (tag, value). The order of the
 *  Vector might be relevant for the one that uses the ParameterSet. Therefore
 *  the class does not mix it up.
 *******************************************************************************/
void tom::tools::ParameterSet::init(const std::vector<std::pair<std::string, std::string> > &m) {
    const std::size_t msize = m.size();
    std::size_t i;

    mtag_.clear();
    midx_.resize(msize);

    btuple p;
    for (i=0; i<msize; i++) {
        const std::pair<std::string, std::string> &mi = m[i];
        boost::tuples::get<0>(p) = i;
        boost::tuples::get<1>(p) = true;
        boost::tuples::get<2>(p) = mi.second;
        mtag_.insert(std::pair<std::string, btuple >(mi.first, p));
        midx_[i] = mi.first;
    }
    mtagit_ = mtag_.begin();
}





// TEMPLATE INSTANTIATIONS
template double            tom::tools::ParameterSet::pop<double            >(const std::string &tag, std::size_t first_idx, const double            *dflt);
template bool              tom::tools::ParameterSet::pop<bool              >(const std::string &tag, std::size_t first_idx, const bool              *dflt);
template unsigned int      tom::tools::ParameterSet::pop<unsigned int      >(const std::string &tag, std::size_t first_idx, const unsigned int      *dflt);
template int               tom::tools::ParameterSet::pop<int               >(const std::string &tag, std::size_t first_idx, const int               *dflt);
template unsigned long int tom::tools::ParameterSet::pop<unsigned long int >(const std::string &tag, std::size_t first_idx, const unsigned long int *dflt);
template std::string       tom::tools::ParameterSet::pop<std::string       >(const std::string &tag, std::size_t first_idx, const std::string       *dflt);







