/****************************************************************************//**
 * \file ParameterSet.hpp
 * \brief
 * \author  Thomas Haller
 * \version 0.2
 * \date    21.11.2008
 *
 * Contains the class tom::tools::ParameterSet to store the program parameters
 * as read from the configuration file.
 * Each parameter is a pair of strings containing the key (tag) and its value.
 * Parameters with the same name can occure several times and for each parameter
 * also the index in the set is relevant.
 *******************************************************************************/
#ifndef ___INCLUDE__TOM__TOOLS__PARAMETERSET_HPP__
#define ___INCLUDE__TOM__TOOLS__PARAMETERSET_HPP__


#include <vector>
#include <map>
#include <string>
#include <utility>
#include <memory>


#include <boost/tuple/tuple.hpp>



namespace tom {
namespace tools {



class ParameterSet {

public:
    typedef boost::tuple<std::size_t, bool, std::string> btuple;

    ParameterSet(const std::string &filename);
    ParameterSet(const std::vector<std::pair<std::string, std::string> > &m);

    /// Gives the zero based INDex of the tuple.
    static std::size_t       &IND(      btuple &b) { return boost::get<0>(b); }
    static const std::size_t &IND(const btuple &b) { return boost::get<0>(b); }

    /// Determines whether the tuple is UNTouched or not.
    static bool              &UNT(      btuple &b) { return boost::get<1>(b); }
    static const bool        &UNT(const btuple &b) { return boost::get<1>(b); }

    /// Gives the VALue of the entry as string.
    static       std::string &VAL(      btuple &b) { return boost::get<2>(b); }
    static const std::string &VAL(const btuple &b) { return boost::get<2>(b); }

    template<typename T> T pop(const std::string &tag, std::size_t first_idx, const T *dflt);
    template<typename T> T get(const std::string &tag, bool search_only_untouched, std::size_t first_idx, const T *dflt);

    std::auto_ptr<tom::tools::ParameterSet::btuple> pop(const std::string &tag, std::size_t first_idx);
    std::auto_ptr<tom::tools::ParameterSet::btuple> get(const std::string &tag, bool search_only_untouched, std::size_t first_idx) const;

    std::auto_ptr<std::pair<std::string, tom::tools::ParameterSet::btuple> > pop(std::size_t first_idx);
    std::auto_ptr<std::pair<std::string, tom::tools::ParameterSet::btuple> > get(bool search_only_untouched, std::size_t first_idx) const;

    static std::vector<std::pair<std::string, std::string> > parse_file(const std::string &filename);

    std::vector<std::pair<std::string, tom::tools::ParameterSet::btuple> > get_all(bool search_only_untouched) const;

private:
    void init(const std::vector<std::pair<std::string, std::string> > &m);

    std::string filename_;

    /// Multimap with all tags. For each entry, its index, whether it was references and the string is saved.
    std::multimap<std::string, boost::tuple<std::size_t, bool, std::string> > mtag_;

    std::multimap<std::string, boost::tuple<std::size_t, bool, std::string> >::iterator mtagit_;


    /// Vector with all the tags, ordered by their index.
    std::vector<std::string> midx_;

}; // class tom::tools::ParameterSet


} // namespace tools
} // namespace tom















#endif
