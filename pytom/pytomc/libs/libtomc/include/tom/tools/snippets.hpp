/****************************************************************************//**
 * \file snippets.hpp
 * \brief Small functions for general purpose.
 * \author  Thomas Haller
 * \version 0.2
 * \date    25.01.2008
 *
 * Contains some useful(?) functions that have not proper place on their own.
 *******************************************************************************/
#ifndef ___INCLUDE__TOM__TOOLS__SNIPPETS_HPP__
#define ___INCLUDE__TOM__TOOLS__SNIPPETS_HPP__


#include <vector>
#include <string>
#include <iosfwd>

#include <cmath>

#include <boost/shared_ptr.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/static_assert.hpp>






namespace tom {
namespace tools {




/// Exception class with an int parameter for exit
class exit_status {
public:
    exit_status(int i)
        :   exit_status_(i) {
    }
    operator int() const {
        return exit_status_;
    }
private:
    int exit_status_;
};


template<typename T> std::size_t decimal_length(T n);


template<typename T>
bool str2range(const std::string &instr, std::vector<T> &range, char delimiter, char delimiter_range);


std::string now2str(const std::string &mode="");


template<typename T, typename _Compare>
void unify_shared_vector(std::vector<boost::shared_ptr<T> > &v, const _Compare &c = _Compare());



template<typename T, typename TIDX>
void subset(const std::vector<T> &v, const std::vector<TIDX> &vidx, std::vector<T> &vsub, bool unique, bool sort);


inline std::vector<const char *> svector2c(std::vector<std::string> &v) {
    std::vector<const char *> r(v.size()+1);
    std::vector<const char *>::iterator rit = r.begin();
    std::vector<std::string>::iterator vit = v.begin();
    for (; vit!=v.end(); vit++, rit++) {
        *rit = vit->c_str();
    }
    *rit = NULL;
    return r;
}


std::string strerror(int errno_);
std::string gethostname();



std::string split_line(std::string &s, std::size_t max_line_length);
void print_parameters(std::ostream &clog, const std::string &name, const std::string &text, const std::string &prefix, std::size_t max_line_length);


std::string prefix_lines(const std::string &s, const std::string &prefix);

std::string getcwd();

void str_trim_spaces(std::string &str);


void nanosleep(double sec);


std::string align_decimal(double x, std::size_t num_int, std::size_t num_prec);



} // namespace tools
} // namespace tom




// INLINE FUNCIONS



/***************************************************************************//**
 * Get the number of characters to write the integral part of the number \c n
 * in decimal notation. The character for a sign is not included.
 ******************************************************************************/
template<typename T>
inline std::size_t tom::tools::decimal_length(T n) {
    if (n < 0) {
        n = -n;
    }
    if (n < 10) {
        return 1;
    }
    return static_cast<std::size_t>(::ceil(std::log(n) / 2.302585092994045901093613792909309268)); /* <= log(10) */
}





#endif


