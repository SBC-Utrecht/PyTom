/***************************************************************************//**
 * \file Precision.hpp
 * \brief Enum class for single/double precision.
 * \author  Thomas Haller
 * \version 0.1
 * \date    29.01.2009
 *
 *
 ******************************************************************************/
#ifndef ___INCLUDE__TOM__TOOLS__PRECISION_HPP__
#define ___INCLUDE__TOM__TOOLS__PRECISION_HPP__



#include <stdexcept>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <limits>


namespace tom {
namespace tools {


/***************************************************************************//**
 *
 ******************************************************************************/
struct Precision {
    typedef enum { FLOAT, DOUBLE } Ptype;

    Precision(Ptype p) { assert(p==FLOAT || p==DOUBLE); prec_ = p; }
    Precision(): prec_(DOUBLE) { }
    Precision(const std::string &s): prec_(str2prec(s)) { }

    operator std::string() const { return prec2str(prec_);  }
    operator Ptype      () const { return prec_;            }
    operator int        () const { return prec_;            }
    Precision &operator=(Ptype             &p) { assert(p==FLOAT || p==DOUBLE); prec_ = p; return *this; };
    Precision &operator=(const std::string &s) { prec_ = str2prec(s); return *this; };
    Precision &operator=(int                p) {
        if (p == static_cast<int>(FLOAT)) {
            prec_ = FLOAT;
        } else {
            assert(p == static_cast<int>(DOUBLE));
            prec_ = DOUBLE;
        }
        return *this;
    }

    static std::string prec2str(Ptype p);
    static Ptype str2prec(const std::string &s);

    bool is_float () const { return prec_ == FLOAT ; }
    bool is_double() const { return prec_ == DOUBLE; }
    void set_float () { prec_ = FLOAT ; }
    void set_double() { prec_ = DOUBLE; }
    void set(Ptype prec) { prec_ = (prec==FLOAT ? FLOAT : DOUBLE); }
    void set(const std::string &s) { prec_ = str2prec(s); }

    bool operator==(Precision &c) const { return (is_float()&&c.is_float()) || (is_double()&&c.is_double()); }

private:
    Ptype prec_;
};

} // namespace tools
} // namespace tom


template <typename T>
inline bool are_same(T a, T b) {
    return std::fabs(a - b) < std::numeric_limits<T>::epsilon();
}


/****************************************************************************//**
 * \brief Convert the precision to a string.
 *******************************************************************************/
inline std::string tom::tools::Precision::prec2str(Ptype p) {
    if (p == FLOAT) {
        return "float";
    } else {
        return "double";
    }
    assert(!"invalid precision type");
}


/****************************************************************************//**
 * \brief Convert the precision to a string.
 *******************************************************************************/
inline tom::tools::Precision::Ptype tom::tools::Precision::str2prec(const std::string &s1) {
    std::string s(s1);
    std::transform(s.begin(), s.end(), s.begin(), &tolower);

    if (s == "single" || s == "float") { return FLOAT ; }
    if (s == "double"                ) { return DOUBLE; }
    throw std::invalid_argument("Invalid value \"" + s1 + "\" for precision. Valid are 'float' ('single') and 'double'.");

}



#endif

