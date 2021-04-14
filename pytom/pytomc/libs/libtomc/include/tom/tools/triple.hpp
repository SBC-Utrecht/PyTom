/****************************************************************************//**
 * \file triple.hpp
 * \brief Class for a triple similar to std::pair.
 * \author  Thomas Haller
 * \version 0.2
 * \date    30.01.2008
 *
 * Contains the class tom::tools::triple which provides similar funtionality
 * as std::pair for three elements.
 *
 * All functions are inline (there is no cpp-file).
 *******************************************************************************/
#ifndef ___INCLUDE__TOM__TOOLS__TRIPLE_HPP__
#define ___INCLUDE__TOM__TOOLS__TRIPLE_HPP__







namespace tom {
namespace tools {



/****************************************************************************//**
 * Template class to store three values of different data types.
 *******************************************************************************/
template<class _T1, class _T2, class _T3>
struct triple {
    typedef _T1 first_type;
    typedef _T2 second_type;
    typedef _T3 third_type;

    union {
        _T1 first;
        _T1 x;
    };
    union {
        _T2 second;
        _T2 y;
    };
    union {
        _T3 third;
        _T3 z;
    };

    triple(): first(), second(), third() { }

    triple(const _T1 &__a, const _T2 &__b, const _T3 &__c)
        : first(__a), second(__b), third(__c) { }

    template<class _U1, class _U2, class _U3> triple(const triple<_U1, _U2, _U3> &__p)
        : first(__p.first), second(__p.second), third(__p.third) { }
}; // class tom::tools::triple



template<class _T1, class _T2, class _T3>
inline bool operator==(const triple<_T1, _T2, _T3> &__x, const triple<_T1, _T2, _T3> &__y) {
    return __x.first==__y.first && __x.second==__y.second && __x.third==__y.third;
}


template<class _T1, class _T2, class _T3>
inline bool operator<(const triple<_T1, _T2, _T3> &__x, const triple<_T1, _T2, _T3> &__y) {
    return __x.first<__y.first || (!(__y.first<__x.first) && (__x.second<__y.second || (!(__y.second<__x.second) && __x.third<__y.third)));
}


template<class _T1, class _T2, class _T3>
inline bool operator!=(const triple<_T1, _T2, _T3> &__x, const triple<_T1, _T2, _T3> &__y) {
    return !(__x == __y);
}


template<class _T1, class _T2, class _T3>
inline bool operator>(const triple<_T1, _T2, _T3> &__x, const triple<_T1, _T2, _T3> &__y) {
    return __y < __x;
}


template<class _T1, class _T2, class _T3>
inline bool operator<=(const triple<_T1, _T2, _T3> &__x, const triple<_T1, _T2, _T3> &__y) {
    return !(__y < __x);
}


template<class _T1, class _T2, class _T3>
inline bool operator>=(const triple<_T1, _T2, _T3> &__x, const triple<_T1, _T2, _T3> &__y) {
    return !(__x < __y);
}


template<class _T1, class _T2, class _T3>
inline triple<_T1, _T2, _T3> make_pair(_T1 __x, _T2 __y, _T3 __z) {
    return triple<_T1, _T2, _T3>(__x, __y, __z);
}



} // namespace tools
} // namespace tom




#endif


