/****************************************************************************//**
 * \file auto_vector.hpp
 * \author  Thomas Haller
 * \version 0.2
 * \date    19.11.2007
 * \brief A template vector containing pointers to objects and its destructor frees the memory.
 *
 * Copied from http://www.relisoft.com/resource/auto_vector.html.
 *******************************************************************************/
#ifndef ___INCLUDE__TOM__TOOLS__AUTO_VECTOR__
#define ___INCLUDE__TOM__TOOLS__AUTO_VECTOR__
//------------------------------------
// Reliable Software (c) 2003
// www.relisoft.com
// Any use, commercial or noncommercial of this code
// is hereby granted, under the condition
// that this copyright notice be not removed.
//------------------------------------

#include <memory>
#include <vector>
#include <cassert>
#include <algorithm>
#include <stdexcept>


//---------------------------------
// Dynamic array of owned pointers.
// Ownership transfer semantics.
//---------------------------------

namespace tom {

/****************************************************************************//**
 * \brief Contains a set of classes and functions which should make life easier :)
 *******************************************************************************/
namespace tools {


/****************************************************************************//**
 * \brief A vector whose destructor destroys the heap objects it contains.
 *
 * Taken from http://www.relisoft.com/resource/auto_vector.html
 *******************************************************************************/
template <class T>
class auto_vector {


public:
    class auto_lvalue {
    public:
        auto_lvalue(T * &p): _p(p) {};
        operator T *()  const { return _p; }
        T *operator->() const { return _p; }
        auto_lvalue &operator=(std::auto_ptr<T> ap) {
            delete this->_p;
            this->_p = ap.release();
            return *this;
        }
    private:
        T * & _p;
    };
public:
    typedef typename std::vector<T*>::size_type size_type;

    explicit auto_vector(size_type capacity=0);
    ~auto_vector();

    // memory management
    size_type size() const     { return this->_arr.size (); }
    size_type capacity() const { return this->_arr.capacity(); }
    void reserve(size_type count);
    void resize(size_type newSize);
    void erase(size_type idx);
    void clear();
    void compact ();
    void swap (auto_vector<T> & other) { this->_arr.swap (other._arr); }


    // array access
    T const *operator[](size_type i) const  { return this->_arr[i]; }
    auto_lvalue operator[](size_type i)     { return auto_lvalue(this->_arr[i]); }
    T const *at(size_type i) const          { return this->_arr.at(i); }
    auto_lvalue at(size_type i)             { return auto_lvalue(this->_arr.at(i)); }
    void assign(size_type i, std::auto_ptr<T> p);
    void assign_direct(size_type i, T * p);
    void insert(size_type idx, std::auto_ptr<T> p);

    // stack access
    void push_back(std::auto_ptr<T> p);
    std::auto_ptr<T> pop_back(); // non-standard
    T *back()               { return _arr.back(); }
    T const *back () const  { return _arr.back(); }
    T *front()              { return _arr.front(); }
    T const *front() const  { return _arr.front(); }

    // iterators
    typedef typename std::vector<T*>::iterator                      iterator;
    typedef typename std::vector<T*>::const_iterator                const_iterator;
    typedef typename std::vector<T*>::reverse_iterator              reverse_iterator;
    typedef typename std::vector<T*>::const_reverse_iterator        const_reverse_iterator;


    iterator                begin()         { return _arr.begin();  }
    iterator                end()           { return _arr.end();    }
    const_iterator          begin()  const  { return _arr.begin();  }
    const_iterator          end()    const  { return _arr.end();    }
    reverse_iterator        rbegin()        { return _arr.rbegin(); }
    reverse_iterator        rend()          { return _arr.rend();   }
    const_reverse_iterator  rbegin() const  { return _arr.rbegin(); }
    const_reverse_iterator  rend()   const  { return _arr.rend();   }

    iterator erase(iterator it);

    // iterator/index conversion
    size_type ToIndex(iterator const &it);
    size_type ToIndex(reverse_iterator const &rit);
    iterator ToIter(size_type idx);
    reverse_iterator ToRIter(size_type idx);

    const std::vector<T *> &get_std_vector() { return this->_arr; }

private:
    std::vector<T*> _arr;
};




template <class T>
inline auto_vector<T>::auto_vector(size_type capacity): _arr(capacity) {
}

template <class T>
inline auto_vector<T>::~auto_vector() {
    this->clear();
}

template <class T>
inline void auto_vector<T>::push_back(std::auto_ptr<T> ptr) {
    this->_arr.push_back(ptr.release());
}

template <class T>
inline std::auto_ptr<T> auto_vector<T>::pop_back() {
    if (this->_arr.empty()) {
        throw std::out_of_range();
    }
    T *p = this->_arr.back();
    this->_arr.pop_back();
    return std::auto_ptr<T>(p);
}

template <class T>
class DeletePtr {
public:
    void operator () (T * p) {
        delete p;
    }
};
template <class T>
void auto_vector<T>::clear() {
    std::for_each(this->begin (), this->end (), DeletePtr<T> ());
    this->_arr.clear ();
}


template <class T>
inline void auto_vector<T>::assign_direct(size_type i, T *p) {
    if (_arr.at(i) == p) {
        return;
    }
    delete _arr[i];
    _arr[i] = p;
}

template <class T>
inline void auto_vector<T>::assign(size_type i, std::auto_ptr<T> p) {
    if (_arr.at(i) != p.get ()) {
        delete _arr[i];
    }
    _arr[i] = p.release();
}

template <class T>
void auto_vector<T>::erase(size_type idx) {
    delete this->_arr.at(idx);
    this->_arr.erase (ToIter(idx));
}

template <class T>
typename auto_vector<T>::iterator auto_vector<T>::erase(typename auto_vector<T>::iterator it) {
    if (!(it < this->end())) {
        throw std::out_of_range();
    }
    delete *it;
    return this->_arr.erase(it);
}

template <class T>
void auto_vector<T>::compact () {
    T *null = 0;
    iterator pos = std::remove(this->begin(), this->end(), null);
    this->_arr.resize(pos - this->begin());
}

template <class T>
typename auto_vector<T>::size_type auto_vector<T>::ToIndex(iterator const &it) {
    if (!(it - begin () >= 0)) {
        throw std::out_of_range();
    }
    return static_cast<size_type> (it - begin ());
}

template <class T>
typename auto_vector<T>::size_type auto_vector<T>::ToIndex(reverse_iterator const &rit) {
    iterator it = rit.base();
    --it;
    if (!(it - begin () >= 0)) {
        throw std::out_of_range();
    }
    return static_cast<size_type> (it - begin ());
}

template <class T>
typename auto_vector<T>::iterator auto_vector<T>::ToIter (size_type idx) {
    return this->begin() + idx;
}

template <class T>
typename auto_vector<T>::reverse_iterator auto_vector<T>::ToRIter (size_type idx) {
    ++idx;
    return reverse_iterator (ToIter (idx));
}


template <class T>
inline void auto_vector <T>::reserve (size_type count) {
    _arr.reserve (count);
}

template <class T>
inline void auto_vector<T>::resize (size_type newSize) {
    if (newSize < size ()) {
        std::for_each (ToIter (newSize), end (), DeletePtr<T> ());
    }
    _arr.resize (newSize);
}

template <class T>
void auto_vector<T>::insert (size_type idx, std::auto_ptr<T> p) {
    _arr.insert(this->ToIter(idx), p.get());
    p.release ();
}


} // namespace tools
} // namespace tom

#endif


