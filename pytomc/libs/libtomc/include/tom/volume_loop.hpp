/***************************************************************************//**
 * \file volume_loop.hpp
 * \brief
 * \author  Thomas Haller
 * \version 0.2
 * \date    17.01.2008
 *
 * Contains template-functions for traversing a tom::volume-object and perform
 * operations on the elements. These functions are not intended to be called
 * directly, instead they help to implement the recuring work of traversing though
 * a volume (considering the different strides), as for example in \c tom::Volume<float>::setValues
 ******************************************************************************/
#ifndef INCLUDE_CORE_VOLUME_LOOP_HPP_
#define INCLUDE_CORE_VOLUME_LOOP_HPP_



// used for std::ptrdiff_t
#include <cstddef>
#include <cassert>

#include <boost/type_traits/is_const.hpp>
#include <boost/static_assert.hpp>


namespace tom {


/// \brief contains main routines for for_each.
namespace loop {

typedef std::size_t index_type;

struct packed_idx {
    index_type x;
    index_type y;
    index_type z;
};


template<typename T> void ptr_add_byte_offset(T *&ptr, std::ptrdiff_t offset);


template<typename TVOL, typename TOP> void for_each           (TVOL &v1, TOP op);
template<typename TVOL, typename TOP> void for_each_no_idx    (TVOL &v1, TOP op);
template<typename TVOL, typename TOP> void for_each_packed_idx(TVOL &v1, TOP op);
template<typename TVOL, typename TOP> void for_each_while     (TVOL &v1, TOP op);
template<typename TVOL, typename TOP> void for_each_step      (TVOL &v1, TOP op);
template<typename TVOL, typename TOP> void for_each_general   (TVOL &v1, TOP op);


template<typename TVOL1, typename TVOL2, typename TOP> void for_each           (TVOL1 &v1, TVOL2 &v2, TOP op);
template<typename TVOL1, typename TVOL2, typename TOP> void for_each_no_idx    (TVOL1 &v1, TVOL2 &v2, TOP op);
template<typename TVOL1, typename TVOL2, typename TOP> void for_each_packed_idx(TVOL1 &v1, TVOL2 &v2, TOP op);
template<typename TVOL1, typename TVOL2, typename TOP> void for_each_while     (TVOL1 &v1, TVOL2 &v2, TOP op);
template<typename TVOL1, typename TVOL2, typename TOP> void for_each_step      (TVOL1 &v1, TVOL2 &v2, TOP op);
template<typename TVOL1, typename TVOL2, typename TOP> void for_each_general   (TVOL1 &v1, TVOL2 &v2, TOP op);

template<typename TVOL1, typename TVOL2, typename TVOL3, typename TOP> void for_each           (TVOL1 &v1, TVOL2 &v2, TVOL3 &v3, TOP op);
template<typename TVOL1, typename TVOL2, typename TVOL3, typename TOP> void for_each_no_idx    (TVOL1 &v1, TVOL2 &v2, TVOL3 &v3, TOP op);
template<typename TVOL1, typename TVOL2, typename TVOL3, typename TOP> void for_each_packed_idx(TVOL1 &v1, TVOL2 &v2, TVOL3 &v3, TOP op);

template<typename T, bool ISCONST> struct const_type;



/***************************************************************************//**
 * contains a typedef type which corresponds to the data-type which the tom::Volume
 * can hold (with proper const/non-const).
 ******************************************************************************/
template<typename T>
struct resolve_volume_type {
    typedef typename boost::remove_reference<T>::type vol_type;
    typedef typename const_type<    typename vol_type::element_type,
                                    boost::is_const<T>::value
                                >::type data_type;
};




}
}





/******************************************************************************/
/* INLINE FUNCTIONS ***********************************************************/
/******************************************************************************/





namespace tom {
namespace loop {
/***************************************************************************//**
 * \brief A template class with a single typedef \c type to create a const datatype.
 *
 * Depending on the second template parameter \c ISCONST, the typedef \c type
 * names either a const or not const version of the first template type.
 *
 * Uses partial template specialization to achieve that.
 ******************************************************************************/
template<typename T, bool ISCONST> struct const_type          { typedef       T type; };
template<typename T>               struct const_type<T, true> { typedef const T type; };
} // namespace loop
} // namespace tom







/***************************************************************************//**
 * \brief Add a byte offset to a pointer of a template-type.
 *
 * \param[in,out] ptr Reference to pointer to which the offset is added.
 * \param[in] offset Offset in bytes (can be negative).
 ******************************************************************************/
template<typename T>
inline void tom::loop::ptr_add_byte_offset(T *&ptr, std::ptrdiff_t offset) {
    ptr = reinterpret_cast<T *>(const_cast<char *>((reinterpret_cast<const char *>(ptr) + offset)));
}







namespace tom { namespace loop { namespace functor { namespace {
// Structure for wrapping the functor to another signature of the function call.
template<typename T, typename TOP>
struct for_each_wrapper_1 {
    for_each_wrapper_1(TOP &op): op_(op) { }
    void stepz(const tom::loop::index_type &z) { }
    void stepy(const tom::loop::index_type &y) { }
    bool operator()(T &val, const tom::loop::packed_idx &idx) {
        op_(val, idx.x, idx.y, idx.z);
        return true;
    }
    TOP &op_;
};
} // namespace
} // namespace functor
} // namespace loop
} // namespace tom

/***************************************************************************//**
 * \brief Template function for traversing a volume and execute an operation
 *
 * \param[in,out] v1 A reference to a volume. Each element is traversed..
 * \param[in,out] op A template which can be called as a function
 *      with 4 parameters: op(T val, index_type x, index_type y, index_type z).
 *      The type of \c val is element type of the volume \c v1.
 *
 * For each element \c val of the volume \c v1 the operation \c op is called with
 * the index of the elements \c x, \c y, \c z. This is usefull for example when
 * looking for the maximum (the parameter \c val actually being a constant reference)
 * or when setting each element to a certain value (\c val being a mutable reference).\\
 * This function considers the stride parameters of the volume.
 *
 * Warning: if you would like to get a result from the functor, you have to
 * explicitly state the type \c TOP to be reference. Otherwise a copy is passed
 * as that is how C++ inferes the template type.
 *
 * The functions \c for_each_no_idx and \c for_each_packed_idx are similar and
 * were introduced, because I was unable ( :-( ) to get boost::lambda swallaw
 * four parameters. I think in the current implemtation, the lambda functors
 * only support up to 3 parameters.
 ******************************************************************************/
template<typename TVOL, typename TOP>
inline void tom::loop::for_each(TVOL &v1, TOP op) {
    for_each_general(v1,
        functor::for_each_wrapper_1<    typename resolve_volume_type<TVOL>::data_type,
                                        typename boost::remove_reference<TOP>::type
        >(op));
}





namespace tom { namespace loop { namespace functor { namespace {
// Structure for wrapping the functor to another signature of the function call.
template<typename T, typename TOP>
struct for_each_wrapper_no_idx_1 {
    for_each_wrapper_no_idx_1(TOP &op): op_(op) { }
    void stepz(const tom::loop::index_type &) { }
    void stepy(const tom::loop::index_type &) { }
    bool operator()(T &val, const tom::loop::packed_idx &) {
        op_(val);
        return true;
    }
    TOP &op_;
};
} } } } // namespace tom::loop::functor::unnamed
/***************************************************************************//**
 * analog to <tt>void tom::loop::for_each(TVOL &v1, TOP op)</tt>, except that
 * the operator can be called in the form op(T val) (i.e. without the index of
 * the current element.
 ******************************************************************************/
template<typename TVOL, typename TOP>
inline void tom::loop::for_each_no_idx(TVOL &v1, TOP op) {
    for_each_general(v1,
        functor::for_each_wrapper_no_idx_1<
            typename resolve_volume_type<TVOL>::data_type,
            typename boost::remove_reference<TOP>::type
        >(op));
}





namespace tom { namespace loop { namespace functor { namespace {
// Structure for wrapping the functor to another signature of the function call.
template<typename T, typename TOP>
struct for_each_wrapper_packed_idx_1 {
    for_each_wrapper_packed_idx_1(TOP &op): op_(op) { }
    void stepz(const tom::loop::index_type &z) { }
    void stepy(const tom::loop::index_type &y) { }
    bool operator()(T &val, const tom::loop::packed_idx &idx) {
        op_(val, idx);
        return true;
    }
    TOP &op_;
};
} } } } // namespace tom::loop::functor::unnamed
/****************************************************************************//**
 * Like <tt>void tom::loop::for_each(TVOL &v1, TOP op)</tt>, but the operator
 * needs to provide the <tt>operator()(T, packed_idx)</tt> instead of
 * <tt>operator()(T, index_type, index_type, index_type)</tt>
 *
 * This is usefull because (currently) boost::lamda does not support AFAIK
 * that much arguments for the functor.
 *******************************************************************************/
template<typename TVOL, typename TOP>
inline void tom::loop::for_each_packed_idx(TVOL &v1, TOP op) {
    for_each_general(v1,
        functor::for_each_wrapper_packed_idx_1<
            typename resolve_volume_type<TVOL>::data_type,
            typename boost::remove_reference<TOP>::type
        >(op));
}






namespace tom { namespace loop { namespace functor { namespace {
// Structure for wrapping the functor to another signature of the function call.
template<typename T, typename TOP>
struct for_each_wrapper_while_1 {
    for_each_wrapper_while_1(TOP &op): op_(op) { }
    void stepz(const tom::loop::index_type &z) { }
    void stepy(const tom::loop::index_type &y) { }
    bool operator()(T &val, const tom::loop::packed_idx &idx) {
        return op_(val, idx.x, idx.y, idx.z);
    }
    TOP &op_;
};
} } } } // namespace tom::loop::functor::unnamed
/***************************************************************************//**
 * \brief Traversing the volume and perform an operation.
 *
 * Similar to <tt>void tom::loop::for_each(TVOL &v1, TOP op)</tt>.
 * The difference is that the operation \c op returns a value
 * which can evaluate to true or false. \c false, breaks the loop over the volume.
 * Usefull for example when comparing an entire volume with one number. When the
 * first mismatch is found, the search can be terminated with a negative result.
 ******************************************************************************/
template<typename TVOL, typename TOP>
inline void tom::loop::for_each_while(TVOL &v1, TOP op) {
    for_each_general(v1,
        functor::for_each_wrapper_while_1<
            typename resolve_volume_type<TVOL>::data_type,
            typename boost::remove_reference<TOP>::type
        >(op));
}




namespace tom { namespace loop { namespace functor { namespace {
// Structure for wrapping the functor to another signature of the function call.
template<typename T, typename TOP>
struct for_each_wrapper_step_1 {
    for_each_wrapper_step_1(TOP &op): op_(op) { }
    void stepz(const tom::loop::index_type &z) {
        op_.stepz(z);
    }
    void stepy(const tom::loop::index_type &y) {
        op_.stepy(y);
    }
    bool operator()(T &val, const tom::loop::packed_idx &idx) {
        op_(val, idx.x);
        return true;
    }
    TOP &op_;
};
} } } } // namespace tom::loop::functor::unnamed
/***************************************************************************//**
 * \brief Traversing the volume and perform an operation.
 *
 * Similar to <tt>void tom::loop::for_each(TVOL &v1, TOP op)</tt>.
 * The differenc is the definition of the operation \c op.
 * \c TOP must now be a functor that supports
 * the methods <tt>stepz(index_type z)</tt>, <tt>stepy(index_type y)</tt> and
 * <tt>operator()(T &, std::size_t x)</c>. In this case \c stepy (\c stepz)
 * are called as the loop reaches a new line (level) while traversing the volume.
 * The object \c op should now remember the current position and can thereby save
 * unnecessary computations.
 ******************************************************************************/
template<typename TVOL, typename TOP>
inline void tom::loop::for_each_step(TVOL &v1, TOP op) {
    for_each_general(v1,
        functor::for_each_wrapper_step_1<   typename resolve_volume_type<TVOL>::data_type,
                                            typename boost::remove_reference<TOP>::type
        >(op));
}




/****************************************************************************//**
 * Like tom::loop::for_each, but the operator needs to provide the
 * operator()(T, packed_idx) instead of operator()(T, index_type, index_type, index_type)
 *
 * This is usefull because (currently) boost::lamda does not support AFAIK
 * that much arguments for the functor.
 *******************************************************************************/
template<typename TVOL, typename TOP>
inline void tom::loop::for_each_general(TVOL &v1, TOP op) {

    typedef typename resolve_volume_type<TVOL>::data_type TPTR1;

    TPTR1 *ptr1 = &v1.get();

    assert( v1.getSizeX()==static_cast<index_type>(v1.getSizeX()) &&
            v1.getSizeY()==static_cast<index_type>(v1.getSizeY()) &&
            v1.getSizeZ()==static_cast<index_type>(v1.getSizeZ()) );

    const index_type sizex = v1.getSizeX();
    const index_type sizey = v1.getSizeY();
    const index_type sizez = v1.getSizeZ();
    std::size_t stridex1 = v1.getStrideX();
    std::size_t stridey1 = v1.getStrideY();
    std::size_t stridez1 = v1.getStrideZ();
    packed_idx idx;
    assert(sizex && sizey && sizez && stridex1 && stridey1 && stridez1);

    // are the volume-strides (in bytes) a multiple of the data-type size?
    if (!(stridex1%sizeof(TPTR1)) &&
        !(stridey1%sizeof(TPTR1)) &&
        !(stridez1%sizeof(TPTR1))) {

        // convert the stride parameters from "bytes" to "number of elements"
        stridex1 /= sizeof(TPTR1);
        stridey1 /= sizeof(TPTR1);
        stridez1 /= sizeof(TPTR1);
        assert(stridex1 && stridey1 && stridez1);

        // Are the elements in the fastest dimension without gap?
        if (stridex1 == 1) {
            if (stridey1 == sizex && stridez1 == sizex*sizey) {
                for (idx.z=0; idx.z<sizez; idx.z++) {
                    op.stepz(idx.z);
                    for (idx.y=0; idx.y<sizey; idx.y++) {
                        op.stepy(idx.y);
                        for (idx.x=0; idx.x<sizex; idx.x++) {
                            if (!op(*ptr1++, const_cast<const packed_idx &>(idx))) {
                                return;
                            }
                        }
                    }
                }
            } else {
                stridez1 -= sizey*stridey1;
                for (idx.z=0; idx.z<sizez; idx.z++) {
                    op.stepz(idx.z);
                    for (idx.y=0; idx.y<sizey; idx.y++) {
                        op.stepy(idx.y);
                        for (idx.x=0; idx.x<sizex; idx.x++) {
                            if (!op(ptr1[idx.x], idx)) {
                                return;
                            }
                        }
                        ptr1 += stridey1;
                    }
                    ptr1 += stridez1;
                }
            }
        } else {
            stridez1 -= sizey*stridey1;
            stridey1 -= sizex*stridex1;
            for (idx.z=0; idx.z<sizez; idx.z++) {
                op.stepz(idx.z);
                for (idx.y=0; idx.y<sizey; idx.y++) {
                    op.stepy(idx.y);
                    for (idx.x=0; idx.x<sizex; idx.x++) {
                        if (!op(ptr1[idx.x], idx)) {
                            return;
                        }
                        ptr1 += stridex1;
                    }
                    ptr1 += stridey1;
                }
                ptr1 += stridez1;
            }
        }
    } else {
        stridez1 -= sizey*stridey1;
        stridey1 -= sizex*stridex1;

        for (idx.z=0; idx.z<sizez; idx.z++) {
            op.stepz(idx.z);
            for (idx.y=0; idx.y<sizey; idx.y++) {
                op.stepy(idx.y);
                for (idx.x=0; idx.x<sizex; idx.x++) {
                    if (!op(ptr1[idx.x], idx)) {
                        return;
                    }
                    ptr_add_byte_offset(ptr1, stridex1);
                }
                ptr_add_byte_offset(ptr1, stridey1);
            }
            ptr_add_byte_offset(ptr1, stridez1);
        }
    }
}






namespace tom { namespace loop { namespace functor { namespace {
// Structure for wrapping the functor to another signature of the function call.
template<typename T1, typename T2, typename TOP>
struct for_each_wrapper_2 {
    for_each_wrapper_2(TOP &op): op_(op) { }
    void stepz(const tom::loop::index_type &z) { }
    void stepy(const tom::loop::index_type &y) { }
    bool operator()(T1 &val1, T2 &val2, const tom::loop::packed_idx &idx) {
        op_(val1, val2, idx.x, idx.y, idx.z);
        return true;
    }
    TOP &op_;
};
} } } } // namespace tom::loop::functor::unnamed
/***************************************************************************//**
 * \brief Traversing two volumes.
 *
 * While in \c tom::loop::for_each only one volume is traversed, now two volumes
 * of the same size (but maybe different types) are traversed elementwise.
 * Usefull for example to copy the content of one volume to the second one
 * (e.g. tom::Volume<T>::setValues(const tom::Volume<T> &v1)).
 *
 * Warning: if you would like to get a result from the functor, you have to
 * explicitly state the type \c TOP as reference. Otherwise a copy is passed.
 * Or use the return value of the loop.
 ******************************************************************************/
template<typename TVOL1, typename TVOL2, typename TOP>
inline void tom::loop::for_each(TVOL1 &v1, TVOL2 &v2, TOP op) {
    typedef functor::for_each_wrapper_2<typename resolve_volume_type<TVOL1>::data_type,
                                        typename resolve_volume_type<TVOL2>::data_type,
                                        typename boost::remove_reference<TOP>::type
                                        > wrapper_type;
    for_each_general(v1, v2, wrapper_type(op));
}





namespace tom { namespace loop { namespace functor { namespace {
// Structure for wrapping the functor to another signature of the function call.
template<typename T1, typename T2, typename TOP>
struct for_each_wrapper_no_idx_2 {
    for_each_wrapper_no_idx_2(TOP &op): op_(op) { }
    void stepz(const tom::loop::index_type &) { }
    void stepy(const tom::loop::index_type &) { }
    bool operator()(T1 &val1, T2 &val2, const tom::loop::packed_idx &) {
        op_(val1, val2);
        return true;
    }
    TOP &op_;
};
} } } } // namespace tom::loop::functor::unnamed
/***************************************************************************//**
 * \brief Traversing two volumes.
 *
 * While in \c tom::loop::for_each only one volume is traversed, now two volumes
 * of the same size (but maybe different types) are traversed elementwise.
 * Usefull for example to copy the content of one volume to the second one
 * (e.g. tom::Volume<T>::setValues(const tom::Volume<T> &v1)).
 *
 * Warning: if you would like to get a result from the functor, you have to
 * explicitly state the type \c TOP as reference. Otherwise a copy is passed.
 * Or use the return value of the loop.
 ******************************************************************************/
template<typename TVOL1, typename TVOL2, typename TOP>
inline void tom::loop::for_each_no_idx(TVOL1 &v1, TVOL2 &v2, TOP op) {
    for_each_general(v1, v2,
        functor::for_each_wrapper_no_idx_2< typename resolve_volume_type<TVOL1>::data_type,
                                            typename resolve_volume_type<TVOL2>::data_type,
                                            typename boost::remove_reference<TOP>::type
        >(op));
}


namespace tom { namespace loop { namespace functor { namespace {
// Structure for wrapping the functor to another signature of the function call.
template<typename T1, typename T2, typename TOP>
struct for_each_wrapper_packed_idx_2 {
    for_each_wrapper_packed_idx_2(TOP &op): op_(op) { }
    void stepz(const tom::loop::index_type &z) { }
    void stepy(const tom::loop::index_type &y) { }
    bool operator()(T1 &val1, T2 &val2, const tom::loop::packed_idx &idx) {
        op_(val1, val2, idx);
        return true;
    }
    TOP &op_;
};
} } } } // namespace tom::loop::functor::unnamed
/****************************************************************************//**
 * \brief Traversing two volumes.
 *
 * While in \c tom::loop::for_each only one volume is traversed, now two volumes
 * of the same size (but maybe different types) are traversed elementwise.
 * Usefull for example to copy the content of one volume to the second one
 * (e.g. tom::Volume<T>::setValues(const tom::Volume<T> &v1)).
 *
 * Warning: if you would like to get a result from the functor, you have to
 * explicitly state the type \c TOP as reference. Otherwise a copy is passed.
 * Or use the return value of the loop.
 *******************************************************************************/
template<typename TVOL1, typename TVOL2, typename TOP>
inline void tom::loop::for_each_packed_idx(TVOL1 &v1, TVOL2 &v2, TOP op) {
    for_each_general(v1, v2,
        functor::for_each_wrapper_packed_idx_2< typename resolve_volume_type<TVOL1>::data_type,
                                                typename resolve_volume_type<TVOL2>::data_type,
                                                typename boost::remove_reference<TOP>::type
        >(op));
}





namespace tom { namespace loop { namespace functor { namespace {
// Structure for wrapping the functor to another signature of the function call.
template<typename T1, typename T2, typename TOP>
struct for_each_wrapper_while_2 {
    for_each_wrapper_while_2(TOP &op): op_(op) { }
    void stepz(const tom::loop::index_type &z) { }
    void stepy(const tom::loop::index_type &y) { }
    bool operator()(T1 &val1, T2 &val2, const tom::loop::packed_idx &idx) {
        return op_(val1, val2, idx.x, idx.y, idx.z);
    }
    TOP &op_;
};
} } } } // namespace tom::loop::functor::unnamed
/****************************************************************************//**
 * \brief Traversing two volumes as long as the operation returns true.
 *
 * Warning: if you would like to get a result from the functor, you have to
 * explicitly state the type \c TOP as reference. Otherwise a copy is passed.
 * Or use the return value of the loop.
 *******************************************************************************/
template<typename TVOL1, typename TVOL2, typename TOP>
inline void tom::loop::for_each_while(TVOL1 &v1, TVOL2 &v2, TOP op) {
    for_each_general(v1, v2,
        functor::for_each_wrapper_while_2<  typename resolve_volume_type<TVOL1>::data_type,
                                            typename resolve_volume_type<TVOL2>::data_type,
                                            typename boost::remove_reference<TOP>::type
        >(op));
}




namespace tom { namespace loop { namespace functor { namespace {
// Structure for wrapping the functor to another signature of the function call.
template<typename T1, typename T2, typename TOP>
struct for_each_wrapper_step_2 {
    for_each_wrapper_step_2(TOP &op): op_(op) { }
    void stepz(const tom::loop::index_type &z) {
        op_.stepz(z);
    }
    void stepy(const tom::loop::index_type &y) {
        op_.stepy(y);
    }
    bool operator()(T1 &val1, T2 &val2, const tom::loop::packed_idx &idx) {
        op_(val1, val2, idx.x);
        return true;
    }
    TOP &op_;
};
} } } } // namespace tom::loop::functor::unnamed
/****************************************************************************//**
 * \brief Traversing two volumes.
 *
 * Warning: if you would like to get a result from the functor, you have to
 * explicitly state the type \c TOP as reference. Otherwise a copy is passed.
 * Or use the return value of the loop.
 *******************************************************************************/
template<typename TVOL1, typename TVOL2, typename TOP>
inline void tom::loop::for_each_step(TVOL1 &v1, TVOL2 &v2, TOP op) {
    for_each_general(v1, v2,
        functor::for_each_wrapper_step_2<   typename resolve_volume_type<TVOL1>::data_type,
                                            typename resolve_volume_type<TVOL2>::data_type,
                                            typename boost::remove_reference<TOP>::type
        >(op));
}





/****************************************************************************//**
 * \brief Traversing two volumes.
 *
 * While in \c tom::loop::for_each only one volume is traversed, now two volumes
 * of the same size (but maybe different types) are traversed elementwise.
 * Usefull for example to copy the content of one volume to the second one
 * (e.g. tom::Volume<T>::setValues(const tom::Volume<T> &v1)).
 *
 * Warning: if you would like to get a result from the functor, you have to
 * explicitly state the type \c TOP as reference. Otherwise a copy is passed.
 * Or use the return value of the loop.
 *******************************************************************************/
template<typename TVOL1, typename TVOL2, typename TOP>
inline void tom::loop::for_each_general(TVOL1 &v1, TVOL2 &v2, TOP op) {
    if (!v1.is_equal_size(v2)) {
        throw std::invalid_argument("Both volumes must have the same size.");
    }
    assert( v1.getSizeX()==static_cast<index_type>(v1.getSizeX()) && v1.getSizeY()==static_cast<index_type>(v1.getSizeY()) && v1.getSizeZ()==static_cast<index_type>(v1.getSizeZ()) );
    assert( v2.getSizeX()==static_cast<index_type>(v2.getSizeX()) && v2.getSizeY()==static_cast<index_type>(v2.getSizeY()) && v2.getSizeZ()==static_cast<index_type>(v2.getSizeZ()) );

    typedef typename resolve_volume_type<TVOL1>::data_type TPTR1;
    typedef typename resolve_volume_type<TVOL2>::data_type TPTR2;

    TPTR1 *ptr1 = &v1.get();
    TPTR2 *ptr2 = &v2.get();
    const index_type sizex = v1.getSizeX();
    const index_type sizey = v1.getSizeY();
    const index_type sizez = v1.getSizeZ();
    std::size_t stridex1 = v1.getStrideX();
    std::size_t stridey1 = v1.getStrideY();
    std::size_t stridez1 = v1.getStrideZ();
    std::size_t stridex2 = v2.getStrideX();
    std::size_t stridey2 = v2.getStrideY();
    std::size_t stridez2 = v2.getStrideZ();
    packed_idx idx;
    assert(sizex && sizey && sizez && stridex1 && stridey1 && stridez1 && stridex2 && stridey2 && stridez2);

    if (!(stridex1%sizeof(TPTR1)) && !(stridey1%sizeof(TPTR1)) && !(stridez1%sizeof(TPTR1)) &&
        !(stridex2%sizeof(TPTR2)) && !(stridey2%sizeof(TPTR2)) && !(stridez2%sizeof(TPTR2))) {

        stridex1 /= sizeof(TPTR1); stridey1 /= sizeof(TPTR1); stridez1 /= sizeof(TPTR1);
        stridex2 /= sizeof(TPTR2); stridey2 /= sizeof(TPTR2); stridez2 /= sizeof(TPTR2);
        assert(stridex1 && stridey1 && stridez1 && stridex2 && stridey2 && stridez2);

        if (stridex1 == 1 && stridex2 == 1) {
            if (stridey1 == sizex &&
                stridez1 == sizex*sizey &&
                stridey1 == stridey2 &&
                stridez1 == stridez2) {
                for (idx.z=0; idx.z<sizez; idx.z++) {
                    op.stepz(idx.z);
                    for (idx.y=0; idx.y<sizey; idx.y++) {
                        op.stepy(idx.y);
                        for (idx.x=0; idx.x<sizex; idx.x++) {
                            if (!op(*ptr1++, *ptr2++, const_cast<const packed_idx &>(idx))) {
                                return;
                            }
                        }
                    }
                }
            } else {
                stridez1 -= sizey*stridey1;
                stridez2 -= sizey*stridey2;

                for (idx.z=0; idx.z<sizez; idx.z++) {
                    op.stepz(idx.z);
                    for (idx.y=0; idx.y<sizey; idx.y++) {
                        op.stepy(idx.y);
                        for (idx.x=0; idx.x<sizex; idx.x++) {
                            if (!op(ptr1[idx.x], ptr2[idx.x], idx)) {
                                return;
                            }
                        }
                        ptr1 += stridey1;
                        ptr2 += stridey2;
                    }
                    ptr1 += stridez1;
                    ptr2 += stridez2;
                }
            }
        } else {
            stridez1 -= sizey*stridey1;
            stridez2 -= sizey*stridey2;
            stridey1 -= sizex*stridex1;
            stridey2 -= sizex*stridex2;
            for (idx.z=0; idx.z<sizez; idx.z++) {
                op.stepz(idx.z);
                for (idx.y=0; idx.y<sizey; idx.y++) {
                    op.stepy(idx.y);
                    for (idx.x=0; idx.x<sizex; idx.x++) {
                        if (!op(ptr1[idx.x], ptr2[idx.x], idx)) {
                            return;
                        }
                        ptr1 += stridex1;
                        ptr2 += stridex2;
                    }
                    ptr1 += stridey1;
                    ptr2 += stridey2;
                }
                ptr1 += stridez1;
                ptr2 += stridez2;
            }
        }
    } else {
        stridez1 -= sizey*stridey1;
        stridez2 -= sizey*stridey2;
        stridey1 -= sizex*stridex1;
        stridey2 -= sizex*stridex2;
        for (idx.z=0; idx.z<sizez; idx.z++) {
            op.stepz(idx.z);
            for (idx.y=0; idx.y<sizey; idx.y++) {
                op.stepy(idx.y);
                for (idx.x=0; idx.x<sizex; idx.x++) {
                    if (!op(ptr1[idx.x], ptr2[idx.x], idx)) {
                        return;
                    }
                    ptr_add_byte_offset(ptr1, stridex1);
                    ptr_add_byte_offset(ptr2, stridex2);
                }
                ptr_add_byte_offset(ptr1, stridey1);
                ptr_add_byte_offset(ptr2, stridey2);
            }
            ptr_add_byte_offset(ptr1, stridez1);
            ptr_add_byte_offset(ptr2, stridez2);
        }
    }
}






namespace tom { namespace loop { namespace functor { namespace {
// Structure for wrapping the functor to another signature of the function call.
template<typename T1, typename T2, typename T3, typename TOP>
struct for_each_wrapper_3 {
    for_each_wrapper_3(TOP &op): op_(op) { }
    void operator()(T1 &val1, T2 &val2, T3 &val3, const tom::loop::packed_idx &idx) {
        op_(val1, val2, val3, idx.x, idx.y, idx.z);
    }
    TOP &op_;
};
} } } } // namespace tom::loop::functor::unnamed
/***************************************************************************//**
 * \brief Traversing three volumes.
 *
 * While in \c tom::loop::for_each only one volume is traversed, now three volumes
 * of the same size (but maybe different types) are traversed elementwise.
 *
 * Warning: if you would like to get a result from the functor, you have to
 * explicitly state the type \c TOP as reference. Otherwise a copy is passed.
 * Or use the return value of the loop.
 ******************************************************************************/
template<typename TVOL1, typename TVOL2, typename TVOL3, typename TOP>
inline void tom::loop::for_each(TVOL1 &v1, TVOL2 &v2, TVOL3 &v3, TOP op) {
    typedef functor::for_each_wrapper_3<typename resolve_volume_type<TVOL1>::data_type,
                                        typename resolve_volume_type<TVOL2>::data_type,
                                        typename resolve_volume_type<TVOL3>::data_type,
                                        typename boost::remove_reference<TOP>::type
                                        > wrapper_type;
    for_each_packed_idx(v1, v2, v3, wrapper_type(op));
}





namespace tom { namespace loop { namespace functor { namespace {
// Structure for wrapping the functor to another signature of the function call.
template<typename T1, typename T2, typename T3, typename TOP>
struct for_each_wrapper_no_idx_3 {
    for_each_wrapper_no_idx_3(TOP &op): op_(op) { }
    void operator()(T1 &val1, T2 &val2, T3 &val3, const tom::loop::packed_idx &) {
        op_(val1, val2, val3);
    }
    TOP &op_;
};
} } } } // namespace tom::loop::functor::unnamed
/***************************************************************************//**
 * \brief Traversing three volumes.
 *
 * While in \c tom::loop::for_each only one volume is traversed, now three volumes
 * of the same size (but maybe different types) are traversed elementwise.
 *
 * Warning: if you would like to get a result from the functor, you have to
 * explicitly state the type \c TOP as reference. Otherwise a copy is passed.
 * Or use the return value of the loop.
 ******************************************************************************/
template<typename TVOL1, typename TVOL2, typename TVOL3, typename TOP>
inline void tom::loop::for_each_no_idx(TVOL1 &v1, TVOL2 &v2, TVOL3 &v3, TOP op) {
    typedef functor::for_each_wrapper_no_idx_3< typename resolve_volume_type<TVOL1>::data_type,
                                                typename resolve_volume_type<TVOL2>::data_type,
                                                typename resolve_volume_type<TVOL3>::data_type,
                                                typename boost::remove_reference<TOP>::type
                                               > wrapper_type;
    for_each_packed_idx(v1, v2, v3, wrapper_type(op));
}




/****************************************************************************//**
 * \brief Traversing three volumes elementwise.
 *
 * Warning: if you would like to get a result from the functor, you have to
 * explicitly state the type \c TOP as reference. Otherwise a copy is passed.
 * Or use the return value of the loop.
 *******************************************************************************/
template<typename TVOL1, typename TVOL2, typename TVOL3, typename TOP>
inline void tom::loop::for_each_packed_idx(TVOL1 &v1, TVOL2 &v2, TVOL3 &v3, TOP op) {

    if (!v1.is_equal_size(v2) || !v1.is_equal_size(v3)) {
        throw std::invalid_argument("All volumes must have the same size.");
    }
    assert( v1.getSizeX()==static_cast<index_type>(v1.getSizeX()) && v1.getSizeY()==static_cast<index_type>(v1.getSizeY()) && v1.getSizeZ()==static_cast<index_type>(v1.getSizeZ()));
    assert( v2.getSizeX()==static_cast<index_type>(v2.getSizeX()) && v2.getSizeY()==static_cast<index_type>(v2.getSizeY()) && v2.getSizeZ()==static_cast<index_type>(v2.getSizeZ()));
    assert( v3.getSizeX()==static_cast<index_type>(v3.getSizeX()) && v3.getSizeY()==static_cast<index_type>(v3.getSizeY()) && v3.getSizeZ()==static_cast<index_type>(v3.getSizeZ()));

    typedef typename resolve_volume_type<TVOL1>::data_type TPTR1;
    typedef typename resolve_volume_type<TVOL2>::data_type TPTR2;
    typedef typename resolve_volume_type<TVOL3>::data_type TPTR3;

    TPTR1 *ptr1 = &v1.get();
    TPTR2 *ptr2 = &v2.get();
    TPTR3 *ptr3 = &v3.get();
    const index_type sizex = v1.getSizeX();
    const index_type sizey = v1.getSizeY();
    const index_type sizez = v1.getSizeZ();
    std::size_t stridex1 = v1.getStrideX();
    std::size_t stridey1 = v1.getStrideY();
    std::size_t stridez1 = v1.getStrideZ();
    std::size_t stridex2 = v2.getStrideX();
    std::size_t stridey2 = v2.getStrideY();
    std::size_t stridez2 = v2.getStrideZ();
    std::size_t stridex3 = v3.getStrideX();
    std::size_t stridey3 = v3.getStrideY();
    std::size_t stridez3 = v3.getStrideZ();
    packed_idx idx;
    assert(sizex && sizey && sizez && stridex1 && stridey1 && stridez1 && stridex2 && stridey2 && stridez2 && stridex3 && stridey3 && stridez3);

    if (!(stridex1%sizeof(TPTR1)) && !(stridey1%sizeof(TPTR1)) && !(stridez1%sizeof(TPTR1)) &&
        !(stridex2%sizeof(TPTR2)) && !(stridey2%sizeof(TPTR2)) && !(stridez2%sizeof(TPTR2)) &&
        !(stridex3%sizeof(TPTR3)) && !(stridey3%sizeof(TPTR3)) && !(stridez3%sizeof(TPTR3))) {

        stridex1 /= sizeof(TPTR1); stridey1 /= sizeof(TPTR1); stridez1 /= sizeof(TPTR1);
        stridex2 /= sizeof(TPTR2); stridey2 /= sizeof(TPTR2); stridez2 /= sizeof(TPTR2);
        stridex3 /= sizeof(TPTR3); stridey3 /= sizeof(TPTR3); stridez3 /= sizeof(TPTR3);
        assert(stridex1 && stridey1 && stridez1 && stridex2 && stridey2 && stridez2 && stridex3 && stridey3 && stridez3);

        if (stridex1 == 1 && stridex2 == 1 && stridex3 == 1) {
            if (stridey1 == sizex &&
                stridez1 == sizex*sizey &&
                stridey1 == stridey2 &&
                stridez1 == stridez2 &&
                stridey1 == stridey3 &&
                stridez1 == stridez3) {
                std::size_t i = 0;
                for (idx.z=0; idx.z<sizez; idx.z++) {
                    for (idx.y=0; idx.y<sizey; idx.y++) {
                        for (idx.x=0; idx.x<sizex; idx.x++, i++) {
                            op(ptr1[i], ptr2[i], ptr3[i], const_cast<const packed_idx &>(idx));
                        }
                    }
                }
            } else {
                stridez1 -= sizey*stridey1;
                stridez2 -= sizey*stridey2;
                stridez3 -= sizey*stridey3;
                for (idx.z=0; idx.z<sizez; idx.z++) {
                    for (idx.y=0; idx.y<sizey; idx.y++) {
                        for (idx.x=0; idx.x<sizex; idx.x++) {
                            op(ptr1[idx.x], ptr2[idx.x], ptr3[idx.x], idx);
                        }
                        ptr1 += stridey1;
                        ptr2 += stridey2;
                        ptr3 += stridey3;
                    }
                    ptr1 += stridez1;
                    ptr2 += stridez2;
                    ptr3 += stridez3;
                }
            }
        } else {
            stridez1 -= sizey*stridey1; stridez2 -= sizey*stridey2; stridez3 -= sizey*stridey3;
            stridey1 -= sizex*stridex1; stridey2 -= sizex*stridex2; stridey3 -= sizex*stridex3;
            for (idx.z=0; idx.z<sizez; idx.z++) {
                for (idx.y=0; idx.y<sizey; idx.y++) {
                    for (idx.x=0; idx.x<sizex; idx.x++) {
                        op(ptr1[idx.x], ptr2[idx.x], ptr3[idx.x], idx);
                        ptr1 += stridex1;
                        ptr2 += stridex2;
                        ptr3 += stridex3;
                    }
                    ptr1 += stridey1;
                    ptr2 += stridey2;
                    ptr3 += stridey3;
                }
                ptr1 += stridez1;
                ptr2 += stridez2;
                ptr3 += stridey3;
            }
        }
    } else {
        stridez1 -= sizey*stridey1;
        stridez2 -= sizey*stridey2;
        stridez3 -= sizey*stridey3;
        stridey1 -= sizex*stridex1;
        stridey2 -= sizex*stridex2;
        stridey3 -= sizex*stridex3;
        for (idx.z=0; idx.z<sizez; idx.z++) {
            for (idx.y=0; idx.y<sizey; idx.y++) {
                for (idx.x=0; idx.x<sizex; idx.x++) {
                    op(ptr1[idx.x], ptr2[idx.x], ptr3[idx.x], idx);
                    ptr_add_byte_offset(ptr1, stridex1);
                    ptr_add_byte_offset(ptr2, stridex2);
                    ptr_add_byte_offset(ptr3, stridex3);
                }
                ptr_add_byte_offset(ptr1, stridey1);
                ptr_add_byte_offset(ptr2, stridey2);
                ptr_add_byte_offset(ptr3, stridey3);
            }
            ptr_add_byte_offset(ptr1, stridez1);
            ptr_add_byte_offset(ptr2, stridez2);
            ptr_add_byte_offset(ptr3, stridez3);
        }
    }
}








#endif


