///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2009                                                 //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file LICENSE.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


//-----------------------------------------------------------------------------
//
// This file defines the Vector class.
//
// A Vector is a mathematical vector whose size is not necessarily
// known at compile time.  (c.f. SmallVector in  TMV_SmallVector.h.)
//
// The Vector class and all associated functions are contained
// in the namespace tmv.  Alse, the Vector class is a template, 
// so for a Vector of doubles, one would write 
// tmv::Vector<double>.  
//
// An optional second template parameter can tell the Vector to
// use Fortran-style indexing instead of C-style.
// ie. the first element is v(1), and the last is v(n).
// So you would write Vector<T,FortranStyle> v(n), for example.
// If you want to use the Fortran-style indexing, you should
// be aware that it only exists for a regular Vector and the two
// view classes, VectorView and ConstVectorView.  
//
// The return type of many of the below methods are given generically
// as value_type, real_type, iterator, etc.
// All of these types are typedefs in Vector.  
// So, you could write, for example:
//
// tmv::Vector<double> v(30);
// typename tmv::Vector<double>::iterator it = v.begin();
//
// The same form can be used for any return type that is not obvious.
// And the same is true of the return types that are vectors.  e.g.:
//
// typedef typename tmv::Vector<double>::subvector_type sub;
// typename sub::reverse_iterator it = v.subVector(10,20).rbegin();
//
//
// Constructors:
//
//    explicit Vector<T>(size_t n)
//        Makes a Vector of size n with _uninitialized_ values
//
//    Vector<T>(size_t n, T x)
//        Makes a Vector of size n with all values = x
//
//    Vector<T>(size_t n, const T* v2)
//    Vector<T>(const std::vector<T>& v2)
//        Makes a vector which copies the elements of v2
//        For the second one, n specifies the length of the vector
//
// Special Creators: 
//
//    VectorView VectorViewOf(T* v, size_t n, int step=1) 
//    ConstVectorView VectorViewOf(const T* v, size_t n, int step=1) 
//        Makes a VectorView of the memory elements stored at vv 
//
// Access
//
//    size_t size() const
//        Returns the size of the Vector
//
//    value_type operator[](int i) const
//    value_type operator()(int i) const
//    value_type cref(int i) const
//        Return the ith element of the vector.
//        The first two respect the index-style of the underlying vector. 
//        The third, cref, always uses CStyle indexing, and does not
//        do any checking of the validity of value of i, so it can
//        be used to speed-up code that has element access in a critical
//        section, when you don't otherwise want to compile with -DNDEBUG
//
//    reference operator[](int i)
//    reference operator()(int i)
//    reference ref(int i)
//        Return a reference to the ith element of the vector.
//        The first two respect the index-style of the underlying vector. 
//        The third, ref, always uses CStyle indexing, and does not
//        do any checking of the validity of value of i.
//        For normal vectors, reference is just value&.  However,
//        we allow for conjugate views of vectors, so in that case,
//        the return type is an object that can be used as a reference,
//        but which applies the conjugation appropriately.
//
//    iterator begin()
//    iterator end()
//    reverse_iterator rbegin()
//    reverse_iterator rend()
//    const_iterator begin() const
//    const_iterator end() const
//    const_reverse_iterator rbegin() const
//    const_reverse_iterator rend() const
//        Return iterators that work in the usual way to traverse the vector
//
// Non-modifying Functions
//
//    Most of these are both member functions and functions of a vector,
//    so Norm(v) and v.norm(), for example, are equivalent.
//
//    real_type norm() const    or Norm(v)
//    real_type norm2() const    or Norm2(v)
//        Returns the 2-norm of a vector = sqrt( sum_i |v_i|^2 )
//
//    real_type normSq() const    or NormSq(v)
//    real_type normSq(real_type scale) const
//        Returns the square of norm().
//        In the method version, you can provide an optional scale, in 
//        which case the output is equal to NormSq(scale*v).
//
//    real_type norm1() const    or Norm1(v) 
//        Returns the 1-norm of a vector = sum_i |v_i|
//
//    real_type normInf() const    or NormInf(v) 
//        Returns the infinity-norm of a vector = max_i |v_i|
//
//    value_type sumElements() const    or SumElements(v) 
//        Returns the sum of all elements in the vector.
//
//    real_type sumAbsElements() const    or SumAbsElements(v) 
//        Returns the sum of absolute values of elements in the vector.
//
//    real_type sumAbs2Elements() const    or SumAbs2Elements(v) 
//        Returns the sum of absolute values of elements in the vector
//        using |real(v(i))| + |imag(v(i))| for complex vectors.
//
//    value_type maxElement() const    or MaxElement(v) 
//    value_type maxElement(int* imax) const
//        Returns the maximum value of any element in the vector.
//        In the method version, you can provide an optional argument imax.
//        On return, *imax holds the index of the element with the 
//        maximum value.
//        The parameter imax can be omitted if it is not desired.
//        As "max" doesn't make sense for complex values, for these
//        we use just the real components.
//
//    value_type minElement() const    or MinElement(v) 
//    value_type minElement(int* imin) const
//        Returns the minimum value of any element in the vector.
//        On return, (in the method version) *imin holds the index of this 
//        element.  Again, the parameter imin can be omitted.
//
//    real_type maxAbsElement() const    or MaxAbsElement(v) 
//    real_type maxAbsElement(int* imax) const
//        The same as MaxElement, except absolute values are used.
//
//    real_type minAbsElement() const    or MinAbsElement(v) 
//    real_type minAbsElement(int* imax) const
//        The same as minElement, except absolute values are used.
//
//    real_type maxAbs2Element() const    or MaxAbs2Element(v) 
//    real_type maxAbs2Element(int* imax) const
//    real_type minAbs2Element() const    or MinAbs2Element(v) 
//    real_type minAbs2Element(int* imax) const
//        The same as the above min/maxAbsElement, but using
//        |real(v(i))| + |imag(v(i))| for complex vectors.
//
//    template <class ret_type>
//    ret_type maxElement(const F& f) const
//    ret_type maxElement(int* imax, const F& f) const
//        Returns the minimum value of f(v[i]).  
//        The funcion f can be any object that has the method:
//        ret_type operator()(value_type x)
//        When calling this, you need to specify ret_type as a template
//        parameter, but F (which is also in fact a template parameter
//        to this method) is inferred from the argument f.
//
//    template <class ret_type>
//    ret_type minElement(int* imin, const F& f) const
//        Returns the minimum value of f(v[i]).  
// 
//    template <class ret_type>
//    ret_type sumElements(const F& f) const
//        Returns the sum of the values f(v[i]).
//
//
// Modifying Functions
//
//    type& zero()
//        Sets all elements to 0
//
//    type& clip(real_type thresh)
//        Set to 0 all elements whose absolute values is < thresh
//
//    type& setAllTo(T x)
//        Sets all elements to x
//        We don't call this v = x, since it doesn't really make sense to
//        think of v as being 'equal' to a scalar.
//        Hence the slightly verbose function name setAllTo.
//
//    type& addToAll(T x)
//        Add x to each element of v
//
//    type& conjugateSelf()
//        Sets all elements to its conjugate
//
//    type& makeBasis(int i)
//        Set all elements to 0, except v(i) = 1
//
//    type& swap(int i1, int i2)
//        Swap elements v(i1) and v(i2)
//
//    type& permute(const int* p)
//    type& permute(const int* p, int i1, int i2)
//        Perform a series of swaps (0,p[0]), (1,p[1])...
//        In the second case, only do (i1,p[i1])...(i2-1,p[i2-1])
//    type& reversePermute(const int* p)
//    type& reversePermute(const int* p, int i1, int i2)
//        The same but perform the swaps in reverse order.
//
//    type& reverseSelf()
//        Reverse the order of the elements of v
//
//    void sort(ad, comp)
//    void sort(int* P, ad, comp)
//        Sorts the vector, returning the swaps required in P.
//        If you do not care about P, you may omit the P parameter.
//        ad = Ascend or Descend (Ascend=default)
//        comp = RealComp, AbsComp, ImagComp, or ArgComp (RealComp=default)
//
//    Swap(v1,v2)
//        Swap vectors v1 and v2
//
// VectorView:
//
//    A VectorView object refers to the elements of a regular Vector
//    or Matrix, so that altering the elements in the view alters the
//    corresponding elements in the original object.  A View can have non-unit
//    steps between elements and can also be a conjugation of the original
//    elements, so that v[3] = z would actually set the original element
//    to z*. 
//
//    Also, since we want to be able to pass these Views around, we have to 
//    keep track of whether we are allowed to alter the original values or
//    just look at them.  Thus, there are two VectorView objects: 
//    ConstVectorView and VectorView.  The first is only allowed to view,
//    not modify, the original elements.  The second is allowed to modify them.
//    This is akin to the const_iterator and iterator types for the
//    standard template library.
//
//    Both of these take several template arguments:
//    T = the underlying data type.
//    S = (int) the step size if known.  (Use UNKNOWN if not known.)
//    C = (bool) whether the view is the conjugate of the data.
//    I = (CStyle or FortranStyle) which indexing style to use.
// 
//    There are default arguments for all but T:
//    S = 1, C = false, I = CStyle
//    so you can omit these in many cases.
//
//    ConstVectorView<T,S,C,I>(const T* p, size_t n, int s=S)
//    VectorView<T,S,C,I>(T* p, size_t n, int s=S)
//        Make a vector view, starting at memory location p, 
//        with n elements, stepping over the data with step size s.
//        Note, that s is generally omitted if you are specifying S
//        as a template argument.  But if S == UNKNOWN, then you would
//        need to specify s as an argument.
//
//    There are also copy constructors that change S (e.g. from UNKNOWN
//    to some compile-time-known value), or change VectorView to 
//    ConstVectorView.  e.g.:
//
//    ConstVectorView<T,S,C,I>(const VectorView<T,S,C,I>& v);
//    ConstVectorView<T,S,C,I>(const VectorView<T,UNKNOWN,C,I>& v);
//    VectorView<T,S,C,I>(const VectorView<T,UNKNOWN,C,I>& v);
// 
// Views:
//
//    All of these methods return some kind of view into the vector
//    data.  The return types are either ConstVectorView or VectorView
//    with the appropriate values for the template parameters S and C.
//    The template parameter I is preserved from the original vector.
//
//    All of these have a const and non-const version.  For the const
//    versions, prepend const_ to the return type name.  
//    e.g. const_subvector_type.  
//   
//    subvector_type subVector(int i1, int i2)
//        This member function returns a subvector view which refers
//        to the same physical elements as the original.
//        i1 is the first element in the subvector.
//        i2 is the usual "one past the end" of the subvector.
//        Thus, for a vector, v, of length 10, you could output the 
//        the first 4 elements of v with:
//        std::cout << v.subVector(0,4);
//      
//    subvector_step_type subVector(int i1, int i2, int istep)
//        This returns a view with non-unit step through the elements.
//        Thus, if you have a vector v of length 10, and you want to
//        set all the even elements to 0, you could write:
//        v.subVector(0,10,2).zero();
//
//    reverse_type reverse()
//        Returns a subvector whose elements are the same as v but 
//        in the reverse order
//
//    conjugate_type conjugate()
//        Returns the conjugate of a Vector (as a view, so it still points
//        to the same physical elements, but modifying this will set the 
//        actual elements to the conjugate of what you set.
//
//    view_type view()
//        Returns a view of a Vector. 
//
//    cview_type cView()
//    fview_type fView()
//        Review the vector with CStyle/FortranStyle respectively.
//
//    xview_type xView()
//        Returns a view of a Vector with the template step parameter 
//        equal to UNKNOWN, rather than any known value.
//        Also, it is always CStyle indexing.
//
//    unitview_type unitView()
//        Returns a view of a Vector with the template step parameter 
//        equal to 1.  This is only useful if the original step parameter
//        is UNKNOWN, since the actual step size must = 1.
//
//    realpart_type realPart()
//    imagpart_type imagPart()
//        Returns a view of the real/imag elements of a complex Vector.
//
//    flatten_type flatten()
//        For complex vectors, this returns a real vector view to 
//        both the real and imag elements in a single view.
//        If you call this on a view, it must have unit step.
//
//    nonconj_type nonConj()
//        Returns a view of the underlying memory elements, removing
//        any vconj that might be set in the current view.
//        This is sometimes useful when an operation has the same
//        effect regardless of vconj, so it is better to just ignore
//        any vconj value.
//
//    nonconst_type nonConst()
//        Returns a mutable view of a const Vector.
//
// Operators:
//    Here we use v for a Vector and x for a Scalar.
//
//    -v
//
//    v = v
//
//    v += v
//    v + v
//
//    v -= v
//    v - v
//
//    v *= x
//    x * v
//    v * x
//    v * v   (inner product)
//
//    v /= x
//    v / x
//
//    v == v
//    v != v
//
//    These all behave in the logical way for dealing with vectors.
//
//
// I/O: 
//
//    void write(std::ostream& os) const    or os << v
//        Writes v to ostream os in the following format:
//          n ( v[0] v[1] v[2] ... v[n] )
//
//    v.write(ostream& os, real_type thresh)
//        Write v to os as above, but if |v[i]| < thresh, write 0 instead
//
//    void read(std::istream& is)    or is >> v
//        Reads the vector from istream is in the same format
//        Note: the vector must already be the correct size
//
//    std::auto_ptr<tmv::Vector<T> > vptr;
//    is >> vptr
//        If you do not know the size of the vector to be read in, you can
//        use this form which will allocate the vector to be the correct
//        size according to the input data.
//
//


#ifndef TMV_Vector_H
#define TMV_Vector_H

#include <vector>
#include "TMV_BaseVector.h"
#include "TMV_VIt.h"
#include "TMV_Array.h"

namespace tmv {

    //
    // Vector
    //

    template <class T, IndexStyle I>
    struct Traits<Vector<T,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef Vector<T,I> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef type copy_type;

        enum { _size = UNKNOWN };
        enum { _fort = (I == FortranStyle) };
        enum { _calc = true };
        enum { _step = 1 };
        enum { _conj = false };
        enum { twoS = isreal ? 1 : 2 };

        typedef ConstVectorView<T,1,false,I> const_subvector_type;
        typedef ConstVectorView<T,UNKNOWN,false,I> const_subvector_step_type;
        typedef ConstVectorView<T,1,false,I> const_view_type;
        typedef ConstVectorView<T,1,false,CStyle> const_cview_type;
        typedef ConstVectorView<T,1,false,FortranStyle> const_fview_type;
        typedef ConstVectorView<T> const_xview_type;
        typedef ConstVectorView<T,1,false,I> const_unitview_type;
        typedef ConstVectorView<T,1,iscomplex,I> const_conjugate_type;
        typedef ConstVectorView<T,-1,false,I> const_reverse_type;
        typedef ConstVectorView<real_type,twoS,false,I> const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstVectorView<real_type,1,false,I> const_flatten_type;
        typedef ConstVectorView<T,1,false,I> const_nonconj_type;
        typedef VectorView<T,1,false,I> nonconst_type;

        typedef CVIt<T,1,false> const_iterator;
        typedef CVIt<T,-1,false> const_reverse_iterator;

        typedef T& reference;

        typedef VectorView<T,1,false,I> subvector_type;
        typedef VectorView<T,UNKNOWN,false,I> subvector_step_type;
        typedef VectorView<T,1,false,I> view_type;
        typedef VectorView<T,1,false,CStyle> cview_type;
        typedef VectorView<T,1,false,FortranStyle> fview_type;
        typedef VectorView<T> xview_type;
        typedef VectorView<T,1,false,I> unitview_type;
        typedef VectorView<T,1,iscomplex,I> conjugate_type;
        typedef VectorView<T,-1,false,I> reverse_type;
        typedef VectorView<real_type,twoS,false,I> realpart_type;
        typedef realpart_type imagpart_type;
        typedef VectorView<real_type,1,false,I> flatten_type;
        typedef VectorView<T,1,false,I> nonconj_type;

        typedef VIt<T,1,false> iterator;
        typedef VIt<T,-1,false> reverse_iterator;
    };

    template <class T, IndexStyle I>
    class Vector : 
        public BaseVector_Mutable<Vector<T,I> >
    {
    public:

        typedef Vector<T,I> type;
        typedef BaseVector_Mutable<type> base_mut;

        enum { _size = Traits<type>::_size };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { isreal = Traits<type>::isreal };
        enum { iscomplex = Traits<type>::iscomplex };
        enum { _step = Traits<type>::_step };
        enum { _conj = Traits<type>::_conj };

        //
        // Constructors
        //

        explicit Vector(size_t n=0) : itssize(n), itsv(n)
        {
            TMVAssert(n>=0);
#ifdef TMV_DEBUG
            this->setAllTo(T(888));
#endif
        }

        Vector(size_t n, T val) : itssize(n), itsv(n)
        {
            TMVAssert(n>=0);
            this->setAllTo(val);
        }

        Vector(size_t n, const T* v2) : itssize(n), itsv(n)
        {
            TMVAssert(n>=0);
            VectorViewOf(v2,n).newAssignTo(*this);
        }

        explicit Vector(const std::vector<T>& v2) : 
            itssize(v2.size()), itsv(itssize)
        {
            TMVAssert(itssize>=0);
            VectorViewOf(&v2[0],itssize).newAssignTo(*this);
        }

        Vector(const type& v2) : itssize(v2.size()), itsv(itssize)
        {
            TMVAssert(itssize>=0);
            v2.newAssignTo(*this);
        }

        template <class V2>
        Vector(const BaseVector<V2>& v2) :
            itssize(v2.size()), itsv(itssize)
        {
            TMVAssert(itssize>=0);
            v2.newAssignTo(*this);
        }

        ~Vector() 
        {
#ifdef TMV_DEBUG
            this->setAllTo(T(999));
#endif
        }


        //
        // Op =
        //

        type& operator=(type& v2)
        {
            TMVAssert(v2.size() == size());
            if (&v2 != this) v2.assignTo(*this);
            return *this; 
        }

        template <class V2>
        type& operator=(const BaseVector<V2>& v2)
        { base_mut::operator=(v2); return *this; }


        //
        // Auxilliary Functions
        //

        const T* cptr() const { return itsv; }
        T* ptr() { return itsv; }
        T cref(int i) const  { return itsv[i]; }
        T& ref(int i) { return itsv[i]; }

        size_t size() const { return itssize; }
        int nElements() const { return itssize; }
        int step() const { return 1; }
        bool isconj() const { return false; }
        void swapWith(type& rhs)
        {
            TMVAssert(rhs.size() == size());
            if (itsv.get() == rhs.itsv.get()) return;
            itsv.swapWith(rhs.itsv);
        }

        void resize(size_t n)
        {
            itssize = n;
            itsv.resize(n);
#ifdef TMV_DEBUG
            this->setAllTo(T(888));
#endif
        }

    private:

        size_t itssize;
        AlignedArray<T> itsv;

    }; // Vector

    template <class T>
    class VectorF : 
        public Vector<T,FortranStyle>
    {
    public:

        typedef VectorF<T> type;
        typedef Vector<T,FortranStyle> vtype;

        explicit VectorF(size_t n) : vtype(n) {}
        VectorF(size_t n, T val) : vtype(n,val) {}
        VectorF(size_t n, const T* v2) : vtype(v2) {}
        explicit VectorF(const std::vector<T>& v2) : vtype(v2) {}
        VectorF(const type& v2) : vtype(v2) {}
        template <class V2>
        VectorF(const BaseVector<V2>& v2) : vtype(v2) {}
        ~VectorF() {}

        type& operator=(type& v2)
        { vtype::operator=(v2); return *this; }
        template <class V2>
        type& operator=(const BaseVector<V2>& v2)
        { vtype::operator=(v2); return *this; }
    }; // VectorF

    //
    // ConstVectorView
    //

    template <class T, int S, bool C, IndexStyle I>
    struct Traits<ConstVectorView<T,S,C,I> >
    {
        typedef T value_type;
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef ConstVectorView<T,S,C,I> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef Vector<T,I> copy_type;

        enum { _size = UNKNOWN };
        enum { _fort = (I==FortranStyle) };
        enum { _calc = true };
        enum { _step = S };
        enum { _conj = C };
        enum { negS = IntTraits<S>::negS };
        enum { twoS = isreal ? S : IntTraits<S>::twoS };
        enum { notC = !C && iscomplex };

        typedef ConstVectorView<T,S,C,I> const_subvector_type;
        typedef ConstVectorView<T,UNKNOWN,C,I> const_subvector_step_type;
        typedef ConstVectorView<T,S,C,I> const_view_type;
        typedef ConstVectorView<T,S,C,CStyle> const_cview_type;
        typedef ConstVectorView<T,S,C,FortranStyle> const_fview_type;
        typedef ConstVectorView<T,UNKNOWN,C> const_xview_type;
        typedef ConstVectorView<T,1,C,I> const_unitview_type;
        typedef ConstVectorView<T,S,notC,I> const_conjugate_type;
        typedef ConstVectorView<T,negS,C,I> const_reverse_type;
        typedef ConstVectorView<real_type,twoS,false,I> const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstVectorView<real_type,1,false,I> const_flatten_type;
        typedef ConstVectorView<T,S,false,I> const_nonconj_type;
        typedef VectorView<T,S,C,I> nonconst_type;

        typedef CVIt<T,S,C> const_iterator;
        typedef CVIt<T,negS,C> const_reverse_iterator;
    };

    template <class T, int S, bool C, IndexStyle I>
    class ConstVectorView : 
        public BaseVector_Calc<ConstVectorView<T,S,C,I> >
    {
    public:

        typedef ConstVectorView<T,S,C,I> type;

        enum { _size = Traits<type>::_size };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { isreal = Traits<type>::isreal };
        enum { iscomplex = Traits<type>::iscomplex };
        enum { _step = Traits<type>::_step };
        enum { _conj = Traits<type>::_conj };

        //
        // Constructors
        //

        ConstVectorView(const T* v, size_t n, int s) : 
            itsv(v), itssize(n), itsstep(s) {}

        ConstVectorView(const T* v, size_t n) : 
            itsv(v), itssize(n), itsstep(S) 
        { TMVStaticAssert(S != UNKNOWN); }

        ConstVectorView(const type& v2) : 
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step()) {}

        ConstVectorView(const VectorView<T,S,C,I>& v2) : 
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step()) {}

        template <int S2, IndexStyle I2>
        ConstVectorView(const ConstVectorView<T,S2,C,I2>& v2) :
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step()) {}

        template <int S2, IndexStyle I2>
        ConstVectorView(const VectorView<T,S2,C,I2>& v2) :
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step()) {}

        template <int N2, int S2, IndexStyle I2>
        ConstVectorView(const ConstSmallVectorView<T,N2,S2,C,I2>& v2) :
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step()) {}

        template <int N2, int S2, IndexStyle I2>
        ConstVectorView(const SmallVectorView<T,N2,S2,C,I2>& v2) :
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step()) {}

        ~ConstVectorView() {
#ifdef TMV_DEBUG
            itsv = 0; 
#endif
        }

    private :
        void operator=(const type& v2);
    public :

        //
        // Auxilliary Functions
        //

        const T* cptr() const { return itsv; }

        T cref(int i) const  { return DoConj<C>(itsv[i*step()]); }

        size_t size() const { return itssize; }
        int nElements() const { return itssize; }
        int step() const { return itsstep; }
        bool isconj() const { return C; }

    protected :

        const T* itsv;
        const size_t itssize;
        const CheckedInt<S> itsstep;

    }; // ConstVectorView

    template <class T, int S, bool C>
    class ConstVectorViewF : 
        public ConstVectorView<T,S,C,FortranStyle>
    {
    public:
        typedef ConstVectorViewF<T,S,C> type;
        typedef ConstVectorView<T,S,C,FortranStyle> vtype;

        ConstVectorViewF(const T* v, size_t n, int s) : vtype(v,n,s) {}
        ConstVectorViewF(const T* v, size_t n) : vtype(v,n) {}
        ConstVectorViewF(const type& v2) : vtype(v2) {}
        ConstVectorViewF(const VectorView<T,S,C,FortranStyle>& v2) : 
            vtype(v2) {}
        template <int S2, IndexStyle I2>
        ConstVectorViewF(const ConstVectorView<T,S2,C,I2>& v2) :
            vtype(v2) {}
        template <int S2, IndexStyle I2>
        ConstVectorViewF(const VectorView<T,S2,C,I2>& v2) : vtype(v2) {}
        template <int N2, int S2, IndexStyle I2>
        ConstVectorViewF(const ConstSmallVectorView<T,N2,S2,C,I2>& v2) :
            vtype(v2) {}
        template <int N2, int S2, IndexStyle I2>
        ConstVectorViewF(const SmallVectorView<T,N2,S2,C,I2>& v2) : 
            vtype(v2) {}
        ~ConstVectorViewF() {}

    private :
        void operator=(const type& v2);
    }; // ConstVectorViewF


    //
    // VectorView
    //

    template <class T, int S, bool C, IndexStyle I>
    struct Traits<VectorView<T,S,C,I> >
    {
        typedef T value_type;
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef VectorView<T,S,C,I> type;
        typedef const ConstVectorView<T,S,C,I> calc_type;
        typedef calc_type eval_type;
        typedef Vector<T,I> copy_type;

        enum { _size = UNKNOWN };
        enum { _fort = (I==FortranStyle) };
        enum { _calc = true };
        enum { _step = S };
        enum { _conj = C };
        enum { negS = IntTraits<S>::negS };
        enum { twoS = isreal ? S : IntTraits<S>::twoS };
        enum { notC = !C && iscomplex };

        typedef ConstVectorView<T,S,C,I> const_subvector_type;
        typedef ConstVectorView<T,UNKNOWN,C,I> const_subvector_step_type;
        typedef ConstVectorView<T,S,C,I> const_view_type;
        typedef ConstVectorView<T,S,C,CStyle> const_cview_type;
        typedef ConstVectorView<T,S,C,FortranStyle> const_fview_type;
        typedef ConstVectorView<T,UNKNOWN,C> const_xview_type;
        typedef ConstVectorView<T,1,C,I> const_unitview_type;
        typedef ConstVectorView<T,S,notC,I> const_conjugate_type;
        typedef ConstVectorView<T,negS,C,I> const_reverse_type;
        typedef ConstVectorView<real_type,twoS,false,I> const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstVectorView<real_type,1,false,I> const_flatten_type;
        typedef ConstVectorView<T,S,false,I> const_nonconj_type;
        typedef VectorView<T,S,C,I> nonconst_type;

        typedef CVIt<T,S,C> const_iterator;
        typedef CVIt<T,negS,C> const_reverse_iterator;

        typedef typename AuxRef<T,C>::reference reference;

        typedef VectorView<T,S,C,I> subvector_type;
        typedef VectorView<T,UNKNOWN,C,I> subvector_step_type;
        typedef VectorView<T,S,C,I> view_type;
        typedef VectorView<T,S,C,CStyle> cview_type;
        typedef VectorView<T,S,C,FortranStyle> fview_type;
        typedef VectorView<T,UNKNOWN,C> xview_type;
        typedef VectorView<T,1,C,I> unitview_type;
        typedef VectorView<T,S,notC,I> conjugate_type;
        typedef VectorView<T,negS,C,I> reverse_type;
        typedef VectorView<real_type,twoS,false,I> realpart_type;
        typedef realpart_type imagpart_type;
        typedef VectorView<real_type,1,false,I> flatten_type;
        typedef VectorView<T,S,false,I> nonconj_type;

        typedef VIt<T,S,C> iterator;
        typedef VIt<T,negS,C> reverse_iterator;
    };

    template <class T, int S, bool C, IndexStyle I>
    class VectorView : 
        public BaseVector_Mutable<VectorView<T,S,C,I> >
    {
    public:

        typedef VectorView<T,S,C,I> type;
        typedef BaseVector_Mutable<type> base_mut;
        typedef typename base_mut::reference reference;

        enum { _size = Traits<type>::_size };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { isreal = Traits<type>::isreal };
        enum { iscomplex = Traits<type>::iscomplex };
        enum { _step = Traits<type>::_step };
        enum { _conj = Traits<type>::_conj };

        //
        // Constructors
        //

        VectorView(T* v, size_t n, int s) : 
            itsv(v), itssize(n), itsstep(s) {}

        VectorView(T* v, size_t n) : 
            itsv(v), itssize(n), itsstep(S) 
        { TMVStaticAssert(S != UNKNOWN); }

        VectorView(const type& v2) : 
            itsv(v2.itsv), itssize(v2.itssize), itsstep(v2.itsstep) {}

        template <int S2, IndexStyle I2>
        VectorView(VectorView<T,S2,C,I2> v2) :
            itsv(v2.ptr()), itssize(v2.size()), itsstep(v2.step()) {}

        template <int N2, int S2, IndexStyle I2>
        VectorView(SmallVectorView<T,N2,S2,C,I2> v2) :
            itsv(v2.ptr()), itssize(v2.size()), itsstep(v2.step()) {}

        ~VectorView() {
#ifdef TMV_DEBUG
            itsv = 0; 
#endif
        }


        //
        // Op =
        //

        template <class V2>
        type& operator=(const BaseVector<V2>& v2)
        { base_mut::operator=(v2); return *this; }

        type& operator=(const type& v2)
        { base_mut::operator=(v2); return *this; }

        //
        // Auxilliary Functions
        //

        const T* cptr() const { return itsv; }
        T* ptr() { return itsv; }

        T cref(int i) const  { return DoConj<C>(itsv[i*step()]); }
        reference ref(int i) { return reference(itsv[i*step()]); }

        size_t size() const { return itssize; }
        int nElements() const { return itssize; }
        int step() const { return itsstep; }
        bool isconj() const { return C; }

    protected :

        T* itsv;
        const size_t itssize;
        const CheckedInt<S> itsstep;

    }; // VectorView

    template <class T, int S, bool C>
    class VectorViewF : 
        public VectorView<T,S,C,FortranStyle>
    {
    public:
        typedef VectorViewF<T,S,C> type;
        typedef VectorView<T,S,C,FortranStyle> vtype;

        VectorViewF(T* v, size_t n, int s) : vtype(v,n,s) {}
        VectorViewF(T* v, size_t n) : vtype(v,n) {}
        VectorViewF(const type& v2) : vtype(v2) {}
        template <int S2, IndexStyle I2>
        VectorViewF(VectorView<T,S2,C,I2> v2) : vtype(v2) {}
        template <int N2, int S2, IndexStyle I2>
        VectorViewF(SmallVectorView<T,N2,S2,C,I2> v2) : vtype(v2) {}
        ~VectorViewF() {}

        template <class V2>
        type& operator=(const BaseVector<V2>& v2)
        { vtype::operator=(v2); return *this; }
        type& operator=(const type& v2)
        { vtype::operator=(v2); return *this; }
    }; // VectorViewF


    //
    // Special Creators: 
    //   VectorViewOf(v,size,step=1) = VectorView of v 
    //

    // VectorView of raw memory:
    template <class T>
    static inline VectorView<T,1> VectorViewOf(T* v, size_t size)
    { return VectorView<T,1>(v,size,1); }

    template <class T>
    static inline VectorView<T,UNKNOWN> VectorViewOf(
        T* v, size_t size, int step)
    { return VectorView<T,UNKNOWN>(v,size,step); }


    //
    // Swap
    //

    template <class T, IndexStyle I>
    static inline void Swap(Vector<T,I>& v1, Vector<T,I>& v2)
    { v1.swapWith(v2); }
    template <class V, class T, int S, bool C, IndexStyle I>
    static inline void Swap(BaseVector_Mutable<V>& v1, VectorView<T,S,C,I> v2)
    { DoSwap(v1,v2); }
    template <class V, class T, int S, bool C, IndexStyle I>
    static inline void Swap(VectorView<T,S,C,I> v1, BaseVector_Mutable<V>& v2)
    { DoSwap(v1,v2); }
    template <class T, int S1, bool C1, IndexStyle I1, int S2, bool C2, IndexStyle I2>
    static inline void Swap(
        VectorView<T,S1,C1,I1> v1, VectorView<T,S2,C2,I2> v2)
    { DoSwap(v1,v2); }


    //
    // Conjugate
    //
    
    template <class T, IndexStyle I>
    static inline VectorView<T,1,true,I> Conjugate(Vector<T,I>& v)
    { return v.conjugate(); }
    template <class T, int S, bool C, IndexStyle I>
    static inline VectorView<T,S,!C,I> Conjugate(VectorView<T,S,C,I> v)
    { return v.conjugate(); }


    //
    // TMV_Text 
    //

    template <class T, IndexStyle I>
    static inline std::string TMV_Text(const Vector<T,I>& )
    {
        std::ostringstream s;
        s << "Vector<"<<TMV_Text(T())<<","<<TMV_Text(I)<<">";
        return s.str();
    }

    template <class T, int S, bool C, IndexStyle I>
    static inline std::string TMV_Text(const ConstVectorView<T,S,C,I>& v)
    {
        std::ostringstream s;
        s << "ConstVectorView<"<<TMV_Text(T());
        s << ","<<IntTraits<S>::text();
        if (S == UNKNOWN) s << "("<<v.step()<<")";
        s <<","<<C<<","<<TMV_Text(I)<<">";
        return s.str();
    }

    template <class T, int S, bool C, IndexStyle I>
    static inline std::string TMV_Text(const VectorView<T,S,C,I>& v)
    {
        std::ostringstream s;
        s << "VectorView<"<<TMV_Text(T());
        s << ","<<IntTraits<S>::text();
        if (S == UNKNOWN) s << "("<<v.step()<<")";
        s <<","<<C<<","<<TMV_Text(I)<<">";
        return s.str();
    }

} // namespace tmv

#endif
