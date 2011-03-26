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
// An optional second template parameter specifies known attributes
// of the vector.  The only one that is valid is
// A = CStyle or FortranStyle.
// The default is CStyle, so you only need to specify if you want
// FortranStyle indexing.
// ie. the first element is v(1), and the last is v(n).
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
//    explicit Vector<T,A>(size_t n)
//        Makes a Vector of size n with _uninitialized_ values
//
//    Vector<T,A>(size_t n, T x)
//        Makes a Vector of size n with all values = x
//
//    Vector<T,A>(size_t n, const T* v2)
//    Vector<T,A>(const std::vector<T>& v2)
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
//    Both of these take two template arguments:
//    T = the underlying data type.
//    A = the known attributes of the vector.
//        Options are:
//        NonConj || Conj
//        NonUnit || Unit
//        CStyle || FortranStyle
// 
//    The default attributes are NonConj, NonUnit, CStyle, and you
//    only need to specify anything that differs from this.
//    So VectorView<T,Conj> means Conj, NonUnit, CStyle.
//    To specify more than one, use the bitwise or: |.
//    So VectorView<T,Unit | FortranStyle> means NonConj, Unit, FortranStyle.
//
//    ConstVectorView<T,A>(const T* p, size_t n, int s)
//    VectorView<T,A>(T* p, size_t n, int s)
//        Make a vector view, starting at memory location p, 
//        with n elements, stepping over the data with step size s.
//
//    You can also convert a VectorView into a ConstVectorView or
//    convert to a view with different attributes, so long as the 
//    conversion is logically valid.  e.g. you cannot just cast away
//    a Conj attribute, but you can change the indexing style or 
//    specify a Unit attribute if you know that the step size really
//    is 1.
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

    template <class T, int A0>
    struct Traits<Vector<T,A0> >
    {
        enum { A = (A0 & ~CheckAlias) | Unit };
        enum { okA = (
                Attrib<A>::vectoronly && 
                !Attrib<A>::conj )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef Vector<T,A0> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef type copy_type;

        enum { _size = UNKNOWN };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _step = 1 };
        enum { _conj = false };
        enum { _checkalias = !Attrib<A>::noalias };
        enum { _unit = true };
        enum { twoS = isreal ? 1 : 2 };

        enum { unitA = A };
        enum { nonunitA = A & ~Unit };
        enum { conjA = iscomplex ? (A ^ Conj) : int(A) };
        enum { nonconjA = A };
        enum { cstyleA = A & ~FortranStyle };
        enum { fstyleA = A | FortranStyle };
        enum { Ac = _checkalias ? (A | CheckAlias) : (A & ~NoAlias) };
        enum { nonunitAc = Ac & ~Unit };
        enum { twosAc = isreal ? int(Ac) : (Ac & ~Conj & ~Unit) };

        typedef ConstVectorView<T,A> const_subvector_type;
        typedef ConstVectorView<T,nonunitA> const_subvector_step_type;
        typedef ConstVectorView<T,A> const_view_type;
        typedef ConstVectorView<T,cstyleA> const_cview_type;
        typedef ConstVectorView<T,fstyleA> const_fview_type;
        typedef ConstVectorView<T> const_xview_type;
        typedef ConstVectorView<T,unitA> const_unitview_type;
        typedef ConstVectorView<T,conjA> const_conjugate_type;
        typedef ConstSmallVectorView<T,UNKNOWN,-1,nonunitAc> 
            const_reverse_type;
        typedef ConstSmallVectorView<real_type,UNKNOWN,twoS,twosAc> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstVectorView<real_type,unitA> const_flatten_type;
        typedef ConstVectorView<T,nonconjA> const_nonconj_type;
        typedef VectorView<T,A> nonconst_type;

        typedef CVIt<T,1,false> const_iterator;
        typedef CVIt<T,-1,false> const_reverse_iterator;

        typedef T& reference;

        typedef VectorView<T,A> subvector_type;
        typedef VectorView<T,nonunitA> subvector_step_type;
        typedef VectorView<T,A> view_type;
        typedef VectorView<T,cstyleA> cview_type;
        typedef VectorView<T,fstyleA> fview_type;
        typedef VectorView<T> xview_type;
        typedef VectorView<T,unitA> unitview_type;
        typedef VectorView<T,conjA> conjugate_type;
        typedef SmallVectorView<T,UNKNOWN,-1,nonunitAc> reverse_type;
        typedef SmallVectorView<real_type,UNKNOWN,twoS,twosAc> realpart_type;
        typedef realpart_type imagpart_type;
        typedef VectorView<real_type,A> flatten_type;
        typedef VectorView<T,A> nonconj_type;

        typedef VIt<T,1,false> iterator;
        typedef VIt<T,-1,false> reverse_iterator;
    };

    template <class T, int A>
    class Vector : 
        public BaseVector_Mutable<Vector<T,A> >
    {
    public:

        typedef Vector<T,A> type;
        typedef BaseVector_Mutable<type> base_mut;

        enum { _size = Traits<type>::_size };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { isreal = Traits<type>::isreal };
        enum { iscomplex = Traits<type>::iscomplex };
        enum { _step = Traits<type>::_step };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _unit = Traits<type>::_unit };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        explicit Vector(size_t n=0) : itssize(n), itsv(n)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVAssert(n>=0);
#ifdef TMV_DEBUG
            this->setAllTo(T(888));
#endif
        }

        Vector(size_t n, T val) : itssize(n), itsv(n)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVAssert(n>=0);
            this->setAllTo(val);
        }

        Vector(size_t n, const T* v2) : itssize(n), itsv(n)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVAssert(n>=0);
            VectorViewOf(v2,n).newAssignTo(*this);
        }

        explicit Vector(const std::vector<T>& v2) : 
            itssize(v2.size()), itsv(itssize)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVAssert(itssize>=0);
            VectorViewOf(&v2[0],itssize).newAssignTo(*this);
        }

        Vector(const type& v2) : itssize(v2.size()), itsv(itssize)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVAssert(itssize>=0);
            v2.newAssignTo(*this);
        }

        template <class V2>
        Vector(const BaseVector<V2>& v2) :
            itssize(v2.size()), itsv(itssize)
        {
            TMVStaticAssert(Traits<type>::okA);
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

    //
    // ConstVectorView
    //

    template <class T, int A0>
    struct Traits<ConstVectorView<T,A0> >
    {
        enum { A = (A0 & ~CheckAlias) };
        enum { okA = (
                Attrib<A>::vectoronly &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef ConstVectorView<T,A0> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        enum { copyA = Attrib<A>::fort ? FortranStyle : CStyle };
        typedef Vector<T,copyA> copy_type;

        enum { _size = UNKNOWN };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _step = Attrib<A>::unit ? 1 : UNKNOWN };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = !Attrib<A>::noalias };
        enum { _unit = Attrib<A>::unit };
        enum { negS = IntTraits<_step>::negS };
        enum { twoS = isreal ? int(_step) : IntTraits<_step>::twoS };

        enum { unitA = A | Unit };
        enum { nonunitA = A & ~Unit };
        enum { conjA = iscomplex ? (A ^ Conj) : int(A) };
        enum { nonconjA = A & ~Conj };
        enum { cstyleA = A & ~FortranStyle };
        enum { fstyleA = A | FortranStyle };
        enum { flatA = isreal ? int(A) : (unitA & ~Conj) };
        enum { twosA = isreal ? int(A) : (A & ~Conj & ~Unit) };
        enum { Ac = _checkalias ? (A | CheckAlias) : (A & ~NoAlias) };
        enum { nonunitAc = Ac & ~Unit };
        enum { twosAc = isreal ? int(Ac) : (Ac & ~Conj & ~Unit) };

        typedef ConstVectorView<T,A> const_subvector_type;
        typedef ConstVectorView<T,nonunitA> const_subvector_step_type;
        typedef ConstVectorView<T,A> const_view_type;
        typedef ConstVectorView<T,cstyleA> const_cview_type;
        typedef ConstVectorView<T,fstyleA> const_fview_type;
        typedef ConstVectorView<T,(_conj ? Conj : NonConj)> const_xview_type;
        typedef ConstVectorView<T,unitA> const_unitview_type;
        typedef ConstVectorView<T,conjA> const_conjugate_type;
        typedef typename TypeSelect<_unit ,
                ConstSmallVectorView<T,UNKNOWN,negS,nonunitAc> ,
                ConstVectorView<T,nonunitA> >::type const_reverse_type;
        typedef typename TypeSelect< (iscomplex && _unit) ,
                ConstSmallVectorView<real_type,UNKNOWN,twoS,twosAc> ,
                ConstVectorView<real_type,twosA> >::type const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstVectorView<real_type,flatA> const_flatten_type;
        typedef ConstVectorView<T,nonconjA> const_nonconj_type;
        typedef VectorView<T,A> nonconst_type;

        typedef CVIt<T,_step,_conj> const_iterator;
        typedef CVIt<T,negS,_conj> const_reverse_iterator;
    };

    template <class T, int A>
    class ConstVectorView : 
        public BaseVector_Calc<ConstVectorView<T,A> >
    {
    public:

        typedef ConstVectorView<T,A> type;

        enum { _size = Traits<type>::_size };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { isreal = Traits<type>::isreal };
        enum { iscomplex = Traits<type>::iscomplex };
        enum { _step = Traits<type>::_step };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _unit = Traits<type>::_unit };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        ConstVectorView(const T* v, size_t n, int s) : 
            itsv(v), itssize(n), itsstep(s) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        ConstVectorView(const T* v, size_t n) : 
            itsv(v), itssize(n), itsstep(_step) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_step != UNKNOWN);
        }

        ConstVectorView(const type& v2) : 
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step())
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        ConstVectorView(const VectorView<T,A>& v2) : 
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step())
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        ConstVectorView(const ConstVectorView<T,A2>& v2) :
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int A2>
        ConstVectorView(const VectorView<T,A2>& v2) :
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int N2, int S2, int A2>
        ConstVectorView(const ConstSmallVectorView<T,N2,S2,A2>& v2) :
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int N2, int S2, int A2>
        ConstVectorView(const SmallVectorView<T,N2,S2,A2>& v2) :
            itsv(v2.cptr()), itssize(v2.size()), itsstep(v2.step()) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

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

        T cref(int i) const  { return DoConj<_conj>(itsv[i*step()]); }

        size_t size() const { return itssize; }
        int nElements() const { return itssize; }
        int step() const { return itsstep; }
        bool isconj() const { return _conj; }

    protected :

        const T* itsv;
        const size_t itssize;
        const CheckedInt<_step> itsstep;

    }; // ConstVectorView


    //
    // VectorView
    //

    template <class T, int A0>
    struct Traits<VectorView<T,A0> >
    {
        enum { A = (A0 & ~CheckAlias) };
        enum { okA = (
                Attrib<A>::vectoronly &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef VectorView<T,A0> type;
        typedef const ConstVectorView<T,A> calc_type;
        typedef calc_type eval_type;
        enum { copyA = Attrib<A>::fort ? FortranStyle : CStyle };
        typedef Vector<T,copyA> copy_type;

        enum { _size = UNKNOWN };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _step = Attrib<A>::unit ? 1 : UNKNOWN };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = !Attrib<A>::noalias };
        enum { _unit = Attrib<A>::unit };
        enum { negS = IntTraits<_step>::negS };
        enum { twoS = isreal ? int(_step) : IntTraits<_step>::twoS };

        enum { unitA = A | Unit };
        enum { nonunitA = A & ~Unit };
        enum { conjA = iscomplex ? (A ^ Conj) : int(A) };
        enum { nonconjA = A & ~Conj };
        enum { cstyleA = A & ~FortranStyle };
        enum { fstyleA = A | FortranStyle };
        enum { flatA = isreal ? int(A) : (unitA & ~Conj) };
        enum { twosA = isreal ? int(A) : (A & ~Conj & ~Unit) };
        enum { Ac = _checkalias ? (A | CheckAlias) : (A & ~NoAlias) };
        enum { nonunitAc = Ac & ~Unit };
        enum { twosAc = isreal ? int(Ac) : (Ac & ~Conj & ~Unit) };

        typedef ConstVectorView<T,A> const_subvector_type;
        typedef ConstVectorView<T,nonunitA> const_subvector_step_type;
        typedef ConstVectorView<T,A> const_view_type;
        typedef ConstVectorView<T,cstyleA> const_cview_type;
        typedef ConstVectorView<T,fstyleA> const_fview_type;
        typedef ConstVectorView<T,(_conj ? Conj : NonConj)> const_xview_type;
        typedef ConstVectorView<T,unitA> const_unitview_type;
        typedef ConstVectorView<T,conjA> const_conjugate_type;
        typedef typename TypeSelect<_unit ,
                ConstSmallVectorView<T,UNKNOWN,negS,nonunitAc> ,
                ConstVectorView<T,nonunitA> >::type const_reverse_type;
        typedef typename TypeSelect< (iscomplex && _unit) ,
                ConstSmallVectorView<real_type,UNKNOWN,twoS,twosAc> ,
                ConstVectorView<real_type,twosA> >::type const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstVectorView<real_type,flatA> const_flatten_type;
        typedef ConstVectorView<T,nonconjA> const_nonconj_type;
        typedef VectorView<T,A> nonconst_type;

        typedef CVIt<T,_step,_conj> const_iterator;
        typedef CVIt<T,negS,_conj> const_reverse_iterator;

        typedef typename AuxRef<T,_conj>::reference reference;

        typedef VectorView<T,A> subvector_type;
        typedef VectorView<T,nonunitA> subvector_step_type;
        typedef VectorView<T,A> view_type;
        typedef VectorView<T,cstyleA> cview_type;
        typedef VectorView<T,fstyleA> fview_type;
        typedef VectorView<T,(_conj ? Conj : NonConj)> xview_type;
        typedef VectorView<T,unitA> unitview_type;
        typedef VectorView<T,conjA> conjugate_type;
        typedef typename TypeSelect<_unit ,
                SmallVectorView<T,UNKNOWN,negS,nonunitAc> ,
                VectorView<T,nonunitA> >::type reverse_type;
        typedef typename TypeSelect< (iscomplex && _unit) ,
                SmallVectorView<real_type,UNKNOWN,twoS,twosAc> ,
                VectorView<real_type,twosA> >::type realpart_type;
        typedef realpart_type imagpart_type;
        typedef VectorView<real_type,flatA> flatten_type;
        typedef VectorView<T,nonconjA> nonconj_type;

        typedef VIt<T,_step,_conj> iterator;
        typedef VIt<T,negS,_conj> reverse_iterator;
    };

    template <class T, int A>
    class VectorView : 
        public BaseVector_Mutable<VectorView<T,A> >
    {
    public:

        typedef VectorView<T,A> type;
        typedef BaseVector_Mutable<type> base_mut;
        typedef typename base_mut::reference reference;

        enum { _size = Traits<type>::_size };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { isreal = Traits<type>::isreal };
        enum { iscomplex = Traits<type>::iscomplex };
        enum { _step = Traits<type>::_step };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _unit = Traits<type>::_unit };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        VectorView(T* v, size_t n, int s) : 
            itsv(v), itssize(n), itsstep(s) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        VectorView(T* v, size_t n) : 
            itsv(v), itssize(n), itsstep(_step) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_step != UNKNOWN);
        }

        VectorView(const type& v2) : 
            itsv(v2.itsv), itssize(v2.size()), itsstep(v2.step())
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        VectorView(VectorView<T,A2> v2) :
            itsv(v2.ptr()), itssize(v2.size()), itsstep(v2.step())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int N2, int S2, int A2>
        VectorView(SmallVectorView<T,N2,S2,A2> v2) :
            itsv(v2.ptr()), itssize(v2.size()), itsstep(v2.step())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

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

        T cref(int i) const  { return DoConj<_conj>(itsv[i*step()]); }
        reference ref(int i) { return reference(itsv[i*step()]); }

        size_t size() const { return itssize; }
        int nElements() const { return itssize; }
        int step() const { return itsstep; }
        bool isconj() const { return _conj; }

    protected :

        T* itsv;
        const size_t itssize;
        const CheckedInt<_step> itsstep;

    }; // VectorView


    //
    // Special Creators: 
    //   VectorViewOf(v,size,step=1) = VectorView of v 
    //

    // VectorView of raw memory:
    template <class T>
    static inline VectorView<T,Unit> VectorViewOf(T* v, size_t size)
    { return VectorView<T,Unit>(v,size); }

    template <class T>
    static inline VectorView<T,NonUnit> VectorViewOf(
        T* v, size_t size, int step)
    { return VectorView<T,NonUnit>(v,size,step); }


    //
    // Swap
    //

    template <class T, int A>
    static inline void Swap(Vector<T,A>& v1, Vector<T,A>& v2)
    { v1.swapWith(v2); }
    template <class V, class T, int A>
    static inline void Swap(BaseVector_Mutable<V>& v1, VectorView<T,A> v2)
    { DoSwap(v1,v2); }
    template <class V, class T, int A>
    static inline void Swap(VectorView<T,A> v1, BaseVector_Mutable<V>& v2)
    { DoSwap(v1,v2); }
    template <class T, int A1, int A2>
    static inline void Swap(
        VectorView<T,A1> v1, VectorView<T,A2> v2)
    { DoSwap(v1,v2); }


    //
    // Conjugate
    //
    
    template <class T, int A>
    static inline typename Vector<T,A>::conjugate_type Conjugate(
        Vector<T,A>& v)
    { return v.conjugate(); }
    template <class T, int A>
    static inline typename VectorView<T,A>::conjugate_type Conjugate(
        VectorView<T,A> v)
    { return v.conjugate(); }


    //
    // TMV_Text 
    //

#ifdef TMV_DEBUG
    template <class T, int A>
    static inline std::string TMV_Text(const Vector<T,A>& v)
    {
        std::ostringstream s;
        s << "Vector<"<<TMV_Text(T());
        s << ","<<Attrib<A>::vtext()<<">";
        s << "("<<v.size()<<","<<v.step()<<")";
        return s.str();
    }

    template <class T, int A>
    static inline std::string TMV_Text(const ConstVectorView<T,A>& v)
    {
        std::ostringstream s;
        s << "ConstVectorView<"<<TMV_Text(T());
        s << ","<<Attrib<A>::vtext()<<">";
        s << "("<<v.size()<<","<<v.step()<<")";
        return s.str();
    }

    template <class T, int A>
    static inline std::string TMV_Text(const VectorView<T,A>& v)
    {
        std::ostringstream s;
        s << "VectorView<"<<TMV_Text(T());
        s << ","<<Attrib<A>::vtext()<<">";
        s << "("<<v.size()<<","<<v.step()<<")";
        return s.str();
    }
#endif

} // namespace tmv

#endif
