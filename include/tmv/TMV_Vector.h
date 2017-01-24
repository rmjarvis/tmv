///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2016                                                 //
// All rights reserved                                                       //
//                                                                           //
// The project is hosted at https://code.google.com/p/tmv-cpp/               //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis17 [at] gmail.                                                 //
//                                                                           //
// This software is licensed under a FreeBSD license.  The file              //
// TMV_LICENSE should have bee included with this distribution.              //
// It not, you can get a copy from https://code.google.com/p/tmv-cpp/.       //
//                                                                           //
// Essentially, you can use this software however you want provided that     //
// you include the TMV_LICENSE file in any distribution that uses it.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


//-----------------------------------------------------------------------------
//
// This file defines the TMV Vector class.
//
// A Vector is a mathematical vector whose size is not necessarily
// known at compile time.  (c.f. SmallVector in TMV_SmallVector.h.)
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
// Constructors:
//
//    explicit Vector<T,A>(int n)
//        Makes a Vector of size n with _uninitialized_ values
//
//    Vector<T,A>(int n, T x)
//        Makes a Vector of size n with all values = x
//
//
// Special Constructors
//
//    makeBasisVector(int n, int i)
//        Makes a Vector whose elements are all 0, except v(i) = 1
//
//    VectorViewOf(T* vv, int n)
//    VectorViewOf(const T* vv, int n)
//        Makes a Vector View (see below) which refers to the exact
//        elements of vv, not copying them to new storage.
//        The first one returns a VectorView, the second a ConstVectorView.
//
// Access Functions
//
//    int size() const
//        Returns the size of the Vector
//
//    T& operator[](int i)
//    T& operator()(int i)
//    T operator[](int i) const
//    T operator()(int i) const
//        Return the ith element of the Vector
//
//    Vector<T,A>::iterator begin()
//    Vector<T,A>::iterator end()
//    Vector<T,A>::reverse_iterator rbegin()
//    Vector<T,A>::reverse_iterator rend()
//    Vector<T,A>::const_iterator begin() const
//    Vector<T,A>::const_iterator end() const
//    Vector<T,A>::const_reverse_iterator rbegin() const
//    Vector<T,A>::const_reverse_iterator rend() const
//        Return iterators that work in the usual way to traverse the Vector
//
// Modifying Functions
//
//    Vector& setZero()
//        Sets all elements to 0
//
//    Vector& clip(RT thresh)
//        Set to 0 all elements whose absolute values is < thresh
//
//    Vector& setAllTo(T x)
//        Sets all elements to x
//        We don't call this v = x, since it doesn't really make sense to
//        think of v as being 'equal' to a scalar.
//        Hence the slightly verbose function name setAllTo.
//
//    Vector& addToAll(T x)
//        Add x to each element of v
//
//    Vector& conjugateSelf()
//        Sets all elements to its conjugate
//
//    Vector& makeBasis(int i)
//        Set all elements to 0, except v(i) = 1
//
//    Vector& swap(int i1, int i2)
//        Swap elements v(i1) and v(i2)
//
//    Vector& permute(const int* p)
//    Vector& permute(const int* p, int i1, int i2)
//        Perform a series of swaps (0,p[0]), (1,p[1])...
//        In the second case, only do (i1,p[i1])...(i2-1,p[i2-1])
//    Vector& reversePermute(const int* p)
//    Vector& reversePermute(const int* p, int i1, int i2)
//        The same but perform the swaps in reverse order.
//
//    Swap(v1,v2)
//        Swap vectors v1 and v2
//
//    Vector& reverseSelf()
//        Reverse the order of the elements of v
//
//    sort(AD, COMP)
//    sort(int* P, AD, COMP)
//    sort(Permutation& P, AD, COMP)
//        Sorts the vector, returning the swaps required in P.
//        If you do not care about P, you may omit the P parameter.
//        AD = Ascend or Descend (Ascend=default)
//        COMP = RealComp, AbsComp, ImagComp, or ArgComp (RealComp=default)
//
// VectorViews:
//
//    A VectorView<T> object refers to the elements of a regular Vector
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
//    Below, VectorView is written for both of these, and I don't bother to
//    indicate the template modifiers.  In general, a ConstVectorView is
//    returned from either a ConstVectorView object, or a const Vector object.
//    A (mutable) VectorView is returned from either a VectorView object or
//    a (non-const) Vector object.
//
//    VectorView subVector(int i1, int i2, int istep=1)
//        This member function will return a subvector view which refers
//        to the same physical elements as the original.
//        i1 is the first element in the subvector.
//        i2 is the usual "one past the end" of the subvector.
//        istep is an optional step size.
//        Thus, if you have a vector v of length 10, and you want to
//        set all the even elements to 0, you could write:
//        v.subVector(0,10,2).Zero();
//        And then to output the first 4 elements of v, you could write:
//        std::cout << v.subVector(0,4);
//
//    VectorView reverse()
//        Returns a subvector whose elements are the same as v but
//        in the reverse order
//
//    VectorView conjugate()
//        Returns the conjugate of a Vector (as a view, so it still points
//        to the same physical elements, but modifying this will set the
//        actual elements to the conjugate of what you set.
//
//    VectorView view()
//        Returns a view of a Vector.
//
//    VectorView cView()
//    VectorView fView()
//        Returns a view of a Vector in either CStyle or FortranStyle
//        form, respectively.
//
//    VectorView realPart()
//    VectorView imagPart()
//        Returns a view of the real/imag elements of a complex Vector.
//    VectorView flatten()
//        For complex vectors with unit step, returns a real vector view to
//        both the real and imag elements in a single view.
//
// Functions of Vectors:
//        (These are all both member functions and functions of a Vector,
//         so Norm(v) and v.norm() for example are equivalent)
//
//    v.norm() or v.norm2() or Norm(v) or Norm2(v)
//        Returns the 2-norm of a Vector
//        = sqrt( sum_i |v_i|^2 )
//
//    v.normSq() or NormSq(v)
//    v.normSq(RT scale=1.)
//        Returns the square of norm()
//        In the method version, you can provide an optional scale, in
//        which case the output is equal to normSq(scale*v).
//
//    v.norm1() or Norm1(v)
//        Returns the 1-norm of a Vector
//        = sum_i |v_i|
//
//    v.norminf() or NormInf(v)
//        Returns the infinity-norm of a Vector
//        = max_i |v_i|
//
//    v.sumElements() or SumElements(v)
//        Returns the sum of all elements in the Vector
//
//    v.sumAbsElements() or SumAbsElements(v)
//        Returns the sum of absolute values of elements in the Vector
//
//    v.sumAbs2Elements() or SumAbs2Elements(v)
//        Returns the sum of absolute values of elements in the Vector
//        using |real(v(i))| + |imag(v(i))| for complex vectors.
//
//    v.maxElement() or MaxElement(v)
//    v.maxElement(int* imax = 0)
//        Returns the maximum value of any element in the Vector
//        In the method version, you can provide an optional argument imax.
//        On return, *imax holds the index of element with the maximum value.
//        The parameter imax can be omitted if it is not desired.
//        As "max" doesn't make sense for complex values, for these
//        we use just the real components.
//
//    v.minElement() or MinElement(v)
//    v.minElement(int* imin = 0)
//        Returns the minimum value of any element in the Vector
//        On return, (in the method version) *imin holds the index of this
//        element.  Again, the parameter imin can be omitted.
//
//    v.maxAbsElement() or MaxAbsElement(v)
//    v.maxAbsElement(int* imax=0)
//        The same as maxElement, except absolute values are used
//
//    v.minAbsElement() or MinAbsElement(v)
//    v.minAbsElement(int* imin=0)
//        The same as minElement, except absolute values are used
//
//    v.maxAbs2Element() or MaxAbs2Element(v)
//    v.maxAbs2Element(int* imax=0)
//    v.minAbs2Element() or MinAbs2Element(v)
//    v.minAbs2Element(int* imin=0)
//        For real v, these are identical to the regular AbsElement routines.
//        But for complex v, instead of calculating |v(i)| for each
//        element, it calculates |v(i).real()| + |v(i).imag()|.
//        This is significantly faster than the full complex abs,
//        and for many purposes, the answer is equally usable.
//        e.g. for scaling a vector by it's "largest" value, this is
//        often just a useful a definition of "largest".
//
// Operators:
//        Here we use v for a Vector and x for a Scalar.
//
//        You can also mix real and complex vectors of the same
//        underlying type.  eg. Vector<double> and Vector<complex<double> >
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
//       These all behave in the logical way for dealing with vectors.
//
//
// I/O:
//
//    os << v
//        Writes v to ostream os in the following format:
//          n ( v[0] v[1] v[2] ... v[n] )
//
//    v.write(ostream& os, TMV_RealType(T) minnonzero)
//        Write v to os as above, but if |v[i]| < minnonzero, write 0 instead
//
//    is >> v
//        Reads v from istream is in the same format
//
//


#ifndef TMV_Vector_H
#define TMV_Vector_H

#include "tmv/TMV_BaseVector.h"
#include "tmv/TMV_ListInit.h"
#include "tmv/TMV_Array.h"
#include "tmv/TMV_IOStyle.h"
#include "tmv/TMV_VIt.h"

namespace tmv {

    template <typename T, typename T1>
    inline void Copy(const GenVector<T1>& v1, VectorView<T> v2);

    class Permutation;

    template <typename T>
    class GenVector : public AssignableToVector<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef T value_type;
        typedef RT real_type;
        typedef CT complex_type;
        typedef GenVector<T> type;
        typedef Vector<T> copy_type;
        typedef ConstVectorView<T> const_view_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_reverse_type;
        typedef ConstVectorView<T,CStyle> const_cview_type;
        typedef ConstVectorView<T,FortranStyle> const_fview_type;
        typedef ConstVectorView<RT> const_real_type;
        typedef VectorView<T> nonconst_type;
        typedef typename RefHelper<T>::const_iterator const_iterator;
        typedef typename RefHelper<T>::const_iterator const_reverse_iterator;

        //
        // Constructor
        //

        inline GenVector() {}
        inline GenVector(const GenVector<T>&) {}
        virtual inline ~GenVector() {}


        //
        // Access Functions
        //

        using AssignableToVector<T>::size;
        void assignToV(VectorView<RT> rhs) const
        {
            TMVAssert(rhs.size() == size());
            TMVAssert(isReal(T()));
            if (!isSameAs(rhs)) Copy(*this,rhs);
        }
        void assignToV(VectorView<CT> rhs) const
        {
            TMVAssert(rhs.size() == size());
            if (!isSameAs(rhs)) Copy(*this,rhs);
        }

        inline T operator()(ptrdiff_t i) const
        {
            TMVAssert(i>=0 && i<size());
            return cref(i);
        }
        inline T operator[](ptrdiff_t i) const
        { return operator()(i); }

        inline const_iterator begin() const
        { return RefHelper<T>::makeIter(cptr(),step(),ct()); }
        inline const_iterator end() const
        { return begin() + size(); }
        inline const_reverse_iterator rbegin() const
        {
            return RefHelper<T>::makeIter(
                cptr()+step()*(size()-1),-step(),ct());
        }
        inline const_reverse_iterator rend() const
        { return rbegin() + size(); }

        template <typename T2>
        inline bool isSameAs(const GenVector<T2>&) const
        { return false; }

        inline bool isSameAs(const GenVector<T>& v2) const
        {
            return (
                this == &v2 ||
                ( cptr()==v2.cptr() && size()==v2.size() &&
                  step()==v2.step() && ct()==v2.ct() ) );
        }

        //
        // subVector
        //

        bool hasSubVector(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const;

        inline const_view_type cSubVector(ptrdiff_t i1, ptrdiff_t i2) const
        { return const_view_type(cptr()+i1*step(),i2-i1,step(),ct()); }

        inline const_view_type subVector(ptrdiff_t i1, ptrdiff_t i2) const
        {
            TMVAssert(hasSubVector(i1,i2,1));
            return cSubVector(i1,i2);
        }

        inline const_view_type cSubVector(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const
        {
            return const_view_type(
                cptr()+i1*step(),(i2-i1)/istep,istep*step(),ct());
        }

        inline const_view_type subVector(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const
        {
            TMVAssert(hasSubVector(i1,i2,istep));
            return cSubVector(i1,i2,istep);
        }

        inline const_view_type reverse() const
        {
            return const_view_type(
                cptr()+(size()-1)*step(),size(),-step(),ct());
        }

        inline const_view_type view() const
        { return const_view_type(cptr(),size(),step(),ct()); }

        inline const_cview_type cView() const
        { return view(); }

        inline const_fview_type fView() const
        { return view(); }

        inline const_view_type conjugate() const
        { return const_view_type(cptr(),size(),step(),TMV_ConjOf(T,ct())); }

        inline const_real_type realPart() const
        {
            return const_real_type(
                reinterpret_cast<const RT*>( cptr()),
                size(), isReal(T()) ? step() : 2*step(), NonConj);
        }

        inline const_real_type imagPart() const
        {
            TMVAssert(isComplex(T()));
            return const_real_type(
                reinterpret_cast<const RT*>( cptr())+1,
                size(), 2*step(), NonConj);
        }

        inline const_real_type flatten() const
        {
            TMVAssert(isComplex(T()));
            TMVAssert(step() == 1);
            return const_real_type(
                reinterpret_cast<const RT*>(cptr()),
                isReal(T()) ? size() : 2*size(), 1, NonConj);
        }

        inline nonconst_type nonConst() const
        {
            return nonconst_type(
                const_cast<T*>(cptr()),
                size(),step(),ct() TMV_FIRSTLAST1(cptr(),end().getP()));
        }

        //
        // Functions of Vector
        //

        inline RT norm() const
        { return norm2(); }

        RT normSq(const RT scale = RT(1)) const;

        inline RT norm1() const
        { return sumAbsElements(); }

        RT norm2() const;

        inline RT normInf() const
        { return size() > 0 ? maxAbsElement() : RT(0); }

        T sumElements() const;

        RT sumAbsElements() const;

        RT sumAbs2Elements() const;

        // Note: if these are just minElement, etc. rather than
        // doMinElement, etc. then the template instantiations
        // get confused by the versions with a template argument.
        T doMinElement(ptrdiff_t* iminout=0) const;

        T doMaxElement(ptrdiff_t* imaxout=0) const;

        RT doMinAbsElement(ptrdiff_t* iminout=0) const;

        RT doMaxAbsElement(ptrdiff_t* imaxout=0) const;

        RT doMinAbs2Element(ptrdiff_t* iminout=0) const;

        RT doMaxAbs2Element(ptrdiff_t* imaxout=0) const;

        inline T minElement(ptrdiff_t* iminout=0) const
        { return doMinElement(iminout); }

        inline T maxElement(ptrdiff_t* imaxout=0) const
        { return doMaxElement(imaxout); }

        inline RT minAbsElement(ptrdiff_t* iminout=0) const
        { return doMinAbsElement(iminout); }

        inline RT maxAbsElement(ptrdiff_t* imaxout=0) const
        { return doMaxAbsElement(imaxout); }

        inline RT minAbs2Element(ptrdiff_t* iminout=0) const
        { return doMinAbs2Element(iminout); }

        inline RT maxAbs2Element(ptrdiff_t* imaxout=0) const
        { return doMaxAbs2Element(imaxout); }

        // Also allow other int types in case ptrdiff_t is not int:
        template <class INT>
        inline T minElement(INT* iminout) const
        { ptrdiff_t i; T temp=minElement(&i); *iminout=i; return temp; }

        template <class INT>
        inline T maxElement(INT* imaxout) const
        { ptrdiff_t i; T temp=maxElement(&i); *imaxout=i; return temp; }

        template <class INT>
        inline RT minAbsElement(INT* iminout) const
        { ptrdiff_t i; RT temp=minAbsElement(&i); *iminout=i; return temp; }

        template <class INT>
        inline RT maxAbsElement(INT* imaxout) const
        { ptrdiff_t i; RT temp=maxAbsElement(&i); *imaxout=i; return temp; }

        template <class INT>
        inline RT minAbs2Element(INT* iminout) const
        { ptrdiff_t i; RT temp=minAbs2Element(&i); *iminout=i; return temp; }

        template <class INT>
        inline RT maxAbs2Element(INT* imaxout) const
        { ptrdiff_t i; RT temp=maxAbs2Element(&i); *imaxout=i; return temp; }


        //
        // I/O
        //

        void write(const TMV_Writer& writer) const;

        virtual const T* cptr() const =0;
        virtual ptrdiff_t step() const =0;
        virtual ConjType ct() const =0;
        inline bool isconj() const
        { return (isComplex(T()) && ct()==Conj); }

        virtual T cref(ptrdiff_t i) const;

    private:

        type& operator=(const GenVector<T>&);

    }; // GenVector

    template <typename T, int A>
    class ConstVectorView : public GenVector<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef ConstVectorView<T,A> type;
        typedef ConstVectorView<RT,A> const_real_type;

        inline ConstVectorView(const ConstVectorView<T,A>& rhs) :
            _v(rhs._v), _size(rhs._size), _step(rhs._step),
            _ct(rhs._ct)
        { TMVAssert(Attrib<A>::vectorok);  }
        inline ConstVectorView(const GenVector<T>& rhs) :
            _v(rhs.cptr()), _size(rhs.size()), _step(rhs.step()),
            _ct(rhs.ct())
        { TMVAssert(Attrib<A>::vectorok);  }
        inline ConstVectorView(
            const T* v, ptrdiff_t size, ptrdiff_t step, ConjType ct) :
            _v(v), _size(size), _step(step), _ct(ct)
        { TMVAssert(Attrib<A>::vectorok);  }
        virtual inline ~ConstVectorView()
        {
#ifdef TMV_EXTRA_DEBUG
            const_cast<const T*&>(_v) = 0;
#endif
        }

        inline ptrdiff_t size() const { return _size; }
        inline const T* cptr() const { return _v; }
        inline ptrdiff_t step() const { return _step; }
        inline ConjType ct() const { return _ct; }

    private:

        const T*const _v;
        const ptrdiff_t _size;
        const ptrdiff_t _step;
        const ConjType _ct;

        type& operator=(const ConstVectorView<T,A>&);

    }; // ConstVectorView

    template <typename T>
    class ConstVectorView<T,FortranStyle> : public ConstVectorView<T,CStyle>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef ConstVectorView<T,FortranStyle> type;
        typedef ConstVectorView<T,FortranStyle> const_view_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_reverse_type;
        typedef ConstVectorView<T,CStyle> const_cview_type;
        typedef ConstVectorView<T,FortranStyle> const_fview_type;
        typedef ConstVectorView<RT,FortranStyle> const_real_type;

        inline ConstVectorView(const ConstVectorView<T,FortranStyle>& rhs) :
            ConstVectorView<T,CStyle>(rhs) {}
        inline ConstVectorView(const ConstVectorView<T,CStyle>& rhs) :
            ConstVectorView<T,CStyle>(rhs) {}
        inline ConstVectorView(const GenVector<T>& rhs) :
            ConstVectorView<T,CStyle>(rhs) {}
        inline ConstVectorView(const T* inv, ptrdiff_t insize, ptrdiff_t instep,
                               ConjType inct) :
            ConstVectorView<T,CStyle>(inv,insize,instep,inct) {}
        virtual inline ~ConstVectorView() {}

        //
        // Access Functions
        //

        inline T operator()(ptrdiff_t i) const
        {
            TMVAssert(i>0 && i<=this->size());
            return GenVector<T>::cref(i-1);
        }
        inline T operator[](ptrdiff_t i) const
        { return operator()(i); }

        //
        // SubVector
        //

        bool hasSubVector(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const;

        inline const_view_type subVector(ptrdiff_t i1, ptrdiff_t i2) const
        {
            TMVAssert(hasSubVector(i1,i2,1));
            return GenVector<T>::cSubVector(i1-1,i2);
        }

        inline const_view_type subVector(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const
        {
            TMVAssert(hasSubVector(i1,i2,istep));
            return GenVector<T>::cSubVector(i1-1,i2-1+istep,istep);
        }

        inline const_view_type reverse() const
        { return GenVector<T>::reverse(); }

        inline const_view_type view() const
        { return GenVector<T>::view(); }

        inline const_cview_type cView() const
        { return view(); }

        inline const_fview_type fView() const
        { return view(); }

        inline const_view_type conjugate() const
        { return GenVector<T>::conjugate(); }

        inline const_real_type realPart() const
        { return GenVector<T>::realPart(); }

        inline const_real_type imagPart() const
        { return GenVector<T>::imagPart(); }

        inline const_real_type flatten() const
        { return GenVector<T>::flatten(); }

        //
        // Functions of Vector
        //

        inline T minElement(ptrdiff_t* iminout=0) const
        {
            T temp = GenVector<T>::minElement(iminout);
            if (iminout) ++(*iminout);
            return temp;
        }

        inline T maxElement(ptrdiff_t* imaxout=0) const
        {
            T temp = GenVector<T>::maxElement(imaxout);
            if (imaxout) ++(*imaxout);
            return temp;
        }

        inline RT minAbsElement(ptrdiff_t* iminout=0) const
        {
            RT temp = GenVector<T>::minAbsElement(iminout);
            if (iminout) ++(*iminout);
            return temp;
        }

        inline RT maxAbsElement(ptrdiff_t* imaxout=0) const
        {
            RT temp = GenVector<T>::maxAbsElement(imaxout);
            if (imaxout) ++(*imaxout);
            return temp;
        }

        inline RT minAbs2Element(ptrdiff_t* iminout=0) const
        {
            RT temp = GenVector<T>::minAbs2Element(iminout);
            if (iminout) ++(*iminout);
            return temp;
        }

        inline RT maxAbs2Element(ptrdiff_t* imaxout=0) const
        {
            RT temp = GenVector<T>::maxAbs2Element(imaxout);
            if (imaxout) ++(*imaxout);
            return temp;
        }

        // Also allow other int types in case ptrdiff_t is not int:
        template <class INT>
        inline T minElement(INT* iminout) const
        { ptrdiff_t i; T temp=minElement(&i); *iminout=i; return temp; }

        template <class INT>
        inline T maxElement(INT* imaxout) const
        { ptrdiff_t i; T temp=maxElement(&i); *imaxout=i; return temp; }

        template <class INT>
        inline RT minAbsElement(INT* iminout) const
        { ptrdiff_t i; RT temp=minAbsElement(&i); *iminout=i; return temp; }

        template <class INT>
        inline RT maxAbsElement(INT* imaxout) const
        { ptrdiff_t i; RT temp=maxAbsElement(&i); *imaxout=i; return temp; }

        template <class INT>
        inline RT minAbs2Element(INT* iminout) const
        { ptrdiff_t i; RT temp=minAbs2Element(&i); *iminout=i; return temp; }

        template <class INT>
        inline RT maxAbs2Element(INT* imaxout) const
        { ptrdiff_t i; RT temp=maxAbs2Element(&i); *imaxout=i; return temp; }



    private:

        const type& operator=(const type&);

    }; // FortranStyle ConstVectorView

    template <typename T, int A>
    class VectorView : public GenVector<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef VectorView<T,A> type;
        typedef GenVector<T> base;
        typedef VectorView<T,A> view_type;
        typedef view_type conjugate_type;
        typedef view_type reverse_type;
        typedef VectorView<T,CStyle> cview_type;
        typedef VectorView<T,FortranStyle> fview_type;
        typedef VectorView<RT,A> real_type;
        typedef ConstVectorView<T,A> const_view_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_reverse_type;
        typedef ConstVectorView<T,CStyle> const_cview_type;
        typedef ConstVectorView<T,FortranStyle> const_fview_type;
        typedef ConstVectorView<RT,A> const_real_type;
        typedef T value_type;
        typedef typename RefHelper<T>::reference reference;
        typedef typename RefHelper<T>::iterator iterator;
        typedef typename RefHelper<T>::iterator reverse_iterator;
        typedef typename RefHelper<T>::const_iterator const_iterator;
        typedef typename RefHelper<T>::const_iterator const_reverse_iterator;

        //
        // Constructors
        //

        inline VectorView(const type& rhs) :
            _v(rhs._v), _size(rhs._size), _step(rhs._step),
            _ct(rhs._ct) TMV_DEFFIRSTLAST(rhs._first,rhs._last)
        { TMVAssert(Attrib<A>::vectorok); }

        inline VectorView(
            T* v, ptrdiff_t size, ptrdiff_t step, ConjType ct
            TMV_PARAMFIRSTLAST(T) ) :
            _v(v), _size(size), _step(step),
            _ct(ct) TMV_DEFFIRSTLAST(_first,_last)
        { TMVAssert(Attrib<A>::vectorok); }

        virtual inline ~VectorView()
        {
            TMV_SETFIRSTLAST(0,0);
#ifdef TMV_EXTRA_DEBUG
            const_cast<T*&>(_v) = 0;
#endif
        }


        //
        // Op =
        //

        inline type& operator=(const type& v2)
        {
            TMVAssert(size() == v2.size());
            v2.assignToV(*this);
            return *this;
        }

        inline type& operator=(const GenVector<RT>& v2)
        {
            TMVAssert(size() == v2.size());
            v2.assignToV(*this);
            return *this;
        }

        inline type& operator=(const GenVector<CT>& v2)
        {
            TMVAssert(size() == v2.size());
            TMVAssert(isComplex(T()));
            v2.assignToV(*this);
            return *this;
        }

        template <typename T2>
        inline type& operator=(const GenVector<T2>& v2)
        {
            TMVAssert(isReal(T2()) || isComplex(T()));
            TMVAssert(size() == v2.size());
            Copy(v2,*this);
            return *this;
        }

        inline type& operator=(const AssignableToVector<RT>& v2)
        {
            TMVAssert(size() == v2.size());
            v2.assignToV(*this);
            return *this;
        }

        inline type& operator=(const AssignableToVector<CT>& v2)
        {
            TMVAssert(size() == v2.size());
            TMVAssert(isComplex(T()));
            v2.assignToV(*this);
            return *this;
        }

        template <typename T2, ptrdiff_t N, int A2>
        inline type& operator=(const SmallVector<T2,N,A2>& v2)
        {
            TMVAssert(size() == v2.size());
            Copy(v2.view(),*this);
            return *this;
        }

        //
        // Access Functions
        //

        inline reference operator()(ptrdiff_t i)
        {
            TMVAssert(i>=0 && i<size());
            return ref(i);
        }
        inline reference operator[](ptrdiff_t i)
        { return operator()(i); }

        inline iterator begin()
        { return RefHelper<T>::makeIter(ptr(),step(),ct() TMV_FIRSTLAST ); }
        inline iterator end()
        { return begin() + size(); }
        inline reverse_iterator rbegin()
        {
            return RefHelper<T>::makeIter(
                ptr()+step()*(size()-1),-step(),ct() TMV_FIRSTLAST );
        }
        inline reverse_iterator rend()
        { return rbegin() + size(); }

        typedef ListAssigner<T,iterator> MyListAssigner;
        MyListAssigner inline operator<<(const T& x)
        { return MyListAssigner(begin(),size(),x); }

        // Repeat const versions
        inline T operator()(ptrdiff_t i) const
        { return base::operator()(i); }
        inline T operator[](ptrdiff_t i) const
        { return base::operator[](i); }
        inline const_iterator begin() const
        { return base::begin(); }
        inline const_iterator end() const
        { return base::end(); }
        inline const_reverse_iterator rbegin() const
        { return base::rbegin(); }
        inline const_reverse_iterator rend() const
        { return base::rend(); }


        //
        // Modifying Functions
        //

        type& setZero();

        type& clip(RT thresh);

        type& setAllTo(const T& x);

        type& addToAll(const T& x);

        type& conjugateSelf();

        type& DoBasis(ptrdiff_t i);
        inline type& makeBasis(ptrdiff_t i)
        { TMVAssert(i>=0 && i<size()); return DoBasis(i); }

        type& DoSwap(ptrdiff_t i1, ptrdiff_t i2);
        inline type& swap(ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(i1>=0 && i1<size());
            TMVAssert(i2>=0 && i2<size());
            return DoSwap(i1,i2);
        }

        type& DoPermute(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2);
        inline type& permute(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(i1>=0 && i1 <= i2 && i2 <= size());
            return DoPermute(p,i1,i2);
        }
        inline type& permute(const ptrdiff_t* p)
        { return DoPermute(p,0,size()); }

        type& DoReversePermute(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2);
        inline type& reversePermute(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(i1>=0 && i1 <= i2 && i2 <= size());
            return DoReversePermute(p,i1,i2);
        }
        inline type& reversePermute(const ptrdiff_t* p)
        { return reversePermute(p,0,size()); }

        type& reverseSelf();

        type& sort(ptrdiff_t* p, ADType ad=Ascend, CompType comp=RealComp);

        inline type& sort(
            Permutation& P, ADType ad=Ascend, CompType comp=RealComp);

        inline type& sort(ADType ad=Ascend, CompType comp=RealComp)
        { sort(0,ad,comp); return *this; }


        //
        // SubVector
        //

        inline view_type cSubVector(ptrdiff_t i1, ptrdiff_t i2)
        {
            return view_type(
                ptr()+i1*step(),(i2-i1),step(),ct() TMV_FIRSTLAST );
        }

        inline view_type subVector(ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(GenVector<T>::hasSubVector(i1,i2,1));
            return cSubVector(i1,i2);
        }

        inline view_type cSubVector(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep)
        {
            return view_type(
                ptr()+i1*step(),(i2-i1)/istep,istep*step(), ct()
                TMV_FIRSTLAST );
        }

        inline view_type subVector(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep)
        {
            TMVAssert(GenVector<T>::hasSubVector(i1,i2,istep));
            return cSubVector(i1,i2,istep);
        }

        inline view_type reverse()
        {
            return view_type(
                ptr()+(size()-1)*step(),size(),-step(), ct()
                TMV_FIRSTLAST );
        }

        inline view_type view()
        { return *this; }

        inline cview_type cView()
        { return view(); }

        inline fview_type fView()
        { return view(); }

        inline view_type conjugate()
        {
            return view_type(
                ptr(),size(),step(), TMV_ConjOf(T,ct()) TMV_FIRSTLAST );
        }

        inline real_type realPart()
        {
            return real_type(
                reinterpret_cast<RT*>(ptr()),
                size(), isReal(T()) ? step() : 2*step(), NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)
                ,reinterpret_cast<const RT*>(_last)
#endif
            );
        }

        inline real_type imagPart()
        {
            TMVAssert(isComplex(T()));
            return real_type(
                reinterpret_cast<RT*>(ptr())+1,size(),2*step(), NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)+1
                ,reinterpret_cast<const RT*>(_last)+1
#endif
            );
        }

        inline real_type flatten()
        {
            TMVAssert(isComplex(T()));
            TMVAssert(step() == 1);
            return real_type(
                reinterpret_cast<RT*>(ptr()),
                isReal(T()) ? size() : 2*size(), 1, NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)
                ,reinterpret_cast<const RT*>(_last)+ (isReal(T())?0:1)
#endif
            );
        }

        inline const_view_type cSubVector(ptrdiff_t i1, ptrdiff_t i2) const
        { return base::cSubVector(i1,i2); }
        inline const_view_type subVector(ptrdiff_t i1, ptrdiff_t i2) const
        { return base::subVector(i1,i2); }
        inline const_view_type cSubVector(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const
        { return base::cSubVector(i1,i2,istep); }
        inline const_view_type subVector(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const
        { return base::subVector(i1,i2,istep); }
        inline const_view_type reverse() const
        { return base::reverse(); }
        inline const_view_type view() const
        { return base::view(); }
        inline const_cview_type cView() const
        { return base::cView(); }
        inline const_fview_type fView() const
        { return base::fView(); }
        inline const_view_type conjugate() const
        { return base::conjugate(); }
        inline const_real_type realPart() const
        { return base::realPart(); }
        inline const_real_type imagPart() const
        { return base::imagPart(); }
        inline const_real_type flatten() const
        { return base::flatten(); }


        //
        // I/O
        //

        void read(const TMV_Reader& reader);


        inline ptrdiff_t size() const { return _size; }
        inline const T* cptr() const { return _v; }
        inline T* ptr() { return _v; }
        inline ptrdiff_t step() const { return _step; }
        inline ConjType ct() const { return _ct; }

        reference ref(ptrdiff_t i);

    private:

        T*const _v;
        const ptrdiff_t _size;
        const ptrdiff_t _step;
        const ConjType _ct;

#ifdef TMVFLDEBUG
    public:
        const T* _first;
        const T* _last;
#endif

    }; // VectorView

    template <typename T>
    class VectorView<T,FortranStyle> : public VectorView<T,CStyle>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef VectorView<T,FortranStyle> type;
        typedef VectorView<T,CStyle> c_type;
        typedef VectorView<T,FortranStyle> view_type;
        typedef view_type conjugate_type;
        typedef view_type reverse_type;
        typedef VectorView<T,CStyle> cview_type;
        typedef VectorView<T,FortranStyle> fview_type;
        typedef VectorView<RT,FortranStyle> real_type;
        typedef ConstVectorView<T,FortranStyle> const_view_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_reverse_type;
        typedef ConstVectorView<T,CStyle> const_cview_type;
        typedef ConstVectorView<T,FortranStyle> const_fview_type;
        typedef ConstVectorView<RT,FortranStyle> const_real_type;
        typedef ConstVectorView<T,FortranStyle> const_type;
        typedef typename RefHelper<T>::reference reference;
        typedef typename RefHelper<T>::iterator iterator;

        //
        // Constructors
        //

        inline VectorView(const type& rhs) : c_type(rhs) {}

        inline VectorView(const c_type& rhs) : c_type(rhs)  {}

        inline VectorView(
            T* inv, ptrdiff_t insize, ptrdiff_t instep, ConjType inct
            TMV_PARAMFIRSTLAST(T) ) :
            c_type(inv,insize,instep,inct TMV_FIRSTLAST1(_first,_last) ) {}

        virtual inline ~VectorView() {}


        //
        // Op =
        //

        inline type& operator=(const type& v2)
        { c_type::operator=(v2); return *this; }

        inline type& operator=(const c_type& v2)
        { c_type::operator=(v2); return *this; }

        inline type& operator=(const GenVector<RT>& v2)
        { c_type::operator=(v2); return *this; }

        inline type& operator=(const GenVector<CT>& v2)
        { c_type::operator=(v2); return *this; }

        template <typename T2>
        inline type& operator=(const GenVector<T2>& v2)
        { c_type::operator=(v2); return *this; }

        inline type& operator=(const AssignableToVector<RT>& v2)
        { c_type::operator=(v2); return *this; }

        inline type& operator=(const AssignableToVector<CT>& v2)
        { c_type::operator=(v2); return *this; }

        template <typename T2, ptrdiff_t N, int A2>
        inline type& operator=(const SmallVector<T2,N,A2>& v2)
        { c_type::operator=(v2); return *this; }


        //
        // Access Functions
        //

        inline reference operator()(ptrdiff_t i)
        {
            TMVAssert(i>0 && i<=this->size());
            return c_type::ref(i-1);
        }
        inline reference operator[](ptrdiff_t i)
        { return operator()(i); }

        typedef ListAssigner<T,iterator> MyListAssigner;
        inline MyListAssigner operator<<(const T& x)
        { return c_type::operator<<(x); }

        inline T operator()(ptrdiff_t i) const
        {
            TMVAssert(i>0 && i<=this->size());
            return c_type::cref(i-1);
        }
        inline T operator[](ptrdiff_t i) const
        { return operator()(i); }

        //
        // Modifying Functions
        //

        inline type& setZero()
        { c_type::setZero(); return *this; }

        inline type& clip(RT thresh)
        { c_type::clip(thresh); return *this; }

        inline type& setAllTo(const T& x)
        { c_type::setAllTo(x); return *this; }

        inline type& addToAll(const T& x)
        { c_type::addToAll(x); return *this; }

        inline type& conjugateSelf()
        { c_type::conjugateSelf(); return *this; }

        inline type& makeBasis(ptrdiff_t i)
        {
            TMVAssert(i>0 && i<=this->size());
            c_type::makeBasis(i-1);
            return *this;
        }

        inline type& swap(ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(i1>0 && i1<=this->size());
            TMVAssert(i2>0 && i2<=this->size());
            if (i1 != i2) c_type::swap(i1-1,i2-1);
            return *this;
        }

        inline type& permute(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(i1>0 && i1 <= i2 && i2 <= this->size());
            c_type::permute(p,i1-1,i2);
            return *this;
        }

        inline type& permute(const ptrdiff_t* p)
        { c_type::permute(p); return *this; }

        inline type& reversePermute(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(i1>0 && i1 <= i2 && i2 <= this->size());
            c_type::reversePermute(p,i1-1,i2);
            return *this;
        }

        inline type& reversePermute(const ptrdiff_t* p)
        { c_type::reversePermute(p); return *this; }

        inline type& reverseSelf()
        { c_type::reverseSelf(); return *this; }

        inline type& sort(ptrdiff_t* p, ADType ad=Ascend, CompType comp=RealComp)
        { c_type::sort(p,ad,comp); return *this; }

        inline type& sort(
            Permutation& P, ADType ad=Ascend, CompType comp=RealComp)
        { c_type::sort(P,ad,comp); return *this; }

        inline type& sort(ADType ad=Ascend, CompType comp=RealComp)
        { c_type::Sort(0,ad,comp); return *this; }

        //
        // SubVector
        //

        inline bool hasSubVector(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const
        {
            return ConstVectorView<T,FortranStyle>(*this).hasSubVector(
                i1,i2,istep);
        }

        inline view_type subVector(ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(hasSubVector(i1,i2,1));
            return c_type::cSubVector(i1-1,i2);
        }

        inline view_type subVector(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep)
        {
            TMVAssert(hasSubVector(i1,i2,istep));
            return c_type::cSubVector(i1-1,i2-1+istep,istep);
        }

        inline view_type reverse()
        { return c_type::reverse(); }

        inline view_type view()
        { return c_type::view(); }

        inline cview_type cView()
        { return view(); }

        inline fview_type fView()
        { return view(); }

        inline view_type conjugate()
        { return c_type::conjugate(); }

        inline real_type realPart()
        { return c_type::realPart(); }

        inline real_type imagPart()
        { return c_type::imagPart(); }

        inline real_type flatten()
        { return c_type::flatten(); }

        inline const_view_type subVector(ptrdiff_t i1, ptrdiff_t i2) const
        {
            TMVAssert(hasSubVector(i1,i2,1));
            return c_type::cSubVector(i1-1,i2);
        }

        inline const_view_type subVector(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const
        {
            TMVAssert(hasSubVector(i1,i2,istep));
            return c_type::cSubVector(i1-1,i2-1+istep,istep);
        }

        inline const_view_type reverse() const
        { return c_type::reverse(); }

        inline const_view_type view() const
        { return c_type::view(); }

        inline const_cview_type cView() const
        { return view(); }

        inline const_fview_type fView() const
        { return view(); }

        inline const_view_type conjugate() const
        { return c_type::conjugate(); }

        inline const_real_type realPart() const
        { return c_type::realPart(); }

        inline const_real_type imagPart() const
        { return c_type::imagPart(); }

        inline const_real_type flatten() const
        { return c_type::flatten(); }


        //
        // Functions of Vector
        //

        inline T minElement(ptrdiff_t* iminout=0) const
        { return const_type(*this).minElement(iminout); }

        inline T maxElement(ptrdiff_t* imaxout=0) const
        { return const_type(*this).maxElement(imaxout); }

        inline RT minAbsElement(ptrdiff_t* iminout=0) const
        { return const_type(*this).minAbsElement(iminout); }

        inline RT maxAbsElement(ptrdiff_t* imaxout=0) const
        { return const_type(*this).maxAbsElement(imaxout); }

        inline RT minAbs2Element(ptrdiff_t* iminout=0) const
        { return const_type(*this).minAbs2Element(iminout); }

        inline RT maxAbs2Element(ptrdiff_t* imaxout=0) const
        { return const_type(*this).maxAbs2Element(imaxout); }

        // Also allow other int types in case ptrdiff_t is not int:
        template <class INT>
        inline T minElement(INT* iminout) const
        { ptrdiff_t i; T temp=minElement(&i); *iminout=i; return temp; }

        template <class INT>
        inline T maxElement(INT* imaxout) const
        { ptrdiff_t i; T temp=maxElement(&i); *imaxout=i; return temp; }

        template <class INT>
        inline RT minAbsElement(INT* iminout) const
        { ptrdiff_t i; RT temp=minAbsElement(&i); *iminout=i; return temp; }

        template <class INT>
        inline RT maxAbsElement(INT* imaxout) const
        { ptrdiff_t i; RT temp=maxAbsElement(&i); *imaxout=i; return temp; }

        template <class INT>
        inline RT minAbs2Element(INT* iminout) const
        { ptrdiff_t i; RT temp=minAbs2Element(&i); *iminout=i; return temp; }

        template <class INT>
        inline RT maxAbs2Element(INT* imaxout) const
        { ptrdiff_t i; RT temp=maxAbs2Element(&i); *imaxout=i; return temp; }


    }; // FortranStyle VectorView


    template <typename T, int A>
    class Vector : public GenVector<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef Vector<T,A> type;
        typedef ConstVectorView<T,A> const_view_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_reverse_type;
        typedef ConstVectorView<T,CStyle> const_cview_type;
        typedef ConstVectorView<T,FortranStyle> const_fview_type;
        typedef ConstVectorView<RT,A> const_real_type;
        typedef VectorView<T,A> view_type;
        typedef view_type conjugate_type;
        typedef view_type reverse_type;
        typedef VectorView<T,CStyle> cview_type;
        typedef VectorView<T,FortranStyle> fview_type;
        typedef VectorView<RT,A> real_type;
        typedef T value_type;
        typedef VIt<T,1,NonConj> iterator;
        typedef CVIt<T,1,NonConj> const_iterator;
        typedef VIt<T,-1,NonConj> reverse_iterator;
        typedef CVIt<T,-1,NonConj> const_reverse_iterator;
        typedef T& reference;

        //
        // Constructors
        //

#define NEW_SIZE(n) \
        _v(n), _size(n) TMV_DEFFIRSTLAST(_v.get(),_v.get()+n)

        explicit inline Vector(ptrdiff_t n=0) : NEW_SIZE(n)
        {
            TMVAssert(Attrib<A>::vectorok);
            TMVAssert(n >= 0);
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(888));
#endif
        }

        inline Vector(ptrdiff_t n, T val) : NEW_SIZE(n)
        {
            TMVAssert(Attrib<A>::vectorok);
            TMVAssert(n >= 0);
            setAllTo(val);
        }

        inline Vector(const Vector<T,A>& rhs) : NEW_SIZE(rhs.size())
        {
            TMVAssert(Attrib<A>::vectorok);
            std::copy(rhs.cptr(),rhs.cptr()+_size,_v.get());
        }

        template <int A2>
        inline Vector(const Vector<T,A2>& rhs) : NEW_SIZE(rhs.size())
        {
            TMVAssert(Attrib<A>::vectorok);
            std::copy(rhs.cptr(),rhs.cptr()+_size,_v.get());
        }

        inline Vector(const GenVector<RT>& rhs) : NEW_SIZE(rhs.size())
        {
            TMVAssert(Attrib<A>::vectorok);
            if (isReal(T()) && rhs.step() == 1 && !rhs.isconj())
                std::copy(rhs.cptr(),rhs.cptr()+_size,_v.get());
            else rhs.assignToV(view());
        }

        inline Vector(const GenVector<CT>& rhs) : NEW_SIZE(rhs.size())
        {
            TMVAssert(Attrib<A>::vectorok);
            TMVAssert(isComplex(T()));
            if (rhs.step() == 1 && !rhs.isconj())
                std::copy(rhs.cptr(),rhs.cptr()+_size,_v.get());
            else rhs.assignToV(view());
        }

        template <typename T2>
        inline Vector(const GenVector<T2>& rhs) : NEW_SIZE(rhs.size())
        {
            TMVAssert(Attrib<A>::vectorok);
            TMVAssert(isReal(T2()) || isComplex(T()));
            Copy(rhs,view());
        }

        inline Vector(const AssignableToVector<RT>& v2) : NEW_SIZE(v2.size())
        {
            TMVAssert(Attrib<A>::vectorok);
            v2.assignToV(view());
        }

        inline Vector(const AssignableToVector<CT>& v2) : NEW_SIZE(v2.size())
        {
            TMVAssert(Attrib<A>::vectorok);
            TMVAssert(isComplex(T()));
            v2.assignToV(view());
        }

        template <typename T2, ptrdiff_t N, int A2>
        inline Vector(const SmallVector<T2,N,A2>& rhs) :
            NEW_SIZE(rhs.size())
        {
            TMVAssert(Attrib<A>::vectorok);
            TMVAssert(isReal(T2()) || isComplex(T()));
            Copy(rhs.view(),view());
        }


#undef NEW_SIZE

        virtual inline ~Vector()
        {
            TMV_SETFIRSTLAST(0,0);
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(999));
#endif
        }


        //
        // Op =
        //

        inline type& operator=(const type& v2)
        {
            TMVAssert(v2.size() == size());
            if (&v2 != this)
                std::copy(v2.cptr(),v2.cptr()+_size,_v.get());
            return *this;
        }

        template <int A2>
        inline type& operator=(const Vector<T,A2>& v2)
        {
            TMVAssert(v2.size() == size());
            std::copy(v2.cptr(),v2.cptr()+_size,_v.get());
            return *this;
        }

        inline type& operator=(const GenVector<RT>& v2)
        {
            TMVAssert(v2.size() == size());
            v2.assignToV(view());
            return *this;
        }

        inline type& operator=(const GenVector<CT>& v2)
        {
            TMVAssert(v2.size() == size());
            TMVAssert(isComplex(T()));
            v2.assignToV(view());
            return *this;
        }

        template <typename T2>
        inline type& operator=(const GenVector<T2>& v2)
        {
            TMVAssert(v2.size() == size());
            view() = v2;
            return *this;
        }

        inline type& operator=(const AssignableToVector<RT>& v2)
        {
            TMVAssert(v2.size() == size());
            v2.assignToV(view());
            return *this;
        }


        inline type& operator=(const AssignableToVector<CT>& v2)
        {
            TMVAssert(v2.size() == size());
            TMVAssert(isComplex(T()));
            v2.assignToV(view());
            return *this;
        }

        template <typename T2, ptrdiff_t N, int A2>
        inline type& operator=(const SmallVector<T2,N,A2>& v2)
        {
            TMVAssert(v2.size() == size());
            view() = v2.view();
            return *this;
        }


        //
        // Access Functions
        //

        inline const_iterator begin() const
        { return const_iterator(cptr(),1); }
        inline const_iterator end() const
        { return begin() + size(); }
        inline const_reverse_iterator rbegin() const
        { return const_reverse_iterator(cptr()+(size()-1),-1); }
        inline const_reverse_iterator rend() const
        { return rbegin() + size(); }
        inline iterator begin()
        { return iterator(ptr(),1); }
        inline iterator end()
        { return begin() + size(); }
        inline reverse_iterator rbegin()
        { return reverse_iterator(ptr()+(size()-1),-1); }
        inline reverse_iterator rend()
        { return rbegin() + size(); }

        inline T operator()(ptrdiff_t i) const
        {
            if (A == CStyle) {
                TMVAssert(i>=0 && i<size()); return cref(i);
            } else {
                TMVAssert(i>0 && i<=size()); return cref(i-1);
            }
        }
        inline T operator[](ptrdiff_t i) const
        { return operator()(i); }

        inline T& operator()(ptrdiff_t i)
        {
            if (A == CStyle) {
                TMVAssert(i>=0 && i<size()); return ref(i);
            } else {
                TMVAssert(i>0 && i<=size()); return ref(i-1);
            }
        }
        inline T& operator[](ptrdiff_t i)
        { return operator()(i); }

        typedef ListAssigner<T,iterator> MyListAssigner;
        inline MyListAssigner operator<<(const T& x)
        { return MyListAssigner(begin(),size(),x); }


        //
        // Modifying Functions
        //

        type& setZero();

        type& clip(RT thresh);

        type& setAllTo(const T& x);

        type& addToAll(const T& x);

        type& conjugateSelf();

        type& DoBasis(ptrdiff_t i);
        inline type& makeBasis(ptrdiff_t i)
        {
            if (A == CStyle) TMVAssert(i>=0 && i<size());
            else TMVAssert(i>0 && i<=size());
            return DoBasis(i);
        }

        type& DoSwap(ptrdiff_t i1, ptrdiff_t i2);
        inline type& swap(ptrdiff_t i1, ptrdiff_t i2)
        {
            if (A == CStyle)
                TMVAssert(i1>=0 && i1<size() && i2>=0 && i2<size());
            else TMVAssert(i1>0 && i1<=size() && i2>0 && i2<=size());
            return DoSwap(i1,i2);
        }

        inline type& permute(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2)
        {
            if (A==FortranStyle) TMVAssert(i1>0);
            TMVAssert(i1>=0 && i1 <= i2 && i2 <= size());
            view().permute(p,i1,i2); return *this;
        }
        inline type& permute(const ptrdiff_t* p)
        { view().permute(p); return *this; }
        inline type& reversePermute(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2)
        {
            if (A==FortranStyle) TMVAssert(i1>0);
            TMVAssert(i1>=0 && i1 <= i2 && i2 <= size());
            view().reversePermute(p,i1,i2); return *this;
        }
        inline type& reversePermute(const ptrdiff_t* p)
        { view().reversePermute(p,0,size()); return *this; }

        inline type& reverseSelf()
        { view().reverseSelf(); return *this; }

        inline type& sort(ptrdiff_t* p, ADType ad=Ascend, CompType comp=RealComp)
        { view().sort(p,ad,comp); return *this; }

        inline type& sort(
            Permutation& P, ADType ad=Ascend, CompType comp=RealComp)
        { view().sort(P,ad,comp); return *this; }

        inline type& sort(ADType ad=Ascend, CompType comp=RealComp)
        { view().sort(0,ad,comp); return *this; }


        //
        // SubVector
        //

        inline const_view_type cSubVector(ptrdiff_t i1, ptrdiff_t i2) const
        { return const_view_type(cptr()+i1,i2-i1,1,NonConj); }

        inline const_view_type subVector(ptrdiff_t i1, ptrdiff_t i2) const
        {
            TMVAssert(view().hasSubVector(i1,i2,1));
            if (A==FortranStyle) --i1;
            return cSubVector(i1,i2);
        }

        inline view_type cSubVector(ptrdiff_t i1, ptrdiff_t i2)
        { return view_type(ptr()+i1,i2-i1,1,NonConj TMV_FIRSTLAST ); }

        inline view_type subVector(ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(view().hasSubVector(i1,i2,1));
            if (A==FortranStyle) --i1;
            return cSubVector(i1,i2);
        }

        inline const_view_type cSubVector(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const
        { return const_view_type(cptr()+i1,(i2-i1)/istep,istep,NonConj); }

        inline const_view_type subVector(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const
        {
            TMVAssert(view().hasSubVector(i1,i2,istep));
            if (A==FortranStyle) { --i1; i2 += istep-1; }
            return cSubVector(i1,i2,istep);
        }

        inline view_type cSubVector(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep)
        {
            return view_type(
                ptr()+i1,(i2-i1)/istep,istep,NonConj TMV_FIRSTLAST );
        }

        inline view_type subVector(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep)
        {
            TMVAssert(view().hasSubVector(i1,i2,istep));
            if (A==FortranStyle) { --i1; i2 += istep-1; }
            return cSubVector(i1,i2,istep);
        }

        inline const_view_type reverse() const
        { return const_view_type(cptr()+size()-1,size(),-1,NonConj); }

        inline view_type reverse()
        { return view_type(ptr()+size()-1,size(),-1,NonConj TMV_FIRSTLAST ); }

        inline const_view_type view() const
        { return const_view_type(cptr(),size(),1,NonConj); }

        inline view_type view()
        { return view_type(ptr(),size(),1,NonConj TMV_FIRSTLAST ); }

        inline const_cview_type cView() const
        { return view(); }

        inline cview_type cView()
        { return view(); }

        inline const_fview_type fView() const
        { return view(); }

        inline fview_type fView()
        { return view(); }

        inline const_view_type conjugate() const
        { return const_view_type(cptr(),size(),1,TMV_ConjOf(T,NonConj)); }

        inline view_type conjugate()
        { return view_type(ptr(),size(),1,TMV_ConjOf(T,NonConj) TMV_FIRSTLAST ); }

        inline const_real_type realPart() const
        { return view().realPart(); }

        inline const_real_type imagPart() const
        { return view().imagPart(); }

        inline const_real_type flatten() const
        { return view().flatten(); }

        inline real_type realPart()
        { return view().realPart(); }

        inline real_type imagPart()
        { return view().imagPart(); }

        inline real_type flatten()
        { return view().flatten(); }


        //
        // Functions of Vector
        //

        inline T minElement(ptrdiff_t* iminout=0) const
        { return view().minElement(iminout); }

        inline T maxElement(ptrdiff_t* imaxout=0) const
        { return view().maxElement(imaxout); }

        inline RT minAbsElement(ptrdiff_t* iminout=0) const
        { return view().minAbsElement(iminout); }

        inline RT maxAbsElement(ptrdiff_t* imaxout=0) const
        { return view().maxAbsElement(imaxout); }

        inline RT minAbs2Element(ptrdiff_t* iminout=0) const
        { return view().minAbs2Element(iminout); }

        inline RT maxAbs2Element(ptrdiff_t* imaxout=0) const
        { return view().maxAbs2Element(imaxout); }

        // Also allow other int types in case ptrdiff_t is not int:
        template <class INT>
        inline T minElement(INT* iminout) const
        { ptrdiff_t i; T temp=minElement(&i); *iminout=i; return temp; }

        template <class INT>
        inline T maxElement(INT* imaxout) const
        { ptrdiff_t i; T temp=maxElement(&i); *imaxout=i; return temp; }

        template <class INT>
        inline RT minAbsElement(INT* iminout) const
        { ptrdiff_t i; RT temp=minAbsElement(&i); *iminout=i; return temp; }

        template <class INT>
        inline RT maxAbsElement(INT* imaxout) const
        { ptrdiff_t i; RT temp=maxAbsElement(&i); *imaxout=i; return temp; }

        template <class INT>
        inline RT minAbs2Element(INT* iminout) const
        { ptrdiff_t i; RT temp=minAbs2Element(&i); *iminout=i; return temp; }

        template <class INT>
        inline RT maxAbs2Element(INT* imaxout) const
        { ptrdiff_t i; RT temp=maxAbs2Element(&i); *imaxout=i; return temp; }


        //
        // I/O
        //

        void read(const TMV_Reader& reader);

        inline ptrdiff_t size() const { return _size; }
        inline const T* cptr() const { return _v.get(); }
        inline T* ptr() { return _v.get(); }
        inline ptrdiff_t step() const { return 1; }
        inline ConjType ct() const { return NonConj; }
        inline bool isconj() const { return false; }

        inline T cref(ptrdiff_t i) const
        { return _v.get()[i]; }

        inline T& ref(ptrdiff_t i)
        { return _v.get()[i]; }

        inline void resize(ptrdiff_t n)
        {
            TMVAssert(n >= 0);
            _v.resize(n);
            _size = n;
#ifdef TMVFLDEBUG
            _first = _v.get();
            _last = _first + n;
#endif
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(888));
#endif
        }

    private:

        AlignedArray<T> _v;
        ptrdiff_t _size;

#ifdef TMVFLDEBUG
    public:
        const T* _first;
        const T* _last;
    protected:
#endif

        template <int A2>
        friend void Swap(Vector<T,A>& v1, Vector<T,A2>& v2)
        {
            TMVAssert(v1.size() == v2.size());
            T* temp = v1._v.release();
            v1._v.reset(v2._v.release());
            v2._v.reset(temp);
#ifdef TMVFLDEBUG
            TMV_SWAP(v1._first,v2._first);
            TMV_SWAP(v1._last,v2._last);
#endif
        }

    }; // Vector


    //
    // Special Constructors
    //

    template <typename T, int A>
    Vector<T,A> DoBasisVector(ptrdiff_t n, ptrdiff_t i);
    template <typename T, int A>
    inline Vector<T,A> BasisVector(ptrdiff_t n, ptrdiff_t i)
    {
        TMVAssert(Attrib<A>::vectorok);
        TMVAssert(n > 0);
        if (A == CStyle) { TMVAssert(i>=0 && i<n); }
        else { TMVAssert(i>0 && i<=n); }
        return DoBasisVector<T,A>(n,i);
    }
    template <typename T>
    inline Vector<T,CStyle> BasisVector(ptrdiff_t n, ptrdiff_t i)
    {
        TMVAssert(n > 0);
        return DoBasisVector<T,CStyle>(n,i);
    }

    template <typename T>
    inline VectorView<T,CStyle> VectorViewOf(T* v, ptrdiff_t size)
    {
        return VectorView<T,CStyle>(
            v,size,1,NonConj TMV_FIRSTLAST1(v,v+size));
    }

    template <typename T>
    inline ConstVectorView<T,CStyle> VectorViewOf(const T* v, ptrdiff_t size)
    { return ConstVectorView<T,CStyle>(v,size,1,NonConj); }

    template <typename T>
    inline VectorView<T,CStyle> VectorViewOf(T* v, ptrdiff_t size, ptrdiff_t step)
    {
        return VectorView<T,CStyle>(
            v,size,step,NonConj TMV_FIRSTLAST1(v,v+size));
    }

    template <typename T>
    inline ConstVectorView<T,CStyle> VectorViewOf(
        const T* v, ptrdiff_t size, ptrdiff_t step)
    { return ConstVectorView<T,CStyle>(v,size,step,NonConj); }


    //
    // Copy Vectors
    //

    inline bool shouldReverse(const ptrdiff_t step1, const ptrdiff_t step2)
    {
        return ( (step2 < 0 && (step1 != 1 || step2 == -1)) ||
                 (step1 == -1 && step2 != 1) );
    }

    template <typename T>
    void DoCopySameType(const GenVector<T>& v1, VectorView<T> v2);

    template <typename T>
    inline void DoCopy(const GenVector<T>& v1, VectorView<T> v2)
    { if (!v1.isSameAs(v2)) DoCopySameType(v1,v2); }

    template <typename T, typename T1>
    inline void DoCopyDiffType(const GenVector<T1>& v1, VectorView<T> v2)
    {
        TMVAssert(isReal(T1()) || isComplex(T()));
        TMVAssert(v1.size()==v2.size());
        TMVAssert(v2.size()>0);
        TMVAssert(v1.ct()==NonConj);
        TMVAssert(v2.ct()==NonConj);
        TMVAssert(v2.step() != -1);
        TMVAssert(v1.step() != -1 || v2.step() == 1);
        TMVAssert(v2.step() > 0 || v1.step() == 1);
        TMVAssert(!v2.isSameAs(v1));

        const T1* v1ptr = v1.cptr();
        T* v2ptr = v2.ptr();
        const ptrdiff_t step1 = v1.step();
        const ptrdiff_t step2 = v2.step();

        if (step1 == 1 && step2 == 1) {
            for(ptrdiff_t i=v2.size();i>0;--i,++v1ptr,++v2ptr) {
#ifdef TMVFLDEBUG
                TMVAssert(v2ptr >= v2._first);
                TMVAssert(v2ptr < v2._last);
#endif
                *v2ptr = *v1ptr;
            }
        } else {
            for(ptrdiff_t i=v2.size();i>0;--i,v1ptr+=step1,v2ptr+=step2) {
#ifdef TMVFLDEBUG
                TMVAssert(v2ptr >= v2._first);
                TMVAssert(v2ptr < v2._last);
#endif
                *v2ptr = *v1ptr;
            }
        }
    }

    template <typename T, typename T1>
    inline void DoCopy(const GenVector<T1>& v1, VectorView<T> v2)
    {
        TMVAssert(isReal(T1()) || isComplex(T()));
        DoCopyDiffType(v1,v2);
    }

    template <typename T>
    inline void DoCopy(const GenVector<std::complex<T> >&, VectorView<T>)
    { TMVAssert(TMV_FALSE); }

    template <typename T, typename T1>
    inline void Copy(const GenVector<T1>& v1, VectorView<T> v2)
    {
        TMVAssert(isReal(T1()) || isComplex(T()));
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v2.step()!=0 || v1.step() == 0);

        if (v1.size() > 0) {
            if (shouldReverse(v1.step(),v2.step())) {
                Copy(v1.reverse(),v2.reverse());
            } else {
                if (v1.isconj()) {
                    if (v2.isconj()) {
                        DoCopy(v1.conjugate(),v2.conjugate());
                    } else {
                        DoCopy(v1.conjugate(),v2);
                        v2.conjugateSelf();
                    }
                } else {
                    if (v2.isconj()) {
                        DoCopy(v1,v2.conjugate());
                        v2.conjugateSelf();
                    }
                    else DoCopy(v1,v2);
                }
            }
        }
    }


    //
    // Swap Vectors
    //

    template <typename T>
    void Swap(VectorView<T> v1, VectorView<T> v2);
    template <typename T, int A>
    inline void Swap(VectorView<T> v1, Vector<T,A>& v2)
    { Swap(v1.view(),v2.view()); }
    template <typename T, int A>
    inline void Swap(Vector<T,A>& v1, VectorView<T> v2)
    { Swap(v1.view(),v2.view()); }


    //
    // Functions of Vectors
    //

    template <typename T>
    inline TMV_RealType(T) Norm(const GenVector<T>& v)
    { return v.norm(); }

    template <typename T>
    inline TMV_RealType(T) Norm1(const GenVector<T>& v)
    { return v.norm1(); }

    template <typename T>
    inline TMV_RealType(T) NormSq(const GenVector<T>& v)
    { return v.normSq(); }

    template <typename T>
    inline TMV_RealType(T) Norm2(const GenVector<T>& v)
    { return v.norm2(); }

    template <typename T>
    inline TMV_RealType(T) NormInf(const GenVector<T>& v)
    { return v.normInf(); }

    template <typename T>
    inline T SumElements(const GenVector<T>& v)
    { return v.sumElements(); }

    template <typename T>
    inline TMV_RealType(T) SumAbsElements(const GenVector<T>& v)
    { return v.sumAbsElements(); }

    template <typename T>
    inline TMV_RealType(T) SumAbs2Elements(const GenVector<T>& v)
    { return v.sumAbs2Elements(); }

    template <typename T>
    inline T MinElement(const GenVector<T>& v)
    { return v.minElement(); }

    template <typename T, int A>
    inline T MinElement(const ConstVectorView<T,A>& v)
    { return v.minElement(); }

    template <typename T, int A>
    inline T MinElement(const VectorView<T,A>& v)
    { return v.minElement(); }

    template <typename T, int A>
    inline T MinElement(const Vector<T,A>& v)
    { return v.minElement(); }

    template <typename T>
    inline T MaxElement(const GenVector<T>& v)
    { return v.maxElement(); }

    template <typename T, int A>
    inline T MaxElement(const ConstVectorView<T,A>& v)
    { return v.maxElement(); }

    template <typename T, int A>
    inline T MaxElement(const VectorView<T,A>& v)
    { return v.maxElement(); }

    template <typename T, int A>
    inline T MaxElement(const Vector<T,A>& v)
    { return v.maxElement(); }

    template <typename T>
    inline TMV_RealType(T) MinAbsElement(const GenVector<T>& v)
    { return v.minAbsElement(); }

    template <typename T, int A>
    inline TMV_RealType(T) MinAbsElement(const ConstVectorView<T,A>& v)
    { return v.minAbsElement(); }

    template <typename T, int A>
    inline TMV_RealType(T) MinAbsElement(const VectorView<T,A>& v)
    { return v.minAbsElement(); }

    template <typename T, int A>
    inline TMV_RealType(T) MinAbsElement(const Vector<T,A>& v)
    { return v.minAbsElement(); }

    template <typename T>
    inline TMV_RealType(T) MinAbs2Element(const GenVector<T>& v)
    { return v.minAbs2Element(); }

    template <typename T, int A>
    inline TMV_RealType(T) MinAbs2Element(const ConstVectorView<T,A>& v)
    { return v.minAbs2Element(); }

    template <typename T, int A>
    inline TMV_RealType(T) MinAbs2Element(const VectorView<T,A>& v)
    { return v.minAbs2Element(); }

    template <typename T, int A>
    inline TMV_RealType(T) MinAbs2Element(const Vector<T,A>& v)
    { return v.minAbs2Element(); }

    template <typename T>
    inline TMV_RealType(T) MaxAbsElement(const GenVector<T>& v)
    { return v.maxAbsElement(); }

    template <typename T, int A>
    inline TMV_RealType(T) MaxAbsElement(const ConstVectorView<T,A>& v)
    { return v.maxAbsElement(); }

    template <typename T, int A>
    inline TMV_RealType(T) MaxAbsElement(const VectorView<T,A>& v)
    { return v.maxAbsElement(); }

    template <typename T, int A>
    inline TMV_RealType(T) MaxAbsElement(const Vector<T,A>& v)
    { return v.maxAbsElement(); }

    template <typename T>
    inline TMV_RealType(T) MaxAbs2Element(const GenVector<T>& v)
    { return v.maxAbs2Element(); }

    template <typename T, int A>
    inline TMV_RealType(T) MaxAbs2Element(const ConstVectorView<T,A>& v)
    { return v.maxAbs2Element(); }

    template <typename T, int A>
    inline TMV_RealType(T) MaxAbs2Element(const VectorView<T,A>& v)
    { return v.maxAbs2Element(); }

    template <typename T, int A>
    inline TMV_RealType(T) MaxAbs2Element(const Vector<T,A>& v)
    { return v.maxAbs2Element(); }

    template <typename T>
    inline ConstVectorView<T> Conjugate(const GenVector<T>& v)
    { return v.conjugate(); }

    template <typename T, int A>
    inline ConstVectorView<T,A> Conjugate(const ConstVectorView<T,A>& v)
    { return v.conjugate(); }

    template <typename T, int A>
    inline ConstVectorView<T,A> Conjugate(const Vector<T,A>& v)
    { return v.conjugate(); }

    template <typename T, int A>
    inline VectorView<T,A> Conjugate(VectorView<T,A> v)
    { return v.conjugate(); }

    template <typename T, int A>
    inline VectorView<T,A> Conjugate(Vector<T,A>& v)
    { return v.conjugate(); }


    //
    // Vector ==, != Vector
    //

    template <typename T1, typename T2>
    bool operator==(const GenVector<T1>& v1, const GenVector<T2>& v2);

    template <typename T1, typename T2>
    inline bool operator!=(const GenVector<T1>& v1, const GenVector<T2>& v2)
    { return !(v1 == v2); }

    //
    // I/O
    //

    template <typename T>
    inline std::ostream& operator<<(
        const TMV_Writer& writer, const GenVector<T>& v)
    { v.write(writer); return writer.getos(); }

    template <typename T>
    inline std::istream& operator>>(
        const TMV_Reader& reader, VectorView<T> v)
    { v.read(reader); return reader.getis(); }

    template <typename T, int A>
    inline std::istream& operator>>(
        const TMV_Reader& reader, Vector<T,A>& v)
    { v.read(reader); return reader.getis(); }


    template <typename T>
    inline std::ostream& operator<<(std::ostream& os, const GenVector<T>& v)
    { return os << IOStyle() << v; }

    template <typename T>
    inline std::istream& operator>>(std::istream& is, VectorView<T> v)
    { return is >> IOStyle() >> v; }

    template <typename T, int A>
    inline std::istream& operator>>(std::istream& is, Vector<T,A>& v)
    { return is >> IOStyle() >> v; }

} // namespace tmv

#endif
