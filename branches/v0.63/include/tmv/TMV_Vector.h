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
// This file defines the TMV Vector class.
//
// The Vector class and all associated functions are contained
// in the namespace tmv.  Alse, the Vector class is a template, 
// so for a Vector of doubles, one would write 
// tmv::Vector<double>.  
//
// A second optional template argument can tell the Vector to
// use Fortran-style indexing instead of C-style.
// ie. the first element is v(1), and the last is v(n).
// So you would write Vector<T,FortranStyle> v(n), for example.
// If you want to use the Fortran-style indexing, you should
// be aware that it only exists for a regular Vector and the two
// view classes, VectorView and ConstVectorView.  
// GenVector and VectorComposite always use C-style indexing,
// so (v1+v2)(1) would return the second element in the sum, even
// if v1 and v2 both use Fortran-style indexing.
// Also, permutation arrays always use the C-style indexing.
//
// Constructors:
//
//    explicit Vector<T>(size_t n)
//        Makes a Vector of size n with _uninitialized_ values
//
//    Vector<T>(size_t n, T x)
//        Makes a Vector of size n with all values = x
//
//    Vector<T>(size_t n, const T* vv)
//    Vector<T>(const vector<T>& vv)
//        Makes a vector which copies the elements of vv
//        For the second one, n specifies the length of the vector
//
// Special Constructors
//
//    makeBasisVector(size_t n, int i)
//        Makes a Vector whose elements are all 0, except v(i) = 1
//
//    VectorViewOf(T* vv, size_t n)
//    VectorViewOf(const T* vv, size_t n)
//        Makes a Vector View (see below) which refers to the exact
//        elements of vv, not copying them to new storage.
//        The first one returns a VectorView, the second a ConstVectorView.
//
// Access Functions
//
//    size_t size() const
//        Returns the size of the Vector
//
//    T& operator[](int i)
//    T& operator()(int i)
//    T operator[](int i) const
//    T operator()(int i) const
//        Return the ith element of the Vector
//
//    Vector<T>::iterator begin()
//    Vector<T>::iterator end()
//    Vector<T>::reverse_iterator rbegin()
//    Vector<T>::reverse_iterator rend()
//    Vector<T>::const_iterator begin() const
//    Vector<T>::const_iterator end() const
//    Vector<T>::const_reverse_iterator rbegin() const
//    Vector<T>::const_reverse_iterator rend() const
//        Return iterators that work in the usual way to traverse the Vector
//
// Modifying Functions
//
//    Vector& zero()
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
//    sort(int* p, AD, COMP)
//        Sorts the vector, returning the swaps required in P.
//        If you do not care about P, you may omit the P parameter.
//        AD = ASCEND or DESCEND (ASCEND=default)
//        COMP = REAL_COMP, ABS_COMP, IMAG_COMP, or ARG_COMP (REAL_COMP=default)
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
//    VectorView real()
//    VectorView imag()
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
//        Returns the sum of absolute values of  elements in the Vector
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
//        Note: v must already be the correct size
//
//    is >> vptr
//        If you do not know the size of the Vector to be read in, you can
//        use this form where vptr is an auto_ptr to an undefined Vector.
//
//


#ifndef TMV_Vector_H
#define TMV_Vector_H

#include "tmv/TMV_BaseVector.h"
#include <vector>
#include "tmv/TMV_ListInit.h"

#ifdef TMVFLDEBUG
#include "tmv/TMV_VIt.h"
#endif

namespace tmv {

    template <class T, class T1> 
    inline void Copy(const GenVector<T1>& v1, const VectorView<T>& v2);

    template <class T> 
    class GenVector : public AssignableToVector<T>
    {
        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef GenVector<T> type;
        typedef ConstVectorView<T> const_view_type;
        typedef ConstVectorView<RT> const_real_type;
        typedef VectorView<T> nonconst_type;

    public:

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
        void assignToV(const VectorView<RT>& rhs) const
        {
            TMVAssert(rhs.size() == size());
            TMVAssert(isReal(T()));
            if (!isSameAs(rhs)) Copy(*this,rhs); 
        }
        void assignToV(const VectorView<CT>& rhs) const
        { 
            TMVAssert(rhs.size() == size());
            if (!isSameAs(rhs)) Copy(*this,rhs); 
        }

        inline T operator[](int i) const 
        { 
            TMVAssert(i>=0 && i<int(size()));
            return cref(i); 
        }
        inline T operator()(int i) const 
        { 
            TMVAssert(i>=0 && i<int(size()));
            return cref(i); 
        }

        typedef T value_type;
        typedef CVIter<T> const_iterator;
        typedef CVIter<T> const_reverse_iterator;

        inline const_iterator begin() const
        { return const_iterator(cptr(),step(),ct()); }
        inline const_iterator end() const
        { return begin() + size(); }
        inline const_reverse_iterator rbegin() const
        { return const_reverse_iterator(cptr(),step(),ct()); }
        inline const_reverse_iterator rend() const
        { return rbegin() + size(); }

        template <class T2> 
        inline bool isSameAs(const GenVector<T2>&) const
        { return false; }

        inline bool isSameAs(const GenVector<T>& v2) const
        {
            return (
                this == &v2 || 
                ( cptr()==v2.cptr() && size()==v2.size() && 
                  step()==v2.step() && ct()==v2.ct() ) );
        }

        template <class T2> 
        TMV_DEPRECATED(bool SameAs(const GenVector<T2>& v2) const);
        TMV_DEPRECATED(bool SameAs(const GenVector<T>& v2) const)
        { return isSameAs(v2); }

        //
        // subVector
        //

        bool hasSubVector(int i1, int i2, int istep) const;

        inline const_view_type subVector(int i1, int i2) const
        {
            TMVAssert(hasSubVector(i1,i2,1));
            return const_view_type(cptr()+i1*step(),i2-i1,step(),ct());
        }

        inline const_view_type subVector(int i1, int i2, int istep) const
        {
            TMVAssert(hasSubVector(i1,i2,istep));
            return const_view_type(
                cptr()+i1*step(),(i2-i1)/istep,istep*step(),ct());
        }

        inline const_view_type reverse() const
        { 
            return const_view_type(
                cptr()+(int(size())-1)*step(),size(),-step(),ct()); 
        }

        inline const_view_type view() const
        { return const_view_type(cptr(),size(),step(),ct()); }

        inline const_view_type conjugate() const
        { return const_view_type(cptr(),size(),step(),TMV_ConjOf(T,ct())); }

        inline const_real_type real() const
        { 
            return const_real_type(
                reinterpret_cast<const RT*>( cptr()),
                size(), isReal(T()) ? step() : 2*step(), NonConj);
        }

        inline const_real_type imag() const
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

        TMV_DEPRECATED(const_view_type SubVector(int i1, int i2) const)
        { return subVector(i1,i2); }
        TMV_DEPRECATED(const_view_type SubVector(
                int i1, int i2, int istep) const)
        { return subVector(i1,i2,istep); }
        TMV_DEPRECATED(const_view_type Reverse() const)
        { return reverse(); }
        TMV_DEPRECATED(const_view_type View() const)
        { return view(); }
        TMV_DEPRECATED(const_view_type Conjugate() const)
        { return conjugate(); }
        TMV_DEPRECATED(const_real_type Real() const)
        { return real(); }
        TMV_DEPRECATED(const_real_type Imag() const)
        { return imag(); }
        TMV_DEPRECATED(const_real_type Flatten() const)
        { return flatten(); }
        TMV_DEPRECATED(nonconst_type NonConst() const)
        { return nonConst(); }


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

        T minElement(int* iminout=0) const;

        T maxElement(int* imaxout=0) const;

        RT minAbsElement(int* iminout=0) const;

        RT maxAbsElement(int* imaxout=0) const;

        TMV_DEPRECATED(RT Norm() const)
        { return norm(); }
        TMV_DEPRECATED(RT NormSq(const RT scale = RT(1)) const)
        { return normSq(scale); }
        TMV_DEPRECATED(RT Norm1() const )
        { return norm1(); }
        TMV_DEPRECATED(RT Norm2() const)
        { return norm2(); }
        TMV_DEPRECATED(RT NormInf() const )
        { return normInf(); }
        TMV_DEPRECATED(T SumElements() const)
        { return sumElements(); }
        TMV_DEPRECATED(RT SumAbsElements() const)
        { return sumAbsElements(); }
        TMV_DEPRECATED(T MinElement(int* iminout=0) const)
        { return minElement(iminout); }
        TMV_DEPRECATED(T MaxElement(int* imaxout=0) const)
        { return maxElement(imaxout); }
        TMV_DEPRECATED(RT MinAbsElement(int* iminout=0) const)
        { return minAbsElement(iminout); }
        TMV_DEPRECATED(RT MaxAbsElement(int* imaxout=0) const)
        { return maxAbsElement(imaxout); }

        //
        // I/O
        //

        void write(std::ostream& os) const;
        void write(std::ostream& os, RT minNonZero) const;

        TMV_DEPRECATED(void Write(std::ostream& os) const)
        { write(os); }
        TMV_DEPRECATED(void Write(std::ostream& os, RT minNonZero) const)
        { write(os,minNonZero); }

        virtual const T* cptr() const =0;
        virtual int step() const =0;
        virtual ConjType ct() const =0;
        inline bool isconj() const 
        { return (isComplex(T()) && ct()==Conj); }

        virtual T cref(int i) const;

    private:

        type& operator=(const GenVector<T>&);

    }; // GenVector

    template <class T> template <class T2> 
    inline bool GenVector<T>::SameAs(const GenVector<T2>& v2) const
    { return isSameAs(v2); }

    template <class T, IndexStyle I> 
    class ConstVectorView : public GenVector<T>
    {
        typedef TMV_RealType(T) RT;
        typedef ConstVectorView<T,I> type;
        typedef ConstVectorView<RT,I> const_real_type;

    public:

        inline ConstVectorView(const ConstVectorView<T,I>& rhs) : 
            _v(rhs._v), _size(rhs._size), _step(rhs._step),
            _ct(rhs._ct) {}
        inline ConstVectorView(const GenVector<T>& rhs) : 
            _v(rhs.cptr()), _size(rhs.size()), _step(rhs.step()),
            _ct(rhs.ct()) {}
        inline ConstVectorView(
            const T* v, size_t size, int step, ConjType ct) : 
            _v(v), _size(size), _step(step), _ct(ct) {}
        virtual inline ~ConstVectorView() 
        {
#ifdef TMVDEBUG
            const_cast<const T*&>(_v) = 0;
#endif
        }

        virtual inline size_t size() const { return _size; }
        virtual inline const T* cptr() const { return _v; }
        virtual inline int step() const { return _step; }
        virtual inline ConjType ct() const { return _ct; }

    private:

        const T*const _v;
        const size_t _size;
        const int _step;
        const ConjType _ct;

        type& operator=(const ConstVectorView<T,I>&);

    }; // ConstVectorView

    template <class T> 
    class ConstVectorView<T,FortranStyle> : public ConstVectorView<T,CStyle>
    {
        typedef TMV_RealType(T) RT;
        typedef ConstVectorView<T,FortranStyle> type;
        typedef ConstVectorView<T,FortranStyle> const_view_type;
        typedef ConstVectorView<RT,FortranStyle> const_real_type;

    public:

        inline ConstVectorView(const ConstVectorView<T,FortranStyle>& rhs) : 
            ConstVectorView<T,CStyle>(rhs) {}
        inline ConstVectorView(const ConstVectorView<T,CStyle>& rhs) : 
            ConstVectorView<T,CStyle>(rhs) {}
        inline ConstVectorView(const GenVector<T>& rhs) : 
            ConstVectorView<T,CStyle>(rhs) {}
        inline ConstVectorView(const T* inv, size_t insize, int instep, 
                               ConjType inct) : 
            ConstVectorView<T,CStyle>(inv,insize,instep,inct) {}
        virtual inline ~ConstVectorView() {}

        //
        // Access Functions
        //

        inline T operator[](int i) const 
        { 
            TMVAssert(i>0 && i<=int(this->size()));
            return GenVector<T>::cref(i-1); 
        }
        inline T operator()(int i) const 
        { 
            TMVAssert(i>0 && i<=int(this->size()));
            return GenVector<T>::cref(i-1); 
        }

        // 
        // SubVector
        //

        bool hasSubVector(int i1, int i2, int istep) const;

        inline const_view_type subVector(int i1, int i2) const
        {
            TMVAssert(hasSubVector(i1,i2,1));
            return GenVector<T>::subVector(i1-1,i2);
        }

        inline const_view_type subVector(
            int i1, int i2, int istep) const
        {
            TMVAssert(hasSubVector(i1,i2,istep));
            return GenVector<T>::subVector(i1-1,i2-1+istep,istep);
        }

        inline const_view_type reverse() const
        { return GenVector<T>::reverse(); }

        inline const_view_type view() const
        { return GenVector<T>::view(); }

        inline const_view_type conjugate() const
        { return GenVector<T>::conjugate(); }

        inline const_real_type real() const
        { return GenVector<T>::real(); }

        inline const_real_type imag() const
        { return GenVector<T>::imag(); }

        inline const_real_type flatten() const
        { return GenVector<T>::flatten(); }

        TMV_DEPRECATED(const_view_type SubVector(int i1, int i2) const)
        { return subVector(i1,i2); }
        TMV_DEPRECATED(const_view_type SubVector(int i1, int i2, int istep) const)
        { return subVector(i1,i2,istep); }
        TMV_DEPRECATED(const_view_type Reverse() const)
        { return reverse(); }
        TMV_DEPRECATED(const_view_type View() const)
        { return view(); }
        TMV_DEPRECATED(const_view_type Conjugate() const)
        { return conjugate(); }
        TMV_DEPRECATED(const_real_type Real() const)
        { return real(); }
        TMV_DEPRECATED(const_real_type Imag() const)
        { return imag(); }
        TMV_DEPRECATED(const_real_type Flatten() const)
        { return flatten(); }

        //
        // Functions of Vector
        //

        inline T minElement(int* iminout=0) const
        { 
            T temp = GenVector<T>::minElement(iminout);
            if (iminout) ++(*iminout);
            return temp;
        }

        inline T maxElement(int* imaxout=0) const
        { 
            T temp = GenVector<T>::maxElement(imaxout);
            if (imaxout) ++(*imaxout);
            return temp;
        }

        inline RT minAbsElement(int* iminout=0) const
        { 
            RT temp = GenVector<T>::minAbsElement(iminout);
            if (iminout) ++(*iminout);
            return temp;
        }

        inline RT maxAbsElement(int* imaxout=0) const
        { 
            RT temp = GenVector<T>::maxAbsElement(imaxout);
            if (imaxout) ++(*imaxout);
            return temp;
        }

        TMV_DEPRECATED(T MinElement(int* iminout=0) const)
        { return minElement(iminout); }
        TMV_DEPRECATED(T MaxElement(int* imaxout=0) const)
        { return maxElement(imaxout); }
        TMV_DEPRECATED(RT MinAbsElement(int* iminout=0) const)
        { return minAbsElement(iminout); }
        TMV_DEPRECATED(RT MaxAbsElement(int* imaxout=0) const)
        { return maxAbsElement(imaxout); }


    private:

        const type& operator=(const type&);

    }; // FortranStyle ConstVectorView

    template <class T, IndexStyle I> 
    class VectorView : public GenVector<T>
    {
        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef VectorView<T,I> type;
        typedef VectorView<T,I> view_type;
        typedef VectorView<RT,I> real_type;

    public:

        //
        // Constructors 
        //

        inline VectorView(const type& rhs) : 
            _v(rhs._v), _size(rhs._size), _step(rhs._step),
            _ct(rhs._ct) TMV_DEFFIRSTLAST(rhs._first,rhs._last) {}

        inline VectorView(
            T* v, size_t size, int step, ConjType ct 
            TMV_PARAMFIRSTLAST(T) ) :
            _v(v), _size(size), _step(step),
            _ct(ct) TMV_DEFFIRSTLAST(_first,_last) {}

        virtual inline ~VectorView() 
        {
#ifdef TMVDEBUG
            const_cast<T*&>(_v) = 0;
#endif
        }


        //
        // Op =
        //

        inline const type& operator=(const type& v2) const
        { 
            TMVAssert(size() == v2.size());
            v2.assignToV(*this);
            return *this; 
        }

        inline const type& operator=(const type& v2) 
        { 
            TMVAssert(size() == v2.size());
            v2.assignToV(*this);
            return *this; 
        }

        inline const type& operator=(const GenVector<RT>& v2) const
        { 
            TMVAssert(size() == v2.size());
            v2.assignToV(*this);
            return *this; 
        }

        inline const type& operator=(const GenVector<CT>& v2) const
        { 
            TMVAssert(size() == v2.size());
            TMVAssert(isComplex(T()));
            v2.assignToV(*this);
            return *this; 
        }

        template <class T2> 
        inline const type& operator=(const GenVector<T2>& v2) const
        { 
            TMVAssert(isReal(T2()) || isComplex(T()));
            TMVAssert(size() == v2.size());
            Copy(v2,*this);
            return *this; 
        }

        inline const type& operator=(const AssignableToVector<RT>& v2) const
        { 
            TMVAssert(size() == v2.size());
            v2.assignToV(*this);
            return *this; 
        }

        inline const type& operator=(const AssignableToVector<CT>& v2) const
        { 
            TMVAssert(size() == v2.size());
            TMVAssert(isComplex(T()));
            v2.assignToV(*this);
            return *this; 
        }

        template <class T2, int N, IndexStyle I2> 
        inline const type& operator=(const SmallVector<T2,N,I2>& v2) const
        { 
            TMVAssert(size() == v2.size());
            Copy(v2.view(),*this);
            return *this; 
        }

        template <int N> 
        inline const type& operator=(
            const SmallVectorComposite<RT,N>& v2) const
        {
            TMVAssert(size() == N);
            v2.assignToV(*this);
            return *this; 
        }

        template <int N> 
        inline const type& operator=(
            const SmallVectorComposite<CT,N>& v2) const
        {
            TMVAssert(size() == N);
            TMVAssert(isComplex(T()));
            v2.assignToV(*this);
            return *this; 
        }

        //
        // Access Functions
        //

        typedef T value_type;
        typedef VIter<T> iterator;
        typedef CVIter<T> const_iterator;
        typedef VIter<T> reverse_iterator;
        typedef CVIter<T> const_reverse_iterator;
        typedef TMV_RefType(T) reference;

        inline reference operator[](int i) const 
        { 
            TMVAssert(i>=0 && i<int(size()));
            return ref(i); 
        }
        inline reference operator()(int i) const 
        { 
            TMVAssert(i>=0 && i<int(size()));
            return ref(i); 
        }

        inline iterator begin() const
        { return iterator(ptr(),step(),ct() TMV_FIRSTLAST ); }
        inline iterator end() const
        { return begin() + size(); }
        inline reverse_iterator rbegin() const
        {
            return reverse_iterator(
                ptr()+step()*(int(size())-1),-step(),ct() TMV_FIRSTLAST ); 
        }
        inline reverse_iterator rend() const
        { return rbegin() + size(); }

        typedef ListAssigner<T,iterator> MyListAssigner;
        MyListAssigner inline operator<<(const T& x)
        { return MyListAssigner(begin(),size(),x); }

        TMV_DEPRECATED(MyListAssigner operator=(ListInitClass))
        { return MyListAssigner(begin(),size()); }


        //
        // Modifying Functions
        //

        const type& zero() const;

        const type& clip(RT thresh) const;

        const type& setAllTo(const T& x) const;

        const type& addToAll(const T& x) const;

        const type& conjugateSelf() const;

        const type& DoBasis(int i) const;
        inline const type& makeBasis(int i) const
        { TMVAssert(i>=0 && i<int(size())); return DoBasis(i); }

        const type& DoSwap(int i1, int i2) const;
        inline const type& swap(int i1, int i2) const
        { 
            TMVAssert(i1>=0 && i1<int(size()));
            TMVAssert(i2>=0 && i2<int(size()));
            return DoSwap(i1,i2);
        }

        const type& DoPermute(const int* p, int i1, int i2) const;
        inline const type& permute(const int* p, int i1, int i2) const
        {
            TMVAssert(i1>=0 && i1 <= i2 && i2 <= int(size()));
            return DoPermute(p,i1,i2);
        }
        inline const type& permute(const int* p) const
        { return DoPermute(p,0,size()); }

        const type& DoReversePermute(const int* p, int i1, int i2) const;
        inline const type& reversePermute(const int* p, int i1, int i2) const
        {
            TMVAssert(i1>=0 && i1 <= i2 && i2 <= int(size()));
            return DoReversePermute(p,i1,i2);
        }
        inline const type& reversePermute(const int* p) const
        { return reversePermute(p,0,size()); }

        const type& reverseSelf() const;

        const type& sort(
            int* p, ADType ad=Ascend, CompType comp=RealComp) const;

        inline const type& sort(ADType ad=Ascend, CompType comp=RealComp) const
        { sort(0,ad,comp); return *this; }


        TMV_DEPRECATED(const type& Zero() const)
        { return zero(); }
        TMV_DEPRECATED(const type& Clip(RT thresh) const)
        { return clip(thresh); }
        TMV_DEPRECATED(const type& SetAllTo(const T& x) const)
        { return setAllTo(x); }
        TMV_DEPRECATED(const type& AddToAll(const T& x) const)
        { return addToAll(x); }
        TMV_DEPRECATED(const type& ConjugateSelf() const)
        { return conjugateSelf(); }
        TMV_DEPRECATED(const type& MakeBasis(int i) const)
        { return makeBasis(i); }
        TMV_DEPRECATED(const type& Swap(int i1, int i2) const)
        { return swap(i1,i2); }
        TMV_DEPRECATED(const type& Permute(
                const int* p, int i1, int i2) const)
        { return permute(p,i1,i2); }
        TMV_DEPRECATED(const type& Permute(const int* p) const)
        { return permute(p); }
        TMV_DEPRECATED(const type& ReversePermute(
                const int* p, int i1, int i2) const)
        { return reversePermute(p,i1,i2); }
        TMV_DEPRECATED(const type& ReversePermute(const int* p) const)
        { return reversePermute(p); }
        TMV_DEPRECATED(const type& ReverseSelf() const)
        { return reverseSelf(); }
        TMV_DEPRECATED(const type& Sort(
                int* p, OldADType ad=ASCEND, OldCompType comp=REAL_COMP) const)
        { return sort(p,ADType(ad),CompType(comp)); }
        TMV_DEPRECATED(const type& Sort(
                OldADType ad=ASCEND, OldCompType comp=REAL_COMP) const)
        { return sort(ADType(ad),CompType(comp)); }

        //
        // SubVector
        //

        inline view_type subVector(int i1, int i2) const
        {
            TMVAssert(GenVector<T>::hasSubVector(i1,i2,1));
            return view_type(ptr()+i1*step(),(i2-i1),step(),
                        ct() TMV_FIRSTLAST );
        }

        inline view_type subVector(int i1, int i2, int istep) const
        {
            TMVAssert(GenVector<T>::hasSubVector(i1,i2,istep));
            return view_type(ptr()+i1*step(),(i2-i1)/istep,istep*step(),
                        ct() TMV_FIRSTLAST );
        }

        inline view_type reverse() const
        { 
            return view_type(ptr()+(int(size())-1)*step(),size(),-step(),
                        ct() TMV_FIRSTLAST ); 
        }

        inline view_type view() const
        { return *this; }

        inline view_type conjugate() const
        {
            return view_type(ptr(),size(),step(),
                        TMV_ConjOf(T,ct()) TMV_FIRSTLAST ); 
        }

        inline real_type real() const
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

        inline real_type imag() const
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

        inline real_type flatten() const
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

        TMV_DEPRECATED(view_type SubVector(int i1, int i2) const)
        { return subVector(i1,i2); }
        TMV_DEPRECATED(view_type SubVector(int i1, int i2, int istep) const)
        { return subVector(i1,i2,istep); }
        TMV_DEPRECATED(view_type Reverse() const)
        { return reverse(); }
        TMV_DEPRECATED(view_type View() const)
        { return view(); }
        TMV_DEPRECATED(view_type Conjugate() const)
        { return conjugate(); }
        TMV_DEPRECATED(real_type Real() const)
        { return real(); }
        TMV_DEPRECATED(real_type Imag() const)
        { return imag(); }
        TMV_DEPRECATED(real_type Flatten() const)
        { return flatten(); }


        //
        // I/O
        //

        void read(std::istream& is) const;

        TMV_DEPRECATED(void Read(std::istream& is) const)
        { read(is); }

        // 
        // Iterator Typedefs
        //

        virtual inline size_t size() const { return _size; }
        virtual inline const T* cptr() const { return _v; }
        virtual inline T* ptr() const { return _v; }
        virtual inline int step() const { return _step; }
        virtual inline ConjType ct() const { return _ct; }

        reference ref(int i) const;

    private:

        T*const _v;
        const size_t _size;
        const int _step;
        const ConjType _ct;

#ifdef TMVFLDEBUG
    public:
        const T*const _first;
        const T*const _last;
#endif

    }; // VectorView

    template <class T> 
    class VectorView<T,FortranStyle> : public VectorView<T,CStyle>
    {
        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef VectorView<T,FortranStyle> type;
        typedef VectorView<T,CStyle> c_type;
        typedef VectorView<T,FortranStyle> view_type;
        typedef VectorView<RT,FortranStyle> real_type;
        typedef ConstVectorView<T,FortranStyle> const_type;

    public:

        //
        // Constructors 
        //

        inline VectorView(const type& rhs) : c_type(rhs) {}

        inline VectorView(const c_type& rhs) : c_type(rhs) {}

        inline VectorView(
            T* inv, size_t insize, int instep, ConjType inct
            TMV_PARAMFIRSTLAST(T) ) :
            c_type(inv,insize,instep,inct TMV_FIRSTLAST1(_first,_last) ) {}

        virtual inline ~VectorView() {}


        //
        // Op =
        //

        inline const type& operator=(const type& v2) const
        { c_type::operator=(v2); return *this; }

        inline const type& operator=(const type& v2)
        { c_type::operator=(v2); return *this; }

        inline const type& operator=(const c_type& v2) const
        { c_type::operator=(v2); return *this; }

        inline const type& operator=(const c_type& v2)
        { c_type::operator=(v2); return *this; }

        inline const type& operator=(const GenVector<RT>& v2) 
        { c_type::operator=(v2); return *this; }

        inline const type& operator=(const GenVector<CT>& v2) const
        { c_type::operator=(v2); return *this; }

        template <class T2> 
        inline const type& operator=(const GenVector<T2>& v2) const
        { c_type::operator=(v2); return *this; }

        inline const type& operator=(
            const AssignableToVector<RT>& v2) const
        { c_type::operator=(v2); return *this; }

        inline const type& operator=(
            const AssignableToVector<CT>& v2) const
        { c_type::operator=(v2); return *this; }

        template <class T2, int N, IndexStyle I2> 
        inline const type& operator=(const SmallVector<T2,N,I2>& v2) const
        { c_type::operator=(v2); return *this; }

        template <int N> 
        inline const type& operator=(const SmallVectorComposite<RT,N>& v2) const
        { c_type::operator=(v2); return *this; }

        template <int N> 
        inline const type& operator=(const SmallVectorComposite<CT,N>& v2) const
        { c_type::operator=(v2); return *this; }


        //
        // Access Functions
        //

        inline TMV_RefType(T) operator[](int i) const 
        { 
            TMVAssert(i>0 && i<=int(this->size()));
            return c_type::ref(i-1); 
        }
        inline TMV_RefType(T) operator()(int i) const 
        { 
            TMVAssert(i>0 && i<=int(this->size()));
            return c_type::ref(i-1); 
        }

        typedef ListAssigner<T,VIter<T> > MyListAssigner;
        inline MyListAssigner operator<<(const T& x)
        { return c_type::operator<<(x); }

        TMV_DEPRECATED(MyListAssigner operator=(ListInitClass li))
        { return c_type::operator=(li); }


        //
        // Modifying Functions
        //

        inline const type& zero() const 
        { c_type::aero(); return *this; }

        inline const type& clip(RT thresh) const
        { c_type::clip(thresh); return *this; }

        inline const type& setAllTo(const T& x) const
        { c_type::setAllTo(x); return *this; }

        inline const type& addToAll(const T& x) const
        { c_type::addToAll(x); return *this; }

        inline const type& conjugateSelf() const
        { c_type::conjugateSelf(); return *this; }

        inline const type& makeBasis(int i) const
        {
            TMVAssert(i>0 && i<=int(this->size()));
            c_type::makeBasis(i-1); 
            return *this; 
        }

        inline const type& swap(int i1, int i2) const
        {
            TMVAssert(i1>0 && i1<=int(this->size()));
            TMVAssert(i2>0 && i2<=int(this->size()));
            if (i1 != i2) c_type::swap(i1-1,i2-1);
            return *this;
        }

        inline const type& permute(const int* p, 
                                   int i1, int i2) const
        { 
            TMVAssert(i1>0 && i1 <= i2 && i2 <= int(this->size()));
            c_type::permute(p,i1-1,i2); 
            return *this; 
        }

        inline const type& permute(const int* p) const
        { c_type::permute(p); return *this; }

        inline const type& reversePermute(
            const int* p, int i1, int i2) const
        { 
            TMVAssert(i1>0 && i1 <= i2 && i2 <= int(this->size()));
            c_type::reversePermute(p,i1-1,i2); 
            return *this; 
        }

        inline const type& reversePermute(
            const int* p) const
        { c_type::reversePermute(p); return *this; }

        inline const type& reverseSelf() const
        { c_type::reverseSelf(); return *this; }

        inline const type& sort(
            int* p, ADType ad=Ascend, CompType comp=RealComp) const
        { c_type::sort(p,ad,comp); return *this; }

        inline const type& sort(ADType ad=Ascend, CompType comp=RealComp) const
        { c_type::Sort(0,ad,comp); return *this; }

        TMV_DEPRECATED(const type& Zero() const)
        { return zero(); }
        TMV_DEPRECATED(const type& Clip(RT thresh) const)
        { return clip(thresh); }
        TMV_DEPRECATED(const type& SetAllTo(const T& x) const)
        { return setAllTo(x); }
        TMV_DEPRECATED(const type& AddToAll(const T& x) const)
        { return addToAll(x); }
        TMV_DEPRECATED(const type& ConjugateSelf() const)
        { return conjugateSelf(); }
        TMV_DEPRECATED(const type& MakeBasis(int i) const)
        { return makeBasis(i); }
        TMV_DEPRECATED(const type& Swap(int i1, int i2) const)
        { return swap(i1,i2); }
        TMV_DEPRECATED(const type& Permute(
                const int* p, int i1, int i2) const)
        { return permute(p,i1,i2); }
        TMV_DEPRECATED(const type& Permute(const int* p) const)
        { return permute(p); }
        TMV_DEPRECATED(const type& ReversePermute(
                const int* p, int i1, int i2) const)
        { return reversePermute(p,i1,i2); }
        TMV_DEPRECATED(const type& ReversePermute(const int* p) const)
        { return reversePermute(p); }
        TMV_DEPRECATED(const type& ReverseSelf() const)
        { return reverseSelf(); }
        TMV_DEPRECATED(const type& Sort(
                int* p, OldADType ad=ASCEND, OldCompType comp=REAL_COMP) const)
        { return sort(p,ADType(ad),CompType(comp)); }
        TMV_DEPRECATED(const type& Sort(
                OldADType ad=ASCEND, OldCompType comp=REAL_COMP) const)
        { return sort(ADType(ad),CompType(comp)); }


        //
        // SubVector
        //

        inline bool hasSubVector(int i1, int i2, int istep) const
        { 
            return ConstVectorView<T,FortranStyle>(*this).hasSubVector(
                i1,i2,istep); 
        }

        inline view_type subVector(int i1, int i2) const
        {
            TMVAssert(hasSubVector(i1,i2,1));
            return c_type::subVector(i1-1,i2);
        }

        inline view_type subVector(int i1, int i2, int istep) const
        {
            TMVAssert(hasSubVector(i1,i2,istep));
            return c_type::subVector(i1-1,i2-1+istep,istep);
        }

        inline view_type reverse() const
        { return c_type::reverse(); }

        inline view_type view() const
        { return c_type::view(); }

        inline view_type conjugate() const
        { return c_type::conjugate(); }

        inline real_type real() const
        { return c_type::real(); }

        inline real_type imag() const
        { return c_type::imag(); }

        inline real_type flatten() const
        { return c_type::flatten(); }


        TMV_DEPRECATED(view_type SubVector(int i1, int i2) const)
        { return subVector(i1,i2); }
        TMV_DEPRECATED(view_type SubVector(int i1, int i2, int istep) const)
        { return subVector(i1,i2,istep); }
        TMV_DEPRECATED(view_type Reverse() const)
        { return reverse(); }
        TMV_DEPRECATED(view_type View() const)
        { return view(); }
        TMV_DEPRECATED(view_type Conjugate() const)
        { return conjugate(); }
        TMV_DEPRECATED(real_type Real() const)
        { return real(); }
        TMV_DEPRECATED(real_type Imag() const)
        { return imag(); }
        TMV_DEPRECATED(real_type Flatten() const)
        { return flatten(); }


        //
        // Functions of Vector
        //

        inline T minElement(int* iminout=0) const
        { return const_type(*this).minElement(iminout); }

        inline T maxElement(int* imaxout=0) const
        { return const_type(*this).maxElement(imaxout); }

        inline RT minAbsElement(int* iminout=0) const
        { return const_type(*this).minAbsElement(iminout); }

        inline RT maxAbsElement(int* imaxout=0) const
        { return const_type(*this).maxAbsElement(imaxout); }

        TMV_DEPRECATED(T MinElement(int* iminout=0) const)
        { return minElement(iminout); }
        TMV_DEPRECATED(T MaxElement(int* imaxout=0) const)
        { return maxElement(imaxout); }
        TMV_DEPRECATED(RT MinAbsElement(int* iminout=0) const)
        { return minAbsElement(iminout); }
        TMV_DEPRECATED(RT MaxAbsElement(int* imaxout=0) const)
        { return maxAbsElement(imaxout); }

    }; // FortranStyle VectorView

#ifdef XTEST
#ifdef TMVDEBUG
#define XTEST_DEBUG
#endif
#endif

    template <class T, IndexStyle I> 
    class Vector : public GenVector<T>
    {
        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef Vector<T,I> type;
        typedef ConstVectorView<T,I> const_view_type;
        typedef ConstVectorView<RT,I> const_real_type;
        typedef VectorView<T,I> view_type;
        typedef VectorView<RT,I> real_type;

    public:

        //
        // Constructors
        //

#define NEW_SIZE(n) \
        _v(new T[(n)]), _size(n) TMV_DEFFIRSTLAST(_v.get(),_v.get()+n)

        explicit inline Vector(size_t n) : NEW_SIZE(n)
        {
#ifdef TMVDEBUG
            setAllTo(T(888));
#endif
        }

        inline Vector(size_t n, T val) : NEW_SIZE(n)
        {
            setAllTo(val); 
        }

        inline Vector(size_t n, const T* vv) : NEW_SIZE(n)
        { 
#ifdef XTEST_DEBUG
            setAllTo(T(888));
#endif
            std::copy(vv,vv+n,_v.get());
        }

        inline explicit Vector(const std::vector<T>& vv) : NEW_SIZE(vv.size())
        {
#ifdef XTEST_DEBUG
            setAllTo(T(888));
#endif
            std::copy(vv.begin(),vv.end(),_v.get());
        }

        inline Vector(const Vector<T,I>& rhs) : NEW_SIZE(rhs.size())
        {
#ifdef XTEST_DEBUG
            setAllTo(T(888));
#endif
            std::copy(rhs.cptr(),rhs.cptr()+_size,_v.get());
        }

        template <IndexStyle I2> 
        inline Vector(const Vector<T,I2>& rhs) : NEW_SIZE(rhs.size())
        {
#ifdef TMVDEBUG
            setAllTo(T(888));
#endif
            std::copy(rhs.cptr(),rhs.cptr()+_size,_v.get());
        }

        inline Vector(const GenVector<RT>& rhs) : NEW_SIZE(rhs.size())
        {
#ifdef XTEST_DEBUG
            setAllTo(T(888));
#endif
            if (isReal(T()) && rhs.step() == 1 && !rhs.isconj())
                std::copy(rhs.cptr(),rhs.cptr()+_size,_v.get());
            else rhs.assignToV(view());
        }

        inline Vector(const GenVector<CT>& rhs) : NEW_SIZE(rhs.size())
        {
#ifdef XTEST_DEBUG
            setAllTo(T(888));
#endif
            TMVAssert(isComplex(T()));
            if (rhs.step() == 1 && !rhs.isconj())
                std::copy(rhs.cptr(),rhs.cptr()+_size,_v.get());
            else rhs.assignToV(view());
        }

        template <class T2> 
        inline Vector(const GenVector<T2>& rhs) : NEW_SIZE(rhs.size())
        {
#ifdef XTEST_DEBUG
            setAllTo(T(888));
#endif
            TMVAssert(isReal(T2()) || isComplex(T()));
            Copy(rhs,view()); 
        }

        inline Vector(const AssignableToVector<RT>& v2) : NEW_SIZE(v2.size())
        { 
#ifdef XTEST_DEBUG
            setAllTo(T(888));
#endif
            v2.assignToV(view()); 
        }

        inline Vector(const AssignableToVector<CT>& v2) : NEW_SIZE(v2.size())
        { 
#ifdef XTEST_DEBUG
            setAllTo(T(888));
#endif
            TMVAssert(isComplex(T()));
            v2.assignToV(view()); 
        }

        template <class T2, int N, IndexStyle I2> 
        inline Vector(const SmallVector<T2,N,I2>& rhs) : 
            NEW_SIZE(rhs.size())
        {
#ifdef XTEST_DEBUG
            setAllTo(T(888));
#endif
            TMVAssert(isReal(T2()) || isComplex(T()));
            Copy(rhs.view(),view()); 
        }

        template <int N> 
        inline Vector(const SmallVectorComposite<RT,N>& v2) : 
            NEW_SIZE(v2.size())
        { 
#ifdef XTEST_DEBUG
            setAllTo(T(888));
#endif
            v2.assignToV(view()); 
        }

        template <int N> 
        inline Vector(const SmallVectorComposite<CT,N>& v2) :
            NEW_SIZE(v2.size())
        { 
#ifdef XTEST_DEBUG
            setAllTo(T(888));
#endif
            TMVAssert(isComplex(T()));
            v2.assignToV(view()); 
        }

#undef NEW_SIZE

        virtual inline ~Vector() 
        {
#ifdef TMVDEBUG
            setAllTo(T(999));
#endif
        }


        //
        // Op =
        //

        inline type& operator=(type& v2)
        { 
            TMVAssert(v2.size() == size());
            if (&v2 != this) 
                std::copy(v2.cptr(),v2.cptr()+_size,_v.get());
            return *this; 
        }

        template <IndexStyle I2> 
        inline type& operator=(Vector<T,I2>& v2)
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

        template <class T2> 
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

        template <class T2, int N, IndexStyle I2> 
        inline type& operator=(const SmallVector<T2,N,I2>& v2) 
        { 
            TMVAssert(v2.size() == size());
            view() = v2.view(); 
            return *this; 
        }

        template <int N> 
        inline type& operator=(const SmallVectorComposite<RT,N>& v2)
        { 
            TMVAssert(N == size());
            v2.assignToV(view()); 
            return *this; 
        }

        template <int N> 
        inline type& operator=(const SmallVectorComposite<CT,N>& v2)
        { 
            TMVAssert(N == size());
            TMVAssert(isComplex(T()));
            v2.assignToV(view()); 
            return *this; 
        }


        //
        // Access Functions
        //

        typedef T value_type;
        typedef VIt<T,Unit,NonConj> iterator;
        typedef CVIt<T,Unit,NonConj> const_iterator;
        typedef VIt<T,Step,NonConj> reverse_iterator;
        typedef CVIt<T,Step,NonConj> const_reverse_iterator;
        typedef T& reference;

        inline const_iterator begin() const
        { return const_iterator(cptr(),1); }
        inline const_iterator end() const
        { return begin() + size(); }
        inline const_reverse_iterator rbegin() const
        { return const_reverse_iterator(cptr()+(int(size())-1),-1); }
        inline const_reverse_iterator rend() const
        { return rbegin() + size(); }
        inline iterator begin()
        { return iterator(ptr(),1); }
        inline iterator end()
        { return begin() + size(); }
        inline reverse_iterator rbegin()
        { return reverse_iterator(ptr()+(int(size())-1),-1); }
        inline reverse_iterator rend()
        { return rbegin() + size(); }

        inline T operator[](int i) const 
        { 
            if (I == CStyle) { 
                TMVAssert(i>=0 && i<int(size())); return cref(i); 
            } else { 
                TMVAssert(i>0 && i<=int(size())); return cref(i-1); 
            }
        }
        inline T operator()(int i) const 
        { 
            if (I == CStyle) {
                TMVAssert(i>=0 && i<int(size())); return cref(i); 
            } else { 
                TMVAssert(i>0 && i<=int(size())); return cref(i-1); 
            }
        }

        inline T& operator[](int i) 
        { 
            if (I == CStyle) {
                TMVAssert(i>=0 && i<int(size())); return ref(i); 
            } else { 
                TMVAssert(i>0 && i<=int(size())); return ref(i-1); 
            }
        }
        inline T& operator()(int i) 
        { 
            if (I == CStyle) {
                TMVAssert(i>=0 && i<int(size())); return ref(i); 
            } else { 
                TMVAssert(i>0 && i<=int(size())); return ref(i-1); 
            }
        }

        typedef ListAssigner<T,iterator> MyListAssigner;
        inline MyListAssigner operator<<(const T& x)
        { return MyListAssigner(begin(),size(),x); }

        TMV_DEPRECATED(MyListAssigner operator=(ListInitClass))
        { return MyListAssigner(begin(),size()); }


        //
        // Modifying Functions
        //

        type& zero();

        type& clip(RT thresh);

        type& setAllTo(const T& x);

        type& addToAll(const T& x);

        type& conjugateSelf();

        type& DoBasis(int i);
        inline type& makeBasis(int i)
        { 
            if (I == CStyle) TMVAssert(i>=0 && i<int(size()));
            else TMVAssert(i>0 && i<=int(size()));
            return DoBasis(i);
        }

        type& DoSwap(int i1, int i2);
        inline type& swap(int i1, int i2)
        {
            if (I == CStyle) 
                TMVAssert(i1>=0 && i1<int(size()) && i2>=0 && i2<int(size()));
            else TMVAssert(i1>0 && i1<=int(size()) && i2>0 && i2<=int(size()));
            return DoSwap(i1,i2);
        }

        inline type& permute(const int* p, int i1, int i2)
        {
            if (I==FortranStyle) TMVAssert(i1>0);
            TMVAssert(i1>=0 && i1 <= i2 && i2 <= int(size()));
            view().permute(p,i1,i2); return *this; 
        }
        inline type& permute(const int* p) 
        { view().permute(p); return *this; }
        inline type& reversePermute(const int* p, int i1, int i2)
        {
            if (I==FortranStyle) TMVAssert(i1>0);
            TMVAssert(i1>=0 && i1 <= i2 && i2 <= int(size()));
            view().reversePermute(p,i1,i2); return *this; 
        }
        inline type& reversePermute(const int* p) 
        { view().reversePermute(p,0,size()); return *this; }

        inline type& reverseSelf()
        { view().reverseSelf(); return *this; }

        inline type& sort(int* p, ADType ad=Ascend, CompType comp=RealComp) 
        { view().sort(p,ad,comp); return *this; }

        inline type& sort(ADType ad=Ascend, CompType comp=RealComp) 
        { view().sort(0,ad,comp); return *this; }


        TMV_DEPRECATED(type& Zero())
        { return zero(); }
        TMV_DEPRECATED(type& Clip(RT thresh))
        { return clip(thresh); }
        TMV_DEPRECATED(type& SetAllTo(const T& x))
        { return setAllTo(x); }
        TMV_DEPRECATED(type& AddToAll(const T& x))
        { return addToAll(x); }
        TMV_DEPRECATED(type& ConjugateSelf())
        { return conjugateSelf(); }
        TMV_DEPRECATED(type& MakeBasis(int i))
        { return makeBasis(i); }
        TMV_DEPRECATED(type& Swap(int i1, int i2))
        { return swap(i1,i2); }
        TMV_DEPRECATED(type& Permute(
                const int* p, int i1, int i2))
        { return permute(p,i1,i2); }
        TMV_DEPRECATED(type& Permute(const int* p))
        { return permute(p); }
        TMV_DEPRECATED(type& ReversePermute(
                const int* p, int i1, int i2))
        { return reversePermute(p,i1,i2); }
        TMV_DEPRECATED(type& ReversePermute(const int* p))
        { return reversePermute(p); }
        TMV_DEPRECATED(type& ReverseSelf())
        { return reverseSelf(); }
        TMV_DEPRECATED(type& Sort(
                int* p, OldADType ad=ASCEND, OldCompType comp=REAL_COMP))
        { return sort(p,ADType(ad),CompType(comp)); }
        TMV_DEPRECATED(type& Sort(
                OldADType ad=ASCEND, OldCompType comp=REAL_COMP))
        { return sort(ADType(ad),CompType(comp)); }


        //
        // SubVector
        //

        inline const_view_type subVector(int i1, int i2) const
        {
            TMVAssert(view().hasSubVector(i1,i2,1));
            if (I==FortranStyle) --i1;
            return const_view_type(cptr()+i1,i2-i1,1,NonConj);
        }

        inline view_type subVector(int i1, int i2)
        {
            TMVAssert(view().hasSubVector(i1,i2,1));
            if (I==FortranStyle) --i1;
            return view_type(ptr()+i1,i2-i1,1,NonConj TMV_FIRSTLAST );
        }

        inline const_view_type subVector(int i1, int i2, int istep) const
        {
            TMVAssert(view().hasSubVector(i1,i2,istep));
            if (I==FortranStyle) { --i1; i2 += istep-1; }
            return const_view_type(cptr()+i1,(i2-i1)/istep,istep,NonConj);
        }

        inline view_type subVector(int i1, int i2, int istep)
        {
            TMVAssert(view().hasSubVector(i1,i2,istep));
            if (I==FortranStyle) { --i1; i2 += istep-1; }
            return view_type(
                ptr()+i1,(i2-i1)/istep,istep,NonConj TMV_FIRSTLAST );
        }

        inline const_view_type reverse() const
        { return const_view_type(cptr()+size()-1,size(),-1,NonConj); }

        inline view_type reverse()
        { return view_type(ptr()+size()-1,size(),-1,NonConj TMV_FIRSTLAST ); }

        inline const_view_type view() const
        { return const_view_type(cptr(),size(),1,NonConj); }

        inline view_type view()
        { return view_type(ptr(),size(),1,NonConj TMV_FIRSTLAST ); }

        inline const_view_type conjugate() const
        { return const_view_type(cptr(),size(),1,TMV_ConjOf(T,NonConj)); }

        inline view_type conjugate()
        { return view_type(ptr(),size(),1,TMV_ConjOf(T,NonConj) TMV_FIRSTLAST ); }

        inline const_real_type real() const
        { return view().real(); }

        inline const_real_type imag() const
        { return view().imag(); }

        inline const_real_type flatten() const
        { return view().flatten(); }

        inline real_type real()
        { return view().real(); }

        inline real_type imag()
        { return view().imag(); }

        inline real_type flatten()
        { return view().flatten(); }


        TMV_DEPRECATED(const_view_type SubVector(int i1, int i2) const)
        { return subVector(i1,i2); }
        TMV_DEPRECATED(const_view_type SubVector(
                int i1, int i2, int istep) const)
        { return subVector(i1,i2,istep); }
        TMV_DEPRECATED(const_view_type Reverse() const)
        { return reverse(); }
        TMV_DEPRECATED(const_view_type View() const)
        { return view(); }
        TMV_DEPRECATED(const_view_type Conjugate() const)
        { return conjugate(); }
        TMV_DEPRECATED(const_real_type Real() const)
        { return real(); }
        TMV_DEPRECATED(const_real_type Imag() const)
        { return imag(); }
        TMV_DEPRECATED(const_real_type Flatten() const)
        { return flatten(); }

        TMV_DEPRECATED(view_type SubVector(int i1, int i2))
        { return subVector(i1,i2); }
        TMV_DEPRECATED(view_type SubVector(int i1, int i2, int istep))
        { return subVector(i1,i2,istep); }
        TMV_DEPRECATED(view_type Reverse())
        { return reverse(); }
        TMV_DEPRECATED(view_type View())
        { return view(); }
        TMV_DEPRECATED(view_type Conjugate())
        { return conjugate(); }
        TMV_DEPRECATED(real_type Real())
        { return real(); }
        TMV_DEPRECATED(real_type Imag())
        { return imag(); }
        TMV_DEPRECATED(real_type Flatten())
        { return flatten(); }


        //
        // Functions of Vector
        //

        inline T minElement(int* iminout=0) const
        { return view().minElement(iminout); }

        inline T maxElement(int* imaxout=0) const
        { return view().maxElement(imaxout); }

        inline RT minAbsElement(int* iminout=0) const
        { return view().minAbsElement(iminout); }

        inline RT maxAbsElement(int* imaxout=0) const
        { return view().maxAbsElement(imaxout); }

        TMV_DEPRECATED(T MinElement(int* iminout=0) const)
        { return minElement(iminout); }
        TMV_DEPRECATED(T MaxElement(int* imaxout=0) const)
        { return maxElement(imaxout); }
        TMV_DEPRECATED(RT MinAbsElement(int* iminout=0) const)
        { return minAbsElement(iminout); }
        TMV_DEPRECATED(RT MaxAbsElement(int* imaxout=0) const)
        { return maxAbsElement(imaxout); }


        // 
        // I/O
        //

        inline void read(std::istream& is)
        { view().read(is); }

        TMV_DEPRECATED(void Read(std::istream& is))
        { read(is); }


        virtual inline size_t size() const { return _size; }
        virtual inline const T* cptr() const { return _v.get(); }
        virtual inline T* ptr() { return _v.get(); }
        virtual inline int step() const { return 1; }
        virtual inline ConjType ct() const { return NonConj; }
        inline bool isconj() const { return false; }

        inline T cref(int i) const
        { return _v.get()[i]; }

        inline T& ref(int i)
        { return _v.get()[i]; }

    private:

        auto_array<T> _v;
        const size_t _size;

#ifdef TMVFLDEBUG
    public:
        const T*const _first;
        const T*const _last;
#endif

    }; // Vector


    //
    // Special Constructors
    //

    template <class T, IndexStyle I> 
    Vector<T,I> DoBasisVector(size_t n, int i);
    template <class T, IndexStyle I> 
    inline Vector<T,I> BasisVector(size_t n, int i)
    { 
        if (I == CStyle) { TMVAssert(i>=0 && i<n); } 
        else { TMVAssert(i>0 && i<=n); }
        return DoBasisVector<T,I>(n,i); 
    }
    template <class T> 
    inline Vector<T,CStyle> BasisVector(size_t n, int i)
    { return DoBasisVector<T,CStyle>(n,i); }

    template <IndexStyle I, class T> 
    inline VectorView<T,I> VectorViewOf(T* v, size_t size)
    { return VectorView<T,I>(v,size,1,NonConj TMV_FIRSTLAST1(v,v+size)); }

    template <IndexStyle I, class T> 
    inline ConstVectorView<T,I> VectorViewOf(const T* v, size_t size)
    { return ConstVectorView<T,I>(v,size,1,NonConj); }

    template <class T> 
    inline VectorView<T,CStyle> VectorViewOf(T* v, size_t size)
    {
        return VectorView<T,CStyle>(
            v,size,1,NonConj TMV_FIRSTLAST1(v,v+size)); 
    }

    template <class T> 
    inline ConstVectorView<T,CStyle> VectorViewOf(const T* v, size_t size)
    { return ConstVectorView<T,CStyle>(v,size,1,NonConj); }


    //
    // Copy Vectors
    //

    inline bool shouldReverse(const int step1, const int step2)
    {
        return ( (step2 < 0 && (step1 != 1 || step2 == -1)) ||
                 (step1 == -1 && step2 != 1) );
    }

    template <class T> 
    void DoCopySameType(const GenVector<T>& v1, const VectorView<T>& v2);

    template <class T> 
    inline void DoCopy(const GenVector<T>& v1, const VectorView<T>& v2)
    { if (!v1.isSameAs(v2)) DoCopySameType(v1,v2); }

    template <class T, class T1> 
    inline void DoCopyDiffType(const GenVector<T1>& v1, const VectorView<T>& v2)
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
        const int step1 = v1.step();
        const int step2 = v2.step();

        if (step1 == 1 && step2 == 1)
            for(int i=v2.size();i>0;--i,++v1ptr,++v2ptr) {
#ifdef TMVFLDEBUG
                TMVAssert(v2ptr >= v2._first);
                TMVAssert(v2ptr < v2._last);
#endif
                *v2ptr = *v1ptr;
            }
        else
            for(int i=v2.size();i>0;--i,v1ptr+=step1,v2ptr+=step2) {
#ifdef TMVFLDEBUG
                TMVAssert(v2ptr >= v2._first);
                TMVAssert(v2ptr < v2._last);
#endif
                *v2ptr = *v1ptr;
            }
    }

    template <class T, class T1> 
    inline void DoCopy(const GenVector<T1>& v1, const VectorView<T>& v2)
    {
        TMVAssert(isReal(T1()) || isComplex(T()));
        DoCopyDiffType(v1,v2);
    }

    template <class T> 
    inline void DoCopy(const GenVector<std::complex<T> >&, const VectorView<T>&)
    { TMVAssert(TMV_FALSE); }

    template <class T, class T1> 
    inline void Copy(const GenVector<T1>& v1, const VectorView<T>& v2)
    { 
        TMVAssert(isReal(T1()) || isComplex(T()));
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v2.step()!=0 || v1.step() == 0);

        if (v1.size() > 0) {
            if (shouldReverse(v1.step(),v2.step()))
                Copy(v1.reverse(),v2.reverse());
            else {
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

    template <class T> 
    void Swap(const VectorView<T>& v1, const VectorView<T>& v2);
    template <class T, IndexStyle I> 
    inline void Swap(const VectorView<T>& v1, Vector<T,I>& v2)
    { Swap(v1.view(),v2.view()); }
    template <class T, IndexStyle I> 
    inline void Swap(Vector<T,I>& v1, const VectorView<T>& v2) 
    { Swap(v1.view(),v2.view()); }
    template <class T, IndexStyle I1, IndexStyle I2> 
    inline void Swap(Vector<T,I1>& v1, Vector<T,I2>& v2)
    { Swap(v1.view(),v2.view()); }


    //
    // Functions of Vectors
    //

    template <class T> 
    inline TMV_RealType(T) Norm(const GenVector<T>& v)
    { return v.norm(); }

    template <class T> 
    inline TMV_RealType(T) Norm1(const GenVector<T>& v)
    { return v.norm1(); }

    template <class T> 
    inline TMV_RealType(T) NormSq(const GenVector<T>& v)
    { return v.normSq(); }

    template <class T> 
    inline TMV_RealType(T) Norm2(const GenVector<T>& v)
    { return v.norm2(); }

    template <class T> 
    inline TMV_RealType(T) NormInf(const GenVector<T>& v)
    { return v.normInf(); }

    template <class T> 
    inline T SumElements(const GenVector<T>& v)
    { return v.sumElements(); }

    template <class T> 
    inline TMV_RealType(T) SumAbsElements(const GenVector<T>& v)
    { return v.sumAbsElements(); }

    template <class T> 
    inline T MinElement(const GenVector<T>& v)
    { return v.minElement(); }

    template <class T, IndexStyle I> 
    inline T MinElement(const ConstVectorView<T,I>& v)
    { return v.minElement(); }

    template <class T, IndexStyle I> 
    inline T MinElement(const VectorView<T,I>& v)
    { return v.minElement(); }

    template <class T, IndexStyle I> 
    inline T MinElement(const Vector<T,I>& v)
    { return v.minElement(); }

    template <class T> 
    inline T MaxElement(const GenVector<T>& v)
    { return v.maxElement(); }

    template <class T, IndexStyle I> 
    inline T MaxElement(const ConstVectorView<T,I>& v)
    { return v.maxElement(); }

    template <class T, IndexStyle I> 
    inline T MaxElement(const VectorView<T,I>& v)
    { return v.maxElement(); }

    template <class T, IndexStyle I> 
    inline T MaxElement(const Vector<T,I>& v)
    { return v.maxElement(); }

    template <class T> 
    inline TMV_RealType(T) MinAbsElement(const GenVector<T>& v)
    { return v.minAbsElement(); }

    template <class T, IndexStyle I> 
    inline TMV_RealType(T) MinAbsElement(const ConstVectorView<T,I>& v)
    { return v.minAbsElement(); }

    template <class T, IndexStyle I> 
    inline TMV_RealType(T) MinAbsElement(const VectorView<T,I>& v)
    { return v.minAbsElement(); }

    template <class T, IndexStyle I> 
    inline TMV_RealType(T) MinAbsElement(const Vector<T,I>& v)
    { return v.minAbsElement(); }

    template <class T> 
    inline TMV_RealType(T) MaxAbsElement(const GenVector<T>& v)
    { return v.maxAbsElement(); }

    template <class T, IndexStyle I> 
    inline TMV_RealType(T) MaxAbsElement(const ConstVectorView<T,I>& v)
    { return v.maxAbsElement(); }

    template <class T, IndexStyle I> 
    inline TMV_RealType(T) MaxAbsElement(const VectorView<T,I>& v)
    { return v.maxAbsElement(); }

    template <class T, IndexStyle I> 
    inline TMV_RealType(T) MaxAbsElement(const Vector<T,I>& v)
    { return v.maxAbsElement(); }

    template <class T> 
    inline ConstVectorView<T> Conjugate(const GenVector<T>& v)
    { return v.conjugate(); }

    template <class T, IndexStyle I> 
    inline ConstVectorView<T,I> Conjugate(const ConstVectorView<T,I>& v)
    { return v.conjugate(); }

    template <class T, IndexStyle I> 
    inline ConstVectorView<T,I> Conjugate(const Vector<T,I>& v)
    { return v.conjugate(); }

    template <class T, IndexStyle I> 
    inline VectorView<T,I> Conjugate(const VectorView<T,I>& v)
    { return v.conjugate(); }

    template <class T, IndexStyle I> 
    inline VectorView<T,I> Conjugate(Vector<T,I>& v)
    { return v.conjugate(); }


    //
    // Vector ==, != Vector
    //

    template <class T1, class T2> 
    bool operator==(const GenVector<T1>& v1, const GenVector<T2>& v2);

    template <class T1, class T2> 
    inline bool operator!=(const GenVector<T1>& v1, const GenVector<T2>& v2)
    { return !(v1 == v2); }

    //
    // I/O
    //

    template <class T> 
    inline std::ostream& operator<<(std::ostream& os, const GenVector<T>& v)
    { v.write(os); return os;}

    template <class T> 
    std::istream& operator>>(std::istream& is, const VectorView<T>& v);

    template <class T, IndexStyle I> 
    inline std::istream& operator>>(std::istream& is, Vector<T,I>& v)
    { return is >> v.view(); }

    template <class T, IndexStyle I> 
    std::istream& operator>>(std::istream& is, auto_ptr<Vector<T,I> >& v);

} // namespace tmv

#endif
