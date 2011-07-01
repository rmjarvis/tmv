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
// This file defines the BaseVector, BaseVector_Calc, and BaseVector_Mutable
// classes.
//
// BaseVector is the base class for all of the various types of vectors.
// It defines all of the functions that you can use on any vector.
//
// There is a single template argument, which is the type of the most
// derived class for the object.
// So, for example, Vector<double> inherits from BaseVector<Vector<double> >.
//
// This is a pretty clever meta-programming trick that allows us to 
// get most of the functionality of virtual methods, but without the 
// expensive vtable lookups.  
// (I didn't invent it, so I'm allowed to say it is clever.  :)  )
//
// Basically the way it works is that BaseVector has a method called
// vec() that returns a down-cast to V:
//
// const V& vec() const { return static_cast<const V&>(*this); }
//
// Then, for any method foo(), where you would normally use virtual 
// inheritance, you instead call vec().foo(). 
//
// The methods that are defined for BaseVector are described in 
// TMV_Vector.h and include all of the constant access methods,
// the non-modifying functions, write, and arithmetic operators.
//
//
// BaseVector_Calc is the base class for all vectors that have their
// data calculated in memory, possibly with a non-unit step size, and
// possibly being the conjugate of the underlying data.
//
// The methods that are defined for BaseVector_Calc are described in
// TMV_Vector.h and include all of the vector views and iterators.
//
//
// BaseVector_Mutable is the base class for all vectors with calculated
// data (like BaseVector_Calc) whose values are allowed to be modified.
// 
// The methods that are defined for BaseVector_Calc are described in
// TMV_Vector.h and include all of the modifying functions, and the
// non-const access, read, vector views, and iterators.
//

#ifndef TMV_BaseVector_H
#define TMV_BaseVector_H

#include "TMV_Base.h"
#include "TMV_Shape.h"
#include "TMV_ListInit.h"

namespace tmv {

    // BaseVector is the base class for all vector classes.
    // All non-modifying functions are defined for BaseVector, such as
    // norm, sumElements, minElement, etc.
    template <class V>
    class BaseVector;

    // BaseVector gets a lot of the relevant information about V 
    // from the Traits class: Traits<V>
    // The types are defined as typedef statements.
    // The integer (or boolean) values are defined as enum statements.
    // See TMV_Vector.h or TMV_SmallVector.h for some concrete examples.
    //
    //
    // value_type = The type of the individual elements
    //
    // type = shorthand for the derived type 
    //
    // calc_type = The type of the calculated version of the vector.
    // i.e. where all of the values are stored in memory somewhere.
    // This is the return type of calc()
    // Use this when you will access elements of a composite vector 
    // multiple times or need all of the elements.
    // Also use this when you will use cptr().
    //
    // eval_type = The type of the evaluated version of the vector
    // This is the return type of eval()
    // Use this when you will only access a few elements of a 
    // composite vector one time each.
    //
    // copy_type = The type of a new copy of the vector
    // This is the return type of copy()
    // Use this type when you need elements in the original vector
    // after overwriting those elements.
    //
    // _size = size of vector (Use UNKNOWN if unknown at compile time)
    //
    // _fort = does the indexing use fortran style?
    //
    // _calc = are the element values already calculated in memory?


    // BaseVector_Calc is derived from BaseVector, and is used
    // for vectors that have their values already calculated somewhere.
    // So composite classes inherit directly from BaseVector, rather
    // than BaseVector_Calc.
    template <class V>
    class BaseVector_Calc;

    // BaseVector_Calc adds some more requirements to the Traits<V> class:
    //
    //  _step = the step size if known (else UNKNOWN)
    //  _conj = is the vector the conjugate of the underlying data?
    //  _checkalias = do we need to check this vector for aliases?
    //
    //  const_subvector_type = return type from subVector(i1,i2) const
    //  const_subvector_step_type = return type from 
    //      subVector(i1,i2,istep) const
    //  const_view_type = return type from view() const
    //  const_cview_type = return type from cView() const
    //  const_fview_type = return type from fView() const
    //  const_xview_type = return type from xView() const
    //  const_unitview_type = return type from unitView() const
    //  const_conjugate_type = return type from conjugate() const
    //  const_reverse_type = return type from reverse() const
    //  const_realpart_type = return type from realPart() const
    //  const_imagpart_type = return type from imagPart() const
    //  const_flatten_type = return type from flatten() const
    //  const_nonconj_type = return type from nonConj() const
    //  nonconst_type = return type from nonConst() const
    //  
    //  const_iterator = return type from begin(), end() const
    //  const_reverse_iterator = return type from rbegin(), rend() const


    // BaseVector_Mutable is derived from BaseVector_Calc, and is used
    // for vectors that are allowed to have their data modified.
    template <class V>
    class BaseVector_Mutable;

    // BaseVector_Mutable adds some more requirements to the Traits<V> class:
    //
    //  reference = return type of v(i)
    //
    //  subvector_type = return type from subVector(i1,i2) 
    //  subvector_step_type = return type from subVector(i1,i2,istep) 
    //  view_type = return type from view() 
    //  cview_type = return type from cView()
    //  fview_type = return type from fView()
    //  xview_type = return type from xView()
    //  unitview_type = return type from unitView()
    //  conjugate_type = return type from conjugate() 
    //  reverse_type = return type from reverse() 
    //  realpart_type = return type from realPart()
    //  imagpart_type = return type from imagPart() 
    //  flatten_type = return type from flatten() 
    //  nonconj_type = return type from nonConj()
    //  
    //  iterator = return type from begin(), end() 
    //  reverse_iterator = return type from rbegin(), rend() 

    // The following all derive from BaseVector_Mutable
    // See TMV_Vector.h and TMV_SmallVector.h for their definitions:
    template <class T, int A=0>
    class Vector;
    template <class T, int A=0>
    class ConstVectorView;
    template <class T, int A=0>
    class VectorView;
    template <class T, int N, int A=0>
    class SmallVector;
    template <class T, int N, int S=UNKNOWN, int A=0>
    class ConstSmallVectorView;
    template <class T, int N, int S=UNKNOWN, int A=0>
    class SmallVectorView;

    // Used by sort(p)
    class Permutation;


    //
    // Helper functions and values:
    //

    // Values used to define how to sort a vector.
    enum ADType { Ascend, Descend };
    enum CompType {
        RealComp, AbsComp, Abs2Comp, ImagComp,
        ArgComp, NormComp, ValueComp, InverseComp, LogComp };

    // A helper class to pick out a component of a possibly complex value
    template <CompType comp, class T>
    struct Component;
    // First the real version:
    template <class T>
    struct Component<RealComp,T>
    {
        typedef T ret_type;
        // return f(x);
        static TMV_INLINE T f(const T& x) { return x; } 
        // apply x = f(x)
        // (usually assigning the result to the real part of x when x is complex
        static TMV_INLINE void applyf(T& x) { }
        // get the applied value of f(x) after doing applyf
        // usually this returns real(x), but sometimes the whole x
        static TMV_INLINE T get(const T& x) { return x; }
    };
    template <class T>
    struct Component<ImagComp,T>
    {
        typedef T ret_type;
        typedef T base_type;
        static TMV_INLINE T f(const T& x) { return T(0); } 
        static TMV_INLINE void applyf(T& x) { x = T(0); }
        static TMV_INLINE T get(const T& x) { return x; }
    };
    template <class T>
    struct Component<AbsComp,T>
    {
        typedef typename Traits<T>::float_type ret_type;
        typedef ret_type base_type;
        static TMV_INLINE ret_type f(const T& x) { return TMV_ABS(x); } 
        static TMV_INLINE void applyf(T& x) { x = TMV_ABS(x); }
        static TMV_INLINE T get(const T& x) { return x; }
    };
    template <class T>
    struct Component<Abs2Comp,T>
    {
        typedef T ret_type;
        typedef T base_type;
        static TMV_INLINE T f(const T& x) { return TMV_ABS2(x); } 
        static TMV_INLINE void applyf(T& x) { x = TMV_ABS2(x); }
        static TMV_INLINE T get(const T& x) { return x; }
    };
    template <class T>
    struct Component<ArgComp,T>
    {
        typedef typename Traits<T>::float_type ret_type;
        typedef ret_type base_type;
        static TMV_INLINE ret_type f(const T& x) { return TMV_ARG(x); } 
        static TMV_INLINE void applyf(T& x) { x = TMV_ARG(x); }
        static TMV_INLINE T get(const T& x) { return x; }
    };
    template <class T>
    struct Component<NormComp,T>
    {
        typedef T ret_type;
        typedef T base_type;
        static TMV_INLINE T f(const T& x) { return TMV_NORM(x); } 
        static TMV_INLINE void applyf(T& x) { x = TMV_NORM(x); }
        static TMV_INLINE T get(const T& x) { return x; }
    };
    template <class T>
    struct Component<ValueComp,T>
    {
        typedef T ret_type;
        typedef T base_type;
        static TMV_INLINE T f(const T& x) { return x; } 
        static TMV_INLINE void applyf(T& x) { }
        static TMV_INLINE T get(const T& x) { return x; }
    };
    template <class T>
    struct Component<InverseComp,T>
    {
        typedef typename Traits<T>::float_type ret_type;
        typedef ret_type base_type;
        static TMV_INLINE ret_type f(const T& x) { return T(1) / x; } 
        static TMV_INLINE void applyf(T& x) { x = T(1) / x; }
        static TMV_INLINE T get(const T& x) { return x; }
    };
    template <class T>
    struct Component<LogComp,T>
    {
        typedef typename Traits<T>::float_type ret_type;
        typedef ret_type base_type;
        static TMV_INLINE ret_type f(const T& x) { return TMV_LOG(x); } 
        static TMV_INLINE void applyf(T& x) { x = TMV_LOG(x); }
        static TMV_INLINE T get(const T& x) { return x; }
    };

    // Now the complex version:
#define CT std::complex<T>
    template <class T>
    struct Component<RealComp,CT>
    {
        typedef T ret_type;
        typedef CT base_type;
        static TMV_INLINE T f(const CT& x) { return real(x); } 
        static TMV_INLINE void applyf(CT& x) { }
        static TMV_INLINE T get(const T& x) { return real(x); }
    };
    template <class T>
    struct Component<ImagComp,CT>
    {
        typedef T ret_type;
        typedef CT base_type;
        static TMV_INLINE T f(const CT& x) { return imag(x); } 
        static TMV_INLINE void applyf(CT& x) { TMV_REAL_PART(x) = imag(x); }
        static TMV_INLINE T get(const CT& x) { return real(x); }
    };
    template <class T>
    struct Component<AbsComp,CT>
    {
        typedef typename Traits<T>::float_type ret_type;
        typedef typename Traits<CT>::float_type base_type;
        static TMV_INLINE ret_type f(const CT& x) { return TMV_ABS(x); } 
        static TMV_INLINE void applyf(CT& x) { TMV_REAL_PART(x) = TMV_ABS(x); }
        static TMV_INLINE T get(const CT& x) { return real(x); }
    };
    template <class T>
    struct Component<Abs2Comp,CT>
    {
        typedef T ret_type;
        typedef CT base_type;
        static TMV_INLINE T f(const CT& x) { return TMV_ABS2(x); } 
        static TMV_INLINE void applyf(CT& x) { TMV_REAL_PART(x) = TMV_ABS2(x); }
        static TMV_INLINE T get(const CT& x) { return real(x); }
    };
    template <class T>
    struct Component<ArgComp,CT>
    {
        typedef typename Traits<T>::float_type ret_type;
        typedef typename Traits<CT>::float_type base_type;
        static TMV_INLINE ret_type f(const CT& x) { return TMV_ARG(x); } 
        static TMV_INLINE void applyf(CT& x) { TMV_REAL_PART(x) = TMV_ARG(x); }
        static TMV_INLINE T get(const CT& x) { return real(x); }
    };
    template <class T>
    struct Component<NormComp,CT>
    {
        typedef T ret_type;
        typedef CT base_type;
        static TMV_INLINE T f(const CT& x) { return TMV_NORM(x); } 
        static TMV_INLINE void applyf(CT& x) { TMV_REAL_PART(x) = TMV_NORM(x); }
        static TMV_INLINE T get(const CT& x) { return real(x); }
    };
    template <class T>
    struct Component<ValueComp,CT>
    {
        typedef CT ret_type;
        typedef CT base_type;
        static TMV_INLINE CT f(const CT& x) { return x; } 
        static TMV_INLINE void applyf(CT& x) { }
        static TMV_INLINE CT get(const CT& x) { return x; }
    };
    template <class T>
    struct Component<InverseComp,CT>
    {
        typedef typename Traits<CT>::float_type ret_type;
        typedef typename Traits<CT>::float_type base_type;
        static TMV_INLINE ret_type f(const CT& x) 
        { return ret_type(std::conj(x)) / TMV_NORM(x); }
        static TMV_INLINE void applyf(CT& x) { x /= TMV_NORM(x); }
        static TMV_INLINE CT get(const CT& x) { return std::conj(x); }
    };
    template <class T>
    struct Component<LogComp,CT>
    {
        typedef typename Traits<T>::float_type ret_type;
        typedef typename Traits<CT>::float_type base_type;
        static TMV_INLINE ret_type f(const CT& x) { return TMV_LOG(real(x)); } 
        static TMV_INLINE void applyf(CT& x) 
        { TMV_REAL_PART(x) = TMV_LOG(real(x)); }
        static TMV_INLINE T get(const CT& x) { return real(x); }
    };
#undef CT

    // These helper functions check the validity of indices according
    // to whether the vector uses CStyle or FortranStyle indexing.
    // They also update the indices to be consistent with CStyle.
    template <bool _fort>
    static TMV_INLINE_ND void CheckIndex(int& i, int n) 
    { TMVAssert(i>=0 && i<n && "index is not valid"); } // CStyle
    template <>
    TMV_INLINE_ND void CheckIndex<true>(int& i, int n) 
    { TMVAssert(i>=1 && i<=n && "index is not valid"); --i; } // FortranStyle
    template <bool _fort>
    static TMV_INLINE_ND void CheckRange(int& i1, int i2, int n)
    { // CStyle
        TMVAssert(i1 >= 0 && "first element must be in range");
        TMVAssert(i2 <= n && "last element must be in range");
        TMVAssert(i2 >= i1 && 
                  "range must have a non-negative number of elements");
    }
    template <>
    TMV_INLINE_ND void CheckRange<true>(int& i1, int i2, int n)
    { // FortranStyle
        TMVAssert(i1 >= 1 && "first element must be in range");
        TMVAssert(i2 <= n && "last element must be in range");
        TMVAssert(i2 >= i1 && "range must have a positive number of elements");
        --i1;
    }
    template <bool _fort>
    static TMV_INLINE_ND void CheckRange(int& i1, int& i2, int istep, int n)
    { // CStyle
        TMVAssert(istep != 0 && "istep cannot be 0");
        TMVAssert(((i1 >= 0 && i1 < n) || i1==i2) && 
                  "first element must be in range");
        TMVAssert(((i2-istep >= 0 && i2-istep < n) || i1==i2) &&
                  "last element must be in range");
        TMVAssert((i2-i1) % istep == 0 &&
                  "range must be an integral multiple of istep");
        TMVAssert((i2-i1) / istep >= 0 &&
                  "must have a non-negative number of elements");
    }
    template <>
    TMV_INLINE_ND void CheckRange<true>(int& i1, int& i2, int istep, int n)
    { // FortranStyle
        TMVAssert(istep != 0 && "istep cannot be 0");
        TMVAssert(i1 >= 1 && i1 <= n && "first element must be in range");
        TMVAssert(i2 >= 1 && i2 <= n && "last element must be in range");
        TMVAssert((i2-i1) % istep == 0 &&
                  "subvector range must be an integral multiple of istep");
        TMVAssert((i2-i1) / istep >= 0 &&
                  "subvector range must have a positive number of elements");
        --i1; i2 += istep-1;
    }

    // This helper class checks two compile time sizes.
    // same = Are they possibly the same size?
    // equal = Are they definitely the same size?
    // size = If knowable, what is that size?
    template <int S1, int S2>
    struct Sizes
    {
        enum { same = (S1==UNKNOWN || S2==UNKNOWN || S1==S2) };
        enum { equal = (S1!=UNKNOWN && S2!=UNKNOWN && S1==S2) };
        enum { size = (S1==UNKNOWN ? S2 : S1) };
    };

    // This helper class helps decide calc_type for composite classes:
    template <class T, int s, bool fort=false>
    struct VCopyHelper
    {
        enum { A2 = (fort ? FortranStyle : CStyle) | NoAlias };
        typedef SmallVector<T,s,A2> type; 
    };
    template <class T, bool fort>
    struct VCopyHelper<T,UNKNOWN,fort>
    {
        enum { A2 = (fort ? FortranStyle : CStyle) | NoAlias };
        typedef Vector<T,A2> type; 
    };

    // This is similar - it defines the right view type when the
    // size or step might be known.
    template <class T, int N, int S, int C=NonConj>
    struct VViewHelper
    { 
        typedef SmallVectorView<T,N,S,C> type; 
        typedef ConstSmallVectorView<T,N,S,C> ctype; 
    };
    template <class T, int S, int C>
    struct VViewHelper<T,UNKNOWN,S,C>
    {
        enum { A = C | (S == 1 ? Unit : NonUnit) };
        typedef VectorView<T,A> type; 
        typedef ConstVectorView<T,A> ctype; 
    };


    // This helper function checks for aliasis in the memory storage of
    // two objects.  We overload it for specific objects that can be 
    // checked to have the same realPart().cptr() values.
    template <class V1, class V2>
    static TMV_INLINE bool SameStorage(
        const BaseVector<V1>& v1, const BaseVector<V2>& v2)
    { return false; }
#ifndef TMV_NO_ALIAS_CHECK
    template <class V1, class V2>
    static TMV_INLINE bool SameStorage(
        const BaseVector_Calc<V1>& v1, const BaseVector_Calc<V2>& v2)
    {
        return static_cast<const void*>(v1.vec().cptr()) == 
            static_cast<const void*>(v2.vec().cptr()); 
    }
#endif

    // This next one is sometime used after SameStorage is found to be true.
    // Here we check whether the elements are exactly the same 
    // by whether the steps are the same.
    // We do not check the _conj values, so that is usually the next 
    // step depending on why we are checking this.
    template <class V1, class V2>
    static TMV_INLINE bool ExactSameStorage(
        const BaseVector<V1>& v1, const BaseVector<V2>& v2)
    { return false; }
    template <class V1, class V2>
    static TMV_INLINE bool ExactSameStorage(
        const BaseVector_Calc<V1>& v1, const BaseVector_Calc<V2>& v2)
    {
        typedef typename V1::value_type T1;
        typedef typename V2::value_type T2;
        return Traits2<T1,T2>::sametype && v1.step() == v2.step(); 
    }

    // This helper class determines if there is ExactSameStorage
    // for two vectors at compile time:
    template <class V1, class V2>
    struct VStepHelper
    {
        enum { known = (
                V1::_step != UNKNOWN &&
                V2::_step != UNKNOWN ) };
        enum { same = (
                known &&
                V1::_step == int(V2::_step) ) };
        enum { noclobber = (
                known &&
                ( same ||
                  (V2::_step > 0 && V1::_step > int(V2::_step)) ||
                  (V2::_step < 0 && V1::_step < int(V2::_step)) ) ) };
    };

    // Defined in TMV_CopyV.h
    template <class V1, class V2>
    static inline void Copy(
        const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2);
    template <class V1, class V2>
    static inline void NoAliasCopy(
        const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2);

    // Defined in TMV_SwapV.h
    template <class V1, class V2>
    static inline void Swap(
        BaseVector_Mutable<V1>& v1, BaseVector_Mutable<V2>& v2);
    template <class V1, class V2>
    static inline void NoAliasSwap(
        BaseVector_Mutable<V1>& v1, BaseVector_Mutable<V2>& v2);
    template <class V>
    static inline void ReverseSelf(BaseVector_Mutable<V>& v);

    // Defined in TMV_ConjugateV.h
    template <class V>
    static inline void ConjugateSelf(BaseVector_Mutable<V>& v);

    // Defined in TMV_NormV.h
    template <class V>
    static inline typename V::real_type DoNormSq(const BaseVector_Calc<V>& v);
    template <class V>
    static inline typename V::float_type DoNormSq(
        const BaseVector_Calc<V>& v, const typename V::float_type scale);
    template <class V>
    static inline typename V::value_type DoSumElements(
        const BaseVector_Calc<V>& v);
    template <class V>
    static inline typename V::float_type DoSumAbsElements(
        const BaseVector_Calc<V>& v);
    template <class V>
    static inline typename V::real_type DoSumAbs2Elements(
        const BaseVector_Calc<V>& v);

    // Defined in TMV_Norm.h
    template <class V>
    static inline typename V::float_type DoNorm2(const BaseVector_Calc<V>& v);

    // Defined in TMV_Det.h
    template <class V>
    static inline typename V::value_type ProdElements(
        const BaseVector_Calc<V>& v);
    template <class V>
    static inline typename V::float_type LogProdElements(
        const BaseVector_Calc<V>& v, typename V::zfloat_type* sign);
    template <class V>
    static inline bool HasZeroElement(const BaseVector_Calc<V>& v);

    // Defined in TMV_MinMax.h
    template <class V>
    static inline typename V::value_type DoMaxElement(
        const BaseVector_Calc<V>& v, int* imax=0);
    template <class V>
    static inline typename V::float_type DoMaxAbsElement(
        const BaseVector_Calc<V>& v, int* imax=0);
    template <class V>
    static inline typename V::real_type DoMaxAbs2Element(
        const BaseVector_Calc<V>& v, int* imax=0);
    template <class V>
    static inline typename V::value_type DoMinElement(
        const BaseVector_Calc<V>& v, int* imin=0);
    template <class V>
    static inline typename V::float_type DoMinAbsElement(
        const BaseVector_Calc<V>& v, int* imin=0);
    template <class V>
    static inline typename V::real_type DoMinAbs2Element(
        const BaseVector_Calc<V>& v, int* imin=0);

    // Defined in TMV_VectorIO.h
    template <class V>
    static inline void Write(std::ostream& os, const BaseVector_Calc<V>& v);
    template <class V>
    static inline void Write(
        std::ostream& os, const BaseVector_Calc<V>& v,
        typename V::float_type thresh) ;
    template <class V>
    static inline void Read(std::istream& is, BaseVector_Mutable<V>& v);

    // Defined in TMV_SortV.h
    template <class V>
    static inline void Sort(BaseVector_Mutable<V>& v, ADType ad, CompType comp);
    template <class V>
    static inline void Sort(
        BaseVector_Mutable<V>& v, int* P, ADType ad, CompType comp);

    // A helper class for returning views without necessarily
    // making a new object.
    template <bool ref, class type, class view_type>
    struct MakeVecView_Helper;

    template <class type, class view_type>
    struct MakeVecView_Helper<true,type,view_type>
    {
        typedef type& ret_type;
        typedef const type& const_ret_type;
        static TMV_INLINE ret_type call(type& v) { return v; }
        static TMV_INLINE const_ret_type call(const type& v) { return v; }
    };

    template <class type, class view_type>
    struct MakeVecView_Helper<false,type,view_type>
    {
        typedef view_type ret_type;
        typedef view_type const_ret_type;
        static TMV_INLINE ret_type call(type& v) 
        { return view_type(v.ptr(),v.size(),v.step()); }
        static TMV_INLINE const_ret_type call(const type& v) 
        { return view_type(v.cptr(),v.size(),v.step()); }
    };

    template <class type, class view_type>
    struct MakeVecView
    {
        enum { ref = Traits2<type,view_type>::sametype };
        typedef MakeVecView_Helper<ref,type,view_type> helper;

        static TMV_INLINE typename helper::ret_type call(type& v)
        { return helper::call(v); }
        static TMV_INLINE typename helper::const_ret_type call(const type& v)
        { return helper::call(v); }
    };

 
    //
    // BaseVector
    //

    template <class V>
    class BaseVector
    {
    public:
        enum { _size = Traits<V>::_size };
        enum { _fort = Traits<V>::_fort };
        enum { _calc = Traits<V>::_calc };
        enum { _shape = Vec };

        typedef V type;

        typedef typename Traits<V>::value_type value_type;
        typedef typename Traits<V>::calc_type calc_type;
        typedef typename Traits<V>::eval_type eval_type;
        typedef typename Traits<V>::copy_type copy_type;

        // Derived values:
        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<real_type>::float_type float_type;
        typedef typename Traits<value_type>::float_type zfloat_type;
        typedef typename Traits<value_type>::complex_type complex_type;
        enum { isreal = Traits<value_type>::isreal };
        enum { iscomplex = Traits<value_type>::iscomplex };

        //
        // Constructor
        //

        TMV_INLINE BaseVector() {}
        TMV_INLINE BaseVector(const BaseVector<V>&) {}
        TMV_INLINE ~BaseVector() {}

        //
        // Access
        // 

        TMV_INLINE value_type operator[](int i) const 
        { return operator()(i); }
        TMV_INLINE value_type operator()(int i) const 
        {
            CheckIndex<_fort>(i,size());
            return cref(i);
        }

        //
        // Functions
        //

        TMV_INLINE value_type sumElements() const
        { return tmv::DoSumElements(calc()); }

        TMV_INLINE float_type sumAbsElements() const 
        { return tmv::DoSumAbsElements(calc()); }

        TMV_INLINE real_type sumAbs2Elements() const
        { return tmv::DoSumAbs2Elements(calc()); }

        TMV_INLINE value_type maxElement(int* imax=0) const
        { return tmv::DoMaxElement(calc(),imax); }

        TMV_INLINE float_type maxAbsElement(int* imax=0) const 
        { return tmv::DoMaxAbsElement(calc(),imax); }

        TMV_INLINE real_type maxAbs2Element(int* imax=0) const
        { return tmv::DoMaxAbs2Element(calc(),imax); }

        TMV_INLINE value_type minElement(int* imin=0) const
        { return tmv::DoMinElement(calc(),imin); }

        TMV_INLINE float_type minAbsElement(int* imin=0) const 
        { return tmv::DoMinAbsElement(calc(),imin); }

        TMV_INLINE real_type minAbs2Element(int* imin=0) const
        { return tmv::DoMinAbs2Element(calc(),imin); }

        TMV_INLINE float_type norm1() const
        { return sumAbsElements(); }

        TMV_INLINE real_type normSq() const
        { return tmv::DoNormSq(calc()); }

        TMV_INLINE float_type normSq(const float_type scale) const
        { return tmv::DoNormSq(calc(),scale); }

        TMV_INLINE float_type norm2() const
        { return tmv::DoNorm2(calc()); }

        TMV_INLINE float_type norm() const
        { return norm2(); }

        TMV_INLINE float_type normInf() const
        { return maxAbsElement(); }

        TMV_INLINE value_type prodElements() const
        { return tmv::ProdElements(calc()); }

        TMV_INLINE float_type logProdElements(zfloat_type* sign=0) const
        { return tmv::LogProdElements(calc(),sign); }

        TMV_INLINE bool hasZeroElement() const
        { return tmv::HasZeroElement(calc()); }



        // 
        // I/O
        //

        TMV_INLINE void write(std::ostream& os) const
        { tmv::Write(os,calc()); }
        TMV_INLINE void write(std::ostream& os, float_type thresh) const
        { tmv::Write(os,calc(),thresh); }



        //
        // Auxilliary routines
        //

        TMV_INLINE const type& vec() const 
        { return static_cast<const type&>(*this); }

        TMV_INLINE calc_type calc() const 
        { return static_cast<calc_type>(vec()); }

        TMV_INLINE eval_type eval() const 
        { return static_cast<eval_type>(vec()); }

        TMV_INLINE copy_type copy() const 
        { return static_cast<copy_type>(vec()); }

        // Note that these last functions need to be defined in a more derived
        // class than this, or an infinite loop will result when compiling.

        TMV_INLINE size_t size() const { return vec().size(); }
        TMV_INLINE int nElements() const { return vec().nElements(); }

        TMV_INLINE value_type cref(int i) const  { return vec().cref(i); }

        template <class V2>
        TMV_INLINE void assignTo(BaseVector_Mutable<V2>& v2) const
        { vec().assignTo(v2); }

        template <class V2>
        TMV_INLINE void newAssignTo(BaseVector_Mutable<V2>& v2) const
        { vec().newAssignTo(v2); }

    private :
        void operator=(const BaseVector<V>& v2);

    }; // BaseVector

    //
    // BaseVector_Calc
    //

    template <class V>
    class BaseVector_Calc : 
        public BaseVector<V>
    {
    public:
        enum { _size = Traits<V>::_size };
        enum { _fort = Traits<V>::_fort };
        enum { _calc = Traits<V>::_calc };
        enum { _step = Traits<V>::_step };
        enum { _conj = Traits<V>::_conj };
        enum { _checkalias = Traits<V>::_checkalias };

        typedef V type;
        typedef BaseVector<V> base;

        typedef typename base::calc_type calc_type;
        typedef typename base::eval_type eval_type;
        typedef typename base::copy_type copy_type;

        typedef typename base::value_type value_type;
        typedef typename base::real_type real_type;
        typedef typename base::complex_type complex_type;
        typedef typename base::float_type float_type;
        typedef typename base::zfloat_type zfloat_type;

        typedef typename Traits<V>::const_view_type const_view_type;
        typedef typename Traits<V>::const_cview_type const_cview_type;
        typedef typename Traits<V>::const_fview_type const_fview_type;
        typedef typename Traits<V>::const_xview_type const_xview_type;
        typedef typename Traits<V>::const_unitview_type const_unitview_type;

        typedef typename Traits<V>::const_subvector_type const_subvector_type;
        typedef typename Traits<V>::const_subvector_step_type 
            const_subvector_step_type;
        typedef typename Traits<V>::const_conjugate_type const_conjugate_type;
        typedef typename Traits<V>::const_reverse_type const_reverse_type;
        typedef typename Traits<V>::const_realpart_type const_realpart_type;
        typedef typename Traits<V>::const_imagpart_type const_imagpart_type;
        typedef typename Traits<V>::const_flatten_type const_flatten_type;
        typedef typename Traits<V>::const_nonconj_type const_nonconj_type;
        typedef typename Traits<V>::nonconst_type nonconst_type;

        typedef typename Traits<V>::const_iterator const_iterator;
        typedef typename Traits<V>::const_reverse_iterator 
            const_reverse_iterator;


        //
        // Constructor
        //

        TMV_INLINE BaseVector_Calc() {}
        TMV_INLINE BaseVector_Calc(const BaseVector_Calc<V>&) {}
        TMV_INLINE ~BaseVector_Calc() {}


        //
        // Access 
        //

        TMV_INLINE const_iterator begin() const
        { return const_iterator(cptr(),step()); }
        TMV_INLINE const_iterator end() const
        { return begin() + size(); }
        TMV_INLINE const_reverse_iterator rbegin() const
        { return const_reverse_iterator(cptr()+(size()-1)*step(),-step()); }
        TMV_INLINE const_reverse_iterator rend() const
        { return rbegin() + size(); }


        //
        // subVector
        //

        // cSubVector always uses CStyle
        const_subvector_type cSubVector(int i1, int i2) const
        { return const_subvector_type(cptr()+i1*step(),i2-i1,step()); }

        const_subvector_step_type cSubVector(
            int i1, int i2, int istep) const
        {
            return const_subvector_step_type(
                cptr()+i1*step(), (i2-i1)/istep, istep*step());
        }

        TMV_INLINE const_subvector_type subVector(int i1, int i2) const
        {
            CheckRange<_fort>(i1,i2,size());
            return cSubVector(i1,i2);
        }

        TMV_INLINE const_subvector_step_type subVector(
            int i1, int i2, int istep) const
        {
            CheckRange<_fort>(i1,i2,istep,size());
            return cSubVector(i1,i2,istep);
        }


        //
        // Views
        //

        TMV_INLINE TMV_MAYBE_CREF(type,const_view_type) view() const
        { return MakeVecView<type,const_view_type>::call(vec()); }

        TMV_INLINE TMV_MAYBE_CREF(type,const_cview_type) cView() const
        { return MakeVecView<type,const_cview_type>::call(vec()); }

        TMV_INLINE TMV_MAYBE_CREF(type,const_fview_type) fView() const
        { return MakeVecView<type,const_fview_type>::call(vec()); }

        TMV_INLINE TMV_MAYBE_CREF(type,const_xview_type) xView() const
        { return MakeVecView<type,const_xview_type>::call(vec()); }

        TMV_INLINE_ND TMV_MAYBE_CREF(type,const_unitview_type) unitView() const
        {
            TMVAssert(step() == 1 && "Called unitView on vector with step!=1");
            return MakeVecView<type,const_unitview_type>::call(vec()); 
        }

        TMV_INLINE TMV_MAYBE_CREF(type,const_view_type) constView() const
        { return MakeVecView<type,const_view_type>::call(vec()); }

        TMV_INLINE TMV_MAYBE_CREF(type,const_conjugate_type) conjugate() const
        { return MakeVecView<type,const_conjugate_type>::call(vec()); }

        const_reverse_type reverse() const
        {
            return const_reverse_type(
                cptr()+(int(size())-1)*step(), size(), -step());
        }

        const_realpart_type realPart() const
        {
            return const_realpart_type(
                reinterpret_cast<const real_type*>(cptr()),
                size(), V::isreal ? step() : 2*step());
        }

        const_imagpart_type imagPart() const
        {
            TMVStaticAssert(V::iscomplex);
            return const_imagpart_type(
                reinterpret_cast<const real_type*>(cptr())+1, size(), 2*step());
        }

        const_flatten_type flatten() const
        {
            TMVStaticAssert(_step == UNKNOWN || _step == 1);
            TMVAssert(step() == 1);
            return const_flatten_type(
                reinterpret_cast<const real_type*>(cptr()),
                V::isreal ? size() : 2*size(), 1);
        }

        TMV_INLINE TMV_MAYBE_CREF(type,const_nonconj_type) nonConj() const
        { return MakeVecView<type,const_nonconj_type>::call(vec()); }

        TMV_INLINE nonconst_type nonConst() const
        {
            return nonconst_type(
                const_cast<value_type*>(cptr()),size(),step());
        }

        //
        // Auxilliary routines
        //

        template <class V2>
        TMV_INLINE void assignTo(BaseVector_Mutable<V2>& v2) const
        { tmv::Copy(vec(),v2); }

        template <class V2>
        TMV_INLINE void newAssignTo(BaseVector_Mutable<V2>& v2) const
        { tmv::NoAliasCopy(vec(),v2); }

        TMV_INLINE const type& vec() const
        { return static_cast<const type&>(*this); }

        TMV_INLINE bool isconj() const { return _conj; }

        // Note that these last functions need to be defined in a more derived
        // class than this, or an infinite loop will result when compiling.

        TMV_INLINE size_t size() const { return vec().size(); }
        TMV_INLINE int step() const { return vec().step(); }
        TMV_INLINE const value_type* cptr() const { return vec().cptr(); }
        TMV_INLINE value_type cref(int i) const  { return vec().cref(i); }

    private:
        void operator=(const BaseVector_Calc<V>&);

    }; // BaseVector_Calc

    //
    // BaseVector_Mutable
    //

    template <class V>
    class BaseVector_Mutable : 
        public BaseVector_Calc<V>
    {
    public:
        enum { _size = Traits<V>::_size };
        enum { _fort = Traits<V>::_fort };
        enum { _calc = Traits<V>::_calc };
        enum { _step = Traits<V>::_step };
        enum { _conj = Traits<V>::_conj };
        enum { _checkalias = Traits<V>::_checkalias };

        typedef V type;
        typedef BaseVector_Calc<V> base;

        typedef typename base::calc_type calc_type;
        typedef typename base::eval_type eval_type;
        typedef typename base::copy_type copy_type;

        typedef typename base::value_type value_type;
        typedef typename base::real_type real_type;
        typedef typename base::complex_type complex_type;
        typedef typename base::float_type float_type;
        typedef typename base::zfloat_type zfloat_type;

        typedef typename base::const_view_type const_view_type;
        typedef typename base::const_cview_type const_cview_type;
        typedef typename base::const_fview_type const_fview_type;
        typedef typename base::const_xview_type const_xview_type;
        typedef typename base::const_unitview_type const_unitview_type;

        typedef typename base::const_subvector_type const_subvector_type;
        typedef typename base::const_subvector_step_type 
            const_subvector_step_type;
        typedef typename base::const_conjugate_type const_conjugate_type;
        typedef typename base::const_reverse_type const_reverse_type;
        typedef typename base::const_realpart_type const_realpart_type;
        typedef typename base::const_imagpart_type const_imagpart_type;
        typedef typename base::const_flatten_type const_flatten_type;
        typedef typename base::const_nonconj_type const_nonconj_type;
        typedef typename base::nonconst_type nonconst_type;

        typedef typename base::const_iterator const_iterator;
        typedef typename base::const_reverse_iterator 
            const_reverse_iterator;

        typedef typename Traits<V>::view_type view_type;
        typedef typename Traits<V>::cview_type cview_type;
        typedef typename Traits<V>::fview_type fview_type;
        typedef typename Traits<V>::xview_type xview_type;
        typedef typename Traits<V>::unitview_type unitview_type;
 
        typedef typename Traits<V>::subvector_type subvector_type;
        typedef typename Traits<V>::subvector_step_type subvector_step_type;
        typedef typename Traits<V>::conjugate_type conjugate_type;
        typedef typename Traits<V>::reverse_type reverse_type;
        typedef typename Traits<V>::realpart_type realpart_type;
        typedef typename Traits<V>::imagpart_type imagpart_type;
        typedef typename Traits<V>::flatten_type flatten_type;
        typedef typename Traits<V>::nonconj_type nonconj_type;

        typedef typename Traits<V>::iterator iterator;
        typedef typename Traits<V>::reverse_iterator reverse_iterator;
        typedef typename Traits<V>::reference reference;

        //
        // Constructor
        //

        TMV_INLINE BaseVector_Mutable() {}
        TMV_INLINE BaseVector_Mutable(const BaseVector_Mutable<V>&) {}
        TMV_INLINE ~BaseVector_Mutable() {}


        //
        // Access 
        //

        TMV_INLINE reference operator[](int i)
        { return operator()(i); }
        TMV_INLINE reference operator()(int i)
        {
            CheckIndex<_fort>(i,size());
            return ref(i);
        }

        TMV_INLINE iterator begin()
        { return iterator(ptr(),step()); }
        TMV_INLINE iterator end()
        { return begin() + size(); }
        TMV_INLINE reverse_iterator rbegin()
        { return reverse_iterator(ptr()+(size()-1)*step(),-step()); }
        TMV_INLINE reverse_iterator rend()
        { return rbegin() + size(); }

        // We need to repeat the const versions so the non-const ones
        // don't clobber them.
        TMV_INLINE value_type operator[](int i) const
        { return base::operator[](i); }
        TMV_INLINE value_type operator()(int i) const
        { return base::operator()(i); }

        TMV_INLINE const_iterator begin() const
        { return base::begin(); }
        TMV_INLINE const_iterator end() const
        { return base::end(); }
        TMV_INLINE const_reverse_iterator rbegin() const
        { return base::rbegin(); }
        TMV_INLINE const_reverse_iterator rend() const
        { return base::rend(); }


        //
        // Op =
        //

        type& operator=(const BaseVector_Mutable<V>& v2) 
        {
            TMVAssert(size() == v2.size());
            v2.assignTo(vec());
            return vec(); 
        }

        template <class V2>
        type& operator=(const BaseVector<V2>& v2) 
        {
            TMVStaticAssert((Sizes<_size,V2::_size>::same));
            TMVAssert(size() == v2.size());
            v2.assignTo(vec());
            return vec(); 
        }

        ListAssigner<value_type,iterator> operator<<(value_type x)
        { return ListAssigner<value_type,iterator>(begin(),size(),x); }


        //
        // Modifying Functions
        //

        type& setZero() 
        {
            const int n=size();
            for(int i=0;i<n;++i) ref(i) = value_type(0);
            return vec();
        }

        type& clip(float_type thresh) 
        {
            const int n=size();
            for(int i=0;i<n;++i) {
                const float_type temp = TMV_ABS(cref(i));
                if (temp < thresh) ref(i) = value_type(0);
            }
            return vec();
        }

        type& setAllTo(value_type x) 
        {
            const int n=size();
            for(int i=0;i<n;++i) ref(i) = x;
            return vec();
        }

        type& addToAll(value_type x) 
        {
            const int n=size();
            for(int i=0;i<n;++i) ref(i) += x;
            return vec();
        }

        TMV_INLINE type& conjugateSelf()
        { tmv::ConjugateSelf(vec()); return vec(); }

        template <class F>
        type& applyToAll(const F& f)
        {
            const int n=size();
            for(int i=0;i<n;++i) ref(i) = f(cref(i));
            return vec();
        }

        type& cMakeBasis(int i, value_type x=value_type(1)) 
        {
            setZero(); ref(i) = x;
            return vec();
        }
        TMV_INLINE type& makeBasis(int i, value_type x=value_type(1)) 
        {
            CheckIndex<_fort>(i,size());
            return cMakeBasis(i,x);
        }

        type& cSwap(int i1, int i2) 
        {
            TMV_SWAP(ref(i1),ref(i2)); 
            return vec();
        }
        TMV_INLINE type& swap(int i1, int i2) 
        {
            CheckIndex<_fort>(i1,size());
            CheckIndex<_fort>(i2,size());
            return cSwap(i1,i2);
        }

        type& cPermute(const int* p, int i1, int i2) 
        {
            for(int i=i1;i<i2;++i) cSwap(i,p[i]); 
            return vec();
        }
        TMV_INLINE type& permute(const int* p, int i1, int i2) 
        {
            CheckRange<_fort>(i1,i2,size());
            return cPermute(p,i1,i2);
        }
        TMV_INLINE type& permute(const int* p) 
        { return cPermute(p,0,size()); }

        type& cReversePermute(const int* p, int i1, int i2) 
        {
            for(int i=i2;i>i1;) { --i; cSwap(i,p[i]); }
            return vec();
        }
        TMV_INLINE type& reversePermute(const int* p, int i1, int i2) 
        {
            CheckRange<_fort>(i1,i2,size());
            return cReversePermute(p,i1,i2);
        }
        TMV_INLINE type& reversePermute(const int* p) 
        { return cReversePermute(p,0,size()); }

        TMV_INLINE type& reverseSelf() 
        { tmv::ReverseSelf(vec()); return vec(); }

        TMV_INLINE type& sort(ADType ad=Ascend, CompType comp=RealComp) 
        { tmv::Sort(vec(),ad,comp); return vec(); }

        // Defined in TMV_Permutation.h
        type& sort(
            Permutation& P, ADType ad=Ascend, CompType comp=RealComp);

        //
        // SubVector
        //

        subvector_type cSubVector(int i1, int i2) 
        { return subvector_type(ptr()+i1*step(),i2-i1,step()); }

        subvector_step_type cSubVector(int i1, int i2, int istep) 
        {
            return subvector_step_type(
                ptr()+i1*step(), (i2-i1)/istep, istep*step());
        }

        TMV_INLINE subvector_type subVector(int i1, int i2) 
        {
            CheckRange<_fort>(i1,i2,size());
            return cSubVector(i1,i2);
        }

        TMV_INLINE subvector_step_type subVector(int i1, int i2, int istep) 
        {
            CheckRange<_fort>(i1,i2,istep,size());
            return cSubVector(i1,i2,istep);
        }

        reverse_type reverse() 
        { return reverse_type(ptr()+(int(size())-1)*step(), size(),-step()); }

        TMV_INLINE TMV_MAYBE_REF(type,view_type) view() 
        { return MakeVecView<type,view_type>::call(vec()); }

        TMV_INLINE TMV_MAYBE_REF(type,cview_type) cView() 
        { return MakeVecView<type,cview_type>::call(vec()); }

        TMV_INLINE TMV_MAYBE_REF(type,fview_type) fView() 
        { return MakeVecView<type,fview_type>::call(vec()); }

        TMV_INLINE TMV_MAYBE_REF(type,xview_type) xView() 
        { return MakeVecView<type,xview_type>::call(vec()); }

        TMV_INLINE TMV_MAYBE_REF(type,unitview_type) unitView() 
        {
            TMVAssert(step() == 1 && "Called unitView on vector with step!=1");
            return MakeVecView<type,unitview_type>::call(vec()); 
        }

        TMV_INLINE TMV_MAYBE_REF(type,conjugate_type) conjugate() 
        { return MakeVecView<type,conjugate_type>::call(vec()); }

        realpart_type realPart() 
        {
            return realpart_type(
                reinterpret_cast<real_type*>(ptr()),
                size(), V::isreal ? step() : 2*step());
        }

        imagpart_type imagPart() 
        {
            TMVStaticAssert(V::iscomplex);
            return imagpart_type(
                reinterpret_cast<real_type*>(ptr())+1, size(), 2*step());
        }

        flatten_type flatten() 
        {
            TMVStaticAssert(_step == UNKNOWN || _step == 1);
            TMVAssert(step() == 1);
            return flatten_type(
                reinterpret_cast<real_type*>(ptr()),
                V::isreal ? size() : 2*size(), 1);
        }

        TMV_INLINE TMV_MAYBE_REF(type,nonconj_type) nonConj()
        { return MakeVecView<type,nonconj_type>::call(vec()); }



        // Repeat the const versions:
        TMV_INLINE const_subvector_type cSubVector(int i1, int i2) const
        { return base::cSubVector(i1,i2); }
        TMV_INLINE const_subvector_step_type cSubVector(
            int i1, int i2, int istep) const
        { return base::cSubVector(i1,i2,istep); }
        TMV_INLINE const_subvector_type subVector(int i1, int i2) const
        { return base::subVector(i1,i2); }
        TMV_INLINE const_subvector_step_type subVector(
            int i1, int i2, int istep) const
        { return base::subVector(i1,i2,istep); }
        TMV_INLINE const_reverse_type reverse() const
        { return base::reverse(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_view_type) view() const
        { return base::view(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_cview_type) cView() const
        { return base::cView(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_fview_type) fView() const
        { return base::fView(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_xview_type) xView() const
        { return base::xView(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_unitview_type) unitView() const
        { return base::unitView(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_conjugate_type) conjugate() const
        { return base::conjugate(); }
        TMV_INLINE const_realpart_type realPart() const
        { return base::realPart(); }
        TMV_INLINE const_imagpart_type imagPart() const
        { return base::imagPart(); }
        TMV_INLINE const_flatten_type flatten() const
        { return base::flatten(); }
        TMV_INLINE TMV_MAYBE_CREF(type,const_nonconj_type) nonConj() const
        { return base::nonConj(); }



        //
        // I/O
        //

        // Defined in TMV_VectorIO.h
        TMV_INLINE void read(std::istream& is)
        { tmv::Read(is,vec()); }


        //
        // Arithmetic
        //

        // These operators need to be here, rather than just defining
        // non-member operator*=, etc., since the argument to a non-member
        // function would have to be BaseVector_Mutable& (ie. a reference).
        // But then you couldn't write something like:
        // m.row(i) *= 3.;
        // since the m.row(i) function returns a view by value, which is not
        // castable to the non-const reference argument of operator*=.
        // So we define all these here with unspecified right hand sides.
        // They just have to be valid arguments to a MultEq, AddEq, etc.
        // function defined as non-member functions for various objects.

        template <class X2>
        TMV_INLINE type& operator+=(const X2& x2)
        { AddEq(vec(),x2); return vec(); }

        template <class X2>
        TMV_INLINE type& operator-=(const X2& x2)
        { SubtractEq(vec(),x2); return vec(); }

        template <class X2>
        TMV_INLINE type& operator*=(const X2& x2)
        { MultEq(vec(),x2); return vec(); }

        template <class X2>
        TMV_INLINE type& operator/=(const X2& x2)
        { LDivEq(vec(),x2); return vec(); }

        template <class X2>
        TMV_INLINE type& operator%=(const X2& x2)
        { RDivEq(vec(),x2); return vec(); }



        //
        // Auxilliary routines
        //

        TMV_INLINE const type& vec() const
        { return static_cast<const type&>(*this); }
        TMV_INLINE type& vec()
        { return static_cast<type&>(*this); }

        // Note that these last functionsneed to be defined in a more derived
        // class than this, or an infinite loop will result when compiling.

        TMV_INLINE size_t size() const { return vec().size(); }
        TMV_INLINE int step() const { return vec().step(); }
        TMV_INLINE value_type* ptr() { return vec().ptr(); }
        TMV_INLINE const value_type* cptr() { return vec().cptr(); }
        TMV_INLINE reference ref(int i) { return vec().ref(i); }
        TMV_INLINE value_type cref(int i) { return vec().cref(i); }

    }; // BaseVector_Mutable

    template <class T>
    class VectorSizer;

    template <class T>
    struct Traits<VectorSizer<T> >
    {
        typedef T value_type;
        typedef VectorSizer<T> type;
        typedef InvalidType calc_type;
        typedef InvalidType eval_type;
        typedef InvalidType copy_type;

        enum { _size = UNKNOWN };
        enum { _shape = Vec };
        enum { _fort = false };
        enum { _calc = false };
    };

    template <class T>
    class VectorSizer : 
        public BaseVector<VectorSizer<T> >
    {
    public:
        TMV_INLINE VectorSizer(const int _s) : s(_s) {}
        TMV_INLINE size_t size() const { return s; }

        TMV_INLINE T cref(int ) const  { return T(0); }

        template <class M2>
        TMV_INLINE void assignTo(BaseVector_Mutable<M2>& ) const {}

        template <class M2>
        TMV_INLINE void newAssignTo(BaseVector_Mutable<M2>& ) const {}

    private :
        const int s;
    }; // VectorSizer


    //
    // Vector ==, != Vector
    //
    
    template <class V1, class V2>
    static bool DoEq(
        const BaseVector_Calc<V1>& v1, const BaseVector_Calc<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same)); 
        TMVAssert(v1.size() == v2.size());
        const int size = Sizes<V1::_size,V2::_size>::size;
        const int n = (size == UNKNOWN ? int(v1.size()) : size);
        for(int i=0;i<n;++i) if (v1.cref(i) != v2.cref(i)) return false;
        return true;
    }

    template <class V1, class V2>
    static TMV_INLINE_ND bool operator==(
        const BaseVector<V1>& v1, const BaseVector<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same)); 
        TMVAssert(v1.size() == v2.size());
        return DoEq(v1.calc(),v2.calc());
    }

    template <class V1, class V2>
    static TMV_INLINE bool operator!=(
        const BaseVector<V1>& v1, const BaseVector<V2>& v2)
    { return !(v1 == v2); }



    //
    // Other Functions of Vectors
    // (These just call the method version.)
    //

    template <class V>
    static TMV_INLINE typename V::value_type SumElements(const BaseVector<V>& v)
    { return v.sumElements(); }
    template <class V>
    static TMV_INLINE typename V::float_type SumAbsElements(
        const BaseVector<V>& v)
    { return v.sumAbsElements(); }
    template <class V>
    static TMV_INLINE typename V::real_type SumAbs2Elements(
        const BaseVector<V>& v)
    { return v.sumAbs2Elements(); }

    template <class V>
    static TMV_INLINE typename V::float_type Norm(const BaseVector<V>& v)
    { return v.norm(); }
    template <class V>
    static TMV_INLINE typename V::float_type Norm1(const BaseVector<V>& v)
    { return v.norm1(); }
    template <class V>
    static TMV_INLINE typename V::float_type Norm2(const BaseVector<V>& v)
    { return v.norm2(); }
    template <class V>
    static TMV_INLINE typename V::real_type NormSq(const BaseVector<V>& v)
    { return v.normSq(); }
    template <class V>
    static TMV_INLINE typename V::float_type NormInf(const BaseVector<V>& v)
    { return v.normInf(); }

    template <class V>
    static TMV_INLINE typename V::value_type MaxElement(const BaseVector<V>& v)
    { return v.maxElement(); }
    template <class V>
    static TMV_INLINE typename V::float_type MaxAbsElement(
        const BaseVector<V>& v)
    { return v.maxAbsElement(); }
    template <class V>
    static TMV_INLINE typename V::real_type MaxAbs2Element(
        const BaseVector<V>& v)
    { return v.maxAbs2Element(); }
    template <class V>
    static TMV_INLINE typename V::value_type MinElement(const BaseVector<V>& v)
    { return v.minElement(); }
    template <class V>
    static TMV_INLINE typename V::float_type MinAbsElement(
        const BaseVector<V>& v)
    { return v.minAbsElement(); }
    template <class V>
    static TMV_INLINE typename V::real_type MinAbs2Element(
        const BaseVector<V>& v)
    { return v.minAbs2Element(); }

    template <class V>
    static TMV_INLINE typename V::const_conjugate_type Conjugate(
        const BaseVector_Calc<V>& v)
    { return v.conjugate(); }

    //
    // TMV_Text 
    //

#ifdef TMV_TEXT
    static inline std::string TMV_Text(ADType ad)
    { return ad == Ascend ? "Ascend" : "Descend"; }

    static inline std::string TMV_Text(CompType comp)
    {
        return comp == RealComp ? "Real" : comp == AbsComp ? "Abs" :
            comp == Abs2Comp ? "Abs2" : comp == ImagComp ? "Imag" : 
            comp == ArgComp ? "Arg" : comp == NormComp ? "Norm" :
            comp == ValueComp ? "Value" : comp == InverseComp ? "Inverse" :
            "UnknownComp";
    }

    template <class V>
    static inline std::string TMV_Text(const BaseVector<V>& v)
    {
        std::ostringstream s;
        s << "BaseVector< "<<TMV_Text(v.vec())<<" >";
        return s.str();
    }

    template <class V>
    static inline std::string TMV_Text(const BaseVector_Calc<V>& v)
    {
        std::ostringstream s;
        s << "BaseVector_Calc< "<<TMV_Text(v.vec())<<" >";
        return s.str();
    }

    template <class V>
    static inline std::string TMV_Text(const BaseVector_Mutable<V>& v)
    {
        std::ostringstream s;
        s << "BaseVector_Mutable< "<<TMV_Text(v.vec())<<" >";
        return s.str();
    }
#endif

} // namespace tmv

#endif
