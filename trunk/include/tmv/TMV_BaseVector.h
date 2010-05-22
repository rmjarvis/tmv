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
// const V& vec() const { return *static_cast<const V*>(this); }
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

#include <sstream>
#include "TMV_Base.h"
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
    template <class T, IndexStyle I=CStyle>
    class Vector;
    template <class T, int S=UNKNOWN, bool C=false, IndexStyle I=CStyle>
    class ConstVectorView;
    template <class T, int S=UNKNOWN, bool C=false, IndexStyle I=CStyle>
    class VectorView;
    template <class T, int N, IndexStyle I=CStyle>
    class SmallVector;
    template <class T, int N, int S=1, bool C=false, IndexStyle I=CStyle>
    class ConstSmallVectorView;
    template <class T, int N, int S=1, bool C=false, IndexStyle I=CStyle>
    class SmallVectorView;

    // These are effectively aliases for I = FortranStyle
    // so you don't have to write all the S,C values if you are using
    // the defaults when you want to specify FortranStyle.
    template <class T>
    class VectorF;
    template <class T, int S=UNKNOWN, bool C=false>
    class ConstVectorViewF;
    template <class T, int S=UNKNOWN, bool C=false>
    class VectorViewF;
    template <class T, int N>
    class SmallVectorF;
    template <class T, int N, int S=1, bool C=false>
    class ConstSmallVectorViewF;
    template <class T, int N, int S=1, bool C=false>
    class SmallVectorViewF;

    // Returned by sort(p)
    class Permutation;


    //
    // Helper functions and values:
    //

    // Values used to define how to sort a vector.
    enum ADType { Ascend, Descend };
    enum CompType { 
        RealComp, AbsComp, Abs2Comp, ImagComp,
        ArgComp, NormComp, ValueComp, InverseComp, LogComp };
    enum OldADType { ASCEND, DESCEND };
    enum OldCOMPType { 
        REAL_COMP, ABS_COMP, ABS2_COMP, IMAG_COMP, ARG_COMP, NORM_COMP, VALUE };

    // A helper class to pick out a component of a (usually) complex value
    template <CompType comp, class T> 
    struct Component;
    // First the real version:
    template <class T> 
    struct Component<RealComp,T>
    { 
        // return f(x);
        static inline T f(const T& x) { return x; } 
        // apply x = f(x)
        // (usually assigning the result to the real part of x when x is complex
        static inline void applyf(T& x) { }
        // get the applied value of f(x) after doing applyf
        // usually this returns real(x), but sometimes the whole x
        static inline T get(const T& x) { return x; }
    };
    template <class T> 
    struct Component<ImagComp,T>
    { 
        static inline T f(const T& x) { return T(0); } 
        static inline void applyf(T& x) { x = T(0); }
        static inline T get(const T& x) { return x; }
    };
    template <class T> 
    struct Component<AbsComp,T>
    {
        static inline T f(const T& x) { return TMV_ABS(x); } 
        static inline void applyf(T& x) { x = TMV_ABS(x); }
        static inline T get(const T& x) { return x; }
    };
    template <class T> 
    struct Component<Abs2Comp,T>
    { 
        static inline T f(const T& x) { return TMV_ABS2(x); } 
        static inline void applyf(T& x) { x = TMV_ABS2(x); }
        static inline T get(const T& x) { return x; }
    };
    template <class T> 
    struct Component<ArgComp,T>
    { 
        static inline T f(const T& x) { return TMV_ARG(x); } 
        static inline void applyf(T& x) { x = TMV_ARG(x); }
        static inline T get(const T& x) { return x; }
    };
    template <class T> 
    struct Component<NormComp,T>
    { 
        static inline T f(const T& x) { return TMV_NORM(x); } 
        static inline void applyf(T& x) { x = TMV_NORM(x); }
        static inline T get(const T& x) { return x; }
    };
    template <class T> 
    struct Component<ValueComp,T>
    {
        static inline T f(const T& x) { return x; } 
        static inline void applyf(T& x) { }
        static inline T get(const T& x) { return x; }
    };
    template <class T> 
    struct Component<InverseComp,T>
    {
        static inline T f(const T& x) { return T(1) / x; } 
        static inline void applyf(T& x) { x = T(1) / x; }
        static inline T get(const T& x) { return x; }
    };
    template <class T> 
    struct Component<LogComp,T>
    {
        static inline T f(const T& x) { return TMV_LOG(x); } 
        static inline void applyf(T& x) { x = TMV_LOG(x); }
        static inline T get(const T& x) { return x; }
    };

    // Now the complex version:
#define CT std::complex<T>
    template <class T> 
    struct Component<RealComp,CT>
    { 
        static inline T f(const CT& x) { return real(x); } 
        static inline void applyf(CT& x) { }
        static inline T get(const T& x) { return real(x); }
    };
    template <class T> 
    struct Component<ImagComp,CT>
    { 
        static inline T f(const CT& x) { return imag(x); } 
        static inline void applyf(CT& x) { real(x) = imag(x); }
        static inline T get(const CT& x) { return real(x); }
    };
    template <class T> 
    struct Component<AbsComp,CT>
    {
        static inline T f(const CT& x) { return TMV_ABS(x); } 
        static inline void applyf(CT& x) { real(x) = TMV_ABS(x); }
        static inline T get(const CT& x) { return real(x); }
    };
    template <class T> 
    struct Component<Abs2Comp,CT>
    { 
        static inline T f(const CT& x) { return TMV_ABS2(x); } 
        static inline void applyf(CT& x) { real(x) = TMV_ABS2(x); }
        static inline T get(const CT& x) { return real(x); }
    };
    template <class T> 
    struct Component<ArgComp,CT>
    { 
        static inline T f(const CT& x) { return TMV_ARG(x); } 
        static inline void applyf(CT& x) { real(x) = TMV_ARG(x); }
        static inline T get(const CT& x) { return real(x); }
    };
    template <class T> 
    struct Component<NormComp,CT>
    { 
        static inline T f(const CT& x) { return TMV_NORM(x); } 
        static inline void applyf(CT& x) { real(x) = TMV_NORM(x); }
        static inline T get(const CT& x) { return real(x); }
    };
    template <class T> 
    struct Component<ValueComp,CT>
    {
        static inline CT f(const CT& x) { return x; } 
        static inline void applyf(CT& x) { }
        static inline CT get(const CT& x) { return x; }
    };
    template <class T> 
    struct Component<InverseComp,CT>
    {
        static inline CT f(const CT& x) { return std::conj(x) / TMV_NORM(x); }
        static inline void applyf(CT& x) { x /= TMV_NORM(x); }
        static inline CT get(const CT& x) { return std::conj(x); }
    };
    template <class T> 
    struct Component<LogComp,CT>
    {
        static inline T f(const CT& x) { return TMV_LOG(real(x)); } 
        static inline void applyf(CT& x) { real(x) = TMV_LOG(real(x)); }
        static inline T get(const CT& x) { return real(x); }
    };
#undef CT

    // These helper functions check the validity of indices according
    // to whether the vector uses CStyle or FortranStyle indexing.
    // They also update the indices to be consistent with CStyle.
    template <bool _fort>
    inline void CheckIndex(int& i, int n) 
    { TMVAssert(i>=0 && i<n && "index is not valid"); } // CStyle
    template <>
    inline void CheckIndex<true>(int& i, int n) 
    { TMVAssert(i>=1 && i<=n && "index is not valid"); --i; } // FortranStyle
    template <bool _fort>
    inline void CheckRange(int& i1, int i2, int n)
    { // CStyle
        TMVAssert(i1 >= 0 && "first element must be in range");
        TMVAssert(i2 <= n && "last element must be in range");
        TMVAssert(i2 >= i1 && 
                  "range must have a non-negative number of elements");
    }
    template <>
    inline void CheckRange<true>(int& i1, int i2, int n)
    { // FortranStyle
        TMVAssert(i1 >= 1 && "first element must be in range");
        TMVAssert(i2 <= n && "last element must be in range");
        TMVAssert(i2 >= i1 && "range must have a positive number of elements");
        --i1;
    }
    template <bool _fort>
    inline void CheckRange(int& i1, int& i2, int istep, int n)
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
    inline void CheckRange<true>(int& i1, int& i2, int istep, int n)
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
    template <class T, int _size, bool _fort>
    struct VCopyHelper
    { typedef SmallVector<T,_size,_fort?FortranStyle:CStyle> type; };
    template <class T, bool _fort>
    struct VCopyHelper<T,UNKNOWN,_fort>
    { typedef Vector<T,_fort?FortranStyle:CStyle> type; };

    // This is a type to use when there is no valid return type for a 
    // particular function. 
    // It will give a compiler error if it is ever used.
    class InvalidType { private: InvalidType(); };

    // This helper function checks for aliasis in the memory storage of
    // two objects.  We overload it for specific objects that can be 
    // checked to have the same realPart().cptr() values.
    template <class V1, class V2>
    inline bool SameStorage(
        const BaseVector<V1>& v1, const BaseVector<V2>& v2)
    { return false; }
#ifndef TMV_NO_ALIAS_CHECK
    template <class V1, class V2>
    inline bool SameStorage(
        const BaseVector_Calc<V1>& v1, const BaseVector_Calc<V2>& v2)
    { return v1.realPart().cptr() == v2.realPart().cptr(); }
#endif

    // This next one is sometime used after SameStorage is found to be true.
    // Here we check whether the elements are exactly the same 
    // by whether the steps are the same.
    // We do not check the _conj values, so that is usually the next 
    // step depending on why we are checking this.
    template <class V1, class V2>
    inline bool ExactSameStorage(
        const BaseVector<V1>& v1, const BaseVector<V2>& v2)
    { return false; }
    template <class V1, class V2>
    inline bool ExactSameStorage(
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
    inline void Copy(const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2);
    template <class V1, class V2>
    inline void NoAliasCopy(
        const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2);

    // Defined in TMV_SwapV.h
    template <class V1, class V2>
    inline void Swap(BaseVector_Mutable<V1>& v1, BaseVector_Mutable<V2>& v2);
    template <class V1, class V2>
    inline void NoAliasSwap(
        BaseVector_Mutable<V1>& v1, BaseVector_Mutable<V2>& v2);
    template <class V>
    inline void ReverseSelf(BaseVector_Mutable<V>& v);
    template <class V>
    inline void ConjugateSelf(BaseVector_Mutable<V>& v);

    // Defined in TMV_NormV.h
    template <class V>
    inline typename V::real_type NormSq(const BaseVector_Calc<V>& v);
    template <class V>
    inline typename V::real_type NormSq(const BaseVector_Calc<V>& v, 
                                        const typename V::real_type scale);
    template <class V>
    inline typename V::real_type Norm2(const BaseVector_Calc<V>& v);
    template <class V>
    inline typename V::value_type SumElements(const BaseVector_Calc<V>& v);
    template <class V>
    inline typename V::real_type SumAbsElements(const BaseVector_Calc<V>& v);
    template <class V>
    inline typename V::real_type SumAbs2Elements(const BaseVector_Calc<V>& v);

    // Defined in TMV_Det.h
    template <class V>
    inline typename V::value_type ProdElements(const BaseVector_Calc<V>& v);
    template <class V>
    inline typename V::real_type LogProdElements(
        const BaseVector_Calc<V>& v, typename V::value_type* sign);
    template <class V>
    inline bool HasZeroElement(const BaseVector_Calc<V>& v);

    // Defined in TMV_MinMax.h
    template <class V>
    inline typename V::value_type MaxElement(
        const BaseVector_Calc<V>& v, int*const imax=0);
    template <class V>
    inline typename V::real_type MaxAbsElement(
        const BaseVector_Calc<V>& v, int*const imax=0);
    template <class V>
    inline typename V::real_type MaxAbs2Element(
        const BaseVector_Calc<V>& v, int*const imax=0);
    template <class V>
    inline typename V::value_type MinElement(
        const BaseVector_Calc<V>& v, int*const imin=0);
    template <class V>
    inline typename V::real_type MinAbsElement(
        const BaseVector_Calc<V>& v, int*const imin=0);
    template <class V>
    inline typename V::real_type MinAbs2Element(
        const BaseVector_Calc<V>& v, int*const imin=0);

    // Defined in TMV_VectorIO.h
    template <class V>
    inline void Write(std::ostream& os, const BaseVector_Calc<V>& v);
    template <class V>
    inline void Write(std::ostream& os, const BaseVector_Calc<V>& v,
                      typename V::real_type thresh) ;
    template <class V>
    inline void Read(std::istream& is, BaseVector_Mutable<V>& v);

    // Defined in TMV_SortV.h
    template <class V>
    inline void Sort(BaseVector_Mutable<V>& v, ADType ad, CompType comp);
    template <class V>
    inline void Sort(
        BaseVector_Mutable<V>& v, int*const P, ADType ad, CompType comp);

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

        typedef V type;

        typedef typename Traits<V>::value_type value_type;
        typedef typename Traits<V>::calc_type calc_type;
        typedef typename Traits<V>::eval_type eval_type;
        typedef typename Traits<V>::copy_type copy_type;

        // Derived values:
        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;
        enum { isreal = Traits<value_type>::isreal };
        enum { iscomplex = Traits<value_type>::iscomplex };
        enum { unknownsizes = _size == UNKNOWN };

        //
        // Constructor
        //

        inline BaseVector() {}
        inline BaseVector(const BaseVector<V>&) {}
        inline ~BaseVector() {}

        //
        // Access
        // 

        inline value_type operator[](int i) const 
        { return operator()(i); }
        inline value_type operator()(int i) const 
        {
            CheckIndex<_fort>(i,size());
            return cref(i);
        }

        //
        // Functions
        //

        inline value_type sumElements() const
        { return tmv::SumElements(calc().cView()); }

        inline real_type sumAbsElements() const 
        { return tmv::SumAbsElements(calc().cView()); }

        inline real_type sumAbs2Elements() const
        { return tmv::SumAbs2Elements(calc().cView()); }

        inline value_type maxElement(int*const imax=0) const
        { return tmv::MaxElement(calc().cView(),imax); }

        inline real_type maxAbsElement(int*const imax=0) const 
        { return tmv::MaxAbsElement(calc().cView(),imax); }

        inline real_type maxAbs2Element(int*const imax=0) const
        { return tmv::MaxAbs2Element(calc().cView(),imax); }

        inline value_type minElement(int*const imin=0) const
        { return tmv::MinElement(calc().cView(),imin); }

        inline real_type minAbsElement(int*const imin=0) const 
        { return tmv::MinAbsElement(calc().cView(),imin); }

        inline real_type minAbs2Element(int*const imin=0) const
        { return tmv::MinAbs2Element(calc().cView(),imin); }

        inline real_type norm1() const
        { return sumAbsElements(); }

        inline real_type normSq() const
        { return tmv::NormSq(calc().cView()); }

        inline real_type normSq(const real_type scale) const
        { return tmv::NormSq(calc().cView(),scale); }

        inline real_type norm2() const
        { return tmv::Norm2(calc().cView()); }

        inline real_type norm() const
        { return norm2(); }

        inline real_type normInf() const
        { return size() > 0 ? maxAbsElement() : real_type(0); }

        inline value_type prodElements() const
        { return tmv::ProdElements(calc()); }

        inline real_type logProdElements(value_type* sign=0) const
        { return tmv::LogProdElements(calc(),sign); }

        inline bool hasZeroElement() const
        { return tmv::HasZeroElement(calc()); }



        // 
        // I/O
        //

        inline void write(std::ostream& os) const
        { tmv::Write(os,calc().cView()); }
        inline void write(std::ostream& os, real_type thresh) const
        { tmv::Write(os,calc().cView(),thresh); }



        //
        // Auxilliary routines
        //

        inline const type& vec() const 
        { return *static_cast<const type*>(this); }

        inline calc_type calc() const 
        { return static_cast<calc_type>(vec()); }

        inline eval_type eval() const 
        { return static_cast<eval_type>(vec()); }

        inline copy_type copy() const 
        { return static_cast<copy_type>(vec()); }

        // Note that these last functions need to be defined in a more derived
        // class than this, or an infinite loop will result when compiling.

        inline size_t size() const { return vec().size(); }

        inline value_type cref(int i) const  { return vec().cref(i); }

        template <class V2>
        inline void assignTo(BaseVector_Mutable<V2>& v2) const
        { vec().assignTo(v2); }

        template <class V2>
        inline void newAssignTo(BaseVector_Mutable<V2>& v2) const
        { vec().newAssignTo(v2); }

    private :
        inline void operator=(const BaseVector<V>& v2);

    }; // BaseVector

    //
    // BaseVector_Calc
    //

    template <class V> 
    class BaseVector_Calc : public BaseVector<V>
    {
    public:
        enum { _size = Traits<V>::_size };
        enum { _fort = Traits<V>::_fort };
        enum { _calc = Traits<V>::_calc };
        enum { _step = Traits<V>::_step };
        enum { _conj = Traits<V>::_conj };

        typedef V type;

        typedef typename Traits<V>::value_type value_type;

        typedef typename Traits<V>::calc_type calc_type;
        typedef typename Traits<V>::eval_type eval_type;
        typedef typename Traits<V>::copy_type copy_type;

        typedef typename Traits<V>::const_subvector_type const_subvector_type;
        typedef typename Traits<V>::const_subvector_step_type 
            const_subvector_step_type;
        typedef typename Traits<V>::const_view_type const_view_type;
        typedef typename Traits<V>::const_cview_type const_cview_type;
        typedef typename Traits<V>::const_fview_type const_fview_type;
        typedef typename Traits<V>::const_xview_type const_xview_type;
        typedef typename Traits<V>::const_unitview_type const_unitview_type;
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

        // Derived values:
        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;

        //
        // Constructor
        //

        inline BaseVector_Calc() {}
        inline BaseVector_Calc(const BaseVector_Calc<V>&) {}
        inline ~BaseVector_Calc() {}


        //
        // Access 
        //

        inline const_iterator begin() const
        { return const_iterator(cptr(),step()); }
        inline const_iterator end() const
        { return begin() + size(); }
        inline const_reverse_iterator rbegin() const
        { return const_reverse_iterator(cptr(),step()); }
        inline const_reverse_iterator rend() const
        { return rbegin() + size(); }


        //
        // subVector
        //

        // cSubVector always uses CStyle
        inline const_subvector_type cSubVector(int i1, int i2) const
        { return const_subvector_type(cptr()+i1*step(),i2-i1,step()); }

        inline const_subvector_step_type cSubVector(
            int i1, int i2, int istep) const
        {
            return const_subvector_step_type(
                cptr()+i1*step(), (i2-i1)/istep, istep*step());
        }

        inline const_subvector_type subVector(int i1, int i2) const
        {
            CheckRange<_fort>(i1,i2,size());
            return cSubVector(i1,i2);
        }

        inline const_subvector_step_type subVector(
            int i1, int i2, int istep) const
        {
            CheckRange<_fort>(i1,i2,istep,size());
            return cSubVector(i1,i2,istep);
        }


        //
        // Views
        //

        inline const_view_type view() const
        { return const_view_type(cptr(),size(),step()); }

        inline const_cview_type cView() const
        { return view(); }

        inline const_fview_type fView() const
        { return view(); }

        inline const_xview_type xView() const
        { return view(); }

        inline const_unitview_type unitView() const
        {
            TMVAssert(step() == 1 && "Called unitView on vector with step!=1");
            return view(); 
        }

        inline const_view_type constView() const
        { return view(); }

        inline const_conjugate_type conjugate() const
        { return const_conjugate_type(cptr(),size(),step()); }

        inline const_reverse_type reverse() const
        {
            return const_reverse_type(
                cptr()+(int(size())-1)*step(), size(), -step());
        }

        inline const_realpart_type realPart() const
        {
            return const_realpart_type(
                reinterpret_cast<const real_type*>(cptr()),
                size(), V::isreal ? step() : 2*step());
        }

        inline const_imagpart_type imagPart() const
        {
            TMVStaticAssert(V::iscomplex);
            return const_imagpart_type(
                reinterpret_cast<const real_type*>(cptr())+1, size(), 2*step());
        }

        inline const_flatten_type flatten() const
        {
            TMVStaticAssert(V::iscomplex);
            TMVStaticAssert(_step == UNKNOWN || _step == 1);
            TMVAssert(step() == 1);
            return const_flatten_type(
                reinterpret_cast<const real_type*>(cptr()),
                V::isreal ? size() : 2*size(), 1);
        }

        inline const_nonconj_type nonConj() const
        { return const_nonconj_type(cptr(),size(),step()); }

        inline nonconst_type nonConst() const
        { 
            return nonconst_type(
                const_cast<value_type*>(cptr()),size(),step());
        }

        //
        // Auxilliary routines
        //

        template <class V2>
        inline void assignTo(BaseVector_Mutable<V2>& v2) const
        { tmv::Copy(*this,v2); }

        template <class V2>
        inline void newAssignTo(BaseVector_Mutable<V2>& v2) const
        { tmv::NoAliasCopy(*this,v2); }

        inline const type& vec() const
        { return *static_cast<const type*>(this); }

        inline bool isconj() const { return _conj; }

        // Note that these last functions need to be defined in a more derived
        // class than this, or an infinite loop will result when compiling.

        inline size_t size() const { return vec().size(); }
        inline int step() const { return vec().step(); }
        inline const value_type* cptr() const { return vec().cptr(); }
        inline value_type cref(int i) const  { return vec().cref(i); }

    private:
        void operator=(const BaseVector_Calc<V>&);

    }; // BaseVector_Calc

    //
    // BaseVector_Mutable
    //

    template <class V> 
    class BaseVector_Mutable : public BaseVector_Calc<V>
    {
    public:
        enum { _size = Traits<V>::_size };
        enum { _fort = Traits<V>::_fort };
        enum { _calc = Traits<V>::_calc };
        enum { _step = Traits<V>::_step };
        enum { _conj = Traits<V>::_conj };

        typedef V type;
        typedef BaseVector_Calc<V> base_inst;
        typedef BaseVector_Mutable<V> base_mut;

        typedef typename Traits<V>::value_type value_type;

        typedef typename Traits<V>::calc_type calc_type;
        typedef typename Traits<V>::eval_type eval_type;
        typedef typename Traits<V>::copy_type copy_type;

        typedef typename Traits<V>::const_subvector_type const_subvector_type;
        typedef typename Traits<V>::const_subvector_step_type 
            const_subvector_step_type;
        typedef typename Traits<V>::const_view_type const_view_type;
        typedef typename Traits<V>::const_cview_type const_cview_type;
        typedef typename Traits<V>::const_fview_type const_fview_type;
        typedef typename Traits<V>::const_xview_type const_xview_type;
        typedef typename Traits<V>::const_unitview_type const_unitview_type;
        typedef typename Traits<V>::const_conjugate_type const_conjugate_type;
        typedef typename Traits<V>::const_reverse_type const_reverse_type;
        typedef typename Traits<V>::const_realpart_type const_realpart_type;
        typedef typename Traits<V>::const_imagpart_type const_imagpart_type;
        typedef typename Traits<V>::const_flatten_type const_flatten_type;
        typedef typename Traits<V>::const_nonconj_type const_nonconj_type;

        typedef typename Traits<V>::const_iterator const_iterator;
        typedef typename Traits<V>::const_reverse_iterator 
            const_reverse_iterator;

        typedef typename Traits<V>::subvector_type subvector_type;
        typedef typename Traits<V>::subvector_step_type subvector_step_type;
        typedef typename Traits<V>::view_type view_type;
        typedef typename Traits<V>::cview_type cview_type;
        typedef typename Traits<V>::fview_type fview_type;
        typedef typename Traits<V>::xview_type xview_type;
        typedef typename Traits<V>::unitview_type unitview_type;
        typedef typename Traits<V>::conjugate_type conjugate_type;
        typedef typename Traits<V>::reverse_type reverse_type;
        typedef typename Traits<V>::realpart_type realpart_type;
        typedef typename Traits<V>::imagpart_type imagpart_type;
        typedef typename Traits<V>::flatten_type flatten_type;
        typedef typename Traits<V>::nonconj_type nonconj_type;

        typedef typename Traits<V>::iterator iterator;
        typedef typename Traits<V>::reverse_iterator reverse_iterator;
        typedef typename Traits<V>::reference reference;

        // Derived values:
        typedef typename Traits<value_type>::real_type real_type;
        typedef typename Traits<value_type>::complex_type complex_type;

        //
        // Constructor
        //

        inline BaseVector_Mutable() {}
        inline BaseVector_Mutable(const base_mut&) {}
        inline ~BaseVector_Mutable() {}


        //
        // Access 
        //

        inline reference operator[](int i)
        { return operator()(i); }
        inline reference operator()(int i)
        {
            CheckIndex<_fort>(i,size());
            return ref(i);
        }

        inline iterator begin()
        { return iterator(ptr(),step()); }
        inline iterator end()
        { return begin() + size(); }
        inline reverse_iterator rbegin()
        { return reverse_iterator(ptr(),step()); }
        inline reverse_iterator rend()
        { return rbegin() + size(); }

        // We need to repeat the const versions so the non-const ones
        // don't clobber them.
        inline value_type operator[](int i) const
        { return base_inst::operator[](i); }
        inline value_type operator()(int i) const
        { return base_inst::operator()(i); }

        inline const_iterator begin() const
        { return base_inst::begin(); }
        inline const_iterator end() const
        { return base_inst::end(); }
        inline const_reverse_iterator rbegin() const
        { return base_inst::rbegin(); }
        inline const_reverse_iterator rend() const
        { return base_inst::rend(); }


        //
        // Op =
        //

        inline base_mut& operator=(base_mut& v2) 
        {
            TMVAssert(size() == v2.size());
            v2.assignTo(*this);
            return *this; 
        }

        template <class V2>
        inline base_mut& operator=(const BaseVector<V2>& v2) 
        {
            TMVStaticAssert((Sizes<_size,V2::_size>::same));
            TMVAssert(size() == v2.size());
            v2.assignTo(*this);
            return *this; 
        }

        inline ListAssigner<value_type,iterator> operator<<(value_type x)
        { return ListAssigner<value_type,iterator>(begin(),size(),x); }


        //
        // Modifying Functions
        //

        inline type& setZero() 
        {
            const int n=size();
            for(int i=0;i<n;++i) ref(i) = value_type(0);
            return vec();
        }

        inline type& clip(real_type thresh) 
        {
            const int n=size();
            for(int i=0;i<n;++i) {
                const real_type temp = TMV_ABS(cref(i));
                if (temp < thresh) ref(i) = value_type(0);
            }
            return vec();
        }

        inline type& setAllTo(value_type x) 
        {
            const int n=size();
            for(int i=0;i<n;++i) ref(i) = x;
            return vec();
        }

        inline type& addToAll(value_type x) 
        {
            const int n=size();
            for(int i=0;i<n;++i) ref(i) += x;
            return vec();
        }

        inline type& conjugateSelf()
        { tmv::ConjugateSelf(vec()); return vec(); }

        template <class F>
        inline type& applyToAll(const F& f)
        {
            const int n=size();
            for(int i=0;i<n;++i) ref(i) = f(cref(i));
            return vec();
        }

        inline base_mut& cMakeBasis(int i, value_type x=value_type(1)) 
        {
            setZero(); ref(i) = x;
            return *this;
        }
        inline base_mut& makeBasis(int i, value_type x=value_type(1)) 
        {
            CheckIndex<_fort>(i,size());
            return cMakeBasis(i,x);
        }

        inline base_mut& cSwap(int i1, int i2) 
        {
            TMV_SWAP(ref(i1),ref(i2)); 
            return *this;
        }
        inline base_mut& swap(int i1, int i2) 
        {
            CheckIndex<_fort>(i1,size());
            CheckIndex<_fort>(i2,size());
            return cSwap(i1,i2);
        }

        inline base_mut& cPermute(const int*const p, int i1, int i2) 
        {
            for(int i=i1;i<i2;++i) cSwap(i,p[i]); 
            return *this;
        }
        inline base_mut& permute(const int*const p, int i1, int i2) 
        {
            CheckRange<_fort>(i1,i2,size());
            return cPermute(p,i1,i2);
        }
        inline base_mut& permute(const int*const p) 
        { return cPermute(p,0,size()); }

        inline base_mut& cReversePermute(const int*const p, int i1, int i2) 
        {
            for(int i=i2;i>i1;) { --i; cSwap(i,p[i]); }
            return *this;
        }
        inline base_mut& reversePermute(const int*const p, int i1, int i2) 
        {
            CheckRange<_fort>(i1,i2,size());
            return cReversePermute(p,i1,i2);
        }
        inline base_mut& reversePermute(const int*const p) 
        { return cReversePermute(p,0,size()); }

        inline base_mut& reverseSelf() 
        { tmv::ReverseSelf(*this); return *this; }

        inline void sort(ADType ad=Ascend, CompType comp=RealComp) 
        { tmv::Sort(*this,ad,comp); }

        // Defined in TMV_Permutation.h
        inline Permutation sort(
            int*const P, ADType ad=Ascend, CompType comp=RealComp);

        //
        // SubVector
        //

        inline subvector_type cSubVector(int i1, int i2) 
        { return subvector_type(ptr()+i1*step(),i2-i1,step()); }

        inline subvector_step_type cSubVector(int i1, int i2, int istep) 
        {
            return subvector_step_type(
                ptr()+i1*step(), (i2-i1)/istep, istep*step());
        }

        inline subvector_type subVector(int i1, int i2) 
        {
            CheckRange<_fort>(i1,i2,size());
            return cSubVector(i1,i2);
        }

        inline subvector_step_type subVector(int i1, int i2, int istep) 
        {
            CheckRange<_fort>(i1,i2,istep,size());
            return cSubVector(i1,i2,istep);
        }

        inline reverse_type reverse() 
        { return reverse_type(ptr()+(int(size())-1)*step(), size(),-step()); }

        inline view_type view() 
        { return view_type(ptr(),size(),step()); }

        inline cview_type cView() 
        { return view(); }

        inline fview_type fView() 
        { return view(); }

        inline xview_type xView() 
        { return view(); }

        inline unitview_type unitView() 
        {
            TMVAssert(step() == 1 && "Called unitView on vector with step!=1");
            return view(); 
        }

        inline conjugate_type conjugate() 
        { return conjugate_type(ptr(),size(),step()); }


        inline realpart_type realPart() 
        {
            return realpart_type(
                reinterpret_cast<real_type*>(ptr()),
                size(), V::isreal ? step() : 2*step());
        }

        inline imagpart_type imagPart() 
        {
            TMVStaticAssert(V::iscomplex);
            return imagpart_type(
                reinterpret_cast<real_type*>(ptr())+1, size(), 2*step());
        }

        inline flatten_type flatten() 
        {
            TMVStaticAssert(V::iscomplex);
            TMVStaticAssert(_step == UNKNOWN || _step == 1);
            TMVAssert(step() == 1);
            return flatten_type(
                reinterpret_cast<real_type*>(ptr()),
                V::isreal ? size() : 2*size(), 1);
        }

        inline nonconj_type nonConj()
        { return nonconj_type(ptr(),size(),step()); }



        // Repeat the const versions:
        inline const_subvector_type cSubVector(int i1, int i2) const
        { return base_inst::cSubVector(i1,i2); }
        inline const_subvector_step_type cSubVector(
            int i1, int i2, int istep) const
        { return base_inst::cSubVector(i1,i2,istep); }
        inline const_subvector_type subVector(int i1, int i2) const
        { return base_inst::subVector(i1,i2); }
        inline const_subvector_step_type subVector(
            int i1, int i2, int istep) const
        { return base_inst::subVector(i1,i2,istep); }
        inline const_reverse_type reverse() const
        { return base_inst::reverse(); }
        inline const_view_type view() const
        { return base_inst::view(); }
        inline const_cview_type cView() const
        { return base_inst::cView(); }
        inline const_fview_type fView() const
        { return base_inst::fView(); }
        inline const_xview_type xView() const
        { return base_inst::xView(); }
        inline const_unitview_type unitView() const
        { return base_inst::unitView(); }
        inline const_conjugate_type conjugate() const
        { return base_inst::conjugate(); }
        inline const_realpart_type realPart() const
        { return base_inst::realPart(); }
        inline const_imagpart_type imagPart() const
        { return base_inst::imagPart(); }
        inline const_flatten_type flatten() const
        { return base_inst::flatten(); }
        inline const_nonconj_type nonConj() const
        { return base_inst::nonConj(); }



        //
        // I/O
        //

        // Defined in TMV_VectorIO.h
        inline void read(std::istream& is)
        {
            cview_type vcv = cView();
            tmv::Read(is,vcv); 
        }


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
        inline base_mut& operator+=(const X2& x2)
        { AddEq(vec(),x2); return *this; }

        template <class X2>
        inline base_mut& operator-=(const X2& x2)
        { SubtractEq(vec(),x2); return *this; }

        template <class X2>
        inline base_mut& operator*=(const X2& x2)
        { MultEq(vec(),x2); return *this; }

        template <class X2>
        inline base_mut& operator/=(const X2& x2)
        { DivEq(vec(),x2); return *this; }

        template <class X2>
        inline base_mut& operator%=(const X2& x2)
        { RDivEq(vec(),x2); return *this; }



        //
        // Auxilliary routines
        //

        inline const type& vec() const
        { return *static_cast<const type*>(this); }
        inline type& vec()
        { return *static_cast<type*>(this); }

        // Note that these last functionsneed to be defined in a more derived
        // class than this, or an infinite loop will result when compiling.

        inline size_t size() const { return vec().size(); }
        inline int step() const { return vec().step(); }
        inline value_type* ptr() { return vec().ptr(); }
        inline const value_type* cptr() { return vec().cptr(); }
        inline reference ref(int i) { return vec().ref(i); }
        inline value_type cref(int i) { return vec().cref(i); }

    }; // BaseVector_Mutable


    //
    // Vector ==, != Vector
    //

    template <class V1, class V2>
    static bool operator==(const BaseVector<V1>& v1, const BaseVector<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same)); 
        TMVAssert(v1.size() == v2.size());
        const int size = Sizes<V1::_size,V2::_size>::size;
        const int n = (size == UNKNOWN ? int(v1.size()) : size);
        for(int i=0;i<n;++i) if (v1.cref(i) != v2.cref(i)) return false;
        return true;
    }

    template <class V1, class V2>
    inline bool operator!=(
        const BaseVector<V1>& v1, const BaseVector<V2>& v2)
    { return !(v1 == v2); }



    //
    // Other Functions of Vectors
    // (These just call the method version.)
    //

    template <class V>
    inline typename V::real_type Norm(const BaseVector<V>& v)
    { return v.norm(); }

    template <class V>
    inline typename V::real_type Norm1(const BaseVector<V>& v)
    { return v.norm1(); }

    template <class V>
    inline typename V::real_type Norm2(const BaseVector<V>& v)
    { return v.norm2(); }

    template <class V>
    inline typename V::real_type NormSq(const BaseVector<V>& v)
    { return v.normSq(); }

    template <class V>
    inline typename V::real_type NormInf(const BaseVector<V>& v)
    { return v.normInf(); }

    template <class V>
    inline typename V::const_conjugate_type Conjugate(
        const BaseVector_Calc<V>& v)
    { return v.conjugate(); }

    template <class V>
    inline typename V::conjugate_type Conjugate(BaseVector_Mutable<V>& v)
    { return v.conjugate(); }

    //
    // TMV_Text 
    //

    inline std::string TMV_Text(ADType ad)
    { return ad == Ascend ? "Ascend" : "Descend"; }

    inline std::string TMV_Text(CompType comp)
    {
        return comp == RealComp ? "Real" : comp == AbsComp ? "Abs" :
            comp == Abs2Comp ? "Abs2" : comp == ImagComp ? "Imag" : 
            comp == ArgComp ? "Arg" : comp == NormComp ? "Norm" :
            comp == ValueComp ? "Value" : comp == InverseComp ? "Inverse" :
            "UnknownComp";
    }

    inline std::string TMV_Text(OldADType ad)
    { return ad == ASCEND ? "Ascend" : "Descend"; }

    inline std::string TMV_Text(OldCOMPType comp)
    {
        return comp == REAL_COMP ? "Real" : comp == ABS_COMP ? "Abs" :
            comp == ABS2_COMP ? "Abs2" : comp == IMAG_COMP ? "Imag" : "Arg";
    }

    template <class V>
    static std::string TMV_Text(const BaseVector<V>& v)
    {
        std::ostringstream s;
        s << "BaseVector< "<<TMV_Text(v.vec())<<" >";
        return s.str();
    }

    template <class V>
    static std::string TMV_Text(const BaseVector_Calc<V>& v)
    {
        std::ostringstream s;
        s << "BaseVector_Calc< "<<TMV_Text(v.vec())<<" >";
        return s.str();
    }

    template <class V>
    static std::string TMV_Text(const BaseVector_Mutable<V>& v)
    {
        std::ostringstream s;
        s << "BaseVector_Mutable< "<<TMV_Text(v.vec())<<" >";
        return s.str();
    }

} // namespace tmv

#endif
