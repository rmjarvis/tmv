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

#ifndef TMV_CopyV_H
#define TMV_CopyV_H

#include "TMV_BaseVector.h"
#include <string.h> // for memmove
#include "TMV_ConjugateV.h"

namespace tmv {

    // Defined below:
    template <class V1, class V2>
    inline void Copy(
        const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2);
    template <class V1, class V2>
    inline void NoAliasCopy(
        const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2);
    template <class V1, class V2>
    inline void InlineCopy(
        const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2);
    template <class V1, class V2>
    inline void AliasCopy(
        const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2);

    // Defined in TMV_Vector.cpp
    template <class T1, bool C1, class T2>
    void InstCopy(
        const ConstVectorView<T1,UNKNOWN,C1>& v1, VectorView<T2> v2); 

    //
    // Copy Vectors
    //

    template <int algo, int size, class V1, class V2> 
    struct CopyV_Helper;

    // algo 0: size = 0, nothing to do
    template <class V1, class V2>
    struct CopyV_Helper<0,0,V1,V2>
    {
        typedef typename V1::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static inline void call(const V1&, V2& ) {}
        static inline void call2(int , IT1 , IT2 ) {}
    };

    // algo 11: simple for loop
    template <int size, class V1, class V2>
    struct CopyV_Helper<11,size,V1,V2>
    {
        typedef typename V1::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static inline void call(const V1& v1, V2& v2)
        {
            const int n = size == UNKNOWN ? int(v1.size()) : size;
            for(int i=0;i<n;++i) v2.ref(i) = v1.cref(i); 
        }
        static inline void call2(int n, IT1 it1, IT2 it2)
        { for(;n;--n) *it2++ = *it1++; }
    };

    // algo 15: fully unroll
    template <int size, class V1, class V2>
    struct CopyV_Helper<15,size,V1,V2>
    {
        template <int I, int N>
        struct Unroller
        {
            static inline void dounroll(const V1& v1, V2& v2)
            {
                Unroller<I,N/2>::dounroll(v1,v2);
                Unroller<I+N/2,N-N/2>::dounroll(v1,v2);
            }
        };
        template <int I>
        struct Unroller<I,1>
        {
            static inline void dounroll(const V1& v1, V2& v2)
            { v2.ref(I) = v1.cref(I); }
        };
        template <int I>
        struct Unroller<I,0>
        { static inline void dounroll(const V1& v1, V2& v2) {} };
        static inline void call(const V1& v1, V2& v2)
        { Unroller<0,size>::dounroll(v1,v2); }
    };

    // algo 21: memmove
    template <int size, class V1, class V2>
    struct CopyV_Helper<21,size,V1,V2>
    {
        static inline void call(const V1& v1, V2& v2)
        {
            const int n = size == UNKNOWN ? int(v1.size()) : size;
            memmove(v2.ptr(),v1.cptr(),n*sizeof(typename V2::value_type));
        }
        static inline void call2(
            const int n, 
            typename V1::const_iterator it1, typename V2::iterator it2)
        { memmove(it2.getP(),it1.getP(),n*sizeof(typename V2::value_type)); }
    };

    // algo 22: std::copy
    template <int size, class V1, class V2>
    struct CopyV_Helper<22,size,V1,V2> 
    {
        static inline void call(const V1& v1, V2& v2)
        { std::copy(v1.begin(),v1.end(),v2.begin()); }
    };


    // algo -4: No branches or copies
    template <int size, class V1, class V2>
    struct CopyV_Helper<-4,size,V1,V2>
    {
        typedef typename V1::const_iterator IT1;
        typedef typename V2::iterator IT2;

        static inline void call(const V1& v1, V2& v2)
        {
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            const int algo = 
                size == 0 ? 0 :
#if TMV_OPT > 0
                ( size != UNKNOWN && size <= 8 ) ? 15 :
                ( Traits2<T1,T2>::sametype && 
                  V1::vconj == int(V2::vconj) &&
                  V1::vstep == 1 && V2::vstep == 1 ) ? 21 :
#endif
                11;
            CopyV_Helper<algo,size,V1,V2>::call(v1,v2);
        }
        static inline void call2(const int n, IT1 it1, IT2 it2)
        {
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            const int algo = 
                size == 0 ? 0 :
#if TMV_OPT > 0
                ( Traits2<T1,T2>::sametype && 
                  V1::vconj == int(V2::vconj) &&
                  V1::vstep == 1 && V2::vstep == 1 ) ? 21 :
#endif
                11;
            CopyV_Helper<algo,size,V1,V2>::call2(n,it1,it2);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int size, class V1, class V2>
    struct CopyV_Helper<-3,size,V1,V2>
    {
        static inline void call(const V1& v1, V2& v2)
        { CopyV_Helper<-4,size,V1,V2>::call(v1,v2); }
    };

    // algo 97: Conjugate
    template <int size, class V1, class V2>
    struct CopyV_Helper<97,size,V1,V2>
    {
        static inline void call(const V1& v1, V2& v2)
        { 
            typedef typename V1::const_conjugate_type V1c;
            typedef typename V2::conjugate_type V2c;
            V1c v1c = v1.conjugate();
            V2c v2c = v2.conjugate();
            CopyV_Helper<-2,size,V1c,V2c>::call(v1c,v2c);
        }
    };

    // algo 98: Call inst
    template <int size, class V1, class V2>
    struct CopyV_Helper<98,size,V1,V2>
    {
        static inline void call(const V1& v1, V2& v2)
        { InstCopy(v1.xView(),v2.xView()); }
    };

    // algo -2: Check for inst
    template <int size, class V1, class V2>
    struct CopyV_Helper<-2,size,V1,V2>
    {
        static inline void call(const V1& v1, V2& v2)
        {
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            const bool inst = 
                V1::vsize == UNKNOWN &&
                V2::vsize == UNKNOWN &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T2>::isinst;
            const bool conj = V2::vconj;
            const int algo = 
                inst ? (conj ? 97 : 98) :
                -4;
            CopyV_Helper<algo,size,V1,V2>::call(v1,v2);
        }
    };

    // algo 99: Check for aliases
    template <int size, class V1, class V2>
    struct CopyV_Helper<99,size,V1,V2>
    {
        static inline void call(const V1& v1, V2& v2)
        {
            if ( !SameStorage(v1,v2) ||
                 v1.step()*v2.step() < 0 || 
                 std::abs(v2.step()) < std::abs(v1.step())) { 
                // No aliasing (or no clobbering)
                CopyV_Helper<-2,size,V1,V2>::call(v1,v2);
            } else if (ExactSameStorage(v1,v2)) {
                // They are already equal modulo a conjugation.
                Maybe<V1::vconj != int(V2::vconj)>::conjself(v2); 
            } else {
                // Need a temporary
                NoAliasCopy(v1.copy(),v2);
            }
        }
    };

    // algo -1: Check for aliases?
    template <int size, class V1, class V2>
    struct CopyV_Helper<-1,size,V1,V2>
    {
        static inline void call(const V1& v1, V2& v2)
        {
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            const bool samestep = VStepHelper<V1,V2>::same;
            const bool noclobber = VStepHelper<V1,V2>::noclobber;
            const bool checkalias =
                V1::vsize == UNKNOWN &&
                V2::vsize == UNKNOWN &&
                (samestep || !noclobber);
            const int algo = 
                checkalias ? 99 :
                -2;
            CopyV_Helper<algo,size,V1,V2>::call(v1,v2);
        }
    };

    template <class V1, class V2> 
    inline void Copy(
        const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same)); 
        TMVAssert(v1.size() == v2.size());
        const int size = Sizes<V1::vsize,V2::vsize>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::cview_type V2v;
        V1v v1v = v1.cView();
        V2v v2v = v2.cView();
        typedef typename V1::value_type T1;
        typedef typename V2::value_type T2;
        CopyV_Helper<-1,size,V1v,V2v>::call(v1v,v2v); 
    }

    template <class V1, class V2>
    inline void NoAliasCopy(
        const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same)); 
        TMVAssert(v1.size() == v2.size());
        const int size = Sizes<V1::vsize,V2::vsize>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::cview_type V2v;
        V1v v1v = v1.cView();
        V2v v2v = v2.cView();
        CopyV_Helper<-2,size,V1v,V2v>::call(v1v,v2v); 
    }

    template <class V1, class V2>
    inline void InlineCopy(
        const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same)); 
        TMVAssert(v1.size() == v2.size());
        const int size = Sizes<V1::vsize,V2::vsize>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::cview_type V2v;
        V1v v1v = v1.cView();
        V2v v2v = v2.cView();
        CopyV_Helper<-4,size,V1v,V2v>::call(v1v,v2v); 
    }

    template <class V1, class V2>
    inline void AliasCopy(
        const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same)); 
        TMVAssert(v1.size() == v2.size());
        const int size = Sizes<V1::vsize,V2::vsize>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::cview_type V2v;
        V1v v1v = v1.cView();
        V2v v2v = v2.cView();
        CopyV_Helper<99,size,V1v,V2v>::call(v1v,v2v); 
    }


} // namespace tmv

#endif
