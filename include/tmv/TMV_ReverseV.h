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

#ifndef TMV_ReverseV_H
#define TMV_ReverseV_H

#include "TMV_BaseVector.h"
#include "TMV_SwapV.h"

namespace tmv {

    // Defined in TMV_Vector.cpp
    template <class T>
    void InstReverseSelf(VectorView<T> v);

    //
    // ReverseSelf
    //

    template <int algo, int size, class V>
    struct ReverseV_Helper;

    // algo 11: simple for loop
    template <int size, class V>
    struct ReverseV_Helper<11,size,V>
    {
        static void call(V& v)
        {
            const int n = size == UNKNOWN ? int(v.size()) : size;
            const int no2 = n/2;
            if (no2) 
                for(int i1=0;i1<no2;++i1) v.cSwap(i1,n-i1-1);
        }
    };

    // algo 12: call swap on two halves
    template <int size, class V>
    struct ReverseV_Helper<12,size,V>
    {
        static void call(V& v)
        {
            typedef typename V::value_type T;
            typedef typename V::subvector_type V1;
            typedef typename V::subvector_type::reverse_type V2;
            const int n = size == UNKNOWN ? int(v.size()) : size;
            const int sizeo2 = size == UNKNOWN ? UNKNOWN : size/2;
            if (n > 1)  {
                V1 v1 = v.cSubVector(0,n/2);
                V2 v2 = v.cSubVector(n-n/2,n).reverse();
                SwapV_Helper<-3,sizeo2,V1,V2>::call(v1,v2);
            }
        }
    };

    // algo 15: fully unroll
    template <int size, class V>
    struct ReverseV_Helper<15,size,V>
    {
        template <int I, int N>
        struct Unroller
        {
            static void unroll(V& v)
            {
                Unroller<I,N-1>::unroll(v);
                v.cSwap(N-1,size-N);
            }
        };
        template <int I>
        struct Unroller<I,0>
        { static void unroll(V& v) {} };
        static void call(V& v)
        { Unroller<0,size/2>::unroll(v); }
    };

    // algo 90: Call inst
    template <int size, class V>
    struct ReverseV_Helper<90,size,V>
    {
        static void call(V& v)
        { InstReverseSelf(v.xView()); }
    };

    // algo 97: Conjugate
    template <int size, class V>
    struct ReverseV_Helper<97,size,V>
    {
        static void call(V& v)
        {
            typedef typename V::conjugate_type Vc;
            Vc vc = v.conjugate();
            ReverseV_Helper<-1,size,Vc>::call(vc);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int size, class V>
    struct ReverseV_Helper<-3,size,V>
    {
        static void call(V& v)
        {
            const int algo = 
                TMV_OPT == 0 ? 11 :
                size != UNKNOWN && size <= 32 ? 15 :
                V::iscomplex ? 12 :
                11;
            ReverseV_Helper<algo,size,V>::call(v);
        }
    };

    // algo -2: Check for inst
    template <int size, class V>
    struct ReverseV_Helper<-2,size,V>
    {
        static void call(V& v)
        {
            typedef typename V::value_type T;
            const bool inst =
                Traits<T>::isinst &&
                V::_size == UNKNOWN;
            const int algo = 
                V::_conj ? 97 :
                inst ? 90 :
                -3;
            ReverseV_Helper<algo,size,V>::call(v);
        }
    };

    template <int size, class V>
    struct ReverseV_Helper<-1,size,V>
    {
        static void call(V& v)
        { ReverseV_Helper<-2,size,V>::call(v); }
    };

    template <class V>
    static inline void ReverseSelf(BaseVector_Mutable<V>& v)
    {
        const int size = V::_size;
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(V,Vv) vv = v.cView();
        ReverseV_Helper<-2,size,Vv>::call(vv);
    }

    template <class V>
    static inline void InlineReverseSelf(BaseVector_Mutable<V>& v)
    {
        const int size = V::_size;
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(V,Vv) vv = v.cView();
        ReverseV_Helper<-3,size,Vv>::call(vv);
    }

} // namespace tmv

#endif
