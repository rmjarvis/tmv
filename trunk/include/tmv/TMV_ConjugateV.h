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

#ifndef TMV_ConjV_H
#define TMV_ConjV_H

#include "TMV_BaseVector.h"
#include "TMV_ScaleV.h"

namespace tmv {

    // Defined in TMV_Vector.cpp
    template <class T>
    void InstConjugateSelf(VectorView<T> v);

    //
    // ConjugateSelf
    //

    template <int algo, int s, class V>
    struct ConjugateV_Helper;

    // algo 0: real, so nothing to do
    template <int s, class V>
    struct ConjugateV_Helper<0,s,V>
    { static void call(V& ) { } };
    
    // algo 11: simple for loop
    template <int s, class V>
    struct ConjugateV_Helper<11,s,V>
    {
        static void call(V& v)
        {
            const int n=v.size();
            for(int i=0;i<n;++i) v.ref(i) = TMV_CONJ(v.cref(i));
        }
    };

    // algo 12: v.imagPart() *= -1
    template <int s, class V>
    struct ConjugateV_Helper<12,s,V>
    {
        static void call(V& v)
        {
            typedef typename V::real_type RT;
            typedef typename V::imagpart_type Vi;
            const Scaling<-1,RT> mone;
            Vi vi = v.imagPart();
            Scale(mone,vi);
        }
    };

    // algo 15: fully unroll
    template <int s, class V>
    struct ConjugateV_Helper<15,s,V>
    {
        template <int I, int N>
        struct Unroller
        {
            static void unroll(V& v)
            {
                Unroller<I,N/2>::unroll(v);
                Unroller<I+N/2,N-N/2>::unroll(v);
            }
        };
        template <int I>
        struct Unroller<I,1>
        {
            static void unroll(V& v)
            { v.ref(I) = TMV_CONJ(v.cref(I)); }
        };
        template <int I>
        struct Unroller<I,0>
        { static void unroll(V& v) {} };
        static void call(V& v)
        { Unroller<0,s>::unroll(v); }
    };

    // algo 90: Call inst
    template <int s, class V>
    struct ConjugateV_Helper<90,s,V>
    {
        static void call(V& v)
        { InstConjugateSelf(v.xView()); }
    };

    // algo 97: Conjugate
    template <int s, class V>
    struct ConjugateV_Helper<97,s,V>
    {
        static void call(V& v)
        {
            typedef typename V::conjugate_type Vc;
            Vc vc = v.conjugate();
            ConjugateV_Helper<-2,s,Vc>::call(vc);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int s, class V>
    struct ConjugateV_Helper<-3,s,V>
    {
        static void call(V& v)
        {
            const int algo = 
                TMV_OPT == 0 ? 12 :
                s != UNKNOWN && s <= 32 ? 15 :
                12;
            ConjugateV_Helper<algo,s,V>::call(v);
        }
    };

    // algo -2: Check for inst
    template <int s, class V>
    struct ConjugateV_Helper<-2,s,V>
    {
        static void call(V& v)
        {
            typedef typename V::value_type T;
            const bool inst =
                (s == UNKNOWN || s > 16) &&
                Traits<T>::isinst;
            const int algo = 
                V::isreal ? 0 :
                V::_conj ? 97 :
                inst ? 90 : 
                -3;
            ConjugateV_Helper<algo,s,V>::call(v);
        }
    };

    template <int s, class V>
    struct ConjugateV_Helper<-1,s,V>
    {
        static void call(V& v)
        { ConjugateV_Helper<-2,s,V>::call(v); }
    };

    template <class V>
    static inline void ConjugateSelf(BaseVector_Mutable<V>& v)
    {
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(V,Vv) vv = v.cView();
        ConjugateV_Helper<-2,V::_size,Vv>::call(vv); 
    }

    template <class V>
    static inline void InlineConjugateSelf(BaseVector_Mutable<V>& v)
    {
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(V,Vv) vv = v.cView();
        ConjugateV_Helper<-3,V::_size,Vv>::call(vv); 
    }

} // namespace tmv

#endif
