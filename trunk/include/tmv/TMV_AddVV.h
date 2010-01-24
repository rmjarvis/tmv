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


#ifndef TMV_AddVV_H
#define TMV_AddVV_H

#include "TMV_BaseVector.h"
#include "TMV_Scaling.h"
#include "TMV_MultXV.h"

namespace tmv {

    //
    // Vector = x1 * Vector + x2 * Vector
    //

    template <int algo, int size, int ix1, class T1, class V1,
              int ix2, class T2, class V2, class V3>
    struct AddVV_Helper;

    // algo 0: size = 0, nothing to do
    template <int ix1, class T1, class V1,
              int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<0,0,ix1,T1,V1,ix2,T2,V2,V3> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix1,T1>& , const V1& ,
            const Scaling<ix2,T2>& , const V2& , V3& ) 
        {}
        static void call2(
            int , const Scaling<ix1,T1>& , IT1 ,
            const Scaling<ix2,T2>& , IT2 , IT3 )
        {}
    };  

    // algo 1: complex vectors with unit step, convert to real version
    template <int size, int ix1, class T1, class V1,
              int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<1,size,ix1,T1,V1,ix2,T2,V2,V3> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        typedef typename V1::const_flatten_type V1f;
        typedef typename V2::const_flatten_type V2f;
        typedef typename V3::flatten_type V3f;
        typedef typename V1f::const_nonconj_type::const_iterator IT1f;
        typedef typename V2f::const_nonconj_type::const_iterator IT2f;
        typedef typename V3f::iterator IT3f;
        enum { size2 = size == UNKNOWN ? UNKNOWN : (size<<1) };
        static inline void call(
            const Scaling<ix1,T1>& x1, const V1& v1,
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            V1f v1f = v1.flatten();
            V2f v2f = v2.flatten();
            V3f v3f = v3.flatten();
            AddVV_Helper<-4,size2,ix1,T1,V1f,ix2,T2,V2f,V3f>::call(
                x1,v1f,x2,v2f,v3f);
        }
        static void call2(
            int n, const Scaling<ix1,T1>& x1, IT1 it1,
            const Scaling<ix2,T2>& x2, IT2 it2, IT3 it3)
        {
            IT1f it1f = it1.flatten();
            IT2f it2f = it2.flatten();
            IT3f it3f = it3.flatten();
            const int n2 = n<<1;
            AddVV_Helper<-4,size2,ix1,T1,V1f,ix2,T2,V2f,V3f>::call(
                n2,x1,it1f,x2,it2f,it3f);
        }
    };

    // algo 11: simple for loop
    template <int size, int ix1, class T1, class V1,
              int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<11,size,ix1,T1,V1,ix2,T2,V2,V3> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix1,T1>& x1, const V1& v1,
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            const int n = size == UNKNOWN ? int(v3.size()) : size;
            for(int i=0;i<n;++i) 
                v3.ref(i) = ZProd<false,false>::prod(x1 , v1.cref(i)) +
                    ZProd<false,false>::prod(x2 , v2.cref(i));
        }
        static void call2(
            int n, const Scaling<ix1,T1>& x1, IT1 it1,
            const Scaling<ix2,T2>& x2, IT2 it2, IT3 it3)
        {
            const bool c1 = V1::vconj;
            const bool c2 = V2::vconj;
            for(;n;--n) 
                *it3++ = ZProd<false,c1>::prod(x1 , *it1++) + 
                    ZProd<false,c2>::prod(x2 , *it2++); 
        }
    };  

    // algo 12: 2 at a time
    template <int size, int ix1, class T1, class V1,
              int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<12,size,ix1,T1,V1,ix2,T2,V2,V3> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static inline void call(
            const Scaling<ix1,T1>& x1, const V1& v1,
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            const int n = size == UNKNOWN ? int(v3.size()) : size;
            call2(n,x1,v1.nonConj().begin(),x2,v2.nonConj().begin(),v3.begin());
        }
        static void call2(
            const int n, const Scaling<ix1,T1>& x1, IT1 it1,
            const Scaling<ix2,T2>& x2, IT2 it2, IT3 it3)
        {
            int n_2 = (n>>1);
            const int nb = n-(n_2<<1);
            const bool c1 = V1::vconj;
            const bool c2 = V2::vconj;

            if (n_2) do {
                it3[0] = ZProd<false,c1>::prod(x1 , it1[0]) +
                    ZProd<false,c2>::prod(x2 , it2[0]);
                it3[1] = (
                    ZProd<false,c1>::prod(x1 , it1[1]) +
                    ZProd<false,c2>::prod(x2 , it2[1]));
                it1 += 2; it2 += 2; it3 += 2;
            } while (--n_2);
            if (nb) {
                *it3 = ZProd<false,c1>::prod(x1 , *it1) +
                    ZProd<false,c2>::prod(x2 , *it2); 
            }
        }
    };

    // algo 13: 4 at a time
    template <int size, int ix1, class T1, class V1,
              int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<13,size,ix1,T1,V1,ix2,T2,V2,V3> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static inline void call(
            const Scaling<ix1,T1>& x1, const V1& v1,
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            const int n = size == UNKNOWN ? int(v3.size()) : size;
            call2(n,x1,v1.nonConj().begin(),x2,v2.nonConj().begin(),v3.begin());
        }
        static void call2(
            const int n, const Scaling<ix1,T1>& x1, IT1 it1,
            const Scaling<ix2,T2>& x2, IT2 it2, IT3 it3)
        {
            int n_4 = (n>>2);
            int nb = n-(n_4<<2);
            const bool c1 = V1::vconj;
            const bool c2 = V2::vconj;

            if (n_4) do {
                it3[0] = ZProd<false,c1>::prod(x1 , it1[0]) +
                    ZProd<false,c2>::prod(x2 , it2[0]);
                it3[1] = ZProd<false,c1>::prod(x1 , it1[1]) +
                    ZProd<false,c2>::prod(x2 , it2[1]);
                it3[2] = ZProd<false,c1>::prod(x1 , it1[2]) +
                    ZProd<false,c2>::prod(x2 , it2[2]);
                it3[3] = ZProd<false,c1>::prod(x1 , it1[3]) +
                    ZProd<false,c2>::prod(x2 , it2[3]);
                it1 += 4; it2 += 4; it3 += 4;
            } while (--n_4);
            if (nb) do {
                *it3++ = ZProd<false,c1>::prod(x1 , *it1++) +
                    ZProd<false,c2>::prod(x2 , *it2++); 
            } while (--nb);
        }
    };

    // algo 14: 8 at a time
    template <int size, int ix1, class T1, class V1,
              int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<14,size,ix1,T1,V1,ix2,T2,V2,V3> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static inline void call(
            const Scaling<ix1,T1>& x1, const V1& v1,
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            const int n = size == UNKNOWN ? int(v3.size()) : size;
            call2(n,x1,v1.nonConj().begin(),x2,v2.nonConj().begin(),v3.begin());
        }
        static void call2(
            const int n, const Scaling<ix1,T1>& x1, IT1 it1,
            const Scaling<ix2,T2>& x2, IT2 it2, IT3 it3)
        {
            int n_8 = (n>>3);
            int nb = n-(n_8<<3);
            const bool c1 = V1::vconj;
            const bool c2 = V2::vconj;

            if (n_8) do {
                it3[0] = ZProd<false,c1>::prod(x1 , it1[0]) +
                    ZProd<false,c2>::prod(x2 , it2[0]);
                it3[1] = ZProd<false,c1>::prod(x1 , it1[1]) +
                    ZProd<false,c2>::prod(x2 , it2[1]);
                it3[2] = ZProd<false,c1>::prod(x1 , it1[2]) +
                    ZProd<false,c2>::prod(x2 , it2[2]);
                it3[3] = ZProd<false,c1>::prod(x1 , it1[3]) +
                    ZProd<false,c2>::prod(x2 , it2[3]);
                it3[4] = ZProd<false,c1>::prod(x1 , it1[4]) +
                    ZProd<false,c2>::prod(x2 , it2[4]);
                it3[5] = ZProd<false,c1>::prod(x1 , it1[5]) +
                    ZProd<false,c2>::prod(x2 , it2[5]);
                it3[6] = ZProd<false,c1>::prod(x1 , it1[6]) +
                    ZProd<false,c2>::prod(x2 , it2[6]);
                it3[7] = ZProd<false,c1>::prod(x1 , it1[7]) +
                    ZProd<false,c2>::prod(x2 , it2[7]);
                it1 += 8; it2 += 8; it3 += 8;
            } while (--n_8);
            if (nb) do {
                *it3++ = ZProd<false,c1>::prod(x1 , *it1++) +
                    ZProd<false,c2>::prod(x2 , *it2++); 
            } while (--nb);
        }
    };

    // algo 15: fully unroll
    template <int size, int ix1, class T1, class V1,
              int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<15,size,ix1,T1,V1,ix2,T2,V2,V3> // known size, unroll
    {
        template <int I, int N>
        struct Unroller
        {
            static inline void unroll(
                const Scaling<ix1,T1>& x1, const V1& v1,
                const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
            {
                Unroller<I,N/2>::unroll(x1,v1,x2,v2,v3);
                Unroller<I+N/2,N-N/2>::unroll(x1,v1,x2,v2,v3);
            }
        };
        template <int I>
        struct Unroller<I,1>
        {
            static inline void unroll(
                const Scaling<ix1,T1>& x1, const V1& v1,
                const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
            {
                v3.ref(I) = ZProd<false,false>::prod(x1 , v1.cref(I)) +
                    ZProd<false,false>::prod(x2 , v2.cref(I)); 
            }
        };
        template <int I>
        struct Unroller<I,0>
        {
            static inline void unroll(
                const Scaling<ix1,T1>& x1, const V1& v1,
                const Scaling<ix2,T2>& x2, const V2& v2, V3& v3) {}
        };
        static inline void call(
            const Scaling<ix1,T1>& x1, const V1& v1,
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        { Unroller<0,size>::unroll(x1,v1,x2,v2,v3); }
    };

    // algo -4: No branches or copies
    template <int size, int ix1, class T1, class V1,
              int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<-4,size,ix1,T1,V1,ix2,T2,V2,V3> 
    {
        typedef typename V3::value_type VT;
        typedef typename V2::real_type RT;
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        enum { allreal = V1::visreal && V2::visreal && V3::visreal };
        enum { allcomplex = (
                V1::viscomplex && V2::viscomplex && V3::viscomplex ) };
        enum { allunit = V1::vstep == 1 && V2::vstep == 1 && V3::vstep == 1 };
        enum { flatten = (
                allunit && allcomplex &&
                Traits<T1>::isreal && Traits<T2>::isreal &&
                V1::vconj == int(V3::vconj) && V2::vconj == int(V3::vconj) ) };
        enum { algo =  (
                size == 0 ? 0 : 
#if TMV_OPT >= 1
                flatten ? 1 :
                (sizeof(RT) == 8 && allunit && allreal) ? 12 :
                (sizeof(RT) == 4 && allunit && allreal) ? 13 :
#endif
                11 ) };
        static inline void call(
            const Scaling<ix1,T1>& x1, const V1& v1,
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            TMVStaticAssert(!V3::vconj);
            const int algo1 = 
                (size != UNKNOWN && size <= int(128/sizeof(VT))) ? 15 :
                algo;
            AddVV_Helper<algo1,size,ix1,T1,V1,ix2,T2,V2,V3>::call(
                x1,v1,x2,v2,v3);
        }
        static inline void call2(
            const int n, const Scaling<ix1,T1>& x1, const IT1& it1,
            const Scaling<ix2,T2>& x2, const IT2& it2, const IT3& it3)
        {
            TMVStaticAssert(!V3::vconj);
            AddVV_Helper<algo,size,ix1,T1,V1,ix2,T2,V2,V3>::call2(
                n,x1,it1,x2,it2,it3);
        } 
    };

    // algo -3: Determine which algorithm to use
    template <int size, int ix1, class T1, class V1,
              int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<-3,size,ix1,T1,V1,ix2,T2,V2,V3> 
    {
        static inline void call(
            const Scaling<ix1,T1>& x1, const V1& v1,
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            AddVV_Helper<-4,size,ix1,T1,V1,ix2,T2,V2,V3>::call(
                x1,v1,x2,v2,v3);
        }
    };

    // algo 97: Call inst
    template <int s, int ix1, class T1, class V1,
              int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<97,s,ix1,T1,V1,ix2,T2,V2,V3>
    {
        static inline void call(
            const Scaling<ix1,T1>& x1, const V1& v1, 
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            NoAliasMultXV<false>(x1,v1,v3);
            NoAliasMultXV<true>(x2,v2,v3);
        }
    };
    // algo -2: Check for inst
    template <int s, int ix1, class T1, class V1,
              int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<-2,s,ix1,T1,V1,ix2,T2,V2,V3>
    {
        static inline void call(
            const Scaling<ix1,T1>& x1, const V1& v1, 
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            typedef typename V1::value_type TV1;
            typedef typename V2::value_type TV2;
            typedef typename V3::value_type TV3;
            const bool inst =
                V1::vsize == UNKNOWN &&
                V2::vsize == UNKNOWN &&
                V3::vsize == UNKNOWN &&
#ifdef TMV_INST_MIX
                Traits2<TV1,TV2>::samebase &&
                Traits2<TV1,TV3>::samebase &&
#else
                Traits2<TV1,TV2>::sametype &&
                Traits2<TV1,TV3>::sametype &&
#endif
                Traits<TV3>::isinst;
            const int algo = 
                inst ? 97 :
                -4;
            AddVV_Helper<algo,s,ix1,T1,V1,ix2,T2,V2,V3>::call(x1,v1,x2,v2,v3);
        }
    };

    // algo 98: Check for aliases when calling Inst functions
    // We don't have a separate Inst function for this.  We just
    // split the operation into two parts: 
    // v3 = x1*v1; 
    // v3 += x2*v2;
    // and let those operations call their Inst functions.
    // However, the alias requirements for this are different than the
    // normal ones, so we need to put some alias checking here.
    template <int s, int ix1, class T1, class V1,
              int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<98,s,ix1,T1,V1,ix2,T2,V2,V3>
    {
        static inline void call(
            const Scaling<ix1,T1>& x1, const V1& v1, 
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            const bool s1 = SameStorage(v1,v3);
            const bool s2 = SameStorage(v2,v3);

            if (!s1 && !s2) {
                // No aliasing
                NoAliasMultXV<false>(x1,v1,v3);
                NoAliasMultXV<true>(x2,v2,v3);
            } else if (!s2) { 
                // Alias with v1 only, do v1 first
                AliasMultXV<false>(x1,v1,v3);
                NoAliasMultXV<true>(x2,v2,v3);
            } else if (!s1) { 
                // Alias with v2 only, do v2 first
                AliasMultXV<false>(x2,v2,v3);
                NoAliasMultXV<true>(x1,v1,v3);
            } else { 
                // Need a temporary
                typename V1::copy_type v1c = v1;
                AliasMultXV<false>(x2,v2,v3);
                NoAliasMultXV<true>(x1,v1c,v3);
            }
        }
    };

    // algo 99: Check for aliases when not calling Inst functions
    template <int s, int ix1, class T1, class V1,
              int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<99,s,ix1,T1,V1,ix2,T2,V2,V3>
    {
        static inline void call(
            const Scaling<ix1,T1>& x1, const V1& v1, 
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            const bool ss1 = SameStorage(v1,v3);
            const bool ss2 = SameStorage(v2,v3);
            const bool s1 = ss1 && !ExactSameStorage(v1,v3);
            const bool s2 = ss2 && !ExactSameStorage(v2,v3);

            if (!s1 && !s2) {
                // No aliasing (or no clobbering)
                AddVV_Helper<-4,s,ix1,T1,V1,ix2,T2,V2,V3>::call(
                    x1,v1,x2,v2,v3);
            } else if (!ss2) { 
                // Alias with v1 only, do v1 first
                AliasMultXV<false>(x1,v1,v3);
                NoAliasMultXV<true>(x2,v2,v3);
            } else if (!ss1) { 
                // Alias with v2 only, do v2 first
                AliasMultXV<false>(x2,v2,v3);
                NoAliasMultXV<true>(x1,v1,v3);
            } else { 
                // Need a temporary
                typename V1::copy_type v1c = v1;
                AliasMultXV<false>(x2,v2,v3);
                NoAliasMultXV<true>(x1,v1c,v3);
            }
        }
    };

    // algo -1: Check for aliases?
    template <int s, int ix1, class T1, class V1,
              int ix2, class T2, class V2, class V3>
    struct AddVV_Helper<-1,s,ix1,T1,V1,ix2,T2,V2,V3>
    {
        static inline void call(
            const Scaling<ix1,T1>& x1, const V1& v1, 
            const Scaling<ix2,T2>& x2, const V2& v2, V3& v3)
        {
            typedef typename V1::value_type TV1;
            typedef typename V2::value_type TV2;
            typedef typename V3::value_type TV3;
            const bool inst =
                V1::vsize == UNKNOWN &&
                V2::vsize == UNKNOWN &&
                V3::vsize == UNKNOWN &&
#ifdef TMV_INST_MIX
                Traits2<TV1,TV2>::samebase &&
                Traits2<TV1,TV3>::samebase &&
#else
                Traits2<TV1,TV2>::sametype &&
                Traits2<TV1,TV3>::sametype &&
#endif
                Traits<TV3>::isinst;
            const bool checkalias =
                V1::vsize == UNKNOWN &&
                V2::vsize == UNKNOWN &&
                V3::vsize == UNKNOWN;
            const int algo = 
                // We do a different check alias with the Inst calls.
                inst ? 98 : 
                checkalias ? 99 : 
                -4;
            AddVV_Helper<algo,s,ix1,T1,V1,ix2,T2,V2,V3>::call(x1,v1,x2,v2,v3);
        }
    };

    template <int ix1, class T1, class V1,
              int ix2, class T2, class V2, class V3>
    inline void AddVV(
        const Scaling<ix1,T1>& x1, const BaseVector_Calc<V1>& v1,
        const Scaling<ix2,T2>& x2, const BaseVector_Calc<V2>& v2,
        BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
        TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
        TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v1.size() == v3.size());
        const int size =
            Sizes<Sizes<V1::vsize,V2::vsize>::size,V3::vsize>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        V1v v1v = v1.cView();
        V2v v2v = v2.cView();
        V3v v3v = v3.cView();
        AddVV_Helper<-1,size,ix1,T1,V1v,ix2,T2,V2v,V3v>::call(
            x1,v1v,x2,v2v,v3v);
    }

    template <int ix1, class T1, class V1,
              int ix2, class T2, class V2, class V3>
    inline void NoAliasAddVV(
        const Scaling<ix1,T1>& x1, const BaseVector_Calc<V1>& v1,
        const Scaling<ix2,T2>& x2, const BaseVector_Calc<V2>& v2,
        BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
        TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
        TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v1.size() == v3.size());
        const int size =
            Sizes<Sizes<V1::vsize,V2::vsize>::size,V3::vsize>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        V1v v1v = v1.cView();
        V2v v2v = v2.cView();
        V3v v3v = v3.cView();
        AddVV_Helper<-2,size,ix1,T1,V1v,ix2,T2,V2v,V3v>::call(
            x1,v1v,x2,v2v,v3v);
    }

    template <int ix1, class T1, class V1,
              int ix2, class T2, class V2, class V3>
    inline void InlineAddVV(
        const Scaling<ix1,T1>& x1, const BaseVector_Calc<V1>& v1,
        const Scaling<ix2,T2>& x2, const BaseVector_Calc<V2>& v2,
        BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
        TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
        TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v1.size() == v3.size());
        const int size =
            Sizes<Sizes<V1::vsize,V2::vsize>::size,V3::vsize>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        V1v v1v = v1.cView();
        V2v v2v = v2.cView();
        V3v v3v = v3.cView();
        AddVV_Helper<-4,size,ix1,T1,V1v,ix2,T2,V2v,V3v>::call(
            x1,v1v,x2,v2v,v3v);
    }

    template <int ix1, class T1, class V1,
              int ix2, class T2, class V2, class V3>
    inline void AliasAddVV(
        const Scaling<ix1,T1>& x1, const BaseVector_Calc<V1>& v1,
        const Scaling<ix2,T2>& x2, const BaseVector_Calc<V2>& v2,
        BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
        TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
        TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v1.size() == v3.size());
        const int size =
            Sizes<Sizes<V1::vsize,V2::vsize>::size,V3::vsize>::size;
        typedef typename V1::value_type TV1;
        typedef typename V2::value_type TV2;
        typedef typename V3::value_type TV3;
        const bool inst =
            V1::vsize == UNKNOWN &&
            V2::vsize == UNKNOWN &&
            V3::vsize == UNKNOWN &&
#ifdef TMV_INST_MIX
            Traits2<TV1,TV2>::samebase &&
            Traits2<TV1,TV3>::samebase &&
#else
            Traits2<TV1,TV2>::sametype &&
            Traits2<TV1,TV3>::sametype &&
#endif
            Traits<TV3>::isinst;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        typedef typename V3::cview_type V3v;
        V1v v1v = v1.cView();
        V2v v2v = v2.cView();
        V3v v3v = v3.cView();
        AddVV_Helper<inst?98:99,size,ix1,T1,V1v,ix2,T2,V2v,V3v>::call(
            x1,v1v,x2,v2v,v3v);
    }

} // namespace tmv

#endif 
