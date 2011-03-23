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


#ifndef TMV_InvertD_H
#define TMV_InvertD_H

#include "TMV_BaseMatrix_Diag.h"


namespace tmv {


    // Defined in TMV_InvertD.cpp
    template <class T>
    void InstElemInvert(VectorView<T> v);

    //
    // DiagMatrix = DiagMatrix^-1
    //

    template <int algo, int s, class V>
    struct ElemInvert_Helper;

    // algo 0: s == 0, nothing to do
    template <class V>
    struct ElemInvert_Helper<0,0,V>
    {
        typedef typename V::iterator IT;
        static void call(V& ) {}
        static void call2(int, IT ) {}
    };

    // algo 11: simple for loop
    template <int s, class V>
    struct ElemInvert_Helper<11,s,V>
    {
        typedef typename V::iterator IT;
        static void call(V& v)
        {
            const int n = s == UNKNOWN ? int(v.size()) : s;
            call2(n,v.begin());
        }
        static void call2(int n, IT A)
        {
            typedef typename V::real_type RT;
            if (n) do {
                *A = ZProd<false,false>::quot(RT(1) , *A); ++A;
            } while (--n);
        }
    };

#ifdef __SSE__
    // algo 21: single precision SSE: v real
    // SSE also has a _mm_rcp_ps command that takes a reciprocal.
    // However, it is only 12 bit accurate, which is not acceptable,
    // so unfortunately we can't use it.
    // Instead, we do the reciprocal as 1 / x with _mm_div_ps.
    template <int s, class V>
    struct ElemInvert_Helper<21,s,V>
    {
        typedef typename V::iterator IT;
        static void call(V& v)
        {
            const int n = s == UNKNOWN ? int(v.size()) : s;
            call2(n,v.begin());
        }
        static void call2(int n, IT A)
        {
            const bool unit = V::_step == 1;
            const typename V::real_type one(1);

            if (unit) {
                while (n && !TMV_Aligned(A.get()) ) {
                    *A = one / *A;
                    ++A; --n;
                }
            }

            int n_4 = (n>>2);
            int nb = n-(n_4<<2);
            
            if (n_4) {
                IT A1 = A+1;
                IT A2 = A+2;
                IT A3 = A+3;

                __m128 xone = _mm_set1_ps(1.F);
                __m128 xA,xB;
                do {
                    Maybe<unit>::sse_load(
                        xA,A.get(),A1.get(),A2.get(),A3.get());
                    xB = _mm_div_ps(xone,xA);
                    Maybe<unit>::sse_store(
                        A.get(),A1.get(),A2.get(),A3.get(),xB);
                    A+=4; A1+=4; A2+=4; A3+=4;
                } while (--n_4);
            }

            if (nb) do { *A = one / *A; ++A; } while (--nb);
        }
    };

    // algo 23: single precision SSE: v complex
    template <int s, class V>
    struct ElemInvert_Helper<23,s,V>
    {
        typedef typename V::iterator IT;
        static void call(V& v)
        {
            const int n = s == UNKNOWN ? int(v.size()) : s;
            call2(n,v.begin());
        }
        static void call2(int n, IT A)
        {
            const bool unit = V::_step == 1;
            const typename V::real_type one(1);

            if (unit) {
                while (n && !TMV_Aligned(A.get()) ) {
                    *A = ZProd<false,false>::quot(one , *A);
                    ++A; --n;
                }
            }

            int n_2 = (n>>1);
            int nb = n-(n_2<<1);
            
            if (n_2) {
                IT A1 = A+1;

                // These look backwards, but order is from hi to lo values.
                __m128 xmone = _mm_set_ps(-1, 1, -1, 1);
                __m128 xA, xB;
                __m128 xAc, xnorm, x1, x2; // temp values
                do {
                    Maybe<unit>::sse_load(xA,A.get(),A1.get());
                    xAc = _mm_mul_ps(xmone,xA); // conj(xA)
                    x1 = _mm_mul_ps(xA,xA);
                    x2 = _mm_shuffle_ps(x1,x1,_MM_SHUFFLE(2,3,0,1));
                    xnorm = _mm_add_ps(x1,x2); // = norm(xA)
                    xB = _mm_div_ps(xAc,xnorm);  // = 1/xA
                    Maybe<unit>::sse_store(A.get(),A1.get(),xB);
                    A+=2; A1+=2;
                } while (--n_2);
            }

            if (nb) *A = ZProd<false,false>::quot(one , *A);
        }
    };
#endif

#ifdef __SSE2__
    // algo 31: double precision SSE2: v real
    template <int s, class V>
    struct ElemInvert_Helper<31,s,V>
    {
        typedef typename V::iterator IT;
        static void call(V& v)
        {
            const int n = s == UNKNOWN ? int(v.size()) : s;
            call2(n,v.begin());
        }
        static void call2(int n, IT A)
        {
            const bool unit = V::_step == 1;
            const typename V::real_type one(1);

            if (unit) {
                while (n && !TMV_Aligned(A.get()) ) {
                    *A = one / *A;
                    ++A; --n;
                }
            }

            int n_2 = (n>>1);
            int nb = n-(n_2<<1);
            
            if (n_2) {
                IT A1 = A+1;

                __m128d xone = _mm_set1_pd(1.);
                __m128d xA,xB;
                do {
                    Maybe<unit>::sse_load(xA,A.get(),A1.get());
                    xB = _mm_div_pd(xone,xA);
                    Maybe<unit>::sse_store(A.get(),A1.get(),xB);
                    A+=2; A1+=2;
                } while (--n_2);
            }

            if (nb) *A = one / *A;
        }
    };

    // algo 33: double precision SSE2: v complex
    template <int s, class V>
    struct ElemInvert_Helper<33,s,V>
    {
        typedef typename V::iterator IT;
        static void call(V& v)
        {
            const int n = s == UNKNOWN ? int(v.size()) : s;
            call2(n,v.begin());
        }
        static void call2(int n, IT A)
        {
            if (n) {
                // These look backwards, but order is from hi to lo values.
                __m128d xmone = _mm_set_pd(-1 , 1);
                __m128d xA,xB;
                __m128d xAc, xnorm, x1, x2; // temp values
                if (TMV_Aligned(A.get()) ) {
                    do {
                        Maybe<true>::sse_load(xA,A.get());
                        xAc = _mm_mul_pd(xmone,xA); // conj(xA)
                        x1 = _mm_mul_pd(xA,xA);
                        x2 = _mm_shuffle_pd(x1,x1,_MM_SHUFFLE2(0,1));
                        xnorm = _mm_add_pd(x1,x2); // = norm(xA)
                        xB = _mm_div_pd(xAc,xnorm);  // = 1/xA
                        Maybe<true>::sse_store(A.get(),xB); ++A;
                    } while (--n);
                } else {
                    do {
                        Maybe<true>::sse_loadu(xA,A.get());
                        xAc = _mm_mul_pd(xmone,xA); // conj(xA)
                        x1 = _mm_mul_pd(xA,xA);
                        x2 = _mm_shuffle_pd(x1,x1,_MM_SHUFFLE2(0,1));
                        xnorm = _mm_add_pd(x1,x2); // = norm(xA)
                        xB = _mm_div_pd(xAc,xnorm);  // = 1/xA
                        Maybe<true>::sse_storeu(A.get(),xB); ++A;
                    } while (--n);
                }
            }
        }
    };
#endif

    // algo 90: Call inst
    template <int s, class V>
    struct ElemInvert_Helper<90,s,V>
    {
        static void call(V& v)
        { InstElemInvert(v.xView()); }
    };

    // algo 97: Conjugate
    template <int s, class V>
    struct ElemInvert_Helper<97,s,V>
    {
        static void call(V& v)
        {
            typedef typename V::conjugate_type Vc;
            Vc vc = v.conjugate();
            ElemInvert_Helper<-2,s,Vc>::call(vc);
        }
    };

    // algo -4: No branches or copies
    template <int s, class V>
    struct ElemInvert_Helper<-4,s,V>
    {
        typedef typename V::iterator IT;
        typedef typename V::real_type RT;
        enum { vfloat = Traits2<RT,float>::sametype };
        enum { vdouble = Traits2<RT,double>::sametype };
        enum { vreal = V::isreal };
        enum { vcomplex = V::iscomplex };
        enum { unit = V::_step == 1 };

#if 0
        enum { algo = 11 };
#else
        enum { algo = (
                s == 0 ? 0 : 
                TMV_OPT == 0? 11 :
#ifdef __SSE__
                ( vfloat && vreal && unit ) ? 21 :
                ( vfloat && vcomplex ) ? 23 :
#endif
#ifdef __SSE2__
                ( vdouble && vreal && unit ) ? 31 :
                ( vdouble && vcomplex ) ? 33 :
#endif
                11 ) };
#endif
        static void call(V& v)
        {
            TMVStaticAssert(!V::_conj);
#ifdef PRINTALGO_InvD
            std::cout<<"InlineElemInvert:  \n";
            std::cout<<"v = "<<TMV_Text(v)<<std::endl;
            std::cout<<"s = "<<s<<" = "<<v.size()<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
            std::cout<<"v = "<<v<<std::endl;
#endif
            ElemInvert_Helper<algo,s,V>::call(v); 
#ifdef PRINTALGO_InvD
            std::cout<<"v => "<<v<<std::endl;
#endif
        }
        static void call2(int n, IT A)
        {
            TMVStaticAssert(!V::_conj);
            ElemInvert_Helper<algo,s,V>::call2(n,A); 
        }
    };

    // algo -3: Determine which algorithm to use
    template <int s, class V>
    struct ElemInvert_Helper<-3,s,V>
    {
        static void call(V& v)
        { ElemInvert_Helper<-4,s,V>::call(v); }
    };

    // algo -2: Check for inst
    template <int s, class V>
    struct ElemInvert_Helper<-2,s,V>
    {
        static void call(V& v)
        {
            typedef typename V::value_type T;
            const bool inst = 
                (s == UNKNOWN || s > 16) &&
                Traits<T>::isinst;
            const int algo = 
                s == 0 ? 0 : 
                V::_conj ? 97 :
                inst ? 90 :
                -4;
            ElemInvert_Helper<algo,s,V>::call(v);
        }
    };

    template <int s, class V>
    struct ElemInvert_Helper<-1,s,V>
    {
        static void call(V& v)
        { ElemInvert_Helper<-2,s,V>::call(v); }
    };

    //
    // Element-wise v(i) = 1/v(i)
    //

    template <class V>
    static inline void ElemInvert(BaseVector_Mutable<V>& v)
    {
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(V,Vv) vv = v.cView();
        ElemInvert_Helper<-2,V::_size,Vv>::call(vv);
    }

    template <class V>
    static inline void InlineElemInvert(BaseVector_Mutable<V>& v)
    {
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(V,Vv) vv = v.cView();
        ElemInvert_Helper<-3,V::_size,Vv>::call(vv);
    }

    template <class M>
    static inline void InvertSelf(BaseMatrix_Diag_Mutable<M>& m)
    {
        if (m.isSingular()) {
#ifdef TMV_NO_THROW
            std::cerr<<"Singular DiagMatrix found\n";
            exit(1);
#else
            throw Singular("DiagMatrix found\n");
#endif
        }
        typename M::diag_type md = m.diag();
        ElemInvert(md);
    }

    template <class M>
    static inline void InlineInvertSelf(BaseMatrix_Diag_Mutable<M>& m)
    {
        if (m.isSingular()) {
#ifdef TMV_NO_THROW
            std::cerr<<"Singular DiagMatrix found\n";
            exit(1);
#else
            throw Singular("DiagMatrix found\n");
#endif
        }
        typename M::diag_type md = m.diag();
        InlineElemInvert(md);
    }


} // namespace tmv

#endif 
