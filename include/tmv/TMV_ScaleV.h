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

#ifndef TMV_ScaleV_H
#define TMV_ScaleV_H

#include "TMV_BaseVector.h"
#include "TMV_Scaling.h"

namespace tmv {

    // Defined below:
    template <int ix, class T, class V>
    inline void Scale(const Scaling<ix,T>& x, BaseVector_Mutable<V>& v);
    template <int ix, class T, class V>
    inline void InlineScale(const Scaling<ix,T>& x, BaseVector_Mutable<V>& v);

    // Defined in TMV_ScaleV.cpp
    template <class T>
    void InstScale(const T x, VectorView<T> v);  

    //
    // Vector *= x
    //

    template <int algo, int size, int ix, class T, class V>
    struct ScaleV_Helper;

    // algo 0: trivial: size == 0 or ix == 1, so nothing to do
    template <int size, int ix, class T, class V>
    struct ScaleV_Helper<0,size,ix,T,V>
    { static inline void call(const Scaling<ix,T>& , V& ) {} };

    // algo 1: complex vector with unit step, convert to real version.
    template <int size, int ix, class T, class V>
    struct ScaleV_Helper<1,size,ix,T,V> 
    {
        typedef typename V::iterator IT;
        typedef typename V::flatten_type Vf;
        typedef typename Vf::iterator ITf;
        typedef typename V::real_type RT;
        enum { size2 = size == UNKNOWN ? UNKNOWN : (size<<1) };
        static inline void call(const Scaling<ix,T>& x, V& v)
        {
            Vf vf = v.flatten();
            ScaleV_Helper<-4,size2,ix,T,Vf>::call(x,vf);
        }
        static inline void call2(int n, const Scaling<ix,T>& x, IT it)
        { 
            ITf itf = it.flatten();
            const int n2 = n<<1;
            ScaleV_Helper<-4,size2,ix,T,Vf>::call2(n2,x,itf);
        }
    };

    // algo 11: simple loop
    template <int size, int ix, class T, class V>
    struct ScaleV_Helper<11,size,ix,T,V>
    {
        typedef typename V::iterator IT;
        static inline void call(const Scaling<ix,T>& x, V& v)
        {
            const int n = size == UNKNOWN ? int(v.size()) : size;
            call2(n,x,v.begin());
        }
        static inline void call2(int n, const Scaling<ix,T>& x, IT it)
        { 
            if (n) do {
                *it = ZProd<false,false>::prod(x,*it); ++it;
            } while (--n);
        }
    };

    // algo 12: 2 at a time
    template <int size, int ix, class T, class V>
    struct ScaleV_Helper<12,size,ix,T,V>
    {
        typedef typename V::iterator IT;
        static inline void call(const Scaling<ix,T>& x, V& v)
        {
            const int n = size == UNKNOWN ? int(v.size()) : size;
            call2(n,x,v.begin());
        }
        static inline void call2(const int n, const Scaling<ix,T>& x, IT it)
        {
            int n_2 = (n>>1);
            const int nb = n-(n_2<<1);

            if (n_2) do {
                it[0] = ZProd<false,false>::prod(x,it[0]);
                it[1] = ZProd<false,false>::prod(x,it[1]); it += 2;
            } while (--n_2);
            if (nb) *it = ZProd<false,false>::prod(x,*it);
        }
    };

    // algo 13: 4 at a time
    template <int size, int ix, class T, class V>
    struct ScaleV_Helper<13,size,ix,T,V>
    {
        typedef typename V::iterator IT;
        static inline void call(const Scaling<ix,T>& x, V& v)
        {
            const int n = size == UNKNOWN ? int(v.size()) : size;
            call2(n,x,v.begin());
        }
        static inline void call2(const int n, const Scaling<ix,T>& x, IT it)
        {
            int n_4 = (n>>2);
            int nb = n-(n_4<<2);

            if (n_4) do {
                it[0] = ZProd<false,false>::prod(x,it[0]);
                it[1] = ZProd<false,false>::prod(x,it[1]);
                it[2] = ZProd<false,false>::prod(x,it[2]);
                it[3] = ZProd<false,false>::prod(x,it[3]); it += 4;
            } while (--n_4);
            if (nb) do {
                *it = ZProd<false,false>::prod(x,*it); ++it;
            } while (--nb);
        }
    };

    // algo 15: fully unroll
    template <int size, int ix, class T, class V>
    struct ScaleV_Helper<15,size,ix,T,V> 
    {
        template <int I, int N>
        struct Unroller
        {
            static inline void unroll(const Scaling<ix,T>& x, V& v)
            {
                Unroller<I,N/2>::unroll(x,v);
                Unroller<I+N/2,N-N/2>::unroll(x,v);
            }
        };
        template <int I>
        struct Unroller<I,1>
        {
            static inline void unroll(const Scaling<ix,T>& x, V& v)
            { v.ref(I) = ZProd<false,false>::prod(x,v.cref(I)); }
        };
        template <int I>
        struct Unroller<I,0>
        { static inline void unroll(const Scaling<ix,T>& , V& ) {} };

        static inline void call(const Scaling<ix,T>& x, V& v)
        { Unroller<0,size>::unroll(x,v); }
    };

#ifdef __SSE__
    // algo 21: single precision SSE: all real
    template <int size, int ix, class T, class V>
    struct ScaleV_Helper<21,size,ix,T,V>
    {
        typedef typename V::iterator IT;
        static inline void call(const Scaling<ix,T>& x, V& v)
        {
            const int n = size == UNKNOWN ? int(v.size()) : size;
            call2(n,x,v.begin());
        }
        static inline void call2(int n, const Scaling<ix,T>& x, IT A)
        {
            const bool unit = V::vstep == 1;
            if (unit ) {
                while (n && (((unsigned int)(A.getP()) & 0xf) != 0) ) {
                    *A++ *= x; 
                    --n;
                }
            }

            int n_4 = (n>>2);
            int nb = n-(n_4<<2);

            if (n_4) {
                IT A1 = A+1;
                IT A2 = A+2;
                IT A3 = A+3;
                float x1 = x;
                __m128 xx = _mm_set1_ps(x1);
                __m128 xA;
                do {
                    Maybe<unit>::sse_load(
                        xA,A.getP(),A1.getP(),A2.getP(),A3.getP());
                    xA = _mm_mul_ps(xA,xx);
                    Maybe<unit>::sse_store(
                        A.getP(),A1.getP(),A2.getP(),A3.getP(),xA);
                    A+=4; A1+=4; A2+=4; A3+=4;
                } while (--n_4);
            }
            if (nb) do { *A++ *= x; } while (--nb);
        }
    };

    // algo 22: single precision SSE: x real, v complex
    template <int size, int ix, class T, class V>
    struct ScaleV_Helper<22,size,ix,T,V>
    {
        typedef typename V::iterator IT;
        static inline void call(const Scaling<ix,T>& x, V& v)
        {
            const int n = size == UNKNOWN ? int(v.size()) : size;
            call2(n,x,v.begin());
        }
        static inline void call2(int n, const Scaling<ix,T>& x, IT A)
        {
            const bool unit = V::vstep == 1;
            if (unit ) {
                while (n && (((unsigned int)(A.getP()) & 0xf) != 0) ) {
                    *A++ *= x; 
                    --n;
                }
            }

            int n_2 = (n>>1);
            int nb = n-(n_2<<1);

            if (n_2) {
                IT A1 = A+1;
                float x1 = x;
                __m128 xx = _mm_set1_ps(x1);
                __m128 xA;
                do {
                    Maybe<unit>::sse_load(xA,A.getP(),A1.getP());
                    xA = _mm_mul_ps(xA,xx);
                    Maybe<unit>::sse_store(A.getP(),A1.getP(),xA);
                    A+=2; A1+=2;
                } while (--n_2);
            }
            if (nb) *A++ *= x; 
        }
    };

    // algo 23: single precision SSE: all complex
    template <int size, int ix, class T, class V>
    struct ScaleV_Helper<23,size,ix,T,V>
    {
        typedef typename V::iterator IT;
        static inline void call(const Scaling<ix,T>& x, V& v)
        {
            const int n = size == UNKNOWN ? int(v.size()) : size;
            call2(n,x,v.begin());
        }
        static inline void call2(int n, const Scaling<ix,T>& x, IT A)
        {
            const bool unit = V::vstep == 1;

            if (unit) {
                while (n && (((unsigned int)(A.getP()) & 0xf) != 0) ) {
                    *A = ZProd<false,false>::prod(x,*A); ++A;
                    --n;
                }
            }

            int n_2 = (n>>1);
            int nb = n-(n_2<<1);

            if (n_2) {
                IT A1 = A+1;
                std::complex<float> xx(x);
                float xr = real(xx);
                float xi = imag(xx);
                __m128 xxr = _mm_set1_ps(xr);
                __m128 xxi = _mm_set_ps(xi , -xi , xi , -xi);
                __m128 xA;
                __m128 x0, x1, x2; // temp vars
                do {
                    // r = xr * Ar - xi * Ai
                    // i = xr * Ai + xi * Ar
                    Maybe<unit>::sse_load(xA,A.getP(),A1.getP());
                    x0 = _mm_shuffle_ps(xA,xA,_MM_SHUFFLE(2,3,0,1));
                    x1 = _mm_mul_ps(xxr,xA); // xr*Ar, xr*Ai
                    x2 = _mm_mul_ps(xxi,x0); // -xi*Ai, xi*Ar
                    xA = _mm_add_ps(x1,x2);
                    Maybe<unit>::sse_store(A.getP(),A1.getP(),xA);
                    A += 2; A1 += 2;
                } while (--n_2);
            }
            if (nb) *A = ZProd<false,false>::prod(x,*A);
        }
    };
#endif

#ifdef __SSE2__
    // algo 31: double precision SSE2: all real
    template <int size, int ix, class T, class V>
    struct ScaleV_Helper<31,size,ix,T,V>
    {
        typedef typename V::iterator IT;
        static inline void call(const Scaling<ix,T>& x, V& v)
        {
            const int n = size == UNKNOWN ? int(v.size()) : size;
            call2(n,x,v.begin());
        }
        static inline void call2(int n, const Scaling<ix,T>& x, IT A)
        {
            const bool unit = V::vstep == 1;
            if (unit) {
                while (n && (((unsigned int)(A.getP()) & 0xf) != 0) ) {
                    *A++ *= x; 
                    --n;
                }
            }

            int n_2 = (n>>1);
            int nb = n-(n_2<<1);

            if (n_2) {
                IT A1 = A+1;
                double x1 = x;
                __m128d xx = _mm_set1_pd(x1);
                __m128d xA;
                do {
                    Maybe<unit>::sse_load(xA,A.getP(),A1.getP());
                    xA = _mm_mul_pd(xA,xx);
                    Maybe<unit>::sse_store(A.getP(),A1.getP(),xA);
                    A+=2; A1+=2;
                } while (--n_2);
            }
            if (nb) *A++ *= x;
        }
    };

    // algo 32: double precision SSE2: x real, v complex
    template <int size, int ix, class T, class V>
    struct ScaleV_Helper<32,size,ix,T,V>
    {
        typedef typename V::iterator IT;
        static inline void call(const Scaling<ix,T>& x, V& v)
        {
            const int n = size == UNKNOWN ? int(v.size()) : size;
            call2(n,x,v.begin());
        }
        static inline void call2(int n, const Scaling<ix,T>& x, IT A)
        {
            if (n) {
                double x1 = x;
                __m128d xx = _mm_set1_pd(x1);
                __m128d xA;
                do {
                    Maybe<true>::sse_load(xA,A.getP());
                    xA = _mm_mul_pd(xA,xx);
                    Maybe<true>::sse_store(A.getP(),xA);
                    ++A;
                } while (--n);
            }
        }
    };

    // algo 33: double precision SSE2: all complex
    template <int size, int ix, class T, class V>
    struct ScaleV_Helper<33,size,ix,T,V>
    {
        typedef typename V::iterator IT;
        static inline void call(const Scaling<ix,T>& x, V& v)
        {
            const int n = size == UNKNOWN ? int(v.size()) : size;
            call2(n,x,v.begin());
        }
        static inline void call2(int n, const Scaling<ix,T>& x, IT A)
        {
            if (n) {
                std::complex<double> xx(x);
                double xr = real(xx);
                double xi = imag(xx);
                __m128d xxr = _mm_set1_pd(xr);
                __m128d xxi = _mm_set_pd(xi , -xi);
                __m128d xA;
                __m128d x0, x1, x2; // temp vars
                if (((unsigned int)(A.getP()) & 0xf) == 0) {
                    do {
                        // r = xr * Ar - xi * Ai
                        // i = xr * Ai + xi * Ar
                        Maybe<true>::sse_load(xA,A.getP());
                        x0 = _mm_shuffle_pd(xA,xA,_MM_SHUFFLE2(0,1));
                        x1 = _mm_mul_pd(xxr,xA);
                        x2 = _mm_mul_pd(xxi,x0);
                        xA = _mm_add_pd(x1,x2);
                        Maybe<true>::sse_store(A.getP(),xA);
                        ++A;
                    } while (--n);
                } else {
                    do {
                        Maybe<true>::sse_loadu(xA,A.getP());
                        x0 = _mm_shuffle_pd(xA,xA,_MM_SHUFFLE2(0,1));
                        x1 = _mm_mul_pd(xxr,xA);
                        x2 = _mm_mul_pd(xxi,x0);
                        xA = _mm_add_pd(x1,x2);
                        Maybe<true>::sse_storeu(A.getP(),xA);
                        ++A;
                    } while (--n);
                }
            }
        }
    };
#endif

    // algo -4: No branches or copies
    template <int size, int ix, class T, class V>
    struct ScaleV_Helper<-4,size,ix,T,V> 
    {
        typedef typename V::value_type VT;
        typedef typename V::real_type RT;
        typedef typename V::iterator IT;
        enum { unit = V::vstep == 1 };
        enum { vfloat = Traits2<RT,float>::sametype };
        enum { vdouble = Traits2<RT,double>::sametype };
        enum { xreal = Traits<T>::isreal };
        enum { xcomplex = Traits<T>::iscomplex };
        enum { vreal = V::visreal };
        enum { vcomplex = V::viscomplex };
        enum { algo = (
                (size == 0 || ix == 1) ? 0 :
#if TMV_OPT >= 1
                (unit && xreal && vcomplex) ? 1 :
#ifdef __SSE__
                (vfloat && xreal && vreal) ? 21 :
                (vfloat && xreal && vcomplex) ? 22 :
                (vfloat && xcomplex && vcomplex) ? 23 :
#endif
#ifdef __SSE2__
                (vdouble && xreal && vreal) ? 31 :
                (vdouble && xreal && vcomplex) ? 32 :
                (vdouble && xcomplex && vcomplex) ? 33 :
#endif
                (vreal && sizeof(RT) == 4) ? 13 :
                (vreal && sizeof(RT) == 8) ? 12 :
#endif
                11 ) };
        static inline void call(const Scaling<ix,T>& x, V& v)
        {
            TMVStaticAssert(!V::vconj);
            ScaleV_Helper<algo,size,ix,T,V>::call(x,v); 
        }
        static inline void call2(
            const int n, const Scaling<ix,T>& x, const IT& it)
        {
            TMVStaticAssert(!V::vconj);
            ScaleV_Helper<algo,size,ix,T,V>::call2(n,x,it); 
        }
    };

    // algo -3: Determine which algorithm to use
    template <int size, int ix, class T, class V>
    struct ScaleV_Helper<-3,size,ix,T,V> 
    {
        static inline void call(const Scaling<ix,T>& x, V& v)
        { ScaleV_Helper<-4,size,ix,T,V>::call(x,v); }
    };

    // algo 97: Conjugate
    template <int size, int ix, class T, class V>
    struct ScaleV_Helper<97,size,ix,T,V> 
    {
        static inline void call(const Scaling<ix,T>& x, V& v)
        {
            typedef typename V::conjugate_type Vc;
            Vc vc = v.conjugate();
            ScaleV_Helper<-1,size,ix,T,Vc>::call(TMV_CONJ(x),vc);
        }
    };

    // algo 98: Call inst
    template <int size, int ix, class T, class V>
    struct ScaleV_Helper<98,size,ix,T,V> 
    {
        static inline void call(const Scaling<ix,T>& x, V& v)
        {
            typename V::value_type xx(x);
            TMVStaticAssert(!V::vconj);
            InstScale(xx,v.xView());
        }
    };

    // algo -2: Check for inst
    template <int size, int ix, class T, class V>
    struct ScaleV_Helper<-2,size,ix,T,V> 
    {
        static inline void call(const Scaling<ix,T>& x, V& v)
        {
            typedef typename V::value_type T1;
            const bool inst =
                Traits<T>::isinst &&
                Traits<T1>::isinst &&
                Traits2<T,T1>::samebase &&
                V::vsize == UNKNOWN;
            const bool conj = V::vconj;
            const int algo = 
                ix == 1 ? 0 :
                conj ? 97 :
                inst ? 98 : 
                -4;
            ScaleV_Helper<algo,size,ix,T,V>::call(x,v);
        }
    };

    // algo -1: Check for aliases?  No.
    template <int size, int ix, class T, class V>
    struct ScaleV_Helper<-1,size,ix,T,V> 
    {
        static inline void call(const Scaling<ix,T>& x, V& v)
        { ScaleV_Helper<-2,size,ix,T,V>::call(x,v); }
    };

    template <int ix, class T, class V>
    inline void Scale(const Scaling<ix,T>& x, BaseVector_Mutable<V>& v)
    {
        typedef typename V::cview_type Vv;
        Vv vv = v.cView();
        ScaleV_Helper<-2,V::vsize,ix,T,Vv>::call(x,vv);
    }

    template <int ix, class T, class V>
    inline void InlineScale(const Scaling<ix,T>& x, BaseVector_Mutable<V>& v)
    {
        typedef typename V::cview_type Vv;
        Vv vv = v.cView();
        ScaleV_Helper<-4,V::vsize,ix,T,Vv>::call(x,vv);
    }

} // namespace tmv

#endif 
