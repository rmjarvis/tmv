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
#include "TMV_Scaling.h"


namespace tmv {

    // Defined below:
    template <int ix, class T, class M1, class M2>
    inline void Invert(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Diag_Mutable<M2>& m2);
    template <int ix, class T, class M1, class M2>
    inline void NoAliasInvert(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Diag_Mutable<M2>& m2);
    template <int ix, class T, class M1, class M2>
    inline void InlineInvert(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Diag_Mutable<M2>& m2);
    template <int ix, class T, class M1, class M2>
    inline void AliasInvert(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Diag_Mutable<M2>& m2);

    template <int ix, class T, class V1, class V2>
    inline void ElemInvert(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2);
    template <int ix, class T, class V1, class V2>
    inline void NoAliasElemInvert(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2);
    template <int ix, class T, class V1, class V2>
    inline void InlineElemInvert(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2);
    template <int ix, class T, class V1, class V2>
    inline void AliasElemInvert(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2);

    // Defined in TMV_InvertD.cpp
    template <class T1, bool C1, class T2>
    void InstInvert(
        const T2 x,
        const ConstDiagMatrixView<T1,UNKNOWN,C1>& m1, DiagMatrixView<T2> m2);
    template <class T1, bool C1, class T2>
    void InstElemInvert(
        const T2 x,
        const ConstVectorView<T1,UNKNOWN,C1>& v1, VectorView<T2> v2);

    //
    // DiagMatrix = DiagMatrix^-1
    //

    template <int algo, int size, int ix, class T, class V1, class V2>
    struct ElemInvert_Helper;

    // algo 0: size == 0, nothing to do
    template <int ix, class T, class V1, class V2>
    struct ElemInvert_Helper<0,0,ix,T,V1,V2> 
    { 
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static inline void call(const Scaling<ix,T>& , const V1& , V2& ) {} 
        static inline void call2(int , const Scaling<ix,T>& , IT1 , IT2 ) {}
    };

    // algo 11: simple for loop
    template <int size, int ix, class T, class V1, class V2>
    struct ElemInvert_Helper<11,size,ix,T,V1,V2> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            const int n = size == UNKNOWN ? int(v2.size()) : size;
            call2(n,x,v1.nonConj().begin(),v2.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B)
        {
            const bool c1 = V1::vconj;
            if (n) do {
                *B++ = ZProd<false,c1>::quot(x , *A++);
            } while (--n);
        }
    };

#ifdef __SSE__
    // algo 21: single precision SSE: all real
    // SSE also has a _mm_rcp_ps command that takes a reciprocal.
    // However, it is only 12 bit accurate, which is not acceptable,
    // so unfortunately we can't use it.
    // Instead, we do the reciprocal as 1 / x with _mm_div_ps.
    template <int size, int ix, class T, class V1, class V2>
    struct ElemInvert_Helper<21,size,ix,T,V1,V2> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            const int n = size == UNKNOWN ? int(v2.size()) : size;
            call2(n,x,v1.nonConj().begin(),v2.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B)
        {
            const bool unit1 = V1::vstep == 1;
            const bool unit2 = V2::vstep == 1;

            if (unit2) {
                while (n && (((unsigned int)(B.getP()) & 0xf) != 0) ) {
                    *B++ = x / *A++;
                    --n;
                }
            } else if (unit1) {
                while (n && (((unsigned int)(A.getP()) & 0xf) != 0) ) {
                    *B++ = x / *A++;
                    --n;
                }
            }

            int n_4 = (n>>2);
            int nb = n-(n_4<<2);
            
            if (n_4) {
                IT1 A1 = A+1;
                IT1 A2 = A+2;
                IT1 A3 = A+3;

                IT2 B1 = B+1;
                IT2 B2 = B+2;
                IT2 B3 = B+3;

                __m128 xx = _mm_set1_ps(float(x));
                __m128 xA,xB;
                do {
                    // If !unit2, then we know that A is aligned, so we 
                    // can use load, otherwise use loadu.
                    Maybe2<!unit2,unit1>::sse_load(
                        xA,A.getP(),A1.getP(),A2.getP(),A3.getP());
                    A+=4; A1+=4; A2+=4; A3+=4;
                    xB = _mm_div_ps(xx,xA);
                    Maybe<unit2>::sse_store(
                        B.getP(),B1.getP(),B2.getP(),B3.getP(),xB);
                    B+=4; B1+=4; B2+=4; B3+=4;
                } while (--n_4);
            }

            if (nb) do { *B++ = x / *A++; } while (--nb);
        }
    };

    // algo 22: single precision SSE: x real v1 real v2 complex
    template <int size, int ix, class T, class V1, class V2>
    struct ElemInvert_Helper<22,size,ix,T,V1,V2> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            const int n = size == UNKNOWN ? int(v2.size()) : size;
            call2(n,x,v1.nonConj().begin(),v2.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B)
        {
            const bool unit1 = V1::vstep == 1;
            const bool unit2 = V2::vstep == 1;

            if (unit1) {
                while (n && (((unsigned int)(A.getP()) & 0xf) != 0) ) {
                    *B++ = x / *A++;
                    --n;
                }
            }

            int n_4 = (n>>2);
            int nb = n-(n_4<<2);
            
            if (n_4) {
                IT1 A1 = A+1;
                IT1 A2 = A+2;
                IT1 A3 = A+3;

                IT2 B1 = B+1;
                IT2 B2 = B+2;
                IT2 B3 = B+3;

                __m128 xx = _mm_set1_ps(float(x));
                __m128 xzero = _mm_set1_ps(0.F);
                __m128 xA,xB;
                do {
                    Maybe<unit1>::sse_load(
                        xA,A.getP(),A1.getP(),A2.getP(),A3.getP());
                    A+=4; A1+=4; A2+=4; A3+=4;
                    xB = _mm_div_ps(xx,xA);
                    Maybe<unit2>::sse_storeu(
                        B.getP(),B1.getP(),_mm_unpacklo_ps(xB,xzero));
                    Maybe<unit2>::sse_storeu(
                        B2.getP(),B3.getP(),_mm_unpackhi_ps(xB,xzero));
                    B+=4; B1+=4; B2+=4; B3+=4;
                } while (--n_4);
            }

            if (nb) do { *B++ = x / *A++; } while (--nb);
        }
    };

    // algo 23: single precision SSE: x real v1 complex v2 complex
    template <int size, int ix, class T, class V1, class V2>
    struct ElemInvert_Helper<23,size,ix,T,V1,V2> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            const int n = size == UNKNOWN ? int(v2.size()) : size;
            call2(n,x,v1.nonConj().begin(),v2.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B)
        {
            const bool unit1 = V1::vstep == 1;
            const bool unit2 = V2::vstep == 1;
            const bool c1 = V1::vconj;

            if (unit2) {
                while (n && (((unsigned int)(B.getP()) & 0xf) != 0) ) {
                    *B++ = ZProd<false,c1>::quot(x , *A++);
                    --n;
                }
            } else if (unit1) {
                while (n && (((unsigned int)(A.getP()) & 0xf) != 0) ) {
                    *B++ = ZProd<false,c1>::quot(x , *A++);
                    --n;
                }
            }

            int n_2 = (n>>1);
            int nb = n-(n_2<<1);
            
            if (n_2) {
                IT1 A1 = A+1;
                IT2 B1 = B+1;

                const float mx = Maybe<c1>::select( float(x) , -float(x) );
                // These look backwards, but order is from hi to lo values.
                __m128 xconj = _mm_set_ps(mx, float(x), mx, float(x));
                __m128 xA,xB;
                __m128 xAc, xnorm, x1, x2; // temp values
                do {
                    Maybe2<!unit2,unit1>::sse_load(xA,A.getP(),A1.getP());
                    A+=2; A1+=2;
                    xAc = _mm_mul_ps(xconj,xA); // conj(xA)
                    x1 = _mm_mul_ps(xA,xA);
                    x2 = _mm_shuffle_ps(x1,x1,_MM_SHUFFLE(2,3,0,1));
                    xnorm = _mm_add_ps(x1,x2); // = norm(xA)
                    xB = _mm_div_ps(xAc,xnorm);  // = x/xA
                    Maybe<unit2>::sse_store(B.getP(),B1.getP(),xB);
                    B+=2; B1+=2;
                } while (--n_2);
            }

            if (nb) *B = ZProd<false,c1>::quot(x , *A);
        }
    };

    // algo 24: single precision SSE: x complex v1 real v2 complex
    template <int size, int ix, class T, class V1, class V2>
    struct ElemInvert_Helper<24,size,ix,T,V1,V2> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            const int n = size == UNKNOWN ? int(v2.size()) : size;
            call2(n,x,v1.nonConj().begin(),v2.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B)
        {
            TMVStaticAssert(ix == 0);
            TMVStaticAssert((Traits2<T,std::complex<float> >::sametype));
            const bool unit1 = V1::vstep == 1;
            const bool unit2 = V2::vstep == 1;

            if (unit1) {
                while (n && (((unsigned int)(A.getP()) & 0xf) != 0) ) {
                    *B++ = x / *A++;
                    --n;
                }
            }

            int n_4 = (n>>2);
            int nb = n-(n_4<<2);
            
            if (n_4) {
                IT1 A1 = A+1;
                IT1 A2 = A+2;
                IT1 A3 = A+3;

                IT2 B1 = B+1;
                IT2 B2 = B+2;
                IT2 B3 = B+3;

                __m128 xone = _mm_set1_ps(1.F);
                __m128 xr = _mm_set1_ps(real(x.x));
                __m128 xi = _mm_set1_ps(imag(x.x));
                __m128 xA,xBr,xBi,xAinv;
                do {
                    Maybe<unit1>::sse_load(
                        xA,A.getP(),A1.getP(),A2.getP(),A3.getP());
                    A+=4; A1+=4; A2+=4; A3+=4;
                    xAinv = _mm_div_ps(xone,xA);
                    xBr = _mm_mul_ps(xr,xAinv);
                    xBi = _mm_mul_ps(xi,xAinv);
                    Maybe<unit2>::sse_storeu(
                        B.getP(),B1.getP(),_mm_unpacklo_ps(xBr,xBi));
                    Maybe<unit2>::sse_storeu(
                        B2.getP(),B3.getP(),_mm_unpackhi_ps(xBr,xBi));
                    B+=4; B1+=4; B2+=4; B3+=4;
                } while (--n_4);
            }

            if (nb) do { *B++ = x / *A++; } while (--nb);
        }
    };

    // algo 25: single precision SSE: all complex
    template <int size, int ix, class T, class V1, class V2>
    struct ElemInvert_Helper<25,size,ix,T,V1,V2> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            const int n = size == UNKNOWN ? int(v2.size()) : size;
            call2(n,x,v1.nonConj().begin(),v2.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B)
        {
            TMVStaticAssert(ix == 0);
            TMVStaticAssert((Traits2<T,std::complex<float> >::sametype));
            const bool unit1 = V1::vstep == 1;
            const bool unit2 = V2::vstep == 1;
            const bool c1 = V1::vconj;

            if (unit2) {
                while (n && (((unsigned int)(B.getP()) & 0xf) != 0) ) {
                    *B++ = ZProd<false,c1>::quot(x , *A++);
                    --n;
                }
            } else if (unit1) {
                while (n && (((unsigned int)(A.getP()) & 0xf) != 0) ) {
                    *B++ = ZProd<false,c1>::quot(x , *A++);
                    --n;
                }
            }

            int n_2 = (n>>1);
            int nb = n-(n_2<<1);
            
            if (n_2) {
                IT1 A1 = A+1;
                IT2 B1 = B+1;

                float xr = real(x.x);
                float mxr = Maybe<c1>::select(xr,-xr);
                float xi = imag(x.x);
                float mxi = Maybe<c1>::select(-xi,xi);
                // B = x * conj(A) / norm(A)
                // Br = xr * Ar + xi * Ai / norm
                // Bi = -xr * Ai + xi * Ar / norm
                __m128 xxr = _mm_set_ps(mxr, xr, mxr, xr);
                __m128 xxi = _mm_set_ps(xi, mxi, xi, mxi);
                __m128 xA,xB;
                __m128 xnorm, x0, x1, x2, x3, x4, x5; // temp values
                do {
                    Maybe2<!unit2,unit1>::sse_load(xA,A.getP(),A1.getP());
                    A+=2; A1+=2;
                    x0 = _mm_shuffle_ps(xA,xA,_MM_SHUFFLE(2,3,0,1));
                    x1 = _mm_mul_ps(xA,xA);
                    // It doesn't seem to matter which way we calculate x2.
                    // They take almost identical time.
                    //x2 = _mm_mul_ps(x0,x0);
                    x2 = _mm_shuffle_ps(x1,x1,_MM_SHUFFLE(2,3,0,1));
                    xnorm = _mm_add_ps(x1,x2); // = norm(xA)
                    x3 = _mm_mul_ps(xxr,xA);
                    x4 = _mm_mul_ps(xxi,x0);
                    x5 = _mm_add_ps(x3,x4);
                    xB = _mm_div_ps(x5,xnorm);  // = x/xA
                    Maybe<unit2>::sse_store(B.getP(),B1.getP(),xB);
                    B+=2; B1+=2;
                } while (--n_2);
            }

            if (nb) *B = ZProd<false,c1>::quot(x , *A);
        }
    };
#endif

#ifdef __SSE2__
    // algo 31: double precision SSE2: all real
    template <int size, int ix, class T, class V1, class V2>
    struct ElemInvert_Helper<31,size,ix,T,V1,V2> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            const int n = size == UNKNOWN ? int(v2.size()) : size;
            call2(n,x,v1.nonConj().begin(),v2.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B)
        {
            const bool unit1 = V1::vstep == 1;
            const bool unit2 = V2::vstep == 1;

            if (unit2) {
                while (n && (((unsigned int)(B.getP()) & 0xf) != 0) ) {
                    *B++ = x / *A++;
                    --n;
                }
            } else if (unit1) {
                while (n && (((unsigned int)(A.getP()) & 0xf) != 0) ) {
                    *B++ = x / *A++;
                    --n;
                }
            }

            int n_2 = (n>>1);
            int nb = n-(n_2<<1);
            
            if (n_2) {
                IT1 A1 = A+1;
                IT2 B1 = B+1;

                __m128d xx = _mm_set1_pd(double(x));
                __m128d xA,xB;
                do {
                    Maybe2<!unit2,unit1>::sse_load(xA,A.getP(),A1.getP());
                    A+=2; A1+=2;
                    xB = _mm_div_pd(xx,xA);
                    Maybe<unit2>::sse_store(B.getP(),B1.getP(),xB);
                    B+=2; B1+=2;
                } while (--n_2);
            }

            if (nb) *B = x / *A;
        }
    };

    // algo 32: double precision SSE2: x real v1 real v2 complex
    template <int size, int ix, class T, class V1, class V2>
    struct ElemInvert_Helper<32,size,ix,T,V1,V2> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            const int n = size == UNKNOWN ? int(v2.size()) : size;
            call2(n,x,v1.nonConj().begin(),v2.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B)
        {
            const bool unit1 = V1::vstep == 1;

            if (unit1) {
                while (n && (((unsigned int)(A.getP()) & 0xf) != 0) ) {
                    *B++ = x / *A++;
                    --n;
                }
            }

            int n_2 = (n>>1);
            int nb = n-(n_2<<1);
            
            if (n_2) {
                IT1 A1 = A+1;
                IT2 B1 = B+1;

                __m128d xx = _mm_set1_pd(double(x));
                __m128d xzero = _mm_set1_pd(0.);
                __m128d xA,xB;
                do {
                    Maybe<unit1>::sse_load(xA,A.getP(),A1.getP());
                    A+=2; A1+=2;
                    xB = _mm_div_pd(xx,xA);
                    Maybe<true>::sse_storeu(
                        B.getP(),_mm_unpacklo_pd(xB,xzero));
                    Maybe<true>::sse_storeu(
                        B1.getP(),_mm_unpackhi_pd(xB,xzero));
                    B+=2; B1+=2;
                } while (--n_2);
            }

            if (nb) *B = x / *A;
        }
    };

    // algo 33: double precision SSE2: x real v1 complex v2 complex
    template <int size, int ix, class T, class V1, class V2>
    struct ElemInvert_Helper<33,size,ix,T,V1,V2> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            const int n = size == UNKNOWN ? int(v2.size()) : size;
            call2(n,x,v1.nonConj().begin(),v2.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B)
        {
            const bool c1 = V1::vconj;
            if (n) {
                const double mone = Maybe<c1>::select( double(x) , -double(x) );
                // These look backwards, but order is from hi to lo values.
                __m128d xconj = _mm_set_pd(mone, double(x));
                __m128d xA,xB;
                __m128d xAc, xnorm, x1, x2; // temp values
                do {
                    Maybe<true>::sse_load(xA,A.getP()); ++A;
                    xAc = _mm_mul_pd(xconj,xA); // x*conj(xA)
                    x1 = _mm_mul_pd(xA,xA);
                    x2 = _mm_shuffle_pd(x1,x1,_MM_SHUFFLE2(0,1));
                    xnorm = _mm_add_pd(x1,x2); // = norm(xA)
                    xB = _mm_div_pd(xAc,xnorm);  // = x/xA
                    Maybe<true>::sse_store(B.getP(),xB); ++B;
                } while (--n);
            }
        }
    };

    // algo 34: double precision SSE2: x complex v1 real v2 complex
    template <int size, int ix, class T, class V1, class V2>
    struct ElemInvert_Helper<34,size,ix,T,V1,V2> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            const int n = size == UNKNOWN ? int(v2.size()) : size;
            call2(n,x,v1.nonConj().begin(),v2.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B)
        {
            TMVStaticAssert(ix == 0);
            TMVStaticAssert((Traits2<T,std::complex<double> >::sametype));
            const bool unit1 = V1::vstep == 1;

            if (unit1) {
                while (n && (((unsigned int)(A.getP()) & 0xf) != 0) ) {
                    *B++ = x / *A++;
                    --n;
                }
            }

            int n_2 = (n>>1);
            int nb = n-(n_2<<1);
            
            if (n_2) {
                IT1 A1 = A+1;
                IT2 B1 = B+1;

                __m128d xone = _mm_set1_pd(1.);
                __m128d xr = _mm_set1_pd(real(x.x));
                __m128d xi = _mm_set1_pd(imag(x.x));
                __m128d xA,xBr,xBi,xAinv;
                do {
                    Maybe<unit1>::sse_load(xA,A.getP(),A1.getP());
                    A+=2; A1+=2;
                    xAinv = _mm_div_pd(xone,xA);
                    xBr = _mm_mul_pd(xr,xAinv);
                    xBi = _mm_mul_pd(xi,xAinv);
                    Maybe<true>::sse_store(B.getP(),_mm_unpacklo_pd(xBr,xBi));
                    Maybe<true>::sse_store(B1.getP(),_mm_unpackhi_pd(xBr,xBi));
                    B+=2; B1+=2;
                } while (--n_2);
            }

            if (nb) *B = x / *A;
        }
    };

    // algo 35: double precision SSE2: all complex
    template <int size, int ix, class T, class V1, class V2>
    struct ElemInvert_Helper<35,size,ix,T,V1,V2> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        static void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            const int n = size == UNKNOWN ? int(v2.size()) : size;
            call2(n,x,v1.nonConj().begin(),v2.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B)
        {
            TMVStaticAssert(ix == 0);
            TMVStaticAssert((Traits2<T,std::complex<double> >::sametype));
            const bool c1 = V1::vconj;
            if (n) {
                double xr = real(x.x);
                double mxr = Maybe<c1>::select(xr,-xr);
                double xi = imag(x.x);
                double mxi = Maybe<c1>::select(-xi,xi);
                // B = x * conj(A) / norm(A)
                // Br = xr * Ar + xi * Ai / norm
                // Bi = -xr * Ai + xi * Ar / norm
                __m128d xxr = _mm_set_pd(mxr, xr);
                __m128d xxi = _mm_set_pd(xi, mxi);
                __m128d xA,xB;
                __m128d xnorm, x0, x1, x2, x3, x4, x5; // temp values
                do {
                    Maybe<true>::sse_load(xA,A.getP()); ++A;
                    x0 = _mm_shuffle_pd(xA,xA,_MM_SHUFFLE2(0,1));
                    x1 = _mm_mul_pd(xA,xA);
                    x2 = _mm_shuffle_pd(x1,x1,_MM_SHUFFLE2(0,1));
                    xnorm = _mm_add_pd(x1,x2); // = norm(xA)
                    x3 = _mm_mul_pd(xxr,xA);
                    x4 = _mm_mul_pd(xxi,x0);
                    x5 = _mm_add_pd(x3,x4);
                    xB = _mm_div_pd(x5,xnorm);  // = x/xA
                    Maybe<true>::sse_store(B.getP(),xB); ++B;
                } while (--n);
            }
        }
    };
#endif

    // algo -4: No branches or copies
    template <int size, int ix, class T, class V1, class V2>
    struct ElemInvert_Helper<-4,size,ix,T,V1,V2> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::iterator IT2;
        typedef typename V1::real_type RT1;
        typedef typename V2::real_type RT2;
        enum { allfloat = 
            Traits2<RT1,float>::sametype && Traits2<RT2,float>::sametype };
        enum { alldouble = 
            Traits2<RT1,double>::sametype && Traits2<RT2,double>::sametype };
        enum { xreal = Traits<T>::isreal };
        enum { xcomplex = Traits<T>::iscomplex };
        enum { v1real = V1::visreal };
        enum { v2real = V2::visreal };
        enum { v1complex = V1::viscomplex };
        enum { v2complex = V2::viscomplex };
        enum { unit = V1::vstep == 1 || V2::vstep == 1 };

        enum { algo = (
                size == 0 ? 0 : 
#if TMV_OPT >= 1
#ifdef __SSE__
                ( allfloat && xreal && v1real && v2real && unit ) ? 21 :
                ( allfloat && xreal && v1real && v2complex ) ? 22 :
                ( allfloat && xreal && v1complex && v2complex ) ? 23 :
                ( allfloat && xcomplex && v1real && v2complex ) ? 24 :
                ( allfloat && xcomplex && v1complex && v2complex ) ? 25 :
#endif
#ifdef __SSE2__
                ( alldouble && xreal && v1real && v2real && unit ) ? 31 :
                ( alldouble && xreal && v1real && v2complex ) ? 32 :
                ( alldouble && xreal && v1complex && v2complex ) ? 33 :
                ( alldouble && xcomplex && v1real && v2complex ) ? 34 :
                ( alldouble && xcomplex && v1complex && v2complex ) ? 35 :
#endif
#endif
                11 ) };
        static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            TMVStaticAssert(!V2::vconj);
#ifdef PRINTALGO_InvD
            if (algo != 11) {
                std::cout<<"InlineElemInvert:  x = "<<ix<<" "<<T(x)<<std::endl;
                std::cout<<"v1 = "<<TMV_Text(v1)<<std::endl;
                std::cout<<"v2 = "<<TMV_Text(v2)<<std::endl;
                std::cout<<"size = "<<size<<" = "<<v2.size()<<std::endl;
                std::cout<<"algo = "<<algo<<std::endl;
                std::cout<<"v1 = "<<v1<<std::endl;
            }
#endif
            ElemInvert_Helper<algo,size,ix,T,V1,V2>::call(x,v1,v2); 
#ifdef PRINTALGO_InvD
            std::cout<<"v2 = "<<v2<<std::endl;
#endif
        }
        static inline void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B)
        {
            TMVStaticAssert(!V2::vconj);
            ElemInvert_Helper<algo,size,ix,T,V1,V2>::call2(n,x,A,B); 
        }
    };

    // algo -3: Determine which algorithm to use
    template <int size, int ix, class T, class V1, class V2>
    struct ElemInvert_Helper<-3,size,ix,T,V1,V2> 
    {
        static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        { ElemInvert_Helper<-4,size,ix,T,V1,V2>::call(x,v1,v2); }
    };

    // algo 97: Conjugate
    template <int size, int ix, class T, class V1, class V2>
    struct ElemInvert_Helper<97,size,ix,T,V1,V2>
    {
        static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        { 
            typedef typename V1::const_conjugate_type V1c;
            typedef typename V2::conjugate_type V2c;
            V1c v1c = v1.conjugate();
            V2c v2c = v2.conjugate();
            ElemInvert_Helper<-2,size,ix,T,V1c,V2c>::call(x,v1c,v2c);
        }
    };

    // algo 98: Call inst
    template <int size, int ix, class T, class V1, class V2>
    struct ElemInvert_Helper<98,size,ix,T,V1,V2>
    {
        static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            typename V2::value_type xx(x);
            InstInvert(xx,v1.xView(),v2.xView()); 
        }
    };

    // algo -2: Check for inst
    template <int size, int ix, class T, class V1, class V2>
    struct ElemInvert_Helper<-2,size,ix,T,V1,V2>
    {
        static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
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
                Traits<T1>::isinst;
            const bool conj = V2::vconj;
            const int algo = 
                size == 0 ? 0 : 
                conj ? 97 :
                inst ? 98 :
                -4;
            ElemInvert_Helper<algo,size,ix,T,V1,V2>::call(x,v1,v2);
        }
    };

    // algo 99: Check for aliases
    template <int size, int ix, class T, class V1, class V2>
    struct ElemInvert_Helper<99,size,ix,T,V1,V2>
    {
        static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            if ( !SameStorage(v1,v2) || 
                 ExactSameStorage(v1,v2) || 
                 v1.step()*v2.step() < 0 || 
                 std::abs(v2.step()) < std::abs(v1.step()) ) { 
                // No aliasing (or no clobering)
                ElemInvert_Helper<-2,size,ix,T,V1,V2>::call(x,v1,v2);
            } else {
                // Need a temporary
                NoAliasElemInvert(x,v1.copy(),v2);
            }
        }
    };

    // algo -1: Check for aliases?
    template <int size, int ix, class T, class V1, class V2>
    struct ElemInvert_Helper<-1,size,ix,T,V1,V2>
    {
        static inline void call(const Scaling<ix,T>& x, const V1& v1, V2& v2)
        {
            const bool noclobber = VStepHelper<V1,V2>::noclobber;
            const bool checkalias =
                V1::vsize == UNKNOWN &&
                V2::vsize == UNKNOWN &&
                !noclobber;
            const int algo = 
                size == 0 ? 0 : 
                checkalias ? 99 : 
                -2;
            ElemInvert_Helper<algo,size,ix,T,V1,V2>::call(x,v1,v2);
        }
    };

    //
    // Element-wise v2(i) = 1/v1(i)
    //

    template <int ix, class T, class V1, class V2>
    inline void ElemInvert(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
        TMVAssert(v1.size() == v2.size());
        const int size = Sizes<V1::vsize,V2::vsize>::size;
        typedef typename V1::const_cview_type V1d;
        typedef typename V2::cview_type V2d;
        V1d v1d = v1.cView();
        V2d v2d = v2.cView();
        ElemInvert_Helper<-1,size,ix,T,V1d,V2d>::call(x,v1d,v2d);
    }

    template <int ix, class T, class V1, class V2>
    inline void NoAliasElemInvert(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
        TMVAssert(v1.size() == v2.size());
        const int size = Sizes<V1::vsize,V2::vsize>::size;
        typedef typename V1::const_cview_type V1d;
        typedef typename V2::cview_type V2d;
        V1d v1d = v1.cView();
        V2d v2d = v2.cView();
        ElemInvert_Helper<-2,size,ix,T,V1d,V2d>::call(x,v1d,v2d);
    }

    template <int ix, class T, class V1, class V2>
    inline void InlineElemInvert(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
        TMVAssert(v1.size() == v2.size());
        const int size = Sizes<V1::vsize,V2::vsize>::size;
        typedef typename V1::const_cview_type V1d;
        typedef typename V2::cview_type V2d;
        V1d v1d = v1.cView();
        V2d v2d = v2.cView();
        ElemInvert_Helper<-4,size,ix,T,V1d,V2d>::call(x,v1d,v2d);
    }

    template <int ix, class T, class V1, class V2>
    inline void AliasElemInvert(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1, BaseVector_Mutable<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
        TMVAssert(v1.size() == v2.size());
        const int size = Sizes<V1::vsize,V2::vsize>::size;
        typedef typename V1::const_cview_type V1d;
        typedef typename V2::cview_type V2d;
        V1d v1d = v1.cView();
        V2d v2d = v2.cView();
        ElemInvert_Helper<99,size,ix,T,V1d,V2d>::call(x,v1d,v2d);
    }

    //
    // D = D^-1
    //
    
    // TODO: I want these to check for singular matrices before inverting.
    // ElemInvert should just invert the elements without the check.
    template <int ix, class T, class M1, class M2>
    inline void Invert(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Diag_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::msize,M2::msize>::same));
        TMVAssert(m1.size() == m2.size());
        const int size = Sizes<M1::msize,M2::msize>::size;
        typedef typename M1::const_diag_type::const_cview_type M1d;
        typedef typename M2::diag_type::cview_type M2d;
        M1d m1d = m1.diag().cView();
        M2d m2d = m2.diag().cView();
        ElemInvert_Helper<-1,size,ix,T,M1d,M2d>::call(x,m1d,m2d);
    }

    template <int ix, class T, class M1, class M2>
    inline void NoAliasInvert(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Diag_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::msize,M2::msize>::same));
        TMVAssert(m1.size() == m2.size());
        const int size = Sizes<M1::msize,M2::msize>::size;
        typedef typename M1::const_diag_type::const_cview_type M1d;
        typedef typename M2::diag_type::cview_type M2d;
        M1d m1d = m1.diag().cView();
        M2d m2d = m2.diag().cView();
        ElemInvert_Helper<-2,size,ix,T,M1d,M2d>::call(x,m1d,m2d);
    }

    template <int ix, class T, class M1, class M2>
    inline void InlineInvert(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Diag_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::msize,M2::msize>::same));
        TMVAssert(m1.size() == m2.size());
        const int size = Sizes<M1::msize,M2::msize>::size;
        typedef typename M1::const_diag_type::const_cview_type M1d;
        typedef typename M2::diag_type::cview_type M2d;
        M1d m1d = m1.diag().cView();
        M2d m2d = m2.diag().cView();
        ElemInvert_Helper<-4,size,ix,T,M1d,M2d>::call(x,m1d,m2d);
    }

    template <int ix, class T, class M1, class M2>
    inline void AliasInvert(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Diag_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::msize,M2::msize>::same));
        TMVAssert(m1.size() == m2.size());
        const int size = Sizes<M1::msize,M2::msize>::size;
        typedef typename M1::const_diag_type::const_cview_type M1d;
        typedef typename M2::diag_type::cview_type M2d;
        M1d m1d = m1.diag().cView();
        M2d m2d = m2.diag().cView();
        ElemInvert_Helper<99,size,ix,T,M1d,M2d>::call(x,m1d,m2d);
    }


    // 
    // M = D^-1
    //
    
    template <int ix, class T, class M1, class M2>
    inline void Invert(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::msize,M2::mcolsize>::same));
        TMVStaticAssert((Sizes<M1::msize,M2::mrowsize>::same));
        TMVAssert(m1.size() == m2.colsize());
        TMVAssert(m1.size() == m2.rowsize());
        const int size2 = Sizes<M2::mcolsize,M2::mrowsize>::size;
        const int size = Sizes<M1::msize,size2>::size;
        typedef typename M1::const_diag_type::const_cview_type M1d;
        typedef typename M2::diag_type::cview_type M2d;
        M1d m1d = m1.diag().cView();
        M2d m2d = m2.diag().cView();
        ElemInvert_Helper<-1,size,ix,T,M1d,M2d>::call(x,m1d,m2d);
        m2.upperTri().offDiag().setZero();
        m2.lowerTri().offDiag().setZero();
    }

    template <int ix, class T, class M1, class M2>
    inline void NoAliasInvert(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::msize,M2::mcolsize>::same));
        TMVStaticAssert((Sizes<M1::msize,M2::mrowsize>::same));
        TMVAssert(m1.size() == m2.colsize());
        TMVAssert(m1.size() == m2.rowsize());
        const int size2 = Sizes<M2::mcolsize,M2::mrowsize>::size;
        const int size = Sizes<M1::msize,size2>::size;
        typedef typename M1::const_diag_type::const_cview_type M1d;
        typedef typename M2::diag_type::cview_type M2d;
        M1d m1d = m1.diag().cView();
        M2d m2d = m2.diag().cView();
        m2.setZero();
        ElemInvert_Helper<-2,size,ix,T,M1d,M2d>::call(x,m1d,m2d);
    }

    template <int ix, class T, class M1, class M2>
    inline void AliasInvert(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::msize,M2::mcolsize>::same));
        TMVStaticAssert((Sizes<M1::msize,M2::mrowsize>::same));
        TMVAssert(m1.size() == m2.colsize());
        TMVAssert(m1.size() == m2.rowsize());
        const int size2 = Sizes<M2::mcolsize,M2::mrowsize>::size;
        const int size = Sizes<M1::msize,size2>::size;
        typedef typename M1::const_diag_type::const_cview_type M1d;
        typedef typename M2::diag_type::cview_type M2d;
        M1d m1d = m1.diag().cView();
        M2d m2d = m2.diag().cView();
        ElemInvert_Helper<99,size,ix,T,M1d,M2d>::call(x,m1d,m2d);
        m2.upperTri().offDiag().setZero();
        m2.lowerTri().offDiag().setZero();
    }



} // namespace tmv

#endif 
