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


#ifndef TMV_DivVD_H
#define TMV_DivVD_H

#include "TMV_BaseMatrix_Diag.h"
#include "TMV_BaseVector.h"

namespace tmv {

    // Defined below:
    template <int ix, class T, class V1, class M2, class V3>
    inline void LDiv(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class M2, class V3>
    inline void NoAliasLDiv(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class M2, class V3>
    inline void InlineLDiv(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class M2, class V3>
    inline void AliasLDiv(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3);

    template <int ix, class T, class V1, class V2, class V3>
    inline void ElemDivVV(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class V2, class V3>
    inline void NoAliasElemDivVV(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class V2, class V3>
    inline void InlineElemDivVV(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3);
    template <int ix, class T, class V1, class V2, class V3>
    inline void AliasElemDivVV(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3);


    // Defined in TMV_DivVD.cpp
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstLDiv(
        const T3 x,
        const ConstVectorView<T1,UNKNOWN,C1>& v1,
        const ConstDiagMatrixView<T2,UNKNOWN,C2>& v2, VectorView<T3> v3);
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstElemDivVV(
        const T3 x,
        const ConstVectorView<T1,UNKNOWN,C1>& v1,
        const ConstVectorView<T2,UNKNOWN,C2>& v2, VectorView<T3> v3);

    //
    // Element-wise division:
    // v3(i) = x * v1(i) / v2(i)
    //

    template <int algo, int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper;

    // algo 0: size == 0, nothing to do
    template <int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<0,0,ix,T,V1,V2,V3> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(const Scaling<ix,T>&, const V1&, const V2&, V3&) {}
        static void call2(int, const Scaling<ix,T>&, IT1, IT2, IT3) {}
    };

    // algo 11: simple for loop
    template <int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<1,size,ix,T,V1,V2,V3> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = size == UNKNOWN ? int(v3.size()) : size;
            call2(n,x,v1.nonConj().begin(),v2.nonConj().begin(),v3.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B, IT3 C)
        {
            const bool c1 = V1::vconj;
            const bool c2 = V2::vconj;
            if (n) do {
                *C++ = ZProd<false,false>::prod(
                    x, ZProd<c1,c2>::quot(*A++,*B++));
            } while (--n);
        }
    };  

#ifdef __SSE__
    // algo 21: single precision SSE: all real
    template <int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<21,size,ix,T,V1,V2,V3> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = size == UNKNOWN ? int(v2.size()) : size;
            call2(n,x,v1.nonConj().begin(),v2.nonConj().begin(),v3.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B, IT3 C)
        {
            const bool unit1 = V1::vstep == 1;
            const bool unit2 = V2::vstep == 1;
            const bool unit3 = V2::vstep == 1;

            if (unit3) {
                while (n && (((unsigned int)(C.getP()) & 0xf) != 0) ) {
                    *C++ = x * *A++ / *B++;
                    --n;
                }
            } else if (unit1) {
                while (n && (((unsigned int)(A.getP()) & 0xf) != 0) ) {
                    *C++ = x * *A++ / *B++;
                    --n;
                }
            } else if (unit2) {
                while (n && (((unsigned int)(B.getP()) & 0xf) != 0) ) {
                    *C++ = x * *A++ / *B++;
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

                IT3 C1 = C+1;
                IT3 C2 = C+2;
                IT3 C3 = C+3;

                __m128 xx = _mm_set1_ps(float(x));
                __m128 xA,xB,xC,x1;
                do {
                    Maybe2<!unit3,unit1>::sse_load(
                        xA,A.getP(),A1.getP(),A2.getP(),A3.getP());
                    A+=4; A1+=4; A2+=4; A3+=4;
                    Maybe2<!(unit3||unit1),unit2>::sse_load(
                        xB,B.getP(),B1.getP(),B2.getP(),B3.getP());
                    B+=4; B1+=4; B2+=4; B3+=4;
                    xC = _mm_div_ps(xA,xB);
                    Maybe<ix!=1>::sse_mult(xx,xC);
                    Maybe<unit2>::sse_store(
                        C.getP(),C1.getP(),C2.getP(),C3.getP(),xC);
                    C+=4; C1+=4; C2+=4; C3+=4;
                } while (--n_4);
            }

            if (nb) do { *C++ = x * *A++ / *B++; } while (--nb);
        }
    };

    // algo 22: single precision SSE: x real v1 real v2 complex
    template <int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<22,size,ix,T,V1,V2,V3> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = size == UNKNOWN ? int(v2.size()) : size;
            call2(n,x,v1.nonConj().begin(),v2.nonConj().begin(),v3.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B, IT3 C)
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
    template <int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<23,size,ix,T,V1,V2,V3> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = size == UNKNOWN ? int(v2.size()) : size;
            call2(n,x,v1.nonConj().begin(),v2.nonConj().begin(),v3.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B, IT3 C)
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
    template <int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<24,size,ix,T,V1,V2,V3> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = size == UNKNOWN ? int(v2.size()) : size;
            call2(n,x,v1.nonConj().begin(),v2.nonConj().begin(),v3.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B, IT3 C)
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
    template <int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<25,size,ix,T,V1,V2,V3> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = size == UNKNOWN ? int(v2.size()) : size;
            call2(n,x,v1.nonConj().begin(),v2.nonConj().begin(),v3.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B, IT3 C)
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
    template <int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<31,size,ix,T,V1,V2,V3> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = size == UNKNOWN ? int(v2.size()) : size;
            call2(n,x,v1.nonConj().begin(),v2.nonConj().begin(),v3.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B, IT3 C)
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
    template <int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<32,size,ix,T,V1,V2,V3> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = size == UNKNOWN ? int(v2.size()) : size;
            call2(n,x,v1.nonConj().begin(),v2.nonConj().begin(),v3.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B, IT3 C)
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
    template <int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<33,size,ix,T,V1,V2,V3> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = size == UNKNOWN ? int(v2.size()) : size;
            call2(n,x,v1.nonConj().begin(),v2.nonConj().begin(),v3.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B, IT3 C)
        {
            const bool c1 = V1::vconj;
            if (n) {
                const double mone = Maybe<c1>::select( double(x) , -double(x) );
                // These look backwards, but order is from hi to lo values.
                __m128d xconj = _mm_set_pd(mone, double(x));
                __m128d xA,xB;
                __m128d xAc, xnorm, x1, x2; // temp values
                if (((unsigned int)(A.getP()) & 0xf) == 0 &&
                    ((unsigned int)(B.getP()) & 0xf) == 0 ) {
                    do {
                        Maybe<true>::sse_load(xA,A.getP()); ++A;
                        xAc = _mm_mul_pd(xconj,xA); // x*conj(xA)
                        x1 = _mm_mul_pd(xA,xA);
                        x2 = _mm_shuffle_pd(x1,x1,_MM_SHUFFLE2(0,1));
                        xnorm = _mm_add_pd(x1,x2); // = norm(xA)
                        xB = _mm_div_pd(xAc,xnorm);  // = x/xA
                        Maybe<true>::sse_store(B.getP(),xB); ++B;
                    } while (--n);
                } else {
                    do {
                        Maybe<true>::sse_loadu(xA,A.getP()); ++A;
                        xAc = _mm_mul_pd(xconj,xA); // x*conj(xA)
                        x1 = _mm_mul_pd(xA,xA);
                        x2 = _mm_shuffle_pd(x1,x1,_MM_SHUFFLE2(0,1));
                        xnorm = _mm_add_pd(x1,x2); // = norm(xA)
                        xB = _mm_div_pd(xAc,xnorm);  // = x/xA
                        Maybe<true>::sse_storeu(B.getP(),xB); ++B;
                    } while (--n);
                }
            }
        }
    };

    // algo 34: double precision SSE2: x complex v1 real v2 complex
    template <int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<34,size,ix,T,V1,V2,V3> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = size == UNKNOWN ? int(v2.size()) : size;
            call2(n,x,v1.nonConj().begin(),v2.nonConj().begin(),v3.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B, IT3 C)
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
    template <int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<35,size,ix,T,V1,V2,V3> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const int n = size == UNKNOWN ? int(v2.size()) : size;
            call2(n,x,v1.nonConj().begin(),v2.nonConj().begin(),v3.begin());
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B, IT3 C)
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
    template <int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<-4,size,ix,T,V1,V2,V3> 
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename V3::iterator IT3;
        enum { algo = (
                1 ) };
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            TMVStaticAssert(!V3::vconj);
            ElemDivVV_Helper<algo,size,ix,T,V1,V2,V3>::call(x,v1,v2,v3); 
        }
        static void call2(int n, const Scaling<ix,T>& x, IT1 A, IT2 B, IT3 C)
        {
            TMVStaticAssert(!V3::vconj);
            ElemDivVV_Helper<algo,size,ix,T,V1,V2,V3>::call2(n,x,A,B,C);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<-3,size,ix,T,V1,V2,V3> 
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        { ElemDivVV_Helper<-4,size,ix,T,V1,V2,V3>::call(x,v1,v2,v3); }
    };

    // algo 97: Conjugate
    template <int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<97,size,ix,T,V1,V2,V3> 
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        { 
            typedef typename V1::const_conjugate_type V1c;
            typedef typename V2::const_conjugate_type V2c;
            typedef typename V3::conjugate_type V3c;
            V1c v1c = v1.conjugate();
            V2c v2c = v2.conjugate();
            V3c v3c = v3.conjugate();
            ElemDivVV_Helper<-2,size,ix,T,V1c,V2c,V3c>::call(
                TMV_CONJ(x),v1c,v2c,v3c);
        }
    };

    // algo 98: Call inst
    template <int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<98,size,ix,T,V1,V2,V3> 
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        { 
            typename V3::value_type xx(x);
            InstElemDivVV(xx,v1.xView(),v2.xView(),v3.xView()); 
        }
    };

    // algo -2: Check for inst
    template <int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<-2,size,ix,T,V1,V2,V3> 
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            typedef typename V3::value_type T3;
            const bool inst =
                V1::vsize == UNKNOWN &&
                V2::vsize == UNKNOWN &&
                V3::vsize == UNKNOWN &&
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const int algo =
                V3::vconj ? 97 : 
                inst ? 98 : 
                -4;
            ElemDivVV_Helper<algo,size,ix,T,V1,V2,V3>::call(x,v1,v2,v3); 
        }
    };

    // algo 99: Check for aliases
    template <int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<99,size,ix,T,V1,V2,V3> 
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const bool noclobber1 =
                !SameStorage(v1,v3) ||
                ExactSameStorage(v1,v3) || 
                v1.step()*v3.step() < 0 || 
                std::abs(v3.step()) < std::abs(v1.step());
            const bool noclobber2 =
                !SameStorage(v2,v3) ||
                ExactSameStorage(v2,v3) || 
                v2.step()*v3.step() < 0 || 
                std::abs(v3.step()) < std::abs(v2.step());
            if (noclobber1) {
                if (noclobber2) {
                    // No aliasing (or no clobering)
                    ElemDivVV_Helper<-2,size,ix,T,V1,V2,V3>::call(x,v1,v2,v3); 
                } else { 
                    // Need a temporary for v2
                    NoAliasElemDivVV(x,v1,v2.copy(),v3);
                }
            } else {
                if (noclobber2) {
                    // Need a temporary for v1
                    NoAliasElemDivVV(v1.copy(),v2,v3);
                } else {
                    // Need a temporary for v3
                    typename V3::copy_type v3c(v3.size());
                    NoAliasElemDivVV(x,v1,v2,v3c);
                    NoAliasCopy(v3c,v3);
                }
            }
        }
    };

    // algo -1: Check for aliases?
    template <int size, int ix, class T, class V1, class V2, class V3>
    struct ElemDivVV_Helper<-1,size,ix,T,V1,V2,V3> 
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const V2& v2, V3& v3)
        {
            const bool noclobber =
                VStepHelper<V1,V3>::noclobber &&
                VStepHelper<V2,V3>::noclobber;
            const bool checkalias =
                V1::vsize == UNKNOWN &&
                V2::vsize == UNKNOWN &&
                V3::vsize == UNKNOWN &&
                !noclobber;
            const int algo =
                checkalias ? 99 : 
                -2;
            ElemDivVV_Helper<algo,size,ix,T,V1,V2,V3>::call(x,v1,v2,v3); 
        }
    };

    template <int ix, class T, class V1, class V2, class V3>
    inline void ElemDivVV(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
        TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v1.size() == v3.size());
        const int s12 = Sizes<V1::vsize,V2::vsize>::size;
        const int size = Sizes<s12,V3::vsize>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2d;
        typedef typename V3::cview_type V3v;
        V1v v1v = v1.cView();
        V2d v2v = v2.cView();
        V3v v3v = v3.cView();
        ElemDivVV_Helper<-1,size,ix,T,V1v,V2d,V3v>::call(x,v1v,v2v,v3v);
    }

    template <int ix, class T, class V1, class V2, class V3>
    inline void NoAliasElemDivVV(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
        TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v1.size() == v3.size());
        const int s12 = Sizes<V1::vsize,V2::vsize>::size;
        const int size = Sizes<s12,V3::vsize>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2d;
        typedef typename V3::cview_type V3v;
        V1v v1v = v1.cView();
        V2d v2v = v2.cView();
        V3v v3v = v3.cView();
        ElemDivVV_Helper<-2,size,ix,T,V1v,V2d,V3v>::call(x,v1v,v2v,v3v);
    }

    template <int ix, class T, class V1, class V2, class V3>
    inline void InlineElemDivVV(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
        TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v1.size() == v3.size());
        const int s12 = Sizes<V1::vsize,V2::vsize>::size;
        const int size = Sizes<s12,V3::vsize>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2d;
        typedef typename V3::cview_type V3v;
        V1v v1v = v1.cView();
        V2d v2v = v2.cView();
        V3v v3v = v3.cView();
        ElemDivVV_Helper<-4,size,ix,T,V1v,V2d,V3v>::call(x,v1v,v2v,v3v);
    }

    template <int ix, class T, class V1, class V2, class V3>
    inline void AliasElemDivVV(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseVector_Calc<V2>& v2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::vsize,V2::vsize>::same));
        TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
        TMVAssert(v1.size() == v2.size());
        TMVAssert(v1.size() == v3.size());
        const int s12 = Sizes<V1::vsize,V2::vsize>::size;
        const int size = Sizes<s12,V3::vsize>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2d;
        typedef typename V3::cview_type V3v;
        V1v v1v = v1.cView();
        V2d v2v = v2.cView();
        V3v v3v = v3.cView();
        ElemDivVV_Helper<99,size,ix,T,V1v,V2d,V3v>::call(x,v1v,v2v,v3v);
    }

    // TODO: Check for singular DiagMatrix
    template <int ix, class T, class V1, class M2, class V3>
    inline void LDiv(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::vsize,M2::msize>::same));
        TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
        TMVAssert(v1.size() == m2.size());
        TMVAssert(v1.size() == v3.size());
        const int s12 = Sizes<V1::vsize,M2::msize>::size;
        const int size = Sizes<s12,V3::vsize>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::const_diag_type::const_cview_type M2d;
        typedef typename V3::cview_type V3v;
        V1v v1v = v1.cView();
        M2d m2d = m2.diag().cView();
        V3v v3v = v3.cView();
        ElemDivVV_Helper<-1,size,ix,T,V1v,M2d,V3v>::call(x,v1v,m2d,v3v);
    }

    template <int ix, class T, class V1, class M2, class V3>
    inline void NoAliasLDiv(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::vsize,M2::msize>::same));
        TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
        TMVAssert(v1.size() == m2.size());
        TMVAssert(v1.size() == v3.size());
        const int s12 = Sizes<V1::vsize,M2::msize>::size;
        const int size = Sizes<s12,V3::vsize>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::const_diag_type::const_cview_type M2d;
        typedef typename V3::cview_type V3v;
        V1v v1v = v1.cView();
        M2d m2d = m2.diag().cView();
        V3v v3v = v3.cView();
        ElemDivVV_Helper<-2,size,ix,T,V1v,M2d,V3v>::call(x,v1v,m2d,v3v);
    }

    template <int ix, class T, class V1, class M2, class V3>
    inline void InlineLDiv(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::vsize,M2::msize>::same));
        TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
        TMVAssert(v1.size() == m2.size());
        TMVAssert(v1.size() == v3.size());
        const int s12 = Sizes<V1::vsize,M2::msize>::size;
        const int size = Sizes<s12,V3::vsize>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::const_diag_type::const_cview_type M2d;
        typedef typename V3::cview_type V3v;
        V1v v1v = v1.cView();
        M2d m2d = m2.diag().cView();
        V3v v3v = v3.cView();
        ElemDivVV_Helper<-4,size,ix,T,V1v,M2d,V3v>::call(x,v1v,m2d,v3v);
    }

    template <int ix, class T, class V1, class M2, class V3>
    inline void AliasLDiv(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::vsize,M2::msize>::same));
        TMVStaticAssert((Sizes<V1::vsize,V3::vsize>::same));
        TMVAssert(v1.size() == m2.size());
        TMVAssert(v1.size() == v3.size());
        const int s12 = Sizes<V1::vsize,M2::msize>::size;
        const int size = Sizes<s12,V3::vsize>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::const_diag_type::const_cview_type M2d;
        typedef typename V3::cview_type V3v;
        V1v v1v = v1.cView();
        M2d m2d = m2.diag().cView();
        V3v v3v = v3.cView();
        ElemDivVV_Helper<99,size,ix,T,V1v,M2d,V3v>::call(x,v1v,m2d,v3v);
    }

    //
    // v1 /= m2
    //

    template <class V1, class M2>
    inline void LDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Diag<M2>& m2)
    { LDiv(v1,m2,v1); }
    template <class V1, class M2>
    inline void NoAliasLDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Diag<M2>& m2)
    { NoAliasLDiv(v1,m2,v1); }
    template <class V1, class M2>
    inline void AliasLDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Diag<M2>& m2)
    { AliasLDiv(v1,m2,v1); }

    //
    // v3 = v1 % m2
    //

    template <int ix, class T, class V1, class M2, class V3>
    inline void RDiv(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3)
    { LDiv(x,v1,m2,v3); }
    template <int ix, class T, class V1, class M2, class V3>
    inline void NoAliasRDiv(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3)
    { NoAliasLDiv(x,v1,m2,v3); }
    template <int ix, class T, class V1, class M2, class V3>
    inline void AliasRDiv(
        const Scaling<ix,T>& x,
        const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Diag<M2>& m2, BaseVector_Mutable<V3>& v3)
    { AliasLDiv(x,v1,m2,v3); }

    //
    // v1 %= m2
    //

    template <class V1, class M2>
    inline void RDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Diag<M2>& m2)
    { LDiv(v1,m2,v1); }
    template <class V1, class M2>
    inline void NoAliasRDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Diag<M2>& m2)
    { NoAliasLDiv(v1,m2,v1); }
    template <class V1, class M2>
    inline void AliasRDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Diag<M2>& m2)
    { AliasLDiv(v1,m2,v1); }

} // namespace tmv

#endif 
