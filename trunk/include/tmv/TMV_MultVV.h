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


#ifndef TMV_MultVV_H
#define TMV_MultVV_H

#include "TMV_Scaling.h"
#include "TMV_BaseVector.h"

namespace tmv {

    // Defined in TMV_MultVV.cpp
    template <class T1, bool C1, class T2>
    T2 InstMultVV(
        const ConstVectorView<T1,UNKNOWN,C1>& v1,
        const ConstVectorView<T2>& v2);


    //
    // Vector * Vector
    //

#if TMV_OPT >= 2
#define TMV_VV_RECURSE
#endif
    const int TMV_MultVV_RecurseSize = 16*1024;

    template <int algo, int s, class V1, class V2>
    struct MultVV_Helper;

    // algo 0: s=0, return 0
    template <class V1, class V2>
    struct MultVV_Helper<0,0,V1,V2>
    {
        typedef typename ProdType<V1,V2>::type PT;
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        static PT call(const V1& , const V2& )
        { return PT(0); }
        static PT call2(int n, IT1 it1, IT2 it2)
        { return PT(0); }
    };

    // algo 1: v1 is complex, v2 is real: swap v1,v2
    template <int s, class V1, class V2>
    struct MultVV_Helper<1,s,V1,V2>
    {
        typedef typename ProdType<V1,V2>::type PT;
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        static PT call(const V1& v1, const V2& v2)
        { return MultVV_Helper<-4,s,V2,V1>::call(v2,v1); }
        static PT call2(int n, IT1 it1, IT2 it2)
        { return MultVV_Helper<-4,s,V2,V1>::call2(n,it2,it1); }
    };

    // algo 2: v2 is conj: conjugate result
    template <int s, class V1, class V2>
    struct MultVV_Helper<2,s,V1,V2>
    {
        typedef typename ProdType<V1,V2>::type PT;
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        static PT call(const V1& v1, const V2& v2)
        {
            typedef typename V1::const_conjugate_type V1c;
            typedef typename V2::const_conjugate_type V2c;
            V1c v1c = v1.conjugate();
            V2c v2c = v2.conjugate();
            return TMV_CONJ(MultVV_Helper<-4,s,V1c,V2c>::call(v1c,v2c));
        }
        static PT call2(int n, IT1 it1, IT2 it2)
        {
            typedef typename V1::const_conjugate_type V1c;
            typedef typename V2::const_conjugate_type V2c;
            return TMV_CONJ(MultVV_Helper<-4,s,V1c,V2c>::call2(n,it1,it2));
        }
    };

    // algo 11: simple for loop
    template <int s, class V1, class V2>
    struct MultVV_Helper<11,s,V1,V2>
    {
        typedef typename ProdType<V1,V2>::type PT;
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        static PT call(const V1& v1, const V2& v2)
        {
            const int n = s == UNKNOWN ? int(v1.size()) : s;
            if (n == 0) return PT(0);
            else {
                PT sum(0);
                for(int i=0;i<n;++i) 
                    sum += ZProd<false,false>::prod(v1.cref(i) , v2.cref(i));
                return sum;
            }
        }
        static PT call2(int n, IT1 it1, IT2 it2)
        {
            const bool c1 = V1::_conj;
            const bool c2 = V2::_conj;
            PT sum(0);
            if (n) do {
                sum += ZProd<c1,c2>::prod(*it1++ , *it2++); 
            } while (--n);
            return sum;
        }
    };

    // algo 12: 2 at once
    template <int s, class V1, class V2>
    struct MultVV_Helper<12,s,V1,V2>
    {
        typedef typename ProdType<V1,V2>::type PT;
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        static PT call(const V1& v1, const V2& v2)
        {
            const int n = s == UNKNOWN ? int(v1.size()) : s;
            return call2(n,v1.begin().nonConj(),v2.begin().nonConj());
        }
        static PT call2(const int n, IT1 it1, IT2 it2)
        {
            PT sum0(0), sum1(0);
            int n_2 = (n>>1);
            const int nb = n-(n_2<<1);
            const bool c1 = V1::_conj;
            const bool c2 = V2::_conj;
            if (n_2) {
                do {
                    sum0 += ZProd<c1,c2>::prod(it1[0] , it2[0]);
                    sum1 += ZProd<c1,c2>::prod(it1[1] , it2[1]);
                    it1 += 2;
                    it2 += 2;
                } while (--n_2);
                sum0 += sum1;
            }
            if (nb) {
                sum0 += ZProd<c1,c2>::prod((*it1) , (*it2));
            }
            return sum0;
        }
    };

    // algo 13: 4 at once
    template <int s, class V1, class V2>
    struct MultVV_Helper<13,s,V1,V2>
    {
        typedef typename ProdType<V1,V2>::type PT;
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        static PT call(const V1& v1, const V2& v2)
        {
            const int n = s == UNKNOWN ? int(v1.size()) : s;
            return call2(n,v1.begin().nonConj(),v2.begin().nonConj());
        }
        static PT call2(const int n, IT1 it1, IT2 it2)
        {
            PT sum0(0), sum1(0);
            int n_4 = (n>>2);
            int nb = n-(n_4<<2);
            const bool c1 = V1::_conj;
            const bool c2 = V2::_conj;
            if (n_4) {
                do {
                    sum0 += ZProd<c1,c2>::prod(it1[0] , it2[0]);
                    sum1 += ZProd<c1,c2>::prod(it1[1] , it2[1]);
                    sum0 += ZProd<c1,c2>::prod(it1[2] , it2[2]);
                    sum1 += ZProd<c1,c2>::prod(it1[3] , it2[3]);
                    it1 += 4;
                    it2 += 4;
                } while (--n_4);
                sum0 += sum1;
            }
            if (nb) do {
                sum0 += ZProd<c1,c2>::prod((*it1++) , (*it2++));
            } while (--nb);
            return sum0;
        }
    };

    // algo 14: 8 at once
    template <int s, class V1, class V2>
    struct MultVV_Helper<14,s,V1,V2>
    {
        typedef typename ProdType<V1,V2>::type PT;
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        static PT call(const V1& v1, const V2& v2)
        {
            const int n = s == UNKNOWN ? int(v1.size()) : s;
            return call2(n,v1.begin().nonConj(),v2.begin().nonConj());
        }
        static PT call2(const int n, IT1 it1, IT2 it2)
        {
            PT sum0(0), sum1(0);
            int n_8 = (n>>3);
            int nb = n-(n_8<<3);
            const bool c1 = V1::_conj;
            const bool c2 = V2::_conj;
            if (n_8) {
                do {
                    sum0 += ZProd<c1,c2>::prod(it1[0] , it2[0]);
                    sum1 += ZProd<c1,c2>::prod(it1[1] , it2[1]);
                    sum0 += ZProd<c1,c2>::prod(it1[2] , it2[2]);
                    sum1 += ZProd<c1,c2>::prod(it1[3] , it2[3]);
                    sum0 += ZProd<c1,c2>::prod(it1[4] , it2[4]);
                    sum1 += ZProd<c1,c2>::prod(it1[5] , it2[5]);
                    sum0 += ZProd<c1,c2>::prod(it1[6] , it2[6]);
                    sum1 += ZProd<c1,c2>::prod(it1[7] , it2[7]);
                    it1 += 8; it2 += 8;
                } while (--n_8);
                sum0 += sum1;
            }
            if (nb) do {
                sum0 += ZProd<c1,c2>::prod((*it1++) , (*it2++));
            } while (--nb);
            return sum0;
        }
    };

    // algo 15: fully unroll
    template <int s, class V1, class V2>
    struct MultVV_Helper<15,s,V1,V2>
    {
        typedef typename ProdType<V1,V2>::type PT;
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        template <int I, int N>
        struct Unroller
        {
            static PT unroll(const V1& v1, const V2& v2)
            {
                return (
                    Unroller<I,N/2>::unroll(v1,v2) +
                    Unroller<I+N/2,N-N/2>::unroll(v1,v2));
            }
            static PT unroll2(const IT1& it1, const IT2& it2)
            {
                return (
                    Unroller<I,N/2>::unroll2(it1,it2) +
                    Unroller<I+N/2,N-N/2>::unroll2(it1,it2));
            }
        };
        template <int I>
        struct Unroller<I,1>
        {
            static PT unroll(const V1& v1, const V2& v2)
            { return ZProd<false,false>::prod(v1.cref(I) , v2.cref(I)); }
            static PT unroll2(const IT1& it1, const IT2& it2)
            {
                const bool c1 = V1::_conj;
                const bool c2 = V2::_conj;
                return ZProd<c1,c2>::prod(it1[I] , it2[I]); 
            }
        };
        template <int I>
        struct Unroller<I,0>
        {
            static PT unroll(const V1& v1, const V2& v2)
            { return PT(0); }
            static PT unroll2(const IT1& it1, const IT2& it2)
            { return PT(0); }
        };

        static PT call(const V1& v1, const V2& v2)
        { return Unroller<0,s>::unroll(v1,v2); }
        static PT call2(const int n, const IT1& it1, const IT2& it2)
        { return Unroller<0,s>::unroll2(it1,it2); }
    };

#ifdef __SSE__
    // algo 21: single precision SSE: all real
    template <int s, class V1, class V2>
    struct MultVV_Helper<21,s,V1,V2>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        static float call(const V1& v1, const V2& v2)
        {
            const int n = s == UNKNOWN ? int(v1.size()) : s;
            return call2(n,v1.begin().nonConj(),v2.begin().nonConj());
        }
        static float call2(int n, IT1 A, IT2 B)
        {
            const bool unit1 = V1::_step == 1;
            const bool unit2 = V2::_step == 1;

            float sum0(0), sum1;

            if (unit2) {
                while (n && !TMV_Aligned(B.get()) ) {
                    sum0 += *A++ * *B++;
                    --n;
                }
            } else if (unit1) {
                while (n && !TMV_Aligned(A.get()) ) {
                    sum0 += *A++ * *B++;
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

                union { __m128 xm; float xf[4]; } xsum;
                xsum.xm = _mm_set1_ps(0.F);
                __m128 xA,xB,x0;
                do {
                    Maybe2<!unit2,unit1>::sse_load(
                        xA,A.get(),A1.get(),A2.get(),A3.get());
                    A+=4; A1+=4; A2+=4; A3+=4;
                    Maybe<unit2>::sse_load(
                        xB,B.get(),B1.get(),B2.get(),B3.get());
                    B+=4; B1+=4; B2+=4; B3+=4;
                    x0 = _mm_mul_ps(xA,xB);
                    xsum.xm = _mm_add_ps(xsum.xm,x0);
                } while (--n_4);
                sum1 = xsum.xf[0] + xsum.xf[1] + xsum.xf[2] + xsum.xf[3];
            } else { sum1 = 0.F; }

            if (nb) do {
                sum0 += *A++ * *B++;
            } while (--nb);
            return sum0 + sum1;
        }
    };

    // algo 22: single precision SSE: v1 real v2 complex
    template <int s, class V1, class V2>
    struct MultVV_Helper<22,s,V1,V2>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        static std::complex<float> call(const V1& v1, const V2& v2)
        {
            const int n = s == UNKNOWN ? int(v1.size()) : s;
            return call2(n,v1.begin().nonConj(),v2.begin().nonConj());
        }
        static std::complex<float> call2(int n, IT1 A, IT2 B)
        {
            TMVStaticAssert(!V2::_conj);
            const bool unit1 = V1::_step == 1;
            const bool unit2 = V2::_step == 1;

            std::complex<float> sum0(0), sum1;

            if (unit2) {
                while (n && !TMV_Aligned(B.get()) ) {
                    sum0 += ZProd<false,false>::prod(*A++ , *B++);
                    --n;
                }
            } else if (unit1) {
                while (n && !TMV_Aligned(A.get()) ) {
                    sum0 += ZProd<false,false>::prod(*A++ , *B++);
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

                union { __m128 xm; float xf[4]; } xsum;
                xsum.xm = _mm_set1_ps(0.F);
                __m128 xsum2 = _mm_set1_ps(0.F);
                __m128 xA,xA1,xA2,xB1,xB2,x1,x2;
                do {
                    Maybe2<!unit2,unit1>::sse_load(
                        xA,A.get(),A1.get(),A2.get(),A3.get());
                    A+=4; A1+=4; A2+=4; A3+=4;
                    Maybe<unit2>::sse_load(xB1,B.get(),B1.get());
                    B+=4; B1+=4;
                    Maybe<unit2>::sse_load(xB2,B2.get(),B3.get());
                    B2+=4; B3+=4;
                    xA1 = _mm_shuffle_ps(xA,xA,_MM_SHUFFLE(1,1,0,0));
                    xA2 = _mm_shuffle_ps(xA,xA,_MM_SHUFFLE(3,3,2,2));
                    x1 = _mm_mul_ps(xA1,xB1);
                    x2 = _mm_mul_ps(xA2,xB2);
                    xsum.xm = _mm_add_ps(xsum.xm,x1);
                    xsum2 = _mm_add_ps(xsum2,x2);
                } while (--n_4);
                xsum.xm = _mm_add_ps(xsum.xm,xsum2);
                sum1 = 
                    std::complex<float>(xsum.xf[0],xsum.xf[1]) + 
                    std::complex<float>(xsum.xf[2],xsum.xf[3]);
            } else { sum1 = 0.F; }
            if (nb) do {
                sum0 += ZProd<false,false>::prod(*A++ , *B++);
            } while (--nb);
            return sum0 + sum1;
        }
    };

    // algo 23: single precision SSE: all complex
    template <int s, class V1, class V2>
    struct MultVV_Helper<23,s,V1,V2>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        static std::complex<float> call(const V1& v1, const V2& v2)
        {
            const int n = s == UNKNOWN ? int(v1.size()) : s;
            return call2(n,v1.begin().nonConj(),v2.begin().nonConj());
        }
        static std::complex<float> call2(int n, IT1 A, IT2 B)
        {
            TMVStaticAssert(!V2::_conj);
            const bool unit1 = V1::_step == 1;
            const bool unit2 = V2::_step == 1;
            const bool c1 = V1::_conj;

            std::complex<float> sum0(0), sum1;

            if (unit2) {
                while (n && !TMV_Aligned(B.get()) ) {
                    sum0 += ZProd<c1,false>::prod(*A++ , *B++);
                    --n;
                }
            } else if (unit1) {
                while (n && !TMV_Aligned(A.get()) ) {
                    sum0 += ZProd<c1,false>::prod(*A++ , *B++);
                    --n;
                }
            }

            int n_2 = (n>>1);
            int nb = n-(n_2<<1);

            if (n_2) {
                IT1 A1 = A+1;
                IT2 B1 = B+1;

                // A*B = ( Ar Br - Ai Bi , Ar Bi + Ai Br )
                // xsum1 = sum ( Ar Br , Ar Bi )
                // xsum2 = sum ( Ai Br , Ai Bi )
                // We combine these appropriately at the end.
                __m128 xsum1 = _mm_set1_ps(0.F);
                __m128 xsum2 = _mm_set1_ps(0.F);
                __m128 xA,xB,xAr,xAi,x1,x2;
                do {
                    Maybe2<!unit2,unit1>::sse_load(xA,A.get(),A1.get());
                    A+=2; A1+=2;
                    Maybe<unit2>::sse_load(xB,B.get(),B1.get());
                    B+=2; B1+=2;
                    xAr = _mm_shuffle_ps(xA,xA,_MM_SHUFFLE(2,2,0,0));
                    xAi = _mm_shuffle_ps(xA,xA,_MM_SHUFFLE(3,3,1,1));
                    x1 = _mm_mul_ps(xAr,xB);
                    x2 = _mm_mul_ps(xAi,xB);
                    xsum1 = _mm_add_ps(xsum1,x1);
                    xsum2 = _mm_add_ps(xsum2,x2);
                } while (--n_2);
                xsum2 = _mm_shuffle_ps(xsum2,xsum2,_MM_SHUFFLE(2,3,0,1));
                const float mone = Maybe<c1>::select(1.F , -1.F);
                const float one = Maybe<c1>::select(-1.F , 1.F);
                __m128 xmone = _mm_set_ps(one,mone,one,mone);
                xsum2 = _mm_mul_ps(xmone,xsum2);
                union { __m128 xm; float xf[4]; } xsum;
                xsum.xm = _mm_add_ps(xsum1,xsum2);
                sum1 = 
                    std::complex<float>(xsum.xf[0],xsum.xf[1]) + 
                    std::complex<float>(xsum.xf[2],xsum.xf[3]);
            } else { sum1 = 0.F; }
            if (nb) do {
                sum0 += ZProd<c1,false>::prod(*A++ , *B++);
            } while (--nb);
            return sum0 + sum1;
        }
    };
#endif

#ifdef __SSE2__
    // algo 31: double precision SSE2: all real
    template <int s, class V1, class V2>
    struct MultVV_Helper<31,s,V1,V2>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        static double call(const V1& v1, const V2& v2)
        {
            const int n = s == UNKNOWN ? int(v1.size()) : s;
            return call2(n,v1.begin().nonConj(),v2.begin().nonConj());
        }
        static double call2(int n, IT1 A, IT2 B)
        {
            const bool unit1 = V1::_step == 1;
            const bool unit2 = V2::_step == 1;

            double sum0(0), sum1;

            if (unit2) {
                while (n && !TMV_Aligned(B.get()) ) {
                    sum0 += *A++ * *B++;
                    --n;
                }
            } else if (unit1) {
                while (n && !TMV_Aligned(A.get()) ) {
                    sum0 += *A++ * *B++;
                    --n;
                }
            }

            int n_2 = (n>>1);
            int nb = n-(n_2<<1);

            if (n_2) {
                IT1 A1 = A+1;
                IT2 B1 = B+1;

                union { __m128d xm; double xd[2]; } xsum;
                xsum.xm = _mm_set1_pd(0.);
                __m128d xA,xB,x0;
                do {
                    Maybe2<!unit2,unit1>::sse_load(xA,A.get(),A1.get());
                    A+=2; A1+=2;
                    Maybe<unit2>::sse_load(xB,B.get(),B1.get());
                    B+=2; B1+=2;
                    x0 = _mm_mul_pd(xA,xB);
                    xsum.xm = _mm_add_pd(xsum.xm,x0);
                } while (--n_2);
                sum1 = xsum.xd[0] + xsum.xd[1];
            } else { sum1 = 0.; }

            if (nb) sum0 += *A * *B;
            return sum0 + sum1;
        }
    };

    // algo 32: double precision SSE2: v1 real v2 complex
    template <int s, class V1, class V2>
    struct MultVV_Helper<32,s,V1,V2>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        static std::complex<double> call(const V1& v1, const V2& v2)
        {
            const int n = s == UNKNOWN ? int(v1.size()) : s;
            return call2(n,v1.begin().nonConj(),v2.begin().nonConj());
        }
        static std::complex<double> call2(int n, IT1 A, IT2 B)
        {
            TMVStaticAssert(!V2::_conj);
            const bool unit1 = V1::_step == 1;
            const bool unit2 = V2::_step == 1;

            std::complex<double> sum0(0), sum1;

            if (unit2) {
                while (n && !TMV_Aligned(B.get()) ) {
                    sum0 += ZProd<false,false>::prod(*A++ , *B++);
                    --n;
                }
            } else if (unit1) {
                while (n && !TMV_Aligned(A.get()) ) {
                    sum0 += ZProd<false,false>::prod(*A++ , *B++);
                    --n;
                }
            }

            int n_2 = (n>>1);
            int nb = n-(n_2<<1);

            if (n_2) {
                IT1 A1 = A+1;
                IT2 B1 = B+1;

                union { __m128d xm; double xd[2]; } xsum;
                xsum.xm = _mm_set1_pd(0.);
                __m128d xsum2 = _mm_set1_pd(0.);
                __m128d xA,xA1,xA2,xB1,xB2,x1,x2;
                do {
                    Maybe2<!unit2,unit1>::sse_load(xA,A.get(),A1.get());
                    A+=2; A1+=2;
                    Maybe<unit2>::sse_load(xB1,B.get());
                    B+=2; 
                    Maybe<unit2>::sse_load(xB2,B1.get());
                    B1+=2;
                    xA1 = _mm_shuffle_pd(xA,xA,_MM_SHUFFLE2(0,0));
                    xA2 = _mm_shuffle_pd(xA,xA,_MM_SHUFFLE2(1,1));
                    x1 = _mm_mul_pd(xA1,xB1);
                    x2 = _mm_mul_pd(xA2,xB2);
                    xsum.xm = _mm_add_pd(xsum.xm,x1);
                    xsum2 = _mm_add_pd(xsum2,x2);
                } while (--n_2);
                xsum.xm = _mm_add_pd(xsum.xm,xsum2);
                sum1 = std::complex<double>(xsum.xd[0],xsum.xd[1]);
            } else { sum1 = 0.; }
            if (nb) do {
                sum0 += ZProd<false,false>::prod(*A++ , *B++);
            } while (--nb);
            return sum0 + sum1;
        }
    };

    // algo 33: double precision SSE2: all complex
    template <int s, class V1, class V2>
    struct MultVV_Helper<33,s,V1,V2>
    {
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        static std::complex<double> call(const V1& v1, const V2& v2)
        {
            const int n = s == UNKNOWN ? int(v1.size()) : s;
            return call2(n,v1.begin().nonConj(),v2.begin().nonConj());
        }
        static std::complex<double> call2(int n, IT1 A, IT2 B)
        {
            TMVStaticAssert(!V2::_conj);
            const bool c1 = V1::_conj;

            std::complex<double> sum(0);

            if (n) {

                // A*B = ( Ar Br - Ai Bi , Ar Bi + Ai Br )
                // xsum1 = sum ( Ar Br , Ar Bi )
                // xsum2 = sum ( Ai Br , Ai Bi )
                // We combine these appropriately at the end.
                __m128d xsum1 = _mm_set1_pd(0.);
                __m128d xsum2 = _mm_set1_pd(0.);
                __m128d xA,xB,xAr,xAi,x1,x2;
                if ( TMV_Aligned(A.get()) && TMV_Aligned(B.get()) ) {
                    do {
                        Maybe<true>::sse_load(xA,A.get()); ++A;
                        Maybe<true>::sse_load(xB,B.get()); ++B;
                        xAr = _mm_shuffle_pd(xA,xA,_MM_SHUFFLE2(0,0));
                        xAi = _mm_shuffle_pd(xA,xA,_MM_SHUFFLE2(1,1));
                        x1 = _mm_mul_pd(xAr,xB);
                        x2 = _mm_mul_pd(xAi,xB);
                        xsum1 = _mm_add_pd(xsum1,x1);
                        xsum2 = _mm_add_pd(xsum2,x2);
                    } while (--n);
                } else {
                    do {
                        Maybe<true>::sse_loadu(xA,A.get()); ++A;
                        Maybe<true>::sse_loadu(xB,B.get()); ++B;
                        xAr = _mm_shuffle_pd(xA,xA,_MM_SHUFFLE2(0,0));
                        xAi = _mm_shuffle_pd(xA,xA,_MM_SHUFFLE2(1,1));
                        x1 = _mm_mul_pd(xAr,xB);
                        x2 = _mm_mul_pd(xAi,xB);
                        xsum1 = _mm_add_pd(xsum1,x1);
                        xsum2 = _mm_add_pd(xsum2,x2);
                    } while (--n);
                }
                xsum2 = _mm_shuffle_pd(xsum2,xsum2,_MM_SHUFFLE2(0,1));
                const double mone = Maybe<c1>::select(1. , -1.);
                const double one = Maybe<c1>::select(-1. , 1.);
                __m128d xmone = _mm_set_pd(one,mone);
                xsum2 = _mm_mul_pd(xmone,xsum2);
                union { __m128d xm; double xd[2]; } xsum;
                xsum.xm = _mm_add_pd(xsum1,xsum2);
                sum += std::complex<double>(xsum.xd[0],xsum.xd[1]);
            } 
            return sum;
        }
    };
#endif

    // algo 41: recurse very large vector product
    // This isn't for speed reasons - it's for increased accuracy.
    // For large vectors, the incremental additions can be much smaller
    // than the running sum, so the relative errors can be huge.
    // With the recursive algorithm, the relative error is generally
    // closer to the expected few * epsilon.
    template <int s, class V1, class V2>
    struct MultVV_Helper<41,s,V1,V2>
    {
        typedef typename ProdType<V1,V2>::type PT;
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename Traits<PT>::real_type RT;
        static PT call(const V1& v1, const V2& v2)
        {
            const int n = s == UNKNOWN ? int(v1.size()) : s;
            return call2(n,v1.begin().nonConj(),v2.begin().nonConj());
        }
        static PT call2(const int n, const IT1& it1, const IT2& it2)
        {
            const int s1 = IntTraits<s>::half_roundup;
            const int s2 = IntTraits2<s,s1>::diff;

            if (n == 0) return PT(0);
            else if (n > TMV_MultVV_RecurseSize) {
                const int no2 = n/2;
                return (
                    MultVV_Helper<41,s1,V1,V2>::call2(no2,it1,it2) + 
                    MultVV_Helper<41,s2,V1,V2>::call2(
                        n-no2,it1+no2,it2+no2));
            } else {
                return MultVV_Helper<-4,s,V1,V2>::call2(n,it1,it2);
            }
        }
    };

    // algo 90: Call inst
    template <int s, class V1, class V2>
    struct MultVV_Helper<90,s,V1,V2>
    {
        typedef typename ProdType<V1,V2>::type PT;
        static PT call(const V1& v1, const V2& v2)
        { return InstMultVV(v1.xView(),v2.xView()); }
    };

    // algo 96: Swap v1,v2
    template <int s, class V1, class V2>
    struct MultVV_Helper<96,s,V1,V2>
    {
        typedef typename ProdType<V1,V2>::type PT;
        static PT call(const V1& v1, const V2& v2)
        { return MultVV_Helper<-2,s,V2,V1>::call(v2,v1); }
    };


    // algo 97: Conjugate
    template <int s, class V1, class V2>
    struct MultVV_Helper<97,s,V1,V2>
    {
        typedef typename ProdType<V1,V2>::type PT;
        static PT call(const V1& v1, const V2& v2)
        {
            typedef typename V1::const_conjugate_type V1c;
            typedef typename V2::const_conjugate_type V2c;
            V1c v1c = v1.conjugate();
            V2c v2c = v2.conjugate();
            return TMV_CONJ(MultVV_Helper<-2,s,V1c,V2c>::call(v1c,v2c));
        }
    };

    // algo -4: No branches or copies
    template <int s, class V1, class V2>
    struct MultVV_Helper<-4,s,V1,V2>
    {
        typedef typename ProdType<V1,V2>::type PT;
        typedef typename V1::const_nonconj_type::const_iterator IT1;
        typedef typename V2::const_nonconj_type::const_iterator IT2;
        typedef typename Traits<V1>::real_type RT1;
        typedef typename Traits<V2>::real_type RT2;
        enum { unit = V1::_step == 1 || V2::_step == 1 };
        enum { allunit = V1::_step == 1 && V2::_step == 1 };
        enum { allfloat =
            Traits2<RT1,float>::sametype && Traits2<RT2,float>::sametype };
        enum { alldouble =
            Traits2<RT1,double>::sametype && Traits2<RT2,double>::sametype };
        enum { v1real = Traits<V1>::isreal };
        enum { v2real = Traits<V2>::isreal };
        enum { v1complex = Traits<V1>::iscomplex };
        enum { v2complex = Traits<V2>::iscomplex };
        enum { algo = (
                s == 0 ? 0 :
                v1complex && v2real ? 1 :
                V2::_conj ? 2 :
                TMV_OPT == 0 ? 11 :
                s != UNKNOWN ? 11 :
                ( s != UNKNOWN && s <= int(128/sizeof(PT)) ) ? 15 :
#ifdef __SSE__
                ( allfloat && v1real && v2real && allunit ) ? 21 :
                ( allfloat && v1real && v2complex ) ? 22 :
                ( allfloat && v1complex && v2complex && unit ) ? 23 :
#endif
#ifdef __SSE2__
                ( alldouble && v1real && v2real && unit ) ? 31 :
                ( alldouble && v1real && v2complex ) ? 32 :
                ( alldouble && v1complex && v2complex ) ? 33 :
#endif
                ( sizeof(RT2) == 4 && allunit && v1real && v2real ) ? 13 :
                ( sizeof(RT2) == 4 && allunit && v1real && v2complex ) ? 12 :
                ( sizeof(RT2) == 8 && allunit && v1real && v2real ) ? 12 :
                11 ) };
        static PT call(const V1& v1, const V2& v2)
        { return MultVV_Helper<algo,s,V1,V2>::call(v1,v2); }
        static PT call2(const int n, const IT1& it1, const IT2& it2)
        { return MultVV_Helper<algo,s,V1,V2>::call2(n,it1,it2); }
    };

    // algo -3: Determine which algorithm to use.
    template <int s, class V1, class V2>
    struct MultVV_Helper<-3,s,V1,V2>
    {
        typedef typename ProdType<V1,V2>::type PT;
        static PT call(const V1& v1, const V2& v2)
        {
            const int algo = 
#ifdef TMV_VV_RECURSE
                (s == UNKNOWN || s > TMV_MultVV_RecurseSize) ? 41 : 
#endif
                -4;
            return MultVV_Helper<algo,s,V1,V2>::call(v1,v2); 
        }
    };

    // algo -2: Check for inst
    template <int s, class V1, class V2>
    struct MultVV_Helper<-2,s,V1,V2>
    {
        typedef typename ProdType<V1,V2>::type PT;
        static PT call(const V1& v1, const V2& v2)
        {
            typedef typename V1::value_type T1;
            typedef typename V2::value_type T2;
            const bool conj = V2::_conj;
            const bool swap = V1::iscomplex && V2::isreal;
            const bool inst = 
                (s == UNKNOWN || s > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T2>::isinst;
            const int algo = 
                swap ? 96 :
                conj ? 97 :
                inst ? 90 :
                -3;
            return MultVV_Helper<algo,s,V1,V2>::call(v1,v2); 
        }
    };

    template <int s, class V1, class V2>
    struct MultVV_Helper<-1,s,V1,V2>
    {
        typedef typename ProdType<V1,V2>::type PT;
        static PT call(const V1& v1, const V2& v2)
        { return MultVV_Helper<-2,s,V1,V2>::call(v1,v2); }
    };

    template <class V1, class V2>
    static inline typename ProdType<V1,V2>::type MultVV(
        const BaseVector_Calc<V1>& v1, const BaseVector_Calc<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same));
        TMVAssert(v1.size() == v2.size());
        const int size = Sizes<V1::_size,V2::_size>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        return MultVV_Helper<-2,size,V1v,V2v>::call(v1v,v2v);
    }

    template <class V1, class V2>
    static inline typename ProdType<V1,V2>::type InlineMultVV(
        const BaseVector_Calc<V1>& v1, const BaseVector_Calc<V2>& v2)
    {
        TMVStaticAssert((Sizes<V1::_size,V2::_size>::same));
        TMVAssert(v1.size() == v2.size());
        const int size = Sizes<V1::_size,V2::_size>::size;
        typedef typename V1::const_cview_type V1v;
        typedef typename V2::const_cview_type V2v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_CREF(V2,V2v) v2v = v2.cView();
        return MultVV_Helper<-3,size,V1v,V2v>::call(v1v,v2v);
    }

} // namespace tmv

#endif 
