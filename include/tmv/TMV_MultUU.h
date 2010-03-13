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


#ifndef TMV_MultUU_H
#define TMV_MultUU_H

#include "TMV_MultMM.h"
#include "TMV_MultUM.h"
#include "TMV_MultUV.h"

#ifdef PRINTALGO_UU
#include <iostream>
#endif

// Use the specialized 1,2,3,4 sized algorithms for the end of the 
// recursive algorithm.
#if TMV_OPT >= 2
#define TMV_OPT_CLEANUP
#endif

// Q1 is the maximum nops to unroll.
#if TMV_OPT >= 3
#define TMV_Q1 200 
#elif TMV_OPT >= 2
#define TMV_Q1 25
#elif TMV_OPT >= 1
#define TMV_Q1 9
#else
#define TMV_Q1 0
#endif

// Q2 is the minimum size to keep recursing.
#define TMV_Q2 8

namespace tmv {

    // Defined below:
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Tri_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void NoAliasMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Tri_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Tri_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void AliasMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Tri_Mutable<M3>& m3);
    template <class M1, int ix, class T, class M2>
    inline void MultEqMM(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2);
    template <class M1, int ix, class T, class M2>
    inline void NoAliasMultEqMM(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2);
    template <class M1, int ix, class T, class M2>
    inline void AliasMultEqMM(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2);

    // Defined in TMV_MultUU.cpp
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstMultMM(
        const T3 x,
        const ConstUpperTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1, 
        const ConstUpperTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2, 
        UpperTriMatrixView<T3,UnknownDiag> m3);
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstAddMultMM(
        const T3 x,
        const ConstUpperTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1, 
        const ConstUpperTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2,
        UpperTriMatrixView<T3,NonUnitDiag> m3);

    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstMultMM(
        const T3 x,
        const ConstLowerTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstLowerTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2, 
        LowerTriMatrixView<T3,UnknownDiag> m3);
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstAddMultMM(
        const T3 x,
        const ConstLowerTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m2,
        const ConstLowerTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m1, 
        LowerTriMatrixView<T3,NonUnitDiag> m3);

    template <int algo, int s, bool add, 
              int ix, class T, class M1, class M2, class M3> 
    struct MultUU_Helper;

    // algo 0: Trivial, nothing to do.
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<0,s,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {}
    };

    // algo 1: s == 1, so reduces to scalar product
    template <bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<1,1,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_UU
            std::cout<<"UU algo 1: N,s,x = "<<1<<','<<1<<','<<T(x)<<std::endl;
#endif
            TMVStaticAssert(!M3::munit);
            Maybe<add>::add(m3.ref(0,0), x*m1.cref(0,0)*m2.cref(0,0));
        }
    };

    // algo 11: UpperTri loop over n
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<11,s,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = s==UNKNOWN ? int(m3.size()) : s;
#ifdef PRINTALGO_UU
            std::cout<<"UU algo 11: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            const bool u1 = M1::munit;
            const bool u2 = M2::munit;
            typedef typename M2::value_type T2;
            typedef typename Maybe<!u2>::template ProdType<T2,T>::type PT2;
            typedef typename M1::const_subtrimatrix_type M1s;
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::const_col_sub_type M2c;
            typedef typename M3::col_sub_type M3c;
            const int ix2 = u2 ? ix : 0;
            const int xx = UNKNOWN;

            for(int j=N-1;j>=0;--j) {
                // m3.col(j,0,j+1) = m1.subTriMatrix(0,j+1) * m2.col(j,0,j+1)
                // ==>
                // m3.col(j,0,j) = m1.col(j,0,j) * m2(j,j) +
                //                 m1.subTriMatrix(0,j) * m2.col(j,0,j)
                // m3(j,j) = m1(j,j) * m2(j,j)
                const Scaling<ix2,PT2> xd(Maybe<!u2>::prod(m2.cref(j,j),x));
                M1s m1s = m1.cSubTriMatrix(0,j);
                M1c m1c = m1.get_col(j,0,j);
                M2c m2c = m2.get_col(j,0,j);
                M3c m3c = m3.get_col(j,0,j);
                MultXV_Helper<-4,xx,add,ix2,PT2,M1c,M3c>::call(xd,m1c,m3c);
                MultUV_Helper<-4,xx,true,ix,T,M1s,M2c,M3c>::call(
                    x,m1s,m2c,m3c);
                Maybe2<!M3::munit,add>::add(
                    m3.ref(j,j), Maybe<!u1>::prod(m1.cref(j,j),xd));
            }
        }
    };

    // algo 12: UpperTri loop over m
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<12,s,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = s==UNKNOWN ? int(m3.size()) : s;
#ifdef PRINTALGO_UU
            std::cout<<"UU algo 12: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            const bool u1 = M1::munit;
            const bool u2 = M2::munit;
            typedef typename M1::value_type T1;
            typedef typename Maybe<!u1>::template ProdType<T1,T>::type PT1;
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M2::const_row_sub_type M2r;
            typedef typename M2::const_subtrimatrix_type M2s;
            typedef typename M2s::const_transpose_type M2t;
            typedef typename M3::row_sub_type M3r;
            const int ix1 = u1 ? ix : 0;
            const int xx = UNKNOWN;

            for(int i=0;i<N;++i) {
                // m3.row(i,i,N) = m1.row(i,i,N) * m2.subTriMatrix(i,N)
                // ==>
                // m3.row(i,i+1,N) = m1.row(i,i+1,N) * m2.subTriMatrix(i+1,N) +
                //                   m1(i,i) * m2.row(i,i+1,N)
                // m3(i,i) = m1(i,i) * m2(i,i)
                const Scaling<ix1,PT1> xd(Maybe<!u1>::prod(m1.cref(i,i),x));
                M1r m1r = m1.get_row(i,i+1,N);
                M2t m2t = m2.cSubTriMatrix(i+1,N).transpose();
                M2r m2r = m2.get_row(i,i+1,N);
                M3r m3r = m3.get_row(i,i+1,N);
                MultUV_Helper<-4,xx,add,ix,T,M2t,M1r,M3r>::call(
                    x,m2t,m1r,m3r);
                MultXV_Helper<-4,xx,true,ix1,PT1,M2r,M3r>::call(xd,m2r,m3r);
                Maybe2<!M3::munit,add>::add(
                    m3.ref(i,i), Maybe<!u2>::prod(m2.cref(i,i),xd));
            }
        }
    };

    // algo 13: UpperTri loop over k
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<13,s,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = s==UNKNOWN ? int(m3.size()) : s;
#ifdef PRINTALGO_UU
            std::cout<<"UU algo 13: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            const bool u1 = M1::munit;
            const bool u2 = M2::munit;
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename Maybe<!u1>::template ProdType<T1,T>::type PT1;
            typedef typename Maybe<!u2>::template ProdType<T2,T>::type PT2;
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::const_row_sub_type M2r;
            typedef typename M3::row_sub_type M3r;
            typedef typename M3::col_sub_type M3c;
            typedef typename M3::submatrix_type M3s;
            const int ix1 = u1 ? ix : 0;
            const int ix2 = u2 ? ix : 0;
            const int xx = UNKNOWN;

            for(int k=N-1;k>=0;--k) {
                // m3.subMatrix(0,k+1,k,N) = m1.col(k,0,k+1) ^ m2.row(k,k,N)
                // ==>
                // m3.subMatrix(0,k,k+1,N) = m1.col(k,0,k) ^ m2.row(k,k+1,N)
                // m3.col(k,0,k) = m1.col(k,0,k) * m2(k,k)
                // m3.row(k,k+1,N) = m1(k,k) * m2.row(k,k+1,N)
                // m3(k,k) = m1(k,k) * m2(k,k)
                const Scaling<ix1,PT1> xd1(Maybe<!u1>::prod(m1.cref(k,k),x));
                const Scaling<ix2,PT2> xd2(Maybe<!u2>::prod(m2.cref(k,k),x));
                M1c m1c = m1.get_col(k,0,k);
                M2r m2r = m2.get_row(k,k+1,N);
                M3c m3c = m3.get_col(k,0,k);
                M3r m3r = m3.get_row(k,k+1,N);
                M3s m3s = m3.cSubMatrix(0,k,k+1,N);
                Rank1VVM_Helper<-4,xx,xx,true,ix,T,M1c,M2r,M3s>::call(
                    x,m1c,m2r,m3s);
                MultXV_Helper<-4,xx,add,ix2,PT2,M1c,M3c>::call(xd2,m1c,m3c);
                MultXV_Helper<-4,xx,true,ix1,PT1,M2r,M3r>::call(xd1,m2r,m3r);
                Maybe2<!M3::munit,add>::add(
                    m3.ref(k,k), Maybe<!u1>::prod(m1.cref(k,k),xd2));
            }
        }
    };

    // algo 16: UpperTri: Unroll small case
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<16,s,add,ix,T,M1,M2,M3> 
    {
        template <int I, int N>
        struct Unroller
        {
            static void unroll(
                const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
            {
                // [ C00 C01 ] = [ A00 A01 ] [ B00 B01 ]
                // [  0  C11 ]   [  0  A11 ] [  0  B11 ]

                const int Nx = N/2;
                const int Ny = N-Nx;
                const int I1 = I+Nx;
                const int I2 = I+N;
                typedef typename M1::const_subtrimatrix_type M1st;
                typedef typename M1::const_submatrix_type M1sm;
                typedef typename M1sm::const_transpose_type M1smt;
                typedef typename M2::const_submatrix_type M2sm;
                typedef typename M2::const_subtrimatrix_type M2st;
                typedef typename M2st::const_transpose_type M2stt;
                typedef typename M3::submatrix_type M3sm;
                typedef typename M3sm::transpose_type M3smt;

                M1st A00 = m1.cSubTriMatrix(I,I1);
                M1smt A01t = m1.cSubMatrix(I,I1,I1,I2).transpose();
                M2sm B01 = m2.cSubMatrix(I,I1,I1,I2);
                M2stt B11t = m2.cSubTriMatrix(I1,I2).transpose();
                M3sm C01 = m3.cSubMatrix(I,I1,I1,I2);
                M3smt C01t = C01.transpose();

                // C01 (+)= A01 B11
                MultUM_Helper<-4,Ny,Nx,add,ix,T,M2stt,M1smt,M3smt>::call(
                    x,B11t,A01t,C01t);

                // C01 += A00 B01
                MultUM_Helper<-4,Nx,Ny,true,ix,T,M1st,M2sm,M3sm>::call(
                    x,A00,B01,C01);

                // C00 (+)= A00 B00
                Unroller<I,Nx>::unroll(x,m1,m2,m3);

                // C11 (+)= A11 B11
                Unroller<I1,Ny>::unroll(x,m1,m2,m3);
            }
        };
        template <int I>
        struct Unroller<I,1>
        {
            static inline void unroll(
                const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
            {
                const bool u1 = M1::munit;
                const bool u2 = M2::munit;
                Maybe2<!M3::munit,add>::add(
                    m3.ref(I,I), Maybe<!u1>::prod(
                            m1.cref(I,I),Maybe<!u2>::prod(m2.cref(I,I),x)));
            }
        };
        template <int I>
        struct Unroller<I,0>
        {
            static inline void unroll(
                const Scaling<ix,T>& , const M1& , const M2& , M3& ) {}
        };
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_UU
            const int N = s==UNKNOWN ? int(m3.size()) : s;
            std::cout<<"UU algo 16: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            Unroller<0,s>::unroll(x,m1,m2,m3); 
        }
    };
    template <bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<16,UNKNOWN,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = m3.size();
#ifdef PRINTALGO_UU
            std::cout<<"UU algo 16: N,s,x = "<<N<<','<<UNKNOWN<<
                ','<<T(x)<<std::endl;
#endif
            const bool rxr = M1::mrowmajor && M3::mrowmajor;
            const bool crx = M1::mcolmajor && M2::mrowmajor;
            const bool xcc = M2::mcolmajor && M3::mcolmajor;
            const int algo2 = rxr ? 12 : crx ? 13 : xcc ? 11 : 13;
            const int algo1 = M3::munit ? 0 : 1;
            switch (N) {
              case 0 :
                   // do nothing
                   break;
              case 1 :
                   MultUU_Helper<algo1,1,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
                   break;
              case 2 :
                   MultUU_Helper<16,2,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
                   break;
              case 3 :
                   MultUU_Helper<16,3,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
                   break;
              default :
                   MultUU_Helper<algo2,UNKNOWN,add,ix,T,M1,M2,M3>::call(
                       x,m1,m2,m3);
            }
        }
    };

    // algo 17: Split the UpperTriMatrix into 3 sections and recurse
    // the calculation on each of them:
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<17,s,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = s==UNKNOWN ? int(m3.size()) : s;
#ifdef PRINTALGO_UU
            std::cout<<"UU algo 17: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif

#if (TMV_Q2 == 1)
            const int algo2 = (s == UNKNOWN || s == 1) ? 1 : 0;
#else
            const int sp1 = IntTraits<s>::Sp1;
            const int sp2 = IntTraits<sp1>::Sp1;
            // nops = 1/6 n(n+1)(n+2)
            const int nops = 
                IntTraits2<IntTraits2<s,sp1>::safeprod,sp2>::safeprod / 6;
            const bool unroll = 
                s > 20 ? false :
                s == UNKNOWN ? false :
                nops > TMV_Q1 ? false :
                s <= 10;
            const bool rxr = M1::mrowmajor && M3::mrowmajor;
            const bool crx = M1::mcolmajor && M2::mrowmajor;
            const bool xcc = M2::mcolmajor && M3::mcolmajor;
            const int algo2 = 
                s == 0 ? 0 :
                s == 1 ? ( M3::munit ? 0 : 1 ) :
                (s != UNKNOWN && s > TMV_Q2) ? 0 :
#ifdef TMV_OPT_CLEANUP
                s == UNKNOWN ? 16 :
#endif
                unroll ? 16 : 
                rxr ? 12 : crx ? 13 : xcc ? 11 : 13;
#endif
            const int algo3 =  // The algorithm for N > Q2
                (s == UNKNOWN || s > TMV_Q2) ? 17 : 0;
            const int algo4 =  // The algorithm for MultUM
                s == UNKNOWN ? -2 : s > 16 ? -3 : s > TMV_Q2 ? -4 : 0;

            if (N > TMV_Q2) {
                // [ C00 C01 ] = [ A00 A01 ] [ B00 B01 ]
                // [  0  C11 ]   [  0  A11 ] [  0  B11 ]

                const int Nx = N > 16 ? ((((N-1)>>5)+1)<<4) : (N>>1);
                // (If N > 16, round N/2 up to a multiple of 16.)
                const int sx = IntTraits<s>::half_roundup;
                const int sy = s == UNKNOWN ? s : s-sx;

                typedef typename M1::const_subtrimatrix_type M1st;
                typedef typename M1::const_submatrix_type M1sm;
                typedef typename M1sm::const_transpose_type M1smt;
                typedef typename M2::const_subtrimatrix_type M2st;
                typedef typename M2::const_submatrix_type M2sm;
                typedef typename M2st::const_transpose_type M2stt;
                typedef typename M3::subtrimatrix_type M3st;
                typedef typename M3::submatrix_type M3sm;
                typedef typename M3sm::transpose_type M3smt;

                M1st A00 = m1.cSubTriMatrix(0,Nx);
                M1smt A01t = m1.cSubMatrix(0,Nx,Nx,N).transpose();
                M1st A11 = m1.cSubTriMatrix(Nx,N);
                M2st B00 = m2.cSubTriMatrix(0,Nx);
                M2sm B01 = m2.cSubMatrix(0,Nx,Nx,N);
                M2st B11 = m2.cSubTriMatrix(Nx,N);
                M2stt B11t = B11.transpose();
                M3st C00 = m3.cSubTriMatrix(0,Nx);
                M3sm C01 = m3.cSubMatrix(0,Nx,Nx,N);
                M3smt C01t = C01.transpose();
                M3st C11 = m3.cSubTriMatrix(Nx,N);

                // C01 (+)= A01 B11
                MultUM_Helper<algo4,sy,sx,add,ix,T,M2stt,M1smt,M3smt>::call(
                    x,B11t,A01t,C01t);

                // C01 += A00 B01
                MultUM_Helper<algo4,sx,sy,true,ix,T,M1st,M2sm,M3sm>::call(
                    x,A00,B01,C01);

                // C00 (+)= x A00 B00
                MultUU_Helper<algo3,sx,add,ix,T,M1st,M2st,M3st>::call(
                    x,A00,B00,C00);

                // C11 (+)= x A11 B11
                MultUU_Helper<algo3,sy,add,ix,T,M1st,M2st,M3st>::call(
                    x,A11,B11,C11);
            }
            else MultUU_Helper<algo2,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo 21: LowerTri loop over n
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<21,s,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = s==UNKNOWN ? int(m3.size()) : s;
#ifdef PRINTALGO_UU
            std::cout<<"UU algo 21: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            const bool u1 = M1::munit;
            const bool u2 = M2::munit;
            typedef typename M2::value_type T2;
            typedef typename Maybe<!u2>::template ProdType<T2,T>::type PT2;
            typedef typename M1::const_subtrimatrix_type M1s;
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::const_col_sub_type M2c;
            typedef typename M3::col_sub_type M3c;
            const int ix2 = u2 ? ix : 0;
            const int xx = UNKNOWN;

            for(int j=0;j<N;++j) {
                // m3.col(j,j,N) = m1.subTriMatrix(j,N) * m2.col(j,j,N)
                // ==>
                // m3.col(j,j+1,N) = m1.col(j,j+1,N) * m2(j,j) +
                //                   m1.subTriMatrix(j+1,N) * m2.col(j,j+1,N)
                // m3(j,j) = m1(j,j) * m2(j,j)
                const Scaling<ix2,PT2> xd(Maybe<!u2>::prod(m2.cref(j,j),x));
                M1s m1s = m1.cSubTriMatrix(j+1,N);
                M1c m1c = m1.get_col(j,j+1,N);
                M2c m2c = m2.get_col(j,j+1,N);
                M3c m3c = m3.get_col(j,j+1,N);
                MultXV_Helper<-4,xx,add,ix2,PT2,M1c,M3c>::call(xd,m1c,m3c);
                MultUV_Helper<-4,xx,true,ix,T,M1s,M2c,M3c>::call(
                    x,m1s,m2c,m3c);
                Maybe2<!M3::munit,add>::add(
                    m3.ref(j,j), Maybe<!u1>::prod(m1.cref(j,j),xd));
            }
        }
    };

    // algo 22: LowerTri loop over m
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<22,s,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = s==UNKNOWN ? int(m3.size()) : s;
#ifdef PRINTALGO_UU
            std::cout<<"UU algo 22: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            const bool u1 = M1::munit;
            const bool u2 = M2::munit;
            typedef typename M1::value_type T1;
            typedef typename Maybe<!u1>::template ProdType<T1,T>::type PT1;
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M2::const_row_sub_type M2r;
            typedef typename M2::const_subtrimatrix_type M2s;
            typedef typename M2s::const_transpose_type M2t;
            typedef typename M3::row_sub_type M3r;
            const int ix1 = u1 ? ix : 0;
            const int xx = UNKNOWN;

            for(int i=0;i<N;++i) {
                // m3.row(i,0,i+1) = m1.row(i,0,i+1) * m2.subTriMatrix(0,i+1)
                // ==>
                // m3.row(i,0,i) = m1.row(i,0,i) * m2.subTriMatrix(0,i) +
                //                 m1(i,i) * m2.row(i,0,i)
                // m3(i,i) = m1(i,i) * m2(i,i)
                const Scaling<ix1,PT1> xd(Maybe<!u1>::prod(m1.cref(i,i),x));
                M1r m1r = m1.get_row(i,0,i);
                M2t m2t = m2.cSubTriMatrix(0,i).transpose();
                M2r m2r = m2.get_row(i,0,i);
                M3r m3r = m3.get_row(i,0,i);
                MultUV_Helper<-4,xx,add,ix,T,M2t,M1r,M3r>::call(
                    x,m2t,m1r,m3r);
                MultXV_Helper<-4,xx,true,ix1,PT1,M2r,M3r>::call(xd,m2r,m3r);
                Maybe2<!M3::munit,add>::add(
                    m3.ref(i,i), Maybe<!u2>::prod(m2.cref(i,i),xd));
            }
        }
    };

    // algo 23: LowerTri loop over k
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<23,s,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = s==UNKNOWN ? int(m3.size()) : s;
#ifdef PRINTALGO_UU
            std::cout<<"UU algo 23: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            const bool u1 = M1::munit;
            const bool u2 = M2::munit;
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename Maybe<!u1>::template ProdType<T1,T>::type PT1;
            typedef typename Maybe<!u2>::template ProdType<T2,T>::type PT2;
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::const_row_sub_type M2r;
            typedef typename M3::row_sub_type M3r;
            typedef typename M3::col_sub_type M3c;
            typedef typename M3::submatrix_type M3s;
            const int ix1 = u1 ? ix : 0;
            const int ix2 = u2 ? ix : 0;
            const int xx = UNKNOWN;

            for(int k=0;k<N;++k) {
                // m3.subMatrix(k,N,0,k+1) = m1.col(k,k,N) ^ m2.row(k,0,k+1)
                // ==>
                // m3.subMatrix(k+1,N,0,k) = m1.col(k,k+1,N) ^ m2.row(k,0,k)
                // m3.col(k,k+1,N) = m1.col(k,k+1,N) * m2(k,k)
                // m3.row(k,0,k) = m1(k,k) * m2.row(k,0,k)
                // m3(k,k) = m1(k,k) * m2(k,k)
                const Scaling<ix1,PT1> xd1(Maybe<!u1>::prod(m1.cref(k,k),x));
                const Scaling<ix2,PT2> xd2(Maybe<!u2>::prod(m2.cref(k,k),x));
                M1c m1c = m1.get_col(k,k+1,N);
                M2r m2r = m2.get_row(k,0,k);
                M3c m3c = m3.get_col(k,k+1,N);
                M3r m3r = m3.get_row(k,0,k);
                M3s m3s = m3.cSubMatrix(k+1,N,0,k);
                Rank1VVM_Helper<-4,xx,xx,true,ix,T,M1c,M2r,M3s>::call(
                    x,m1c,m2r,m3s);
                MultXV_Helper<-4,xx,add,ix2,PT2,M1c,M3c>::call(xd2,m1c,m3c);
                MultXV_Helper<-4,xx,true,ix1,PT1,M2r,M3r>::call(xd1,m2r,m3r);
                Maybe2<!M3::munit,add>::add(
                    m3.ref(k,k), Maybe<!u1>::prod(m1.cref(k,k),xd2));
            }
        }
    };

    // algo 26: LowerTri: Unroll small case
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<26,s,add,ix,T,M1,M2,M3> 
    {
        template <int I, int N>
        struct Unroller
        {
            static void unroll(
                const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
            {
                // [ C00  0  ] = [ A00  0  ] [ B00  0  ]
                // [ C10 C11 ]   [ A10 A11 ] [ B10 B11 ]

                const int Nx = N/2;
                const int Ny = N-Nx;
                const int I1 = I+Nx;
                const int I2 = I+N;
                typedef typename M1::const_subtrimatrix_type M1st;
                typedef typename M1::const_submatrix_type M1sm;
                typedef typename M1sm::const_transpose_type M1smt;
                typedef typename M2::const_submatrix_type M2sm;
                typedef typename M2::const_subtrimatrix_type M2st;
                typedef typename M2st::const_transpose_type M2stt;
                typedef typename M3::submatrix_type M3sm;
                typedef typename M3sm::transpose_type M3smt;

                M1smt A10t = m1.cSubMatrix(I1,I2,I,I1).transpose();
                M1st A11 = m1.cSubTriMatrix(I1,I2);
                M2sm B10 = m2.cSubMatrix(I1,I2,I,I1);
                M2stt B00t = m2.cSubTriMatrix(I,I1).transpose();
                M3sm C10 = m3.cSubMatrix(I1,I2,I,I1);
                M3smt C10t = C10.transpose();

                // C10 (+)= A10 B00
                MultUM_Helper<-4,Nx,Ny,add,ix,T,M2stt,M1smt,M3smt>::call(
                    x,B00t,A10t,C10t);

                // C10 += A11 B10
                MultUM_Helper<-4,Ny,Nx,true,ix,T,M1st,M2sm,M3sm>::call(
                    x,A11,B10,C10);

                // C00 (+)= A00 B00
                Unroller<I,Nx>::unroll(x,m1,m2,m3);

                // C11 (+)= A11 B11
                Unroller<I1,Ny>::unroll(x,m1,m2,m3);
            }
        };
        template <int I>
        struct Unroller<I,1>
        {
            static inline void unroll(
                const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
            {
                const bool u1 = M1::munit;
                const bool u2 = M2::munit;
                Maybe2<!M3::munit,add>::add(
                    m3.ref(I,I), Maybe<!u1>::prod(
                            m1.cref(I,I),Maybe<!u2>::prod(m2.cref(I,I),x)));
            }
        };
        template <int I>
        struct Unroller<I,0>
        {
            static inline void unroll(
                const Scaling<ix,T>& , const M1& , const M2& , M3& ) {}
        };
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_UU
            std::cout<<"UU algo 26: N,s,x = "<<2<<','<<2<<','<<T(x)<<std::endl;
#endif
            Unroller<0,s>::unroll(x,m1,m2,m3); 
        }
    };
    template <bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<26,UNKNOWN,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = m3.size();
#ifdef PRINTALGO_UU
            std::cout<<"UU algo 26: N,s,x = "<<N<<','<<UNKNOWN<<
                ','<<T(x)<<std::endl;
#endif
            const bool rxr = M1::mrowmajor && M3::mrowmajor;
            const bool crx = M1::mcolmajor && M2::mrowmajor;
            const bool xcc = M2::mcolmajor && M3::mcolmajor;
            const int algo2 = rxr ? 22 : crx ? 23 : xcc ? 21 : 23;
            const int algo1 = M3::munit ? 0 : 1;
            switch (N) {
              case 0 :
                   // do nothing
                   break;
              case 1 :
                   MultUU_Helper<algo1,1,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
                   break;
              case 2 :
                   MultUU_Helper<26,2,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
                   break;
              case 3 :
                   MultUU_Helper<26,3,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
                   break;
              default :
                   MultUU_Helper<algo2,UNKNOWN,add,ix,T,M1,M2,M3>::call(
                       x,m1,m2,m3);
            }
        }
    };

    // algo 27: Split the LowerTriMatrix into 3 sections and recurse
    // the calculation on each of them:
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<27,s,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = s==UNKNOWN ? int(m3.size()) : s;
#ifdef PRINTALGO_UU
            std::cout<<"UU algo 27: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif

#if (TMV_Q2 == 1)
            const int algo2 = (s == UNKNOWN || s == 1) ? 1 : 0;
#else
            const int sp1 = IntTraits<s>::Sp1;
            const int sp2 = IntTraits<sp1>::Sp1;
            // nops = 1/6 n(n+1)(n+2)
            const int nops = 
                IntTraits2<IntTraits2<s,sp1>::safeprod,sp2>::safeprod / 6;
            const bool unroll = 
                s > 20 ? false :
                s == UNKNOWN ? false :
                nops > TMV_Q1 ? false :
                s <= 10;
            const bool rxr = M1::mrowmajor && M3::mrowmajor;
            const bool crx = M1::mcolmajor && M2::mrowmajor;
            const bool xcc = M2::mcolmajor && M3::mcolmajor;
            const int algo2 = 
                s == 0 ? 0 :
                s == 1 ? ( M3::munit ? 0 : 1 ) :
                (s != UNKNOWN && s > TMV_Q2) ? 0 :
#ifdef TMV_OPT_CLEANUP
                s == UNKNOWN ? 26 :
#endif
                unroll ? 26 : 
                rxr ? 22 : crx ? 23 : xcc ? 21 : 23;
#endif
            const int algo3 =  // The algorithm for N > Q2
                (s == UNKNOWN || s > TMV_Q2) ? 27 : 0;
            const int algo4 =  // The algorithm for MultUM
                s == UNKNOWN ? -2 : s > 16 ? -3 : s > TMV_Q2 ? -4 : 0;

            if (N > TMV_Q2) {
                // [ C00  0  ] = [ A00  0  ] [ B00  0  ]
                // [ C10 C11 ]   [ A10 A11 ] [ B10 B11 ]

                const int Nx = N > 16 ? ((((N-1)>>5)+1)<<4) : (N>>1);
                // (If N > 16, round N/2 up to a multiple of 16.)
                const int sx = IntTraits<s>::half_roundup;
                const int sy = s == UNKNOWN ? s : s-sx;

                typedef typename M1::const_subtrimatrix_type M1st;
                typedef typename M1::const_submatrix_type M1sm;
                typedef typename M1sm::const_transpose_type M1smt;
                typedef typename M2::const_subtrimatrix_type M2st;
                typedef typename M2::const_submatrix_type M2sm;
                typedef typename M2st::const_transpose_type M2stt;
                typedef typename M3::subtrimatrix_type M3st;
                typedef typename M3::submatrix_type M3sm;
                typedef typename M3sm::transpose_type M3smt;

                M1st A00 = m1.cSubTriMatrix(0,Nx);
                M1smt A10t = m1.cSubMatrix(Nx,N,0,Nx).transpose();
                M1st A11 = m1.cSubTriMatrix(Nx,N);
                M2st B00 = m2.cSubTriMatrix(0,Nx);
                M2stt B00t = B00.transpose();
                M2sm B10 = m2.cSubMatrix(Nx,N,0,Nx);
                M2st B11 = m2.cSubTriMatrix(Nx,N);
                M3st C00 = m3.cSubTriMatrix(0,Nx);
                M3sm C10 = m3.cSubMatrix(Nx,N,0,Nx);
                M3smt C10t = C10.transpose();
                M3st C11 = m3.cSubTriMatrix(Nx,N);

                // C10 (+)= A10 B00
                MultUM_Helper<algo4,sx,sy,add,ix,T,M2stt,M1smt,M3smt>::call(
                    x,B00t,A10t,C10t);

                // C10 += A11 B10
                MultUM_Helper<algo4,sy,sx,true,ix,T,M1st,M2sm,M3sm>::call(
                    x,A11,B10,C10);

                // C00 (+)= x A00 B00
                MultUU_Helper<algo3,sx,add,ix,T,M1st,M2st,M3st>::call(
                    x,A00,B00,C00);

                // C11 (+)= x A11 B11
                MultUU_Helper<algo3,sy,add,ix,T,M1st,M2st,M3st>::call(
                    x,A11,B11,C11);
            }
            else MultUU_Helper<algo2,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo 87: Use temporary for m1*m2
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<87,s,add,ix,T,M1,M2,M3>
    {
        // This algo is slightly complicated by the fact that we allow
        // UnknownDiag.  If m3 is UnknownDiag, then its mshape dictates
        // NonUnitDiag.  But that means at m3c is not copyable back to m3 
        // directly (can't copy NonUnitDiag -> UnitDiag).
        // So we have another layer of indirection at the end to make sure 
        // that an UnknownDiag m3 is copied correctly.
        template <bool unknowndiag, class M3c> 
        struct copyBack
        { // unknowndiag = false
            static inline void call(
                const Scaling<ix,T>& x, const M3c& m3c, M3& m3)
            { MultXU_Helper<-2,s,add,ix,T,M3c,M3>::call(x,m3c,m3); }
        };
        template <class M3c>
        struct copyBack<true,M3c>
        {
            static inline void call(
                const Scaling<1,T>& , const M3c& m3c, M3& m3)
            {
                TMVStaticAssert(ix == 1);
                TMVStaticAssert(!add);
                if (m3.isunit()) {
                    typedef typename M3c::const_unitdiag_type M3cu;
                    M3cu m3cu = m3c.viewAsUnitDiag();
                    CopyU_Helper<-2,s,M3cu,M3>::call(m3cu,m3);
                } else {
                    CopyU_Helper<-2,s,M3c,M3>::call(m3c,m3);
                }
            }
        };
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = s == UNKNOWN ? int(m3.size()) : s;
#ifdef PRINTALGO_UU
            std::cout<<"UU algo 87: N,s,x = "<<N<< ','<<s<<','<<T(x)<<std::endl;
#endif
            typedef typename Traits<T>::real_type RT;
            const Scaling<1,RT> one;
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename Traits2<T1,T2>::type PT3;
            const bool rm = M1::mrowmajor && M2::mrowmajor;
            const int s3 = M3::mshape;
            typedef typename MCopyHelper<PT3,s3,s,s,rm,false>::type M3c;
            M3c m3c(N);
            typedef typename M3c::view_type M3cv;
            typedef typename M3c::const_view_type M3ccv;
            M3cv m3cv = m3c.view();
            M3ccv m3ccv = m3c.view();
            MultUU_Helper<-2,s,false,1,RT,M1,M2,M3cv>::call(
                one,m1,m2,m3cv);
            // can't be unitdiag unless ix == 1 and !add
            const bool unknowndiag = M3::munknowndiag && ix == 1 && !add;
            copyBack<unknowndiag,M3ccv>::call(x,m3ccv,m3);
        }
    };

    // algo -4: No branches or copies
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<-4,s,add,ix,T,M1,M2,M3> 
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        { MultUU_Helper<-3,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3); }
    };

    // algo -3: Determine which algorithm to use
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<-3,s,add,ix,T,M1,M2,M3> 
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            // Possible algorithms are:
            //
            // Trivial and special for small TriMatrix sizes.
            //  0 = s == 0, so nothing to do
            //  1 = s == 1: reduces to trivial scalar product function
            //
            // UpperTri:
            // 11 = loop over n: MultUV
            // 12 = loop over m: MultUV 
            // 13 = loop over k: Rank1
            // 16 = Unroll small case
            // 17 = split each trimatrix into 3 submatrices and recurse
            // 
            // LowerTri:
            // 21 = loop over n: MultUV
            // 22 = loop over m: MultUV 
            // 23 = loop over k: Rank1
            // 26 = Unroll small case
            // 27 = split each trimatrix into 3 submatrices and recurse

            const bool upper = M1::mupper;
#if TMV_OPT == 0 
            const bool rxr = M1::mrowmajor && M3::mrowmajor;
            const bool crx = M1::mcolmajor && M2::mrowmajor;
            const bool xcc = M2::mcolmajor && M3::mcolmajor;
            const int algo =
                upper ? ( rxr ? 12 : crx ? 13 : xcc ? 11 : 13 ) :
                ( rxr ? 22 : crx ? 23 : xcc ? 21 : 23 );
#else
            const int s2 = s > 20 ? UNKNOWN : s;
            const int s2p1 = IntTraits<s2>::Sp1;
            const int s2p2 = IntTraits<s2p1>::Sp1;
            // nops = 1/6 n(n+1)(n+2)
            const int nops = 
                IntTraits2<IntTraits2<s2,s2p1>::safeprod,s2p2>::safeprod / 6;
            const bool unroll = 
                s == UNKNOWN ? false :
                nops > TMV_Q1 ? false :
                s <= 10;
            const int algo = 
                s == 0 ? 0 :
                s == 1 ? ( M3::munit ? 0 : 1 ) :
                unroll ? ( upper ? 16 : 26 ) :
                upper ? 17 : 27;
#endif
#ifdef PRINTALGO_UU
            const int N = s==UNKNOWN ? int(m3.size()) : s;
            std::cout<<"InlineMultUU: x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"N = "<<N<<std::endl;
            std::cout<<"s = "<<s<<std::endl;
            std::cout<<"add = "<<add<<", algo = "<<algo<<std::endl;
            //std::cout<<"m1 = "<<m1<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
            //std::cout<<"m3 = "<<m3<<std::endl;
#endif
            MultUU_Helper<algo,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo 96: Transpose
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<96,s,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::const_transpose_type M1t;
            typedef typename M2::const_transpose_type M2t;
            typedef typename M3::transpose_type M3t;
            M1t m1t = m1.transpose();
            M2t m2t = m2.transpose();
            M3t m3t = m3.transpose();
            MultUU_Helper<-2,s,add,ix,T,M2t,M1t,M3t>::call(x,m2t,m1t,m3t);
        }
    };

    // algo 97: Conjugate
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<97,s,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        { 
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            typedef typename M3::conjugate_type M3c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            M3c m3c = m3.conjugate();
            MultUU_Helper<-2,s,add,ix,T,M1c,M2c,M3c>::call(
                TMV_CONJ(x),m1c,m2c,m3c);
        }
    };

    // algo 98: call InstMultMM
    template <int s, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<98,s,false,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typename M3::value_type xx(x);
            InstMultMM(xx,m1.xdView(),m2.xdView(),m3.xdView());
        }
    };
    template <int s, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<98,s,true,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typename M3::value_type xx(x);
            InstAddMultMM(
                xx,m1.xdView(),m2.xdView(),m3.xView().viewAsNonUnitDiag());
        }
    };

    // algo -2: Check for inst
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<-2,s,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            const bool inst = 
                M1::mcolsize == UNKNOWN && M1::mrowsize == UNKNOWN &&
                M2::mcolsize == UNKNOWN && M2::mrowsize == UNKNOWN &&
                M3::mcolsize == UNKNOWN && M3::mrowsize == UNKNOWN &&
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const int algo = 
                s == 0 ? 0 :
                s == 1 ? ( M3::munit ? 0 : 1 ) :
                M3::mconj ? 97 :
                inst ? 98 : 
                -3;
            MultUU_Helper<algo,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo 99: Check for aliases
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<99,s,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            // We set up the algorithm so that m1,m3 can be in the 
            // same storage.
            const bool s1 = 
                SameStorage(m1,m3) && !OppositeStorage(m1,m3);
            const bool s2 = 
                SameStorage(m2,m3) && !OppositeStorage(m2,m3);
#ifdef PRINTALGO_UU
            std::cout<<"UU algo 99:\n";
            std::cout<<"s1 = "<<s1<<std::endl;
            std::cout<<"s2 = "<<s2<<std::endl;
            std::cout<<"Exact(13) = "<<ExactSameStorage(m1,m3);
            std::cout<<"Exact(23) = "<<ExactSameStorage(m2,m3);
#endif
            if ((!s1 || ExactSameStorage(m1,m3)) && !s2) {
#ifdef PRINTALGO_UU
                std::cout<<"No alias\n";
#endif
                // No aliasing (or no clobbering)
                MultUU_Helper<-2,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else if ((!s2 || ExactSameStorage(m2,m3)) && !s1) {
#ifdef PRINTALGO_UU
                std::cout<<"No alias for transpose\n";
#endif
                // Can transpose to get no clobbering storage
                MultUU_Helper<96,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else if (!s2) {
#ifdef PRINTALGO_UU
                std::cout<<"Copy m1 -> m3\n";
#endif
                // Copy m1 to m3, becomes MultEq op
                AliasCopy(m1,m3);
                NoAliasMultMM<add>(x,m3,m2,m3);
            } else if (!s1) {
#ifdef PRINTALGO_UU
                std::cout<<"Copy m2 -> m3\n";
#endif
                // Copy m2 to m3, becomes MultEq op
                AliasCopy(m2,m3);
                typename M3::transpose_type m3t = m3.transpose();
                NoAliasMultMM<add>(x,m3t,m1.transpose(),m3t);
            } else { 
#ifdef PRINTALGO_UU
                std::cout<<"Temp m3\n";
#endif
                // Use temporary for m1*m2
                MultUU_Helper<87,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            }
        }
    };

    // algo -1: Check for aliases?
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUU_Helper<-1,s,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const bool checkalias =
                M1::msize == UNKNOWN && 
                M2::msize == UNKNOWN && 
                M3::msize == UNKNOWN;
            const int algo = 
                s == 0 ? 0 :
                s == 1 ? ( M3::munit ? 0 : 1 ) :
                checkalias ? 99 : 
                -2;
            MultUU_Helper<algo,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    template <int algo, bool add, int ix, class T, class M1, class M2, class M3>
    inline void DoMultUU(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Tri_Mutable<M3>& m3)
    {
        TMVStaticAssert(M1::mupper == int(M2::mupper));
        TMVStaticAssert(M1::mupper == int(M3::mupper));
        TMVStaticAssert((Sizes<M1::msize,M3::msize>::same));
        TMVStaticAssert((Sizes<M1::msize,M2::msize>::same));
        TMVStaticAssert(!M3::munit || ix == 1);
        TMVAssert(m1.size() == m3.size());
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m3.isunit() || ix == 1);

        const int s = Sizes<Sizes<M1::msize,M2::msize>::size,M3::msize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        M1v m1v = m1.cView();
        M2v m2v = m2.cView();
        M3v m3v = m3.cView();
        MultUU_Helper<algo,s,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Tri_Mutable<M3>& m3)
    { DoMultUU<-1,add>(x,m1,m2,m3); }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void NoAliasMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Tri_Mutable<M3>& m3)
    { DoMultUU<-2,add>(x,m1,m2,m3); }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Tri_Mutable<M3>& m3)
    { DoMultUU<-3,add>(x,m1,m2,m3); }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void AliasMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Tri_Mutable<M3>& m3)
    { DoMultUU<99,add>(x,m1,m2,m3); }

    template <class M1, int ix, class T, class M2>
    inline void MultEqMM(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2)
    {
        TMVStaticAssert(M1::mupper == int(M2::mupper));
        TMVStaticAssert(!M1::munit || ix == 1);
        TMVAssert(!m1.isunit() || ix == 1);
        MultMM<false>(x,m1.mat(),m2.mat(),m1.mat());
    }

    template <class M1, int ix, class T, class M2>
    inline void NoAliasMultEqMM(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2)
    {
        TMVStaticAssert(M1::mupper == int(M2::mupper));
        TMVStaticAssert(!M1::munit || ix == 1);
        TMVAssert(!m1.isunit() || ix == 1);
        NoAliasMultMM<false>(x,m1.mat(),m2.mat(),m1.mat());
    }

    template <class M1, int ix, class T, class M2>
    inline void AliasMultEqMM(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2)
    {
        TMVStaticAssert(M1::mupper == int(M2::mupper));
        TMVStaticAssert(!M1::munit || ix == 1);
        TMVAssert(!m1.isunit() || ix == 1);
        AliasMultMM<false>(x,m1.mat(),m2.mat(),m1.mat());
    }

#undef TMV_Q1
#undef TMV_OPT_CLEANUP
#undef TMV_Q2

} // namespace tmv

#endif 
