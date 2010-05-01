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


#ifndef TMV_MultUL_H
#define TMV_MultUL_H

#include "TMV_MultMM.h"
#include "TMV_MultXV.h"
#include "TMV_MultVV.h"
#include "TMV_MultUV.h"
#include "TMV_MultUM.h"

#ifdef PRINTALGO_UL
#include <iostream>
#include "TMV_MatrixIO.h"
#include "TMV_TriMatrixIO.h"
#endif

// Use the specialized 1,2,3 sized algorithms for the end of the 
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

// Q2 is the size to stop recursing.
#define TMV_Q2 8

namespace tmv {

    // Defined below:
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void NoAliasMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void AliasMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);

    // Defined in TMV_MultUL.cpp
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstMultMM(
        const T3 x,
        const ConstUpperTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1, 
        const ConstLowerTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2, 
        MatrixView<T3> m3);
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstAddMultMM(
        const T3 x,
        const ConstUpperTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1, 
        const ConstLowerTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2,
        MatrixView<T3> m3);

    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstMultMM(
        const T3 x,
        const ConstLowerTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1, 
        const ConstUpperTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2, 
        MatrixView<T3> m3);
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstAddMultMM(
        const T3 x,
        const ConstLowerTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1, 
        const ConstUpperTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2,
        MatrixView<T3> m3);

    // Defined in TMV_MultUU.h
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


    //
    // Matrix * Matrix
    //

    template <int algo, int s, bool add, 
              int ix, class T, class M1, class M2, class M3> 
    struct MultUL_Helper;

    // algo 0: Trivial, nothing to do.
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<0,s,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {}
    };

    // algo 1: s == 1, so reduces to scalar product
    template <bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<1,1,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_UL
            std::cout<<"UL algo 1: N,s,x = "<<1<<','<<1<<','<<T(x)<<std::endl;
#endif
            Maybe<add>::add(m3.ref(0,0),x*m1.cref(0,0)*m2.cref(0,0));
        }
    };

    // algo 11: U*L loop over n
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<11,s,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = s==UNKNOWN ? int(m1.size()) : s;
#ifdef PRINTALGO_UL
            std::cout<<"UL algo 11: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            const bool u2 = M2::_unit;
            typedef typename M2::value_type T2;
            typedef typename Maybe<!u2>::template ProdType<T2,T>::type PT2;
            typedef typename M1::const_submatrix_type M1sm;
            typedef typename M1::const_subtrimatrix_type M1st;
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M2::const_col_sub_type M2c;
            typedef typename M3::col_sub_type M3c;
            const int ix2 = u2 ? ix : 0;
            const int xx = UNKNOWN;

            for(int j=0;j<N;++j) {
                // m3.col(j,0,j) = m1.subMatrix(0,j,j,N) * m2.col(j,j,N)
                // m3.col(j,j,N) = m1.subTriMatrix(j,N) * m2.col(j,j,N)
                // ==>
                // m3.col(j,0,j) = m1.col(j,0,j) * m2(j,j) +
                //                 m1.subMatrix(0,j,j+1,N) * m2.col(j,j+1,N)
                // m3(j,j) = m1.row(j,j+1,N) * m2.col(j,j+1,N) + 
                //           m1(j,j) * m2(j,j)
                // m3.col(j,j+1,N) = m1.subTriMatrix(j+1,N) * m2.col(j,j+1,N)
                
                const Scaling<ix2,PT2> xd(Maybe<!u2>::prod(m2.cref(j,j),x));
                M1sm m1m = m1.cSubMatrix(0,j,j+1,N);
                M1st m1t = m1.cSubTriMatrix(j+1,N);
                M1c m1a = m1.get_col(j,0,j);
                M1r m1b = m1.get_row(j,j+1,N);
                M2c m2b = m2.get_col(j,j+1,N);
                M3c m3a = m3.get_col(j,0,j);
                M3c m3b = m3.get_col(j,j+1,N);

                MultXV_Helper<-4,xx,add,ix2,PT2,M1c,M3c>::call(xd,m1a,m3a);
                MultMV_Helper<-4,xx,xx,true,ix,T,M1sm,M2c,M3c>::call(
                    x,m1m,m2b,m3a);
                Maybe<add>::add(
                    m3.ref(j,j) ,
                    x * (MultVV_Helper<-4,xx,M1r,M2c>::call(m1b,m2b) +
                         Maybe<!u2>::prod(m2.cref(j,j) , m1.cref(j,j))));
                MultUV_Helper<-4,xx,add,ix,T,M1st,M2c,M3c>::call(
                    x,m1t,m2b,m3b);
            }
        }
    };

    // algo 12: U*L loop over m
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<12,s,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = s==UNKNOWN ? int(m1.size()) : s;
#ifdef PRINTALGO_UL
            std::cout<<"UL algo 12: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            const bool u1 = M1::_unit;
            typedef typename M1::value_type T1;
            typedef typename Maybe<!u1>::template ProdType<T1,T>::type PT1;
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M2::const_submatrix_type M2sm;
            typedef typename M2sm::const_transpose_type M2smt;
            typedef typename M2::const_subtrimatrix_type M2st;
            typedef typename M2st::const_transpose_type M2stt;
            typedef typename M2::const_row_sub_type M2r;
            typedef typename M2::const_col_sub_type M2c;
            typedef typename M3::row_sub_type M3r;
            const int ix1 = u1 ? ix : 0;
            const int xx = UNKNOWN;

            for(int i=0;i<N;++i) {
                // m3.row(i,0,i) = m1.row(i,i,N) * m2.subMatrix(i,N,0,i)
                // m3.row(i,i,N) = m1.row(i,i,N) * m2.subTriMatrix(i,N)
                // ==>
                // m3.row(i,0,i) = m1(i,i) * m2.row(i,0,i) +
                //                 m1.row(i,i+1,N) * m2.subMatrix(i+1,N,0,i) 
                // m3(i,i) = m1.row(i,i+1,N) * m2.col(i,i+1,N) +
                //           m1(i,i) * m2(i,i)
                // m3.row(i,i+1,N) = m1.row(i,i+1,N) * m2.subTriMatrix(i+1,N) 

                const Scaling<ix1,PT1> xd(Maybe<!u1>::prod(m1.cref(i,i),x));

                M1r m1b = m1.get_row(i,i+1,N);
                M2smt m2m = m2.cSubMatrix(i+1,N,0,i).transpose();
                M2stt m2t = m2.cSubTriMatrix(i+1,N).transpose();
                M2r m2a = m2.get_row(i,0,i);
                M2c m2b = m2.get_col(i,i+1,N);
                M3r m3a = m3.get_row(i,0,i);
                M3r m3b = m3.get_row(i,i+1,N);

                MultXV_Helper<-4,xx,add,ix1,PT1,M2r,M3r>::call(xd,m2a,m3a);
                MultMV_Helper<-4,xx,xx,true,ix,T,M2smt,M1r,M3r>::call(
                    x,m2m,m1b,m3a);
                Maybe<add>::add(
                    m3.ref(i,i) ,
                    x * (MultVV_Helper<-4,xx,M1r,M2c>::call(m1b,m2b) +
                         Maybe<!u1>::prod(m1.cref(i,i) , m2.cref(i,i))));
                MultUV_Helper<-4,xx,add,ix,T,M2stt,M1r,M3r>::call(
                    x,m2t,m1b,m3b);
            }
        }
    };

    // algo 13: U*L loop over k
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<13,s,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = s==UNKNOWN ? int(m1.size()) : s;
#ifdef PRINTALGO_UL
            std::cout<<"UL algo 13: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            const bool u1 = M1::_unit;
            const bool u2 = M2::_unit;
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
                // m3.subMatrix(0,k+1,0,k+1) = 
                //                     m1.col(k,0,k+1) ^ m2.row(k,0,k+1)
                // ==>
                // m3.subMatrix(0,k,0,k) = m1.col(k,0,k) ^ m2.row(k,0,k)
                // m3.col(k,0,k) = m1.col(k,0,k) * m2(k,k)
                // m3.row(k,0,k) = m1(k,k) * m2.row(k,0,k)
                // m3(k,k) = m1(k,k) * m2(k,k)
                M1c m1c = m1.get_col(k,0,k);
                M2r m2r = m2.get_row(k,0,k);
                M3c m3c = m3.get_col(k,0,k);
                M3r m3r = m3.get_row(k,0,k);
                M3s m3s = m3.cSubMatrix(0,k,0,k);
                const Scaling<ix1,PT1> xd1(Maybe<!u1>::prod(m1.cref(k,k),x));
                const Scaling<ix2,PT2> xd2(Maybe<!u2>::prod(m2.cref(k,k),x));
                Rank1VVM_Helper<-4,xx,xx,true,ix,T,M1c,M2r,M3s>::call(
                    x,m1c,m2r,m3s);
                MultXV_Helper<-4,xx,add,ix1,PT1,M2r,M3r>::call(xd1,m2r,m3r);
                MultXV_Helper<-4,xx,add,ix2,PT2,M1c,M3c>::call(xd2,m1c,m3c);
                Maybe<add>::add(
                    m3.ref(k,k),Maybe<!u2>::prod(m2.cref(k,k),xd1));
            }
        }
    };

    // algo 16: Unroll small case
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<16,s,add,ix,T,M1,M2,M3> 
    {
        template <int I, int N>
        struct Unroller
        {
            static void unroll(
                const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
            {
                // [ C00 C01 ] = [ A00 A01 ] [ B00  0  ]
                // [ C10 C11 ]   [  0  A11 ] [ B10 B11 ]
                
                const int Nx = N/2;
                const int Ny = N-Nx;
                const int I1 = I+Nx;
                const int I2 = I+N;
                typedef typename M1::const_submatrix_type M1sm;
                typedef typename M1sm::const_transpose_type M1smt;
                typedef typename M1::const_subtrimatrix_type M1st;
                typedef typename M2::const_submatrix_type M2sm;
                typedef typename M2::const_subtrimatrix_type M2st;
                typedef typename M2st::const_transpose_type M2stt;
                typedef typename M3::submatrix_type M3sm;
                typedef typename M3sm::transpose_type M3smt;

                M1sm A01 = m1.cSubMatrix(I,I1,I1,I2);
                M1smt A01t = A01.transpose();
                M1st A11 = m1.cSubTriMatrix(I1,I2);
                M2sm B10 = m2.cSubMatrix(I1,I2,I,I1);
                M2stt B11t = m2.cSubTriMatrix(I1,I2).transpose();
                M3sm C00 = m3.cSubMatrix(I,I1,I,I1);
                M3smt C01t = m3.cSubMatrix(I,I1,I1,I2).transpose();
                M3sm C10 = m3.cSubMatrix(I1,I2,I,I1);

                // C00 (+)= A00 B00
                Unroller<I,Nx>::unroll(x,m1,m2,m3);

                // C00 += A01 B10
                MultMM_Helper<-4,Nx,Nx,Ny,true,ix,T,M1sm,M2sm,M3sm>::call(
                    x,A01,B10,C00);

                // C10 (+)= A11 B10
                MultUM_Helper<-4,Ny,Nx,add,ix,T,M1st,M2sm,M3sm>::call(
                    x,A11,B10,C10);

                // C01 (+)= A01 B11
                MultUM_Helper<-4,Ny,Nx,add,ix,T,M2stt,M1smt,M3smt>::call(
                    x,B11t,A01t,C01t);

                // C11 = A11 B11
                Unroller<I1,Ny>::unroll(x,m1,m2,m3);
            }
        };
        template <int I>
        struct Unroller<I,1>
        {
            static inline void unroll(
                const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
            {
                const bool u1 = M1::_unit;
                const bool u2 = M2::_unit;
                Maybe<add>::add(
                    m3.ref(I,I), Maybe<!u1>::prod(
                        m1.cref(I,I),Maybe<!u2>::prod(
                            m2.cref(I,I),x)));
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
#ifdef PRINTALGO_UL
            std::cout<<"UL algo 16: N,s,x = "<<2<<','<<2<<','<<T(x)<<std::endl;
#endif
            Unroller<0,s>::unroll(x,m1,m2,m3); 
        }
    };

    template <bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<16,UNKNOWN,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = m1.size();
#ifdef PRINTALGO_UL
            std::cout<<"UL algo 16: N,s,x = "<<N<<','<<UNKNOWN<<
                ','<<T(x)<<std::endl;
#endif
            const bool rxr = M1::_rowmajor && M3::_rowmajor;
            const bool crx = M1::_colmajor && M2::_rowmajor;
            const bool xcc = M2::_colmajor && M3::_colmajor;
            const int algo2 = rxr ? 12 : crx ? 13 : xcc ? 11 : 13;
            switch (N) {
              case 0 :
                   // do nothing
                   break;
              case 1 : 
                   MultUL_Helper<1,1,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
                   break;
              case 2 :
                   MultUL_Helper<16,2,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
                   break;
              case 3 :
                   MultUL_Helper<16,3,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
                   break;
              default :
                   MultUL_Helper<algo2,UNKNOWN,add,ix,T,M1,M2,M3>::call(
                       x,m1,m2,m3);
            }
        }
    };

    // algo 17: Split the TriMatrixes into 3 sections and recurse
    // the calculation on each of them:
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<17,s,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = s==UNKNOWN ? int(m1.size()) : s;
#ifdef PRINTALGO_UL
            std::cout<<"UL algo 17: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif

#if (TMV_Q2 == 1)
            const int algo2 = (s == UNKNOWN || s == 1) ? 1 : 0;
#else
            const int sp1 = IntTraits<s>::Sp1;
            const int twosp1 = IntTraits<IntTraits<s>::twoS>::Sp1;
            // nops = 1/6 n(n+1)(2n+1)
            const int nops = 
                IntTraits2<IntTraits2<s,sp1>::safeprod,twosp1>::safeprod / 6;
            const bool unroll = 
                s > 20 ? false :
                s == UNKNOWN ? false :
                nops > TMV_Q1 ? false :
                s <= 10;
            const bool rxr = M1::_rowmajor && M3::_rowmajor;
            const bool crx = M1::_colmajor && M2::_rowmajor;
            const bool xcc = M2::_colmajor && M3::_colmajor;
            const int algo2 = 
                s == 0 ? 0 :
                s == 1 ? 1 :
                (s != UNKNOWN && s > TMV_Q2) ? 0 :
#ifdef TMV_OPT_CLEANUP
                s == UNKNOWN ? 16 :
#endif
                unroll ? 16 : 
                rxr ? 12 : crx ? 13 : xcc ? 11 : 13;
#endif
            const int algo3 =  // The algorithm for N > Q2
                (s == UNKNOWN || s > TMV_Q2) ? 17 : 0;
            const int algo4 =  // The algorithm for MultMM, MultUM
                s == UNKNOWN ? -2 : s > 16 ? -3 : s > TMV_Q2 ? -4 : 0;

            if (N > TMV_Q2) {
                // [ C00 C01 ] = [ A00 A01 ] [ B00  0  ]
                // [ C10 C11 ]   [  0  A11 ] [ B10 B11 ]

                const int Nx = N > 16 ? ((((N-1)>>5)+1)<<4) : (N>>1);
                // (If N > 16, round N/2 up to a multiple of 16.)
                const int sx = IntTraits<s>::half_roundup;
                const int sy = IntTraits2<s,sx>::diff;

                typedef typename M1::const_subtrimatrix_type M1st;
                typedef typename M1::const_submatrix_type M1sm;
                typedef typename M1sm::const_transpose_type M1smt;
                typedef typename M2::const_subtrimatrix_type M2st;
                typedef typename M2st::const_transpose_type M2stt;
                typedef typename M2::const_submatrix_type M2sm;
                typedef typename M3::submatrix_type M3sm;
                typedef typename M3sm::transpose_type M3smt;

                M1st A00 = m1.cSubTriMatrix(0,Nx);
                M1sm A01 = m1.cSubMatrix(0,Nx,Nx,N);
                M1smt A01t = A01.transpose();
                M1st A11 = m1.cSubTriMatrix(Nx,N);
                M2st B00 = m2.cSubTriMatrix(0,Nx);
                M2sm B10 = m2.cSubMatrix(Nx,N,0,Nx);
                M2st B11 = m2.cSubTriMatrix(Nx,N);
                M2stt B11t = B11.transpose();
                M3sm C00 = m3.cSubMatrix(0,Nx,0,Nx);
                M3sm C01 = m3.cSubMatrix(0,Nx,Nx,N);
                M3smt C01t = C01.transpose();
                M3sm C10 = m3.cSubMatrix(Nx,N,0,Nx);
                M3sm C11 = m3.cSubMatrix(Nx,N,Nx,N);

                // C00 (+)= x A00 B00
                MultUL_Helper<algo3,sx,add,ix,T,M1st,M2st,M3sm>::call(
                    x,A00,B00,C00);

                // C00 += A01 B10
                MultMM_Helper<algo4,sx,sx,sy,true,ix,T,M1sm,M2sm,M3sm>::call(
                    x,A01,B10,C00);

                // We set up the storage so that it's ok if all three
                // matrices are in the same place with the same steps.
                // Or if B is the opposite storage and A is somewhere else.
                // So in the next lines, it is guaranteed that the 
                // C10 assignment doesn't clobber the A01 that is needed
                // next.

                // C10 (+)= A11 B10
                MultUM_Helper<algo4,sy,sx,add,ix,T,M1st,M2sm,M3sm>::call(
                    x,A11,B10,C10);

                // C01 (+)= A01 B11
                MultUM_Helper<algo4,sy,sx,add,ix,T,M2stt,M1smt,M3smt>::
                    call(x,B11t,A01t,C01t);

                // C11 (+)= x A11 B11
                MultUL_Helper<algo3,sy,add,ix,T,M1st,M2st,M3sm>::call(
                    x,A11,B11,C11);
            } else {
                MultUL_Helper<algo2,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            }
        }
    };

    // algo 21: L*U loop over n
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<21,s,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = s==UNKNOWN ? int(m1.size()) : s;
#ifdef PRINTALGO_UL
            std::cout<<"UL algo 21: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            const bool u2 = M2::_unit;
            typedef typename M2::value_type T2;
            typedef typename Maybe<!u2>::template ProdType<T2,T>::type PT2;
            typedef typename M1::const_submatrix_type M1sm;
            typedef typename M1::const_subtrimatrix_type M1st;
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M2::const_col_sub_type M2c;
            typedef typename M3::col_sub_type M3c;
            const int ix2 = u2 ? ix : 0;
            const int xx = UNKNOWN;

            for(int j=N-1;j>=0;--j) {
                // m3.col(j,0,j+1) = m1.subTriMatrix(0,j+1) * m2.col(j,0,j+1)
                // m3.col(j,j+1,N) = m1.subMatrix(j+1,N,0,j+1) * m2.col(j,0,j+1)
                // ==>
                // m3.col(j,j+1,N) = m1.col(j,j+1,N) * m2(j,j) +
                //                   m1.subMatrix(j+1,N,0,j) * m2.col(j,0,j)
                // m3(j,j) = m1.row(j,0,j) * m2.col(j,0,j) + 
                //           m1(j,j) * m2(j,j)
                // m3.col(j,0,j) = m1.subTriMatrix(0,j) * m2.col(j,0,j)
                
                const Scaling<ix2,PT2> xd(Maybe<!u2>::prod(m2.cref(j,j),x));
                M1st m1t = m1.cSubTriMatrix(0,j);
                M1sm m1m = m1.cSubMatrix(j+1,N,0,j);
                M1r m1a = m1.get_row(j,0,j);
                M1c m1b = m1.get_col(j,j+1,N);
                M2c m2a = m2.get_col(j,0,j);
                M3c m3a = m3.get_col(j,0,j);
                M3c m3b = m3.get_col(j,j+1,N);

                MultXV_Helper<-4,xx,add,ix2,PT2,M1c,M3c>::call(xd,m1b,m3b);
                MultMV_Helper<-4,xx,xx,true,ix,T,M1sm,M2c,M3c>::call(
                    x,m1m,m2a,m3b);
                Maybe<add>::add(
                    m3.ref(j,j) ,
                    x * (MultVV_Helper<-4,xx,M1r,M2c>::call(m1a,m2a) +
                         Maybe<!u2>::prod(m2.cref(j,j) , m1.cref(j,j))));
                MultUV_Helper<-4,xx,add,ix,T,M1st,M2c,M3c>::call(
                    x,m1t,m2a,m3a);
            }
        }
    };

    // algo 22: L*U loop over m
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<22,s,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = s==UNKNOWN ? int(m1.size()) : s;
#ifdef PRINTALGO_UL
            std::cout<<"UL algo 22: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            const bool u1 = M1::_unit;
            typedef typename M1::value_type T1;
            typedef typename Maybe<!u1>::template ProdType<T1,T>::type PT1;
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M2::const_submatrix_type M2sm;
            typedef typename M2sm::const_transpose_type M2smt;
            typedef typename M2::const_subtrimatrix_type M2st;
            typedef typename M2st::const_transpose_type M2stt;
            typedef typename M2::const_row_sub_type M2r;
            typedef typename M2::const_col_sub_type M2c;
            typedef typename M3::row_sub_type M3r;
            const int ix1 = u1 ? ix : 0;
            const int xx = UNKNOWN;

            for(int i=N-1;i>=0;--i) {
                // m3.row(i,0,i+1) = m1.row(i,0,i+1) * m2.subTriMatrix(0,i+1)
                // m3.row(i,i+1,N) = m1.row(i,0,i+1) * m2.subMatrix(0,i+1,i+1,N)
                // ==>
                // m3.row(i,i+1,N) = m1(i,i) * m2.row(i,i+1,N) +
                //                   m1.row(i,0,i) * m2.subMatrix(0,i,i+1,N) 
                // m3(i,i) = m1.row(i,0,i) * m2.col(i,0,i) +
                //           m1(i,i) * m2(i,i)
                // m3.row(i,0,i) = m1.row(i,0,i) * m2.subTriMatrix(0,i) 

                const Scaling<ix1,PT1> xd(Maybe<!u1>::prod(m1.cref(i,i),x));

                M1r m1a = m1.get_row(i,0,i);
                M2stt m2t = m2.cSubTriMatrix(0,i).transpose();
                M2smt m2m = m2.cSubMatrix(0,i,i+1,N).transpose();
                M2r m2b = m2.get_row(i,i+1,N);
                M2c m2a = m2.get_col(i,0,i);
                M3r m3a = m3.get_row(i,0,i);
                M3r m3b = m3.get_row(i,i+1,N);

                MultXV_Helper<-4,xx,add,ix1,PT1,M2r,M3r>::call(xd,m2b,m3b);
                MultMV_Helper<-4,xx,xx,true,ix,T,M2smt,M1r,M3r>::call(
                    x,m2m,m1a,m3b);
                Maybe<add>::add(
                    m3.ref(i,i) ,
                    x * (MultVV_Helper<-4,xx,M1r,M2c>::call(m1a,m2a) +
                         Maybe<!u1>::prod(m1.cref(i,i) , m2.cref(i,i))));
                MultUV_Helper<-4,xx,add,ix,T,M2stt,M1r,M3r>::call(
                    x,m2t,m1a,m3a);
            }
        }
    };

    // algo 23: L*U loop over k
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<23,s,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = s==UNKNOWN ? int(m1.size()) : s;
#ifdef PRINTALGO_UL
            std::cout<<"UL algo 23: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            const bool u1 = M1::_unit;
            const bool u2 = M2::_unit;
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
                // m3.subMatrix(k,N,k,N) = m1.col(k,k,N) ^ m2.row(k,k,N)
                // ==>
                // m3.subMatrix(k+1,N,k+1,N) = m1.col(k,k+1,N) ^ m2.row(k,k+1,N)
                // m3.col(k,k+1,N) = m1.col(k,k+1,N) * m2(k,k)
                // m3.row(k,k+1,N) = m1(k,k) * m2.row(k,k+1,N)
                // m3(k,k) = m1(k,k) * m2(k,k)
                M1c m1c = m1.get_col(k,k+1,N);
                M2r m2r = m2.get_row(k,k+1,N);
                M3c m3c = m3.get_col(k,k+1,N);
                M3r m3r = m3.get_row(k,k+1,N);
                M3s m3s = m3.cSubMatrix(k+1,N,k+1,N);
                const Scaling<ix1,PT1> xd1(Maybe<!u1>::prod(m1.cref(k,k),x));
                const Scaling<ix2,PT2> xd2(Maybe<!u2>::prod(m2.cref(k,k),x));
                Rank1VVM_Helper<-4,xx,xx,true,ix,T,M1c,M2r,M3s>::call(
                    x,m1c,m2r,m3s);
                MultXV_Helper<-4,xx,add,ix1,PT1,M2r,M3r>::call(xd1,m2r,m3r);
                MultXV_Helper<-4,xx,add,ix2,PT2,M1c,M3c>::call(xd2,m1c,m3c);
                Maybe<add>::add(
                    m3.ref(k,k),Maybe<!u2>::prod(m2.cref(k,k),xd1));
            }
        }
    };

    // algo 26: Unroll small case
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<26,s,add,ix,T,M1,M2,M3> 
    {
        template <int I, int N>
        struct Unroller
        {
            static void unroll(
                const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
            {
                // [ C00 C01 ] = [ A00  0  ] [ B00 B01 ]
                // [ C10 C11 ]   [ A10 A11 ] [  0  B11 ]
                
                const int Nx = N/2;
                const int Ny = N-Nx;
                const int I1 = I+Nx;
                const int I2 = I+N;
                typedef typename M1::const_subtrimatrix_type M1st;
                typedef typename M1::const_submatrix_type M1sm;
                typedef typename M1sm::const_transpose_type M1smt;
                typedef typename M2::const_subtrimatrix_type M2st;
                typedef typename M2::const_submatrix_type M2sm;
                typedef typename M2st::const_transpose_type M2stt;
                typedef typename M3::submatrix_type M3sm;
                typedef typename M3sm::transpose_type M3smt;

                M1st A00 = m1.cSubTriMatrix(I,I1);
                M1sm A10 = m1.cSubMatrix(I1,I2,I,I1);
                M1smt A10t = A10.transpose();
                M2stt B00t = m2.cSubTriMatrix(I,I1).transpose();
                M2sm B01 = m2.cSubMatrix(I,I1,I1,I2);
                M3sm C01 = m3.cSubMatrix(I,I1,I1,I2);
                M3smt C10t = m3.cSubMatrix(I1,I2,I,I1).transpose();
                M3sm C11 = m3.cSubMatrix(I1,I2,I1,I2);

                // C11 = A11 B11
                Unroller<I1,Ny>::unroll(x,m1,m2,m3);

                // C11 += A10 B01
                MultMM_Helper<-4,Ny,Ny,Nx,true,ix,T,M1sm,M2sm,M3sm>::call(
                    x,A10,B01,C11);

                // C01 (+)= A00 B01
                MultUM_Helper<-4,Nx,Ny,add,ix,T,M1st,M2sm,M3sm>::call(
                    x,A00,B01,C01);

                // C10 (+)= A10 B00
                MultUM_Helper<-4,Nx,Ny,add,ix,T,M2stt,M1smt,M3smt>::call(
                    x,B00t,A10t,C10t);

                // C00 (+)= A00 B00
                Unroller<I,Nx>::unroll(x,m1,m2,m3);
            }
        };
        template <int I>
        struct Unroller<I,1>
        {
            static inline void unroll(
                const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
            {
                const bool u1 = M1::_unit;
                const bool u2 = M2::_unit;
                Maybe<add>::add(
                    m3.ref(I,I), Maybe<!u1>::prod(
                        m1.cref(I,I),Maybe<!u2>::prod(
                            m2.cref(I,I),x)));
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
#ifdef PRINTALGO_UL
            std::cout<<"UL algo 26: N,s,x = "<<2<<','<<2<<','<<T(x)<<std::endl;
#endif
            Unroller<0,s>::unroll(x,m1,m2,m3); 
        }
    };
    template <bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<26,UNKNOWN,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = m1.size();
#ifdef PRINTALGO_UL
            std::cout<<"UL algo 26: N,s,x = "<<N<<','<<UNKNOWN<<
                ','<<T(x)<<std::endl;
#endif
            const bool rxr = M1::_rowmajor && M3::_rowmajor;
            const bool crx = M1::_colmajor && M2::_rowmajor;
            const bool xcc = M2::_colmajor && M3::_colmajor;
            const int algo2 = rxr ? 22 : crx ? 23 : xcc ? 21 : 23;
            switch (N) {
              case 0 :
                   // do nothing
                   break;
              case 1 :
                   MultUL_Helper<1,1,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
                   break;
              case 2 :
                   MultUL_Helper<26,2,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
                   break;
              case 3 :
                   MultUL_Helper<26,3,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
                   break;
              default :
                   MultUL_Helper<algo2,UNKNOWN,add,ix,T,M1,M2,M3>::call(
                       x,m1,m2,m3);
            }
        }
    };

    // algo 27: Split the TriMatrixes into 3 sections and recurse
    // the calculation on each of them:
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<27,s,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = s==UNKNOWN ? int(m1.size()) : s;
#ifdef PRINTALGO_UL
            std::cout<<"UL algo 27: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif

#if (TMV_Q2 == 1)
            const int algo2 = (s == UNKNOWN || s == 1) ? 1 : 0;
#else
            const int sp1 = IntTraits<s>::Sp1;
            const int twosp1 = IntTraits<IntTraits<s>::twoS>::Sp1;
            // nops = 1/6 n(n+1)(2n+1)
            const int nops = 
                IntTraits2<IntTraits2<s,sp1>::safeprod,twosp1>::safeprod / 6;
            const bool unroll = 
                s > 20 ? false :
                nops > TMV_Q1 ? false :
                s <= 10;
            const bool rxr = M1::_rowmajor && M3::_rowmajor;
            const bool crx = M1::_colmajor && M2::_rowmajor;
            const bool xcc = M2::_colmajor && M3::_colmajor;
            const int algo2 = 
                s == 0 ? 0 :
                s == 1 ? 1 :
                (s != UNKNOWN && s > TMV_Q2) ? 0 :
#ifdef TMV_OPT_CLEANUP
                s == UNKNOWN ? 26 :
#endif
                unroll ? 26 : 
                rxr ? 22 : crx ? 23 : xcc ? 21 : 23;
#endif
            const int algo3 =  // The algorithm for N > Q2
                (s == UNKNOWN || s > TMV_Q2) ? 27 : 0;
            const int algo4 =  // The algorithm for MultMM, MultUM
                s == UNKNOWN ? -2 : s > 16 ? -3 : s > TMV_Q2 ? -4 : 0;

            if (N > TMV_Q2) {
                // [ C00 C01 ] = [ A00  0  ] [ B00 B01 ]
                // [ C10 C11 ]   [ A10 A11 ] [  0  B11 ]

                const int Nx = N > 16 ? ((((N-1)>>5)+1)<<4) : (N>>1);
                // (If N > 16, round N/2 up to a multiple of 16.)
                const int sx = IntTraits<s>::half_roundup;
                const int sy = IntTraits2<s,sx>::diff;

                typedef typename M1::const_subtrimatrix_type M1st;
                typedef typename M1::const_submatrix_type M1sm;
                typedef typename M1sm::const_transpose_type M1smt;
                typedef typename M2::const_subtrimatrix_type M2st;
                typedef typename M2st::const_transpose_type M2stt;
                typedef typename M2::const_submatrix_type M2sm;
                typedef typename M3::submatrix_type M3sm;
                typedef typename M3sm::transpose_type M3smt;

                M1st A00 = m1.cSubTriMatrix(0,Nx);
                M1sm A10 = m1.cSubMatrix(Nx,N,0,Nx);
                M1smt A10t = A10.transpose();
                M1st A11 = m1.cSubTriMatrix(Nx,N);
                M2st B00 = m2.cSubTriMatrix(0,Nx);
                M2stt B00t = B00.transpose();
                M2sm B01 = m2.cSubMatrix(0,Nx,Nx,N);
                M2st B11 = m2.cSubTriMatrix(Nx,N);
                M3sm C00 = m3.cSubMatrix(0,Nx,0,Nx);
                M3sm C01 = m3.cSubMatrix(0,Nx,Nx,N);
                M3smt C10t = m3.cSubMatrix(Nx,N,0,Nx).transpose();
                M3sm C11 = m3.cSubMatrix(Nx,N,Nx,N);

                // C11 (+)= x A11 B11
                MultUL_Helper<algo3,sy,add,ix,T,M1st,M2st,M3sm>::call(
                    x,A11,B11,C11);

                // C11 += A10 B01
                MultMM_Helper<algo4,sy,sy,sx,true,ix,T,M1sm,M2sm,M3sm>::call(
                    x,A10,B01,C11);

                // C01 (+)= A00 B01
                MultUM_Helper<algo4,sx,sy,add,ix,T,M1st,M2sm,M3sm>::call(
                    x,A00,B01,C01);

                // C10 (+)= A10 B00
                MultUM_Helper<algo4,sx,sy,add,ix,T,M2stt,M1smt,M3smt>::
                    call(x,B00t,A10t,C10t);

                // C00 (+)= x A00 B00
                MultUL_Helper<algo3,sx,add,ix,T,M1st,M2st,M3sm>::call(
                    x,A00,B00,C00);
            } else {
                MultUL_Helper<algo2,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            }
        }
    };

    // algo 87: Use temporary for m1*m2
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<87,s,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = s==UNKNOWN ? int(m1.size()) : s;
#ifdef PRINTALGO_UL
            std::cout<<"UL algo 87: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            typedef typename Traits<T>::real_type RT;
            const Scaling<1,RT> one;
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename Traits2<T1,T2>::type PT3;
            const bool rm = M1::_rowmajor && M2::_rowmajor;
            typedef typename MCopyHelper<PT3,Rec,s,s,rm,false>::type M3c;
            M3c m3c(N,N);
            typedef typename M3c::view_type M3cv;
            typedef typename M3c::const_view_type M3ccv;
            M3cv m3cv = m3c.view();
            M3ccv m3ccv = m3c.view();
            MultUL_Helper<-2,s,false,1,RT,M1,M2,M3cv>::call(one,m1,m2,m3cv);
            MultXM_Helper<-2,s,s,add,ix,T,M3ccv,M3>::call(x,m3ccv,m3);
        }
    };

    // algo 88: Copy m2 to the same storage as m3
    template <int s, int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<88,s,false,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = s==UNKNOWN ? int(m1.size()) : s;
#ifdef PRINTALGO_UL
            std::cout<<"UL algo 88: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            typedef typename TypeSelect<M2::_upper,
                    typename M3::uppertri_type,
                    typename M3::lowertri_type>::type M3u;
            M3u m3u = Maybe<M2::_upper>::uppertri(m3);
            typedef typename M3u::const_view_type M3ucv;
            M3ucv m3ucv = m3u.view();
            CopyU_Helper<-2,s,M2,M3u>::call(m2,m3u);
            MultUL_Helper<-2,s,false,ix,T,M1,M3ucv,M3>::call(x,m1,m3ucv,m3);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<-3,s,add,ix,T,M1,M2,M3> 
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            TMVStaticAssert(M1::_upper != int(M2::_upper));
            
            // Possible algorithms are:
            //
            // Trivial and special for small TriMatrix sizes.
            //  0 = s == 0, so nothing to do
            //  1 = s == 1: reduces to trivial scalar product function
            //
            // U * L:
            // 11 = loop over n: MultUV
            // 12 = loop over m: MultUV
            // 13 = loop over k: Rank1
            // 16 = Unroll small case
            // 17 = split each trimatrix into 3 submatrices and recurse
            //
            // L * U:
            // 21 = loop over n: MultUV
            // 22 = loop over m: MultUV
            // 23 = loop over k: Rank1
            // 26 = Unroll small case
            // 27 = split each trimatrix into 3 submatrices and recurse
            //
            // 87 = Use a temporary for m1*m2
            // 88 = Copy m2 to the same storage as m3

            const bool ul = M1::_upper;
#if 0
            const int algo = 27;
#else
#if TMV_OPT == 0 
            const bool xcc = M2::_colmajor && M3::_colmajor;
            const bool rxr = M1::_rowmajor && M3::_rowmajor;
            const bool crx = M1::_colmajor && M2::_rowmajor;
            const int algo =
                ul ? ( rxr ? 12 : crx ? 13 : xcc ? 11 : 13 ) :
                ( rxr ? 22 : crx ? 23 : xcc ? 22 : 23 );
#else
            const int s2 = s > 20 ? UNKNOWN : s;
            const int s2p1 = IntTraits<s2>::Sp1;
            const int twos2p1 = IntTraits<IntTraits<s2>::twoS>::Sp1;
            // nops = 1/6 n(n+1)(2n+1)
            const int nops = 
                IntTraits2<IntTraits2<s2,s2p1>::safeprod,twos2p1>::safeprod / 6;
            const bool unroll = 
                s == UNKNOWN ? false :
                nops > TMV_Q1 ? false :
                s <= 10;
            const int algo = 
                s == 0 ? 0 :
                s == 1 ? 1 :
                unroll ? ( ul ? 16 : 26 ) :
                ul ? 17 : 27;
#endif
#endif
#ifdef PRINTALGO_UL
            const int N = s==UNKNOWN ? int(m1.size()) : s;
            std::cout<<"InlineMultUL: x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"N = "<<N<<std::endl;
            std::cout<<"s = "<<s<<std::endl;
            std::cout<<"add = "<<add<<", algo = "<<algo<<std::endl;
#endif
            MultUL_Helper<algo,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo 94: U*U, so convert to MultUU call
    template <int s, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<94,s,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_UL
            std::cout<<"UL algo 94\n";
#endif
            typename M3::uppertri_type m3u = m3.upperTri();
            NoAliasMultMM<add>(x,m1,m2,m3u);
            Maybe<!add>::zero2(m3.lowerTri());
        }
    };

    // algo 95: L*L, so convert to MultUU call
    template <int s, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<95,s,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_UL
            std::cout<<"UL algo 95\n";
#endif
            typename M3::lowertri_type m3l = m3.lowerTri();
            NoAliasMultMM<add>(x,m1,m2,m3l);
            Maybe<!add>::zero2(m3.upperTri());
        }
    };

    // algo 96: Transpose
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<96,s,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        { 
#ifdef PRINTALGO_UL
            std::cout<<"UL algo 96\n";
#endif
            typedef typename M1::const_transpose_type M1t;
            typedef typename M2::const_transpose_type M2t;
            typedef typename M3::transpose_type M3t;
            M1t m1t = m1.transpose();
            M2t m2t = m2.transpose();
            M3t m3t = m3.transpose();
            MultUL_Helper<-2,s,add,ix,T,M1t,M2t,M3t>::call(x,m1t,m2t,m3t);
        }
    };

    // algo 97: Conjugate
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<97,s,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        { 
#ifdef PRINTALGO_UL
            std::cout<<"UL algo 97\n";
#endif
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            typedef typename M3::conjugate_type M3c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            M3c m3c = m3.conjugate();
            MultUL_Helper<-2,s,add,ix,T,M1c,M2c,M3c>::call(
                TMV_CONJ(x),m1c,m2c,m3c);
        }
    };

    // algo 98: call InstMultMM
    template <int s, int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<98,s,false,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_UL
            std::cout<<"UL algo 98\n";
#endif
            typename M3::value_type xx(x);
            InstMultMM(xx,m1.xdView(),m2.xdView(),m3.xView());
        }
    };
    template <int s, int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<98,s,true,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_UL
            std::cout<<"UL algo 98\n";
#endif
            typename M3::value_type xx(x);
            InstAddMultMM(xx,m1.xdView(),m2.xdView(),m3.xView());
        }
    };

    // algo -2: Check for inst
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<-2,s,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            const bool inst = 
                M1::unknownsizes &&
                M2::unknownsizes &&
                M3::unknownsizes &&
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
                s == 1 ? 1 :
                M1::_upper && M2::_upper ? 94 :
                !M1::_upper && !M2::_upper ? 95 :
                M3::_conj ? 97 :
                inst ? 98 : 
                -3;
#ifdef PRINTALGO_UL
            std::cout<<"UL algo 98: algo = "<<algo<<std::endl;
#endif
            MultUL_Helper<algo,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo 940: U*U, so convert to MultUU call (with alias check)
    template <int s, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<940,s,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_UL
            std::cout<<"UL algo 940\n";
#endif
            typename M3::uppertri_type m3u = m3.upperTri();
            MultMM<add>(x,m1,m2,m3u);
            typename M3::lowertri_type::offdiag_type m3lo = 
                m3.lowerTri().offDiag();
            Maybe<!add>::zero(m3lo);
        }
    };

    // algo 950: L*L, so convert to MultUU call (with alias check)
    template <int s, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<950,s,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_UL
            std::cout<<"UL algo 950\n";
#endif
            typename M3::lowertri_type m3l = m3.lowerTri();
            MultMM<add>(x,m1,m2,m3l);
            typename M3::uppertri_type::offdiag_type m3uo = 
                m3.upperTri().offDiag();
            Maybe<!add>::zero(m3uo);
        }
    };

    // algo 99: Check for aliases
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<99,s,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_UL
            std::cout<<"UL algo 99\n";
#endif
            // We set up the algorithm so that it's ok if all three
            // matrices are in the same place with the same steps.
            // Or if m1 is the same storage and m2 is the opposite storage
            // (which means that it is really in the same place as m1).
            // Or if m2 is the opposite storage and m1 is somewhere else.
            const bool s1 = SameStorage(m1,m3);
            const bool s2 = SameStorage(m2,m3);
            const bool ok1 = !s1 || ExactSameStorage(m1,m3);
            const bool ok2 = !s2 || 
                ExactSameStorage(m2,m3) || OppositeStorage(m2,m3);
            if ( ok1 && ok2 ) {
                // No aliasing (or no clobbering)
                MultUL_Helper<-2,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else if ( ok2 && OppositeStorage(m1,m3) ) {
                // Can transpose to get no clobbering storage
                MultUL_Helper<96,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else { 
                // Use temporary for m1*m2
                MultUL_Helper<87,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            }
        }
    };

    // algo -1: Check for aliases?
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUL_Helper<-1,s,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const bool checkalias =
                M1::_size == UNKNOWN &&
                M2::_size == UNKNOWN &&
                M3::_colsize == UNKNOWN && M3::_rowsize == UNKNOWN;
            const int algo = 
                s == 0 ? 0 :
                s == 1 ? 1 :
                M1::_upper && M2::_upper ? 940 :
                !M1::_upper && !M2::_upper ? 950 :
                checkalias ? 99 : 
                -2;
#ifdef PRINTALGO_UL
            std::cout<<"UL algo -1: algo = "<<algo<<std::endl;
#endif
            MultUL_Helper<algo,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    template <int algo, bool add, int ix, class T, class M1, class M2, class M3>
    inline void DoMultUL(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
        TMVStaticAssert((Sizes<M1::_size,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_size,M3::_rowsize>::same));
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.size() == m3.colsize());
        TMVAssert(m1.size() == m3.rowsize());

        const int s1 = Sizes<M1::_size,M2::_size>::size;
        const int s2 = Sizes<M3::_colsize,M3::_rowsize>::size;
        const int s = Sizes<s1,s2>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        M1v m1v = m1.cView();
        M2v m2v = m2.cView();
        M3v m3v = m3.cView();
        MultUL_Helper<algo,s,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    { DoMultUL<-1,add>(x,m1,m2,m3); }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void NoAliasMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    { DoMultUL<-2,add>(x,m1,m2,m3); }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    { DoMultUL<-3,add>(x,m1,m2,m3); }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void AliasMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    { DoMultUL<99,add>(x,m1,m2,m3); }


} // namespace tmv

#undef TMV_OPT_CLEANUP
#undef TMV_Q1
#undef TMV_Q2

#endif 
