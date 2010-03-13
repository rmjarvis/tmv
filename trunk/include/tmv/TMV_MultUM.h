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


#ifndef TMV_MultUM_H
#define TMV_MultUM_H

#include "TMV_MultMM.h"
#include "TMV_MultUV.h"
#include "TMV_CopyU.h"
#include "TMV_MultXU.h"

#ifdef _OPENMP
#include "omp.h"
#ifdef PRINTALGO_UM_OMP
#include <fstream>
#endif
#endif

#ifdef PRINTALGO_UM
#include <iostream>
#endif

// Use the specialized 1,2,3 sized algorithms for the end of the 
// recursive algorithm.
#if TMV_OPT >= 2
#define TMV_OPT_CLEANUP
#endif

// Check for small (<=5) values of cs or rs
// This leads to significant speed improvements for such matrices at little
// cost to larger matrices, but the large increase in code size and the fact
// that it only benefits a few particular sizes of matrices mean that
// we require TMV_OPT = 3 for it.
#if TMV_OPT >= 3
#define TMV_OPT_SMALL
#endif

// Q2 is the minimum size to keep recursing
#define TMV_Q2 8

// Q6 is the minimum value of (M^2*N / 16^3) to use multiple threads.
// (We also require that N > 64 so we have something to split.)
#define TMV_Q6 64

namespace tmv {

    // Defined below:
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Rec<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void NoAliasMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Rec<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Rec<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void AliasMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Rec<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void NoAliasMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void AliasMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <class M1, int ix, class T, class M2>
    inline void MultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2);
    template <class M1, int ix, class T, class M2>
    inline void NoAliasMultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2);
    template <class M1, int ix, class T, class M2>
    inline void AliasMultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2);

    // Defined in TMV_MultUM.cpp
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstMultMM(
        const T3 x,
        const ConstUpperTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1, 
        const ConstMatrixView<T2,UNKNOWN,UNKNOWN,C2>& m2, MatrixView<T3> m3);
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstAddMultMM(
        const T3 x,
        const ConstUpperTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1, 
        const ConstMatrixView<T2,UNKNOWN,UNKNOWN,C2>& m2, MatrixView<T3> m3);

    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstMultMM(
        const T3 x,
        const ConstLowerTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1, 
        const ConstMatrixView<T2,UNKNOWN,UNKNOWN,C2>& m2, MatrixView<T3> m3);
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstAddMultMM(
        const T3 x,
        const ConstLowerTriMatrixView<T1,UnknownDiag,UNKNOWN,UNKNOWN,C1>& m1, 
        const ConstMatrixView<T2,UNKNOWN,UNKNOWN,C2>& m2, MatrixView<T3> m3);


    //
    // Matrix * Matrix
    //

    template <int algo, int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3> 
    struct MultUM_Helper;

    // algo 0: Trivial, nothing to do.
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<0,cs,rs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {}
    };

    // algo 1: cs == 1, so reduces to MultXV
    template <int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<1,1,rs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_UM
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
            std::cout<<"UM algo 1: M,N,cs,rs,x = "<<1<<','<<N<<
                ','<<1<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename Traits2<T,typename M1::value_type>::type PT1;
            typedef typename M2::const_row_type M2r;
            typedef typename M3::row_type M3r;

            M2r m2r = m2.get_row(0);
            M3r m3r = m3.get_row(0);
            MultXV_Helper<-1,rs,add,0,PT1,M2r,M3r>::call(
                Scaling<0,PT1>(x*m1.cref(0,0)), m2r, m3r);
        }
    };

    // algo 2: rs == 1, so reduces to MultUV
    template <int cs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<2,cs,1,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_UM
            const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
            std::cout<<"UM algo 2: M,N,cs,rs,x = "<<M<<','<<1<<
                ','<<cs<<','<<1<<','<<T(x)<<std::endl;
#endif
            typedef typename M2::const_col_type M2c;
            typedef typename M3::col_type M3c;

            M2c m2c = m2.get_col(0);
            M3c m3c = m3.get_col(0);
            MultUV_Helper<-1,cs,add,ix,T,M1,M2c,M3c>::call(x,m1,m2c,m3c);
        }
    };

    // algo 201: same as 1, but use -2 algo
    template <int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<201,1,rs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_UM
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
            std::cout<<"UM algo 201: M,N,cs,rs,x = "<<1<<','<<N<<
                ','<<1<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename Traits2<T,typename M1::value_type>::type PT1;
            typedef typename M2::const_row_type M2r;
            typedef typename M3::row_type M3r;

            M2r m2r = m2.get_row(0);
            M3r m3r = m3.get_row(0);
            MultXV_Helper<-2,rs,add,0,PT1,M2r,M3r>::call(
                Scaling<0,PT1>(x*m1.cref(0,0)), m2r, m3r);
        }
    };

    // algo 202: same as 2, but use -2 algo
    template <int cs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<202,cs,1,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_UM
            const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
            std::cout<<"UM algo 202: M,N,cs,rs,x = "<<M<<','<<1<<
                ','<<cs<<','<<1<<','<<T(x)<<std::endl;
#endif
            typedef typename M2::const_col_type M2c;
            typedef typename M3::col_type M3c;

            M2c m2c = m2.get_col(0);
            M3c m3c = m3.get_col(0);
            MultUV_Helper<-2,cs,add,ix,T,M1,M2c,M3c>::call(x,m1,m2c,m3c);
        }
    };

    // algo 401: same as 1, but use -4 algo
    template <int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<401,1,rs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_UM
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
            std::cout<<"UM algo 401: M,N,cs,rs,x = "<<1<<','<<N<<
                ','<<1<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename Traits2<T,typename M1::value_type>::type PT1;
            typedef typename M2::const_row_type M2r;
            typedef typename M3::row_type M3r;

            M2r m2r = m2.get_row(0);
            M3r m3r = m3.get_row(0);
            MultXV_Helper<-4,rs,add,0,PT1,M2r,M3r>::call(
                Scaling<0,PT1>(x*m1.cref(0,0)), m2r, m3r);
        }
    };

    // algo 402: same as 2, but use -4 algo
    template <int cs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<402,cs,1,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_UM
            const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
            std::cout<<"UM algo 402: M,N,cs,rs,x = "<<M<<','<<1<<
                ','<<cs<<','<<1<<','<<T(x)<<std::endl;
#endif
            typedef typename M2::const_col_type M2c;
            typedef typename M3::col_type M3c;

            M2c m2c = m2.get_col(0);
            M3c m3c = m3.get_col(0);
            MultUV_Helper<-4,cs,add,ix,T,M1,M2c,M3c>::call(x,m1,m2c,m3c);
        }
    };

    // algo 11: UpperTri loop over n
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<11,cs,rs,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
#ifdef PRINTALGO_UM
            const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
            std::cout<<"UM algo 11: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename M2::const_col_type M2c;
            typedef typename M3::col_type M3c;

            for(int j=0;j<N;++j) {
                // m3.col(j) = m1 * m2.col(j)
                M2c m2j = m2.get_col(j);
                M3c m3j = m3.get_col(j);
                MultUV_Helper<-4,cs,add,ix,T,M1,M2c,M3c>::call(x,m1,m2j,m3j);
            }
        }
    };

    // algo 12: UpperTri loop over m
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<12,cs,rs,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
#ifdef PRINTALGO_UM
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
            std::cout<<"UM algo 12: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            const bool unit = M1::munit;
            typedef typename M1::value_type T1;
            typedef typename Maybe<!unit>::template ProdType<T1,T>::type PT1;
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M2::const_row_type M2r;
            typedef typename M2::const_rowrange_type M2rr;
            typedef typename M2rr::const_transpose_type M2t;
            typedef typename M3::row_type M3r;
            const int ix1 = unit ? ix : 0;
            const int xx = UNKNOWN;

            for(int i=0;i<M;++i) {
                // m3.row(i) = m1.row(i,i,M) * m2.rowRange(i,M)
                //           = m1(i,i) * m2.row(i) +
                //             m1.row(i,i+1,M) * m2.rowRange(i+1,M)
                M1r m1i = m1.get_row(i,i+1,M);
                M2t m2t = m2.cRowRange(i+1,M).transpose();
                M2r m2i = m2.get_row(i);
                M3r m3i = m3.get_row(i);
                const Scaling<ix1,PT1> xdi(Maybe<!unit>::prod(m1.cref(i,i),x));
                MultXV_Helper<-4,rs,add,ix1,PT1,M2r,M3r>::call(xdi,m2i,m3i);
                MultMV_Helper<-4,rs,xx,true,ix,T,M2t,M1r,M3r>::call(
                    x,m2t,m1i,m3i);
            }
        }
    };

    // algo 13: UpperTri loop over k
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<13,cs,rs,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
#ifdef PRINTALGO_UM
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
            std::cout<<"UM algo 13: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            const bool unit = M1::munit;
            typedef typename M1::value_type T1;
            typedef typename Maybe<!unit>::template ProdType<T1,T>::type PT1;
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::const_row_type M2r;
            typedef typename M3::row_type M3r;
            typedef typename M3::rowrange_type M3rr;
            const int ix1 = unit ? ix : 0;
            const int xx = UNKNOWN;

            for(int k=0;k<M;++k) {
                // m3.rowRange(0,k+1) = m1.col(k,0,k+1) ^ m2.row(k)
                // ==> 
                // m3.rowRange(0,k) += m1.col(k,0,k) ^ m2.row(k)
                // m3.row(k) = m1(k,k) * m2.row(k)
                M1c m1k = m1.get_col(k,0,k);
                M2r m2k = m2.get_row(k);
                M3r m3k = m3.get_row(k);
                M3rr m3r = m3.cRowRange(0,k);
                const Scaling<ix1,PT1> xdk(Maybe<!unit>::prod(m1.cref(k,k),x));
                Rank1VVM_Helper<-4,xx,rs,true,ix,T,M1c,M2r,M3rr>::call(
                    x,m1k,m2k,m3r);
                MultXV_Helper<-4,rs,add,ix1,PT1,M2r,M3r>::call(xdk,m2k,m3k);
            }
        }
    };

    // algo 16: Specialization for cs = 2,3,4,5
    template <int rs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<16,2,rs,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
#ifdef PRINTALGO_UM
            std::cout<<"UM algo 16: M,N,cs,rs,x = "<<2<<','<<N<<
                ','<<2<<','<<rs<<','<<T(x)<<std::endl;
#endif
            const bool unit = M1::munit;
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            typedef typename Maybe<!unit>::template ProdType<T1,T>::type PT1d;
            typedef typename Traits2<T,T1>::type PT1;
            typedef typename M2::const_row_type M2r;
            typedef typename M2r::const_iterator IT2;
            typedef typename M3::row_type M3r;
            typedef typename M3r::iterator IT3;

            T2 b0, b1;
            T3 c0, c1;

            const int ix1 = unit ? ix : 0;
            const Scaling<ix1,PT1d> a00(Maybe<!unit>::prod(m1.cref(0,0),x));
            const Scaling<0,PT1> a01(x * m1.cref(0,1));
            const Scaling<ix1,PT1d> a11(Maybe<!unit>::prod(m1.cref(1,1),x));

            const int Bstepi = m2.stepi();
            const int Cstepi = m3.stepi();
            IT2 B0 = m2.get_row(0).begin();;
            IT2 B1 = B0; B1.shiftP(Bstepi);
            IT3 C0 = m3.get_row(0).begin();;
            IT3 C1 = C0; C1.shiftP(Cstepi);

            Prefetch_Read(B0.getP());
            Prefetch_Read(B1.getP());
            Prefetch_Write(C0.getP());
            Prefetch_Write(C1.getP());

            int j;

            j = N; do {
                b0 = *B0++; b1 = *B1++;
                c0 = a00 * b0 + a01 * b1;
                Maybe<add>::add(*C0++ , c0);
                c1 = a11 * b1;
                Maybe<add>::add(*C1++ , c1);
            } while (--j);
        }
    };
    template <int rs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<16,3,rs,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
#ifdef PRINTALGO_UM
            std::cout<<"UM algo 16: M,N,cs,rs,x = "<<3<<','<<N<<
                ','<<3<<','<<rs<<','<<T(x)<<std::endl;
#endif
            const bool unit = M1::munit;
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            typedef typename Maybe<!unit>::template ProdType<T1,T>::type PT1d;
            typedef typename Traits2<T,T1>::type PT1;
            typedef typename M2::const_row_type M2r;
            typedef typename M2r::const_iterator IT2;
            typedef typename M3::row_type M3r;
            typedef typename M3r::iterator IT3;

            T2 b0, b1, b2;
            T3 c0, c1, c2;

            const int ix1 = unit ? ix : 0;
            const Scaling<ix1,PT1d> a00(Maybe<!unit>::prod(m1.cref(0,0),x));
            const Scaling<0,PT1> a01(x * m1.cref(0,1));
            const Scaling<0,PT1> a02(x * m1.cref(0,2));
            const Scaling<ix1,PT1d> a11(Maybe<!unit>::prod(m1.cref(1,1),x));
            const Scaling<0,PT1> a12(x * m1.cref(1,2));
            const Scaling<ix1,PT1d> a22(Maybe<!unit>::prod(m1.cref(2,2),x));

            const int Bstepi = m2.stepi();
            const int Cstepi = m3.stepi();
            IT2 B0 = m2.get_row(0).begin();;
            IT2 B1 = B0; B1.shiftP(Bstepi);
            IT2 B2 = B1; B2.shiftP(Bstepi);
            IT3 C0 = m3.get_row(0).begin();;
            IT3 C1 = C0; C1.shiftP(Cstepi);
            IT3 C2 = C1; C2.shiftP(Cstepi);

            Prefetch_Read(B0.getP());
            Prefetch_Read(B1.getP());
            Prefetch_Read(B2.getP());
            Prefetch_Write(C0.getP());
            Prefetch_Write(C1.getP());
            Prefetch_Write(C2.getP());

            int j;

            j = N; do {
                b0 = *B0++; b1 = *B1++; b2 = *B2++;
                c0 = a00 * b0 + a01 * b1 + a02 * b2;
                Maybe<add>::add(*C0++ , c0);
                c1 = a11 * b1 + a12 * b2;
                Maybe<add>::add(*C1++ , c1);
                c2 = a22 * b2;
                Maybe<add>::add(*C2++ , c2);
            } while (--j);
        }
    };
    template <int rs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<16,4,rs,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
#ifdef PRINTALGO_UM
            std::cout<<"UM algo 16: M,N,cs,rs,x = "<<4<<','<<N<<
                ','<<4<<','<<rs<<','<<T(x)<<std::endl;
#endif
            const bool unit = M1::munit;
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            typedef typename Maybe<!unit>::template ProdType<T1,T>::type PT1d;
            typedef typename Traits2<T,T1>::type PT1;
            typedef typename M2::const_row_type M2r;
            typedef typename M2r::const_iterator IT2;
            typedef typename M3::row_type M3r;
            typedef typename M3r::iterator IT3;

            T2 b0, b1, b2, b3;
            T3 c0, c1, c2, c3;

            const int ix1 = unit ? ix : 0;
            const Scaling<ix1,PT1d> a00(Maybe<!unit>::prod(m1.cref(0,0),x));
            const Scaling<0,PT1> a01(x * m1.cref(0,1));
            const Scaling<0,PT1> a02(x * m1.cref(0,2));
            const Scaling<0,PT1> a03(x * m1.cref(0,3));
            const Scaling<ix1,PT1d> a11(Maybe<!unit>::prod(m1.cref(1,1),x));
            const Scaling<0,PT1> a12(x * m1.cref(1,2));
            const Scaling<0,PT1> a13(x * m1.cref(1,3));
            const Scaling<ix1,PT1d> a22(Maybe<!unit>::prod(m1.cref(2,2),x));
            const Scaling<0,PT1> a23(x * m1.cref(2,3));
            const Scaling<ix1,PT1d> a33(Maybe<!unit>::prod(m1.cref(3,3),x));

            const int Bstepi = m2.stepi();
            const int Cstepi = m3.stepi();
            IT2 B0 = m2.get_row(0).begin();;
            IT2 B1 = B0; B1.shiftP(Bstepi);
            IT2 B2 = B1; B2.shiftP(Bstepi);
            IT2 B3 = B2; B3.shiftP(Bstepi);
            IT3 C0 = m3.get_row(0).begin();;
            IT3 C1 = C0; C1.shiftP(Cstepi);
            IT3 C2 = C1; C2.shiftP(Cstepi);
            IT3 C3 = C2; C3.shiftP(Cstepi);

            Prefetch_Read(B0.getP());
            Prefetch_Read(B1.getP());
            Prefetch_Read(B2.getP());
            Prefetch_Read(B3.getP());
            Prefetch_Write(C0.getP());
            Prefetch_Write(C1.getP());
            Prefetch_Write(C2.getP());
            Prefetch_Write(C3.getP());

            int j;

            j = N; do {
                b0 = *B0++; b1 = *B1++; b2 = *B2++; b3 = *B3++;
                c0 = a00 * b0 + a01 * b1 + a02 * b2 + a03 * b3;
                Maybe<add>::add(*C0++ , c0);
                c1 = a11 * b1 + a12 * b2 + a13 * b3;
                Maybe<add>::add(*C1++ , c1);
                c2 = a22 * b2 + a23 * b3;
                Maybe<add>::add(*C2++ , c2);
                c3 = a33 * b3;
                Maybe<add>::add(*C3++ , c3);
            } while (--j);
        }
    };
    template <int rs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<16,5,rs,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
#ifdef PRINTALGO_UM
            std::cout<<"UM algo 16: M,N,cs,rs,x = "<<5<<','<<N<<
                ','<<5<<','<<rs<<','<<T(x)<<std::endl;
#endif
            const bool unit = M1::munit;
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            typedef typename Maybe<!unit>::template ProdType<T1,T>::type PT1d;
            typedef typename Traits2<T,T1>::type PT1;
            typedef typename M2::const_row_type M2r;
            typedef typename M2r::const_iterator IT2;
            typedef typename M3::row_type M3r;
            typedef typename M3r::iterator IT3;

            T2 b0, b1, b2, b3, b4;
            T3 c0, c1, c2, c3, c4;

            const int ix1 = unit ? ix : 0;
            const Scaling<ix1,PT1d> a00(Maybe<!unit>::prod(m1.cref(0,0),x));
            const Scaling<0,PT1> a01(x * m1.cref(0,1));
            const Scaling<0,PT1> a02(x * m1.cref(0,2));
            const Scaling<0,PT1> a03(x * m1.cref(0,3));
            const Scaling<0,PT1> a04(x * m1.cref(0,4));
            const Scaling<ix1,PT1d> a11(Maybe<!unit>::prod(m1.cref(1,1),x));
            const Scaling<0,PT1> a12(x * m1.cref(1,2));
            const Scaling<0,PT1> a13(x * m1.cref(1,3));
            const Scaling<0,PT1> a14(x * m1.cref(1,4));
            const Scaling<ix1,PT1d> a22(Maybe<!unit>::prod(m1.cref(2,2),x));
            const Scaling<0,PT1> a23(x * m1.cref(2,3));
            const Scaling<0,PT1> a24(x * m1.cref(2,4));
            const Scaling<ix1,PT1d> a33(Maybe<!unit>::prod(m1.cref(3,3),x));
            const Scaling<0,PT1> a34(x * m1.cref(3,4));
            const Scaling<ix1,PT1d> a44(Maybe<!unit>::prod(m1.cref(4,4),x));

            const int Bstepi = m2.stepi();
            const int Cstepi = m3.stepi();
            IT2 B0 = m2.get_row(0).begin();;
            IT2 B1 = B0; B1.shiftP(Bstepi);
            IT2 B2 = B1; B2.shiftP(Bstepi);
            IT2 B3 = B2; B3.shiftP(Bstepi);
            IT2 B4 = B3; B4.shiftP(Bstepi);
            IT3 C0 = m3.get_row(0).begin();;
            IT3 C1 = C0; C1.shiftP(Cstepi);
            IT3 C2 = C1; C2.shiftP(Cstepi);
            IT3 C3 = C2; C3.shiftP(Cstepi);
            IT3 C4 = C3; C4.shiftP(Cstepi);

            Prefetch_Read(B0.getP());
            Prefetch_Read(B1.getP());
            Prefetch_Read(B2.getP());
            Prefetch_Read(B3.getP());
            Prefetch_Read(B4.getP());
            Prefetch_Write(C0.getP());
            Prefetch_Write(C1.getP());
            Prefetch_Write(C2.getP());
            Prefetch_Write(C3.getP());
            Prefetch_Write(C4.getP());

            int j;

            j = N; do {
                b0 = *B0++; b1 = *B1++; b2 = *B2++; b3 = *B3++; b4 = *B4++;
                c0 = a00 * b0 + a01 * b1 + a02 * b2 + a03 * b3 + a04 * b4;
                Maybe<add>::add(*C0++ , c0);
                c1 = a11 * b1 + a12 * b2 + a13 * b3 + a14 * b4;
                Maybe<add>::add(*C1++ , c1);
                c2 = a22 * b2 + a23 * b3 + a24 * b4;
                Maybe<add>::add(*C2++ , c2);
                c3 = a33 * b3 + a34 * b4;
                Maybe<add>::add(*C3++ , c3);
                c4 = a44 * b4;
                Maybe<add>::add(*C4++ , c4);
            } while (--j);
        }
    };
    template <int rs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<16,UNKNOWN,rs,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = m3.colsize();
#ifdef PRINTALGO_UM
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
            std::cout<<"UM algo 16: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<UNKNOWN<<','<<rs<<','<<T(x)<<std::endl;
#endif
            const bool rxr = M1::mrowmajor && M3::mrowmajor;
            const bool crx = M1::mcolmajor && M2::mrowmajor;
            const bool xcc = M2::mcolmajor && M3::mcolmajor;
            const int algo2 = rxr ? 12 : crx ? 13 : xcc ? 11 : 13;
            switch (M) {
              case 0 :
                   // do nothing
                   break;
              case 1 :
                   MultUM_Helper<1,1,rs,add,ix,T,M1,M2,M3>::call(
                       x,m1,m2,m3);
                   break;
              case 2 :
                   MultUM_Helper<16,2,rs,add,ix,T,M1,M2,M3>::call(
                       x,m1,m2,m3);
                   break;
              case 3 :
                   MultUM_Helper<16,3,rs,add,ix,T,M1,M2,M3>::call(
                       x,m1,m2,m3);
                   break;
              default :
                   MultUM_Helper<algo2,UNKNOWN,rs,add,ix,T,M1,M2,M3>::call(
                       x,m1,m2,m3);
            }
        }
    };

    // algo 17: Split the UpperTriMatrix into 3 sections and recurse
    // the calculation on each of them:
    // ( A B ) ( D ) = ( F )
    // ( 0 C ) ( E )   ( G )
    // F = AD + BE
    // G = CE
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<17,cs,rs,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
#ifdef PRINTALGO_UM
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
            std::cout<<"UM algo 17: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif

#if (TMV_Q2 == 1)
            const int algo2 = (cs == UNKNOWN || cs == 1) ? 1 : 0;
#else
            const bool rxr = M1::mrowmajor && M3::mrowmajor;
            const bool crx = M1::mcolmajor && M2::mrowmajor;
            const bool xcc = M2::mcolmajor && M3::mcolmajor;
            const int algo2 = 
                cs == 0 ? 0 :
                cs == 1 ? 1 :
                (cs != UNKNOWN && cs > TMV_Q2) ? 0 :
#ifdef TMV_OPT_CLEANUP
                cs == UNKNOWN ? 16 :
#endif
                cs <= 5 ? 16 :
                rxr ? 12 : crx ? 13 : xcc ? 11 : 13;
#endif
            const int algo3 =  // The algorithm for M > Q2
                (cs == UNKNOWN || cs > TMV_Q2) ? 17 : 0;
            const int algo4 =  // The algorithm for MultMM
                cs == UNKNOWN ? -2 : cs > 16 ? -3 : cs > TMV_Q2 ? -4 : 0;

            if (M > TMV_Q2) {
                const int Mx = M > 16 ? ((((M-1)>>5)+1)<<4) : (M>>1);
                // (If M > 16, round M/2 up to a multiple of 16.)
                const int csx = IntTraits<cs>::half_roundup;
                const int csy = cs == UNKNOWN ? cs : cs-csx;

                typedef typename M1::const_subtrimatrix_type M1a;
                typedef typename M1::const_submatrix_type M1b;
                typedef typename M2::const_rowrange_type M2r;
                typedef typename M3::rowrange_type M3r;

                M1a A = m1.cSubTriMatrix(0,Mx);
                M1b B = m1.cSubMatrix(0,Mx,Mx,M);
                M1a C = m1.cSubTriMatrix(Mx,M);
                M2r D = m2.cRowRange(0,Mx);
                M2r E = m2.cRowRange(Mx,M);
                M3r F = m3.cRowRange(0,Mx);
                M3r G = m3.cRowRange(Mx,M);

                MultUM_Helper<algo3,csx,rs,add,ix,T,M1a,M2r,M3r>::call(x,A,D,F);
                MultMM_Helper<algo4,csx,rs,csy,true,ix,T,M1b,M2r,M3r>::call(
                    x,B,E,F);
                MultUM_Helper<algo3,csy,rs,add,ix,T,M1a,M2r,M3r>::call(x,C,E,G);
            } else {
                MultUM_Helper<algo2,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            }
        }
    };

    // algo 21: LowerTri loop over n
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<21,cs,rs,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
#ifdef PRINTALGO_UM
            const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
            std::cout<<"UM algo 21: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename M2::const_col_type M2c;
            typedef typename M3::col_type M3c;

            for(int j=0;j<N;++j) {
                // m3.col(j) = m1 * m2.col(j)
                M2c m2j = m2.get_col(j);
                M3c m3j = m3.get_col(j);
                MultUV_Helper<-4,cs,add,ix,T,M1,M2c,M3c>::call(x,m1,m2j,m3j);
            }
        }
    };

    // algo 22: LowerTri loop over m
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<22,cs,rs,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
#ifdef PRINTALGO_UM
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
            std::cout<<"UM algo 22: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            const bool unit = M1::munit;
            typedef typename M1::value_type T1;
            typedef typename Maybe<!unit>::template ProdType<T1,T>::type PT1;
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M2::const_row_type M2r;
            typedef typename M2::const_rowrange_type M2rr;
            typedef typename M2rr::const_transpose_type M2t;
            typedef typename M3::row_type M3r;
            const int ix1 = unit ? ix : 0;
            const int xx = UNKNOWN;

            for(int i=M-1;i>=0;--i) {
                // m3.row(i) = m1.row(i,0,i+1) * m2.rowRange(0,i+1)
                //           = m1(i,i) * m2.row(i) +
                //             m1.row(i,0,i) * m2.rowRange(0,i)
                M1r m1i = m1.get_row(i,0,i);
                M2t m2t = m2.cRowRange(0,i).transpose();
                M2r m2i = m2.get_row(i);
                M3r m3i = m3.get_row(i);
                const Scaling<ix1,PT1> xdi(Maybe<!unit>::prod(m1.cref(i,i),x));
                MultXV_Helper<-4,rs,add,ix1,PT1,M2r,M3r>::call(xdi,m2i,m3i);
                MultMV_Helper<-4,rs,xx,true,ix,T,M2t,M1r,M3r>::call(
                    x,m2t,m1i,m3i);
            }
        }
    };

    // algo 23: LowerTri loop over k
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<23,cs,rs,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
#ifdef PRINTALGO_UM
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
            std::cout<<"UM algo 23: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            const bool unit = M1::munit;
            typedef typename M1::value_type T1;
            typedef typename Maybe<!unit>::template ProdType<T1,T>::type PT1;
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::const_row_type M2r;
            typedef typename M3::row_type M3r;
            typedef typename M3::rowrange_type M3rr;
            const int ix1 = unit ? ix : 0;
            const int xx = UNKNOWN;

            for(int k=M-1;k>=0;--k) {
                // m3.rowRange(k,M) = m1.col(k,k,M) ^ m2.row(k)
                // ==> m3.row(k) = m1(k,k) * m2.row(k)
                //     m3.rowRange(k+1,M) += m1.col(k,k+1,M) ^ m2.row(k)
                M1c m1k = m1.get_col(k,k+1,M);
                M2r m2k = m2.get_row(k);
                M3r m3k = m3.get_row(k);
                M3rr m3r = m3.cRowRange(k+1,M);
                const Scaling<ix1,PT1> xdk(Maybe<!unit>::prod(m1.cref(k,k),x));
                Rank1VVM_Helper<-4,xx,rs,true,ix,T,M1c,M2r,M3rr>::call(
                    x,m1k,m2k,m3r);
                MultXV_Helper<-4,rs,add,ix1,PT1,M2r,M3r>::call(xdk,m2k,m3k);
            }
        }
    };

    // algo 26: Specialization for cs = 2,3,4,5
    template <int rs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<26,2,rs,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
#ifdef PRINTALGO_UM
            std::cout<<"UM algo 26: M,N,cs,rs,x = "<<2<<','<<N<<
                ','<<2<<','<<rs<<','<<T(x)<<std::endl;
#endif
            const bool unit = M1::munit;
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            typedef typename Maybe<!unit>::template ProdType<T1,T>::type PT1d;
            typedef typename Traits2<T,T1>::type PT1;
            typedef typename M2::const_row_type M2r;
            typedef typename M2r::const_iterator IT2;
            typedef typename M3::row_type M3r;
            typedef typename M3r::iterator IT3;

            T2 b0, b1;
            T3 c0, c1;

            const int ix1 = unit ? ix : 0;
            const Scaling<ix1,PT1d> a00(Maybe<!unit>::prod(m1.cref(0,0),x));
            const Scaling<0,PT1> a10(x * m1.cref(1,0));
            const Scaling<ix1,PT1d> a11(Maybe<!unit>::prod(m1.cref(1,1),x));

            const int Bstepi = m2.stepi();
            const int Cstepi = m3.stepi();
            IT2 B0 = m2.get_row(0).begin();;
            IT2 B1 = B0; B1.shiftP(Bstepi);
            IT3 C0 = m3.get_row(0).begin();;
            IT3 C1 = C0; C1.shiftP(Cstepi);

            Prefetch_Read(B0.getP());
            Prefetch_Read(B1.getP());
            Prefetch_Write(C0.getP());
            Prefetch_Write(C1.getP());

            int j;

            j = N; do {
                b0 = *B0++; b1 = *B1++;
                c1 = a10 * b0 + a11 * b1;
                Maybe<add>::add(*C1++ , c1);
                c0 = a00 * b0;
                Maybe<add>::add(*C0++ , c0);
            } while (--j);
        }
    };
    template <int rs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<26,3,rs,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
#ifdef PRINTALGO_UM
            std::cout<<"UM algo 26: M,N,cs,rs,x = "<<3<<','<<N<<
                ','<<3<<','<<rs<<','<<T(x)<<std::endl;
#endif
            const bool unit = M1::munit;
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            typedef typename Maybe<!unit>::template ProdType<T1,T>::type PT1d;
            typedef typename Traits2<T,T1>::type PT1;
            typedef typename M2::const_row_type M2r;
            typedef typename M2r::const_iterator IT2;
            typedef typename M3::row_type M3r;
            typedef typename M3r::iterator IT3;

            T2 b0, b1, b2;
            T3 c0, c1, c2;

            const int ix1 = unit ? ix : 0;
            const Scaling<ix1,PT1d> a00(Maybe<!unit>::prod(m1.cref(0,0),x));
            const Scaling<0,PT1> a10(x * m1.cref(1,0));
            const Scaling<0,PT1> a20(x * m1.cref(2,0));
            const Scaling<ix1,PT1d> a11(Maybe<!unit>::prod(m1.cref(1,1),x));
            const Scaling<0,PT1> a21(x * m1.cref(2,1));
            const Scaling<ix1,PT1d> a22(Maybe<!unit>::prod(m1.cref(2,2),x));

            const int Bstepi = m2.stepi();
            const int Cstepi = m3.stepi();
            IT2 B0 = m2.get_row(0).begin();;
            IT2 B1 = B0; B1.shiftP(Bstepi);
            IT2 B2 = B1; B2.shiftP(Bstepi);
            IT3 C0 = m3.get_row(0).begin();;
            IT3 C1 = C0; C1.shiftP(Cstepi);
            IT3 C2 = C1; C2.shiftP(Cstepi);

            Prefetch_Read(B0.getP());
            Prefetch_Read(B1.getP());
            Prefetch_Read(B2.getP());
            Prefetch_Write(C0.getP());
            Prefetch_Write(C1.getP());
            Prefetch_Write(C2.getP());

            int j;

            j = N; do {
                b0 = *B0++; b1 = *B1++; b2 = *B2++;
                c2 = a20 * b0 + a21 * b1 + a22 * b2;
                Maybe<add>::add(*C2++ , c2);
                c1 = a10 * b0 + a11 * b1;
                Maybe<add>::add(*C1++ , c1);
                c0 = a00 * b0;
                Maybe<add>::add(*C0++ , c0);
            } while (--j);
        }
    };
    template <int rs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<26,4,rs,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
#ifdef PRINTALGO_UM
            std::cout<<"UM algo 26: M,N,cs,rs,x = "<<4<<','<<N<<
                ','<<4<<','<<rs<<','<<T(x)<<std::endl;
#endif
            const bool unit = M1::munit;
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            typedef typename Maybe<!unit>::template ProdType<T1,T>::type PT1d;
            typedef typename Traits2<T,T1>::type PT1;
            typedef typename M2::const_row_type M2r;
            typedef typename M2r::const_iterator IT2;
            typedef typename M3::row_type M3r;
            typedef typename M3r::iterator IT3;

            T2 b0, b1, b2, b3;
            T3 c0, c1, c2, c3;

            const int ix1 = unit ? ix : 0;
            const Scaling<ix1,PT1d> a00(Maybe<!unit>::prod(m1.cref(0,0),x));
            const Scaling<0,PT1> a10(x * m1.cref(1,0));
            const Scaling<0,PT1> a20(x * m1.cref(2,0));
            const Scaling<0,PT1> a30(x * m1.cref(3,0));
            const Scaling<ix1,PT1d> a11(Maybe<!unit>::prod(m1.cref(1,1),x));
            const Scaling<0,PT1> a21(x * m1.cref(2,1));
            const Scaling<0,PT1> a31(x * m1.cref(3,1));
            const Scaling<ix1,PT1d> a22(Maybe<!unit>::prod(m1.cref(2,2),x));
            const Scaling<0,PT1> a32(x * m1.cref(3,2));
            const Scaling<ix1,PT1d> a33(Maybe<!unit>::prod(m1.cref(3,3),x));

            const int Bstepi = m2.stepi();
            const int Cstepi = m3.stepi();
            IT2 B0 = m2.get_row(0).begin();;
            IT2 B1 = B0; B1.shiftP(Bstepi);
            IT2 B2 = B1; B2.shiftP(Bstepi);
            IT2 B3 = B2; B3.shiftP(Bstepi);
            IT3 C0 = m3.get_row(0).begin();;
            IT3 C1 = C0; C1.shiftP(Cstepi);
            IT3 C2 = C1; C2.shiftP(Cstepi);
            IT3 C3 = C2; C3.shiftP(Cstepi);

            Prefetch_Read(B0.getP());
            Prefetch_Read(B1.getP());
            Prefetch_Read(B2.getP());
            Prefetch_Read(B3.getP());
            Prefetch_Write(C0.getP());
            Prefetch_Write(C1.getP());
            Prefetch_Write(C2.getP());
            Prefetch_Write(C3.getP());

            int j;

            j = N; do {
                b0 = *B0++; b1 = *B1++; b2 = *B2++; b3 = *B3++;
                c3 = a30 * b0 + a31 * b1 + a32 * b2 + a33 * b3;
                Maybe<add>::add(*C3++ , c3);
                c2 = a20 * b0 + a21 * b1 + a22 * b2;
                Maybe<add>::add(*C2++ , c2);
                c1 = a10 * b0 + a11 * b1;
                Maybe<add>::add(*C1++ , c1);
                c0 = a00 * b0;
                Maybe<add>::add(*C0++ , c0);
            } while (--j);
        }
    };
    template <int rs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<26,5,rs,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
#ifdef PRINTALGO_UM
            std::cout<<"UM algo 26: M,N,cs,rs,x = "<<5<<','<<N<<
                ','<<5<<','<<rs<<','<<T(x)<<std::endl;
#endif
            const bool unit = M1::munit;
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            typedef typename Maybe<!unit>::template ProdType<T1,T>::type PT1d;
            typedef typename Traits2<T,T1>::type PT1;
            typedef typename M2::const_row_type M2r;
            typedef typename M2r::const_iterator IT2;
            typedef typename M3::row_type M3r;
            typedef typename M3r::iterator IT3;

            T2 b0, b1, b2, b3, b4;
            T3 c0, c1, c2, c3, c4;

            const int ix1 = unit ? ix : 0;
            const Scaling<ix1,PT1d> a00(Maybe<!unit>::prod(m1.cref(0,0),x));
            const Scaling<0,PT1> a10(x * m1.cref(1,0));
            const Scaling<0,PT1> a20(x * m1.cref(2,0));
            const Scaling<0,PT1> a30(x * m1.cref(3,0));
            const Scaling<0,PT1> a40(x * m1.cref(4,0));
            const Scaling<ix1,PT1d> a11(Maybe<!unit>::prod(m1.cref(1,1),x));
            const Scaling<0,PT1> a21(x * m1.cref(2,1));
            const Scaling<0,PT1> a31(x * m1.cref(3,1));
            const Scaling<0,PT1> a41(x * m1.cref(4,1));
            const Scaling<ix1,PT1d> a22(Maybe<!unit>::prod(m1.cref(2,2),x));
            const Scaling<0,PT1> a32(x * m1.cref(3,2));
            const Scaling<0,PT1> a42(x * m1.cref(4,2));
            const Scaling<ix1,PT1d> a33(Maybe<!unit>::prod(m1.cref(3,3),x));
            const Scaling<0,PT1> a43(x * m1.cref(4,3));
            const Scaling<ix1,PT1d> a44(Maybe<!unit>::prod(m1.cref(4,4),x));

            const int Bstepi = m2.stepi();
            const int Cstepi = m3.stepi();
            IT2 B0 = m2.get_row(0).begin();;
            IT2 B1 = B0; B1.shiftP(Bstepi);
            IT2 B2 = B1; B2.shiftP(Bstepi);
            IT2 B3 = B2; B3.shiftP(Bstepi);
            IT2 B4 = B3; B4.shiftP(Bstepi);
            IT3 C0 = m3.get_row(0).begin();;
            IT3 C1 = C0; C1.shiftP(Cstepi);
            IT3 C2 = C1; C2.shiftP(Cstepi);
            IT3 C3 = C2; C3.shiftP(Cstepi);
            IT3 C4 = C3; C4.shiftP(Cstepi);

            Prefetch_Read(B0.getP());
            Prefetch_Read(B1.getP());
            Prefetch_Read(B2.getP());
            Prefetch_Read(B3.getP());
            Prefetch_Read(B4.getP());
            Prefetch_Write(C0.getP());
            Prefetch_Write(C1.getP());
            Prefetch_Write(C2.getP());
            Prefetch_Write(C3.getP());
            Prefetch_Write(C4.getP());

            int j;

            j = N; do {
                b0 = *B0++; b1 = *B1++; b2 = *B2++; b3 = *B3++; b4 = *B4++;
                c4 = a40 * b0 + a41 * b1 + a42 * b2 + a43 * b3 + a44 * b4;
                Maybe<add>::add(*C4++ , c4);
                c3 = a30 * b0 + a31 * b1 + a32 * b2 + a33 * b3;
                Maybe<add>::add(*C3++ , c3);
                c2 = a20 * b0 + a21 * b1 + a22 * b2;
                Maybe<add>::add(*C2++ , c2);
                c1 = a10 * b0 + a11 * b1;
                Maybe<add>::add(*C1++ , c1);
                c0 = a00 * b0;
                Maybe<add>::add(*C0++ , c0);
            } while (--j);
        }
    };
    template <int rs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<26,UNKNOWN,rs,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = m3.colsize();
#ifdef PRINTALGO_UM
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
            std::cout<<"UM algo 26: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<UNKNOWN<<','<<rs<<','<<T(x)<<std::endl;
#endif
            const bool rxr = M1::mrowmajor && M3::mrowmajor;
            const bool crx = M1::mcolmajor && M2::mrowmajor;
            const bool xcc = M2::mcolmajor && M3::mcolmajor;
            const int algo2 = rxr ? 22 : crx ? 23 : xcc ? 21 : 23;
            switch (M) {
              case 0 :
                   // do nothing
                   break;
              case 1 :
                   MultUM_Helper<1,1,rs,add,ix,T,M1,M2,M3>::call(
                       x,m1,m2,m3);
                   break;
              case 2 :
                   MultUM_Helper<26,2,rs,add,ix,T,M1,M2,M3>::call(
                       x,m1,m2,m3);
                   break;
              case 3 :
                   MultUM_Helper<26,3,rs,add,ix,T,M1,M2,M3>::call(
                       x,m1,m2,m3);
                   break;
              default :
                   MultUM_Helper<algo2,UNKNOWN,rs,add,ix,T,M1,M2,M3>::call(
                       x,m1,m2,m3);
            }
        }
    };

    // algo 27: Split the LowerTriMatrix into 3 sections and recurse
    // the calculation on each of them:
    // ( A 0 ) ( D ) = ( F )
    // ( B C ) ( E )   ( G )
    // G = BD + CE
    // F = AD
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<27,cs,rs,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
#ifdef PRINTALGO_UM
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
            std::cout<<"UM algo 27: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif

#if (TMV_Q2 == 1)
            const int algo2 = (cs == UNKNOWN || cs == 1) ? 1 : 0;
#else
            const bool rxr = M1::mrowmajor && M3::mrowmajor;
            const bool crx = M1::mcolmajor && M2::mrowmajor;
            const bool xcc = M2::mcolmajor && M3::mcolmajor;
            const int algo2 = 
                cs == 0 ? 0 :
                cs == 1 ? 1 :
                (cs != UNKNOWN && cs > TMV_Q2) ? 0 :
#ifdef TMV_OPT_CLEANUP
                cs == UNKNOWN ? 26 :
#endif
                cs <= 5 ? 26 :
                rxr ? 22 : crx ? 23 : xcc ? 21 : 23;
#endif
            const int algo3 =  // The algorithm for M > Q2
                (cs == UNKNOWN || cs > TMV_Q2) ? 27 : 0;
            const int algo4 =  // The algorithm for MultMM
                cs == UNKNOWN ? -2 : cs > 16 ? -3 : cs > TMV_Q2 ? -4 : 0;

            if (M > TMV_Q2) {
                const int Mx = M > 16 ? ((((M-1)>>5)+1)<<4) : (M>>1);
                // (If M > 16, round M/2 up to a multiple of 16.)
                const int csx = IntTraits<cs>::half_roundup;
                const int csy = cs == UNKNOWN ? cs : cs-csx;

                typedef typename M1::const_subtrimatrix_type M1a;
                typedef typename M1::const_submatrix_type M1b;
                typedef typename M2::const_rowrange_type M2r;
                typedef typename M3::rowrange_type M3r;

                M1a A = m1.cSubTriMatrix(0,Mx);
                M1b B = m1.cSubMatrix(Mx,M,0,Mx);
                M1a C = m1.cSubTriMatrix(Mx,M);
                M2r D = m2.cRowRange(0,Mx);
                M2r E = m2.cRowRange(Mx,M);
                M3r F = m3.cRowRange(0,Mx);
                M3r G = m3.cRowRange(Mx,M);

                MultUM_Helper<algo3,csy,rs,add,ix,T,M1a,M2r,M3r>::call(x,C,E,G);
                MultMM_Helper<algo4,csy,rs,csx,true,ix,T,M1b,M2r,M3r>::call(
                    x,B,D,G);
                MultUM_Helper<algo3,csx,rs,add,ix,T,M1a,M2r,M3r>::call(x,A,D,F);
            } else {
                MultUM_Helper<algo2,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            }
        }
    };

    // algo 31: Determine which algorithm to use based on the runtime
    // knowledge of the sizes.
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<31,cs,rs,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
#ifdef PRINTALGO_UM
            std::cout<<"UM algo 31: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            const bool rxr = M1::mrowmajor && M3::mrowmajor;
            const bool crx = M1::mcolmajor && M2::mrowmajor;
            const bool xcc = M2::mcolmajor && M3::mcolmajor;
            const bool nxx = !(M1::mcolmajor || M1::mrowmajor);
            const bool xxn = !(M3::mcolmajor || M3::mrowmajor);
            // xcc, rxr, crx have fast algorithms already
            // add can't be turned into a MultEq op
            // nxx, xxn don't benefit from the change
            //   (and would go into an infinite loop)
            const bool multeq = !(xcc || rxr || crx || add || nxx || xxn);
            const bool maybesmall = 
                cs==UNKNOWN || (cs<=5 && (rs==UNKNOWN || rs<16)) || cs<=3;
            const bool upper1 = M1::mupper;
#ifdef TMV_OPT_CLEANUP
            const int algo2 = maybesmall ? ( upper1 ? 16 : 26 ) : 0;
#else
            const int algo2 = maybesmall ? (
                upper1 ?  ( rxr ? 12 : crx ? 13 : xcc ? 11 : 13 ) :
                ( rxr ? 22 : crx ? 23 : xcc ? 21 : 23 ) ) :
                0; 
#endif
            // This next trick just avoids compiling algo 33 for 
            // cases where multeq = false, since algo3 is only called
            // when multeq is true.
            const int algo3 = multeq ? 33 : 0;
            const int algo4 = 
                (rs == UNKNOWN || rs <= 3) ? (
                    upper1 ? 11 : 21 ) :
                0;
            const int algo5 = upper1 ? 17 : 27;

#ifdef _OPENMP
            const int Mc = M < 16 ? 1 : M>>4; // M/16
            const int Nc = N < 16 ? 1 : N>>4; // N/16
#endif

            // Put the small matrix option first, so it doesn't have to 
            // go through a bunch of if/else statements.  For large matrices,
            // all these if/else's don't matter for the total time.
            if ((M <= 5 && N < 16) || (M <= 3))
                MultUM_Helper<algo2,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            else if (multeq) 
                MultUM_Helper<algo3,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            else if (N <= 3)
                MultUM_Helper<algo4,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
#ifdef _OPENMP
            else if (!omp_in_parallel() && N >= 64 && Mc*Mc*Nc > TMV_Q6)
                MultUM_Helper<36,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
#endif
            else
                MultUM_Helper<algo5,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo 33: Bad majority, copy m2 and turn into a MultEq operation.
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<33,cs,rs,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_UM
            const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
            std::cout<<"UM algo 33: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            NoAliasCopy(m2,m3);
            typedef typename M3::const_view_type M3cv;
            MultUM_Helper<-2,cs,rs,add,ix,T,M1,M3cv,M3>::call(
                x,m1,m3.view(),m3);
        }
    };

#ifdef _OPENMP
    // algo 36: Split problem into smaller parts with OpenMP for 
    // parallelization.
    // We take a pretty simple approach here, and just split up m2 and m3
    // matrices by columns and let each thread do a single matrix.
    // Then each thread calls algo 17 or 27 to calculate its product.
    // Also, we require that all but the last thread has a column width
    // that is a multiple of 16.  This way we get the maximum advantage from
    // our blocking structure while keeping the threads as balanced as 
    // possible.
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<36,cs,rs,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
#ifdef PRINTALGO_UM
            const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
            std::cout<<"UM algo 36: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            const bool upper1 = M1::mupper;
            const int algo2 = upper1 ? 17 : 27;
            bool bad_alloc = false;
#ifdef PRINTALGO_UM_OMP
            std::ofstream fout("omp.out");
#endif
#pragma omp parallel
            {
                try {
                    int num_threads = omp_get_num_threads();
                    int mythread = omp_get_thread_num();
#ifdef PRINTALGO_UM_OMP
#pragma omp critical
                    {
                        fout<<"thread "<<mythread<<"/"<<num_threads<<std::endl;
                    }
#endif
                    if (num_threads == 1) {
#ifdef PRINTALGO_UM_OMP
#pragma omp critical
                        {
                            fout<<"thread "<<mythread<<"/"<<num_threads;
                            fout<<"\nonly 1 thread"<<std::endl;
                        }
#endif
                        MultUM_Helper<algo2,cs,rs,add,ix,T,M1,M2,M3>::call(
                            x,m1,m2,m3);
                    } else {
                        int Nx = N / num_threads;
                        Nx = ((((Nx-1)>>4)+1)<<4); 
                        int j1 = mythread * Nx;
                        int j2 = (mythread+1) * Nx;
                        if (j2 > N || mythread == num_threads-1) j2 = N;
#ifdef PRINTALGO_UM_OMP
#pragma omp critical
                        {
                            fout<<"thread "<<mythread<<"/"<<num_threads;
                            fout<<"Nx = "<<Nx<<std::endl;
                            fout<<"j1 = "<<j1<<std::endl;
                            fout<<"j2 = "<<j2<<std::endl;
                        }
#endif
                        if (j1 < N)  {
                            typedef typename M2::const_colrange_type M2c;
                            typedef typename M3::colrange_type M3c;
                            const int rsx = UNKNOWN; 
                            M2c m2c = m2.cColRange(j1,j2);
                            M3c m3c = m3.cColRange(j1,j2);
#ifdef PRINTALGO_UM_OMP
#pragma omp critical
                            {
                                fout<<"thread "<<mythread<<"/"<<num_threads;
                                fout<<"\nm1 = "<<m1<<std::endl;
                                fout<<"m2c = "<<m2c<<std::endl;
                                fout<<"m3c = "<<m3c<<std::endl;
                            }
#endif
                            MultUM_Helper<
                                algo2,cs,rsx,add,ix,T,M1,M2c,M3c>::call(
                                    x,m1,m2c,m3c);
#ifdef PRINTALGO_UM_OMP
#pragma omp critical
                            {
                                fout<<"thread "<<mythread<<"/"<<num_threads;
                                fout<<"\nm3c => "<<m3c<<std::endl;
                            }
#endif
                        }
                    }
                } catch (...) {
                    // should only be std::bad_alloc, but it's good form to 
                    // catch everything inside a parallel block
                    bad_alloc = true;
                }
            }
            if (bad_alloc)
                MultUM_Helper<algo2,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };
#endif

    // algo 38: Unknown cs, check if M is small
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<38,cs,rs,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef TMV_OPT_SMALL
            TMVStaticAssert(cs == UNKNOWN);
            const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
#ifdef PRINTALGO_UM
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
            std::cout<<"UM algo 38: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            const bool upper1 = M1::mupper;
            const int algo2 = upper1 ? 16 : 26;

            if (M <= 5) {
                // then it is worth figuring out what M is.
                switch (M) {
                  case 0 :
                       // do nothing
                       break;
                  case 1 :
                       MultUM_Helper<1,1,rs,add,ix,T,M1,M2,M3>::call(
                           x,m1,m2,m3);
                       break;
                  case 2 :
                       MultUM_Helper<algo2,2,rs,add,ix,T,M1,M2,M3>::call(
                           x,m1,m2,m3);
                       break;
                  case 3 :
                       MultUM_Helper<algo2,3,rs,add,ix,T,M1,M2,M3>::call(
                           x,m1,m2,m3);
                       break;
                  case 4 :
                       MultUM_Helper<algo2,4,rs,add,ix,T,M1,M2,M3>::call(
                           x,m1,m2,m3);
                       break;
                  case 5 :
                       MultUM_Helper<algo2,5,rs,add,ix,T,M1,M2,M3>::call(
                           x,m1,m2,m3);
                }
            } else 
#endif
                MultUM_Helper<31,cs,rs,add,ix,T,M1,M2,M3>::call( x,m1,m2,m3); 
        }
    };

    // algo 81: copy m1
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<81,cs,rs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = cs == UNKNOWN ? int(m3.colsize()) : cs;
#ifdef PRINTALGO_MV_MM
            const int N = rs == UNKNOWN ? int(m3.rowsize()) : rs;
            std::cout<<"UM algo 81: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::value_type T1;
            const bool rm = M2::mcolmajor || M3::mrowmajor; 
            const int s1 = M1::mshape;
            typedef typename MCopyHelper<T1,s1,cs,cs,rm,false>::type M1c;
            M1c m1c(M);
            typedef typename M1c::view_type M1cv;
            typedef typename M1c::const_view_type M1ccv;
            M1cv m1cv = m1c.view();
            M1ccv m1ccv = m1c.view();
            CopyU_Helper<-2,cs,M1,M1cv>::call(m1,m1cv);
            MultUM_Helper<-2,cs,rs,add,ix,T,M1ccv,M2,M3>::call(
                x,m1ccv,m2,m3);
        }
    };

    // algo 82: copy x*m1
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<82,cs,rs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = cs == UNKNOWN ? int(m3.colsize()) : cs;
#ifdef PRINTALGO_UM
            const int N = rs == UNKNOWN ? int(m3.rowsize()) : rs;
            std::cout<<"UM algo 82: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::value_type T1;
            typedef typename Traits2<T,T1>::type PT1;
            const bool rm = M2::mcolmajor || M3::mrowmajor; 
            const int s1 = ShapeTraits<M1::mshape>::nonunit_shape;
            typedef typename MCopyHelper<PT1,s1,cs,cs,rm,false>::type M1c;
            M1c m1c(M);
            typedef typename M1c::view_type M1cv;
            typedef typename M1c::const_view_type M1ccv;
            M1cv m1cv = m1c.view();
            M1ccv m1ccv = m1c.view();
            typedef typename Traits<T>::real_type RT;
            const Scaling<1,RT> one;
            MultXU_Helper<-2,cs,false,ix,T,M1,M1cv>::call(x,m1,m1cv);
            MultUM_Helper<-2,cs,rs,add,1,RT,M1ccv,M2,M3>::call(
                one,m1ccv,m2,m3);
        }
    };

    // algo 83: Copy m1, figure out where to put x
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<83,cs,rs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = cs == UNKNOWN ? int(m3.colsize()) : cs;
            const int N = rs == UNKNOWN ? int(m3.rowsize()) : rs;
            if (M > N) {
                MultUM_Helper<81,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else {
                MultUM_Helper<82,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            }
        }
    };
    // If ix == 1, don't need the branch - just go to 81
    template <int cs, int rs, 
              bool add, class T, class M1, class M2, class M3>
    struct MultUM_Helper<83,cs,rs,add,1,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<1,T>& x, const M1& m1, const M2& m2, M3& m3)
        { MultUM_Helper<81,cs,rs,add,1,T,M1,M2,M3>::call(x,m1,m2,m3); }
    };

    // algo 84: copy m2
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<84,cs,rs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = cs == UNKNOWN ? int(m3.colsize()) : cs;
            const int N = rs == UNKNOWN ? int(m3.rowsize()) : rs;
#ifdef PRINTALGO_UM
            std::cout<<"UM algo 84: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename M2::value_type T2;
            const bool rm = M1::mcolmajor && M3::mrowmajor;
            typedef typename MCopyHelper<T2,Rec,cs,rs,rm,false>::type M2c;
            M2c m2c(M,N);
            typedef typename M2c::view_type M2cv;
            typedef typename M2c::const_view_type M2ccv;
            M2cv m2cv = m2c.view();
            M2ccv m2ccv = m2c.view();
            CopyM_Helper<-2,cs,rs,M2,M2cv>::call(m2,m2cv);
            MultUM_Helper<-2,cs,rs,add,ix,T,M1,M2ccv,M3>::call(
                x,m1,m2ccv,m3);
        }
    };

    // algo 85: copy x*m2
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<85,cs,rs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = cs == UNKNOWN ? int(m3.colsize()) : cs;
            const int N = rs == UNKNOWN ? int(m3.rowsize()) : rs;
#ifdef PRINTALGO_UM
            std::cout<<"UM algo 85: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename M2::value_type T2;
            typedef typename Traits2<T,T2>::type PT2;
            const bool rm = M1::mcolmajor && M3::mrowmajor;
            typedef typename MCopyHelper<PT2,Rec,cs,rs,rm,false>::type M2c;
            M2c m2c(M,N);
            typedef typename M2c::view_type M2cv;
            typedef typename M2c::const_view_type M2ccv;
            M2cv m2cv = m2c.view();
            M2ccv m2ccv = m2c.view();
            typedef typename Traits<T>::real_type RT;
            const Scaling<1,RT> one;
            MultXM_Helper<-2,cs,rs,false,ix,T,M2,M2cv>::call(x,m2,m2cv);
            MultUM_Helper<-2,cs,rs,add,1,RT,M1,M2ccv,M3>::call(
                one,m1,m2ccv,m3);
        }
    };

    // algo 86: Copy m2, figure out where to put x
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<86,cs,rs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = cs == UNKNOWN ? int(m3.colsize()) : cs;
            const int N = rs == UNKNOWN ? int(m3.rowsize()) : rs;
            if (N > M) {
                MultUM_Helper<84,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else {
                MultUM_Helper<85,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            }
        }
    };
    // If ix == 1, don't need the branch - just go to 84
    template <int cs, int rs, bool add, class T, class M1, class M2, class M3>
    struct MultUM_Helper<86,cs,rs,add,1,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<1,T>& x, const M1& m1, const M2& m2, M3& m3)
        { MultUM_Helper<84,cs,rs,add,1,T,M1,M2,M3>::call(x,m1,m2,m3); }
    };

    // algo 87: Use temporary for m1*m2
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<87,cs,rs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = cs == UNKNOWN ? int(m3.colsize()) : cs;
            const int N = rs == UNKNOWN ? int(m3.rowsize()) : rs;
#ifdef PRINTALGO_UM
            std::cout<<"UM algo 87: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename Traits<T>::real_type RT;
            const Scaling<1,RT> one;
            if (N > M) {
                typedef typename M1::value_type T1;
                typedef typename M2::value_type T2;
                typedef typename Traits2<T1,T2>::type PT3;
                const bool rm = M1::mrowmajor && M2::mrowmajor;
                typedef typename MCopyHelper<PT3,Rec,cs,rs,rm,false>::type M3c;
                M3c m3c(M,N);
                typedef typename M3c::view_type M3cv;
                typedef typename M3c::const_view_type M3ccv;
                M3cv m3cv = m3c.view();
                M3ccv m3ccv = m3c.view();
                MultUM_Helper<-2,cs,rs,false,1,RT,M1,M2,M3cv>::call(
                    one,m1,m2,m3cv);
                MultXM_Helper<-2,cs,rs,add,ix,T,M3ccv,M3>::call(x,m3ccv,m3);
            } else {
                typedef typename M3::value_type T3;
                const bool rm = M1::mrowmajor && M2::mrowmajor;
                typedef typename MCopyHelper<T3,Rec,cs,rs,rm,false>::type M3c;
                M3c m3c(M,N);
                typedef typename M3c::view_type M3cv;
                typedef typename M3c::const_view_type M3ccv;
                M3cv m3cv = m3c.view();
                M3ccv m3ccv = m3c.view();
                MultUM_Helper<-2,cs,rs,false,ix,T,M1,M2,M3cv>::call(
                    x,m1,m2,m3cv);
                MultXM_Helper<-2,cs,rs,add,1,RT,M3ccv,M3>::call(one,m3ccv,m3);
            }
        }
    };
    // If ix == 1, don't need the branch
    template <int cs, int rs, bool add, 
              class T, class M1, class M2, class M3>
    struct MultUM_Helper<87,cs,rs,add,1,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<1,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = cs == UNKNOWN ? int(m3.colsize()) : cs;
            const int N = rs == UNKNOWN ? int(m3.rowsize()) : rs;
#ifdef PRINTALGO_UM
            std::cout<<"UM algo 87: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename Traits<T>::real_type RT;
            const Scaling<1,RT> one;
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename Traits2<T1,T2>::type PT3;
            const bool rm = M1::mrowmajor && M2::mrowmajor;
            typedef typename MCopyHelper<PT3,Rec,cs,rs,rm,false>::type M3c;
            M3c m3c(M,N);
            typedef typename M3c::view_type M3cv;
            typedef typename M3c::const_view_type M3ccv;
            M3cv m3cv = m3c.view();
            M3ccv m3ccv = m3c.view();
            MultUM_Helper<-2,cs,rs,false,1,RT,M1,M2,M3cv>::call(
                one,m1,m2,m3cv);
            MultXM_Helper<-2,cs,rs,add,1,T,M3ccv,M3>::call(x,m3ccv,m3);
        }
    };

    // algo -4: No branches or copies
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<-4,cs,rs,add,ix,T,M1,M2,M3> 
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const bool upper1 = M1::mupper;
            const int algo = 
                ( cs == 0 || rs == 0 ) ? 0 :
                cs == 1 ? 401 :
                rs == 1 ? 402 :
                ( cs != UNKNOWN ) ? (
                    upper1 ? (cs <= 5 ? 16 : 17 ) :
                    (cs <= 5 ? 26 : 27 ) ) :
                upper1 ? 17 : 27;
            MultUM_Helper<algo,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<-3,cs,rs,add,ix,T,M1,M2,M3> 
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            // Possible algorithms are:
            //
            // Trivial and special for small TriMatrix sizes.
            //  0 = cs or rs == 0, so nothing to do
            //  1 = cs == 1: reduces to trivial MultXV function
            //  2 = rs == 1: reduces to trivial MultMV function
            //
            // UpperTri:
            // 11 = loop over n: MultUV
            // 12 = loop over m: MultMV
            // 13 = loop over k: Rank1 
            // 16 = specialization for cs = 2,3,4,5
            // 17 = split trimatrix into 3 submatrices
            //
            // LowerTri:
            // 21 = loop over n: MultUV
            // 22 = loop over m: MultMV
            // 23 = loop over k: Rank1
            // 26 = specialization for cs = 2,3,4,5
            // 27 = split trimatrix into 3 submatrices
            // 
            // Overall drivers:
            // 31 = Choose algorithm based on the (runtime) size
            // 33 = Bad majority, copy m2 and turn into a MultEq operation.
            // 36 = Parallelize using openmp 
            // 38 = Check if M is small

            const bool xcc = M2::mcolmajor && M3::mcolmajor;
            const bool rxr = M1::mrowmajor && M3::mrowmajor;
            const bool crx = M1::mcolmajor && M2::mrowmajor;
            const bool upper1 = M1::mupper;

#if 0
            const int algo = 33;
#else
#if TMV_OPT == 0 
            const int algo =
                upper1 ? 
                ( rxr ? 12 : crx ? 13 : xcc ? 11 : 13 ) :
                ( rxr ? 22 : crx ? 23 : xcc ? 21 : 23 );
#else
            const bool nxx = !(M1::mcolmajor || M1::mrowmajor);
            const bool xxn = !(M3::mcolmajor || M3::mrowmajor);
            // xcc, rxr, crx have fast algorithms already
            // add can't be turned into a MultEq op
            // nxx, xxn don't benefit from the change
            //   (and would go into an infinite loop)
            const bool multeq = !(xcc || rxr || crx || add || nxx || xxn);
#ifdef _OPENMP
            const int Mc = cs == UNKNOWN ? UNKNOWN : (cs < 16) ? 1 : (cs>>4);
            const int Nc = rs == UNKNOWN ? UNKNOWN : (rs < 16) ? 1 : (rs>>4);
            const int McMcNc = IntTraits2<IntTraits2<Mc,Mc>::prod,Nc>::prod;
#endif
            const int algo = 
                ( cs == 0 || rs == 0 ) ? 0 :
                cs == 1 ? 201 :
                rs == 1 ? 202 :
                upper1 ? (
                    cs == UNKNOWN ? 38 : 
                    cs <= 5 ? 16 :
                    rs == UNKNOWN ? 31 :
                    multeq ? 33 :
#ifdef _OPENMP
                    (rs >= 64 && McMcNc >= TMV_Q6) ? 36 :
#endif
                    17 ) :
                ( // lowertri
                    cs == UNKNOWN ? 38 : 
                    cs <= 5 ? 26 :
                    rs == UNKNOWN ? 31 :
                    multeq ? 33 :
#ifdef _OPENMP
                    (rs >= 64 && McMcNc >= TMV_Q6) ? 36 :
#endif
                    27 );
#endif
#endif
#ifdef PRINTALGO_UM
            const int M = cs==UNKNOWN ? int(m3.colsize()) : cs;
            const int N = rs==UNKNOWN ? int(m3.rowsize()) : rs;
            std::cout<<"InlineMultUM: x = "<<ix<<"  "<<T(x)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"M = "<<M<<"  N = "<<N<<std::endl;
            std::cout<<"cs = "<<cs<<"  rs = "<<rs<<std::endl;
            std::cout<<"add = "<<add<<", algo = "<<algo<<std::endl;
#endif
            MultUM_Helper<algo,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo 97: Conjugate
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<97,cs,rs,add,ix,T,M1,M2,M3>
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
            MultUM_Helper<-2,cs,rs,add,ix,T,M1c,M2c,M3c>::call(
                TMV_CONJ(x),m1c,m2c,m3c);
        }
    };

    // algo 98: call InstMultMM
    template <int cs, int rs, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<98,cs,rs,false,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typename M3::value_type xx(x);
            InstMultMM(xx,m1.xdView(),m2.xView(),m3.xView());
        }
    };
    template <int cs, int rs, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<98,cs,rs,true,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typename M3::value_type xx(x);
            InstAddMultMM(xx,m1.xdView(),m2.xView(),m3.xView());
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<-2,cs,rs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            const bool inst = 
                M1::msize == UNKNOWN &&
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
                ( cs == 0 || rs == 0 ) ? 0 :
                cs == 1 ? 201 :
                rs == 1 ? 202 :
                M3::mconj ? 97 :
                inst ? 98 : 
                -3;
            MultUM_Helper<algo,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo 99: Check for aliases
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<99,cs,rs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const bool s1 = SameStorage(m1,m3);
            const bool s2 = SameStorage(m2,m3) && !ExactSameStorage(m2,m3);
            if (!s1 && !s2) {
                // No aliasing (or no clobber)
                MultUM_Helper<-2,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else if (!add && !s1) {
                // MultEq: copy m2
                AliasCopy(m2,m3);
                NoAliasMultMM<add>(x,m1,m3,m3);
            } else if (s1 && !s2) {
                // copy m1
                MultUM_Helper<83,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else { 
                // Use temporary for m1*m2
                MultUM_Helper<87,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            }
        }
    };

    // algo -1: Check for aliases?
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultUM_Helper<-1,cs,rs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const bool checkalias =
                M1::msize == UNKNOWN &&
                M2::mcolsize == UNKNOWN && M2::mrowsize == UNKNOWN &&
                M3::mcolsize == UNKNOWN && M3::mrowsize == UNKNOWN;
            const int algo = 
                ( cs == 0 || rs == 0 ) ? 0 :
                cs == 1 ? 1 :
                rs == 1 ? 2 :
                checkalias ? 99 : 
                -2;
            MultUM_Helper<algo,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    template <int algo, bool add, int ix, class T, class M1, class M2, class M3>
    inline void DoMultUM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Rec<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::msize,M3::mcolsize>::same));
        TMVStaticAssert((Sizes<M1::msize,M2::mcolsize>::same));
        TMVStaticAssert((Sizes<M2::mrowsize,M3::mrowsize>::same));
        TMVAssert(m1.size() == m3.colsize());
        TMVAssert(m1.size() == m2.colsize());
        TMVAssert(m2.rowsize() == m3.rowsize());

        const int cs = Sizes<M3::mcolsize,M1::msize>::size;
        const int rs = Sizes<M3::mrowsize,M2::mrowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        M1v m1v = m1.cView();
        M2v m2v = m2.cView();
        M3v m3v = m3.cView();
        MultUM_Helper<algo,cs,rs,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Rec<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    { DoMultUM<-1,add>(x,m1,m2,m3); }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void NoAliasMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Rec<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    { DoMultUM<-2,add>(x,m1,m2,m3); }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Rec<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    { DoMultUM<-3,add>(x,m1,m2,m3); }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void AliasMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Tri<M1>& m1, const BaseMatrix_Rec<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    { DoMultUM<99,add>(x,m1,m2,m3); }

    template <int algo, bool add, int ix, class T, class M1, class M2, class M3>
    inline void DoMultMU(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::mcolsize,M3::mcolsize>::same));
        TMVStaticAssert((Sizes<M1::mrowsize,M2::msize>::same));
        TMVStaticAssert((Sizes<M2::msize,M3::mrowsize>::same));
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.size());
        TMVAssert(m2.size() == m3.rowsize());

        const int cs = Sizes<M3::mcolsize,M1::mcolsize>::size;
        const int rs = Sizes<M3::mrowsize,M2::msize>::size;
        typedef typename M1::const_transpose_type M1t;
        typedef typename M2::const_transpose_type M2t;
        typedef typename M3::transpose_type M3t;
        M1t m1t = m1.transpose();
        M2t m2t = m2.transpose();
        M3t m3t = m3.transpose();
        MultUM_Helper<algo,rs,cs,add,ix,T,M2t,M1t,M3t>::call(x,m2t,m1t,m3t);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    { DoMultMU<-1,add>(x,m1,m2,m3); }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void NoAliasMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    { DoMultMU<-2,add>(x,m1,m2,m3); }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    { DoMultMU<-3,add>(x,m1,m2,m3); }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void AliasMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    { DoMultMU<99,add>(x,m1,m2,m3); }

    template <class M1, int ix, class T, class M2>
    inline void MultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2)
    { MultMM<false>(x,m1.mat(),m2.mat(),m1.mat()); }

    template <class M1, int ix, class T, class M2>
    inline void NoAliasMultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2)
    { NoAliasMultMM<false>(x,m1.mat(),m2.mat(),m1.mat()); }

    template <class M1, int ix, class T, class M2>
    inline void AliasMultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M2>& m2)
    { AliasMultMM<false>(x,m1.mat(),m2.mat(),m1.mat()); }

#undef TMV_OPT_SMALL
#undef TMV_OPT_CLEANUP
#undef TMV_Q2
#undef TMV_Q6

} // namespace tmv

#endif 
