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


#ifndef TMV_DivMU_H
#define TMV_DivMU_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_ScaleV.h"
#include "TMV_MultMV.h"
#include "TMV_Rank1VVM.h"
#include "TMV_DivVU.h"
#include "TMV_InvertU.h"
#include "TMV_MultXM_Funcs.h"
#include "TMV_DivMM_Funcs.h"

#ifdef _OPENMP
#include "omp.h"
#ifdef PRINTALGO_DivU_OMP
#include <fstream>
#endif
#endif

#ifdef PRINTALGO_DivU
#include <iostream>
#endif

// Use the specialized 1,2,3,4 sized algorithms for the end of the 
// recursive algorithm.
#if TMV_OPT >= 2
#define TMV_DIVMU_CLEANUP
#endif

// Check for small (<=5) values of cs or rs
// This leads to significant speed improvements for such matrices at little
// cost to larger matrices, but the large increase in code size and the fact
// that it only benefits a few particular sizes of matrices mean that
// we require TMV_OPT = 3 for it.
#if TMV_OPT >= 3
#define TMV_DIVMU_SMALL
#endif

// The maximum size to stop recursing
#if TMV_OPT >= 3
#define TMV_DIVMU_RECURSE 4
#else
#define TMV_DIVMU_RECURSE 1
#endif

// Inline the MV (MultUV, Rank1VVM) calls.
#if TMV_OPT >= 1
#define TMV_DIVMU_INLINE_MV
#endif

// Inline the small-sized MM (MultMM) calls.
// (Only relevant if RECURSE <= 16.)
#if TMV_OPT >= 2
#define TMV_DIVMU_INLINE_MM
#endif

// The minimum value of (M^2*N / 16^3) to use multiple threads.
// (We also require that N > 64 so we have something to split.)
#define TMV_DIVMU_OMP_THRESH 64


namespace tmv {

    // Defined in TMV_DivMU.cpp
    template <class T1, class T2, int C2>
    void InstTriLDivEq(
        MatrixView<T1> m1, const ConstUpperTriMatrixView<T2,C2>& m2);
    template <class T1, class T2, int C2>
    void InstTriLDivEq(
        MatrixView<T1> m1, const ConstLowerTriMatrixView<T2,C2>& m2);

    template <class T1, class T2, int C2>
    void InstAliasTriLDivEq(
        MatrixView<T1> m1, const ConstUpperTriMatrixView<T2,C2>& m2);
    template <class T1, class T2, int C2>
    void InstAliasTriLDivEq(
        MatrixView<T1> m1, const ConstLowerTriMatrixView<T2,C2>& m2);


    //
    // m1 /= m2
    // m2 is UpperTri or LowerTri
    //

    template <int algo, int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper;

    // algo 0: Trivial, nothing to do.
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<0,cs,rs,M1,M2>
    { static TMV_INLINE void call(M1& , const M2& ) {} };

    // algo 1: M == 1, so reduces to ScaleV
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<1,cs,rs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DivU
            const int N = rs==UNKNOWN ? int(m1.rowsize()) : rs;
            std::cout<<"LDivEqMU algo 1: M,N,cs,rs = "<<1<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            const bool u2 = M2::_unit;
            typedef typename M1::row_type M1r;
            typedef typename M2::real_type RT;
            typedef typename M2::value_type T2;
            typedef typename TypeSelect<u2,RT,T2>::type XT2;
            const int ix2 = u2 ? 1 : 0;

            M1r m1r = m1.get_row(0);
            const Scaling<ix2,XT2> inv00(
                Maybe<!u2>::invprod( m2.cref(0,0) , RT(1) ));
            ScaleV_Helper<-3,rs,ix2,XT2,M1r>::call(inv00, m1r);
        }
    };

    // algo 101: same as 1, but use -1 algo
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<101,cs,rs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DivU
            const int N = rs==UNKNOWN ? int(m1.rowsize()) : rs;
            std::cout<<"LDivEqMU algo 101: M,N,cs,rs = "<<1<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            const bool u2 = M2::_unit;
            typedef typename M1::row_type M1r;
            typedef typename M2::real_type RT;
            typedef typename M2::value_type T2;
            typedef typename TypeSelect<u2,RT,T2>::type XT2;
            const int ix2 = u2 ? 1 : 0;

            M1r m1r = m1.get_row(0);
            const Scaling<ix2,XT2> inv00(
                Maybe<!u2>::invprod( m2.cref(0,0) , RT(1) ));
            ScaleV_Helper<-1,rs,ix2,XT2,M1r>::call(inv00, m1r);
        }
    };

    // algo 201: same as 1, but use -2 algo
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<201,cs,rs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DivU
            const int N = rs==UNKNOWN ? int(m1.rowsize()) : rs;
            std::cout<<"LDivEqMU algo 201: M,N,cs,rs = "<<1<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            const bool u2 = M2::_unit;
            typedef typename M1::row_type M1r;
            typedef typename M2::real_type RT;
            typedef typename M2::value_type T2;
            typedef typename TypeSelect<u2,RT,T2>::type XT2;
            const int ix2 = u2 ? 1 : 0;

            M1r m1r = m1.get_row(0);
            const Scaling<ix2,XT2> inv00(
                Maybe<!u2>::invprod( m2.cref(0,0) , RT(1) ));
            ScaleV_Helper<-2,rs,ix2,XT2,M1r>::call(inv00, m1r);
        }
    };

    // algo 2: N == 1, so reduces to LDivVU
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<2,cs,rs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DivU
            const int M = cs==UNKNOWN ? int(m1.colsize()) : cs;
            std::cout<<"LDivEqMU algo 2: M,N,cs,rs = "<<M<<','<<1<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::col_type M1c;

            M1c m1c = m1.get_col(0);
            LDivEqVU_Helper<-4,cs,M1c,M2>::call(m1c,m2);
        }
    };

    // algo 102: same as 2, but use -1 algo
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<102,cs,rs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DivU
            const int M = cs==UNKNOWN ? int(m1.colsize()) : cs;
            std::cout<<"LDivEqMU algo 102: M,N,cs,rs = "<<M<<','<<1<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::col_type M1c;

            M1c m1c = m1.get_col(0);
            LDivEqVU_Helper<-1,cs,M1c,M2>::call(m1c,m2);
        }
    };

    // algo 202: same as 2, but use -2 algo
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<202,cs,rs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DivU
            const int M = cs==UNKNOWN ? int(m1.colsize()) : cs;
            std::cout<<"LDivEqMU algo 202: M,N,cs,rs = "<<M<<','<<1<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::col_type M1c;

            M1c m1c = m1.get_col(0);
            LDivEqVU_Helper<-2,cs,M1c,M2>::call(m1c,m2);
        }
    };

    // algo 11: UpperTri loop over n
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<11,cs,rs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            const int N = rs==UNKNOWN ? int(m1.rowsize()) : rs;
#ifdef PRINTALGO_DivU
            const int M = cs==UNKNOWN ? int(m1.colsize()) : cs;
            std::cout<<"LDivEqMU algo 11: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif

            typedef typename M1::col_type M1c;
            for(int j=0;j<N;++j) {
                // m1.col(j) /= m2
                M1c m1j = m1.get_col(j);
                LDivEqVU_Helper<-4,cs,M1c,M2>::call(m1j,m2);
            }
        }
    };

    // algo 12: UpperTri loop over m
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<12,cs,rs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            const int M = cs==UNKNOWN ? int(m1.colsize()) : cs;
#ifdef PRINTALGO_DivU
            const int N = rs==UNKNOWN ? int(m1.rowsize()) : rs;
            std::cout<<"LDivEqMU algo 12: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            const bool u2 = M2::_unit;
            typedef typename M2::real_type RT;
            typedef typename M2::value_type T2;
            typedef typename TypeSelect<u2,RT,T2>::type XT2;
            typedef typename M2::const_row_sub_type M2r;
            typedef typename M1::row_type M1r;
            typedef typename M1::rowrange_type::const_transpose_type M1rrt;
            const int ix2 = u2 ? 1 : 0;
            const int xx = UNKNOWN;
            const Scaling<-1,RT> mone;
#ifdef TMV_DIVMU_INLINE_MV
            const int algo2 = -4;
#else
            const int algo2 = -2;
#endif

            for(int i=M;i--;) {
                // m1.row(i) -= m2.row(i,i+1,M) * m1.rowRange(i+1,M)
                // m1.row(i) /= m2(i,i)
                M1r m1i = m1.get_row(i);
                M1rrt m1rrt = m1.cRowRange(i+1,M).transpose();
                M2r m2i = m2.get_row(i,i+1,M);
                MultMV_Helper<algo2,rs,xx,true,-1,RT,M1rrt,M2r,M1r>::call(
                    mone,m1rrt,m2i,m1i);
                const Scaling<ix2,XT2> invii(
                    Maybe<!u2>::invprod( m2.cref(i,i) , RT(1) ));
                ScaleV_Helper<-3,rs,ix2,XT2,M1r>::call(invii,m1i);
            }
        }
    };

    // algo 13: UpperTri loop over k
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<13,cs,rs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            const int M = cs==UNKNOWN ? int(m1.colsize()) : cs;
#ifdef PRINTALGO_DivU
            const int N = rs==UNKNOWN ? int(m1.rowsize()) : rs;
            std::cout<<"LDivEqMU algo 13: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            const bool u2 = M2::_unit;
            typedef typename M2::real_type RT;
            typedef typename M2::value_type T2;
            typedef typename TypeSelect<u2,RT,T2>::type XT2;
            typedef typename M1::row_type M1r;
            typedef typename M1::rowrange_type M1rr;
            typedef typename M1::const_row_type M1rc;
            typedef typename M2::const_col_sub_type M2c;
            const int ix2 = u2 ? 1 : 0;
            const int xx = UNKNOWN;
            const Scaling<-1,RT> mone;
#ifdef TMV_DIVMU_INLINE_MV
            const int algo2 = -4;
#else
            const int algo2 = -2;
#endif

            for(int k=M;k--;) {
                // m1.row(k) /= m2(k,k)
                // m1.rowRange(0,k) -= m2.col(k,0,k) ^ m1.row(k)
                M1r m1k = m1.get_row(k);
                M1rr m1rr = m1.cRowRange(0,k);
                M2c m2k = m2.get_col(k,0,k);
                const Scaling<ix2,XT2> invkk(
                    Maybe<!u2>::invprod( m2.cref(k,k) , RT(1) ));
                ScaleV_Helper<-3,rs,ix2,XT2,M1r>::call(invkk,m1k);
                Rank1VVM_Helper<algo2,xx,rs,true,-1,RT,M2c,M1rc,M1rr>::call(
                    mone,m2k,m1k,m1rr);
            }
        }
    };

    // algo 16: For small cs, create m2.inverse() explicitly
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<16,cs,rs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DivU
            const int M = cs==UNKNOWN ? int(m1.colsize()) : cs;
            const int N = rs==UNKNOWN ? int(m1.rowsize()) : rs;
            std::cout<<"LDivEqMU algo 16: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M2::real_type RT;
            typedef typename M2::value_type T2;
            typedef SmallUpperTriMatrix<T2,cs,M2::_dt> Uinv;
            Uinv m2inv = m2;
            InvertU_Helper<-3,cs,Uinv>::call(m2inv);
            const Scaling<1,RT> one;
            NoAliasMultMM<false>(one,m2inv,m1,m1);
        }
    };
    template <int rs, class M1, class M2>
    struct LDivEqMU_Helper<16,UNKNOWN,rs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            const int M = m1.colsize();
#ifdef PRINTALGO_DivU
            const int N = rs==UNKNOWN ? int(m1.rowsize()) : rs;
            std::cout<<"LDivEqMU algo 16: M,N,cs,rs = "<<M<<','<<N<<
                ','<<UNKNOWN<<','<<rs<<std::endl;
#endif
#if TMV_DIVMU_RECURSE > 4
            const bool rr = M1::_rowmajor && M2::_rowmajor;
            const bool rc = M1::_rowmajor && M2::_colmajor;
            const bool cx = M1::_colmajor;
            const int algo2 = rr ? 12 : rc ? 13 : cx ? 11 : 13;
#endif
            switch (M) {
              case 0 :
                   // do nothing
                   break;
              case 1 :
                   LDivEqMU_Helper<1,1,rs,M1,M2>::call(m1,m2);
                   break;
#if TMV_DIVMU_RECURSE > 1
              case 2 :
                   LDivEqMU_Helper<16,2,rs,M1,M2>::call(m1,m2);
                   break;
#endif
#if TMV_DIVMU_RECURSE > 2
              case 3 :
                   LDivEqMU_Helper<16,3,rs,M1,M2>::call(m1,m2);
                   break;
#endif
#if TMV_DIVMU_RECURSE > 3
              case 4 :
                   LDivEqMU_Helper<16,4,rs,M1,M2>::call(m1,m2);
                   break;
#endif
#if TMV_DIVMU_RECURSE > 4
              default :
                   LDivEqMU_Helper<algo2,UNKNOWN,rs,M1,M2>::call(m1,m2);
#endif
            }
        }
    };

    // algo 17: Split the UpperTriMatrix into 3 sections and recurse
    // the calculation on each of them:
    // ( A B ) ( D ) = ( F )
    // ( 0 C ) ( E )   ( G )
    // F = AD + BE
    // G = CE
    // Solving for D,E yields:
    // E /= C
    // D -= BE
    // D /= A
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<17,cs,rs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            const int M = cs==UNKNOWN ? int(m1.colsize()) : cs;
#ifdef PRINTALGO_DivU
            const int N = rs==UNKNOWN ? int(m1.rowsize()) : rs;
            std::cout<<"LDivEqMU algo 17: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif

            const bool rr = M1::_rowmajor && M2::_rowmajor;
            const bool rc = M1::_rowmajor && M2::_colmajor;
            const bool cx = M1::_colmajor;
            const int algo2 =  // The algorithm for M <= RECURSE
                cs == 0 ? 0 :
                cs == 1 ? 1 :
                (cs != UNKNOWN && cs > TMV_DIVMU_RECURSE) ? 0 :
                TMV_DIVMU_RECURSE == 1 ? 1 :
#ifdef TMV_DIVMU_CLEANUP
                cs == UNKNOWN ? 16 :
#endif
                (cs != UNKNOWN && cs <= 3) ? 16 :
                rr ? 12 : rc ? 13 : cx ? 11 : 13;
            const int algo3 =  // The algorithm for M > 16
                cs == UNKNOWN ? 17 : 
#ifdef TMV_DIVMU_INLINE_MM
                cs <= 16 ? 0 : 
#endif
                cs > TMV_DIVMU_RECURSE ? 17 :
                0;
            const int algo4 =  // The algorithm for MultMM
                cs == UNKNOWN ? -2 : 
#ifdef TMV_DIVMU_INLINE_MM
                cs <= 16 ? 0 : 
#endif
                cs > TMV_DIVMU_RECURSE ? -3 :
                0;
#ifdef TMV_DIVMU_INLINE_MM
            const int algo3b =  // The algorithm for M > RECURSE
                cs == UNKNOWN || 
                (cs > TMV_DIVMU_RECURSE && cs <= 16) ? 17 : 0;
            const int algo4b =  // The algorithm for MultMM
                cs == UNKNOWN || 
                (cs > TMV_DIVMU_RECURSE && cs <= 16) ? -4 : 0;
#endif
#ifdef PRINTALGO_DivU
            std::cout<<"algo2,3,4 = "<<algo2<<"  "<<algo3<<"  "<<algo4<<std::endl;
#ifdef TMV_DIVMU_INLINE_MM
            std::cout<<"algo3b,4b = "<<algo3b<<"  "<<algo4b<<std::endl;
#endif
#endif

            typedef typename M2::real_type RT;
            typedef typename M2::const_subtrimatrix_type M2a;
            typedef typename M2::const_submatrix_type M2b;
            typedef typename M1::const_rowrange_type M1rc;
            typedef typename M1::rowrange_type M1r;
            const Scaling<-1,RT> mone;

            if (M > TMV_DIVMU_RECURSE) {
                const int Mx = M > 16 ? ((((M-1)>>5)+1)<<4) : (M>>1);
                // (If M > 16, round M/2 up to a multiple of 16.)
                const int csx = IntTraits<cs>::half_roundup;
                const int csy = IntTraits2<cs,csx>::diff;

                M2a A = m2.cSubTriMatrix(0,Mx);
                M2b B = m2.cSubMatrix(0,Mx,Mx,M);
                M2a C = m2.cSubTriMatrix(Mx,M);
                M1r D = m1.cRowRange(0,Mx);
                M1r E = m1.cRowRange(Mx,M);

#ifdef TMV_DIVMU_INLINE_MM
                if (M > 16) {
                    // For large M, make sure to use good MultMM algo
#endif
                    LDivEqMU_Helper<algo3,csy,rs,M1r,M2a>::call(E,C);
                    MultMM_Helper<algo4,csx,rs,csy,true,-1,RT,M2b,M1rc,M1r>::
                        call(mone,B,E,D);
                    LDivEqMU_Helper<algo3,csx,rs,M1r,M2a>::call(D,A);
#ifdef TMV_DIVMU_INLINE_MM
                } else {
                    // For smaller M, do the no branching algorithm
                    LDivEqMU_Helper<algo3b,csy,rs,M1r,M2a>::call(E,C);
                    MultMM_Helper<algo4b,csx,rs,csy,true,-1,RT,M2b,M1rc,M1r>::
                        call(mone,B,E,D);
                    LDivEqMU_Helper<algo3b,csx,rs,M1r,M2a>::call(D,A);
                }
#endif
            } else {
                LDivEqMU_Helper<algo2,cs,rs,M1,M2>::call(m1,m2);
            }
        }
    };

    // algo 21: LowerTri loop over n
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<21,cs,rs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            const int N = rs==UNKNOWN ? int(m1.rowsize()) : rs;
#ifdef PRINTALGO_DivU
            const int M = cs==UNKNOWN ? int(m1.colsize()) : cs;
            std::cout<<"LDivEqMU algo 21: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
#ifdef TMV_DIVMU_INLINE_MV
            const int algo2 = -4;
#else
            const int algo2 = -2;
#endif

            typedef typename M1::col_type M1c;
            for(int j=0;j<N;++j) {
                // m1.col(j) /= m2
                M1c m1j = m1.get_col(j);
                LDivEqVU_Helper<algo2,cs,M1c,M2>::call(m1j,m2);
            }
        }
    };

    // algo 22: LowerTri loop over m
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<22,cs,rs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            const int M = cs==UNKNOWN ? int(m1.colsize()) : cs;
#ifdef PRINTALGO_DivU
            const int N = rs==UNKNOWN ? int(m1.rowsize()) : rs;
            std::cout<<"LDivEqMU algo 22: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            const bool u2 = M2::_unit;
            typedef typename M2::real_type RT;
            typedef typename M2::value_type T2;
            typedef typename TypeSelect<u2,RT,T2>::type XT2;
            typedef typename M2::const_row_sub_type M2r;
            typedef typename M1::rowrange_type::const_transpose_type M1rrt;
            typedef typename M1::row_type M1r;
            const int ix2 = u2 ? 1 : 0;
            const int xx = UNKNOWN;
            const Scaling<-1,RT> mone;
#ifdef TMV_DIVMU_INLINE_MV
            const int algo2 = -4;
#else
            const int algo2 = -2;
#endif

            for(int i=0;i<M;++i) {
                // m1.row(i) -= m2.row(i,0,i) * m1.rowRange(0,i)
                // m1.row(i) /= m2(i,i)
                M1r m1i = m1.get_row(i);
                M2r m2i = m2.get_row(i,0,i);
                M1rrt m1rrt = m1.cRowRange(0,i).transpose();
                MultMV_Helper<algo2,rs,xx,true,-1,RT,M1rrt,M2r,M1r>::call(
                    mone,m1rrt,m2i,m1i);
                const Scaling<ix2,XT2> invii(
                    Maybe<!u2>::invprod( m2.cref(i,i) , RT(1) ));
                ScaleV_Helper<-3,rs,ix2,XT2,M1r>::call(invii,m1i);
            }
        }
    };

    // algo 23: LowerTri loop over k
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<23,cs,rs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            const int M = cs==UNKNOWN ? int(m1.colsize()) : cs;
#ifdef PRINTALGO_DivU
            const int N = rs==UNKNOWN ? int(m1.rowsize()) : rs;
            std::cout<<"LDivEqMU algo 23: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            const bool u2 = M2::_unit;
            typedef typename M2::real_type RT;
            typedef typename M2::value_type T2;
            typedef typename TypeSelect<u2,RT,T2>::type XT2;
            typedef typename M1::row_type M1r;
            typedef typename M1::const_row_type M1rc;
            typedef typename M1::rowrange_type M1rr;
            typedef typename M2::const_col_sub_type M2c;
            const int ix2 = u2 ? 1 : 0;
            const int xx = UNKNOWN;
            const Scaling<-1,RT> mone;
#ifdef TMV_DIVMU_INLINE_MV
            const int algo2 = -4;
#else
            const int algo2 = -2;
#endif

            for(int k=0;k<M;++k) {
                // m1.row(k) /= m2(k,k)
                // m1.rowRange(k+1,M) -= m2.col(k,k+1,M) ^ m1.row(k)
                M1r m1k = m1.get_row(k);
                M1rr m1rr = m1.cRowRange(k+1,M);
                M2c m2k = m2.get_col(k,k+1,M);
                const Scaling<ix2,XT2> invkk(
                    Maybe<!u2>::invprod( m2.cref(k,k) , RT(1) ));
                ScaleV_Helper<-3,rs,ix2,XT2,M1r>::call(invkk,m1k);
                Rank1VVM_Helper<algo2,xx,rs,true,-1,RT,M2c,M1rc,M1rr>::call(
                    mone,m2k,m1k,m1rr);
            }
        }
    };

    // algo 26: For small cs, create m2.inverse() explicitly
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<26,cs,rs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DivU
            const int M = cs==UNKNOWN ? int(m1.colsize()) : cs;
            const int N = rs==UNKNOWN ? int(m1.rowsize()) : rs;
            std::cout<<"LDivEqMU algo 26: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M2::real_type RT;
            typedef typename M2::value_type T2;
            typedef SmallLowerTriMatrix<T2,cs,M2::_dt> Linv;
            Linv m2inv = m2;
            InvertU_Helper<-3,cs,Linv>::call(m2inv);
            const Scaling<1,RT> one;
            NoAliasMultMM<false>(one,m2inv,m1,m1);
        }
    };
    template <int rs, class M1, class M2>
    struct LDivEqMU_Helper<26,UNKNOWN,rs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            const int M = m1.colsize();
#ifdef PRINTALGO_DivU
            const int N = rs==UNKNOWN ? int(m1.rowsize()) : rs;
            std::cout<<"LDivEqMU algo 26: M,N,cs,rs = "<<M<<','<<N<<
                ','<<UNKNOWN<<','<<rs<<std::endl;
#endif
#if TMV_DIVMU_RECURSE > 4
            const bool rr = M1::_rowmajor && M2::_rowmajor;
            const bool rc = M1::_rowmajor && M2::_colmajor;
            const bool cx = M1::_colmajor;
            const int algo2 = rr ? 22 : rc ? 23 : cx ? 21 : 23;
#endif
            switch (M) {
              case 0 :
                   // do nothing
                   break;
              case 1 :
                   LDivEqMU_Helper<1,1,rs,M1,M2>::call(m1,m2);
                   break;
#if TMV_DIVMU_RECURSE > 1
              case 2 :
                   LDivEqMU_Helper<26,2,rs,M1,M2>::call(m1,m2);
                   break;
#endif
#if TMV_DIVMU_RECURSE > 2
              case 3 :
                   LDivEqMU_Helper<26,3,rs,M1,M2>::call(m1,m2);
                   break;
#endif
#if TMV_DIVMU_RECURSE > 3
              case 4 :
                   LDivEqMU_Helper<26,4,rs,M1,M2>::call(m1,m2);
                   break;
#endif
#if TMV_DIVMU_RECURSE > 4
              default :
                   LDivEqMU_Helper<algo2,UNKNOWN,rs,M1,M2>::call(m1,m2);
#endif
            }
        }
    };

    // algo 27: Split the LowerTriMatrix into 3 sections and recurse
    // the calculation on each of them:
    // ( A 0 ) ( D ) = ( F )
    // ( B C ) ( E )   ( G )
    // G = BD + CE
    // F = AD
    // Solving for D,E yields:
    // D /= A
    // E -= BD
    // E /= C
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<27,cs,rs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            const int M = cs==UNKNOWN ? int(m1.colsize()) : cs;
#ifdef PRINTALGO_DivU
            const int N = rs==UNKNOWN ? int(m1.rowsize()) : rs;
            std::cout<<"LDivEqMU algo 27: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif

            const bool rr = M1::_rowmajor && M2::_rowmajor;
            const bool rc = M1::_rowmajor && M2::_colmajor;
            const bool cx = M1::_colmajor;
            const int algo2 = 
                cs == 0 ? 0 :
                cs == 1 ? 1 :
                (cs != UNKNOWN && cs > TMV_DIVMU_RECURSE) ? 0 :
                TMV_DIVMU_RECURSE == 1 ? 1 :
#ifdef TMV_DIVMU_CLEANUP
                cs == UNKNOWN ? 26 :
#endif
                (cs != UNKNOWN && cs <= TMV_DIVMU_RECURSE) ? 26 :
                rr ? 22 : rc ? 23 : cx ? 21 : 23;
            const int algo3 =  // The algorithm for M > 16
                cs == UNKNOWN ? 27 : 
#ifdef TMV_DIVMU_INLINE_MM
                cs <= 16 ? 0 : 
#endif
                cs > TMV_DIVMU_RECURSE ? 27 :
                0;
            const int algo4 =  // The algorithm for MultMM
                cs == UNKNOWN ? -2 : 
#ifdef TMV_DIVMU_INLINE_MM
                cs <= 16 ? 0 : 
#endif
                cs > TMV_DIVMU_RECURSE ? -3 :
                0;
#ifdef TMV_DIVMU_INLINE_MM
            const int algo3b =  // The algorithm for M > RECURSE
                cs == UNKNOWN || 
                (cs > TMV_DIVMU_RECURSE && cs <= 16) ? 27 : 0;
            const int algo4b =  // The algorithm for MultMM
                cs == UNKNOWN || 
                (cs > TMV_DIVMU_RECURSE && cs <= 16) ? -4 : 0;
#endif
#ifdef PRINTALGO_DivU
            std::cout<<"algo2,3,4 = "<<algo2<<"  "<<algo3<<"  "<<algo4<<std::endl;
#ifdef TMV_DIVMU_INLINE_MM
            std::cout<<"algo3b,4b = "<<algo3b<<"  "<<algo4b<<std::endl;
#endif
#endif

            typedef typename M2::real_type RT;
            typedef typename M2::const_subtrimatrix_type M2a;
            typedef typename M2::const_submatrix_type M2b;
            typedef typename M1::const_rowrange_type M1rc;
            typedef typename M1::rowrange_type M1r;
            const Scaling<-1,RT> mone;

            if (M > TMV_DIVMU_RECURSE) {
                const int Mx = M > 16 ? ((((M-1)>>5)+1)<<4) : (M>>1);
                // (If M > 16, round M/2 up to a multiple of 16.)
                const int csx = IntTraits<cs>::half_roundup;
                const int csy = IntTraits2<cs,csx>::diff;

                M2a A = m2.cSubTriMatrix(0,Mx);
                M2b B = m2.cSubMatrix(Mx,M,0,Mx);
                M2a C = m2.cSubTriMatrix(Mx,M);
                M1r D = m1.cRowRange(0,Mx);
                M1r E = m1.cRowRange(Mx,M);

#ifdef TMV_DIVMU_INLINE_MM
                if (M > 16) {
                    // For large M, make sure to use good MultMM algo
#endif
                    LDivEqMU_Helper<algo3,csx,rs,M1r,M2a>::call(D,A);
                    MultMM_Helper<algo4,csy,rs,csx,true,-1,RT,M2b,M1rc,M1r>::
                        call(mone,B,D,E);
                    LDivEqMU_Helper<algo3,csy,rs,M1r,M2a>::call(E,C);
#ifdef TMV_DIVMU_INLINE_MM
                } else {
                    // For smaller M, do the no branching algorithm
                    LDivEqMU_Helper<algo3b,csx,rs,M1r,M2a>::call(D,A);
                    MultMM_Helper<algo4b,csy,rs,csx,true,-1,RT,M2b,M1rc,M1r>::
                        call(mone,B,D,E);
                    LDivEqMU_Helper<algo3b,csy,rs,M1r,M2a>::call(E,C);
                }
#endif
            } else {
                LDivEqMU_Helper<algo2,cs,rs,M1,M2>::call(m1,m2);
            }
        }
    };

    // algo 31: Determine which algorithm to use based on the runtime
    // knowledge of the sizes.
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<31,cs,rs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            const int M = cs==UNKNOWN ? int(m1.colsize()) : cs;
#if defined(PRINTALGO_DivU) || defined(_OPENMP)
            const int N = rs==UNKNOWN ? int(m1.rowsize()) : rs;
#endif
#ifdef PRINTALGO_DivU
            std::cout<<"LDivEqMU algo 31: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            const bool upper = M2::_upper;
            const bool rr = M1::_rowmajor && M2::_rowmajor;
            const bool rc = M1::_rowmajor && M2::_colmajor;
            const bool cx = M1::_colmajor;
            const int algo2 = 
                TMV_DIVMU_RECURSE == 1 ? 1 :
#ifdef TMV_DIVMU_CLEANUP
                true ? ( upper ? 16 : 26 ) :
#endif
                upper ? ( rr ? 12 : rc ? 13 : cx ? 11 : 13 ) :
                ( rr ? 22 : rc ? 23 : cx ? 21 : 23 );
            const int algo3 = upper ? 17 : 27;

#ifdef _OPENMP
            const int Mc = M < 16 ? 1 : M>>4; // M/16
            const int Nc = N < 16 ? 1 : N>>4; // N/16
#endif

            // Put the small matrix option first, so it doesn't have to 
            // go through a bunch of if/else statements.  For large matrices,
            // all these if/else's don't matter for the total time.
            if (M <= TMV_DIVMU_RECURSE)
                LDivEqMU_Helper<algo2,cs,rs,M1,M2>::call(m1,m2);
#ifdef _OPENMP
            else if (N >= 64 && Mc*Mc*Nc > TMV_DIVMU_OMP_THRESH)
                LDivEqMU_Helper<36,cs,rs,M1,M2>::call(m1,m2);
#endif
            else
                LDivEqMU_Helper<algo3,cs,rs,M1,M2>::call(m1,m2);
        }
    };

#ifdef _OPENMP
    // algo 36: Split problem into smaller parts with OpenMP for 
    // parallelization.
    // We take a pretty simple approach here, and just split up m2 and m1
    // matrices by columns and let each thread do a single matrix.
    // Then each thread calls algo 17 or 27 to calculate its product.
    // Also, we require that all but the last thread has a column width
    // that is a multiple of 16.  This way we get the maximum advantage from
    // our blocking structure while keeping the threads as balanced as 
    // possible.
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<36,cs,rs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            const int N = rs==UNKNOWN ? int(m1.rowsize()) : rs;
#ifdef PRINTALGO_DivU
            const int M = cs==UNKNOWN ? int(m1.colsize()) : cs;
            std::cout<<"LDivEqMU algo 36: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            const bool upper = M2::_upper;
            const int algo2 = upper ? 17 : 27;
#ifdef PRINTALGO_DivU_OMP
            std::ofstream fout("omp.out");
#endif
#pragma omp parallel
            {
                int num_threads = omp_get_num_threads();
                int mythread = omp_get_thread_num();
#ifdef PRINTALGO_DivU_OMP
#pragma omp critical
                {
                    fout<<"thread "<<mythread<<"/"<<num_threads<<std::endl;
                }
#endif
                if (num_threads == 1) {
#ifdef PRINTALGO_DivU_OMP
#pragma omp critical
                    {
                        fout<<"thread "<<mythread<<"/"<<num_threads;
                        fout<<"\nonly 1 thread"<<std::endl;
                    }
#endif
                    LDivEqMU_Helper<algo2,cs,rs,M1,M2>::call(m1,m2);
                } else {
                    int Nx = N / num_threads;
                    Nx = ((((Nx-1)>>4)+1)<<4); 
                    int j1 = mythread * Nx;
                    int j2 = (mythread+1) * Nx;
                    if (j2 > N || mythread == num_threads-1) j2 = N;
#ifdef PRINTALGO_DivU_OMP
#pragma omp critical
                    {
                        fout<<"thread "<<mythread<<"/"<<num_threads;
                        fout<<"Nx = "<<Nx<<std::endl;
                        fout<<"j1 = "<<j1<<std::endl;
                        fout<<"j2 = "<<j2<<std::endl;
                    }
#endif
                    if (j1 < N)  {
                        typedef typename M1::colrange_type M1c;
                        const int rsx = UNKNOWN; 
                        M1c m1c = m1.cColRange(j1,j2);
#ifdef PRINTALGO_DivU_OMP
#pragma omp critical
                        {
                            fout<<"thread "<<mythread<<"/"<<num_threads;
                            fout<<"\nm2 = "<<m2<<std::endl;
                            fout<<"m1c = "<<m1c<<std::endl;
                        }
#endif
                        LDivEqMU_Helper<algo2,cs,rsx,M1c,M2>::call(m1c,m2);
#ifdef PRINTALGO_DivU_OMP
#pragma omp critical
                        {
                            fout<<"thread "<<mythread<<"/"<<num_threads;
                            fout<<"\nm1c => "<<m1c<<std::endl;
                        }
#endif
                    }
                }
            }
        }
    };
#endif

    // algo 38: Unknown cs, check if M is small
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<38,cs,rs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            TMVStaticAssert(cs == UNKNOWN);
            const int M = cs==UNKNOWN ? int(m1.colsize()) : cs;
#ifdef PRINTALGO_DivU
            const int N = rs==UNKNOWN ? int(m1.rowsize()) : rs;
            std::cout<<"LDivEqMU algo 38: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            const bool upper = M2::_upper;
            const int algo2 = upper ? 16 : 26;

            if (M <= 3) {
                // then it is worth figuring out what M is.
                switch (M) {
                  case 0 :
                       // do nothing
                       break;
                  case 1 :
                       LDivEqMU_Helper<1,1,rs,M1,M2>::call(m1,m2);
                       break;
                  case 2 :
                       LDivEqMU_Helper<algo2,2,rs,M1,M2>::call(m1,m2);
                       break;
                  case 3 :
                       LDivEqMU_Helper<algo2,3,rs,M1,M2>::call(m1,m2);
                       break;
                }
            } else 
                LDivEqMU_Helper<31,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo 81: copy m2
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<81,cs,rs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DivU
            const int M = cs == UNKNOWN ? int(m1.colsize()) : cs;
            const int N = rs == UNKNOWN ? int(m1.rowsize()) : rs;
            std::cout<<"LDivEqMU algo 81: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            NoAliasTriLDivEq(m1,m2.copy());
        }
    };

    // algo 84: copy m1
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<84,cs,rs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DivU
            const int M = cs == UNKNOWN ? int(m1.colsize()) : cs;
            const int N = rs == UNKNOWN ? int(m1.rowsize()) : rs;
            std::cout<<"LDivEqMU algo 84: M,N,cs,rs = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<std::endl;
#endif
            typedef typename M1::value_type T1;
            // I think rowmajor is usually better, since it only inverts
            // the diagonal elements of m2 once each.
            const bool rm = true; 
            typedef typename MCopyHelper<T1,Rec,cs,rs,rm,false>::type M1c;
            M1c m1c = m1;
            NoAliasTriLDivEq(m1c,m2);
            NoAliasCopy(m1c,m1);
        }
    };

    // algo 90: call inst
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<90,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, const M2& m2)
        { InstTriLDivEq(m1.xView(),m2.xView()); }
    };

    // algo 91: call inst alias
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<91,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, const M2& m2)
        { InstAliasTriLDivEq(m1.xView(),m2.xView()); }
    };

    // algo 97: Conjugate
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<97,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, const M2& m2)
        {
            typedef typename M1::conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            LDivEqMU_Helper<-2,cs,rs,M1c,M2c>::call(m1c,m2c);
        }
    };

    // algo 197: Conjugate
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<197,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, const M2& m2)
        {
            typedef typename M1::conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            LDivEqMU_Helper<99,cs,rs,M1c,M2c>::call(m1c,m2c);
        }
    };

    // algo 98: Inline check for aliases
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<98,cs,rs,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            if (!SameStorage(m1,m2)) {
                // No aliasing 
                LDivEqMU_Helper<-2,cs,rs,M1,M2>::call(m1,m2);
            } else if (m1.colsize() <= 2*m1.rowsize()) {
                // copy m2
                LDivEqMU_Helper<81,cs,rs,M1,M2>::call(m1,m2);
            } else {
                // Use temporary for m1
                LDivEqMU_Helper<84,cs,rs,M1,M2>::call(m1,m2);
            }
        }
    };

    // algo 99: Check for aliases
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<99,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, const M2& m2)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst = 
                (cs == UNKNOWN || cs > 16) &&
                (rs == UNKNOWN || rs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T1>::isinst;
            const int algo = 
                ( cs == 0 || rs == 0 ) ? 0 :
                M1::_conj ? 197 :
                inst ? 91 : 
                98;
            LDivEqMU_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo -4: No branches or copies
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<-4,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, const M2& m2)
        {
            const bool upper = M2::_upper;
            const int algo = 
                ( cs == 0 || rs == 0 ) ? 0 :
                cs == 1 ? 1 :
                rs == 1 ? 2 :
                upper ? 

                cs != UNKNOWN && cs <= 5 ? 16 :
                17  :

                cs != UNKNOWN && cs <= 5 ? 26 :
                27;

            LDivEqMU_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<-3,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, const M2& m2)
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
            // 16 = invert m2 explicitly and do MultEq op
            // 17 = split trimatrix into 3 submatrices
            //
            // LowerTri:
            // 21 = loop over n: MultUV
            // 22 = loop over m: MultMV
            // 23 = loop over k: Rank1
            // 26 = invert m2 explicitly and do MultEq op
            // 27 = split trimatrix into 3 submatrices
            // 
            // Overall drivers:
            // 31 = Choose algorithm based on the (runtime) size
            // 36 = Parallelize using openmp 
            // 38 = Check if M is small

            const bool upper = M2::_upper;
#if 0
            //const int algo = upper ? 17 : 27;
            const int algo = 36;
#else
            const bool rr = M1::_rowmajor && M2::_rowmajor;
            const bool rc = M1::_rowmajor && M2::_colmajor;
            const bool cx = M1::_colmajor;
#ifdef _OPENMP
            const int Mc = cs == UNKNOWN ? UNKNOWN : (cs < 16) ? 1 : (cs>>4);
            const int Nc = rs == UNKNOWN ? UNKNOWN : (rs < 16) ? 1 : (rs>>4);
            const int McMcNc = IntTraits2<IntTraits2<Mc,Mc>::prod,Nc>::prod;
#endif
            const int algo = 
                ( cs == 0 || rs == 0 ) ? 0 :
                cs == 1 ? 1 :
                rs == 1 ? 2 :
                upper ? 

                TMV_OPT == 0 ? ( rr ? 12 : rc ? 13 : cx ? 11 : 13 ) :
#ifdef TMV_DIVMU_SMALL
                cs == UNKNOWN ? 38 : 
#endif
                cs == UNKNOWN ? 31 : 
                cs <= 3 && (rs == UNKNOWN || rs > 3) ? 16 :
                rs == UNKNOWN ? 31 :
#ifdef _OPENMP
                (rs >= 64 && McMcNc >= TMV_DIVMU_OMP_THRESH) ? 36 :
#endif
                17  :

                // lowertri
                TMV_OPT == 0 ? ( rr ? 22 : rc ? 23 : cx ? 21 : 23 ) :
#ifdef TMV_DIVMU_SMALL
                cs == UNKNOWN ? 38 : 
#endif
                cs == UNKNOWN ? 31 : 
                cs <= 3 && (rs == UNKNOWN || rs > 3) ? 26 :
                rs == UNKNOWN ? 31 :
#ifdef _OPENMP
                (rs >= 64 && McMcNc >= TMV_DIVMU_OMP_THRESH) ? 36 :
#endif
                27 ;
#endif
#ifdef PRINTALGO_DivU
            const int M = cs==UNKNOWN ? int(m1.colsize()) : cs;
            const int N = rs==UNKNOWN ? int(m1.rowsize()) : rs;
            std::cout<<"InlineLDivEqMU: \n";
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"M = "<<M<<"  N = "<<N<<std::endl;
            std::cout<<"cs = "<<cs<<"  rs = "<<rs<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
            //std::cout<<"m1 = "<<m1<<std::endl;
            //std::cout<<"m2 = "<<m2<<std::endl;
#endif
            if (m2.isSingular()) ThrowSingular("TriMatrix");
            LDivEqMU_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
#ifdef PRINTALGO_DivU
            //std::cout<<"m1 => "<<m1<<std::endl;
#endif
        }
    };

    // algo -402: Same as algo -2, but use -4 if no Inst version
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<-402,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, const M2& m2)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst = 
                (cs == UNKNOWN || cs > 16) &&
                (rs == UNKNOWN || rs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T1>::isinst;
            const int algo = 
                ( cs == 0 || rs == 0 ) ? 0 :
                cs == 1 ? 1 :
                rs == 1 ? 2 :
                inst ? 90 : 
                -4;
            LDivEqMU_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<-2,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, const M2& m2)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst = 
                (cs == UNKNOWN || cs > 16) &&
                (rs == UNKNOWN || rs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T1>::isinst;
            const int algo = 
                ( cs == 0 || rs == 0 ) ? 0 :
                cs == 1 ? 201 :
                rs == 1 ? 202 :
                M1::_conj ? 97 :
                inst ? 90 : 
                -3;
            LDivEqMU_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo -1: Check for aliases?
    template <int cs, int rs, class M1, class M2>
    struct LDivEqMU_Helper<-1,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, const M2& m2)
        {
            const int algo = 
                ( cs == 0 || rs == 0 ) ? 0 :
                cs == 1 ? 101 :
                rs == 1 ? 102 :
                M1::_checkalias ? 99 : 
                -2;
            LDivEqMU_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    template <class M1, class M2>
    static inline void TriLDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const BaseMatrix_Tri<M2>& m2)
    {
        TMVStaticAssert((Sizes<M2::_size,M1::_colsize>::same));
        TMVAssert(m2.size() == m1.colsize());

        const int cs = Sizes<M1::_colsize,M2::_size>::size;
        const int rs = M1::_rowsize;
        typedef typename M1::cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        TMV_MAYBE_REF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        LDivEqMU_Helper<-1,cs,rs,M1v,M2v>::call(m1v,m2v);
    }
    template <class M1, class M2>
    static inline void NoAliasTriLDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const BaseMatrix_Tri<M2>& m2)
    {
        TMVStaticAssert((Sizes<M2::_size,M1::_colsize>::same));
        TMVAssert(m2.size() == m1.colsize());

        const int cs = Sizes<M1::_colsize,M2::_size>::size;
        const int rs = M1::_rowsize;
        typedef typename M1::cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        TMV_MAYBE_REF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        LDivEqMU_Helper<-2,cs,rs,M1v,M2v>::call(m1v,m2v);
    }
    template <class M1, class M2>
    static inline void InlineTriLDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const BaseMatrix_Tri<M2>& m2)
    {
        TMVStaticAssert((Sizes<M2::_size,M1::_colsize>::same));
        TMVAssert(m2.size() == m1.colsize());

        const int cs = Sizes<M1::_colsize,M2::_size>::size;
        const int rs = M1::_rowsize;
        typedef typename M1::cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        TMV_MAYBE_REF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        LDivEqMU_Helper<-3,cs,rs,M1v,M2v>::call(m1v,m2v);
    }
    template <class M1, class M2>
    static inline void InlineAliasTriLDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const BaseMatrix_Tri<M2>& m2)
    {
        TMVStaticAssert((Sizes<M2::_size,M1::_colsize>::same));
        TMVAssert(m2.size() == m1.colsize());

        const int cs = Sizes<M1::_colsize,M2::_size>::size;
        const int rs = M1::_rowsize;
        typedef typename M1::cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        TMV_MAYBE_REF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        LDivEqMU_Helper<98,cs,rs,M1v,M2v>::call(m1v,m2v);
    }
    template <class M1, class M2>
    static inline void AliasTriLDivEq(
        BaseMatrix_Rec_Mutable<M1>& m1, const BaseMatrix_Tri<M2>& m2)
    {
        TMVStaticAssert((Sizes<M2::_size,M1::_colsize>::same));
        TMVAssert(m2.size() == m1.colsize());

        const int cs = Sizes<M1::_colsize,M2::_size>::size;
        const int rs = M1::_rowsize;
        typedef typename M1::cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        TMV_MAYBE_REF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        LDivEqMU_Helper<99,cs,rs,M1v,M2v>::call(m1v,m2v);
    }

} // namespace tmv

#endif 
