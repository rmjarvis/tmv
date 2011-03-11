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


#ifndef TMV_DivUU_H
#define TMV_DivUU_H

#include "TMV_MultMM.h"
#include "TMV_MultUV.h"
#include "TMV_MultUM.h"
#include "TMV_CopyU.h"
#include "TMV_DivVU.h"
#include "TMV_DivMU.h"
#include "TMV_Rank1VVM.h"

#ifdef PRINTALGO_DivU
#include <iostream>
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

// Q2 is the maximum size to stop recursing
#define TMV_Q2 1

namespace tmv {

    // Defined below:
    template <class M1, class M2>
    static void LDivEq(
        BaseMatrix_Tri_Mutable<M1>& m1, const BaseMatrix_Tri<M2>& m2);
    template <class M1, class M2>
    static void NoAliasLDivEq(
        BaseMatrix_Tri_Mutable<M1>& m1, const BaseMatrix_Tri<M2>& m2);
    template <class M1, class M2>
    static void InlineLDivEq(
        BaseMatrix_Tri_Mutable<M1>& m1, const BaseMatrix_Tri<M2>& m2);
    template <class M1, class M2>
    static void AliasLDivEq(
        BaseMatrix_Tri_Mutable<M1>& m1, const BaseMatrix_Tri<M2>& m2);

    // Defined in TMV_DivUU.cpp
    template <class T1, class T2, bool C2>
    void InstLDivEq(
        UpperTriMatrixView<T1> m1,
        const ConstUpperTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2);
    template <class T1, class T2, bool C2>
    void InstLDivEq(
        LowerTriMatrixView<T1> m1,
        const ConstLowerTriMatrixView<T2,UnknownDiag,UNKNOWN,UNKNOWN,C2>& m2);


    template <int algo, int s, class M1, class M2>
    struct LDivEqUU_Helper;

    // algo 0: Trivial, nothing to do.
    template <int s, class M1, class M2>
    struct LDivEqUU_Helper<0,s,M1,M2>
    { static void call(M1& , const M2& ) {} };

    // algo 1: N == 1, so reduces to scalar quotient
    template <int s, class M1, class M2>
    struct LDivEqUU_Helper<1,s,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DivU
            std::cout<<"LDivEqUU algo 1: N,s = "<<1<<','<<s<<std::endl;
#endif
            Maybe<!M1::_unit && !M2::_unit>::invscale(
                m1.ref(0,0) , m2.cref(0,0));
        }
    };

    // algo 11: UpperTri loop over n
    template <int s, class M1, class M2>
    struct LDivEqUU_Helper<11,s,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            const int N = s==UNKNOWN ? int(m1.size()) : s;
#ifdef PRINTALGO_DivU
            std::cout<<"LDivEqUU algo 11: N,s = "<<N<<','<<s<<std::endl;
#endif
            typedef typename M1::col_sub_type M1c;
            typedef typename M2::const_subtrimatrix_type M2s;
            const int xx = UNKNOWN;

            for(int j=0;j<N;++j) {
                // m1.col(j,0,j) /= m2.subTriMatrix(0,j)
                M1c m1j = m1.get_col(j,0,j);
                M2s m2s = m2.cSubTriMatrix(0,j);
                LDivEqVU_Helper<-4,xx,M1c,M2s>::call(m1j,m2s);
            }
        }
    };

    // algo 12: UpperTri loop over m
    template <int s, class M1, class M2>
    struct LDivEqUU_Helper<12,s,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            const int N = s==UNKNOWN ? int(m1.size()) : s;
#ifdef PRINTALGO_DivU
            std::cout<<"LDivEqUU algo 12: N,s = "<<N<<','<<s<<std::endl;
#endif
            const bool u2 = M2::_unit;
            typedef typename M2::real_type RT;
            typedef typename M2::value_type T2;
            typedef typename TypeSelect<u2,RT,T2>::type XT2;
            typedef typename M2::const_row_sub_type M2r;
            typedef typename M1::row_sub_type M1r;
            typedef typename M1::subtrimatrix_type::const_transpose_type M1st;
            const int ix2 = u2 ? 1 : 0;
            const int xx = UNKNOWN;
            const Scaling<-1,RT> mone;

            for(int i=N;i--;) {
                // m1.row(i,i+1,N) -= m2.row(i,i+1,N) * m1.subTriMatrix(i+1,N)
                // m1.row(i,i,N) /= m2(i,i)
                M1r m1a = m1.get_row(i,i+1,N);
                M1r m1b = m1.get_row(i,i,N);
                M1st m1st = m1.cSubTriMatrix(i+1,N).transpose();
                M2r m2i = m2.get_row(i,i+1,N);
                MultUV_Helper<-4,xx,true,-1,RT,M1st,M2r,M1r>::call(
                    mone,m1st,m2i,m1a);
                const Scaling<ix2,XT2> invii(
                    Maybe<!u2>::invprod( m2.cref(i,i) , RT(1) ));
                ScaleV_Helper<-4,xx,ix2,T2,M1r>::call(invii,m1b);
            }
        }
    };

    // algo 13: UpperTri loop over k
    template <int s, class M1, class M2>
    struct LDivEqUU_Helper<13,s,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            const int N = s==UNKNOWN ? int(m1.size()) : s;
#ifdef PRINTALGO_DivU
            std::cout<<"LDivEqUU algo 13: N,s = "<<N<<','<<s<<std::endl;
#endif
            const bool u2 = M2::_unit;
            typedef typename M2::real_type RT;
            typedef typename M2::value_type T2;
            typedef typename TypeSelect<u2,RT,T2>::type XT2;
            typedef typename M1::row_sub_type M1r;
            typedef typename M1::submatrix_type M1s;
            typedef typename M1::const_row_sub_type M1rc;
            typedef typename M2::const_col_sub_type M2c;
            const int ix2 = u2 ? 1 : 0;
            const int xx = UNKNOWN;
            const Scaling<-1,RT> mone;

            for(int k=N;k--;) {
                // m1.row(k,k,N) /= m2(k,k)
                // m1.subMatrix(0,k,k,N) -= m2.col(k,0,k) ^ m1.row(k,k,N)
                M1r m1k = m1.get_row(k,k,N);
                M1s m1s = m1.cSubMatrix(0,k,k,N);
                M2c m2k = m2.get_col(k,0,k);
                const Scaling<ix2,XT2> invkk(
                    Maybe<!u2>::invprod( m2.cref(k,k) , RT(1) ));
                ScaleV_Helper<-4,xx,ix2,XT2,M1r>::call(invkk,m1k);
                Rank1VVM_Helper<-4,xx,xx,true,-1,RT,M2c,M1rc,M1s>::call(
                    mone,m2k,m1k,m1s);
            }
        }
    };

    // algo 17: Split the UpperTriMatrix into 3 sections and recurse
    // the calculation on each of them:
    // ( A B ) ( D E ) = ( G H )
    // ( 0 C ) ( 0 F )   ( 0 I )
    // G = AD 
    // H = AE + BF
    // I = CF
    // Solving for D,E,F yields:
    // F /= C
    // E -= BF
    // E /= A
    // D /= A
    template <int s, class M1, class M2>
    struct LDivEqUU_Helper<17,s,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            const int N = s==UNKNOWN ? int(m1.size()) : s;
#ifdef PRINTALGO_DivU
            std::cout<<"LDivEqUU algo 17: N,s = "<<N<<','<<s<<std::endl;
#endif

            const bool rr = M1::_rowmajor && M2::_rowmajor;
            const bool rc = M1::_rowmajor && M2::_colmajor;
            const bool cx = M1::_colmajor;
            const int algo2 = 
                s == 0 ? 0 :
                s == 1 ? 1 :
                (s != UNKNOWN && s > TMV_Q2) ? 0 :
                TMV_Q2 == 1 ? 1 :
                rr ? 12 : rc ? 13 : cx ? 11 : 13;
            const int algo3 =  // The algorithm for N > 32
                s == UNKNOWN || s > 32 ? 17 : 0;
            const int algo4 =  // The algorithm for MultMM
                s == UNKNOWN ? -2 : s > 32 ? -3 : 0;
#if TMV_Q2 < 32
            const int algo3b =  // The algorithm for N > Q2
                s == UNKNOWN || s > TMV_Q2 ? 17 : 0;
            const int algo4b =  // The algorithm for MultMM
                s == UNKNOWN || s > TMV_Q2 ? -4 : 0;
#endif

            typedef typename M2::real_type RT;
            typedef typename M2::const_subtrimatrix_type M2a;
            typedef typename M2::const_submatrix_type M2b;
            typedef typename M2b::const_transpose_type M2bt;
            typedef typename M1::subtrimatrix_type M1a;
            typedef typename M1a::const_transpose_type M1at;
            typedef typename M1::submatrix_type M1b;
            typedef typename M1b::transpose_type M1bt;
            const Scaling<-1,RT> mone;

            if (N > TMV_Q2) { 
                const int Nx = N > 16 ? ((((N-1)>>5)+1)<<4) : (N>>1);
                // (If N > 16, round N/2 up to a multiple of 16.)
                const int sx = IntTraits<s>::half_roundup;
                const int sy = IntTraits2<s,sx>::diff;

                M2a A = m2.cSubTriMatrix(0,Nx);
                M2bt Bt = m2.cSubMatrix(0,Nx,Nx,N).transpose();
                M2a C = m2.cSubTriMatrix(Nx,N);
                M1a D = m1.cSubTriMatrix(0,Nx);
                M1b E = m1.cSubMatrix(0,Nx,Nx,N);
                M1bt Et = E.transpose();
                M1a F = m1.cSubTriMatrix(Nx,N);
                M1at Ft = F.transpose();

#if TMV_Q2 < 32
                if (N > 32) {
                    // For large N, make sure to use good MultMM algo
#endif
                    LDivEqUU_Helper<algo3,sy,M1a,M2a>::call(F,C);
                    MultUM_Helper<algo4,sy,sx,true,-1,RT,M1at,M2bt,M1bt>::
                        call(mone,Ft,Bt,Et);
                    LDivEqMU_Helper<algo4,sx,sy,M1b,M2a>::call(E,A);
                    LDivEqUU_Helper<algo3,sx,M1a,M2a>::call(D,A);
#if TMV_Q2 < 32
                } else {
                    // For smaller N, do the no branching algorithms
                    LDivEqUU_Helper<algo3b,sy,M1a,M2a>::call(F,C);
                    MultUM_Helper<algo4b,sy,sx,true,-1,RT,M1at,M2bt,M1bt>::
                        call(mone,Ft,Bt,Et);
                    LDivEqMU_Helper<algo4b,sx,sy,M1b,M2a>::call(E,A);
                    LDivEqUU_Helper<algo3b,sx,M1a,M2a>::call(D,A);
                }
#endif
            }
            else LDivEqUU_Helper<algo2,s,M1,M2>::call(m1,m2);
        }
    };

    // algo 21: LowerTri loop over n
    template <int s, class M1, class M2>
    struct LDivEqUU_Helper<21,s,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            const int N = s==UNKNOWN ? int(m1.size()) : s;
#ifdef PRINTALGO_DivU
            std::cout<<"LDivEqUU algo 21: N,s = "<<N<<','<<s<<std::endl;
#endif
            typedef typename M1::col_sub_type M1c;
            typedef typename M2::const_subtrimatrix_type M2s;
            for(int j=0;j<N;++j) {
                // m1.col(j,j,N) /= m2.subTriMatrix(j,N)
                M1c m1j = m1.get_col(j,j,N);
                M2s m2s = m2.cSubTriMatrix(j,N);
                LDivEqVU_Helper<-4,s,M1c,M2s>::call(m1j,m2s);
            }
        }
    };

    // algo 22: LowerTri loop over m
    template <int s, class M1, class M2>
    struct LDivEqUU_Helper<22,s,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            const int N = s==UNKNOWN ? int(m1.size()) : s;
#ifdef PRINTALGO_DivU
            std::cout<<"LDivEqUU algo 22: N,s = "<<N<<','<<s<<std::endl;
#endif
            const bool u2 = M2::_unit;
            typedef typename M2::real_type RT;
            typedef typename M2::value_type T2;
            typedef typename TypeSelect<u2,RT,T2>::type XT2;
            typedef typename M2::const_row_sub_type M2r;
            typedef typename M1::subtrimatrix_type::const_transpose_type M1st;
            typedef typename M1::row_sub_type M1r;
            const int ix2 = u2 ? 1 : 0;
            const int xx = UNKNOWN;
            const Scaling<-1,RT> mone;

            for(int i=0;i<N;++i) {
                // m1.row(i,0,i) -= m2.row(i,0,i) * m1.subTriMatrix(0,i)
                // m1.row(i,0,i+1) /= m2(i,i)
                M1r m1a = m1.get_row(i,0,i);
                M1r m1b = m1.get_row(i,0,i+1);
                M2r m2i = m2.get_row(i,0,i);
                M1st m1st = m1.cSubTriMatrix(0,i).transpose();
                MultUV_Helper<-4,xx,true,-1,RT,M1st,M2r,M1r>::call(
                    mone,m1st,m2i,m1a);
                const Scaling<ix2,XT2> invii(
                    Maybe<!u2>::invprod( m2.cref(i,i) , RT(1) ));
                ScaleV_Helper<-4,xx,ix2,XT2,M1r>::call(invii,m1b);
            }
        }
    };

    // algo 23: LowerTri loop over k
    template <int s, class M1, class M2>
    struct LDivEqUU_Helper<23,s,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            const int N = s==UNKNOWN ? int(m1.size()) : s;
#ifdef PRINTALGO_DivU
            std::cout<<"LDivEqUU algo 23: N,s = "<<N<<','<<s<<std::endl;
#endif
            const bool u2 = M2::_unit;
            typedef typename M2::real_type RT;
            typedef typename M2::value_type T2;
            typedef typename TypeSelect<u2,RT,T2>::type XT2;
            typedef typename M1::row_sub_type M1r;
            typedef typename M1::const_row_sub_type M1rc;
            typedef typename M1::submatrix_type M1s;
            typedef typename M2::const_col_sub_type M2c;
            const int ix2 = u2 ? 1 : 0;
            const int xx = UNKNOWN;
            const Scaling<-1,RT> mone;

            for(int k=0;k<N;++k) {
                // m1.row(k,0,k+1) /= m2(k,k)
                // m1.subMatrix(k+1,N,0,k+1) -= m2.col(k,k+1,N) ^ m1.row(k)
                M1r m1k = m1.get_row(k,0,k+1);
                M1s m1s = m1.cSubMatrix(k+1,N,0,k+1);
                M2c m2k = m2.get_col(k,k+1,N);
                const Scaling<ix2,XT2> invkk(
                    Maybe<!u2>::invprod( m2.cref(k,k) , RT(1) ));
                ScaleV_Helper<-4,xx,ix2,XT2,M1r>::call(invkk,m1k);
                Rank1VVM_Helper<-4,xx,xx,true,-1,RT,M2c,M1rc,M1s>::call(
                    mone,m2k,m1k,m1s);
            }
        }
    };

    // algo 27: Split the LowerTriMatrix into 3 sections and recurse
    // the calculation on each of them:
    // ( A 0 ) ( D 0 ) = ( G 0 )
    // ( B C ) ( E F )   ( H I )
    // G = AD
    // H = BD + CE
    // I = CF
    // Solving for D,E,F yields:
    // D /= A
    // E -= BD
    // E /= C
    // F /= C
    template <int s, class M1, class M2>
    struct LDivEqUU_Helper<27,s,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            const int N = s==UNKNOWN ? int(m1.size()) : s;
#ifdef PRINTALGO_DivU
            std::cout<<"LDivEqUU algo 27: N,s = "<<N<<','<<s<<std::endl;
#endif

            const bool rr = M1::_rowmajor && M2::_rowmajor;
            const bool rc = M1::_rowmajor && M2::_colmajor;
            const bool cx = M1::_colmajor;
            const int algo2 = 
                s == 0 ? 0 :
                s == 1 ? 1 :
                (s != UNKNOWN && s > TMV_Q2) ? 0 :
                TMV_Q2 == 1 ? 1 :
                rr ? 22 : rc ? 23 : cx ? 21 : 23;
            const int algo3 =  // The algorithm for N > Q2
                s == UNKNOWN || s > 32 ? 27 : 0;
            const int algo4 =  // The algorithm for MultMM
                s == UNKNOWN ? -2 : s > 32 ? -3 : 0;
#if TMV_Q2 < 32
            const int algo3b =  // The algorithm for N > Q2
                s == UNKNOWN || s > TMV_Q2 ? 27 : 0;
            const int algo4b =  // The algorithm for MultMM
                s == UNKNOWN || s > TMV_Q2 ? -4 : 0;
#endif

            typedef typename M2::real_type RT;
            typedef typename M2::const_subtrimatrix_type M2a;
            typedef typename M2::const_submatrix_type M2b;
            typedef typename M2b::const_transpose_type M2bt;
            typedef typename M1::subtrimatrix_type M1a;
            typedef typename M1a::const_transpose_type M1at;
            typedef typename M1::submatrix_type M1b;
            typedef typename M1b::transpose_type M1bt;
            const Scaling<-1,RT> mone;

            if (N > TMV_Q2) {
                const int Nx = N > 16 ? ((((N-1)>>5)+1)<<4) : (N>>1);
                // (If N > 16, round M/2 up to a multiple of 16.)
                const int sx = IntTraits<s>::half_roundup;
                const int sy = IntTraits2<s,sx>::diff;

                M2a A = m2.cSubTriMatrix(0,Nx);
                M2b B = m2.cSubMatrix(Nx,N,0,Nx);
                M2bt Bt = B.transpose();
                M2a C = m2.cSubTriMatrix(Nx,N);
                M1a D = m1.cSubTriMatrix(0,Nx);
                M1at Dt = D.transpose();
                M1b E = m1.cSubMatrix(Nx,N,0,Nx);
                M1bt Et = E.transpose();
                M1a F = m1.cSubTriMatrix(Nx,N);

#if TMV_Q2 < 32
                if (N > 32) {
                    // For large N, make sure to use good MultMM algo
#endif
                    LDivEqUU_Helper<algo3,sx,M1a,M2a>::call(D,A);
                    MultUM_Helper<algo4,sx,sy,true,-1,RT,M1at,M2bt,M1bt>::
                        call(mone,Dt,Bt,Et);
                    LDivEqMU_Helper<algo3,sy,sx,M1b,M2a>::call(E,C);
                    LDivEqUU_Helper<algo3,sy,M1a,M2a>::call(F,C);
#if TMV_Q2 < 32
                } else {
                    // For smaller N, do the no branching algorithms
                    LDivEqUU_Helper<algo3b,sx,M1a,M2a>::call(D,A);
                    MultUM_Helper<algo4b,sx,sy,true,-1,RT,M1at,M2bt,M1bt>::
                        call(mone,Dt,Bt,Et);
                    LDivEqMU_Helper<algo3b,sy,sx,M1b,M2a>::call(E,C);
                    LDivEqUU_Helper<algo3b,sy,M1a,M2a>::call(F,C);
                }
#endif
            }
            else LDivEqUU_Helper<algo2,s,M1,M2>::call(m1,m2);
        }
    };

    // algo 87: Use temporary for m1/m2
    template <int s, class M1, class M2>
    struct LDivEqUU_Helper<87,s,M1,M2>
    {
        template <bool unknowndiag, int dummy>
        struct copyBack
        { // unknowndiag = false
            template <class M1c>
            static void call(const M1c& m1c, M1& m1)
            { NoAliasCopy(m1c,m1); }
        };
        template <int dummy>
        struct copyBack<true,dummy>
        {
            template <class M1c>
            static void call(const M1c& m1c, M1& m1)
            {
                if (m1.isunit()) NoAliasCopy(m1c.viewAsUnitDiag(),m1); 
                else NoAliasCopy(m1c,m1); 
            }
        };
        static void call(M1& m1, const M2& m2)
        {
#ifdef PRINTALGO_DivU
            const int N = s == UNKNOWN ? int(m1.size()) : s;
            std::cout<<"LDivEqUU algo 87: N,s = "<<N<<','<<s<<std::endl;
#endif
            typename M1::copy_type m1c = m1;
            NoAliasLDivEq(m1c,m2);
            const bool unknowndiag = M1::_unknowndiag && 
                (M2::_unit || M2::_unknowndiag);
            copyBack<unknowndiag,1>::call(m1c,m1);
        }
    };

    // algo -4: No branches or copies
    template <int s, class M1, class M2>
    struct LDivEqUU_Helper<-4,s,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            const bool upper2 = M2::_upper;
            const int algo = 
                ( s == 0 ) ? 0 :
                s == 1 ? 1 :
                //( s != UNKNOWN ) ? (
                    //upper2 ? (s <= 5 ? 16 : 17 ) :
                    //(s <= 5 ? 26 : 27 ) ) :
                upper2 ? 17 : 27;
            LDivEqUU_Helper<algo,s,M1,M2>::call(m1,m2);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int s, class M1, class M2>
    struct LDivEqUU_Helper<-3,s,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            // Possible algorithms are:
            //
            // Trivial and special for small TriMatrix sizes.
            //  0 = s == 0, so nothing to do
            //  1 = s == 1: reduces to scalar quotient
            //
            // UpperTri:
            // 11 = loop over n: MultUV
            // 12 = loop over m: MultMV
            // 13 = loop over k: Rank1 
            // 16 = unroll
            // 17 = split trimatrix into 3 submatrices
            //
            // LowerTri:
            // 21 = loop over n: MultUV
            // 22 = loop over m: MultMV
            // 23 = loop over k: Rank1
            // 26 = unroll
            // 27 = split trimatrix into 3 submatrices

            const bool upper2 = M2::_upper;
#if 0
            //const int algo = upper2 ? 17 : 27;
            const int algo = 36;
#else
#if TMV_OPT == 0 
            const bool rr = M1::_rowmajor && M2::_rowmajor;
            const bool rc = M1::_rowmajor && M2::_colmajor;
            const bool cx = M1::_colmajor;
            const int algo = 
                upper2 ? 
                ( rr ? 12 : rc ? 13 : cx ? 11 : 13 ) :
                ( rr ? 22 : rc ? 23 : cx ? 21 : 23 );
#else
            const int algo = 
                ( s == 0 ) ? 0 :
                s == 1 ? 201 :
                upper2 ? (
                    //s <= 5 ? 16 :
                    17 ) :
                ( // lowertri
                    //s <= 5 ? 26 :
                    27 );
#endif
#endif
#ifdef PRINTALGO_DivU
            const int N = s==UNKNOWN ? int(m1.size()) : s;
            std::cout<<"InlineLDivEqUU: \n";
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"N = "<<N<<std::endl;
            std::cout<<"s = "<<s<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            LDivEqUU_Helper<algo,s,M1,M2>::call(m1,m2);
        }
    };

    // algo 97: Conjugate
    template <int s, class M1, class M2>
    struct LDivEqUU_Helper<97,s,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        { 
            typedef typename M1::conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            LDivEqUU_Helper<-2,s,M1c,M2c>::call(m1c,m2c);
        }
    };

    // algo 98: call InstMultMM
    template <int s, class M1, class M2>
    struct LDivEqUU_Helper<98,s,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        { InstLDivEq(m1.xdView(),m2.xdView()); }
    };

    // algo -2: Check for inst
    template <int s, class M1, class M2>
    struct LDivEqUU_Helper<-2,s,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst = 
                (s == UNKNOWN || s > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T1>::isinst;
            const int algo = 
                ( s == 0 ) ? 0 :
                s == 1 ? 1 :
                M1::_conj ? 97 :
                inst ? 98 : 
                -3;
            LDivEqUU_Helper<algo,s,M1,M2>::call(m1,m2);
        }
    };

    // algo 99: Check for aliases
    template <int s, class M1, class M2>
    struct LDivEqUU_Helper<99,s,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            if (!SameStorage(m1,m2)) {
                // No aliasing 
                LDivEqUU_Helper<-2,s,M1,M2>::call(m1,m2);
            } else { 
                // Use temporary for m1/m2
                LDivEqUU_Helper<87,s,M1,M2>::call(m1,m2);
            }
        }
    };

    // algo -1: Check for aliases?
    template <int s, class M1, class M2>
    struct LDivEqUU_Helper<-1,s,M1,M2>
    {
        static void call(M1& m1, const M2& m2)
        {
            const bool checkalias =
                M1::_size == UNKNOWN && 
                M2::_size == UNKNOWN;
            const int algo = 
                ( s == 0 ) ? 0 :
                s == 1 ? 1 :
                checkalias ? 99 : 
                -2;
            LDivEqUU_Helper<algo,s,M1,M2>::call(m1,m2);
        }
    };

    template <int algo, class M1, class M2>
    static void DoLDivEq(
        BaseMatrix_Tri_Mutable<M1>& m1, const BaseMatrix_Tri<M2>& m2)
    {
        TMVStaticAssert((Sizes<M2::_size,M1::_size>::same));
        TMVAssert(m2.size() == m1.size());

        const int s = Sizes<M1::_size,M2::_size>::size;
        typedef typename M1::cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        M1v m1v = m1.cView();
        M2v m2v = m2.cView();
        if (m2.isSingular()) {
#ifdef TMV_NO_THROW
            std::cerr<<"Singular TriMatrix found in Division\n";
            exit(1);
#else
            throw SingularMatrix<M2>(m2.mat());
#endif
        }
        LDivEqUU_Helper<algo,s,M1v,M2v>::call(m1v,m2v);
    }
    template <class M1, class M2>
    static void LDivEq(
        BaseMatrix_Tri_Mutable<M1>& m1, const BaseMatrix_Tri<M2>& m2)
    { DoLDivEq<-1>(m1,m2); }
    template <class M1, class M2>
    static void NoAliasLDivEq(
        BaseMatrix_Tri_Mutable<M1>& m1, const BaseMatrix_Tri<M2>& m2)
    { DoLDivEq<-2>(m1,m2); }
    template <class M1, class M2>
    static void InlineLDivEq(
        BaseMatrix_Tri_Mutable<M1>& m1, const BaseMatrix_Tri<M2>& m2)
    { DoLDivEq<-3>(m1,m2); }
    template <class M1, class M2>
    static void AliasLDivEq(
        BaseMatrix_Tri_Mutable<M1>& m1, const BaseMatrix_Tri<M2>& m2)
    { DoLDivEq<99>(m1,m2); }

    //
    // m3 = m1 / m2
    //

    template <int algo, int ix, class T, class M1, class M2, class M3>
    struct LDivUU_Helper;

    // algo -2: NoAlias: Move along to LDivEq
    template <int ix, class T, class M1, class M2, class M3>
    struct LDivUU_Helper<-2,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            NoAliasMultXM<false>(x,m1,m3);
            NoAliasLDivEq(m3,m2);
        }
    };

    // algo 99: Check for aliases
    template <int ix, class T, class M1, class M2, class M3>
    struct LDivUU_Helper<99,ix,T,M1,M2,M3>
    {
        template <bool unknowndiag, class M3c>
        struct copyBack
        { // unknowndiag = false
            static void call(const M3c& m3c, M3& m3)
            { NoAliasCopy(m3c,m3); }
        };
        template <class M3c>
        struct copyBack<true,M3c>
        {
            static void call(const M3c& m3c, M3& m3)
            {
                if (m3.isunit()) {
                    typedef typename M3c::const_unitdiag_type M3cu;
                    M3cu m3cu = m3c.viewAsUnitDiag();
                    NoAliasCopy(m3cu,m3);
                } else {
                    NoAliasCopy(m3c,m3);
                }
            }
        };

        // The copyBack call is only needed if M3 is a TriMatrix, not if it 
        // is a regular Matrix.  So break out two options here.
        template <class M3x>
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, 
            BaseMatrix_Tri_Mutable<M3x>& m3)
        {
            if ( !SameStorage(m2,m3.mat()) ) {
                AliasMultXM<false>(x,m1,m3.mat());
                NoAliasLDivEq(m3.mat(),m2);
            } else {
                typedef typename M3::copy_type M3c;
                M3c m3c(m3.size());
                NoAliasMultXM<false>(x,m1,m3c);
                NoAliasLDivEq(m3c,m2);
                const bool unknowndiag = M3::_unknowndiag && ix == 1 && 
                    (M1::_unit || M1::_unknowndiag) &&
                    (M2::_unit || M2::_unknowndiag);
                copyBack<unknowndiag,M3c>::call(m3c,m3.mat());
            }
        }
        template <class M3x>
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2,
            BaseMatrix_Mutable<M3x>& m3)
        {
            if ( !SameStorage(m2,m3) ) {
                AliasMultXM<false>(x,m1,m3.mat());
                NoAliasLDivEq(m3.mat(),m2);
            } else {
                typedef typename M3::copy_type M3c;
                M3c m3c(m3.colsize(),m3.rowsize());
                NoAliasMultXM<false>(x,m1,m3c);
                NoAliasLDivEq(m3c,m2);
                NoAliasCopy(m3c,m3.mat());
            }
        }
    };

    // algo -1: Check for aliases?
    template <int ix, class T, class M1, class M2, class M3>
    struct LDivUU_Helper<-1,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const bool checkalias =
                M1::_size == UNKNOWN &&
                M2::_size == UNKNOWN &&
                M3::_size == UNKNOWN;
            const int algo =
                checkalias ? 99 :
                -2;
            LDivUU_Helper<algo,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

#if 0
    template <int algo, int ix, class T, class M1, class M2, class M3>
    static void DoLDiv(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Calc<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_size,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_size,M3::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
        TMVAssert(m1.size() == m3.colsize());
        TMVAssert(m1.size() == m3.rowsize());
        TMVAssert(m1.size() == m2.size());
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        M1v m1v = m1.cView();
        M2v m2v = m2.cView();
        M3v m3v = m3.cView();
        LDivUU_Helper<algo,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }
    template <int ix, class T, class M1, class M2, class M3>
    static void LDiv(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Calc<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Mutable<M3>& m3)
    { DoLDiv<-1>(x,m1,m2,m3); }
    template <int ix, class T, class M1, class M2, class M3>
    static void NoAliasLDiv(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Calc<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Mutable<M3>& m3)
    { DoLDiv<-2>(x,m1,m2,m3); }
    template <int ix, class T, class M1, class M2, class M3>
    static void AliasLDiv(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Calc<M1>& m1, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Mutable<M3>& m3)
    { DoLDiv<99>(x,m1,m2,m3); }

    //
    // m1 %= m2
    //

    template <class M1, class M2>
    static void RDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix_Tri<M2>& m2)
    { 
        typename M1::transpose_type m1t = m1.transpose();
        LDivEq(m1t,m2.transpose()); 
    }
    template <class M1, class M2>
    static void NoAliasRDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix_Tri<M2>& m2)
    { 
        typename M1::transpose_type m1t = m1.transpose();
        NoAliasLDivEq(m1t,m2.transpose()); 
    }
    template <class M1, class M2>
    static void AliasRDivEq(
        BaseMatrix_Mutable<M1>& m1, const BaseMatrix_Tri<M2>& m2)
    {
        typename M1::transpose_type m1t = m1.transpose();
        AliasLDivEq(m1t,m2.transpose()); 
    }

    //
    // m3 = m1 % m2
    //

    template <int ix, class T, class M1, class M2, class M3>
    static void RDiv(
        const Scaling<ix,T>& x,
        const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3)
    { 
        typename M3::transpose_type m3t = m3.transpose();
        LDiv(x,m1.transpose(),m2.transpose(),m3t); 
    }
    template <int ix, class T, class M1, class M2, class M3>
    static void NoAliasRDiv(
        const Scaling<ix,T>& x,
        const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3)
    { 
        typename M3::transpose_type m3t = m3.transpose();
        NoAliasLDiv(x,m1.transpose(),m2.transpose(),m3t); 
    }
    template <int ix, class T, class M1, class M2, class M3>
    static void AliasRDiv(
        const Scaling<ix,T>& x,
        const BaseMatrix_Calc<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3)
    {
        typename M3::transpose_type m3t = m3.transpose();
        AliasLDiv(x,m1.transpose(),m2.transpose(),m3t); 
    }
#endif

#undef TMV_Q1
#undef TMV_Q2

} // namespace tmv

#endif 
