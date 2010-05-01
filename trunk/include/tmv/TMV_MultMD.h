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
// Boston, MA  02320-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#ifndef TMV_MultMD_H
#define TMV_MultMD_H

#include "TMV_ElemMultVV.h"
#include "TMV_MultXV.h"
#include "TMV_Vector.h"
#include "TMV_SmallVector.h"
#include "TMV_BaseMatrix_Rec.h"
#include "TMV_Prefetch.h"
#include "TMV_DiagMatrix.h"
#include "TMV_CopyM.h"
#include "TMV_MultXM.h"

//#define PRINTALGO_MD

#ifdef PRINTALGO_MD
#include <iostream>
#endif

namespace tmv {

    // Defined below:
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void NoAliasMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void AliasMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, const BaseMatrix_Rec<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void NoAliasMultMM(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, const BaseMatrix_Rec<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void AliasMultMM(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, const BaseMatrix_Rec<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3);
    template <class M1, int ix, class T, class M2>
    inline void MultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2);
    template <class M1, int ix, class T, class M2>
    inline void NoAliasMultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2);
    template <class M1, int ix, class T, class M2>
    inline void AliasMultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2);


    // Defined in TMV_MultMD.cpp
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstMultMM(
        const T3 x,
        const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstDiagMatrixView<T2,UNKNOWN,C2>& m2, MatrixView<T3> m3);
    template <class T1, bool C1, class T2, bool C2, class T3>
    void InstAddMultMM(
        const T3 x,
        const ConstMatrixView<T1,UNKNOWN,UNKNOWN,C1>& m1,
        const ConstDiagMatrixView<T2,UNKNOWN,C2>& m2, MatrixView<T3> m3);

    //
    // Matrix * DiagMatrix
    //

    // Q3 is the crossover memory size to start using prefetch commands.
    // This is undoubtedly a function of the L1 (and L2?) cache size,
    // but 2KBytes is probably not too bad for most machines.
    // (That's an empirical value for my Intel Core 2 Duo.)
#define TMV_Q3 2048

    // ZeroIX controls whether ix = -1 should act like ix = 1 or ix = 0.
    // It doesn't really seem to matter much either way.
#define TMV_ZeroIX (ix == 0)
    //#define TMV_ZeroIX (ix != 1)

    template <int algo, int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultMD_Helper;

    // algo 0: cs or rs = 0, so nothing to do
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultMD_Helper<0,cs,rs,add,ix,T,M1,M2,M3>
    {
        static void call(const Scaling<ix,T>& , const M1& , const M2& , M3& ) 
        {} 
    };

    // algo 1: cs == 1, so simplifies to a MultDV function
    template <int rs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMD_Helper<1,1,rs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::const_row_type M1r;
            typedef typename M2::const_diag_type M2d;
            typedef typename M3::row_type M3r;
            M1r m10 = m1.get_row(0);
            M2d m2d = m2.diag();
            M3r m30 = m3.get_row(0);
            ElemMultVV_Helper<-1,rs,add,ix,T,M1r,M2d,M3r>::call(x,m10,m2d,m30);
        }
    };

    // algo 2: rs == 1, so simplifies to a MultXV function
    template <int cs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMD_Helper<2,cs,1,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::const_col_type M1c;
            typedef typename M3::col_type M3c;
            typedef typename Traits2<T,typename M2::value_type>::type PT2;
            M1c m10 = m1.get_col(0);
            PT2 xm2 = ZProd<false,false>::prod(x , m2.cref(0));
            M3c m30 = m3.get_col(0);
            MultXV_Helper<-1,cs,add,0,PT2,M1c,M3c>::call(xm2,m10,m30);
        }
    };

    // algo 201: same as 1, but use -2 algo
    template <int rs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMD_Helper<201,1,rs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::const_row_type M1r;
            typedef typename M2::const_diag_type M2d;
            typedef typename M3::row_type M3r;
            M1r m10 = m1.get_row(0);
            M2d m2d = m2.diag();
            M3r m30 = m3.get_row(0);
            ElemMultVV_Helper<-2,rs,add,ix,T,M1r,M2d,M3r>::call(x,m10,m2d,m30);
        }
    };

    // algo 202: same as 2, but use -2 algo
    template <int cs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMD_Helper<202,cs,1,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::const_col_type M1c;
            typedef typename M3::col_type M3c;
            typedef typename Traits2<T,typename M2::value_type>::type PT2;
            M1c m10 = m1.get_col(0);
            PT2 xm2 = ZProd<false,false>::prod(x , m2.cref(0));
            M3c m30 = m3.get_col(0);
            MultXV_Helper<-2,cs,add,0,PT2,M1c,M3c>::call(xm2,m10,m30);
        }
    };

    // algo 401: same as 1, but use -4 algo
    template <int rs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMD_Helper<401,1,rs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::const_row_type M1r;
            typedef typename M2::const_diag_type M2d;
            typedef typename M3::row_type M3r;
            M1r m10 = m1.get_row(0);
            M2d m2d = m2.diag();
            M3r m30 = m3.get_row(0);
            ElemMultVV_Helper<-4,rs,add,ix,T,M1r,M2d,M3r>::call(x,m10,m2d,m30);
        }
    };

    // algo 402: same as 2, but use -4 algo
    template <int cs, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultMD_Helper<402,cs,1,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::const_col_type M1c;
            typedef typename M3::col_type M3c;
            typedef typename Traits2<T,typename M2::value_type>::type PT2;
            M1c m10 = m1.get_col(0);
            PT2 xm2 = ZProd<false,false>::prod(x , m2.cref(0));
            M3c m30 = m3.get_col(0);
            MultXV_Helper<-4,cs,add,0,PT2,M1c,M3c>::call(xm2,m10,m30);
        }
    };

    // algo 11: The basic column major loop
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultMD_Helper<11,cs,rs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = (rs == UNKNOWN ? m3.rowsize() : rs);
#ifdef PRINTALGO_MD
            const int M = (cs == UNKNOWN ? m3.colsize() : cs);
            std::cout<<"MD algo 11: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif

            typedef typename M1::const_col_type M1c;
            typedef typename M2::value_type T2;
            typedef typename Traits2<T,T2>::type PT2;
            typedef typename M3::col_type M3c;

            for(int j=0;j<N;++j) {
                PT2 dj = x * m2.cref(j);
                M1c m1j = m1.get_col(j);
                M3c m3j = m3.get_col(j);
                MultXV_Helper<-4,cs,add,0,PT2,M1c,M3c>::call(dj,m1j,m3j);
            } 
        }
    };

    // algo 12: Column major loop with iterators
    template <int cs, int rs, bool add,
              int ix, class T, class M1, class M2, class M3>
    struct MultMD_Helper<12,cs,rs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = (cs == UNKNOWN ? m3.colsize() : cs);
            int N = (rs == UNKNOWN ? m3.rowsize() : rs);
#ifdef PRINTALGO_MD
            std::cout<<"MD algo 12: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            const bool c2 = M2::_conj;

            typedef typename M1::value_type T1;
            typedef typename M1::const_col_type M1c;
            typedef typename M1c::const_nonconj_type::const_iterator IT1;
            const int Astepj = m1.stepj();
            IT1 A = m1.get_col(0).nonConj().begin();

            typedef typename M2::const_diag_type M2d;
            typedef typename M2d::const_nonconj_type::const_iterator IT2;
            typedef typename M2::value_type T2;
            typedef typename Traits2<T,T2>::type PT2;
            IT2 D = m2.diag().nonConj().begin();
            PT2 dj;

            typedef typename M3::col_type M3c;
            typedef typename M3c::iterator IT3;
            const int Bstepj = m3.stepj();
            IT3 B = m3.get_col(0).begin();

            const bool dopref = M * sizeof(T1) >= TMV_Q3;

            Prefetch_Read(D.get());
            Prefetch_Read(A.get());
            Prefetch_Write(B.get());

            if (N) do {
                dj = ZProd<false,c2>::prod(x , *D++);
                MultXV_Helper<-4,UNKNOWN,add,0,PT2,M1c,M3c>::call2(
                    M,dj,A,B);
                A.shiftP(Astepj);
                B.shiftP(Bstepj);
                if (dopref) {
                    Prefetch_Read(A.get());
                    Prefetch_Write(B.get());
                }
            } while (--N);
        }
    };

    // algo 14: The basic row major loop
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultMD_Helper<14,cs,rs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int M = (cs == UNKNOWN ? m3.colsize() : cs);
#ifdef PRINTALGO_MD
            const int N = (rs == UNKNOWN ? m3.rowsize() : rs);
            std::cout<<"MD algo 14: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif

            typedef typename M1::const_row_type M1r;
            typedef typename M2::const_diag_type M2d;
            typedef typename M3::row_type M3r;
            M2d m2d = m2.diag();

            for (int i=0;i<M;++i) {
                M1r m1i = m1.get_row(i);
                M3r m3i = m3.get_row(i);
                ElemMultVV_Helper<-4,rs,add,ix,T,M1r,M2d,M3r>::call(
                    x,m1i,m2d,m3i);
            }
        }
    };

    // algo 15: Row major loop with iterators
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultMD_Helper<15,cs,rs,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            int M = (cs == UNKNOWN ? m3.colsize() : cs);
            const int N = (rs == UNKNOWN ? m3.rowsize() : rs);
#ifdef PRINTALGO_MD
            std::cout<<"MD algo 15: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::value_type T1;
            typedef typename M1::const_row_type M1r;
            typedef typename M1r::const_nonconj_type::const_iterator IT1;
            const int Astepi = m1.stepi();
            IT1 A = m1.get_row(0).nonConj().begin();

            typedef typename M2::const_diag_type M2d;
            typedef typename M2d::const_nonconj_type::const_iterator IT2;
            IT2 D = m2.diag().nonConj().begin();

            typedef typename M3::row_type M3r;
            typedef typename M3r::iterator IT3;
            const int Bstepi = m3.stepi();
            IT3 B = m3.get_row(0).begin();

            const bool dopref = N * sizeof(T1) >= TMV_Q3;

            Prefetch_MultiRead(D.get());
            Prefetch_Read(A.get());
            Prefetch_Write(B.get());

            if (M) do {
                ElemMultVV_Helper<-4,UNKNOWN,add,ix,T,M1r,M2d,M3r>::call2(
                    N,x,A,D,B);
                A.shiftP(Astepi);
                B.shiftP(Bstepi);
                if (dopref) {
                    Prefetch_Read(A.get());
                    Prefetch_Write(B.get());
                }
            } while (--M);
        }
    };

    // algo 82: copy x*m2
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultMD_Helper<82,cs,rs,add,ix,T,M1,M2,M3> 
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = (rs == UNKNOWN ? m3.rowsize() : rs);
#ifdef PRINTALGO_MD
            const int M = (cs == UNKNOWN ? m3.colsize() : cs);
            std::cout<<"MD algo 82: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            typedef typename M2::value_type T2;
            typedef typename Traits2<T,T2>::type PT2;
            typedef typename MCopyHelper<PT2,Diag,rs,rs,false,false>::type M2c;
            M2c m2c(N);
            typedef typename M2::const_diag_type M2d;
            typedef typename M2c::diag_type M2cd;
            typedef typename M2c::const_view_type M2cv;
            M2d m2d = m2.diag();
            M2cd m2cd = m2c.diag();
            M2cv m2cv = m2c.view();
            typedef typename Traits<T>::real_type RT;
            const Scaling<1,RT> one;
            MultXV_Helper<-2,rs,false,ix,T,M2d,M2cd>::call(x,m2d,m2cd);
            MultMD_Helper<-2,cs,rs,add,1,RT,M1,M2cv,M3>::call(one,m1,m2cv,m3);
        }
    };

    // algo -4: No branches or copies
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultMD_Helper<-4,cs,rs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const bool bothrm = M1::_rowmajor && M3::_rowmajor;
            const bool bothcm = M1::_colmajor && M3::_colmajor;
            const int algo = 
                ( rs == 0 || cs == 0 ) ? 0 :
                ( cs == 1 ) ? 401 :
                ( rs == 1 ) ? 402 :
                bothrm ? ( 
                    ( rs != UNKNOWN && rs <= 5 && cs > rs ) ? 11 : 
                    ( rs != UNKNOWN && rs <= 10 ) ? 14 : 15 ) :
                bothcm ? ( 
                    ( cs != UNKNOWN && cs <= 10 ) ? 11 : 12 ) :
                ( cs != UNKNOWN && cs <= 10 ) ? 11 : 12;
            MultMD_Helper<algo,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultMD_Helper<-3,cs,rs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            // Possible algorithms to choose from:
            //  0 = cs or rs == 0, so nothing to do
            //  1 = cs == 1: reduces to trivial MultDV function
            //  2 = rs == 1: reduces to trivial MultXV function
            //
            // 11 = column major loop: MultXV on each column
            // 14 = row major loop: MultDV on each row
            //
            // 82 = copy x*m2

            const bool bothrm = M1::_rowmajor && M3::_rowmajor;
            const bool bothcm = M1::_colmajor && M3::_colmajor;
#if TMV_OPT >= 1
            const bool docopy = TMV_ZeroIX || M2::_diagstep != 1;
#else
            const bool docopy = false;
#endif

            const int algo = 
                ( rs == 0 || cs == 0 ) ? 0 :
                ( cs == 1 ) ? 401 :
                ( rs == 1 ) ? 402 :
                bothrm ? ( 
                    ( rs != UNKNOWN && rs <= 5 && cs > rs ) ? 11 : 
                    docopy ? 82 : 
                    ( rs != UNKNOWN && rs <= 10 ) ? 14 : 15 ) :
                bothcm ? ( 
                    ( cs != UNKNOWN && cs <= 10 ) ? 11 : 12 ) :
                ( cs != UNKNOWN && cs <= 10 ) ? 11 : 12;
#ifdef PRINTALGO_MD
            std::cout<<"InlineMultMD: \n";
            std::cout<<"x = "<<ix<<"  "<<T(x)<<"  add = "<<add<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"cs,rs,algo = "<<cs<<"  "<<rs<<"  "<<algo<<std::endl;
#endif
            MultMD_Helper<algo,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo 97: Conjugate
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultMD_Helper<97,cs,rs,add,ix,T,M1,M2,M3>
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
            MultMD_Helper<-2,cs,rs,add,ix,T,M1c,M2c,M3c>::call(
                TMV_CONJ(x),m1c,m2c,m3c);
        }
    };

    // algo 98: call InstMultMM
    template <int cs, int rs, int ix, class T, class M1, class M2, class M3>
    struct MultMD_Helper<98,cs,rs,false,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typename M3::value_type xx(x);
            InstMultMM(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };
    template <int cs, int rs, int ix, class T, class M1, class M2, class M3>
    struct MultMD_Helper<98,cs,rs,true,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typename M3::value_type xx(x);
            InstAddMultMM(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultMD_Helper<-2,cs,rs,add,ix,T,M1,M2,M3>
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
                ( rs == 0 || cs == 0 ) ? 0 :
                ( cs == 1 ) ? 201 :
                ( rs == 1 ) ? 202 :
                M3::_conj ? 97 :
                inst ? 98 : 
                -3;
            MultMD_Helper<algo,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo 99: Check for aliases
    template <int cs, int rs, int ix, class T, class M1, class M2, class M3>
    struct MultMD_Helper<99,cs,rs,true,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            if ( ( !SameStorage(m1,m3) ||
                   ExactSameStorage(m1,m3) ) &&
                 !SameStorage(m2,m3) ) {
                // No aliasing (or no clobbering)
                MultMD_Helper<-2,cs,rs,true,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else if (SameStorage(m2,m3)) {
                // Use temporary for m2
                MultMD_Helper<82,cs,rs,true,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else { 
                // SameStorage(m1,m3)
                // Need a temporary
                NoAliasMultMM<true>(x,m1.copy(),m2,m3);
            }
        }
    };
    template <int cs, int rs, int ix, class T, class M1, class M2, class M3>
    struct MultMD_Helper<99,cs,rs,false,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            if ( ( !SameStorage(m1,m3) ||
                   ExactSameStorage(m1,m3) ) &&
                 !SameStorage(m2,m3) ) {
                // No aliasing (or no clobbering)
                MultMD_Helper<-2,cs,rs,false,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else if (SameStorage(m2,m3)) {
                // Use temporary for m2
                MultMD_Helper<82,cs,rs,false,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else { 
                // SameStorage(m1,m3)
                // Let Copy handle the aliasing
                AliasCopy(m1,m3);
                NoAliasMultMM<false>(x,m3,m2,m3);
            }
        }
    };

    // algo -1: Check for aliases?
    template <int cs, int rs, bool add, 
              int ix, class T, class M1, class M2, class M3>
    struct MultMD_Helper<-1,cs,rs,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const bool checkalias =
                M1::_colsize == UNKNOWN && M1::_rowsize == UNKNOWN &&
                M2::_size == UNKNOWN && 
                M3::_colsize == UNKNOWN && M3::_rowsize == UNKNOWN;
            const int algo = 
                ( rs == 0 || cs == 0 ) ? 0 :
                ( cs == 1 ) ? 1 :
                ( rs == 1 ) ? 2 :
                checkalias ? 99 : 
                -2;
            MultMD_Helper<algo,cs,rs,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    template <int algo, bool add, int ix, class T, class M1, class M2, class M3>
    inline void DoMultMD(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M3::_colsize,M1::_colsize>::same));
        TMVStaticAssert((Sizes<M3::_rowsize,M1::_rowsize>::same));
        TMVStaticAssert((Sizes<M3::_rowsize,M2::_size>::same));
        TMVAssert(m3.colsize() == m1.colsize());
        TMVAssert(m3.rowsize() == m1.rowsize());
        TMVAssert(m3.rowsize() == m2.size());
        const int cs = Sizes<M3::_colsize,M1::_colsize>::size;
        const int rs1 = Sizes<M3::_rowsize,M1::_rowsize>::size;
        const int rs = Sizes<rs1,M2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        M1v m1v = m1.cView();
        M2v m2v = m2.cView();
        M3v m3v = m3.cView();
        MultMD_Helper<algo,cs,rs,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    { DoMultMD<-1,add>(x,m1,m2,m3); }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void NoAliasMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    { DoMultMD<-2,add>(x,m1,m2,m3); }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    { DoMultMD<-3,add>(x,m1,m2,m3); }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void AliasMultMM(
        const Scaling<ix,T>& x, 
        const BaseMatrix_Rec<M1>& m1, const BaseMatrix_Diag<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    { DoMultMD<99,add>(x,m1,m2,m3); }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void MultMM(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, const BaseMatrix_Rec<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        typename M3::transpose_type m3t = m3.transpose();
        MultMM<add>(x,m2.transpose(),m1.transpose(),m3t);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void NoAliasMultMM(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, const BaseMatrix_Rec<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        typename M3::transpose_type m3t = m3.transpose();
        NoAliasMultMM<add>(x,m2.transpose(),m1.transpose(),m3t);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void InlineMultMM(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, const BaseMatrix_Rec<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        typename M3::transpose_type m3t = m3.transpose();
        InlineMultMM<add>(x,m2.transpose(),m1.transpose(),m3t);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    inline void AliasMultMM(
        const Scaling<ix,T>& x,
        const BaseMatrix_Diag<M1>& m1, const BaseMatrix_Rec<M2>& m2,
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        typename M3::transpose_type m3t = m3.transpose();
        AliasMultMM<add>(x,m2.transpose(),m1.transpose(),m3t);
    }

    template <class M1, int ix, class T, class M2>
    inline void MultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    { MultMM<false>(x,m1.mat(),m2.mat(),m1.mat()); }

    template <class M1, int ix, class T, class M2>
    inline void NoAliasMultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    { NoAliasMultMM<false>(x,m1.mat(),m2.mat(),m1.mat()); }

    template <class M1, int ix, class T, class M2>
    inline void AliasMultEqMM(
        BaseMatrix_Rec_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    { AliasMultMM<false>(x,m1.mat(),m2.mat(),m1.mat()); }

} // namespace tmv

#undef TMV_ZeroIX

#endif 
