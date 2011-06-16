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


#ifndef TMV_MultUD_H
#define TMV_MultUD_H

#include "TMV_BaseMatrix_Tri.h"
#include "TMV_BaseMatrix_Diag.h"
#include "TMV_MultXV.h"
#include "TMV_ElemMultVV.h"
#include "TMV_CopyU.h"

#ifdef PRINTALGO_UD
#include <iostream>
#endif

namespace tmv {

    // Defined in TMV_MultUD.cpp
    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1,
        const ConstDiagMatrixView<T2,C2>& m2, UpperTriMatrixView<T3> m3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1,
        const ConstDiagMatrixView<T2,C2>& m2, UpperTriMatrixView<T3> m3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstMultMM(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m1,
        const ConstDiagMatrixView<T2,C2>& m2, LowerTriMatrixView<T3> m3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMultMM(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m1,
        const ConstDiagMatrixView<T2,C2>& m2, LowerTriMatrixView<T3> m3);

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasMultMM(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1,
        const ConstDiagMatrixView<T2,C2>& m2, UpperTriMatrixView<T3> m3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMultMM(
        const T3 x, const ConstUpperTriMatrixView<T1,C1>& m1,
        const ConstDiagMatrixView<T2,C2>& m2, UpperTriMatrixView<T3> m3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasMultMM(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m1,
        const ConstDiagMatrixView<T2,C2>& m2, LowerTriMatrixView<T3> m3);
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMultMM(
        const T3 x, const ConstLowerTriMatrixView<T1,C1>& m1,
        const ConstDiagMatrixView<T2,C2>& m2, LowerTriMatrixView<T3> m3);

    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3);

    //
    // Matrix * DiagMatrix
    //

    // The maximum nops to unroll.
#if TMV_OPT >= 3
#define TMV_UD_UNROLL 200 
#elif TMV_OPT >= 2
#define TMV_UD_UNROLL 25
#elif TMV_OPT >= 1
#define TMV_UD_UNROLL 9
#else
#define TMV_UD_UNROLL 0
#endif

    // ZeroIX controls whether ix = -1 should act like ix = 1 or ix = 0.
    // It doesn't really seem to matter much either way.
#define TMV_ZeroIX (ix == 0)
    //#define TMV_ZeroIX (ix != 1)

    template <int algo, int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper;

    // algo 0: s = 0, so nothing to do
    template <bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<0,0,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& , const M1& , const M2& , M3& ) 
        {} 
    };

    // algo 1: s == 1, so simplifies to a scalar product
    template <bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<1,1,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_UD
            std::cout<<"UD algo 1: N,s,x = "<<1<<','<<1<<','<<T(x)<<std::endl;
#endif
            Maybe<add>::add(
                m3.ref(0,0),
                ZProd<false,false>::prod(
                    x,ZProd<false,false>::prod(m1.cref(0,0) , m2.cref(0,0)) ));
        }
    };

    // algo 11: UpperTri: Basic column major loop 
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<11,s,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            int N = (s == UNKNOWN ? m3.size() : s);
#ifdef PRINTALGO_UD
            std::cout<<"UD algo 11: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            const bool unit = m1.isunit();
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::value_type T2;
            typedef typename Traits2<T,T2>::type PT2;
            typedef typename M2::const_diag_type M2d;
            typedef typename M3::col_sub_type M3c;
            typedef typename M3::diag_type M3d;

            for(int j=0,jj=(unit?0:1);j<N;++j,++jj) {
                PT2 dj = x * m2.cref(j);
                M1c m1j = m1.get_col(j,0,jj);
                M3c m3j = m3.get_col(j,0,jj);
                MultXV_Helper<-4,UNKNOWN,add,0,PT2,M1c,M3c>::call(
                    dj,m1j,m3j);
            }
            if (unit) {
                M2d m2d = m2.diag();
                M3d m3d = m3.diag();
                MultXV_Helper<-4,s,add,ix,T,M2d,M3d>::call(x,m2d,m3d);
            }
        }
    };

    // algo 12: UpperTri: Column major loop with iterators
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<12,s,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            int N = (s == UNKNOWN ? m3.size() : s);
#ifdef PRINTALGO_UD
            std::cout<<"UD algo 12: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            const bool unit = m1.isunit();
            const bool c2 = M2::_conj;

            typedef typename M1::const_col_sub_type M1c;
            typedef typename M1c::const_nonconj_type::const_iterator IT1;
            const int Astepj = m1.stepj();
            IT1 A = m1.get_col(0,0,1).begin().nonConj();

            typedef typename M2::const_diag_type M2d;
            typedef typename M2d::const_nonconj_type::const_iterator IT2;
            typedef typename M2::value_type T2;
            typedef typename Traits2<T,T2>::type PT2;
            IT2 D = m2.diag().begin().nonConj();
            PT2 dj;

            typedef typename M3::col_sub_type M3c;
            typedef typename M3c::iterator IT3;
            const int Bstepj = m3.stepj();
            IT3 B = m3.get_col(0,0,1).begin();

            int len = unit?0:1;
            if (N) do {
                dj = ZProd<false,c2>::prod(x , *D++);
                MultXV_Helper<-4,UNKNOWN,add,0,PT2,M1c,M3c>::call2(
                    len,dj,A,B);
                if (unit) Maybe<add>::add(B[len], dj);
                ++len;
                A.shiftP(Astepj);
                B.shiftP(Bstepj);
            } while (--N);
        }
    };

    // algo 14: UpperTri: Basic row major loop
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<14,s,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            int N = (s == UNKNOWN ? m3.size() : s);
#ifdef PRINTALGO_UD
            std::cout<<"UD algo 14: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            const bool unit = m1.isunit();
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M2::const_diag_type M2d;
            typedef typename M2d::const_subvector_type M2ds;
            typedef typename M3::row_sub_type M3r;
            typedef typename M3::diag_type M3d;
            M2d m2d = m2.diag();

            for(int i=0,ii=(unit?1:0);i<N;++i,++ii) {
                M1r m1i = m1.row(i,ii,N);
                M2ds m2di = m2d.cSubVector(ii,N);
                M3r m3i = m3.row(i,ii,N);
                ElemMultVV_Helper<-4,UNKNOWN,add,ix,T,M1r,M2ds,M3r>::call(
                    x,m1i,m2di,m3i);
            } 
            if (unit) {
                M3d m3d = m3.diag();
                MultXV_Helper<-4,s,add,ix,T,M2d,M3d>::call(x,m2d,m3d);
            }
        }
    };

    // algo 15: UpperTri: Row major loop with iterators
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<15,s,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            int N = (s == UNKNOWN ? m3.size() : s);
#ifdef PRINTALGO_UD
            std::cout<<"UD algo 15: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            const bool unit = m1.isunit();
            const bool c2 = M2::_conj;

            typedef typename M1::const_row_sub_type M1r;
            typedef typename M1r::const_nonconj_type::const_iterator IT1;
            const int Adiagstep = m1.diagstep();
            IT1 A = m1.get_row(0,unit?1:0,N).begin().nonConj();

            typedef typename M2::const_diag_type M2d;
            typedef typename M2d::const_nonconj_type::const_iterator IT2;
            IT2 D = m2.diag().begin().nonConj();

            typedef typename M3::row_sub_type M3r;
            typedef typename M3r::iterator IT3;
            const int Bdiagstep = m3.diagstep();
            IT3 B = m3.get_row(0,unit?1:0,N).begin();

            int len = unit ? N-1 : N;
            if (unit) ++D;
            if (N) do {
                if (unit) Maybe<add>::add(
                    *(B-1), ZProd<false,c2>::prod(x , *(D-1)));
                ElemMultVV_Helper<-4,UNKNOWN,add,ix,T,M1r,M2d,M3r>::call2(
                    len--,x,A,D++,B);
                A.shiftP(Adiagstep);
                B.shiftP(Bdiagstep);
            } while (--N);
        }
    };

    // algo 16: UpperTri: Unroll by columns
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<16,s,add,ix,T,M1,M2,M3>
    {
        template <int I, int N>
        struct Unroller
        {
            static inline void unroll(
                const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
            {
                Unroller<I,N/2>::unroll(x,m1,m2,m3);
                Unroller<I+N/2,N-N/2>::unroll(x,m1,m2,m3);
            }
        };
        template <int I>
        struct Unroller<I,1>
        {
            static inline void unroll(
                const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
            {
                typedef typename M2::value_type T2;
                typedef typename Traits2<T,T2>::type PT2;
                typedef typename M1::const_col_sub_type M1c;
                typedef typename M3::col_sub_type M3c;
                const bool u1 = M1::_unit;

                PT2 xd = x * m2.cref(I,I);
                M1c m1c = m1.get_col(I,0,I);
                M3c m3c = m3.get_col(I,0,I);
                MultXV_Helper<-4,I,add,0,PT2,M1c,M3c>::call(xd,m1c,m3c);
                Maybe<add>::add(
                    m3.ref(I,I), Maybe<!u1>::prod(m1.cref(I,I),xd));
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
#ifdef PRINTALGO_UD
            std::cout<<"UD algo 16: N,s,x = "<<2<<','<<2<<','<<T(x)<<std::endl;
#endif
            Unroller<0,s>::unroll(x,m1,m2,m3); 
        }
    };

    // algo 19: UpperTri: m1 is unit diag, do calculation on offdiag part.
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<19,s,add,ix,T,M1,M2,M3> // s known
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = (s == UNKNOWN ? m3.size() : s);
            TMVStaticAssert(M1::_unit);
#ifdef PRINTALGO_UD
            std::cout<<"UD algo 19: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_offdiag_type M1o;
            typedef typename M2::const_subdiagmatrix_type M2s;
            typedef typename M2::const_diag_type M2d;
            typedef typename M3::offdiag_type M3o;
            typedef typename M3::diag_type M3d;
            const int sm1 = IntTraits2<s,-1>::sum;
            if (N > 1) {
                M1o m1o = m1.offDiag();
                M2s m2s = m2.subDiagMatrix(1,N);
                M3o m3o = m3.offDiag();
                MultUD_Helper<-2,sm1,add,ix,T,M1o,M2s,M3o>::call(x,m1o,m2s,m3o);
            }
            M2d m2d = m2.diag();
            M3d m3d = m3.diag();
            MultXV_Helper<-2,s,add,ix,T,M2d,M3d>::call(x,m2d,m3d);
        }
    };

    // algo 21: LowerTri: Basic column major loop
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<21,s,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            int N = (s == UNKNOWN ? m3.size() : s);
#ifdef PRINTALGO_UD
            std::cout<<"UD algo 21: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            const bool unit = m1.isunit();
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::value_type T2;
            typedef typename Traits2<T,T2>::type PT2;
            typedef typename M2::const_diag_type M2d;
            typedef typename M3::col_sub_type M3c;
            typedef typename M3::diag_type M3d;

            for(int j=0,jj=(unit?1:0);j<N;++j,++jj) {
                PT2 dj = x * m2.cref(j);
                M1c m1j = m1.get_col(j,jj,N);
                M3c m3j = m3.get_col(j,jj,N);
                MultXV_Helper<-4,UNKNOWN,add,0,PT2,M1c,M3c>::call(
                    dj,m1j,m3j);
            }
            if (unit) {
                M2d m2d = m2.diag();
                M3d m3d = m3.diag();
                MultXV_Helper<-4,s,add,ix,T,M2d,M3d>::call(x,m2d,m3d);
            }
        }
    };

    // algo 22: LowerTri: Column major with iterators
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<22,s,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            int N = (s == UNKNOWN ? m3.size() : s);
#ifdef PRINTALGO_UD
            std::cout<<"UD algo 22: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            const bool unit = m1.isunit();
            const bool c2 = M2::_conj;

            typedef typename M1::const_col_sub_type M1c;
            typedef typename M1c::const_nonconj_type::const_iterator IT1;
            const int Adiagstep = m1.diagstep();
            IT1 A = m1.get_col(0,unit?1:0,N).begin().nonConj();

            typedef typename M2::const_diag_type M2d;
            typedef typename M2d::const_nonconj_type::const_iterator IT2;
            typedef typename M2::value_type T2;
            typedef typename Traits2<T,T2>::type PT2;
            IT2 D = m2.diag().begin().nonConj();
            PT2 dj;

            typedef typename M3::col_sub_type M3c;
            typedef typename M3c::iterator IT3;
            const int Bdiagstep = m3.diagstep();
            IT3 B = m3.get_col(0,unit?1:0,N).begin();

            int len = unit ? N-1 : N;
            if (N) do {
                dj = ZProd<false,c2>::prod(x , *D++);
                if (unit) Maybe<add>::add(*(B-1),dj);
                MultXV_Helper<-4,UNKNOWN,add,0,PT2,M1c,M3c>::call2(
                    len--,dj,A,B);
                A.shiftP(Adiagstep);
                B.shiftP(Bdiagstep);
            } while (--N);
        }
    };

    // algo 24: LowerTri: Basic row major loop
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<24,s,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            int N = (s == UNKNOWN ? m3.size() : s);
#ifdef PRINTALGO_UD
            std::cout<<"UD algo 24: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            const bool unit = m1.isunit();
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M2::const_diag_type M2d;
            typedef typename M2d::const_subvector_type M2ds;
            typedef typename M3::row_sub_type M3r;
            typedef typename M3::diag_type M3d;
            M2d m2d = m2.diag();

            for(int i=0,ii=(unit?0:1);i<N;++i,++ii) {
                M1r m1i = m1.row(i,0,ii);
                M2ds m2di = m2d.cSubVector(0,ii);
                M3r m3i = m3.row(i,0,ii);
                ElemMultVV_Helper<-4,UNKNOWN,add,ix,T,M1r,M2ds,M3r>::call(
                    x,m1i,m2di,m3i);
            } 
            if (unit) {
                M3d m3d = m3.diag();
                MultXV_Helper<-4,s,add,ix,T,M2d,M3d>::call(x,m2d,m3d);
            }
        }
    };

    // algo 25: LowerTri: Row major loop with iterators
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<25,s,add,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            int N = (s == UNKNOWN ? m3.size() : s);
#ifdef PRINTALGO_UD
            std::cout<<"UD algo 25: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            const bool unit = m1.isunit();
            const bool c2 = M2::_conj;

            typedef typename M1::const_row_sub_type M1r;
            typedef typename M1r::const_nonconj_type::const_iterator IT1;
            const int Astepi = m1.stepi();
            IT1 A = m1.get_row(0,0,unit?0:1).begin().nonConj();

            typedef typename M2::const_diag_type M2d;
            typedef typename M2d::const_nonconj_type::const_iterator IT2;
            const IT2 D = m2.diag().begin().nonConj();

            typedef typename M3::row_sub_type M3r;
            typedef typename M3r::iterator IT3;
            const int Bstepi = m3.stepi();
            IT3 B = m3.get_row(0,0,unit?0:1).begin();

            int len = unit?0:1;
            if (N) do {
                ElemMultVV_Helper<-4,UNKNOWN,add,ix,T,M1r,M2d,M3r>::call2(
                    len,x,A,D,B);
                if (unit) Maybe<add>::add(
                    B[len], ZProd<false,c2>::prod(x , D[len]));
                ++len;
                A.shiftP(Astepi);
                B.shiftP(Bstepi);
            } while (--N);
        }
    };

    // algo 26: LowerTri: Unroll by columns
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<26,s,add,ix,T,M1,M2,M3>
    {
        // [ C00  0  ] = [ A00  0  ] [ B00  0  ]
        // [ C10 C11 ]   [ A10 A11 ] [  0  B11 ]
        template <int I, int N>
        struct Unroller
        {
            static inline void unroll(
                const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
            {
                Unroller<I,N/2>::unroll(x,m1,m2,m3);
                Unroller<I+N/2,N-N/2>::unroll(x,m1,m2,m3);
            }
        };
        template <int I>
        struct Unroller<I,1>
        {
            static inline void unroll(
                const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
            {
                typedef typename M2::value_type T2;
                typedef typename Traits2<T,T2>::type PT2;
                typedef typename M1::const_col_sub_type M1c;
                typedef typename M3::col_sub_type M3c;
                const bool u1 = M1::_unit;

                PT2 xd = x * m2.cref(I,I);
                M1c m1c = m1.get_col(I,I+1,s);
                M3c m3c = m3.get_col(I,I+1,s);
                MultXV_Helper<-4,s-I-1,add,0,PT2,M1c,M3c>::call(xd,m1c,m3c);
                Maybe<add>::add(
                    m3.ref(I,I), Maybe<!u1>::prod(m1.cref(I,I),xd));
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
#ifdef PRINTALGO_UD
            std::cout<<"UD algo 26: N,s,x = "<<2<<','<<2<<','<<T(x)<<std::endl;
#endif
            Unroller<0,s>::unroll(x,m1,m2,m3); 
        }
    };

    // algo 29: LowerTri: m1 is unit diag, do calculation on offdiag part.
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<29,s,add,ix,T,M1,M2,M3> // s known
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int N = (s == UNKNOWN ? m3.size() : s);
            TMVStaticAssert(M1::_unit);
#ifdef PRINTALGO_UD
            std::cout<<"UD algo 29: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            typedef typename M1::const_offdiag_type M1o;
            typedef typename M2::const_subdiagmatrix_type M2s;
            typedef typename M2::const_diag_type M2d;
            typedef typename M3::offdiag_type M3o;
            typedef typename M3::diag_type M3d;
            const int sm1 = IntTraits2<s,-1>::sum;
            if (N > 1) {
                M1o m1o = m1.offDiag();
                M2s m2s = m2.subDiagMatrix(0,N-1);
                M3o m3o = m3.offDiag();
                MultUD_Helper<-2,sm1,add,ix,T,M1o,M2s,M3o>::call(x,m1o,m2s,m3o);
            }
            M2d m2d = m2.diag();
            M3d m3d = m3.diag();
            MultXV_Helper<-2,s,add,ix,T,M2d,M3d>::call(x,m2d,m3d);
        }
    };

    // algo 82: Copy x*m2
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<82,s,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_UD
            const int N = (s == UNKNOWN ? m3.size() : s);
            std::cout<<"UD algo 82: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            typedef typename M3::real_type RT;
            const Scaling<1,RT> one;
            NoAliasMultMM<add>(one,m1,(x*m2).calc(),m3);
        }
    };

    // algo 88: m3 = m1, then MultEq
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<88,s,add,ix,T,M1,M2,M3>
    {
        static inline void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_UD
            const int N = (s == UNKNOWN ? m3.size() : s);
            std::cout<<"UD algo 88: N,s,x = "<<N<<','<<s<<','<<T(x)<<std::endl;
#endif
            CopyU_Helper<-2,s,M1,M3>::call(m1,m3);
            typedef typename M3::const_view_type M3cv;
            M3cv m3cv = m3.view();
            MultUD_Helper<-2,s,add,ix,T,M3cv,M2,M3>::call(x,m3cv,m2,m3);
        }
    };

    // algo 90: call inst
    template <int s, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<90,s,false,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstMultMM(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };
    template <int s, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<90,s,true,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAddMultMM(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };

    // algo 91: call inst alias
    template <int s, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<91,s,false,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasMultMM(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };
    template <int s, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<91,s,true,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasAddMultMM(xx,m1.xView(),m2.xView(),m3.xView());
        }
    };

    // algo 97: Conjugate
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<97,s,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            typedef typename M3::conjugate_type M3c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            M3c m3c = m3.conjugate();
            MultUD_Helper<-2,s,add,ix,T,M1c,M2c,M3c>::call(
                TMV_CONJ(x),m1c,m2c,m3c);
        }
    };

    // algo 197: Conjugate
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<197,s,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            typedef typename M3::conjugate_type M3c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            M3c m3c = m3.conjugate();
            MultUD_Helper<99,s,add,ix,T,M1c,M2c,M3c>::call(
                TMV_CONJ(x),m1c,m2c,m3c);
        }
    };

    // algo 98: Inline check for aliases
    template <int s, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<98,s,true,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            if ( ( !SameStorage(m1,m3) ||
                   ExactSameStorage(m1,m3) ||
                   OppositeStorage(m1,m3) ) &&
                 !SameStorage(m2,m3) ) {
                // No aliasing (or no clobbering)
                MultUD_Helper<-2,s,true,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else if (SameStorage(m2,m3)) {
                // Use temporary for m2
                MultUD_Helper<82,s,true,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else {
                // SameStorage(m1,m3)
                // Need a temporary
                NoAliasMultMM<true>(x,m1.copy(),m2,m3);
            }
        }
    };
    template <int s, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<98,s,false,ix,T,M1,M2,M3>
    {
        static void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            if ( ( !SameStorage(m1,m3) ||
                   ExactSameStorage(m1,m3) ||
                   OppositeStorage(m1,m3) ) &&
                 !SameStorage(m2,m3) ) {
                // No aliasing (or no clobbering)
                MultUD_Helper<-2,s,false,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else if (SameStorage(m2,m3)) {
                // Use temporary for m2
                MultUD_Helper<82,s,false,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
            } else {
                // SameStorage(m1,m3)
                // Let Copy handle the aliasing
                AliasCopy(m1,m3);
                NoAliasMultMM<false>(x,m3,m2,m3);
            }
        }
    };

    // algo 99: Check for aliases
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<99,s,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            const bool inst = 
                (s == UNKNOWN || s > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const int algo = 
                ( s == 0 ) ? 0 :
                ( s == 1 ) ? 1 :
                M3::_conj ? 197 :
                inst ? 91 : 
                98;
            MultUD_Helper<algo,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo -4: No branches or copies
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<-4,s,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            const bool bothrm = M1::_rowmajor && M3::_rowmajor;
            const bool bothcm = M1::_colmajor && M3::_colmajor;
            const int s2 = s > 20 ? UNKNOWN : s;
            const int s2p1 = IntTraits<s2>::Sp1;
            // nops = n(n+1)/2
            const int nops = IntTraits2<s2,s2p1>::safeprod / 2;
            const bool unroll = 
                s > 10 ? false :
                s == UNKNOWN ? false :
                nops <= TMV_UD_UNROLL;
            const bool docopy =
                TMV_OPT == 0 ? false :
                TMV_ZeroIX || M2::_diagstep != 1;
            const int algo = 
                ( s == 0 ) ? 0 :
                ( s == 1 ) ? 1 :
                !(Traits2<T1,T2>::samebase && Traits2<T1,T3>::samebase) ?
                ( M3::_upper ? bothrm ? 14 : 11 : bothrm ? 24 : 21 ) :
                M3::_upper ? ( // UpperTri
                    unroll ? 16 :
                    bothrm ? ( ( s != UNKNOWN && s <= 10 ) ? 14 : 15 ) :
                    bothcm ? ( ( s != UNKNOWN && s <= 10 ) ? 11 : 12 ) :
                    ( s != UNKNOWN && s <= 10 ) ? 11 : 12 ) :
                ( // LowerTri
                    unroll ? 26 :
                    bothrm ? ( ( s != UNKNOWN && s <= 10 ) ? 24 : 25 ) :
                    bothcm ? ( ( s != UNKNOWN && s <= 10 ) ? 21 : 22 ) :
                    ( M1::_rowmajor && !docopy ) ? ( 
                        ( s != UNKNOWN && s <= 10 ) ? 24 : 25 ) :
                    ( s != UNKNOWN && s <= 10 ) ? 21 : 22 ); 
            MultUD_Helper<algo,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<-3,s,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            // Possible algorithms to choose from:
            //
            // Trivial:
            //  0 = s == 0, so nothing to do
            //  1 = s == 1: reduces to trivial scalar product
            //
            // UpperTri:
            // 11 = column major - basic loop
            // 12 = column major - use iterators
            // 14 = row major - basic loop
            // 15 = row major - use iterators
            // 16 = unroll
            // 19 = m1 is unitdiag
            //
            // LowerTri:
            // 21 = column major - basic loop
            // 22 = column major - use iterators
            // 24 = row major - basic loop
            // 25 = row major - use iterators
            // 26 = unroll
            // 29 = m1 is unitdiag
            // 
            // 82 = copy x*m2
            // 88 = copy m3 = m1, to make a MultEq op

            const bool bothrm = M1::_rowmajor && M3::_rowmajor;
            const bool bothcm = M1::_colmajor && M3::_colmajor;
            const bool docopy = 
                TMV_OPT == 0 ? false :
                TMV_ZeroIX || M2::_diagstep != 1;
            const int s2 = s > 20 ? UNKNOWN : s;
            const int s2p1 = IntTraits<s2>::Sp1;
            // nops = n(n+1)/2
            const int nops = IntTraits2<s2,s2p1>::safeprod / 2;
            const bool unroll = 
                s > 10 ? false :
                s == UNKNOWN ? false :
                nops <= TMV_UD_UNROLL;
            const int algo = 
                ( s == 0 ) ? 0 :
                ( s == 1 ) ? 1 :
                !(Traits2<T1,T2>::samebase && Traits2<T1,T3>::samebase) ?
                ( M3::_upper ? bothrm ? 14 : 11 : bothrm ? 24 : 21 ) :
                M3::_upper ? ( // UpperTri
                    unroll ? 16 :
                    bothrm ? (
                        docopy ? 82 :
                        ( s != UNKNOWN && s <= 10 ) ? 14 : 15 ) :
                    bothcm ? (
                        ( s != UNKNOWN && s <= 10 ) ? 11 : 12 ) :
                    ( s != UNKNOWN && s <= 10 ) ? 11 : 12 ) :
                ( // LowerTri
                    unroll ? 26 :
                    bothrm ? (
                        docopy ? 82 :
                        ( s != UNKNOWN && s <= 10 ) ? 24 : 25 ) :
                    bothcm ? (
                        ( s != UNKNOWN && s <= 10 ) ? 21 : 22 ) :
                    // For some reason, this seems to be faster than 22:
                    M1::_rowmajor ? ( 
                        docopy ? 82 :
                        ( s != UNKNOWN && s <= 10 ) ? 24 : 25 ) :
                    ( s != UNKNOWN && s <= 10 ) ? 21 : 22 ); 
#ifdef PRINTALGO_UD
            std::cout<<"InlineMultUD: \n";
            std::cout<<"x = "<<ix<<"  "<<T(x)<<"  add = "<<add<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"s,algo = "<<s<<"  "<<algo<<std::endl;
#endif
#ifdef XDEBUG_UD
            typedef typename M3::real_type RT;
            typedef typename M3::value_type T3;
            Matrix<T3> m1c = m1;
            Matrix<T3> m2c = m2;
            Matrix<T3> m3i = m3;
            Matrix<T3> m3c = m3;
            NoAliasMultMM<add>(x,m1c,m2c,m3c);
#endif
            MultUD_Helper<algo,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
#ifdef XDEBUG_UD
            if (Norm(m3-m3c) > 1.e-3*(Norm(m1c)*Norm(m2c)+(add?Norm(m3i):RT(0)))) {
                std::cout<<"m1 = "<<m1c<<std::endl;
                std::cout<<"m2 = "<<m2c<<std::endl;
                std::cout<<"m3 = "<<m3i<<std::endl;
                std::cout<<"m3 => "<<m3<<std::endl;
                std::cout<<"Correct m3 = "<<m3c<<std::endl;
                abort();
            }
#endif
        }
    };

    // algo -2: Check for inst
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<-2,s,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M3::value_type T3;
            const bool inst = 
                (s == UNKNOWN || s > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T3>::samebase &&
                Traits2<T2,T3>::samebase &&
#else
                Traits2<T1,T3>::sametype &&
                Traits2<T2,T3>::sametype &&
#endif
                Traits<T3>::isinst;
            const int algo = 
                ( s == 0 ) ? 0 :
                ( s == 1 ) ? 1 :
                M3::_conj ? 97 :
                inst ? 90 : 
                -3;
            MultUD_Helper<algo,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    // algo -1: Check for aliases?
    template <int s, bool add, int ix, class T, class M1, class M2, class M3>
    struct MultUD_Helper<-1,s,add,ix,T,M1,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, const M2& m2, M3& m3)
        {
            const int algo = 
                ( s == 0 ) ? 0 :
                ( s == 1 ) ? 1 :
                M3::_checkalias ? 99 : 
                -2;
            MultUD_Helper<algo,s,add,ix,T,M1,M2,M3>::call(x,m1,m2,m3);
        }
    };

    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3)
    {
        TMVStaticAssert(!M3::_unit);
        TMVStaticAssert(M1::_upper == int(M3::_upper));
        TMVStaticAssert((Sizes<M3::_size,M1::_size>::same));
        TMVStaticAssert((Sizes<M2::_size,M1::_size>::same));
        TMVStaticAssert((Sizes<M3::_size,M2::_size>::same));
        TMVAssert(m3.size() == m1.size());
        TMVAssert(m3.size() == m2.size());
        TMVAssert(!m3.isunit());
        const int s = Sizes<Sizes<M3::_size,M1::_size>::size,M2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        MultUD_Helper<-1,s,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3)
    {
        TMVStaticAssert(!M3::_unit);
        TMVStaticAssert(M1::_upper == int(M3::_upper));
        TMVStaticAssert((Sizes<M3::_size,M1::_size>::same));
        TMVStaticAssert((Sizes<M2::_size,M1::_size>::same));
        TMVStaticAssert((Sizes<M3::_size,M2::_size>::same));
        TMVAssert(m3.size() == m1.size());
        TMVAssert(m3.size() == m2.size());
        TMVAssert(!m3.isunit());
        const int s = Sizes<Sizes<M3::_size,M1::_size>::size,M2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        MultUD_Helper<-2,s,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void InlineMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3)
    {
        TMVStaticAssert(!M3::_unit);
        TMVStaticAssert(M1::_upper == int(M3::_upper));
        TMVStaticAssert((Sizes<M3::_size,M1::_size>::same));
        TMVStaticAssert((Sizes<M2::_size,M1::_size>::same));
        TMVStaticAssert((Sizes<M3::_size,M2::_size>::same));
        TMVAssert(m3.size() == m1.size());
        TMVAssert(m3.size() == m2.size());
        TMVAssert(!m3.isunit());
        const int s = Sizes<Sizes<M3::_size,M1::_size>::size,M2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        MultUD_Helper<-3,s,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void InlineAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3)
    {
        TMVStaticAssert(!M3::_unit);
        TMVStaticAssert(M1::_upper == int(M3::_upper));
        TMVStaticAssert((Sizes<M3::_size,M1::_size>::same));
        TMVStaticAssert((Sizes<M2::_size,M1::_size>::same));
        TMVStaticAssert((Sizes<M3::_size,M2::_size>::same));
        TMVAssert(m3.size() == m1.size());
        TMVAssert(m3.size() == m2.size());
        TMVAssert(!m3.isunit());
        const int s = Sizes<Sizes<M3::_size,M1::_size>::size,M2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        MultUD_Helper<98,s,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Tri_Mutable<M3>& m3)
    {
        TMVStaticAssert(!M3::_unit);
        TMVStaticAssert(M1::_upper == int(M3::_upper));
        TMVStaticAssert((Sizes<M3::_size,M1::_size>::same));
        TMVStaticAssert((Sizes<M2::_size,M1::_size>::same));
        TMVStaticAssert((Sizes<M3::_size,M2::_size>::same));
        TMVAssert(m3.size() == m1.size());
        TMVAssert(m3.size() == m2.size());
        TMVAssert(!m3.isunit());
        const int s = Sizes<Sizes<M3::_size,M1::_size>::size,M2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        MultUD_Helper<99,s,add,ix,T,M1v,M2v,M3v>::call(x,m1v,m2v,m3v);
    }

    template <class M1, int ix, class T, class M2>
    static TMV_INLINE void MultEqMM(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    { MultMM<false>(x,m1.mat(),m2.mat(),m1.mat()); }

    template <class M1, int ix, class T, class M2>
    static TMV_INLINE void NoAliasMultEqMM(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    { NoAliasMultMM<false>(x,m1.mat(),m2.mat(),m1.mat()); }

    template <class M1, int ix, class T, class M2>
    static TMV_INLINE void AliasMultEqMM(
        BaseMatrix_Tri_Mutable<M1>& m1,
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M2>& m2)
    { AliasMultMM<false>(x,m1.mat(),m2.mat(),m1.mat()); }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        const bool upper = M1::_upper;
        typedef typename TypeSelect<upper,
                typename M3::uppertri_type,
                typename M3::lowertri_type>::type M3u;
        M3u m3u = Maybe<upper>::uppertri(m3);
        MultMM<add>(x,m1,m2,m3u);
        Maybe<!add>::zero_offdiag2(Maybe<!upper>::uppertri(m3));
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        const bool upper = M1::_upper;
        typedef typename TypeSelect<upper,
                typename M3::uppertri_type,
                typename M3::lowertri_type>::type M3u;
        M3u m3u = Maybe<upper>::uppertri(m3);
        Maybe<!add>::zero(m3);
        NoAliasMultMM<add>(x,m1,m2,m3u);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        const BaseMatrix_Diag<M2>& m2, BaseMatrix_Rec_Mutable<M3>& m3)
    {
        const bool upper = M1::_upper;
        typedef typename TypeSelect<upper,
                typename M3::uppertri_type,
                typename M3::lowertri_type>::type M3u;
        M3u m3u = Maybe<upper>::uppertri(m3);
        AliasMultMM<add>(x,m1,m2,m3u);
        Maybe<!add>::zero_offdiag2(Maybe<!upper>::uppertri(m3));
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void MultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    {
        typename M3::transpose_type m3t = m3.transpose();
        MultMM<add>(x,m2.transpose(),m1.transpose(),m3t);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void NoAliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    {
        typename M3::transpose_type m3t = m3.transpose();
        NoAliasMultMM<add>(x,m2.transpose(),m1.transpose(),m3t);
    }

    template <bool add, int ix, class T, class M1, class M2, class M3>
    static inline void AliasMultMM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        const BaseMatrix_Tri<M2>& m2, BaseMatrix_Mutable<M3>& m3)
    {
        typename M3::transpose_type m3t = m3.transpose();
        AliasMultMM<add>(x,m2.transpose(),m1.transpose(),m3t);
    }


#undef TMV_UD_UNROLL
#undef TMV_ZeroIX

} // namespace tmv

#endif 
