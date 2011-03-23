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


#ifndef TMV_DivVU_H
#define TMV_DivVU_H

#include "TMV_BaseMatrix_Tri.h"
#include "TMV_BaseVector.h"
#include "TMV_MultXV.h"
#include "TMV_MultVV.h"
#include "TMV_Prefetch.h"
#include "TMV_MultXV_Funcs.h"
#include "TMV_DivVM_Funcs.h"


#ifdef PRINTALGO_DivU
#include <iostream>
#include "TMV_MatrixIO.h"
#include "TMV_VectorIO.h"
#include "TMV_TriMatrixIO.h"
#endif

namespace tmv {

    // Defined in TMV_DivVU.cpp
    template <class T1, class T2, int C2>
    void InstLDivEq(
        VectorView<T1> v1, const ConstUpperTriMatrixView<T2,C2>& m2);
    template <class T1, class T2, int C2>
    void InstLDivEq(
        VectorView<T1> v1, const ConstLowerTriMatrixView<T2,C2>& m2);

    //
    // v1 /= m2
    // m2 is UpperTri or LowerTri
    //

    // The maximum nops to unroll.
#if TMV_OPT >= 3
#define TMV_DIVVU_UNROLL 200 
#elif TMV_OPT >= 2
#define TMV_DIVVU_UNROLL 25
#elif TMV_OPT >= 1
#define TMV_DIVVU_UNROLL 9
#else
#define TMV_DIVVU_UNROLL 0
#endif

    // The minimum size to copy a vector if its step == UNKNOWN.
#define TMV_DIVVU_COPYSIZE 4

    // The crossover memory size to start using prefetch commands.
    // This is undoubtedly a function of the L1 (and L2?) cache size,
    // but 2KBytes is probably not too bad for most machines.
    // (That's an empirical value for my Intel Core 2 Duo.)
#define TMV_DIVVU_PREFETCH 2048


    template <int algo, int s, class V1, class M2>
    struct LDivEqVU_Helper;

    template <int algo, int s, class V1, class M2>
    struct LDivEqVU_Helper;

    // algo 0: s = 0, so nothing to do
    template <class V1, class M2>
    struct LDivEqVU_Helper<0,0,V1,M2>
    { static void call(V1& , const M2& ) {} };

    // algo 1: s == 1, so simplifies to a scalar quotient
    template <class V1, class M2>
    struct LDivEqVU_Helper<1,1,V1,M2>
    {
        static void call(V1& v1, const M2& m2)
        {
#ifdef PRINTALGO_DivU
            std::cout<<"LDivEqVU algo 1: N,s = "<<1<<','<<1<<std::endl;
#endif
            Maybe<!M2::_unit>::invscale(v1.ref(0) , m2.cref(0,0)); 
        }
    };

    // algo 11: The basic column major loop for UpperTri
    template <int s, class V1, class M2>
    struct LDivEqVU_Helper<11,s,V1,M2>
    {
        static void call(V1& v1, const M2& m2)
        {
            const int N = (s == UNKNOWN ? m2.size() : s);
            if (N == 0) return;
#ifdef PRINTALGO_DivU
            std::cout<<"LDivEqVU algo 11: N,s = "<<N<<','<<s<<std::endl;
#endif
            typedef typename V1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M2::const_col_sub_type M2c;
            typedef typename M2::const_diag_type M2d;
            T1 Xj;

            const bool unit = M2::_unit;

            typedef typename V1::iterator IT1;
            typedef typename M2c::const_nonconj_type::const_iterator IT2;

            IT1 X0 = v1.begin();
            IT1 X = X0 + N-1;
            IT2 A0j = m2.get_col(N-1,0,N).begin().nonConj();
            const int Astepj = m2.stepj();

            const bool dopref = N * Astepj * sizeof(T1) >= TMV_DIVVU_PREFETCH;

            Prefetch_MultiWrite(v1.ptr());
            if (dopref) Prefetch_Read(A0j.get());
            else Prefetch_Read(m2.cptr());

            for(int j=N;j--;) {
                // loop from j = N-1 .. 0
                if (*X != T2(0)) {
                    // x(j) = x(j) / A(j,j);
                    Xj = Maybe<!unit>::invprod(m2.cref(j,j),*X);
                    *X-- = Xj;
                    // x.subVector(0,j) -= x(j) * A.col(j,0,j);
                    MultXV_Helper<-4,UNKNOWN,true,0,T1,M2c,V1>::call2(
                        j,Scaling<0,T1>(-Xj),A0j,X0);
                    A0j.shiftP(-Astepj);
                    if (dopref) Prefetch_Read(A0j.get());
                } else {
                    A0j.shiftP(-Astepj);
                    if (dopref) Prefetch_Read(A0j.get());
                    --X; 
                }
            }
        }
    };

    // algo 12: The basic row major loop for UpperTri
    template <int s, class V1, class M2>
    struct LDivEqVU_Helper<12,s,V1,M2>
    {
        static void call(V1& v1, const M2& m2)
        {
            int N = (s == UNKNOWN ? m2.size() : s);
            if (N == 0) return;
#ifdef PRINTALGO_DivU
            std::cout<<"LDivEqVU algo 12: N,s = "<<N<<','<<s<<std::endl;
#endif
            typedef typename V1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M2::const_row_sub_type M2r;
            T1 Xi;

            typedef typename V1::iterator IT1;
            typedef typename M2r::const_nonconj_type::const_iterator IT2;

            const bool unit = M2::_unit;

            IT1 X = v1.begin() + N-1;
            //IT2 A00 = m2.get_row(0,0,N).begin().nonConj();
            // Actually Aii is the address of A(i,i+1)
            IT2 Aii = m2.get_row(N-1,N,N).begin().nonConj();
            const int Adiagstep = m2.diagstep();

            const bool dopref = N * m2.stepi() * sizeof(T1) >= TMV_DIVVU_PREFETCH;

            Prefetch_Write(v1.ptr());
            if (dopref) Prefetch_Read(Aii.get());
            else Prefetch_Read(m2.cptr());

            // [  A11 A12 ] [ X1 ] = [ A11 X1 ]
            // [   0  A22 ] [  0 ]   [    0   ]
            // Unlike with MultUV, there is no gain from checking for 0's
            // at the start of v1.
            int i=N-1;
            while (N > 0 && *X == T2(0)) {
                --N; --X; --i; Aii.shiftP(-Adiagstep); 
            }

            for(int len=0;N--;--i) {
                // loop from i = N-1..0
                // x(i) = A.row(i,i,N) * x.subVector(i,N)
                //      = A(i,i)^-1 * (
                //        x(i) - A.row(i,i+1,N) * x.subVector(i+1,N) )
                Xi = -MultVV_Helper<-4,UNKNOWN,M2r,V1>::call2(len++,Aii,X+1);
                Xi += *X;
                Maybe<!unit>::invscale(Xi,m2.cref(i,i));
                Aii.shiftP(-Adiagstep);
                if (dopref) Prefetch_Read(Aii.get()-1);
                *X-- = Xi;
            }
        } 
    };

    // algo 15: UpperTri: unroll by rows
    template <int s, class V1, class M2>
    struct LDivEqVU_Helper<15,s,V1,M2>
    {
        template <int I, int N>
        struct Unroller
        {
            static void unroll(V1& v1, const M2& m2)
            {
                Unroller<I+N/2,N-N/2>::unroll(v1,m2);
                Unroller<I,N/2>::unroll(v1,m2);
            }
        };
        template <int I>
        struct Unroller<I,1>
        {
            static void unroll(V1& v1, const M2& m2)
            {
                typedef typename V1::const_subvector_type V1s;
                typedef typename M2::const_row_sub_type M2r;
                const bool unit = M2::_unit;
                v1.ref(I) = Maybe<!unit>::invprod(
                    m2.cref(I,I) , v1.cref(I) -
                    MultVV_Helper<-4,s-I-1,M2r,V1s>::call(
                        m2.get_row(I,I+1,s),v1.cSubVector(I+1,s)));
            }
        };
        template <int I>
        struct Unroller<I,0>
        { static void unroll(V1& , const M2& ) {} };
        static void call(V1& v1, const M2& m2)
        {
#ifdef PRINTALGO_DivU
            std::cout<<"LDivEqVU algo 15: N,s = "<<s<<','<<s<<std::endl;
#endif
            Unroller<0,s>::unroll(v1,m2); 
        }
    };

    // algo 21: The basic column major loop for LowerTri
    template <int s, class V1, class M2>
    struct LDivEqVU_Helper<21,s,V1,M2>
    {
        static void call(V1& v1, const M2& m2)
        {
            int N = (s == UNKNOWN ? m2.size() : s);
            if (N == 0) return;
#ifdef PRINTALGO_DivU
            std::cout<<"LDivEqVU algo 21: N,s = "<<N<<','<<s<<std::endl;
#endif
            typedef typename V1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M2::const_col_sub_type M2c;
            typedef typename M2::const_diag_type M2d;
            T1 Xj;

            const bool unit = M2::_unit;

            typedef typename V1::iterator IT1;
            typedef typename M2c::const_nonconj_type::const_iterator IT2;

            // Actually Ajj is the address of A(j,j+1)
            IT1 X = v1.begin();
            IT2 Ajj = m2.get_col(0,1,N).begin().nonConj();
            const int Adiagstep = m2.diagstep();

            const bool dopref = N * sizeof(T1) >= TMV_DIVVU_PREFETCH;

            Prefetch_MultiWrite(v1.ptr());
            Prefetch_Read(m2.cptr());

            for(int j=0;N--;++j) {
                // loop from j = 0 .. N-1
                if (*X != T2(0)) {
                    // x(j) = x(j) / A(j,j);
                    Xj = Maybe<!unit>::invprod(m2.cref(j,j),*X);
                    *X++ = Xj;
                    // y.subVector(j+1,N) -= x(j) * A.col(j,j+1,N);
                    MultXV_Helper<-4,UNKNOWN,true,0,T1,M2c,V1>::call2(
                        N,Scaling<0,T1>(-Xj),Ajj,X);
                    Ajj.shiftP(Adiagstep);
                    if (dopref) Prefetch_Read(Ajj.get()-1);
                } else {
                    ++X;
                    Ajj.shiftP(Adiagstep);
                    if (dopref) Prefetch_Read(Ajj.get()-1);
                }
            }
        }
    };

    // algo 22: The basic row major loop for LowerTri
    template <int s, class V1, class M2>
    struct LDivEqVU_Helper<22,s,V1,M2>
    {
        static void call(V1& v1, const M2& m2)
        {
            int N = (s == UNKNOWN ? m2.size() : s);
            if (N == 0) return;
#ifdef PRINTALGO_DivU
            std::cout<<"LDivEqVU algo 22: N,s = "<<N<<','<<s<<std::endl;
#endif
            typedef typename V1::value_type T1;
            typedef typename M2::value_type T2;
            typedef typename M2::const_row_sub_type M2r;
            T1 Xi;

            typedef typename V1::iterator IT1;
            typedef typename M2r::const_nonconj_type::const_iterator IT2;

            const bool unit = M2::_unit;

            IT1 X0 = v1.begin();
            IT1 X = X0;
            IT2 Ai0 = m2.get_row(0,0,1).begin().nonConj();
            const int Astepi = m2.stepi();
            const int Adiagstep = m2.diagstep();

            const bool dopref = N * sizeof(T1) >= TMV_DIVVU_PREFETCH;

            Prefetch_Write(v1.ptr());
            Prefetch_Read(m2.cptr());

            // [ A00  0  ] [  0 ]   [    0   ]
            // [ A10 A11 ] [ B1 ] = [ A11 B1 ]
            int i=0;
            while (N > 0 && *X == T2(0)) {
                --N; ++X; ++X0; ++i; Ai0.shiftP(Adiagstep);
            }
            for(int len=0;N--;++i) {
                // loop from i = 0 .. N-1
                // Yi = A.row(i,0,i+1) * x.subVector(0,i+1)
                //    = A(i,i)^-1 * (
                //        x(i) - A.row(i,0,i) * x.subVector(0,i))
                Xi = -MultVV_Helper<-4,UNKNOWN,M2r,V1>::call2(len++,Ai0,X0);
                Ai0.shiftP(Astepi);
                if (dopref) Prefetch_Read(Ai0.get());
                Xi += *X;
                Maybe<!unit>::invscale(Xi,m2.cref(i,i));
                *X++ = Xi;
            } 
        }
    };

    // algo 25: LowerTri: unroll by rows
    template <int s, class V1, class M2>
    struct LDivEqVU_Helper<25,s,V1,M2>
    {
        template <int I, int N>
        struct Unroller
        {
            static void unroll(V1& v1, const M2& m2)
            {
                Unroller<I,N/2>::unroll(v1,m2);
                Unroller<I+N/2,N-N/2>::unroll(v1,m2);
            }
        };
        template <int I>
        struct Unroller<I,1>
        {
            static void unroll(V1& v1, const M2& m2)
            {
                typedef typename M2::const_row_sub_type M2r;
                typedef typename V1::const_subvector_type V1s;
                const bool unit = M2::_unit;
                v1.ref(I) = Maybe<!unit>::invprod(
                    m2.cref(I,I) , v1.cref(I) -
                    MultVV_Helper<-4,I,M2r,V1s>::call(
                        m2.get_row(I,0,I),v1.cSubVector(0,I)));
            }
        };
        template <int I>
        struct Unroller<I,0>
        { static void unroll(V1& , const M2& ) {} };
        static void call(V1& v1, const M2& m2)
        {
#ifdef PRINTALGO_DivU
            std::cout<<"LDivEqVU algo 25: N,s = "<<s<<','<<s<<std::endl;
#endif
            Unroller<0,s>::unroll(v1,m2); 
        }
    };

    // algo 43: v1.step == UNKNOWN, so maybe copy v1
    template <int s, class V1, class M2>
    struct LDivEqVU_Helper<43,s,V1,M2>
    {
        static void call(V1& v1, const M2& m2)
        {
            const int N = s == UNKNOWN ? int(m2.size()) : s;
#ifdef PRINTALGO_DivU
            std::cout<<"LDivEqVU algo 43: N,s = "<<N<<','<<s<<std::endl;
#endif
            if (N > TMV_DIVVU_COPYSIZE) {
                LDivEqVU_Helper<85,s,V1,M2>::call(v1,m2);
            } else 
                LDivEqVU_Helper<-4,s,V1,M2>::call(v1,m2);
        }
    };

    // algo 85: v1c = v1/m2, v1 = v1c
    template <int s, class V1, class M2>
    struct LDivEqVU_Helper<85,s,V1,M2>
    {
        static void call(V1& v1, const M2& m2)
        {
#ifdef PRINTALGO_DivU
            const int N = s == UNKNOWN ? int(m2.size()) : s;
            std::cout<<"LDivEqVU algo 85: N,s = "<<N<<','<<s<<std::endl;
#endif
            typename V1::copy_type v1c(v1);
            NoAliasLDivEq(v1c,m2);
            NoAliasCopy(v1c,v1);
        }
    };

    // algo 90: Call inst
    template <int s, class V1, class M2>
    struct LDivEqVU_Helper<90,s,V1,M2>
    {
        static void call(V1& v1, const M2& m2)
        { InstLDivEq(v1.xView(),m2.xView()); }
    };

    // algo 97: Conjugate
    template <int s, class V1, class M2>
    struct LDivEqVU_Helper<97,s,V1,M2>
    {
        static void call(V1& v1, const M2& m2)
        {
            typedef typename V1::conjugate_type V1c;
            typedef typename M2::const_conjugate_type M2c;
            V1c v1c = v1.conjugate();
            M2c m2c = m2.conjugate();
            LDivEqVU_Helper<-2,s,V1c,M2c>::call(v1c,m2c);
        }
    };

    // algo 99: Check for aliases
    template <int s, class V1, class M2>
    struct LDivEqVU_Helper<99,s,V1,M2>
    {
        static void call(V1& v1, const M2& m2)
        {
            if ( !SameStorage(v1,m2) ) {
                // No aliasing
                LDivEqVU_Helper<-2,s,V1,M2>::call(v1,m2);
            } else {
                // Use temporary for v1
                LDivEqVU_Helper<85,s,V1,M2>::call(v1,m2);
            }
        }
    };

    // algo -4: No branches or copies
    template <int s, class V1, class M2>
    struct LDivEqVU_Helper<-4,s,V1,M2>
    {
        static void call(V1& v1, const M2& m2)
        {
            TMVStaticAssert(!V1::_conj);
            const int s2 = s > 20 ? UNKNOWN : s;
            const int s2p1 = IntTraits<s2>::Sp1;
            // nops = n(n+1)/2
            const int nops = IntTraits2<s2,s2p1>::safeprod / 2;
            const bool unroll = 
                s > 10 ? false :
                s == UNKNOWN ? false :
                nops <= TMV_DIVVU_UNROLL;
            const int algo = 
                ( s == 0 ) ? 0 : 
                ( s == 1 ) ? 1 : 
                M2::_upper ? (
                    unroll ? 15 :
                    M2::_colmajor ? 11 : 12 ) :
                ( // lowertri
                    unroll ? 25 :
                    M2::_colmajor ? 21 : 22 );
            LDivEqVU_Helper<algo,s,V1,M2>::call(v1,m2); 
        }
    };

    // algo -3: Determine which algorithm to use
    template <int s, class V1, class M2>
    struct LDivEqVU_Helper<-3,s,V1,M2>
    {
        static void call(V1& v1, const M2& m2)
        {
            typedef typename V1::value_type T1;
            typedef typename M2::value_type T2;
            TMVStaticAssert(!V1::_conj);
            // Possible algorithms to choose from:
            //
            // Trivial:
            //  0 = s == 0, so nothing to do
            //  1 = s == 1: reduces to trivial scalar division
            //
            // UpperTri:
            // 11 = column major, simple for loop
            // 12 = row major, simple for loop
            // 15 = fully unroll by rows
            //
            // LowerTri:
            // 21 = column major, simple for loop
            // 22 = row major, simple for loop
            // 25 = fully unroll by rows
            //
            // Copy vector to new storage:
            // 85 = temp v1 = v1/m2

            const int s2 = s > 20 ? UNKNOWN : s;
            const int s2p1 = IntTraits<s2>::Sp1;
            // nops = n(n+1)/2
            const int nops = IntTraits2<s2,s2p1>::safeprod / 2;
            const bool unroll = 
                s > 10 ? false :
                s == UNKNOWN ? false :
                nops <= TMV_DIVVU_UNROLL;
            const int algo = 
                ( s == 0 ) ? 0 : // trivial - nothing to do
                ( s == 1 ) ? 1 : // trivial - s = 1
                !Traits2<T1,T2>::samebase ?
                ( M2::_upper ? M2::_colmajor?11:12 : M2::_colmajor?21:22 ) :

                M2::_upper ? 
                TMV_OPT == 0 ? ( M2::_colmajor ? 11 : 12 ) :
                unroll ? 15 :
                M2::_colmajor ? (
                    V1::_step == UNKNOWN ? (
#ifdef TMV_OPT_SCALE
                        s == UNKNOWN ? 43 :
#endif
                        s == UNKNOWN ? 85 :
                        s > TMV_DIVVU_COPYSIZE ? 85 : 11 ) :
                    11 ) :
                M2::_rowmajor ? (
                    V1::_step == UNKNOWN ? (
#ifdef TMV_OPT_SCALE
                        s == UNKNOWN ? 43 :
#endif
                        s == UNKNOWN ? 85 :
                        s > TMV_DIVVU_COPYSIZE ? 85 : 12 ) :
                    12 ) :
                12 :

                // lowertri
                TMV_OPT == 0 ? ( M2::_colmajor ? 21 : 22 ) : 
                unroll ? 25 :
                M2::_colmajor ? (
                    V1::_step == UNKNOWN ? (
#ifdef TMV_OPT_SCALE
                        s == UNKNOWN ? 43 :
#endif
                        s == UNKNOWN ? 85 :
                        s > TMV_DIVVU_COPYSIZE ? 85 : 21 ) :
                    21 ) :
                M2::_rowmajor ? (
                    V1::_step == UNKNOWN ? (
#ifdef TMV_OPT_SCALE
                        s == UNKNOWN ? 43 :
#endif
                        s == UNKNOWN ? 85 :
                        s > TMV_DIVVU_COPYSIZE ? 85 : 22 ) :
                    22 ) :
                22;
#ifdef PRINTALGO_DivU
            std::cout<<"InlineLDivEqVU: \n";
            std::cout<<"v1 = "<<TMV_Text(v1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"s,algo = "<<s<<"  "<<algo<<std::endl;
#endif
            if (m2.isSingular()) {
#ifdef TMV_NO_THROW
                std::cerr<<"Singular TriMatrix found\n";
                exit(1);
#else
                throw Singular("TriMatrix found\n");
#endif
            }
            LDivEqVU_Helper<algo,s,V1,M2>::call(v1,m2);
        }
    };

    // algo -2: Check for inst
    template <int s, class V1, class M2>
    struct LDivEqVU_Helper<-2,s,V1,M2>
    {
        static void call(V1& v1, const M2& m2)
        {
            typedef typename V1::value_type T1;
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
                ( s == 1 ) ? 1 :
                V1::_conj ? 97 : 
                inst ? 90 : 
                -3;
            LDivEqVU_Helper<algo,s,V1,M2>::call(v1,m2); 
        }
    };

    // algo -1: Check for aliases?
    template <int s, class V1, class M2>
    struct LDivEqVU_Helper<-1,s,V1,M2>
    {
        static void call(V1& v1, const M2& m2)
        {
            const bool checkalias =
                V1::_size == UNKNOWN &&
                M2::_size == UNKNOWN;
            const int algo =
                ( s == 0 ) ? 0 : 
                ( s == 1 ) ? 1 :
                checkalias ? 99 : 
                -2;
            LDivEqVU_Helper<algo,s,V1,M2>::call(v1,m2); 
        }
    };

    template <class V1, class M2>
    static inline void LDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Tri<M2>& m2)
    {
        TMVStaticAssert((Sizes<V1::_size,M2::_size>::same));
        TMVAssert(v1.size() == m2.size());
        const int s = Sizes<V1::_size,M2::_size>::size;
        typedef typename V1::cview_type V1v;
        typedef typename M2::const_cview_type M2v;
        TMV_MAYBE_REF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        LDivEqVU_Helper<-1,s,V1v,M2v>::call(v1v,m2v);
    }
    template <class V1, class M2>
    static inline void NoAliasLDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Tri<M2>& m2)
    {
        TMVStaticAssert((Sizes<V1::_size,M2::_size>::same));
        TMVAssert(v1.size() == m2.size());
        const int s = Sizes<V1::_size,M2::_size>::size;
        typedef typename V1::cview_type V1v;
        typedef typename M2::const_cview_type M2v;
        TMV_MAYBE_REF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        LDivEqVU_Helper<-2,s,V1v,M2v>::call(v1v,m2v);
    }
    template <class V1, class M2>
    static inline void InlineLDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Tri<M2>& m2)
    {
        TMVStaticAssert((Sizes<V1::_size,M2::_size>::same));
        TMVAssert(v1.size() == m2.size());
        const int s = Sizes<V1::_size,M2::_size>::size;
        typedef typename V1::cview_type V1v;
        typedef typename M2::const_cview_type M2v;
        TMV_MAYBE_REF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        LDivEqVU_Helper<-3,s,V1v,M2v>::call(v1v,m2v);
    }
    template <class V1, class M2>
    static inline void AliasLDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Tri<M2>& m2)
    {
        TMVStaticAssert((Sizes<V1::_size,M2::_size>::same));
        TMVAssert(v1.size() == m2.size());
        const int s = Sizes<V1::_size,M2::_size>::size;
        typedef typename V1::cview_type V1v;
        typedef typename M2::const_cview_type M2v;
        TMV_MAYBE_REF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        LDivEqVU_Helper<99,s,V1v,M2v>::call(v1v,m2v);
    }

    //
    // v3 = v1 / m2
    // 

    // This is an abbreviated set that just checks for aliases
    // in the right situations.  Then it calls LDivEq once everything is 
    // known to be ok alias-wise.
    template <int algo, int ix, class T, class V1, class M2, class V3>
    struct LDivVU_Helper;

    // algo 99: Check for aliases
    template <int ix, class T, class V1, class M2, class V3>
    struct LDivVU_Helper<99,ix,T,V1,M2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const M2& m2, V3& v3)
        {
            if ( !SameStorage(m2,v3) ) {
                AliasMultXV<false>(x,v1,v3);
                NoAliasLDivEq(v3,m2);
            } else {
                typename V3::copy_type v3c(v3.size());
                NoAliasMultXV<false>(x,v1,v3c);
                NoAliasLDivEq(v3c,m2);
                NoAliasCopy(v3c,v3);
            }
        }
    };

    // algo -2: NoAlias: Move along to LDivEq
    template <int ix, class T, class V1, class M2, class V3>
    struct LDivVU_Helper<-2,ix,T,V1,M2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const M2& m2, V3& v3)
        {
            NoAliasMultXV<false>(x,v1,v3);
            NoAliasLDivEq(v3,m2);
        }
    };

    // algo -1: Check for aliases?
    template <int ix, class T, class V1, class M2, class V3>
    struct LDivVU_Helper<-1,ix,T,V1,M2,V3>
    {
        static void call(
            const Scaling<ix,T>& x, const V1& v1, const M2& m2, V3& v3)
        {
            const bool checkalias =
                V1::_size == UNKNOWN &&
                M2::_size == UNKNOWN &&
                V3::_size == UNKNOWN;
            const int algo =
                checkalias ? 99 : 
                -2;
            LDivVU_Helper<algo,ix,T,V1,M2,V3>::call(x,v1,m2,v3); 
        }
    };

    template <int ix, class T, class V1, class M2, class V3>
    static inline void LDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Tri<M2>& m2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::_size,M2::_size>::same));
        TMVStaticAssert((Sizes<V1::_size,V3::_size>::same));
        TMVAssert(v1.size() == m2.size());
        TMVAssert(v1.size() == v3.size());
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        LDivVU_Helper<-1,ix,T,V1v,M2v,V3v>::call(x,v1v,m2v,v3v);
    }
    template <int ix, class T, class V1, class M2, class V3>
    static inline void NoAliasLDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Tri<M2>& m2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::_size,M2::_size>::same));
        TMVStaticAssert((Sizes<V1::_size,V3::_size>::same));
        TMVAssert(v1.size() == m2.size());
        TMVAssert(v1.size() == v3.size());
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        LDivVU_Helper<-2,ix,T,V1v,M2v,V3v>::call(x,v1v,m2v,v3v);
    }
    template <int ix, class T, class V1, class M2, class V3>
    static inline void InlineLDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Tri<M2>& m2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::_size,M2::_size>::same));
        TMVStaticAssert((Sizes<V1::_size,V3::_size>::same));
        TMVAssert(v1.size() == m2.size());
        TMVAssert(v1.size() == v3.size());
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        LDivVU_Helper<-3,ix,T,V1v,M2v,V3v>::call(x,v1v,m2v,v3v);
    }
    template <int ix, class T, class V1, class M2, class V3>
    static inline void AliasLDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Tri<M2>& m2, BaseVector_Mutable<V3>& v3)
    {
        TMVStaticAssert((Sizes<V1::_size,M2::_size>::same));
        TMVStaticAssert((Sizes<V1::_size,V3::_size>::same));
        TMVAssert(v1.size() == m2.size());
        TMVAssert(v1.size() == v3.size());
        typedef typename V1::const_cview_type V1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename V3::cview_type V3v;
        TMV_MAYBE_CREF(V1,V1v) v1v = v1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(V3,V3v) v3v = v3.cView();
        LDivVU_Helper<99,ix,T,V1v,M2v,V3v>::call(x,v1v,m2v,v3v);
    }

    //
    // v1 %= m2
    //

    template <class V1, class M2>
    static inline void RDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Tri<M2>& m2)
    { LDivEq(v1,m2.transpose()); }
    template <class V1, class M2>
    static inline void NoAliasRDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Tri<M2>& m2)
    { NoAliasLDivEq(v1,m2.transpose()); }
    template <class V1, class M2>
    static inline void AliasRDivEq(
        BaseVector_Mutable<V1>& v1, const BaseMatrix_Tri<M2>& m2)
    { AliasLDivEq(v1,m2.transpose()); }

    //
    // v3 = v1 % m2
    //

    template <int ix, class T, class V1, class M2, class V3>
    static inline void RDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Tri<M2>& m2, BaseVector_Mutable<V3>& v3)
    { LDiv(x,v1,m2.transpose(),v3); }
    template <int ix, class T, class V1, class M2, class V3>
    static inline void NoAliasRDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Tri<M2>& m2, BaseVector_Mutable<V3>& v3)
    { NoAliasLDiv(x,v1,m2.transpose(),v3); }
    template <int ix, class T, class V1, class M2, class V3>
    static inline void AliasRDiv(
        const Scaling<ix,T>& x, const BaseVector_Calc<V1>& v1,
        const BaseMatrix_Tri<M2>& m2, BaseVector_Mutable<V3>& v3)
    { AliasLDiv(x,v1,m2.transpose(),v3); }

#undef TMV_DIVVU_UNROLL
#undef TMV_DIVVU_COPYSIZE
#undef TMV_DIVVU_PREFETCH

} // namespace tmv

#endif 
