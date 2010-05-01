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


#ifndef TMV_AddMM_H
#define TMV_AddMM_H

#include "TMV_MultXV.h"
#include "TMV_AddVV.h"
#include "TMV_CopyM.h"

namespace tmv {

    //
    // Matrix + Matrix
    //

    template <int algo, int cs, int rs, int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    struct AddMM_Helper;

    // algo 0: size == 0, nothing to do
    template <int cs, int rs, int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<0,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static inline void call(
            const Scaling<ix1,T1>& , const M1& , 
            const Scaling<ix2,T2>& , const M2& , M3& )
        {}
    };

    // algo 1: Linearize to vector version
    template <int cs, int rs, int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<1,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static inline void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            typedef typename M1::const_linearview_type M1l;
            typedef typename M2::const_linearview_type M2l;
            typedef typename M3::linearview_type M3l;
            M1l m1l = m1.linearView();
            M2l m2l = m2.linearView();
            M3l m3l = m3.linearView();
            InlineAddVV(x1,m1l,x2,m2l,m3l);
        }
    };

    // algo 11: Loop over columns
    template <int cs, int rs, int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<11,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            TMVStaticAssert(!M3::_conj);
            const int M = cs == UNKNOWN ? int(m2.colsize()) : cs;
            int N = rs == UNKNOWN ? int(m2.rowsize()) : rs;
            typedef typename M1::const_col_type M1c;
            typedef typename M2::const_col_type M2c;
            typedef typename M3::col_type M3c;
            typedef typename M1c::const_nonconj_type::const_iterator IT1;
            typedef typename M2c::const_nonconj_type::const_iterator IT2;
            typedef typename M3c::iterator IT3;
            TMVStaticAssert(!M3c::_conj);
            const int step1 = m1.stepj();
            const int step2 = m2.stepj();
            const int step3 = m3.stepj();
            IT1 it1 = m1.get_col(0).nonConj().begin();
            IT2 it2 = m2.get_col(0).nonConj().begin();
            IT3 it3 = m3.get_col(0).begin();
            for(;N;--N) {
                AddVV_Helper<-4,UNKNOWN,ix1,T1,M1c,ix2,T2,M2c,M3c>::call2(
                    M,x1,it1,x2,it2,it3);
                it1.shiftP(step1);
                it2.shiftP(step2);
                it3.shiftP(step3);
            }
        }
    };

    // algo 15: Fully unroll by columns
    template <int cs, int rs, int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<15,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        template <int I, int M, int J, int N>
        struct Unroller
        {
            static void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
            {
                Unroller<I,M,J,N/2>::unroll(x1,m1,x2,m2,m3);
                Unroller<I,M,J+N/2,N-N/2>::unroll(x1,m1,x2,m2,m3);
            }
        };
        template <int I, int M, int J>
        struct Unroller<I,M,J,1>
        {
            static void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
            {
                Unroller<I,M/2,J,1>::unroll(x1,m1,x2,m2,m3);
                Unroller<I+M/2,M-M/2,J,1>::unroll(x1,m1,x2,m2,m3);
            }
        };
        template <int I, int M, int J>
        struct Unroller<I,M,J,0>
        {
            static inline void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3) {}
        };
        template <int I, int J>
        struct Unroller<I,1,J,1>
        {
            static inline void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
            { m3.ref(I,J) = x1 * m1.cref(I,J) + x2 * m2.cref(I,J); }
        };
        template <int I, int J>
        struct Unroller<I,0,J,1>
        {
            static inline void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3) {}
        };
        static inline void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        { Unroller<0,cs,0,rs>::unroll(x1,m1,x2,m2,m3); }
    };

    // algo 21: Loop over rows
    template <int cs, int rs, int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<21,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            int M = cs == UNKNOWN ? int(m2.colsize()) : cs;
            const int N = rs == UNKNOWN ? int(m2.rowsize()) : rs;
            typedef typename M1::const_row_type M1r;
            typedef typename M2::const_row_type M2r;
            typedef typename M3::row_type M3r;
            typedef typename M1r::const_nonconj_type::const_iterator IT1;
            typedef typename M2r::const_nonconj_type::const_iterator IT2;
            typedef typename M3r::iterator IT3;
            const int step1 = m1.stepi();
            const int step2 = m2.stepi();
            const int step3 = m3.stepi();
            IT1 it1 = m1.get_row(0).nonConj().begin();
            IT2 it2 = m2.get_row(0).nonConj().begin();
            IT3 it3 = m3.get_row(0).begin();
            for(;M;--M) {
                AddVV_Helper<-4,UNKNOWN,ix1,T1,M1r,ix2,T2,M2r,M3r>::call2(
                    N,x1,it1,x2,it2,it3);
                it1.shiftP(step1);
                it2.shiftP(step2);
                it3.shiftP(step3);
            }
        }
    };

    // algo 25: Fully unroll by rows
    template <int cs, int rs, int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<25,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        template <int I, int M, int J, int N>
        struct Unroller
        {
            static void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
            {
                Unroller<I,M/2,J,N>::unroll(x1,m1,x2,m2,m3);
                Unroller<I+M/2,M-M/2,J,N>::unroll(x1,m1,x2,m2,m3);
            }
        };
        template <int I, int J, int N>
        struct Unroller<I,1,J,N>
        {
            static void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
            {
                Unroller<I,1,J,N/2>::unroll(x1,m1,x2,m2,m3);
                Unroller<I,1,J+N/2,N-N/2>::unroll(x1,m1,x2,m2,m3);
            }
        };
        template <int I, int J, int N>
        struct Unroller<I,0,J,N>
        {
            static inline void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3) {}
        };
        template <int I, int J>
        struct Unroller<I,1,J,1>
        {
            static inline void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
            { m3.ref(I,J) = x1 * m1.cref(I,J) + x2 * m2.cref(I,J); }
        };
        template <int I, int J>
        struct Unroller<I,1,J,0>
        {
            static inline void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3) {}
        };
        static inline void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        { Unroller<0,cs,0,rs>::unroll(x1,m1,x2,m2,m3); }
    };

    // algo 41: Unknown sizes, determine which algorithm to use
    template <int cs, int rs, int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<41,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            if (m1.canLinearize() && m2.canLinearize() && m3.canLinearize() &&
                m1.stepi() == m2.stepi() && m1.stepj() == m2.stepj() && 
                m1.stepi() == m3.stepi() && m1.stepj() == m3.stepj())
                AddMM_Helper<1,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(
                    x1,m1,x2,m2,m3);
            else if  (
                (m1.isrm() && m2.isrm()) || 
                (m1.isrm() && m3.isrm()) ||
                (m2.isrm() && m3.isrm()) ||
                ( !( (m1.iscm() && m2.iscm()) ||
                     (m1.iscm() && m3.iscm()) ||
                     (m2.iscm() && m3.iscm()) ) && 
                  m1.colsize() > m1.rowsize() ) ) 
                AddMM_Helper<21,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(
                        x1,m1,x2,m2,m3);
            else 
                AddMM_Helper<11,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(
                    x1,m1,x2,m2,m3);
        }
    };

    // algo -4: No branches or copies
    template <int cs, int rs, int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<-4,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static inline void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            TMVStaticAssert(!M3::_conj);
            typedef typename M3::value_type T3;
            const bool allrm = M1::_rowmajor && M2::_rowmajor && M3::_rowmajor;
#if TMV_OPT == 0 
            const int algo = allrm ? 21 : 11;
#else
            const bool allcm = M1::_colmajor && M2::_colmajor && M3::_colmajor;
            const bool tworm = ( 
                (M1::_rowmajor && M2::_rowmajor) ||
                (M1::_rowmajor && M3::_rowmajor) ||
                (M2::_rowmajor && M3::_rowmajor) );
            const bool twocm = ( 
                (M1::_colmajor && M2::_colmajor) ||
                (M1::_colmajor && M3::_colmajor) ||
                (M2::_colmajor && M3::_colmajor) );
            const bool canlin = 
                M1::_canlin && M2::_canlin && M3::_canlin &&
                ( allrm || allcm );
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                canlin ? 1 :
                ( cs != UNKNOWN && rs != UNKNOWN ) ? (
                    ( IntTraits2<cs,rs>::prod <= int(128/sizeof(T2)) ) ? (
                        tworm ? 25 : 15 ) :
                    tworm ? 21 : 
                    twocm ? 11 :
                    ( cs > rs ) ? 21 : 11 ) :
                tworm ? 21 : 11;
#endif
            AddMM_Helper<algo,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(
                x1,m1,x2,m2,m3);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<-3,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static inline void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            TMVStaticAssert(!M3::_conj);
#if TMV_OPT <= 1
            const int algo = -4;
#else
            const bool allrm = M1::_rowmajor && M2::_rowmajor && M3::_rowmajor;
            const bool allcm = M1::_colmajor && M2::_colmajor && M3::_colmajor;
            const bool canlin = 
                M1::_canlin && M2::_canlin && M3::_canlin &&
                ( allrm || allcm );
            const int algo = 
                cs == 0 || rs == 0 ? 0 :
                canlin ? 1 :
                ( cs == UNKNOWN || rs == UNKNOWN ) ? 41 :
                -4;
#endif
            AddMM_Helper<algo,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(
                x1,m1,x2,m2,m3);
        }
    };

    // algo 96: Conjugate
    template <int cs, int rs, int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<96,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static inline void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            typedef typename M3::conjugate_type M3c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            M3c m3c = m3.conjugate();
            AddMM_Helper<-2,cs,rs,ix1,T1,M1c,ix2,T2,M2c,M3c>::call(
                TMV_CONJ(x1),m1c,TMV_CONJ(x2),m2c,m3c);
        }
    };

    // algo 97: Call inst
    template <int cs, int rs, int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<97,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static inline void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            NoAliasMultXM<false>(x1,m1,m3);
            NoAliasMultXM<true>(x2,m2,m3);
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<-2,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static inline void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type TM1;
            typedef typename M2::value_type TM2;
            typedef typename M3::value_type TM3;
            const bool inst =
                M1::unknownsizes &&
                M2::unknownsizes &&
                M3::unknownsizes &&
#ifdef TMV_INST_MIX
                Traits2<TM1,TM2>::samebase &&
                Traits2<TM1,TM3>::samebase &&
#else
                Traits2<TM1,TM2>::sametype &&
                Traits2<TM1,TM3>::sametype &&
#endif
                Traits<TM3>::isinst;
            const bool conj = M3::_conj;
            const int algo = 
                conj ? 96 :
                inst ? 97 :
                -3;
            AddMM_Helper<algo,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(
                x1,m1,x2,m2,m3);
        }
    };

    // algo 98: Check for aliases when calling Inst functions
    // We don't have a separate Inst function for this.  We just
    // split the operation into two parts: 
    // m3 = x1*m1; 
    // m3 += x2*m2;
    // and let those operations call their Inst functions.
    // However, the alias requirements for this are different than the
    // normal ones, so we need to put some alias checking here.
    template <int cs, int rs, int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<98,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static inline void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            const bool s1 = SameStorage(m1,m3);
            const bool s2 = SameStorage(m2,m3);

            if (!s1 && !s2) {
                // No aliasing
                NoAliasMultXM<false>(x1,m1,m3);
                NoAliasMultXM<true>(x2,m2,m3);
            } else if (!s2) { 
                // Alias with m1 only, do m1 first
                AliasMultXM<false>(x1,m1,m3);
                NoAliasMultXM<true>(x2,m2,m3);
            } else if (!s1) { 
                // Alias with m2 only, do m2 first
                AliasMultXM<false>(x2,m2,m3);
                NoAliasMultXM<true>(x1,m1,m3);
            } else { 
                // Need a temporary
                typename M1::copy_type m1c = m1;
                AliasMultXM<false>(x2,m2,m3);
                NoAliasMultXM<true>(x1,m1c,m3);
            }
        }
    };

    // algo 99: Check for aliases when not calling Inst functions
    template <int cs, int rs, int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<99,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static inline void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            const bool ss1 = SameStorage(m1,m3);
            const bool ss2 = SameStorage(m2,m3);
            const bool s1 = ss1 && !ExactSameStorage(m1,m3);
            const bool s2 = ss2 && !ExactSameStorage(m2,m3);

            if (!s1 && !s2) {
                // No aliasing (or no clobbering)
                AddMM_Helper<-2,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(
                    x1,m1,x2,m2,m3);
            } else if (!ss2) { 
                // Alias with m1 only, do m1 first
                AliasMultXM<false>(x1,m1,m3);
                NoAliasMultXM<true>(x2,m2,m3);
            } else if (!ss1) { 
                // Alias with m2 only, do m2 first
                AliasMultXM<false>(x2,m2,m3);
                NoAliasMultXM<true>(x1,m1,m3);
            } else { 
                // Need a temporary
                typename M1::copy_type m1c = m1;
                AliasMultXM<false>(x2,m2,m3);
                NoAliasMultXM<true>(x1,m1c,m3);
            }
        }
    };

    // algo -1: Check for aliases?
    template <int cs, int rs, int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    struct AddMM_Helper<-1,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static inline void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            const bool noclobber = 
                MStepHelper<M1,M3>::same && MStepHelper<M2,M3>::same;
            const bool checkalias =
                M1::_colsize == UNKNOWN && M1::_rowsize == UNKNOWN &&
                M2::_colsize == UNKNOWN && M2::_rowsize == UNKNOWN &&
                M3::_colsize == UNKNOWN && M3::_rowsize == UNKNOWN &&
                !noclobber;
            typedef typename M1::value_type TM1;
            typedef typename M2::value_type TM2;
            typedef typename M3::value_type TM3;
            const bool inst =
                M1::unknownsizes &&
                M2::unknownsizes &&
                M3::unknownsizes &&
#ifdef TMV_INST_MIX
                Traits2<TM1,TM2>::samebase &&
                Traits2<TM1,TM3>::samebase &&
#else
                Traits2<TM1,TM2>::sametype &&
                Traits2<TM1,TM3>::sametype &&
#endif
                Traits<TM3>::isinst;
            const int algo = 
                // We do a different alias check for the Inst calls.
                inst ? 98 :
                checkalias ? 99 : 
                -2;
            AddMM_Helper<algo,cs,rs,ix1,T1,M1,ix2,T2,M2,M3>::call(
                x1,m1,x2,m2,m3);
        }
    };

    template <int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    inline void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Rec<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Rec<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.rowsize() == m3.rowsize());
        const int cs = 
            Sizes<Sizes<M1::_colsize,M2::_colsize>::size,M3::_colsize>::size;
        const int rs = 
            Sizes<Sizes<M1::_rowsize,M2::_rowsize>::size,M3::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        M1v m1v = m1.cView();
        M2v m2v = m2.cView();
        M3v m3v = m3.cView();
        AddMM_Helper<-1,cs,rs,ix1,T1,M1v,ix2,T2,M2v,M3v>::call(
            x1,m1v,x2,m2v,m3v);
    }

    template <int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    inline void NoAliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Rec<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Rec<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.rowsize() == m3.rowsize());
        const int cs = 
            Sizes<Sizes<M1::_colsize,M2::_colsize>::size,M3::_colsize>::size;
        const int rs = 
            Sizes<Sizes<M1::_rowsize,M2::_rowsize>::size,M3::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        M1v m1v = m1.cView();
        M2v m2v = m2.cView();
        M3v m3v = m3.cView();
        AddMM_Helper<-2,cs,rs,ix1,T1,M1v,ix2,T2,M2v,M3v>::call(
            x1,m1v,x2,m2v,m3v);
    }

    template <int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    inline void InlineAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Rec<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Rec<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.rowsize() == m3.rowsize());
        const int cs = 
            Sizes<Sizes<M1::_colsize,M2::_colsize>::size,M3::_colsize>::size;
        const int rs = 
            Sizes<Sizes<M1::_rowsize,M2::_rowsize>::size,M3::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        M1v m1v = m1.cView();
        M2v m2v = m2.cView();
        M3v m3v = m3.cView();
        AddMM_Helper<-3,cs,rs,ix1,T1,M1v,ix2,T2,M2v,M3v>::call(
            x1,m1v,x2,m2v,m3v);
    }

    template <int ix1, class T1, class M1,
              int ix2, class T2, class M2, class M3>
    inline void AliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Rec<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Rec<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_colsize,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M3::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.colsize() == m3.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.rowsize() == m3.rowsize());
        const int cs = 
            Sizes<Sizes<M1::_colsize,M2::_colsize>::size,M3::_colsize>::size;
        const int rs = 
            Sizes<Sizes<M1::_rowsize,M2::_rowsize>::size,M3::_rowsize>::size;
        typedef typename M1::value_type TM1;
        typedef typename M2::value_type TM2;
        typedef typename M3::value_type TM3;
        const bool inst =
            M1::unknownsizes &&
            M2::unknownsizes &&
            M3::unknownsizes &&
#ifdef TMV_INST_MIX
            Traits2<TM1,TM2>::samebase &&
            Traits2<TM1,TM3>::samebase &&
#else
            Traits2<TM1,TM2>::sametype &&
            Traits2<TM1,TM3>::sametype &&
#endif
            Traits<TM3>::isinst;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        M1v m1v = m1.cView();
        M2v m2v = m2.cView();
        M3v m3v = m3.cView();
        AddMM_Helper<inst?98:99,cs,rs,ix1,T1,M1v,ix2,T2,M2v,M3v>::call(
            x1,m1v,x2,m2v,m3v);
    }

} // namespace tmv

#endif 
