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

#ifndef TMV_MultXU_H
#define TMV_MultXU_H

#include "TMV_BaseMatrix_Tri.h"
#include "TMV_CopyU.h"
#include "TMV_MultXV.h"

namespace tmv {

    // Defined in TMV_MultXU.cpp
    template <class T1, int C1, class T2>
    void InstMultXM(
        const T2 x, const ConstUpperTriMatrixView<T1,C1>& m1,
        UpperTriMatrixView<T2> m2);
    template <class T1,  int C1, class T2>
    void InstAddMultXM(
        const T2 x, const ConstUpperTriMatrixView<T1,C1>& m1, 
        UpperTriMatrixView<T2> m2);

    template <class T1, int C1, class T2>
    void InstAliasMultXM(
        const T2 x, const ConstUpperTriMatrixView<T1,C1>& m1,
        UpperTriMatrixView<T2> m2);
    template <class T1,  int C1, class T2>
    void InstAliasAddMultXM(
        const T2 x, const ConstUpperTriMatrixView<T1,C1>& m1, 
        UpperTriMatrixView<T2> m2);


    //
    // U (+)= x * U
    //

    template <bool add, int ix, class T, class M1, class M2>
    static inline void NoAliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1, 
        BaseMatrix_Tri_Mutable<M2>& m2);

    // The maximum nops to unroll.
#if TMV_OPT >= 3
#define TMV_XU_UNROLL 200 
#elif TMV_OPT >= 2
#define TMV_XU_UNROLL 25
#elif TMV_OPT >= 1
#define TMV_XU_UNROLL 9
#else
#define TMV_XU_UNROLL 0
#endif

    template <int algo, int s, bool add, int ix, class T, class M1, class M2>
    struct MultXU_Helper;

    // algo 0: trivial: ix == 1, !add, so call Copy
    template <int s, class T, class M1, class M2>
    struct MultXU_Helper<0,s,false,1,T,M1,M2>
    {
        static void call(const Scaling<1,T>& , const M1& m1, M2& m2)
        { CopyU_Helper<-1,s,M1,M2>::call(m1,m2); }
    };

    // algo 1: Transpose (and go back to -3, rather than -2)
    template <int s, bool add, int ix, class T, class M1, class M2>
    struct MultXU_Helper<1,s,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M1::const_transpose_type M1t;
            typedef typename M2::transpose_type M2t;
            M1t m1t = m1.transpose();
            M2t m2t = m2.transpose();
            MultXU_Helper<-3,s,add,ix,T,M1t,M2t>::call(x,m1t,m2t);
        }
    };

    // algo 2: UnknownDiag
    template <int s, bool add, int ix, class T, class M1, class M2>
    struct MultXU_Helper<2,s,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            if (m1.isunit())
                MultXU_Helper<3,s,add,ix,T,M1,M2>::call(x,m1,m2);
            else
                MultXU_Helper<-4,s,add,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    // algo 3: m1 is unitdiag
    template <int s, bool add, int ix, class T, class M1, class M2>
    struct MultXU_Helper<3,s,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            if (m2.size() > 1) {
                typedef typename M1::const_offdiag_type M1o;
                typedef typename M2::offdiag_type M2o;
                M1o m1o = m1.offDiag();
                M2o m2o = m2.offDiag();
                const int sm1 = IntTraits2<s,-1>::sum;
                MultXU_Helper<-2,sm1,add,ix,T,M1o,M2o>::call(x,m1o,m2o);
            }
            typedef typename M2::diag_type M2d;
            M2d m2d = m2.diag();
            Maybe<add>::addtoall(m2d,T(x));
        }
    };

    // algo 11: Loop over columns
    template <int s, bool add, int ix, class T, class M1, class M2>
    struct MultXU_Helper<11,s,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            int N = (s == UNKNOWN ? m2.size() : s);
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::col_sub_type M2c;
            typedef typename M1c::const_nonconj_type::const_iterator IT1;
            typedef typename M2c::iterator IT2;
            const int step1 = m1.stepj();
            const int step2 = m2.stepj();
            IT1 it1 = m1.get_col(0,0,1).begin().nonConj();
            IT2 it2 = m2.get_col(0,0,1).begin();
            int M=1;
            for(;N;--N) {
                MultXV_Helper<-4,UNKNOWN,add,ix,T,M1c,M2c>::call2(
                    M++,x,it1,it2);
                it1.shiftP(step1);
                it2.shiftP(step2);
            }
        }
    };

    // algo 12: Loop over rows
    template <int s, bool add, int ix, class T, class M1, class M2>
    struct MultXU_Helper<12,s,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            int N = (s == UNKNOWN ? m2.size() : s);
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M2::row_sub_type M2r;
            typedef typename M1r::const_nonconj_type::const_iterator IT1;
            typedef typename M2r::iterator IT2;
            const int step1 = m1.diagstep();
            const int step2 = m2.diagstep();
            IT1 it1 = m1.get_row(0,0,N).begin().nonConj();
            IT2 it2 = m2.get_row(0,0,N).begin();
            for(;N;--N) {
                MultXV_Helper<-4,UNKNOWN,add,ix,T,M1r,M2r>::call2(
                    N,x,it1,it2);
                it1.shiftP(step1);
                it2.shiftP(step2);
            }
        }
    };

    // algo 15: Fully unroll by columns
    template <int s, bool add, int ix, class T, class M1, class M2>
    struct MultXU_Helper<15,s,add,ix,T,M1,M2>
    {
        template <int I, int M, int J, int N>
        struct Unroller
        {
            static void unroll(
                const Scaling<ix,T>& x, const M1& m1, M2& m2)
            {
                Unroller<I,M-(N-N/2),J,N/2>::unroll(x,m1,m2);
                Unroller<I,M,J+N/2,N-N/2>::unroll(x,m1,m2);
            }
        };
        template <int I, int M, int J>
        struct Unroller<I,M,J,1>
        {
            static void unroll(
                const Scaling<ix,T>& x, const M1& m1, M2& m2)
            {
                Unroller<I,M/2,J,1>::unroll(x,m1,m2);
                Unroller<I+M/2,M-M/2,J,1>::unroll(x,m1,m2);
            }
        };
        template <int I, int M, int J>
        struct Unroller<I,M,J,0>
        {
            static void unroll(const Scaling<ix,T>& , const M1& , M2& ) 
            {} 
        };
        template <int I, int J>
        struct Unroller<I,1,J,1>
        {
            static void unroll(
                const Scaling<ix,T>& x, const M1& m1, M2& m2)
            {
                Maybe<add>::add( 
                    m2.ref(I,J), ZProd<false,false>::prod(x,m1.cref(I,J)) );
            }
        };
        template <int I, int J>
        struct Unroller<I,0,J,1>
        {
            static void unroll(const Scaling<ix,T>& , const M1& , M2& ) 
            {} 
        };

        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        { Unroller<0,s,0,s>::unroll(x,m1,m2); }
    };

    // algo 16: Fully unroll by rows
    template <int s, bool add, int ix, class T, class M1, class M2>
    struct MultXU_Helper<16,s,add,ix,T,M1,M2>
    {
        template <int I, int M, int J, int N>
        struct Unroller
        {
            static void unroll(
                const Scaling<ix,T>& x, const M1& m1, M2& m2)
            {
                Unroller<I,M/2,J,N>::unroll(x,m1,m2);
                Unroller<I+M/2,M-M/2,J+M/2,N-M/2>::unroll(x,m1,m2);
            }
        };
        template <int I, int J, int N>
        struct Unroller<I,1,J,N>
        {
            static void unroll(
                const Scaling<ix,T>& x, const M1& m1, M2& m2)
            {
                Unroller<I,1,J,N/2>::unroll(x,m1,m2);
                Unroller<I,1,J+N/2,N-N/2>::unroll(x,m1,m2);
            }
        };
        template <int I, int J, int N>
        struct Unroller<I,0,J,N>
        {
            static void unroll(const Scaling<ix,T>& , const M1& , M2& ) 
            {} 
        };
        template <int I, int J>
        struct Unroller<I,1,J,1>
        {
            static void unroll(
                const Scaling<ix,T>& x, const M1& m1, M2& m2)
            {
                Maybe<add>::add( 
                    m2.ref(I,J), ZProd<false,false>::prod(x,m1.cref(I,J)) );
            }
        };
        template <int I, int J>
        struct Unroller<I,1,J,0>
        {
            static void unroll(const Scaling<ix,T>& , const M1& , M2& ) 
            {} 
        };

        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        { Unroller<0,s,0,s>::unroll(x,m1,m2); }
    };

    // algo 90: Call inst
    template <int s, int ix, class T, class M1, class M2>
    struct MultXU_Helper<90,s,true,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M2::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAddMultXM(xx,m1.xView(),m2.xView());
        }
    };
    template <int s, int ix, class T, class M1, class M2>
    struct MultXU_Helper<90,s,false,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M2::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstMultXM(xx,m1.xView(),m2.xView());
        }
    };

    // algo 91: Call inst alias
    template <int s, int ix, class T, class M1, class M2>
    struct MultXU_Helper<91,s,true,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M2::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasAddMultXM(xx,m1.xView(),m2.xView());
        }
    };
    template <int s, int ix, class T, class M1, class M2>
    struct MultXU_Helper<91,s,false,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M2::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasMultXM(xx,m1.xView(),m2.xView());
        }
    };

    // algo 96: Transpose
    template <int s, bool add, int ix, class T, class M1, class M2>
    struct MultXU_Helper<96,s,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M1::const_transpose_type M1t;
            typedef typename M2::transpose_type M2t;
            M1t m1t = m1.transpose();
            M2t m2t = m2.transpose();
            MultXU_Helper<-2,s,add,ix,T,M1t,M2t>::call(x,m1t,m2t);
        }
    };

    // algo 196: Transpose
    template <int s, bool add, int ix, class T, class M1, class M2>
    struct MultXU_Helper<196,s,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M1::const_transpose_type M1t;
            typedef typename M2::transpose_type M2t;
            M1t m1t = m1.transpose();
            M2t m2t = m2.transpose();
            MultXU_Helper<99,s,add,ix,T,M1t,M2t>::call(x,m1t,m2t);
        }
    };

    // algo 97: Conjugate
    template <int s, bool add, int ix, class T, class M1, class M2>
    struct MultXU_Helper<97,s,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            MultXU_Helper<-2,s,add,ix,T,M1c,M2c>::call(TMV_CONJ(x),m1c,m2c);
        }
    };

    // algo 197: Conjugate
    template <int s, bool add, int ix, class T, class M1, class M2>
    struct MultXU_Helper<197,s,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            MultXU_Helper<99,s,add,ix,T,M1c,M2c>::call(TMV_CONJ(x),m1c,m2c);
        }
    };

    // algo 98: Inline check for aliases
    template <int s, int ix, class T, class M1, class M2>
    struct MultXU_Helper<98,s,false,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            if ( !SameStorage(m1,m2) ||
                 ExactSameStorage(m1,m2) ||
                 OppositeStorage(m1,m2) ) {
                // No aliasing (or no clobbering)
                MultXU_Helper<-2,s,false,ix,T,M1,M2>::call(x,m1,m2);
            } else {
                // Let Copy handle the aliasing
                AliasCopy(m1,m2);
                Scale(x,m2);
            }
        }
    };
    template <int s, int ix, class T, class M1, class M2>
    struct MultXU_Helper<98,s,true,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            if ( !SameStorage(m1,m2) ||
                 ExactSameStorage(m1,m2) ||
                 OppositeStorage(m1,m2) ) {
                // No aliasing (or no clobbering)
                MultXU_Helper<-2,s,true,ix,T,M1,M2>::call(x,m1,m2);
            } else {
                // Need a temporary
                NoAliasMultXM<true>(x,m1.copy(),m2);
            }
        }
    };

    // algo 99: Check for aliases
    template <int s, bool add, int ix, class T, class M1, class M2>
    struct MultXU_Helper<99,s,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            TMVStaticAssert(M1::_upper == int(M2::_upper));
            TMVStaticAssert(!M2::_unit || (!add && ix == 1 && M1::_unit));
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
                (ix == 1 && !add) ? 0 :
                M2::_lower ? 196 :
                M2::_conj ? 197 :
                inst ? 91 :
                98;
            MultXU_Helper<algo,s,add,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    // algo -4: No branches or copies
    // And m1 is NonUnitDiag
    template <int s, bool add, int ix, class T, class M1, class M2>
    struct MultXU_Helper<-4,s,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            TMVStaticAssert(M1::_upper);
            TMVStaticAssert(M2::_upper);
            TMVStaticAssert(!M1::_unit);
            TMVStaticAssert(!M2::_unit);
            TMVAssert(!m1.isunit());
            TMVAssert(!m2.isunit());
            typedef typename M2::value_type T2;
            const int s2 = s > 20 ? UNKNOWN : s;
            const int s2p1 = IntTraits<s2>::Sp1;
            // nops = n(n+1)/2
            const int nops = IntTraits2<s2,s2p1>::safeprod / 2;
            const bool unroll = 
                s > 10 ? false :
                s == UNKNOWN ? false :
                nops <= TMV_XU_UNROLL;
            const int algo = 
                unroll ? ( M2::_rowmajor ? 16 : 15 ) :
                M2::_rowmajor ? 12 : 11;
            MultXU_Helper<algo,s,add,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int s, bool add, int ix, class T, class M1, class M2>
    struct MultXU_Helper<-3,s,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            TMVStaticAssert(M1::_upper == int(M2::_upper));
            TMVStaticAssert(!M2::_unit || (!add && ix == 1 && M1::_unit));
            const int algo = 
                (ix == 1 && !add) ? 0 :
                M2::_lower ? 1 :
                M1::_unknowndiag ? 2 :
                M1::_unit ? 3 :
                -4;
            MultXU_Helper<algo,s,add,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    // algo -2: Check for inst
    template <int s, bool add, int ix, class T, class M1, class M2>
    struct MultXU_Helper<-2,s,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            TMVStaticAssert(M1::_upper == int(M2::_upper));
            TMVStaticAssert(!M2::_unit || (!add && ix == 1 && M1::_unit));
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
                (ix == 1 && !add) ? 0 :
                M2::_lower ? 96 :
                M2::_conj ? 97 :
                inst ? 90 :
                -3;
            MultXU_Helper<algo,s,add,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    // algo -1: Check for aliases?
    template <int s, bool add, int ix, class T, class M1, class M2>
    struct MultXU_Helper<-1,s,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            TMVStaticAssert(M1::_upper == int(M2::_upper));
            TMVStaticAssert(!M2::_unit || (!add && ix == 1 && M1::_unit));
            const bool noclobber =
                MStepHelper<M1,M2>::same ||
                MStepHelper<M1,M2>::opp;
            const bool checkalias =
                M2::_checkalias && !noclobber;
            const int algo =
                (ix == 1 && !add) ? 0 :
                checkalias ? 99 : 
                -2;
            MultXU_Helper<algo,s,add,ix,T,M1,M2>::call(x,m1,m2);
        }
    };


    template <bool add, int ix, class T, class M1, class M2>
    static inline void MultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1, 
        BaseMatrix_Tri_Mutable<M2>& m2)
    {
        TMVStaticAssert(M1::_upper == int(M2::_upper));
        TMVStaticAssert(!M2::_unit || (!add && ix == 1 && M1::_unit));
        TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m2.isunit() || (!add && ix == 1 && m1.isunit()));
        const int s = Sizes<M1::_size,M2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        MultXU_Helper<-1,s,add,ix,T,M1v,M2v>::call(x,m1v,m2v);
    }

    template <bool add, int ix, class T, class M1, class M2>
    static inline void NoAliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1, 
        BaseMatrix_Tri_Mutable<M2>& m2)
    {
        TMVStaticAssert(M1::_upper == int(M2::_upper));
        TMVStaticAssert(!M2::_unit || (!add && ix == 1 && M1::_unit));
        TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m2.isunit() || (!add && ix == 1 && m1.isunit()));
        const int s = Sizes<M1::_size,M2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        MultXU_Helper<-2,s,add,ix,T,M1v,M2v>::call(x,m1v,m2v);
    }

    template <bool add, int ix, class T, class M1, class M2>
    static inline void InlineMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1, 
        BaseMatrix_Tri_Mutable<M2>& m2)
    {
        TMVStaticAssert(M1::_upper == int(M2::_upper));
        TMVStaticAssert(!M2::_unit || (!add && ix == 1 && M1::_unit));
        TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m2.isunit() || (!add && ix == 1 && m1.isunit()));
        const int s = Sizes<M1::_size,M2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        MultXU_Helper<-3,s,add,ix,T,M1v,M2v>::call(x,m1v,m2v);
    }

    template <bool add, int ix, class T, class M1, class M2>
    static inline void InlineAliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1, 
        BaseMatrix_Tri_Mutable<M2>& m2)
    {
        TMVStaticAssert(M1::_upper == int(M2::_upper));
        TMVStaticAssert(!M2::_unit || (!add && ix == 1 && M1::_unit));
        TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m2.isunit() || (!add && ix == 1 && m1.isunit()));
        const int s = Sizes<M1::_size,M2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        MultXU_Helper<98,s,add,ix,T,M1v,M2v>::call(x,m1v,m2v);
    }

    template <bool add, int ix, class T, class M1, class M2>
    static inline void AliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1, 
        BaseMatrix_Tri_Mutable<M2>& m2)
    {
        TMVStaticAssert(M1::_upper == int(M2::_upper));
        TMVStaticAssert(!M2::_unit || (!add && ix == 1 && M1::_unit));
        TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
        TMVAssert(m1.size() == m2.size());
        TMVAssert(!m2.isunit() || (!add && ix == 1 && m1.isunit()));
        const int s = Sizes<M1::_size,M2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        MultXU_Helper<99,s,add,ix,T,M1v,M2v>::call(x,m1v,m2v);
    }

    template <class T1, int C1, class T2>
    static inline void InstMultXM(
        const T2 x, const ConstLowerTriMatrixView<T1,C1>& m1,
        LowerTriMatrixView<T2> m2)
    { InstMultXM(x,m1.transpose(),m2.transpose()); }
    template <class T1,  int C1, class T2>
    static inline void InstAddMultXM(
        const T2 x, const ConstLowerTriMatrixView<T1,C1>& m1, 
        LowerTriMatrixView<T2> m2)
    { InstAddMultXM(x,m1.transpose(),m2.transpose()); }

    //
    // M (+)= x * U
    //

    template <bool add, int ix, class T, class M1, class M2>
    static inline void MultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        BaseMatrix_Rec_Mutable<M2>& m2)
    {
        const bool upper = M1::_upper;
        typedef typename TypeSelect<upper,
                typename M2::uppertri_type,
                typename M2::lowertri_type>::type M2u;
        M2u m2u = Maybe<upper>::uppertri(m2);
        MultXM<add>(x,m1,m2u);
        Maybe<!add>::zero_offdiag2(Maybe<!upper>::uppertri(m2));
    }

    template <bool add, int ix, class T, class M1, class M2>
    static inline void NoAliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        BaseMatrix_Rec_Mutable<M2>& m2)
    {
        const bool upper = M1::_upper;
        typedef typename TypeSelect<upper,
                typename M2::uppertri_type,
                typename M2::lowertri_type>::type M2u;
        M2u m2u = Maybe<upper>::uppertri(m2);
        Maybe<!add>::zero(m2);
        NoAliasMultXM<add>(x,m1,m2u);
    }

    template <bool add, int ix, class T, class M1, class M2>
    static inline void AliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        BaseMatrix_Rec_Mutable<M2>& m2)
    {
        const bool upper = M1::_upper;
        typedef typename TypeSelect<upper,
                typename M2::uppertri_type,
                typename M2::lowertri_type>::type M2u;
        M2u m2u = Maybe<upper>::uppertri(m2);
        AliasMultXM<add>(x,m1,m2u);
        Maybe<!add>::zero_offdiag2(Maybe<!upper>::uppertri(m2));
    }


    //
    // U (+)= x * D
    //

    template <bool add, int ix, class T, class M1, class M2>
    static inline void MultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        BaseMatrix_Tri_Mutable<M2>& m2)
    {
        typename M1::const_diag_type m1d = m1.diag();
        typename M2::diag_type m2d = m2.diag();
        MultXV<add>(x,m1d,m2d);
        Maybe<!add>::zero_offdiag(m2);
    }

    template <bool add, int ix, class T, class M1, class M2>
    static inline void NoAliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        BaseMatrix_Tri_Mutable<M2>& m2)
    {
        typename M1::const_diag_type m1d = m1.diag();
        typename M2::diag_type m2d = m2.diag();
        Maybe<!add>::zero(m2);
        NoAliasMultXV<add>(x,m1d,m2d);
    }

    template <bool add, int ix, class T, class M1, class M2>
    static inline void AliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        BaseMatrix_Tri_Mutable<M2>& m2)
    {
        typename M1::const_diag_type m1d = m1.diag();
        typename M2::diag_type m2d = m2.diag();
        AliasMultXV<add>(x,m1d,m2d);
        Maybe<!add>::zero_offdiag(m2);
    }

#undef TMV_XU_UNROLL

} // namespace tmv

#endif
