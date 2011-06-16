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

#ifndef TMV_CopyU_H
#define TMV_CopyU_H

#include "TMV_BaseMatrix_Tri.h"
#include "TMV_CopyV.h"

namespace tmv {

    // Defined in TMV_TriMatrix.cpp
    template <class T1, int C1, class T2>
    void InstCopy(
        const ConstUpperTriMatrixView<T1,C1>& m1, UpperTriMatrixView<T2> m2);
    template <class T1, int C1, class T2>
    void InstAliasCopy(
        const ConstUpperTriMatrixView<T1,C1>& m1, UpperTriMatrixView<T2> m2);

    //
    // Copy Matrices
    //
    
    // The maximum nops to unroll.
#if TMV_OPT >= 3
#define TMV_COPYU_UNROLL 200 
#elif TMV_OPT >= 2
#define TMV_COPYU_UNROLL 25
#elif TMV_OPT >= 1
#define TMV_COPYU_UNROLL 9
#else
#define TMV_COPYU_UNROLL 0
#endif

    template <int algo, int s, class M1, class M2>
    struct CopyU_Helper;

    // algo 0: s == 0, nothing to do
    template <class M1, class M2>
    struct CopyU_Helper<0,0,M1,M2>
    {
        static TMV_INLINE void call(const M1& , M2& ) {}
    };

    // algo 1: transpose 
    template <int s, class M1, class M2>
    struct CopyU_Helper<1,s,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            typedef typename M1::const_transpose_type M1t;
            typedef typename M2::transpose_type M2t;
            M1t m1t = m1.transpose();
            M2t m2t = m2.transpose();
            CopyU_Helper<-3,s,M1t,M2t>::call(m1t,m2t);
        }
    };

    // algo 2: m1 is unitdiag, break out the offdiag part
    template <int s, class M1, class M2>
    struct CopyU_Helper<2,s,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
            const int N = (s == UNKNOWN ? m2.size() : s);
            if (!m2.isunit()) m2.diag().setAllTo(1);
            if (N > 1) {
                typedef typename M1::const_offdiag_type M1o;
                typedef typename M2::offdiag_type M2o;
                M1o m1o = m1.offDiag();
                M2o m2o = m2.offDiag();
                const int sm1 = IntTraits2<s,-1>::sum;
                CopyU_Helper<-4,sm1,M1o,M2o>::call(m1o,m2o);
            }
        }
    };

    // algo 3: UnknownDiag
    template <int s, class M1, class M2>
    struct CopyU_Helper<3,s,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
            if (m1.isunit())
                CopyU_Helper<2,s,M1,M2>::call(m1,m2);
            else
                CopyU_Helper<-4,s,M1,M2>::call(m1,m2);
        }
    };

    // algo 11: Loop over columns
    template <int s, class M1, class M2>
    struct CopyU_Helper<11,s,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
            int N = (s == UNKNOWN ? m2.size() : s);
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::col_sub_type M2c;
            typedef typename M1c::const_iterator IT1;
            typedef typename M2c::iterator IT2;
            const int step1 = m1.stepj();
            const int step2 = m2.stepj();
            IT1 it1 = m1.get_col(0,0,1).begin();
            IT2 it2 = m2.get_col(0,0,1).begin();
            int M=1;
            for(;N;--N) {
                CopyV_Helper<-3,UNKNOWN,M1c,M2c>::call2(M++,it1,it2);
                it1.shiftP(step1);
                it2.shiftP(step2);
            } 
        }
    };

    // algo 15: Fully unroll by columns
    template <int s, class M1, class M2>
    struct CopyU_Helper<15,s,M1,M2>
    {
        template <int I, int M, int J, int N>
        struct Unroller
        {
            static inline void unroll(const M1& m1, M2& m2)
            {
                Unroller<I,M-(N-N/2),J,N/2>::unroll(m1,m2);
                Unroller<I,M,J+N/2,N-N/2>::unroll(m1,m2);
            }
        };
        template <int I, int M, int J>
        struct Unroller<I,M,J,1>
        {
            static inline void unroll(const M1& m1, M2& m2)
            {
                Unroller<I,M/2,J,1>::unroll(m1,m2);
                Unroller<I+M/2,M-M/2,J,1>::unroll(m1,m2);
            }
        };
        template <int I, int M, int J>
        struct Unroller<I,M,J,0>
        { static inline void unroll(const M1& , M2& ) {} };
        template <int I, int J>
        struct Unroller<I,1,J,1>
        {
            static inline void unroll(const M1& m1, M2& m2)
            { m2.ref(I,J) = m1.cref(I,J); }
        };
        template <int I, int J>
        struct Unroller<I,0,J,1>
        { static inline void unroll(const M1& , M2& ) {} };

        static inline void call(const M1& m1, M2& m2)
        { Unroller<0,s,0,s>::unroll(m1,m2); }
    };

    // algo 21: Loop over rows
    template <int s, class M1, class M2>
    struct CopyU_Helper<21,s,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            int N = (s == UNKNOWN ? m2.size() : s);
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M2::row_sub_type M2r;
            typedef typename M1r::const_iterator IT1;
            typedef typename M2r::iterator IT2;
            const int step1 = m1.diagstep();
            const int step2 = m2.diagstep();
            IT1 it1 = m1.get_row(0,0,N).begin();
            IT2 it2 = m2.get_row(0,0,N).begin();
            for(;N;--N) {
                CopyV_Helper<-3,UNKNOWN,M1r,M2r>::call2(N,it1,it2);
                it1.shiftP(step1);
                it2.shiftP(step2);
            }
        }
    };

    // algo 25: Fully unroll by rows
    template <int s, class M1, class M2>
    struct CopyU_Helper<25,s,M1,M2>
    {
        template <int I, int M, int J, int N>
        struct Unroller
        {
            static inline void unroll(const M1& m1, M2& m2)
            {
                Unroller<I,M/2,J,N>::unroll(m1,m2);
                Unroller<I+M/2,M-M/2,J+M/2,N-M/2>::unroll(m1,m2);
            }
        };
        template <int I, int J, int N>
        struct Unroller<I,1,J,N>
        {
            static inline void unroll(const M1& m1, M2& m2)
            {
                Unroller<I,1,J,N/2>::unroll(m1,m2);
                Unroller<I,1,J+N/2,N-N/2>::unroll(m1,m2);
            }
        };
        template <int I, int J, int N>
        struct Unroller<I,0,J,N>
        { static inline void unroll(const M1& , M2& ) {} };
        template <int I, int J>
        struct Unroller<I,1,J,1>
        {
            static inline void unroll(const M1& m1, M2& m2)
            { m2.ref(I,J) = m1.cref(I,J); }
        };
        template <int I, int J>
        struct Unroller<I,1,J,0>
        { static inline void unroll(const M1& , M2& ) {} };

        static inline void call(const M1& m1, M2& m2)
        { Unroller<0,s,0,s>::unroll(m1,m2); }
    };

    // algo 90: Call inst
    template <int size, class M1, class M2>
    struct CopyU_Helper<90,size,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
#ifndef TMV_INST_MIX
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            TMVStaticAssert((Traits2<T1,T2>::sametype));
#endif
            InstCopy(m1.xView(),m2.xView()); 
        }
    };

    // algo 91: Call inst alias
    template <int size, class M1, class M2>
    struct CopyU_Helper<91,size,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
#ifndef TMV_INST_MIX
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            TMVStaticAssert((Traits2<T1,T2>::sametype));
#endif
            InstAliasCopy(m1.xView(),m2.xView()); 
        }
    };

    // algo 96: Transpose
    template <int size, class M1, class M2>
    struct CopyU_Helper<96,size,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            typedef typename M1::const_transpose_type M1t;
            typedef typename M2::transpose_type M2t;
            M1t m1t = m1.transpose();
            M2t m2t = m2.transpose();
            CopyU_Helper<-2,size,M1t,M2t>::call(m1t,m2t);
        }
    };

    // algo 196: Transpose
    template <int size, class M1, class M2>
    struct CopyU_Helper<196,size,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            typedef typename M1::const_transpose_type M1t;
            typedef typename M2::transpose_type M2t;
            M1t m1t = m1.transpose();
            M2t m2t = m2.transpose();
            CopyU_Helper<99,size,M1t,M2t>::call(m1t,m2t);
        }
    };

    // algo 97: Conjugate
    template <int size, class M1, class M2>
    struct CopyU_Helper<97,size,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            CopyU_Helper<-2,size,M1c,M2c>::call(m1c,m2c);
        }
    };

    // algo 197: Conjugate
    template <int size, class M1, class M2>
    struct CopyU_Helper<197,size,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            CopyU_Helper<99,size,M1c,M2c>::call(m1c,m2c);
        }
    };

    // algo 98: Inline Check for aliases
    template <int size, class M1, class M2>
    struct CopyU_Helper<98,size,M1,M2>
    {
        // When we need a temporary below, there is a complication by
        // the possibility of m1 being UnknownDiag.  
        // In this case, m1.copy() becomes NonUnitDiag, which can't be
        // copied back to a UnitDiag m2.
        template <int algo, int dummy>
        struct CopyWithTemp;
        template <int dummy>
        struct CopyWithTemp<0,dummy> // Normal case.
        {
            static TMV_INLINE void call(const M1& m1, M2& m2)
            { NoAliasCopy(m1.copy(),m2); }
        };
        template <int dummy>
        struct CopyWithTemp<1,dummy> // Use UnitDiag
        {
            static TMV_INLINE void call(const M1& m1, M2& m2)
            {
                typename M2::unitdiag_type m2u = m2.viewAsUnitDiag();
                NoAliasCopy(m1.viewAsUnitDiag().copy(),m2u); 
            }
        };
        template <int dummy>
        struct CopyWithTemp<2,dummy> // Maybe use UnitDiag
        {
            static void call(const M1& m1, M2& m2)
            {
                if (m2.isunit()) {
                    typename M2::unitdiag_type m2u = m2.viewAsUnitDiag();
                    NoAliasCopy(m1.viewAsUnitDiag().copy(),m2u); 
                } else {
                    typename M2::nonunitdiag_type m2n = m2.viewAsUnitDiag();
                    NoAliasCopy(m1.copy(),m2n); 
                }
            }
        };
        static void call(const M1& m1, M2& m2)
        {
            TMVStaticAssert(M1::_upper == int(M2::_upper));
            TMVStaticAssert(M1::_unit || M1::_unknowndiag || !M2::_unit);
            TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
            TMVAssert(m1.size() == m2.size());
            TMVAssert(m1.isunit() || !m2.isunit());
            if ( !SameStorage(m1,m2) ||
                 OppositeStorage(m1,m2) ) {
                // No aliasing (or no clobbering)
                CopyU_Helper<-2,size,M1,M2>::call(m1,m2);
            } else if (ExactSameStorage(m1,m2)) {
                // They are already equal modulo a conjugation
                Maybe<M1::_conj != int(M2::_conj)>::conjself(m2);
                // Except possibly for the diagonal...
                if (m1.isunit() && !m2.isunit()) m2.diag().setAllTo(1);
            } else {
                // Need a temporary
                const int algo = 
                    M2::_unit ? 1 :
                    M1::_unknowndiag && M2::_unknowndiag ? 2 :
                    0;
                CopyWithTemp<algo,1>::call(m1,m2);
            }
        }
    };

    // algo 99: Check for aliases
    template <int s, class M1, class M2>
    struct CopyU_Helper<99,s,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            TMVStaticAssert(M1::_upper == int(M2::_upper));
            TMVStaticAssert(M1::_unit || M1::_unknowndiag || !M2::_unit);
            TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
            TMVAssert(m1.size() == m2.size());
            TMVAssert(m1.isunit() || !m2.isunit());
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst = 
                (s == UNKNOWN || s > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T2>::isinst;
            const int algo = 
                M2::_lower ? 196 :
                M2::_conj ? 197 :
                inst ? 91 :
                98;
            CopyU_Helper<algo,s,M1,M2>::call(m1,m2);
        }
    };

    // algo -4: No branches or copies
    // Also requires that M1,M2 are both NonUnitDiag and UpperTri
    template <int s, class M1, class M2>
    struct CopyU_Helper<-4,s,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            TMVStaticAssert(M1::_upper);
            TMVStaticAssert(M2::_upper);
            TMVStaticAssert(!M1::_unit && !M2::_unit);
            TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
            TMVAssert(m1.size() == m2.size());
            TMVAssert(!m1.isunit() && !m2.isunit());
            typedef typename M2::value_type T2;
            const int s2 = s > 20 ? UNKNOWN : s;
            const int s2p1 = IntTraits<s2>::Sp1;
            // nops = n(n+1)/2
            const int nops = IntTraits2<s2,s2p1>::safeprod / 2;
            const bool unroll = 
                s > 10 ? false :
                s == UNKNOWN ? false :
                nops <= TMV_COPYU_UNROLL;
            const int algo = 
                unroll ? ( M2::_rowmajor ? 25 : 15 ) :
                M2::_rowmajor ? 21 : 11;
            CopyU_Helper<algo,s,M1,M2>::call(m1,m2);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int s, class M1, class M2>
    struct CopyU_Helper<-3,s,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            TMVStaticAssert(M1::_upper == int(M2::_upper));
            TMVStaticAssert(M1::_unit || M1::_unknowndiag || !M2::_unit);
            TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
            TMVAssert(m1.size() == m2.size());
            TMVAssert(m1.isunit() || !m2.isunit());
            typedef typename M2::value_type T2;
            const int algo = 
                s == 0 ? 0 : 
                M2::_lower ? 1 :
                M1::_unit ? 2 :
                M1::_unknowndiag ? 3 :
                -4;
            CopyU_Helper<algo,s,M1,M2>::call(m1,m2);
        }
    };

    // algo -2: Check for inst
    template <int s, class M1, class M2>
    struct CopyU_Helper<-2,s,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            TMVStaticAssert(M1::_upper == int(M2::_upper));
            TMVStaticAssert(M1::_unit || M1::_unknowndiag || !M2::_unit);
            TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
            TMVAssert(m1.size() == m2.size());
            TMVAssert(m1.isunit() || !m2.isunit());
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst = 
                (s == UNKNOWN || s > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T2>::isinst;
            const int algo = 
                M2::_lower ? 96 :
                M2::_conj ? 97 :
                inst ? 90 :
                -3;
            CopyU_Helper<algo,s,M1,M2>::call(m1,m2);
        }
    };

    // algo -1: Check for aliases?
    template <int s, class M1, class M2>
    struct CopyU_Helper<-1,s,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            TMVStaticAssert(M1::_upper == int(M2::_upper));
            TMVStaticAssert(M1::_unit || M1::_unknowndiag || !M2::_unit);
            TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
            TMVAssert(m1.size() == m2.size());
            TMVAssert(m1.isunit() || !m2.isunit());
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool noclobber = MStepHelper<M1,M2>::opp;
            const bool checkalias = 
                M2::_checkalias && !noclobber;
            const int algo = 
                checkalias ? 99 : 
                -2;
            CopyU_Helper<algo,s,M1,M2>::call(m1,m2);
        }
    };

    template <class M1, class M2>
    static inline void Copy(
        const BaseMatrix_Tri<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
    {
        TMVStaticAssert(M1::_upper == int(M2::_upper));
        TMVStaticAssert(M1::_unit || M1::_unknowndiag || !M2::_unit);
        TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.isunit() || !m2.isunit());
        const int s = Sizes<M1::_size,M2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        CopyU_Helper<-1,s,M1v,M2v>::call(m1v,m2v);
    }

    template <class M1, class M2>
    static inline void NoAliasCopy(
        const BaseMatrix_Tri<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
    {
        TMVStaticAssert(M1::_upper == int(M2::_upper));
        TMVStaticAssert(M1::_unit || M1::_unknowndiag || !M2::_unit);
        TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.isunit() || !m2.isunit());
        const int s = Sizes<M1::_size,M2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        CopyU_Helper<-2,s,M1v,M2v>::call(m1v,m2v);
    }

    template <class M1, class M2>
    static inline void InlineCopy(
        const BaseMatrix_Tri<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
    {
        TMVStaticAssert(M1::_upper == int(M2::_upper));
        TMVStaticAssert(M1::_unit || M1::_unknowndiag || !M2::_unit);
        TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.isunit() || !m2.isunit());
        const int s = Sizes<M1::_size,M2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        CopyU_Helper<-3,s,M1v,M2v>::call(m1v,m2v);
    }

    template <class M1, class M2>
    static inline void InlineAliasCopy(
        const BaseMatrix_Tri<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
    {
        TMVStaticAssert(M1::_upper == int(M2::_upper));
        TMVStaticAssert(M1::_unit || M1::_unknowndiag || !M2::_unit);
        TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.isunit() || !m2.isunit());
        const int s = Sizes<M1::_size,M2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        CopyU_Helper<98,s,M1v,M2v>::call(m1v,m2v);
    }

    template <class M1, class M2>
    static inline void AliasCopy(
        const BaseMatrix_Tri<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
    {
        TMVStaticAssert(M1::_upper == int(M2::_upper));
        TMVStaticAssert(M1::_unit || M1::_unknowndiag || !M2::_unit);
        TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.isunit() || !m2.isunit());
        const int s = Sizes<M1::_size,M2::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        CopyU_Helper<99,s,M1v,M2v>::call(m1v,m2v);
    }


    // 
    // M = U
    //

    template <class M1, class M2>
    static inline void Copy(
        const BaseMatrix_Tri<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_size,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_size,M2::_rowsize>::same));
        TMVAssert(m1.size() == m2.colsize());
        TMVAssert(m1.size() == m2.rowsize());
        const bool upper = M1::_upper;
        typedef typename TypeSelect<upper,
                typename M2::uppertri_type,
                typename M2::lowertri_type>::type M2u;
        M2u m2u = Maybe<upper>::uppertri(m2);
        Copy(m1,m2u);
        if (m1.size() > 1) 
            Maybe<!upper>::uppertri(m2).offDiag().setZero();
    }

    template <class M1, class M2>
    static inline void NoAliasCopy(
        const BaseMatrix_Tri<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_size,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_size,M2::_rowsize>::same));
        TMVAssert(m1.size() == m2.colsize());
        TMVAssert(m1.size() == m2.rowsize());
        const bool upper = M1::_upper;
        typedef typename TypeSelect<upper,
                typename M2::uppertri_type,
                typename M2::lowertri_type>::type M2u;
        M2u m2u = Maybe<upper>::uppertri(m2);
        m2.setZero();
        NoAliasCopy(m1,m2u);
    }

    template <class M1, class M2>
    static inline void AliasCopy(
        const BaseMatrix_Tri<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_size,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_size,M2::_rowsize>::same));
        TMVAssert(m1.size() == m2.colsize());
        TMVAssert(m1.size() == m2.rowsize());
        const bool upper = M1::_upper;
        typedef typename TypeSelect<upper,
                typename M2::uppertri_type,
                typename M2::lowertri_type>::type M2u;
        M2u m2u = Maybe<upper>::uppertri(m2);
        AliasCopy(m1,m2u);
        if (m1.size() > 1) 
            Maybe<!upper>::uppertri(m2).offDiag().setZero();
    }


    //
    // U = D
    //

    template <class M1, class M2>
    static inline void Copy(
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
        TMVAssert(m1.size() == m2.size());
        typename M1::const_diag_type m1d = m1.diag();
        typename M2::diag_type m2d = m2.diag();
        Copy(m1d,m2d);
        if (m2.size() > 1) m2.offDiag().setZero();
    }

    template <class M1, class M2>
    static inline void NoAliasCopy(
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
        TMVAssert(m1.size() == m2.size());
        typename M1::const_diag_type m1d = m1.diag();
        typename M2::diag_type m2d = m2.diag();
        m2.setZero();
        NoAliasCopy(m1d,m2d);
    }

    template <class M1, class M2>
    static TMV_INLINE void AliasCopy(
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
    { Copy(m1,m2); }


    //
    // L = L
    //

    template <class T1, int C1, class T2>
    static TMV_INLINE void InstCopy(
        const ConstLowerTriMatrixView<T1,C1>& m1, LowerTriMatrixView<T2> m2)
    { InstCopy(m1.transpose(),m2.transpose()); }


#undef TMV_COPYU_UNROLL

} // namespace tmv

#endif
