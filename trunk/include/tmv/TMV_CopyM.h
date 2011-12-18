
#ifndef TMV_CopyM_H
#define TMV_CopyM_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_CopyV.h"

namespace tmv {

    // Defined in TMV_Matrix.cpp
    template <class T1, int C1, class T2>
    void InstCopy(const ConstMatrixView<T1,C1>& m1, MatrixView<T2> m2); 
    template <class T1, int C1, class T2>
    void InstAliasCopy(const ConstMatrixView<T1,C1>& m1, MatrixView<T2> m2); 

    //
    // Copy Matrices
    //

    template <int algo, int cs, int rs, class M1, class M2>
    struct CopyM_Helper;

    // algo 0: size = 0, nothing to do
    template <int cs, int rs, class M1, class M2>
    struct CopyM_Helper<0,cs,rs,M1,M2>
    { static TMV_INLINE void call(const M1& , M2& ) {} };

    // algo 1: Linearize to vector version
    template <int cs, int rs, class M1, class M2>
    struct CopyM_Helper<1,cs,rs,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            typedef typename M1::const_linearview_type M1l;
            typedef typename M2::linearview_type M2l;
            M1l m1l = m1.linearView();
            M2l m2l = m2.linearView();
            const int cs_rs = IntTraits2<cs,rs>::prod;
            CopyV_Helper<-3,cs_rs,M1l,M2l>::call(m1l,m2l);
        }
    };

    // algo 2: Only one column
    template <int cs, int rs, class M1, class M2>
    struct CopyM_Helper<2,cs,rs,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            typedef typename M1::const_col_type M1c;
            typedef typename M2::col_type M2c;
            M1c m1c = m1.col(0);
            M2c m2c = m2.col(0);
            CopyV_Helper<-3,cs,M1c,M2c>::call(m1c,m2c);
        }
    };

    // algo 102: Only one column, use algo -1
    template <int cs, int rs, class M1, class M2>
    struct CopyM_Helper<102,cs,rs,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            typedef typename M1::const_col_type M1c;
            typedef typename M2::col_type M2c;
            M1c m1c = m1.col(0);
            M2c m2c = m2.col(0);
            CopyV_Helper<-1,cs,M1c,M2c>::call(m1c,m2c);
        }
    };

    // algo 202: Only one column, use algo -2
    template <int cs, int rs, class M1, class M2>
    struct CopyM_Helper<202,cs,rs,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            typedef typename M1::const_col_type M1c;
            typedef typename M2::col_type M2c;
            M1c m1c = m1.col(0);
            M2c m2c = m2.col(0);
            CopyV_Helper<-2,cs,M1c,M2c>::call(m1c,m2c);
        }
    };

    // algo 3: Only one row
    template <int cs, int rs, class M1, class M2>
    struct CopyM_Helper<3,cs,rs,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            typedef typename M1::const_row_type M1r;
            typedef typename M2::row_type M2r;
            M1r m1r = m1.row(0);
            M2r m2r = m2.row(0);
            CopyV_Helper<-3,rs,M1r,M2r>::call(m1r,m2r);
        }
    };

    // algo 103: Only one row, use algo -1
    template <int cs, int rs, class M1, class M2>
    struct CopyM_Helper<103,cs,rs,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            typedef typename M1::const_row_type M1r;
            typedef typename M2::row_type M2r;
            M1r m1r = m1.row(0);
            M2r m2r = m2.row(0);
            CopyV_Helper<-1,rs,M1r,M2r>::call(m1r,m2r);
        }
    };

    // algo 203: Only one row, use algo -2
    template <int cs, int rs, class M1, class M2>
    struct CopyM_Helper<203,cs,rs,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            typedef typename M1::const_row_type M1r;
            typedef typename M2::row_type M2r;
            M1r m1r = m1.row(0);
            M2r m2r = m2.row(0);
            CopyV_Helper<-2,rs,M1r,M2r>::call(m1r,m2r);
        }
    };

    // algo 11: Loop over columns
    template <int cs, int rs, class M1, class M2>
    struct CopyM_Helper<11,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
            const int M = cs == TMV_UNKNOWN ? int(m2.colsize()) : cs;
            int N = rs == TMV_UNKNOWN ? int(m2.rowsize()) : rs;
            typedef typename M1::const_col_type M1c;
            typedef typename M2::col_type M2c;
            typedef typename M1c::const_iterator IT1;
            typedef typename M2c::iterator IT2;
            const int step1 = m1.stepj();
            const int step2 = m2.stepj();
            IT1 it1 = m1.get_col(0).begin();
            IT2 it2 = m2.get_col(0).begin();
            for(;N;--N) {
                CopyV_Helper<-3,cs,M1c,M2c>::call2(M,it1,it2);
                it1.shiftP(step1);
                it2.shiftP(step2);
            }
        }
    };

    // algo 15: Fully unroll by columns
    template <int cs, int rs, class M1, class M2>
    struct CopyM_Helper<15,cs,rs,M1,M2>
    {
        template <int I, int M, int J, int N>
        struct Unroller
        {
            static inline void unroll(const M1& m1, M2& m2)
            {
                Unroller<I,M,J,N/2>::unroll(m1,m2);
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
        { Unroller<0,cs,0,rs>::unroll(m1,m2); }
    };

    // algo 21: Loop over rows
    template <int cs, int rs, class M1, class M2>
    struct CopyM_Helper<21,cs,rs,M1,M2>
    {
        static inline void call(const M1& m1, M2& m2)
        {
            int M = cs == TMV_UNKNOWN ? int(m2.colsize()) : cs;
            const int N = rs == TMV_UNKNOWN ? int(m2.rowsize()) : rs;
            typedef typename M1::const_row_type M1r;
            typedef typename M2::row_type M2r;
            typedef typename M1r::const_iterator IT1;
            typedef typename M2r::iterator IT2;
            const int step1 = m1.stepi();
            const int step2 = m2.stepi();
            IT1 it1 = m1.get_row(0).begin();
            IT2 it2 = m2.get_row(0).begin();
            for(;M;--M) {
                CopyV_Helper<-3,rs,M1r,M2r>::call2(N,it1,it2);
                it1.shiftP(step1);
                it2.shiftP(step2);
            }
        }
    };

    // algo 25: Fully unroll by rows
    template <int cs, int rs, class M1, class M2>
    struct CopyM_Helper<25,cs,rs,M1,M2>
    {
        template <int I, int M, int J, int N>
        struct Unroller
        {
            static TMV_INLINE void unroll(const M1& m1, M2& m2)
            {
                Unroller<I,M/2,J,N>::unroll(m1,m2);
                Unroller<I+M/2,M-M/2,J,N>::unroll(m1,m2);
            }
        };
        template <int I, int J, int N>
        struct Unroller<I,1,J,N>
        {
            static TMV_INLINE void unroll(const M1& m1, M2& m2)
            {
                Unroller<I,1,J,N/2>::unroll(m1,m2);
                Unroller<I,1,J+N/2,N-N/2>::unroll(m1,m2);
            }
        };
        template <int I, int J, int N>
        struct Unroller<I,0,J,N>
        { static TMV_INLINE void unroll(const M1& , M2& ) {} };
        template <int I, int J>
        struct Unroller<I,1,J,1>
        {
            static TMV_INLINE void unroll(const M1& m1, M2& m2)
            { m2.ref(I,J) = m1.cref(I,J); }
        };
        template <int I, int J>
        struct Unroller<I,1,J,0>
        { static TMV_INLINE void unroll(const M1& , M2& ) {} };

        static inline void call(const M1& m1, M2& m2)
        { Unroller<0,cs,0,rs>::unroll(m1,m2); }
    };

    // algo 31: Unknown sizes, determine which algorithm to use
    template <int cs, int rs, class M1, class M2>
    struct CopyM_Helper<31,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
            if (m1.canLinearize() && m2.canLinearize() &&
                m1.stepi() == m2.stepi() && m1.stepj() == m2.stepj()) 
                CopyM_Helper<1,cs,rs,M1,M2>::call(m1,m2);
            else if ((m1.isrm() && m2.isrm()) || 
                     ( !(m1.iscm() && m2.iscm()) &&
                       (m1.colsize() > m1.rowsize()) ) )
                CopyM_Helper<21,cs,rs,M1,M2>::call(m1,m2);
            else 
                CopyM_Helper<11,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo 90: Call inst
    template <int cs, int rs, class M1, class M2>
    struct CopyM_Helper<90,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        { InstCopy(m1.xView(),m2.xView()); }
    };

    // algo 91: Call inst alias
    template <int cs, int rs, class M1, class M2>
    struct CopyM_Helper<91,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        { InstAliasCopy(m1.xView(),m2.xView()); }
    };

    // algo 97: Conjugate
    template <int cs, int rs, class M1, class M2>
    struct CopyM_Helper<97,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            CopyM_Helper<-2,cs,rs,M1c,M2c>::call(m1c,m2c);
        }
    };

    // algo 197: Conjugate
    template <int cs, int rs, class M1, class M2>
    struct CopyM_Helper<197,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            CopyM_Helper<99,cs,rs,M1c,M2c>::call(m1c,m2c);
        }
    };

    // algo 98: Inline check for aliases
    template <int cs, int rs, class M1, class M2>
    struct CopyM_Helper<98,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
            if (!SameStorage(m1,m2)) {
                // No aliasing
                CopyM_Helper<-2,cs,rs,M1,M2>::call(m1,m2);
            } else if (ExactSameStorage(m1,m2)) {
                // They are already equal modulo a conjugation.
                Maybe<M1::_conj != int(M2::_conj)>::conjself(m2);
            } else if (m2.colsize() == m2.rowsize() &&
                       OppositeStorage(m1,m2)) {
                // Then transpose
                const bool sq = cs == rs && cs != TMV_UNKNOWN;
                Maybe<sq>::tranself(m2);
                // And maybe conjugate
                Maybe<sq && M1::_conj != int(M2::_conj)>::conjself(m2);
            } else if (m2.ptr()) {
                // Need a temporary
                m2.noAlias() = m1.copy();
            } else {
                // else m2.ptr == 0, so don't need to do anything.
                // m1 and m2 are degenerate
                TMVAssert(m1.cptr() == 0);
                TMVAssert(m2.cptr() == 0);
                TMVAssert(m2.colsize() == 0 || m2.rowsize() == 0);
            }
        }
    };

    // algo 99: Check for aliases
    template <int cs, int rs, class M1, class M2>
    struct CopyM_Helper<99,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst = 
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T1>::isinst;
            const int algo = 
                M2::_conj ? 197 :
                inst ? 91 :
                98;
            CopyM_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo -4: No branches or copies
    template <int cs, int rs, class M1, class M2>
    struct CopyM_Helper<-4,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            typedef typename M2::value_type T2;
            const bool allrm = M1::_rowmajor && M2::_rowmajor;
            const bool allcm = M1::_colmajor && M2::_colmajor;
            const bool canlin = 
                M1::_canlin && M2::_canlin && (allrm || allcm);
            const int algo = 
                (cs == 0 || rs == 0) ? 0 :
                canlin ? 1 :
                rs == 1 ? 2 : 
                cs == 1 ? 3 :
                TMV_OPT == 0 ? (allrm ? 21 : 11) :
                ( cs != TMV_UNKNOWN && rs != TMV_UNKNOWN ) ? (
                    ( IntTraits2<cs,rs>::prod <= int(128/sizeof(T2)) ) ? (
                        ( M1::_rowmajor && M2::_rowmajor ) ? 25 : 15 ) :
                    allrm ? 21 : 
                    allcm ? 11 :
                    ( cs > rs ) ? 21 : 11 ) :
                allrm ? 21 : 11;
            CopyM_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, class M1, class M2>
    struct CopyM_Helper<-3,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            const int algo = 
                rs == 1 ? 2 : 
                cs == 1 ? 3 :
                TMV_OPT <= 1 ? -4 : 
                cs == TMV_UNKNOWN || rs == TMV_UNKNOWN ? 31 :
                -4;
            CopyM_Helper<algo,cs,rs,M1,M2>::call(m1,m2); 
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, class M1, class M2>
    struct CopyM_Helper<-2,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst = 
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T1>::isinst;
            const int algo = 
                rs == 1 ? 202 : 
                cs == 1 ? 203 :
                M2::_conj ? 97 :
                inst ? 90 :
                -3;
            CopyM_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo -1: Check for aliases?
    template <int cs, int rs, class M1, class M2>
    struct CopyM_Helper<-1,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const int algo = 
                rs == 1 ? 102 : 
                cs == 1 ? 103 :
                M2::_checkalias ? 99 : 
                -2;
            CopyM_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    template <class M1, class M2>
    inline void Copy(
        const BaseMatrix_Rec<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        const int cs = Sizes<M1::_colsize,M2::_colsize>::size;
        const int rs = Sizes<M1::_rowsize,M2::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        CopyM_Helper<-1,cs,rs,M1v,M2v>::call(m1v,m2v);
    }

    template <class M1, class M2>
    inline void InlineCopy(
        const BaseMatrix_Rec<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        const int cs = Sizes<M1::_colsize,M2::_colsize>::size;
        const int rs = Sizes<M1::_rowsize,M2::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        CopyM_Helper<-3,cs,rs,M1v,M2v>::call(m1v,m2v);
    }

    template <class M1, class M2>
    inline void InlineAliasCopy(
        const BaseMatrix_Rec<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        const int cs = Sizes<M1::_colsize,M2::_colsize>::size;
        const int rs = Sizes<M1::_rowsize,M2::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        CopyM_Helper<98,cs,rs,M1v,M2v>::call(m1v,m2v);
    }

} // namespace tmv

#endif
