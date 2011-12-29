
#ifndef TMV_SwapM_H
#define TMV_SwapM_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_SwapV.h"

namespace tmv {

    // Defined in TMV_Matrix.cpp
    template <class T, int C>
    void InstSwap(MatrixView<T,C> m1, MatrixView<T> m2); 
    template <class T, int C>
    void InstAliasSwap(MatrixView<T,C> m1, MatrixView<T> m2); 

    //
    // Swap Matrices
    //

    template <int algo, int cs, int rs, class M1, class M2>
    struct SwapM_Helper;

    // algo 1: Linearize to vector version
    template <int cs, int rs, class M1, class M2>
    struct SwapM_Helper<1,cs,rs,M1,M2>
    {
        static inline void call(M1& m1, M2& m2)
        {
            typedef typename M1::linearview_type M1l;
            typedef typename M2::linearview_type M2l;
            M1l m1l = m1.linearView();
            M2l m2l = m2.linearView();
            const int cs_rs = IntTraits2<cs,rs>::prod;
            SwapV_Helper<-3,cs_rs,M1l,M2l>::call(m1l,m2l);
        }
    };

    // algo 11: Loop over columns
    template <int cs, int rs, class M1, class M2>
    struct SwapM_Helper<11,cs,rs,M1,M2>
    {
        static void call(M1& m1, M2& m2)
        {
            const int M = cs == TMV_UNKNOWN ? m2.colsize() : cs;
            int N = rs == TMV_UNKNOWN ? m2.rowsize() : rs;
            typedef typename M1::col_type M1c;
            typedef typename M2::col_type M2c;
            typedef typename M1c::iterator IT1;
            typedef typename M2c::iterator IT2;
            const int step1 = m1.stepj();
            const int step2 = m2.stepj();
            IT1 it1 = m1.get_col(0).begin();
            IT2 it2 = m2.get_col(0).begin();
            for(;N;--N) {
                SwapV_Helper<-3,cs,M1c,M2c>::call2(M,it1,it2);
                it1.shiftP(step1);
                it2.shiftP(step2);
            } 
        }
    };

    // algo 12: Loop over rows
    template <int cs, int rs, class M1, class M2>
    struct SwapM_Helper<12,cs,rs,M1,M2>
    {
        static void call(M1& m1, M2& m2)
        {
            int M = cs == TMV_UNKNOWN ? m2.colsize() : cs;
            const int N = rs == TMV_UNKNOWN ? m2.rowsize() : rs;
            typedef typename M1::row_type M1r;
            typedef typename M2::row_type M2r;
            typedef typename M1r::iterator IT1;
            typedef typename M2r::iterator IT2;
            const int step1 = m1.stepi();
            const int step2 = m2.stepi();
            IT1 it1 = m1.get_row(0).begin();
            IT2 it2 = m2.get_row(0).begin();
            for(;M;--M) {
                SwapV_Helper<-3,rs,M1r,M2r>::call2(N,it1,it2);
                it1.shiftP(step1);
                it2.shiftP(step2);
            } 
        }
    };

    // algo 30: Unknown sizes, determine which algorithm to use
    template <int cs, int rs, class M1, class M2>
    struct SwapM_Helper<30,cs,rs,M1,M2>
    {
        static void call(M1& m1, M2& m2)
        {
            if (m1.canLinearize() && m2.canLinearize() &&
                m1.stepi() == m2.stepi() && m1.stepj() == m2.stepj()) 
                SwapM_Helper<1,cs,rs,M1,M2>::call(m1,m2);
            else if ( ( m1.isrm() && m2.isrm() ) || 
                      ( !(m1.iscm() && m2.iscm()) &&
                        (m1.colsize() > m1.rowsize()) ) )
                SwapM_Helper<12,cs,rs,M1,M2>::call(m1,m2);
            else 
                SwapM_Helper<11,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo 15: Fully unroll by columns
    template <int cs, int rs, class M1, class M2>
    struct SwapM_Helper<15,cs,rs,M1,M2>
    {
        template <int I, int M, int J, int N>
        struct Unroller
        {
            static TMV_INLINE void unroll(M1& m1, M2& m2)
            {
                Unroller<I,M,J,N/2>::unroll(m1,m2);
                Unroller<I,M,J+N/2,N-N/2>::unroll(m1,m2);
            }
        };
        template <int I, int M, int J>
        struct Unroller<I,M,J,1>
        {
            static TMV_INLINE void unroll(M1& m1, M2& m2)
            {
                Unroller<I,M/2,J,1>::unroll(m1,m2);
                Unroller<I+M/2,M-M/2,J,1>::unroll(m1,m2);
            }
        };
        template <int I, int M, int J>
        struct Unroller<I,M,J,0>
        { static TMV_INLINE void unroll(M1& , M2& ) {} };
        template <int I, int J>
        struct Unroller<I,1,J,1>
        {
            static TMV_INLINE void unroll(M1& m1, M2& m2)
            { TMV_SWAP(m2.ref(I,J) , m1.ref(I,J)); }
        };
        template <int I, int J>
        struct Unroller<I,0,J,1>
        { static TMV_INLINE void unroll(M1& , M2& ) {} };

        static inline void call(M1& m1, M2& m2)
        { Unroller<0,cs,0,rs>::unroll(m1,m2); }
    };

    // algo 16: Fully unroll by rows
    template <int cs, int rs, class M1, class M2>
    struct SwapM_Helper<16,cs,rs,M1,M2>
    {
        template <int I, int M, int J, int N>
        struct Unroller
        {
            static TMV_INLINE void unroll(M1& m1, M2& m2)
            {
                Unroller<I,M/2,J,N>::unroll(m1,m2);
                Unroller<I+M/2,M-M/2,J,N>::unroll(m1,m2);
            }
        };
        template <int I, int J, int N>
        struct Unroller<I,1,J,N>
        {
            static TMV_INLINE void unroll(M1& m1, M2& m2)
            {
                Unroller<I,1,J,N/2>::unroll(m1,m2);
                Unroller<I,1,J+N/2,N-N/2>::unroll(m1,m2);
            }
        };
        template <int I, int J, int N>
        struct Unroller<I,0,J,N>
        { static TMV_INLINE void unroll(M1& , M2& ) {} };
        template <int I, int J>
        struct Unroller<I,1,J,1>
        {
            static TMV_INLINE void unroll(M1& m1, M2& m2)
            { TMV_SWAP(m2.ref(I,J) , m1.ref(I,J) ); }
        };
        template <int I, int J>
        struct Unroller<I,1,J,0>
        { static TMV_INLINE void unroll(M1& , M2& ) {} };

        static inline void call(M1& m1, M2& m2)
        { Unroller<0,cs,0,rs>::unroll(m1,m2); }
    };

    // algo 90: Call inst
    template <int cs, int rs, class M1, class M2>
    struct SwapM_Helper<90,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        { InstSwap(m1.xView(),m2.xView()); }
    };

    // algo 91: Call inst alias
    template <int cs, int rs, class M1, class M2>
    struct SwapM_Helper<91,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        { InstAliasSwap(m1.xView(),m2.xView()); }
    };

    // algo 97: Conjugate
    template <int cs, int rs, class M1, class M2>
    struct SwapM_Helper<97,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        {
            typedef typename M1::conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            SwapM_Helper<-2,cs,rs,M1c,M2c>::call(m1c,m2c);
        }
    };

    // algo 197: Conjugate
    template <int cs, int rs, class M1, class M2>
    struct SwapM_Helper<197,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        {
            typedef typename M1::conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            SwapM_Helper<99,cs,rs,M1c,M2c>::call(m1c,m2c);
        }
    };

    // algo 98: Inline check for aliases
    template <int cs, int rs, class M1, class M2>
    struct SwapM_Helper<98,cs,rs,M1,M2>
    {
        static void call(M1& m1, M2& m2)
        {
            if (!SameStorage(m1,m2)) {
                // No aliasing
                SwapM_Helper<-2,cs,rs,M1,M2>::call(m1,m2);
            } else if (ExactSameStorage(m1,m2)) {
                // They are equal modulo a conjugation
                Maybe<M1::_conj != int(M2::_conj)>::conjself(m1);
            } else if (m2.colsize() == m2.rowsize() &&
                       OppositeStorage(m1,m2)) {
                // Then transpose
                m2.transposeSelf();
                // And maybe conjugate
                Maybe<M1::_conj != int(M2::_conj)>::conjself(m2);
            } else {
                // Need a temporary
                typename M1::copy_type m1c = m1;
                m1.noAlias() = m2;
                m2.noAlias() = m1c;
            }
        }
    };

    // algo 99: Check for aliases
    template <int cs, int rs, class M1, class M2>
    struct SwapM_Helper<99,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst = 
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
                Traits2<T1,T2>::sametype &&
                Traits<T1>::isinst;
            const int algo =
                M2::_conj ? 197 :
                inst ? 91 :
                98;
            SwapM_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, class M1, class M2>
    struct SwapM_Helper<-3,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        {
            typedef typename M2::value_type T2;
            const bool canlin = 
                M1::_canlin && M2::_canlin &&
                ( (M1::_rowmajor && M2::_rowmajor) ||
                  (M1::_colmajor && M2::_colmajor) );
            const int algo = 
                canlin ? 1 :
                TMV_OPT == 0 ? (( M1::_rowmajor && M2::_rowmajor ) ? 2 : 11 ) :
                ( cs != TMV_UNKNOWN && rs != TMV_UNKNOWN ) ? (
                    ( IntTraits2<cs,rs>::prod <= int(128/sizeof(T2)) ) ? (
                        ( M1::_rowmajor && M2::_rowmajor ) ? 16 : 15 ) :
                    ( M1::_rowmajor && M2::_rowmajor ) ? 12 :
                    ( M1::_colmajor && M2::_colmajor ) ? 11 :
                    ( cs > rs ) ? 12 : 11 ) :
                TMV_OPT >= 2 ? 30 :
                ( M1::_rowmajor && M2::_rowmajor ) ? 12 :
                ( M1::_colmajor && M2::_colmajor ) ? 11 :
                ( cs == TMV_UNKNOWN || rs == TMV_UNKNOWN ) ? 11 :
                ( cs > rs ) ? 12 : 
                11;
            SwapM_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, class M1, class M2>
    struct SwapM_Helper<-2,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst = 
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
                Traits2<T1,T2>::sametype &&
                Traits<T1>::isinst;
            const int algo =
                M2::_conj ? 97 :
                inst ? 90 :
                -3;
            SwapM_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo -1: Check for aliases?
    template <int cs, int rs, class M1, class M2>
    struct SwapM_Helper<-1,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        {
            const bool checkalias = 
                M1::_checkalias || M2::_checkalias;
            const int algo =
                checkalias ? 99 : 
                -2;
            SwapM_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    template <class M1, class M2>
    inline void DoSwap(
        BaseMatrix_Rec_Mutable<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        const int cs = Sizes<M1::_colsize,M2::_colsize>::size;
        const int rs = Sizes<M1::_rowsize,M2::_rowsize>::size;
        typedef typename M1::cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_REF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        SwapM_Helper<-1,cs,rs,M1v,M2v>::call(m1v,m2v);
    }

    template <class M1, class M2>
    inline void InlineSwap(
        BaseMatrix_Rec_Mutable<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        const int cs = Sizes<M1::_colsize,M2::_colsize>::size;
        const int rs = Sizes<M1::_rowsize,M2::_rowsize>::size;
        typedef typename M1::cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_REF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        SwapM_Helper<-3,cs,rs,M1v,M2v>::call(m1v,m2v);
    }

    template <class M1, class M2>
    inline void InlineAliasSwap(
        BaseMatrix_Rec_Mutable<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        const int cs = Sizes<M1::_colsize,M2::_colsize>::size;
        const int rs = Sizes<M1::_rowsize,M2::_rowsize>::size;
        typedef typename M1::cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_REF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        SwapM_Helper<98,cs,rs,M1v,M2v>::call(m1v,m2v);
    }

    template <class M1, class M2>
    TMV_INLINE void Swap(
        BaseMatrix_Rec_Mutable<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
    { DoSwap(m1,m2); }

} // namespace tmv

#endif
