
#ifndef TMV_SwapM_H
#define TMV_SwapM_H

#include "TMV_BaseMatrix_Band.h"
#include "TMV_SwapV.h"

namespace tmv {

    // Defined in TMV_BandMatrix.cpp
    template <class T, int C>
    void InstSwap(BandMatrixView<T,C> m1, BandMatrixView<T> m2); 
    template <class T, int C>
    void InstAliasSwap(BandMatrixView<T,C> m1, BandMatrixView<T> m2); 

    //
    // Swap Matrices
    //

    template <int algo, ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct SwapB_Helper;

    // algo 11: Loop over columns
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct SwapB_Helper<11,cs,rs,M1,M2>
    {
        static void call(M1& m1, M2& m2)
        {
            const ptrdiff_t M = cs == Unknown ? m2.colsize() : cs;
            const ptrdiff_t N = rs == Unknown ? m2.rowsize() : rs;
            const ptrdiff_t xx = Unknown;
            typedef typename M1::col_sub_type M1c;
            typedef typename M2::col_sub_type M2c;
            typedef typename M1c::iterator IT1;
            typedef typename M2c::iterator IT2;
            const ptrdiff_t rowstep1 = m1.stepj();
            const ptrdiff_t rowstep2 = m2.stepj();
            const ptrdiff_t diagstep1 = m1.diagstep();
            const ptrdiff_t diagstep2 = m2.diagstep();

            const ptrdiff_t lh = IntTraits<IntTraits2<M1::_nlo,M1::_nhi>::sum>::Sp1;
            const ptrdiff_t j1 = m1.nhi();
            const ptrdiff_t j2 = TMV_MIN(N,M-m1.nlo());
            const ptrdiff_t j3 = TMV_MIN(N,M+m1.nhi());
            ptrdiff_t len = m1.nlo()+1;
            IT1 it1 = m1.get_col(0,0,len).begin();
            IT2 it2 = m2.get_col(0,0,len).begin();
            ptrdiff_t j=0;
            for(;j<j1;++j) {
                SwapV_Helper<-3,xx,M1c,M2c>::call2(len,it1,it2);
                it1.shiftP(rowstep1);
                it2.shiftP(rowstep2);
                if (len < M) ++len;
            }
            for(;j<j2;++j) {
                SwapV_Helper<-3,lh,M1c,M2c>::call2(len,it1,it2);
                it1.shiftP(diagstep1);
                it2.shiftP(diagstep2);
            }
            if (j1 >= j2) ++len;
            for(;j<j3;++j) {
                SwapV_Helper<-3,xx,M1c,M2c>::call2(--len,it1,it2);
                it1.shiftP(diagstep1);
                it2.shiftP(diagstep2);
            }
        }
    };

    // algo 12: Loop over rows
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct SwapB_Helper<12,cs,rs,M1,M2>
    {
        static void call(M1& m1, M2& m2)
        {
            const ptrdiff_t M = cs == Unknown ? m2.colsize() : cs;
            const ptrdiff_t N = rs == Unknown ? m2.rowsize() : rs;
            const ptrdiff_t xx = Unknown;
            typedef typename M1::row_sub_type M1r;
            typedef typename M2::row_sub_type M2r;
            typedef typename M1r::iterator IT1;
            typedef typename M2r::iterator IT2;
            const ptrdiff_t colstep1 = m1.stepi();
            const ptrdiff_t colstep2 = m2.stepi();
            const ptrdiff_t diagstep1 = m1.diagstep();
            const ptrdiff_t diagstep2 = m2.diagstep();

            const ptrdiff_t lh = IntTraits<IntTraits2<M1::_nlo,M1::_nhi>::sum>::Sp1;
            const ptrdiff_t i1 = m1.nlo();
            const ptrdiff_t i2 = TMV_MIN(M,N-m1.nhi());
            const ptrdiff_t i3 = TMV_MIN(M,N+m1.nlo());
            ptrdiff_t len = m1.nhi()+1;
            IT1 it1 = m1.get_row(0,0,len).begin();
            IT2 it2 = m2.get_row(0,0,len).begin();
            ptrdiff_t i=0;
            for(;i<i1;++i) {
                SwapV_Helper<-3,xx,M1r,M2r>::call2(len,it1,it2);
                it1.shiftP(colstep1);
                it2.shiftP(colstep2);
                if (len < N) ++len;
            }
            for(;i<i2;++i) {
                SwapV_Helper<-3,lh,M1r,M2r>::call2(len,it1,it2);
                it1.shiftP(diagstep1);
                it2.shiftP(diagstep2);
            }
            if (i1 >= i2) ++len;
            for(;i<i3;++i) {
                SwapV_Helper<-3,xx,M1r,M2r>::call2(--len,it1,it2);
                it1.shiftP(diagstep1);
                it2.shiftP(diagstep2);
            }
        }
    };

    // algo 13: Loop over diagonals
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct SwapB_Helper<13,cs,rs,M1,M2>
    {
        static void call(M1& m1, M2& m2)
        {
            const ptrdiff_t M = cs == Unknown ? m2.colsize() : cs;
            const ptrdiff_t N = rs == Unknown ? m2.rowsize() : rs;
            const ptrdiff_t xx = Unknown;
            typedef typename M1::diag_sub_type M1d;
            typedef typename M2::diag_sub_type M2d;
            typedef typename M1d::iterator IT1;
            typedef typename M2d::iterator IT2;
            const ptrdiff_t colstep1 = m1.stepi();
            const ptrdiff_t colstep2 = m2.stepi();
            const ptrdiff_t rowstep1 = m1.stepj();
            const ptrdiff_t rowstep2 = m2.stepj();
            IT1 it1 = m1.get_diag(-m1.nlo()).begin();
            IT2 it2 = m2.get_diag(-m1.nlo()).begin();
            ptrdiff_t len = TMV_MIN(M-m1.nlo(),N);
            for(ptrdiff_t k=m1.nlo();k;--k) {
                SwapV_Helper<-3,xx,M1d,M2d>::call2(len,it1,it2);
                it1.shiftP(-colstep1);
                it2.shiftP(-colstep2);
                if (len < N) ++len;
            }
            TMVAssert(len == TMV_MIN(M,N));
            const ptrdiff_t ds = IntTraits2<cs,rs>::min;
            SwapV_Helper<-3,ds,M1d,M2d>::call2(len,it1,it2);
            for(ptrdiff_t k=1;k<=m1.nhi();++k) {
                it1.shiftP(rowstep1);
                it2.shiftP(rowstep2);
                if (k+len > N) --len;
                SwapV_Helper<-3,xx,M1d,M2d>::call2(len,it1,it2);
            }
        }
    };

    // algo 90: Call inst
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct SwapB_Helper<90,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        { InstSwap(m1.xView(),m2.xView()); }
    };

    // algo 91: Call inst alias
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct SwapB_Helper<91,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        { InstAliasSwap(m1.xView(),m2.xView()); }
    };

    // algo 97: Conjugate
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct SwapB_Helper<97,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        {
            typedef typename M1::conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            SwapB_Helper<-2,cs,rs,M1c,M2c>::call(m1c,m2c);
        }
    };

    // algo 197: Conjugate
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct SwapB_Helper<197,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        {
            typedef typename M1::conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            SwapB_Helper<99,cs,rs,M1c,M2c>::call(m1c,m2c);
        }
    };

    // algo 98: Inline check for aliases
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct SwapB_Helper<98,cs,rs,M1,M2>
    {
        static void call(M1& m1, M2& m2)
        {
            if (!SameStorage(m1,m2)) {
                // No aliasing
                SwapB_Helper<-2,cs,rs,M1,M2>::call(m1,m2);
            } else if (ExactSameStorage(m1,m2)) {
                // They are equal modulo a conjugation
                Maybe<M1::_conj != int(M2::_conj)>::conjself(m1);
            } else if (m2.colsize() == m2.rowsize() &&
                       m2.nlo() == m2.nhi() &&
                       m1.nlo() == m2.nhi() &&
                       m1.nhi() == m2.nhi() &&
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
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct SwapB_Helper<99,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst = 
                (cs == Unknown || cs > 16) &&
                (rs == Unknown || rs > 16) &&
                Traits2<T1,T2>::sametype &&
                Traits<T1>::isinst;
            const int algo =
                M2::_conj ? 197 :
                inst ? 91 :
                98;
            SwapB_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo -3: Determine which algorithm to use
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct SwapB_Helper<-3,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        {
            const int algo = 
                TMV_OPT == 0 ? 13 :
                ( M1::_colmajor && M2::_colmajor ) ? 11 :
                ( M1::_rowmajor && M2::_rowmajor ) ? 12 :
                13;
            SwapB_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo -2: Check for inst
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct SwapB_Helper<-2,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst = 
                (cs == Unknown || cs > 16) &&
                (rs == Unknown || rs > 16) &&
                Traits2<T1,T2>::sametype &&
                Traits<T1>::isinst;
            const int algo =
                M2::_conj ? 97 :
                inst ? 90 :
                -3;
            SwapB_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo -1: Check for aliases?
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct SwapB_Helper<-1,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        {
            const bool checkalias = 
                M1::_checkalias || M2::_checkalias;
            const int algo =
                checkalias ? 99 : 
                -2;
            SwapB_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    template <class M1, class M2>
    inline void DoSwap(
        BaseMatrix_Band_Mutable<M1>& m1, BaseMatrix_Band_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        const ptrdiff_t cs = Sizes<M1::_colsize,M2::_colsize>::size;
        const ptrdiff_t rs = Sizes<M1::_rowsize,M2::_rowsize>::size;
        typedef typename M1::cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_REF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        SwapB_Helper<-1,cs,rs,M1v,M2v>::call(m1v,m2v);
    }

    template <class M1, class M2>
    inline void InlineSwap(
        BaseMatrix_Band_Mutable<M1>& m1, BaseMatrix_Band_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        const ptrdiff_t cs = Sizes<M1::_colsize,M2::_colsize>::size;
        const ptrdiff_t rs = Sizes<M1::_rowsize,M2::_rowsize>::size;
        typedef typename M1::cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_REF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        SwapB_Helper<-3,cs,rs,M1v,M2v>::call(m1v,m2v);
    }

    template <class M1, class M2>
    inline void InlineAliasSwap(
        BaseMatrix_Band_Mutable<M1>& m1, BaseMatrix_Band_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        const ptrdiff_t cs = Sizes<M1::_colsize,M2::_colsize>::size;
        const ptrdiff_t rs = Sizes<M1::_rowsize,M2::_rowsize>::size;
        typedef typename M1::cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_REF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        SwapB_Helper<98,cs,rs,M1v,M2v>::call(m1v,m2v);
    }

    template <class M1, class M2>
    TMV_INLINE void Swap(
        BaseMatrix_Band_Mutable<M1>& m1, BaseMatrix_Band_Mutable<M2>& m2)
    { DoSwap(m1,m2); }

} // namespace tmv

#endif
