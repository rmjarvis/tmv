
#ifndef TMV_CopyB_H
#define TMV_CopyB_H

#include "TMV_BaseMatrix_Band.h"
#include "TMV_CopyV.h"
#include "TMV_SmallBandMatrix.h"

namespace tmv {

    // Defined in TMV_Matrix.cpp
    template <class T1, int C1, class T2>
    void InstCopy(
        const ConstBandMatrixView<T1,C1>& m1, BandMatrixView<T2> m2); 
    template <class T1, int C1, class T2>
    void InstAliasCopy(
        const ConstBandMatrixView<T1,C1>& m1, BandMatrixView<T2> m2); 

    //
    // Copy Matrices
    //

    template <int algo, ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct CopyB_Helper;

    // algo 0: size = 0, nothing to do
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct CopyB_Helper<0,cs,rs,M1,M2>
    { static TMV_INLINE void call(const M1& , M2& ) {} };

    // algo 11: Loop over columns
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct CopyB_Helper<11,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
            const ptrdiff_t M = cs == Unknown ? m2.colsize() : cs;
            const ptrdiff_t N = rs == Unknown ? m2.rowsize() : rs;
            const ptrdiff_t xx = Unknown;
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::col_sub_type M2c;
            typedef typename M1c::const_iterator IT1;
            typedef typename M2c::iterator IT2;
            const ptrdiff_t rowstep1 = m1.stepj();
            const ptrdiff_t rowstep2 = m2.stepj();
            const ptrdiff_t diagstep1 = m1.diagstep();
            const ptrdiff_t diagstep2 = m2.diagstep();

            // Split calculation up into 3 sections.
            // The first and third require care about the number of 
            // elements in each column, but the middle section has the 
            // same number of elements in each: m1.nlo()+m1.nhi()+1
            // so we can use the compile-time knowledge if possible
            const ptrdiff_t lh = IntTraits<IntTraits2<M1::_nlo,M1::_nhi>::sum>::Sp1;
            const ptrdiff_t j1 = m1.nhi();
            const ptrdiff_t j2 = TMV_MIN(N,M-m1.nlo());
            const ptrdiff_t j3 = TMV_MIN(N,M+m1.nhi());
            ptrdiff_t len = m1.nlo()+1;
            IT1 it1 = m1.get_col(0,0,len).begin();
            IT2 it2 = m2.get_col(0,0,len).begin();
            ptrdiff_t j=0;
            for(;j<j1;++j) {
                CopyV_Helper<-3,xx,M1c,M2c>::call2(len,it1,it2);
                it1.shiftP(rowstep1);
                it2.shiftP(rowstep2);
                if (len < M) ++len;
            }
            if (j1 < j2) TMVAssert(len == m1.nlo()+m1.nhi()+1);
            for(;j<j2;++j) {
                CopyV_Helper<-3,lh,M1c,M2c>::call2(len,it1,it2);
                it1.shiftP(diagstep1);
                it2.shiftP(diagstep2);
            }
            if (j1 >= j2) ++len;
            for(;j<j3;++j) {                
                CopyV_Helper<-3,xx,M1c,M2c>::call2(--len,it1,it2);
                it1.shiftP(diagstep1);
                it2.shiftP(diagstep2);
            }
        }
    };

    // algo 12: Loop over rows
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct CopyB_Helper<12,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
            const ptrdiff_t M = cs == Unknown ? m2.colsize() : cs;
            const ptrdiff_t N = rs == Unknown ? m2.rowsize() : rs;
            const ptrdiff_t xx = Unknown;
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M2::row_sub_type M2r;
            typedef typename M1r::const_iterator IT1;
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
                CopyV_Helper<-3,xx,M1r,M2r>::call2(len,it1,it2);
                it1.shiftP(colstep1);
                it2.shiftP(colstep2);
                if (len < N) ++len;
            }
            if (i1 < i2) TMVAssert(len == m1.nlo()+m1.nhi()+1);
            for(;i<i2;++i) {
                CopyV_Helper<-3,lh,M1r,M2r>::call2(len,it1,it2);
                it1.shiftP(diagstep1);
                it2.shiftP(diagstep2);
            }
            if (i1 >= i2) ++len;
            for(;i<i3;++i) {                
                CopyV_Helper<-3,xx,M1r,M2r>::call2(--len,it1,it2);
                it1.shiftP(diagstep1);
                it2.shiftP(diagstep2);
            }
        }
    };

    // algo 13: Loop over diagonals
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct CopyB_Helper<13,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
            const ptrdiff_t M = cs == Unknown ? m2.colsize() : cs;
            const ptrdiff_t N = rs == Unknown ? m2.rowsize() : rs;
            const ptrdiff_t xx = Unknown;
            typedef typename M1::const_diag_sub_type M1d;
            typedef typename M2::diag_sub_type M2d;
            typedef typename M1d::const_iterator IT1;
            typedef typename M2d::iterator IT2;
            const ptrdiff_t colstep1 = m1.stepi();
            const ptrdiff_t colstep2 = m2.stepi();
            const ptrdiff_t rowstep1 = m1.stepj();
            const ptrdiff_t rowstep2 = m2.stepj();
            IT1 it1 = m1.get_diag(-m1.nlo()).begin();
            IT2 it2 = m2.get_diag(-m1.nlo()).begin();
            ptrdiff_t len = TMV_MIN(M-m1.nlo(),N);
            for(ptrdiff_t k=m1.nlo();k;--k) {
                CopyV_Helper<-3,xx,M1d,M2d>::call2(len,it1,it2);
                it1.shiftP(-colstep1);
                it2.shiftP(-colstep2);
                if (len < N) ++len;
            }
            TMVAssert(len == TMV_MIN(M,N));
            const ptrdiff_t ds = IntTraits2<cs,rs>::min;
            CopyV_Helper<-3,ds,M1d,M2d>::call2(len,it1,it2);
            for(ptrdiff_t k=1;k<=m1.nhi();++k) {
                it1.shiftP(rowstep1);
                it2.shiftP(rowstep2);
                if (k+len > N) --len;
                CopyV_Helper<-3,xx,M1d,M2d>::call2(len,it1,it2);
            }
        }
    };

    // algo 90: Call inst
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct CopyB_Helper<90,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        { InstCopy(m1.xView(),m2.xView()); }
    };

    // algo 91: Call inst alias
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct CopyB_Helper<91,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        { InstAliasCopy(m1.xView(),m2.xView()); }
    };

    // algo 97: Conjugate
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct CopyB_Helper<97,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            CopyB_Helper<-2,cs,rs,M1c,M2c>::call(m1c,m2c);
        }
    };

    // algo 197: Conjugate
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct CopyB_Helper<197,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            CopyB_Helper<99,cs,rs,M1c,M2c>::call(m1c,m2c);
        }
    };

    // algo 98: Inline check for aliases
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct CopyB_Helper<98,cs,rs,M1,M2>
    {
        static void call(const M1& m1, M2& m2)
        {
            if (!SameStorage(m1,m2)) {
                // No aliasing
                CopyB_Helper<-2,cs,rs,M1,M2>::call(m1,m2);
            } else if (ExactSameStorage(m1,m2)) {
                // They are already equal modulo a conjugation.
                Maybe<M1::_conj != int(M2::_conj)>::conjself(m2);
            } else if (m2.colsize() == m2.rowsize() &&
                       m2.nlo() == m2.nhi() &&
                       m1.nlo() == m2.nhi() &&
                       m1.nhi() == m2.nhi() &&
                       OppositeStorage(m1,m2)) {
                // Then transpose
                m2.transposeSelf();
                // And maybe conjugate
                Maybe<M1::_conj != int(M2::_conj)>::conjself(m2);
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
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct CopyB_Helper<99,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst = 
                (cs == Unknown || cs > 16) &&
                (rs == Unknown || rs > 16) &&
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
            CopyB_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo -4: No branches or copies
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct CopyB_Helper<-4,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            typedef typename M2::value_type T2;
            const int algo = 
                (cs == 0 || rs == 0) ? 0 :
                TMV_OPT == 0 ? 13 :
                (M1::_colmajor && M2::_colmajor) ? 11 :
                (M1::_rowmajor && M2::_rowmajor) ? 12 :
                13;
            CopyB_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo -3: Determine which algorithm to use
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct CopyB_Helper<-3,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            CopyB_Helper<-4,cs,rs,M1,M2>::call(m1,m2); 
            if (m2.nlo() > m1.nlo()) 
                m2.cDiagRange(-m2.nlo(),-m1.nlo()).setZero();
            if (m2.nhi() > m1.nhi()) 
                m2.cDiagRange(m1.nhi()+1,m2.nhi()+1).setZero();
        }
    };

    // algo -2: Check for inst
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct CopyB_Helper<-2,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst = 
                (cs == Unknown || cs > 16) &&
                (rs == Unknown || rs > 16) &&
#ifdef TMV_INST_MIX
                Traits2<T1,T2>::samebase &&
#else
                Traits2<T1,T2>::sametype &&
#endif
                Traits<T1>::isinst;
            const int algo = 
                M2::_conj ? 97 :
                inst ? 90 :
                -3;
            CopyB_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    // algo -1: Check for aliases?
    template <ptrdiff_t cs, ptrdiff_t rs, class M1, class M2>
    struct CopyB_Helper<-1,cs,rs,M1,M2>
    {
        static TMV_INLINE void call(const M1& m1, M2& m2)
        {
            const int algo = 
                M2::_checkalias ? 99 : 
                -2;
            CopyB_Helper<algo,cs,rs,M1,M2>::call(m1,m2);
        }
    };

    template <class M1, class M2>
    inline void Copy(
        const BaseMatrix_Band<M1>& m1, BaseMatrix_Band_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.nlo() <= m2.nlo());
        TMVAssert(m1.nhi() <= m2.nhi());
        const ptrdiff_t cs = Sizes<M1::_colsize,M2::_colsize>::size;
        const ptrdiff_t rs = Sizes<M1::_rowsize,M2::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        CopyB_Helper<-1,cs,rs,M1v,M2v>::call(m1v,m2v);
    }

    template <class M1, class M2>
    inline void InlineCopy(
        const BaseMatrix_Band<M1>& m1, BaseMatrix_Band_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.nlo() <= m2.nlo());
        TMVAssert(m1.nhi() <= m2.nhi());
        const ptrdiff_t cs = Sizes<M1::_colsize,M2::_colsize>::size;
        const ptrdiff_t rs = Sizes<M1::_rowsize,M2::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        CopyB_Helper<-3,cs,rs,M1v,M2v>::call(m1v,m2v);
    }

    template <class M1, class M2>
    inline void InlineAliasCopy(
        const BaseMatrix_Band<M1>& m1, BaseMatrix_Band_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.nlo() <= m2.nlo());
        TMVAssert(m1.nhi() <= m2.nhi());
        const ptrdiff_t cs = Sizes<M1::_colsize,M2::_colsize>::size;
        const ptrdiff_t rs = Sizes<M1::_rowsize,M2::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        CopyB_Helper<98,cs,rs,M1v,M2v>::call(m1v,m2v);
    }

    //
    // M = B
    //

    template <class M1, class M2>
    inline void Copy(
        const BaseMatrix_Band<M1>& m1, BaseMatrix_Rec_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        typename BMVO<M2>::b b2 = BandMatrixViewOf(m2,m1.nlo(),m1.nhi());
        Copy(m1,b2);
        if (m1.nlo() < m1.colsize()-1)
            BandMatrixViewOf(
                m2.cRowRange(m1.nlo()+1,m1.colsize()),
                m1.colsize()-m1.nlo()-2,0).setZero();
        if (m1.nhi() < m1.rowsize()-1)
            BandMatrixViewOf(
                m2.cColRange(m1.nhi()+1,m1.rowsize()),
                0,m1.rowsize()-m1.nhi()-2).setZero();
    }


    //
    // U = B
    //

    template <class M1, class M2>
    inline void Copy(
        const BaseMatrix_Band<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_size>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_size>::same));
        TMVStaticAssert(!M2::_unit);
        TMVAssert(m1.colsize() == m2.size());
        TMVAssert(m1.rowsize() == m2.size());
        const ptrdiff_t lo = Maybe<M2::_upper>::select(m1.nlo(),m1.nhi());
        const ptrdiff_t hi = Maybe<M2::_upper>::select(m1.nhi(),m1.nlo());
        TMVAssert(lo == 0);
        typename BMVOTri<M2>::b b2 = BandMatrixViewOf(m2,hi);
        Copy(m1,b2);
        m2.offDiag(hi+1).setZero();
    }


    //
    // B = D
    //

    template <class M1, class M2>
    inline void Copy(
        const BaseMatrix_Diag<M1>& m1, BaseMatrix_Band_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_size,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_size,M2::_rowsize>::same));
        TMVAssert(m1.size() == m2.colsize());
        TMVAssert(m1.size() == m2.rowsize());
        typename M1::const_diag_type d1 = m1.diag();
        typename M2::diag_type d2 = m2.diag();
        Copy(d1,d2);
        if (m2.nlo() > 0) m2.cDiagRange(-m2.nlo(),0).setZero();
        if (m2.nhi() > 0) m2.cDiagRange(1,m2.nhi()+1).setZero();
    }

    //
    // B = U
    //

    template <class M1, class M2>
    inline void Copy(
        const BaseMatrix_Tri<M1>& m1, BaseMatrix_Band_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_size,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_size,M2::_rowsize>::same));
        TMVAssert(m1.size() == m2.colsize());
        TMVAssert(m1.size() == m2.rowsize());
        const ptrdiff_t k1 = Maybe<M1::_upper>::select(0,-m1.size()+1);
        const ptrdiff_t k2 = Maybe<M1::_upper>::select(m1.size(),1);
        const ptrdiff_t k3 = Maybe<M1::_upper>::select(m2.nlo(),1);
        const ptrdiff_t k4 = Maybe<M1::_upper>::select(0,m2.nhi()+1);
        TMVAssert(Maybe<M1::_upper>::select(m2.nhi(),m2.nlo()) == m1.size()-1);
        typename M2::diagrange_type u2 = m2.cDiagRange(k1,k2);
        Copy(m1,u2);
        if (k4 > k3) m2.cDiagRange(k3,k4).setZero();
    }

} // namespace tmv

#endif
