
#ifndef TMV_MultXB_H
#define TMV_MultXB_H

#include "TMV_BaseMatrix_Band.h"
#include "TMV_Scaling.h"
#include "TMV_MultXV.h"
#include "TMV_CopyB.h"

namespace tmv {

    // Defined in TMV_MultXB.cpp
    template <class T1, int C1, class T2>
    void InstMultXM(
        const T2 x, const ConstBandMatrixView<T1,C1>& m1,
        BandMatrixView<T2> m2);
    template <class T1, int C1, class T2>
    void InstAddMultXM(
        const T2 x, const ConstBandMatrixView<T1,C1>& m1,
        BandMatrixView<T2> m2);

    template <class T1, int C1, class T2>
    void InstAliasMultXM(
        const T2 x, const ConstBandMatrixView<T1,C1>& m1,
        BandMatrixView<T2> m2);
    template <class T1, int C1, class T2>
    void InstAliasAddMultXM(
        const T2 x, const ConstBandMatrixView<T1,C1>& m1,
        BandMatrixView<T2> m2);

    //
    // BandMatrix += x * BandMatrix
    //

    template <bool add, int ix, class T, class M1, class M2>
    static inline void NoAliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1, 
        BaseMatrix_Band_Mutable<M2>& m2);

    template <int algo, int cs, int rs, bool add, int ix, class T, class M1, class M2>
    struct MultXB_Helper;

    // algo 1: trivial: ix == 1, !add, so call Copy (with alias check)
    template <int cs, int rs, class T, class M1, class M2>
    struct MultXB_Helper<1,cs,rs,false,1,T,M1,M2>
    {
        static TMV_INLINE void call(const Scaling<1,T>& , const M1& m1, M2& m2)
        { CopyB_Helper<-3,cs,rs,M1,M2>::call(m1,m2); }
    };

    // algo 101: Same as algo 1, but use algo -3 for Copy
    template <int cs, int rs, class T, class M1, class M2>
    struct MultXB_Helper<101,cs,rs,false,1,T,M1,M2>
    {
        static TMV_INLINE void call(const Scaling<1,T>& , const M1& m1, M2& m2)
        { CopyB_Helper<-1,cs,rs,M1,M2>::call(m1,m2); }
    };

    // algo 201: Same as algo 1, but use algo -2 for Copy
    template <int cs, int rs, class T, class M1, class M2>
    struct MultXB_Helper<201,cs,rs,false,1,T,M1,M2>
    {
        static void call(const Scaling<1,T>& , const M1& m1, M2& m2)
        {
            // Need a quick alias check here, since ExactSameStorage
            // is allowed for MultXB, but not Copy
            if (!SameStorage(m1,m2)) {
                CopyB_Helper<-2,cs,rs,M1,M2>::call(m1,m2); 
            } else {
                TMVAssert(ExactSameStorage(m1,m2));
                Maybe<(M2::_conj != int(M1::_conj))>::conjself(m2);
            }
        }
    };

    // algo 11: Loop over columns
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
    struct MultXB_Helper<11,cs,rs,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            const int M = cs == TMV_UNKNOWN ? int(m2.colsize()) : cs;
            const int N = rs == TMV_UNKNOWN ? int(m2.rowsize()) : rs;
            const int xx = TMV_UNKNOWN;
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::col_sub_type M2c;
            typedef typename M1c::const_nonconj_type::const_iterator IT1;
            typedef typename M2c::iterator IT2;
            const int rowstep1 = m1.stepj();
            const int rowstep2 = m2.stepj();
            const int diagstep1 = m1.diagstep();
            const int diagstep2 = m2.diagstep();

            const int lh = IntTraits<IntTraits2<M1::_nlo,M1::_nhi>::sum>::Sp1;
            const int j1 = m1.nhi();
            const int j2 = TMV_MIN(N,M-m1.nlo());
            const int j3 = TMV_MIN(N,M+m1.nhi());
            int len = m1.nlo()+1;
            IT1 it1 = m1.get_col(0,0,len).begin().nonConj();
            IT2 it2 = m2.get_col(0,0,len).begin();
            for(int j=0;j<j1;++j) {
                MultXV_Helper<-4,xx,add,ix,T,M1c,M2c>::call2(len,x,it1,it2);
                it1.shiftP(rowstep1);
                it2.shiftP(rowstep2);
                if (len < M) ++len;
            }
            if (j1 < j2) TMVAssert(len == m1.nlo()+m1.nhi()+1);
            for(int j=j1;j<j2;++j) {
                MultXV_Helper<-4,lh,add,ix,T,M1c,M2c>::call2(len,x,it1,it2);
                it1.shiftP(diagstep1);
                it2.shiftP(diagstep2);
            }
            for(int j=j2;j<j3;++j) {
                MultXV_Helper<-4,xx,add,ix,T,M1c,M2c>::call2(--len,x,it1,it2);
                it1.shiftP(diagstep1);
                it2.shiftP(diagstep2);
            }
        }
    };

    // algo 12: Loop over rows
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
    struct MultXB_Helper<12,cs,rs,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            const int M = cs == TMV_UNKNOWN ? int(m2.colsize()) : cs;
            const int N = rs == TMV_UNKNOWN ? int(m2.rowsize()) : rs;
            const int xx = TMV_UNKNOWN;
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M2::row_sub_type M2r;
            typedef typename M1r::const_nonconj_type::const_iterator IT1;
            typedef typename M2r::iterator IT2;
            const int colstep1 = m1.stepi();
            const int colstep2 = m2.stepi();
            const int diagstep1 = m1.diagstep();
            const int diagstep2 = m2.diagstep();

            const int lh = IntTraits<IntTraits2<M1::_nlo,M1::_nhi>::sum>::Sp1;
            const int i1 = m1.nlo();
            const int i2 = TMV_MIN(M,N-m1.nhi());
            const int i3 = TMV_MIN(M,N+m1.nlo());
            int len = m1.nhi()+1;
            IT1 it1 = m1.get_row(0,0,len).begin().nonConj();
            IT2 it2 = m2.get_row(0,0,len).begin();
            for(int i=0;i<i1;++i) {
                MultXV_Helper<-4,xx,add,ix,T,M1r,M2r>::call2(len,x,it1,it2);
                it1.shiftP(colstep1);
                it2.shiftP(colstep2);
                if (len < N) ++len;
            }
            if (i1 < i2) TMVAssert(len == m1.nlo()+m1.nhi()+1);
            for(int i=i1;i<i2;++i) {
                MultXV_Helper<-4,lh,add,ix,T,M1r,M2r>::call2(len,x,it1,it2);
                it1.shiftP(diagstep1);
                it2.shiftP(diagstep2);
            }
            for(int i=i2;i<i3;++i) {
                MultXV_Helper<-4,xx,add,ix,T,M1r,M2r>::call2(--len,x,it1,it2);
                it1.shiftP(diagstep1);
                it2.shiftP(diagstep2);
            }
        }
    };

    // algo 13: Loop over diagonals
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
    struct MultXB_Helper<13,cs,rs,add,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            const int M = cs == TMV_UNKNOWN ? int(m2.colsize()) : cs;
            const int N = rs == TMV_UNKNOWN ? int(m2.rowsize()) : rs;
            const int xx = TMV_UNKNOWN;
            typedef typename M1::const_diag_sub_type M1d;
            typedef typename M2::diag_sub_type M2d;
            typedef typename M1d::const_nonconj_type::const_iterator IT1;
            typedef typename M2d::iterator IT2;
            const int colstep1 = m1.stepi();
            const int colstep2 = m2.stepi();
            const int rowstep1 = m1.stepj();
            const int rowstep2 = m2.stepj();
            IT1 it1 = m1.get_diag(-m1.nlo()).begin().nonConj();
            IT2 it2 = m2.get_diag(-m1.nlo()).begin();
            int len = TMV_MIN(M-m1.nlo(),N);
            for(int k=m1.nlo();k;--k) {
                MultXV_Helper<-4,xx,add,ix,T,M1d,M2d>::call2(len,x,it1,it2);
                it1.shiftP(-colstep1);
                it2.shiftP(-colstep2);
                if (len < N) ++len;
            }
            TMVAssert(len == TMV_MIN(M,N));
            const int ds = IntTraits2<cs,rs>::min;
            MultXV_Helper<-4,ds,add,ix,T,M1d,M2d>::call2(len,x,it1,it2);
            for(int k=0;k<m1.nhi();++k) {
                it1.shiftP(rowstep1);
                it2.shiftP(rowstep2);
                if (k+len >= N) --len;
                MultXV_Helper<-4,xx,add,ix,T,M1d,M2d>::call2(len,x,it1,it2);
            }
        }
    };

    // algo 90: Call inst
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct MultXB_Helper<90,cs,rs,true,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M2::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAddMultXM(xx,m1.xView(),m2.xView()); 
        }
    };
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct MultXB_Helper<90,cs,rs,false,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M2::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstMultXM(xx,m1.xView(),m2.xView()); 
        }
    };

    // algo 91: Call inst alias
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct MultXB_Helper<91,cs,rs,true,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M2::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasAddMultXM(xx,m1.xView(),m2.xView()); 
        }
    };
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct MultXB_Helper<91,cs,rs,false,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M2::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstAliasMultXM(xx,m1.xView(),m2.xView()); 
        }
    };

    // algo 97: Conjugate
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
    struct MultXB_Helper<97,cs,rs,add,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            MultXB_Helper<-2,cs,rs,add,ix,T,M1c,M2c>::call(
                TMV_CONJ(x),m1c,m2c);
        }
    };

    // algo 197: Conjugate
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
    struct MultXB_Helper<197,cs,rs,add,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            MultXB_Helper<99,cs,rs,add,ix,T,M1c,M2c>::call(
                TMV_CONJ(x),m1c,m2c);
        }
    };

    // algo 98: Inline check for aliases
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct MultXB_Helper<98,cs,rs,true,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            // More than just these are non-clobbering.
            // TODO: Add in the steps that work.
            if ( !SameStorage(m1,m2) ||
                 ExactSameStorage(m1,m2) ) {
                // No aliasing (or no clobbering)
                MultXB_Helper<-2,cs,rs,true,ix,T,M1,M2>::call(x,m1,m2);
            } else {
                // Need a temporary
                NoAliasMultXM<true>(x,m1.copy(),m2);
            }
        }
    };
    template <int cs, int rs, int ix, class T, class M1, class M2>
    struct MultXB_Helper<98,cs,rs,false,ix,T,M1,M2>
    {
        static void call(const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            if ( !SameStorage(m1,m2) ||
                 ExactSameStorage(m1,m2) ) {
                // No aliasing (or no clobbering)
                MultXB_Helper<-2,cs,rs,false,ix,T,M1,M2>::call(x,m1,m2);
            } else {
                // Let Copy handle the aliasing
                AliasCopy(m1,m2);
                Scale(x,m2);
            }
        }
    };

    // algo 99: Check for aliases
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
    struct MultXB_Helper<99,cs,rs,add,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
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
            MultXB_Helper<algo,cs,rs,add,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    // algo -4: No branches or copies 
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
    struct MultXB_Helper<-4,cs,rs,add,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            typedef typename M2::value_type T2;
            const int algo = 
                ( ix == 1 && !add ) ? 1 :
                TMV_OPT == 0 ? 13 :
                ( M1::_colmajor && M2::_colmajor ) ? 11 : 
                ( M1::_rowmajor && M2::_rowmajor ) ? 12 : 
                13;
#ifdef PRINTALGO_XB
            std::cout<<"XB: algo = "<<algo<<std::endl;
#endif
            MultXB_Helper<algo,cs,rs,add,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
    struct MultXB_Helper<-3,cs,rs,add,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            MultXB_Helper<-4,cs,rs,add,ix,T,M1,M2>::call(x,m1,m2); 
            if (m2.nlo() > m1.nlo())
                Maybe<!add>::zero2(m2.diagRange(-m2.nlo(),-m1.nlo()));
            if (m2.nhi() > m1.nhi())
                Maybe<!add>::zero2(m2.diagRange(m1.nhi()+1,m2.nhi()+1));
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
    struct MultXB_Helper<-2,cs,rs,add,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
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
                ( ix == 1 && !add ) ? 201 :
                M2::_conj ? 97 :
                inst ? 90 :
                -3;
            MultXB_Helper<algo,cs,rs,add,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    // algo -1: Check for aliases?
    template <int cs, int rs, bool add, int ix, class T, class M1, class M2>
    struct MultXB_Helper<-1,cs,rs,add,ix,T,M1,M2>
    {
        static TMV_INLINE void call(
            const Scaling<ix,T>& x, const M1& m1, M2& m2)
        {
            const int algo = 
                ( ix == 1 && !add ) ? 101 :
                M2::_checkalias ? 99 : 
                -2;
            MultXB_Helper<algo,cs,rs,add,ix,T,M1,M2>::call(x,m1,m2);
        }
    };

    template <bool add, int ix, class T, class M1, class M2>
    static inline void MultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1, 
        BaseMatrix_Band_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.nlo() <= m2.nlo());
        TMVAssert(m1.nhi() <= m2.nhi());
        const int cs = Sizes<M1::_colsize,M2::_colsize>::size;
        const int rs = Sizes<M1::_rowsize,M2::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        MultXB_Helper<-1,cs,rs,add,ix,T,M1v,M2v>::call(x,m1v,m2v);
    }

    template <bool add, int ix, class T, class M1, class M2>
    static inline void NoAliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1, 
        BaseMatrix_Band_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.nlo() <= m2.nlo());
        TMVAssert(m1.nhi() <= m2.nhi());
        const int cs = Sizes<M1::_colsize,M2::_colsize>::size;
        const int rs = Sizes<M1::_rowsize,M2::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        MultXB_Helper<-2,cs,rs,add,ix,T,M1v,M2v>::call(x,m1v,m2v);
    }

    template <bool add, int ix, class T, class M1, class M2>
    static inline void InlineMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1, 
        BaseMatrix_Band_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.nlo() <= m2.nlo());
        TMVAssert(m1.nhi() <= m2.nhi());
        const int cs = Sizes<M1::_colsize,M2::_colsize>::size;
        const int rs = Sizes<M1::_rowsize,M2::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        MultXB_Helper<-3,cs,rs,add,ix,T,M1v,M2v>::call(x,m1v,m2v);
    }

    template <bool add, int ix, class T, class M1, class M2>
    static inline void InlineAliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1, 
        BaseMatrix_Band_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.nlo() <= m2.nlo());
        TMVAssert(m1.nhi() <= m2.nhi());
        const int cs = Sizes<M1::_colsize,M2::_colsize>::size;
        const int rs = Sizes<M1::_rowsize,M2::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        MultXB_Helper<98,cs,rs,add,ix,T,M1v,M2v>::call(x,m1v,m2v);
    }

    template <bool add, int ix, class T, class M1, class M2>
    static inline void AliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1, 
        BaseMatrix_Band_Mutable<M2>& m2)
    {
        TMVStaticAssert((Sizes<M1::_colsize,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_rowsize,M2::_rowsize>::same));
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        TMVAssert(m1.nlo() <= m2.nlo());
        TMVAssert(m1.nhi() <= m2.nhi());
        const int cs = Sizes<M1::_colsize,M2::_colsize>::size;
        const int rs = Sizes<M1::_rowsize,M2::_rowsize>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        MultXB_Helper<99,cs,rs,add,ix,T,M1v,M2v>::call(x,m1v,m2v);
    }


    //
    // M (+)= x * B
    //

    template <bool add, int ix, class T, class M1, class M2>
    static inline void MultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        BaseMatrix_Rec_Mutable<M2>& m2)
    {
        typename BMVO<M2>::b b2 = BandMatrixViewOf(m2,m1.nlo(),m1.nhi());
        MultXM<add>(x,m1,b2);
        if (m1.nlo() < m2.colsize()-1)
            Maybe<!add>::zero2(
                BandMatrixViewOf(
                    m2.cRowRange(m1.nlo()+1,m1.colsize()),
                    m1.colsize()-m1.nlo()-2,0));
        if (m1.nhi() < m2.rowsize()-1)
            Maybe<!add>::zero2(
                BandMatrixViewOf(
                    m2.cColRange(m1.nhi()+1,m1.rowsize()),
                    0,m1.rowsize()-m1.nhi()-2));
    }

    template <bool add, int ix, class T, class M1, class M2>
    static inline void NoAliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        BaseMatrix_Rec_Mutable<M2>& m2)
    {
        m2.setZero();
        typename BMVO<M2>::b b2 = BandMatrixViewOf(m2,m1.nlo(),m1.nhi());
        NoAliasMultXM<add>(x,m1,b2);
    }

    template <bool add, int ix, class T, class M1, class M2>
    static inline void AliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        BaseMatrix_Rec_Mutable<M2>& m2)
    {
        typename BMVO<M2>::b b2 = BandMatrixViewOf(m2,m1.nlo(),m1.nhi());
        AliasMultXM<add>(x,m1,b2);
        if (m1.nlo() < m1.colsize()-1)
            Maybe<!add>::zero2(
                BandMatrixViewOf(
                    m2.cRowRange(m1.nlo()+1,m1.colsize()),
                    m1.colsize()-m1.nlo()-2,0));
        if (m1.nhi() < m1.rowsize()-1)
            Maybe<!add>::zero2(
                BandMatrixViewOf(
                    m2.cColRange(m1.nhi()+1,m1.rowsize()),
                    0,m1.rowsize()-m1.nhi()-2));
    }


    //
    // U (+)= x * B
    //

    template <bool add, int ix, class T, class M1, class M2>
    static inline void MultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        BaseMatrix_Tri_Mutable<M2>& m2)
    {
        const int lo = Maybe<M2::_upper>::select(m1.nlo(),m1.nhi());
        const int hi = Maybe<M2::_upper>::select(m1.nhi(),m1.nlo());
        TMVAssert(lo == 0);
        typename BMVOTri<M2>::b b2 = BandMatrixViewOf(m2,hi);
        MultXM<add>(x,m1,b2);
        Maybe<!add>::zero2(m2.offDiag(hi+1));
    }

    template <bool add, int ix, class T, class M1, class M2>
    static inline void NoAliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        BaseMatrix_Tri_Mutable<M2>& m2)
    {
        const int lo = Maybe<M2::_upper>::select(m1.nlo(),m1.nhi());
        const int hi = Maybe<M2::_upper>::select(m1.nhi(),m1.nlo());
        TMVAssert(lo == 0);
        typename BMVOTri<M2>::b b2 = BandMatrixViewOf(m2,hi);
        m2.setZero();
        NoAliasMultXM<add>(x,m1,b2);
    }

    template <bool add, int ix, class T, class M1, class M2>
    static inline void AliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Band<M1>& m1,
        BaseMatrix_Tri_Mutable<M2>& m2)
    {
        const int lo = Maybe<M2::_upper>::select(m1.nlo(),m1.nhi());
        const int hi = Maybe<M2::_upper>::select(m1.nhi(),m1.nlo());
        TMVAssert(lo == 0);
        typename BMVOTri<M2>::b b2 = BandMatrixViewOf(m2,hi);
        AliasMultXM<add>(x,m1,b2);
        Maybe<!add>::zero2(m2.offDiag(hi+1));
    }


    //
    // B (+)= x * D
    //

    template <bool add, int ix, class T, class M1, class M2>
    static inline void MultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        BaseMatrix_Band_Mutable<M2>& m2)
    {
        typename M1::const_diag_type d1 = m1.diag();
        typename M2::diag_type d2 = m2.diag();
        MultXV<add>(x,d1,d2);
        if (m2.nlo() > 0) Maybe<!add>::zero2(m2.cDiagRange(-m2.nlo(),0));
        if (m2.nhi() > 0) Maybe<!add>::zero2(m2.cDiagRange(1,m2.nhi()+1));
    }

    template <bool add, int ix, class T, class M1, class M2>
    static inline void NoAliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        BaseMatrix_Band_Mutable<M2>& m2)
    {
        typename M1::const_diag_type d1 = m1.diag();
        typename M2::diag_type d2 = m2.diag();
        m2.setZero();
        NoAliasMultXV<add>(x,d1,d2);
    }

    template <bool add, int ix, class T, class M1, class M2>
    static inline void AliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Diag<M1>& m1,
        BaseMatrix_Band_Mutable<M2>& m2)
    {
        typename M1::const_diag_type d1 = m1.diag();
        typename M2::diag_type d2 = m2.diag();
        AliasMultXV<add>(x,d1,d2);
        if (m2.nlo() > 0) Maybe<!add>::zero2(m2.cDiagRange(-m2.nlo(),0));
        if (m2.nhi() > 0) Maybe<!add>::zero2(m2.cDiagRange(1,m2.nhi()+1));
    }

    //
    // B (+)= x * U
    //

    template <bool add, int ix, class T, class M1, class M2>
    static inline void MultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        BaseMatrix_Band_Mutable<M2>& m2)
    {
        const int k1 = Maybe<M1::_upper>::select(0,-m1.size()+1);
        const int k2 = Maybe<M1::_upper>::select(m1.size(),1);
        const int k3 = Maybe<M1::_upper>::select(m2.nlo(),1);
        const int k4 = Maybe<M1::_upper>::select(0,m2.nhi()+1);
        TMVAssert(Maybe<M1::_upper>::select(m2.nhi(),m2.nlo()) == m1.size()-1);
        typename M2::diagrange_type u2 = m2.cDiagRange(k1,k2);
        MultXV<add>(x,m1,u2);
        if (k4 > k3) Maybe<!add>::zero2(m2.cDiagRange(k3,k4));
    }

    template <bool add, int ix, class T, class M1, class M2>
    static inline void NoAliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        BaseMatrix_Band_Mutable<M2>& m2)
    {
        const int k1 = Maybe<M1::_upper>::select(0,-m1.size()+1);
        const int k2 = Maybe<M1::_upper>::select(m1.size(),1);
        const int k3 = Maybe<M1::_upper>::select(m2.nlo(),1);
        const int k4 = Maybe<M1::_upper>::select(0,m2.nhi()+1);
        TMVAssert(Maybe<M1::_upper>::select(m2.nhi(),m2.nlo()) == m1.size()-1);
        typename M2::diagrange_type u2 = m2.cDiagRange(k1,k2);
        NoAliasMultXV<add>(x,m1,u2);
        if (k4 > k3) Maybe<!add>::zero2(m2.cDiagRange(k3,k4));
    }

    template <bool add, int ix, class T, class M1, class M2>
    static inline void AliasMultXM(
        const Scaling<ix,T>& x, const BaseMatrix_Tri<M1>& m1,
        BaseMatrix_Band_Mutable<M2>& m2)
    {
        const int k1 = Maybe<M1::_upper>::select(0,-m1.size()+1);
        const int k2 = Maybe<M1::_upper>::select(m1.size(),1);
        const int k3 = Maybe<M1::_upper>::select(m2.nlo(),1);
        const int k4 = Maybe<M1::_upper>::select(0,m2.nhi()+1);
        TMVAssert(Maybe<M1::_upper>::select(m2.nhi(),m2.nlo()) == m1.size()-1);
        typename M2::diagrange_type u2 = m2.cDiagRange(k1,k2);
        AliasMultXV<add>(x,m1,u2);
        if (k4 > k3) Maybe<!add>::zero2(m2.cDiagRange(k3,k4));
    }


} // namespace tmv

#endif 
