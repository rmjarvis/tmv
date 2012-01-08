
//
// This file defines the basic composite type for a product of a 
// matrix and a scalar.  It also implements the calculation for 
// dense rectangular matrices.

#ifndef TMV_ScaleB_H
#define TMV_ScaleB_H

#include "TMV_BaseMatrix_Band.h"
#include "TMV_ScaleV.h"

#ifdef PRINTALGO_XB
#include <iostream>
#include "TMV_MatrixIO.h"
#include "TMV_BandMatrixIO.h"
#endif

namespace tmv {

    // Defined in TMV_ScaleB.cpp
    template <class T>
    void InstScale(const T x, BandMatrixView<T> m);

    //
    // BandMatrix *= x
    //

    template <int algo, int cs, int rs, int ix, class T, class M1>
    struct ScaleB_Helper;

    // algo 0: trivial: cs == 0, rs == 0 or ix == 1, so nothing to do
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleB_Helper<0,cs,rs,ix,T,M1>
    { static TMV_INLINE void call(const Scaling<ix,T>& , M1& ) {} };

    // algo 1: Linearize to vector version
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleB_Helper<1,cs,rs,ix,T,M1>
    {
        static TMV_INLINE void call(const Scaling<ix,T>& x, M1& m)
        {
#ifdef PRINTALGO_XB
            std::cout<<"ScaleB algo 1\n";
#endif
            TMVStaticAssert(M1::_canlin);
            typedef typename M1::linearview_type Ml;
            Ml ml = m.linearView();
            ScaleV_Helper<-3,Ml::_size,ix,T,Ml>::call(x,ml);
        }
    };

    // algo 11: Loop over columns
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleB_Helper<11,cs,rs,ix,T,M1>
    {
        static void call(const Scaling<ix,T>& x, M1& m)
        {
            const int M = cs == Unknown ? m.colsize() : cs;
            const int N = rs == Unknown ? m.rowsize() : rs;
#ifdef PRINTALGO_XB
            std::cout<<"ScaleB algo 11: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            const int xx = Unknown;
            typedef typename M1::col_sub_type M1c;
            typedef typename M1c::iterator IT;
            const int rowstep = m.stepj();
            const int diagstep = m.diagstep();

            const int lh = IntTraits<IntTraits2<M1::_nlo,M1::_nhi>::sum>::Sp1;
            const int j1 = m.nhi();
            const int j2 = TMV_MIN(N,M-m.nlo());
            const int j3 = TMV_MIN(N,M+m.nhi());
            //std::cout<<"j1,j2,j3 = "<<j1<<','<<j2<<','<<j3<<std::endl;
            int len = m.nlo()+1;
            IT it = m.get_col(0,0,len).begin();
            int j=0;
            for(;j<j1;++j) {
                //std::cout<<"A j = "<<j<<", len  = "<<len<<std::endl;
                ScaleV_Helper<-3,xx,ix,T,M1c>::call2(len,x,it);
                it.shiftP(rowstep);
                if (len < M) ++len;
            }
            for(;j<j2;++j) {
                //std::cout<<"B j = "<<j<<", len  = "<<len<<std::endl;
                ScaleV_Helper<-3,lh,ix,T,M1c>::call2(len,x,it);
                it.shiftP(diagstep);
            }
            if (j1 >= j2) ++len;
            for(;j<j3;++j) {
                //std::cout<<"C j = "<<j<<", len  = "<<len<<std::endl;
                ScaleV_Helper<-3,xx,ix,T,M1c>::call2(--len,x,it);
                it.shiftP(diagstep);
            }
        }
    };

    // algo 12: Loop over rows
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleB_Helper<12,cs,rs,ix,T,M1>
    {
        static void call(const Scaling<ix,T>& x, M1& m)
        {
            const int M = cs == Unknown ? m.colsize() : cs;
            const int N = rs == Unknown ? m.rowsize() : rs;
#ifdef PRINTALGO_XB
            std::cout<<"ScaleB algo 11: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            const int xx = Unknown;
            typedef typename M1::row_sub_type M1r;
            typedef typename M1r::iterator IT;
            const int colstep = m.stepi();
            const int diagstep = m.diagstep();

            const int lh = IntTraits<IntTraits2<M1::_nlo,M1::_nhi>::sum>::Sp1;
            const int i1 = m.nlo();
            const int i2 = TMV_MIN(M,N-m.nhi());
            const int i3 = TMV_MIN(M,N+m.nlo());
            int len = m.nhi()+1;
            IT it = m.get_row(0,0,len).begin();
            int i=0;
            for(;i<i1;++i) {
                ScaleV_Helper<-3,xx,ix,T,M1r>::call2(len,x,it);
                it.shiftP(colstep);
                if (len < N) ++len;
            }
            for(;i<i2;++i) {
                ScaleV_Helper<-3,lh,ix,T,M1r>::call2(len,x,it);
                it.shiftP(diagstep);
            }
            if (i1 >= i2) ++len;
            for(;i<i3;++i) {
                ScaleV_Helper<-3,xx,ix,T,M1r>::call2(--len,x,it);
                it.shiftP(diagstep);
            }
        }
    };

    // algo 13: Loop over diagonals
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleB_Helper<13,cs,rs,ix,T,M1>
    {
        static void call(const Scaling<ix,T>& x, M1& m)
        {
            const int M = cs == Unknown ? m.colsize() : cs;
            const int N = rs == Unknown ? m.rowsize() : rs;
#ifdef PRINTALGO_XB
            std::cout<<"ScaleB algo 11: M,N,cs,rs,x = "<<M<<','<<N<<
                ','<<cs<<','<<rs<<','<<T(x)<<std::endl;
#endif
            const int xx = Unknown;
            typedef typename M1::diag_sub_type M1d;
            typedef typename M1d::iterator IT;
            const int colstep = m.stepi();
            const int rowstep = m.stepj();
            IT it = m.get_diag(-m.nlo()).begin();
            int len = TMV_MIN(M-m.nlo(),N);
            for(int k=m.nlo();k;--k) {
                ScaleV_Helper<-3,xx,ix,T,M1d>::call2(len,x,it);
                it.shiftP(-colstep);
                if (len < N) ++len;
            }
            TMVAssert(len == TMV_MIN(M,N));
            const int ds = IntTraits2<cs,rs>::min;
            ScaleV_Helper<-3,ds,ix,T,M1d>::call2(len,x,it);
            for(int k=1;k<=m.nhi();++k) {
                it.shiftP(rowstep);
                if (k+len > N) --len;
                ScaleV_Helper<-3,xx,ix,T,M1d>::call2(len,x,it);
            }
        }
    };

    // algo 90: Call inst
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleB_Helper<90,cs,rs,ix,T,M1>
    {
        static TMV_INLINE void call(const Scaling<ix,T>& x, M1& m)
        {
            typedef typename M1::value_type VT;
            VT xx = Traits<VT>::convert(T(x));
            InstScale(xx,m.xView());
        }
    };

    // algo 97: Conjugate
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleB_Helper<97,cs,rs,ix,T,M1>
    {
        static TMV_INLINE void call(const Scaling<ix,T>& x, M1& m)
        {
            typedef typename M1::conjugate_type Mc;
            Mc mc = m.conjugate();
            ScaleB_Helper<-2,cs,rs,ix,T,Mc>::call(TMV_CONJ(x),mc);
        }
    };

    // algo -4: No copies or branches
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleB_Helper<-4,cs,rs,ix,T,M1>
    {
        static TMV_INLINE void call(const Scaling<ix,T>& x, M1& m)
        {
            typedef typename M1::value_type T1;
            const int algo = 
                (cs == 0 || rs == 0) ? 0 :
                (ix == 1) ? 1 :
                M1::_canlin ? 2 :
                TMV_OPT == 0 ? 13 :
                M1::_colmajor ? 11 :
                M1::_rowmajor ? 12 :
                13;
#ifdef PRINTALGO_XB
            std::cout<<"ScaleB: algo = "<<algo<<std::endl;
            std::cout<<"x = "<<ix<<" "<<T(x)<<std::endl;
            std::cout<<"m = "<<TMV_Text(m)<<std::endl;
            //std::cout<<"m = "<<m<<std::endl;
#endif
            ScaleB_Helper<algo,cs,rs,ix,T,M1>::call(x,m);
#ifdef PRINTALGO_XB
            //std::cout<<"m => "<<m<<std::endl;
#endif
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleB_Helper<-3,cs,rs,ix,T,M1>
    {
        static TMV_INLINE void call(const Scaling<ix,T>& x, M1& m)
        { ScaleB_Helper<-4,cs,rs,ix,T,M1>::call(x,m); }
    };

    // algo -2: Check for inst
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleB_Helper<-2,cs,rs,ix,T,M1>
    {
        static TMV_INLINE void call(const Scaling<ix,T>& x, M1& m)
        {
            typedef typename M1::value_type T1;
            const bool inst =
                (cs == Unknown || cs > 16) &&
                (rs == Unknown || rs > 16) &&
                Traits<T1>::isinst;
            const int algo = 
                ix == 1 ? 0 :
                M1::_conj ? 97 :
                inst ? 90 : 
                -3;
            ScaleB_Helper<algo,cs,rs,ix,T,M1>::call(x,m);
        }
    };

    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleB_Helper<-1,cs,rs,ix,T,M1>
    {
        static TMV_INLINE void call(const Scaling<ix,T>& x, M1& m)
        { ScaleB_Helper<-2,cs,rs,ix,T,M1>::call(x,m); }
    };

    template <int ix, class T, class M>
    inline void Scale(
        const Scaling<ix,T>& x, BaseMatrix_Band_Mutable<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        ScaleB_Helper<-2,cs,rs,ix,T,Mv>::call(x,mv);
    }

    template <int ix, class T, class M>
    inline void InlineScale(
        const Scaling<ix,T>& x, BaseMatrix_Band_Mutable<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        ScaleB_Helper<-3,cs,rs,ix,T,Mv>::call(x,mv);
    }

} // namespace tmv

#endif
