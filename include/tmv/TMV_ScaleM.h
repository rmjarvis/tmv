
//
// This file defines the basic composite type for a product of a 
// matrix and a scalar.  It also implements the calculation for 
// dense rectangular matrices.

#ifndef TMV_ScaleM_H
#define TMV_ScaleM_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_ScaleV.h"

namespace tmv {

    // Defined in TMV_ScaleM.cpp
    template <class T>
    void InstScale(const T x, MatrixView<T> m);

    //
    // Matrix *= x
    //

    template <int algo, int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper;

    // algo 0: trivial: cs == 0, rs == 0 or ix == 1, so nothing to do
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper<0,cs,rs,ix,T,M1>
    { static TMV_INLINE void call(const Scaling<ix,T>& , M1& ) {} };

    // algo 1: Linearize to vector version
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper<1,cs,rs,ix,T,M1>
    {
        static TMV_INLINE void call(const Scaling<ix,T>& x, M1& m)
        {
            typedef typename M1::linearview_type Ml;
            Ml ml = m.linearView();
            ScaleV_Helper<-3,Ml::_size,ix,T,Ml>::call(x,ml);
        }
    };

    // algo 11: Loop over columns
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper<11,cs,rs,ix,T,M1>
    {
        static void call(const Scaling<ix,T>& x, M1& m)
        {
            const int M = cs == TMV_UNKNOWN ? m.colsize() : cs;
            int N = rs == TMV_UNKNOWN ? m.rowsize() : rs;
            typedef typename M1::col_type Mc;
            typedef typename Mc::iterator IT;
            const int step = m.stepj();
            IT it = m.get_col(0).begin();
            for(;N;--N) {
                ScaleV_Helper<-3,cs,ix,T,Mc>::call2(M,x,it);
                it.shiftP(step);
            }
        }
    };

    // algo 21: Loop over rows
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper<21,cs,rs,ix,T,M1>
    {
        static void call(const Scaling<ix,T>& x, M1& m)
        {
            int M = cs == TMV_UNKNOWN ? m.colsize() : cs;
            const int N = rs == TMV_UNKNOWN ? m.rowsize() : rs;
            typedef typename M1::row_type Mr;
            typedef typename Mr::iterator IT;
            const int step = m.stepi();
            IT it = m.get_row(0).begin();
            for(;M;--M) {
                ScaleV_Helper<-3,rs,ix,T,Mr>::call2(N,x,it);
                it.shiftP(step);
            }
        }
    };

    // algo 31: Unknown sizes, determine which algorithm to use
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper<31,cs,rs,ix,T,M1>
    {
        static void call(const Scaling<ix,T>& x, M1& m)
        {
            if (m.canLinearize()) 
                ScaleM_Helper<1,cs,rs,ix,T,M1>::call(x,m);
            else if ( m.isrm() ||
                      (!m.iscm() && (m.colsize() > m.rowsize()) ) )
                ScaleM_Helper<21,cs,rs,ix,T,M1>::call(x,m);
            else
                ScaleM_Helper<11,cs,rs,ix,T,M1>::call(x,m);
        }
    };

    // algo 25: Fully unroll by rows
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper<25,cs,rs,ix,T,M1>
    {
        template <int I, int M, int J, int N, bool iscomplex>
        struct Unroller
        {
            static TMV_INLINE void unroll(const Scaling<ix,T>& x, M1& m)
            {
                Unroller<I,M/2,J,N,iscomplex>::unroll(x,m);
                Unroller<I+M/2,M-M/2,J,N,iscomplex>::unroll(x,m);
            }
        };
        template <int I, int J, int N, bool iscomplex>
        struct Unroller<I,1,J,N,iscomplex>
        {
            static TMV_INLINE void unroll(const Scaling<ix,T>& x, M1& m)
            {
                Unroller<I,1,J,N/2,iscomplex>::unroll(x,m);
                Unroller<I,1,J+N/2,N-N/2,iscomplex>::unroll(x,m);
            }
        };
        template <int I, int J, int N, bool iscomplex>
        struct Unroller<I,0,J,N,iscomplex>
        { static TMV_INLINE void unroll(const Scaling<ix,T>& , M1& ) {} };
        template <int I, int J>
        struct Unroller<I,1,J,1,false>
        {
            static TMV_INLINE void unroll(const Scaling<ix,T>& x, M1& m)
            { m.ref(I,J) *= x; }
        };
        template <int I, int J>
        struct Unroller<I,1,J,1,true>
        {
            static TMV_INLINE void unroll(const Scaling<ix,T>& x, M1& m)
            {
                typedef typename M1::real_type RT;
                typedef typename M1::value_type VT;
                const RT rm = 
                    ZProd<false,false>::rprod(x,m.cref(I,J));
                const RT im = 
                    ZProd<false,false>::iprod(x,m.cref(I,J));
                m.ref(I,J) = VT(rm,im);
            }
        };
        template <int I, int J, bool iscomplex>
        struct Unroller<I,1,J,0,iscomplex>
        { static TMV_INLINE void unroll(const Scaling<ix,T>& , M1& ) {} };

        static inline void call(const Scaling<ix,T>& x, M1& m)
        { Unroller<0,cs,0,rs,M1::iscomplex>::unroll(x,m); }
    };

    // algo 15: Fully unroll by columns
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper<15,cs,rs,ix,T,M1>
    {
        template <int I, int M, int J, int N, bool iscomplex>
        struct Unroller
        {
            static TMV_INLINE void unroll(const Scaling<ix,T>& x, M1& m)
            {
                Unroller<I,M,J,N/2,iscomplex>::unroll(x,m);
                Unroller<I,M,J+N/2,N-N/2,iscomplex>::unroll(x,m);
            }
        };
        template <int I, int M, int J, bool iscomplex>
        struct Unroller<I,M,J,1,iscomplex>
        {
            static TMV_INLINE void unroll(const Scaling<ix,T>& x, M1& m)
            {
                Unroller<I,M/2,J,1,iscomplex>::unroll(x,m);
                Unroller<I+M/2,M-M/2,J,1,iscomplex>::unroll(x,m);
            }
        };
        template <int I, int M, int J, bool iscomplex>
        struct Unroller<I,M,J,0,iscomplex>
        { static TMV_INLINE void unroll(const Scaling<ix,T>& , M1& ) {} };
        template <int I, int J>
        struct Unroller<I,1,J,1,false>
        {
            static TMV_INLINE void unroll(const Scaling<ix,T>& x, M1& m)
            { m.ref(I,J) *= x; }
        };
        template <int I, int J>
        struct Unroller<I,1,J,1,true>
        {
            static TMV_INLINE void unroll(const Scaling<ix,T>& x, M1& m)
            {
                typedef typename M1::real_type RT;
                typedef typename M1::value_type VT;
                const RT rm =
                    ZProd<false,false>::rprod(x,m.cref(I,J));
                const RT im =
                    ZProd<false,false>::iprod(x,m.cref(I,J));
                m.ref(I,J) = VT(rm,im);
            }
        };
        template <int I, int J, bool iscomplex>
        struct Unroller<I,0,J,1,iscomplex>
        { static TMV_INLINE void unroll(const Scaling<ix,T>& , M1& ) {} };

        static inline void call(const Scaling<ix,T>& x, M1& m)
        { Unroller<0,cs,0,rs,M1::iscomplex>::unroll(x,m); }
    };

    // algo 90: Call inst
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper<90,cs,rs,ix,T,M1>
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
    struct ScaleM_Helper<97,cs,rs,ix,T,M1>
    {
        static TMV_INLINE void call(const Scaling<ix,T>& x, M1& m)
        {
            typedef typename M1::conjugate_type Mc;
            Mc mc = m.conjugate();
            ScaleM_Helper<-2,cs,rs,ix,T,Mc>::call(TMV_CONJ(x),mc);
        }
    };

    // algo -4: No copies or branches
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper<-4,cs,rs,ix,T,M1>
    {
        static TMV_INLINE void call(const Scaling<ix,T>& x, M1& m)
        {
            typedef typename M1::value_type T1;
            const int algo = 
                (cs == 0 || rs == 0) ? 0 :
                (ix == 1) ? 1 :
                M1::_canlin ? 2 :
                TMV_OPT == 0 ? 11 :
                ( cs != TMV_UNKNOWN && rs != TMV_UNKNOWN ) ? (
                    ( IntTraits2<cs,rs>::prod <= int(128/sizeof(T1)) ) ? (
                        ( M1::_rowmajor ? 25 : 15 ) ) :
                    M1::_rowmajor ? 21 :
                    M1::_colmajor ? 11 :
                    ( cs > rs ) ? 21 : 11 ) :
                M1::_rowmajor ? 21 :
                11;
            ScaleM_Helper<algo,cs,rs,ix,T,M1>::call(x,m);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper<-3,cs,rs,ix,T,M1>
    {
        static TMV_INLINE void call(const Scaling<ix,T>& x, M1& m)
        {
            const int algo = 
                (cs == 0 || rs == 0 || ix == 1) ? 0 :
                M1::_canlin ? 1 :
                TMV_OPT >= 2 && (cs == TMV_UNKNOWN || rs == TMV_UNKNOWN) ? 31 :
                -4;
            ScaleM_Helper<algo,cs,rs,ix,T,M1>::call(x,m);
        }
    };

    // algo -2: Check for inst
    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper<-2,cs,rs,ix,T,M1>
    {
        static TMV_INLINE void call(const Scaling<ix,T>& x, M1& m)
        {
            typedef typename M1::value_type T1;
            const bool inst =
                (cs == TMV_UNKNOWN || cs > 16) &&
                (rs == TMV_UNKNOWN || rs > 16) &&
                Traits<T1>::isinst;
            const int algo = 
                ix == 1 ? 0 :
                M1::_conj ? 97 :
                inst ? 90 : 
                -3;
            ScaleM_Helper<algo,cs,rs,ix,T,M1>::call(x,m);
        }
    };

    template <int cs, int rs, int ix, class T, class M1>
    struct ScaleM_Helper<-1,cs,rs,ix,T,M1>
    {
        static TMV_INLINE void call(const Scaling<ix,T>& x, M1& m)
        { ScaleM_Helper<-2,cs,rs,ix,T,M1>::call(x,m); }
    };

    template <int ix, class T, class M>
    inline void Scale(
        const Scaling<ix,T>& x, BaseMatrix_Rec_Mutable<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        ScaleM_Helper<-2,cs,rs,ix,T,Mv>::call(x,mv);
    }

    template <int ix, class T, class M>
    inline void InlineScale(
        const Scaling<ix,T>& x, BaseMatrix_Rec_Mutable<M>& m)
    {
        const int cs = M::_colsize;
        const int rs = M::_rowsize;
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        ScaleM_Helper<-3,cs,rs,ix,T,Mv>::call(x,mv);
    }

} // namespace tmv

#endif
