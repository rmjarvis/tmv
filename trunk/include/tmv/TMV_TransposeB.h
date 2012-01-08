
#ifndef TMV_TransposeB_H
#define TMV_TransposeB_H

#include "TMV_BaseMatrix_Band.h"
#include "TMV_SwapV.h"

namespace tmv {

    // Defined in TMV_BandMatrix.cpp
    template <class T>
    void InstTransposeSelf(BandMatrixView<T> m);

    //
    // TransposeSelf
    //

    template <int algo, int hi, class M1>
    struct TransposeSelfB_Helper;

    // algo 0: nhi == 0, nothing to do.
    template <int hi, class M1>
    struct TransposeSelfB_Helper<0,hi,M1>
    { static TMV_INLINE void call(M1& m) {} };

    // algo 11: Simple for loop
    template <int hi, class M1>
    struct TransposeSelfB_Helper<11,hi,M1>
    {
        static void call(M1& m)
        {
            const int nhi = (hi == Unknown ? m.nhi() : hi);
            for(int k=1;k<=nhi;++k) {
                typename M1::diag_sub_type::noalias_type v1 = 
                    m.get_diag(-k).noAlias();
                typename M1::diag_sub_type::noalias_type v2 =
                    m.get_diag(k).noAlias();
                Swap(v1,v2);
            }
        }
    };

    // algo 90: Call inst
    template <int hi, class M1>
    struct TransposeSelfB_Helper<90,hi,M1>
    {
        static TMV_INLINE void call(M1& m)
        { InstTransposeSelf(m.xView()); }
    };

    // algo 97: Conjugate
    template <int hi, class M1>
    struct TransposeSelfB_Helper<97,hi,M1>
    {
        static TMV_INLINE void call(M1& m)
        {
            typedef typename M1::nonconj_type Mnc;
            Mnc mnc = m.nonConj();
            TransposeSelfB_Helper<-2,hi,Mnc>::call(mnc);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int hi, class M1>
    struct TransposeSelfB_Helper<-3,hi,M1>
    {
        static TMV_INLINE void call(M1& m)
        {
            TMVStaticAssert(!M1::_conj);
            const int algo = 
                hi == 0 ? 0 :
                11;
            TransposeSelfB_Helper<algo,hi,M1>::call(m);
        }
    };

    // algo -2: Check for inst
    template <int hi, class M1>
    struct TransposeSelfB_Helper<-2,hi,M1>
    {
        static TMV_INLINE void call(M1& m)
        {
            typedef typename M1::value_type T;
            const int s = Sizes<M1::_colsize,M1::_rowsize>::size;
            const bool inst = 
                (s == Unknown || s > 16) &&
                Traits<T>::isinst;
            const int algo = 
                M1::_conj ? 97 :
                inst ? 90 :
                -3;
            TransposeSelfB_Helper<algo,hi,M1>::call(m);
        }
    };

    template <int hi, class M1>
    struct TransposeSelfB_Helper<-1,hi,M1>
    {
        static TMV_INLINE void call(M1& m)
        { TransposeSelfB_Helper<-2,hi,M1>::call(m); }
    };

    template <class M>
    inline void TransposeSelf(BaseMatrix_Band_Mutable<M>& m)
    {
        TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same)); 
        TMVStaticAssert((Sizes<M::_nlo,M::_nhi>::same)); 
        TMVAssert(m.colsize() == m.rowsize());
        TMVAssert(m.nlo() == m.nhi());
        const int hi = Sizes<M::_nlo,M::_nhi>::size;
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        TransposeSelfB_Helper<-2,hi,Mv>::call(mv);
    }

    template <class M>
    inline void InlineTransposeSelf(BaseMatrix_Band_Mutable<M>& m)
    {
        TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same)); 
        TMVStaticAssert((Sizes<M::_nlo,M::_nhi>::same)); 
        TMVAssert(m.colsize() == m.rowsize());
        TMVAssert(m.nlo() == m.nhi());
        const int hi = Sizes<M::_nlo,M::_nhi>::size;
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        TransposeSelfB_Helper<-3,hi,Mv>::call(mv);
    }

} // namespace tmv

#endif
