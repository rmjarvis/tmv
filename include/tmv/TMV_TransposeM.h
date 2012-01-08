
#ifndef TMV_TransposeM_H
#define TMV_TransposeM_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_SwapV.h"

namespace tmv {

    // Defined in TMV_Matrix.cpp
    template <class T>
    void InstTransposeSelf(MatrixView<T> m);

    //
    // TransposeSelf
    //

// UNROLL is the maximum nops to unroll.
#if TMV_OPT >= 3
#define TMV_TRANSPOSEM_UNROLL 200 
#elif TMV_OPT >= 2
#define TMV_TRANSPOSEM_UNROLL 25
#elif TMV_OPT >= 1
#define TMV_TRANSPOSEM_UNROLL 9
#else
#define TMV_TRANSPOSEM_UNROLL 0
#endif

    template <int algo, int s, class M1>
    struct TransposeSelf_Helper;

    // algo 0: s == 0 or s == 1, nothing to do.
    template <int s, class M1>
    struct TransposeSelf_Helper<0,s,M1>
    { static TMV_INLINE void call(M1& m) {} };

    // algo 11: Simple for loop
    template <int s, class M1>
    struct TransposeSelf_Helper<11,s,M1>
    {
        static void call(M1& m)
        {
            const int n = (s == Unknown ? m.colsize() : s);
            for(int i=1;i<n;++i) {
                typename M1::row_sub_type::noalias_type v1 =
                    m.get_row(i,0,i).noAlias();
                typename M1::col_sub_type::noalias_type v2 =
                    m.get_col(i,0,i).noAlias();
                Swap(v1,v2);
            }
        }
    };

    // algo 12: Same as algo 11, but with iterators.
    template <int s, class M1>
    struct TransposeSelf_Helper<12,s,M1>
    {
        static void call(M1& m)
        {
            const int n = s == Unknown ? m.colsize() : s;
            if (n <= 1) return;
            typedef typename M1::row_type Mr;
            typedef typename M1::col_type Mc;
            typedef typename Mr::iterator IT1;
            typedef typename Mc::iterator IT2;
            IT1 it1 = m.row(1).begin();
            IT2 it2 = m.col(1).begin();
            const int step1 = m.stepi();
            const int step2 = m.stepj();
            for(int i=1;i<n;++i) {
                SwapV_Helper<-3,Unknown,Mr,Mc>::call2(i,it1,it2);
                it1.shiftP(step1);
                it2.shiftP(step2);
            }
        }
    };

    // algo 13: The other way to do the loop.  
    // This way seems to be a little bit slower.
    template <int s, class M1>
    struct TransposeSelf_Helper<13,s,M1>
    {
        static void call(M1& m)
        {
            const int n = (s == Unknown ? m.colsize() : s);
            for(int i=0;i<n-1;++i) {
                typename M1::row_sub_type::noalias_type v1 =
                    m.get_row(i,i+1,n).noAlias();
                typename M1::col_sub_type::noalias_type v2 =
                    m.get_col(i,i+1,n).noAlias();
                Swap(v1,v2);
            }
        }
    };

    // algo 14: Same as algo 13, but with iterators.
    template <int s, class M1>
    struct TransposeSelf_Helper<14,s,M1>
    {
        static void call(M1& m)
        {
            int n = s == Unknown ? m.colsize() : s;
            if (n <= 1) return;
            typedef typename M1::row_type Mr;
            typedef typename M1::col_type Mc;
            typedef typename Mr::iterator IT1;
            typedef typename Mc::iterator IT2;
            IT1 it1 = m.row(0).begin(); ++it1;
            IT2 it2 = m.col(0).begin(); ++it2;
            const int step = m.stepi() + m.stepj();
            for(--n;n;--n) {
                SwapV_Helper<-3,Unknown,Mr,Mc>::call2(n,it1,it2);
                it1.shiftP(step);
                it2.shiftP(step);
            }
        }
    };

    // algo 15: Fully unroll
    template <int s, class M1>
    struct TransposeSelf_Helper<15,s,M1>
    {
        template <int I, int M, int J, int N>
        struct Unroller
        {
            static TMV_INLINE void unroll(M1& m)
            {
                Unroller<I,M/2,0,I>::unroll(m);
                Unroller<I+M/2,M-M/2,0,I+M/2>::unroll(m);
            }
        };
        template <int I, int J, int N>
        struct Unroller<I,1,J,N>
        {
            static TMV_INLINE void unroll(M1& m)
            {
                Unroller<I,1,J,N/2>::unroll(m);
                Unroller<I,1,J+N/2,N-N/2>::unroll(m);
            }
        };
        template <int I, int J, int N>
        struct Unroller<I,0,J,N>
        { static TMV_INLINE void unroll(M1& ) {} };
        template <int I, int J>
        struct Unroller<I,1,J,1>
        {
            static TMV_INLINE void unroll(M1& m)
            { TMV_SWAP(m.ref(I,J) , m.ref(J,I) ); }
        };
        template <int I>
        struct Unroller<I,1,I,1>
        { static TMV_INLINE void unroll(M1& m) {} };
        template <int I, int J>
        struct Unroller<I,1,J,0>
        { static TMV_INLINE void unroll(M1& ) {} };

        static inline void call(M1& m)
        { Unroller<0,s,0,s>::unroll(m); }
    };

    // algo 90: Call inst
    template <int s, class M1>
    struct TransposeSelf_Helper<90,s,M1>
    {
        static TMV_INLINE void call(M1& m)
        { InstTransposeSelf(m.xView()); }
    };

    // algo 97: Conjugate
    template <int s, class M1>
    struct TransposeSelf_Helper<97,s,M1>
    {
        static TMV_INLINE void call(M1& m)
        {
            typedef typename M1::nonconj_type Mnc;
            Mnc mnc = m.nonConj();
            TransposeSelf_Helper<-2,s,Mnc>::call(mnc);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int s, class M1>
    struct TransposeSelf_Helper<-3,s,M1>
    {
        static TMV_INLINE void call(M1& m)
        {
            TMVStaticAssert(!M1::_conj);
            // nops = n(n-1)/2
            const int sm1 = IntTraits<s>::Sm1;
            const int nops = IntTraits2<s,sm1>::safeprod / 2;
            const bool unroll = 
                s >= 8 ? false :
                s == Unknown ? false :
                nops <= TMV_TRANSPOSEM_UNROLL;
            const int algo = 
                s == 0 || s == 1 ? 0 :
                unroll ? 15 :
                12;
            TransposeSelf_Helper<algo,s,M1>::call(m);
        }
    };

    // algo -2: Check for inst
    template <int s, class M1>
    struct TransposeSelf_Helper<-2,s,M1>
    {
        static TMV_INLINE void call(M1& m)
        {
            typedef typename M1::value_type T;
            const bool inst = 
                (s == Unknown || s > 16) &&
                Traits<T>::isinst;
            const int algo = 
                M1::_conj ? 97 :
                inst ? 90 :
                -3;
            TransposeSelf_Helper<algo,s,M1>::call(m);
        }
    };

    template <int s, class M1>
    struct TransposeSelf_Helper<-1,s,M1>
    {
        static TMV_INLINE void call(M1& m)
        { TransposeSelf_Helper<-2,s,M1>::call(m); }
    };

    template <class M>
    inline void TransposeSelf(BaseMatrix_Rec_Mutable<M>& m)
    {
        TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same)); 
        TMVAssert(m.colsize() == m.rowsize());
        const int s = Sizes<M::_colsize,M::_rowsize>::size;
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        TransposeSelf_Helper<-2,s,Mv>::call(mv);
    }

    template <class M>
    inline void InlineTransposeSelf(BaseMatrix_Rec_Mutable<M>& m)
    {
        TMVStaticAssert((Sizes<M::_colsize,M::_rowsize>::same)); 
        TMVAssert(m.colsize() == m.rowsize());
        const int s = Sizes<M::_colsize,M::_rowsize>::size;
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        TransposeSelf_Helper<-3,s,Mv>::call(mv);
    }

} // namespace tmv

#endif
