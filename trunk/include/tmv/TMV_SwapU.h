
#ifndef TMV_SwapU_H
#define TMV_SwapU_H

#include "TMV_BaseMatrix_Tri.h"
#include "TMV_SwapV.h"

namespace tmv {

    // Defined in TMV_TriMatrix.cpp
    template <class T, int C1>
    void InstSwap(UpperTriMatrixView<T,C1> m1, UpperTriMatrixView<T> m2); 
    template <class T, int C1>
    void InstAliasSwap(UpperTriMatrixView<T,C1> m1, UpperTriMatrixView<T> m2); 

    //
    // Swap Matrices
    //

    // The maximum nops to unroll.
#if TMV_OPT >= 3
#define TMV_SWAPU_UNROLL 200
#elif TMV_OPT >= 2
#define TMV_SWAPU_UNROLL 25
#elif TMV_OPT >= 1
#define TMV_SWAPU_UNROLL 9
#else
#define TMV_SWAPU_UNROLL 0
#endif

    template <int algo, int s, class M1, class M2>
    struct SwapU_Helper;

    // algo 1: m1,m2 are unitdiag
    template <int s, class M1, class M2>
    struct SwapU_Helper<1,s,M1,M2>
    {
        static inline void call(M1& m1, M2& m2)
        {
            if (m1.size() > 1) {
                typedef typename M1::offdiag_type M1o;
                typedef typename M2::offdiag_type M2o;
                M1o m1o = m1.offDiag();
                M2o m2o = m2.offDiag();
                const int sm1 = IntTraits2<s,-1>::sum;
                SwapU_Helper<-3,sm1,M1o,M2o>::call(m1o,m2o);
            }
        }
    };

    // algo 2: UnknownDiag
    template <int s, class M1, class M2>
    struct SwapU_Helper<2,s,M1,M2>
    {
        static void call(M1& m1, M2& m2)
        {
            const int s2 = s > 20 ? Unknown : s;
            const int s2p1 = IntTraits<s2>::Sp1;
            const int nops = IntTraits2<s2,s2p1>::safeprod / 2;
            const bool unroll = 
                s > 10 ? false :
                s == Unknown ? false :
                nops <= TMV_SWAPU_UNROLL;
            const int algo2 = 
                unroll ? ( M2::_rowmajor ? 16 : 15 ) :
                M2::_rowmajor ? 12 : 11;
            if (m1.isunit()) 
                SwapU_Helper<1,s,M1,M2>::call(m1,m2);
            else
                SwapU_Helper<algo2,s,M1,M2>::call(m1,m2);
        }       
    };          

    // algo 3: Transpose (and go back to -3, rather than -2)
    template <int s, class M1, class M2>
    struct SwapU_Helper<3,s,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        {
            typedef typename M1::transpose_type M1t;
            typedef typename M2::transpose_type M2t;
            M1t m1t = m1.transpose();
            M2t m2t = m2.transpose();
            SwapU_Helper<-3,s,M1t,M2t>::call(m1t,m2t);
        }       
    };          

    // algo 11: Loop over columns
    template <int s, class M1, class M2>
    struct SwapU_Helper<11,s,M1,M2>
    {
        static void call(M1& m1, M2& m2)
        {
            const int N = (s == Unknown ? m2.size() : s);
            typedef typename M1::col_sub_type M1c;
            typedef typename M2::col_sub_type M2c;
            typedef typename M1c::iterator IT1;
            typedef typename M2c::iterator IT2;
            const int step1 = m1.stepj();
            const int step2 = m2.stepj();
            IT1 it1 = m1.get_col(0,0,1).begin();
            IT2 it2 = m2.get_col(0,0,1).begin();
            for(int j=0;j<N;++j) {
                SwapV_Helper<-3,Unknown,M1c,M2c>::call2(j+1,it1,it2);
                it1.shiftP(step1);
                it2.shiftP(step2);
            } 
        }
    };

    // algo 12: Loop over rows
    template <int s, class M1, class M2>
    struct SwapU_Helper<12,s,M1,M2>
    {
        static void call(M1& m1, M2& m2)
        {
            int N = (s == Unknown ? m2.size() : s);
            typedef typename M1::row_sub_type M1r;
            typedef typename M2::row_sub_type M2r;
            typedef typename M1r::iterator IT1;
            typedef typename M2r::iterator IT2;
            const int step1 = m1.diagstep();
            const int step2 = m2.diagstep();
            IT1 it1 = m1.get_row(0,0,N).begin();
            IT2 it2 = m2.get_row(0,0,N).begin();
            for(;N;--N) {
                SwapV_Helper<-3,Unknown,M1r,M2r>::call2(N,it1,it2);
                it1.shiftP(step1);
                it2.shiftP(step2);
            } 
        }
    };

    // algo 15: Fully unroll by columns
    template <int s, class M1, class M2>
    struct SwapU_Helper<15,s,M1,M2>
    {
        template <int I, int M, int J, int N>
        struct Unroller
        {
            static TMV_INLINE void unroll(M1& m1, M2& m2)
            {
                Unroller<I,M-(N-N/2),J,N/2>::unroll(m1,m2);
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
        { Unroller<0,s,0,s>::unroll(m1,m2); }
    };

    // algo 16: Fully unroll by rows
    template <int s, class M1, class M2>
    struct SwapU_Helper<16,s,M1,M2>
    {
        template <int I, int M, int J, int N>
        struct Unroller
        {
            static TMV_INLINE void unroll(M1& m1, M2& m2)
            {
                Unroller<I,M/2,J,N>::unroll(m1,m2);
                Unroller<I+M/2,M-M/2,J+M/2,N-M/2>::unroll(m1,m2);
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
        { Unroller<0,s,0,s>::unroll(m1,m2); }
    };

    // algo 90: Call inst
    template <int s, class M1, class M2>
    struct SwapU_Helper<90,s,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        { InstSwap(m1.xView(),m2.xView()); }
    };

    // algo 91: Call inst alias
    template <int s, class M1, class M2>
    struct SwapU_Helper<91,s,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        { InstAliasSwap(m1.xView(),m2.xView()); }
    };

    // algo 96: Transpose
    template <int s, class M1, class M2>
    struct SwapU_Helper<96,s,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        {
            typedef typename M1::transpose_type M1t;
            typedef typename M2::transpose_type M2t;
            M1t m1t = m1.transpose();
            M2t m2t = m2.transpose();
            SwapU_Helper<-2,s,M1t,M2t>::call(m1t,m2t);
        }       
    };          

    // algo 196: Transpose
    template <int s, class M1, class M2>
    struct SwapU_Helper<196,s,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        {
            typedef typename M1::transpose_type M1t;
            typedef typename M2::transpose_type M2t;
            M1t m1t = m1.transpose();
            M2t m2t = m2.transpose();
            SwapU_Helper<99,s,M1t,M2t>::call(m1t,m2t);
        }       
    };          

    // algo 97: Conjugate
    template <int s, class M1, class M2>
    struct SwapU_Helper<97,s,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        {
            typedef typename M1::conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            SwapU_Helper<-2,s,M1c,M2c>::call(m1c,m2c);
        }       
    };          

    // algo 197: Conjugate
    template <int s, class M1, class M2>
    struct SwapU_Helper<197,s,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        {
            typedef typename M1::conjugate_type M1c;
            typedef typename M2::conjugate_type M2c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            SwapU_Helper<99,s,M1c,M2c>::call(m1c,m2c);
        }       
    };          

    // algo 98: Inline check for aliases
    template <int s, class M1, class M2>
    struct SwapU_Helper<98,s,M1,M2>
    {
        static void call(M1& m1, M2& m2)
        {
            if (!SameStorage(m1,m2)) {
                // No aliasing
                SwapU_Helper<-2,s,M1,M2>::call(m1,m2);
            } else if (ExactSameStorage(m1,m2)) {
                // They are equal modulo a conjugation
                Maybe<M1::_conj != int(M2::_conj)>::conjself(m1);
            } else if (OppositeStorage(m1,m2)) {
                // Then only need to swap the offdiag
                if (m2.size() > 1) {
                    typedef typename M1::offdiag_type M1o;
                    typedef typename M2::offdiag_type M2o;
                    M1o m1o = m1.offDiag();
                    M2o m2o = m2.offDiag();
                    const int sm1 = IntTraits<s>::Sm1;
                    SwapU_Helper<-2,sm1,M1o,M2o>::call(m1o,m2o);
                }
                // And maybe conjugate the diag
                typename M2::diag_type m2d = m2.diag();
                Maybe<M1::_conj != int(M2::_conj)>::conjself(m2d);
            } else {
                // Need a temporary
                typename M1::copy_type m1c = m1;
                m1.noAlias() = m2;
                m2.noAlias() = m1c;
            }
        }
    };

    // algo 99: Check for aliases
    template <int s, class M1, class M2>
    struct SwapU_Helper<99,s,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst =
                (s == Unknown || s > 16) &&
                Traits2<T1,T2>::sametype &&
                Traits<T1>::isinst;
            const int algo =
                M2::_lower ? 196 :
                M2::_conj ? 197 :
                inst ? 91 :
                98;
            SwapU_Helper<algo,s,M1,M2>::call(m1,m2);
        }
    };

    // algo -4: No branches or copies, and Upper, NonUnitDiag
    template <int s, class M1, class M2>
    struct SwapU_Helper<-4,s,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        {
            TMVStaticAssert(!M1::_unit);
            TMVStaticAssert(!M2::_unit);
            TMVStaticAssert(M1::_upper);
            TMVStaticAssert(M2::_upper);
            typedef typename M2::value_type T2;
            const int s2 = s > 20 ? Unknown : s;
            const int s2p1 = IntTraits<s2>::Sp1;
            // nops = n(n+1)/2
            const int nops = IntTraits2<s2,s2p1>::safeprod / 2;
            const bool unroll = 
                s == Unknown ? false :
                nops > TMV_SWAPU_UNROLL ? false :
                s <= 10;
            const int algo = 
                unroll ? ( M2::_rowmajor ? 16 : 15 ) :
                M2::_rowmajor ? 12 : 11;
            SwapU_Helper<algo,s,M1,M2>::call(m1,m2);
        }
    };

    // algo -3: Determine which algorithm to use
    template <int s, class M1, class M2>
    struct SwapU_Helper<-3,s,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        {
            TMVStaticAssert(M1::_unit == int(M2::_unit));
            TMVStaticAssert(M1::_upper == int(M2::_upper));
            const int algo = 
                M1::_unit || M2::_unit ? 1 :
                M1::_lower ? 3 :
                M1::_unknowndiag && M2::_unknowndiag ? 2 :
                -4;
            SwapU_Helper<algo,s,M1,M2>::call(m1,m2);
        }
    };

    // algo -2: Check for inst
    template <int s, class M1, class M2>
    struct SwapU_Helper<-2,s,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        {
            typedef typename M1::value_type T1;
            typedef typename M2::value_type T2;
            const bool inst =
                (s == Unknown || s > 16) &&
                Traits2<T1,T2>::sametype &&
                Traits<T1>::isinst;
            const int algo =
                M2::_lower ? 96 :
                M2::_conj ? 97 :
                inst ? 90 :
                -3;
            SwapU_Helper<algo,s,M1,M2>::call(m1,m2);
        }
    };

    // algo -1: Check for aliases?
    template <int s, class M1, class M2>
    struct SwapU_Helper<-1,s,M1,M2>
    {
        static TMV_INLINE void call(M1& m1, M2& m2)
        {
            const bool noclobber =
                MStepHelper<M1,M2>::opp && M1::_unit;
            const bool checkalias =
                (M1::_checkalias || M2::_checkalias) && !noclobber;
            const int algo =
                checkalias ? 99 :
                -2;
            SwapU_Helper<algo,s,M1,M2>::call(m1,m2);
        }
    };


    template <class M1, class M2>
    inline void DoSwap(
        BaseMatrix_Tri_Mutable<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
    {
        TMVStaticAssert(M1::_upper == int(M2::_upper));
        TMVStaticAssert(M1::_unit == int(M2::_unit));
        TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.isunit() == m2.isunit());
        const int s = Sizes<M1::_size,M2::_size>::size;
        typedef typename M1::cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_REF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        SwapU_Helper<-1,s,M1v,M2v>::call(m1v,m2v);
    }

    template <class M1, class M2>
    inline void InlineSwap(
        BaseMatrix_Tri_Mutable<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
    {
        TMVStaticAssert(M1::_upper == int(M2::_upper));
        TMVStaticAssert(M1::_unit == int(M2::_unit));
        TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.isunit() == m2.isunit());
        const int s = Sizes<M1::_size,M2::_size>::size;
        typedef typename M1::cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_REF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        SwapU_Helper<-3,s,M1v,M2v>::call(m1v,m2v);
    }

    template <class M1, class M2>
    inline void InlineAliasSwap(
        BaseMatrix_Tri_Mutable<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
    {
        TMVStaticAssert(M1::_upper == int(M2::_upper));
        TMVStaticAssert(M1::_unit == int(M2::_unit));
        TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.isunit() == m2.isunit());
        const int s = Sizes<M1::_size,M2::_size>::size;
        typedef typename M1::cview_type M1v;
        typedef typename M2::cview_type M2v;
        TMV_MAYBE_REF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_REF(M2,M2v) m2v = m2.cView();
        SwapU_Helper<98,s,M1v,M2v>::call(m1v,m2v);
    }

    template <class M1, class M2>
    TMV_INLINE void Swap(
        BaseMatrix_Tri_Mutable<M1>& m1, BaseMatrix_Tri_Mutable<M2>& m2)
    { DoSwap(m1,m2); }

} // namespace tmv

#endif
