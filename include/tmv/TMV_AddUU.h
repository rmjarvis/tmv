

#ifndef TMV_AddUU_H
#define TMV_AddUU_H

#include "TMV_BaseMatrix_Tri.h"
#include "TMV_Scaling.h"
#include "TMV_AddVV.h"
#include "TMV_MultXM_Funcs.h"
#include "TMV_MultXV_Funcs.h"

#ifdef PRINTALGO_AddUU
#include <iostream>
#include "tmv/TMV_MatrixIO.h"
#include "tmv/TMV_TriMatrixIO.h"
#endif

namespace tmv {

    // Defined in TMV_AddUU.cpp
    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddMM(
        const T3 x1, const ConstUpperTriMatrixView<T1,C1>& m1,
        const T3 x2, const ConstUpperTriMatrixView<T2,C2>& m2, 
        UpperTriMatrixView<T3> m3);

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddMM(
        const T3 x1, const ConstUpperTriMatrixView<T1,C1>& m1,
        const T3 x2, const ConstUpperTriMatrixView<T2,C2>& m2, 
        UpperTriMatrixView<T3> m3);

    //
    // U = x * U + x * U
    //

    // The maximum nops to unroll.
#if TMV_OPT >= 3
#define TMV_ADDUU_UNROLL 200 
#elif TMV_OPT >= 2
#define TMV_ADDUU_UNROLL 25
#elif TMV_OPT >= 1
#define TMV_ADDUU_UNROLL 9
#else
#define TMV_ADDUU_UNROLL 0
#endif

    template <int algo, ptrdiff_t s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUU_Helper;

    // algo 2: m1 and/or m2 is unitdiag
    template <ptrdiff_t s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUU_Helper<2,s,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            const ptrdiff_t N = m3.size();
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUU algo 2: N,s = "<<N<<','<<s<<std::endl;
#endif
            if (N > 1) {
                typedef typename M1::const_offdiag_type M1o;
                typedef typename M2::const_offdiag_type M2o;
                typedef typename M3::offdiag_type M3o;
                M1o m1o = m1.offDiag();
                M2o m2o = m2.offDiag();
                M3o m3o = m3.offDiag();
                const ptrdiff_t sm1 = IntTraits2<s,-1>::sum;
                AddUU_Helper<-2,sm1,ix1,T1,M1o,ix2,T2,M2o,M3o>::call(
                    x1,m1o,x2,m2o,m3o);
            }
            typename M3::diag_type m3d = m3.diag();
            if (m1.isunit()) m3d.setAllTo(T1(x1));
            else MultXV<false>(x1,m1.diag(),m3d);
            if (m2.isunit()) m3d.addToAll(T2(x2));
            else MultXV<true>(x2,m2.diag(),m3d);
        }
    };

    // algo 3: UnknownDiag
    template <ptrdiff_t s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUU_Helper<3,s,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUU algo 3: N,s = "<<m3.size()<<','<<s<<std::endl;
#endif
            if (m1.isunit() || m2.isunit())
                AddUU_Helper<2,s,ix1,T1,M1,ix2,T2,M2,M3>::call(
                    x1,m1,x2,m2,m3);
            else
                AddUU_Helper<-4,s,ix1,T1,M1,ix2,T2,M2,M3>::call(
                    x1,m1,x2,m2,m3);
        }
    };


    // algo 11: Loop over columns
    template <ptrdiff_t s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUU_Helper<11,s,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUU algo 11: N,s = "<<m3.size()<<','<<s<<std::endl;
#endif
            ptrdiff_t N = (s == Unknown ? m2.size() : s);
            typedef typename M1::const_col_sub_type M1c;
            typedef typename M2::const_col_sub_type M2c;
            typedef typename M3::col_sub_type M3c;
            typedef typename M1c::const_nonconj_type::const_iterator IT1;
            typedef typename M2c::const_nonconj_type::const_iterator IT2;
            typedef typename M3c::iterator IT3;
            const ptrdiff_t step1 = m1.stepj();
            const ptrdiff_t step2 = m2.stepj();
            const ptrdiff_t step3 = m3.stepj();
            IT1 it1 = m1.get_col(0,0,1).begin().nonConj();
            IT2 it2 = m2.get_col(0,0,1).begin().nonConj();
            IT3 it3 = m3.get_col(0,0,1).begin();
            ptrdiff_t M=1;
            for(;N;--N) {
                AddVV_Helper<-4,Unknown,ix1,T1,M1c,ix2,T2,M2c,M3c>::call2(
                    M++,x1,it1,x2,it2,it3);
                it1.shiftP(step1);
                it2.shiftP(step2);
                it3.shiftP(step3);
            }
        }
    };

    // algo 12: Loop over rows
    template <ptrdiff_t s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUU_Helper<12,s,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUU algo 12: N,s = "<<m3.size()<<','<<s<<std::endl;
#endif
            ptrdiff_t N = (s == Unknown ? m2.size() : s);
            typedef typename M1::const_row_sub_type M1r;
            typedef typename M2::const_row_sub_type M2r;
            typedef typename M3::row_sub_type M3r;
            typedef typename M1r::const_nonconj_type::const_iterator IT1;
            typedef typename M2r::const_nonconj_type::const_iterator IT2;
            typedef typename M3r::iterator IT3;
            const ptrdiff_t step1 = m1.diagstep();
            const ptrdiff_t step2 = m2.diagstep();
            const ptrdiff_t step3 = m3.diagstep();
            IT1 it1 = m1.get_row(0,0,N).begin().nonConj();
            IT2 it2 = m2.get_row(0,0,N).begin().nonConj();
            IT3 it3 = m3.get_row(0,0,N).begin();
            for(;N;--N) {
                AddVV_Helper<-4,Unknown,ix1,T1,M1r,ix2,T2,M2r,M3r>::call2(
                    N,x1,it1,x2,it2,it3);
                it1.shiftP(step1);
                it2.shiftP(step2);
                it3.shiftP(step3);
            }
        }
    };

    // algo 15: Fully unroll by rows
    template <ptrdiff_t s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUU_Helper<15,s,ix1,T1,M1,ix2,T2,M2,M3>
    {
        template <ptrdiff_t I, ptrdiff_t M, ptrdiff_t J, ptrdiff_t N>
        struct Unroller
        {
            static TMV_INLINE void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
            {
                Unroller<I,M/2,J,N>::unroll(x1,m1,x2,m2,m3);
                Unroller<I+M/2,M-M/2,J+M/2,N-M/2>::unroll(x1,m1,x2,m2,m3);
            }
        };
        template <ptrdiff_t I, ptrdiff_t J, ptrdiff_t N>
        struct Unroller<I,1,J,N>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
            {
                Unroller<I,1,J,N/2>::unroll(x1,m1,x2,m2,m3);
                Unroller<I,1,J+N/2,N-N/2>::unroll(x1,m1,x2,m2,m3);
            }
        };
        template <ptrdiff_t I, ptrdiff_t J, ptrdiff_t N>
        struct Unroller<I,0,J,N>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3) {}
        };
        template <ptrdiff_t I, ptrdiff_t J>
        struct Unroller<I,1,J,1>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
            { m3.ref(I,J) = x1 * m1.cref(I,J) + x2 * m2.cref(I,J); }
        };
        template <ptrdiff_t I, ptrdiff_t J>
        struct Unroller<I,1,J,0>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3) {}
        };
        static inline void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUU algo 15: N,s = "<<m3.size()<<','<<s<<std::endl;
#endif
            Unroller<0,s,0,s>::unroll(x1,m1,x2,m2,m3); 
        }
    };

    // algo 16: Fully unroll by columns
    template <ptrdiff_t s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUU_Helper<16,s,ix1,T1,M1,ix2,T2,M2,M3>
    {
        template <ptrdiff_t I, ptrdiff_t M, ptrdiff_t J, ptrdiff_t N>
        struct Unroller
        {
            static TMV_INLINE void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
            {
                Unroller<I,M-(N-N/2),J,N/2>::unroll(x1,m1,x2,m2,m3);
                Unroller<I,M,J+N/2,N-N/2>::unroll(x1,m1,x2,m2,m3);
            }
        };
        template <ptrdiff_t I, ptrdiff_t M, ptrdiff_t J>
        struct Unroller<I,M,J,1>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
            {
                Unroller<I,M/2,J,1>::unroll(x1,m1,x2,m2,m3);
                Unroller<I+M/2,M-M/2,J,1>::unroll(x1,m1,x2,m2,m3);
            }
        };
        template <ptrdiff_t I, ptrdiff_t M, ptrdiff_t J>
        struct Unroller<I,M,J,0>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3) {}
        };
        template <ptrdiff_t I, ptrdiff_t J>
        struct Unroller<I,1,J,1>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
            { m3.ref(I,J) = x1 * m1.cref(I,J) + x2 * m2.cref(I,J); }
        };
        template <ptrdiff_t I, ptrdiff_t J>
        struct Unroller<I,0,J,1>
        {
            static TMV_INLINE void unroll(
                const Scaling<ix1,T1>& x1, const M1& m1, 
                const Scaling<ix2,T2>& x2, const M2& m2, M3& m3) {}
        };
        static inline void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUU algo 16: N,s = "<<m3.size()<<','<<s<<std::endl;
#endif
            Unroller<0,s,0,s>::unroll(x1,m1,x2,m2,m3); 
        }
    };

    // algo 90: Call inst
    template <ptrdiff_t s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUU_Helper<90,s,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type T3;
            T3 xx1 = Traits<T3>::convert(T1(x1));
            T3 xx2 = Traits<T3>::convert(T2(x2));
            InstAddMM(xx1,m1.xView(),xx2,m2.xView(),m3.xView());
        }
    };

    // algo 91: Call inst alias
    template <ptrdiff_t s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUU_Helper<91,s,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            typedef typename M3::value_type T3;
            T3 xx1 = Traits<T3>::convert(T1(x1));
            T3 xx2 = Traits<T3>::convert(T2(x2));
            InstAliasAddMM(xx1,m1.xView(),xx2,m2.xView(),m3.xView());
        }
    };

    // algo 96: LowerTri, transpose:
    template <ptrdiff_t s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUU_Helper<96,s,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUU algo 201: N,s = "<<m3.size()<<','<<s<<std::endl;
#endif
            typedef typename M1::const_transpose_type M1t;
            typedef typename M2::const_transpose_type M2t;
            typedef typename M3::transpose_type M3t;
            M1t m1t = m1.transpose();
            M2t m2t = m2.transpose();
            M3t m3t = m3.transpose();
            AddUU_Helper<-2,s,ix1,T1,M1t,ix2,T2,M2t,M3t>::call(
                x1,m1t,x2,m2t,m3t);
        }
    };

    // algo 196: LowerTri, transpose:
    template <ptrdiff_t s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUU_Helper<196,s,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUU algo 196: N,s = "<<m3.size()<<','<<s<<std::endl;
#endif
            typedef typename M1::const_transpose_type M1t;
            typedef typename M2::const_transpose_type M2t;
            typedef typename M3::transpose_type M3t;
            M1t m1t = m1.transpose();
            M2t m2t = m2.transpose();
            M3t m3t = m3.transpose();
            AddUU_Helper<99,s,ix1,T1,M1t,ix2,T2,M2t,M3t>::call(
                x1,m1t,x2,m2t,m3t);
        }
    };

    // algo 97: Conjugate
    template <ptrdiff_t s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUU_Helper<97,s,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            typedef typename M3::conjugate_type M3c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            M3c m3c = m3.conjugate();
            AddUU_Helper<-2,s,ix1,T1,M1c,ix2,T2,M2c,M3c>::call(
                TMV_CONJ(x1),m1c,TMV_CONJ(x2),m2c,m3c);
        }
    };

    // algo 197: Conjugate
    template <ptrdiff_t s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUU_Helper<197,s,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            typedef typename M1::const_conjugate_type M1c;
            typedef typename M2::const_conjugate_type M2c;
            typedef typename M3::conjugate_type M3c;
            M1c m1c = m1.conjugate();
            M2c m2c = m2.conjugate();
            M3c m3c = m3.conjugate();
            AddUU_Helper<99,s,ix1,T1,M1c,ix2,T2,M2c,M3c>::call(
                TMV_CONJ(x1),m1c,TMV_CONJ(x2),m2c,m3c);
        }
    };

    // algo 98: Inlinst check for aliases 
    template <ptrdiff_t s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUU_Helper<98,s,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            const bool s1 = SameStorage(m1,m3);
            const bool s2 = SameStorage(m2,m3);
            const bool ss1 = s1 && !(M1::_unit && OppositeStorage(m1,m3));
            const bool ss2 = s2 && !(M2::_unit && OppositeStorage(m2,m3));
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUU Check aliases:\n";
            std::cout<<"s1,s2 = "<<s1<<','<<s2<<std::endl;
            std::cout<<"ss1,ss2 = "<<ss1<<','<<ss2<<std::endl;
#endif

            if (!ss1 && !ss2) {
                // No aliasing (or no clobbering)
                AddUU_Helper<-2,s,ix1,T1,M1,ix2,T2,M2,M3>::call(
                    x1,m1,x2,m2,m3);
            } else if (!ss2) {
                // Alias with m1 only, do m1 first
                MultXM<false>(x1,m1,m3);
                typename M3::noalias_type m3na = m3.noAlias();
                MultXM<true>(x2,m2,m3na);
            } else if (!ss1) {
                // Alias with m2 only, do m2 first
                MultXM<false>(x2,m2,m3);
                typename M3::noalias_type m3na = m3.noAlias();
                MultXM<true>(x1,m1,m3na);
            } else {
                // Need a temporary
                typename M1::copy_type m1c = m1;
                MultXM<false>(x2,m2,m3);
                typename M3::noalias_type m3na = m3.noAlias();
                MultXM<true>(x1,m1c,m3na);
            }
        }
    };

    // algo 99: Check for aliases
    template <ptrdiff_t s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUU_Helper<99,s,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            typedef typename M1::value_type TM1;
            typedef typename M2::value_type TM2;
            typedef typename M3::value_type TM3;
            const bool inst =
                (s == Unknown || s > 16) &&
#ifdef TMV_INST_MIX
                Traits2<TM1,TM2>::samebase &&
                Traits2<TM1,TM3>::samebase &&
#else
                Traits2<TM1,TM2>::sametype &&
                Traits2<TM1,TM3>::sametype &&
#endif
                Traits<TM3>::isinst;
            const bool conj = M3::_conj;
            const int algo = 
                M3::_lower ? 196 :
                conj ? 197 :
                inst ? 91 :
                98;
            AddUU_Helper<algo,s,ix1,T1,M1,ix2,T2,M2,M3>::call(x1,m1,x2,m2,m3);
        }
    };
    // algo -4: No branches or copies
    template <ptrdiff_t s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUU_Helper<-4,s,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            TMVStaticAssert(!M3::_conj);
            TMVStaticAssert(!M1::_unit);
            TMVStaticAssert(!M2::_unit);
            TMVStaticAssert(!M3::_unit);
            TMVStaticAssert(M1::_upper);
            TMVStaticAssert(M2::_upper);
            TMVStaticAssert(M3::_upper);
            TMVAssert(!m1.isunit());
            TMVAssert(!m2.isunit());
            TMVAssert(!m3.isunit());
            const ptrdiff_t s2 = s > 20 ? Unknown : s;
            const ptrdiff_t s2p1 = IntTraits<s2>::Sp1;
            // nops = n(n+1)
            const ptrdiff_t nops = IntTraits2<s2,s2p1>::safeprod;
            const bool unroll = 
                s > 10 ? false :
                s == Unknown ? false :
                nops <= TMV_ADDUU_UNROLL;
            const int algo = 
                unroll ? (
                    ( (M1::_rowmajor && M2::_rowmajor) ||
                      (M1::_rowmajor && M3::_rowmajor) ||
                      (M2::_rowmajor && M3::_rowmajor) ) ? 15 : 16 ) :
                ( (M1::_rowmajor && M2::_rowmajor) ||
                  (M1::_rowmajor && M3::_rowmajor) ||
                  (M2::_rowmajor && M3::_rowmajor) ) ? 12 : 11;
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUU algo -4: N,s = "<<m3.size()<<','<<s<<std::endl;
            std::cout<<"x1 = "<<ix1<<"  "<<T1(x1)<<std::endl;
            std::cout<<"x2 = "<<ix2<<"  "<<T2(x2)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            AddUU_Helper<algo,s,ix1,T1,M1,ix2,T2,M2,M3>::call(
                x1,m1,x2,m2,m3);
        }
    };

    // algo -3: Determine which algorithm to use
    template <ptrdiff_t s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUU_Helper<-3,s,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            TMVStaticAssert(!M3::_conj);
            TMVStaticAssert(M3::_upper);
            const ptrdiff_t s2 = s > 20 ? Unknown : s;
            const ptrdiff_t s2p1 = IntTraits<s2>::Sp1;
            // nops = n(n+1)
            const ptrdiff_t nops = IntTraits2<s2,s2p1>::safeprod;
            const bool unroll = 
                s > 10 ? false :
                s == Unknown ? false :
                nops <= TMV_ADDUU_UNROLL;
            const int algo = 
                ( M1::_unit || M2::_unit ) ? 2 :
                unroll ? (
                    ( (M1::_rowmajor && M2::_rowmajor) ||
                      (M1::_rowmajor && M3::_rowmajor) ||
                      (M2::_rowmajor && M3::_rowmajor) ) ? 15 : 16 ) :
                (M1::_unknowndiag || M2::_unknowndiag) ? 3 :
                ( (M1::_rowmajor && M2::_rowmajor) ||
                  (M1::_rowmajor && M3::_rowmajor) ||
                  (M2::_rowmajor && M3::_rowmajor) ) ? 12 : 11;
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUU algo -3: N,s = "<<m3.size()<<','<<s<<std::endl;
            std::cout<<"x1 = "<<ix1<<"  "<<T1(x1)<<std::endl;
            std::cout<<"x2 = "<<ix2<<"  "<<T2(x2)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            AddUU_Helper<algo,s,ix1,T1,M1,ix2,T2,M2,M3>::call(
                x1,m1,x2,m2,m3);
        }
    };

    // algo -2: Check for inst
    template <ptrdiff_t s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUU_Helper<-2,s,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUU algo -2: N,s = "<<m3.size()<<','<<s<<std::endl;
#endif
            typedef typename M1::value_type TM1;
            typedef typename M2::value_type TM2;
            typedef typename M3::value_type TM3;
            const bool inst =
                (s == Unknown || s > 16) &&
#ifdef TMV_INST_MIX
                Traits2<TM1,TM2>::samebase &&
                Traits2<TM1,TM3>::samebase &&
#else
                Traits2<TM1,TM2>::sametype &&
                Traits2<TM1,TM3>::sametype &&
#endif
                Traits<TM3>::isinst;
            const bool conj = M3::_conj;
            const int algo = 
                M3::_lower ? 96 :
                conj ? 97 :
                inst ? 90 :
                -3;
            AddUU_Helper<algo,s,ix1,T1,M1,ix2,T2,M2,M3>::call(x1,m1,x2,m2,m3);
        }
    };

    // algo -1: Check for aliases?
    template <ptrdiff_t s, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUU_Helper<-1,s,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUU algo -1: N,s = "<<m3.size()<<','<<s<<std::endl;
#endif
            const int algo = 
                M3::_checkalias ? 99 : 
                -2;
            AddUU_Helper<algo,s,ix1,T1,M1,ix2,T2,M2,M3>::call(x1,m1,x2,m2,m3);
        }
    };

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    inline void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Tri<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Tri_Mutable<M3>& m3)
    {
        TMVStaticAssert(M1::_upper == int(M3::_upper));
        TMVStaticAssert(M2::_upper == int(M3::_upper));
        TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
        TMVStaticAssert((Sizes<M1::_size,M3::_size>::same));
        TMVStaticAssert(!M3::_unit);
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.size() == m3.size());
        TMVAssert(!m3.isunit());
        const ptrdiff_t s = Sizes<Sizes<M1::_size,M2::_size>::size,M3::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        AddUU_Helper<-1,s,ix1,T1,M1v,ix2,T2,M2v,M3v>::call(x1,m1v,x2,m2v,m3v);
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    inline void InlineAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Tri<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Tri_Mutable<M3>& m3)
    {
        TMVStaticAssert(M1::_upper == int(M3::_upper));
        TMVStaticAssert(M2::_upper == int(M3::_upper));
        TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
        TMVStaticAssert((Sizes<M1::_size,M3::_size>::same));
        TMVStaticAssert(!M3::_unit);
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.size() == m3.size());
        TMVAssert(!m3.isunit());
        const ptrdiff_t s = Sizes<Sizes<M1::_size,M2::_size>::size,M3::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        AddUU_Helper<-3,s,ix1,T1,M1v,ix2,T2,M2v,M3v>::call(x1,m1v,x2,m2v,m3v);
    }

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    inline void InlineAliasAddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Tri<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Tri_Mutable<M3>& m3)
    {
        TMVStaticAssert(M1::_upper == int(M3::_upper));
        TMVStaticAssert(M2::_upper == int(M3::_upper));
        TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
        TMVStaticAssert((Sizes<M1::_size,M3::_size>::same));
        TMVStaticAssert(!M3::_unit);
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.size() == m3.size());
        TMVAssert(!m3.isunit());
        const ptrdiff_t s = Sizes<Sizes<M1::_size,M2::_size>::size,M3::_size>::size;
        typedef typename M1::const_cview_type M1v;
        typedef typename M2::const_cview_type M2v;
        typedef typename M3::cview_type M3v;
        TMV_MAYBE_CREF(M1,M1v) m1v = m1.cView();
        TMV_MAYBE_CREF(M2,M2v) m2v = m2.cView();
        TMV_MAYBE_REF(M3,M3v) m3v = m3.cView();
        AddUU_Helper<98,s,ix1,T1,M1v,ix2,T2,M2v,M3v>::call(x1,m1v,x2,m2v,m3v);
    }


    //
    // M = x * U + x * U
    //

    template <int algo, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUUM_Helper;

    // algo 11: M = x * U + x * U  (with alias checking)
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUUM_Helper<11,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUUM algo 11\n";
#endif
            typename M3::uppertri_type m3u = m3.upperTri();
            AddMM(x1,m1,x2,m2,m3u);
            if (m1.size() > 1) m3.lowerTri().offDiag().setZero();
        }
    };

    // algo 12: M = x * U + x * L  (with alias checking)
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUUM_Helper<12,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            const bool s1 = SameStorage(m1,m3);
            const bool s2 = SameStorage(m2,m3);
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUUM algo 12\n";
            std::cout<<"s1,s22 = "<<s1<<','<<s2<<std::endl;
#endif

            if (!s1 && !s2) {
                // No aliasing
                AddUUM_Helper<22,ix1,T1,M1,ix2,T2,M2,M3>::call(x1,m1,x2,m2,m3);
            } else if (!s2) {
                // Alias with m1 only, do m1 first
                typename M3::uppertri_type m3u = m3.upperTri();
                typename M3::lowertri_type::noalias_type m3l =
                    m3.lowerTri().noAlias();
                MultXM<false>(x1,m1,m3u);
                MultXM<true>(x2,m2,m3l);
            } else if (!s1) {
                // Alias with m2 only, do m2 first
                typename M3::uppertri_type::noalias_type m3u = 
                    m3.upperTri().noAlias();
                typename M3::lowertri_type m3l = m3.lowerTri();
                MultXM<false>(x2,m2,m3l);
                MultXM<true>(x1,m1,m3u);
            } else {
                // Need a temporary
                typename M3::uppertri_type::noalias_type m3u =
                    m3.upperTri().noAlias();
                typename M3::lowertri_type m3l = m3.lowerTri();
                typename M1::copy_type m1c = m1;
                MultXM<false>(x2,m2,m3l);
                MultXM<true>(x1,m1c,m3u);
            }
        }
    };

    // algo 13: M = x * L + x * U  (with alias checking)
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUUM_Helper<13,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            const bool s1 = SameStorage(m1,m3);
            const bool s2 = SameStorage(m2,m3);
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUUM algo 13\n";
            std::cout<<"s1,s22 = "<<s1<<','<<s2<<std::endl;
#endif

            if (!s1 && !s2) {
                // No aliasing
                AddUUM_Helper<23,ix1,T1,M1,ix2,T2,M2,M3>::call(x1,m1,x2,m2,m3);
            } else if (!s2) {
                // Alias with m1 only, do m1 first
                typename M3::uppertri_type::noalias_type m3u =
                    m3.upperTri().noAlias();
                typename M3::lowertri_type m3l = m3.lowerTri();
                MultXM<false>(x1,m1,m3l);
                MultXM<true>(x2,m2,m3u);
            } else if (!s1) {
                // Alias with m2 only, do m2 first
                typename M3::uppertri_type m3u = m3.upperTri();
                typename M3::lowertri_type::noalias_type m3l =
                    m3.lowerTri().noAlias();
                MultXM<false>(x2,m2,m3u);
                MultXM<true>(x1,m1,m3l);
            } else {
                // Need a temporary
                typename M3::uppertri_type m3u = m3.upperTri();
                typename M3::lowertri_type::noalias_type m3l =
                    m3.lowerTri().noAlias();
                typename M1::copy_type m1c = m1;
                MultXM<false>(x2,m2,m3u);
                MultXM<true>(x1,m1c,m3l);
            }
        }
    };

    // algo 14: M = x * L + x * L  (with alias checking)
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUUM_Helper<14,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUUM algo 14\n";
#endif
            typename M3::lowertri_type m3l = m3.lowerTri();
            AddMM(x1,m1,x2,m2,m3l);
            if (m1.size() > 1) m3.upperTri().offDiag().setZero();
        }
    };

    // algo 21: M = x * U + x * U  (no alias checking)
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUUM_Helper<21,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUUM algo 21\n";
#endif
            typename M3::uppertri_type::noalias_type m3u =
                m3.upperTri().noAlias();
            m3.setZero();
            AddMM(x1,m1,x2,m2,m3u);
            if (m1.size() > 1) m3.lowerTri().offDiag().setZero();
        }
    };

    // algo 22: M = x * U + x * L  (no alias checking)
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUUM_Helper<22,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUUM algo 22\n";
#endif
            if (m1.size() > 0) {
                typename M3::uppertri_type::offdiag_type::noalias_type m3u = 
                    m3.upperTri().offDiag().noAlias();
                typename M3::lowertri_type::offdiag_type::noalias_type m3l = 
                    m3.lowerTri().offDiag().noAlias();
                MultXM<false>(x1,m1.offDiag(),m3u);
                MultXM<false>(x2,m2.offDiag(),m3l);
            }
            AddUUM_Helper<31,ix1,T1,M1,ix2,T2,M2,M3>::call(x1,m1,x2,m2,m3);
        }
    };

    // algo 23: M = x * L + x * U  (no alias checking)
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUUM_Helper<23,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUUM algo 23\n";
#endif
            if (m1.size() > 0) {
                typename M3::uppertri_type::offdiag_type::noalias_type m3u = 
                    m3.upperTri().offDiag().noAlias();
                typename M3::lowertri_type::offdiag_type::noalias_type m3l = 
                    m3.lowerTri().offDiag().noAlias();
                MultXM<false>(x1,m1.offDiag(),m3l);
                MultXM<false>(x2,m2.offDiag(),m3u);
            }
            AddUUM_Helper<31,ix1,T1,M1,ix2,T2,M2,M3>::call(x1,m1,x2,m2,m3);
        }
    };

    // algo 24: M = x * L + x * L  (no alias checking)
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUUM_Helper<24,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUUM algo 24\n";
#endif
            typename M3::lowertri_type::noalias_type m3l =
                m3.lowerTri().noAlias();
            m3.setZero();
            AddMM(x1,m1,x2,m2,m3l);
        }
    };

    // 31-34 just handle the diagonals after the offdiagonals have been 
    // done already.
    
    // algo 31: check for unit and unknown diags
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUUM_Helper<31,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUUM algo 31\n";
#endif
            const int algo = 
                (M1::_unit || M2::_unit) ? 33 :
                (M1::_unknowndiag || M2::_unknowndiag) ? 32 :
                34;
            AddUUM_Helper<algo,ix1,T1,M1,ix2,T2,M2,M3>::call(x1,m1,x2,m2,m3);
        }
    };

    // algo 32: UnknownDiag
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUUM_Helper<32,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUUM algo 32\n";
#endif
            if (m1.isunit() || m2.isunit())
                AddUUM_Helper<33,ix1,T1,M1,ix2,T2,M2,M3>::call(
                    x1,m1,x2,m2,m3);
            else
                AddUUM_Helper<34,ix1,T1,M1,ix2,T2,M2,M3>::call(
                    x1,m1,x2,m2,m3);
        }
    };

    // algo 33: m1 and/or m2 is unitdiag
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUUM_Helper<33,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUUM algo 33\n";
#endif
            typename M3::diag_type m3d = m3.diag();
            if (m1.isunit()) {
                if (m2.isunit()) {
                    m3d.setAllTo(T1(x1)+T2(x2));
                } else {
                    MultXV<false>(x2,m2.diag(),m3d);
                    m3d.addToAll(T1(x1));
                }
            } else {
                MultXV<false>(x1,m1.diag(),m3d);
                m3d.addToAll(T2(x2));
            }
        }
    };

    // algo 34: no unitdiag's
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUUM_Helper<34,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUUM algo 34\n";
#endif
            typename M3::diag_type m3d = m3.diag();
            AddVV(x1,m1.diag(),x2,m2.diag(),m3d);
        }
    };

    // algo 99: Check for aliases
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUUM_Helper<99,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            const int algo = 
                (M1::_upper && M2::_upper) ? 11 :
                (M1::_upper && M2::_lower) ? 12 :
                (M1::_lower && M2::_upper) ? 13 :
                14;
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUUM algo 99\n";
            std::cout<<"x1 = "<<ix1<<"  "<<T1(x1)<<std::endl;
            std::cout<<"x2 = "<<ix2<<"  "<<T2(x2)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            AddUUM_Helper<algo,ix1,T1,M1,ix2,T2,M2,M3>::call(x1,m1,x2,m2,m3);
        }
    };
    // algo -2: No alias checking
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUUM_Helper<-2,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            const int algo = 
                (M1::_upper && M2::_upper) ? 21 :
                (M1::_upper && M2::_lower) ? 22 :
                (M1::_lower && M2::_upper) ? 23 :
                24;
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUUM algo -2\n";
            std::cout<<"x1 = "<<ix1<<"  "<<T1(x1)<<std::endl;
            std::cout<<"x2 = "<<ix2<<"  "<<T2(x2)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            AddUUM_Helper<algo,ix1,T1,M1,ix2,T2,M2,M3>::call(x1,m1,x2,m2,m3);
        }
    };

    // algo -1: Check for aliases?
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUUM_Helper<-1,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUUM algo -1\n";
#endif
            const int algo = 
                M3::_checkalias ? 99 :
                -2;
            AddUUM_Helper<algo,ix1,T1,M1,ix2,T2,M2,M3>::call(x1,m1,x2,m2,m3);
        }
    };

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    inline void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Tri<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_size,M2::_size>::same));
        TMVStaticAssert((Sizes<M1::_size,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_size,M3::_rowsize>::same));
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.size() == m3.colsize());
        TMVAssert(m1.size() == m3.rowsize());
        AddUUM_Helper<-1,ix1,T1,M1,ix2,T2,M2,M3>::call(
            x1,m1.mat(),x2,m2.mat(),m3.mat());
    }


    //
    // M = x * U + x * M
    //

    template <int algo, int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUMM_Helper;

    // algo 11: M = x * U + x * M  (with alias checking)
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUMM_Helper<11,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            TMVStaticAssert(M1::_upper);
            const bool s1 = SameStorage(m1,m3);
            const bool s2 = SameStorage(m2,m3);
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUMM algo 11\n";
            std::cout<<"s1,s2 = "<<s1<<','<<s2<<std::endl;
#endif

            if (!s1 && !s2) {
                // No aliasing
                AddUMM_Helper<21,ix1,T1,M1,ix2,T2,M2,M3>::call(x1,m1,x2,m2,m3);
            } else if (!s2) {
                // Alias with m1 only, do m1 first
                typename M3::uppertri_type m3u = m3.upperTri();
                MultXM<false>(x1,m1,m3u);
                typename M3::noalias_type m3na = m3.noAlias();
                MultXM<true>(x2,m2,m3na);
            } else if (!s1) {
                // Alias with m2 only, do m2 first
                typename M3::uppertri_type::noalias_type m3u =
                    m3.upperTri().noAlias();
                MultXM<false>(x2,m2,m3);
                MultXM<true>(x1,m1,m3u);
            } else {
                // Need a temporary
                typename M3::uppertri_type::noalias_type m3u =
                    m3.upperTri().noAlias();
                typename M1::copy_type m1c = m1;
                MultXM<false>(x2,m2,m3);
                MultXM<true>(x1,m1c,m3u);
            }
        }
    };

    // algo 12: M = x * L + x * M  (with alias checking)
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUMM_Helper<12,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            TMVStaticAssert(M1::_lower);
            const bool s1 = SameStorage(m1,m3);
            const bool s2 = SameStorage(m2,m3);
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUMM algo 12\n";
            std::cout<<"s1,s2 = "<<s1<<','<<s2<<std::endl;
#endif

            if (!s1 && !s2) {
                // No aliasing
                AddUMM_Helper<22,ix1,T1,M1,ix2,T2,M2,M3>::call(x1,m1,x2,m2,m3);
            } else if (!s2) {
                // Alias with m1 only, do m1 first
                typename M3::lowertri_type m3l = m3.lowerTri();
                MultXM<false>(x1,m1,m3l);
                typename M3::noalias_type m3na = m3.noAlias();
                MultXM<true>(x2,m2,m3na);
            } else if (!s1) {
                // Alias with m2 only, do m2 first
                typename M3::lowertri_type::noalias_type m3l =
                    m3.lowerTri().noAlias();
                MultXM<false>(x2,m2,m3);
                MultXM<true>(x1,m1,m3l);
            } else {
                // Need a temporary
                typename M3::lowertri_type::noalias_type m3l =
                    m3.lowerTri().noAlias();
                typename M1::copy_type m1c = m1;
                MultXM<false>(x2,m2,m3);
                MultXM<true>(x1,m1c,m3l);
            }
        }
    };

    // algo 21: M = x * U + x * M  (no alias checking)
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUMM_Helper<21,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUMM algo 21\n";
#endif
            TMVStaticAssert(M1::_upper);
            typename M3::uppertri_type::noalias_type m3u =
                m3.upperTri().noAlias();
            typename M3::noalias_type m3na = m3.noAlias();
            MultXM<false>(x2,m2,m3na);
            MultXM<true>(x1,m1,m3u);
        }
    };

    // algo 22: M = x * L + x * M  (no alias checking)
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUMM_Helper<22,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUMM algo 22\n";
#endif
            TMVStaticAssert(M1::_lower);
            typename M3::lowertri_type::noalias_type m3l =
                m3.lowerTri().noAlias();
            typename M3::noalias_type m3na = m3.noAlias();
            MultXM<false>(x2,m2,m3na);
            MultXM<true>(x1,m1,m3l);
        }
    };

    // 31-34 just handle the diagonals after the offdiagonals have been 
    // done already.
    
    // algo 31: check for unit and unknown diags
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUMM_Helper<31,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUMM algo 31\n";
#endif
            const int algo = 
                M1::_unit ? 33 :
                M1::_unknowndiag ? 32 :
                34;
            AddUMM_Helper<algo,ix1,T1,M1,ix2,T2,M2,M3>::call(x1,m1,x2,m2,m3);
        }
    };

    // algo 32: UnknownDiag
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUMM_Helper<32,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUMM algo 32\n";
#endif
            if (m1.isunit())
                AddUMM_Helper<33,ix1,T1,M1,ix2,T2,M2,M3>::call(
                    x1,m1,x2,m2,m3);
            else
                AddUMM_Helper<34,ix1,T1,M1,ix2,T2,M2,M3>::call(
                    x1,m1,x2,m2,m3);
        }
    };

    // algo 33: m1 is UnitDiag
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUMM_Helper<33,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUMM algo 33\n";
#endif
            typename M3::diag_type m3d = m3.diag();
            MultXV<false>(x2,m2.diag(),m3d);
            m3d.addToAll(T1(x1));
        }
    };

    // algo 34: m1 is NonUnitDiag
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUMM_Helper<34,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUMM algo 34\n";
#endif
            typename M3::diag_type m3d = m3.diag();
            AddVV(x1,m1.diag(),x2,m2.diag(),m3d);
        }
    };

    // algo -2: No alias checking
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUMM_Helper<-2,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            const int algo = M1::_upper ? 21 : 22;
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUMM algo -2\n";
            std::cout<<"x1 = "<<ix1<<"  "<<T1(x1)<<std::endl;
            std::cout<<"x2 = "<<ix2<<"  "<<T2(x2)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            AddUMM_Helper<algo,ix1,T1,M1,ix2,T2,M2,M3>::call(x1,m1,x2,m2,m3);
        }
    };

    // algo 99: Check for aliases
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUMM_Helper<99,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
            const int algo = M1::_upper ? 11 : 12;
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUMM algo 99\n";
            std::cout<<"x1 = "<<ix1<<"  "<<T1(x1)<<std::endl;
            std::cout<<"x2 = "<<ix2<<"  "<<T2(x2)<<std::endl;
            std::cout<<"m1 = "<<TMV_Text(m1)<<std::endl;
            std::cout<<"m2 = "<<TMV_Text(m2)<<std::endl;
            std::cout<<"m3 = "<<TMV_Text(m3)<<std::endl;
            std::cout<<"algo = "<<algo<<std::endl;
#endif
            AddUMM_Helper<algo,ix1,T1,M1,ix2,T2,M2,M3>::call(x1,m1,x2,m2,m3);
        }
    };

    // algo -1: Check for aliases?
    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    struct AddUMM_Helper<-1,ix1,T1,M1,ix2,T2,M2,M3>
    {
        static TMV_INLINE void call(
            const Scaling<ix1,T1>& x1, const M1& m1, 
            const Scaling<ix2,T2>& x2, const M2& m2, M3& m3)
        {
#ifdef PRINTALGO_AddUU
            std::cout<<"AddUMM algo -1\n";
#endif
            const int algo = 
                M3::_checkalias ? 99 :
                -2;
            AddUMM_Helper<algo,ix1,T1,M1,ix2,T2,M2,M3>::call(x1,m1,x2,m2,m3);
        }
    };

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    inline void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Tri<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Rec<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    {
        TMVStaticAssert((Sizes<M1::_size,M2::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_size,M2::_rowsize>::same));
        TMVStaticAssert((Sizes<M1::_size,M3::_colsize>::same));
        TMVStaticAssert((Sizes<M1::_size,M3::_rowsize>::same));
        TMVAssert(m1.size() == m2.colsize());
        TMVAssert(m1.size() == m2.rowsize());
        TMVAssert(m1.size() == m3.colsize());
        TMVAssert(m1.size() == m3.rowsize());
        AddUMM_Helper<-1,ix1,T1,M1,ix2,T2,M2,M3>::call(
            x1,m1.mat(),x2,m2.mat(),m3.mat());
    }


    //
    // M = x * M + x * U
    //

    template <int ix1, class T1, class M1, int ix2, class T2, class M2, class M3>
    TMV_INLINE void AddMM(
        const Scaling<ix1,T1>& x1, const BaseMatrix_Rec<M1>& m1, 
        const Scaling<ix2,T2>& x2, const BaseMatrix_Tri<M2>& m2, 
        BaseMatrix_Rec_Mutable<M3>& m3)
    { AddMM(x2,m2,x1,m1,m3); }


#undef TMV_ADDUU_UNROLL

} // namespace tmv

#endif 
