
//#define PRINTALGO_LU

#include "TMV_Blas.h"
#include "tmv/TMV_LUDiv.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_SmallMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_CopyV.h"
#include "tmv/TMV_DivVU.h"
#include "tmv/TMV_DivMU.h"
#include "tmv/TMV_PermuteM.h"
#include "tmv/TMV_ScaleM.h"

namespace tmv {

    template <bool trans>
    struct NonLapLUSolve // trans == false here
    {
        template <class M1, class M2>
        static void call(const M1& m1, const Permutation& P, M2& m2)
        { InlineLU_SolveInPlace(m1,P,m2); }
    };
    template <>
    struct NonLapLUSolve<true>
    {
        template <class M1, class M2>
        static void call(const M1& m1, const Permutation& P, M2& m2)
        { InlineLU_SolveTransposeInPlace(m1,P,m2); }
    };

#ifdef ALAP
    // ALAP, not LAP, since ATLAS has these routines
    template <bool trans, class M1, class T2, class T1> 
    static inline void LapLUSolve(
        const M1& m1, const Permutation& P, MatrixView<T2> m2, T1)
    { NonLapLUSolve<trans>::call(m1,P,m2); }
    template <bool trans, class M1, class T2, class T1> 
    static inline void LapLUSolve(
        const M1& m1, const Permutation& P, VectorView<T2,Unit> v2, T1)
    { NonLapLUSolve<trans>::call(m1,P,v2); }
#ifdef TMV_INST_DOUBLE
    template <bool trans, class M1> 
    static void LapLUSolve(
        const M1& m1, const Permutation& P, MatrixView<double> m2, double)
    {
        TMVAssert(P.isInverse());
        int n = m1.colsize();
        int nrhs = m2.rowsize();
        int lda = m1.stepj();
        int ldb = m2.stepj();
#ifdef CLAP
        const int* ipiv = P.getValues();
#else
        AlignedArray<int> ipiv1(n);
        for(int i=0;i<n;++i) ipiv1[i] = P.getValues()[i] LAPPLUS1;
        const int* ipiv = ipiv1.get();
#endif
        LAPNAME(dgetrs) (
            LAPCM trans?LAPCH_T:LAPCH_NT,
            LAPV(n),LAPV(nrhs),LAPP(m1.cptr()),LAPV(lda),
            LAPP(ipiv),LAPP(m2.ptr()),LAPV(ldb) LAPINFO LAP1);
        LAP_Results("dgetrs");
    }
    template <bool trans, class M1> 
    static void LapLUSolve(
        const M1& m1, const Permutation& P,
        MatrixView<std::complex<double> > m2, std::complex<double> )
    {
        TMVAssert(P.isInverse());
        int n = m1.colsize();
        int nrhs = m2.rowsize();
        int lda = m1.stepj();
        int ldb = m2.stepj();
#ifdef CLAP
        const int* ipiv = P.getValues();
#else
        AlignedArray<int> ipiv1(n);
        for(int i=0;i<n;++i) ipiv1[i] = P.getValues()[i] LAPPLUS1;
        const int* ipiv = ipiv1.get();
#endif
        if (!trans && m1.isconj()) m2.conjugateSelf();
        LAPNAME(zgetrs) (
            LAPCM trans?(m1.isconj()?LAPCH_CT:LAPCH_T):LAPCH_NT,
            LAPV(n),LAPV(nrhs),LAPP(m1.cptr()),LAPV(lda),
            LAPP(ipiv),LAPP(m2.ptr()),LAPV(ldb) LAPINFO LAP1);
        if (!trans && m1.isconj()) m2.conjugateSelf();
        LAP_Results("zgetrs");
    }
    template <bool trans, class M1> 
    static void LapLUSolve(
        const M1& m1, const Permutation& P,
        MatrixView<std::complex<double> > m2, double t1)
    {
        Matrix<double,ColMajor> temp = m2.realPart();
        LapLUSolve<trans>(m1,P,temp.xView(),t1);
        m2.realPart() = temp;
        temp = m2.imagPart();
        LapLUSolve<trans>(m1,P,temp.xView(),t1);
        m2.imagPart() = temp;
    }
    template <bool trans, class M1> 
    static inline void LapLUSolve(
        const M1& m1, const Permutation& P, VectorView<double,Unit> v2, double t1)
    { LapLUSolve<trans>(m1,P,ColVectorViewOf(v2).xView(),t1); }
    template <bool trans, class M1> 
    static inline void LapLUSolve(
        const M1& m1, const Permutation& P, 
        VectorView<std::complex<double>,Unit> v2, std::complex<double> t1)
    { LapLUSolve<trans>(m1,P,ColVectorViewOf(v2).xView(),t1); }
    template <bool trans, class M1> 
    static inline void LapLUSolve(
        const M1& m1, const Permutation& P, 
        VectorView<std::complex<double>,Unit> v2, double t1)
    { LapLUSolve<trans>(m1,P,ColVectorViewOf(v2).xView(),t1); }
#endif
#ifdef TMV_INST_FLOAT
#ifndef MKL // TODO: Check if this is still required...
    template <bool trans, class M1> 
    static void LapLUSolve(
        const M1& m1, const Permutation& P, MatrixView<float> m2, float)
    {
        TMVAssert(P.isInverse());
        int n = m1.colsize();
        int nrhs = m2.rowsize();
        int lda = m1.stepj();
        int ldb = m2.stepj();
#ifdef CLAP
        const int* ipiv = P.getValues();
#else
        AlignedArray<int> ipiv1(n);
        for(int i=0;i<n;++i) ipiv1[i] = P.getValues()[i] LAPPLUS1;
        const int* ipiv = ipiv1.get();
#endif
        LAPNAME(sgetrs) (
            LAPCM trans?LAPCH_T:LAPCH_NT,
            LAPV(n),LAPV(nrhs),LAPP(m1.cptr()),LAPV(lda),
            LAPP(ipiv),LAPP(m2.ptr()),LAPV(ldb) LAPINFO LAP1);
        LAP_Results("sgetrs");
    }
    template <bool trans, class M1> 
    static void LapLUSolve(
        const M1& m1, const Permutation& P,
        MatrixView<std::complex<float> > m2, std::complex<float> )
    {
        TMVAssert(P.isInverse());
        int n = m1.colsize();
        int nrhs = m2.rowsize();
        int lda = m1.stepj();
        int ldb = m2.stepj();
#ifdef CLAP
        const int* ipiv = P.getValues();
#else
        AlignedArray<int> ipiv1(n);
        for(int i=0;i<n;++i) ipiv1[i] = P.getValues()[i] LAPPLUS1;
        const int* ipiv = ipiv1.get();
#endif
        if (!trans && m1.isconj()) m2.conjugateSelf();
        LAPNAME(cgetrs) (
            LAPCM trans?(m1.isconj()?LAPCH_CT:LAPCH_T):LAPCH_NT,
            LAPV(n),LAPV(nrhs),LAPP(m1.cptr()),LAPV(lda),
            LAPP(ipiv),LAPP(m2.ptr()),LAPV(ldb) LAPINFO LAP1);
        if (!trans && m1.isconj()) m2.conjugateSelf();
        LAP_Results("cgetrs");
    }
    template <bool trans, class M1> 
    static void LapLUSolve(
        const M1& m1, const Permutation& P,
        MatrixView<std::complex<float> > m2, float t1)
    {
        Matrix<float,ColMajor> temp = m2.realPart();
        LapLUSolve<trans>(m1,P,temp.xView(),t1);
        m2.realPart() = temp;
        temp = m2.imagPart();
        LapLUSolve<trans>(m1,P,temp.xView(),t1);
        m2.imagPart() = temp;
    }
    template <bool trans, class M1> 
    static inline void LapLUSolve(
        const M1& m1, const Permutation& P, VectorView<float,Unit> v2, float t1)
    { LapLUSolve<trans>(m1,P,ColVectorViewOf(v2).xView(),t1); }
    template <bool trans, class M1> 
    static inline void LapLUSolve(
        const M1& m1, const Permutation& P, 
        VectorView<std::complex<float>,Unit> v2, std::complex<float> t1)
    { LapLUSolve<trans>(m1,P,ColVectorViewOf(v2).xView(),t1); }
    template <bool trans, class M1> 
    static inline void LapLUSolve(
        const M1& m1, const Permutation& P,
        VectorView<std::complex<float>,Unit> v2, float t1)
    { LapLUSolve<trans>(m1,P,ColVectorViewOf(v2).xView(),t1); }
#endif // MKL
#endif // FLOAT
#endif // ALAP

    template <bool trans, class M1, class T2>
    static void DoLUSolve(
        const M1& m1, const Permutation& P, MatrixView<T2> m2)
    {
#ifdef ALAP
        if (m1.stepj()>0) {
            const typename M1::value_type t1(0);
            if ( (m2.iscm() && m2.stepj()>0) ) {
                LapLUSolve<trans>(m1,P,m2,t1);
            } else {
                Matrix<T2,ColMajor|NoDivider> m2c(m2);
                LapLUSolve<trans>(m1,P,m2c.xView(),t1);
                InstCopy(m2c.constView().xView(),m2);
            }
        } else {
            Matrix<typename M1::value_type,ColMajor|NoDivider> m1c(m1);
            DoLUSolve<trans>(m1c.constView().xView(),P,m2);
        }
#else
        NonLapLUSolve<trans>::call(m1,P,m2);
#endif
    }

    template <bool trans, class M1, class T2>
    static void DoLUSolve(
        const M1& m1, const Permutation& P, VectorView<T2> v2)
    {
#ifdef ALAP
        if (m1.stepj()>0) {
            const typename M1::value_type t1(0);
            if ( v2.step() == 1 ) {
                LapLUSolve<trans>(m1,P,v2.unitView(),t1);
            } else  {
                Vector<T2> v2c(v2);
                LapLUSolve<trans>(m1,P,v2c.xView().unitView(),t1);
                InstCopy(v2c.constView().xView(),v2);
            }
        } else {
            Matrix<typename M1::value_type,ColMajor|NoDivider> m1c(m1);
            DoLUSolve<trans>(m1c.constView().xView(),P,v2);
        }
#else
        NonLapLUSolve<trans>::call(m1,P,v2);
#endif
    }


    template <class T1, int C1, class T2>
    void InstLU_SolveInPlace(
        const ConstMatrixView<T1,C1>& m1, const Permutation& P,
        MatrixView<T2> m2)
    { DoLUSolve<false>(m1,P,m2); }
    template <class T1, int C1, class T2>
    void InstLU_SolveTransposeInPlace(
        const ConstMatrixView<T1,C1>& m1, const Permutation& P,
        MatrixView<T2> m2)
    { DoLUSolve<true>(m1,P,m2); }
    template <class T1, int C1, class T2>
    void InstLU_SolveInPlace(
        const ConstMatrixView<T1,C1>& m1, const Permutation& P,
        VectorView<T2> v2)
    { DoLUSolve<false>(m1,P,v2); }
    template <class T1, int C1, class T2>
    void InstLU_SolveTransposeInPlace(
        const ConstMatrixView<T1,C1>& m1, const Permutation& P,
        VectorView<T2> v2)
    { DoLUSolve<true>(m1,P,v2); }

#define InstFile "TMV_LUDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


