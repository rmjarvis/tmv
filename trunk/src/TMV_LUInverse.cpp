
//#define PRINTALGO_LU

#include "TMV_Blas.h"
#include "tmv/TMV_LUInverse.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_DivMU.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_MultUL.h"
#include "tmv/TMV_PermuteM.h"

namespace tmv {

    template <class T1>
    static inline void NonLapLUInverse(MatrixView<T1> m1, const Permutation& P)
    { InlineLU_Inverse(m1,P); }

#ifdef ALAP
    // ALAP, not LAP, since ATLAS has these routines
    template <class T1> 
    static inline void LapLUInverse(
        MatrixView<T1,ColMajor> m1, const Permutation& P)
    { NonLapLUInverse(m1,P); }
#ifdef INST_DOUBLE
    template <>
    static void LapLUInverse(
        MatrixView<double,ColMajor> m1, const Permutation& P)
    {
        TMVAssert(P.isInverse());
        int n = m1.colsize();
        int lda = m1.stepj();
#ifdef CLAP
        const int* ipiv = P.getValues();
#else
        auto_array<int> ipiv1(new int[n]);
        const int* ipiv = ipiv1.get();
        for(int i=0;i<n;++i) ipiv1[i] = P.getValues()[i]+1;
#endif
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = n*LAP_BLOCKSIZE;
        auto_array<double> work(new double[lwork]);
#else
        int lwork = -1;
        auto_array<double> work(new double[1]);
        LAPNAME(dgetri) (
            LAPCM LAPV(n),LAPP(m1.ptr()),LAPV(lda),
            LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(work[0]);
        work.reset(new double[lwork]);
#endif
#endif
        LAPNAME(dgetri) (
            LAPCM LAPV(n),LAPP(m1.ptr()),LAPV(lda),
            LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results("dgetri");
#else
        LAP_Results(int(work[0]),n,n,lwork,"dgetri");
#endif
    }
    template <>
    static void LapLUInverse(
        MatrixView<std::complex<double>,ColMajor> m1, const Permutation& P)
    {
        TMVAssert(P.isInverse());
        int n = m1.colsize();
        int lda = m1.stepj();
#ifdef CLAP
        const int* ipiv = P.getValues();
#else
        auto_array<int> ipiv1(new int[n]);
        const int* ipiv = ipiv1.get();
        for(int i=0;i<n;++i) ipiv1[i] = P.getValues()[i]+1;
#endif
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = n*LAP_BLOCKSIZE;
        auto_array<std::complex<double> > work(new std::complex<double>[lwork]);
#else
        int lwork = -1;
        auto_array<std::complex<double> > work(new std::complex<double>[1]);
        LAPNAME(zgetri) (
            LAPCM LAPV(n),LAPP(m1.ptr()),LAPV(lda),
            LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(std::real(work[0]));
        work.reset(new std::complex<double>[lwork]);
#endif
#endif
        LAPNAME(zgetri) (
            LAPCM LAPV(n),LAPP(m1.ptr()),LAPV(lda),
            LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results("zgetri");
#else
        LAP_Results(int(std::real(work[0])),n,n,lwork,"zgetri");
#endif
    }
#endif
#ifdef INST_FLOAT
    template <>
    static void LapLUInverse(
        MatrixView<float,ColMajor> m1, const Permutation& P)
    {
        TMVAssert(P.isInverse());
        int n = m1.colsize();
        int lda = m1.stepj();
#ifdef CLAP
        const int* ipiv = P.getValues();
#else
        auto_array<int> ipiv1(new int[n]);
        const int* ipiv = ipiv1.get();
        for(int i=0;i<n;++i) ipiv1[i] = P.getValues()[i]+1;
#endif
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = n*LAP_BLOCKSIZE;
        auto_array<float> work(new float[lwork]);
#else
        int lwork = -1;
        auto_array<float> work(new float[1]);
        LAPNAME(sgetri) (
            LAPCM LAPV(n),LAPP(m1.ptr()),LAPV(lda),
            LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(work[0]);
        work.reset(new float[lwork]);
#endif
#endif
        LAPNAME(sgetri) (
            LAPCM LAPV(n),LAPP(m1.ptr()),LAPV(lda),
            LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results("sgetri");
#else
        LAP_Results(int(work[0]),n,n,lwork,"sgetri");
#endif
    }
    template <>
    static void LapLUInverse(
        MatrixView<std::complex<float>,ColMajor> m1, const Permutation& P)
    {
        TMVAssert(P.isInverse());
        int n = m1.colsize();
        int lda = m1.stepj();
#ifdef CLAP
        const int* ipiv = P.getValues();
#else
        auto_array<int> ipiv1(new int[n]);
        const int* ipiv = ipiv1.get();
        for(int i=0;i<n;++i) ipiv1[i] = P.getValues()[i]+1;
#endif
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = n*LAP_BLOCKSIZE;
        auto_array<std::complex<float> > work(new std::complex<float>[lwork]);
#else
        int lwork = -1;
        auto_array<std::complex<float> > work(new std::complex<float>[1]);
        LAPNAME(cgetri) (
            LAPCM LAPV(n),LAPP(m1.ptr()),LAPV(lda),
            LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(std::real(work[0]));
        work.reset(new std::complex<float>[lwork]);
#endif
#endif
        LAPNAME(cgetri) (
            LAPCM LAPV(n),LAPP(m1.ptr()),LAPV(lda),
            LAPP(ipiv.get()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results("cgetri");
#else
        LAP_Results(int(std::real(work[0])),n,n,lwork,"cgetri");
#endif
    }
#endif // FLOAT
#endif // ALAP

    template <class T1>
    void InstLU_Inverse(MatrixView<T1> m1, const Permutation& P)
    {
#ifdef ALAP
        if (m1.iscm() && m1.stepj() > 0) {
            LapLUInverse(m1.cmView(),P);
        } else {
            Matrix<T1,ColMajor|NoDivider> m1c(m1);
            LapLUInverse(m1c.view(),P);
            InstCopy(m1c.xView().constView(),m1);
        }
#else
        NonLapLUInverse(m1,P);
#endif
    }

    template <class T1, class T2, int C1>
    void InstLU_InverseATA(
        const ConstMatrixView<T1,C1>& m1, const Permutation& P,
        const bool trans, MatrixView<T2> m2)
    {
        InlineLU_InverseATA(m1,P,trans,m2);
    }

#define InstFile "TMV_LUInverse.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


