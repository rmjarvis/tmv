
//#define PRINTALGO_LU
//#define XDEBUG_LU

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

    template <class M1>
    static inline void DoLUInverse(M1& m1, const Permutation& P)
    { InlineLU_Inverse(m1,P); }

#ifdef ALAP
    // ALAP, not LAP, since ATLAS has these routines
#ifdef TMV_INST_DOUBLE
    static void DoLUInverse(
        MatrixView<double,ColMajor> m1, const Permutation& P)
    {
        TMVAssert(P.isInverse());
        int n = m1.colsize();
        int lda = m1.stepj();
#ifdef CLAP
        const int* ipiv = P.getValues();
#else
        AlignedArray<int> ipiv1(n);
        for(int i=0;i<n;++i) ipiv1[i] = P.getValues()[i]+1;
        const int* ipiv = ipiv1.get();
#endif
        int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = n*LAP_BLOCKSIZE;
        AlignedArray<double> work(lwork);
#else
        int lwork = -1;
        AlignedArray<double> work(1);
        LAPNAME(dgetri) (
            LAPCM LAPV(n),LAPP(m1.ptr()),LAPV(lda),
            LAPP(ipiv) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(work[0]);
        work.resize(lwork);
#endif
#endif
        LAPNAME(dgetri) (
            LAPCM LAPV(n),LAPP(m1.ptr()),LAPV(lda),
            LAPP(ipiv) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results(Lap_info,"dgetri");
#else
        LAP_Results(Lap_info,int(work[0]),n,n,lwork,"dgetri");
#endif
    }
    static void DoLUInverse(
        MatrixView<std::complex<double>,ColMajor> m1, const Permutation& P)
    {
        TMVAssert(P.isInverse());
        int n = m1.colsize();
        int lda = m1.stepj();
#ifdef CLAP
        const int* ipiv = P.getValues();
#else
        AlignedArray<int> ipiv1(n);
        for(int i=0;i<n;++i) ipiv1[i] = P.getValues()[i]+1;
        const int* ipiv = ipiv1.get();
#endif
        int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = n*LAP_BLOCKSIZE;
        AlignedArray<std::complex<double> > work(lwork);
#else
        int lwork = -1;
        AlignedArray<std::complex<double> > work(1);
        LAPNAME(zgetri) (
            LAPCM LAPV(n),LAPP(m1.ptr()),LAPV(lda),
            LAPP(ipiv) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(std::real(work[0]));
        work.resize(lwork);
#endif
#endif
        LAPNAME(zgetri) (
            LAPCM LAPV(n),LAPP(m1.ptr()),LAPV(lda),
            LAPP(ipiv) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results(Lap_info,"zgetri");
#else
        LAP_Results(Lap_info,int(std::real(work[0])),n,n,lwork,"zgetri");
#endif
    }
#endif
#ifdef TMV_INST_FLOAT
    static void DoLUInverse(
        MatrixView<float,ColMajor> m1, const Permutation& P)
    {
        TMVAssert(P.isInverse());
        int n = m1.colsize();
        int lda = m1.stepj();
#ifdef CLAP
        const int* ipiv = P.getValues();
#else
        AlignedArray<int> ipiv1(n);
        for(int i=0;i<n;++i) ipiv1[i] = P.getValues()[i]+1;
        const int* ipiv = ipiv1.get();
#endif
        int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = n*LAP_BLOCKSIZE;
        AlignedArray<float> work(lwork);
#else
        int lwork = -1;
        AlignedArray<float> work(1);
        LAPNAME(sgetri) (
            LAPCM LAPV(n),LAPP(m1.ptr()),LAPV(lda),
            LAPP(ipiv) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(work[0]);
        work.resize(lwork);
#endif
#endif
        LAPNAME(sgetri) (
            LAPCM LAPV(n),LAPP(m1.ptr()),LAPV(lda),
            LAPP(ipiv) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results(Lap_info,"sgetri");
#else
        LAP_Results(Lap_info,int(work[0]),n,n,lwork,"sgetri");
#endif
    }
    static void DoLUInverse(
        MatrixView<std::complex<float>,ColMajor> m1, const Permutation& P)
    {
        TMVAssert(P.isInverse());
        int n = m1.colsize();
        int lda = m1.stepj();
#ifdef CLAP
        const int* ipiv = P.getValues();
#else
        AlignedArray<int> ipiv1(n);
        for(int i=0;i<n;++i) ipiv1[i] = P.getValues()[i]+1;
        const int* ipiv = ipiv1.get();
#endif
        int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = n*LAP_BLOCKSIZE;
        AlignedArray<std::complex<float> > work(lwork);
#else
        int lwork = -1;
        AlignedArray<std::complex<float> > work(1);
        LAPNAME(cgetri) (
            LAPCM LAPV(n),LAPP(m1.ptr()),LAPV(lda),
            LAPP(ipiv) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(std::real(work[0]));
        work.resize(lwork);
#endif
#endif
        LAPNAME(cgetri) (
            LAPCM LAPV(n),LAPP(m1.ptr()),LAPV(lda),
            LAPP(ipiv) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results(Lap_info,"cgetri");
#else
        LAP_Results(Lap_info,int(std::real(work[0])),n,n,lwork,"cgetri");
#endif
    }
#endif // FLOAT
#endif // ALAP

    template <class T1>
    void InstLU_Inverse(MatrixView<T1> m1, const Permutation& P)
    {
#ifdef ALAP
        if (m1.iscm() && m1.stepj() > 0) {
            DoLUInverse(m1.cmView(),P);
        } else {
            Matrix<T1,ColMajor|NoDivider> m1c(m1);
            DoLUInverse(m1c.view(),P);
            InstCopy(m1c.xView().constView(),m1);
        }
#else
        DoLUInverse(m1,P);
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


