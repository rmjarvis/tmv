
//#define PRINTALGO_LU

#include "TMV_Blas.h"
#include "tmv/TMV_LUDecompose.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_SmallTriMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_SmallVector.h"
#include "tmv/TMV_CopyV.h"
#include "tmv/TMV_SwapV.h"
#include "tmv/TMV_MinMax.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_PermuteM.h"
#include "tmv/TMV_Det.h"

namespace tmv {

#ifdef ALAP
    template <class T> 
    static inline void LapLU_Decompose(
        MatrixView<T,ColMajor>& A, int* P)
    { InlineLU_Decompose(A,P); }
#ifdef TMV_INST_DOUBLE
    static void LapLU_Decompose(
        MatrixView<double,ColMajor>& A, int* P)
    {
        TMVAssert(A.iscm());

        int m = A.colsize();
        int n = A.rowsize();
        int lda = A.stepj();
        int* lap_p = new int[n];
        LAPNAME(dgetrf) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p) LAPINFO);
        LAP_Results("dgetrf");
        const int M = A.colsize();
        for(int i=0;i<M;i++) {
            P[i] = lap_p[i] LAPMINUS1;
        }
        delete [] lap_p;
    }
    static void LapLU_Decompose(
        MatrixView<std::complex<double>,ColMajor>& A, int* P)
    {
        TMVAssert(A.iscm());

        int m = A.colsize();
        int n = A.rowsize();
        int lda = A.stepj();
        int* lap_p = new int[n];
        LAPNAME(zgetrf) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p) LAPINFO);
        LAP_Results("zgetrf");
        const int M = A.colsize();
        for(int i=0;i<M;i++) {
            P[i] = lap_p[i] LAPMINUS1;
        }
        delete [] lap_p;
    }
#endif
#ifdef TMV_INST_FLOAT
#ifndef MKL
    // This is giving me a weird runtime error sometimes with MKL:
    //   OMP abort: Unable to set worker thread stack size to 2098176 bytes
    //   Try reducing KMP_STACKSIZE or increasing the shell stack limit.
    // So I'm cutting it out for MKL compilations
    static void LapLU_Decompose(
        MatrixView<float,ColMajor>& A, int* P)
    {
        int m = A.colsize();
        int n = A.rowsize();
        int lda = A.stepj();
        int* lap_p = new int[n];
        LAPNAME(sgetrf) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p) LAPINFO);
        LAP_Results("sgetrf");
        const int M = A.colsize();
        for(int i=0;i<M;i++) {
            P[i] = lap_p[i] LAPMINUS1;
        }
        delete [] lap_p;
    }
    static void LapLU_Decompose(
        MatrixView<std::complex<float>,ColMajor>& A, int* P)
    {
        int m = A.colsize();
        int n = A.rowsize();
        int lda = A.stepj();
        int* lap_p = new int[n];
        LAPNAME(cgetrf) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p) LAPINFO);
        LAP_Results("cgetrf");
        const int M = A.colsize();
        for(int i=0;i<M;i++) {
            P[i] = lap_p[i] LAPMINUS1;
        }
        delete [] lap_p;
    }
#endif // MKL
#endif // FLOAT
#endif // ALAP

    template <class T> 
    void InstLU_Decompose(MatrixView<T> A, int* P)
    {
        if (A.colsize() > 0 && A.rowsize() > 0) {
            if (A.iscm()) {
                MatrixView<T,ColMajor> Acm = A;
#ifdef ALAP
                LapLU_Decompose(Acm,P);
#else
                InlineLU_Decompose(Acm,P);
#endif
            } else {
                Matrix<T,ColMajor|NoDivider> Ac = A;
                InstLU_Decompose(Ac.xView(),P);
                InstCopy(Ac.constView().xView(),A);
            }
        }
    }

#define InstFile "TMV_LUDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


