
//#define PRINTALGO_QR
//#define XDEBUG_QR

#include "TMV_Blas.h"
#include "tmv/TMV_QRDecompose.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_SmallMatrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_SmallTriMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_SmallVector.h"
#include "tmv/TMV_CopyV.h"
#include "tmv/TMV_SwapV.h"
#include "tmv/TMV_MinMax.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_Rank1VVM.h"
#include "tmv/TMV_MultMV.h"
#include "tmv/TMV_NormM.h"
#include "tmv/TMV_MultMM.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_MultUM.h"

#include "tmv/TMV_ProdVV.h"
#include "tmv/TMV_ProdXM.h"
#include "tmv/TMV_ProdMM.h"
#include "tmv/TMV_ProdMV.h"

namespace tmv {

    template <class M, class V> 
    static inline void DoQR_Decompose(M& A, V& beta)
    {
        typename M::cmview_type Acm = A.cmView();
        typename V::unitview_type betau = beta.unitView();
        InlineQR_Decompose(Acm,betau); 
    }

#ifdef LAP
#ifdef TMV_INST_DOUBLE
    void DoQR_Decompose(MatrixView<double> A, VectorView<double> beta)
    {
        int m = A.colsize();
        int n = A.rowsize();
        beta.setZero();
        if (A.iscm()) {
            int lda = A.stepj();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = 2*n*LAP_BLOCKSIZE;
            AlignedArray<double> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<double> work(1);
            work.get()[0] = 0.;
            LAPNAME(dgeqrf) (
                LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
                LAPP(beta.ptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(work[0]);
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(dgeqrf) (
                LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
                LAPP(beta.ptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"dgeqrf");
#else
            LAP_Results(Lap_info,int(work[0]),m,n,lwork,"dgeqrf");
#endif
        } else {
            int lda = A.stepi();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = 2*n*LAP_BLOCKSIZE;
            AlignedArray<double> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<double> work(1);
            work.get()[0] = 0.;
            LAPNAME(dgelqf) (
                LAPCM LAPV(n),LAPV(m),LAPP(A.ptr()),LAPV(lda),
                LAPP(beta.ptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(work[0]);
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(dgelqf) (
                LAPCM LAPV(n),LAPV(m),LAPP(A.ptr()),LAPV(lda),
                LAPP(beta.ptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"dgelqf");
#else
            LAP_Results(Lap_info,int(work[0]),m,n,lwork,"dgelqf");
#endif
        }
        double* bi = beta.ptr();
        for(int i=0;i<n;++i,++bi)  {
            if (TMV_ABS(*bi-1.) > 1.01) *bi = 0.;
        }
    }
    // Note: No complex decomposition, since LAPACK uses a complex
    // beta vector, which we do not.  So the complex decomposition is
    // always done with the TMV code.
#endif
#ifdef TMV_INST_FLOAT
    void DoQR_Decompose(MatrixView<float> A, VectorView<float> beta)
    {
        int m = A.colsize();
        int n = A.rowsize();
        beta.setZero();
        if (A.iscm()) {
            int lda = A.stepj();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = 2*n*LAP_BLOCKSIZE;
            AlignedArray<float> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<float> work(1);
            work.get()[0] = 0.F;
            LAPNAME(sgeqrf) (
                LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
                LAPP(beta.ptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(work[0]);
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(sgeqrf) (
                LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
                LAPP(beta.ptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"sgeqrf");
#else
            LAP_Results(Lap_info,int(work[0]),m,n,lwork,"sgeqrf");
#endif
        } else {
            int lda = A.stepi();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = 2*n*LAP_BLOCKSIZE;
            AlignedArray<float> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<float> work(1);
            work.get()[0] = 0.F;
            LAPNAME(sgelqf) (
                LAPCM LAPV(n),LAPV(m),LAPP(A.ptr()),LAPV(lda),
                LAPP(beta.ptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(work[0]);
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(sgelqf) (
                LAPCM LAPV(n),LAPV(m),LAPP(A.ptr()),LAPV(lda),
                LAPP(beta.ptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"sgelqf");
#else
            LAP_Results(Lap_info,int(work[0]),m,n,lwork,"sgelqf");
#endif
        }
        float* bi = beta.ptr();
        for(int i=0;i<n;++i,++bi) {
            if (TMV_ABS(*bi-1.F) > 1.01F) *bi = 0.F;
        }
    }
#endif
#endif // LAP

    template <class T, class RT> 
    void InstQR_Decompose(MatrixView<T> A, VectorView<RT> beta)
    {
        if (beta.step() == 1) {
            if (A.iscm()) {
                MatrixView<T,ColMajor> Acm = A;
                DoQR_Decompose(Acm,beta);
            } else {
                Matrix<T,ColMajor|NoDivider> Ac = A;
                InstQR_Decompose(Ac.xView(),beta);
                InstCopy(Ac.constView().xView(),A);
            }
        } else {
            Vector<RT> betac = beta;
            InstQR_Decompose(A,betac.xView());
            InstCopy(betac.constView().xView(),beta);
        }
    }

#define InstFile "TMV_QRDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


