
//#define PRINTALGO_QR
//#define XDEBUG_QR

#ifdef NOGEQP3
#ifdef LAP
#undef LAP
#endif
#endif

#include "TMV_Blas.h"
#include "tmv/TMV_QRPDecompose.h"
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
#include "tmv/TMV_Norm.h"
#include "tmv/TMV_MultMM.h"
#include "tmv/TMV_CopyU.h"
#include "tmv/TMV_MultUM.h"

#include "tmv/TMV_SumVV.h"
#include "tmv/TMV_ProdMM.h"
#include "tmv/TMV_ProdMV.h"

#ifdef LAP
#include "tmv/TMV_SortV.h"
#endif

namespace tmv {

    template <class M, class V> 
    static inline void DoQRP_Decompose(M& A, V& beta, ptrdiff_t* P, bool strict)
    { 
        typename M::cmview_type Acm = A.cmView();
        typename V::unitview_type betau = beta.unitView();
        InlineQRP_Decompose(Acm,betau,P,strict);
    }

#ifdef LAP
#ifdef TMV_INST_DOUBLE
    void DoQRP_Decompose(
        MatrixView<double> A, VectorView<double> beta, ptrdiff_t* P, bool strict)
    {
        int m = A.colsize();
        int n = A.rowsize();
        AlignedArray<int> lap_p(n);
        for(int i=0;i<n;++i) (lap_p.get())[i] = 0;
        beta.setZero();
        int lda = A.stepj();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = 3*n+1;
        AlignedArray<double> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#else
        int lwork = -1;
        AlignedArray<double> work(1);
        work.get()[0] = 0.;
        LAPNAME(dgeqp3) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(beta.ptr()) 
            LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(work[0]);
        work.resize(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
        LAPNAME(dgeqp3) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(beta.ptr()) 
            LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results("dgeqp3");
#else
        LAP_Results(int(work[0]),m,n,lwork,"dgeqp3");
#endif
        double thresh = TMV_Epsilon<double>()*A.normF();
        for(int i=0;i<n;++i) {
            if (TMV_ABS(A(i,i)) < thresh) {
                TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
                A.subMatrix(i,A.colsize(),i,A.rowsize()).setZero();
                beta.subVector(i,n).setZero();
                break;
            }
        }
#ifndef CLAP
        for(int i=0;i<n;++i) --((lap_p.get())[i]);
#endif
        ConvertIndexToPermute(n,lap_p.get(),P);
        double* bi = beta.ptr();
        for(int i=0;i<n;++i,++bi) {
            if (TMV_ABS(*bi-1.) > 1.01) *bi = 0.;
        }
    }
#endif
#ifdef TMV_INST_FLOAT
    void DoQRP_Decompose(
        MatrixView<float> A, VectorView<float> beta, ptrdiff_t* P, bool strict)
    {
        int m = A.colsize();
        int n = A.rowsize();
        AlignedArray<int> lap_p(n);
        for(int i=0;i<n;++i) (lap_p.get())[i] = 0;
        beta.setZero();
        int lda = A.stepj();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = 3*n+1;
        AlignedArray<float> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#else
        int lwork = -1;
        AlignedArray<float> work(1);
        work.get()[0] = 0.;
        LAPNAME(sgeqp3) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(beta.ptr()) 
            LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(work[0]);
        work.resize(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
        LAPNAME(sgeqp3) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
            LAPP(lap_p.get()),LAPP(beta.ptr()) 
            LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
        LAP_Results("sgeqp3");
#else
        LAP_Results(int(work[0]),m,n,lwork,"sgeqp3");
#endif
        float thresh = TMV_Epsilon<float>()*A.normF();
        for(int i=0;i<n;++i) {
            if (TMV_ABS(A(i,i)) < thresh) {
                TMVAssert(NormInf(A.diag(0,i+1,A.rowsize())) < thresh);
                A.subMatrix(i,A.colsize(),i,A.rowsize()).setZero();
                beta.subVector(i,n).setZero();
                break;
            }
        }
#ifndef CLAP
        for(int i=0;i<n;++i) --((lap_p.get())[i]);
#endif
        ConvertIndexToPermute(n,lap_p.get(),P);
        float* bi = beta.ptr();
        for(int i=0;i<n;++i,++bi) {
            if (TMV_ABS(*bi-1.F) > 1.01F) *bi = 0.F;
        }
    }
#endif // FLOAT
#endif // LAP

    template <class T, class RT> 
    void InstQRP_Decompose(
        MatrixView<T> A, VectorView<RT> beta, ptrdiff_t* P, bool strict)
    {
        if (A.rowsize() > 0) {
            if (beta.step() == 1) {
                VectorView<RT,Unit> beta1 = beta;
                if (A.iscm()) {
                    MatrixView<T,ColMajor> Acm = A;
                    DoQRP_Decompose(Acm,beta1,P,strict);
                } else {
                    Matrix<T,ColMajor|NoDivider> Ac = A;
                    InstQRP_Decompose(Ac.xView(),beta,P,strict);
                    InstCopy(Ac.constView().xView(),A);
                }
            } else {
                Vector<RT> betac = beta;
                InstQRP_Decompose(A,betac.xView(),P,strict);
                InstCopy(betac.constView().xView(),beta);
            }
        }
    }

#define InstFile "TMV_QRPDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


