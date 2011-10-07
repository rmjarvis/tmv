
//#undef NDEBUG
//#define PRINTALGO_SVD
//#define XDEBUG_SVD
//#include "TMV.h"

#include "TMV_Blas.h"
#include "tmv/TMV_SVDecompose_Bidiag.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_MultMM.h"

#include "tmv/TMV_ProdMM.h"
#include "tmv/TMV_ProdMV.h"
#include "tmv/TMV_SumMM.h"

namespace tmv {

    template <class T, class RT>
    void DoBidiagonalize(
        MatrixView<T> A, VectorView<RT> Ubeta, VectorView<RT> Vbeta,
        VectorView<T> D, VectorView<T> E)
    {
        MatrixView<T,ColMajor> Acm = A;
        VectorView<RT,Unit> Ubu = Ubeta;
        VectorView<RT,Unit> Vbu = Vbeta;
        VectorView<T,Unit> Du = D;
        VectorView<T,Unit> Eu = E;
        InlineBidiagonalize(Acm,Ubu,Vbu,Du,Eu);
    }

#ifdef LAP
#ifdef TMV_INST_DOUBLE
    void DoBidiagonalize(
        MatrixView<double>& A, VectorView<double>& Ubeta,
        VectorView<double>& Vbeta, VectorView<double>& D,
        VectorView<double>& E)
    {
        int m = A.colsize();
        int n = A.rowsize();
        int ldu = A.stepj();
        D.setZero();
        E.setZero();
        Ubeta.setZero();
        Vector<double> Vbeta2(n,0.);  
        // Stupid LAPACK requires an extra element in the Vbeta vector
        // which it sets to 0 (!!!) rather than ignores.
        // So we need to create a temporary Vector which is size n.
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = (m+n)*LAP_BLOCKSIZE;
        AlignedArray<double> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#else
        int lwork = -1;
        AlignedArray<double> work(1);
        work.get()[0] = 0.;
        LAPNAME(dgebrd) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(ldu),
            LAPP(D.ptr()),LAPP(E.ptr()),LAPP(Ubeta.ptr()),LAPP(Vbeta2.ptr())
            LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(work[0]);
        work.resize(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
        //std::cout<<"Before dgebrd"<<std::endl;
        LAPNAME(dgebrd) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(ldu),
            LAPP(D.ptr()),LAPP(E.ptr()),LAPP(Ubeta.ptr()),LAPP(Vbeta2.ptr())
            LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        //std::cout<<"After dgebrd"<<std::endl;
        Vbeta = Vbeta2.subVector(0,n-1);
#ifdef LAPNOWORK
        LAP_Results("dgebrd");
#else
        LAP_Results(int(work[0]),m,n,lwork,"dgebrd");
#endif
    }
#endif
#ifdef TMV_INST_FLOAT
    void DoBidiagonalize(
        MatrixView<float>& A, VectorView<float>& Ubeta,
        VectorView<float>& Vbeta, VectorView<float>& D,
        VectorView<float>& E)
    {
        int m = A.colsize();
        int n = A.rowsize();
        int ldu = A.stepj();
        D.setZero();
        E.setZero();
        Ubeta.setZero();
        Vector<float> Vbeta2(n,0.F);  
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
        int lwork = (m+n)*LAP_BLOCKSIZE;
        AlignedArray<float> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#else
        int lwork = -1;
        AlignedArray<float> work(1);
        work.get()[0] = 0.F;
        LAPNAME(sgebrd) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(ldu),
            LAPP(D.ptr()),LAPP(E.ptr()),LAPP(Ubeta.ptr()),LAPP(Vbeta2.ptr())
            LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        lwork = int(work[0]);
        work.resize(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
        //std::cout<<"Before sgebrd"<<std::endl;
        LAPNAME(sgebrd) (
            LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(ldu),
            LAPP(D.ptr()),LAPP(E.ptr()),LAPP(Ubeta.ptr()),LAPP(Vbeta2.ptr())
            LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
        //std::cout<<"After sgebrd"<<std::endl;
        Vbeta = Vbeta2.subVector(0,n-1);
#ifdef LAPNOWORK
        LAP_Results("sgebrd");
#else
        LAP_Results(int(work[0]),m,n,lwork,"sgebrd");
#endif
    }
#endif 
#endif // LAP

    template <class T, class RT>
    void InstBidiagonalize(
        MatrixView<T> A, VectorView<RT> Ubeta, VectorView<RT> Vbeta,
        VectorView<T> D, VectorView<T> E)
    {
        TMVAssert(A.iscm());
        TMVAssert(Ubeta.step() == 1);
        TMVAssert(Vbeta.step() == 1);
        TMVAssert(D.step() == 1);
        TMVAssert(E.step() == 1);
        if (D.size() > 0) {
            DoBidiagonalize(A,Ubeta,Vbeta,D,E);
        }
    }

#define InstFile "TMV_SVDecompose_Bidiag.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


