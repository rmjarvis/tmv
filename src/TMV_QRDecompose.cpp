///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2016                                                 //
// All rights reserved                                                       //
//                                                                           //
// The project is hosted at https://code.google.com/p/tmv-cpp/               //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis17 [at] gmail.                                                 //
//                                                                           //
// This software is licensed under a FreeBSD license.  The file              //
// TMV_LICENSE should have bee included with this distribution.              //
// It not, you can get a copy from https://code.google.com/p/tmv-cpp/.       //
//                                                                           //
// Essentially, you can use this software however you want provided that     //
// you include the TMV_LICENSE file in any distribution that uses it.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


//#define XDEBUG


#include "TMV_Blas.h"
#include "TMV_QRDiv.h"
#include "tmv/TMV_QRD.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_Householder.h"
#include "tmv/TMV_TriMatrixArith.h"
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define QR_BLOCKSIZE TMV_BLOCKSIZE
#else
#define QR_BLOCKSIZE 64
#endif

    //
    // QR Decompose
    //

    template <class T> 
    static void NonBlockQRDecompose(
        MatrixView<T> A, VectorView<T> beta, T& det)
    {
#ifdef XDEBUG
        cout<<"Start NonBlockQRDecompose\n";
        cout<<"Norm(A) = "<<Norm(A)<<endl;
        cout<<"Norm(beta) = "<<Norm(A)<<endl;
        Matrix<T> A0(A);
#endif

        TMVAssert(A.colsize() >= A.rowsize());
        TMVAssert(A.rowsize() == beta.size());
        TMVAssert(A.rowsize() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(!beta.isconj());
        TMVAssert(beta.step()==1);
        // Decompose A into A = Q R 
        // where Q is unitary, and R is upper triangular
        // Q and R are stored in the same matrix (output of A), 
        // with the beta's for the Householder matrices returned in beta.
        const ptrdiff_t M = A.colsize();
        const ptrdiff_t N = A.rowsize();

        T* bj = beta.ptr();
        for(ptrdiff_t j=0;j<N;++j,++bj) {
            // Apply the Householder Reflection for this column
#ifdef TMVFLDEBUG
            TMVAssert(bj >= beta._first);
            TMVAssert(bj < beta._last);
#endif
            *bj = HouseholderReflect(A.subMatrix(j,M,j,N),det);
        }
#ifdef XDEBUG
        Matrix<T> R = A.upperTri();
        Matrix<T> Q(A);
        GetQFromQR(Q.view(),beta);
        Matrix<T> AA = Q*R;
        cout<<"Done NonBlockQRDecompose\n";
        cout<<"Norm(A) = "<<Norm(A)<<endl;
        cout<<"Norm(beta) = "<<Norm(A)<<endl;
        cout<<"Norm(AA-A0) = "<<Norm(AA-A0)<<endl;
        if (!(Norm(AA-A0) < 0.0001*Norm(Q)*Norm(R))) {
            cerr<<"NonBlockQRDecompose: A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"-> "<<A<<endl;
            cerr<<"beta = "<<beta<<endl;
            cerr<<"Q = "<<Q<<endl;
            cerr<<"R = "<<R<<endl;
            cerr<<"Q*R = "<<Q*R<<endl;
            cerr<<"Norm(diff) = "<<Norm(AA-A0)<<endl;
            abort(); 
        }
#endif
    }

    template <class T> 
    static void RecursiveQRDecompose(
        MatrixView<T> A, UpperTriMatrixView<T> Z, T& det,
        bool makeZ)
    {
#ifdef XDEBUG
        cout<<"Start RecursiveQRDecompose\n";
        cout<<"A = "<<A<<endl;
        cout<<"Norm(A) = "<<Norm(A)<<endl;
        cout<<"Norm(beta) = "<<Norm(A)<<endl;
        Matrix<T> A0(A);
#endif
        TMVAssert(A.colsize() >= A.rowsize());
        TMVAssert(A.rowsize() == Z.size());
        TMVAssert(A.rowsize() > 0);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(Z.iscm());
        TMVAssert(!A.isconj());
        TMVAssert(!Z.isconj());
        // This is very similar to the BlockHouseholder_MakeZ function
        // in Householder.cpp.  The difference is the addition of the 
        // Householder_Reflects.
        // The makeZ parameter should be set to true if you want the Z
        // matrix to be correct on output.  If you don't need the Z matrix
        // after this call, setting makeZ to false will speed it up slightly.
        // (In either case, the diagonal of Z is correctly set to be 
        // beta.conjugate().)
        const ptrdiff_t M = A.colsize();
        const ptrdiff_t N = A.rowsize();

        if (N==1) {
            T b = HouseholderReflect(A.col(0),det);
            *Z.ptr() = TMV_CONJ(b);
        } else if (N==2) {
            T* Z00 = Z.ptr();
            T* Z01 = Z00 + Z.stepj();
            T* Z11 = Z01 + 1;

            T b0 = HouseholderReflect(A,det);
            *Z00 = TMV_CONJ(b0);
            T b1 = HouseholderReflect(A.col(1,1,M),det);
            *Z11 = TMV_CONJ(b1);

            if (makeZ) {
                const T* A10 = A.cptr()+A.stepi();
                T temp = A.col(0,2,M).conjugate()*A.col(1,2,M);
                temp += TMV_CONJ(*A10);
                *Z01 = -TMV_CONJ(b0*b1)*temp;
            }
        } else {
            ptrdiff_t j1 = (N+1)/2;
            MatrixView<T> A1 = A.colRange(0,j1);
            UpperTriMatrixView<T> Z1 = Z.subTriMatrix(0,j1);
            RecursiveQRDecompose(A1,Z1,det,true);

            BlockHouseholderLDiv(A1,Z1,A.colRange(j1,N));

            MatrixView<T> A2 = A.subMatrix(j1,M,j1,N);
            UpperTriMatrixView<T> Z2 = Z.subTriMatrix(j1,N);
            RecursiveQRDecompose(A2,Z2,det,makeZ);

            if (makeZ) {
                MatrixView<T> Z3 = Z.subMatrix(0,j1,j1,N);
                Z3 = A1.rowRange(j1,N).adjoint() *
                    A.subMatrix(j1,N,j1,N).lowerTri(UnitDiag);
                Z3 += A1.rowRange(N,M).adjoint() * A.subMatrix(N,M,j1,N);
                Z3 = -Z1*Z3;
                Z3 *= Z2;
            }
        }
    }

    template <class T> 
    static void BlockQRDecompose(
        MatrixView<T> A, VectorView<T> beta, T& det)
    {
        TMVAssert(A.colsize() >= A.rowsize());
        TMVAssert(A.rowsize() == beta.size());
        TMVAssert(A.rowsize() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(beta.ct() == NonConj);
        TMVAssert(beta.step()==1);

#ifdef XDEBUG
        cout<<"Start BlockQRDecompose\n";
        cout<<"Norm(A) = "<<Norm(A)<<endl;
        cout<<"Norm(beta) = "<<Norm(A)<<endl;
        Matrix<T> A0(A);
#endif
        const ptrdiff_t M = A.colsize();
        const ptrdiff_t N = A.rowsize();

        UpperTriMatrix<T,NonUnitDiag|ColMajor> BaseZ(
            TMV_MIN(QR_BLOCKSIZE,int(N)));
        for(ptrdiff_t j1=0;j1<N;) {
            ptrdiff_t j2 = TMV_MIN(N,j1+QR_BLOCKSIZE);
            MatrixView<T> A1 = A.subMatrix(j1,M,j1,j2);
            UpperTriMatrixView<T> Z = BaseZ.subTriMatrix(0,j2-j1);

            RecursiveQRDecompose(A1,Z,det,j2<N);
            beta.subVector(j1,j2) = Z.diag().conjugate();

            if (j2 < N) 
                BlockHouseholderLDiv(A1,Z,A.subMatrix(j1,M,j2,N));
            j1 = j2;
        }
#ifdef XDEBUG
        Matrix<T> R = A.upperTri();
        Matrix<T> Q(A);
        GetQFromQR(Q.view(),beta);
        Matrix<T> AA = Q*R;
        cout<<"Done BlockQRDecompose\n";
        cout<<"Norm(A) = "<<Norm(A)<<endl;
        cout<<"Norm(beta) = "<<Norm(A)<<endl;
        cout<<"Norm(AA-A0) = "<<Norm(AA-A0)<<endl;
        cout<<"Norm(BaseZ) = "<<Norm(BaseZ)<<endl;
        if (!(Norm(AA-A0) < 0.0001*Norm(Q)*Norm(R))) {
            cerr<<"BlockQRDecompose: A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"-> "<<A<<endl;
            cerr<<"beta = "<<beta<<endl;
            cerr<<"Q*R = "<<Q*R<<endl;
            cerr<<"Norm(diff) = "<<Norm(AA-A0)<<endl;
            Matrix<T> A2 = A0;
            Vector<T> beta2(beta.size());
            T det2(0);
            NonBlockQRDecompose(A2.view(),beta2.view(),det2);
            cerr<<"NonBlock "<<A2<<endl;
            cerr<<"beta = "<<beta2<<endl;
            abort(); 
        }
#endif
    }

    template <class T> 
    static void NonLapQRDecompose(
        MatrixView<T> A, VectorView<T> beta, T& det)
    {
        TMVAssert(A.colsize() >= A.rowsize());
        TMVAssert(A.rowsize() == beta.size());
        TMVAssert(A.rowsize() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(beta.ct() == NonConj);
        TMVAssert(beta.step()==1);

        if (A.rowsize() > QR_BLOCKSIZE)
            BlockQRDecompose(A,beta,det);
        else {
            UpperTriMatrix<T,NonUnitDiag|ColMajor> Z(A.rowsize());
            RecursiveQRDecompose(A,Z.view(),det,false);
            beta = Z.diag().conjugate();
        }
    }

#ifdef LAP
    template <class T> 
    static inline void LapQRDecompose(
        MatrixView<T> A, VectorView<T> beta, T& det)
    { NonLapQRDecompose(A,beta,det); }
#ifdef INST_DOUBLE
    template <> 
    void LapQRDecompose(
        MatrixView<double> A,
        VectorView<double> beta, double& det)
    {
        TMVAssert(A.colsize() >= A.rowsize());
        TMVAssert(A.rowsize() == beta.size());
        TMVAssert(A.rowsize() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(beta.ct() == NonConj);
        TMVAssert(beta.step()==1);
        int m = A.colsize();
        int n = A.rowsize();
        beta.setZero();
        if (BlasIsCM(A)) {
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
            if (det && *bi != 0.) det = -det;
        }
    }
    template <> 
    void LapQRDecompose(
        MatrixView<std::complex<double> > A,
        VectorView<std::complex<double> > beta,
        std::complex<double>& det)
    {
        TMVAssert(A.colsize() >= A.rowsize());
        TMVAssert(A.rowsize() == beta.size());
        TMVAssert(A.rowsize() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(beta.ct() == NonConj);
        TMVAssert(beta.step()==1);
        int m = A.colsize();
        int n = A.rowsize();
        beta.setZero();
        if (BlasIsCM(A)) {
            int lda = A.stepj();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = 2*n*LAP_BLOCKSIZE;
            AlignedArray<std::complex<double> > work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<std::complex<double> > work(1);
            work.get()[0] = 0.;
            LAPNAME(zgeqrf) (
                LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
                LAPP(beta.ptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(TMV_REAL(work[0]));
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(zgeqrf) (
                LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
                LAPP(beta.ptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"zgeqrf");
#else
            LAP_Results(Lap_info,int(TMV_REAL(work[0])),m,n,lwork,"zgeqrf");
#endif
            beta.conjugateSelf();
        } else {
            int lda = A.stepi();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = 2*n*LAP_BLOCKSIZE;
            AlignedArray<std::complex<double> > work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<std::complex<double> > work(1);
            work.get()[0] = 0.;
            LAPNAME(zgelqf) (
                LAPCM LAPV(n),LAPV(m),LAPP(A.ptr()),LAPV(lda),
                LAPP(beta.ptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(TMV_REAL(work[0]));
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(zgelqf) (
                LAPCM LAPV(n),LAPV(m),LAPP(A.ptr()),LAPV(lda),
                LAPP(beta.ptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"zgelqf");
#else
            LAP_Results(Lap_info,int(TMV_REAL(work[0])),m,n,lwork,"zgelqf");
#endif
        }
        std::complex<double>* bi = beta.ptr();
        for(int i=0;i<n;++i,++bi) {
            // Sometimes (very rarely!) LAPACK comes back with an 
            // impossible beta value with |beta - 1| > 1.  e.g. (2,-1) (?!?)
            // This seems to only happen when beta should really be 0.  
            // So just set it to 0 now.
            if (TMV_ABS(*bi-1.) > 1.01) *bi = 0.;
            if (TMV_IMAG(*bi) != 0.) 
                det *= -TMV_CONJ(*bi * *bi)/norm(*bi);
            else if (TMV_REAL(*bi) != 0.)
                det = -det;
        }
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void LapQRDecompose(
        MatrixView<float> A,
        VectorView<float> beta, float& det)
    {
        TMVAssert(A.colsize() >= A.rowsize());
        TMVAssert(A.rowsize() == beta.size());
        TMVAssert(A.rowsize() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(beta.ct() == NonConj);
        TMVAssert(beta.step()==1);
        int m = A.colsize();
        int n = A.rowsize();
        beta.setZero();
        if (BlasIsCM(A)) {
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
            if (*bi != 0.F) det = -det;
        }
    }
    template <> 
    void LapQRDecompose(
        MatrixView<std::complex<float> > A,
        VectorView<std::complex<float> > beta, std::complex<float>& det)
    {
        TMVAssert(A.colsize() >= A.rowsize());
        TMVAssert(A.rowsize() == beta.size());
        TMVAssert(A.rowsize() > 0);
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(beta.ct() == NonConj);
        TMVAssert(beta.step()==1);
        int m = A.colsize();
        int n = A.rowsize();
        beta.setZero();
        if (BlasIsCM(A)) {
            int lda = A.stepj();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = 2*n*LAP_BLOCKSIZE;
            AlignedArray<std::complex<float> > work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<std::complex<float> > work(1);
            work.get()[0] = 0.F;
            LAPNAME(cgeqrf) (
                LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
                LAPP(beta.ptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(TMV_REAL(work[0]));
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(cgeqrf) (
                LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
                LAPP(beta.ptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"cgeqrf");
#else
            LAP_Results(Lap_info,int(TMV_REAL(work[0])),m,n,lwork,"cgeqrf");
#endif
            beta.conjugateSelf();
        } else {
            int lda = A.stepi();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = 2*n*LAP_BLOCKSIZE;
            AlignedArray<std::complex<float> > work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<std::complex<float> > work(1);
            work.get()[0] = 0.F;
            LAPNAME(cgelqf) (
                LAPCM LAPV(n),LAPV(m),LAPP(A.ptr()),LAPV(lda),
                LAPP(beta.ptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(TMV_REAL(work[0]));
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(cgelqf) (
                LAPCM LAPV(n),LAPV(m),LAPP(A.ptr()),LAPV(lda),
                LAPP(beta.ptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"cgelqf");
#else
            LAP_Results(Lap_info,int(TMV_REAL(work[0])),m,n,lwork,"cgelqf");
#endif
        }
        std::complex<float>* bi = beta.ptr();
        for(int i=0;i<n;++i,++bi) {
            if (TMV_ABS(*bi-1.F) > 1.01F) *bi = 0.F;
            if (TMV_IMAG(*bi) != 0.F) 
                det *= -TMV_CONJ(*bi * *bi)/norm(*bi);
            else if (TMV_REAL(*bi) != 0.F)
                det = -det;
        }
    }
#endif
#endif

    template <class T> 
    void QR_Decompose(MatrixView<T> A, VectorView<T> beta, T& det)
    {
#ifdef XDEBUG
        cout<<"Start QR_Decompose\n";
        cout<<"A = "<<TMV_Text(A)<<endl;
        cout<<"Norm(A) = "<<Norm(A)<<endl;
        cout<<"beta = "<<TMV_Text(beta)<<endl;
        cout<<"Norm(beta) = "<<Norm(beta)<<endl;
        Matrix<T> A0(A);
#endif

        TMVAssert(A.colsize() >= A.rowsize());
        TMVAssert(A.rowsize() == beta.size());
        TMVAssert(A.ct() == NonConj);
        TMVAssert(A.isrm() || A.iscm());
        TMVAssert(beta.ct() == NonConj);
        TMVAssert(beta.step()==1);
        if (A.rowsize() > 0) {
#ifdef LAP
            TMVAssert(BlasIsRM(A) || BlasIsCM(A));
            LapQRDecompose(A,beta,det);
#else
            NonLapQRDecompose(A,beta,det);
#endif
        }
#ifdef XDEBUG
        cout<<"Done QR_Decompose:\n";
        cout<<"Norm(A) = "<<Norm(A)<<endl;
        cout<<"beta = "<<beta<<endl;
        cout<<"Norm(beta) = "<<Norm(beta)<<endl;
        Matrix<T> R = A.upperTri();
        Matrix<T> Q(A);
        GetQFromQR(Q.view(),beta);
        Matrix<T> AA = Q*R;
        cout<<"Norm(AA-A0) = "<<Norm(AA-A0)<<endl;
        if (!(Norm(AA-A0) <= 0.0001*Norm(Q)*Norm(R))) {
            cerr<<"QRDecompose: A = "<<TMV_Text(A)<<"  "<<A0<<endl;
            cerr<<"-> "<<A<<endl;
            cerr<<"beta = "<<beta<<endl;
            cerr<<"Q*R = "<<Q*R<<endl;
            cerr<<"Norm(diff) = "<<Norm(AA-A0)<<endl;
            Matrix<T> A2 = A0;
            Vector<T> beta2(beta.size());
            T det2(0);
            NonBlockQRDecompose(A2.view(),beta2.view(),det2);
            cerr<<"NonBlock "<<A2<<endl;
            cerr<<"beta = "<<beta2<<endl;
            abort(); 
        }
#endif
    }


    //
    // QR Decompose - Unpacked
    //

    template <class T> 
    void QR_Decompose(
        MatrixView<T> Q, UpperTriMatrixView<T> R, T& det)
    {
        // Decompose A (input as Q) into A = Q R 
        // where Q is unitary and R is upper triangular

        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(R.colsize() == Q.rowsize());
        TMVAssert(R.rowsize() == Q.rowsize());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(R.ct() == NonConj);

        Vector<T> beta(Q.rowsize());
        QR_Decompose(Q,beta.view(),det);
        R = Q.upperTri();
        GetQFromQR(Q,beta.view());
    }

    template <class T> 
    void QR_Decompose(MatrixView<T> Q, UpperTriMatrixView<T> R)
    {
        T d(0);
        if (Q.isconj()) {
            if (R.isconj()) {
                QR_Decompose(Q.conjugate(),R.conjugate(),d);
            } else {
                QR_Decompose(Q.conjugate(),R,d);
                R.conjugateSelf();
            }
        } else {
            if (R.isconj()) {
                QR_Decompose(Q,R.conjugate(),d);
                R.conjugateSelf();
            } else {
                QR_Decompose(Q,R,d);
            }
        }
    }

    template <class T> 
    void QR_Decompose(MatrixView<T> A)
    {
        // Decompose A into Q R, but ignore Q.
        // On output, R is A.upperTri()

        TMVAssert(A.colsize() >= A.rowsize());

        Vector<T> beta(A.rowsize());
        T d(0);
        if (A.isconj()) 
            QR_Decompose(A.conjugate(),beta.view(),d);
        else
            QR_Decompose(A,beta.view(),d);
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_QRDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


