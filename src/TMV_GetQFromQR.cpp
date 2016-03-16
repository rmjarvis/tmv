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
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_Householder.h"

#ifdef XDEBUG
#include "tmv/TMV_MatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define QR_BLOCKSIZE TMV_BLOCKSIZE
#else
#define QR_BLOCKSIZE 64
#endif

    // 
    // Get Q from QR
    // 

    template <class T> 
    static void NonBlockGetQFromQR(
        MatrixView<T> Q, const GenVector<T>& beta)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(beta.ct() == NonConj);
        const ptrdiff_t M = Q.colsize();
        const ptrdiff_t N = Q.rowsize();
        Q.upperTri().setZero();
        const ptrdiff_t sb = beta.step();
        if (sb == 1) {
            const T* bi = beta.cptr()+(N-1);
            for(ptrdiff_t i=N-1;i>=0;--i,--bi) {
                HouseholderUnpack(Q.subMatrix(i,M,i,N),*bi);
            }
        } else {
            const T* bi = beta.cptr()+(N-1)*sb;
            for(ptrdiff_t i=N-1;i>=0;--i,bi-=sb) {
                HouseholderUnpack(Q.subMatrix(i,M,i,N),*bi);
            }
        }
    }

    template <class T> 
    static void BlockGetQFromQR(
        MatrixView<T> Q, const GenVector<T>& beta)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(beta.ct() == NonConj);
#ifdef XDEBUG
        Matrix<T> Q0(Q);
        Matrix<T> Q2(Q);
        NonBlockGetQFromQR(Q2.view(),beta);
        std::cout<<"Start Block GetQFromQR"<<std::endl;
#endif
        const ptrdiff_t M = Q.colsize();
        const ptrdiff_t N = Q.rowsize();
        Q.upperTri().setZero();
        UpperTriMatrix<T,NonUnitDiag|ColMajor> BaseZ(
            TMV_MIN(QR_BLOCKSIZE,int(N)));
        for(ptrdiff_t j2=N;j2>0;) {
            ptrdiff_t j1 = j2 > QR_BLOCKSIZE ? j2-QR_BLOCKSIZE : 0;
            MatrixView<T> Y = Q.subMatrix(j1,M,j1,j2);
            UpperTriMatrixView<T> Z = BaseZ.subTriMatrix(0,Y.rowsize());
            BlockHouseholderMakeZ(Y,Z,beta.subVector(j1,j2));
            BlockHouseholderUnpack(Y,Z,Q.subMatrix(j1,M,j2,N));
            j2 = j1;
        }
#ifdef XDEBUG
        if (Norm(Q-Q2) > 0.001*Norm(Q)*Norm(beta)) {
            cerr<<"BlockGetQ: Q = "<<TMV_Text(Q)<<"  "<<Q0<<endl;
            cerr<<"beta = "<<beta<<endl;
            cerr<<"-> Q = "<<Q<<endl;
            cerr<<"NonBlock Q = "<<Q2<<endl;
            abort(); 
        }
#endif
    }

    template <class T> 
    static inline void NonLapGetQFromQR(
        MatrixView<T> Q, const GenVector<T>& beta)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(beta.ct() == NonConj);
        if (Q.rowsize() >= QR_BLOCKSIZE)
            BlockGetQFromQR(Q,beta);
        else
            NonBlockGetQFromQR(Q,beta);
    }
#ifdef LAP
    template <class T> 
    static inline void LapGetQFromQR(
        MatrixView<T> Q, const GenVector<T>& beta)
    { NonLapGetQFromQR(Q,beta); }
#ifdef INST_DOUBLE
    template <> 
    void LapGetQFromQR(
        MatrixView<double> Q, const GenVector<double>& beta)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);
        TMVAssert(beta.step() == 1);
        int m = Q.colsize();
        int n = Q.rowsize();
        if (BlasIsCM(Q)) {
            int ldq = Q.stepj();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<double> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<double> work(1);
            work.get()[0] = 0.;
            LAPNAME(dorgqr) (
                LAPCM LAPV(m),LAPV(n),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
                LAPP(beta.cptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(work[0]);
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(dorgqr) (
                LAPCM LAPV(m),LAPV(n),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
                LAPP(beta.cptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"dorgqr");
#else
            LAP_Results(Lap_info,int(work[0]),m,n,lwork,"dorgqr");
#endif
        } else {
            int ldq = Q.stepi();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<double> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<double> work(1);
            work.get()[0] = 0.;
            LAPNAME(dorglq) (
                LAPCM LAPV(n),LAPV(m),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
                LAPP(beta.cptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(work[0]);
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(dorglq) (
                LAPCM LAPV(n),LAPV(m),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
                LAPP(beta.cptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"dorglq");
#else
            LAP_Results(Lap_info,int(work[0]),m,n,lwork,"dorglq");
#endif
        }
    }
    template <> 
    void LapGetQFromQR(
        MatrixView<std::complex<double> > Q,
        const GenVector<std::complex<double> >& beta)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);
        TMVAssert(beta.step() == 1);
        int m = Q.colsize();
        int n = Q.rowsize();
        if (BlasIsCM(Q)) {
            int ldq = Q.stepj();
            Vector<std::complex<double> > conjbeta = beta.conjugate();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<std::complex<double> > work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<std::complex<double> > work(1);
            work.get()[0] = 0.;
            LAPNAME(zungqr) (
                LAPCM LAPV(m),LAPV(n),LAPV(n),
                LAPP(Q.ptr()),LAPV(ldq),LAPP(conjbeta.cptr())
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(TMV_REAL(work[0]));
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(zungqr) (
                LAPCM LAPV(m),LAPV(n),LAPV(n),
                LAPP(Q.ptr()),LAPV(ldq),LAPP(conjbeta.cptr())
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"zungqr");
#else
            LAP_Results(Lap_info,int(TMV_REAL(work[0])),m,n,lwork,"zungqr");
#endif
        } else {
            int ldq = Q.stepi();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<std::complex<double> > work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<std::complex<double> > work(1);
            work.get()[0] = 0.;
            LAPNAME(zunglq) (
                LAPCM LAPV(n),LAPV(m),LAPV(n),
                LAPP(Q.ptr()),LAPV(ldq),LAPP(beta.cptr())
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(TMV_REAL(work[0]));
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(zunglq) (
                LAPCM LAPV(n),LAPV(m),LAPV(n),
                LAPP(Q.ptr()),LAPV(ldq),LAPP(beta.cptr())
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"zunglq");
#else
            LAP_Results(Lap_info,int(TMV_REAL(work[0])),m,n,lwork,"zunglq");
#endif
        }
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void LapGetQFromQR(
        MatrixView<float> Q, const GenVector<float>& beta)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);
        TMVAssert(beta.step() == 1);
        int m = Q.colsize();
        int n = Q.rowsize();
        if (BlasIsCM(Q)) {
            int ldq = Q.stepj();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<float> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<float> work(1);
            work.get()[0] = 0.F;
            LAPNAME(sorgqr) (
                LAPCM LAPV(m),LAPV(n),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
                LAPP(beta.cptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(work[0]);
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(sorgqr) (
                LAPCM LAPV(m),LAPV(n),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
                LAPP(beta.cptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"sorgqr");
#else
            LAP_Results(Lap_info,int(work[0]),m,n,lwork,"sorgqr");
#endif
        } else {
            int ldq = Q.stepi();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<float> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<float> work(1);
            work.get()[0] = 0.F;
            LAPNAME(sorglq) (
                LAPCM LAPV(n),LAPV(m),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
                LAPP(beta.cptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(work[0]);
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(sorglq) (
                LAPCM LAPV(n),LAPV(m),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
                LAPP(beta.cptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"sorglq");
#else
            LAP_Results(Lap_info,int(work[0]),m,n,lwork,"sorglq");
#endif
        }
    }
    template <> 
    void LapGetQFromQR(
        MatrixView<std::complex<float> > Q,
        const GenVector<std::complex<float> >& beta)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);
        TMVAssert(beta.step() == 1);
        int m = Q.colsize();
        int n = Q.rowsize();
        if (BlasIsCM(Q)) {
            int ldq = Q.stepj();
            Vector<std::complex<float> > conjbeta = beta.conjugate();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<std::complex<float> > work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<std::complex<float> > work(1);
            work.get()[0] = 0.F;
            LAPNAME(cungqr) (
                LAPCM LAPV(m),LAPV(n),LAPV(n),
                LAPP(Q.ptr()),LAPV(ldq),LAPP(conjbeta.cptr())
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(TMV_REAL(work[0]));
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(cungqr) (
                LAPCM LAPV(m),LAPV(n),LAPV(n),
                LAPP(Q.ptr()),LAPV(ldq),LAPP(conjbeta.cptr())
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"cungqr");
#else
            LAP_Results(Lap_info,int(TMV_REAL(work[0])),m,n,lwork,"cungqr");
#endif
        } else {
            int ldq = Q.stepi();
            int Lap_info=0;
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<std::complex<float> > work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#else
            int lwork = -1;
            AlignedArray<std::complex<float> > work(1);
            work.get()[0] = 0.F;
            LAPNAME(cunglq) (
                LAPCM LAPV(n),LAPV(m),LAPV(n),
                LAPP(Q.ptr()),LAPV(ldq),LAPP(beta.cptr())
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(TMV_REAL(work[0]));
            work.resize(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
#endif
            LAPNAME(cunglq) (
                LAPCM LAPV(n),LAPV(m),LAPV(n),
                LAPP(Q.ptr()),LAPV(ldq),LAPP(beta.cptr())
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results(Lap_info,"cunglq");
#else
            LAP_Results(Lap_info,int(TMV_REAL(work[0])),m,n,lwork,"cunglq");
#endif
        }
    }
#endif
#endif

    template <class T> 
    void GetQFromQR(MatrixView<T> Q, const GenVector<T>& beta)
    {
        // Extract the Q matrix from a packed Q matrix
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(beta.size() == Q.rowsize());
        TMVAssert(Q.isrm() || Q.iscm());
        TMVAssert(Q.ct() == NonConj);
        if (beta.isconj()) {
            Vector<T> b2 = beta;
            GetQFromQR(Q,b2.view());
        } else {
#ifdef XDEBUG
            std::cout<<"Start GetQFromQR: Q = "<<Q<<"beta = "<<beta<<std::endl;
            Matrix<T> Q0(Q);
            Matrix<T> Q2(Q);
            NonBlockGetQFromQR(Q2.view(),beta);
#endif
#ifdef LAP
            LapGetQFromQR(Q,beta);
#else
            NonLapGetQFromQR(Q,beta);
#endif
#ifdef XDEBUG
            std::cout<<"Done GetQFromQR: Q = "<<Q<<std::endl;
            std::cout<<"(QtQ-1).diag = "<<Matrix<T>(Q.adjoint()*Q-T(1)).diag()<<std::endl;
            std::cout<<"(Q2tQ2-1).diag = "<<Matrix<T>(Q2.adjoint()*Q2-T(1)).diag()<<std::endl;
            std::cout<<"Norm(QtQ-1) = "<<Norm(Q.adjoint()*Q-T(1))<<std::endl;
            std::cout<<"Norm(Q2tQ2-1) = "<<Norm(Q2.adjoint()*Q2-T(1))<<std::endl;
            if (Norm(Q-Q2) > 0.001*Norm(Q)*Norm(beta)) {
                cerr<<"LapGetQ: Q = "<<TMV_Text(Q)<<"  "<<Q0<<endl;
                cerr<<"beta = "<<beta<<endl;
                cerr<<"-> Q = "<<Q<<endl;
                cerr<<"NonBlock Q = "<<Q2<<endl;
                cerr<<"diff = "<<tmv::Matrix<T>(Q-Q2).clip(1.e-6)<<endl;
                cerr<<"Norm(diff) = "<<Norm(Q-Q2)<<endl;
                abort(); 
            }
#endif
        }
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_GetQFromQR.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


