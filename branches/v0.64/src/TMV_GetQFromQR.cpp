///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2009                                                 //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file LICENSE.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


//#define XDEBUG

#include "TMV_Blas.h"
#include "TMV_QRDiv.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "TMV_Householder.h"

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
        const MatrixView<T>& Q, const GenVector<T>& beta)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(beta.ct() == NonConj);
        const int M = Q.colsize();
        const int N = Q.rowsize();
        Q.upperTri().setZero();
        const int sb = beta.step();
        if (sb == 1) {
            const T* bi = beta.cptr()+(N-1);
            for(int i=N-1;i>=0;--i,--bi) {
                HouseholderUnpack(Q.subMatrix(i,M,i,N),*bi);
            }
        } else {
            const T* bi = beta.cptr()+(N-1)*sb;
            for(int i=N-1;i>=0;--i,bi-=sb) {
                HouseholderUnpack(Q.subMatrix(i,M,i,N),*bi);
            }
        }
    }

    template <class T> 
    static void BlockGetQFromQR(
        const MatrixView<T>& Q, const GenVector<T>& beta)
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
        const int M = Q.colsize();
        const int N = Q.rowsize();
        Q.upperTri().setZero();
        UpperTriMatrix<T,NonUnitDiag,ColMajor> BaseZ(TMV_MIN(QR_BLOCKSIZE,N));
        for(int j2=N;j2>0;) {
            int j1 = j2 > QR_BLOCKSIZE ? j2-QR_BLOCKSIZE : 0;
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
        const MatrixView<T>& Q, const GenVector<T>& beta)
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
        const MatrixView<T>& Q, const GenVector<T>& beta)
    { NonLapGetQFromQR(Q,beta); }
#ifdef INST_DOUBLE
    template <> 
    void LapGetQFromQR(
        const MatrixView<double>& Q, const GenVector<double>& beta)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);
        TMVAssert(beta.step() == 1);
        int m = Q.colsize();
        int n = Q.rowsize();
        if (Q.isrm()) {
            int ldq = Q.stepi();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<double> work(lwork);
#else
            int lwork = -1;
            AlignedArray<double> work(1);
            LAPNAME(dorglq) (
                LAPCM LAPV(n),LAPV(m),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
                LAPP(beta.cptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(work[0]);
            work.resize(lwork);
#endif
#endif
            LAPNAME(dorglq) (
                LAPCM LAPV(n),LAPV(m),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
                LAPP(beta.cptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results("dorglq");
#else
            LAP_Results(int(work[0]),m,n,lwork,"dorglq");
#endif
        } else {
            int ldq = Q.stepj();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<double> work(lwork);
#else
            int lwork = -1;
            AlignedArray<double> work(1);
            LAPNAME(dorgqr) (
                LAPCM LAPV(m),LAPV(n),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
                LAPP(beta.cptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(work[0]);
            work.resize(lwork);
#endif
#endif
            LAPNAME(dorgqr) (
                LAPCM LAPV(m),LAPV(n),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
                LAPP(beta.cptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results("dorgqr");
#else
            LAP_Results(int(work[0]),m,n,lwork,"dorgqr");
#endif
        }
    }
    template <> 
    void LapGetQFromQR(
        const MatrixView<std::complex<double> >& Q,
        const GenVector<std::complex<double> >& beta)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);
        TMVAssert(beta.step() == 1);
        int m = Q.colsize();
        int n = Q.rowsize();
        if (Q.isrm()) {
            int ldq = Q.stepi();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<std::complex<double> > work(lwork);
#else
            int lwork = -1;
            AlignedArray<std::complex<double> > work(1);
            LAPNAME(zunglq) (
                LAPCM LAPV(n),LAPV(m),LAPV(n),
                LAPP(Q.ptr()),LAPV(ldq),LAPP(beta.cptr())
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(std::real(work[0]));
            work.resize(lwork);
#endif
#endif
            LAPNAME(zunglq) (
                LAPCM LAPV(n),LAPV(m),LAPV(n),
                LAPP(Q.ptr()),LAPV(ldq),LAPP(beta.cptr())
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results("zunglq");
#else
            LAP_Results(int(std::real(work[0])),m,n,lwork,"zunglq");
#endif
        } else {
            int ldq = Q.stepj();
            Vector<std::complex<double> > conjbeta = beta.conjugate();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<std::complex<double> > work(lwork);
#else
            int lwork = -1;
            AlignedArray<std::complex<double> > work(1);
            LAPNAME(zungqr) (
                LAPCM LAPV(m),LAPV(n),LAPV(n),
                LAPP(Q.ptr()),LAPV(ldq),LAPP(conjbeta.cptr())
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(std::real(work[0]));
            work.resize(lwork);
#endif
#endif
            LAPNAME(zungqr) (
                LAPCM LAPV(m),LAPV(n),LAPV(n),
                LAPP(Q.ptr()),LAPV(ldq),LAPP(conjbeta.cptr())
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results("zungqr");
#else
            LAP_Results(int(std::real(work[0])),m,n,lwork,"zungqr");
#endif
        }
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void LapGetQFromQR(
        const MatrixView<float>& Q, const GenVector<float>& beta)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);
        TMVAssert(beta.step() == 1);
        int m = Q.colsize();
        int n = Q.rowsize();
        if (Q.isrm()) {
            int ldq = Q.stepi();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<float> work(lwork);
#else
            int lwork = -1;
            AlignedArray<float> work(1);
            LAPNAME(sorglq) (
                LAPCM LAPV(n),LAPV(m),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
                LAPP(beta.cptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(work[0]);
            work.resize(lwork);
#endif
#endif
            LAPNAME(sorglq) (
                LAPCM LAPV(n),LAPV(m),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
                LAPP(beta.cptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results("sorglq");
#else
            LAP_Results(int(work[0]),m,n,lwork,"sorglq");
#endif
        } else {
            int ldq = Q.stepj();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<float> work(lwork);
#else
            int lwork = -1;
            AlignedArray<float> work(1);
            LAPNAME(sorgqr) (
                LAPCM LAPV(m),LAPV(n),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
                LAPP(beta.cptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(work[0]);
            work.resize(lwork);
#endif
#endif
            LAPNAME(sorgqr) (
                LAPCM LAPV(m),LAPV(n),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
                LAPP(beta.cptr()) LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results("sorgqr");
#else
            LAP_Results(int(work[0]),m,n,lwork,"sorgqr");
#endif
        }
    }
    template <> 
    void LapGetQFromQR(
        const MatrixView<std::complex<float> >& Q,
        const GenVector<std::complex<float> >& beta)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);
        TMVAssert(beta.step() == 1);
        int m = Q.colsize();
        int n = Q.rowsize();
        if (Q.isrm()) {
            int ldq = Q.stepi();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<std::complex<float> > work(lwork);
#else
            int lwork = -1;
            AlignedArray<std::complex<float> > work(1);
            LAPNAME(cunglq) (
                LAPCM LAPV(n),LAPV(m),LAPV(n),
                LAPP(Q.ptr()),LAPV(ldq),LAPP(beta.cptr())
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(std::real(work[0]));
            work.resize(lwork);
#endif
#endif
            LAPNAME(cunglq) (
                LAPCM LAPV(n),LAPV(m),LAPV(n),
                LAPP(Q.ptr()),LAPV(ldq),LAPP(beta.cptr())
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results("cunglq");
#else
            LAP_Results(int(std::real(work[0])),m,n,lwork,"cunglq");
#endif
        } else {
            int ldq = Q.stepj();
            Vector<std::complex<float> > conjbeta = beta.conjugate();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = n*LAP_BLOCKSIZE;
            AlignedArray<std::complex<float> > work(lwork);
#else
            int lwork = -1;
            AlignedArray<std::complex<float> > work(1);
            LAPNAME(cungqr) (
                LAPCM LAPV(m),LAPV(n),LAPV(n),
                LAPP(Q.ptr()),LAPV(ldq),LAPP(conjbeta.cptr())
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
            lwork = int(std::real(work[0]));
            work.resize(lwork);
#endif
#endif
            LAPNAME(cungqr) (
                LAPCM LAPV(m),LAPV(n),LAPV(n),
                LAPP(Q.ptr()),LAPV(ldq),LAPP(conjbeta.cptr())
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
            LAP_Results("cungqr");
#else
            LAP_Results(int(std::real(work[0])),m,n,lwork,"cungqr");
#endif
        }
    }
#endif
#endif

    template <class T> 
    void GetQFromQR(const MatrixView<T>& Q, const GenVector<T>& beta)
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

#define InstFile "TMV_GetQFromQR.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


