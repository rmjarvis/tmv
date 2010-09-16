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
#include "tmv/TMV_QRD.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "TMV_Householder.h"
#include "tmv/TMV_PackedQ.h"

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
    // Packed Q - LDivEq
    //

    template <class T, class T1> 
    static void NonBlockQLDivEq(
        const GenMatrix<T1>& Q, const GenVector<T1>& beta,
        const MatrixView<T>& m)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.colsize() == m.colsize());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);
        // Solve Q x = m in place 
        // where Q is stored as Householder vectors along with beta
        //
        // Q is H0t H1t ... H_N-1t
        // So x = H_N-1 .. H1 H0 m
        //
        const int M = Q.colsize();
        const int N = Q.rowsize();
        for(int j=0;j<N;++j) if (beta(j) != T1(0)) {
            HouseholderLMult(Q.col(j,j+1,M),beta(j),m.rowRange(j,M));
        }
    }

    template <class T, class T1> 
    static void BlockQLDivEq(
        const GenMatrix<T1>& Q, const GenVector<T1>& beta,
        const MatrixView<T>& m)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.colsize() == m.colsize());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);

#ifdef XDEBUG
        Matrix<T> m0(m);
        Matrix<T> m2(m);
        NonBlockQLDivEq(Q,beta,m2.view());
#endif

        // x = H_N-1 .. H1 H0 m
        // In first block step:
        // m = (Hr .. H0) m
        //   = (H0t .. Hrt)t m
        // So form Y,Z from Ht's, rather than H's, and then call LDiv

        const int M = Q.colsize();
        const int N = Q.rowsize();
        UpperTriMatrix<T1,NonUnitDiag,ColMajor> BaseZ(TMV_MIN(QR_BLOCKSIZE,N));
        for(int j1=0;j1<N;) {
            int j2 = TMV_MIN(N,j1+QR_BLOCKSIZE);
            ConstMatrixView<T1> Y = Q.subMatrix(j1,M,j1,j2);
            UpperTriMatrixView<T1> Z = BaseZ.subTriMatrix(0,Y.rowsize());
            BlockHouseholderMakeZ(Y,Z,beta.subVector(j1,j2));
            BlockHouseholderLDiv(Y,Z,m.rowRange(j1,M));
            j1 = j2;
        }
#ifdef XDEBUG
        if (Norm(m-m2) > 0.001*Norm(Q)*Norm(m0)) {
            cerr<<"BlockQLDivEq: Q = "<<TMV_Text(Q)<<"  "<<Q<<endl;
            cerr<<"beta = "<<beta<<endl;
            cerr<<"m = "<<TMV_Text(m)<<"  "<<m0<<endl;
            cerr<<"-> m = "<<m<<endl;
            cerr<<"NonBlock m = "<<m2<<endl;
            abort(); 
        }
#endif
    }

    template <class T, class T1> 
    static void NonLapQLDivEq(
        const GenMatrix<T1>& Q, const GenVector<T1>& beta,
        const MatrixView<T>& m)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.colsize() == m.colsize());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);

        if (Q.rowsize() > QR_BLOCKSIZE && m.rowsize() > QR_BLOCKSIZE)
            BlockQLDivEq(Q,beta,m);
        else
            NonBlockQLDivEq(Q,beta,m);
    }
#ifdef LAP
    template <class T, class T1> 
    static inline void LapQLDivEq(
        const GenMatrix<T1>& Q, const GenVector<T1>& beta,
        const MatrixView<T>& m)
    { NonLapQLDivEq(Q,beta,m); }
#ifdef INST_DOUBLE
    template <> 
    void LapQLDivEq(
        const GenMatrix<double>& Q,
        const GenVector<double>& beta, const MatrixView<double>& x)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.colsize() == x.colsize());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);
        int ldx = x.isrm() ? x.stepi() : x.stepj();
        int m = x.isrm() ? x.rowsize() : x.colsize();
        int n = x.isrm() ? x.colsize() : x.rowsize();
        int k = Q.rowsize();
        if (Q.isrm()) {
            int ldq = Q.stepi();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.rowsize()*LAP_BLOCKSIZE;
            AlignedArray<double> work(lwork);
#else
            int lwork = -1;
            AlignedArray<double> work(1);
            LAPNAME(dormlq) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_T:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(work[0]);
            work.resize(lwork);
#endif
#endif
            LAPNAME(dormlq) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_T:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
#ifdef LAPNOWORK
            LAP_Results("dormlq");
#else
            LAP_Results(int(work[0]),m,n,lwork,"dormlq");
#endif
        } else {
            int ldq = Q.stepj();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.rowsize()*LAP_BLOCKSIZE;
            AlignedArray<double> work(lwork);
#else
            int lwork = -1;
            AlignedArray<double> work(1);
            LAPNAME(dormqr) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L, 
                x.isrm()?LAPCH_NT:LAPCH_T, LAPV(m),LAPV(n),LAPV(k),
                LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(work[0]);
            work.resize(lwork);
#endif
#endif
            LAPNAME(dormqr) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L, 
                x.isrm()?LAPCH_NT:LAPCH_T, LAPV(m),LAPV(n),LAPV(k),
                LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
#ifdef LAPNOWORK
            LAP_Results("dormqr");
#else
            LAP_Results(int(work[0]),m,n,lwork,"dormqr");
#endif
        }
    }
    template <> 
    void LapQLDivEq(
        const GenMatrix<std::complex<double> >& Q,
        const GenVector<std::complex<double> >& beta,
        const MatrixView<std::complex<double> >& x)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.colsize() == x.colsize());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);
        int k = Q.rowsize();
        if (Q.isrm()) {
            int ldx = x.isrm() ? x.stepi() : x.stepj();
            int m = x.isrm() ? x.rowsize() : x.colsize();
            int n = x.isrm() ? x.colsize() : x.rowsize();
            int ldq = Q.stepi();
            if (x.isrm() == x.isconj()) x.conjugateSelf();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.rowsize()*LAP_BLOCKSIZE;
            AlignedArray<std::complex<double> > work(lwork);
#else
            int lwork = -1;
            AlignedArray<std::complex<double> > work(1);
            LAPNAME(zunmlq) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_CT:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(std::real(work[0]));
            work.resize(lwork);
#endif
#endif
            LAPNAME(zunmlq) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_CT:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            if (x.isrm() == x.isconj()) x.conjugateSelf();
#ifdef LAPNOWORK
            LAP_Results("zunmlq");
#else
            LAP_Results(int(TMV_REAL(work[0])),m,n,lwork,"zunmlq");
#endif
        } else {
            int ldx = x.isrm() ? x.stepi() : x.stepj();
            int m = x.isrm() ? x.rowsize() : x.colsize();
            int n = x.isrm() ? x.colsize() : x.rowsize();
            int ldq = Q.stepj();
            if (x.isrm() != x.isconj()) x.conjugateSelf();
            Vector<std::complex<double> > conjbeta = beta.conjugate();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.rowsize()*LAP_BLOCKSIZE;
            AlignedArray<std::complex<double> > work(lwork);
#else
            int lwork = -1;
            AlignedArray<std::complex<double> > work(1);
            LAPNAME(zunmqr) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_NT:LAPCH_CT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(conjbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(std::real(work[0]));
            work.resize(lwork);
#endif
#endif
            LAPNAME(zunmqr) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_NT:LAPCH_CT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(conjbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            if (x.isrm() != x.isconj()) x.conjugateSelf();
#ifdef LAPNOWORK
            LAP_Results("zunmqr");
#else
            LAP_Results(int(TMV_REAL(work[0])),m,n,lwork,"zunmqr");
#endif
        }
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void LapQLDivEq(
        const GenMatrix<float>& Q,
        const GenVector<float>& beta, const MatrixView<float>& x)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.colsize() == x.colsize());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);
        int ldx = x.isrm() ? x.stepi() : x.stepj();
        int m = x.isrm() ? x.rowsize() : x.colsize();
        int n = x.isrm() ? x.colsize() : x.rowsize();
        int k = Q.rowsize();
        if (Q.isrm()) {
            int ldq = Q.stepi();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.rowsize()*LAP_BLOCKSIZE;
            AlignedArray<float> work(lwork);
#else
            int lwork = -1;
            AlignedArray<float> work(1);
            LAPNAME(sormlq) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_T:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(work[0]);
            work.resize(lwork);
#endif
#endif
            LAPNAME(sormlq) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_T:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
#ifdef LAPNOWORK
            LAP_Results("sormlq");
#else
            LAP_Results(int(work[0]),m,n,lwork,"sormlq");
#endif
        } else {
            int ldq = Q.stepj();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.rowsize()*LAP_BLOCKSIZE;
            AlignedArray<float> work(lwork);
#else
            int lwork = -1;
            AlignedArray<float> work(1);
            LAPNAME(sormqr) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L, 
                x.isrm()?LAPCH_NT:LAPCH_T, LAPV(m),LAPV(n),LAPV(k),
                LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(work[0]);
            work.resize(lwork);
#endif
#endif
            LAPNAME(sormqr) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L, 
                x.isrm()?LAPCH_NT:LAPCH_T, LAPV(m),LAPV(n),LAPV(k),
                LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
#ifdef LAPNOWORK
            LAP_Results("sormqr");
#else
            LAP_Results(int(work[0]),m,n,lwork,"sormqr");
#endif
        }
    }
    template <> 
    void LapQLDivEq(
        const GenMatrix<std::complex<float> >& Q,
        const GenVector<std::complex<float> >& beta,
        const MatrixView<std::complex<float> >& x)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.colsize() == x.colsize());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);
        int k = Q.rowsize();
        if (Q.isrm()) {
            int ldx = x.isrm() ? x.stepi() : x.stepj();
            int m = x.isrm() ? x.rowsize() : x.colsize();
            int n = x.isrm() ? x.colsize() : x.rowsize();
            int ldq = Q.stepi();
            if (x.isrm() == x.isconj()) x.conjugateSelf();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.rowsize()*LAP_BLOCKSIZE;
            AlignedArray<std::complex<float> > work(lwork);
#else
            int lwork = -1;
            AlignedArray<std::complex<float> > work(1);
            LAPNAME(cunmlq) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_CT:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(std::real(work[0]));
            work.resize(lwork);
#endif
#endif
            LAPNAME(cunmlq) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_CT:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            if (x.isrm() == x.isconj()) x.conjugateSelf();
#ifdef LAPNOWORK
            LAP_Results("cunmlq");
#else
            LAP_Results(int(TMV_REAL(work[0])),m,n,lwork,"cunmlq");
#endif
        } else {
            int ldx = x.isrm() ? x.stepi() : x.stepj();
            int m = x.isrm() ? x.rowsize() : x.colsize();
            int n = x.isrm() ? x.colsize() : x.rowsize();
            int ldq = Q.stepj();
            if (x.isrm() != x.isconj()) x.conjugateSelf();
            Vector<std::complex<float> > conjbeta = beta.conjugate();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.rowsize()*LAP_BLOCKSIZE;
            AlignedArray<std::complex<float> > work(lwork);
#else
            int lwork = -1;
            AlignedArray<std::complex<float> > work(1);
            LAPNAME(cunmqr) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_NT:LAPCH_CT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(conjbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(std::real(work[0]));
            work.resize(lwork);
#endif
#endif
            LAPNAME(cunmqr) (
                LAPCM x.isrm()?LAPCH_R:LAPCH_L,
                x.isrm()?LAPCH_NT:LAPCH_CT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(conjbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            if (x.isrm() != x.isconj()) x.conjugateSelf();
#ifdef LAPNOWORK
            LAP_Results("cunmqr");
#else
            LAP_Results(int(TMV_REAL(work[0])),m,n,lwork,"cunmqr");
#endif
        }
    }
#endif
#endif

    template <class T, class T1> 
    void Q_LDivEq(
        const GenMatrix<T1>& Q, const GenVector<T1>& beta,
        const MatrixView<T>& m)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(beta.size() == Q.rowsize());
        TMVAssert(m.colsize() == Q.colsize());
        TMVAssert(Q.isrm() || Q.iscm());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);

#ifdef XDEBUG
        Matrix<T1> QQ(Q.colsize(),Q.colsize(),0.);
        QQ.setToIdentity();
        QQ.colRange(0,Q.rowsize()) = Q;
        Vector<T1> bb(Q.colsize(),0.);
        bb.subVector(0,beta.size()) = beta;
        GetQFromQR(QQ.view(),bb);
        Matrix<T1> Qx = Q;
        Vector<T1> bx = beta;
        GetQFromQR(Qx.view(),bx);
        Matrix<T> m0 = m;
        Matrix<T> m2 = QQ.adjoint() * m;
        std::cout<<"Start Q_LDivEq: \n";
        std::cout<<"Q = "<<TMV_Text(Q)<<std::endl;
        std::cout<<"beta = "<<TMV_Text(beta)<<std::endl;
        std::cout<<"m = "<<TMV_Text(m)<<std::endl;
        std::cout<<"Q = "<<Q<<std::endl;
        std::cout<<"beta = "<<beta<<std::endl;
        std::cout<<"QQ = "<<QQ<<std::endl;
        std::cout<<"Qx = "<<Qx<<std::endl;
        std::cout<<"Norm(Qx-QQx) = "<<Norm(Qx-QQ.colRange(0,Q.rowsize()))<<std::endl;
        std::cout<<"m0 = "<<m<<std::endl;
        std::cout<<"m2 = "<<m2<<std::endl;
        std::cout<<"m2' = "<<(Qx.adjoint() * m)<<std::endl;
#ifdef LAP
        Matrix<T> m3 = m;
        NonLapQLDivEq(Q,beta,m3.view());
        std::cout<<"m3 = "<<m3<<std::endl;
#endif
#endif

        if (m.colsize() > 0 && m.rowsize() > 0) {
#ifdef LAP
            if ( m.isrm() || m.iscm() )
                LapQLDivEq(Q,beta,m);
            else
#endif
                NonLapQLDivEq(Q,beta,m);
        }

#ifdef XDEBUG
        std::cout<<"m = "<<m<<std::endl;
        std::cout<<"m2 = "<<m2<<std::endl;
        std::cout<<"m-m2 = "<<m-m2<<std::endl;
        std::cout<<"Norm(m-m2) = "<<Norm(m-m2)<<std::endl;
        if (Norm(m-m2) > 0.001*Norm(m0)) {
            cerr<<"Q_LDivEq\n";
            cerr<<"Q = "<<Q<<endl;
            cerr<<"beta = "<<beta<<endl;
            cerr<<"m = "<<TMV_Text(m)<<"  "<<m0<<endl;
            cerr<<"m => "<<m<<endl;
            cerr<<"m2 = "<<m2<<endl;
            cerr<<"Norm(m-m2) = "<<Norm(m-m2)<<endl;
            abort();
        }
#endif
    }

    //
    // Packed Q - RDivEq
    //

    template <class T, class T1> 
    static void NonBlockQRDivEq(
        const GenMatrix<T1>& Q, const GenVector<T1>& beta,
        const MatrixView<T>& m)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(beta.size() == Q.rowsize());
        TMVAssert(m.rowsize() == Q.colsize());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);
        // Solve x Q = m in place
        // where Q is stored as Householder vectors along with beta
        //
        // x = m Qt 
        // Qt is H_N-1 H_N-2 ... H1 H0
        const int M = Q.colsize();
        const int N = Q.rowsize();
        for(int j=N-1;j>=0;--j) if (beta(j) != T1(0)) {
            HouseholderLMult(
                Q.col(j,j+1,M).conjugate(),beta(j),
                m.colRange(j,M).transpose());
        }
    }

    template <class T, class T1> 
    static void BlockQRDivEq(
        const GenMatrix<T1>& Q, const GenVector<T1>& beta,
        const MatrixView<T>& m)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.colsize() == m.rowsize());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);

#ifdef XDEBUG
        Matrix<T> m0(m);
        Matrix<T> m2(m);
        NonBlockQRDivEq(Q,beta,m2.view());
#endif

        // x = m Qt 
        // x = m H_N-1 H_N-2 ... H1 H0
        // Again form Y,Z from Ht's, rather than H's, and then call RDiv

        const int M = Q.colsize();
        const int N = Q.rowsize();
        UpperTriMatrix<T1,NonUnitDiag,ColMajor> BaseZ(TMV_MIN(QR_BLOCKSIZE,N));
        for(int j2=N;j2>0;) {
            int j1 = j2 > QR_BLOCKSIZE ? j2-QR_BLOCKSIZE : 0;
            ConstMatrixView<T1> Y = Q.subMatrix(j1,M,j1,j2);
            UpperTriMatrixView<T1> Z = BaseZ.subTriMatrix(0,Y.rowsize());
            BlockHouseholderMakeZ(Y,Z,beta.subVector(j1,j2));
            BlockHouseholderLMult(Y,Z,m.colRange(j1,M).adjoint());
            j2 = j1;
        }
#ifdef XDEBUG
        if (Norm(m-m2) > 0.001*Norm(Q)*Norm(m0)) {
            cerr<<"BlockQRDivEq: Q = "<<TMV_Text(Q)<<"  "<<Q<<endl;
            cerr<<"beta = "<<beta<<endl;
            cerr<<"m = "<<TMV_Text(m)<<"  "<<m0<<endl;
            cerr<<"-> m = "<<m<<endl;
            cerr<<"NonBlock m = "<<m2<<endl;
            abort(); 
        }
#endif
    }

    template <class T, class T1> 
    static void NonLapQRDivEq(
        const GenMatrix<T1>& Q, const GenVector<T1>& beta,
        const MatrixView<T>& m)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.colsize() == m.rowsize());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);

        if (Q.rowsize() > QR_BLOCKSIZE && m.colsize() > QR_BLOCKSIZE)
            BlockQRDivEq(Q,beta,m);
        else
            NonBlockQRDivEq(Q,beta,m);
    }

#ifdef LAP
    template <class T, class T1> 
    static inline void LapQRDivEq(
        const GenMatrix<T1>& Q, const GenVector<T1>& beta,
        const MatrixView<T>& m)
    { NonLapQRDivEq(Q,beta,m); }
#ifdef INST_DOUBLE
    template <> 
    void LapQRDivEq(
        const GenMatrix<double>& Q,
        const GenVector<double>& beta, const MatrixView<double>& x)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.colsize() == x.rowsize());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);

        int ldx = x.isrm() ? x.stepi() : x.stepj();
        int m = x.isrm() ? x.rowsize() : x.colsize();
        int n = x.isrm() ? x.colsize() : x.rowsize();
        int k = Q.rowsize();
        if (Q.isrm()) {
            int ldq = Q.stepi();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.colsize()*LAP_BLOCKSIZE;
            AlignedArray<double> work(lwork);
#else
            int lwork = -1;
            AlignedArray<double> work(1);
            LAPNAME(dormlq) (
                LAPCM x.isrm()?LAPCH_L:LAPCH_R,
                x.isrm()?LAPCH_T:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(work[0]);
            work.resize(lwork);
#endif
#endif
            LAPNAME(dormlq) (
                LAPCM x.isrm()?LAPCH_L:LAPCH_R,
                x.isrm()?LAPCH_T:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
#ifdef LAPNOWORK
            LAP_Results("dormlq");
#else
            LAP_Results(int(work[0]),m,n,lwork,"dormlq");
#endif
        } else {
            int ldq = Q.stepj();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.colsize()*LAP_BLOCKSIZE;
            AlignedArray<double> work(lwork);
#else
            int lwork = -1;
            AlignedArray<double> work(1);
            LAPNAME(dormqr) (
                LAPCM x.isrm()?LAPCH_L:LAPCH_R, 
                x.isrm()?LAPCH_NT:LAPCH_T, LAPV(m),LAPV(n),LAPV(k),
                LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(work[0]);
            work.resize(lwork);
#endif
#endif
            LAPNAME(dormqr) (
                LAPCM x.isrm()?LAPCH_L:LAPCH_R, 
                x.isrm()?LAPCH_NT:LAPCH_T, LAPV(m),LAPV(n),LAPV(k),
                LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
#ifdef LAPNOWORK
            LAP_Results("dormqr");
#else
            LAP_Results(int(work[0]),m,n,lwork,"dormqr");
#endif
        }
    }
    template <> 
    void LapQRDivEq(
        const GenMatrix<std::complex<double> >& Q,
        const GenVector<std::complex<double> >& beta,
        const MatrixView<std::complex<double> >& x)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.colsize() == x.rowsize());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);
        int k = Q.rowsize();
        if (Q.isrm()) {
            int ldx = x.isrm() ? x.stepi() : x.stepj();
            int m = x.isrm() ? x.rowsize() : x.colsize();
            int n = x.isrm() ? x.colsize() : x.rowsize();
            int ldq = Q.stepi();
            if (x.isrm() == x.isconj()) x.conjugateSelf();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.colsize()*LAP_BLOCKSIZE;
            AlignedArray<std::complex<double> > work(lwork);
#else
            int lwork = -1;
            AlignedArray<std::complex<double> > work(1);
            LAPNAME(zunmlq) (
                LAPCM x.isrm()?LAPCH_L:LAPCH_R,
                x.isrm()?LAPCH_CT:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(std::real(work[0]));
            work.resize(lwork);
#endif
#endif
            LAPNAME(zunmlq) (
                LAPCM x.isrm()?LAPCH_L:LAPCH_R,
                x.isrm()?LAPCH_CT:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            if (x.isrm() == x.isconj()) x.conjugateSelf();
#ifdef LAPNOWORK
            LAP_Results("zunmlq");
#else
            LAP_Results(int(TMV_REAL(work[0])),m,n,lwork,"zunmlq");
#endif
        } else {
            int ldx = x.isrm() ? x.stepi() : x.stepj();
            int m = x.isrm() ? x.rowsize() : x.colsize();
            int n = x.isrm() ? x.colsize() : x.rowsize();
            int ldq = Q.stepj();
            if (x.isrm() != x.isconj()) x.conjugateSelf();
            Vector<std::complex<double> > conjbeta = beta.conjugate();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.colsize()*LAP_BLOCKSIZE;
            AlignedArray<std::complex<double> > work(lwork);
#else
            int lwork = -1;
            AlignedArray<std::complex<double> > work(1);
            LAPNAME(zunmqr) (
                LAPCM x.isrm()?LAPCH_L:LAPCH_R,
                x.isrm()?LAPCH_NT:LAPCH_CT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(conjbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(std::real(work[0]));
            work.resize(lwork);
#endif
#endif
            LAPNAME(zunmqr) (
                LAPCM x.isrm()?LAPCH_L:LAPCH_R,
                x.isrm()?LAPCH_NT:LAPCH_CT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(conjbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            if (x.isrm() != x.isconj()) x.conjugateSelf();
#ifdef LAPNOWORK
            LAP_Results("zunmqr");
#else
            LAP_Results(int(TMV_REAL(work[0])),m,n,lwork,"zunmqr");
#endif
        }
    }
#endif
#ifdef INST_FLOAT
    template <> 
    void LapQRDivEq(
        const GenMatrix<float>& Q,
        const GenVector<float>& beta, const MatrixView<float>& x)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.colsize() == x.rowsize());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);

        int ldx = x.isrm() ? x.stepi() : x.stepj();
        int m = x.isrm() ? x.rowsize() : x.colsize();
        int n = x.isrm() ? x.colsize() : x.rowsize();
        int k = Q.rowsize();
        if (Q.isrm()) {
            int ldq = Q.stepi();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.colsize()*LAP_BLOCKSIZE;
            AlignedArray<float> work(lwork);
#else
            int lwork = -1;
            AlignedArray<float> work(1);
            LAPNAME(sormlq) (
                LAPCM x.isrm()?LAPCH_L:LAPCH_R,
                x.isrm()?LAPCH_T:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(work[0]);
            work.resize(lwork);
#endif
#endif
            LAPNAME(sormlq) (
                LAPCM x.isrm()?LAPCH_L:LAPCH_R,
                x.isrm()?LAPCH_T:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
#ifdef LAPNOWORK
            LAP_Results("sormlq");
#else
            LAP_Results(int(work[0]),m,n,lwork,"sormlq");
#endif
        } else {
            int ldq = Q.stepj();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.colsize()*LAP_BLOCKSIZE;
            AlignedArray<float> work(lwork);
#else
            int lwork = -1;
            AlignedArray<float> work(1);
            LAPNAME(sormqr) (
                LAPCM x.isrm()?LAPCH_L:LAPCH_R, 
                x.isrm()?LAPCH_NT:LAPCH_T, LAPV(m),LAPV(n),LAPV(k),
                LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(work[0]);
            work.resize(lwork);
#endif
#endif
            LAPNAME(sormqr) (
                LAPCM x.isrm()?LAPCH_L:LAPCH_R, 
                x.isrm()?LAPCH_NT:LAPCH_T, LAPV(m),LAPV(n),LAPV(k),
                LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
#ifdef LAPNOWORK
            LAP_Results("sormqr");
#else
            LAP_Results(int(work[0]),m,n,lwork,"sormqr");
#endif
        }
    }
    template <> 
    void LapQRDivEq(
        const GenMatrix<std::complex<float> >& Q,
        const GenVector<std::complex<float> >& beta,
        const MatrixView<std::complex<float> >& x)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(Q.rowsize() == beta.size());
        TMVAssert(Q.colsize() == x.rowsize());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);
        int k = Q.rowsize();
        if (Q.isrm()) {
            int ldx = x.isrm() ? x.stepi() : x.stepj();
            int m = x.isrm() ? x.rowsize() : x.colsize();
            int n = x.isrm() ? x.colsize() : x.rowsize();
            int ldq = Q.stepi();
            if (x.isrm() == x.isconj()) x.conjugateSelf();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.colsize()*LAP_BLOCKSIZE;
            AlignedArray<std::complex<float> > work(lwork);
#else
            int lwork = -1;
            AlignedArray<std::complex<float> > work(1);
            LAPNAME(cunmlq) (
                LAPCM x.isrm()?LAPCH_L:LAPCH_R,
                x.isrm()?LAPCH_CT:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(std::real(work[0]));
            work.resize(lwork);
#endif
#endif
            LAPNAME(cunmlq) (
                LAPCM x.isrm()?LAPCH_L:LAPCH_R,
                x.isrm()?LAPCH_CT:LAPCH_NT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            if (x.isrm() == x.isconj()) x.conjugateSelf();
#ifdef LAPNOWORK
            LAP_Results("cunmlq");
#else
            LAP_Results(int(TMV_REAL(work[0])),m,n,lwork,"cunmlq");
#endif
        } else {
            int ldx = x.isrm() ? x.stepi() : x.stepj();
            int m = x.isrm() ? x.rowsize() : x.colsize();
            int n = x.isrm() ? x.colsize() : x.rowsize();
            int ldq = Q.stepj();
            if (x.isrm() != x.isconj()) x.conjugateSelf();
            Vector<std::complex<float> > conjbeta = beta.conjugate();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
            int lwork = x.colsize()*LAP_BLOCKSIZE;
            AlignedArray<std::complex<float> > work(lwork);
#else
            int lwork = -1;
            AlignedArray<std::complex<float> > work(1);
            LAPNAME(cunmqr) (
                LAPCM x.isrm()?LAPCH_L:LAPCH_R,
                x.isrm()?LAPCH_NT:LAPCH_CT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(conjbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            lwork = int(std::real(work[0]));
            work.resize(lwork);
#endif
#endif
            LAPNAME(cunmqr) (
                LAPCM x.isrm()?LAPCH_L:LAPCH_R,
                x.isrm()?LAPCH_NT:LAPCH_CT,
                LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
                LAPP(conjbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
                LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
            if (x.isrm() != x.isconj()) x.conjugateSelf();
#ifdef LAPNOWORK
            LAP_Results("cunmqr");
#else
            LAP_Results(int(TMV_REAL(work[0])),m,n,lwork,"cunmqr");
#endif
        }
    }
#endif // FLOAT
#endif // LAP

    template <class T, class T1> 
    void Q_RDivEq(
        const GenMatrix<T1>& Q, const GenVector<T1>& beta,
        const MatrixView<T>& m)
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(beta.size() == Q.rowsize());
        TMVAssert(m.rowsize() == Q.colsize());
        TMVAssert(Q.isrm() || Q.iscm());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);

#ifdef XDEBUG
        Matrix<T1> QQ(Q.colsize(),Q.colsize(),0.);
        QQ.setToIdentity();
        QQ.colRange(0,Q.rowsize()) = Q;
        Vector<T1> bb(Q.colsize(),0.);
        bb.subVector(0,beta.size()) = beta;
        GetQFromQR(QQ.view(),bb);
        Matrix<T> m0 = m;
        Matrix<T> m2 = m * QQ.adjoint();
        std::cout<<"Start Q_RDivEq: \n";
        std::cout<<"Q = "<<TMV_Text(Q)<<std::endl;
        std::cout<<"beta = "<<TMV_Text(beta)<<std::endl;
        std::cout<<"m = "<<TMV_Text(m)<<std::endl;
        std::cout<<"Q = "<<Q<<std::endl;
        std::cout<<"beta = "<<beta<<std::endl;
        std::cout<<"m = "<<m<<std::endl;
#endif

        if (m.colsize() > 0 && m.rowsize() > 0) {
#ifdef LAP
            if ( m.isrm() || m.iscm() )
                LapQRDivEq(Q,beta,m);
            else
#endif
                NonLapQRDivEq(Q,beta,m);
        }

#ifdef XDEBUG
        std::cout<<"Norm(m-m2) = "<<Norm(m-m2)<<std::endl;
        if (Norm(m-m2) > 0.001*Norm(m0)) {
            cerr<<"Q_RDivEq\n";
            cerr<<"Q = "<<Q<<endl;
            cerr<<"beta = "<<beta<<endl;
            cerr<<"m = "<<TMV_Text(m)<<"  "<<m0<<endl;
            cerr<<"m => "<<m<<endl;
            cerr<<"m2 = "<<m2<<endl;
            cerr<<"Norm(m-m2) = "<<Norm(m-m2)<<endl;
            abort();
        }
#endif
    }

    //
    // PackedQ routines
    //

    template <class T> 
    static void unpack(
        const GenMatrix<T>& Q, const GenVector<T>& beta,
        const MatrixView<T>& m)
    {
        if (m.isrm() || m.iscm()) {
            m = Q;
            GetQFromQR(m,beta);
        } else {
            Matrix<T> m1 = Q;
            GetQFromQR(m1.view(),beta);
            m = m1;
        }
    }

    template <class T, class T1> 
    static void unpack(
        const GenMatrix<T>& Q, const GenVector<T>& beta,
        const MatrixView<T1>& m)
    {
        Matrix<T> m1 = Q;
        GetQFromQR(m1.view(),beta);
        m = m1;
    }

    template <class T> template <class T1> 
    void PackedQ<T>::doAssignToM(const MatrixView<T1>& m) const
    { unpack(Q,beta,m); }

    template <class T> template <class T1> 
    void PackedQ<T>::LDivEq(const VectorView<T1>& v) const
    {
        TMVAssert(Q.colsize() == Q.rowsize());
        TMVAssert(beta.size() == Q.rowsize());
        TMVAssert(v.size() == Q.colsize());
        TMVAssert(Q.isrm() || Q.iscm());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);

        Q_LDivEq(Q,beta,ColVectorViewOf(v));
    }

    template <class T> template <class T1> 
    void PackedQ<T>::RDivEq(const VectorView<T1>& v) const
    {
        TMVAssert(Q.colsize() == Q.rowsize());
        TMVAssert(beta.size() == Q.rowsize());
        TMVAssert(v.size() == Q.colsize());
        TMVAssert(Q.isrm() || Q.iscm());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);

        Q_RDivEq(Q,beta,RowVectorViewOf(v));
    }

    template <class T> template <class T1, class T2> 
    void PackedQ<T>::LDiv(
        const GenVector<T1>& v, const VectorView<T2>& x) const
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(beta.size() == Q.rowsize());
        TMVAssert(v.size() == Q.colsize());
        TMVAssert(x.size() == Q.rowsize());
        TMVAssert(Q.isrm() || Q.iscm());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);

        // Solve Q x = v
        if (Q.isSquare()) {
            x = v;
            Q_LDivEq(Q,beta,ColVectorViewOf(x));
        } else {
            Vector<T1> v1 = v;
            Q_LDivEq(Q,beta,ColVectorViewOf(v1));
            x = v1.subVector(0,x.size());
        }
    }

    template <class T> template <class T1, class T2> 
    void PackedQ<T>::RDiv(
        const GenVector<T1>& v, const VectorView<T2>& x) const
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(beta.size() == Q.rowsize());
        TMVAssert(x.size() == Q.colsize());
        TMVAssert(v.size() == Q.rowsize());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);

        x.subVector(0,v.size()) = v;
        x.subVector(v.size(),x.size()).setZero();
        Q_RDivEq(Q,beta,RowVectorViewOf(x));
    }

    template <class T> template <class T1> 
    void PackedQ<T>::LDivEq(const MatrixView<T1>& m) const
    {
        TMVAssert(Q.colsize() == Q.rowsize());
        TMVAssert(beta.size() == Q.rowsize());
        TMVAssert(m.colsize() == Q.colsize());
        TMVAssert(Q.isrm() || Q.iscm());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);

        Q_LDivEq(Q,beta,m);
    }

    template <class T> template <class T1> 
    void PackedQ<T>::RDivEq(const MatrixView<T1>& m) const
    {
        TMVAssert(Q.colsize() == Q.rowsize());
        TMVAssert(beta.size() == Q.rowsize());
        TMVAssert(m.rowsize() == Q.rowsize());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);

        Q_RDivEq(Q,beta,m);
    }

    template <class T> template <class T1, class T2> 
    void PackedQ<T>::LDiv(
        const GenMatrix<T1>& m, const MatrixView<T2>& x) const
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(beta.size() == Q.rowsize());
        TMVAssert(m.colsize() == Q.colsize());
        TMVAssert(x.colsize() == Q.rowsize());
        TMVAssert(x.rowsize() == m.rowsize());
        TMVAssert(Q.isrm() || Q.iscm());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);

        // Solve Q x = m
        if (Q.isSquare()) {
            x = m;
            Q_LDivEq(Q,beta,x);
        } else {
            Matrix<T,ColMajor> m1 = m;
            Q_LDivEq(Q,beta,m1.view());
            x = m1.rowRange(0,x.colsize());
        }
    }

    template <class T> template <class T1, class T2> 
    void PackedQ<T>::RDiv(
        const GenMatrix<T1>& m, const MatrixView<T2>& x) const
    {
        TMVAssert(Q.colsize() >= Q.rowsize());
        TMVAssert(beta.size() == Q.rowsize());
        TMVAssert(x.rowsize() == Q.colsize());
        TMVAssert(m.rowsize() == Q.rowsize());
        TMVAssert(x.colsize() == m.colsize());
        TMVAssert(Q.ct() == NonConj);
        TMVAssert(beta.ct() == NonConj);

        x.colRange(0,m.rowsize()) = m;
        x.colRange(m.rowsize(),x.rowsize()).setZero();
        Q_RDivEq(Q,beta,x);
    }

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_PackedQ.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


