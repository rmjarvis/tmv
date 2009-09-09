///////////////////////////////////////////////////////////////////////////////
// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
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

  template <class T, class T1> static void NonBlockQ_LDivEq(
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
      Householder_LMult(Q.col(j,j+1,M),beta(j),m.Rows(j,M));
    }
  }

  template <class T, class T1> static void BlockQ_LDivEq(
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
    NonBlockQ_LDivEq(Q,beta,m2.View());
#endif

    // x = H_N-1 .. H1 H0 m
    // In first block step:
    // m = (Hr .. H0) m
    //   = (H0t .. Hrt)t m
    // So form Y,Z from Ht's, rather than H's, and then call LDiv

    const int M = Q.colsize();
    const int N = Q.rowsize();
    UpperTriMatrix<T1,NonUnitDiag,ColMajor> BaseZ(MIN(QR_BLOCKSIZE,N));
    for(int j1=0;j1<N;) {
      int j2 = MIN(N,j1+QR_BLOCKSIZE);
      ConstMatrixView<T1> Y = Q.SubMatrix(j1,M,j1,j2);
      UpperTriMatrixView<T1> Z = BaseZ.SubTriMatrix(0,Y.rowsize());
      BlockHouseholder_MakeZ(Y,Z,beta.SubVector(j1,j2));
      BlockHouseholder_LDiv(Y,Z,m.Rows(j1,M));
      j1 = j2;
    }
#ifdef XDEBUG
    if (Norm(m-m2) > 0.001*Norm(Q)*Norm(m0)) {
      cerr<<"BlockQ_LDivEq: Q = "<<TypeText(Q)<<"  "<<Q<<endl;
      cerr<<"beta = "<<beta<<endl;
      cerr<<"m = "<<TypeText(m)<<"  "<<m0<<endl;
      cerr<<"-> m = "<<m<<endl;
      cerr<<"NonBlock m = "<<m2<<endl;
      abort(); 
    }
#endif
  }

  template <class T, class T1> static void NonLapQ_LDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& beta,
      const MatrixView<T>& m)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.colsize() == m.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

    if (Q.rowsize() > QR_BLOCKSIZE && m.rowsize() > QR_BLOCKSIZE)
      BlockQ_LDivEq(Q,beta,m);
    else
      NonBlockQ_LDivEq(Q,beta,m);
  }
#ifdef LAP
  template <class T, class T1> static inline void LapQ_LDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& beta,
      const MatrixView<T>& m)
  { NonLapQ_LDivEq(Q,beta,m); }
#ifdef INST_DOUBLE
  template <> void LapQ_LDivEq(const GenMatrix<double>& Q,
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
      auto_array<double> work(new double[lwork]);
#else
      int lwork = -1;
      auto_array<double> work(new double[1]);
      LAPNAME(dormlq) (LAPCM x.isrm()?LAPCH_R:LAPCH_L,
          x.isrm()?LAPCH_T:LAPCH_NT,
          LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
          LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
          LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      lwork = int(work[0]);
      work.reset(new double[lwork]);
#endif
#endif
      LAPNAME(dormlq) (LAPCM x.isrm()?LAPCH_R:LAPCH_L,
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
      auto_array<double> work(new double[lwork]);
#else
      int lwork = -1;
      auto_array<double> work(new double[1]);
      LAPNAME(dormqr) (LAPCM x.isrm()?LAPCH_R:LAPCH_L, 
          x.isrm()?LAPCH_NT:LAPCH_T, LAPV(m),LAPV(n),LAPV(k),
          LAPP(Q.cptr()),LAPV(ldq),
          LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
          LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      lwork = int(work[0]);
      work.reset(new double[lwork]);
#endif
#endif
      LAPNAME(dormqr) (LAPCM x.isrm()?LAPCH_R:LAPCH_L, 
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
  template <> void LapQ_LDivEq(
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
      if (x.isrm() == x.isconj()) x.ConjugateSelf();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
      int lwork = x.rowsize()*LAP_BLOCKSIZE;
      auto_array<std::complex<double> > work(new std::complex<double>[lwork]);
#else
      int lwork = -1;
      auto_array<std::complex<double> > work(new std::complex<double>[1]);
      LAPNAME(zunmlq) (LAPCM x.isrm()?LAPCH_R:LAPCH_L,
          x.isrm()?LAPCH_CT:LAPCH_NT,
          LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
          LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
          LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      lwork = int(std::real(work[0]));
      work.reset(new std::complex<double>[lwork]);
#endif
#endif
      LAPNAME(zunmlq) (LAPCM x.isrm()?LAPCH_R:LAPCH_L,
          x.isrm()?LAPCH_CT:LAPCH_NT,
          LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
          LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
          LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      if (x.isrm() == x.isconj()) x.ConjugateSelf();
#ifdef LAPNOWORK
      LAP_Results("zunmlq");
#else
      LAP_Results(int(REAL(work[0])),m,n,lwork,"zunmlq");
#endif
    } else {
      int ldx = x.isrm() ? x.stepi() : x.stepj();
      int m = x.isrm() ? x.rowsize() : x.colsize();
      int n = x.isrm() ? x.colsize() : x.rowsize();
      int ldq = Q.stepj();
      if (x.isrm() != x.isconj()) x.ConjugateSelf();
      Vector<std::complex<double> > conjbeta = beta.Conjugate();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
      int lwork = x.rowsize()*LAP_BLOCKSIZE;
      auto_array<std::complex<double> > work(new std::complex<double>[lwork]);
#else
      int lwork = -1;
      auto_array<std::complex<double> > work(new std::complex<double>[1]);
      LAPNAME(zunmqr) (LAPCM x.isrm()?LAPCH_R:LAPCH_L,
          x.isrm()?LAPCH_NT:LAPCH_CT,
          LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
          LAPP(conjbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
          LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      lwork = int(std::real(work[0]));
      work.reset(new std::complex<double>[lwork]);
#endif
#endif
      LAPNAME(zunmqr) (LAPCM x.isrm()?LAPCH_R:LAPCH_L,
          x.isrm()?LAPCH_NT:LAPCH_CT,
          LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
          LAPP(conjbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
          LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      if (x.isrm() != x.isconj()) x.ConjugateSelf();
#ifdef LAPNOWORK
      LAP_Results("zunmqr");
#else
      LAP_Results(int(REAL(work[0])),m,n,lwork,"zunmqr");
#endif
    }
  }
#endif
#ifdef INST_FLOAT
  template <> void LapQ_LDivEq(const GenMatrix<float>& Q,
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
      auto_array<float> work(new float[lwork]);
#else
      int lwork = -1;
      auto_array<float> work(new float[1]);
      LAPNAME(sormlq) (LAPCM x.isrm()?LAPCH_R:LAPCH_L,
          x.isrm()?LAPCH_T:LAPCH_NT,
          LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
          LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
          LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      lwork = int(work[0]);
      work.reset(new float[lwork]);
#endif
#endif
      LAPNAME(sormlq) (LAPCM x.isrm()?LAPCH_R:LAPCH_L,
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
      auto_array<float> work(new float[lwork]);
#else
      int lwork = -1;
      auto_array<float> work(new float[1]);
      LAPNAME(sormqr) (LAPCM x.isrm()?LAPCH_R:LAPCH_L, 
          x.isrm()?LAPCH_NT:LAPCH_T, LAPV(m),LAPV(n),LAPV(k),
          LAPP(Q.cptr()),LAPV(ldq),
          LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
          LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      lwork = int(work[0]);
      work.reset(new float[lwork]);
#endif
#endif
      LAPNAME(sormqr) (LAPCM x.isrm()?LAPCH_R:LAPCH_L, 
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
  template <> void LapQ_LDivEq(
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
      if (x.isrm() == x.isconj()) x.ConjugateSelf();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
      int lwork = x.rowsize()*LAP_BLOCKSIZE;
      auto_array<std::complex<float> > work(new std::complex<float>[lwork]);
#else
      int lwork = -1;
      auto_array<std::complex<float> > work(new std::complex<float>[1]);
      LAPNAME(cunmlq) (LAPCM x.isrm()?LAPCH_R:LAPCH_L,
          x.isrm()?LAPCH_CT:LAPCH_NT,
          LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
          LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
          LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      lwork = int(std::real(work[0]));
      work.reset(new std::complex<float>[lwork]);
#endif
#endif
      LAPNAME(cunmlq) (LAPCM x.isrm()?LAPCH_R:LAPCH_L,
          x.isrm()?LAPCH_CT:LAPCH_NT,
          LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
          LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
          LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      if (x.isrm() == x.isconj()) x.ConjugateSelf();
#ifdef LAPNOWORK
      LAP_Results("cunmlq");
#else
      LAP_Results(int(REAL(work[0])),m,n,lwork,"cunmlq");
#endif
    } else {
      int ldx = x.isrm() ? x.stepi() : x.stepj();
      int m = x.isrm() ? x.rowsize() : x.colsize();
      int n = x.isrm() ? x.colsize() : x.rowsize();
      int ldq = Q.stepj();
      if (x.isrm() != x.isconj()) x.ConjugateSelf();
      Vector<std::complex<float> > conjbeta = beta.Conjugate();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
      int lwork = x.rowsize()*LAP_BLOCKSIZE;
      auto_array<std::complex<float> > work(new std::complex<float>[lwork]);
#else
      int lwork = -1;
      auto_array<std::complex<float> > work(new std::complex<float>[1]);
      LAPNAME(cunmqr) (LAPCM x.isrm()?LAPCH_R:LAPCH_L,
          x.isrm()?LAPCH_NT:LAPCH_CT,
          LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
          LAPP(conjbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
          LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      lwork = int(std::real(work[0]));
      work.reset(new std::complex<float>[lwork]);
#endif
#endif
      LAPNAME(cunmqr) (LAPCM x.isrm()?LAPCH_R:LAPCH_L,
          x.isrm()?LAPCH_NT:LAPCH_CT,
          LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
          LAPP(conjbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
          LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      if (x.isrm() != x.isconj()) x.ConjugateSelf();
#ifdef LAPNOWORK
      LAP_Results("cunmqr");
#else
      LAP_Results(int(REAL(work[0])),m,n,lwork,"cunmqr");
#endif
    }
  }
#endif
#endif

  template <class T, class T1> void Q_LDivEq(
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
    QQ.SetToIdentity();
    QQ.Cols(0,Q.rowsize()) = Q;
    Vector<T1> bb(Q.colsize(),0.);
    bb.SubVector(0,beta.size()) = beta;
    GetQFromQR(QQ.View(),bb);
    Matrix<T> m0 = m;
    Matrix<T> m2 = QQ.Adjoint() * m;
#endif

    if (m.colsize() > 0 && m.rowsize() > 0) {
#ifdef LAP
      if ( m.isrm() || m.iscm() )
        LapQ_LDivEq(Q,beta,m);
      else
#endif
        NonLapQ_LDivEq(Q,beta,m);
    }
#ifdef XDEBUG
    if (Norm(m-m2) > 0.001*Norm(m0)) {
      cerr<<"Q_LDivEq\n";
      cerr<<"Q = "<<Q<<endl;
      cerr<<"beta = "<<beta<<endl;
      cerr<<"m = "<<TypeText(m)<<"  "<<m0<<endl;
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

  template <class T, class T1> static void NonBlockQ_RDivEq(
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
      Householder_LMult(Q.col(j,j+1,M).Conjugate(),beta(j),
          m.Cols(j,M).Transpose());
    }
  }

  template <class T, class T1> static void BlockQ_RDivEq(
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
    NonBlockQ_RDivEq(Q,beta,m2.View());
#endif

    // x = m Qt 
    // x = m H_N-1 H_N-2 ... H1 H0
    // Again form Y,Z from Ht's, rather than H's, and then call RDiv

    const int M = Q.colsize();
    const int N = Q.rowsize();
    UpperTriMatrix<T1,NonUnitDiag,ColMajor> BaseZ(MIN(QR_BLOCKSIZE,N));
    for(int j2=N;j2>0;) {
      int j1 = j2 > QR_BLOCKSIZE ? j2-QR_BLOCKSIZE : 0;
      ConstMatrixView<T1> Y = Q.SubMatrix(j1,M,j1,j2);
      UpperTriMatrixView<T1> Z = BaseZ.SubTriMatrix(0,Y.rowsize());
      BlockHouseholder_MakeZ(Y,Z,beta.SubVector(j1,j2));
      BlockHouseholder_LMult(Y,Z,m.Cols(j1,M).Adjoint());
      j2 = j1;
    }
#ifdef XDEBUG
    if (Norm(m-m2) > 0.001*Norm(Q)*Norm(m0)) {
      cerr<<"BlockQ_RDivEq: Q = "<<TypeText(Q)<<"  "<<Q<<endl;
      cerr<<"beta = "<<beta<<endl;
      cerr<<"m = "<<TypeText(m)<<"  "<<m0<<endl;
      cerr<<"-> m = "<<m<<endl;
      cerr<<"NonBlock m = "<<m2<<endl;
      abort(); 
    }
#endif
  }

  template <class T, class T1> static void NonLapQ_RDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& beta,
      const MatrixView<T>& m)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.colsize() == m.rowsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

    if (Q.rowsize() > QR_BLOCKSIZE && m.colsize() > QR_BLOCKSIZE)
      BlockQ_RDivEq(Q,beta,m);
    else
      NonBlockQ_RDivEq(Q,beta,m);
  }

#ifdef LAP
  template <class T, class T1> static inline void LapQ_RDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& beta,
      const MatrixView<T>& m)
  { NonLapQ_RDivEq(Q,beta,m); }
#ifdef INST_DOUBLE
  template <> void LapQ_RDivEq(const GenMatrix<double>& Q,
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
      auto_array<double> work(new double[lwork]);
#else
      int lwork = -1;
      auto_array<double> work(new double[1]);
      LAPNAME(dormlq) (LAPCM x.isrm()?LAPCH_L:LAPCH_R,
          x.isrm()?LAPCH_T:LAPCH_NT,
          LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
          LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
          LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      lwork = int(work[0]);
      work.reset(new double[lwork]);
#endif
#endif
      LAPNAME(dormlq) (LAPCM x.isrm()?LAPCH_L:LAPCH_R,
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
      auto_array<double> work(new double[lwork]);
#else
      int lwork = -1;
      auto_array<double> work(new double[1]);
      LAPNAME(dormqr) (LAPCM x.isrm()?LAPCH_L:LAPCH_R, 
          x.isrm()?LAPCH_NT:LAPCH_T, LAPV(m),LAPV(n),LAPV(k),
          LAPP(Q.cptr()),LAPV(ldq),
          LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
          LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      lwork = int(work[0]);
      work.reset(new double[lwork]);
#endif
#endif
      LAPNAME(dormqr) (LAPCM x.isrm()?LAPCH_L:LAPCH_R, 
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
  template <> void LapQ_RDivEq(
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
      if (x.isrm() == x.isconj()) x.ConjugateSelf();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
      int lwork = x.colsize()*LAP_BLOCKSIZE;
      auto_array<std::complex<double> > work(new std::complex<double>[lwork]);
#else
      int lwork = -1;
      auto_array<std::complex<double> > work(new std::complex<double>[1]);
      LAPNAME(zunmlq) (LAPCM x.isrm()?LAPCH_L:LAPCH_R,
          x.isrm()?LAPCH_CT:LAPCH_NT,
          LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
          LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
          LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      lwork = int(std::real(work[0]));
      work.reset(new std::complex<double>[lwork]);
#endif
#endif
      LAPNAME(zunmlq) (LAPCM x.isrm()?LAPCH_L:LAPCH_R,
          x.isrm()?LAPCH_CT:LAPCH_NT,
          LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
          LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
          LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      if (x.isrm() == x.isconj()) x.ConjugateSelf();
#ifdef LAPNOWORK
      LAP_Results("zunmlq");
#else
      LAP_Results(int(REAL(work[0])),m,n,lwork,"zunmlq");
#endif
    } else {
      int ldx = x.isrm() ? x.stepi() : x.stepj();
      int m = x.isrm() ? x.rowsize() : x.colsize();
      int n = x.isrm() ? x.colsize() : x.rowsize();
      int ldq = Q.stepj();
      if (x.isrm() != x.isconj()) x.ConjugateSelf();
      Vector<std::complex<double> > conjbeta = beta.Conjugate();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
      int lwork = x.colsize()*LAP_BLOCKSIZE;
      auto_array<std::complex<double> > work(new std::complex<double>[lwork]);
#else
      int lwork = -1;
      auto_array<std::complex<double> > work(new std::complex<double>[1]);
      LAPNAME(zunmqr) (LAPCM x.isrm()?LAPCH_L:LAPCH_R,
          x.isrm()?LAPCH_NT:LAPCH_CT,
          LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
          LAPP(conjbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
          LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      lwork = int(std::real(work[0]));
      work.reset(new std::complex<double>[lwork]);
#endif
#endif
      LAPNAME(zunmqr) (LAPCM x.isrm()?LAPCH_L:LAPCH_R,
          x.isrm()?LAPCH_NT:LAPCH_CT,
          LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
          LAPP(conjbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
          LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      if (x.isrm() != x.isconj()) x.ConjugateSelf();
#ifdef LAPNOWORK
      LAP_Results("zunmqr");
#else
      LAP_Results(int(REAL(work[0])),m,n,lwork,"zunmqr");
#endif
    }
  }
#endif
#ifdef INST_FLOAT
  template <> void LapQ_RDivEq(const GenMatrix<float>& Q,
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
      auto_array<float> work(new float[lwork]);
#else
      int lwork = -1;
      auto_array<float> work(new float[1]);
      LAPNAME(sormlq) (LAPCM x.isrm()?LAPCH_L:LAPCH_R,
          x.isrm()?LAPCH_T:LAPCH_NT,
          LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
          LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
          LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      lwork = int(work[0]);
      work.reset(new float[lwork]);
#endif
#endif
      LAPNAME(sormlq) (LAPCM x.isrm()?LAPCH_L:LAPCH_R,
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
      auto_array<float> work(new float[lwork]);
#else
      int lwork = -1;
      auto_array<float> work(new float[1]);
      LAPNAME(sormqr) (LAPCM x.isrm()?LAPCH_L:LAPCH_R, 
          x.isrm()?LAPCH_NT:LAPCH_T, LAPV(m),LAPV(n),LAPV(k),
          LAPP(Q.cptr()),LAPV(ldq),
          LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
          LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      lwork = int(work[0]);
      work.reset(new float[lwork]);
#endif
#endif
      LAPNAME(sormqr) (LAPCM x.isrm()?LAPCH_L:LAPCH_R, 
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
  template <> void LapQ_RDivEq(
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
      if (x.isrm() == x.isconj()) x.ConjugateSelf();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
      int lwork = x.colsize()*LAP_BLOCKSIZE;
      auto_array<std::complex<float> > work(new std::complex<float>[lwork]);
#else
      int lwork = -1;
      auto_array<std::complex<float> > work(new std::complex<float>[1]);
      LAPNAME(cunmlq) (LAPCM x.isrm()?LAPCH_L:LAPCH_R,
          x.isrm()?LAPCH_CT:LAPCH_NT,
          LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
          LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
          LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      lwork = int(std::real(work[0]));
      work.reset(new std::complex<float>[lwork]);
#endif
#endif
      LAPNAME(cunmlq) (LAPCM x.isrm()?LAPCH_L:LAPCH_R,
          x.isrm()?LAPCH_CT:LAPCH_NT,
          LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
          LAPP(beta.cptr()),LAPP(x.ptr()),LAPV(ldx)
          LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      if (x.isrm() == x.isconj()) x.ConjugateSelf();
#ifdef LAPNOWORK
      LAP_Results("cunmlq");
#else
      LAP_Results(int(REAL(work[0])),m,n,lwork,"cunmlq");
#endif
    } else {
      int ldx = x.isrm() ? x.stepi() : x.stepj();
      int m = x.isrm() ? x.rowsize() : x.colsize();
      int n = x.isrm() ? x.colsize() : x.rowsize();
      int ldq = Q.stepj();
      if (x.isrm() != x.isconj()) x.ConjugateSelf();
      Vector<std::complex<float> > conjbeta = beta.Conjugate();
#ifndef LAPNOWORK
#ifdef NOWORKQUERY
      int lwork = x.colsize()*LAP_BLOCKSIZE;
      auto_array<std::complex<float> > work(new std::complex<float>[lwork]);
#else
      int lwork = -1;
      auto_array<std::complex<float> > work(new std::complex<float>[1]);
      LAPNAME(cunmqr) (LAPCM x.isrm()?LAPCH_L:LAPCH_R,
          x.isrm()?LAPCH_NT:LAPCH_CT,
          LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
          LAPP(conjbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
          LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      lwork = int(std::real(work[0]));
      work.reset(new std::complex<float>[lwork]);
#endif
#endif
      LAPNAME(cunmqr) (LAPCM x.isrm()?LAPCH_L:LAPCH_R,
          x.isrm()?LAPCH_NT:LAPCH_CT,
          LAPV(m),LAPV(n),LAPV(k),LAPP(Q.cptr()),LAPV(ldq),
          LAPP(conjbeta.cptr()),LAPP(x.ptr()),LAPV(ldx)
          LAPWK(work.get()) LAPVWK(lwork) LAPINFO LAP1 LAP1);
      if (x.isrm() != x.isconj()) x.ConjugateSelf();
#ifdef LAPNOWORK
      LAP_Results("cunmqr");
#else
      LAP_Results(int(REAL(work[0])),m,n,lwork,"cunmqr");
#endif
    }
  }
#endif // FLOAT
#endif // LAP

  template <class T, class T1> void Q_RDivEq(
      const GenMatrix<T1>& Q, const GenVector<T1>& beta,
      const MatrixView<T>& m)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(beta.size() == Q.rowsize());
    TMVAssert(m.rowsize() == Q.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

#ifdef XDEBUG
    Matrix<T1> QQ(Q.colsize(),Q.colsize(),0.);
    QQ.SetToIdentity();
    QQ.Cols(0,Q.rowsize()) = Q;
    Vector<T1> bb(Q.colsize(),0.);
    bb.SubVector(0,beta.size()) = beta;
    GetQFromQR(QQ.View(),bb);
    Matrix<T> m0 = m;
    Matrix<T> m2 = m * QQ.Adjoint();
#endif

    if (m.colsize() > 0 && m.rowsize() > 0) {
#ifdef LAP
      if (m.isrm() || m.iscm())
        LapQ_RDivEq(Q,beta,m);
      else
#endif
        NonLapQ_RDivEq(Q,beta,m);
    }

#ifdef XDEBUG
    if (Norm(m-m2) > 0.001*Norm(m0)) {
      cerr<<"Q_RDivEq\n";
      cerr<<"Q = "<<Q<<endl;
      cerr<<"beta = "<<beta<<endl;
      cerr<<"m = "<<TypeText(m)<<"  "<<m0<<endl;
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

  template <class T> static void UnPack(
      const GenMatrix<T>& Q, const GenVector<T>& beta,
      const MatrixView<T>& m)
  {
    if (m.isrm() || m.iscm()) {
      m = Q;
      GetQFromQR(m,beta);
    } else {
      Matrix<T> m1 = Q;
      GetQFromQR(m1.View(),beta);
      m = m1;
    }
  }

  template <class T, class T1> static void UnPack(
      const GenMatrix<T>& Q, const GenVector<T>& beta,
      const MatrixView<T1>& m)
  {
    Matrix<T> m1 = Q;
    GetQFromQR(m1.View(),beta);
    m = m1;
  }

  template <class T> template <class T1> void PackedQ<T>::DoAssignToM(
      const MatrixView<T1>& m) const
  { UnPack(Q,beta,m); }

  template <class T> template <class T1> void PackedQ<T>::LDivEq(
      const VectorView<T1>& v) const
  {
    TMVAssert(Q.colsize() == Q.rowsize());
    TMVAssert(beta.size() == Q.rowsize());
    TMVAssert(v.size() == Q.colsize());
    TMVAssert(Q.isrm() || Q.iscm());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

    Q_LDivEq(Q,beta,ColVectorViewOf(v));
  }

  template <class T> template <class T1> void PackedQ<T>::RDivEq(
      const VectorView<T1>& v) const
  {
    TMVAssert(Q.colsize() == Q.rowsize());
    TMVAssert(beta.size() == Q.rowsize());
    TMVAssert(v.size() == Q.colsize());
    TMVAssert(Q.isrm() || Q.iscm());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

    Q_RDivEq(Q,beta,RowVectorViewOf(v));
  }

  template <class T> template <class T1, class T2> void PackedQ<T>::LDiv(
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
    if (Q.IsSquare()) {
      x = v;
      Q_LDivEq(Q,beta,ColVectorViewOf(x));
    } else {
      Vector<T1> v1 = v;
      Q_LDivEq(Q,beta,ColVectorViewOf(v1));
      x = v1.SubVector(0,x.size());
    }
  }

  template <class T> template <class T1, class T2> void PackedQ<T>::RDiv(
      const GenVector<T1>& v, const VectorView<T2>& x) const
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(beta.size() == Q.rowsize());
    TMVAssert(x.size() == Q.colsize());
    TMVAssert(v.size() == Q.rowsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

    x.SubVector(0,v.size()) = v;
    x.SubVector(v.size(),x.size()).Zero();
    Q_RDivEq(Q,beta,RowVectorViewOf(x));
  }

  template <class T> template <class T1> void PackedQ<T>::LDivEq(
      const MatrixView<T1>& m) const
  {
    TMVAssert(Q.colsize() == Q.rowsize());
    TMVAssert(beta.size() == Q.rowsize());
    TMVAssert(m.colsize() == Q.colsize());
    TMVAssert(Q.isrm() || Q.iscm());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

    Q_LDivEq(Q,beta,m);
  }

  template <class T> template <class T1> void PackedQ<T>::RDivEq(
      const MatrixView<T1>& m) const
  {
    TMVAssert(Q.colsize() == Q.rowsize());
    TMVAssert(beta.size() == Q.rowsize());
    TMVAssert(m.rowsize() == Q.rowsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

    Q_RDivEq(Q,beta,m);
  }

  template <class T> template <class T1, class T2> void PackedQ<T>::LDiv(
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
    if (Q.IsSquare()) {
      x = m;
      Q_LDivEq(Q,beta,x);
    } else {
      Matrix<T,ColMajor> m1 = m;
      Q_LDivEq(Q,beta,m1.View());
      x = m1.Rows(0,x.colsize());
    }
  }

  template <class T> template <class T1, class T2> void PackedQ<T>::RDiv(
      const GenMatrix<T1>& m, const MatrixView<T2>& x) const
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(beta.size() == Q.rowsize());
    TMVAssert(x.rowsize() == Q.colsize());
    TMVAssert(m.rowsize() == Q.rowsize());
    TMVAssert(x.colsize() == m.colsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);

    x.Cols(0,m.rowsize()) = m;
    x.Cols(m.rowsize(),x.rowsize()).Zero();
    Q_RDivEq(Q,beta,x);
  }

#define InstFile "TMV_PackedQ.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


