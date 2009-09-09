///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2007                                                        //
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
// along with this program in the file gpl.txt.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



#include "TMV_Blas.h"
#include "TMV_Matrix.h"
#include "TMV_TriMatrix.h"
#include "TMV_DiagMatrix.h"
#include "TMV_QRDiv.h"
#include "TMV_Householder.h"

//#define XDEBUG

#ifdef XDEBUG
#include "TMV_MatrixArith.h"
#include <iostream>
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
 
  template <class T> inline void NonBlockGetQFromQR(
      const MatrixView<T>& Q, const GenVector<T>& beta)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(beta.ct() == NonConj);
    const size_t M = Q.colsize();
    const size_t N = Q.rowsize();
    UpperTriMatrixViewOf(Q).Zero();
    const int sb = beta.step();
    if (sb == 1) {
      const T* bi = beta.cptr()+(N-1);
      for(int i=N-1;i>=0;--i,--bi) 
	Householder_Unpack(Q.SubMatrix(i,M,i,N),*bi);
    }
    else {
      const T* bi = beta.cptr()+(N-1)*sb;
      for(int i=N-1;i>=0;--i,bi-=sb) 
	Householder_Unpack(Q.SubMatrix(i,M,i,N),*bi);
    }
  }

  template <class T> inline void BlockGetQFromQR(
      const MatrixView<T>& Q, const GenVector<T>& beta)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(beta.ct() == NonConj);
#ifdef XDEBUG
    Matrix<T> Q0(Q);
    Matrix<T> Q2(Q);
    NonBlockGetQFromQR(Q2.View(),beta);
#endif
    const size_t M = Q.colsize();
    const size_t N = Q.rowsize();
    UpperTriMatrixViewOf(Q).Zero();
    UpperTriMatrix<T,NonUnitDiag,ColMajor> BaseZ(
	std::min(size_t(QR_BLOCKSIZE),N));
    for(size_t j2=N;j2>0;) {
      size_t j1 = j2 > QR_BLOCKSIZE ? j2-QR_BLOCKSIZE : 0;
      MatrixView<T> Y = Q.SubMatrix(j1,M,j1,j2);
      UpperTriMatrixView<T> Z = BaseZ.SubTriMatrix(0,Y.rowsize());
      BlockHouseholder_MakeZ(Y,Z,beta.SubVector(j1,j2));
      BlockHouseholder_Unpack(Y,Z,Q.SubMatrix(j1,M,j2,N));
      j2 = j1;
    }
#ifdef XDEBUG
    if (Norm(Q-Q2) > 0.001*Norm(Q)*Norm(beta)) {
      cerr<<"BlockGetQ: Q = "<<Type(Q)<<"  "<<Q0<<endl;
      cerr<<"beta = "<<beta<<endl;
      cerr<<"-> Q = "<<Q<<endl;
      cerr<<"NonBlock Q = "<<Q2<<endl;
      abort(); 
    }
#endif
  }

  template <class T> inline void NonLapGetQFromQR(
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
  template <class T> inline void LapGetQFromQR(
      const MatrixView<T>& Q, const GenVector<T>& beta)
  { NonLapGetQFromQR(Q,beta); }
#ifdef INST_DOUBLE
  template <> inline void LapGetQFromQR(const MatrixView<double>& Q,
      const GenVector<double>& beta)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int m = Q.colsize();
    int n = Q.rowsize();
#ifndef LAPNOWORK
    int lwork = n*LAP_BLOCKSIZE;
    double* work = LAP_DWork(lwork);
#endif
    if (Q.isrm()) {
      int ldq = Q.stepi();
      LAPNAME(dorglq) (LAPCM LAPV(n),LAPV(m),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
	  LAPP(beta.cptr()) LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("dorglq");
#else
      LAP_Results(int(work[0]),m,n,lwork,"dorglq");
#endif
    } else {
      int ldq = Q.stepj();
      LAPNAME(dorgqr) (LAPCM LAPV(m),LAPV(n),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
	  LAPP(beta.cptr()) LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("dorgqr");
#else
      LAP_Results(int(work[0]),m,n,lwork,"dorgqr");
#endif
    }
  }
  template <> inline void LapGetQFromQR(
      const MatrixView<std::complex<double> >& Q,
      const GenVector<std::complex<double> >& beta)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int m = Q.colsize();
    int n = Q.rowsize();
#ifndef LAPNOWORK
    int lwork = n*LAP_BLOCKSIZE;
    std::complex<double>* work = LAP_ZWork(lwork);
#endif
    if (Q.isrm()) {
      int ldq = Q.stepi();
      LAPNAME(zunglq) (LAPCM LAPV(n),LAPV(m),LAPV(n),
	  LAPP(Q.ptr()),LAPV(ldq),LAPP(beta.cptr())
	  LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("zunglq");
#else
      LAP_Results(int(REAL(work[0])),m,n,lwork,"zunglq");
#endif
    } else {
      int ldq = Q.stepj();
      Vector<std::complex<double> > conjbeta = beta.Conjugate();
      LAPNAME(zungqr) (LAPCM LAPV(m),LAPV(n),LAPV(n),
	  LAPP(Q.ptr()),LAPV(ldq),LAPP(conjbeta.cptr())
	  LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("zungqr");
#else
      LAP_Results(int(REAL(work[0])),m,n,lwork,"zungqr");
#endif
    }
  }
#endif
#ifdef INST_FLOAT
  template <> inline void LapGetQFromQR(const MatrixView<float>& Q,
      const GenVector<float>& beta)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int m = Q.colsize();
    int n = Q.rowsize();
#ifndef LAPNOWORK
    int lwork = n*LAP_BLOCKSIZE;
    float* work = LAP_SWork(lwork);
#endif
    if (Q.isrm()) {
      int ldq = Q.stepi();
      LAPNAME(sorglq) (LAPCM LAPV(n),LAPV(m),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
	  LAPP(beta.cptr()) LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("sorglq");
#else
      LAP_Results(int(work[0]),m,n,lwork,"sorglq");
#endif
    } else {
      int ldq = Q.stepj();
      LAPNAME(sorgqr) (LAPCM LAPV(m),LAPV(n),LAPV(n),LAPP(Q.ptr()),LAPV(ldq),
	  LAPP(beta.cptr()) LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("sorgqr");
#else
      LAP_Results(int(work[0]),m,n,lwork,"sorgqr");
#endif
    }
  }
  template <> inline void LapGetQFromQR(
      const MatrixView<std::complex<float> >& Q,
      const GenVector<std::complex<float> >& beta)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == beta.size());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(beta.ct() == NonConj);
    int m = Q.colsize();
    int n = Q.rowsize();
#ifndef LAPNOWORK
    int lwork = n*LAP_BLOCKSIZE;
    std::complex<float>* work = LAP_CWork(lwork);
#endif
    if (Q.isrm()) {
      int ldq = Q.stepi();
      LAPNAME(cunglq) (LAPCM LAPV(n),LAPV(m),LAPV(n),
	  LAPP(Q.ptr()),LAPV(ldq),LAPP(beta.cptr())
	  LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("cunglq");
#else
      LAP_Results(int(REAL(work[0])),m,n,lwork,"cunglq");
#endif
    } else {
      int ldq = Q.stepj();
      Vector<std::complex<float> > conjbeta = beta.Conjugate();
      LAPNAME(cungqr) (LAPCM LAPV(m),LAPV(n),LAPV(n),
	  LAPP(Q.ptr()),LAPV(ldq),LAPP(conjbeta.cptr())
	  LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("cungqr");
#else
      LAP_Results(int(REAL(work[0])),m,n,lwork,"cungqr");
#endif
    }
  }
#endif
#endif
  template <class T> void GetQFromQR(
      const MatrixView<T>& Q, const GenVector<T>& beta)
  {
    // Extract the Q matrix from a combined QRx matrix
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(beta.size() == Q.rowsize());
    TMVAssert(Q.isrm() || Q.iscm());
    TMVAssert(Q.ct() == NonConj);
    if (beta.isconj()) {
      Vector<T> b2 = beta;
      GetQFromQR(Q,b2.View());
    } else {
#ifdef LAP
      LapGetQFromQR(Q,beta);
#else
      NonLapGetQFromQR(Q,beta);
#endif
    }
  }

#define InstFile "TMV_QRDiv_D.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


