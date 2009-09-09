///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2008                                                        //
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



#include "TMV_Blas.h"
#include "TMV_QRDiv.h"
#include "TMV_QRD.h"
#include "TMV_Matrix.h"
#include "TMV_TriMatrix.h"
#include "TMV_Householder.h"
#include "TMV_TriMatrixArith.h"
#include "TMV_VectorArith.h"
#include "TMV_MatrixArith.h"

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
  // QR Decompose
  //

  template <class T> static void NonBlockQR_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& beta, T& det)
  {
#ifdef XDEBUG
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
    const int M = A.colsize();
    const int N = A.rowsize();

    T* bj = beta.ptr();
    for(int j=0;j<N;++j,++bj) {
      // Apply the Householder Reflection for this column
#ifdef TMVFLDEBUG
      TMVAssert(bj >= beta.first);
      TMVAssert(bj < beta.last);
#endif
      *bj = Householder_Reflect(A.SubMatrix(j,M,j,N),det);
    }
#ifdef XDEBUG
    Matrix<T> R(UpperTriMatrixViewOf(A));
    Matrix<T> Q(A);
    GetQFromQR(Q.View(),beta);
    Matrix<T> AA = Q*R;
    if (Norm(AA-A0) > 0.0001*Norm(Q)*Norm(R)) {
      cerr<<"NonBlockQR_Decompose: A = "<<Type(A)<<"  "<<A0<<endl;
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

  template <class T> static void RecursiveQR_Decompose(
      const MatrixView<T>& A, const UpperTriMatrixView<T>& Z, T& det,
      bool makeZ)
  {
#ifdef XDEBUG
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
    // beta.Conjugate().)
    const int M = A.colsize();
    const int N = A.rowsize();

    if (N==1) {
      T b = Householder_Reflect(A.col(0),det);
#ifdef TMVFLDEBUG
      TMVAssert(Z.ptr() >= Z.first);
      TMVAssert(Z.ptr() < Z.last);
#endif
      *Z.ptr() = CONJ(b);
    } else if (N==2) {
      T* Z00 = Z.ptr();
      T* Z01 = Z00 + Z.stepj();
      T* Z11 = Z01 + 1;

      T b0 = Householder_Reflect(A,det);
#ifdef TMVFLDEBUG
      TMVAssert(Z00 >= Z.first);
      TMVAssert(Z00 < Z.last);
      TMVAssert(Z01 >= Z.first);
      TMVAssert(Z01 < Z.last);
      TMVAssert(Z11 >= Z.first);
      TMVAssert(Z11 < Z.last);
#endif
      *Z00 = CONJ(b0);
      T b1 = Householder_Reflect(A.col(1,1,M),det);
      *Z11 = CONJ(b1);

      if (makeZ) {
	const T* A10 = A.cptr()+A.stepi();
	T temp = A.col(0,2,M).Conjugate()*A.col(1,2,M);
	temp += CONJ(*A10);
	*Z01 = -CONJ(b0*b1)*temp;
      }
    } else {
      int j1 = (N+1)/2;
      MatrixView<T> A1 = A.Cols(0,j1);
      UpperTriMatrixView<T> Z1 = Z.SubTriMatrix(0,j1);
      RecursiveQR_Decompose(A1,Z1,det,true);

      BlockHouseholder_LDiv(A1,Z1,A.Cols(j1,N));

      MatrixView<T> A2 = A.SubMatrix(j1,M,j1,N);
      UpperTriMatrixView<T> Z2 = Z.SubTriMatrix(j1,N);
      RecursiveQR_Decompose(A2,Z2,det,makeZ);

      if (makeZ) {
	MatrixView<T> Z3 = Z.SubMatrix(0,j1,j1,N);
	Z3 = A1.Rows(j1,N).Adjoint() *
	  LowerTriMatrixViewOf(A.SubMatrix(j1,N,j1,N),UnitDiag);
	Z3 += A1.Rows(N,M).Adjoint() * A.SubMatrix(N,M,j1,N);
	Z3 = -Z1*Z3;
	Z3 *= Z2;
      }
    }
#ifdef XDEBUG
    Matrix<T> R(UpperTriMatrixViewOf(A));
    Matrix<T> Q(A);
    GetQFromQR(Q.View(),Z.diag().Conjugate());
    Matrix<T> AA = Q*R;
    if (Norm(AA-A0) > 0.0001*Norm(Q)*Norm(R)) {
      cerr<<"RecursiveQR_Decompose: A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> "<<A<<endl;
      cerr<<"Z = "<<Z<<endl;
      cerr<<"beta = "<<Z.diag().Conjugate()<<endl;
      cerr<<"Q = "<<Q<<endl;
      cerr<<"R = "<<R<<endl;
      cerr<<"Q*R = "<<Q*R<<endl;
      cerr<<"Norm(diff) = "<<Norm(AA-A0)<<endl;
      Matrix<T> A2 = A0;
      Vector<T> beta2(Z.size());
      T det2(0);
      NonBlockQR_Decompose(A2.View(),beta2.View(),det2);
      cerr<<"NonBlock "<<A2<<endl;
      cerr<<"beta = "<<beta2<<endl;
      abort(); 
    }
#endif
  }

  template <class T> static void BlockQR_Decompose(
      const MatrixView<T>& A, const VectorView<T>& beta, T& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);

#ifdef XDEBUG
    Matrix<T> A0(A);
#endif
    const int M = A.colsize();
    const int N = A.rowsize();

    UpperTriMatrix<T,NonUnitDiag,ColMajor> BaseZ(
	MIN(QR_BLOCKSIZE,N));
    for(int j1=0;j1<N;) {
      int j2 = MIN(N,j1+QR_BLOCKSIZE);
      MatrixView<T> A1 = A.SubMatrix(j1,M,j1,j2);
      UpperTriMatrixView<T> Z = BaseZ.SubTriMatrix(0,j2-j1);

      RecursiveQR_Decompose(A1,Z,det,j2<N);
      beta.SubVector(j1,j2) = Z.diag().Conjugate();

      if (j2 < N) 
	BlockHouseholder_LDiv(A1,Z,A.SubMatrix(j1,M,j2,N));
      j1 = j2;
    }
#ifdef XDEBUG
    Matrix<T> R(UpperTriMatrixViewOf(A));
    Matrix<T> Q(A);
    GetQFromQR(Q.View(),beta);
    Matrix<T> AA = Q*R;
    if (Norm(AA-A0) > 0.0001*Norm(Q)*Norm(R)) {
      cerr<<"BlockQR_Decompose: A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> "<<A<<endl;
      cerr<<"beta = "<<beta<<endl;
      cerr<<"Q*R = "<<Q*R<<endl;
      cerr<<"Norm(diff) = "<<Norm(AA-A0)<<endl;
      Matrix<T> A2 = A0;
      Vector<T> beta2(beta.size());
      T det2(0);
      NonBlockQR_Decompose(A2.View(),beta2.View(),det2);
      cerr<<"NonBlock "<<A2<<endl;
      cerr<<"beta = "<<beta2<<endl;
      abort(); 
    }
#endif
  }

  template <class T> static void NonLapQR_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& beta, T& det)
  {
    TMVAssert(A.colsize() >= A.rowsize());
    TMVAssert(A.rowsize() == beta.size());
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(beta.ct() == NonConj);
    TMVAssert(beta.step()==1);

    if (A.rowsize() > QR_BLOCKSIZE)
      BlockQR_Decompose(A,beta,det);
    else {

      UpperTriMatrix<T,NonUnitDiag,ColMajor> Z(A.rowsize());
      RecursiveQR_Decompose(A,Z.View(),det,false);
      beta = Z.diag().Conjugate();
    }
  }

#ifdef LAP
  template <class T> static inline void LapQR_Decompose(
      const MatrixView<T>& A,
      const VectorView<T>& beta, T& det)
  { NonLapQR_Decompose(A,beta,det); }
#ifdef INST_DOUBLE
  template <> void LapQR_Decompose(const MatrixView<double>& A,
      const VectorView<double>& beta, double& det)
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
#ifndef LAPNOWORK
    int lwork = 2*n*LAP_BLOCKSIZE;
    double* work = LAP_DWork(lwork);
#endif
    if (A.isrm()) {
      int lda = A.stepi();
      LAPNAME(dgelqf) (LAPCM LAPV(n),LAPV(m),LAPP(A.ptr()),LAPV(lda),
	  LAPP(beta.ptr()) LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("dgelqf");
#else
      LAP_Results(int(work[0]),m,n,lwork,"dgelqf");
#endif
    } else {
      int lda = A.stepj();
      LAPNAME(dgeqrf) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
	  LAPP(beta.ptr()) LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("dgeqrf");
#else
      LAP_Results(int(work[0]),m,n,lwork,"dgeqrf");
#endif
    }
    const double* bi = beta.cptr();
    if (det) for(int i=0;i<n;++i,++bi) 
      if (*bi != 0.) det = -det;
  }
  template <> void LapQR_Decompose(
      const MatrixView<std::complex<double> >& A,
      const VectorView<std::complex<double> >& beta, std::complex<double>& det)
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
#ifndef LAPNOWORK
    int lwork = 2*n*LAP_BLOCKSIZE;
    std::complex<double>* work = LAP_ZWork(lwork);
#endif
    if (A.isrm()) {
      int lda = A.stepi();
      LAPNAME(zgelqf) (LAPCM LAPV(n),LAPV(m),LAPP(A.ptr()),LAPV(lda),
	  LAPP(beta.ptr()) LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("zgelqf");
#else
      LAP_Results(int(REAL(work[0])),m,n,lwork,"zgelqf");
#endif
    } else {
      int lda = A.stepj();
      LAPNAME(zgeqrf) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
	  LAPP(beta.ptr()) LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("zgeqrf");
#else
      LAP_Results(int(REAL(work[0])),m,n,lwork,"zgeqrf");
#endif
      beta.ConjugateSelf();
    }
    if (det!=double(0)) {
      const std::complex<double>* bi = beta.cptr();
      for(int i=0;i<n;++i,++bi)
	if (IMAG(*bi) != 0.) 
	  det *= -CONJ(*bi * *bi)/norm(*bi);
	else if (REAL(*bi) != 0.)
	  det = -det;
    }
  }
#endif
#ifdef INST_FLOAT
  template <> void LapQR_Decompose(const MatrixView<float>& A,
      const VectorView<float>& beta, float& det)
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
#ifndef LAPNOWORK
    int lwork = 2*n*LAP_BLOCKSIZE;
    float* work = LAP_SWork(lwork);
#endif
    if (A.isrm()) {
      int lda = A.stepi();
      LAPNAME(sgelqf) (LAPCM LAPV(n),LAPV(m),LAPP(A.ptr()),LAPV(lda),
	  LAPP(beta.ptr()) LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("sgelqf");
#else
      LAP_Results(int(work[0]),m,n,lwork,"sgelqf");
#endif
    } else {
      int lda = A.stepj();
      LAPNAME(sgeqrf) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
	  LAPP(beta.ptr()) LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("sgeqrf");
#else
      LAP_Results(int(work[0]),m,n,lwork,"sgeqrf");
#endif
    }
    const float* bi = beta.cptr();
    if (det) for(int i=0;i<n;++i,++bi) 
      if (*bi != 0.) det = -det;
  }
  template <> void LapQR_Decompose(
      const MatrixView<std::complex<float> >& A,
      const VectorView<std::complex<float> >& beta, std::complex<float>& det)
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
#ifndef LAPNOWORK
    int lwork = 2*n*LAP_BLOCKSIZE;
    std::complex<float>* work = LAP_CWork(lwork);
#endif
    if (A.isrm()) {
      int lda = A.stepi();
      LAPNAME(cgelqf) (LAPCM LAPV(n),LAPV(m),LAPP(A.ptr()),LAPV(lda),
	  LAPP(beta.ptr()) LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("cgelqf");
#else
      LAP_Results(int(REAL(work[0])),m,n,lwork,"cgelqf");
#endif
    } else {
      int lda = A.stepj();
      LAPNAME(cgeqrf) (LAPCM LAPV(m),LAPV(n),LAPP(A.ptr()),LAPV(lda),
	  LAPP(beta.ptr()) LAPWK(work) LAPVWK(lwork) LAPINFO);
#ifdef LAPNOWORK
      LAP_Results("cgeqrf");
#else
      LAP_Results(int(REAL(work[0])),m,n,lwork,"cgeqrf");
#endif
      beta.ConjugateSelf();
    }
    if (det!=float(0)) {
      const std::complex<float>* bi = beta.cptr();
      for(int i=0;i<n;++i,++bi)
	if (IMAG(*bi) != 0.) 
	  det *= -CONJ(*bi * *bi)/norm(*bi);
	else if (REAL(*bi) != 0.)
	  det = -det;
    }
  }
#endif
#endif
  template <class T> void QR_Decompose(
      const MatrixView<T>& A, const VectorView<T>& beta, T& det)
  {
#ifdef XDEBUG
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
      LapQR_Decompose(A,beta,det);
#else
      NonLapQR_Decompose(A,beta,det);
#endif
    }
#ifdef XDEBUG
    Matrix<T> R(UpperTriMatrixViewOf(A));
    Matrix<T> Q(A);
    GetQFromQR(Q.View(),beta);
    Matrix<T> AA = Q*R;
    if (Norm(AA-A0) > 0.0001*Norm(Q)*Norm(R)) {
      cerr<<"BlockQR_Decompose: A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"-> "<<A<<endl;
      cerr<<"beta = "<<beta<<endl;
      cerr<<"Q*R = "<<Q*R<<endl;
      cerr<<"Norm(diff) = "<<Norm(AA-A0)<<endl;
      Matrix<T> A2 = A0;
      Vector<T> beta2(beta.size());
      T det2(0);
      NonBlockQR_Decompose(A2.View(),beta2.View(),det2);
      cerr<<"NonBlock "<<A2<<endl;
      cerr<<"beta = "<<beta2<<endl;
      abort(); 
    }
#endif
  }


  //
  // QR Decompose - Unpacked
  //

  template <class T> void QR_Decompose(
      const MatrixView<T>& Q, const UpperTriMatrixView<T>& R, T& det)
  {
    // Decompose A (input as Q) into A = Q R 
    // where Q is unitary and R is upper triangular

    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(R.colsize() == Q.rowsize());
    TMVAssert(R.rowsize() == Q.rowsize());
    TMVAssert(Q.ct() == NonConj);
    TMVAssert(R.ct() == NonConj);

    Vector<T> beta(Q.rowsize());
    QR_Decompose(Q,beta.View(),det);
    R = UpperTriMatrixViewOf(Q);
    GetQFromQR(Q,beta.View());
  }

  template <class T> void QR_Decompose(
      const MatrixView<T>& Q, const UpperTriMatrixView<T>& R)
  {
    T d(0);
    if (Q.isconj()) {
      if (R.isconj()) {
	QR_Decompose(Q.Conjugate(),R.Conjugate(),d);
      } else {
	QR_Decompose(Q.Conjugate(),R,d);
	R.ConjugateSelf();
      }
    } else {
      if (R.isconj()) {
	QR_Decompose(Q,R.Conjugate(),d);
	R.ConjugateSelf();
      } else {
	QR_Decompose(Q,R,d);
      }
    }
  }

  template <class T> void QR_Decompose(const MatrixView<T>& A)
  {
    // Decompose A into Q R, but ignore Q.
    // On output, R is UpperTriMatrixViewOf(A).

    TMVAssert(A.colsize() >= A.rowsize());

    Vector<T> beta(A.rowsize());
    T d(0);
    if (A.isconj()) 
      QR_Decompose(A.Conjugate(),beta.View(),d);
    else
      QR_Decompose(A,beta.View(),d);
  }

#define InstFile "TMV_QRDecompose.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


