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
#include "TMV_MatrixArithFunc.h"
#include "TMV_Matrix.h"
#include "TMV_MatrixArith.h"
#include "TMV_MultMM.h"

//#define XDEBUG

#ifdef XDEBUG
#include "TMV_VectorArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

  template <bool add, class T, class Ta, class Tb> static void NonBlasMultMM(
      const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.iscm());

    if (A.iscm()) 
      if (B.iscm()) 
	CCCMultMM<add>(alpha,A,B,C);
      else 
	CRCMultMM<add>(alpha,A,B,C);
    else
      if (B.iscm()) 
	RCCMultMM<add>(alpha,A,B,C);
      else {
	// With RRC, there is no way to make the innermost loop have
	// dual unit-stride vectors.  So it is always faster to just
	// copy one matrix to the opposite storage.  The fastest
	// algorithm is RCC, so it is best to copy the B matrix.
	int N=B.rowsize();
	if (N > MM_BLOCKSIZE) {
	  int j1=0;
	  int K=B.colsize();
	  Matrix<T,ColMajor> B1(K,MM_BLOCKSIZE);
	  for (int j2=MM_BLOCKSIZE;j2<N;j1=j2,j2+=MM_BLOCKSIZE) {
	    B1 = B.Cols(j1,j2);
	    RCCMultMM<add>(alpha,A,B1,C.Cols(j1,j2));
	  }
	  B1.Cols(0,N-j1) = B.Cols(j1,N);
	  RCCMultMM<add>(alpha,A,B1.Cols(0,N-j1),C.Cols(j1,N));
	} else {
	  Matrix<T,ColMajor> B1 = B;
	  RCCMultMM<add>(alpha,A,B1,C);
	}
      }
  }

#ifdef BLAS
  template <class T, class Ta, class Tb> static inline void BlasMultMM(
      const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const int beta, const MatrixView<T>& C)
  {
    if (beta == 0) NonBlasMultMM<false>(alpha,A,B,C); 
    else NonBlasMultMM<true>(alpha,A,B,C); 
  }
#ifdef INST_DOUBLE
  template <> void BlasMultMM(
      const double alpha, const GenMatrix<double>& A,
      const GenMatrix<double>& B, const int beta, const MatrixView<double>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != double(0));
    TMVAssert(A.ct()==NonConj);
    TMVAssert(B.ct()==NonConj);
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.iscm());

    int m = C.colsize();
    int n = C.rowsize();
    int k = A.rowsize();
    int lda = A.isrm()?A.stepi():A.stepj();
    int ldb = B.isrm()?B.stepi():B.stepj();
    int ldc = C.stepj();
    double xbeta(beta);
    BLASNAME(dgemm) (BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
	B.iscm()?BLASCH_NT:BLASCH_T,
	BLASV(m),BLASV(n),BLASV(k),BLASV(alpha),
	BLASP(A.cptr()),BLASV(lda),BLASP(B.cptr()),BLASV(ldb),
	BLASV(xbeta),BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
  }
  template <> void BlasMultMM(
      const std::complex<double> alpha,
      const GenMatrix<std::complex<double> >& A,
      const GenMatrix<std::complex<double> >& B,
      const int beta, const MatrixView<std::complex<double> >& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != double(0));
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.iscm());

    if (A.iscm() && A.isconj()) {
      Matrix<std::complex<double>,ColMajor> AA = alpha*A;
      return BlasMultMM(std::complex<double>(1),AA,B,beta,C);
    } else if (B.iscm() && B.isconj()) {
      Matrix<std::complex<double>,ColMajor> BB = alpha*B;
      return BlasMultMM(std::complex<double>(1),A,BB,beta,C);
    } else {
      int m = C.colsize();
      int n = C.rowsize();
      int k = A.rowsize();
      int lda = A.isrm()?A.stepi():A.stepj();
      int ldb = B.isrm()?B.stepi():B.stepj();
      int ldc = C.stepj();
      std::complex<double> xbeta(beta);
      BLASNAME(zgemm) (BLASCM A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
	  B.iscm()?BLASCH_NT:B.isconj()?BLASCH_CT:BLASCH_T,
	  BLASV(m),BLASV(n),BLASV(k),BLASP(&alpha),
	  BLASP(A.cptr()),BLASV(lda),BLASP(B.cptr()),BLASV(ldb),
	  BLASP(&xbeta),BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
    }
  }
  template <> void BlasMultMM(
      const std::complex<double> alpha,
      const GenMatrix<std::complex<double> >& A,
      const GenMatrix<double>& B,
      const int beta, const MatrixView<std::complex<double> >& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != double(0));
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.iscm());

    if (A.iscm() && ((!A.isconj() && IMAG(alpha)==double(0)) || beta == 0)) {
      int m = 2*C.colsize();
      int n = C.rowsize();
      int k = A.rowsize();
      int lda = 2*A.stepj();
      int ldb = B.isrm()?B.stepi():B.stepj();
      int ldc = 2*C.stepj();
      if (beta == 0) {
	double xalpha(1);
	double xbeta(0);
	BLASNAME(dgemm) (BLASCM BLASCH_NT, B.iscm()?BLASCH_NT:BLASCH_T,
	    BLASV(m),BLASV(n),BLASV(k),BLASV(xalpha),
	    BLASP((const double*)(A.cptr())),BLASV(lda),
	    BLASP(B.cptr()),BLASV(ldb),
	    BLASV(xbeta),BLASP((double*)(C.ptr())),BLASV(ldc) BLAS1 BLAS1);
	if (A.isconj()) C.ConjugateSelf();
	C *= alpha;
      } else {
	double xalpha(REAL(alpha));
	double xbeta(beta);
	BLASNAME(dgemm) (BLASCM BLASCH_NT, B.iscm()?BLASCH_NT:BLASCH_T,
	    BLASV(m),BLASV(n),BLASV(k),BLASV(xalpha),
	    BLASP((const double*)(A.cptr())),BLASV(lda),
	    BLASP(B.cptr()),BLASV(ldb),
	    BLASV(xbeta),BLASP((double*)(C.ptr())),BLASV(ldc) BLAS1 BLAS1);
      } 
    } else {
      if (IMAG(alpha) == double(0)) {
	Matrix<double,ColMajor> A1 = A.Real();
	Matrix<double,ColMajor> C1 = REAL(alpha)*A1*B;
	if (beta == 0) C.Real() = C1;
	else C.Real() += C1;
	if (A.isconj()) C1 = -REAL(alpha)*(A1=A.Conjugate().Imag())*B;
	else C1 = REAL(alpha)*(A1=A.Imag())*B;
	if (beta == 0) C.Imag() = C1;
	else C.Imag() += C1;
      } else {
	Matrix<double,ColMajor> Ar = A.Real();
	Matrix<double,ColMajor> Ai = A.isconj()?A.Conjugate().Imag():A.Imag();
	Matrix<double,ColMajor> C1 = REAL(alpha)*Ar*B;
	if (A.isconj()) C1 += IMAG(alpha)*Ai*B;
	else C1 -= IMAG(alpha)*Ai*B;
	if (beta == 0) C.Real() = C1;
	else C.Real() += C1;

	if (A.isconj()) C1 = -REAL(alpha)*Ai*B;
	else C1 = REAL(alpha)*Ai*B;
	C1 += IMAG(alpha)*Ar*B;
	if (beta == 0) C.Imag() = C1;
	else C.Imag() += C1;
      }
    }
  }
  template <> void BlasMultMM(
      const std::complex<double> alpha,
      const GenMatrix<double>& A,
      const GenMatrix<std::complex<double> >& B,
      const int beta, const MatrixView<std::complex<double> >& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != double(0));
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.iscm());

    if (IMAG(alpha) == double(0)) {
      Matrix<double,ColMajor> B1 = B.Real();
      Matrix<double,ColMajor> C1 = REAL(alpha)*A*B1;
      if (beta == 0) C.Real() = C1;
      else C.Real() += C1;
      B1 = B.Imag();
      if (B.isconj()) C1 = -REAL(alpha)*A*B1;
      else C1 = REAL(alpha)*A*B1;
      if (beta == 0) C.Imag() = C1;
      else C.Imag() += C1;
    } else {
      Matrix<double,ColMajor> Br = B.Real();
      Matrix<double,ColMajor> Bi = B.Imag();
      Matrix<double,ColMajor> C1 = REAL(alpha)*A*Br;
      if (B.isconj()) C1 += IMAG(alpha)*A*Bi;
      else C1 -= IMAG(alpha)*A*Bi;
      if (beta == 0) C.Real() = C1;
      else C.Real() += C1;

      if (B.isconj()) C1 = -REAL(alpha)*A*Bi;
      else C1 = REAL(alpha)*A*Bi;
      C1 += IMAG(alpha)*A*Br;
      if (beta == 0) C.Imag() = C1;
      else C.Imag() += C1;
    }
  }
  template <> void BlasMultMM(
      const std::complex<double> alpha,
      const GenMatrix<double>& A,
      const GenMatrix<double>& B,
      const int beta, const MatrixView<std::complex<double> >& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != double(0));
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.iscm());

    Matrix<double,ColMajor> C1 = A*B;
    if (beta == 0) C = alpha*C1;
    else C += alpha*C1;
  }
#endif // INST_DOUBLE
#ifdef INST_FLOAT
  template <> void BlasMultMM(
      const float alpha, const GenMatrix<float>& A,
      const GenMatrix<float>& B, const int beta, const MatrixView<float>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != float(0));
    TMVAssert(A.ct()==NonConj);
    TMVAssert(B.ct()==NonConj);
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.iscm());

    int m = C.colsize();
    int n = C.rowsize();
    int k = A.rowsize();
    int lda = A.isrm()?A.stepi():A.stepj();
    int ldb = B.isrm()?B.stepi():B.stepj();
    int ldc = C.stepj();
    float xbeta(beta);
    BLASNAME(sgemm) (BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
	B.iscm()?BLASCH_NT:BLASCH_T,
	BLASV(m),BLASV(n),BLASV(k),BLASV(alpha),
	BLASP(A.cptr()),BLASV(lda),BLASP(B.cptr()),BLASV(ldb),
	BLASV(xbeta),BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
  }
  template <> void BlasMultMM(
      const std::complex<float> alpha,
      const GenMatrix<std::complex<float> >& A,
      const GenMatrix<std::complex<float> >& B,
      const int beta, const MatrixView<std::complex<float> >& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != float(0));
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.iscm());

    if (A.iscm() && A.isconj()) {
      Matrix<std::complex<float> > AA = alpha*A;
      return BlasMultMM(std::complex<float>(1),AA,B,beta,C);
    } else if (B.iscm() && B.isconj()) {
      Matrix<std::complex<float> > BB = alpha*B;
      return BlasMultMM(std::complex<float>(1),A,BB,beta,C);
    } else {
      int m = C.colsize();
      int n = C.rowsize();
      int k = A.rowsize();
      int lda = A.isrm()?A.stepi():A.stepj();
      int ldb = B.isrm()?B.stepi():B.stepj();
      int ldc = C.stepj();
      std::complex<float> xbeta(beta);
      BLASNAME(cgemm) (BLASCM A.iscm()?BLASCH_NT:A.isconj()?BLASCH_CT:BLASCH_T,
	  B.iscm()?BLASCH_NT:B.isconj()?BLASCH_CT:BLASCH_T,
	  BLASV(m),BLASV(n),BLASV(k),BLASP(&alpha),
	  BLASP(A.cptr()),BLASV(lda),BLASP(B.cptr()),BLASV(ldb),
	  BLASP(&xbeta),BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
    }
  }
  template <> void BlasMultMM(
      const std::complex<float> alpha,
      const GenMatrix<std::complex<float> >& A,
      const GenMatrix<float>& B,
      const int beta, const MatrixView<std::complex<float> >& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != float(0));
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.iscm());

    if (A.iscm() && !A.isconj() && (IMAG(alpha)==float(0) || beta == 0)) {
      int m = 2*C.colsize();
      int n = C.rowsize();
      int k = A.rowsize();
      int lda = 2*A.stepj();
      int ldb = B.isrm()?B.stepi():B.stepj();
      int ldc = 2*C.stepj();
      if (IMAG(alpha)==float(0)) {
	float xalpha(REAL(alpha));
	float xbeta(beta);
	BLASNAME(sgemm) (BLASCM BLASCH_NT, B.iscm()?BLASCH_NT:BLASCH_T,
	    BLASV(m),BLASV(n),BLASV(k),BLASV(xalpha),
	    BLASP((const float*)(A.cptr())),BLASV(lda),
	    BLASP(B.cptr()),BLASV(ldb),
	    BLASV(xbeta),BLASP((float*)(C.ptr())),BLASV(ldc) BLAS1 BLAS1);
      } else { // beta == 0
	float xalpha(1);
	float xbeta(0);
	BLASNAME(sgemm) (BLASCM BLASCH_NT, B.iscm()?BLASCH_NT:BLASCH_T,
	    BLASV(m),BLASV(n),BLASV(k),BLASV(xalpha),
	    BLASP((const float*)(A.cptr())),BLASV(lda),
	    BLASP(B.cptr()),BLASV(ldb),
	    BLASV(xbeta),BLASP((float*)(C.ptr())),BLASV(ldc) BLAS1 BLAS1);
	C *= alpha;
      } 
    } else {
      if (IMAG(alpha) == float(0)) {
	Matrix<float,ColMajor> A1 = A.Real();
	Matrix<float,ColMajor> C1 = REAL(alpha)*A1*B;
	if (beta == 0) C.Real() = C1;
	else C.Real() += C1;
	if (A.isconj()) C1 = -REAL(alpha)*(A1=A.Conjugate().Imag())*B;
	else C1 = REAL(alpha)*(A1=A.Imag())*B;
	if (beta == 0) C.Imag() = C1;
	else C.Imag() += C1;
      } else {
	Matrix<float,ColMajor> Ar = A.Real();
	Matrix<float,ColMajor> Ai = A.isconj()?A.Conjugate().Imag():A.Imag();
	Matrix<float,ColMajor> C1 = REAL(alpha)*Ar*B;
	if (A.isconj()) C1 += IMAG(alpha)*Ai*B;
	else C1 -= IMAG(alpha)*Ai*B;
	if (beta == 0) C.Real() = C1;
	else C.Real() += C1;

	if (A.isconj()) C1 = -REAL(alpha)*Ai*B;
	else C1 = REAL(alpha)*Ai*B;
	C1 += IMAG(alpha)*Ar*B;
	if (beta == 0) C.Imag() = C1;
	else C.Imag() += C1;
      }
    }
  }
  template <> void BlasMultMM(
      const std::complex<float> alpha,
      const GenMatrix<float>& A,
      const GenMatrix<std::complex<float> >& B,
      const int beta, const MatrixView<std::complex<float> >& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != float(0));
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.iscm());

    if (IMAG(alpha) == float(0)) {
      Matrix<float,ColMajor> B1 = B.Real();
      Matrix<float,ColMajor> C1 = REAL(alpha)*A*B1;
      if (beta == 0) C.Real() = C1;
      else C.Real() += C1;
      B1 = B.Imag();
      if (B.isconj()) C1 = -REAL(alpha)*A*B1;
      else C1 = REAL(alpha)*A*B1;
      if (beta == 0) C.Imag() = C1;
      else C.Imag() += C1;
    } else {
      Matrix<float,ColMajor> Br = B.Real();
      Matrix<float,ColMajor> Bi = B.Imag();
      Matrix<float,ColMajor> C1 = REAL(alpha)*A*Br;
      if (B.isconj()) C1 += IMAG(alpha)*A*Bi;
      else C1 -= IMAG(alpha)*A*Bi;
      if (beta == 0) C.Real() = C1;
      else C.Real() += C1;

      if (B.isconj()) C1 = -REAL(alpha)*A*Bi;
      else C1 = REAL(alpha)*A*Bi;
      C1 += IMAG(alpha)*A*Br;
      if (beta == 0) C.Imag() = C1;
      else C.Imag() += C1;
    }
  }
  template <> void BlasMultMM(
      const std::complex<float> alpha,
      const GenMatrix<float>& A,
      const GenMatrix<float>& B,
      const int beta, const MatrixView<std::complex<float> >& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != float(0));
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.iscm());

    Matrix<float,ColMajor> C1 = A*B;
    if (beta == 0) C = alpha*C1;
    else C += alpha*C1;
  }
#endif // INST_FLOAT
#endif // BLAS

  template <bool add, class T, class Ta, class Tb> static void DoMultMM(
      const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct() == NonConj);
    TMVAssert(A.isrm() || A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.iscm());

#ifdef BLAS
    if (IsComplex(T()) && (IsReal(Ta()) || IsReal(Tb())))
      BlasMultMM(alpha,A,B,add?1:0,C);
    else if (!((A.isrm() && A.stepi()>0) || (A.iscm() && A.stepj()>0))) {
      Matrix<T,ColMajor> A2 = alpha*A;
      DoMultMM<add>(T(1),A2,B,C);
    } else if (!((B.isrm() && B.stepi()>0) || (B.iscm() && B.stepj()>0))) {
      Matrix<T,ColMajor> B2 = alpha*B;
      DoMultMM<add>(T(1),A,B2,C);
    } else {
      BlasMultMM(alpha,A,B,add?1:0,C);
    }
#else
    NonBlasMultMM<add>(alpha,A,B,C);
#endif
  }

  template <bool add, class T, class Ta, class Tb> static void FullTempMultMM(
      const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
    // C (+)= alpha * A * B  via a temporary
  {
    Matrix<T,ColMajor> C2(C.colsize(),C.rowsize());
    DoMultMM<false>(T(1),A,B,C2.View());

    if (add) C += alpha*C2;
    else C = alpha*C2;
  }

  // Block Temp allows B to be the same storage as C
  template <bool add, class T, class Ta, class Tb> static void BlockTempMultMM(
      const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
    // C (+)= alpha * A * C
  {
    const int N = C.rowsize();
    for (int j=0,j2;j<N;j=j2) {
      j2 = MIN(N,j+MM_BLOCKSIZE);
      if (C.isrm()) {
	Matrix<T,ColMajor> B2 = alpha * B.Cols(j,j2);
	DoMultMM<add>(T(1),B2.Transpose(),A.Transpose(),
	    C.Cols(j,j2).Transpose());
      } else {
	Matrix<T,ColMajor> B2 = alpha * B.Cols(j,j2);
	DoMultMM<add>(T(1),A,B2,C.Cols(j,j2));
      }
    }
  }


  template <bool add, class T, class Ta, class Tb> void MultMM(
      const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
    // C (+)= alpha * A * B
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
#ifdef XDEBUG
    Matrix<T> A0 = A;
    Matrix<T> B0 = B;
    Matrix<T> C0 = C;
    Matrix<T> C2 = C;
    for(size_t i=0;i<C.colsize();i++)
      for(size_t j=0;j<C.rowsize();j++)
	C2(i,j) = A0.row(i) * B0.col(j);
    C2 *= alpha;
    if (add) C2 += C0;
    //cout<<"MultMM: add = "<<add<<", alpha = "<<alpha<<endl;
    //cout<<"A = "<<Type(A)<<" "<<A0<<endl;
    //cout<<"B = "<<Type(B)<<" "<<B0<<endl;
    //cout<<"C = "<<Type(C)<<" "<<C0<<endl;
#endif

    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (A.rowsize() == 0 || alpha == T(0))  {
	if (!add) C.Zero();
      }
      else if (C.isconj()) 
	MultMM<add>(CONJ(alpha),A.Conjugate(),B.Conjugate(),C.Conjugate());
      else if (C.isrm()) 
	MultMM<add>(alpha,B.Transpose(),A.Transpose(),C.Transpose());
      else if (!(B.isrm() || B.iscm())) 
	MultMM<add>(alpha,A,Matrix<Tb,ColMajor>(B),C);
      else if (!(A.isrm() || A.iscm())) 
	if (B.iscm())
	  MultMM<add>(alpha,Matrix<Ta,RowMajor>(A),B,C);
	else
	  MultMM<add>(alpha,Matrix<Ta,ColMajor>(A),B,C);
      else if (!C.iscm()) 
	FullTempMultMM<add>(alpha,A,B,C);
      else if (SameStorage(A,C)) 
	if (SameStorage(B,C)) 
	  FullTempMultMM<add>(alpha,A,B,C);
	else if (C.stepi() == A.stepi() && C.stepj() == A.stepj())
	  BlockTempMultMM<add>(alpha,B.Transpose(),A.Transpose(),
	      C.Transpose());
	else
	  FullTempMultMM<add>(alpha,A,B,C);
      else if (SameStorage(B,C))
	if (C.stepi() == B.stepi() && C.stepj() == B.stepj())
	  BlockTempMultMM<add>(alpha,A,B,C);
	else
	  FullTempMultMM<add>(alpha,A,B,C);
      else
	DoMultMM<add>(alpha,A,B,C);
    }
      
#ifdef XDEBUG
    if (Norm(C2-C) > 0.001*(ABS(alpha)*Norm(A0)*Norm(B0)+
	  (add?Norm(C0):RealType(T)(0)))) {
      cerr<<"MultMM: alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"A = "<<Type(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"--> C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_MultMM.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


