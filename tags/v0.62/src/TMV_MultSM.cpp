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
#include "tmv/TMV_SymMatrixArithFunc.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_SymMatrixArith.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_VectorArith.h"
#ifdef BLAS
#include "tmv/TMV_TriMatrixArith.h"
#endif

#ifdef XDEBUG
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
#define SYM_MM_BLOCKSIZE TMV_BLOCKSIZE
#define SYM_MM_BLOCKSIZE2 (TMV_BLOCKSIZE/2)
#else
#define SYM_MM_BLOCKSIZE 64
#define SYM_MM_BLOCKSIZE2 32
#endif

  //
  // MultMM
  //

  template <bool add, class T, class Ta, class Tb> static void RRowMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.uplo() == Lower);

    const int N = A.size();
    for(int j=0;j<N;++j) {
      if (add) C.row(j) += alpha * A.row(j,0,j+1) * B.Rows(0,j+1);
      else C.row(j) = alpha * A.row(j,0,j+1) * B.Rows(0,j+1);
      C.Rows(0,j) += alpha * A.col(j,0,j) ^ B.row(j);
    }
  }

  template <bool add, class T, class Ta, class Tb> static void CRowMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.uplo() == Lower);

    const int N = A.size();
    for(int j=N-1;j>=0;--j) {
      if (add) C.row(j) += alpha * A.row(j,j,N) * B.Rows(j,N);
      else C.row(j) = alpha * A.row(j,j,N) * B.Rows(j,N);
      C.Rows(j+1,N) += alpha * A.col(j,j+1,N) ^ B.row(j);
    }
  }

  template <bool add, class T, class Ta, class Tb> static inline void RowMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
    if (A.iscm()) CRowMultMM<add>(alpha,A,B,C);
    else RRowMultMM<add>(alpha,A,B,C);
  }

  template <bool add, class T, class Ta, class Tb> static void ColMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.uplo() == Lower);

    const int N = C.rowsize();
    for(int j=0;j<N;++j) 
      if (add) C.col(j) += alpha * A * B.col(j);
      else C.col(j) = alpha * A * B.col(j);
  }

  template <bool add, class T, class Ta, class Tb> static void RecursiveMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(C.ct()==NonConj);
    TMVAssert(A.uplo() == Lower);

    const int N = A.size();
    if (N <= SYM_MM_BLOCKSIZE2) {
      if (B.isrm() && C.isrm()) RowMultMM<add>(alpha,A,B,C);
      else if (B.iscm() && C.iscm()) ColMultMM<add>(alpha,A,B,C);
      else if (C.colsize() < C.rowsize()) RowMultMM<add>(alpha,A,B,C);
      else ColMultMM<add>(alpha,A,B,C);
    } else {
      int k = N/2;
      const int nb = SYM_MM_BLOCKSIZE;
      if (k > nb) k = k/nb*nb;

      // [ A00 A10t ] [ B0 ] = [ A00 B0 + A10t B1 ]
      // [ A10 A11  ] [ B1 ]   [ A10 B0 + A11 B1  ]

      ConstSymMatrixView<Ta> A00 = A.SubSymMatrix(0,k);
      ConstSymMatrixView<Ta> A11 = A.SubSymMatrix(k,N);
      ConstMatrixView<Ta> A10 = A.SubMatrix(k,N,0,k);
      ConstMatrixView<Tb> B0 = B.Rows(0,k);
      ConstMatrixView<Tb> B1 = B.Rows(k,N);
      MatrixView<T> C0 = C.Rows(0,k);
      MatrixView<T> C1 = C.Rows(k,N);

      RecursiveMultMM<add>(alpha,A00,B0,C0);
      RecursiveMultMM<add>(alpha,A11,B1,C1);
      C1 += alpha * A10 * B0;
      if (A.issym())
        C0 += alpha * A10.Transpose() * B1;
      else
        C0 += alpha * A10.Adjoint() * B1;
    }
  }

  template <bool add, class T, class Ta, class Tb> static void NonBlasMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));

    if (A.uplo() == Upper)
      if (A.isherm()) NonBlasMultMM<add>(alpha,A.Adjoint(),B,C);
      else NonBlasMultMM<add>(alpha,A.Transpose(),B,C);
    else if (C.isconj())
      NonBlasMultMM<add>(CONJ(alpha),A.Conjugate(),B.Conjugate(),
          C.Conjugate());
    else RecursiveMultMM<add>(alpha,A,B,C);
  }

#ifdef BLAS
  template <class T, class Ta, class Tb> static inline void BlasMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const int beta, const MatrixView<T>& C)
  { 
    if (beta == 1) NonBlasMultMM<true>(alpha,A,B,C); 
    else NonBlasMultMM<false>(alpha,A,B,C); 
  }
#ifdef INST_DOUBLE
  template <> void BlasMultMM(
      const double alpha, const GenSymMatrix<double>& A,
      const GenMatrix<double>& B, const int beta, const MatrixView<double>& C)
  {
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != 0.);
    TMVAssert(A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.isrm() || C.iscm());
    TMVAssert(B.stor() == C.stor());

    int m = C.iscm() ? C.colsize() : C.rowsize();
    int n = C.iscm() ? C.rowsize() : C.colsize();
    int lda = A.stepj();
    int ldb = B.iscm()?B.stepj():B.stepi();
    int ldc = C.iscm()?C.stepj():C.stepi();
    double xbeta(beta);
    BLASNAME(dsymm) (BLASCM C.isrm()?BLASCH_R:BLASCH_L,
        A.uplo() == Upper ? BLASCH_UP : BLASCH_LO,
        BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
        BLASP(B.cptr()),BLASV(ldb),BLASV(xbeta),
        BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
  }
  template <> void BlasMultMM(
      std::complex<double> alpha,
      const GenSymMatrix<std::complex<double> >& A,
      const GenMatrix<std::complex<double> >& B,
      const int beta, const MatrixView<std::complex<double> >& C)
  {
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != 0.);
    TMVAssert(A.ct()==NonConj);
    TMVAssert(B.ct()==C.ct());
    TMVAssert(A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.isrm() || C.iscm());
    TMVAssert(B.stor() == C.stor());
    TMVAssert(A.issym() || C.iscm()!=C.isconj());
    TMVAssert(A.isherm() || !C.isconj());

    int m = C.iscm() ? C.colsize() : C.rowsize();
    int n = C.iscm() ? C.rowsize() : C.colsize();
    int lda = A.stepj();
    int ldb = B.iscm()?B.stepj():B.stepi();
    int ldc = C.iscm()?C.stepj():C.stepi();
    std::complex<double> xbeta(beta);
    if (A.issym())
      BLASNAME(zsymm) (BLASCM C.isrm()?BLASCH_R:BLASCH_L,
          A.uplo() == Upper ? BLASCH_UP : BLASCH_LO,
          BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
          BLASP(B.cptr()),BLASV(ldb),BLASP(&xbeta),
          BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
    else {
      if (C.isrm()) alpha = CONJ(alpha);
      BLASNAME(zhemm) (BLASCM C.isrm()?BLASCH_R:BLASCH_L,
          A.uplo() == Upper ? BLASCH_UP : BLASCH_LO,
          BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
          BLASP(B.cptr()),BLASV(ldb),BLASP(&xbeta),
          BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
    }
  }
#endif
#ifdef INST_FLOAT
  template <> void BlasMultMM(
      const float alpha, const GenSymMatrix<float>& A,
      const GenMatrix<float>& B, const int beta, const MatrixView<float>& C)
  {
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != 0.F);
    TMVAssert(A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.isrm() || C.iscm());
    TMVAssert(B.stor() == C.stor());

    int m = C.iscm() ? C.colsize() : C.rowsize();
    int n = C.iscm() ? C.rowsize() : C.colsize();
    int lda = A.stepj();
    int ldb = B.iscm()?B.stepj():B.stepi();
    int ldc = C.iscm()?C.stepj():C.stepi();
    float xbeta(beta);
    BLASNAME(ssymm) (BLASCM C.isrm()?BLASCH_R:BLASCH_L,
        A.uplo() == Upper ? BLASCH_UP : BLASCH_LO,
        BLASV(m),BLASV(n),BLASV(alpha),BLASP(A.cptr()),BLASV(lda),
        BLASP(B.cptr()),BLASV(ldb),BLASV(xbeta),
        BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
  }
  template <> void BlasMultMM(
      std::complex<float> alpha,
      const GenSymMatrix<std::complex<float> >& A,
      const GenMatrix<std::complex<float> >& B,
      const int beta, const MatrixView<std::complex<float> >& C)
  {
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != 0.F);
    TMVAssert(A.ct()==NonConj);
    TMVAssert(B.ct()==C.ct());
    TMVAssert(A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.isrm() || C.iscm());
    TMVAssert(B.stor() == C.stor());
    TMVAssert(A.issym() || C.iscm()!=C.isconj());
    TMVAssert(A.isherm() || !C.isconj());

    int m = C.iscm() ? C.colsize() : C.rowsize();
    int n = C.iscm() ? C.rowsize() : C.colsize();
    int lda = A.stepj();
    int ldb = B.iscm()?B.stepj():B.stepi();
    int ldc = C.iscm()?C.stepj():C.stepi();
    std::complex<float> xbeta(beta);
    if (A.issym())
      BLASNAME(csymm) (BLASCM C.isrm()?BLASCH_R:BLASCH_L,
          A.uplo() == Upper ? BLASCH_UP : BLASCH_LO,
          BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
          BLASP(B.cptr()),BLASV(ldb),BLASP(&xbeta),
          BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
    else {
      if (C.isrm()) alpha = CONJ(alpha);
      BLASNAME(chemm) (BLASCM C.isrm()?BLASCH_R:BLASCH_L,
          A.uplo() == Upper ? BLASCH_UP : BLASCH_LO,
          BLASV(m),BLASV(n),BLASP(&alpha),BLASP(A.cptr()),BLASV(lda),
          BLASP(B.cptr()),BLASV(ldb),BLASP(&xbeta),
          BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
    }
  }
#endif
  template <class T> static void BlasMultMM(
      const std::complex<T> alpha,
      const GenSymMatrix<std::complex<T> >& A, const GenMatrix<T>& B,
      const int beta, const MatrixView<std::complex<T> >& C)
  {
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(A.ct()==NonConj);
    TMVAssert(A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.isrm() || C.iscm());
    TMVAssert(B.stor() == C.stor());
    TMVAssert(A.issym() || C.iscm()!=C.isconj());
    TMVAssert(A.isherm() || !C.isconj());

    if (IMAG(alpha) == T(0)) {
      SymMatrix<T,Lower,ColMajor> A1 = A.Real();
      Matrix<T,ColMajor> C1 = REAL(alpha)*A1*B;
      if (beta == 0) C.Real() = C1;
      else C.Real() += C1;
      if (A.issym()) {
        A1 = A.Imag();
        if (C.isconj()) C1 = -REAL(alpha)*A1*B;
        else C1 = REAL(alpha)*A1*B;
      } else {
        LowerTriMatrixView<T> L = A1.LowerTri();
        L = A.LowerTri().Imag();
        // A.Imag() = L - LT
        if (A.LowerTri().isconj() != C.isconj()) {
          C1 = -REAL(alpha)*L*B;
          C1 += REAL(alpha)*L.Transpose()*B;
        } else {
          C1 = REAL(alpha)*L*B;
          C1 -= REAL(alpha)*L.Transpose()*B;
        }
      }
      if (beta == 0) C.Imag() = C1;
      else C.Imag() += C1;
    } else {
      SymMatrix<T,Lower,ColMajor> Ar = A.Real();
      SymMatrix<T,Lower,ColMajor> Ai(A.size());
      LowerTriMatrixView<T> L = Ai.LowerTri();
      Matrix<T,ColMajor> C1 = REAL(alpha)*Ar*B;
      if (A.issym()) {
        Ai = A.Imag();
        C1 -= IMAG(alpha)*Ai*B;
      } else {
        L = A.LowerTri().Imag();
        if (A.LowerTri().isconj()) L *= T(-1);
        C1 -= IMAG(alpha)*L*B;
        C1 += IMAG(alpha)*L.Transpose()*B;
      }
      if (beta == 0) C.Real() = C1;
      else C.Real() += C1;
      C1 = IMAG(alpha)*Ar*B;
      if (A.issym()) {
        C1 += REAL(alpha)*Ai*B;
      } else {
        C1 += REAL(alpha)*L*B;
        C1 -= REAL(alpha)*L.Transpose()*B;
      }
      if (C.isconj()) C1 *= T(-1);
      if (beta == 0) C.Imag() = C1;
      else C.Imag() += C1;
    }
  }
  template <class T> static void BlasMultMM(
      const std::complex<T> alpha,
      const GenSymMatrix<T>& A, const GenMatrix<std::complex<T> >& B,
      const int beta, const MatrixView<std::complex<T> >& C)
  {
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(B.ct()==C.ct());
    TMVAssert(A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.isrm() || C.iscm());
    TMVAssert(B.stor() == C.stor());
    TMVAssert(A.issym() || C.iscm()!=C.isconj());
    TMVAssert(A.isherm() || !C.isconj());

    if (IMAG(alpha) == T(0)) {
      Matrix<T,ColMajor> B1 = B.Real();
      Matrix<T,ColMajor> C1 = REAL(alpha)*A*B1;
      if (beta == 0) C.Real() = C1;
      else C.Real() += C1;
      B1 = B.Imag();
      if (B.isconj()) C1 = -REAL(alpha)*A*B1;
      else C1 = REAL(alpha)*A*B1;
      if (beta == 0) C.Imag() = C1;
      else C.Imag() += C1;
    } else {
      Matrix<T,ColMajor> Br = B.Real();
      Matrix<T,ColMajor> Bi = B.Imag();
      Matrix<T,ColMajor> C1 = REAL(alpha)*A*Br;
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
  template <class T> static void BlasMultMM(
      const std::complex<T> alpha,
      const GenSymMatrix<T>& A, const GenMatrix<T>& B,
      const int beta, const MatrixView<std::complex<T> >& C)
  {
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));
    TMVAssert(A.iscm());
    TMVAssert(B.isrm() || B.iscm());
    TMVAssert(C.isrm() || C.iscm());
    TMVAssert(B.stor() == C.stor());
    TMVAssert(A.issym() || C.iscm()!=C.isconj());
    TMVAssert(A.isherm() || !C.isconj());

    Matrix<T,ColMajor> C1 = A*B;
    if (beta == 0) C = alpha*C1;
    else C += alpha*C1;
  }
#endif // BLAS

  template <bool add, class T, class Ta, class Tb> static void DoMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.rowsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(alpha != T(0));

#ifdef BLAS
    if (A.isrm())
      DoMultMM<add>(alpha,A.issym()?A.Transpose():A.Adjoint(),B,C);
    else if (A.isconj())
      DoMultMM<add>(CONJ(alpha),A.Conjugate(),B.Conjugate(),C.Conjugate());
    else if (!((C.isrm() && C.stepi()>0) || (C.iscm() && C.stepj()>0)) ||
        (C.iscm() && C.isconj()) || (C.isrm() && C.isconj()==A.issym())) {
      Matrix<T,ColMajor> C2(C.colsize(),C.rowsize());
      DoMultMM<false>(T(1),A,B,C2.View());
      if (add) C += alpha*C2;
      else C = alpha*C2;
    } else if (!(A.iscm() && A.stepj()>0)) {
      if (IMAG(alpha) == RealType(T)(0)) {
        if (A.isherm()) {
          if (A.uplo() == Upper) {
            HermMatrix<Ta,Upper,ColMajor> A2 = REAL(alpha)*A;
            DoMultMM<add>(T(1),A2,B,C);
          } else {
            HermMatrix<Ta,Lower,ColMajor> A2 = REAL(alpha)*A;
            DoMultMM<add>(T(1),A2,B,C);
          }
        } else {
          if (A.uplo() == Upper) {
            SymMatrix<Ta,Upper,ColMajor> A2 = REAL(alpha)*A;
            DoMultMM<add>(T(1),A2,B,C);
          } else {
            SymMatrix<Ta,Lower,ColMajor> A2 = REAL(alpha)*A;
            DoMultMM<add>(T(1),A2,B,C);
          }
        }
      } else {
        if (!A.issym()) {
          if (A.uplo() == Upper) {
            // alpha * A is not Hermitian, so can't do A2 = alpha * A
            HermMatrix<Ta,Upper,ColMajor> A2 = A;
            DoMultMM<add>(alpha,A2,B,C);
          } else {
            HermMatrix<Ta,Lower,ColMajor> A2 = A;
            DoMultMM<add>(alpha,A2,B,C);
          }
        } else {
          if (A.uplo() == Upper) {
            SymMatrix<T,Upper,ColMajor> A2 = alpha*A;
            DoMultMM<add>(T(1),A2,B,C);
          } else {
            SymMatrix<T,Lower,ColMajor> A2 = alpha*A;
            DoMultMM<add>(T(1),A2,B,C);
          }
        }
      }
    } else if ((B.stor() != C.stor()) || 
        (IsComplex(Tb()) && B.isconj() != C.isconj()) || 
        !((B.isrm() && B.stepi()>0) || (B.iscm() && B.stepj()>0))) {
      if (IMAG(alpha) == RealType(T)(0)) {
        if (C.isconj()) {
          if (C.iscm()) {
            Matrix<Tb,ColMajor> B2 = REAL(alpha)*B.Conjugate();
            DoMultMM<add>(T(1),A,B2.Conjugate(),C);
          } else {
            Matrix<Tb,RowMajor> B2 = REAL(alpha)*B.Conjugate();
            DoMultMM<add>(T(1),A,B2.Conjugate(),C);
          }
        } else {
          if (C.iscm()) {
            Matrix<Tb,ColMajor> B2 = REAL(alpha)*B;
            DoMultMM<add>(T(1),A,B2,C);
          } else {
            Matrix<Tb,RowMajor> B2 = REAL(alpha)*B;
            DoMultMM<add>(T(1),A,B2,C);
          }
        }
      } else {
        if (C.isconj()) {
          if (C.iscm()) {
            Matrix<T,ColMajor> B2 = CONJ(alpha)*B.Conjugate();
            DoMultMM<add>(T(1),A,B2.Conjugate(),C);
          } else {
            Matrix<T,RowMajor> B2 = CONJ(alpha)*B.Conjugate();
            DoMultMM<add>(T(1),A,B2.Conjugate(),C);
          }
        } else {
          if (C.iscm()) {
            Matrix<T,ColMajor> B2 = alpha*B;
            DoMultMM<add>(T(1),A,B2,C);
          } else {
            Matrix<T,RowMajor> B2 = alpha*B;
            DoMultMM<add>(T(1),A,B2,C);
          }
        }
      }
    } else {
      BlasMultMM(alpha,A,B,add?1:0,C);
    }
#else
    NonBlasMultMM<add>(alpha,A,B,C);
#endif
  }

  template <bool add, class T, class Ta, class Tb> static void FullTempMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
    if (C.isrm()) {
      Matrix<T,RowMajor> C2(C.colsize(),C.rowsize());
      DoMultMM<false>(T(1),A,B,C2.View());
      if (add) C += alpha*C2;
      else C = alpha*C2;
    } else {
      Matrix<T,ColMajor> C2(C.colsize(),C.rowsize());
      DoMultMM<false>(T(1),A,B,C2.View());
      if (add) C += alpha*C2;
      else C = alpha*C2;
    }
  }

  template <bool add, class T, class Ta, class Tb> static void BlockTempMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
    const int N = C.rowsize();
    for(int j=0;j<N;) {
      int j2 = MIN(N,j+SYM_MM_BLOCKSIZE);
      if (IMAG(alpha) == RealType(T)(0)) {
        if (C.isrm()) {
          Matrix<Tb,RowMajor> B2 = REAL(alpha) * B.Cols(j,j2);
          DoMultMM<add>(T(1),A,B2,C.Cols(j,j2));
        } else {
          Matrix<Tb,ColMajor> B2 = REAL(alpha) * B.Cols(j,j2);
          DoMultMM<add>(T(1),A,B2,C.Cols(j,j2));
        }
      } else {
        if (C.isrm()) {
          Matrix<T,RowMajor> B2 = alpha * B.Cols(j,j2);
          DoMultMM<add>(T(1),A,B2,C.Cols(j,j2));
        } else {
          Matrix<T,ColMajor> B2 = alpha * B.Cols(j,j2);
          DoMultMM<add>(T(1),A,B2,C.Cols(j,j2));
        }
      }
      j = j2;
    }
  }

  template <bool add, class T, class Ta, class Tb> void MultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  // C (+)= alpha * A * B
  {
#ifdef XTEST
    TMVAssert(A.HermOK());
#endif
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
#ifdef XDEBUG
    //cout<<"Start MultMM: alpha = "<<alpha<<endl;
    //cout<<"A = "<<A.cptr()<<"  "<<TypeText(A)<<"  "<<A<<endl;
    //cout<<"B = "<<B.cptr()<<"  "<<TypeText(B)<<"  "<<B<<endl;
    //cout<<"C = "<<C.cptr()<<"  "<<TypeText(C)<<"  "<<C<<endl;
    Matrix<Ta> A0 = A;
    Matrix<Tb> B0 = B;
    Matrix<T> C0 = C;
    Matrix<T> C2 = alpha*A0*B0;
    if (add) C2 += C0;
#endif

    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (alpha == T(0)) {
        if (!add) C.Zero();
      }
      else if (SameStorage(A,C)) 
        FullTempMultMM<add>(alpha,A,B,C);
      else if (SameStorage(B,C)) 
        if (C.stepi() == B.stepi() && C.stepj() == B.stepj())
          BlockTempMultMM<add>(alpha,A,B,C);
        else
          FullTempMultMM<add>(alpha,A,B,C);
      else DoMultMM<add>(alpha, A, B, C);
    }

#ifdef XDEBUG
    //cout<<"Done: C = "<<C<<endl;
    if (Norm(C-C2) > 0.001*(ABS(alpha)*Norm(A0)*Norm(B0)+
          (add?Norm(C0):RealType(T)(0)))) {
      cerr<<"MultMM: alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"A = "<<TypeText(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<TypeText(B)<<"  "<<B0<<endl;
      cerr<<"C = "<<TypeText(C)<<"  "<<C0<<endl;
      cerr<<"--> C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

  template <bool add, class T, class Ta, class Tb> static void BlockTempMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenSymMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == C.rowsize());
    TMVAssert(A.size() > 0);
    TMVAssert(alpha != T(0));

    const int N = A.size();

    for(int j=0;j<N;) {
      int j2 = MIN(N,j+SYM_MM_BLOCKSIZE);
      if (IMAG(alpha) == RealType(T)(0)) {
        if (C.isrm()) {
          Matrix<Tb,RowMajor> B2(N,j2-j);
          B2.Rows(0,j) = REAL(alpha) * B.SubMatrix(0,j,j,j2);
          B2.Rows(j,j2) = REAL(alpha) * B.SubSymMatrix(j,j2);
          B2.Rows(j2,N) = REAL(alpha) * B.SubMatrix(j2,N,j,j2);
          DoMultMM<add>(T(1),A,B2.View(),C.Cols(j,j2));
        } else {
          Matrix<Tb,ColMajor> B2(N,j2-j);
          B2.Rows(0,j) = REAL(alpha) * B.SubMatrix(0,j,j,j2);
          B2.Rows(j,j2) = REAL(alpha) * B.SubSymMatrix(j,j2);
          B2.Rows(j2,N) = REAL(alpha) * B.SubMatrix(j2,N,j,j2);
          DoMultMM<add>(T(1),A,B2.View(),C.Cols(j,j2));
        }
      } else {
        if (C.isrm()) {
          Matrix<T,RowMajor> B2(N,j2-j);
          B2.Rows(0,j) = alpha * B.SubMatrix(0,j,j,j2);
          B2.Rows(j,j2) = alpha * B.SubSymMatrix(j,j2);
          B2.Rows(j2,N) = alpha * B.SubMatrix(j2,N,j,j2);
          DoMultMM<add>(T(1),A,B2.View(),C.Cols(j,j2));
        } else {
          Matrix<T,ColMajor> B2(N,j2-j);
          B2.Rows(0,j) = alpha * B.SubMatrix(0,j,j,j2);
          B2.Rows(j,j2) = alpha * B.SubSymMatrix(j,j2);
          B2.Rows(j2,N) = alpha * B.SubMatrix(j2,N,j,j2);
          DoMultMM<add>(T(1),A,B2.View(),C.Cols(j,j2));
        }
      }
      j = j2;
    }
  }

  template <bool add, class T, class Ta, class Tb> static void FullTempMultMM(
      const T alpha, const GenSymMatrix<Ta>& A, const GenSymMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
    if (C.isrm()) {
      Matrix<T,RowMajor> C2(C.colsize(),C.rowsize());
      BlockTempMultMM<false>(T(1),A,B,C2.View());
      if (add) C += alpha*C2;
      else C = alpha*C2;
    } else {
      Matrix<T,ColMajor> C2(C.colsize(),C.rowsize());
      BlockTempMultMM<false>(T(1),A,B,C2.View());
      if (add) C += alpha*C2;
      else C = alpha*C2;
    }
  }

  template <bool add, class T, class Ta, class Tb> void MultMM(const T alpha,
      const GenSymMatrix<Ta>& A, const GenSymMatrix<Tb>& B,
      const MatrixView<T>& C)
  // C (+)= alpha * A * B
  {
#ifdef XTEST
    TMVAssert(A.HermOK());
    TMVAssert(B.HermOK());
#endif
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == C.rowsize());
#ifdef XDEBUG
    //cout<<"Start MultMM: alpha = "<<alpha<<endl;
    //cout<<"A = "<<A.cptr()<<"  "<<TypeText(A)<<"  "<<A<<endl;
    //cout<<"B = "<<B.cptr()<<"  "<<TypeText(B)<<"  "<<B<<endl;
    //cout<<"C = "<<C.cptr()<<"  "<<TypeText(C)<<"  "<<C<<endl;
    Matrix<Ta> A0 = A;
    Matrix<Tb> B0 = B;
    Matrix<T> C0 = C;
    Matrix<T> C2 = alpha*A0*B0;
    if (add) C2 += C0;
#endif

    if (A.size() > 0) {
      if (SameStorage(A,C) || SameStorage(B,C))
        FullTempMultMM<add>(alpha,A,B,C);
      else BlockTempMultMM<add>(alpha, A, B, C);
    }

#ifdef XDEBUG
    //cout<<"Done: C = "<<C<<endl;
    if (Norm(C-C2) > 0.001*(ABS(alpha)*Norm(A0)*Norm(B0)+
          (add?Norm(C0):RealType(T)(0)))) {
      cerr<<"MultMM: alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"A = "<<TypeText(A)<<"  "<<A0<<endl;
      cerr<<"B = "<<TypeText(B)<<"  "<<B0<<endl;
      cerr<<"C = "<<TypeText(C)<<"  "<<C0<<endl;
      cerr<<"--> C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_MultSM.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


