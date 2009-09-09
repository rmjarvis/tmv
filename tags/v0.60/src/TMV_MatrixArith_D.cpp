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
#include "TMV_VectorArith.h"
#include "TMV_MatrixArith.h"

//#define XDEBUG

#ifdef XDEBUG
#include <iostream>
#include "TMV_VectorArith.h"
using std::cerr;
using std::endl;
#endif

namespace tmv {

#ifdef TMV_BLOCKSIZE
  const size_t MM_BLOCKSIZE = TMV_BLOCKSIZE;
#else
  const size_t MM_BLOCKSIZE = 64;
#endif

  // MJ: Look at Atlas code, and try to mimic structure to make this faster.
 
  //
  // MultMM
  //

  template <bool add, class T, class Ta, class Tb> inline void RowMultMM(
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

    if (add)
      for(size_t i=0;i<C.colsize();++i) 
	C.row(i) += alpha * A.row(i) * B;
    else
      for(size_t i=0;i<C.colsize();++i) 
	C.row(i) = alpha * A.row(i) * B;
  }

  template <bool add, class T, class Ta, class Tb> inline void OPMultMM(
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

    if (!add) C.Zero();
    for(size_t k=0;k<A.rowsize();++k) 
      C += alpha * A.col(k) ^ B.row(k);
  }

  template <bool add, class T, class Ta, class Tb> inline void NonBlasMultMM(
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
    TMVAssert(!C.isrm());
    // We want to make the inner loops as efficient as possible.
    // We want the inner loops to be unit stride or with long vectors.

    if (B.iscm() && C.iscm()) 
      RowMultMM<add>(alpha,B.Transpose(),A.Transpose(),C.Transpose());
    else if (A.iscm() && B.isrm()) OPMultMM<add>(alpha,A,B,C);
    else {
      const size_t M = C.colsize();
      const size_t N = C.rowsize();
      const size_t K = A.rowsize();
      if (M < N && M < K) RowMultMM<add>(alpha,A,B,C);
      else if (N < M && N < K) 
	RowMultMM<add>(alpha,B.Transpose(),A.Transpose(),C.Transpose());
      else if (K < M && K < N) OPMultMM<add>(alpha,A,B,C);
      else if (M < N) RowMultMM<add>(alpha,A,B,C);
      else if (A.isrm()) RowMultMM<add>(alpha,A,B,C);
      else if (B.iscm() || C.iscm()) 
	RowMultMM<add>(alpha,B.Transpose(),A.Transpose(),C.Transpose());
      else if (A.iscm() || B.isrm()) OPMultMM<add>(alpha,A,B,C);
      else RowMultMM<add>(alpha,A,B,C);
    }
  }

#ifdef BLAS
  template <class T, class Ta, class Tb> inline void BlasMultMM(
      const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const int beta, const MatrixView<T>& C)
  {
    if (beta == 0) NonBlasMultMM<false>(alpha,A,B,C); 
    else NonBlasMultMM<true>(alpha,A,B,C); 
  }
#ifdef INST_DOUBLE
  template <> inline void BlasMultMM(
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
  template <> inline void BlasMultMM(
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
      Matrix<std::complex<double> > AA = alpha*A;
      return BlasMultMM(std::complex<double>(1),AA,B,beta,C);
    } else if (B.iscm() && B.isconj()) {
      Matrix<std::complex<double> > BB = alpha*B;
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
#endif
#ifdef INST_FLOAT
  template <> inline void BlasMultMM(
      const float alpha, const GenMatrix<float>& A,
      const GenMatrix<float>& B, const int beta, const MatrixView<float>& C)
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
    float xbeta(beta);
    BLASNAME(sgemm) (BLASCM A.iscm()?BLASCH_NT:BLASCH_T,
	B.iscm()?BLASCH_NT:BLASCH_T,
	BLASV(m),BLASV(n),BLASV(k),BLASV(alpha),
	BLASP(A.cptr()),BLASV(lda),BLASP(B.cptr()),BLASV(ldb),
	BLASV(xbeta),BLASP(C.ptr()),BLASV(ldc) BLAS1 BLAS1);
  }
  template <> inline void BlasMultMM(
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
#endif 
#ifdef ELAP
#ifdef INST_DOUBLE
  template <> inline void BlasMultMM(
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
    TMVAssert(C.iscm());

    if (beta == 0 && A.iscm() && B.iscm() && B.IsSquare() && !A.isconj()) {
      int m = C.colsize();
      int n = C.rowsize();
      int lda = A.stepj();
      int ldb = B.stepj();
      int ldc = C.stepj();
#ifndef LAPNOWORK
      int lwork = 2*m*n;
      double* rwork = LAP_DWork(lwork);
#endif
      LAPNAMEX(zlacrm) (LAPCM LAPV(m),LAPV(n),LAPP(A.cptr()),LAPV(lda),
	  LAPP(B.cptr()),LAPV(ldb),LAPP(C.ptr()),LAPV(ldc) LAPWK(rwork));
      C *= alpha;
    } else 
      if (beta == 0) NonBlasMultMM<false>(alpha,A,B,C);
      else NonBlasMultMM<true>(alpha,A,B,C);
  }
#endif
#ifdef INST_FLOAT
  template <> inline void BlasMultMM(
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
    TMVAssert(C.iscm());

    if (beta == 0 && A.iscm() && B.iscm() && B.IsSquare() && !A.isconj()) {
      int m = C.colsize();
      int n = C.rowsize();
      int lda = A.stepj();
      int ldb = B.stepj();
      int ldc = C.stepj();
#ifndef LAPNOWORK
      int lwork = 2*m*n;
      float* rwork = LAP_SWork(lwork);
#endif
      LAPNAMEX(clacrm) (LAPCM LAPV(m),LAPV(n),LAPP(A.cptr()),LAPV(lda),
	  LAPP(B.cptr()),LAPV(ldb),LAPP(C.ptr()),LAPV(ldc) LAPWK(rwork));
      C *= alpha;
    } else 
      if (beta == 0) NonBlasMultMM<false>(alpha,A,B,C);
      else NonBlasMultMM<true>(alpha,A,B,C);
  }
#endif 
#endif // ELAP
#endif // BLAS

  template <bool add, class T, class Ta, class Tb> inline void DoMultMM(
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
    TMVAssert(!C.isrm());

#ifdef BLAS
    if (IsComplex(T()) && (IsReal(Ta()) || IsReal(Tb())))
      BlasMultMM(alpha,A,B,add?1:0,C);
    else if (!(C.iscm() && C.stepj()>0)) {
      Matrix<T,ColMajor> C2(C.colsize(),C.rowsize(),T(0));
      DoMultMM<false>(alpha,A,B,C2.View());
      if (add) C += C2;
      else C = C2;
    } else if (!((A.isrm() && A.stepi()>0) || (A.iscm() && A.stepj()>0))) {
      if (IMAG(alpha) == RealType(T)(0)) {
	Matrix<Ta,ColMajor> A2 = REAL(alpha)*A;
	DoMultMM<add>(T(1),A2,B,C);
      } else {
	Matrix<T,ColMajor> A2 = alpha*A;
	DoMultMM<add>(T(1),A2,B,C);
      }
    } else if (!((B.isrm() && B.stepi()>0) || (B.iscm() && B.stepj()>0))) {
      if (IMAG(alpha) == RealType(T)(0)) {
	Matrix<Tb,ColMajor> B2 = REAL(alpha)*B;
	DoMultMM<add>(T(1),A,B2,C);
      } else {
	Matrix<T,ColMajor> B2 = alpha*B;
	DoMultMM<add>(T(1),A,B2,C);
      }
    } else {
      BlasMultMM(alpha,A,B,add?1:0,C);
    }
#else
    NonBlasMultMM<add>(alpha,A,B,C);
#endif
  }

  template <bool add, class T, class Ta, class Tb> inline void FullTempMultMM(
      const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
    // C (+)= alpha * A * B  via a temporary
  {
    Matrix<T,ColMajor> C2(C.colsize(),C.rowsize());
    DoMultMM<false>(T(1),A,B,C2.View());

    if (add) C += alpha*C2;
    else C = alpha*C2;
  }

  template <bool add, class T, class Ta, class Tb> inline void BlockTempMultMM(
      const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
    // C (+)= alpha * A * C
  {
    for (size_t j=0;j<C.rowsize();) {
      size_t j2 = std::min(C.rowsize(),j+MM_BLOCKSIZE);
      if (C.isrm()) {
	if (IMAG(alpha) == RealType(T)(0)) {
	  Matrix<Tb,ColMajor> B2 = REAL(alpha) * B.Cols(j,j2);
	  DoMultMM<add>(T(1),B2.Transpose(),A.Transpose(),
	      C.Cols(j,j2).Transpose());
	} else {
	  Matrix<T,ColMajor> B2 = alpha * B.Cols(j,j2);
	  DoMultMM<add>(T(1),B2.Transpose(),A.Transpose(),
	      C.Cols(j,j2).Transpose());
	}
      } else {
	if (IMAG(alpha) == RealType(T)(0)) {
	  Matrix<Tb,ColMajor> B2 = REAL(alpha) * B.Cols(j,j2);
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
    //cerr<<"MultMM: add = "<<add<<", alpha = "<<alpha<<endl;
    //cerr<<"A = "<<Type(A)<<" "<<A0<<endl;
    //cerr<<"B = "<<Type(B)<<" "<<B0<<endl;
    //cerr<<"C = "<<Type(C)<<" "<<C0<<endl;
#endif

    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (A.rowsize() == 0 || alpha == T(0))  {
	if (!add) C.Zero();
      }
      else if (C.isconj()) 
	MultMM<add>(CONJ(alpha),A.Conjugate(),B.Conjugate(),C.Conjugate());
      else if (C.isrm()) 
	MultMM<add>(alpha,B.Transpose(),A.Transpose(),C.Transpose());
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

#define InstFile "TMV_MatrixArith_D.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


