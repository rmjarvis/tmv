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



#include "TMV_DiagTriArithFunc.h"
#include "TMV_DiagMatrix.h"
#include "TMV_TriMatrix.h"
#include "TMV_TriMatrixArith.h"
#include "TMV_DiagMatrixArith.h"

//#define XDEBUG

#ifdef XDEBUG
#include "TMV_MatrixArith.h"
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#endif

namespace tmv {

  template <bool a1, bool ca, class T, class T1,  class Ta> 
    static void DoMultEqMM(const T1 alpha, const GenDiagMatrix<Ta>& A,
	const UpperTriMatrixView<T>& B)
    {
      TMVAssert(A.size() == B.size());
      TMVAssert(A.size() > 0);
      TMVAssert(B.ct()==NonConj);
      TMVAssert(ca == A.diag().isconj());
      TMVAssert(B.dt() == NonUnitDiag);
      TMVAssert(a1 == (alpha == T(1)));

      if (A.size() == 1) {
	const Ta Ax = ca ? CONJ(*A.diag().cptr()) : *A.diag().cptr();
#ifdef TMVFLDEBUG
	TMVAssert(B.ptr() >= B.first);
	TMVAssert(B.ptr() < B.last);
#endif
	if (a1) *B.ptr() *= Ax;
	else *B.ptr() *= alpha * Ax;
      } else {
	const int N = A.size();
	const int k = N/2;
	// [ B00 B01 ] = [ A00  0  ] * [ B00 B01 ]
	// [  0  B11 ]   [  0  A11 ]   [  0  B11 ]
	// B00 = A00 * B00
	// B01 = A00 * B01
	// B11 = A11 * B11
	ConstDiagMatrixView<Ta> A00 = A.SubDiagMatrix(0,k);
	ConstDiagMatrixView<Ta> A11 = A.SubDiagMatrix(k,N);
	UpperTriMatrixView<T> B00 = B.SubTriMatrix(0,k);
	UpperTriMatrixView<T> B11 = B.SubTriMatrix(k,N);
	MatrixView<T> B01 = B.SubMatrix(0,k,k,N);

	DoMultEqMM<a1,ca>(alpha,A00,B00);
	B01 = alpha * A00 * B01;
	DoMultEqMM<a1,ca>(alpha,A11,B11);
      }
    }

  template <class T, class Ta> static void MultEqMM(const T alpha,
      const GenDiagMatrix<Ta>& A, const UpperTriMatrixView<T>& B)
    // B = alpha * A * B
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(alpha != T(0));
    TMVAssert(B.dt() == NonUnitDiag);
#ifdef XDEBUG
    //cout<<"Start MultEqMM: alpha = "<<alpha<<endl;
    //cout<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cout<<"B = "<<Type(B)<<"  "<<B<<endl;
    Matrix<T> B0 = B;
    Matrix<Ta> A0 = A;
    Matrix<T> B2 = alpha*A0*B0;
#endif

    if (B.isconj()) MultEqMM(CONJ(alpha),A.Conjugate(),B.Conjugate());
    else if (A.size() > 0) {
      if (alpha == T(1)) 
	if (A.diag().isconj())
	  DoMultEqMM<true,true>(REAL(alpha),A,B);
	else
	  DoMultEqMM<true,false>(REAL(alpha),A,B);
      else if (IMAG(alpha) == RealType(T)(0)) 
	if (A.diag().isconj())
	  DoMultEqMM<false,true>(REAL(alpha),A,B);
	else
	  DoMultEqMM<false,false>(REAL(alpha),A,B);
      else 
	if (A.diag().isconj())
	  DoMultEqMM<false,true>(alpha,A,B);
	else
	  DoMultEqMM<false,false>(alpha,A,B);
    }

#ifdef XDEBUG
    if (Norm(Matrix<T>(B)-B2) > 0.001*(ABS(alpha)*Norm(A0)*Norm(B0))) {
      cerr<<"MultEqMM: alpha = "<<alpha<<endl;
      cerr<<"A = "<<Type(A)<<" step "<<A.diag().step()<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"-> B = "<<B<<endl;
      cerr<<"B2 = "<<B2<<endl;
      abort();
    }
#endif
  }

  template <bool a1, bool ca, bool ub, bool cb, class T, class T1, class Ta, class Tb> 
    static void DoAddMultMM(const T1 alpha,
	const GenDiagMatrix<Ta>& A, const GenUpperTriMatrix<Tb>& B,
	const UpperTriMatrixView<T>& C)
    {
      TMVAssert(A.size() == B.size());
      TMVAssert(A.size() == C.size());
      TMVAssert(A.size() > 0);
      TMVAssert(C.ct()==NonConj);
      TMVAssert(ca == A.diag().isconj());
      TMVAssert(cb == B.isconj());
      TMVAssert(C.dt() == NonUnitDiag);
      TMVAssert(ub == B.isunit());
      TMVAssert(a1 == (alpha == T(1)));

      if (A.size() == 1) {
	const Ta Ax = ca ? CONJ(*A.diag().cptr()) : *A.diag().cptr();
#ifdef TMVFLDEBUG
	TMVAssert(C.ptr() >= C.first);
	TMVAssert(C.ptr() < C.last);
#endif
	if (a1)
	  if (ub) *C.ptr() += Ax;
	  else *C.ptr() += Ax * (cb?CONJ(*B.cptr()):*B.cptr());
	else
	  if (ub) *C.ptr() += alpha * Ax;
	  else *C.ptr() += alpha * Ax * (cb?CONJ(*B.cptr()):*B.cptr());
      } else {
	const int N = A.size();
	const int k = N/2;
	// [ C00 C01 ] += [ A00  0  ] * [ B00 B01 ]
	// [  0  C11 ]    [  0  A11 ]   [  0  B11 ]
	// C00 += A00 * B00
	// C01 += A00 * B01
	// C11 += A11 * B11
	ConstDiagMatrixView<Ta> A00 = A.SubDiagMatrix(0,k);
	ConstDiagMatrixView<Ta> A11 = A.SubDiagMatrix(k,N);
	ConstUpperTriMatrixView<Tb> B00 = B.SubTriMatrix(0,k);
	ConstUpperTriMatrixView<Tb> B11 = B.SubTriMatrix(k,N);
	ConstMatrixView<Tb> B01 = B.SubMatrix(0,k,k,N);
	UpperTriMatrixView<T> C00 = C.SubTriMatrix(0,k);
	UpperTriMatrixView<T> C11 = C.SubTriMatrix(k,N);
	MatrixView<T> C01 = C.SubMatrix(0,k,k,N);

	DoAddMultMM<a1,ca,ub,cb>(alpha,A00,B00,C00);
	C01 += alpha * A00 * B01;
	DoAddMultMM<a1,ca,ub,cb>(alpha,A11,B11,C11);
      }
    } 

  template <bool a1, bool ca, class T, class T1, class Ta, class Tb> 
    static inline void DoAddMultMMa(const T1 alpha,
	const GenDiagMatrix<Ta>& A, const GenUpperTriMatrix<Tb>& B,
	const UpperTriMatrixView<T>& C)
    {
      if (B.isunit())
	if (B.isconj()) DoAddMultMM<a1,ca,true,true>(alpha,A,B,C);
	else DoAddMultMM<a1,ca,false,true>(alpha,A,B,C);
      else
	if (B.isconj()) DoAddMultMM<a1,ca,true,false>(alpha,A,B,C);
	else DoAddMultMM<a1,ca,false,false>(alpha,A,B,C);
    }

  template <class T, class Ta, class Tb> static void AddMultMM(const T alpha,
      const GenDiagMatrix<Ta>& A, const GenUpperTriMatrix<Tb>& B,
      const UpperTriMatrixView<T>& C)
    // C += alpha * A * B
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.size());
    TMVAssert(alpha != T(0));
    TMVAssert(C.dt() == NonUnitDiag);
#ifdef XDEBUG
    //cout<<"Start AddMultMM: alpha = "<<alpha<<endl;
    //cout<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cout<<"B = "<<Type(B)<<"  "<<B<<endl;
    //cout<<"C = "<<Type(C)<<"  "<<C<<endl;
    Matrix<Ta> A0 = A;
    Matrix<Tb> B0 = B;
    Matrix<T> C0 = C;
    Matrix<T> C2 = C0+alpha*A0*B0;
#endif

    if (C.isconj()) AddMultMM(CONJ(alpha),A.Conjugate(),
	B.Conjugate(),C.Conjugate());
    else if (A.size() > 0) {
      if (alpha == T(1)) 
	if (A.diag().isconj())
	  DoAddMultMMa<true,true>(REAL(alpha),A,B,C);
	else
	  DoAddMultMMa<true,false>(REAL(alpha),A,B,C);
      else if (IMAG(alpha) == RealType(T)(0)) 
	if (A.diag().isconj())
	  DoAddMultMMa<false,true>(REAL(alpha),A,B,C);
	else
	  DoAddMultMMa<false,false>(REAL(alpha),A,B,C);
      else 
	if (A.diag().isconj())
	  DoAddMultMMa<false,true>(alpha,A,B,C);
	else
	  DoAddMultMMa<false,false>(alpha,A,B,C);
    }

#ifdef XDEBUG
    if (Norm(Matrix<T>(C)-C2) > 0.001*(ABS(alpha)*Norm(A0)*Norm(B0)+Norm(C0))) {
      cerr<<"AddMultMM: alpha = "<<alpha<<endl;
      cerr<<"A = "<<Type(A)<<" step "<<A.diag().step()<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"-> C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

  template <bool add, class T, class Ta, class Tb> void MultMM(const T alpha,
      const GenDiagMatrix<Ta>& A, const GenUpperTriMatrix<Tb>& B,
      const UpperTriMatrixView<T>& C)
    // C (+)= alpha * A * B
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.size());
    TMVAssert(C.dt() == NonUnitDiag);
#ifdef XDEBUG
    //cout<<"Start MultMM: alpha = "<<alpha<<endl;
    //cout<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cout<<"B = "<<Type(B)<<"  "<<B<<endl;
    //cout<<"C = "<<Type(C)<<"  "<<C<<endl;
    Matrix<Ta> A0 = A;
    Matrix<Tb> B0 = B;
    Matrix<T> C0 = C;
    Matrix<T> C2 = alpha*A0*B0;
    if (add) C2 += C0;
#endif

    if (A.size() > 0) {
      if (alpha==T(0)) {
	if (!add) C.Zero();
      }
      else if (SameStorage(A.diag(),C)) {
	DiagMatrix<T> tempA = A;
	MultMM<add>(alpha,tempA,B,C);
      }
      else if (!add) {
	C = B;
	MultEqMM(alpha,A,C);
      }
      else if (SameStorage(B,C)) {
	if (B.isrm()) {
	  UpperTriMatrix<T,NonUnitDiag,RowMajor> tempB = B;
	  MultEqMM(alpha,A,tempB.View());
	  C += tempB;
	}
	else {
	  UpperTriMatrix<T,NonUnitDiag,ColMajor> tempB = B;
	  MultEqMM(alpha,A,tempB.View());
	  C += tempB;
	}
      }
      else {
	AddMultMM(alpha,A,B,C);
      }
    }
#ifdef XDEBUG
    if (Norm(Matrix<T>(C)-C2) > 0.001*(ABS(alpha)*Norm(A0)*Norm(B0)+
	  add?Norm(C0):RealType(T)(0))) {
      cerr<<"MultMM: alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"A = "<<Type(A)<<" step "<<A.diag().step()<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"-> C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

  template <bool a1, bool ca, class T, class T1,  class Ta> 
    static void DoMultEqMM(const T1 alpha, const GenDiagMatrix<Ta>& A,
	const LowerTriMatrixView<T>& B)
    {
      TMVAssert(A.size() == B.size());
      TMVAssert(A.size() > 0);
      TMVAssert(B.ct()==NonConj);
      TMVAssert(ca == A.diag().isconj());
      TMVAssert(B.dt() == NonUnitDiag);
      TMVAssert(a1 == (alpha == T(1)));

      if (A.size() == 1) {
	const Ta Ax = ca ? CONJ(*A.diag().cptr()) : *A.diag().cptr();
#ifdef TMVFLDEBUG
	TMVAssert(B.ptr() >= B.first);
	TMVAssert(B.ptr() < B.last);
#endif
	if (a1) *B.ptr() *= Ax;
	else *B.ptr() *= alpha * Ax;
      } else {
	const int N = A.size();
	const int k = N/2;
	// [ B00  0  ] = [ A00  0  ] * [ B00  0  ]
	// [ B10 B11 ]   [  0  A11 ]   [ B10 B11 ]
	// B00 = A00 * B00
	// B10 = A11 * B10
	// B11 = A11 * B11
	ConstDiagMatrixView<Ta> A00 = A.SubDiagMatrix(0,k);
	ConstDiagMatrixView<Ta> A11 = A.SubDiagMatrix(k,N);
	LowerTriMatrixView<T> B00 = B.SubTriMatrix(0,k);
	LowerTriMatrixView<T> B11 = B.SubTriMatrix(k,N);
	MatrixView<T> B10 = B.SubMatrix(k,N,0,k);

	DoMultEqMM<a1,ca>(alpha,A00,B00);
	B10 = alpha * A11 * B10;
	DoMultEqMM<a1,ca>(alpha,A11,B11);
      }
    }

  template <class T, class Ta> static void MultEqMM(const T alpha,
      const GenDiagMatrix<Ta>& A, const LowerTriMatrixView<T>& B)
    // B = alpha * A * B
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(alpha != T(0));
    TMVAssert(B.dt() == NonUnitDiag);
#ifdef XDEBUG
    //cout<<"Start MultEqMM: alpha = "<<alpha<<endl;
    //cout<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cout<<"B = "<<Type(B)<<"  "<<B<<endl;
    Matrix<T> B0 = B;
    Matrix<Ta> A0 = A;
    Matrix<T> B2 = alpha*A0*B0;
#endif

    if (B.isconj()) MultEqMM(CONJ(alpha),A.Conjugate(),B.Conjugate());
    else if (A.size() > 0) {
      if (alpha == T(1)) 
	if (A.diag().isconj())
	  DoMultEqMM<true,true>(REAL(alpha),A,B);
	else
	  DoMultEqMM<true,false>(REAL(alpha),A,B);
      else if (IMAG(alpha) == RealType(T)(0)) 
	if (A.diag().isconj())
	  DoMultEqMM<false,true>(REAL(alpha),A,B);
	else
	  DoMultEqMM<false,false>(REAL(alpha),A,B);
      else 
	if (A.diag().isconj())
	  DoMultEqMM<false,true>(alpha,A,B);
	else
	  DoMultEqMM<false,false>(alpha,A,B);
    }

#ifdef XDEBUG
    if (Norm(Matrix<T>(B)-B2) > 0.001*(ABS(alpha)*Norm(A0)*Norm(B0))) {
      cerr<<"MultEqMM: alpha = "<<alpha<<endl;
      cerr<<"A = "<<Type(A)<<" step "<<A.diag().step()<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"-> B = "<<B<<endl;
      cerr<<"B2 = "<<B2<<endl;
      abort();
    }
#endif
  }

  template <bool a1, bool ca, bool ub, bool cb, class T, class T1, class Ta, class Tb> 
    static void DoAddMultMM(const T1 alpha,
	const GenDiagMatrix<Ta>& A, const GenLowerTriMatrix<Tb>& B,
	const LowerTriMatrixView<T>& C)
    {
      TMVAssert(A.size() == B.size());
      TMVAssert(A.size() == C.size());
      TMVAssert(A.size() > 0);
      TMVAssert(C.ct()==NonConj);
      TMVAssert(ca == A.diag().isconj());
      TMVAssert(cb == B.isconj());
      TMVAssert(C.dt() == NonUnitDiag);
      TMVAssert(ub == B.isunit());
      TMVAssert(a1 == (alpha == T(1)));

      if (A.size() == 1) {
	const Ta Ax = ca ? CONJ(*A.diag().cptr()) : *A.diag().cptr();
#ifdef TMVFLDEBUG
	TMVAssert(C.ptr() >= C.first);
	TMVAssert(C.ptr() < C.last);
#endif
	if (a1)
	  if (ub) *C.ptr() += Ax;
	  else *C.ptr() += Ax * (cb?CONJ(*B.cptr()):*B.cptr());
	else
	  if (ub) *C.ptr() += alpha * Ax;
	  else *C.ptr() += alpha * Ax * (cb?CONJ(*B.cptr()):*B.cptr());
      } else {
	const int N = A.size();
	const int k = N/2;
	// [ C00  0  ] += [ A00  0  ] * [ B00  0  ]
	// [ C10 C11 ]    [  0  A11 ]   [ B10 B11 ]
	// C00 += A00 * B00
	// C10 += A11 * B10
	// C11 += A11 * B11
	ConstDiagMatrixView<Ta> A00 = A.SubDiagMatrix(0,k);
	ConstDiagMatrixView<Ta> A11 = A.SubDiagMatrix(k,N);
	ConstLowerTriMatrixView<Tb> B00 = B.SubTriMatrix(0,k);
	ConstLowerTriMatrixView<Tb> B11 = B.SubTriMatrix(k,N);
	ConstMatrixView<Tb> B10 = B.SubMatrix(k,N,0,k);
	LowerTriMatrixView<T> C00 = C.SubTriMatrix(0,k);
	LowerTriMatrixView<T> C11 = C.SubTriMatrix(k,N);
	MatrixView<T> C10 = C.SubMatrix(k,N,0,k);

	DoAddMultMM<a1,ca,ub,cb>(alpha,A00,B00,C00);
	C10 += alpha * A11 * B10;
	DoAddMultMM<a1,ca,ub,cb>(alpha,A11,B11,C11);
      }
    } 

  template <bool a1, bool ca, class T, class T1, class Ta, class Tb> 
    static inline void DoAddMultMMa(const T1 alpha,
	const GenDiagMatrix<Ta>& A, const GenLowerTriMatrix<Tb>& B,
	const LowerTriMatrixView<T>& C)
    {
      if (B.isunit())
	if (B.isconj()) DoAddMultMM<a1,ca,true,true>(alpha,A,B,C);
	else DoAddMultMM<a1,ca,false,true>(alpha,A,B,C);
      else
	if (B.isconj()) DoAddMultMM<a1,ca,true,false>(alpha,A,B,C);
	else DoAddMultMM<a1,ca,false,false>(alpha,A,B,C);
    }

  template <class T, class Ta, class Tb> static void AddMultMM(const T alpha,
      const GenDiagMatrix<Ta>& A, const GenLowerTriMatrix<Tb>& B,
      const LowerTriMatrixView<T>& C)
    // C += alpha * A * B
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.size());
    TMVAssert(alpha != T(0));
    TMVAssert(C.dt() == NonUnitDiag);
#ifdef XDEBUG
    //cout<<"Start AddMultMM: alpha = "<<alpha<<endl;
    //cout<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cout<<"B = "<<Type(B)<<"  "<<B<<endl;
    //cout<<"C = "<<Type(C)<<"  "<<C<<endl;
    Matrix<Ta> A0 = A;
    Matrix<Tb> B0 = B;
    Matrix<T> C0 = C;
    Matrix<T> C2 = C0+alpha*A0*B0;
#endif

    if (C.isconj()) AddMultMM(CONJ(alpha),A.Conjugate(),
	B.Conjugate(),C.Conjugate());
    else if (A.size() > 0) {
      if (alpha == T(1)) 
	if (A.diag().isconj())
	  DoAddMultMMa<true,true>(REAL(alpha),A,B,C);
	else
	  DoAddMultMMa<true,false>(REAL(alpha),A,B,C);
      else if (IMAG(alpha) == RealType(T)(0)) 
	if (A.diag().isconj())
	  DoAddMultMMa<false,true>(REAL(alpha),A,B,C);
	else
	  DoAddMultMMa<false,false>(REAL(alpha),A,B,C);
      else 
	if (A.diag().isconj())
	  DoAddMultMMa<false,true>(alpha,A,B,C);
	else
	  DoAddMultMMa<false,false>(alpha,A,B,C);
    }

#ifdef XDEBUG
    if (Norm(Matrix<T>(C)-C2) > 0.001*(ABS(alpha)*Norm(A0)*Norm(B0)+Norm(C0))) {
      cerr<<"AddMultMM: alpha = "<<alpha<<endl;
      cerr<<"A = "<<Type(A)<<" step "<<A.diag().step()<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"-> C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

  template <bool add, class T, class Ta, class Tb> void MultMM(const T alpha,
      const GenDiagMatrix<Ta>& A, const GenLowerTriMatrix<Tb>& B,
      const LowerTriMatrixView<T>& C)
    // C (+)= alpha * A * B
  {
    TMVAssert(A.size() == B.size());
    TMVAssert(A.size() == C.size());
    TMVAssert(C.dt() == NonUnitDiag);
#ifdef XDEBUG
    //cout<<"Start MultMM: alpha = "<<alpha<<endl;
    //cout<<"A = "<<Type(A)<<"  "<<A<<endl;
    //cout<<"B = "<<Type(B)<<"  "<<B<<endl;
    //cout<<"C = "<<Type(C)<<"  "<<C<<endl;
    Matrix<Ta> A0 = A;
    Matrix<Tb> B0 = B;
    Matrix<T> C0 = C;
    Matrix<T> C2 = alpha*A0*B0;
    if (add) C2 += C0;
#endif

    if (A.size() > 0) {
      if (alpha==T(0)) {
	if (!add) C.Zero();
      }
      else if (SameStorage(A.diag(),C)) {
	DiagMatrix<T> tempA = A;
	MultMM<add>(alpha,tempA,B,C);
      }
      else if (!add) {
	C = B;
	MultEqMM(alpha,A,C);
      }
      else if (SameStorage(B,C)) {
	if (B.isrm()) {
	  LowerTriMatrix<T,NonUnitDiag,RowMajor> tempB = B;
	  MultEqMM(alpha,A,tempB.View());
	  C += tempB;
	}
	else {
	  LowerTriMatrix<T,NonUnitDiag,ColMajor> tempB = B;
	  MultEqMM(alpha,A,tempB.View());
	  C += tempB;
	}
      }
      else {
	AddMultMM(alpha,A,B,C);
      }
    }
#ifdef XDEBUG
    if (Norm(Matrix<T>(C)-C2) > 0.001*(ABS(alpha)*Norm(A0)*Norm(B0)+
	  add?Norm(C0):RealType(T)(0))) {
      cerr<<"MultMM: alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"A = "<<Type(A)<<" step "<<A.diag().step()<<"  "<<A0<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"C = "<<Type(C)<<"  "<<C0<<endl;
      cerr<<"-> C = "<<C<<endl;
      cerr<<"C2 = "<<C2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_MultDU.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv
