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



#include "TMV_Matrix.h"
#include "TMV_DiagMatrix.h"
#include "TMV_VectorArith.h"
#include "TMV_MatrixArith.h"
#include "TMV_DiagMatrixArith.h"

//#define XDEBUG

#ifdef XDEBUG
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

  template <class T> ConstVectorView<T> DiagMatrixComposite<T>::cdiag() const
  {
    if (!inst.get()) inst.reset(new DiagMatrix<T>(*this));
    return inst->diag();
  }

  template <class T, class Ta, class Tb> void AddMM(
      const T alpha, const GenDiagMatrix<Ta>& A,
      const T beta, const GenMatrix<Tb>& B, const MatrixView<T>& C)
    // C = alpha*A + beta*B
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == B.rowsize());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(A.size() == C.rowsize());

    if (A.size() > 0) {
      if (SameStorage(A.diag(),C)) {
	DiagMatrix<Ta> tempA = A;
	C = beta * B;
	AddMM(alpha,tempA,C);
      } else {
	C = beta * B;
	AddMM(alpha,A,C);
      }
    }
  }

  template <bool add, class T, class Ta, class Tx> void MultMV(
      const T alpha, const GenDiagMatrix<Ta>& A, const GenVector<Tx>& x,
      const VectorView<T>& y)
    // y (+)= alpha * A * x 
    // yi (+)= alpha * Ai * xi
  {
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());
#ifdef XDEBUG
    //std::cerr<<"MultMV: \n";
    //std::cerr<<"alpha = "<<alpha<<std::endl;
    //std::cerr<<"A = "<<Type(A)<<"  "<<A<<std::endl;
    //std::cerr<<"x = "<<Type(x)<<"  "<<x<<std::endl;
    //std::cerr<<"y = "<<Type(y)<<"  "<<y<<std::endl;
    Vector<T> y0 = y;
    Vector<Tx> x0 = x;
    Matrix<Ta> A0 = A;
    Vector<T> y2 = alpha*A0*x0;
    if (add) y2 += y0;
#endif

    if (y.size() > 0) {
      if (alpha == T(0)) {
	if (!add) y.Zero();
      } 
      else if (!add) {
	if (SameStorage(A.diag(),y)) {
	  if (y.SameAs(A.diag())) ElementProd(alpha,x,y);
	  else if (SameStorage(x,y)) {
	    if (y.SameAs(x)) ElementProd(alpha,A.diag(),y);
	    else {
	      Vector<T> y2 = x;
	      ElementProd(alpha,A.diag(),y2.View());
	      y = y2;
	    }
	  } 
	  else {
	    y = A.diag();
	    ElementProd(alpha,x,y);
	  }
	}
	else {
	  y = x;
	  ElementProd(alpha,A.diag(),y);
	}
      }
      else AddElementProd(alpha,A.diag(),x,y);
    }
#ifdef XDEBUG
    if (Norm(y-y2) > 0.001*(ABS(alpha)*Norm(A0)*Norm(x0)+
	  (add?Norm(y0):RealType(T)(0)))) {
      cerr<<"MultMV: alpha = "<<alpha<<endl;
      cerr<<"add = "<<add<<endl;
      cerr<<"A = "<<Type(A)<<" step "<<A.diag().step()<<"  "<<A0<<endl;
      cerr<<"x = "<<Type(x)<<" step "<<x.step()<<"  "<<x0<<endl;
      cerr<<"y = "<<Type(y)<<" step "<<y.step()<<"  "<<y0<<endl;
      cerr<<"-> y = "<<y<<endl;
      cerr<<"y2 = "<<y2<<endl;
      abort();
    }
#endif
  }

  template <bool rm, bool ca, class T, class Ta> inline void RowMultEqMM(
      const GenDiagMatrix<Ta>& A, const MatrixView<T>& B)
  {
    // B = A * B
    // Bij = Ai * Bij
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(B.ct() == NonConj);
    TMVAssert(rm == B.isrm());
    TMVAssert(ca == A.diag().isconj());

    const Ta* Ai = A.diag().cptr();
    T* Browi = B.ptr();
    const int Astep = A.diag().step();
    const int stepj = B.stepj();
    const int stepi = B.stepi();
    const size_t M = B.colsize();
    const size_t N = B.rowsize();

    for(size_t i=M;i>0;--i,Ai+=Astep,Browi+=stepi) {
      T* Bij = Browi;
      if (*Ai == Ta(0)) 
	for(size_t j=N;j>0;--j,(rm?++Bij:Bij+=stepj))
	  *Bij = T(0);
      else if (IMAG(*Ai) == RealType(Ta)(0))
	for(size_t j=N;j>0;--j,(rm?++Bij:Bij+=stepj))
	  *Bij *= REAL(*Ai);
      else
	for(size_t j=N;j>0;--j,(rm?++Bij:Bij+=stepj))
	  *Bij *= (ca?CONJ(*Ai):*Ai);
    }
  }

  template <bool cm, bool ca, class T, class Ta> inline void DoColMultEqMM(
      const GenDiagMatrix<Ta>& A, const MatrixView<T>& B)
  {
    // B = A * B 
    // Bij = Ai * Bij
    TMVAssert(A.size() == B.colsize());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(B.ct()==NonConj);
    TMVAssert(A.diag().step() == 1);
    TMVAssert(cm == B.iscm());
    TMVAssert(ca == A.diag().isconj());

    const Ta*const Aptr = A.diag().cptr();
    T* Bcolj = B.ptr();
    const int stepj = B.stepj();
    const int stepi = B.stepi();
    const size_t M = B.colsize();
    const size_t N = B.rowsize();

    for(size_t j=N;j>0;--j,Bcolj+=stepj) {
      T* Bij = Bcolj;
      const Ta* Ai = Aptr;
      for(size_t i=M;i>0;--i,++Ai,(cm?++Bij:Bij+=stepi))
	*Bij *= (ca ? CONJ(*Ai) : *Ai);
    }
  }

  template <bool cm, bool ca, class T, class Ta> inline void ColMultEqMM(
      const GenDiagMatrix<Ta>& A, const MatrixView<T>& B)
  {
    if (A.diag().step() == 1)
      DoColMultEqMM<cm,ca>(A,B);
    else {
      DiagMatrix<Ta> AA = A;
      DoColMultEqMM<cm,false>(AA,B);
    }
  }

  template <class T, class Ta> inline void MultEqMM(const T alpha,
      const GenDiagMatrix<Ta>& A, const MatrixView<T>& B)
    // B = alpha * A * B
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(alpha != T(0));
#ifdef XDEBUG
    Matrix<T> B0 = B;
    Matrix<Ta> A0 = A;
    Matrix<T> B2 = alpha*A0*B0;
#endif

    if (B.isconj()) MultEqMM(CONJ(alpha),A.Conjugate(),B.Conjugate());
    else if (B.colsize() > 0 && B.rowsize() > 0) {
      if (alpha != T(1)) {
	if (IMAG(alpha) == RealType(T)(0)) {
	  DiagMatrix<Ta> AA = REAL(alpha) * A;
	  if (B.isrm()) RowMultEqMM<true,false>(AA,B);
	  else if (B.iscm()) DoColMultEqMM<true,false>(AA,B);
	  else if (B.colsize() > B.rowsize()) DoColMultEqMM<false,false>(AA,B);
	  else RowMultEqMM<false,false>(AA,B);
	}
	else {
	  // AA = alpha * A;
	  DiagMatrix<T> AA = alpha * A;
	  if (B.isrm()) RowMultEqMM<true,false>(AA,B);
	  else if (B.iscm()) DoColMultEqMM<true,false>(AA,B);
	  else if (B.colsize() > B.rowsize()) DoColMultEqMM<false,false>(AA,B);
	  else RowMultEqMM<false,false>(AA,B);
	}
      }
      else if (A.diag().isconj()) {
        if (B.isrm()) RowMultEqMM<true,true>(A,B);
	else if (B.iscm()) ColMultEqMM<true,true>(A,B);
	else if (B.colsize() > B.rowsize()) ColMultEqMM<false,true>(A,B);
	else RowMultEqMM<false,true>(A,B);
      }
      else {
        if (B.isrm()) RowMultEqMM<true,false>(A,B);
	else if (B.iscm()) ColMultEqMM<true,false>(A,B);
	else if (B.colsize() > B.rowsize()) ColMultEqMM<false,false>(A,B);
	else RowMultEqMM<false,false>(A,B);
      }
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

  template <bool rm, bool ca, bool cb, class T, class Ta, class Tb>
    inline void DoRowAddMultMM(
	const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
	const MatrixView<T>& C)
    {
      // C += A * B
      // Cij += Ai * Bij
      TMVAssert(A.size() == B.colsize());
      TMVAssert(A.size() == C.colsize());
      TMVAssert(B.rowsize() == C.rowsize());
      TMVAssert(C.rowsize() > 0);
      TMVAssert(C.colsize() > 0);
      TMVAssert(C.ct() == NonConj);
      TMVAssert(rm == (B.isrm() && C.isrm()));
      TMVAssert(ca == A.diag().isconj());
      TMVAssert(cb == B.isconj());

      const Ta* Ai = A.diag().cptr();
      const Tb* Browi = B.cptr();
      T* Crowi = C.ptr();
      const int Astep = A.diag().step();
      const int Bstepj = B.stepj();
      const int Bstepi = B.stepi();
      const int Cstepj = C.stepj();
      const int Cstepi = C.stepi();
      const size_t M = C.colsize();
      const size_t N = C.rowsize();

      for(size_t i=M;i>0;--i,Ai+=Astep,Browi+=Bstepi,Crowi+=Cstepi) {
	const Tb* Bij = Browi;
	T* Cij = Crowi;
	if (IMAG(*Ai) == RealType(Ta)(0)) {
	  if (REAL(*Ai) != RealType(Ta)(0))
	    for(size_t j=N;j>0;--j,(rm?++Bij:Bij+=Bstepj),
		(rm?++Cij:Cij+=Cstepj))
	      *Cij += REAL(*Ai)*(cb?CONJ(*Bij):*Bij);
	}
	else
	  for(size_t j=N;j>0;--j,(rm?++Bij:Bij+=Bstepj),
	      (rm?++Cij:Cij+=Cstepj))
	    *Cij += (ca?CONJ(*Ai):*Ai)*(cb?CONJ(*Bij):*Bij);
      }
    }

  template <bool rm, bool ca, class T, class Ta, class Tb>
    inline void RowAddMultMM(
	const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
	const MatrixView<T>& C)
    {
      if (B.isconj()) DoRowAddMultMM<rm,ca,true>(A,B,C);
      else DoRowAddMultMM<rm,ca,false>(A,B,C);
    }

  template <bool cm, bool ca, bool cb, class T, class Ta, class Tb> 
    inline void DoColAddMultMM(
	const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
	const MatrixView<T>& C)
    {
      // C += A * B 
      // Cij = Ai * Bij
      TMVAssert(A.size() == B.colsize());
      TMVAssert(A.size() == C.colsize());
      TMVAssert(B.rowsize() == C.rowsize());
      TMVAssert(C.rowsize() > 0);
      TMVAssert(C.colsize() > 0);
      TMVAssert(A.diag().step() == 1);
      TMVAssert(cm == (B.iscm() && C.iscm()));
      TMVAssert(ca == A.diag().isconj());
      TMVAssert(cb == B.isconj());

      const Ta*const Aptr = A.diag().cptr();
      const Tb* Bcolj = B.cptr();
      T* Ccolj = C.ptr();
      const int Cstepj = C.stepj();
      const int Cstepi = C.stepi();
      const int Bstepj = B.stepj();
      const int Bstepi = B.stepi();
      const size_t M = C.colsize();
      const size_t N = C.rowsize();

      for(size_t j=N;j>0;--j,Bcolj+=Bstepj,Ccolj+=Cstepj) {
	const Tb* Bij = Bcolj;
	T* Cij = Ccolj;
	const Ta* Ai = Aptr;
	for(size_t i=M;i>0;--i,++Ai,(cm?++Bij:Bij+=Bstepi),
	    (cm?++Cij:Cij+=Cstepi))
	  *Cij += (ca ? CONJ(*Ai) : *Ai) * (cb ? CONJ(*Bij) : *Bij);
      }
    }

  template <bool cm, bool ca, class T, class Ta, class Tb> 
    inline void ColAddMultMM(
	const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
	const MatrixView<T>& C)
    { 
      if (A.diag().step() == 1)
	if (B.isconj())
	  DoColAddMultMM<cm,ca,true>(A,B,C);
	else
	  DoColAddMultMM<cm,ca,false>(A,B,C);
      else {
	DiagMatrix<Ta> AA = A;
	if (B.isconj())
	  DoColAddMultMM<cm,ca,true>(AA,B,C);
	else
	  DoColAddMultMM<cm,ca,false>(AA,B,C);
      }
    }

  template <class T, class Ta, class Tb> inline void AddMultMM(const T alpha,
      const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
    // C += alpha * A * B
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha != T(0));
#ifdef XDEBUG
    Matrix<Ta> A0 = A;
    Matrix<Tb> B0 = B;
    Matrix<T> C0 = C;
    Matrix<T> C2 = C0+alpha*A0*B0;
#endif

    if (C.isconj()) AddMultMM(CONJ(alpha),A.Conjugate(),
	B.Conjugate(),C.Conjugate());
    else if (C.colsize() > 0 && C.rowsize() > 0) {
      if (alpha != T(1)) {
	if (IMAG(alpha) == RealType(T)(0)) {
	  DiagMatrix<Ta> AA = REAL(alpha) * A;
	  if (B.isrm() && C.isrm()) 
	    RowAddMultMM<true,false>(AA,B,C);
	  else if (B.iscm() && C.iscm()) 
	    ColAddMultMM<true,false>(AA,B,C);
	  else if (B.colsize() > B.rowsize()) 
	    ColAddMultMM<false,false>(AA,B,C);
	  else 
	    RowAddMultMM<false,false>(AA,B,C);
	}
	else {
	  DiagMatrix<T> AA = alpha * A;
	  if (B.isrm() && C.isrm()) 
	    RowAddMultMM<true,false>(AA,B,C);
	  else if (B.iscm() && C.iscm()) 
	    ColAddMultMM<true,false>(AA,B,C);
	  else if (B.colsize() > B.rowsize()) 
	    ColAddMultMM<false,false>(AA,B,C);
	  else 
	    RowAddMultMM<false,false>(AA,B,C);
	}
      }
      else if (A.diag().isconj()) {
        if (B.isrm() && C.isrm()) 
	  RowAddMultMM<true,true>(A,B,C);
	else if (B.iscm() && C.iscm()) 
	  ColAddMultMM<true,true>(A,B,C);
	else if (B.colsize() > B.rowsize()) 
	  ColAddMultMM<false,true>(A,B,C);
	else 
	  RowAddMultMM<false,true>(A,B,C);
      }
      else {
        if (B.isrm() && C.isrm()) 
	  RowAddMultMM<true,false>(A,B,C);
	else if (B.iscm() && C.iscm()) 
	  ColAddMultMM<true,false>(A,B,C);
	else if (B.colsize() > B.rowsize()) 
	  ColAddMultMM<false,false>(A,B,C);
	else 
	  RowAddMultMM<false,false>(A,B,C);
      }
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
      const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
    // C (+)= alpha * A * B
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
#ifdef XDEBUG
    Matrix<Ta> A0 = A;
    Matrix<Tb> B0 = B;
    Matrix<T> C0 = C;
    Matrix<T> C2 = alpha*A0*B0;
    if (add) C2 += C0;
#endif

    if (C.colsize() > 0 && C.rowsize() > 0) {
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
	  Matrix<T,RowMajor> tempB = B;
	  MultEqMM(alpha,A,tempB.View());
	  C += tempB;
	}
	else {
	  Matrix<T,ColMajor> tempB = B;
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

#define InstFile "TMV_DiagMatrixArith.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


