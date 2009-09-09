
#include "TMV.h"
#include "TMV_Diag.h"

//#define XDEBUG

namespace tmv {

  template <class T, class Ta> void AddMM(
      const T alpha, const GenDiagMatrix<Ta>& A,
      const T beta, const MatrixView<T>& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == B.rowsize());

    if (beta == T(1)) 
      B.diag() += alpha * A.diag();
    else if (A.SameStorageAs(B)) {
      Vector<Ta> A2 = A.diag();
      B *= beta;
      B.diag() += alpha * A2;
    } else {
      B *= beta;
      B.diag() += alpha * A.diag();
    }
  }

  template <class T, class Ta, class Tx> void MultMV(
      const T alpha, const GenDiagMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
    // y = alpha * A * x + beta * y
    // yi = alpha * Ai * xi + beta * yi
  {
#ifdef XDEBUG
    Vector<T> y0 = y;
    Vector<T> x0 = x;
    Matrix<T> A0 = A;
    Vector<T> y2 = beta*y0+alpha*A0*x0;
#endif
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());

    if (y.size() > 0) {
      if (alpha == T(0)) {
	y *= beta;
      } else if (beta == T(0)) {
	if (y.SameStorageAs(A.diag())) {
	  if (y.SameAs(A.diag())) ElementProd(alpha,x,y);
	  else {
	    if (y.SameStorageAs(x)) {
	      if (y.SameAs(x)) ElementProd(alpha,A.diag(),y);
	      else {
		Vector<T> y2 = x;
		ElementProd(alpha,A.diag(),y2.View());
	      }
	    } else {
	      y = A.diag();
	      ElementProd(alpha,x,y);
	    }
	  }
	} else {
	  y = x;
	  ElementProd(alpha,A.diag(),y);
	}
      } else if (beta == T(1)) {
	AddElementProd(alpha,A.diag(),x,y);
      } else if (y.SameStorageAs(A.diag()) || y.SameStorageAs(x)) {
	Vector<T> temp = x;
	ElementProd(alpha,A.diag(),temp.View());
	AddVV(T(1),temp,beta,y);
      } else {
	y *= beta;
	AddElementProd(alpha,A.diag(),x,y);
      }
    }
#ifdef XDEBUG
    if (Norm(y-y2) > 0.001*(abs(alpha)*Norm(A0)*Norm(x0)+
	  (beta==T(0)?RealType(T)(0):abs(beta)*Norm(y0)))) {
      cerr<<"MultMV: alpha, beta = "<<alpha<<"  "<<beta<<endl;
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
#ifdef XDEBUG
    Matrix<T> B0 = B;
    Matrix<T> A0 = A;
    Matrix<T> B2 = alpha*A0*B0;
#endif
    TMVAssert(A.size() == B.colsize());
    TMVAssert(alpha != T(0));

    if (B.isconj()) MultEqMM(CONJ(alpha),A.Conjugate(),B.Conjugate());
    else if (B.colsize() > 0 && B.rowsize() > 0) {
      if (alpha != T(1)) {
	if (IMAG(alpha) == RealType(T)(0)) {
	  DiagMatrix<Ta> AA = REAL(alpha)*A;
	  if (B.isrm()) RowMultEqMM<true,false>(AA,B);
	  else if (B.iscm()) DoColMultEqMM<true,false>(AA,B);
	  else if (B.colsize() > B.rowsize()) DoColMultEqMM<false,false>(AA,B);
	  else RowMultEqMM<false,false>(AA,B);
	} else {
	  DiagMatrix<T> AA = alpha*A;
	  if (B.isrm()) RowMultEqMM<true,false>(AA,B);
	  else if (B.iscm()) DoColMultEqMM<true,false>(AA,B);
	  else if (B.colsize() > B.rowsize()) DoColMultEqMM<false,false>(AA,B);
	  else RowMultEqMM<false,false>(AA,B);
	}
      } else if (A.diag().isconj()) {
        if (B.isrm()) RowMultEqMM<true,true>(A,B);
	else if (B.iscm()) ColMultEqMM<true,true>(A,B);
	else if (B.colsize() > B.rowsize()) ColMultEqMM<false,true>(A,B);
	else RowMultEqMM<false,true>(A,B);
      } else {
        if (B.isrm()) RowMultEqMM<true,false>(A,B);
	else if (B.iscm()) ColMultEqMM<true,false>(A,B);
	else if (B.colsize() > B.rowsize()) ColMultEqMM<false,false>(A,B);
	else RowMultEqMM<false,false>(A,B);
      }
    }

#ifdef XDEBUG
    if (Norm(B-B2) > 0.001*(abs(alpha)*Norm(A0)*Norm(B0))) {
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
#ifdef XDEBUG
    Matrix<T> A0 = A;
    Matrix<T> B0 = B;
    Matrix<T> C0 = C;
    Matrix<T> C2 = C0+alpha*A0*B0;
#endif
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha != T(0));

    if (C.isconj()) AddMultMM(CONJ(alpha),A.Conjugate(),
	B.Conjugate(),C.Conjugate());
    else if (C.colsize() > 0 && C.rowsize() > 0) {
      if (alpha != T(1)) {
	if (IMAG(alpha) == RealType(T)(0)) {
	  DiagMatrix<Ta> AA = REAL(alpha)*A;
	  if (B.isrm() && C.isrm()) 
	    RowAddMultMM<true,false>(AA,B,C);
	  else if (B.iscm() && C.iscm()) 
	    ColAddMultMM<true,false>(AA,B,C);
	  else if (B.colsize() > B.rowsize()) 
	    ColAddMultMM<false,false>(AA,B,C);
	  else 
	    RowAddMultMM<false,false>(AA,B,C);
	} else {
	  DiagMatrix<T> AA = alpha*A;
	  if (B.isrm() && C.isrm()) 
	    RowAddMultMM<true,false>(AA,B,C);
	  else if (B.iscm() && C.iscm()) 
	    ColAddMultMM<true,false>(AA,B,C);
	  else if (B.colsize() > B.rowsize()) 
	    ColAddMultMM<false,false>(AA,B,C);
	  else 
	    RowAddMultMM<false,false>(AA,B,C);
	}
      } else if (A.diag().isconj()) {
        if (B.isrm() && C.isrm()) 
	  RowAddMultMM<true,true>(A,B,C);
	else if (B.iscm() && C.iscm()) 
	  ColAddMultMM<true,true>(A,B,C);
	else if (B.colsize() > B.rowsize()) 
	  ColAddMultMM<false,true>(A,B,C);
	else 
	  RowAddMultMM<false,true>(A,B,C);
      } else {
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
    if (Norm(C-C2) > 0.001*(abs(alpha)*Norm(A0)*Norm(B0)+Norm(C0))) {
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

  template <class T, class Ta, class Tb> void MultMM(const T alpha,
      const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C
  {
#ifdef XDEBUG
    Matrix<T> A0 = A;
    Matrix<T> B0 = B;
    Matrix<T> C0 = C;
    Matrix<T> C2 = beta*C0+alpha*A0*B0;
#endif
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());

    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (alpha==T(0)) {
	C *= beta;
      } else if (beta == T(0)) {
	C = B;
	MultEqMM(alpha,A,C);
      } else if (C.SameStorageAs(B)) {
	if (B.isrm()) {
	  Matrix<T,RowMajor> tempB = B;
	  MultEqMM(alpha,A,tempB.View());
	  C *= beta;
	  C += tempB;
	} else {
	  Matrix<T,ColMajor> tempB = B;
	  MultEqMM(alpha,A,tempB.View());
	  C *= beta;
	  C += tempB;
	}
      } else {
	C *= beta;
	AddMultMM(alpha,A,B,C);
      }
    }
#ifdef XDEBUG
    if (Norm(C-C2) > 0.001*(abs(alpha)*Norm(A0)*Norm(B0)+
	  beta==T(0)?RealType(T)(0):abs(beta)*Norm(C0))) {
      cerr<<"MultMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
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


