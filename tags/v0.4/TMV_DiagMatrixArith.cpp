
#include "TMV.h"
#include "TMV_Diag.h"
#include "TMV_VectorArith_Inline.h"

//#define XDEBUG

namespace tmv {

  template <class T, class Ta> void AddMM(
      const T alpha, const GenDiagMatrix<Ta>& A,
      const T beta, const MatrixView<T>& B)
  {
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == B.rowsize());

    if (beta == T(1)) AddVV(alpha,A.diag(),B.diag());
    else if (A.SameStorageAs(B)) {
      Vector<Ta> A2 = A.diag();
      MultXM(beta,B);
      AddVV(alpha,A2,B.diag());
    } else {
      MultXM(beta,B);
      AddVV(alpha,A.diag(),B.diag());
    }
  }

  template <class T, class Ta, class Tx> void MultMV(
      const T alpha, const GenDiagMatrix<Ta>& A, const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
    // y = alpha * A * x + beta * y
    // yi = alpha * Ai * xi + beta * yi
  {
#ifdef XDEBUG
    Vector<T> y2 = beta*Vector<T>(y)+alpha*Matrix<T>(A)*Vector<T>(x);
    Vector<T> y0 = y;
#endif
    TMVAssert(A.size() == x.size());
    TMVAssert(A.size() == y.size());

    if (y.size() > 0) {
      if (alpha == T(0)) {
	MultXV(beta,y);
      } else {
	if (beta == T(0)) {
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
	  MultXV(beta,y);
	  AddVV(T(1),temp,y);
	} else {
	  MultXV(beta,y);
	  AddElementProd(alpha,A.diag(),x,y);
	}
      }
    }
#ifdef XDEBUG
    if (Norm(y-y2) > 1.e-4) {
      cerr<<"MultMV: alpha, beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<" step "<<A.diag().step()<<"  "<<A<<endl;
      cerr<<"x = "<<Type(x)<<" step "<<x.step()<<"  "<<x<<endl;
      cerr<<"y = "<<Type(y)<<" step "<<y.step()<<"  "<<y0<<endl;
      cerr<<"-> y = "<<y<<endl;
      cerr<<"y2 = "<<y2<<endl;
      abort();
    }
#endif

  }

  template <class T1, class Ta, class Tb> void RowMajorMultEqMM(
      const T1 alpha, const GenDiagMatrix<Ta>& A, const MatrixView<Tb>& B)
  {
    // B = alpha * A * B
    // Bij = alpha * Ai * Bij
    TMVAssert(A.size() == B.colsize());
    TMVAssert(alpha!=T1(0));
    TMVAssert(B.isrm());
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(B.ct() == NonConj);

    CVIter<Ta> Ait = A.diag().begin();
    Tb* Bptr = B.ptr();
    if (alpha == T1(1)) 
      for(size_t i=0;i<B.colsize();++i,++Ait,Bptr+=B.stepi()) {
	if (*Ait != Ta(0)) 
	  DoMultXV(*Ait,VIt<Tb,Unit,NonConj>(Bptr,1 
		FIRSTLAST1(B.first,B.last) ),B.rowsize());
      }
    else 
      for(size_t i=0;i<B.colsize();++i,++Ait,Bptr+=B.stepi()) {
	if (*Ait != Ta(0)) 
	  DoMultXV(*Ait*alpha,VIt<Tb,Unit,NonConj>(Bptr,1 
		FIRSTLAST1(B.first,B.last) ),B.rowsize());
      }
  }

  template <class T, class Ta> void DoRowMajorMultEqMM(
      const T alpha, const GenDiagMatrix<Ta>& A, const MatrixView<T>& B)
  {
    if (IMAG(alpha) == RealType(T)(0)) RowMajorMultEqMM(REAL(alpha),A,B);
    else RowMajorMultEqMM(alpha,A,B);
  }
      
  template <class T1, class Ta, ConjItType Ca, class Tb> void ColMajorMultEqMM(
      const T1 alpha, const CVIt<Ta,Unit,Ca>& Ait, const MatrixView<Tb>& B)
  {
    // B = alpha * A * B 
    // Bij = alpha * Ai * Bij
    TMVAssert(alpha!=T1(0));
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(B.ct()==NonConj);
    TMVAssert(B.iscm());

    Tb* Bptr = B.ptr();
    for(size_t j=0;j<B.rowsize();++j,Bptr+=B.stepj())
      DoElementProd(alpha,Ait,VIt<Tb,Unit,NonConj>(Bptr,1 
	    FIRSTLAST1(B.first,B.last) ),B.colsize());
  }

  template <class T, class Ta> void DoColMajorMultEqMM(
      const T alpha, const GenDiagMatrix<Ta>& A, const MatrixView<T>& B)
  {
    if (IMAG(alpha) == RealType(T)(0)) 
      if (A.diag().isconj())
	ColMajorMultEqMM(REAL(alpha),CVIt<Ta,Unit,Conj>(A.diag().begin()),B);
      else
	ColMajorMultEqMM(REAL(alpha),CVIt<Ta,Unit,NonConj>(A.diag().begin()),B);
    else 
      if (A.diag().isconj())
	ColMajorMultEqMM(alpha,CVIt<Ta,Unit,Conj>(A.diag().begin()),B);
      else
	ColMajorMultEqMM(alpha,CVIt<Ta,Unit,NonConj>(A.diag().begin()),B);
  }
      
  template <class T, class Ta> void ColMultEqMM(const T alpha,
      const GenDiagMatrix<Ta>& A, const MatrixView<T>& B)
  {
    // B = alpha * A * B 
    // Bij = alpha * Ai * Bij
    TMVAssert(A.size() == B.colsize());
    TMVAssert(alpha!=T(0));
    TMVAssert(B.rowsize() > 0);
    TMVAssert(B.colsize() > 0);
    TMVAssert(B.ct()==NonConj);

    for(size_t j=0;j<B.rowsize();++j) ElementProd(alpha,A.diag(),B.col(j));
  }

  template <class T, class Ta> void MultEqMM(const T alpha,
      const GenDiagMatrix<Ta>& A, const MatrixView<T>& B)
    // B = alpha * A * B
  {
#ifdef XDEBUG
    Matrix<T> B2 = alpha*Matrix<T>(A)*Matrix<T>(B);
    Matrix<T> B0 = B;
#endif
    TMVAssert(A.size() == B.colsize());
    TMVAssert(alpha != T(0));

    if (B.isconj()) MultEqMM(CONJ(alpha),A.QuickConjugate(),B.QuickConjugate());
    else if (B.isrm()) DoRowMajorMultEqMM(alpha,A,B);
    else if (B.iscm() && A.diag().step()==1) DoColMajorMultEqMM(alpha,A,B);
    else ColMultEqMM(alpha,A,B);
#ifdef XDEBUG
    if (Norm(B-B2) > 1.e-4) {
      cerr<<"MultEqMM: alpha = "<<alpha<<endl;
      cerr<<"A = "<<Type(A)<<" step "<<A.diag().step()<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B0<<endl;
      cerr<<"-> B = "<<B<<endl;
      cerr<<"B2 = "<<B2<<endl;
      abort();
    }
#endif
  }

  template <class T1, class Ta, class Tb, class Tc> void RowMajorAddMultEqMM(
      const T1 alpha, const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<Tc>& C)
  {
    // C = alpha * A * B + beta * C
    // Cij = alpha * Ai * Bij + beta Cij
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha!=T1(0));
    TMVAssert(B.isrm());
    TMVAssert(C.isrm());
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct() == NonConj);

    CVIter<Ta> Ait = A.diag().begin();
    Tc* Cptr = C.ptr();
    const Tb* Bptr = B.cptr();
    if (alpha == T1(1)) 
      for(size_t i=0;i<C.colsize();++i,++Ait,Bptr+=B.stepi(),
	  Cptr+=C.stepi()) {
	if (*Ait != Ta(0)) 
	  if (B.isconj())
	    DoAddVV(*Ait,CVIt<Tb,Unit,Conj>(Bptr,1),
		VIt<Tc,Unit,NonConj>(Cptr,1 FIRSTLAST1(C.first,C.last) ),
		C.rowsize());
	  else
	    DoAddVV(*Ait,CVIt<Tb,Unit,NonConj>(Bptr,1),
		VIt<Tc,Unit,NonConj>(Cptr,1 FIRSTLAST1(C.first,C.last) ),
		C.rowsize());
      }
    else 
      for(size_t i=0;i<C.colsize();++i,++Ait,Bptr+=B.stepi(),
	  Cptr+=C.stepi()) {
	if (*Ait != Ta(0)) 
	  if (B.isconj())
	    DoAddVV(*Ait*alpha,CVIt<Tb,Unit,Conj>(Bptr,1),
		VIt<Tc,Unit,NonConj>(Cptr,1 FIRSTLAST1(C.first,C.last) ),
		C.rowsize());
	  else
	    DoAddVV(*Ait*alpha,CVIt<Tb,Unit,NonConj>(Bptr,1),
		VIt<Tc,Unit,NonConj>(Cptr,1 FIRSTLAST1(C.first,C.last) ),
		C.rowsize());
      }
  }

  template <class T, class Ta, class Tb> void DoRowMajorAddMultEqMM(
      const T alpha, const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
    if (IMAG(alpha) == RealType(T)(0)) RowMajorAddMultEqMM(REAL(alpha),A,B,C);
    else RowMajorAddMultEqMM(alpha,A,B,C);
  }
      
  template <class T1, class Ta, ConjItType Ca, class Tb, class Tc> 
    void ColMajorAddMultEqMM(
	const T1 alpha, CVIt<Ta,Unit,Ca> Ait,
	const GenMatrix<Tb>& B, const MatrixView<Tc>& C)
    {
      // C += alpha * A * B 
      // Cij += alpha * Ai * Bij
      TMVAssert(B.colsize() == C.colsize());
      TMVAssert(B.rowsize() == C.rowsize());
      TMVAssert(alpha!=T1(0));
      TMVAssert(C.rowsize() > 0);
      TMVAssert(C.colsize() > 0);
      TMVAssert(B.iscm());
      TMVAssert(C.iscm());
      TMVAssert(C.ct()==NonConj);

      Tc* Cptr = C.ptr();
      const Tb* Bptr = B.cptr();
      for(size_t j=0;j<C.rowsize();++j,Bptr+=B.stepj(),Cptr+=C.stepj()) {
	if (B.isconj())
	  DoAddElementProd(alpha,Ait,CVIt<Tb,Unit,Conj>(Bptr,1),
	      VIt<Tc,Unit,NonConj>(Cptr,1 FIRSTLAST1(C.first,C.last) ),
	      C.colsize());
	else
	  DoAddElementProd(alpha,Ait,CVIt<Tb,Unit,NonConj>(Bptr,1),
	      VIt<Tc,Unit,NonConj>(Cptr,1 FIRSTLAST1(C.first,C.last) ),
	      C.colsize());
      }
    }

  template <class T, class Ta, class Tb> void DoColMajorAddMultEqMM(
      const T alpha, const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
    if (IMAG(alpha) == RealType(T)(0)) 
      if (A.diag().isconj())
	ColMajorAddMultEqMM(REAL(alpha),CVIt<Ta,Unit,Conj>(A.diag().begin()),
	    B,C);
      else
	ColMajorAddMultEqMM(REAL(alpha),CVIt<Ta,Unit,NonConj>(A.diag().begin()),
	    B,C);
    else
      if (A.diag().isconj())
	ColMajorAddMultEqMM(alpha,CVIt<Ta,Unit,Conj>(A.diag().begin()),B,C);
      else
	ColMajorAddMultEqMM(alpha,CVIt<Ta,Unit,NonConj>(A.diag().begin()),B,C);
  }

  template <class T, class Ta, class Tb> void ColAddMultEqMM(const T alpha,
      const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
    // C += alpha * A * B 
    // Cij += alpha * Ai * Bij
    TMVAssert(A.size() == B.colsize());
    TMVAssert(alpha!=T(0));
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct()==NonConj);

    for(size_t j=0;j<C.rowsize();++j) 
      AddElementProd(alpha,A.diag(),B.col(j),C.col(j));
  }

  template <class T, class Ta, class Tb> void AddMultEqMM(const T alpha,
      const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
    // C += alpha * A * B
  {
#ifdef XDEBUG
    Matrix<T> C2 = Matrix<T>(C)+alpha*Matrix<T>(A)*Matrix<T>(B);
    Matrix<T> C0 = C;
#endif
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha != T(0));

    if (C.isconj()) AddMultEqMM(CONJ(alpha),A.QuickConjugate(),
	B.QuickConjugate(),C.QuickConjugate());
    else if (C.isrm() && B.isrm()) DoRowMajorAddMultEqMM(alpha,A,B,C);
    else if (C.iscm() && B.iscm() && A.diag().step()==1) 
      DoColMajorAddMultEqMM(alpha,A,B,C);
    else ColAddMultEqMM(alpha,A,B,C);
#ifdef XDEBUG
    if (Norm(C-C2) > 1.e-4) {
      cerr<<"AddMultEqMM: alpha = "<<alpha<<endl;
      cerr<<"A = "<<Type(A)<<" step "<<A.diag().step()<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
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
    Matrix<T> C2 = beta*Matrix<T>(C)+alpha*Matrix<T>(A)*Matrix<T>(B);
    Matrix<T> C0 = C;
#endif
    TMVAssert(A.size() == B.colsize());
    TMVAssert(A.size() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());

    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (alpha==T(0)) {
	MultXM(beta,C);
      } else if (beta == T(0)) {
	C = B;
	MultEqMM(alpha,A,C);
      } else if (C.SameStorageAs(B)) {
	if (B.isrm()) {
	  Matrix<T,RowMajor> tempB = B;
	  MultEqMM(alpha,A,tempB.QuickView());
	  MultXM(beta,C);
	  AddMM(T(1),tempB,T(1),C);
	} else {
	  Matrix<T,ColMajor> tempB = B;
	  MultEqMM(alpha,A,tempB.QuickView());
	  MultXM(beta,C);
	  AddMM(T(1),tempB,T(1),C);
	}
      } else {
	MultXM(beta,C);
	AddMultEqMM(alpha,A,B,C);
      }
    }
#ifdef XDEBUG
    if (Norm(C-C2) > 1.e-4) {
      cerr<<"MultMM: alpha,beta = "<<alpha<<"  "<<beta<<endl;
      cerr<<"A = "<<Type(A)<<" step "<<A.diag().step()<<"  "<<A<<endl;
      cerr<<"B = "<<Type(B)<<"  "<<B<<endl;
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


