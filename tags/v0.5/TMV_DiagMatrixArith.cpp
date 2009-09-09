
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
    Vector<T> y2 = beta*Vector<T>(y)+alpha*Matrix<T>(A)*Vector<T>(x);
    Vector<T> y0 = y;
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
	y *= beta;
	y += temp;
      } else {
	y *= beta;
	AddElementProd(alpha,A.diag(),x,y);
      }
    }
#ifdef XDEBUG
    if (Norm(y-y2) > 0.001*Norm(y)) {
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
    for(size_t i=0;i<B.colsize();++i,++Ait,Bptr+=B.stepi()) {
      if (*Ait == Ta(0)) {
	std::fill(Bptr,Bptr+B.rowsize(),Tb(0));
      } else {
	Tb aa = *Ait;
	if (alpha != T1(1)) aa *= alpha;
	if (IMAG(aa) == RealType(Tb)(0)) {
	  // B.row(i) *= alpha * (*Ait);
	  if (REAL(aa) != RealType(Tb)(1))
	    DoMultXV(REAL(aa),VIt<Tb,Unit,NonConj>(Bptr,1 
		  FIRSTLAST1(B.first,B.last) ),B.rowsize());
	}
	else
	  DoMultXV(aa,VIt<Tb,Unit,NonConj>(Bptr,1 
		FIRSTLAST1(B.first,B.last) ),B.rowsize());
      }
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

    if (B.isconj()) MultEqMM(CONJ(alpha),A.Conjugate(),B.Conjugate());
    else if (B.isrm()) DoRowMajorMultEqMM(alpha,A,B);
    else if (B.iscm() && A.diag().step()==1) DoColMajorMultEqMM(alpha,A,B);
    else ColMultEqMM(alpha,A,B);
#ifdef XDEBUG
    if (Norm(B-B2) > 0.001*Norm(B)) {
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
    for(size_t i=0;i<C.colsize();++i,++Ait,Bptr+=B.stepi(),
	Cptr+=C.stepi()) {
      if (*Ait != Ta(0)) {
	Tc aa = *Ait;
	if (alpha != T1(1)) aa *= alpha;
	// C.row(i) += alpha * (*Ait) * B.row(i);
	if (IMAG(aa) == RealType(Tc)(0)) 
	  if (B.isconj())
	    DoAddVV(REAL(aa),CVIt<Tb,Unit,Conj>(Bptr,1),
		VIt<Tc,Unit,NonConj>(Cptr,1 FIRSTLAST1(C.first,C.last) ),
		C.rowsize());
	  else
	    DoAddVV(REAL(aa),CVIt<Tb,Unit,NonConj>(Bptr,1),
		VIt<Tc,Unit,NonConj>(Cptr,1 FIRSTLAST1(C.first,C.last) ),
		C.rowsize());
	else
	  if (B.isconj())
	    DoAddVV(aa,CVIt<Tb,Unit,Conj>(Bptr,1),
		VIt<Tc,Unit,NonConj>(Cptr,1 FIRSTLAST1(C.first,C.last) ),
		C.rowsize());
	  else
	    DoAddVV(aa,CVIt<Tb,Unit,NonConj>(Bptr,1),
		VIt<Tc,Unit,NonConj>(Cptr,1 FIRSTLAST1(C.first,C.last) ),
		C.rowsize());
      }
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

    if (C.isconj()) AddMultEqMM(CONJ(alpha),A.Conjugate(),
	B.Conjugate(),C.Conjugate());
    else if (C.isrm() && B.isrm()) DoRowMajorAddMultEqMM(alpha,A,B,C);
    else if (C.iscm() && B.iscm() && A.diag().step()==1) 
      DoColMajorAddMultEqMM(alpha,A,B,C);
    else ColAddMultEqMM(alpha,A,B,C);
#ifdef XDEBUG
    if (Norm(C-C2) > 0.001*Norm(C)) {
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
	AddMultEqMM(alpha,A,B,C);
      }
    }
#ifdef XDEBUG
    if (Norm(C-C2) > 0.001*Norm(C)) {
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


