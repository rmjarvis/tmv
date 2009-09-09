
#include "TMV.h"
#include "TMV_Diag.h"
#include "TMV_VectorArith_Inline.h"

namespace tmv {

  template <class T, class Ta, class Tx> void MultMV(
      const T alpha, const GenDiagMatrix<Ta>& A,
      const GenVector<Tx>& x,
      const T beta, const VectorView<T>& y)
    // y = alpha * A * x + beta * y
    // yi = alpha * Ai * xi + beta * yi
  {
    TMVAssert(A.rowsize() == x.size());
    TMVAssert(A.colsize() == y.size());

    if (y.size() > 0) {
      if (alpha == T(0)) y *= beta;
      else AddElementProd(alpha,A.diag(),x,beta,y);
    }
  }

  template <class T, class Ta> void AddMM(const T alpha,
      const GenDiagMatrix<Ta>& A, 
      const T beta, const MatrixView<T>& B)
    // B = alpha * A + beta * B
  {
    TMVAssert(A.colsize() == B.colsize());
    TMVAssert(A.rowsize() == B.rowsize());

    if (B.colsize() > 0 && B.rowsize() > 0) {
      if (beta == T(0)) {
	B.Zero();
	if (alpha == T(1)) B.diag() = A.diag();
	else B.diag() = alpha * A.diag();
      } else {
	if (beta != T(1)) B *= beta;
	if (alpha == T(1)) B.diag() += A.diag();
	else B.diag() += alpha * A.diag();
      }
    }
  }

  // MJ: convert to Aptr version
  template <class T1, class T2, class T3, class T4, class T5>
    void RowMajorMultMM(const T1 alpha,
	const GenDiagMatrix<T2>& A, const GenMatrix<T3>& B,
	const T4 beta, const MatrixView<T5>& C)
    {
      // C = alpha * A * B + beta * C
      // Cij = alpha * Ai * Bij + beta Cij
      TMVAssert(B.isrm());
      TMVAssert(C.isrm());
      TMVAssert(B.rowsize() == C.rowsize());
      TMVAssert(alpha!=T1(0));
      TMVAssert(C.rowsize() > 0);
      TMVAssert(C.colsize() > 0);
      TMVAssert(C.ct() == NonConj);

      CVIter<T2> Ait = A.diag().begin();
      if (beta == T4(0)) {
	if (C.SameStorageAs(B) && !C.SameAs(B)) {
	  Matrix<T5,RowMajor> Cx(C.colsize(),C.rowsize());
	  RowMajorMultMM(alpha,A,B,T4(0),Cx.QuickView());
	  C = Cx;
	} else {
	  C = B;
	  if (alpha == T1(1))
	    for(size_t i=0;i<C.colsize();++i,++Ait) 
	      DoMultXV(*Ait,VIt<T5,Unit,NonConj>(C.row(i).begin()),
		  C.rowsize());
	  else
	    for(size_t i=0;i<C.colsize();++i,++Ait) 
	      DoMultXV(*Ait*alpha,VIt<T5,Unit,NonConj>(C.row(i).begin()),
		  C.rowsize());
	}
      } else if (beta == T4(1)) {
	if (alpha == T1(1)) 
	  for(size_t i=0;i<C.colsize();++i,++Ait) {
	    if (*Ait != T2(0)) {
	      if (B.isconj())
		DoAddVV(*Ait,CVIt<T3,Unit,Conj>(B.row(i).begin()),
		    VIt<T5,Unit,NonConj>(C.row(i).begin()),C.rowsize());
	      else
		DoAddVV(*Ait,CVIt<T3,Unit,NonConj>(B.row(i).begin()),
		    VIt<T5,Unit,NonConj>(C.row(i).begin()),C.rowsize());
	    }
	  }
	else 
	  for(size_t i=0;i<C.colsize();++i,++Ait) {
	    if (*Ait != T2(0)) {
	      if (B.isconj()) 
		DoAddVV(*Ait*alpha,CVIt<T3,Unit,Conj>(B.row(i).begin()),
		    VIt<T5,Unit,NonConj>(C.row(i).begin()),C.rowsize());
	      else
		DoAddVV(*Ait*alpha,CVIt<T3,Unit,NonConj>(B.row(i).begin()),
		    VIt<T5,Unit,NonConj>(C.row(i).begin()),C.rowsize());
	    }
	  }
      } else {
	if (C.SameStorageAs(B)) {
	  Matrix<T5,RowMajor> Cx = C;;
	  RowMajorMultMM(alpha,A,B,beta,Cx.QuickView());
	  C = Cx;
	} else if (alpha == T1(1)) {
	  for(size_t i=0;i<C.colsize();++i,++Ait) {
	    DoMultXV(beta,VIt<T5,Unit,NonConj>(C.row(i).begin()),C.rowsize());
	    if (*Ait != T2(0)) {
	      if (B.isconj())
		DoAddVV(*Ait,CVIt<T3,Unit,Conj>(B.row(i).begin()),
		    VIt<T5,Unit,NonConj>(C.row(i).begin()),C.rowsize());
	      else
		DoAddVV(*Ait,CVIt<T3,Unit,NonConj>(B.row(i).begin()),
		    VIt<T5,Unit,NonConj>(C.row(i).begin()),C.rowsize());
	    }
	  }
	} else {
	  for(size_t i=0;i<C.colsize();++i,++Ait) {
	    DoMultXV(beta,VIt<T5,Unit,NonConj>(C.row(i).begin()),C.rowsize());
	    if (*Ait != T2(0)) {
	      if (B.isconj())
		DoAddVV(*Ait*alpha,CVIt<T3,Unit,Conj>(B.row(i).begin()),
		    VIt<T5,Unit,NonConj>(C.row(i).begin()),C.rowsize());
	      else
		DoAddVV(*Ait*alpha,CVIt<T3,Unit,NonConj>(B.row(i).begin()),
		    VIt<T5,Unit,NonConj>(C.row(i).begin()),C.rowsize());
	    }
	  }
	}
      }
    }

  template <class T, class Ta, class Tb> void DoRowMajorMultMM(
      const T alpha, const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha!=T(0));
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(B.isrm());
    TMVAssert(C.isrm());
    TMVAssert(C.ct()==NonConj);

    if (IMAG(alpha) == RealType(T)(0))
      if (IMAG(beta) == RealType(T)(0))
	RowMajorMultMM(REAL(alpha),A,B,REAL(beta),C);
      else
	RowMajorMultMM(REAL(alpha),A,B,beta,C);
    else
      if (IMAG(beta) == RealType(T)(0))
	RowMajorMultMM(alpha,A,B,REAL(beta),C);
      else
	RowMajorMultMM(alpha,A,B,beta,C);
  }
      
  template <class T, class Ta, class Tb> void RowMultMM(const T alpha,
      const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    // C = alpha * A * B + beta * C
    // Cij = alpha * Ai * Bij + beta Cij
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha!=T(0));
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct()==NonConj);

    if (beta == T(0)) {
      if (C.SameStorageAs(B) && !C.SameAs(B)) {
	Matrix<T,RowMajor> C2(C.colsize(),C.rowsize());
	RowMultMM(alpha,A,B,T(0),C2.QuickView());
	C = C2;
      } else {
	C = B;
	CVIter<Ta> Ait = A.diag().begin();
	if (alpha == T(1))
	  for(size_t i=0;i<C.colsize();++i,++Ait) C.row(i) *= *Ait;
	else
	  for(size_t i=0;i<C.colsize();++i,++Ait) C.row(i) *= *Ait*alpha;
      }
    } else if (beta == T(1)) {
      CVIter<Ta> Ait = A.diag().begin();
      if (alpha == T(1)) 
	for(size_t i=0;i<C.colsize();++i,++Ait) 
	  C.row(i) += (*Ait) * B.row(i);
      else 
	for(size_t i=0;i<C.colsize();++i,++Ait) 
	  C.row(i) += (alpha * (*Ait)) * B.row(i);
    } else {
      CVIter<Ta> Ait = A.diag().begin();
      if (alpha == T(1)) 
	for(size_t i=0;i<C.colsize();++i,++Ait) 
	  C.row(i) = (*Ait) * B.row(i) + beta*C.row(i);
      else 
	for(size_t i=0;i<C.colsize();++i,++Ait) 
	  C.row(i) = (alpha * (*Ait)) * B.row(i) + beta*C.row(i);
    }
  }

  template <class T, class Ta, class Tb> void ColMultMM(const T alpha,
      const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
  {
    // C = alpha * A * B + beta * C
    // Cij = alpha * Ai * Bij + beta Cij
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(alpha!=T(0));
    TMVAssert(C.rowsize() > 0);
    TMVAssert(C.colsize() > 0);
    TMVAssert(C.ct()==NonConj);

    for(size_t j=0;j<C.rowsize();++j) {
      MultMV(alpha,A,B.col(j),beta,C.col(j));
    }
  }

  template <class T, class Ta, class Tb> void MultMM(const T alpha,
      const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const T beta, const MatrixView<T>& C)
    // C = alpha * A * B + beta * C
  {
    TMVAssert(A.rowsize() == B.colsize());
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());

    if (C.colsize() > 0 && C.rowsize() > 0) {
      if (C.isconj()) MultMM(CONJ(alpha),A.QuickConjugate(),
	  B.QuickConjugate(),CONJ(beta),C.QuickConjugate());
      else if (alpha == T(0)) C *= beta;
      else if (C.isrm() && B.isrm()) DoRowMajorMultMM(alpha,A,B,beta,C);
      else if (C.iscm() && B.iscm()) ColMultMM(alpha,A,B,beta,C);
      else if ((C.iscm() || B.iscm()) && C.colsize() >= B.colsize()) 
	ColMultMM(alpha,A,B,beta,C);
      else if ((C.isrm() || B.isrm()) && C.colsize() <= B.colsize()) 
	RowMultMM(alpha,A,B,beta,C);
      else if (C.iscm() || B.iscm()) 
	ColMultMM(alpha,A,B,beta,C);
      else RowMultMM(alpha,A,B,beta,C);
    }
  }

  template <class T, class Ta, class Tb> void MultMM(const T alpha,
	const GenDiagMatrix<Ta>& A, const GenDiagMatrix<Tb>& B,
	const T beta, const DiagMatrixView<T>& C)
    // C = alpha * A * B + beta * C
  { 
    TMVAssert(A.size() == C.size());
    TMVAssert(B.size() == C.size());
    if (C.diag().isconj())
      MultMM(CONJ(alpha),A.Conjugate(),B.Conjugate(),CONJ(beta),C.Conjugate());
    else
      AddElementProd(alpha,A.diag(),B.diag(),beta,C.diag()); 
  }

#define InstFile "TMV_DiagMatrixArith.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


