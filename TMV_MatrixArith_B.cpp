
#include "TMV.h"

namespace tmv {

  //
  // MultXM
  //

  template <class T, class Ta> inline void RowMajorMultXM(
      const Ta alpha, const MatrixView<T>& A)
  {
    TMVAssert(A.isrm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.colsize() > 0 && A.rowsize() > 0);
    TMVAssert(alpha != Ta(1));
    TMVAssert(alpha != Ta(0));

    const size_t M = A.colsize();
    const size_t N = A.rowsize();

    T* Ai0 = A.ptr();

    for(size_t i=M;i>0;--i,Ai0+=A.stepi()) {
      // A.row(i) *= alpha;
      T* Aij = Ai0;
      for(size_t j=N;j>0;--j,++Aij) *Aij *= alpha;
    }
  }

  template <class T> void MultXM(const T alpha, const MatrixView<T>& A)
    // A = alpha * A
  {
    if (A.colsize() > 0 && A.rowsize() > 0 && alpha != T(1)) {
      if (A.isconj()) MultXM(CONJ(alpha),A.Conjugate());
      else if (alpha == T(0)) A.Zero();
      else if (A.CanLinearize()) A.LinearView() *= alpha;
      else if (A.isrm())
	if (IsComplex(T()) && IMAG(alpha) == RealType(T)(0))
	  RowMajorMultXM(REAL(alpha),A);
	else
	  RowMajorMultXM(alpha,A);
      else if (A.iscm())
	if (IsComplex(T()) && IMAG(alpha) == RealType(T)(0))
	  RowMajorMultXM(REAL(alpha),A.Transpose());
	else 
	  RowMajorMultXM(alpha,A.Transpose());
      else 
	if (A.colsize() < A.rowsize())
	  for(size_t i=0;i<A.colsize();++i) A.row(i) *= alpha;
	else 
	  for(size_t j=0;j<A.colsize();++j) A.col(j) *= alpha;
    }
  }

  template <class T, class Ta> void ElementProd(
      const T alpha, const GenMatrix<Ta>& A, const MatrixView<T>& B)
  {
    TMVAssert(A.colsize() == B.colsize());
    TMVAssert(A.rowsize() == B.rowsize());
    if (A.stor() == B.stor() && A.CanLinearize() && B.CanLinearize()) {
      TMVAssert(A.stepi() == B.stepi() && A.stepj() == B.stepj());
      ElementProd(alpha,A.ConstLinearView(),B.LinearView());
    } else if (B.isrm()) {
      for(size_t i=0;i<B.colsize();i++)
	ElementProd(alpha,A.row(i),B.row(i));
    } else {
      for(size_t j=0;j<B.rowsize();j++)
	ElementProd(alpha,A.col(j),B.col(j));
    }
  }

  template <class T, class Ta, class Tb> void AddElementProd(
      const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
      const MatrixView<T>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == C.rowsize());
    TMVAssert(B.colsize() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    if (A.stor() == C.stor() && B.stor() == C.stor() &&
	A.CanLinearize() && B.CanLinearize() && C.CanLinearize()) {
      TMVAssert(A.stepi() == C.stepi() && A.stepj() == C.stepj());
      TMVAssert(B.stepi() == C.stepi() && B.stepj() == C.stepj());
      AddElementProd(alpha,A.ConstLinearView(),B.ConstLinearView(),
	  C.LinearView());
    } else if (C.isrm()) {
      for(size_t i=0;i<C.colsize();i++)
	AddElementProd(alpha,A.row(i),B.row(i),C.row(i));
    } else {
      for(size_t j=0;j<C.rowsize();j++)
	AddElementProd(alpha,A.col(j),B.col(j),C.col(j));
    }
  }

#define InstFile "TMV_MatrixArith_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


