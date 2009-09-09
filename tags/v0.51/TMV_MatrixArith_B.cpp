
#include "TMV.h"

namespace tmv {

  //
  // MultXM
  //

  template <class T, class Ta> void RowMajorMultXM(
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

#define InstFile "TMV_MatrixArith_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


