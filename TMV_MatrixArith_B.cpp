
#include "TMV_VectorArith_Inline.h"
#include "TMV.h"

namespace tmv {

  //
  // MultXM
  //

  template <class T> inline void RowMajorMultXM(const T alpha,
      const MatrixView<T>& A)
  {
    TMVAssert(A.isrm());
    TMVAssert(A.ct() == NonConj);
    TMVAssert(A.colsize() > 0 && A.rowsize() > 0);
    TMVAssert(alpha != T(1));
    TMVAssert(alpha != T(0));

    T* Aptr = A.ptr();
    for(size_t i=0;i<A.colsize();++i,Aptr+=A.stepi()) 
      // A.row(i) *= alpha;
      if (IMAG(alpha) == RealType(T)(0)) 
	DoMultXV(REAL(alpha),VIt<T,Unit,NonConj>(Aptr,1
	      FIRSTLAST1(A.first,A.last) ),A.rowsize());
      else 
	DoMultXV(alpha,VIt<T,Unit,NonConj>(Aptr,1
	      FIRSTLAST1(A.first,A.last) ),A.rowsize());
  }

  template <class T> inline void MultXM(const T alpha, const MatrixView<T>& A)
    // A = alpha * A
  {
    if (A.isconj()) MultXM(CONJ(alpha),A.Conjugate());
    else if (A.colsize() > 0 && A.rowsize() > 0 && alpha != T(1)) {
      if (alpha == T(0)) A.Zero();
      else if ( (A.stepj()==int(1) && A.stepi()==int(A.rowsize())) || 
	  (A.stepi()==int(1) && A.stepj()==int(A.colsize())) || 
	  (A.stepj()<=A.stepi() && A.stepi()==A.stepj()*int(A.rowsize())) ||
	  (A.stepi()<=A.stepj() && A.stepj()==A.stepi()*int(A.colsize())) ) {
	A.LinearView() *= alpha;
      } else {
	if (A.isrm()) RowMajorMultXM(alpha,A);
	else if (A.iscm()) RowMajorMultXM(alpha,A.Transpose());
	else if (A.colsize() < A.rowsize())
	  for(size_t i=0;i<A.colsize();++i) A.row(i) *= alpha;
	else 
	  for(size_t j=0;j<A.colsize();++j) A.col(j) *= alpha;
      }
    }
  }

#define InstFile "TMV_MatrixArith_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


