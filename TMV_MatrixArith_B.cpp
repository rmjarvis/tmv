
#include "TMV_VectorArith_Inline.h"
#include "TMV.h"

namespace tmv {

  //
  // MultXM
  //

  template <class T> inline void MultXM(
      const T alpha, const MatrixView<T>& A)
    // A = alpha * A
  {
    if (A.isconj()) MultXM(CONJ(alpha),A.Conjugate());
    else if (A.colsize() > 0 && A.rowsize() > 0 && alpha != T(1)) {
      if (alpha == T(0)) A.Zero();
      else if ( (A.stepj()==int(1) && A.stepi()==int(A.rowsize())) || 
	  (A.stepi()==int(1) && A.stepj()==int(A.colsize())) || 
	  (A.stepj()<=A.stepi() && A.stepi()==A.stepj()*int(A.rowsize())) ||
	  (A.stepi()<=A.stepj() && A.stepj()==A.stepi()*int(A.colsize())) ) {
	MultXV(alpha,A.LinearView());
      } else {
	if (A.isrm()) {
	  T* Aptr = A.ptr();
	  if (IMAG(alpha) == RealType(T)(0)) 
	    for(size_t i=0;i<A.colsize();++i,Aptr+=A.stepi()) 
	      DoMultXV(REAL(alpha),VIt<T,Unit,NonConj>(Aptr,1
		    FIRSTLAST1(A.first,A.last) ),A.rowsize());
	  else 
	    for(size_t i=0;i<A.colsize();++i,Aptr+=A.stepi()) 
	      DoMultXV(alpha,VIt<T,Unit,NonConj>(Aptr,1
		    FIRSTLAST1(A.first,A.last) ),A.rowsize());
	}
	else if (A.iscm()) {
	  T* Aptr = A.ptr();
	  if (IMAG(alpha) == RealType(T)(0)) 
	    for(size_t j=0;j<A.rowsize();++j,Aptr+=A.stepj()) 
	      DoMultXV(REAL(alpha),VIt<T,Unit,NonConj>(Aptr,1
		    FIRSTLAST1(A.first,A.last) ),A.colsize());
	  else 
	    for(size_t j=0;j<A.rowsize();++j,Aptr+=A.stepj()) 
	      DoMultXV(alpha,VIt<T,Unit,NonConj>(Aptr,1
		    FIRSTLAST1(A.first,A.last) ),A.colsize());
	}
	else if (A.colsize() < A.rowsize()) {
	  for(size_t i=0;i<A.colsize();++i) MultXV(alpha,A.row(i));
	}
	else {
	  for(size_t i=0;i<A.colsize();++i) MultXV(alpha,A.col(i));
	}
      }
    }
  }

#define InstFile "TMV_MatrixArith_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


