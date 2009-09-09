
#include "TMV_VectorArith_Inline.h"
#include "TMV.h"
#include "TMV_Tri.h"

namespace tmv {

  //
  // MultXM
  //

  template <class T, class T1> void RowMajorMultXM(
      const T1 alpha, const UpperTriMatrixView<T>& A)
  {
    TMVAssert(A.isrm());
    TMVAssert(!A.isunit());
    TMVAssert(alpha != T1(1));
    TMVAssert(A.size() > 0);
    TMVAssert(A.ct() == NonConj);

    const int ds = A.stepi()+1;
    T* Aptr = A.ptr();
    for(size_t i=0,len=A.rowsize();len>0;++i,--len,Aptr+=ds) 
      DoMultXV(alpha,VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),
	  len);
  }

  template <class T, class T1> void ColMajorMultXM(
      const T1 alpha, const UpperTriMatrixView<T>& A)
  {
    TMVAssert(A.iscm());
    TMVAssert(!A.isunit());
    TMVAssert(alpha != T1(1));
    TMVAssert(A.size() > 0);
    TMVAssert(A.ct() == NonConj);

    T* Aptr = A.ptr();
    for(size_t j=0;j<A.rowsize();++j,Aptr+=A.stepj()) 
      DoMultXV(alpha,VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),
	  j+1);
  }

  template <class T> void MultXM(
      const T alpha, const UpperTriMatrixView<T>& A)
    // A = alpha * A
  {
    TMVAssert(!A.isunit());

    if (A.size() > 0 && alpha != T(1)) {
      if (alpha == T(0)) A.Zero();
      else {
	if (A.isconj()) MultXM(CONJ(alpha),A.Conjugate());
	else if (A.isrm())
	  if (IMAG(alpha) == RealType(T)(0))
	    RowMajorMultXM(REAL(alpha),A);
	  else
	    RowMajorMultXM(alpha,A);
	else if (A.iscm())
	  if (IMAG(alpha) == RealType(T)(0))
	    ColMajorMultXM(REAL(alpha),A);
	  else
	    ColMajorMultXM(alpha,A);
	else 
	  for(size_t i=0;i<A.colsize();++i) 
	    MultXV(alpha,A.row(i,i,A.rowsize()));
      }
    }
  }

#define InstFile "TMV_TriMatrixArith_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv

