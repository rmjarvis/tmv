
#include "TMV_VectorArith_Inline.h"
#include "TMV.h"
#include "TMV_Band.h"

namespace tmv {

  //
  // MultXM
  //

  template <class T, class T1> void RowMajorMultXM(
      const T1 alpha, const BandMatrixView<T>& A)
  {
    TMVAssert(alpha != T1(1));
    TMVAssert(alpha != T1(0));
    TMVAssert(A.isrm());
    TMVAssert(A.colsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);

    size_t j1=0;
    size_t j2=A.nhi()+1;
    size_t k=A.nlo();
    size_t len = A.nhi()+1;
    T* Aptr = A.ptr();
    for(size_t i=0;i<A.colsize();++i) {
      // A.row(i,j1,j2) *= alpha;
      DoMultXV(alpha,VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),
	  len);
      if (k>0) {--k; ++len; Aptr+=A.stepi(); } 
      else { ++j1; Aptr+=A.diagstep(); }
      if (j2<A.rowsize()) ++j2;
      else if (j1==A.rowsize()) break;
      else --len;
    }
  }

  template <class T, class T1> void DiagMajorMultXM(
      const T1 alpha, const BandMatrixView<T>& A)
  {
    TMVAssert(alpha != T1(1));
    TMVAssert(alpha != T1(0));
    TMVAssert(A.isdm());
    TMVAssert(A.colsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);

    size_t len = min(A.colsize()-A.nlo(),A.rowsize());
    size_t j2 = len;
    T* Aptr = A.ptr()+A.nlo()*A.stepi();
    for(int i=-A.nlo();i<=A.nhi();++i) {
      // A.diag(i) *= alpha;
      DoMultXV(alpha,VIt<T,Unit,NonConj>(Aptr,1 FIRSTLAST1(A.first,A.last) ),
	  len);
      if (i<0) { 
	Aptr -= A.stepi();
	if (j2<A.rowsize()) { ++len; ++j2; } 
      }
      else { 
	Aptr += A.stepj();
	if (j2<A.rowsize()) ++j2; 
	else --len; 
      }
    }
  }

  template <class T> void MultXM(const T alpha, const BandMatrixView<T>& A)
    // A = alpha * A
  {
    if (A.rowsize() > 0 && A.colsize() > 0 && alpha != T(1)) {
      if (alpha == T(0)) A.Zero();
      else {
	if (A.isconj()) MultXM(CONJ(alpha),A.Conjugate());
	else if (A.isrm()) {
	  if (IMAG(alpha)==RealType(T)(0)) RowMajorMultXM(REAL(alpha),A);
	  else RowMajorMultXM(alpha,A);
	} else if (A.iscm()) {
	  if (IMAG(alpha)==RealType(T)(0)) 
	    RowMajorMultXM(REAL(alpha),A.Transpose());
	  else RowMajorMultXM(alpha,A.Transpose());
	} else if (A.isdm()) {
	  if (IMAG(alpha)==RealType(T)(0)) DiagMajorMultXM(REAL(alpha),A);
	  else DiagMajorMultXM(alpha,A);
	} else {
	  for(int i=-A.nlo();i<=A.nhi();++i) A.diag(i) *= alpha;
	}
      }
    }
  }

#define InstFile "TMV_BandMatrixArith_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


