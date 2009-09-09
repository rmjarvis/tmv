
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
    for(size_t i=0;i<A.colsize();++i) {
      DoMultXV(alpha,VIt<T,Unit,NonConj>(A.row(i,j1,j2).begin()),len);
      if (k>0) {--k;++len;} else ++j1;
      if (j2<A.rowsize()) ++j2;
      else if (j1==A.rowsize()) break;
      else --len;
    }
  }

  template <class T, class T1> void ColMajorMultXM(
      const T1 alpha, const BandMatrixView<T>& A)
  {
    TMVAssert(alpha != T1(1));
    TMVAssert(alpha != T1(0));
    TMVAssert(A.iscm());
    TMVAssert(A.colsize() > 0);
    TMVAssert(A.rowsize() > 0);
    TMVAssert(A.ct() == NonConj);

    size_t i1=0;
    size_t i2=A.nlo()+1;
    size_t k=A.nhi();
    size_t len = A.nlo()+1;
    for(size_t j=0;j<A.rowsize();++j) {
      DoMultXV(alpha,VIt<T,Unit,NonConj>(A.col(j,i1,i2).begin()),len);
      if (k>0) {--k;++len;} else ++i1;
      if (i2<A.colsize()) ++i2;
      else if (i1==A.colsize()) break;
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
    for(int i=-A.nlo();i<=A.nhi();++i) {
      DoMultXV(alpha,VIt<T,Unit,NonConj>(A.diag(i).begin()),len);
      if (i<0) { if (j2<A.rowsize()) { ++len; ++j2; } }
      else { if (j2<A.rowsize()) ++j2; else --len; }
    }
  }

  template <class T> void MultXM(const T alpha, const BandMatrixView<T>& A)
    // A = alpha * A
  {
    if (A.rowsize() > 0 && A.colsize() > 0 && alpha != T(1)) {
      if (alpha == T(0)) A.Zero();
      else {
	if (A.isconj()) MultXM(CONJ(alpha),A.QuickConjugate());
	else if (A.isrm()) {
	  if (IMAG(alpha)==RealType(T)(0)) RowMajorMultXM(REAL(alpha),A);
	  else RowMajorMultXM(alpha,A);
	} else if (A.iscm()) {
	  if (IMAG(alpha)==RealType(T)(0)) ColMajorMultXM(REAL(alpha),A);
	  else ColMajorMultXM(alpha,A);
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


