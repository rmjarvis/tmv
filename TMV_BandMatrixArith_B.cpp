
#include "TMV.h"
#include "TMV_Band.h"

namespace tmv {

  //
  // MultXM
  //

  template <class T> void MultXM(const T alpha, const BandMatrixView<T>& A)
    // A = alpha * A
  {
    if (A.rowsize() > 0 && A.colsize() > 0 && alpha != T(1)) {
      if (A.isconj()) MultXM(CONJ(alpha),A.Conjugate());
      else if (alpha == T(0)) A.Zero();
      else if (A.CanLinearize()) A.LinearView() *= alpha;
      else 
	for(int i=-A.nlo();i<=A.nhi();++i) A.diag(i) *= alpha;
    }
  }

  template <class T, class Ta> void ElementProd(
      const T alpha, const GenBandMatrix<Ta>& A, const BandMatrixView<T>& B)
  {
    TMVAssert(A.colsize() == B.colsize());
    TMVAssert(A.rowsize() == B.rowsize());
    TMVAssert(A.nlo() == B.nlo());
    TMVAssert(A.nhi() == B.nhi());
    if (A.stor() == B.stor() && A.CanLinearize() && B.CanLinearize()) {
      TMVAssert(A.stepi() == B.stepi() && A.stepj() == B.stepj());
      ElementProd(alpha,A.ConstLinearView(),B.LinearView());
    } else {
      for(int i=-A.nlo();i<=A.nhi();++i) 
	ElementProd(alpha,A.diag(i),B.diag(i));
    }
  }

  template <class T, class Ta, class Tb> void AddElementProd(
      const T alpha, const GenBandMatrix<Ta>& A, const GenBandMatrix<Tb>& B,
      const BandMatrixView<T>& C)
  {
    TMVAssert(A.colsize() == C.colsize());
    TMVAssert(A.rowsize() == C.rowsize());
    TMVAssert(B.colsize() == C.colsize());
    TMVAssert(B.rowsize() == C.rowsize());
    TMVAssert(A.nlo() == C.nlo());
    TMVAssert(A.nhi() == C.nhi());
    TMVAssert(B.nlo() == C.nlo());
    TMVAssert(B.nhi() == C.nhi());
    if (A.stor() == C.stor() && B.stor() == C.stor() &&
	A.CanLinearize() && B.CanLinearize() && C.CanLinearize()) {
      TMVAssert(A.stepi() == C.stepi() && A.stepj() == C.stepj());
      TMVAssert(B.stepi() == C.stepi() && B.stepj() == C.stepj());
      AddElementProd(alpha,A.ConstLinearView(),B.ConstLinearView(),
	  C.LinearView());
    } else {
      for(int i=-A.nlo();i<=A.nhi();++i) 
	AddElementProd(alpha,A.diag(i),B.diag(i),C.diag(i));
    }
  }

#define InstFile "TMV_BandMatrixArith_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


