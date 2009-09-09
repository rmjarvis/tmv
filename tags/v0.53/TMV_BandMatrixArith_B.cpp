
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

#define InstFile "TMV_BandMatrixArith_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


