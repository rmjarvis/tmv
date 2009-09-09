
#include "TMV.h"
#include "TMV_Diag.h"

namespace tmv {

  template <class T> void SVDiv<T>::Inverse(const MatrixView<T>& minv) const
  { 
    TMVAssert(U.get() && V.get());
    if (istrans) {
      // A^-1 = (Vt S^-1 Ut)T = U* S^-1 V*
      Matrix<T,ColMajor> SinvV = V->Conjugate().Rows(0,kmax) /
	DiagMatrixViewOf(S.SubVector(0,kmax));
      minv = U->Conjugate().Cols(0,kmax) * SinvV;
    } else {
      // A^-1 = Vt S^-1 Ut
      Matrix<T,ColMajor> SinvUt = U->Adjoint().Rows(0,kmax) /
	DiagMatrixViewOf(S.SubVector(0,kmax));
      minv = V->Adjoint().Cols(0,kmax) * SinvUt;
    }
  }

  template <class T> void SVDiv<T>::DoInverseATA(
      const MatrixView<T>& minv) const
  {
    TMVAssert(V.get());
    // A = U S V
    // At = Vt S Ut
    // AtA = Vt S^2 V
    // (AtA)^-1 = Vt S^-2 V
    //
    // if istrans:
    // AT = U S V
    // At = U* S V*
    // AAt = VT S^2 V*
    // (AAt)^-1 = VT S^-2 V*
    //
    Matrix<T,ColMajor> SinvV = V->Rows(0,kmax) /
      DiagMatrixViewOf(S.SubVector(0,kmax));
    if (istrans)
      minv = SinvV.Transpose() * SinvV.Conjugate();
    else
      minv = SinvV.Adjoint() * SinvV;
  }

#define InstFile "TMV_SVDiv_C.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


