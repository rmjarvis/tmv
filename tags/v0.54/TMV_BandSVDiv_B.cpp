
#include "TMV_Band.h"
#include "TMV_Diag.h"


namespace tmv {

  template <class T> void BandSVDiv<T>::Inverse(
      const MatrixView<T>& minv) const
  { 
    TMVAssert(U.get() && V.get());
    if (istrans) {
      Matrix<T,ColMajor> SinvV = V->Conjugate().Rows(0,kmax) /
	DiagMatrixViewOf(S.SubVector(0,kmax));
      minv = U->Conjugate().Cols(0,kmax) * SinvV;
    } else  {
      Matrix<T,ColMajor> SinvUt = U->Adjoint().Rows(0,kmax) /
	DiagMatrixViewOf(S.SubVector(0,kmax));
      minv = V->Adjoint().Cols(0,kmax) * SinvUt;
    }
  }

  template <class T> void BandSVDiv<T>::DoInverseATA(
      const MatrixView<T>& minv) const
  {
    TMVAssert(V.get());
    Matrix<T,ColMajor> SinvV = V->Rows(0,kmax) /
      DiagMatrixViewOf(S.SubVector(0,kmax));
    if (istrans)
      minv = SinvV.Transpose() * SinvV.Conjugate();
    else
      minv = SinvV.Adjoint() * SinvV;
  }

#define InstFile "TMV_BandSVDiv_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


