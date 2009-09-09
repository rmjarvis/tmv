
#include "TMV_Band.h"

namespace tmv {

  template <class T> void BandLUDiv<T>::Inverse(
      const MatrixView<T>& minv) const
  {
    TMVAssert(minv.colsize() == LUx.colsize());
    TMVAssert(minv.rowsize() == LUx.colsize());
    minv.SetToIdentity();
    LDivEq(minv);
  }

  template <class T> void BandLUDiv<T>::InverseATA(
      const MatrixView<T>& minv) const
  {
    Matrix<T,ColMajor> temp(minv.colsize(),minv.rowsize());
    Inverse(temp.View());
    minv = temp*Transpose(temp);
  }

#define InstFile "TMV_BandLUDiv_C.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


