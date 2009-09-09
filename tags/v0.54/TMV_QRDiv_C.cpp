
#include "TMV.h"
#include "TMV_Tri.h"

namespace tmv {

  template <class T> void QRDiv<T>::Inverse(const MatrixView<T>& minv) const
  {
    MatrixView<T> minv2 = istrans ? minv.Transpose() : minv;
    // minv = R^-1 Q^-1
    TMVAssert(minv2.colsize() == QRx.rowsize());
    TMVAssert(minv2.rowsize() == QRx.colsize());

    const size_t N = QRx.rowsize();
    const size_t M = QRx.colsize();

    minv2.Cols(N,M).Zero();
    minv2.Cols(0,N) = UpperTriMatrixViewOf(QRx).Inverse();
    Q_RDivEq(QRx,beta,minv2);
  }

  template <class T> void QRDiv<T>::DoInverseATA(
      const MatrixView<T>& ata) const
  {
    // At A = Rt Qt Q R = Rt R
    // (At A)^-1 = (Rt R)^-1 = R^-1 * Rt^-1

    UpperTriMatrixView<T> rinv = UpperTriMatrixViewOf(ata);
    rinv = UpperTriMatrixViewOf(QRx).Inverse();
    ata = rinv * rinv.Adjoint();
  }

  template <class T> Matrix<T> QRDiv<T>::GetQ() const
  {
    Matrix<T> Q = QRx;
    GetQFromQR(Q.View(),beta);
    return Q;
  }

#define InstFile "TMV_QRDiv_C.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


