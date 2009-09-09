
#include "TMV.h"
#include "TMV_Band.h"

namespace tmv {

  template <class T> void BandQRDiv<T>::Inverse(const MatrixView<T>& minv) const
  {
    MatrixView<T> minv2 = istrans ? minv.Transpose() : minv;

    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(minv2.colsize() == QRx.rowsize());
    TMVAssert(minv2.rowsize() == QRx.colsize());
    const size_t M = QRx.colsize();
    const size_t N = QRx.rowsize();
    minv2.Cols(N,M).Zero();
    minv2.Cols(0,N).SetToIdentity();
    BandTriLDivEq(QRx.SubBandMatrix(0,N,0,N,0,QRx.nhi()).Transpose(),
	minv2.Cols(0,N).Transpose(),NonUnitDiag);
    BandQ_RDivEq(QRx,Qbeta,minv2);
  }

  template <class T> void BandQRDiv<T>::DoInverseATA(
      const MatrixView<T>& minv) const
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    const size_t N = QRx.rowsize();
    // At A = Rt R
    // (At A)^-1 = (Rt R)^-1 = R^-1 Rt^-1
    Matrix<T,ColMajor> Rinv(N,N);
    Rinv.SetToIdentity();

    BandTriLDivEq(QRx.SubBandMatrix(0,N,0,N,0,QRx.nhi()),
	Rinv.View(),NonUnitDiag);
    minv = Rinv * Rinv.Transpose();
  }

#define InstFile "TMV_BandQRDiv_C.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


