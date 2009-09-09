
#include "TMV.h"
#include "TMV_Tri.h"

namespace tmv {

  template <class T> void QRPDiv<T>::Inverse(const MatrixView<T>& minv) const
  {
    MatrixView<T> minv2 = istrans ? minv.Transpose() : minv;
    // minv2 = I Pt R^-1 Q^-1
    TMVAssert(minv2.colsize() == QRx.rowsize());
    TMVAssert(minv2.rowsize() == QRx.colsize());

    const size_t N = QRx.rowsize();
    const size_t M = QRx.colsize();

    minv2.Cols(0,N).SetToIdentity();
    minv2.Cols(0,N).PermuteCols(P.get());
    minv2.Cols(N1,M).Zero();
    minv2.Cols(0,N1) %= UpperTriMatrixViewOf(QRx).SubTriMatrix(0,N1);
    Q_RDivEq(QRx,beta,minv2);
  }

  template <class T> void QRPDiv<T>::DoInverseATA(
      const MatrixView<T>& ata) const
  {
    // At A = Pt Rt R P
    // (At A)^-1 = (Pt Rt R P)^-1 = Pt R^-1 R^-1t P
    UpperTriMatrixView<T> rinv = UpperTriMatrixViewOf(ata);
    rinv = UpperTriMatrixViewOf(QRx).Inverse();
    ata = rinv * rinv.Adjoint();
    ata.ReversePermuteRows(P.get());
    ata.ReversePermuteCols(P.get());
  }

#define InstFile "TMV_QRPDiv_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


