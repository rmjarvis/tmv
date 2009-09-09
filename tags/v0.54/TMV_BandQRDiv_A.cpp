
#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_Householder.h"

//#define XDEBUG

namespace tmv {

  template <class T> void BandQR_Decompose(
      const BandMatrixView<T>& QRx,
      const VectorView<T>& Qbeta, T& det)
  {
    // Decompose A (input as QRx) into A = Q R 
    // where Q is unitary, and R is upper triangular
    // Q and R are stored in the same matrix (QRx)
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(Qbeta.size() == QRx.rowsize());

#ifdef XDEBUG
    Matrix<T> A0 = QRx;
#endif
    size_t endcol = QRx.nlo()+1;
    size_t endrow = QRx.nhi()+1;
    for(size_t j=0;j<QRx.rowsize();++j) {
      // Apply the Householder Reflection for this column
      Qbeta(j) = Householder_Reflect(QRx.SubMatrix(j,endcol,j,endrow),det);
      if (endcol < QRx.colsize()) ++endcol;
      if (endrow < QRx.rowsize()) ++endrow;
    }
#ifdef XDEBUG
    Matrix<T> Q = GetQFromBandQR(QRx,Qbeta);
    Matrix<T> R = QRx.Diags(0,QRx.nhi()+1);
    Matrix<T> QR = Q*R;
    if (Norm(QR-A0) > 0.00001*Norm(A0)) {
      cerr<<"BandQR_Decompose: \n";
      cerr<<"A = "<<Type(QRx)<<"  "<<A0<<endl;
      cerr<<"QRx = "<<QRx<<endl;
      cerr<<"Q = "<<Q<<endl;
      cerr<<"R = "<<R<<endl;
      cerr<<"QR = "<<QR<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_BandQRDiv_A.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


