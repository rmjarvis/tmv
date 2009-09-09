
#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_Householder.h"
#include "TMV_Diag.h"

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

    size_t endcol = QRx.nlo()+1;
    size_t endrow = QRx.nhi()+1;
    for(size_t j=0;j<QRx.rowsize();++j) {
      // Apply the Householder Reflection for this column
      Qbeta(j) = Householder_Reflect(QRx.SubMatrix(j,endcol,j,endrow),det);
      if (endcol < QRx.colsize()) ++endcol;
      if (endrow < QRx.rowsize()) ++endrow;
    }
  }

  template <class T> inline Matrix<T,ColMajor> GetQFromBandQR(
      const GenBandMatrix<T>& QRx, const GenVector<T>& Qbeta) 
  {
    // Extract the Q matrix from a combined QRx matrix
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    const size_t N = QRx.rowsize();
    const size_t M = QRx.colsize();
    Matrix<T,ColMajor> Q = BandMatrixViewOf(QRx,QRx.nlo(),0);

    for(int j=N-1;j>=0;j--) {
      Householder_Unpack(Q.SubMatrix(j,M,j,N),Qbeta(j));
    }
    return Q;
  }

  //
  // Packed BandQ - LDivEq/RDivEq
  //

  template <class T1, class T2> void BandQ_LDivEq(
      const GenBandMatrix<T1>& Q, const GenVector<T1>& Qbeta,
      const MatrixView<T2>& m)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Q.rowsize() == Qbeta.size());
    TMVAssert(Q.colsize() == m.colsize());

    size_t i2 = Q.nlo()+1;
    for(size_t j=0,i1=1;j<Q.rowsize();++j,++i1) {
      if (Qbeta(j) != T1(0)) 
	Householder_Mult(Q.col(j,i1,i2),Qbeta(j),m.Rows(j,i2));
      if (i2<Q.colsize()) ++i2;
    }
  }

  template <class T1, class T2> void BandQ_RDivEq(
      const GenBandMatrix<T1>& Q, const GenVector<T1>& Qbeta,
      const MatrixView<T2>& m)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Qbeta.size() == Q.rowsize());
    TMVAssert(m.rowsize() == Q.colsize());

    const size_t N = Q.rowsize();
    size_t i1 = N;
    size_t i2 = Q.IsSquare() ? N : min(N+Q.nlo(),Q.colsize());
    size_t k=Q.IsSquare() ? Q.nlo() : N+Q.nlo()-i2;
    for(size_t j=N-1;i1>0;--j,--i1) {
      if (Qbeta(j) != T1(0)) 
	Householder_ConjMult(Q.col(j,i1,i2),CONJ(Qbeta(j)),
	    m.Transpose().Rows(j,i2));
      if (k>0) --k; else --i2;
    }
  }

  //
  // LDiv
  //

  template <class T1, class T2, class T3> void BandQR_LDiv(
      const GenBandMatrix<T1>& QRx, const GenVector<T1>& Qbeta,
      const GenMatrix<T2>& m, const MatrixView<T3>& x)
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(Qbeta.size() == QRx.rowsize());
    TMVAssert(m.colsize() == QRx.colsize());
    TMVAssert(x.colsize() == QRx.rowsize());
    TMVAssert(x.rowsize() == m.rowsize());

    const size_t N = QRx.rowsize();

    if (QRx.IsSquare()) {
      x = m;
      BandQ_LDivEq(QRx,Qbeta,x);
    } else {
      if (m.isrm()) {
	Matrix<T3,RowMajor> m1 = m;
	BandQ_LDivEq(QRx,Qbeta,m1.QuickView());
	x = m1.Rows(0,N);
      } else {
	Matrix<T3,ColMajor> m1 = m;
	BandQ_LDivEq(QRx,Qbeta,m1.QuickView());
	x = m1.Rows(0,N);
      }
    }

    BandTriLDivEq(QRx.SubBandMatrix(0,N,0,N,0,QRx.nhi()),x,NonUnitDiag);
  }

  //
  // LDivEq
  //

  template <class T1, class T2> void BandQR_LDivEq(
      const GenBandMatrix<T1>& QRx, const GenVector<T1>& Qbeta, 
      const MatrixView<T2>& m)
  {
    TMVAssert(QRx.colsize() == QRx.rowsize());
    TMVAssert(Qbeta.size() == QRx.rowsize());
    TMVAssert(m.colsize() == QRx.colsize());

    BandQ_LDivEq(QRx,Qbeta,m);
    BandTriLDivEq(BandMatrixViewOf(QRx,0,QRx.nhi()),m,NonUnitDiag);
  }

  //
  // RDiv
  //

  template <class T1, class T2, class T3> void BandQR_RDiv(
      const GenBandMatrix<T1>& QRx, const GenVector<T1>& Qbeta, 
      const GenMatrix<T2>& m, const MatrixView<T3>& x)
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    TMVAssert(Qbeta.size() == QRx.rowsize());
    TMVAssert(x.rowsize() == QRx.colsize());
    TMVAssert(m.rowsize() == QRx.rowsize());
    TMVAssert(x.colsize() == m.colsize());

    const size_t M = QRx.colsize();
    const size_t N = QRx.rowsize();

    x.Cols(N,M).Zero();
    x.Cols(0,N) = m;
    BandTriLDivEq(QRx.SubBandMatrix(0,N,0,N,0,QRx.nhi()).Transpose(),
	x.Cols(0,N).Transpose(),NonUnitDiag);

    BandQ_RDivEq(QRx,Qbeta,x);
  }

  //
  // RDivEq
  //

  template <class T1, class T2> void BandQR_RDivEq(
      const GenBandMatrix<T1>& QRx, const GenVector<T1>& Qbeta, 
      const MatrixView<T2>& m)
  {
    TMVAssert(QRx.colsize() == QRx.rowsize());
    TMVAssert(Qbeta.size() == QRx.rowsize());
    TMVAssert(m.rowsize() == QRx.colsize());

    // Solve x Q R = m in place (m <- x)
    BandTriLDivEq(BandMatrixViewOf(QRx,0,QRx.nhi()).Transpose(),
	m.Transpose(),NonUnitDiag);
    BandQ_RDivEq(QRx,Qbeta,m);
  }

  template <class T> T BandQRDiv<T>::Det() const
  {
    if (!donedet) {
      det *= DiagMatrixViewOf(QRx.diag()).Det();
      donedet = true;
    }
    return det;
  }

  template <class T> Matrix<T,ColMajor> BandQRDiv<T>::Inverse() const
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    if (istrans) {
      Matrix<T,ColMajor> inv(QRx.colsize(),QRx.rowsize());
      LDiv(Eye<RealType(T),ColMajor>(QRx.rowsize()),inv.QuickView());
      return inv;
    } else {
      Matrix<T,ColMajor> inv(QRx.rowsize(),QRx.colsize());
      RDiv(Eye<RealType(T),ColMajor>(QRx.rowsize()),inv.QuickView());
      return inv;
    }
  }

  template <class T> Matrix<T,ColMajor> BandQRDiv<T>::DoInverseATA() const
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    const size_t N = QRx.rowsize();
    // At A = Rt R
    // (At A)^-1 = (Rt R)^-1 = R^-1 Rt^-1
    Matrix<T,ColMajor> Rinv(N,N);
    Rinv.SetToIdentity();

    BandTriLDivEq(QRx.SubBandMatrix(0,N,0,N,0,QRx.nhi()),
	Rinv.QuickView(),NonUnitDiag);
    return Rinv * Rinv.Transpose();
  }

  template <class T> Matrix<T,ColMajor> BandQRDiv<T>::QR_GetQ() const 
  {
    // Extract the Q matrix from a combined QRx matrix
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    return GetQFromBandQR(QRx,Qbeta);
  }

#define InstFile "TMV_BandQRDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


