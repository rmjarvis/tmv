
#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_Householder.h"

//#define XDEBUG

namespace tmv {

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

    if (Q.nlo() > 0) {
      size_t i2 = Q.nlo()+1;
      for(size_t j=0,i1=1;j<Q.rowsize();++j,++i1) {
	if (Qbeta(j) != T1(0)) 
	  Householder_LMult(Q.col(j,i1,i2),Qbeta(j),m.Rows(j,i2));
	if (i2<Q.colsize()) ++i2;
      }
    }
  }

  template <class T1, class T2> void BandQ_RDivEq(
      const GenBandMatrix<T1>& Q, const GenVector<T1>& Qbeta,
      const MatrixView<T2>& m)
  {
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Qbeta.size() == Q.rowsize());
    TMVAssert(m.rowsize() == Q.colsize());

    if (Q.nlo() > 0) {
      const size_t N = Q.rowsize();
      size_t i1 = N;
      size_t i2 = Q.IsSquare() ? N : min(N+Q.nlo(),Q.colsize());
      size_t k=Q.IsSquare() ? Q.nlo() : N+Q.nlo()-i2;
      for(size_t j=N-1;i1>0;--j,--i1) {
	if (Qbeta(j) != T1(0)) 
	  Householder_LMult(Q.col(j,i1,i2).Conjugate(),Qbeta(j),
	      m.Cols(j,i2).Transpose());
	if (k>0) --k; else --i2;
      }
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

#ifdef XDEBUG
    Matrix<T1> QR2 = GetQFromBandQR(QRx,Qbeta)*QRx.Diags(0,QRx.nhi()+1);
    QR2.DivideUsing(tmv::QR);
    Matrix<T3> x2 = m / QR2;
#endif
    const size_t N = QRx.rowsize();

    if (QRx.IsSquare()) {
      x = m;
      BandQ_LDivEq(QRx,Qbeta,x);
    } else if (QRx.nlo() > 0) {
      if (m.isrm()) {
	Matrix<T3,RowMajor> m1 = m;
	BandQ_LDivEq(QRx,Qbeta,m1.View());
	x = m1.Rows(0,N);
      } else {
	Matrix<T3,ColMajor> m1 = m;
	BandQ_LDivEq(QRx,Qbeta,m1.View());
	x = m1.Rows(0,N);
      }
    } else {
      x = m.Rows(0,N);
    }

    BandTriLDivEq(QRx.SubBandMatrix(0,N,0,N,0,QRx.nhi()),x,NonUnitDiag);

#ifdef XDEBUG
    if (Norm(x2-x) > 0.001*Norm(x)) {
      cerr<<"BandQR_LDiv: \n";
      cerr<<"m = "<<Type(m)<<"  "<<m<<endl;
      cerr<<"x = "<<Type(x)<<endl;
      cerr<<"-> x = "<<x<<endl;
      cerr<<"x2 = "<<x2<<endl;
      cerr<<"QR = "<<QR2<<endl;
      cerr<<"QR x = "<<QR2*x<<endl;
      cerr<<"QR x2 = "<<QR2*x2<<endl;
      abort();
    }
#endif
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

#ifdef XDEBUG
    Matrix<T2> m0 = m;
    Matrix<T1> QR2 = GetQFromBandQR(QRx,Qbeta)*QRx.Diags(0,QRx.nhi()+1);
    QR2.DivideUsing(tmv::QR);
    Matrix<T2> m2 = m / QR2;
#endif
    const size_t N = QRx.rowsize();

    BandQ_LDivEq(QRx,Qbeta,m);
    BandTriLDivEq(QRx.SubBandMatrix(0,N,0,N,0,QRx.nhi()),m,NonUnitDiag);

#ifdef XDEBUG
    if (Norm(m2-m) > 0.001*Norm(m)) {
      cerr<<"BandQR_LDivEq: \n";
      cerr<<"m = "<<Type(m)<<"  "<<m0<<endl;
      cerr<<"-> m = "<<m<<endl;
      cerr<<"m2 = "<<m2<<endl;
      cerr<<"QR = "<<QR2<<endl;
      cerr<<"QR m = "<<QR2*m<<endl;
      cerr<<"QR m2 = "<<QR2*m2<<endl;
      abort();
    }
#endif
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

#ifdef XDEBUG
    Matrix<T1> QR2 = GetQFromBandQR(QRx,Qbeta)*QRx.Diags(0,QRx.nhi()+1);
    QR2.DivideUsing(tmv::QR);
    Matrix<T3> x2 = m % QR2;
#endif
    const size_t M = QRx.colsize();
    const size_t N = QRx.rowsize();

    x.Cols(N,M).Zero();
    x.Cols(0,N) = m;
    BandTriLDivEq(QRx.SubBandMatrix(0,N,0,N,0,QRx.nhi()).Transpose(),
	x.Cols(0,N).Transpose(),NonUnitDiag);

    BandQ_RDivEq(QRx,Qbeta,x);

#ifdef XDEBUG
    if (Norm(x2-x) > 0.001*Norm(x)) {
      cerr<<"BandQR_RDiv: \n";
      cerr<<"m = "<<Type(m)<<"  "<<m<<endl;
      cerr<<"x = "<<Type(x)<<endl;
      cerr<<"-> x = "<<x<<endl;
      cerr<<"x2 = "<<x2<<endl;
      cerr<<"QR = "<<QR2<<endl;
      cerr<<"x QR = "<<x*QR2<<endl;
      cerr<<"x2 QR = "<<x2*QR2<<endl;
      abort();
    }
#endif
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

#ifdef XDEBUG
    Matrix<T2> m0 = m;
    Matrix<T1> QR2 = GetQFromBandQR(QRx,Qbeta)*QRx.Diags(0,QRx.nhi()+1);
    QR2.DivideUsing(tmv::QR);
    Matrix<T2> m2 = m % QR2;
#endif

    const size_t N = QRx.rowsize();

    // Solve x Q R = m in place (m <- x)
    BandTriLDivEq(QRx.SubBandMatrix(0,N,0,N,0,QRx.nhi()).Transpose(),
	m.Transpose(),NonUnitDiag);
    BandQ_RDivEq(QRx,Qbeta,m);

#ifdef XDEBUG
    if (Norm(m2-m) > 0.001*Norm(m)) {
      cerr<<"BandQR_RDivEq: \n";
      cerr<<"m = "<<Type(m)<<"  "<<m0<<endl;
      cerr<<"-> m = "<<m<<endl;
      cerr<<"m2 = "<<m2<<endl;
      cerr<<"QR = "<<QR2<<endl;
      cerr<<"m QR = "<<m*QR2<<endl;
      cerr<<"m2 QR = "<<m2*QR2<<endl;
      abort();
    }
#endif
  }

#define InstFile "TMV_BandQRDiv_B.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


