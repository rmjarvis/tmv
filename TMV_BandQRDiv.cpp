
#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_Householder.h"
#include "TMV_Diag.h"

//#define XDEBUG

namespace tmv {

#define NEWLO (istrans ? A.nhi() : A.nlo())
#define NEWHI min(A.nlo()+A.nhi(),int(A.colsize())-1)
#define APTR (inplace ? A.NonConst().ptr() : \
  new T[BandStorageLength(ColMajor, \
      istrans ? A.rowsize() : A.colsize(), \
      istrans ? A.colsize() : A.rowsize(), NEWLO, NEWHI)])
#define QRX (istrans ? \
  (inplace ? A.NonConst().Transpose() : \
  BandMatrixViewOf(Aptr,A.rowsize(),A.colsize(),A.nhi(), NEWHI, ColMajor)) : \
  (inplace ? A.NonConst().View() : \
  BandMatrixViewOf(Aptr,A.colsize(),A.rowsize(),A.nlo(), NEWHI, ColMajor)))

  template <class T> BandQRDiv<T>::BandQRDiv(const GenBandMatrix<T>& A,
      bool _inplace) :
    istrans(A.colsize()<A.rowsize() || (A.IsSquare() && A.nhi()<A.nlo())
	|| (A.IsSquare() && A.nhi()==A.nlo() && A.isrm())),
    inplace(_inplace || NEWLO == 0),
    Aptr(APTR), QRx(QRX), Qbeta(QRx.rowsize()),
    det(T(1)), donedet(false)
  {
    if (inplace) {
      // For inplace decomposition, make sure the original band matrix
      // has room for the extra upper diagonals...
      // if isrm stepi >= (2*A.nlo()+A.nhi())
      // if iscm stepj >= (2*A.nlo()+A.nhi())
      // if isdm extra diags appear at end, so can't really check
      TMVAssert(!QRx.isrm() || QRx.stepi()>=QRx.nlo()+QRx.nhi());
      TMVAssert(!QRx.iscm() || QRx.stepj()>=QRx.nlo()+QRx.nhi());
      TMVAssert(QRx == (istrans ? A.Transpose() : A.View()));
      if (QRx.nlo() > 0)
	QRx.Diags(QRx.nhi()-QRx.nlo()+1,QRx.nhi()+1).Zero();
    } else {
      if (istrans) QRx = A.Transpose();
      else QRx = A;
    }
    if (QRx.nlo() > 0) 
      BandQR_Decompose(QRx,Qbeta.View(),det);
  }

#undef NEWLO
#undef NEWHI
#undef APTR
#undef QRX

  template <class T> BandQRDiv<T>::~BandQRDiv() 
  { if (!inplace) delete[] Aptr; }

  template <class T> bool BandQRDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenBandMatrix<T>* bm = dynamic_cast<const GenBandMatrix<T>*>(&m);
    TMVAssert(bm);
    if (fout) {
      *fout << "M = "<<tmv::Type(*bm)<<"  ";
      *fout << (istrans ? bm->Transpose() : bm->View()) <<endl;
      *fout << "Q = "<<GetQ()<<endl;
      *fout << "R = "<<GetR()<<endl;
    }
    Matrix<T> qr = GetQ()*GetR();
    RealType(T) nm = Norm(qr- (istrans ? bm->Transpose() : bm->View()) );
    nm /= Norm(GetQ())*Norm(GetR());
    if (fout) {
      *fout << "QR = "<<qr<<endl;
      *fout << "Norm(M-QR)/Norm(QR) = "<<nm<<endl;
    }
    BandMatrix<T> bm2 = *bm;
    bm2.DivideUsing(SVS);
    bm2.SetDiv();
    return nm < bm2.Condition()*bm2.colsize()*Epsilon<T>();
  }

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
    BandMatrix<T> A0 = QRx.Diags(-QRx.nlo(),QRx.nhi()-QRx.nlo()+1);
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
    BandMatrix<T> R = QRx.Diags(0,QRx.nhi()+1);
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

  template <class T> inline Matrix<T> GetQFromBandQR(
      const GenBandMatrix<T>& QRx, const GenVector<T>& Qbeta) 
  {
    // Extract the Q matrix from a combined QRx matrix
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    const size_t N = QRx.rowsize();
    const size_t M = QRx.colsize();
    Matrix<T> Q = BandMatrixViewOf(QRx,QRx.nlo(),0);
    if (QRx.nlo() == 0) Q.SetToIdentity();
    else {
      for(int j=N-1;j>=0;j--) {
	Householder_Unpack(Q.SubMatrix(j,M,j,N),Qbeta(j));
      }
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
    Matrix<T1> QR = GetQFromBandQR(QRx,Qbeta)*QRx.Diags(0,QRx.nhi()+1);
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
    Matrix<T3> mm = QR * x;
    if (Norm(mm-m) > 0.001*Norm(m)) {
      cerr<<"BandQR_LDiv: \n";
      cerr<<"m = "<<Type(m)<<"  "<<m<<endl;
      cerr<<"x = "<<Type(x)<<endl;
      cerr<<"-> x = "<<x<<endl;
      cerr<<"QR = "<<QR<<endl;
      cerr<<"x QR = "<<mm<<endl;
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
    Matrix<T1> QR = GetQFromBandQR(QRx,Qbeta)*QRx.Diags(0,QRx.nhi()+1);
    Matrix<T2> m0 = m;
#endif
    const size_t N = QRx.rowsize();

    BandQ_LDivEq(QRx,Qbeta,m);
    BandTriLDivEq(QRx.SubBandMatrix(0,N,0,N,0,QRx.nhi()),m,NonUnitDiag);

#ifdef XDEBUG
    Matrix<T2> mm = QR*m;
    if (Norm(mm-m0) > 0.001*Norm(m0)) {
      cerr<<"BandQR_LDivEq: \n";
      cerr<<"m = "<<Type(m)<<"  "<<m0<<endl;
      cerr<<"-> m = "<<m<<endl;
      cerr<<"QR = "<<QR<<endl;
      cerr<<"QR m = "<<mm<<endl;
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
    Matrix<T1> Q = GetQFromBandQR(QRx,Qbeta);
    Matrix<T1> R = QRx.Diags(0,QRx.nhi()+1);
    Matrix<T1> QR = Q*R;
#endif
    const size_t M = QRx.colsize();
    const size_t N = QRx.rowsize();

    x.Cols(N,M).Zero();
    x.Cols(0,N) = m;
    BandTriLDivEq(QRx.SubBandMatrix(0,N,0,N,0,QRx.nhi()).Transpose(),
	x.Cols(0,N).Transpose(),NonUnitDiag);

    BandQ_RDivEq(QRx,Qbeta,x);

#ifdef XDEBUG
    Matrix<T3> mm = x * QR;
    if (Norm(mm-m) > 0.001*Norm(m)) {
      cerr<<"BandQR_RDiv: \n";
      cerr<<"QRx = "<<QRx<<endl;
      cerr<<"Q = "<<Q<<endl;
      cerr<<"R = "<<R<<endl;
      cerr<<"m = "<<Type(m)<<"  "<<m<<endl;
      cerr<<"x = "<<Type(x)<<endl;
      cerr<<"-> x = "<<x<<endl;
      cerr<<"QR = "<<QR<<endl;
      cerr<<"x QR = "<<mm<<endl;
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
    Matrix<T1> QR = GetQFromBandQR(QRx,Qbeta)*QRx.Diags(0,QRx.nhi()+1);
    Matrix<T2> m0 = m;
#endif

    const size_t N = QRx.rowsize();

    // Solve x Q R = m in place (m <- x)
    BandTriLDivEq(QRx.SubBandMatrix(0,N,0,N,0,QRx.nhi()).Transpose(),
	m.Transpose(),NonUnitDiag);
    BandQ_RDivEq(QRx,Qbeta,m);

#ifdef XDEBUG
    Matrix<T2> mm = m * QR;
    if (Norm(mm-m0) > 0.001*Norm(m0)) {
      cerr<<"BandQR_RDivEq: \n";
      cerr<<"m = "<<Type(m)<<"  "<<m0<<endl;
      cerr<<"-> m = "<<m<<endl;
      cerr<<"QR = "<<QR<<endl;
      cerr<<"m QR = "<<mm<<endl;
      abort();
    }
#endif
  }

  template <class T> T BandQRDiv<T>::Det() const
  {
    if (!donedet) {
      det *= DiagMatrixViewOf(QRx.diag()).Det();
      donedet = true;
    }
    return det;
  }

  template <class T> void BandQRDiv<T>::Inverse(const MatrixView<T>& minv) const
  {
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    if (istrans) {
      TMVAssert(minv.colsize() == QRx.colsize());
      TMVAssert(minv.rowsize() == QRx.rowsize());
      LDiv(Eye<RealType(T),ColMajor>(QRx.rowsize()),minv);
    } else {
      TMVAssert(minv.colsize() == QRx.rowsize());
      TMVAssert(minv.rowsize() == QRx.colsize());
      RDiv(Eye<RealType(T),ColMajor>(QRx.rowsize()),minv);
    }
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

  template <class T> Matrix<T> BandQRDiv<T>::GetQ() const 
  {
    // Extract the Q matrix from a combined QRx matrix
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    return GetQFromBandQR(QRx,Qbeta);
  }

#define InstFile "TMV_BandQRDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


