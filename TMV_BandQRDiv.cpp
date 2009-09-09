
#include "TMV.h"
#include "TMV_Band.h"
#include "TMV_Householder.h"
#include "TMV_Diag.h"

namespace tmv {

#define NEWLO (istrans ? A.nhi() : A.nlo())
#define NEWHI min(A.nlo()+A.nhi(),int(A.colsize())-1)
#define APTR1 inplace ? 0 : \
  new T[BandStorageLength(ColMajor, \
      istrans ? A.rowsize() : A.colsize(), \
      istrans ? A.colsize() : A.rowsize(), NEWLO, NEWHI)]
#define APTR inplace ? A.NonConst().ptr() : Aptr1.get()
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
    Aptr1(APTR1), Aptr(APTR), QRx(QRX), Qbeta(QRx.rowsize()),
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
    if (QRx.nlo() > 0) BandQR_Decompose(QRx,Qbeta.View(),det);
  }

#undef NEWLO
#undef NEWHI
#undef APTR
#undef QRX

} // namespace tmv

namespace tmv {
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

  template <class T> Matrix<T> GetQFromBandQR(
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

  template <class T> T BandQRDiv<T>::Det() const
  {
    if (!donedet) {
      det *= DiagMatrixViewOf(QRx.diag()).Det();
      donedet = true;
    }
    return det;
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
