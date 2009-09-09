
#include "TMV.h"
#include "TMV_Diag.h"
#include "TMV_Tri.h"

//#define XDEBUG

namespace tmv {

#define APTR1 inplace ? 0 : new T[A.colsize()*A.rowsize()]
#define APTR inplace ? A.NonConst().ptr() : Aptr1.get()
#define QRX istrans ? \
  inplace ? A.NonConst().Transpose() : \
  MatrixViewOf(Aptr,A.rowsize(),A.colsize(),ColMajor) : \
  inplace ? A.NonConst().View() : \
  MatrixViewOf(Aptr,A.colsize(),A.rowsize(),ColMajor)

  template <class T> QRDiv<T>::QRDiv(const GenMatrix<T>& A,
      bool _inplace) :
    istrans(A.colsize()<A.rowsize()),
    inplace(_inplace && (A.iscm() || A.isrm())), 
    Aptr1(APTR1), Aptr(APTR), QRx(QRX), 
    beta(QRx.rowsize()), det(T(1)), donedet(false)
  {
    if (istrans) {
      if (inplace) TMVAssert(A.Transpose() == QRx);
      else QRx = A.Transpose();
    }
    else {
      if (inplace) TMVAssert(A == QRx); 
      else QRx = A;
    }
    QR_Decompose(QRx,beta.View(),det);
  }

#undef QRX
#undef APTR

  template <class T> bool QRDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenMatrix<T>* mm = dynamic_cast<const GenMatrix<T>*>(&m);
    TMVAssert(mm);
#ifdef XDEBUG
    bool printmat = fout && m.colsize() < 100 && m.rowsize() < 100;
    if (printmat) {
      *fout << "M = "<<tmv::Type(*mm)<<"  ";
      *fout << (istrans ? mm->Transpose() : mm->View()) <<endl;
      *fout << "Q = "<<GetQ()<<endl;
      *fout << "R = "<<GetR()<<endl;
    }
#endif
    Matrix<T> qr = GetQ()*GetR();
    RealType(T) nm = Norm(qr- (istrans ? mm->Transpose() : mm->View()) );
    nm /= Norm(GetQ())*Norm(GetR());
#ifdef XDEBUG
    if (printmat) {
      *fout << "QR = "<<qr<<endl;
    }
#endif
    if (fout) {
      *fout << "Norm(M-QR)/Norm(QR) = "<<nm<<"  "<<QRx.rowsize()*Epsilon<T>()<<endl;
    }
    Matrix<T> m2 = *mm;
    m2.DivideUsing(SVS);
    m2.SetDiv();
    return nm < m2.Condition()*m2.colsize()*Epsilon<T>();
  }

  template <class T> T QRDiv<T>::Det() const
  {
    if (!donedet) {
      det *= DiagMatrixViewOf(QRx.diag()).Det();
      donedet = true;
    }
    return det;
  }

  template <class T> Matrix<T> QRDiv<T>::GetQ() const
  {
    Matrix<T> Q = QRx;
    GetQFromQR(Q.View(),beta);
    return Q;
  }

#define InstFile "TMV_QRDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


