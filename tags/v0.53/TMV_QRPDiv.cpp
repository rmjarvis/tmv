
#include "TMV.h"
#include "TMV_Householder.h"
#include "TMV_Tri.h"
#include "TMV_Diag.h"

//#define XDEBUG

namespace tmv {

#ifdef TMV_BLOCKSIZE
  const size_t QRP_BLOCKSIZE = TMV_BLOCKSIZE;
#else
  const size_t QRP_BLOCKSIZE = 32;
#endif

#define APTR inplace ? A.NonConst().ptr() : new T[A.colsize()*A.rowsize()]
#define QRX istrans ? \
  inplace ? A.NonConst().Transpose() : \
  MatrixViewOf(Aptr,A.rowsize(),A.colsize(),ColMajor) : \
  inplace ? A.NonConst().View() : \
  MatrixViewOf(Aptr,A.colsize(),A.rowsize(),ColMajor)

  template <class T> QRPDiv<T>::QRPDiv(const GenMatrix<T>& A, bool _inplace) :
    istrans(A.colsize()<A.rowsize()),
    inplace(_inplace && (A.iscm() || A.isrm())), Aptr(APTR), QRx(QRX), 
    beta(QRx.rowsize()), P(new size_t[beta.size()]),
    det(T(1)), donedet(false), N1(beta.size())
  {
    if (istrans) {
      if (inplace) { TMVAssert(A.Transpose() == QRx); }
      else QRx = A.Transpose();
    }
    else {
      if (inplace) { TMVAssert(A == QRx); }
      else QRx = A;
    }
    QRP_Decompose(QRx,beta.View(),P,det);
    while(N1>0 && QRx.diag()(N1-1)==T(0)) --N1;
  }
#undef QRX
#undef APTR

  template <class T> QRPDiv<T>::~QRPDiv<T>()
  { delete[] P; if (!inplace) delete[] Aptr; }

  template <class T> bool QRPDiv<T>::CheckDecomp(
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
      *fout << "P = ";
      for(size_t i=0;i<QRx.rowsize();i++) *fout<<P[i]<<" ";
      *fout<<endl;
    }
#endif
    Matrix<T> qr = GetQ()*GetR();
    qr.ReversePermuteCols(GetP());
    RealType(T) nm = Norm(qr- (istrans ? mm->Transpose() : mm->View()) );
    nm /= Norm(GetQ())*Norm(GetR());
#ifdef XDEBUG
    if (printmat) {
      *fout << "QRP = "<<qr<<endl;
    }
#endif
    if (fout) {
      *fout << "Norm(M-QR)/Norm(QR) = "<<nm<<"  "<<QRx.rowsize()*Epsilon<T>()<<endl;
    }
    Matrix<T> m2 = *mm;
    m2.DivideUsing(SVS);
    m2.SetDiv();
    return nm < m2.Condition()*Epsilon<T>();
  }

  template <class T> T QRPDiv<T>::Det() const
  {
    if (!donedet) {
      det *= DiagMatrixViewOf(QRx.diag()).Det();
      donedet = true;
    }
    return det;
  }

  template <class T> Matrix<T> QRPDiv<T>::GetQ() const
  {
    Matrix<T> Q = QRx;
    GetQFromQR(Q.View(),beta);
    return Q;
  }

#define InstFile "TMV_QRPDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


