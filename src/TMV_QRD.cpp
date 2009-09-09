///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2008                                                        //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file LICENSE.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



#include "TMV_QRD.h"
#include "TMV_QRDiv.h"
#include "TMV_Matrix.h"
#include "TMV_TriMatrix.h"
#include "TMV_DiagMatrix.h"
#include "TMV_TriMatrixArith.h"
#include "TMV_MatrixArith.h"
#include <ostream>

namespace tmv {

  template <class T> struct QRDiv<T>::QRDiv_Impl
  {
    public :
      QRDiv_Impl(const GenMatrix<T>& A, bool inplace);

      const bool istrans;
      const bool inplace;
      auto_array<T> Aptr1;
      T* Aptr;
      MatrixView<T> QRx;
      Vector<T> beta;
      mutable RealType(T) logdet;
      mutable T signdet;
      mutable bool donedet;

    private :
      QRDiv_Impl(const QRDiv_Impl&) :
	istrans(false), inplace(false), Aptr1(0), Aptr(0), 
	QRx(MatrixViewOf(Aptr,0,0,ColMajor)), 
	beta(0), logdet(0), signdet(0), donedet(false)
      { TMVAssert(FALSE); }
      QRDiv_Impl& operator=(const QRDiv_Impl&)
      { TMVAssert(FALSE); return *this; }
  };

#define APTR1 (inplace ? 0 : new T[A.colsize()*A.rowsize()])
#define APTR (inplace ? A.NonConst().ptr() : Aptr1.get())
#define QRX (istrans ? \
    (inplace ? A.NonConst().Transpose() : \
     MatrixViewOf(Aptr,A.rowsize(),A.colsize(),ColMajor)) : \
    (inplace ? A.NonConst().View() : \
     MatrixViewOf(Aptr,A.colsize(),A.rowsize(),ColMajor)))

  template <class T> QRDiv<T>::QRDiv_Impl::QRDiv_Impl(
      const GenMatrix<T>& A, bool _inplace) :
    istrans(A.colsize()<A.rowsize()),
    inplace(_inplace && (A.iscm() || A.isrm())), 
    Aptr1(APTR1), Aptr(APTR), QRx(QRX), 
    beta(QRx.rowsize()), logdet(0), signdet(1), donedet(false) {}

#undef QRX
#undef APTR

  template <class T> QRDiv<T>::QRDiv(const GenMatrix<T>& A,
      bool inplace) :
    pimpl(new QRDiv_Impl(A,inplace))
  {
    if (pimpl->istrans) {
      if (inplace) TMVAssert(A.Transpose() == pimpl->QRx);
      else pimpl->QRx = A.Transpose();
    }
    else {
      if (inplace) TMVAssert(A == pimpl->QRx); 
      else pimpl->QRx = A;
    }
    QR_Decompose(pimpl->QRx,pimpl->beta.View(),pimpl->signdet);
  }

  template <class T> QRDiv<T>::~QRDiv() { delete pimpl; pimpl=0; }

  template <class T> template <class T1> void QRDiv<T>::DoLDivEq(
      const MatrixView<T1>& m) const
  {
    if (pimpl->istrans) 
      QR_LDivEq(pimpl->QRx,pimpl->beta,0,m.Transpose(),pimpl->QRx.rowsize());
    else 
      QR_LDivEq(pimpl->QRx,pimpl->beta,0,m,pimpl->QRx.rowsize());
  }

  template <class T> template <class T1> void QRDiv<T>::DoRDivEq(
      const MatrixView<T1>& m) const
  {
    if (pimpl->istrans) 
      QR_RDivEq(pimpl->QRx,pimpl->beta,0,m.Transpose(),pimpl->QRx.rowsize());
    else 
      QR_RDivEq(pimpl->QRx,pimpl->beta,0,m,pimpl->QRx.rowsize());
  }

  template <class T> template <class T1, class T2> void QRDiv<T>::DoLDiv(
      const GenMatrix<T1>& m, const MatrixView<T2>& x) const
  {
    TMVAssert(m.colsize() == (pimpl->istrans ? pimpl->QRx.rowsize() :
	pimpl->QRx.colsize()));
    TMVAssert(x.colsize() == (pimpl->istrans ? pimpl->QRx.colsize() :
	pimpl->QRx.rowsize()));
    if (pimpl->istrans) 
      QR_RDiv(pimpl->QRx,pimpl->beta,0,m.Transpose(),x.Transpose(),
	  pimpl->QRx.rowsize());
    else QR_LDiv(pimpl->QRx,pimpl->beta,0,m,x,pimpl->QRx.rowsize());
  }

  template <class T> template <class T1, class T2> void QRDiv<T>::DoRDiv(
      const GenMatrix<T1>& m, const MatrixView<T2>& x) const
  {
    TMVAssert(m.colsize() == x.colsize());
    TMVAssert(m.rowsize() == (pimpl->istrans ? pimpl->QRx.colsize() :
	pimpl->QRx.rowsize()));
    TMVAssert(x.rowsize() == (pimpl->istrans ? pimpl->QRx.rowsize() :
	pimpl->QRx.colsize()));
    if (pimpl->istrans) 
      QR_LDiv(pimpl->QRx,pimpl->beta,0,m.Transpose(),x.Transpose(),
	  pimpl->QRx.rowsize());
    else 
      QR_RDiv(pimpl->QRx,pimpl->beta,0,m,x,pimpl->QRx.rowsize());
  }

  template <class T> T QRDiv<T>::Det() const
  {
    if (!pimpl->donedet) {
      T s;
      pimpl->logdet = DiagMatrixViewOf(pimpl->QRx.diag()).LogDet(&s);
      pimpl->signdet *= s;
      pimpl->donedet = true;
    }         
    if (pimpl->signdet == T(0)) return T(0);
    else return pimpl->signdet * std::exp(pimpl->logdet);  
  }                  

  template <class T> RealType(T) QRDiv<T>::LogDet(T* sign) const
  {
    if (!pimpl->donedet) {
      T s;
      pimpl->logdet = DiagMatrixViewOf(pimpl->QRx.diag()).LogDet(&s);
      pimpl->signdet *= s;
      pimpl->donedet = true;
    }
    if (sign) *sign = pimpl->signdet;
    return pimpl->logdet;  
  }                  
  template <class T> template <class T1> void QRDiv<T>::DoInverse(
      const MatrixView<T1>& minv) const
  {
    if (pimpl->istrans)
      QR_Inverse(pimpl->QRx,pimpl->beta,minv.Transpose()); 
    else
      QR_Inverse(pimpl->QRx,pimpl->beta,minv); 
  }

  template <class T> void QRDiv<T>::DoInverseATA(
      const MatrixView<T>& ata) const
  {
    // At A = Rt Qt Q R = Rt R
    // (At A)^-1 = (Rt R)^-1 = R^-1 * Rt^-1

    UpperTriMatrixView<T> rinv = UpperTriMatrixViewOf(ata);
    rinv = UpperTriMatrixViewOf(pimpl->QRx).Inverse();
    ata = rinv * rinv.Adjoint();
  }

  template <class T> bool QRDiv<T>::Singular() const 
  { return Det() == T(0); }

  template <class T> bool QRDiv<T>::IsTrans() const 
  { return pimpl->istrans; }

  template <class T> Matrix<T> QRDiv<T>::GetQ() const
  {
    Matrix<T> Q = pimpl->QRx;
    GetQFromQR(Q.View(),pimpl->beta);
    return Q;
  }

  template <class T> ConstUpperTriMatrixView<T> QRDiv<T>::GetR() const
  { return UpperTriMatrixViewOf(pimpl->QRx,NonUnitDiag); }

  template <class T> const GenMatrix<T>& QRDiv<T>::GetQR() const 
  { return pimpl->QRx; }

  template <class T> const GenVector<T>& QRDiv<T>::Getbeta() const 
  { return pimpl->beta; }

  template <class T> bool QRDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, std::ostream* fout) const
  {
    Matrix<T> mm = m;
    bool printmat = fout && m.colsize() < 100 && m.rowsize() < 100;
    if (printmat) {
      *fout << "QRDiv:\n";
      *fout << "M = "<<(pimpl->istrans?mm.Transpose():mm.View())<<std::endl;
      *fout << "Q = "<<GetQ()<<std::endl;
      *fout << "R = "<<GetR()<<std::endl;
    }
    Matrix<T> qr = GetQ()*GetR();
    RealType(T) nm = Norm(qr-(pimpl->istrans ? mm.Transpose() : mm.View()));
    nm /= Norm(GetQ())*Norm(GetR());
    if (printmat) {
      *fout << "QR = "<<qr<<std::endl;
    }
    if (fout) {
      *fout << "Norm(M-QR)/Norm(QR) = "<<nm<<"  ";
      *fout << pimpl->QRx.rowsize()*Epsilon<T>()<<std::endl;
    }
    return nm < mm.DoCondition()*mm.colsize()*Epsilon<T>();
  }

  template <class T> size_t QRDiv<T>::colsize() const
  { return pimpl->istrans ? pimpl->QRx.rowsize() : pimpl->QRx.colsize(); }

  template <class T> size_t QRDiv<T>::rowsize() const
  { return pimpl->istrans ? pimpl->QRx.colsize() : pimpl->QRx.rowsize(); }

#define InstFile "TMV_QRD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


