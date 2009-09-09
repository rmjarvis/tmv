///////////////////////////////////////////////////////////////////////////////
// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2009                                                 //
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



#include "tmv/TMV_BandQRD.h"
#include "TMV_BandQRDiv.h"
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_DiagMatrix.h"
#include "TMV_Householder.h"
#include "TMV_BandLUDiv.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_TriMatrixArith.h"
#include "tmv/TMV_BandMatrixArith.h"
#include "tmv/TMV_MatrixArith.h"
#include <ostream>

namespace tmv {

  template <class T> struct BandQRDiv<T>::BandQRDiv_Impl
  {
  public :
    BandQRDiv_Impl(const GenBandMatrix<T>& A, bool _inplace);

    const bool istrans;
    const bool inplace;
    auto_array<T> Aptr1;
    T* Aptr;
    BandMatrixView<T> QRx;
    Vector<T> Qbeta;
    mutable RealType(T) logdet;
    mutable T signdet;
    mutable bool donedet;
  };

#define NEWLO (istrans ? A.nhi() : A.nlo())
#define NEWHI MIN(A.nlo()+A.nhi(),int(istrans?A.colsize():A.rowsize())-1)
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

  template <class T> BandQRDiv<T>::BandQRDiv_Impl::BandQRDiv_Impl(
      const GenBandMatrix<T>& A, bool _inplace) :
    istrans(A.colsize()<A.rowsize() || (A.IsSquare() && A.nhi()<A.nlo())
        || (A.IsSquare() && A.nhi()==A.nlo() && A.isrm())),
    inplace(_inplace || NEWLO == 0),
    Aptr1(APTR1), Aptr(APTR), QRx(QRX), Qbeta(QRx.rowsize()),
    logdet(0), signdet(1), donedet(false) {}

#undef NEWLO
#undef NEWHI
#undef APTR
#undef QRX

  template <class T> BandQRDiv<T>::BandQRDiv(const GenBandMatrix<T>& A,
      bool inplace) : pimpl(new BandQRDiv_Impl(A,inplace))
  {
    if (inplace) {
      // For inplace decomposition, make sure the original band matrix
      // has room for the extra upper diagonals...
      // if isrm stepi >= (2*A.nlo()+A.nhi())
      // if iscm stepj >= (2*A.nlo()+A.nhi())
      // if isdm extra diags appear at end, so can't really check
      TMVAssert(!pimpl->QRx.isrm() || 
          pimpl->QRx.stepi()>=pimpl->QRx.nlo()+pimpl->QRx.nhi());
      TMVAssert(!pimpl->QRx.iscm() || 
          pimpl->QRx.stepj()>=pimpl->QRx.nlo()+pimpl->QRx.nhi());
      TMVAssert(pimpl->QRx == (pimpl->istrans ? A.Transpose() : A.View()));
      if (pimpl->QRx.nlo() > 0)
        pimpl->QRx.Diags(pimpl->QRx.nhi()-pimpl->QRx.nlo()+1,
            pimpl->QRx.nhi()+1).Zero();
    } else {
      if (pimpl->istrans) pimpl->QRx = A.Transpose();
      else pimpl->QRx = A;
    }
    TMVAssert(pimpl->QRx.colsize() == (pimpl->istrans?A.rowsize():A.colsize()));
    TMVAssert(pimpl->QRx.rowsize() == (pimpl->istrans?A.colsize():A.rowsize()));
    if (pimpl->QRx.nlo() > 0) 
      QR_Decompose(pimpl->QRx,pimpl->Qbeta.View(),pimpl->signdet);
  }

  template <class T> BandQRDiv<T>::~BandQRDiv() {}

  template <class T> template <class T1> void BandQRDiv<T>::DoLDivEq(
      const MatrixView<T1>& m) const
  {
    if (pimpl->istrans)
      QR_RDivEq(pimpl->QRx,pimpl->Qbeta,m.Transpose());
    else 
      QR_LDivEq(pimpl->QRx,pimpl->Qbeta,m);
  }

  template <class T> template <class T1> void BandQRDiv<T>::DoRDivEq(
      const MatrixView<T1>& m) const
  {
    if (pimpl->istrans) QR_LDivEq(pimpl->QRx,pimpl->Qbeta,m.Transpose());
    else QR_RDivEq(pimpl->QRx,pimpl->Qbeta,m);
  }

  template <class T> template <class T1, class T2> void BandQRDiv<T>::DoLDiv(
      const GenMatrix<T1>& m, const MatrixView<T2>& x) const
  {
    TMVAssert(m.rowsize() == x.rowsize());
    TMVAssert(m.colsize() == colsize());
    TMVAssert(x.colsize() == rowsize());
    if (pimpl->istrans) 
      QR_RDiv(pimpl->QRx,pimpl->Qbeta,m.Transpose(),x.Transpose());
    else 
      QR_LDiv(pimpl->QRx,pimpl->Qbeta,m,x);
  }

  template <class T> template <class T1, class T2> void BandQRDiv<T>::DoRDiv(
      const GenMatrix<T1>& m, const MatrixView<T2>& x) const
  {
    TMVAssert(m.colsize() == x.colsize());
    TMVAssert(m.rowsize() == rowsize());
    TMVAssert(x.rowsize() == colsize());
    if (pimpl->istrans) 
      QR_LDiv(pimpl->QRx,pimpl->Qbeta,m.Transpose(),x.Transpose());
    else 
      QR_RDiv(pimpl->QRx,pimpl->Qbeta,m,x);
  }

  template <class T> T BandQRDiv<T>::Det() const
  {
    if (!pimpl->donedet) {
      T s;
      pimpl->logdet = DiagMatrixViewOf(pimpl->QRx.diag()).LogDet(&s);
      pimpl->signdet *= s;
      pimpl->donedet = true;
    }         
    if (pimpl->signdet == T(0)) return T(0);
    else return pimpl->signdet * EXP(pimpl->logdet);  
  }                  

  template <class T> RealType(T) BandQRDiv<T>::LogDet(T* sign) const
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

  template <class T> template <class T1> void BandQRDiv<T>::DoInverse(
      const MatrixView<T1>& minv) const
  {
    MatrixView<T1> minv2 = pimpl->istrans ? minv.Transpose() : minv;
    TMVAssert(pimpl->QRx.colsize() >= pimpl->QRx.rowsize());
    TMVAssert(minv2.colsize() == pimpl->QRx.rowsize());
    TMVAssert(minv2.rowsize() == pimpl->QRx.colsize());
    QR_Inverse(pimpl->QRx,pimpl->Qbeta,minv2);
  }

  template <class T> void BandQRDiv<T>::DoInverseATA(
      const MatrixView<T>& minv) const
  {
    TMVAssert(minv.colsize() == pimpl->QRx.rowsize());
    TMVAssert(minv.rowsize() == pimpl->QRx.rowsize());
    TMVAssert(pimpl->QRx.colsize() >= pimpl->QRx.rowsize());
    const int N = pimpl->QRx.rowsize();
    // At A = Rt R
    // (At A)^-1 = (Rt R)^-1 = R^-1 Rt^-1
    UpperTriMatrixView<T> Rinv = minv.UpperTri();
    Rinv = pimpl->QRx.SubBandMatrix(0,N,0,N,0,pimpl->QRx.nhi());
    Tri_Inverse(Rinv,pimpl->QRx.nhi());
    minv = Rinv * Rinv.Transpose();
  }

  template <class T> bool BandQRDiv<T>::Singular() const
  { return Det() == T(0); }

  template <class T> bool BandQRDiv<T>::IsTrans() const
  { return pimpl->istrans; }

  template <class T> void GetQFromBandQR(
      const MatrixView<T>& Q, const GenVector<T>& Qbeta, const int nlo) 
  {
    // Extract the Q matrix from a combined QRx matrix
    TMVAssert(Q.colsize() >= Q.rowsize());
    TMVAssert(Qbeta.size() == Q.rowsize());
    const int M = Q.colsize();
    const int N = Q.rowsize();
    Q.UpperTri().Zero();
    for(int j=N-1;j>=0;j--) {
      if (j+nlo+1 > M)
        Householder_Unpack(Q.SubMatrix(j,M,j,N),Qbeta(j));
      else
        Householder_Unpack(Q.SubMatrix(j,j+nlo+1,j,N),Qbeta(j));
    }
  }

  template <class T> Matrix<T> BandQRDiv<T>::GetQ() const 
  {
    // Extract the Q matrix from a combined QRx matrix
    TMVAssert(pimpl->QRx.colsize() >= pimpl->QRx.rowsize());
    Matrix<T> temp(pimpl->QRx.colsize(),pimpl->QRx.rowsize());
    temp = pimpl->QRx;
    if (pimpl->QRx.nlo() == 0) temp.SetToIdentity();
    else GetQFromBandQR(temp.View(),pimpl->Qbeta.View(),pimpl->QRx.nlo());
    return temp;
  }

  template <class T> ConstBandMatrixView<T> BandQRDiv<T>::GetR() const
  { 
    return pimpl->QRx.SubBandMatrix(
        0,pimpl->QRx.rowsize(),0,pimpl->QRx.rowsize(),0,pimpl->QRx.nhi()); 
  }
  template <class T> const GenBandMatrix<T>& BandQRDiv<T>::GetQR() const 
  { return pimpl->QRx; }
  template <class T> const GenVector<T>& BandQRDiv<T>::GetQbeta() const 
  { return pimpl->Qbeta; }

  template <class T> bool BandQRDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, std::ostream* fout) const
  {
    Matrix<T> mm = m;
    Matrix<T> Q = GetQ();
    if (fout) {
      *fout <<"BandQRDiv:\n";
      *fout << "M = "<<(pimpl->istrans?mm.Transpose():mm.View())<<std::endl;
      *fout << "Q = "<<Q<<std::endl;
      *fout << "R = "<<GetR()<<std::endl;
    }
    Matrix<T> qr = Q*GetR();
    RealType(T) nm = Norm(qr-(pimpl->istrans ? mm.Transpose() : mm.View()));
    nm /= Norm(Q)*Norm(GetR());
    if (fout) {
      *fout << "QR = "<<qr<<std::endl;
      *fout << "Norm(M-QR)/Norm(QR) = "<<nm<<std::endl;
    }
    return nm < mm.DoCondition()*mm.colsize()*Epsilon<T>();
  }

  template <class T> size_t BandQRDiv<T>::colsize() const
  { return pimpl->istrans ? pimpl->QRx.rowsize() : pimpl->QRx.colsize(); }

  template <class T> size_t BandQRDiv<T>::rowsize() const
  { return pimpl->istrans ? pimpl->QRx.colsize() : pimpl->QRx.rowsize(); }

#define InstFile "TMV_BandQRD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv
