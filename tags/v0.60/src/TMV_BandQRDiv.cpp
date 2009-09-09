///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 2007                                                        //
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
// along with this program in the file gpl.txt.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



#include "TMV_BandMatrix.h"
#include "TMV_MatrixArith.h"
#include "TMV_BandMatrixArith.h"
#include "TMV_Householder.h"
#include "TMV_BandQRDiv.h"
#include "TMV_BandLUDiv.h"
#include "TMV_DiagMatrix.h"
#include <ostream>

namespace tmv {

  template <class T> struct BandQRDiv<T>::BandQRDiv_Impl
  {
    BandQRDiv_Impl(const GenBandMatrix<T>& A, bool _inplace);

    const bool istrans;
    const bool inplace;
    auto_array<T> Aptr1;
    T* Aptr;
    BandMatrixView<T> QRx;
    Vector<T> Qbeta;
    mutable T det;
    mutable bool donedet;
  };

#define NEWLO (istrans ? A.nhi() : A.nlo())
#define NEWHI std::min(A.nlo()+A.nhi(),int(istrans?A.colsize():A.rowsize())-1)
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
    det(T(1)), donedet(false) {}

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
    if (pimpl->QRx.nlo() > 0) 
      BandQR_Decompose(pimpl->QRx,pimpl->Qbeta.View(),pimpl->det);
  }

  template <class T> BandQRDiv<T>::~BandQRDiv() { delete pimpl; }

#define CT std::complex<T>
  template <class T> inline void BandQR_LDivEq(
      const GenBandMatrix<CT>& , const GenVector<CT>& ,
      const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void BandQR_RDivEq(
      const GenBandMatrix<CT>& , const GenVector<CT>& ,
      const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void BandQR_LDiv(
      const GenBandMatrix<CT>& , const GenVector<CT>& ,
      const GenMatrix<T>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void BandQR_LDiv(
      const GenBandMatrix<CT>& , const GenVector<CT>& ,
      const GenMatrix<CT>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void BandQR_LDiv(
      const GenBandMatrix<T>& , const GenVector<T>& ,
      const GenMatrix<CT>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void BandQR_RDiv(
      const GenBandMatrix<CT>& , const GenVector<CT>& ,
      const GenMatrix<T>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void BandQR_RDiv(
      const GenBandMatrix<CT>& , const GenVector<CT>& ,
      const GenMatrix<CT>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void BandQR_RDiv(
      const GenBandMatrix<T>& , const GenVector<T>& ,
      const GenMatrix<CT>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void BandQ_LDivEq(
      const GenBandMatrix<CT>& , const GenVector<CT>& ,
      const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void BandQ_RDivEq(
      const GenBandMatrix<CT>& , const GenVector<CT>& ,
      const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void BandTriLDivEq(
      const GenBandMatrix<std::complex<T> >& , const MatrixView<T>& , 
      DiagType )
  { TMVAssert(FALSE); }
  template <class T> inline void BandTriLDivEq(
      const GenBandMatrix<std::complex<T> >& , const VectorView<T>& , 
      DiagType )
  { TMVAssert(FALSE); }
#undef CT

  template <class T> template <class T1> void BandQRDiv<T>::DoLDivEq(
      const MatrixView<T1>& m) const
  {
    if (pimpl->istrans)
      BandQR_RDivEq(pimpl->QRx,pimpl->Qbeta,m.Transpose());
    else 
      BandQR_LDivEq(pimpl->QRx,pimpl->Qbeta,m);
  }

  template <class T> template <class T1> void BandQRDiv<T>::DoRDivEq(
      const MatrixView<T1>& m) const
  {
    if (pimpl->istrans) BandQR_LDivEq(pimpl->QRx,pimpl->Qbeta,m.Transpose());
    else BandQR_RDivEq(pimpl->QRx,pimpl->Qbeta,m);
  }

  template <class T> template <class T1, class T2> void BandQRDiv<T>::DoLDiv(
      const GenMatrix<T1>& m, const MatrixView<T2>& x) const
  {
    TMVAssert(m.rowsize() == x.rowsize());
    TMVAssert(m.colsize() == colsize());
    TMVAssert(x.colsize() == rowsize());
    if (pimpl->istrans) 
      BandQR_RDiv(pimpl->QRx,pimpl->Qbeta,m.Transpose(),x.Transpose());
    else 
      BandQR_LDiv(pimpl->QRx,pimpl->Qbeta,m,x);
  }

  template <class T> template <class T1, class T2> void BandQRDiv<T>::DoRDiv(
      const GenMatrix<T1>& m, const MatrixView<T2>& x) const
  {
    TMVAssert(m.colsize() == x.colsize());
    TMVAssert(m.rowsize() == rowsize());
    TMVAssert(x.rowsize() == colsize());
    if (pimpl->istrans) 
      BandQR_LDiv(pimpl->QRx,pimpl->Qbeta,m.Transpose(),x.Transpose());
    else 
      BandQR_RDiv(pimpl->QRx,pimpl->Qbeta,m,x);
  }

  template <class T> T BandQRDiv<T>::Det() const
  {
    if (!pimpl->donedet) {
      pimpl->det *= DiagMatrixViewOf(pimpl->QRx.diag()).Det();
      pimpl->donedet = true;
    }
    return pimpl->det;
  }

  template <class T> template <class T1> void BandQRDiv<T>::DoInverse(
      const MatrixView<T1>& minv) const
  {
    MatrixView<T1> minv2 = pimpl->istrans ? minv.Transpose() : minv;

    TMVAssert(pimpl->QRx.colsize() >= pimpl->QRx.rowsize());
    TMVAssert(minv2.colsize() == pimpl->QRx.rowsize());
    TMVAssert(minv2.rowsize() == pimpl->QRx.colsize());
    const size_t M = pimpl->QRx.colsize();
    const size_t N = pimpl->QRx.rowsize();
    minv2.Cols(N,M).Zero();
    minv2.Cols(0,N).SetToIdentity();
    BandTriLDivEq(
	pimpl->QRx.SubBandMatrix(0,N,0,N,0,pimpl->QRx.nhi()).Transpose(),
	minv2.Cols(0,N).Transpose(),NonUnitDiag);
    BandQ_RDivEq(pimpl->QRx,pimpl->Qbeta,minv2);
  }

  template <class T> void BandQRDiv<T>::DoInverseATA(
      const MatrixView<T>& minv) const
  {
    TMVAssert(pimpl->QRx.colsize() >= pimpl->QRx.rowsize());
    const size_t N = pimpl->QRx.rowsize();
    // At A = Rt R
    // (At A)^-1 = (Rt R)^-1 = R^-1 Rt^-1
    Matrix<T,ColMajor> Rinv(N,N);
    Rinv.SetToIdentity();

    BandTriLDivEq(pimpl->QRx.SubBandMatrix(0,N,0,N,0,pimpl->QRx.nhi()),
	Rinv.View(),NonUnitDiag);
    minv = Rinv * Rinv.Transpose();
  }

  template <class T> bool BandQRDiv<T>::Singular() const
  { return Det() == T(0); }

  template <class T> bool BandQRDiv<T>::IsTrans() const
  { return pimpl->istrans; }

  template <class T> Matrix<T> GetQFromBandQR(
      const GenBandMatrix<T>& QRx, const GenVector<T>& Qbeta) 
  {
    // Extract the Q matrix from a combined QRx matrix
    TMVAssert(QRx.colsize() >= QRx.rowsize());
    const size_t N = QRx.rowsize();
    const size_t M = QRx.colsize();
    Matrix<T> Q(BandMatrixViewOf(QRx,QRx.nlo(),0));
    if (QRx.nlo() == 0) Q.SetToIdentity();
    else {
      for(int j=N-1;j>=0;j--) {
	Householder_Unpack(Q.SubMatrix(j,M,j,N),Qbeta(j));
      }
    }
    return Q;
  }

  template <class T> Matrix<T> BandQRDiv<T>::GetQ() const 
  {
    // Extract the Q matrix from a combined QRx matrix
    TMVAssert(pimpl->QRx.colsize() >= pimpl->QRx.rowsize());
    return GetQFromBandQR(pimpl->QRx,pimpl->Qbeta);
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

  template <class T> std::string BandQRDiv<T>::Type() const
  { return std::string("BandQRDiv<") + tmv::Type(T()) + ">"; }
  template <class T> DivType BandQRDiv<T>::GetDivType() const 
  { return QR; }

  template <class T> bool BandQRDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, std::ostream* fout) const
  {
    Matrix<T> mm = m;
    if (fout) {
      *fout <<"BandQRDiv:\n";
      *fout << "M = "<<(pimpl->istrans?mm.Transpose():mm.View())<<std::endl;
      *fout << "Q = "<<GetQ()<<std::endl;
      *fout << "R = "<<GetR()<<std::endl;
    }
    Matrix<T> qr = GetQ()*GetR();
    RealType(T) nm = Norm(qr-(pimpl->istrans ? mm.Transpose() : mm.View()));
    nm /= Norm(GetQ())*Norm(GetR());
    if (fout) {
      *fout << "QR = "<<qr<<std::endl;
      *fout << "Norm(M-QR)/Norm(QR) = "<<nm<<std::endl;
    }
    mm.DivideUsing(SVS);
    mm.SetDiv();
    return nm < mm.Condition()*mm.colsize()*Epsilon<T>();
  }

  template <class T> size_t BandQRDiv<T>::colsize() const
  { return pimpl->istrans ? pimpl->QRx.rowsize() : pimpl->QRx.colsize(); }

  template <class T> size_t BandQRDiv<T>::rowsize() const
  { return pimpl->istrans ? pimpl->QRx.colsize() : pimpl->QRx.rowsize(); }

#define InstFile "TMV_BandQRDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv
