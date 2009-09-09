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



#include "TMV_BandLUD.h"
#include "TMV_BandLUDiv.h"
#include "TMV_BandMatrix.h"
#include "TMV_TriMatrix.h"
#include "TMV_DiagMatrix.h"
#include "TMV_MatrixArith.h"
#include "TMV_TriMatrixArith.h"
#include <ostream>

namespace tmv {

  template <class T> struct BandLUDiv<T>::BandLUDiv_Impl
  {
    public :
      BandLUDiv_Impl(const GenBandMatrix<T>& A, bool _inplace);
      BandLUDiv_Impl(const AssignableToBandMatrix<T>& A);

      const bool istrans;
      const bool inplace;
      auto_array<T> Aptr1;
      T* Aptr;
      BandMatrixView<T> LUx;
      auto_array<int> P;
      mutable RealType(T) logdet;
      mutable T signdet;
      mutable bool donedet;

    private :
      BandLUDiv_Impl(const BandLUDiv_Impl&) : 
	istrans(false), inplace(false), Aptr1(0), Aptr(0),
	LUx(BandMatrixViewOf(Aptr,0,0,0,0,ColMajor)),
	P(0), logdet(0), signdet(0), donedet(false)
      { TMVAssert(FALSE); }
      BandLUDiv_Impl& operator=(const BandLUDiv_Impl&)
      { TMVAssert(FALSE); return *this; }
  };

#define NEWLO MIN(A.nlo(),A.nhi())
#define NEWHI MIN(A.nlo()+A.nhi(),int(A.colsize())-1)
#define APTR1 (inplace ? 0 : \
  new T[BandStorageLength(ColMajor,A.colsize(),A.colsize(),NEWLO,NEWHI)])
#define APTR (inplace ? A.NonConst().ptr() : Aptr1.get())

#define LUX (istrans ? \
    (inplace ? \
     BandMatrixView<T>(A.NonConst().ptr(),A.colsize(),A.colsize(),A.nhi(), \
       NEWHI,A.stepj(),A.stepi(),A.diagstep(),TransOf(A.stor()),A.ct() \
       FIRSTLAST1(A.NonConst().first,A.NonConst().last) ) : \
     BandMatrixViewOf(Aptr,A.colsize(),A.colsize(),A.nhi(), NEWHI,  \
       (A.nlo() == 1 && A.nhi() == 1) ? DiagMajor : ColMajor)) : \
    (inplace ? \
     BandMatrixView<T>(A.NonConst().ptr(),A.colsize(),A.colsize(),A.nlo(), \
       NEWHI,A.stepi(),A.stepj(),A.diagstep(),A.stor(),A.ct() \
       FIRSTLAST1(A.NonConst().first,A.NonConst().last) ) : \
     BandMatrixViewOf(Aptr,A.colsize(),A.colsize(),A.nlo(), NEWHI, \
       (A.nlo() == 1 && A.nhi() == 1) ? DiagMajor : ColMajor)))

  template <class T> BandLUDiv<T>::BandLUDiv_Impl::BandLUDiv_Impl(
      const GenBandMatrix<T>& A, bool _inplace) :
    istrans(A.nhi()<A.nlo() || (A.nhi()==A.nlo() && A.isrm())),
    inplace(NEWLO == 0 || (_inplace && 
	  ((A.isrm() && istrans) || (A.iscm() && !istrans) || 
	   (A.isdm() && A.nlo()==1 && A.nhi()==1)))),
    Aptr1(APTR1), Aptr(APTR), LUx(LUX),
    P(new int[A.colsize()]), logdet(0), signdet(1), donedet(false) {}

#undef LUX
#undef APTR
#undef APTR1
#undef NEWLO
#undef NEWHI

  template <class T> BandLUDiv<T>::BandLUDiv(const GenBandMatrix<T>& A,
      bool inplace) : pimpl(new BandLUDiv_Impl(A,inplace))
  {
    TMVAssert(A.IsSquare());
    if (inplace) {
      // For inplace decomposition, make sure the original band matrix
      // has room for the extra upper diagonals...
      // if iscm stepj >= (2*A.nlo()+A.nhi())
      // if isdm extra diags appear at end, so can't really check
      // if isrm stepi >= (2*A.nhi()+A.nlo())
      if (A.iscm()) {
	TMVAssert(!pimpl->istrans);
	TMVAssert(A.stepj() >= MIN(int(A.colsize()),2*A.nlo()+A.nhi()));
	TMVAssert(pimpl->LUx.Diags(-A.nlo(),A.nhi()+1) == A);
      } else if (A.isrm()) {
	TMVAssert(pimpl->istrans);
	TMVAssert(A.stepi() >= MIN(int(A.colsize()),2*A.nhi()+A.nlo()));
	TMVAssert(pimpl->LUx.Diags(-A.nhi(),A.nlo()+1).Transpose() == A);
      } else {
	TMVAssert(A.isdm());
	if (pimpl->istrans)
	  TMVAssert(pimpl->LUx.Diags(-A.nhi(),A.nlo()+1).Transpose() == A);
	else
	  TMVAssert(pimpl->LUx.Diags(-A.nlo(),A.nhi()+1) == A);
      }
    } else {
      if (pimpl->istrans) 
	BandMatrixViewOf(pimpl->LUx,A.nhi(),A.nlo()) = A.Transpose();
      else BandMatrixViewOf(pimpl->LUx,A.nlo(),A.nhi()) = A;
    }

    if (pimpl->LUx.nlo() > 0) {
      int Anhi = pimpl->istrans ? A.nlo() : A.nhi();
      if (Anhi < pimpl->LUx.nhi())
	pimpl->LUx.Diags(Anhi+1,pimpl->LUx.nhi()+1).Zero();
      LU_Decompose(pimpl->LUx,pimpl->P.get(),pimpl->signdet,Anhi);
    } else {
      pimpl->P.get()[0] = 1; // A tag to indicate that P is not set yet.
    }
  }

#define NEWLO MIN(A.nlo(),A.nhi())
#define NEWHI MIN(A.nlo()+A.nhi(),int(A.colsize())-1)
#define APTR1 \
  new T[BandStorageLength(ColMajor,A.colsize(),A.colsize(),NEWLO,NEWHI)]
#define APTR Aptr1.get()

#define LUX \
  BandMatrixViewOf(Aptr,A.colsize(),A.colsize(),NEWLO,NEWHI, \
      (A.nlo() == 1 && A.nhi() == 1) ? DiagMajor : ColMajor)

  template <class T> BandLUDiv<T>::BandLUDiv_Impl::BandLUDiv_Impl(
      const AssignableToBandMatrix<T>& A) :
    istrans(A.nhi()<A.nlo()), inplace(false),
    Aptr1(APTR1), Aptr(APTR), LUx(LUX),
    P(new int[A.colsize()]), logdet(0), signdet(1), donedet(false) {}

#undef LUX
#undef APTR
#undef APTR1
#undef NEWLO
#undef NEWHI

  template <class T> BandLUDiv<T>::BandLUDiv(
      const AssignableToBandMatrix<T>& A) : pimpl(new BandLUDiv_Impl(A))
  {
    TMVAssert(A.IsSquare());
    if (pimpl->istrans) 
      BandMatrixViewOf(pimpl->LUx,A.nhi(),A.nlo()).Transpose() = A;
    else BandMatrixViewOf(pimpl->LUx,A.nlo(),A.nhi()) = A;

    if (pimpl->LUx.nlo() > 0) {
      int Anhi = pimpl->istrans ? A.nlo() : A.nhi();
      if (Anhi < pimpl->LUx.nhi())
	pimpl->LUx.Diags(Anhi+1,pimpl->LUx.nhi()+1).Zero();
      LU_Decompose(pimpl->LUx,pimpl->P.get(),pimpl->signdet,Anhi);
    } else {
      pimpl->P.get()[0] = 1; // A tag to indicate that P is not set yet.
    }
  }

  template <class T> BandLUDiv<T>::~BandLUDiv() { delete pimpl; pimpl=0; }

  template <class T> template <class T1> void BandLUDiv<T>::DoLDivEq(
      const MatrixView<T1>& m) const
  {
    if (pimpl->istrans) LU_RDivEq(pimpl->LUx,pimpl->P.get(),m.Transpose());
    else LU_LDivEq(pimpl->LUx,pimpl->P.get(),m);
  }

  template <class T> template <class T1> void BandLUDiv<T>::DoRDivEq(
      const MatrixView<T1>& m) const
  {
    if (pimpl->istrans) LU_LDivEq(pimpl->LUx,pimpl->P.get(),m.Transpose());
    else LU_RDivEq(pimpl->LUx,pimpl->P.get(),m);
  }

  template <class T> template <class T1, class T2> void BandLUDiv<T>::DoLDiv(
      const GenMatrix<T1>& m, const MatrixView<T2>& x) const
  {
    if (pimpl->istrans) 
      LU_RDivEq(pimpl->LUx,pimpl->P.get(),(x=m).Transpose());
    else LU_LDivEq(pimpl->LUx,pimpl->P.get(),x=m);
  }

  template <class T> template <class T1, class T2> void BandLUDiv<T>::DoRDiv(
      const GenMatrix<T1>& m, const MatrixView<T2>& x) const
  {
    if (pimpl->istrans) 
      LU_LDivEq(pimpl->LUx,pimpl->P.get(),(x=m).Transpose());
    else LU_RDivEq(pimpl->LUx,pimpl->P.get(),x=m);
  }

  template <class T> T BandLUDiv<T>::Det() const
  {
    if (!pimpl->donedet) {
      T s;
      pimpl->logdet = DiagMatrixViewOf(pimpl->LUx.diag()).LogDet(&s);
      pimpl->signdet *= s;
      pimpl->donedet = true;
    }         
    if (pimpl->signdet == T(0)) return T(0);
    else return pimpl->signdet * std::exp(pimpl->logdet);  
  }                  

  template <class T> RealType(T) BandLUDiv<T>::LogDet(T* sign) const
  {
    if (!pimpl->donedet) {
      T s;
      pimpl->logdet = DiagMatrixViewOf(pimpl->LUx.diag()).LogDet(&s);
      pimpl->signdet *= s;
      pimpl->donedet = true;
    }
    if (sign) *sign = pimpl->signdet;
    return pimpl->logdet;  
  }                  

  template <class T> template <class T1> void BandLUDiv<T>::DoInverse(
      const MatrixView<T1>& minv) const
  {
    if (pimpl->istrans)
      LU_Inverse(pimpl->LUx,pimpl->P.get(),minv.Transpose());
    else
      LU_Inverse(pimpl->LUx,pimpl->P.get(),minv);
  }

  template <class T> void BandLUDiv<T>::DoInverseATA(
      const MatrixView<T>& ata) const
  {
    // See corresponding routine in TMV_LUD.cpp
    if (pimpl->istrans) {
      UpperTriMatrixView<T> uinv = UpperTriMatrixViewOf(ata);
      uinv = GetU();
      Tri_Inverse(uinv,pimpl->LUx.nhi());
      ata = uinv.Transpose() * uinv.Conjugate();
      // ata /= L.Transpose()
      // ata = LT^-1 ata
      // ataT = ataT L^-1
      LU_PackedPL_RDivEq(pimpl->LUx,pimpl->P.get(),ata.Transpose());
      // ata %= L.Conjugate()
      // ata* = ata* L^-1
      LU_PackedPL_RDivEq(pimpl->LUx,pimpl->P.get(),ata.Conjugate());
    } else {
      LowerTriMatrixView<T> linv = LowerTriMatrixViewOf(ata,UnitDiag);
      // linv = L.Inverse()
      ata.SetToIdentity();
      LU_PackedPL_LDivEq(pimpl->LUx,pimpl->P.get(),ata);
      ata = linv * linv.Adjoint();
      ConstBandMatrixView<T> U = GetU();
      // ata /= U;
      TriLDivEq(U,ata,NonUnitDiag);
      // ata %= U.Adjoint();
      // ata = ata Ut^-1
      // atat = U^-1 atat
      TriLDivEq(U,ata.Adjoint(),NonUnitDiag);
    }
  }

  template <class T> bool BandLUDiv<T>::Singular() const 
  { return Det() == T(0); }

  template <class T> bool BandLUDiv<T>::IsTrans() const 
  { return pimpl->istrans; }

  template <class T> ConstBandMatrixView<T> BandLUDiv<T>::GetU() const
  { return BandMatrixViewOf(pimpl->LUx,0,pimpl->LUx.nhi()); }

  template <class T> LowerTriMatrix<T,UnitDiag> BandLUDiv<T>::GetL() const
  {
    LowerTriMatrix<T,UnitDiag> L(pimpl->LUx.colsize());
    LU_PackedPL_Unpack(pimpl->LUx,pimpl->P.get(),L.View());
    return L;
  }

  template <class T> const GenBandMatrix<T>& BandLUDiv<T>::GetLU() const 
  { return pimpl->LUx; }

  template <class T> const int* BandLUDiv<T>::GetP() const
  {
    if (pimpl->LUx.nlo() == 0 && pimpl->LUx.colsize() > 0 && 
	pimpl->P.get()[0] == 1) {
      // Then P hasn't been set up yet.
      const int N = pimpl->LUx.colsize();
      for(int i=0;i<N;++i) pimpl->P.get()[i] = i;
    }
    return pimpl->P.get();
  }

  template <class T> bool BandLUDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, std::ostream* fout) const
  {
    Matrix<T> mm = m;
    if (fout) {
      *fout << "BandLUDiv:\n";
      *fout << "M = "<<(pimpl->istrans?mm.Transpose():mm.View())<<std::endl;
      *fout << "L = "<<GetL()<<std::endl;
      *fout << "U = "<<GetU()<<std::endl;
    }
    Matrix<T> lu = GetL() * Matrix<T>(GetU());
    lu.ReversePermuteRows(GetP());
    RealType(T) nm = Norm(lu-(pimpl->istrans ? mm.Transpose() : mm.View()));
    nm /= Norm(GetL())*Norm(GetU());
    if (fout) {
      *fout << "PLU = "<<lu<<std::endl;
      *fout << "Norm(M-PLU)/Norm(PLU) = "<<nm<<std::endl;
    }
    return nm < mm.DoCondition()*mm.colsize()*Epsilon<T>();
  }

  template <class T> size_t BandLUDiv<T>::colsize() const
  { return pimpl->LUx.colsize(); }

  template <class T> size_t BandLUDiv<T>::rowsize() const
  { return pimpl->LUx.rowsize(); }

#define InstFile "TMV_BandLUD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


