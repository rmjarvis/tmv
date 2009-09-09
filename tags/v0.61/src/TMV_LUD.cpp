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



#include "TMV_LUD.h"
#include "TMV_LUDiv.h"
#include "TMV_Matrix.h"
#include "TMV_TriMatrix.h"
#include "TMV_DiagMatrix.h"
#include "TMV_TriMatrixArith.h"
#include "TMV_MatrixArith.h"
#include <ostream>

namespace tmv {

  template <class T> struct LUDiv<T>::LUDiv_Impl
  {
    public :
      LUDiv_Impl(const GenMatrix<T>& m, bool inplace);

      const bool istrans;
      const bool inplace;
      auto_array<T> Aptr1;
      T* Aptr;
      MatrixView<T> LUx;
      auto_array<int> P;
      mutable RealType(T) logdet;
      mutable T signdet;
      mutable bool donedet;

    private :
      LUDiv_Impl(const LUDiv_Impl&) :
	istrans(false), inplace(false), Aptr1(0), Aptr(0), 
	LUx(MatrixViewOf(Aptr,0,0,ColMajor)),
	P(0), logdet(0), signdet(0), donedet(false)
      { TMVAssert(FALSE); }
      LUDiv_Impl& operator=(const LUDiv_Impl&) 
      { TMVAssert(FALSE); return *this; }

  };

#define APTR1 (inplace ? 0 : new T[A.colsize()*A.rowsize()])
#define APTR (inplace ? A.NonConst().ptr() : Aptr1.get())
#define LUX (istrans ? \
    (inplace ? A.NonConst().Transpose() : \
     MatrixViewOf(Aptr,A.rowsize(),A.colsize(),ColMajor) ): \
    (inplace ? A.NonConst().View() : \
     MatrixViewOf(Aptr,A.colsize(),A.rowsize(),ColMajor)))

  template <class T> LUDiv<T>::LUDiv_Impl::LUDiv_Impl(
      const GenMatrix<T>& A, bool _inplace) :
    istrans(A.isrm()), inplace(_inplace && (A.iscm() || A.isrm())),
    Aptr1(APTR1), Aptr(APTR), LUx(LUX), P(new int[A.colsize()]),
    logdet(0), signdet(1), donedet(false) {}

#undef LUX
#undef APTR

  template <class T> LUDiv<T>::LUDiv(const GenMatrix<T>& A,
      bool inplace) :
    pimpl(new LUDiv_Impl(A,inplace)) 
  {
    TMVAssert(A.IsSquare());
    if (pimpl->istrans) {
      if (pimpl->inplace) TMVAssert(A.Transpose() == pimpl->LUx);
      else pimpl->LUx = A.Transpose();
    }
    else {
      if (inplace) TMVAssert(A == pimpl->LUx);
      else pimpl->LUx = A;
    }
    LU_Decompose(pimpl->LUx,pimpl->P.get(),pimpl->signdet);
  }

  template <class T> LUDiv<T>::~LUDiv() { delete pimpl; pimpl=0; }

  template <class T> template <class T1> void LUDiv<T>::DoLDivEq(
      const MatrixView<T1>& m) const
  {
    TMVAssert(pimpl->LUx.colsize() == m.colsize());
    if (pimpl->istrans) LU_RDivEq(pimpl->LUx,pimpl->P.get(),m.Transpose());
    else LU_LDivEq(pimpl->LUx,pimpl->P.get(),m);
  }

  template <class T> template <class T1> void LUDiv<T>::DoRDivEq(
      const MatrixView<T1>& m) const
  {
    TMVAssert(pimpl->LUx.colsize() == m.rowsize());
    if (pimpl->istrans) LU_LDivEq(pimpl->LUx,pimpl->P.get(),m.Transpose());
    else LU_RDivEq(pimpl->LUx,pimpl->P.get(),m);
  }

  template <class T> template <class T1, class T2> void LUDiv<T>::DoLDiv(
      const GenMatrix<T1>& m1, const MatrixView<T2>& m0) const
  {
    TMVAssert(m1.colsize() == m0.colsize());
    TMVAssert(m1.rowsize() == m0.rowsize());
    TMVAssert(pimpl->LUx.colsize() == m1.colsize());
    DoLDivEq(m0=m1);
  }

  template <class T> template <class T1, class T2> void LUDiv<T>::DoRDiv(
      const GenMatrix<T1>& m1, const MatrixView<T2>& m0) const
  {
    TMVAssert(m1.colsize() == m0.colsize());
    TMVAssert(m1.rowsize() == m0.rowsize());
    TMVAssert(pimpl->LUx.colsize() == m1.rowsize());
    DoRDivEq(m0=m1);
  }

  template <class T> T LUDiv<T>::Det() const
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

  template <class T> RealType(T) LUDiv<T>::LogDet(T* sign) const
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
  template <class T> template <class T1> void LUDiv<T>::DoInverse(
      const MatrixView<T1>& minv) const
  {
    TMVAssert(minv.colsize() == pimpl->LUx.colsize());
    TMVAssert(minv.rowsize() == pimpl->LUx.colsize());
    // m = P L U
    // m^-1 = U^-1 L^-1 Pt
    if (pimpl->istrans) LU_Inverse(pimpl->LUx,pimpl->P.get(),minv.Transpose());
    else LU_Inverse(pimpl->LUx,pimpl->P.get(),minv);
  }

  template <class T> void LUDiv<T>::DoInverseATA(
      const MatrixView<T>& ata) const
  {
    TMVAssert(ata.colsize() == pimpl->LUx.colsize());
    TMVAssert(ata.rowsize() == pimpl->LUx.colsize());
    // (At A)^-1 = A^-1 (A^-1)t
    // = (U^-1 L^-1 Pt) (P L^-1t U^-1t)
    // = U^-1 L^-1 L^-1t U^-1t
    //
    // if PLU is really AT, then
    // A^-1 = P L^-1T U^-1T
    // (At A)^-1 = P L^-1T U^-1T U^-1* L^-1* Pt

    LowerTriMatrixView<T> L = LowerTriMatrixViewOf(pimpl->LUx,UnitDiag);
    UpperTriMatrixView<T> U = UpperTriMatrixViewOf(pimpl->LUx);

    if (pimpl->istrans) {
      UpperTriMatrixView<T> uinv = UpperTriMatrixViewOf(ata);
      uinv = U.Inverse();
      ata = uinv.Transpose() * uinv.Conjugate();
      ata /= L.Transpose();
      ata %= L.Conjugate();
      ata.ReversePermuteCols(pimpl->P.get());
      ata.ReversePermuteRows(pimpl->P.get());
    } else {
      LowerTriMatrixView<T> linv = LowerTriMatrixViewOf(ata,UnitDiag);
      linv = L.Inverse();
      ata = linv * linv.Adjoint();
      ata /= U;
      ata %= U.Adjoint();
    }
  }

  template <class T> bool LUDiv<T>::Singular() const 
  { return Det() == T(0); }

  template <class T> bool LUDiv<T>::IsTrans() const 
  { return pimpl->istrans; }

  template <class T> ConstLowerTriMatrixView<T> LUDiv<T>::GetL() const 
  { return LowerTriMatrixViewOf(pimpl->LUx,UnitDiag); }

  template <class T> ConstUpperTriMatrixView<T> LUDiv<T>::GetU() const 
  { return UpperTriMatrixViewOf(pimpl->LUx,NonUnitDiag); }

  template <class T> const GenMatrix<T>& LUDiv<T>::GetLU() const 
  { return pimpl->LUx; }

  template <class T> const int* LUDiv<T>::GetP() const 
  { return pimpl->P.get(); }

  template <class T> bool LUDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, std::ostream* fout) const
  {
    Matrix<T> mm = m;
    if (fout) {
      *fout << "LUDiv:\n";
      *fout << "M = "<<(pimpl->istrans?mm.Transpose():mm.View())<<std::endl;
      *fout << "L = "<<GetL()<<std::endl;
      *fout << "U = "<<GetU()<<std::endl;
    }
    Matrix<T> lu = GetL()*GetU();
    lu.ReversePermuteRows(GetP());
    RealType(T) nm = Norm(lu-(pimpl->istrans ? mm.Transpose() : mm.View()));
    nm /= Norm(GetL())*Norm(GetU());
    if (fout) {
      *fout << "PLU = "<<lu<<std::endl;
      *fout << "Norm(M-PLU)/Norm(PLU) = "<<nm<<std::endl;
    }
    return nm < mm.DoCondition()*mm.colsize()*Epsilon<T>();
  }

  template <class T> size_t LUDiv<T>::colsize() const
  { return pimpl->LUx.colsize(); }

  template <class T> size_t LUDiv<T>::rowsize() const
  { return pimpl->LUx.rowsize(); }

#define InstFile "TMV_LUD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


