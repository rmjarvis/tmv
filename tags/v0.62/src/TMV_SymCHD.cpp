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



#include "tmv/TMV_SymCHD.h"
#include "TMV_SymCHDiv.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "TMV_SymSquare.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_TriMatrixArith.h"
#include "tmv/TMV_MatrixArith.h"
#include <ostream>

namespace tmv {

  template <class T> struct HermCHDiv<T>::HermCHDiv_Impl
  {
  public :
    HermCHDiv_Impl(const GenSymMatrix<T>& m, bool inplace);

    const bool inplace;
    auto_array<T> Aptr1;
    T* Aptr;
    SymMatrixView<T> LLx;
    mutable bool zerodet;
    mutable RealType(T) logdet;
    mutable bool donedet;
  };

#define APTR1 (inplace ? 0 : new T[A.size()*A.size()])
#define APTR (inplace ? A.NonConst().ptr() : Aptr1.get())
#define LLX \
  (inplace ? (A.uplo()==Upper ? A.NonConst().Adjoint() : A.NonConst()) : \
   HermMatrixViewOf(Aptr,A.size(),Lower,BaseStorOf(A.LowerTri())))

  template <class T> HermCHDiv<T>::HermCHDiv_Impl::HermCHDiv_Impl(
      const GenSymMatrix<T>& A, bool _inplace) :
    inplace(_inplace && (A.iscm() || A.isrm())), 
    Aptr1(APTR1), Aptr(APTR), LLx(LLX), 
    zerodet(false), logdet(0),donedet(false) {}

#undef APTR
#undef APTR1
#undef LLX

  template <class T> HermCHDiv<T>::HermCHDiv(const GenSymMatrix<T>& A,
      bool inplace) : pimpl(new HermCHDiv_Impl(A,inplace))
  {
#ifdef XTEST
    TMVAssert(A.HermOK());
#endif
    TMVAssert(IsReal(T()) || A.isherm());
    if (inplace) TMVAssert(A==pimpl->LLx); 
    else pimpl->LLx = A;
#ifndef NOTHROW
    try {
#endif
      CH_Decompose(pimpl->LLx);
#ifndef NOTHROW
    }
    catch (NonPosDef)
    {
      if (inplace) throw NonPosDefHermMatrix<T>(pimpl->LLx);
      else throw NonPosDefHermMatrix2<T>(pimpl->LLx,A);
    }
#endif
  }

  template <class T> HermCHDiv<T>::~HermCHDiv() {}

  template <class T> template <class T1> void HermCHDiv<T>::DoLDivEq(
      const MatrixView<T1>& m) const
  {
    TMVAssert(pimpl->LLx.size() == m.colsize());
    CH_LDivEq(pimpl->LLx,m);
  }

  template <class T> template <class T1> void HermCHDiv<T>::DoRDivEq(
      const MatrixView<T1>& m) const
  {
    TMVAssert(pimpl->LLx.size() == m.rowsize());
    CH_RDivEq(pimpl->LLx,m);
  }

  template <class T> template <class T1, class T2> void HermCHDiv<T>::DoLDiv(
      const GenMatrix<T1>& m1, const MatrixView<T2>& m0) const
  {
    TMVAssert(m1.colsize() == m0.colsize());
    TMVAssert(m1.rowsize() == m0.rowsize());
    TMVAssert(pimpl->LLx.size() == m1.colsize());
    CH_LDivEq(pimpl->LLx,m0=m1);
  }

  template <class T> template <class T1, class T2> void HermCHDiv<T>::DoRDiv(
      const GenMatrix<T1>& m1, const MatrixView<T2>& m0) const
  {
    TMVAssert(m1.colsize() == m0.colsize());
    TMVAssert(m1.rowsize() == m0.rowsize());
    TMVAssert(pimpl->LLx.size() == m1.rowsize());
    CH_RDivEq(pimpl->LLx,m0=m1);
  }

  template <class T> T HermCHDiv<T>::Det() const
  {
    if (!pimpl->donedet) {
      T s;
      pimpl->logdet = DiagMatrixViewOf(pimpl->LLx.diag()).LogDet(&s);
      pimpl->logdet *= RealType(T)(2);
      pimpl->zerodet = s == T(0);
      pimpl->donedet = true;
    }         
    if (pimpl->zerodet) return T(0);
    else return EXP(pimpl->logdet);  
  }                  

  template <class T> RealType(T) HermCHDiv<T>::LogDet(T* sign) const
  {
    if (!pimpl->donedet) {
      T s;
      pimpl->logdet = DiagMatrixViewOf(pimpl->LLx.diag()).LogDet(&s);
      pimpl->logdet *= RealType(T)(2);
      pimpl->zerodet = s == T(0);
      pimpl->donedet = true;
    }
    if (sign) {
      if (pimpl->zerodet) *sign = T(0);
      else *sign = T(1);
    }
    return pimpl->logdet;  
  }                  

  template <class T> template <class T1> void HermCHDiv<T>::DoInverse(
      const SymMatrixView<T1>& sinv) const
  {
    TMVAssert(sinv.size() == pimpl->LLx.size());
    TMVAssert(sinv.isherm());
    CH_Inverse(pimpl->LLx,sinv);
  }

  template <class T> template <class T1> void HermCHDiv<T>::DoInverse(
      const MatrixView<T1>& minv) const
  {
    TMVAssert(minv.colsize() == pimpl->LLx.size());
    TMVAssert(minv.rowsize() == pimpl->LLx.size());

    if (IsComplex(T1())) minv.diag().Imag().Zero();
    Inverse(HermMatrixViewOf(minv,Lower));
    if (minv.colsize() > 1)
      minv.UpperTri().OffDiag() = minv.LowerTri().OffDiag().Adjoint();
  }

  template <class T> void HermCHDiv<T>::DoInverseATA(
      const MatrixView<T>& ata) const
  {
    // ata = (At A)^-1 = A^-1 (A^-1)t
    //     = A^-1 A^-1
    SymMatrixView<T> hermata = HermMatrixViewOf(ata,Lower);
    DoInverse(hermata);
    SymSquare<true>(ata);
  }

  template <class T> bool HermCHDiv<T>::Singular() const 
  { return Det() == T(0); }

  template <class T> 
  const ConstLowerTriMatrixView<T> HermCHDiv<T>::GetL() const 
  { return pimpl->LLx.LowerTri(); }

  template <class T> const GenSymMatrix<T>& HermCHDiv<T>::GetLL() const 
  { return pimpl->LLx; }

  template <class T> bool HermCHDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, std::ostream* fout) const
  {
    Matrix<T> mm = m;
    if (fout) {
      *fout << "HermCHDiv:\n";
      *fout << "M = "<<mm<<std::endl;
      *fout << "L = "<<GetL()<<std::endl;
    }
    Matrix<T> lu = GetL()*GetL().Adjoint();
    RealType(T) nm = Norm(lu-mm);
    nm /= SQR(Norm(GetL()));
    if (fout) {
      *fout << "LLt = "<<lu<<std::endl;
      *fout << "Norm(M-LLt)/Norm(LLt) = "<<nm<<std::endl;
    }
    return nm < mm.DoCondition()*mm.colsize()*Epsilon<T>();
  }

  template <class T> size_t HermCHDiv<T>::colsize() const
  { return pimpl->LLx.size(); }

  template <class T> size_t HermCHDiv<T>::rowsize() const
  { return pimpl->LLx.size(); }

#define InstFile "TMV_SymCHD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv
