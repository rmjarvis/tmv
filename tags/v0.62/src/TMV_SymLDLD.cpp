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



#include "tmv/TMV_SymLDLD.h"
#include "TMV_SymLDLDiv.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "TMV_SymSquare.h"
#include "tmv/TMV_TriMatrixArith.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_BandMatrixArith.h"
#include <ostream>

namespace tmv {

  template <class T> struct SymLDLDiv<T>::SymLDLDiv_Impl
  {
  public :
    SymLDLDiv_Impl(const GenSymMatrix<T>& m, bool inplace);

    const bool inplace;
    auto_array<T> Aptr1;
    T* Aptr;
    SymMatrixView<T> LLx;
    Vector<T> xD;
    auto_array<int> P;
    mutable RealType(T) logdet;
    mutable T signdet;
  };

#define APTR1 (inplace ? 0 : new T[A.size()*A.size()])
#define APTR (inplace ? A.NonConst().ptr() : Aptr1.get())
#define LLX \
  (inplace ? \
   (A.uplo()==Upper ? \
    (A.isherm() ? A.NonConst().Adjoint() : A.NonConst().Transpose()) : \
    A.NonConst()) : \
   (A.isherm() ? \
    HermMatrixViewOf(Aptr,A.size(),Lower,BaseStorOf(A.LowerTri())) : \
    SymMatrixViewOf(Aptr,A.size(),Lower,BaseStorOf(A.LowerTri()))))

  template <class T> SymLDLDiv<T>::SymLDLDiv_Impl::SymLDLDiv_Impl(
      const GenSymMatrix<T>& A, bool _inplace) :
    inplace(_inplace && (A.isrm() || A.iscm())), 
    Aptr1(APTR1), Aptr(APTR), LLx(LLX), xD(A.size()-1),
    P(new int[A.colsize()]), logdet(0), signdet(1) {}

#undef APTR1
#undef APTR
#undef LLX

  template <class T> SymLDLDiv<T>::SymLDLDiv(const GenSymMatrix<T>& A,
      bool inplace) : pimpl(new SymLDLDiv_Impl(A,inplace))
  {
#ifdef XTEST
    TMVAssert(A.HermOK());
#endif
    if (inplace) TMVAssert(A == pimpl->LLx); 
    else pimpl->LLx = A;

    LDL_Decompose(pimpl->LLx,pimpl->xD.View(),pimpl->P.get(),
        pimpl->logdet,pimpl->signdet);
#ifdef XTEST
    TMVAssert(pimpl->LLx.HermOK());
#endif
  }

  template <class T> SymLDLDiv<T>::~SymLDLDiv() {}

  template <class T> template <class T1> void SymLDLDiv<T>::DoLDivEq(
      const MatrixView<T1>& m) const
  { LDL_LDivEq(pimpl->LLx,pimpl->xD,pimpl->P.get(),m); }

  template <class T> template <class T1> void SymLDLDiv<T>::DoRDivEq(
      const MatrixView<T1>& m) const
  { LDL_RDivEq(pimpl->LLx,pimpl->xD,pimpl->P.get(),m); }

  template <class T> template <class T1, class T2> void SymLDLDiv<T>::DoLDiv(
      const GenMatrix<T1>& m1, const MatrixView<T2>& m0) const
  { LDL_LDivEq(pimpl->LLx,pimpl->xD,pimpl->P.get(),m0=m1); }

  template <class T> template <class T1, class T2> void SymLDLDiv<T>::DoRDiv(
      const GenMatrix<T1>& m1, const MatrixView<T2>& m0) const
  { LDL_RDivEq(pimpl->LLx,pimpl->xD,pimpl->P.get(),m0=m1); }

  template <class T> T SymLDLDiv<T>::Det() const 
  { 
    if (pimpl->signdet == T(0)) return T(0);
    else return pimpl->signdet * EXP(pimpl->logdet);  
  }

  template <class T> RealType(T) SymLDLDiv<T>::LogDet(T* sign) const 
  { 
    if (sign) *sign = pimpl->signdet;
    return pimpl->logdet;
  }

  template <class T> template <class T1> void SymLDLDiv<T>::DoInverse(
      const SymMatrixView<T1>& sinv) const
  { 
    TMVAssert(IsReal(T()) || issym() == sinv.issym());
    TMVAssert(IsReal(T()) || isherm() == sinv.isherm());
    LDL_Inverse(pimpl->LLx,pimpl->xD,pimpl->P.get(),sinv); 
#ifdef XTEST
    TMVAssert(sinv.HermOK());
#endif
  }

  template <class T> template <class T1> void SymLDLDiv<T>::DoInverse(
      const MatrixView<T1>& minv) const
  {
    if (pimpl->LLx.isherm()) {
      DoInverse(HermMatrixViewOf(minv,Lower));
      if (minv.colsize() > 1)
        minv.UpperTri().OffDiag() =
        minv.LowerTri().OffDiag().Adjoint();
    } else {
      DoInverse(SymMatrixViewOf(minv,Lower));
      if (minv.colsize() > 1)
        minv.UpperTri().OffDiag() =
        minv.LowerTri().OffDiag().Transpose();
    }
  }

  template <class T> void SymLDLDiv<T>::DoInverseATA(
      const MatrixView<T>& ata) const
  {
    TMVAssert(ata.colsize() == pimpl->LLx.size());
    TMVAssert(ata.rowsize() == pimpl->LLx.size());

    if (pimpl->LLx.isherm()) {
      SymMatrixView<T> symata = HermMatrixViewOf(ata,Lower);
      Inverse(symata);
      SymSquare<true>(ata);
    } else {
      SymMatrixView<T> symata = SymMatrixViewOf(ata,Lower);
      Inverse(symata);
      SymSquare<false>(ata);
    }
  }

  template <class T> bool SymLDLDiv<T>::Singular() const 
  { return Det() == T(0); }

  template <class T> const ConstLowerTriMatrixView<T> SymLDLDiv<T>::GetL() const
  { return pimpl->LLx.LowerTri().ViewAsUnitDiag(); }

  template <class T> const BandMatrix<T> SymLDLDiv<T>::GetD() const 
  { 
    BandMatrix<T> temp(pimpl->LLx.size(),pimpl->LLx.size(),1,1);
    temp.diag() = pimpl->LLx.diag();
    temp.diag(-1) = pimpl->xD;
    temp.diag(1) = (pimpl->LLx.isherm() ? pimpl->xD.Conjugate() : 
        pimpl->xD.View());
    return temp;
  }

  template <class T> const int* SymLDLDiv<T>::GetP() const 
  { return pimpl->P.get(); }

  template <class T> const GenSymMatrix<T>& SymLDLDiv<T>::GetLL() const 
  { return pimpl->LLx; }

  template <class T> const GenVector<T>& SymLDLDiv<T>::GetxD() const 
  { return pimpl->xD; }

  template <class T> bool SymLDLDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, std::ostream* fout) const
  {
    Matrix<T> mm = m;
    if (fout) {
      *fout << "SymLDLDiv:\n";
      *fout << "M = "<<mm<<std::endl;
      *fout << "L = "<<GetL()<<std::endl;
      *fout << "D = "<<GetD()<<std::endl;
      *fout << "P = ";
      for(int i=0;i<int(mm.colsize());i++) *fout<<pimpl->P.get()[i]<<" ";
      *fout <<std::endl;
    }
    Matrix<T> lu = GetL()*GetD()*
    (pimpl->LLx.isherm()?GetL().Adjoint():GetL().Transpose());
    lu.ReversePermuteRows(pimpl->P.get());
    lu.ReversePermuteCols(pimpl->P.get());
    RealType(T) nm = Norm(lu-mm);
    nm /= SQR(Norm(GetL()))*Norm(GetD());
    if (fout) {
      *fout << "LDLt = "<<lu<<std::endl;
      *fout << "Norm(M-LDLt)/Norm(LDLt) = "<<nm<<std::endl;
    }
    return nm < mm.DoCondition()*mm.colsize()*Epsilon<T>();
  }

  template <class T> size_t SymLDLDiv<T>::colsize() const
  { return pimpl->LLx.size(); }

  template <class T> size_t SymLDLDiv<T>::rowsize() const
  { return pimpl->LLx.size(); }

  template <class T> bool SymLDLDiv<T>::isherm() const
  { return pimpl->LLx.isherm(); }

  template <class T> bool SymLDLDiv<T>::issym() const
  { return pimpl->LLx.issym(); }

#define InstFile "TMV_SymLDLD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


