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



#include "TMV_Matrix.h"
#include "TMV_SymMatrix.h"
#include "TMV_SymLUDiv.h"
#include "TMV_SymCHDiv_B.h"
#include "TMV_TriMatrixArith.h"
#include "TMV_SymMatrixArith.h"
#include "TMV_MatrixArith.h"
#include <ostream>

namespace tmv {

  template <class T> struct SymLUDiv<T>::SymLUDiv_Impl
  {
    SymLUDiv_Impl(const GenSymMatrix<T>& m, bool inplace);

    const bool inplace;
    auto_array<T> Aptr1;
    T* Aptr;
    SymMatrixView<T> LLx;
    Vector<T> xD;
    auto_array<size_t> P;
    mutable T det;
  };

#define APTR1 inplace ? 0 : new T[A.size()*A.size()]
#define APTR inplace ? A.NonConst().ptr() : Aptr1.get()
#define LLX \
  inplace ? \
    A.uplo()==Upper ? \
      (A.isherm() ? A.NonConst().Adjoint() : A.NonConst().Transpose()) : \
      A.NonConst() : \
    A.isherm() ? \
      HermMatrixViewOf(Aptr,A.size(),Lower,BaseStorOf(A.LowerTri())) : \
      SymMatrixViewOf(Aptr,A.size(),Lower,BaseStorOf(A.LowerTri()))

  template <class T> SymLUDiv<T>::SymLUDiv_Impl::SymLUDiv_Impl(
      const GenSymMatrix<T>& A, bool _inplace) :
    inplace(_inplace && (A.isrm() || A.iscm())), 
    Aptr1(APTR1), Aptr(APTR), LLx(LLX), 
    xD(Vector<T>(A.size()-1)),
    P(new size_t[A.colsize()]), det(T(1)) {}

#undef APTR1
#undef APTR
#undef LLX

  template <class T> SymLUDiv<T>::SymLUDiv(const GenSymMatrix<T>& A,
      bool inplace) : pimpl(new SymLUDiv_Impl(A,inplace))
  {
#ifdef XTEST
    TMVAssert(A.HermOK());
#endif
    if (inplace) TMVAssert(A == pimpl->LLx); 
    else pimpl->LLx = A;

    SymLU_Decompose(pimpl->LLx,pimpl->xD.View(),pimpl->P.get(),pimpl->det);
#ifdef XTEST
    TMVAssert(pimpl->LLx.HermOK());
#endif
  }

  template <class T> SymLUDiv<T>::~SymLUDiv() { delete pimpl; }

#define CT std::complex<T>
  template <class T> inline void SymLU_LDivEq(
      const GenSymMatrix<CT>& , const GenVector<CT>& , const size_t* ,
      const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void SymLU_RDivEq(
      const GenSymMatrix<CT>& , const GenVector<CT>& , const size_t* ,
      const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void SymLU_Inverse(
      const GenSymMatrix<CT>& , const GenVector<CT>& , const size_t* ,
      const SymMatrixView<T>& )
  { TMVAssert(FALSE); }
#undef CT

  template <class T> template <class T1> void SymLUDiv<T>::DoLDivEq(
      const MatrixView<T1>& m) const
  { SymLU_LDivEq(pimpl->LLx,pimpl->xD,pimpl->P.get(),m); }

  template <class T> template <class T1> void SymLUDiv<T>::DoRDivEq(
      const MatrixView<T1>& m) const
  { SymLU_RDivEq(pimpl->LLx,pimpl->xD,pimpl->P.get(),m); }

  template <class T> template <class T1, class T2> void SymLUDiv<T>::DoLDiv(
      const GenMatrix<T1>& m1, const MatrixView<T2>& m0) const
  { SymLU_LDivEq(pimpl->LLx,pimpl->xD,pimpl->P.get(),m0=m1); }

  template <class T> template <class T1, class T2> void SymLUDiv<T>::DoRDiv(
      const GenMatrix<T1>& m1, const MatrixView<T2>& m0) const
  { SymLU_RDivEq(pimpl->LLx,pimpl->xD,pimpl->P.get(),m0=m1); }

  template <class T> T SymLUDiv<T>::Det() const 
  { return pimpl->det; }

  template <class T> template <class T1> void SymLUDiv<T>::DoInverse(
      const SymMatrixView<T1>& sinv) const
  { 
    TMVAssert(IsReal(T()) || issym() == sinv.issym());
    TMVAssert(IsReal(T()) || isherm() == sinv.isherm());
    SymLU_Inverse(pimpl->LLx,pimpl->xD,pimpl->P.get(),sinv); 
#ifdef XTEST
    TMVAssert(sinv.HermOK());
#endif
  }

  template <class T> template <class T1> void SymLUDiv<T>::DoInverse(
      const MatrixView<T1>& minv) const
  {
    if (pimpl->LLx.isherm()) {
      DoInverse(HermMatrixViewOf(minv,Lower));
      if (minv.colsize() > 1)
	UpperTriMatrixViewOf(minv).OffDiag() =
	  LowerTriMatrixViewOf(minv).OffDiag().Adjoint();
    } else {
      DoInverse(SymMatrixViewOf(minv,Lower));
      if (minv.colsize() > 1)
	UpperTriMatrixViewOf(minv).OffDiag() =
	  LowerTriMatrixViewOf(minv).OffDiag().Transpose();
    }
  }

  template <class T> void SymLUDiv<T>::DoInverseATA(
      const MatrixView<T>& ata) const
  {
    TMVAssert(ata.colsize() == pimpl->LLx.size());
    TMVAssert(ata.rowsize() == pimpl->LLx.size());

    if (pimpl->LLx.isherm()) {
      SymMatrixView<T> symata = HermMatrixViewOf(ata,Lower);
      Inverse(symata);
      SymATASquare<true>(ata);
    } else {
      SymMatrixView<T> symata = SymMatrixViewOf(ata,Lower);
      Inverse(symata);
      SymATASquare<false>(ata);
    }
  }

  template <class T> bool SymLUDiv<T>::Singular() const 
  { return Det() == T(0); }

  template <class T> const ConstLowerTriMatrixView<T> SymLUDiv<T>::GetL() const
  { return pimpl->LLx.LowerTri().MakeUnitDiag(); }

  template <class T> const Matrix<T> SymLUDiv<T>::GetD() const 
  { 
    Matrix<T> temp(pimpl->LLx.size(),pimpl->LLx.size(),T(0));
    temp.diag() = pimpl->LLx.diag();
    temp.diag(-1) = pimpl->xD;
    temp.diag(1) = pimpl->LLx.isherm() ? pimpl->xD.Conjugate() : 
      pimpl->xD.View();
    return temp;
  }

  template <class T> const size_t* SymLUDiv<T>::GetP() const 
  { return pimpl->P.get(); }

  template <class T> const GenSymMatrix<T>& SymLUDiv<T>::GetLL() const 
  { return pimpl->LLx; }

  template <class T> const GenVector<T>& SymLUDiv<T>::GetxD() const 
  { return pimpl->xD; }

  template <class T> std::string SymLUDiv<T>::Type() const
  { return std::string("SymLUDiv<") + tmv::Type(T()) + ">"; }

  template <class T> DivType SymLUDiv<T>::GetDivType() const 
  { return LU; }

  template <class T> bool SymLUDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, std::ostream* fout) const
  {
    Matrix<T> mm = m;
    if (fout) {
      *fout << "SymLUDiv:\n";
      *fout << "M = "<<mm<<std::endl;
      *fout << "L = "<<GetL()<<std::endl;
      *fout << "D = "<<GetD()<<std::endl;
      *fout << "P = ";
      for(size_t i=0;i<mm.colsize();i++) *fout<<pimpl->P.get()[i]<<" ";
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
    mm.DivideUsing(SVS);
    return nm < mm.Condition()*mm.colsize()*Epsilon<T>();
  }

  template <class T> size_t SymLUDiv<T>::colsize() const
  { return pimpl->LLx.size(); }

  template <class T> size_t SymLUDiv<T>::rowsize() const
  { return pimpl->LLx.size(); }

  template <class T> bool SymLUDiv<T>::isherm() const
  { return pimpl->LLx.isherm(); }

  template <class T> bool SymLUDiv<T>::issym() const
  { return pimpl->LLx.issym(); }

#define InstFile "TMV_SymLUDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


