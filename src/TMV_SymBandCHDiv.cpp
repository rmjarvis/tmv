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



#include "TMV_SymBandMatrix.h"
#include "TMV_SymBandCHDiv.h"
#include "TMV_SymBandCHDiv_A.h"
#include "TMV_SymCHDiv_B.h"
#include "TMV_DiagMatrix.h"
#include "TMV_BandMatrixArith.h"
#include "TMV_DiagMatrixArith.h"
#include "TMV_SymMatrixArith.h"
#include "TMV_SymBandMatrixArith.h"
#include "TMV_MatrixArith.h"
#include "TMV_BandLUDiv.h"
#include <ostream>

namespace tmv {

  template <class T> struct HermBandCHDiv<T>::HermBandCHDiv_Impl
  {
    HermBandCHDiv_Impl(const GenSymBandMatrix<T>& m, bool inplace);

    const bool inplace;
    auto_array<T> Aptr1;
    T* Aptr;
    SymBandMatrixView<T> LLx;
    mutable T det;
    mutable bool donedet;
  };

#define APTR1 inplace ? 0 : \
  new T[BandStorageLength(ColMajor,A.size(),A.size(),A.nlo(),A.nlo())]
#define APTR inplace ? A.NonConst().ptr() : Aptr1.get()
#define LLX \
  inplace ? A.uplo()==Upper ? A.NonConst().Adjoint() : A.NonConst() : \
  HermBandMatrixViewOf(Aptr,A.size(),A.nlo(),Lower, \
      A.nlo() == 1 ? DiagMajor : ColMajor)

  template <class T> HermBandCHDiv<T>::HermBandCHDiv_Impl::HermBandCHDiv_Impl(
      const GenSymBandMatrix<T>& A, bool _inplace) :
    inplace((_inplace && (((A.iscm() || A.isrm()) && A.nlo()>1) || 
	  (A.isdm() && A.nlo()==1))) || A.nlo()==0 ),
    Aptr1(APTR1), Aptr(APTR), LLx(LLX), det(T(0)), donedet(false) {}

#undef APTR
#undef APTR1
#undef LLX

  template <class T> HermBandCHDiv<T>::HermBandCHDiv(
      const GenSymBandMatrix<T>& A, bool inplace) :
    pimpl(new HermBandCHDiv_Impl(A,inplace))
  {
#ifdef XTEST
    TMVAssert(A.HermOK());
#endif
    TMVAssert(IsReal(T()) || A.isherm());
    if (inplace) TMVAssert(A==pimpl->LLx); 
    else pimpl->LLx = A;
    try {
      if (A.nlo() > 0)
	HermBandCH_Decompose(pimpl->LLx);
      else {
	if (A.Real().diag().MinElement() <= RealType(T)(0)) throw NonPosDef();
      }
    }
    catch (NonPosDef)
    {
      if (inplace) throw NonPosDefHermBandMatrix<T>(pimpl->LLx);
      else throw NonPosDefHermBandMatrix2<T>(pimpl->LLx,A);
    }
#ifdef XTEST
    TMVAssert(pimpl->LLx.HermOK());
#endif
  }

  template <class T> HermBandCHDiv<T>::~HermBandCHDiv() { delete pimpl; }

#define CT std::complex<T>
  template <class T> inline void HermBandCH_LDivEq(
      const GenSymBandMatrix<CT>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void HermBandCH_RDivEq(
      const GenSymBandMatrix<CT>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void HermBandCH_Inverse(
      const GenSymBandMatrix<CT>& , const SymMatrixView<T>& )
  { TMVAssert(FALSE); }
#undef CT

  template <class T> template <class T1> void HermBandCHDiv<T>::DoLDivEq(
      const MatrixView<T1>& m) const
  {
    TMVAssert(pimpl->LLx.size() == m.colsize());
    HermBandCH_LDivEq(pimpl->LLx,m);
  }

  template <class T> template <class T1> void HermBandCHDiv<T>::DoRDivEq(
      const MatrixView<T1>& m) const
  {
    TMVAssert(pimpl->LLx.size() == m.rowsize());
    HermBandCH_RDivEq(pimpl->LLx,m);
  }

  template <class T> template <class T1, class T2> void HermBandCHDiv<T>::DoLDiv(
      const GenMatrix<T1>& m1, const MatrixView<T2>& m0) const
  {
    TMVAssert(m1.colsize() == m0.colsize());
    TMVAssert(m1.rowsize() == m0.rowsize());
    TMVAssert(pimpl->LLx.size() == m1.colsize());
    HermBandCH_LDivEq(pimpl->LLx,m0=m1);
  }

  template <class T> template <class T1, class T2> void HermBandCHDiv<T>::DoRDiv(
      const GenMatrix<T1>& m1, const MatrixView<T2>& m0) const
  {
    TMVAssert(m1.colsize() == m0.colsize());
    TMVAssert(m1.rowsize() == m0.rowsize());
    TMVAssert(pimpl->LLx.size() == m1.rowsize());
    HermBandCH_RDivEq(pimpl->LLx,m0=m1);
  }

  template <class T1, class T2> void HermBandCH_LDivEq(
      const GenSymBandMatrix<T1>& LL, const MatrixView<T2>& m)
  {
    TMVAssert(LL.size() == m.colsize());
    if (LL.nlo() == 0) {
      m /= DiagMatrixViewOf(LL.diag().Real());
    } else if (LL.nlo() == 1) {
      // m = (LDLt)^-1 m
      //   = Lt^-1 D^-1 L^-1 m

      // m /= LL.LowerBand(UnitDiag)
      BandTriLDivEq(LL.LowerBand(),m,UnitDiag);
      m /= DiagMatrixViewOf(LL.diag().Real());
      // m /= LL.UpperBand(UnitDiag)
      BandTriLDivEq(LL.UpperBand(),m,UnitDiag);
    } else {
      // m = (LLt)^-1 m
      //   = Lt^-1 L^-1 m
      m /= LL.LowerBand();
      m /= LL.UpperBand();
    }
  }

  template <class T1, class T2> void HermBandCH_RDivEq(
      const GenSymBandMatrix<T1>& LL, const MatrixView<T2>& m)
  {
    TMVAssert(LL.size() == m.rowsize());
    if (LL.nlo() == 0) {
      m %= DiagMatrixViewOf(LL.diag().Real());
    } else if (LL.nlo() == 1) {
      // m = m (LDLt)^-1
      //   = m Lt^-1 D^-1 L^-1
      // mt = Lt^-1 D^-1 L^-1 mt

      // mt /= LL.LowerBand(UnitDiag)
      BandTriLDivEq(LL.LowerBand(),m.Adjoint(),UnitDiag);
      m %= DiagMatrixViewOf(LL.diag().Real());
      // mt /= LL.UpperBand(UnitDiag)
      BandTriLDivEq(LL.UpperBand(),m.Adjoint(),UnitDiag);
    } else {
      // m = m (LLt)^-1
      //   = m Lt^-1 L^-1
      // mt = Lt^-1 L^-1 mt
      m.Adjoint() /= LL.LowerBand();
      m.Adjoint() /= LL.UpperBand();
    }
  }

  template <class T> T HermBandCHDiv<T>::Det() const
  {
    if (!pimpl->donedet) {
      pimpl->det = DiagMatrixViewOf(pimpl->LLx.diag()).Det();
      if (pimpl->LLx.nlo() > 1) pimpl->det *= pimpl->det;
      pimpl->donedet = true;
    }
    return pimpl->det;
  }

  template <class T> template <class T1> void HermBandCHDiv<T>::DoInverse(
      const SymMatrixView<T1>& sinv) const
  {
    TMVAssert(sinv.size() == pimpl->LLx.size());
    TMVAssert(sinv.isherm());
    HermBandCH_Inverse(pimpl->LLx,sinv);
  }

  template <class T> template <class T1> void HermBandCHDiv<T>::DoInverse(
      const MatrixView<T1>& minv) const
  {
    TMVAssert(minv.colsize() == pimpl->LLx.size());
    TMVAssert(minv.rowsize() == pimpl->LLx.size());

    if (IsComplex(T1())) minv.diag().Imag().Zero();
    Inverse(HermMatrixViewOf(minv,Lower));
    if (minv.colsize() > 1)
      UpperTriMatrixViewOf(minv).OffDiag() =
	LowerTriMatrixViewOf(minv).OffDiag().Adjoint();
  }

  template <class T> void HermBandCHDiv<T>::DoInverseATA(
      const MatrixView<T>& ata) const
  {
    // ata = (At A)^-1 = A^-1 (A^-1)t
    //     = A^-1 A^-1
    SymMatrixView<T> hermata = HermMatrixViewOf(ata,Lower);
    DoInverse(hermata);
    SymATASquare<true>(ata);
  }

  template <class T> bool HermBandCHDiv<T>::Singular() const 
  { return Det() == T(0); }

  template <class T> const BandMatrix<T> HermBandCHDiv<T>::GetL() const 
  { 
    BandMatrix<T> L = pimpl->LLx.LowerBand();
    if (pimpl->LLx.nlo() <= 1) L.diag().SetAllTo(T(1));
    return L;
  }

  template <class T> const DiagMatrix<T> HermBandCHDiv<T>::GetD() const
  { 
    if (pimpl->LLx.nlo() <= 1) return DiagMatrix<T>(pimpl->LLx.diag());
    else return DiagMatrix<T>(pimpl->LLx.size(),T(1));
  }

  template <class T> const GenSymBandMatrix<T>& HermBandCHDiv<T>::GetLL() const 
  { return pimpl->LLx; }

  template <class T> std::string HermBandCHDiv<T>::Type() const
  { return std::string("HermBandCHDiv<") + tmv::Type(T()) + ">"; }

  template <class T> DivType HermBandCHDiv<T>::GetDivType() const 
  { return CH; }

  template <class T> bool HermBandCHDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, std::ostream* fout) const
  {
    Matrix<T> mm = m;
    if (fout) {
      *fout << "HermBandCHDiv:\n";
      *fout << "M = "<<mm<<std::endl;
      *fout << "L = "<<GetL()<<std::endl;
      *fout << "D = "<<GetD()<<std::endl;
    }
    BandMatrix<T> lu = GetL()*GetD()*GetL().Adjoint();
    RealType(T) nm = Norm(lu-mm);
    nm /= SQR(Norm(GetL())) * Norm(GetD());
    if (fout) {
      *fout << "LDLt = "<<lu<<std::endl;
      *fout << "M-LDLt = "<<(mm-lu)<<std::endl;
      *fout << "Norm(M-LDLt)/Norm(LDLt) = "<<nm<<std::endl;
    }
    mm.DivideUsing(SVS);
    mm.SetDiv();
    return nm < mm.Condition()*mm.colsize()*Epsilon<T>();
  }

  template <class T> size_t HermBandCHDiv<T>::colsize() const
  { return pimpl->LLx.size(); }

  template <class T> size_t HermBandCHDiv<T>::rowsize() const
  { return pimpl->LLx.size(); }

#define InstFile "TMV_SymBandCHDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


