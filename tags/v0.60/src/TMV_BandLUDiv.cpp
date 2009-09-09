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
#include "TMV_VectorArith.h"
#include "TMV_MatrixArith.h"
#include "TMV_TriMatrixArith.h"
#include "TMV_BandMatrixArith.h"
#include "TMV_BandLUDiv.h"
#include "TMV_DiagMatrix.h"
#include "TMV_TriMatrix.h"
#include <ostream>

namespace tmv {

  // MJ: Maybe I should add this AssignableTo / Gen flexibility to the 
  // other Div's, but this is the only one I really needed (for 
  // SymBandMatrix), so I haven't bothered to add it to the others.

  template <class T> struct BandLUDiv<T>::BandLUDiv_Impl
  {
    BandLUDiv_Impl(const GenBandMatrix<T>& A, bool _inplace);
    BandLUDiv_Impl(const AssignableToBandMatrix<T>& A);

    const bool istrans;
    const bool inplace;
    auto_array<T> Aptr1;
    T* Aptr;
    BandMatrixView<T> LUx;
    auto_array<size_t> P;
    mutable T det;
    mutable bool donedet;
  };

#define NEWLO std::min(A.nlo(),A.nhi())
#define NEWHI std::min(A.nlo()+A.nhi(),int(A.colsize())-1)
#define APTR1 inplace ? 0 : \
  new T[BandStorageLength(ColMajor,A.colsize(),A.colsize(),NEWLO,NEWHI)]
#define APTR inplace ? A.NonConst().ptr() : Aptr1.get()

#define LUX istrans ? \
  inplace ? \
  BandMatrixView<T>(A.NonConst().ptr(),A.colsize(),A.colsize(),A.nhi(), \
      NEWHI,A.stepj(),A.stepi(),A.diagstep(),TransOf(A.stor()),A.ct() \
      FIRSTLAST1(A.NonConst().first,A.NonConst().last) ) : \
  BandMatrixViewOf(Aptr,A.colsize(),A.colsize(),A.nhi(), NEWHI,  \
      (A.nlo() == 1 && A.nhi() == 1) ? DiagMajor : ColMajor) : \
  inplace ? \
  BandMatrixView<T>(A.NonConst().ptr(),A.colsize(),A.colsize(),A.nlo(), \
      NEWHI,A.stepi(),A.stepj(),A.diagstep(),A.stor(),A.ct() \
      FIRSTLAST1(A.NonConst().first,A.NonConst().last) ) : \
  BandMatrixViewOf(Aptr,A.colsize(),A.colsize(),A.nlo(), NEWHI, \
      (A.nlo() == 1 && A.nhi() == 1) ? DiagMajor : ColMajor)

  template <class T> BandLUDiv<T>::BandLUDiv_Impl::BandLUDiv_Impl(
      const GenBandMatrix<T>& A, bool _inplace) :
    istrans(A.nhi()<A.nlo() || (A.nhi()==A.nlo() && A.isrm())),
    inplace(NEWLO == 0 || (_inplace && 
	  ((A.isrm() && istrans) || (A.iscm() && !istrans) || 
	   (A.isdm() && A.nlo()==1 && A.nhi()==1)))),
    Aptr1(APTR1), Aptr(APTR), LUx(LUX),
    P(new size_t[A.colsize()]), det(T(1)), donedet(false) {}

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
	TMVAssert(A.stepj() >= std::min(int(A.colsize()),2*A.nlo()+A.nhi()));
	TMVAssert(pimpl->LUx.Diags(-A.nlo(),A.nhi()+1) == A);
      } else if (A.isrm()) {
	TMVAssert(pimpl->istrans);
	TMVAssert(A.stepi() >= std::min(int(A.colsize()),2*A.nhi()+A.nlo()));
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
      BandLU_Decompose(pimpl->LUx,pimpl->P.get(),pimpl->det,Anhi);
    } else {
      pimpl->P.get()[0] = 1; // A tag to indicate that P is not set yet.
    }
  }

#define NEWLO std::min(A.nlo(),A.nhi())
#define NEWHI std::min(A.nlo()+A.nhi(),int(A.colsize())-1)
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
    P(new size_t[A.colsize()]), det(T(1)), donedet(false) {}

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
      BandLU_Decompose(pimpl->LUx,pimpl->P.get(),pimpl->det,Anhi);
    } else {
      pimpl->P.get()[0] = 1; // A tag to indicate that P is not set yet.
    }
  }

  template <class T> BandLUDiv<T>::~BandLUDiv() { delete pimpl; }

#define CT std::complex<T>
  template <class T> inline void BandLU_LDivEq(const GenBandMatrix<CT>& ,
      const size_t* , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void BandLU_RDivEq(const GenBandMatrix<CT>& ,
      const size_t* , const MatrixView<T>& )
  { TMVAssert(FALSE); }
#undef CT

  template <class T> template <class T1> void BandLUDiv<T>::DoLDivEq(
      const MatrixView<T1>& m) const
  {
    if (pimpl->istrans) BandLU_RDivEq(pimpl->LUx,pimpl->P.get(),m.Transpose());
    else BandLU_LDivEq(pimpl->LUx,pimpl->P.get(),m);
  }

  template <class T> template <class T1> void BandLUDiv<T>::DoRDivEq(
      const MatrixView<T1>& m) const
  {
    if (pimpl->istrans) BandLU_LDivEq(pimpl->LUx,pimpl->P.get(),m.Transpose());
    else BandLU_RDivEq(pimpl->LUx,pimpl->P.get(),m);
  }

  template <class T> template <class T1, class T2> void BandLUDiv<T>::DoLDiv(
      const GenMatrix<T1>& m, const MatrixView<T2>& x) const
  {
    if (pimpl->istrans) 
      BandLU_RDivEq(pimpl->LUx,pimpl->P.get(),(x=m).Transpose());
    else BandLU_LDivEq(pimpl->LUx,pimpl->P.get(),x=m);
  }

  template <class T> template <class T1, class T2> void BandLUDiv<T>::DoRDiv(
      const GenMatrix<T1>& m, const MatrixView<T2>& x) const
  {
    if (pimpl->istrans) 
      BandLU_LDivEq(pimpl->LUx,pimpl->P.get(),(x=m).Transpose());
    else BandLU_RDivEq(pimpl->LUx,pimpl->P.get(),x=m);
  }

  template <class T> T BandLUDiv<T>::Det() const
  {
    if (!pimpl->donedet) {
      pimpl->det *= DiagMatrixViewOf(pimpl->LUx.diag()).Det();
      pimpl->donedet = true;
    }         
    return pimpl->det;  
  }                  

  template <class T> template <class T1> void BandLUDiv<T>::DoInverse(
      const MatrixView<T1>& minv) const
  {
    minv.SetToIdentity();
    LDivEq(minv);
  }

  template <class T> void BandLUDiv<T>::DoInverseATA(
      const MatrixView<T>& minv) const
  {
    Matrix<T,ColMajor> temp(minv.colsize(),minv.rowsize());
    Inverse(temp.View());
    minv = temp*Transpose(temp);
  }

  template <class T> bool BandLUDiv<T>::Singular() const 
  { return Det() == T(0); }

  template <class T> bool BandLUDiv<T>::IsTrans() const 
  { return pimpl->istrans; }

  template <class T> ConstBandMatrixView<T> BandLUDiv<T>::GetU() const
  { return BandMatrixViewOf(pimpl->LUx,0,pimpl->LUx.nhi()); }

  template <class T> LowerTriMatrix<T,UnitDiag> BandLUDiv<T>::GetL() const
  {
    size_t N = pimpl->LUx.colsize();
    int nlo = pimpl->LUx.nlo();
    LowerTriMatrix<T,UnitDiag> L(N,T(0));
    if (nlo == 0) 
      L.SetToIdentity();
    else {
      for(size_t i=0;i<N;++i) {
	Swap(L.row(i,0,i),L.row(pimpl->P.get()[i],0,i));
	size_t end = std::min(i+nlo+1,N);
	L.col(i,i+1,end) = pimpl->LUx.col(i,i+1,end);
      }
    }
    return L;
  }

  template <class T> const GenBandMatrix<T>& BandLUDiv<T>::GetLU() const 
  { return pimpl->LUx; }

  template <class T> const size_t* BandLUDiv<T>::GetP() const
  {
    if (pimpl->LUx.nlo() == 0 && pimpl->LUx.colsize() > 0 && 
	pimpl->P.get()[0] == 1) {
      for(size_t i=0;i<pimpl->LUx.colsize();++i) pimpl->P.get()[i] = i;
    }
    return pimpl->P.get();
  }

  template <class T> std::string BandLUDiv<T>::Type() const
  { return std::string("BandLUDiv<") + tmv::Type(T()) + ">"; }

  template <class T> DivType BandLUDiv<T>::GetDivType() const 
  { return LU; }

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
    mm.DivideUsing(SVS);
    mm.SetDiv();
    return nm < mm.Condition()*mm.colsize()*Epsilon<T>();
  }

  template <class T> size_t BandLUDiv<T>::colsize() const
  { return pimpl->LUx.colsize(); }

  template <class T> size_t BandLUDiv<T>::rowsize() const
  { return pimpl->LUx.rowsize(); }

#define InstFile "TMV_BandLUDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


