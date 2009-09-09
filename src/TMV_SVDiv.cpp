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
#include "TMV_DiagMatrix.h"
#include "TMV_SVDiv.h"
#include "TMV_MatrixArith.h"
#include "TMV_DiagMatrixArith.h"
#include <algorithm>
#include <ostream>

namespace tmv {

  template <class T> struct SVDiv<T>::SVDiv_Impl
  {
    SVDiv_Impl(const GenMatrix<T>& m, bool inplace); 

    const bool istrans;
    const bool inplace;
    auto_array<T> Aptr1;
    T* Aptr;
    auto_ptr<MatrixView<T> > U;
    Vector<RealType(T)> S;
    auto_ptr<Matrix<T,ColMajor> > V;
    T det;
    mutable size_t kmax;
  };

#define APTR1 inplace ? 0 : new T[A.colsize()*A.rowsize()]     
#define APTR inplace ? A.NonConst().ptr() : Aptr1.get()
#define UX istrans ? \
  inplace ? A.NonConst().Transpose() : \
  MatrixViewOf(Aptr,A.rowsize(),A.colsize(),ColMajor) : \
  inplace ? A.NonConst().View() : \
  MatrixViewOf(Aptr,A.colsize(),A.rowsize(), \
      A.IsSquare() ? BaseStorOf(A) : ColMajor)

  template <class T> SVDiv<T>::SVDiv_Impl::SVDiv_Impl(
      const GenMatrix<T>& A, bool _inplace) :
    istrans(A.colsize() < A.rowsize()), 
    inplace(_inplace && (A.isrm() || A.iscm())), Aptr1(APTR1), Aptr(APTR),
    U(new MatrixView<T>(UX)), S(U->rowsize()), V(0), det(T(1)) {}

#undef UX
#undef APTR

  template <class T> SVDiv<T>::SVDiv(const GenMatrix<T>& A,
      bool inplace, bool StoreU, bool StoreV) :
    pimpl(new SVDiv_Impl(A,inplace))
  {
    if (pimpl->istrans) {
      if (inplace) TMVAssert(A.Transpose() == *pimpl->U); 
      else *pimpl->U = A.Transpose();
    } else {
      if (inplace) TMVAssert(A == *pimpl->U); 
      else *pimpl->U = A;
    }

    if (pimpl->istrans) std::swap(StoreU,StoreV);

    if (StoreV) {
      pimpl->V.reset(
	  new Matrix<T,ColMajor>(pimpl->U->rowsize(),pimpl->U->rowsize()));
      SV_Decompose(*pimpl->U,pimpl->S.View(),pimpl->V->View(),
	  pimpl->det,StoreU);
    }
    else SV_Decompose(*pimpl->U,pimpl->S.View(),pimpl->det,StoreU);

    if (!StoreU && !inplace) { pimpl->Aptr1.reset(0); }
    if (!StoreU) { pimpl->U.reset(0); }

    // Set kmax for actual 0 elements (to within machine precision).
    // Any further cut in the number of singular values to use
    // should be done by the user.
    Thresh(Epsilon<T>());
  }

  template <class T> SVDiv<T>::~SVDiv() { delete pimpl; }

#define CT std::complex<T>
  template <class T> inline void SV_LDiv(
      const GenMatrix<CT>& , const GenVector<T>& ,
      const GenMatrix<CT>& , size_t ,
      const GenMatrix<T>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void SV_LDiv(
      const GenMatrix<CT>& , const GenVector<T>& ,
      const GenMatrix<CT>& , size_t ,
      const GenMatrix<CT>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void SV_LDiv(
      const GenMatrix<T>& , const GenVector<T>& ,
      const GenMatrix<T>& , size_t ,
      const GenMatrix<CT>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void SV_RDiv(
      const GenMatrix<CT>& , const GenVector<T>& ,
      const GenMatrix<CT>& , size_t ,
      const GenMatrix<T>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void SV_RDiv(
      const GenMatrix<CT>& , const GenVector<T>& ,
      const GenMatrix<CT>& , size_t ,
      const GenMatrix<CT>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
  template <class T> inline void SV_RDiv(
      const GenMatrix<T>& , const GenVector<T>& ,
      const GenMatrix<T>& , size_t ,
      const GenMatrix<CT>& , const MatrixView<T>& )
  { TMVAssert(FALSE); }
#undef CT

  template <class T> template <class T1> void SVDiv<T>::DoLDivEq2(
      const MatrixView<T1>& m) const
  {
    TMVAssert(pimpl->U.get() && pimpl->V.get());
    TMVAssert(m.colsize() == colsize());
    TMVAssert(m.colsize() == rowsize());
    if (pimpl->istrans) 
      SV_RDiv(*pimpl->U,pimpl->S,*pimpl->V,pimpl->kmax,
	  m.Transpose(),m.Transpose());
    else 
      SV_LDiv(*pimpl->U,pimpl->S,*pimpl->V,pimpl->kmax,m,m);
  }

  template <class T> template <class T1> void SVDiv<T>::DoRDivEq2(
      const MatrixView<T1>& m) const
  {
    TMVAssert(pimpl->U.get() && pimpl->V.get());
    TMVAssert(m.rowsize() == colsize());
    TMVAssert(m.rowsize() == rowsize());

    if (pimpl->istrans) 
      SV_LDiv(*pimpl->U,pimpl->S,*pimpl->V,pimpl->kmax,
	  m.Transpose(),m.Transpose());
    else 
      SV_RDiv(*pimpl->U,pimpl->S,*pimpl->V,pimpl->kmax,m,m);
  }

  template <class T> template <class T1, class T2> void SVDiv<T>::DoLDiv2(
      const GenMatrix<T1>& m, const MatrixView<T2>& x) const
  {
    TMVAssert(pimpl->U.get() && pimpl->V.get());
    TMVAssert(m.rowsize() == x.rowsize());
    TMVAssert(m.colsize() == colsize());
    TMVAssert(x.colsize() == rowsize());
    if (pimpl->istrans) 
      SV_RDiv(*pimpl->U,pimpl->S,*pimpl->V,pimpl->kmax,
	m.Transpose(),x.Transpose());
    else 
      SV_LDiv(*pimpl->U,pimpl->S,*pimpl->V,pimpl->kmax,m,x);
  }

  template <class T> template <class T1, class T2> void SVDiv<T>::DoRDiv2(
      const GenMatrix<T1>& m, const MatrixView<T2>& x) const
  {
    TMVAssert(pimpl->U.get() && pimpl->V.get());
    TMVAssert(m.colsize() == x.colsize());
    TMVAssert(m.rowsize() == rowsize());
    TMVAssert(x.rowsize() == colsize());

    if (pimpl->istrans) 
      SV_LDiv(*pimpl->U,pimpl->S,*pimpl->V,pimpl->kmax,
	  m.Transpose(),x.Transpose());
    else 
      SV_RDiv(*pimpl->U,pimpl->S,*pimpl->V,pimpl->kmax,m,x);
  }

  template <class T> T SVDiv<T>::Det() const 
  { return pimpl->det; }

  template <class T> template <class T1> void SVDiv<T>::DoInverse2(
      const MatrixView<T1>& minv) const
  { 
    TMVAssert(pimpl->U.get() && pimpl->V.get());
    if (pimpl->istrans) {
      // A^-1 = (Vt S^-1 Ut)T = U* S^-1 V*
      Matrix<T,ColMajor> SinvV = pimpl->V->Conjugate().Rows(0,pimpl->kmax) /
	DiagMatrixViewOf(pimpl->S.SubVector(0,pimpl->kmax));
      minv = pimpl->U->Conjugate().Cols(0,pimpl->kmax) * SinvV;
    } else {
      // A^-1 = Vt S^-1 Ut
      Matrix<T,ColMajor> SinvUt = pimpl->U->Adjoint().Rows(0,pimpl->kmax) /
	DiagMatrixViewOf(pimpl->S.SubVector(0,pimpl->kmax));
      minv = pimpl->V->Adjoint().Cols(0,pimpl->kmax) * SinvUt;
    }
  }

  template <class T> void SVDiv<T>::DoInverseATA2(
      const MatrixView<T>& minv) const
  {
    TMVAssert(pimpl->V.get());
    // A = U S V
    // At = Vt S Ut
    // AtA = Vt S^2 V
    // (AtA)^-1 = Vt S^-2 V
    //
    // if istrans:
    // AT = U S V
    // At = U* S V*
    // AAt = VT S^2 V*
    // (AAt)^-1 = VT S^-2 V*
    //
    Matrix<T,ColMajor> SinvV = pimpl->V->Rows(0,pimpl->kmax) /
      DiagMatrixViewOf(pimpl->S.SubVector(0,pimpl->kmax));
    if (pimpl->istrans)
      minv = SinvV.Transpose() * SinvV.Conjugate();
    else
      minv = SinvV.Adjoint() * SinvV;
  }

  template <class T> bool SVDiv<T>::Singular() const
  { return pimpl->kmax < pimpl->S.size(); }
  
  template <class T> RealType(T) SVDiv<T>::Norm2() const 
  { return pimpl->S(0); }

  template <class T> RealType(T) SVDiv<T>::Condition() const 
  { return pimpl->S(0)/pimpl->S(pimpl->S.size()-1); }

  template <class T> void SVDiv<T>::Thresh(RealType(T) toler,
      std::ostream* debugout) const
  {
    TMVAssert(toler < RealType(T)(1));
    RealType(T) thresh = pimpl->S(0)*toler;
    for(pimpl->kmax=pimpl->S.size(); 
	pimpl->kmax>0 && pimpl->S(pimpl->kmax-1)<=thresh; --pimpl->kmax);
    if(debugout) {
      (*debugout)<<"S = "<<pimpl->S<<std::endl;
      (*debugout)<<"Smax = "<<pimpl->S(0)<<", thresh = "<<thresh<<std::endl;
      (*debugout)<<"kmax = "<<pimpl->kmax;
      (*debugout)<<" (S.size = "<<pimpl->S.size()<<")"<<std::endl;
    }
  }

  template <class T> void SVDiv<T>::Top(size_t neigen,
      std::ostream* debugout) const
  {
    TMVAssert(neigen > 0);
    TMVAssert(neigen <= pimpl->S.size());
    pimpl->kmax = neigen;
    if(debugout) {
      (*debugout)<<"S = "<<pimpl->S<<std::endl;
      (*debugout)<<"kmax = "<<pimpl->kmax;
      (*debugout)<<" (S.size = "<<pimpl->S.size()<<")"<<std::endl;
    }
  }

  template <class T> size_t SVDiv<T>::GetKMax() const
  { return pimpl->kmax; }

  template <class T> ConstMatrixView<T> SVDiv<T>::GetU() const
  {
    if (pimpl->istrans) { 
      TMVAssert(pimpl->V.get());
      return pimpl->V->Transpose(); 
    }
    else {
      TMVAssert(pimpl->U.get());
      return pimpl->U->View(); 
    }
  }
  template <class T> ConstVectorView<RealType(T)> SVDiv<T>::GetS() const
  { return pimpl->S.View(); }

  template <class T> ConstMatrixView<T> SVDiv<T>::GetV() const
  {
    if (pimpl->istrans) {
      TMVAssert(pimpl->U.get());
      return pimpl->U->Transpose(); 
    }
    else {
      TMVAssert(pimpl->V.get());
      return pimpl->V->View(); 
    }
  }

  template <class T> std::string SVDiv<T>::Type() const
  { return std::string("SVDiv<") + tmv::Type(T()) + ">"; }

  template <class T> DivType SVDiv<T>::GetDivType() const
  { 
    return pimpl->U.get() ? pimpl->V.get() ? SV : SVU : 
      pimpl->V.get() ? SVV : SVS; 
  }

  template <class T> bool SVDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, std::ostream* fout) const
  {
    Matrix<T> mm = m;
    if (fout) {
      *fout << "SVDiv:\n";
      *fout << "M = "<<mm<<std::endl;
      *fout << "U = "<<GetU()<<std::endl;
      *fout << "S = "<<GetS()<<std::endl;
      *fout << "V = "<<GetV()<<std::endl;
    }
    Matrix<T> usv = GetU()*DiagMatrixViewOf(GetS())*GetV();
    RealType(T) nm = Norm(usv-mm);
    nm /= Norm(GetU())*Norm(GetS())*Norm(GetV());
    RealType(T) cond = GetS()(0) / GetS()(GetKMax()-1);
    if (fout) {
      *fout << "USV = "<<usv<<std::endl;
      *fout << "Norm(M-USV)/Norm(USV) = "<<nm;
      *fout <<"  "<<cond<<" * "<<Epsilon<T>()<<std::endl;
    }
    return nm < cond*Epsilon<T>();
  }

  template <class T> size_t SVDiv<T>::colsize() const
  { return pimpl->istrans ? pimpl->U->rowsize() : pimpl->U->colsize(); }

  template <class T> size_t SVDiv<T>::rowsize() const
  { return pimpl->istrans ? pimpl->U->colsize() : pimpl->U->rowsize(); }

  template <class T> bool SVDiv<T>::Uisset() const
  { return pimpl->U.get(); }

  template <class T> bool SVDiv<T>::Visset() const
  { return pimpl->V.get(); }

#define InstFile "TMV_SVDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


