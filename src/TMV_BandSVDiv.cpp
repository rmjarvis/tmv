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
#include "TMV_DiagMatrix.h"
#include "TMV_DiagMatrixArith.h"
#include "TMV_BandMatrixArith.h"
#include "TMV_SVDiv.h"
#include "TMV_BandSVDiv.h"
#include <algorithm>
#include <ostream>

namespace tmv {

  template <class T> struct BandSVDiv<T>::BandSVDiv_Impl
  {
    BandSVDiv_Impl(const GenBandMatrix<T>& A);

    const bool istrans;
    auto_ptr<Matrix<T,ColMajor> > U;
    Vector<RealType(T)> S;
    auto_ptr<Matrix<T,ColMajor> > V;
    T det;
    mutable size_t kmax;
  }; 

  template <class T> BandSVDiv<T>::BandSVDiv_Impl::BandSVDiv_Impl(
      const GenBandMatrix<T>& A) :
    istrans(A.colsize() < A.rowsize()),
    U(0), S(std::min(A.rowsize(),A.colsize())), V(0), det(T(1)) {}

  template <class T> BandSVDiv<T>::BandSVDiv(const GenBandMatrix<T>& A,
      bool StoreU, bool StoreV) : pimpl(new BandSVDiv_Impl(A))
  {
    size_t M = std::max(A.colsize(),A.rowsize());
    size_t N = std::min(A.colsize(),A.rowsize());

    if (pimpl->istrans) std::swap(StoreU,StoreV);

    auto_ptr<MatrixView<T> > Uv(0);
    auto_ptr<MatrixView<T> > Vv(0);
    if (StoreU) {
      pimpl->U.reset(new Matrix<T,ColMajor>(M,N));
      Uv.reset(new MatrixView<T>(pimpl->U->View()));
    }
    if (StoreV) {
      pimpl->V.reset(new Matrix<T,ColMajor>(N,N));
      Vv.reset(new MatrixView<T>(pimpl->V->View()));
    }

    if (pimpl->istrans)
      BandSV_Decompose(A.Transpose(),Uv.get(),pimpl->S.View(),Vv.get(),
	  pimpl->det);
    else
      BandSV_Decompose(A,Uv.get(),pimpl->S.View(),Vv.get(),pimpl->det);

    Thresh(Epsilon<T>());
  }

  template <class T> BandSVDiv<T>::~BandSVDiv() { delete pimpl; }

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

  template <class T> template <class T1> void BandSVDiv<T>::DoLDivEq2(
      const MatrixView<T1>& m) const
  {
    TMVAssert(pimpl->U.get() && pimpl->V.get());
    TMVAssert(m.colsize() == rowsize());
    TMVAssert(m.colsize() == colsize());
    if (pimpl->istrans) SV_RDiv(*pimpl->U,pimpl->S,*pimpl->V,pimpl->kmax,
	m.Transpose(),m.Transpose());
    else SV_LDiv(*pimpl->U,pimpl->S,*pimpl->V,pimpl->kmax,m,m);
  }

  template <class T> template <class T1> void BandSVDiv<T>::DoRDivEq2(
      const MatrixView<T1>& m) const
  {
    TMVAssert(pimpl->U.get() && pimpl->V.get());
    TMVAssert(m.rowsize() == rowsize());
    TMVAssert(m.rowsize() == colsize());
    if (pimpl->istrans) SV_LDiv(*pimpl->U,pimpl->S,*pimpl->V,pimpl->kmax,
	m.Transpose(),m.Transpose());
    else SV_RDiv(*pimpl->U,pimpl->S,*pimpl->V,pimpl->kmax,m,m);
  }

  template <class T> template <class T1, class T2> void BandSVDiv<T>::DoLDiv2(
      const GenMatrix<T1>& m, const MatrixView<T2>& x) const
  {
    TMVAssert(pimpl->U.get() && pimpl->V.get());
    TMVAssert(m.rowsize() == x.rowsize());
    TMVAssert(m.colsize() == colsize());
    TMVAssert(x.colsize() == rowsize());
    if (pimpl->istrans) SV_RDiv(*pimpl->U,pimpl->S,*pimpl->V,pimpl->kmax,
	m.Transpose(),x.Transpose());
    else SV_LDiv(*pimpl->U,pimpl->S,*pimpl->V,pimpl->kmax,m,x);
  }

  template <class T> template <class T1, class T2> void BandSVDiv<T>::DoRDiv2(
      const GenMatrix<T1>& m, const MatrixView<T2>& x) const
  {
    TMVAssert(pimpl->U.get() && pimpl->V.get());
    TMVAssert(m.colsize() == x.colsize());
    TMVAssert(m.rowsize() == rowsize());
    TMVAssert(x.rowsize() == colsize());
    if (pimpl->istrans) SV_LDiv(*pimpl->U,pimpl->S,*pimpl->V,pimpl->kmax,
	m.Transpose(),x.Transpose());
    else SV_RDiv(*pimpl->U,pimpl->S,*pimpl->V,pimpl->kmax,m,x);
  }

  template <class T> T BandSVDiv<T>::Det() const 
  { return pimpl->det; }

  template <class T> template <class T1> void BandSVDiv<T>::DoInverse2(
      const MatrixView<T1>& minv) const
  { 
    TMVAssert(pimpl->U.get() && pimpl->V.get());
    if (pimpl->istrans) {
      Matrix<T,ColMajor> SinvV = pimpl->V->Conjugate().Rows(0,pimpl->kmax) /
	DiagMatrixViewOf(pimpl->S.SubVector(0,pimpl->kmax));
      minv = pimpl->U->Conjugate().Cols(0,pimpl->kmax) * SinvV;
    } else  {
      Matrix<T,ColMajor> SinvUt = pimpl->U->Adjoint().Rows(0,pimpl->kmax) /
	DiagMatrixViewOf(pimpl->S.SubVector(0,pimpl->kmax));
      minv = pimpl->V->Adjoint().Cols(0,pimpl->kmax) * SinvUt;
    }
  }

  template <class T> void BandSVDiv<T>::DoInverseATA2(
      const MatrixView<T>& minv) const
  {
    TMVAssert(pimpl->V.get());
    Matrix<T,ColMajor> SinvV = pimpl->V->Rows(0,pimpl->kmax) /
      DiagMatrixViewOf(pimpl->S.SubVector(0,pimpl->kmax));
    if (pimpl->istrans)
      minv = SinvV.Transpose() * SinvV.Conjugate();
    else
      minv = SinvV.Adjoint() * SinvV;
  }

  template <class T> bool BandSVDiv<T>::Singular() const 
  { return pimpl->kmax < pimpl->S.size(); }

  template <class T> RealType(T) BandSVDiv<T>::Norm2() const 
  { return pimpl->S(0); }

  template <class T> RealType(T) BandSVDiv<T>::Condition() const 
  { return pimpl->S(0)/pimpl->S(pimpl->S.size()-1); }

  template <class T> void BandSVDiv<T>::Thresh(RealType(T) toler,
      std::ostream* debugout) const
  {
    TMVAssert(toler < RealType(T)(1));
    RealType(T) thresh = pimpl->S(0)*toler;
    for(pimpl->kmax=pimpl->S.size(); 
	pimpl->kmax>0 && pimpl->S(pimpl->kmax-1)<=thresh; 
	--pimpl->kmax);
    if(debugout) {
      (*debugout)<<"S = "<<pimpl->S<<std::endl;
      (*debugout)<<"Smax = "<<pimpl->S(0)<<", thresh = "<<thresh<<std::endl;
      (*debugout)<<"kmax = "<<pimpl->kmax;
      (*debugout)<<" (S.size = "<<pimpl->S.size()<<")"<<std::endl;
    }
  }

  template <class T> void BandSVDiv<T>::Top(size_t neigen,
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

  template <class T> size_t BandSVDiv<T>::GetKMax() const 
  { return pimpl->kmax; }

  template <class T> ConstMatrixView<T> BandSVDiv<T>::GetU() const
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

  template <class T> ConstVectorView<RealType(T)> BandSVDiv<T>::GetS() const
  { return pimpl->S.View(); }

  template <class T> ConstMatrixView<T> BandSVDiv<T>::GetV() const
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

  template <class T> std::string BandSVDiv<T>::Type() const
  { return std::string("BandSVDiv<") + tmv::Type(T()) + ">"; }

  template <class T> DivType BandSVDiv<T>::GetDivType() const 
  { return SV; }

  template <class T> bool BandSVDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, std::ostream* fout) const
  {
    Matrix<T> mm = m;
    if (fout) {
      *fout << "BandSVDiv:\n";
      *fout << "M = "<<m<<std::endl;
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
      *fout << "Norm(M-USV) = "<<Norm(nm-usv)<<std::endl;
      *fout << "Norm(M-USV)/Norm(USV) = "<<nm<<std::endl;
    }
    return nm < cond*mm.colsize()*Epsilon<T>();
  }

  template <class T> size_t BandSVDiv<T>::colsize() const
  { return pimpl->istrans ? pimpl->U->rowsize() : pimpl->U->colsize(); }

  template <class T> size_t BandSVDiv<T>::rowsize() const
  { return pimpl->istrans ? pimpl->U->colsize() : pimpl->U->rowsize(); }

  template <class T> bool BandSVDiv<T>::Uisset() const
  { return pimpl->U.get(); }

  template <class T> bool BandSVDiv<T>::Visset() const
  { return pimpl->V.get(); }

#define InstFile "TMV_BandSVDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


