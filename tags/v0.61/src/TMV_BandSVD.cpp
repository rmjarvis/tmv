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



#include "TMV_BandSVD.h"
#include "TMV_BandSVDiv.h"
#include "TMV_BandMatrix.h"
#include "TMV_Matrix.h"
#include "TMV_Vector.h"
#include "TMV_DiagMatrix.h"
#include "TMV_MatrixArith.h"
#include "TMV_DiagMatrixArith.h"
#include <ostream>

namespace tmv {

#define RT RealType(T)

  template <class T> struct BandSVDiv<T>::BandSVDiv_Impl
  {
    public:
      BandSVDiv_Impl(const GenBandMatrix<T>& A);

      const bool istrans;
      Matrix<T,ColMajor> U;
      DiagMatrix<RT> S;
      Matrix<T,ColMajor> V;
      RT logdet;
      T signdet;
      mutable int kmax;

    private :
      BandSVDiv_Impl(const BandSVDiv_Impl&) :
	istrans(false), U(0,0), S(0), V(0,0), logdet(0), signdet(0), kmax(0)
      { TMVAssert(FALSE); }
      BandSVDiv_Impl& operator=(const BandSVDiv_Impl&)
      { TMVAssert(FALSE); return *this; }
  }; 

#define M MAX(A.colsize(),A.rowsize())
#define N MIN(A.colsize(),A.rowsize())

  template <class T> BandSVDiv<T>::BandSVDiv_Impl::BandSVDiv_Impl(
      const GenBandMatrix<T>& A) :
    istrans(A.colsize() < A.rowsize()),
    U(M,N), S(N), V(N,N), logdet(0), signdet(1), kmax(0) {}

  template <class T> BandSVDiv<T>::BandSVDiv(const GenBandMatrix<T>& A) :
    pimpl(new BandSVDiv_Impl(A))
  {
    SV_Decompose<T>(pimpl->istrans ? A.Transpose() : A.View(),
	pimpl->U.View(),pimpl->S.View(),pimpl->V.View(),
	pimpl->logdet,pimpl->signdet);
    Thresh(Epsilon<T>());
  }

  template <class T> BandSVDiv<T>::~BandSVDiv() { delete pimpl; pimpl=0; }

  template <class T> template <class T1> void BandSVDiv<T>::DoLDivEq(
      const MatrixView<T1>& m) const
  {
    TMVAssert(m.colsize() == rowsize());
    TMVAssert(m.colsize() == colsize());
    if (pimpl->istrans) SV_RDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,
	m.Transpose(),m.Transpose());
    else SV_LDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,m,m);
  }

  template <class T> template <class T1> void BandSVDiv<T>::DoRDivEq(
      const MatrixView<T1>& m) const
  {
    TMVAssert(m.rowsize() == rowsize());
    TMVAssert(m.rowsize() == colsize());
    if (pimpl->istrans) SV_LDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,
	m.Transpose(),m.Transpose());
    else SV_RDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,m,m);
  }

  template <class T> template <class T1, class T2> void BandSVDiv<T>::DoLDiv(
      const GenMatrix<T1>& m, const MatrixView<T2>& x) const
  {
    TMVAssert(m.rowsize() == x.rowsize());
    TMVAssert(m.colsize() == colsize());
    TMVAssert(x.colsize() == rowsize());
    if (pimpl->istrans) SV_RDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,
	m.Transpose(),x.Transpose());
    else SV_LDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,m,x);
  }

  template <class T> template <class T1, class T2> void BandSVDiv<T>::DoRDiv(
      const GenMatrix<T1>& m, const MatrixView<T2>& x) const
  {
    TMVAssert(m.colsize() == x.colsize());
    TMVAssert(m.rowsize() == rowsize());
    TMVAssert(x.rowsize() == colsize());
    if (pimpl->istrans) SV_LDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,
	m.Transpose(),x.Transpose());
    else SV_RDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,m,x);
  }

  template <class T> T BandSVDiv<T>::Det() const 
  { 
    if (pimpl->signdet == T(0)) return T(0);
    else return pimpl->signdet * std::exp(pimpl->logdet);  
  }

  template <class T> RT BandSVDiv<T>::LogDet(T* sign) const 
  { 
    if (sign) *sign = pimpl->signdet;
    return pimpl->logdet;
  }

  template <class T> template <class T1> void BandSVDiv<T>::DoInverse(
      const MatrixView<T1>& minv) const
  { 
    if (pimpl->istrans) {
      Matrix<T,ColMajor> SinvV = pimpl->V.Conjugate().Rows(0,pimpl->kmax) /
	pimpl->S.SubDiagMatrix(0,pimpl->kmax);
      minv = pimpl->U.Conjugate().Cols(0,pimpl->kmax) * SinvV;
    } else  {
      Matrix<T,ColMajor> SinvUt = pimpl->U.Adjoint().Rows(0,pimpl->kmax) /
	pimpl->S.SubDiagMatrix(0,pimpl->kmax);
      minv = pimpl->V.Adjoint().Cols(0,pimpl->kmax) * SinvUt;
    }
  }

  template <class T> void BandSVDiv<T>::DoInverseATA(
      const MatrixView<T>& minv) const
  {
    Matrix<T,ColMajor> SinvV = pimpl->V.Rows(0,pimpl->kmax) /
      pimpl->S.SubDiagMatrix(0,pimpl->kmax);
    if (pimpl->istrans)
      minv = SinvV.Transpose() * SinvV.Conjugate();
    else
      minv = SinvV.Adjoint() * SinvV;
  }

  template <class T> bool BandSVDiv<T>::Singular() const 
  { return pimpl->kmax < int(pimpl->S.size()); }

  template <class T> RT BandSVDiv<T>::Norm2() const 
  { return pimpl->S.size() > 0 ? pimpl->S(0) : RT(0); }

  template <class T> RT BandSVDiv<T>::Condition() const 
  {
    return pimpl->S.size() > 0 ? 
      pimpl->S(0)/pimpl->S(pimpl->S.size()-1) : RT(1); 
  }

  template <class T> void BandSVDiv<T>::Thresh(RT toler,
      std::ostream* debugout) const
  {
    if (pimpl->S.size() == 0) pimpl->kmax = 0;
    else {
      TMVAssert(toler < RT(1));
      RT thresh = pimpl->S(0)*toler;
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
  }

  template <class T> void BandSVDiv<T>::Top(int neigen,
      std::ostream* debugout) const
  {
    TMVAssert(neigen > 0);
    if (neigen < int(pimpl->S.size())) pimpl->kmax = neigen;
    else pimpl->kmax = pimpl->S.size();
    if(debugout) {
      (*debugout)<<"S = "<<pimpl->S<<std::endl;
      (*debugout)<<"kmax = "<<pimpl->kmax;
      (*debugout)<<" (S.size = "<<pimpl->S.size()<<")"<<std::endl;
    }
  }

  template <class T> int BandSVDiv<T>::GetKMax() const 
  { return pimpl->kmax; }

  template <class T> ConstMatrixView<T> BandSVDiv<T>::GetU() const
  {
    if (pimpl->istrans) {
      return pimpl->V.Transpose(); 
    }
    else { 
      return pimpl->U.View(); 
    }
  }

  template <class T> ConstDiagMatrixView<RT> BandSVDiv<T>::GetS() const
  { return pimpl->S.View(); }

  template <class T> ConstMatrixView<T> BandSVDiv<T>::GetV() const
  {
    if (pimpl->istrans) {
      return pimpl->U.Transpose(); 
    }
    else { 
      return pimpl->V.View();
    }
  }

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
    Matrix<T> usv = GetU()*GetS()*GetV();
    RT nm = Norm(usv-mm);
    nm /= Norm(GetU())*Norm(GetS())*Norm(GetV());
    RT cond = GetS()(0) / GetS()(GetKMax()-1);
    if (fout) {
      *fout << "USV = "<<usv<<std::endl;
      *fout << "Norm(M-USV) = "<<Norm(mm-usv)<<std::endl;
      *fout << "Norm(M-USV)/Norm(USV) = "<<nm<<std::endl;
    }
    return nm < cond*mm.colsize()*Epsilon<T>();
  }

  template <class T> size_t BandSVDiv<T>::colsize() const
  { return pimpl->istrans ? pimpl->U.rowsize() : pimpl->U.colsize(); }

  template <class T> size_t BandSVDiv<T>::rowsize() const
  { return pimpl->istrans ? pimpl->U.colsize() : pimpl->U.rowsize(); }

#undef RT

#define InstFile "TMV_BandSVD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


