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



#include "TMV_SVD.h"
#include "TMV_SVDiv.h"
#include "TMV_Matrix.h"
#include "TMV_DiagMatrix.h"
#include "TMV_DiagMatrixArith.h"
#include "TMV_MatrixArith.h"
#include <ostream>

namespace tmv {

#define RT RealType(T)

  template <class T> struct SVDiv<T>::SVDiv_Impl
  {
    public :
      SVDiv_Impl(const GenMatrix<T>& m, bool inplace); 

      const bool istrans;
      const bool inplace;
      auto_array<T> Aptr1;
      T* Aptr;
      MatrixView<T> U;
      DiagMatrix<RT> S;
      Matrix<T,ColMajor> V;
      RT logdet;
      T signdet;
      mutable int kmax;

    private :
      SVDiv_Impl(const SVDiv_Impl&) :
	istrans(false), inplace(false), Aptr1(0), Aptr(0), 
	U(MatrixViewOf((T*)(0),0,0,ColMajor)), S(0), V(0,0), 
	logdet(0), signdet(0), kmax(0)
      { TMVAssert(FALSE); }
      SVDiv_Impl& operator=(const SVDiv_Impl&)
      { TMVAssert(FALSE); return *this; }
  };

#define APTR1 (inplace ? 0 : new T[A.colsize()*A.rowsize()])
#define APTR (inplace ? A.NonConst().ptr() : Aptr1.get())
#define UX (istrans ? \
  (inplace ? A.NonConst().Transpose() : \
   MatrixViewOf(Aptr,A.rowsize(),A.colsize(),ColMajor)) : \
  (inplace ? A.NonConst().View() : \
   MatrixViewOf(Aptr,A.colsize(),A.rowsize(), \
     A.IsSquare() ? BaseStorOf(A) : ColMajor)))

  template <class T> SVDiv<T>::SVDiv_Impl::SVDiv_Impl(
      const GenMatrix<T>& A, bool _inplace) :
    istrans(A.colsize() < A.rowsize()), 
    inplace(_inplace && (A.isrm() || A.iscm())), Aptr1(APTR1), Aptr(APTR),
    U(UX), S(U.rowsize()), V(U.rowsize(),U.rowsize()), 
    logdet(0), signdet(1), kmax(0) {}

#undef UX
#undef APTR

  template <class T> SVDiv<T>::SVDiv(const GenMatrix<T>& A,
      bool inplace) :
    pimpl(new SVDiv_Impl(A,inplace))
  {
    if (pimpl->istrans) {
      if (inplace) TMVAssert(A.Transpose() == pimpl->U); 
      else pimpl->U = A.Transpose();
    } else {
      if (inplace) TMVAssert(A == pimpl->U); 
      else pimpl->U = A;
    }

    SV_Decompose<T>(pimpl->U,pimpl->S.View(),pimpl->V.View(),
	pimpl->logdet, pimpl->signdet,true);

    // Set kmax for actual 0 elements (to within machine precision).
    // Any further cut in the number of singular values to use
    // should be done by the user.
    Thresh(Epsilon<T>());
  }

  template <class T> SVDiv<T>::~SVDiv() { delete pimpl; pimpl=0; }

  template <class T> template <class T1> void SVDiv<T>::DoLDivEq(
      const MatrixView<T1>& m) const
  {
    TMVAssert(m.colsize() == colsize());
    TMVAssert(m.colsize() == rowsize());
    if (pimpl->istrans) 
      SV_RDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,
	  m.Transpose(),m.Transpose());
    else 
      SV_LDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,m,m);
  }

  template <class T> template <class T1> void SVDiv<T>::DoRDivEq(
      const MatrixView<T1>& m) const
  {
    TMVAssert(m.rowsize() == colsize());
    TMVAssert(m.rowsize() == rowsize());

    if (pimpl->istrans) 
      SV_LDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,
	  m.Transpose(),m.Transpose());
    else 
      SV_RDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,m,m);
  }

  template <class T> template <class T1, class T2> void SVDiv<T>::DoLDiv(
      const GenMatrix<T1>& m, const MatrixView<T2>& x) const
  {
    TMVAssert(m.rowsize() == x.rowsize());
    TMVAssert(m.colsize() == colsize());
    TMVAssert(x.colsize() == rowsize());
    if (pimpl->istrans) 
      SV_RDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,
	m.Transpose(),x.Transpose());
    else 
      SV_LDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,m,x);
  }

  template <class T> template <class T1, class T2> void SVDiv<T>::DoRDiv(
      const GenMatrix<T1>& m, const MatrixView<T2>& x) const
  {
    TMVAssert(m.colsize() == x.colsize());
    TMVAssert(m.rowsize() == rowsize());
    TMVAssert(x.rowsize() == colsize());

    if (pimpl->istrans) 
      SV_LDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,
	  m.Transpose(),x.Transpose());
    else 
      SV_RDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,m,x);
  }

  template <class T> T SVDiv<T>::Det() const 
  { 
    if (pimpl->signdet == T(0)) return T(0);
    else return pimpl->signdet * std::exp(pimpl->logdet);  
  }

  template <class T> RT SVDiv<T>::LogDet(T* sign) const 
  { 
    if (sign) *sign = pimpl->signdet;
    return pimpl->logdet;
  }

  template <class T> template <class T1> void SVDiv<T>::DoInverse(
      const MatrixView<T1>& minv) const
  { 
    if (pimpl->istrans) {
      // A^-1 = (Vt S^-1 Ut)T = U* S^-1 V*
      Matrix<T,ColMajor> SinvV = pimpl->V.Conjugate().Rows(0,pimpl->kmax) /
	pimpl->S.SubDiagMatrix(0,pimpl->kmax);
      minv = pimpl->U.Conjugate().Cols(0,pimpl->kmax) * SinvV;
    } else {
      // A^-1 = Vt S^-1 Ut
      Matrix<T,ColMajor> SinvUt = pimpl->U.Adjoint().Rows(0,pimpl->kmax) /
	pimpl->S.SubDiagMatrix(0,pimpl->kmax);
      minv = pimpl->V.Adjoint().Cols(0,pimpl->kmax) * SinvUt;
    }
  }

  template <class T> void SVDiv<T>::DoInverseATA(
      const MatrixView<T>& minv) const
  {
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
    Matrix<T,ColMajor> SinvV = pimpl->V.Rows(0,pimpl->kmax) /
      pimpl->S.SubDiagMatrix(0,pimpl->kmax);
    if (pimpl->istrans)
      minv = SinvV.Transpose() * SinvV.Conjugate();
    else
      minv = SinvV.Adjoint() * SinvV;
  }

  template <class T> bool SVDiv<T>::Singular() const
  { return pimpl->kmax < int(pimpl->S.size()); }
  
  template <class T> RT SVDiv<T>::Norm2() const 
  { return pimpl->S.size() > 0 ? pimpl->S(0) : RT(0) ; }

  template <class T> RT SVDiv<T>::Condition() const 
  { 
    return pimpl->S.size() > 0 ?
      pimpl->S(0)/pimpl->S(pimpl->S.size()-1) : RT(1); 
  }

  template <class T> void SVDiv<T>::Thresh(RT toler,
      std::ostream* debugout) const
  {
    if (pimpl->S.size() == 0) pimpl->kmax = 0;
    else {
      TMVAssert(toler < RT(1));
      RT thresh = pimpl->S(0)*toler;
      for(pimpl->kmax=pimpl->S.size(); 
	  pimpl->kmax>0 && pimpl->S(pimpl->kmax-1)<=thresh; --pimpl->kmax);
      if(debugout) {
	(*debugout)<<"S = "<<pimpl->S<<std::endl;
	(*debugout)<<"Smax = "<<pimpl->S(0)<<", thresh = "<<thresh<<std::endl;
	(*debugout)<<"kmax = "<<pimpl->kmax;
	(*debugout)<<" (S.size = "<<pimpl->S.size()<<")"<<std::endl;
      }
    }
  }

  template <class T> void SVDiv<T>::Top(int neigen,
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

  template <class T> int SVDiv<T>::GetKMax() const
  { return pimpl->kmax; }

  template <class T> ConstMatrixView<T> SVDiv<T>::GetU() const
  {
    if (pimpl->istrans) return pimpl->V.Transpose(); 
    else return pimpl->U.View(); 
  }
  template <class T> ConstDiagMatrixView<RT> SVDiv<T>::GetS() const
  { return pimpl->S.View(); }

  template <class T> ConstMatrixView<T> SVDiv<T>::GetV() const
  {
    if (pimpl->istrans) return pimpl->U.Transpose(); 
    else return pimpl->V.View(); 
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
    Matrix<T> usv = GetU()*GetS()*GetV();
    RT nm = Norm(usv-mm);
    nm /= Norm(GetU())*Norm(GetS())*Norm(GetV());
    RT cond = GetS()(0) / GetS()(GetKMax()-1);
    if (fout) {
      *fout << "USV = "<<usv<<std::endl;
      *fout << "Norm(M-USV)/Norm(USV) = "<<nm;
      *fout <<"  "<<cond<<" * "<<Epsilon<T>()<<std::endl;
    }
    return nm < cond*Epsilon<T>();
  }

  template <class T> size_t SVDiv<T>::colsize() const
  { return pimpl->istrans ? pimpl->U.rowsize() : pimpl->U.colsize(); }

  template <class T> size_t SVDiv<T>::rowsize() const
  { return pimpl->istrans ? pimpl->U.colsize() : pimpl->U.rowsize(); }

#undef RT

#define InstFile "TMV_SVD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


