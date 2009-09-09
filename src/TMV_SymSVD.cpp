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



#include "TMV_SymSVD.h"
#include "TMV_SymSVDiv.h"
#include "TMV_SymMatrix.h"
#include "TMV_Matrix.h"
#include "TMV_DiagMatrix.h"
#include "TMV_DiagMatrixArith.h"
#include "TMV_MatrixArith.h"
#include "TMV_VectorArith.h"
#include <ostream>

namespace tmv {

#define RT RealType(T)

  template <class T> struct HermSVDiv<T>::HermSVDiv_Impl
  {
    public :
      HermSVDiv_Impl(const GenSymMatrix<T>& m) :
	U(m.size(),m.size()), S(m.size()), 
	logdet(0), signdet(0), donedet(false), kmax(0) {}

      Matrix<T,ColMajor> U;
      DiagMatrix<RT> S;
      mutable RT logdet;
      mutable RT signdet;
      mutable bool donedet;
      mutable int kmax;

    private :
      HermSVDiv_Impl(const HermSVDiv_Impl&) :
	U(0,0), S(0), logdet(0), signdet(1), donedet(false), kmax(0)
      { TMVAssert(FALSE); }
      HermSVDiv_Impl& operator=(const HermSVDiv_Impl&)
      { TMVAssert(FALSE); return *this; }
  }; // HermSVDiv

  template <class T> HermSVDiv<T>::HermSVDiv(const GenSymMatrix<T>& A) :
      pimpl(new HermSVDiv_Impl(A))
  {
#ifdef XTEST
    TMVAssert(A.HermOK());
#endif
    TMVAssert(A.isherm());
    LowerTriMatrixViewOf(pimpl->U) = A.LowerTri();
    HermSV_Decompose<T>(pimpl->U.View(),pimpl->S.View());

    Thresh(Epsilon<T>());
  }

  template <class T> HermSVDiv<T>::~HermSVDiv() 
  { delete pimpl; pimpl=0; }

  template <class T> template <class T1> void HermSVDiv<T>::DoLDivEq(
      const MatrixView<T1>& m) const
  {
    TMVAssert(m.colsize() == colsize());
    SV_LDiv(pimpl->U,pimpl->S,pimpl->U.Adjoint(),pimpl->kmax,m,m);
  }

  template <class T> template <class T1> void HermSVDiv<T>::DoRDivEq(
      const MatrixView<T1>& m) const
  {
    TMVAssert(m.rowsize() == rowsize());
    SV_RDiv(pimpl->U,pimpl->S,pimpl->U.Adjoint(),pimpl->kmax,m,m);
  }

  template <class T> template <class T1, class T2> void HermSVDiv<T>::DoLDiv(
      const GenMatrix<T1>& m, const MatrixView<T2>& x) const
  {
    TMVAssert(m.rowsize() == x.rowsize());
    TMVAssert(m.colsize() == colsize());
    TMVAssert(x.colsize() == rowsize());
    SV_LDiv(pimpl->U,pimpl->S,pimpl->U.Adjoint(),pimpl->kmax,m,x);
  }

  template <class T> template <class T1, class T2> void HermSVDiv<T>::DoRDiv(
      const GenMatrix<T1>& m, const MatrixView<T2>& x) const
  {
    TMVAssert(m.colsize() == x.colsize());
    TMVAssert(m.rowsize() == rowsize());
    TMVAssert(x.rowsize() == colsize());
    SV_RDiv(pimpl->U,pimpl->S,pimpl->U.Adjoint(),pimpl->kmax,m,x);
  }

  template <class T> T HermSVDiv<T>::Det() const
  {
    if (!pimpl->donedet) {
      pimpl->logdet = pimpl->S.LogDet(&pimpl->signdet);
      pimpl->donedet = true;
    }         
    if (pimpl->signdet == T(0)) return T(0);
    else return pimpl->signdet * std::exp(pimpl->logdet);  
  }                  

  template <class T> RealType(T) HermSVDiv<T>::LogDet(T* sign) const
  {
    if (!pimpl->donedet) {
      pimpl->logdet = pimpl->S.LogDet(&pimpl->signdet);
      pimpl->donedet = true;
    }
    if (sign) *sign = pimpl->signdet;
    return pimpl->logdet;  
  }                  
  template <class T> template <class T1> void HermSVDiv<T>::DoInverse(
      const SymMatrixView<T1>& sinv) const
  {
    TMVAssert(sinv.size() == pimpl->S.size());
    TMVAssert(sinv.isherm());
    HermSV_Inverse(pimpl->U,pimpl->S,pimpl->kmax,sinv);
#ifdef XTEST
    TMVAssert(sinv.HermOK());
#endif
  }

  template <class T> template <class T1>  void HermSVDiv<T>::DoInverse(
      const MatrixView<T1>& minv) const
  { 
    // A^-1 = U S^-1 Ut
    HermSV_Inverse(pimpl->U,pimpl->S,pimpl->kmax,HermMatrixViewOf(minv,Upper));
    if (pimpl->S.size() > 1)
      LowerTriMatrixViewOf(minv).OffDiag() = 
	UpperTriMatrixViewOf(minv).OffDiag().Adjoint();
  }

  template <class T> void HermSVDiv<T>::DoInverseATA(
      const MatrixView<T>& minv) const
  {
    // A = U S Ut
    // At = U S Ut
    // AtA = U S^2 Ut
    // (AtA)^-1 = U S^-2 Ut
    //
    Matrix<T,RowMajor> SinvUt = pimpl->U.Adjoint().Rows(0,pimpl->kmax) /
      pimpl->S.SubDiagMatrix(0,pimpl->kmax);
    minv = SinvUt.Adjoint() * SinvUt;
  }

  template <class T> bool HermSVDiv<T>::Singular() const 
  { return pimpl->kmax < int(pimpl->S.size()); }

  template <class T> RT HermSVDiv<T>::Norm2() const 
  { return pimpl->S.size() > 0 ? ABS(pimpl->S(0)) : RT(0); }

  template <class T> RT HermSVDiv<T>::Condition() const 
  { 
    return pimpl->S.size() > 0 ? 
      ABS(pimpl->S(0)/pimpl->S(pimpl->S.size()-1)) : RT(1);
  }

  template <class T> void HermSVDiv<T>::Thresh(RT toler,
      std::ostream* debugout) const
  {
    if (pimpl->S.size() == 0) pimpl->kmax = 0;
    else {
      TMVAssert(toler < RT(1));
      RT thresh = ABS(pimpl->S(0))*toler;
      for(pimpl->kmax=pimpl->S.size(); 
	  pimpl->kmax>0 && ABS(pimpl->S(pimpl->kmax-1))<=thresh; --pimpl->kmax);
      if(debugout) {
	(*debugout)<<"S = "<<pimpl->S<<std::endl;
	(*debugout)<<"Smax = "<<ABS(pimpl->S(0));
	(*debugout)<<", thresh = "<<thresh<<std::endl;
	(*debugout)<<"kmax = "<<pimpl->kmax;
	(*debugout)<<" (S.size = "<<pimpl->S.size()<<")"<<std::endl;
      }
    }
  }

  template <class T> void HermSVDiv<T>::Top(int neigen,
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

  template <class T> int HermSVDiv<T>::GetKMax() const 
  { return pimpl->kmax; }

  template <class T> ConstMatrixView<T> HermSVDiv<T>::GetU() const
  { return pimpl->U.View(); }

  template <class T> DiagMatrix<RT> HermSVDiv<T>::GetS() const
  { 
    DiagMatrix<RT> temp = pimpl->S;
    const size_t N = pimpl->S.size();
    for(size_t i=0;i<N;i++) 
      if (temp(i) < 0) temp(i) = -temp(i);
    return temp;
  }

  template <class T> Matrix<T> HermSVDiv<T>::GetV() const
  {
    Matrix<T> temp = pimpl->U.Adjoint();
    const int N = pimpl->S.size();
    for(int i=0;i<N;i++) 
      if (pimpl->S(i) < 0) temp.row(i) *= RT(-1);
    return temp;
  }

  template <class T> bool HermSVDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, std::ostream* fout) const
  {
    Matrix<T> mm = m;
    if (fout) {
      *fout << "HermSVDiv:\n";
      *fout << "M = "<<mm<<std::endl;
      *fout << "U = "<<GetU()<<std::endl;
      *fout << "S = "<<GetS()<<std::endl;
      *fout << "V = "<<GetV()<<std::endl;
    }
    Matrix<T> usv = GetU()*GetS()*GetV();
    RT nm = Norm(usv-mm);
    nm /= Norm(GetU())*Norm(GetS())*Norm(GetV());
    RT cond = Condition();
    if (fout) {
      *fout << "USV = "<<usv<<std::endl;
      *fout << "Norm(M-USV)/Norm(USV) = "<<nm;
      *fout <<"  "<<cond<<" * "<<Epsilon<T>()<<std::endl;
    }
    return nm < cond*mm.colsize()*Epsilon<T>();
  }

  template <class T> size_t HermSVDiv<T>::colsize() const
  { return pimpl->S.size(); }

  template <class T> size_t HermSVDiv<T>::rowsize() const
  { return pimpl->S.size(); }

  template <class T> struct SymSVDiv<T>::SymSVDiv_Impl
  {
    public :
      SymSVDiv_Impl(const GenSymMatrix<T>& m) :
	U(m.size(),m.size()),
	S(m.size()), V(m.size(),m.size()), logdet(0), signdet(1), kmax(0) {}

      Matrix<T,ColMajor> U;
      DiagMatrix<RT> S;
      Matrix<T,ColMajor> V;
      RT logdet;
      T signdet;
      mutable int kmax;

    private :
      SymSVDiv_Impl(const SymSVDiv_Impl&) :
	U(0,0), S(0), V(0,0), logdet(0), signdet(0), kmax(0)
      { TMVAssert(FALSE); }
      SymSVDiv_Impl& operator=(const SymSVDiv_Impl&)
      { TMVAssert(FALSE); return *this; }
  }; // SymSVDiv

  template <class T> SymSVDiv<T>::SymSVDiv(const GenSymMatrix<T>& A) :
      pimpl(new SymSVDiv_Impl(A))
  {
    TMVAssert(IsComplex(T()));
    LowerTriMatrixViewOf(pimpl->U) = A.LowerTri();
    SymSV_Decompose<T>(pimpl->U.View(),pimpl->S.View(),
	pimpl->V.View(),pimpl->logdet,pimpl->signdet);

    Thresh(Epsilon<T>());
  }

  template <class T> SymSVDiv<T>::~SymSVDiv() 
  { delete pimpl; pimpl=0; }

  template <class T> template <class T1> void SymSVDiv<T>::DoLDivEq(
      const MatrixView<T1>& m) const
  {
    TMVAssert(m.colsize() == pimpl->U.colsize());
    SV_LDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,m,m);
  }

  template <class T> template <class T1> void SymSVDiv<T>::DoRDivEq(
      const MatrixView<T1>& m) const
  {
    TMVAssert(m.rowsize() == pimpl->U.rowsize());
    SV_RDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,m,m);
  }

  template <class T> template <class T1, class T2> void SymSVDiv<T>::DoLDiv(
      const GenMatrix<T1>& m, const MatrixView<T2>& x) const
  {
    TMVAssert(m.rowsize() == x.rowsize());
    TMVAssert(m.colsize() == colsize());
    TMVAssert(x.colsize() == rowsize());
    SV_LDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,m,x);
  }

  template <class T> template <class T1, class T2> void SymSVDiv<T>::DoRDiv(
      const GenMatrix<T1>& m, const MatrixView<T2>& x) const
  {
    TMVAssert(m.colsize() == x.colsize());
    TMVAssert(m.rowsize() == rowsize());
    TMVAssert(x.rowsize() == colsize());
    SV_RDiv(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,m,x);
  }

  template <class T> T SymSVDiv<T>::Det() const 
  { 
    if (pimpl->signdet == T(0)) return T(0);
    else return pimpl->signdet * std::exp(pimpl->logdet);  
  }

  template <class T> RT SymSVDiv<T>::LogDet(T* sign) const 
  { 
    if (sign) *sign = pimpl->signdet;
    return pimpl->logdet;
  }

  template <class T> template <class T1> void SymSVDiv<T>::DoInverse(
      const SymMatrixView<T1>& sinv) const
  {
    TMVAssert(sinv.size() == pimpl->S.size());
    TMVAssert(sinv.issym());
    SymSV_Inverse(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,sinv);
  }

  template <class T> template <class T1>  void SymSVDiv<T>::DoInverse(
      const MatrixView<T1>& minv) const
  { 
    SymSV_Inverse(pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,
	SymMatrixViewOf(minv,Upper));
    if (pimpl->S.size() > 1)
      LowerTriMatrixViewOf(minv).OffDiag() = 
	UpperTriMatrixViewOf(minv).OffDiag().Transpose();
  }

  template <class T> void SymSVDiv<T>::DoInverseATA(
      const MatrixView<T>& minv) const
  {
    // A = U S V
    // At = Vt S Ut
    // AtA = Vt S^2 V
    // (AtA)^-1 = Vt S^-2 V
    //
    Matrix<T,RowMajor> SinvV = pimpl->V.Rows(0,pimpl->kmax) /
      pimpl->S.SubDiagMatrix(0,pimpl->kmax);
    minv = SinvV.Adjoint() * SinvV;
  }

  template <class T> bool SymSVDiv<T>::Singular() const 
  { return pimpl->kmax < int(pimpl->S.size()); }

  template <class T> RT SymSVDiv<T>::Norm2() const 
  { return pimpl->S.size() > 0 ? pimpl->S(0) : RT(0); }

  template <class T> RT SymSVDiv<T>::Condition() const 
  { 
    return pimpl->S.size() > 0 ? 
      pimpl->S(0)/pimpl->S(pimpl->S.size()-1) : RT(1);
  }

  template <class T> void SymSVDiv<T>::Thresh(RT toler,
      std::ostream* debugout) const
  {
    TMVAssert(toler < RT(1));
    RT thresh = pimpl->S(0)*toler;
    for(pimpl->kmax=pimpl->S.size(); 
	pimpl->kmax>0 && ABS(pimpl->S(pimpl->kmax-1))<=thresh; 
	--pimpl->kmax);
    if(debugout) {
      (*debugout)<<"S = "<<pimpl->S<<std::endl;
      (*debugout)<<"Smax = "<<pimpl->S(0)<<", thresh = "<<thresh<<std::endl;
      (*debugout)<<"kmax = "<<pimpl->kmax;
      (*debugout)<<" (S.size = "<<pimpl->S.size()<<")"<<std::endl;
    }
  }

  template <class T> void SymSVDiv<T>::Top(int neigen,
      std::ostream* debugout) const
  {
    TMVAssert(neigen > 0);
    TMVAssert(neigen <= int(pimpl->S.size()));
    pimpl->kmax = neigen;
    if(debugout) {
      (*debugout)<<"S = "<<pimpl->S<<std::endl;
      (*debugout)<<"kmax = "<<pimpl->kmax;
      (*debugout)<<" (S.size = "<<pimpl->S.size()<<")"<<std::endl;
    }
  }

  template <class T> int SymSVDiv<T>::GetKMax() const 
  { return pimpl->kmax; }

  template <class T> ConstMatrixView<T> SymSVDiv<T>::GetU() const
  { return pimpl->U.View(); }

  template <class T> ConstDiagMatrixView<RT> SymSVDiv<T>::GetS() const
  { return pimpl->S.View(); }

  template <class T> ConstMatrixView<T> SymSVDiv<T>::GetV() const
  { return pimpl->V.View(); }

  template <class T> bool SymSVDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, std::ostream* fout) const
  {
    Matrix<T> mm = m;
    if (fout) {
      *fout << "SymSVDiv:\n";
      *fout << "M = "<<mm<<std::endl;
      *fout << "U = "<<GetU()<<std::endl;
      *fout << "S = "<<GetS()<<std::endl;
      *fout << "V = "<<GetV()<<std::endl;
    }
    Matrix<T> usv = GetU()*GetS()*GetV();
    RT nm = Norm(usv-mm);
    nm /= Norm(GetU())*Norm(GetS())*Norm(GetV());
    RT cond = Condition();
    if (fout) {
      *fout << "USV = "<<usv<<std::endl;
      *fout << "Norm(M-USV)/Norm(USV) = "<<nm;
      *fout <<"  "<<cond<<" * "<<Epsilon<T>()<<std::endl;
    }
    return nm < cond*mm.colsize()*Epsilon<T>();
  }

  template <class T> size_t SymSVDiv<T>::colsize() const
  { return pimpl->S.size(); }

  template <class T> size_t SymSVDiv<T>::rowsize() const
  { return pimpl->S.size(); }

#undef RT

#define InstFile "TMV_SymSVD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


