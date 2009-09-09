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
#include "TMV_DiagMatrix.h"
#include "TMV_SymBandMatrix.h"
#include "TMV_SymBandSVDiv.h"
#include "TMV_SVDiv.h"
#include "TMV_MatrixArith.h"
#include "TMV_DiagMatrixArith.h"
#include "TMV_SymMatrixArith.h"
#include "TMV_SymBandMatrixArith.h"
#include "TMV_SymCHDiv_B.h"
#include <ostream>

namespace tmv {

  template <class T> struct HermBandSVDiv<T>::HermBandSVDiv_Impl
  {
    HermBandSVDiv_Impl(const GenSymBandMatrix<T>& m) :
      U(0), S(m.size()), det(T(1))
    {}

    auto_ptr<Matrix<T,ColMajor> > U;
    Vector<RealType(T)> S;
    T det;
    mutable size_t kmax;
  }; // HermBandSVDiv

  template <class T> inline const MatrixView<T>* ZMV() 
  { return (const MatrixView<T>*)(0); }

  template <class T> HermBandSVDiv<T>::HermBandSVDiv(
      const GenSymBandMatrix<T>& A, bool StoreU) :
    pimpl(new HermBandSVDiv_Impl(A))
  {
#ifdef XTEST
    TMVAssert(A.HermOK());
#endif
    TMVAssert(A.isherm());
    TMVAssert(pimpl);

    const size_t N = A.size();
    if (StoreU) {
      pimpl->U.reset(new Matrix<T,ColMajor>(N,N));
      MatrixView<T> Uv = pimpl->U->View();
      HermBandSV_Decompose(A,&Uv,pimpl->S.View(),pimpl->det);
    }
    else {
      HermBandSV_Decompose(A,ZMV<T>(),pimpl->S.View(),pimpl->det);
    }

    Thresh(Epsilon<T>());
  }

  template <class T> HermBandSVDiv<T>::~HermBandSVDiv() { delete pimpl; }

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

  template <class T> template <class T1> void HermBandSVDiv<T>::DoLDivEq2(
      const MatrixView<T1>& m) const
  {
    TMVAssert(pimpl);
    TMVAssert(pimpl->U.get());
    TMVAssert(m.colsize() == colsize());
    SV_LDiv(*pimpl->U,pimpl->S,pimpl->U->Adjoint(),pimpl->kmax,m,m);
  }

  template <class T> template <class T1> void HermBandSVDiv<T>::DoRDivEq2(
      const MatrixView<T1>& m) const
  {
    TMVAssert(pimpl);
    TMVAssert(pimpl->U.get());
    TMVAssert(m.rowsize() == rowsize());
    SV_RDiv(*pimpl->U,pimpl->S,pimpl->U->Adjoint(),pimpl->kmax,m,m);
  }

  template <class T> template <class T1, class T2> 
    void HermBandSVDiv<T>::DoLDiv2(
	const GenMatrix<T1>& m, const MatrixView<T2>& x) const
    {
      TMVAssert(pimpl);
      TMVAssert(pimpl->U.get());
      TMVAssert(m.rowsize() == x.rowsize());
      TMVAssert(m.colsize() == colsize());
      TMVAssert(x.colsize() == rowsize());
      SV_LDiv(*pimpl->U,pimpl->S,pimpl->U->Adjoint(),pimpl->kmax,m,x);
    }

  template <class T> template <class T1, class T2> 
    void HermBandSVDiv<T>::DoRDiv2(
	const GenMatrix<T1>& m, const MatrixView<T2>& x) const
    {
      TMVAssert(pimpl);
      TMVAssert(pimpl->U.get());
      TMVAssert(m.colsize() == x.colsize());
      TMVAssert(m.rowsize() == rowsize());
      TMVAssert(x.rowsize() == colsize());
      SV_RDiv(*pimpl->U,pimpl->S,pimpl->U->Adjoint(),pimpl->kmax,m,x);
    }

  template <class T> T HermBandSVDiv<T>::Det() const
  {
    TMVAssert(pimpl);
    return pimpl->det; 
  }

  template <class T, class T1> inline void HermBandSV_Inverse(
      const GenMatrix<T>& U, const GenVector<RealType(T)>& S, size_t kmax,
      const SymMatrixView<T1>& sinv)
  {
    TMVAssert(sinv.isherm());
    Matrix<T,ColMajor> SinvUt = U.Adjoint().Rows(0,kmax) /
      DiagMatrixViewOf(S.SubVector(0,kmax));
    SymMultMM<false>(T1(1),U.Cols(0,kmax),SinvUt,sinv);
#ifdef XTEST
    TMVAssert(sinv.HermOK());
#endif
  }
  template <class T> inline void HermBandSV_Inverse(
      const GenMatrix<std::complex<T> >&, const GenVector<T>&, size_t,
      const SymMatrixView<T>&)
  { TMVAssert(FALSE); }

  template <class T> template <class T1> void HermBandSVDiv<T>::DoInverse(
      const SymMatrixView<T1>& sinv) const
  {
    TMVAssert(pimpl);
    TMVAssert(sinv.size() == pimpl->S.size());
    TMVAssert(sinv.isherm());
    HermBandSV_Inverse(*pimpl->U,pimpl->S,pimpl->kmax,sinv);
#ifdef XTEST
    TMVAssert(sinv.HermOK());
#endif
  }

  template <class T> template <class T1>  void HermBandSVDiv<T>::DoInverse2(
      const MatrixView<T1>& minv) const
  { 
    TMVAssert(pimpl);
    TMVAssert(pimpl->U.get());
    // A^-1 = U S^-1 Ut
    HermBandSV_Inverse(*pimpl->U,pimpl->S,pimpl->kmax,
	HermMatrixViewOf(minv,Upper));
    if (pimpl->S.size() > 1)
      LowerTriMatrixViewOf(minv).OffDiag() = 
	UpperTriMatrixViewOf(minv).OffDiag().Adjoint();
  }

  template <class T> void HermBandSVDiv<T>::DoInverseATA2(
      const MatrixView<T>& minv) const
  {
    TMVAssert(pimpl);
    TMVAssert(pimpl->U.get());
    // A = U S Ut
    // At = U S Ut
    // AtA = U S^2 Ut
    // (AtA)^-1 = U S^-2 Ut
    //
    Matrix<T,RowMajor> SinvUt = pimpl->U->Adjoint().Rows(0,pimpl->kmax) /
      DiagMatrixViewOf(pimpl->S.SubVector(0,pimpl->kmax));
    minv = SinvUt.Adjoint() * SinvUt;
  }

  template <class T> bool HermBandSVDiv<T>::Singular() const 
  { 
    TMVAssert(pimpl);
    return pimpl->kmax < pimpl->S.size(); 
  }

  template <class T> RealType(T) HermBandSVDiv<T>::Norm2() const 
  { 
    TMVAssert(pimpl);
    return ABS(pimpl->S(0)); 
  }

  template <class T> RealType(T) HermBandSVDiv<T>::Condition() const 
  { 
    TMVAssert(pimpl);
    return ABS(pimpl->S(0)/pimpl->S(pimpl->S.size()-1)); 
  }


  template <class T> void HermBandSVDiv<T>::Thresh(RealType(T) toler,
      std::ostream* debugout) const
  {
    TMVAssert(pimpl);
    TMVAssert(toler < RealType(T)(1));
    RealType(T) thresh = ABS(pimpl->S(0))*toler;
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

  template <class T> void HermBandSVDiv<T>::Top(size_t neigen,
      std::ostream* debugout) const
  {
    TMVAssert(pimpl);
    TMVAssert(neigen > 0);
    TMVAssert(neigen <= pimpl->S.size());
    pimpl->kmax = neigen;
    if(debugout) {
      (*debugout)<<"S = "<<pimpl->S<<std::endl;
      (*debugout)<<"kmax = "<<pimpl->kmax;
      (*debugout)<<" (S.size = "<<pimpl->S.size()<<")"<<std::endl;
    }
  }

  template <class T> size_t HermBandSVDiv<T>::GetKMax() const 
  { 
    TMVAssert(pimpl);
    return pimpl->kmax; 
  }

  template <class T> ConstMatrixView<T> HermBandSVDiv<T>::GetU() const
  {
    TMVAssert(pimpl);
    TMVAssert(pimpl->U.get());
    return pimpl->U->View(); 
  }

  template <class T> ConstVectorView<RealType(T)> HermBandSVDiv<T>::GetS() const
  {
    TMVAssert(pimpl);
    return pimpl->S.View(); 
  }

  template <class T> ConstMatrixView<T> HermBandSVDiv<T>::GetV() const
  {
    TMVAssert(pimpl);
    TMVAssert(pimpl->U.get());
    return pimpl->U->Adjoint(); 
  }

  template <class T> std::string HermBandSVDiv<T>::Type() const
  { return std::string("HermBandSVDiv<") + tmv::Type(T()) + ">"; }

  template <class T> DivType HermBandSVDiv<T>::GetDivType() const 
  { return SV; }

  template <class T> bool HermBandSVDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, std::ostream* fout) const
  {
    Matrix<T> mm = m;
    if (fout) {
      *fout << "HermBandSVDiv:\n";
      *fout << "M = "<<mm<<std::endl;
      *fout << "U = "<<GetU()<<std::endl;
      *fout << "S = "<<GetS()<<std::endl;
      *fout << "V = "<<GetV()<<std::endl;
    }
    Matrix<T> usv = GetU()*DiagMatrixViewOf(GetS())*GetV();
    RealType(T) nm = Norm(usv-mm);
    nm /= Norm(GetU())*Norm(GetS())*Norm(GetV());
    RealType(T) cond = Condition();
    if (fout) {
      *fout << "USV = "<<usv<<std::endl;
      *fout << "Norm(M-USV)/Norm(USV) = "<<nm;
      *fout <<"  "<<cond<<" * "<<Epsilon<T>()<<std::endl;
    }
    return nm < cond*mm.colsize()*Epsilon<T>();
  }

  template <class T> size_t HermBandSVDiv<T>::colsize() const
  { 
    TMVAssert(pimpl);
    return pimpl->S.size(); 
  }

  template <class T> size_t HermBandSVDiv<T>::rowsize() const
  {
    TMVAssert(pimpl);
    return pimpl->S.size(); 
  }

  template <class T> bool HermBandSVDiv<T>::Uisset() const 
  { 
    TMVAssert(pimpl);
    return pimpl->U.get(); 
  }

  template <class T> struct SymBandSVDiv<T>::SymBandSVDiv_Impl
  {
    SymBandSVDiv_Impl(const GenSymBandMatrix<T>& m) :
      U(0), S(m.size()), V(0), det(T(1)) {}

    auto_ptr<Matrix<T,ColMajor> > U;
    Vector<RealType(T)> S;
    auto_ptr<Matrix<T,ColMajor> > V;
    T det;
    mutable size_t kmax;
  }; // SymBandSVDiv

  template <class T> SymBandSVDiv<T>::SymBandSVDiv(
      const GenSymBandMatrix<T>& A, bool StoreU, bool StoreV) :
    pimpl(new SymBandSVDiv_Impl(A))
  {
    TMVAssert(IsComplex(T()));
    const size_t N = A.size();

    auto_ptr<MatrixView<T> > Uv(0);
    auto_ptr<MatrixView<T> > Vv(0);

    if (StoreU) {
      pimpl->U.reset(new Matrix<T,ColMajor>(N,N));
      Uv.reset(new MatrixView<T>(pimpl->U->View()));
    }
    if (StoreV) {
      pimpl->V.reset(new Matrix<T,ColMajor>(N,N));
      Vv.reset(new MatrixView<T>(pimpl->V->View()));
    }

    SymBandSV_Decompose(A,Uv.get(),pimpl->S.View(),Vv.get(),pimpl->det);

    Thresh(Epsilon<T>());
  }

  template <class T> SymBandSVDiv<T>::~SymBandSVDiv() { delete pimpl; }

  template <class T> template <class T1> void SymBandSVDiv<T>::DoLDivEq2(
      const MatrixView<T1>& m) const
  {
    TMVAssert(pimpl);
    TMVAssert(pimpl->U.get() && pimpl->V.get());
    TMVAssert(m.colsize() == pimpl->U->colsize());
    SV_LDiv(*pimpl->U,pimpl->S,*pimpl->V,pimpl->kmax,m,m);
  }

  template <class T> template <class T1> void SymBandSVDiv<T>::DoRDivEq2(
      const MatrixView<T1>& m) const
  {
    TMVAssert(pimpl);
    TMVAssert(pimpl->U.get() && pimpl->V.get());
    TMVAssert(m.rowsize() == pimpl->U->rowsize());
    SV_RDiv(*pimpl->U,pimpl->S,*pimpl->V,pimpl->kmax,m,m);
  }

  template <class T> template <class T1, class T2> void SymBandSVDiv<T>::DoLDiv2(
      const GenMatrix<T1>& m, const MatrixView<T2>& x) const
  {
    TMVAssert(pimpl);
    TMVAssert(pimpl->U.get() && pimpl->V.get());
    TMVAssert(m.rowsize() == x.rowsize());
    TMVAssert(m.colsize() == colsize());
    TMVAssert(x.colsize() == rowsize());
    SV_LDiv(*pimpl->U,pimpl->S,*pimpl->V,pimpl->kmax,m,x);
  }

  template <class T> template <class T1, class T2> void SymBandSVDiv<T>::DoRDiv2(
      const GenMatrix<T1>& m, const MatrixView<T2>& x) const
  {
    TMVAssert(pimpl);
    TMVAssert(pimpl->U.get() && pimpl->V.get());
    TMVAssert(m.colsize() == x.colsize());
    TMVAssert(m.rowsize() == rowsize());
    TMVAssert(x.rowsize() == colsize());
    SV_RDiv(*pimpl->U,pimpl->S,*pimpl->V,pimpl->kmax,m,x);
  }

  template <class T> T SymBandSVDiv<T>::Det() const 
  { 
    TMVAssert(pimpl);
    return pimpl->det; 
  }

  template <class T, class T1> inline void SymBandSV_Inverse(
      const GenMatrix<T>& U, const GenVector<RealType(T)>& S, 
      const GenMatrix<T>& V, size_t kmax, const SymMatrixView<T1>& sinv)
  {
    // A = U S V
    // A^-1 = Vt S^-1 Ut
    Matrix<T,ColMajor> SinvUt = U.Adjoint().Rows(0,kmax) /
      DiagMatrixViewOf(S.SubVector(0,kmax));
    SymMultMM<false>(T1(1),V.Adjoint().Cols(0,kmax),SinvUt,sinv);
  }
  template <class T> inline void SymBandSV_Inverse(
      const GenMatrix<std::complex<T> >&, 
      const GenVector<T>&, const GenMatrix<std::complex<T> >&,
      size_t , const SymMatrixView<T>&)
  { TMVAssert(FALSE); }

  template <class T> template <class T1> void SymBandSVDiv<T>::DoInverse(
      const SymMatrixView<T1>& sinv) const
  {
    TMVAssert(pimpl);
    TMVAssert(sinv.size() == pimpl->S.size());
    TMVAssert(sinv.issym());
    SymBandSV_Inverse(*pimpl->U,pimpl->S,*pimpl->V,pimpl->kmax,sinv);
  }

  template <class T> template <class T1>  void SymBandSVDiv<T>::DoInverse2(
      const MatrixView<T1>& minv) const
  { 
    TMVAssert(pimpl);
    TMVAssert(pimpl->U.get() && pimpl->V.get());
    SymBandSV_Inverse(*pimpl->U,pimpl->S,*pimpl->V,pimpl->kmax,
	SymMatrixViewOf(minv,Upper));
    if (pimpl->S.size() > 1)
      LowerTriMatrixViewOf(minv).OffDiag() = 
	UpperTriMatrixViewOf(minv).OffDiag().Transpose();
  }

  template <class T> void SymBandSVDiv<T>::DoInverseATA2(
      const MatrixView<T>& minv) const
  {
    TMVAssert(pimpl);
    TMVAssert(pimpl->U.get() && pimpl->V.get());
    // A = U S V
    // At = Vt S Ut
    // AtA = Vt S^2 V
    // (AtA)^-1 = Vt S^-2 V
    //
    Matrix<T,RowMajor> SinvV = pimpl->V->Rows(0,pimpl->kmax) /
      DiagMatrixViewOf(pimpl->S.SubVector(0,pimpl->kmax));
    minv = SinvV.Adjoint() * SinvV;
  }

  template <class T> bool SymBandSVDiv<T>::Singular() const 
  { 
    TMVAssert(pimpl);
    return pimpl->kmax < pimpl->S.size(); 
  }

  template <class T> RealType(T) SymBandSVDiv<T>::Norm2() const 
  { 
    TMVAssert(pimpl);
    return pimpl->S(0); 
  }

  template <class T> RealType(T) SymBandSVDiv<T>::Condition() const 
  { 
    TMVAssert(pimpl);
    return pimpl->S(0)/pimpl->S(pimpl->S.size()-1); 
  }

  template <class T> void SymBandSVDiv<T>::Thresh(RealType(T) toler,
      std::ostream* debugout) const
  {
    TMVAssert(pimpl);
    TMVAssert(toler < RealType(T)(1));
    RealType(T) thresh = pimpl->S(0)*toler;
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

  template <class T> void SymBandSVDiv<T>::Top(size_t neigen,
      std::ostream* debugout) const
  {
    TMVAssert(pimpl);
    TMVAssert(neigen > 0);
    TMVAssert(neigen <= pimpl->S.size());
    pimpl->kmax = neigen;
    if(debugout) {
      (*debugout)<<"S = "<<pimpl->S<<std::endl;
      (*debugout)<<"kmax = "<<pimpl->kmax;
      (*debugout)<<" (S.size = "<<pimpl->S.size()<<")"<<std::endl;
    }
  }

  template <class T> size_t SymBandSVDiv<T>::GetKMax() const 
  { 
    TMVAssert(pimpl);
    return pimpl->kmax; 
  }

  template <class T> ConstMatrixView<T> SymBandSVDiv<T>::GetU() const
  { 
    TMVAssert(pimpl);
    TMVAssert(pimpl->U.get());
    return pimpl->U->View(); 
  }

  template <class T> ConstVectorView<RealType(T)> SymBandSVDiv<T>::GetS() const
  { 
    TMVAssert(pimpl);
    return pimpl->S.View(); 
  }

  template <class T> ConstMatrixView<T> SymBandSVDiv<T>::GetV() const
  {
    TMVAssert(pimpl);
    TMVAssert(pimpl->V.get());
    return pimpl->V->View(); 
  }

  template <class T> std::string SymBandSVDiv<T>::Type() const
  { return std::string("SymBandSVDiv<") + tmv::Type(T()) + ">"; }

  template <class T> DivType SymBandSVDiv<T>::GetDivType() const 
  { return SV; }

  template <class T> bool SymBandSVDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, std::ostream* fout) const
  {
    Matrix<T> mm = m;
    if (fout) {
      *fout << "SymBandSVDiv:\n";
      *fout << "M = "<<mm<<std::endl;
      *fout << "U = "<<GetU()<<std::endl;
      *fout << "S = "<<GetS()<<std::endl;
      *fout << "V = "<<GetV()<<std::endl;
    }
    Matrix<T> usv = GetU()*DiagMatrixViewOf(GetS())*GetV();
    RealType(T) nm = Norm(usv-mm);
    nm /= Norm(GetU())*Norm(GetS())*Norm(GetV());
    RealType(T) cond = Condition();
    if (fout) {
      *fout << "USV = "<<usv<<std::endl;
      *fout << "Norm(M-USV)/Norm(USV) = "<<nm;
      *fout <<"  "<<cond<<" * "<<Epsilon<T>()<<std::endl;
    }
    return nm < cond*mm.colsize()*Epsilon<T>();
  }

  template <class T> size_t SymBandSVDiv<T>::colsize() const
  { 
    TMVAssert(pimpl);
    return pimpl->S.size(); 
  }

  template <class T> size_t SymBandSVDiv<T>::rowsize() const
  { 
    TMVAssert(pimpl);
    return pimpl->S.size(); 
  }

  template <class T> bool SymBandSVDiv<T>::Uisset() const 
  { 
    TMVAssert(pimpl);
    return pimpl->U.get(); 
  }

  template <class T> bool SymBandSVDiv<T>::Visset() const 
  {
    TMVAssert(pimpl);
    return pimpl->V.get(); 
  }

#define InstFile "TMV_SymBandSVDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


