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



#include "TMV_Blas.h"
#include "tmv/TMV_SymBandMatrix.h"
#include "tmv/TMV_SymBandCHD.h"
#include "tmv/TMV_BandLUD.h"
#include "tmv/TMV_SymBandSVD.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_VIt.h"
#include "tmv/TMV_SymBandMatrixArith.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_BandMatrix.h"
#include <iostream>
#include <string>
#include "portable_platform.h"

namespace tmv {

#define RT RealType(T)
#define CT ComplexType(T)

  //
  // Access
  //

  template <class T> T GenSymBandMatrix<T>::cref(int i, int j) const
  {
    if ((uplo() == Upper && i<=j) || (uplo() == Lower && i>=j)) {
      const T* mi = cptr() + int(i)*stepi() + int(j)*stepj();
      return isconj() ? CONJ(*mi) : *mi;
    } else {
      const T* mi = cptr() + int(j)*stepi() + int(i)*stepj();
      return issym() != isconj() ? *mi : CONJ(*mi);
    }
  }

  template <class T, IndexStyle I> RefType(T) SymBandMatrixView<T,I>::ref(
      int i, int j) const
  {
    if ((uplo() == Upper && i<=j) || (uplo() == Lower && i>=j)) {
      T* mi = ptr() + int(i)*stepi() + int(j)*stepj();
      return REF(mi,ct());
    } else {
      T* mi = ptr() + int(j)*stepi() + int(i)*stepj();
      return this->issym() != this->isconj() ? REF(mi,NonConj) : REF(mi,Conj);
    }
  }

  template <class T> 
  auto_ptr<BaseMatrix<T> > GenSymBandMatrix<T>::NewCopy() const
  {
    auto_ptr<BaseMatrix<T> > a;
    if (issym()) {
      if (uplo() == Upper) {
        if (isrm()) a.reset(new SymBandMatrix<T,Upper,RowMajor>(*this));
        else if (iscm()) a.reset(new SymBandMatrix<T,Upper,ColMajor>(*this));
        else a.reset(new SymBandMatrix<T,Upper,DiagMajor>(*this));
      } else {
        if (isrm()) a.reset(new SymBandMatrix<T,Lower,RowMajor>(*this));
        else if (iscm()) a.reset(new SymBandMatrix<T,Lower,ColMajor>(*this));
        else a.reset(new SymBandMatrix<T,Lower,DiagMajor>(*this));
      }
    } else {
      if (uplo() == Upper) {
        if (isrm()) a.reset(new HermBandMatrix<T,Upper,RowMajor>(*this));
        else if (iscm()) a.reset(new HermBandMatrix<T,Upper,ColMajor>(*this));
        else a.reset(new HermBandMatrix<T,Upper,DiagMajor>(*this));
      } else {
        if (isrm()) a.reset(new HermBandMatrix<T,Lower,RowMajor>(*this));
        else if (iscm()) a.reset(new HermBandMatrix<T,Lower,ColMajor>(*this));
        else a.reset(new HermBandMatrix<T,Lower,DiagMajor>(*this));
      }
    }
    return a;
  }

  template <class T> 
  auto_ptr<BaseMatrix<T> > GenSymBandMatrix<T>::NewView() const
  {
    auto_ptr<BaseMatrix<T> > a(new ConstSymBandMatrixView<T>(View()));
    return a;
  }

  template <class T> 
  auto_ptr<BaseMatrix<T> > GenSymBandMatrix<T>::NewTranspose() const
  {
    auto_ptr<BaseMatrix<T> > a(new ConstSymBandMatrixView<T>(Transpose()));
    return a;
  }

  template <class T> 
  auto_ptr<BaseMatrix<T> > GenSymBandMatrix<T>::NewConjugate() const
  {
    auto_ptr<BaseMatrix<T> > a(new ConstSymBandMatrixView<T>(Conjugate()));
    return a;
  }

  template <class T> 
  auto_ptr<BaseMatrix<T> > GenSymBandMatrix<T>::NewAdjoint() const
  {
    auto_ptr<BaseMatrix<T> > a(new ConstSymBandMatrixView<T>(Adjoint()));
    return a;
  }

  template <class T> 
  auto_ptr<BaseMatrix<T> > GenSymBandMatrix<T>::NewInverse() const
  {
    if (issym()) {
      auto_ptr<SymMatrix<T,Upper,ColMajor> > a(
          new SymMatrix<T,Upper,ColMajor>(size()));
      Inverse(a->View());
      BaseMatrix<T>* ret1 = a.release();
      auto_ptr<BaseMatrix<T> > ret(ret1);
      return ret;
    } else {
      auto_ptr<HermMatrix<T,Upper,ColMajor> > a(
          new HermMatrix<T,Upper,ColMajor>(size()));
      Inverse(a->View());
      BaseMatrix<T>* ret1 = a.release();
      auto_ptr<BaseMatrix<T> > ret(ret1);
      return ret;
    }
  }


  template <class T> void GenSymBandMatrix<T>::NewDivider() const
  { 
    switch(this->GetDivType()) {
      case LU : {
                  this->SetDiv(new BandLUDiv<T>(*this));
                  break; 
                }
      case SV : {
                  if (isherm()) this->SetDiv(new HermBandSVDiv<T>(*this)); 
                  else this->SetDiv(new SymBandSVDiv<T>(*this)); 
                  break; 
                }
      case CH : {
                  this->SetDiv(
                      new HermBandCHDiv<T>(*this,this->IsDivInPlace()));
                  break;
                }
      default : TMVAssert(FALSE); 
    }
  }

  template <class T> QuotXsB<T,T> GenSymBandMatrix<T>::QInverse() const
  { return QuotXsB<T,T>(T(1),*this); }

  template <class T> template <class T1> void GenSymBandMatrix<T>::DoInverse(
      const SymMatrixView<T1>& sinv) const
  {
    TMVAssert(issym() == sinv.issym());
    TMVAssert(isherm() == sinv.isherm());

    this->SetDiv();
    const SymDivider<T>* sdiv = dynamic_cast<const SymDivider<T>*>(
        this->GetDiv());
    TMVAssert(sdiv);
    sdiv->Inverse(sinv);
  }

  //
  // OK? (SubBandMatrix, etc.)
  //

  template <class T> bool GenSymBandMatrix<T>::OKSubMatrix(
      int i1, int i2, int j1, int j2, int istep, int jstep) const
  { 
    if (i1==i2 || j1==j2) return true; // no elements, so whatever...
    bool ok = true;
    int i2x = i2-istep;
    int j2x = j2-jstep;
    if (istep == 0) {
      ok = false;
      std::cout<<"istep ("<<istep<<") can not be 0\n";
    }
    if (i1 < 0 || i1 >= int(size())) {
      ok = false;
      std::cout<<"first col element ("<<i1<<") must be in 0 -- ";
      std::cout<<size()-1<<std::endl;
    }
    if (i2x < 0 || i2x >= int(size())) {
      ok = false;
      std::cout<<"last col element ("<<i2x<<") must be in 0 -- ";
      std::cout<<size()-1<<std::endl;
    }
    if ((i2-i1)%istep != 0) {
      ok = false;
      std::cout<<"col range ("<<i2-i1<<") must be multiple of istep (";
      std::cout<<istep<<")\n";
    }
    if ((i2-i1)/istep < 0) {
      ok = false;
      std::cout<<"n col elements ("<<(i2-i1)/istep<<") must be nonnegative\n";
    }
    if (jstep == 0) {
      ok = false;
      std::cout<<"jstep ("<<jstep<<") can not be 0\n";
    }
    if (j1 < 0 || j1 >= int(size())) {
      ok = false;
      std::cout<<"first row element ("<<j1<<") must be in 0 -- ";
      std::cout<<size()-1<<std::endl;
    }
    if (j2-jstep < 0 || j2-jstep >= int(size())) {
      ok = false;
      std::cout<<"last row element ("<<j2-jstep<<") must be in 0 -- ";
      std::cout<<size()-1<<std::endl;
    }
    if ((j2-j1)%jstep != 0) {
      ok = false;
      std::cout<<"row range ("<<j2-j1<<") must be multiple of jstep (";
      std::cout<<jstep<<")\n";
    }
    if ((j2-j1)/jstep < 0) {
      ok = false;
      std::cout<<"n row elements ("<<(j2-j1)/jstep<<") must be nonnegative\n";
    }
    if ((i1<j1 && i2x>j2x) || (i1>j1 && i2x<j2x)) {
      ok = false;
      std::cout<<"Upper left ("<<i1<<','<<j1<<") and lower right (";
      std::cout<<i2x<<','<<j2x<<") corners must be in same triangle\n";
    }
    if ((i2x<j1 && i1>j2x) || (i2x>j1 && i1<j2x)) {
      ok = false;
      std::cout<<"Upper right ("<<i1<<','<<j2x<<") and lower left (";
      std::cout<<i2x<<','<<j1<<") corners must be in same triangle\n";
    }
    if (!okij(i1,j2x)) {
      ok = false;
      std::cout<<"Upper right ("<<i1<<','<<j2x<<") corner must be in band.\n";
    }
    if (!okij(i2x,j1)) {
      ok = false;
      std::cout<<"Lower left ("<<i2x<<','<<j1<<") corner must be in band.\n";
    }
    return ok;
  }

  template <class T> bool GenSymBandMatrix<T>::OKSubBandMatrix(
      int i1, int i2, int j1, int j2, int newnlo, int newnhi,
      int istep, int jstep) const
  {
    if (i1==i2 || j1==j2) return true; // no elements, so whatever...
    bool ok = true;
    if (istep == 0) {
      ok = false;
      std::cout<<"istep ("<<istep<<") can not be 0\n";
    }
    if (i1 < 0 || i1 >= int(colsize())) {
      ok = false;
      std::cout<<"first col element ("<<i1<<") must be in 0 -- ";
      std::cout<<colsize()-1<<std::endl;
    }
    if (i2-istep < 0 || i2-istep >= int(colsize())) {
      ok = false;
      std::cout<<"last col element ("<<i2-istep<<") must be in 0 -- ";
      std::cout<<colsize()-1<<std::endl;
    }
    if ((i2-i1)%istep != 0) {
      ok = false;
      std::cout<<"col range ("<<i2-i1<<") must be multiple of istep (";
      std::cout<<istep<<")\n";
    }
    if ((i2-i1)/istep < 0) {
      ok = false;
      std::cout<<"n col elements ("<<(i2-i1)/istep<<") must be nonnegative\n";
    }
    if (jstep == 0) {
      ok = false;
      std::cout<<"jstep ("<<jstep<<") can not be 0\n";
    }
    if (j1 < 0 || j1 >= int(rowsize())) {
      ok = false;
      std::cout<<"first row element ("<<j1<<") must be in 0 -- ";
      std::cout<<rowsize()-1<<std::endl;
    }
    if (j2-jstep < 0 || j2-jstep >= int(rowsize())) {
      ok = false;
      std::cout<<"last row element ("<<j2-jstep<<") must be in 0 -- ";
      std::cout<<rowsize()-1<<std::endl;
    }
    if ((j2-j1)%jstep != 0) {
      ok = false;
      std::cout<<"row range ("<<j2-j1<<") must be multiple of istep (";
      std::cout<<jstep<<")\n";
    }
    if ((j2-j1)/jstep < 0) {
      ok = false;
      std::cout<<"n row elements ("<<(j2-j1)/jstep<<") must be nonnegative\n";
    }
    if (!okij(i1,j1)) {
      ok = false;
      std::cout<<"Upper left corner ("<<i1<<','<<j1<<") must be in band\n";
    }
    if (!okij(i1,j1+newnhi)) {
      ok = false;
      std::cout<<"Start of top diagonal ("<<i1<<','<<j1+newnhi;
      std::cout<<") must be in band\n";
    }
    if (!okij(i1+newnlo,j1)) {
      ok = false;
      std::cout<<"Start of bottom diagonal ("<<i1+newnlo<<','<<j1;
      std::cout<<") must be in band\n";
    }
    if (newnhi >= j2-j1) {
      ok = false;
      std::cout<<"new nhi ("<<newnhi<<") must be less than the new rowsize (";
      std::cout<<j2-j1<<")\n";
    }
    if (newnlo >= i2-i1) {
      ok = false;
      std::cout<<"new nlo ("<<newnlo<<") must be less than the new colsize (";
      std::cout<<i2-i1<<")\n";
    }
    if ((i1<j1+newnhi && i1+newnlo>j1) || (i1>j1+newnhi && i1+newnlo<j1)) {
      ok = false;
      std::cout<<"Start of top ("<<i1<<','<<j1+newnhi<<") and bottom (";
      std::cout<<i1+newnlo<<','<<j1<<") diagonals must be in same triangle\n";
    }
    return ok;
  }

  template <class T> bool GenSymBandMatrix<T>::OKSubVector(
      int i, int j, int istep, int jstep, int n) const 
  {
    if (n==0) return true;
    bool ok = true;
    if (istep == 0 && jstep == 0) {
      ok = false;
      std::cout<<"istep ("<<istep<<") and jstep ("<<jstep;
      std::cout<<") can not both be 0\n";
    }
    if (i<0 || i >= int(size())) {
      ok = false;
      std::cout<<"i ("<<i<<") must be in 0 -- "<<size()-1<<std::endl;
    }
    if (j<0 || j >= int(size())) {
      ok = false;
      std::cout<<"j ("<<j<<") must be in 0 -- "<<size()-1<<std::endl;
    }
    int i2 = int(i)+istep*int(n-1);
    int j2 = int(j)+jstep*int(n-1);
    if (i2 < 0 || i2 >= int(size())) {
      ok = false;
      std::cout<<"last element's i ("<<i2<<") must be in 0 -- ";
      std::cout<<size()-1<<std::endl;
    }
    if (j2 < 0 || j2 >= int(size())) {
      ok = false;
      std::cout<<"last element's j ("<<j2<<") must be in 0 -- ";
      std::cout<<size()-1<<std::endl;
    }
    if ((i<j && i2>j2) || (i>j && i2<j2)) {
      ok = false;
      std::cout<<"First ("<<i<<','<<j<<") and last ("<<i2<<','<<j2;
      std::cout<<") elements must be in same triangle\n";
    }
    if (!okij(i,j)) {
      ok = false;
      std::cout<<"First ("<<i<<','<<j<<") element must be in band\n";
    }
    if (!okij(i2,j2)) {
      ok = false;
      std::cout<<"Last ("<<i2<<','<<j2<<") element must be in band\n";
    }
    return ok;
  }

  template <class T> bool GenSymBandMatrix<T>::OKSubSymMatrix(
      int i1, int i2, int istep) const 
  {
    if (i1==i2) return true;
    bool ok=true;
    if (istep == 0) {
      ok = false; 
      std::cout<<"istep ("<<istep<<") can not be 0\n";
    }
    if (i1<0 || i1 >= int(size())) {
      ok = false;
      std::cout<<"first diag element ("<<i1<<") must be in 0 -- ";
      std::cout<<size()-1<<std::endl;
    }
    if (i2-istep<0 || i2-istep >= int(size())) {
      ok = false;
      std::cout<<"last diag element ("<<i2-istep<<") must be in 0 -- ";
      std::cout<<size()-1<<std::endl;
    }
    if ((i2-i1)%istep != 0) {
      ok = false;
      std::cout<<"range ("<<i2-i1<<") must be multiple of istep (";
      std::cout<<istep<<")\n";
    }
    if ((i2-i1)/istep < 0) {
      ok = false;
      std::cout<<"n diag elements ("<<(i2-i1)/istep<<") must be nonnegative\n";
    }
    if (!okij(i1,i2-istep)) {
      ok = false;
      std::cout<<"Upper right ("<<i1<<','<<i2-istep;
      std::cout<<") corner must be in band\n";
    }
    return ok;
  }

  template <class T> bool GenSymBandMatrix<T>::OKSubSymBandMatrix(
      int i1, int i2, int newnlo, int istep) const 
  {
    if (i1==i2) return true;
    bool ok=true;
    if (istep == 0) {
      ok = false; 
      std::cout<<"istep ("<<istep<<") can not be 0\n";
    }
    if (i1<0 || i1 >= int(size())) {
      ok = false;
      std::cout<<"first diag element ("<<i1<<") must be in 0 -- ";
      std::cout<<size()-1<<std::endl;
    }
    if (i2-istep<0 || i2-istep >= int(size())) {
      ok = false;
      std::cout<<"last diag element ("<<i2-istep<<") must be in 0 -- ";
      std::cout<<size()-1<<std::endl;
    }
    if ((i2-i1)%istep != 0) {
      ok = false;
      std::cout<<"range ("<<i2-i1<<") must be multiple of istep (";
      std::cout<<istep<<")\n";
    }
    if ((i2-i1)/istep < 0) {
      ok = false;
      std::cout<<"n diag elements ("<<(i2-i1)/istep<<") must be nonnegative\n";
    }
    if (newnlo > nlo()) {
      ok = false;
      std::cout<<"new number of off-diagonals ("<<newnlo<<") must be less ";
      std::cout<<"than or equal to the current value ("<<nlo()<<")\n";
    }
    return ok;
  }

  template <class T> bool ConstSymBandMatrixView<T,FortranStyle>::OKSubMatrix(
      int i1, int i2, int j1, int j2, int istep, int jstep) const
  { 
    if (i1==i2 || j1==j2) return true; // no elements, so whatever...
    bool ok = true;
    if (istep == 0) {
      ok = false;
      std::cout<<"istep ("<<istep<<") can not be 0\n";
    }
    if (i1 < 1 || i1 > int(this->size())) {
      ok = false;
      std::cout<<"first col element ("<<i1<<") must be in 1 -- ";
      std::cout<<this->size()<<std::endl;
    }
    if (i2 < 1 || i2 > int(this->size())) {
      ok = false;
      std::cout<<"last col element ("<<i2-istep<<") must be in 1 -- ";
      std::cout<<this->size()<<std::endl;
    }
    if ((i2-i1)%istep != 0) {
      ok = false;
      std::cout<<"col range ("<<i2-i1<<") must be multiple of istep (";
      std::cout<<istep<<")\n";
    }
    if ((i2-i1)/istep < 0) {
      ok = false;
      std::cout<<"n col elements ("<<(i2-i1)/istep+1<<") must be positive\n";
    }
    if (jstep == 0) {
      ok = false;
      std::cout<<"jstep ("<<jstep<<") can not be 0\n";
    }
    if (j1 < 1 || j1 > int(this->size())) {
      ok = false;
      std::cout<<"first row element ("<<j1<<") must be in 1 -- ";
      std::cout<<this->size()<<std::endl;
    }
    if (j2 < 1 || j2 > int(this->size())) {
      ok = false;
      std::cout<<"last row element ("<<j2-jstep<<") must be in 1 -- ";
      std::cout<<this->size()<<std::endl;
    }
    if ((j2-j1)%jstep != 0) {
      ok = false;
      std::cout<<"row range ("<<j2-j1<<") must be multiple of istep (";
      std::cout<<jstep<<")\n";
    }
    if ((j2-j1)/jstep < 0) {
      ok = false;
      std::cout<<"n row elements ("<<(j2-j1)/jstep+1<<") must be positive\n";
    }
    if ((i1<j1 && i2>j2) || (i1>j1 && i2<j2)) {
      ok = false;
      std::cout<<"Upper left ("<<i1<<','<<j1<<") and lower right (";
      std::cout<<i2<<','<<j2<<") corners must be in same triangle\n";
    }
    if ((i2<j1 && i1>j2) || (i2>j1 && i1<j2)) {
      ok = false;
      std::cout<<"Upper right ("<<i1<<','<<j2<<") and lower left (";
      std::cout<<i2<<','<<j1<<") corners must be in same triangle\n";
    }
    if (!this->okij(i1-1,j2-1)) {
      ok = false;
      std::cout<<"Upper right ("<<i1<<','<<j2<<") corner must be in band.\n";
    }
    if (!this->okij(i2-1,j1-1)) {
      ok = false;
      std::cout<<"Lower left ("<<i2<<','<<j1<<") corner must be in band.\n";
    }
    return ok;
  }

  template <class T> 
  bool ConstSymBandMatrixView<T,FortranStyle>::OKSubBandMatrix(
      int i1, int i2, int j1, int j2, int newnlo, int newnhi,
      int istep, int jstep) const
  {
    if (i1==i2 || j1==j2) return true; // no elements, so whatever...
    bool ok = true;
    if (istep == 0) {
      ok = false;
      std::cout<<"istep ("<<istep<<") can not be 0\n";
    }
    if (i1 < 1 || i1 > int(this->colsize())) {
      ok = false;
      std::cout<<"first col element ("<<i1<<") must be in 1 -- ";
      std::cout<<this->colsize()<<std::endl;
    }
    if (i2 < 1 || i2 > int(this->colsize())) {
      ok = false;
      std::cout<<"last col element ("<<i2<<") must be in 1 -- ";
      std::cout<<this->colsize()<<std::endl;
    }
    if ((i2-i1)%istep != 0) {
      ok = false;
      std::cout<<"col range ("<<i2-i1<<") must be multiple of istep (";
      std::cout<<istep<<")\n";
    }
    if ((i2-i1)/istep < 0) {
      ok = false;
      std::cout<<"n col elements ("<<(i2-i1)/istep+1<<") must be positive\n";
    }
    if (jstep == 0) {
      ok = false;
      std::cout<<"jstep ("<<jstep<<") can not be 0\n";
    }
    if (j1 < 1 || j1 > int(this->rowsize())) {
      ok = false;
      std::cout<<"first row element ("<<j1<<") must be in 1 -- ";
      std::cout<<this->rowsize()<<std::endl;
    }
    if (j2 < 1 || j2 > int(this->rowsize())) {
      ok = false;
      std::cout<<"last row element ("<<j2<<") must be in 1 -- ";
      std::cout<<this->rowsize()<<std::endl;
    }
    if ((j2-j1)%jstep != 0) {
      ok = false;
      std::cout<<"row range ("<<j2-j1<<") must be multiple of istep (";
      std::cout<<jstep<<")\n";
    }
    if ((j2-j1)/jstep < 0) {
      ok = false;
      std::cout<<"n row elements ("<<(j2-j1)/jstep+1<<") must be positive\n";
    }
    if (!this->okij(i1-1,j1-1)) {
      ok = false;
      std::cout<<"Upper left corner ("<<i1<<','<<j1<<") must be in band\n";
    }
    if (!this->okij(i1-1,j1-1+newnhi)) {
      ok = false;
      std::cout<<"Start of top diagonal ("<<i1<<','<<j1+newnhi;
      std::cout<<") must be in band\n";
    }
    if (!this->okij(i1-1+newnlo,j1-1)) {
      ok = false;
      std::cout<<"Start of bottom diagonal ("<<i1+newnlo<<','<<j1;
      std::cout<<") must be in band\n";
    }
    if (newnhi >= j2-j1+1) {
      ok = false;
      std::cout<<"new nhi ("<<newnhi<<") must be less than the new rowsize (";
      std::cout<<j2-j1+1<<")\n";
    }
    if (newnlo >= i2-i1+1) {
      ok = false;
      std::cout<<"new nlo ("<<newnlo<<") must be less than the new colsize (";
      std::cout<<i2-i1+1<<")\n";
    }
    if ((i1<j1+newnhi && i1+newnlo>j1) || (i1>j1+newnhi && i1+newnlo<j1)) {
      ok = false;
      std::cout<<"Start of top ("<<i1<<','<<j1+newnhi<<") and bottom (";
      std::cout<<i1+newnlo<<','<<j1<<") diagonals must be in same triangle\n";
    }
    return ok;
  }

  template <class T> bool ConstSymBandMatrixView<T,FortranStyle>::OKSubVector(
      int i, int j, int istep, int jstep, int n) const 
  {
    if (n==0) return true;
    bool ok = true;
    if (istep == 0 && jstep == 0) {
      ok = false;
      std::cout<<"istep ("<<istep<<") and jstep ("<<jstep;
      std::cout<<") can not both be 0\n";
    }
    if (i < 1 || i > int(this->size())) {
      ok = false;
      std::cout<<"i ("<<i<<") must be in 1 -- "<<this->size()<<std::endl;
    }
    if (i < 1 || j > int(this->size())) {
      ok = false;
      std::cout<<"j ("<<j<<") must be in 1 -- "<<this->size()<<std::endl;
    }
    int i2 = int(i)+istep*int(n-1);
    int j2 = int(j)+jstep*int(n-1);
    if (i2 < 1 || i2 > int(this->size())) {
      ok = false;
      std::cout<<"last element's i ("<<i2<<") must be in 1 -- ";
      std::cout<<this->size()<<std::endl;
    }
    if (j2 < 1 || j2 > int(this->size())) {
      ok = false;
      std::cout<<"last element's j ("<<j2<<") must be in 1 -- ";
      std::cout<<this->size()<<std::endl;
    }
    if ((i<j && i2>j2) || (i>j && i2<j2)) {
      ok = false;
      std::cout<<"First ("<<i<<','<<j<<") and last ("<<i2<<','<<j2;
      std::cout<<") elements must be in same triangle\n";
    }
    if (!this->okij(i-1,j-1)) {
      ok = false;
      std::cout<<"First ("<<i<<','<<j<<") element must be in band\n";
    }
    if (!this->okij(i2-1,j2-1)) {
      ok = false;
      std::cout<<"Last ("<<i2<<','<<j2<<") element must be in band\n";
    }
    return ok;
  }

  template <class T> 
  bool ConstSymBandMatrixView<T,FortranStyle>::OKSubSymMatrix(
      int i1, int i2, int istep) const 
  {
    if (i1==i2) return true;
    bool ok=true;
    if (istep == 0) {
      ok = false; 
      std::cout<<"istep ("<<istep<<") can not be 0\n";
    }
    if (i1<1 || i1 > int(this->size())) {
      ok = false;
      std::cout<<"first diag element ("<<i1<<") must be in 1 -- ";
      std::cout<<this->size()<<std::endl;
    }
    if (i2-istep<1 || i2-istep > int(this->size())) {
      ok = false;
      std::cout<<"last diag element ("<<i2-istep<<") must be in 1 -- ";
      std::cout<<this->size()<<std::endl;
    }
    if ((i2-i1)%istep != 0) {
      ok = false;
      std::cout<<"range ("<<i2-i1<<") must be multiple of istep ("<<istep;
      std::cout<<")\n";
    }
    if ((i2-i1)/istep < 0) {
      ok = false;
      std::cout<<"n diag elements ("<<(i2-i1)/istep+1<<") must be positive\n";
    }
    if (!this->okij(i1-1,i2-1)) {
      ok = false;
      std::cout<<"Upper right ("<<i1<<','<<i2<<") corner must be in band\n";
    }
    return ok;
  }

  template <class T> 
  bool ConstSymBandMatrixView<T,FortranStyle>::OKSubSymBandMatrix(
      int i1, int i2, int newnlo, int istep) const 
  {
    if (i1==i2) return true;
    bool ok=true;
    if (istep == 0) {
      ok = false; 
      std::cout<<"istep ("<<istep<<") can not be 0\n";
    }
    if (i1<1 || i1 > int(this->size())) {
      ok = false;
      std::cout<<"first diag element ("<<i1<<") must be in 1 -- ";
      std::cout<<this->size()<<std::endl;
    }
    if (i2-istep<1 || i2-istep > int(this->size())) {
      ok = false;
      std::cout<<"last diag element ("<<i2-istep<<") must be in 1 -- ";
      std::cout<<this->size()<<std::endl;
    }
    if ((i2-i1)%istep != 0) {
      ok = false;
      std::cout<<"range ("<<i2-i1<<") must be multiple of istep ("<<istep;
      std::cout<<")\n";
    }
    if ((i2-i1)/istep < 0) {
      ok = false;
      std::cout<<"n diag elements ("<<(i2-i1)/istep+1<<") must be positive\n";
    }
    if (newnlo > this->nlo()) {
      ok = false;
      std::cout<<"new number of off-diagonals ("<<newnlo<<") must be less ";
      std::cout<<"than or equal to the current value ("<<this->nlo()<<")\n";
    }
    return ok;
  }


  //
  // Norms
  //

  template <class T> RT GenSymBandMatrix<T>::NormSq(const RT scale) const
  { 
    RT ans = diag().NormSq(scale);
    if (size() > 0 && nlo() > 0) 
      ans += RT(2) * UpperBandOff().NormSq(scale);
    return ans;
  }

  template <class T> static RT NonLapNorm1(
      const GenSymBandMatrix<T>& m) 
  { 
    if (m.nlo() > 0) {
      RT max(0);
      const int N = m.size();
      if (N > 0) {
        int i1=0;
        int i2=m.nlo()+1;
        int k=m.nlo();
        for(int j=0;j<N;++j) {
          RT temp = m.col(j,i1,j).Norm1();
          temp += m.col(j,j,i2).Norm1();
          if (temp > max) max = temp;
          if (k>0) --k; else ++i1;
          if (i2 < N) ++i2;
        }
      }
      return max;
    } else {
      return m.diag().NormInf();
    }
  } 

  template <class T> static RT NonLapNormF(
      const GenSymBandMatrix<T>& m)
  {
    const RT eps = Epsilon<T>();
    const RT inveps = RT(1)/eps;

    RT mmax = m.MaxAbsElement();
    if (mmax == RT(0)) return RT(0);
    else if (mmax * mmax * eps == RT(0)) {
      // Then we need to rescale, since underflow has caused rounding errors
      // Epsilon is a pure power of 2, so no rounding errors from rescaling.
      RT scale = inveps;
      mmax *= scale;
      while (mmax < eps) { scale *= inveps; mmax *= inveps; }
      return SQRT(m.NormSq(scale))/scale;
    } else if (RT(1) / (mmax*mmax) == RT(0)) {
      // Then we have overflow, so we need to rescale:
      RT scale = eps;
      mmax *= scale;
      while (mmax > RT(1)) { scale *= eps; mmax *= eps; }
      return SQRT(m.NormSq(scale))/scale;
    } 
    return SQRT(m.NormSq());
  }

  template <class T> static inline RT NonLapMaxAbsElement(
      const GenSymBandMatrix<T>& m)
  { return m.UpperBand().MaxAbsElement(); }

#ifdef XLAP
  template <class T> static RT LapNorm(const char c, 
      const GenSymBandMatrix<T>& m)
  {
    switch(c) {
      case 'M' : return NonLapMaxAbsElement(m);
      case '1' : return NonLapNorm1(m);
      case 'F' : return NonLapNormF(m);
      default : TMVAssert(FALSE);
    }
    return RT(0);
  }
#ifdef INST_DOUBLE
  template <> double LapNorm(const char c, const GenSymBandMatrix<double>& m)
  {
    TMVAssert(m.iscm() || m.isrm());
    if (m.isrm()) return LapNorm(c,m.Transpose());
    else {
      char cc = c;
      int N = m.size();
      int noff = m.nlo();
      int lda = m.diagstep();
#ifndef LAPNOWORK
      auto_array<double> work(cc == '1' ? new double[N] : 0);
#endif
      const double* mp = m.cptr();
      if (m.uplo()==Upper) mp -= m.nhi();
      return LAPNAME(dlansb) (LAPCM LAPV(cc),
          (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO,
          LAPV(N),LAPV(noff),LAPP(mp),LAPV(lda) LAPWK(work.get()) 
          LAP1 LAP1);
    }
  }
  template <> double LapNorm(const char c, 
      const GenSymBandMatrix<std::complex<double> >& m)
  {
    TMVAssert(m.iscm() || m.isrm());
    if (m.isrm()) return LapNorm(c,m.Transpose());
    else {
      char cc = c;
      int N = m.size();
      int noff = m.nlo();
      int lda = m.diagstep();
#ifndef LAPNOWORK
      auto_array<double> work(cc == '1' ? new double[N] : 0);
#endif
      const std::complex<double>* mp = m.cptr();
      if (m.uplo()==Upper) mp -= m.nhi();
      if (m.isherm()) 
        return LAPNAME(zlanhb) (LAPCM LAPV(cc),
            (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO,
            LAPV(N),LAPV(noff),LAPP(mp),LAPV(lda) LAPWK(work.get()) 
            LAP1 LAP1);
      else 
        return LAPNAME(zlansb) (LAPCM LAPV(cc),
            (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO,
            LAPV(N),LAPV(noff),LAPP(mp),LAPV(lda) LAPWK(work.get()) 
            LAP1 LAP1);
    }
  }
#endif
#ifdef INST_FLOAT
  template <> float LapNorm(const char c, const GenSymBandMatrix<float>& m)
  {
    TMVAssert(m.iscm() || m.isrm());
    if (m.isrm()) return LapNorm(c,m.Transpose());
    else {
      char cc = c;
      int N = m.size();
      int noff = m.nlo();
      int lda = m.diagstep();
#ifndef LAPNOWORK
      auto_array<float> work(cc == '1' ? new float[N] : 0);
#endif
      const float* mp = m.cptr();
      if (m.uplo()==Upper) mp -= m.nhi();
      return LAPNAME(slansb) (LAPCM LAPV(cc),
          (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO,
          LAPV(N),LAPV(noff),LAPP(mp),LAPV(lda) LAPWK(work.get()) 
          LAP1 LAP1);
    }
  }
  template <> float LapNorm(const char c, 
      const GenSymBandMatrix<std::complex<float> >& m)
  {
    TMVAssert(m.iscm() || m.isrm());
    if (m.isrm()) return LapNorm(c,m.Transpose());
    else {
      char cc = c;
      int N = m.size();
      int noff = m.nlo();
      int lda = m.diagstep();
#ifndef LAPNOWORK
      auto_array<float> work(cc == '1' ? new float[N] : 0);
#endif
      const std::complex<float>* mp = m.cptr();
      if (m.uplo()==Upper) mp -= m.nhi();
      if (m.isherm()) 
        return LAPNAME(clanhb) (LAPCM LAPV(cc),
            (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO,
            LAPV(N),LAPV(noff),LAPP(mp),LAPV(lda) LAPWK(work.get()) 
            LAP1 LAP1);
      else 
        return LAPNAME(clansb) (LAPCM LAPV(cc),
            (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO,
            LAPV(N),LAPV(noff),LAPP(mp),LAPV(lda) LAPWK(work.get()) 
            LAP1 LAP1);
    }
  }
#endif
#endif // XLAP

  template <class T> RT GenSymBandMatrix<T>::MaxAbsElement() const
  {
#ifdef XLAP
    if ((isrm() && stepi()>0) || (iscm() && stepj()>0))
      return LapNorm('M',*this);
    else
#endif
      return NonLapMaxAbsElement(*this);
  }
  template <class T> RT GenSymBandMatrix<T>::Norm1() const
  {
#ifdef XLAP
    if ((isrm() && stepi()>0) || (iscm() && stepj()>0))
      return LapNorm('1',*this);
    else
#endif
      return NonLapNorm1(*this);
  }
  template <class T> RT GenSymBandMatrix<T>::NormF() const
  {
#ifdef XLAP
    if ((isrm() && stepi()>0) || (iscm() && stepj()>0))
      return LapNorm('F',*this);
    else
#endif
      return NonLapNormF(*this);
  }

  template <class T> RT GenSymBandMatrix<T>::DoNorm2() const
  {
    if (this->colsize() < this->rowsize()) return Transpose().DoNorm2();
    if (size() == 0) return RT(0);
    DiagMatrix<RT> S(this->size());
    SV_Decompose(*this,S.View());
    return std::abs(S(0));
  }

  template <class T> RT GenSymBandMatrix<T>::DoCondition() const
  {
    if (this->colsize() < this->rowsize()) return Transpose().DoNorm2();
    if (size() == 0) return RT(1);
    DiagMatrix<RT> S(this->size());
    SV_Decompose(*this,S.View());
    return std::abs(S(0)/S(S.size()-1));
  }

  template <class T> 
  SymBandMatrix<T,Upper,DiagMajor> SymTriDiagMatrix(
      const GenVector<T>& v1, const GenVector<T>& v2)
  {
    TMVAssert2(v2.size() == v1.size()-1);
    SymBandMatrix<T,Upper,DiagMajor> temp(v1.size(),1);
    temp.diag() = v1;
    temp.diag(1) = v2;
    return temp;
  }

  template <class T> 
  HermBandMatrix<T,Upper,DiagMajor> HermTriDiagMatrix(
      const GenVector<T>& v1, const GenVector<T>& v2, UpLoType uplo)
  {
    TMVAssert2(v2.size() == v1.size()-1);
    HermBandMatrix<T,Upper,DiagMajor> temp(v1.size(),1);
    temp.diag() = v1;
    if (uplo == Upper)
      temp.diag(1) = v2;
    else
      temp.diag(-1) = v2;
    return temp;
  }

  template <class T> 
  HermBandMatrix<std::complex<T>,Upper,DiagMajor> HermTriDiagMatrix(
      const GenVector<T>& v1, const GenVector<std::complex<T> >& v2,
      UpLoType uplo)
  {
    TMVAssert2(v2.size() == v1.size()-1);
    HermBandMatrix<std::complex<T>,Upper,DiagMajor> temp(v1.size(),1);
    temp.diag() = v1;
    if (uplo == Upper)
      temp.diag(1) = v2;
    else
      temp.diag(-1) = v2;
    return temp;
  }

  template <class T> ConstSymBandMatrixView<T> SymBandMatrixViewOf(
      const T* m, size_t s, int nlo, UpLoType uplo, StorageType stor)
  {
    TMVAssert2(stor==RowMajor || stor==ColMajor || stor==DiagMajor);
    TMVAssert2(s>0);
    TMVAssert2(nlo<int(s));
    if (stor == DiagMajor) 
      return ConstSymBandMatrixView<T>(m + (uplo==Upper ? 0 : nlo*(s-1)),
          s,nlo,-int(s)+1,int(s),1,Sym,uplo,DiagMajor,NonConj);
    else if (stor == RowMajor)
      return ConstSymBandMatrixView<T>(m,s,nlo,nlo,1,nlo+1,
          Sym,uplo,RowMajor,NonConj);
    else
      return ConstSymBandMatrixView<T>(m,s,nlo,1,nlo,nlo+1,
          Sym,uplo,ColMajor,NonConj);
  }

  template <class T> SymBandMatrixView<T> SymBandMatrixViewOf(
      T* m, size_t s, int nlo, UpLoType uplo, StorageType stor)
  {
    TMVAssert2(stor==RowMajor || stor==ColMajor || stor==DiagMajor);
    TMVAssert2(s>0);
    TMVAssert2(nlo<int(s));
    if (stor == DiagMajor) 
      return SymBandMatrixView<T>(m + (uplo==Upper ? 0 : nlo*(s-1)),
          s,nlo,-int(s)+1,int(s),1,Sym,uplo,DiagMajor,NonConj
          FIRSTLAST1(m,m+BandStorageLength(DiagMajor,s,s,nlo,nlo)));
    else if (stor == RowMajor)
      return SymBandMatrixView<T>(m,s,nlo,nlo,1,nlo+1,
          Sym,uplo,RowMajor,NonConj
          FIRSTLAST1(m,m+BandStorageLength(RowMajor,s,s,nlo,nlo)));
    else
      return SymBandMatrixView<T>(m,s,nlo,1,nlo,nlo+1,
          Sym,uplo,ColMajor,NonConj
          FIRSTLAST1(m,m+BandStorageLength(ColMajor,s,s,nlo,nlo)));
  }

  template <class T> ConstSymBandMatrixView<T> HermBandMatrixViewOf(
      const T* m, size_t s, int nlo, UpLoType uplo, StorageType stor)
  {
    TMVAssert2(stor==RowMajor || stor==ColMajor || stor==DiagMajor);
    TMVAssert2(s>0);
    TMVAssert2(nlo<int(s));
    if (stor == DiagMajor) 
      return ConstSymBandMatrixView<T>(m + (uplo==Upper ? 0 : nlo*(s-1)),
          s,nlo,-int(s)+1,int(s),1,Herm,uplo,DiagMajor,NonConj);
    else if (stor == RowMajor)
      return ConstSymBandMatrixView<T>(m,s,nlo,nlo,1,nlo+1,
          Herm,uplo,RowMajor,NonConj);
    else
      return ConstSymBandMatrixView<T>(m,s,nlo,1,nlo,nlo+1,
          Herm,uplo,ColMajor,NonConj);
  }

  template <class T> SymBandMatrixView<T> HermBandMatrixViewOf(
      T* m, size_t s, int nlo, UpLoType uplo, StorageType stor)
  {
    TMVAssert2(stor==RowMajor || stor==ColMajor || stor==DiagMajor);
    TMVAssert2(s>0);
    TMVAssert2(nlo<int(s));
    if (stor == DiagMajor) 
      return SymBandMatrixView<T>(m + (uplo==Upper ? 0 : nlo*(s-1)),
          s,nlo,-int(s)+1,int(s),1,Herm,uplo,DiagMajor,NonConj
          FIRSTLAST1(m,m+BandStorageLength(DiagMajor,s,s,nlo,nlo)));
    else if (stor == RowMajor)
      return SymBandMatrixView<T>(m,s,nlo,nlo,1,nlo+1,
          Herm,uplo,RowMajor,NonConj
          FIRSTLAST1(m,m+BandStorageLength(RowMajor,s,s,nlo,nlo)));
    else
      return SymBandMatrixView<T>(m,s,nlo,1,nlo,nlo+1,
          Herm,uplo,ColMajor,NonConj
          FIRSTLAST1(m,m+BandStorageLength(ColMajor,s,s,nlo,nlo)));
  }


  //
  // I/O
  //

  // This bit is to workaround a bug in pgCC that was fixed in version 7.
  // I don't know if versions earlier than 6.1 had the bug, but 
  // I apply the workaround to all version before 7.
  template <class T> inline T Value(T x) { return x; }
#ifdef PLATFORM_COMPILER_PGI
#if PLATFORM_COMPILER_VERSION < 0x070000
  inline double Value(long double x) { return double(x); }
  inline std::complex<double> Value(std::complex<long double> x) 
  { return std::complex<double>(x); }
#endif
#endif

  template <bool sym, bool conj, bool rm, bool cm, bool compact, bool th, class T>
  static void DoWrite(std::ostream& os, const GenSymBandMatrix<T>& m, 
      RT thresh)
  {
    TMVAssert(m.uplo() == Lower);
    int gaplen1 = 0;
    int rowlen1 = compact ? 1 : 0;
    int rowlen2 = m.nlo()+1;
    int gaplen2 = m.size()-rowlen2;
    const int si = m.stepi();
    const int sj = rm?1:m.stepj();
    const int sd = m.diagstep();
    const T* mrowi = m.cptr();
    int k=m.nlo();

    if (compact) 
      os << (sym?"sB ":"hB ")<< m.size()<<' '<<m.nlo() << std::endl;
    else
      os << m.size() <<' '<< m.size() << std::endl;

    while (rowlen2>0) {
      os << "( ";
      if (!compact) {
        for(int j=gaplen1;j>0;--j) os << ' '<<Value(T(0))<<' ';
      }

      const T* mij = mrowi;
      for(int j=rowlen1;j>0;--j,rm?++mij:mij+=sj) {
        if (conj)
          if (th)
            os << ' '<<Value(ABS(*mij)<thresh ? T(0) : CONJ(*mij))<<' ';
          else
            os << ' '<<Value(CONJ(*mij))<<' ';
        else
          if (th)
            os << ' '<<Value(ABS(*mij)<thresh ? T(0) : *mij)<<' ';
          else
            os << ' '<<Value(*mij)<<' ';
      }

      if (!compact) {
        for(int j=rowlen2;j>0;--j,cm?++mij:mij+=si) {
          if (sym == conj)
            if (th)
              os << ' '<<Value(ABS(*mij)<thresh ? T(0) : CONJ(*mij))<<' ';
            else
              os << ' '<<Value(CONJ(*mij))<<' ';
          else
            if (th)
              os << ' '<<Value(ABS(*mij)<thresh ? T(0) : *mij)<<' ';
            else
              os << ' '<<Value(*mij)<<' ';
        }

        for(int j=gaplen2;j>0;--j) os << ' '<<Value(T(0))<<' ';
      }
      os << " )\n";

      if (k>0) { --k; ++rowlen1; mrowi += si; }
      else { if(!compact) ++gaplen1; mrowi += sd; }
      if (gaplen2 > 0) --gaplen2; 
      else --rowlen2;
    }
  }

  template <bool rm, bool cm, bool compact, bool th, class T> static inline void DoWrite1(
      std::ostream& os, const GenSymBandMatrix<T>& m, T thresh)
  { DoWrite<true,false,rm,cm,compact,th>(os,m,thresh); }

  template <bool rm, bool cm, bool compact, bool th, class T> static void DoWrite1(
      std::ostream& os, const GenSymBandMatrix<std::complex<T> >& m, T thresh)
  {
    if (m.issym())
      if (m.isconj())
        DoWrite<true,true,rm,cm,compact,th>(os,m,thresh);
      else
        DoWrite<true,false,rm,cm,compact,th>(os,m,thresh);
    else
      if (m.isconj())
        DoWrite<false,true,rm,cm,compact,th>(os,m,thresh);
      else
        DoWrite<false,false,rm,cm,compact,th>(os,m,thresh);
  }

  template <class T> void GenSymBandMatrix<T>::Write(std::ostream& os) const
  {
    if (uplo() == Upper) 
      if (issym()) Transpose().Write(os);
      else Adjoint().Write(os);
    else {
      if (isrm())
        DoWrite1<true,false,false,false>(os,*this,RT(0));
      else if (iscm()) 
        DoWrite1<false,true,false,false>(os,*this,RT(0));
      else
        DoWrite1<false,false,false,false>(os,*this,RT(0));
    }
  }

  template <class T> void GenSymBandMatrix<T>::Write(std::ostream& os,
      RT thresh) const
  {
    if (uplo() == Upper) 
      if (issym()) Transpose().Write(os,thresh);
      else Adjoint().Write(os,thresh);
    else 
      if (isrm())
        DoWrite1<true,false,false,true>(os,*this,thresh);
      else if (iscm())
        DoWrite1<false,true,false,true>(os,*this,thresh);
      else
        DoWrite1<false,false,false,true>(os,*this,thresh);
  }

  template <class T> void GenSymBandMatrix<T>::WriteCompact(
      std::ostream& os) const
  {
    if (uplo() == Upper) 
      if (issym()) Transpose().WriteCompact(os);
      else Adjoint().WriteCompact(os);
    else 
      if (isrm())
        DoWrite1<true,false,true,false>(os,*this,RT(0));
      else if (iscm())
        DoWrite1<false,true,true,false>(os,*this,RT(0));
      else
        DoWrite1<false,false,true,false>(os,*this,RT(0));
  }

  template <class T> void GenSymBandMatrix<T>::WriteCompact(std::ostream& os,
      RT thresh) const
  {
    if (uplo() == Upper) 
      if (issym()) Transpose().WriteCompact(os);
      else Adjoint().WriteCompact(os);
    else 
      if (isrm())
        DoWrite1<true,false,true,true>(os,*this,thresh);
      else if (iscm())
        DoWrite1<false,true,true,true>(os,*this,thresh);
      else
        DoWrite1<false,false,true,true>(os,*this,thresh);
  }

#ifndef NOTHROW
  template <class T> class SymBandMatrixReadError :
    public ReadError
  {
  public :
    int i,j;
    mutable auto_ptr<SymBandMatrix<T> > m;
    std::string exp,got;
    size_t s;
    int lo;
    bool is, iseof, isbad;

    SymBandMatrixReadError(
        int _i, int _j, const GenSymBandMatrix<T>& _m,
        std::istream& _is) throw() :
      ReadError("SymBandMatrix."),
      i(_i), j(_j), m(new SymBandMatrix<T>(_m)), exp(""), got(""), 
      s(_m.size()), lo(_m.nlo()),
      is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
    SymBandMatrixReadError(std::istream& _is) throw() :
      ReadError("SymBandMatrix."),
      i(0), j(0), m(0), exp(""), got(""), s(0), lo(0),
      is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
    SymBandMatrixReadError(
        int _i, int _j, const GenSymBandMatrix<T>& _m,
        std::istream& _is, std::string _e, std::string _g) throw() :
      ReadError("SymBandMatrix."),
      i(_i), j(_j), m(new SymBandMatrix<T>(_m)), 
      exp(_e), got(_g), s(_m.size()), lo(_m.nlo()),
      is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
    SymBandMatrixReadError(std::istream& _is, 
        std::string _e, std::string _g) throw() :
      ReadError("SymBandMatrix."),
      i(0), j(0), m(0), exp(_e), got(_g), s(0), lo(0),
      is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
    SymBandMatrixReadError(
        const GenSymBandMatrix<T>& _m, std::istream& _is, 
        size_t _s, int _lo) throw() :
      ReadError("SymBandMatrix."),
      i(0), j(0), m(new SymBandMatrix<T>(_m)), exp(""), got(""), 
      s(_s), lo(_lo),
      is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

    SymBandMatrixReadError(const SymBandMatrixReadError<T>& rhs) :
      i(rhs.i), j(rhs.j), m(rhs.m), exp(rhs.exp), got(rhs.got), 
      s(rhs.s), lo(rhs.lo),
      is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}

    virtual ~SymBandMatrixReadError() throw() {}

    virtual void Write(std::ostream& os) const throw()
    {
      os<<"TMV Read Error: Reading istream input for SymBandMatrix\n";
      if (exp != got) {
        os<<"Wrong format: expected '"<<exp<<"'";
        if (IsReal(T()) && exp == "sB") os<<" (or 'hB')";
        os<<", got '"<<got<<"'.\n";
      }
      if (m.get() && s != m->size()) {
        os<<"Wrong size: expected "<<m->size()<<", got "<<s<<".\n";
      }
      if (m.get() && lo != m->nlo()) {
        os<<"Write nlo: expected "<<m->nlo()<<", got "<<lo<<".\n";
      }
      if (!is) {
        if (iseof) {
          os<<"Input stream reached end-of-file prematurely.\n";
        } else if (isbad) {
          os<<"Input stream is corrupted.\n";
        } else {
          os<<"Input stream cannot read next character.\n";
        }
      }
      if (m.get()) {
        os<<"The portion of the SymBandMatrix which was successfully read is: \n";
        ConstSymBandMatrixView<T> mm = m->View();
        int k = mm.nlo();
        int rowlen = 1;
        int rowstart = 0;
        for(int ii=0;ii<i;++ii) {
          os<<"( ";
          for(int jj=0;jj<rowlen;++jj)
            os<<' '<<mm(ii,jj+rowstart)<<' ';
          os<<" )\n";
          if (k > 0) { ++rowlen; --k; }
          else ++rowstart;
        }
        os<<"( ";
        for(int jj=0;jj<j;++jj)
          os<<' '<<mm(i,jj+rowstart)<<' ';      
        os<<" )\n";
      }
    }
  };

  template <class T> class HermBandMatrixReadError :
    public ReadError
  {
  public :
    int i,j;
    mutable auto_ptr<HermBandMatrix<T> > m;
    std::string exp,got;
    size_t s;
    int lo;
    T dv;
    bool is, iseof, isbad;

    HermBandMatrixReadError(
        int _i, int _j, const GenSymBandMatrix<T>& _m,
        std::istream& _is) throw() :
      ReadError("HermBandMatrix."),
      i(_i), j(_j), m(new HermBandMatrix<T>(_m)), exp(""), got(""), 
      s(_m.size()), lo(_m.nlo()),
      dv(0), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
    HermBandMatrixReadError(std::istream& _is) throw() :
      ReadError("HermBandMatrix."),
      i(0), j(0), m(0), exp(""), got(""), s(0), lo(0),
      dv(0), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
    HermBandMatrixReadError(
        int _i, int _j, const GenSymBandMatrix<T>& _m,
        std::istream& _is, std::string _e, std::string _g) throw() :
      ReadError("HermBandMatrix."),
      i(_i), j(_j), m(new HermBandMatrix<T>(_m)), exp(_e), got(_g), 
      s(_m.size()), lo(_m.nlo()),
      dv(0), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
    HermBandMatrixReadError(
        int _i, int _j, const GenSymBandMatrix<T>& _m,
        std::istream& _is, T _dv) throw() :
      ReadError("HermBandMatrix."),
      i(_i), j(_j), m(new HermBandMatrix<T>(_m)), exp(""), got(""), 
      s(_m.size()), lo(_m.nlo()),
      dv(_dv), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
    HermBandMatrixReadError(std::istream& _is, 
        std::string _e, std::string _g) throw() :
      ReadError("HermBandMatrix."),
      i(0), j(0), m(0), exp(_e), got(_g), s(0), lo(0),
      dv(0), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
    HermBandMatrixReadError(
        const GenSymBandMatrix<T>& _m, std::istream& _is, 
        size_t _s, int _lo) throw() :
      ReadError("HermBandMatrix."),
      i(0), j(0), m(new HermBandMatrix<T>(_m)), exp(""), got(""), 
      s(_s), lo(_lo),
      dv(0), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

    HermBandMatrixReadError(const HermBandMatrixReadError<T>& rhs) :
      i(rhs.i), j(rhs.j), m(rhs.m), exp(rhs.exp), got(rhs.got),
      s(rhs.s), lo(rhs.lo),
      dv(rhs.dv), is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}

    virtual ~HermBandMatrixReadError() throw() {}

    virtual void Write(std::ostream& os) const throw()
    {
      os<<"TMV Read Error: Reading istream input for HermBandMatrix\n";
      if (exp != got) {
        os<<"Wrong format: expected '"<<exp<<"'";
        if (IsReal(T()) && exp == "hB") os<<" (or 'sB')";
        os<<", got '"<<got<<"'.\n";
      }
      if (dv != T(0)) {
        os<<"Non-real value found on diagonal: "<<dv<<std::endl;
      }
      if (m.get() && s != m->size()) {
        os<<"Wrong size: expected "<<m->size()<<", got "<<s<<".\n";
      }
      if (m.get() && s != m->size()) {
        os<<"Wrong nlo: expected "<<m->nlo()<<", got "<<lo<<".\n";
      }
      if (!is) {
        if (iseof) {
          os<<"Input stream reached end-of-file prematurely.\n";
        } else if (isbad) {
          os<<"Input stream is corrupted.\n";
        } else {
          os<<"Input stream cannot read next character.\n";
        }
      }
      if (m.get()) {
        os<<"The portion of the HermBandMatrix which was successfully read is: \n";
        ConstSymBandMatrixView<T> mm = m->View();
        for(int ii=0;ii<i;++ii) {
          os<<"( ";
          for(int jj=0;jj<(ii<j?i+1:i);++jj)
            os<<' '<<mm(ii,jj)<<' ';
          os<<" )\n";
        }
        os<<"( ";
        for(int jj=0;jj<j;++jj)
          os<<' '<<mm(i,jj)<<' ';      
        os<<" )\n";
      }
    }
  };
#endif

  template <class T, IndexStyle I> void SymBandMatrixView<T,I>::Read(
      std::istream& is) const
  {
    if (uplo() == Upper) 
      if (this->issym()) Transpose().Read(is);
      else Adjoint().Read(is);
    else {
      TMVAssert(uplo() == Lower);
      std::string paren;
      T* mrowi = ptr();
      int rowlen = 1;
      int k = nlo();
      const int si = stepi();
      const int sj = stepj();
      const int sd = si+sj;
      for(int i=size();i>0;--i) {
#ifdef NOSTL
        char p;
        is >> p;
        paren = p;
#else
        is >> paren;
#endif
        if (!is || paren != "(") {
          if (this->issym())
#ifdef NOTHROW
          { std::cerr<<"SymBandMatrix ReadError: "<<paren<<" != (\n"; exit(1); }
#else
          throw SymBandMatrixReadError<T>(size()-i,0,*this,is,
              "(",is?paren:"(");
#endif
          else
#ifdef NOTHROW
          { std::cerr<<"HermBandMatrix ReadError: "<<paren<<" != (\n"; exit(1); }
#else
          throw HermBandMatrixReadError<T>(size()-i,0,*this,is,
              "(",is?paren:"(");
#endif
        }
        T* mij = mrowi;
        if (this->isrm()) {
          for(int j=rowlen;j>0;--j,++mij) {
            is >> *mij;
            if (!is) {
              if (this->issym())
#ifdef NOTHROW
              { std::cerr<<"SymBandMatrix ReadError: !is \n"; exit(1); }
#else
              throw SymBandMatrixReadError<T>(size()-i,rowlen-j,*this,is);
#endif
              else
#ifdef NOTHROW
              { std::cerr<<"HermBandMatrix ReadError: !is \n"; exit(1); }
#else
              throw HermBandMatrixReadError<T>(size()-i,rowlen-j,*this,is);
#endif
            }
            if (this->isherm() && j==1 && IMAG(*mij) != RT(0))
#ifdef NOTHROW
            { std::cerr<<"HermBandMatrix ReadError: "<<*mij<" not real\n"; exit(1); }
#else
            throw HermBandMatrixReadError<T>(size()-i,rowlen-j,*this,is,*mij);
#endif
          }
        }
        else {
          for(int j=rowlen;j>0;--j,mij+=sj) {
            is >> *mij;
            if (!is) {
              if (this->issym())
#ifdef NOTHROW
              { std::cerr<<"SymBandMatrix ReadError: !is \n"; exit(1); }
#else
              throw SymBandMatrixReadError<T>(size()-i,rowlen-j,*this,is);
#endif
              else
#ifdef NOTHROW
              { std::cerr<<"HermBandMatrix ReadError: !is \n"; exit(1); }
#else
              throw HermBandMatrixReadError<T>(size()-i,rowlen-j,*this,is);
#endif
            }
            if (this->isherm() && j==1 && IMAG(*mij) != RT(0))
#ifdef NOTHROW
            { std::cerr<<"HermBandMatrix ReadError: "<<*mij<<" not real \n"; exit(1); }
#else
            throw HermBandMatrixReadError<T>(size()-i,rowlen-j,*this,is,*mij);
#endif
          }
        }
#ifdef NOSTL
        is >> p;
        paren = p;
#else
        is >> paren;
#endif
        if ((!is && i>1)  || paren != ")") {
          if (this->issym())
#ifdef NOTHROW
          { std::cerr<<"SymBandMatrix ReadError: "<<paren<<" != )\n"; exit(1); }
#else
          throw SymBandMatrixReadError<T>(size()-i,size()-i+1,*this,is,
              ")",is?paren:")");
#endif
          else
#ifdef NOTHROW
          { std::cerr<<"HermBandMatrix ReadError: "<<paren<<" != ) \n"; exit(1); }
#else
          throw HermBandMatrixReadError<T>(size()-i,size()-i+1,*this,is,
              ")",is?paren:")");
#endif
        }
        if (k > 0) { ++rowlen; --k; mrowi+=si; }
        else { mrowi+=sd; }
      }
    }
    if (this->isconj()) ConjugateSelf();
  }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  std::istream& operator>>(std::istream& is, 
      auto_ptr<SymBandMatrix<T,U,S,I> >& m)
  { 
    std::string sh;
#ifdef NOSTL
    sh = "  ";
    is >> sh[0];
    if (is) is >> sh[1];
#else
    is >> sh;
#endif
    if (!is)
#ifdef NOTHROW
    { std::cerr<<"SymBandMatrix ReadError: !is \n"; exit(1); }
#else
    throw SymBandMatrixReadError<T>(is);
#endif
    if (IsReal(T())) {
      if (sh != "sB" && sh != "hB") 
#ifdef NOTHROW
      { std::cerr<<"SymBandMatrix ReadError: "<<sh<<" != sB\n"; exit(1); }
#else
      throw SymBandMatrixReadError<T>(is,"sB",sh);
#endif
    } else {
      if (sh != "sB") 
#ifdef NOTHROW
      { std::cerr<<"SymBandMatrix ReadError: "<<sh<<" != sB\n"; exit(1); }
#else
      throw SymBandMatrixReadError<T>(is,"sB",sh);
#endif
    }
    size_t size;
    int nlo;
    is >> size >> nlo;
    if (!is) 
#ifdef NOTHROW
    { std::cerr<<"SymBandMatrix ReadError: !is \n"; exit(1); }
#else
    throw SymBandMatrixReadError<T>(is);
#endif
    m.reset(new SymBandMatrix<T,U,S,I>(size,nlo));
    m->View().Read(is); 
    return is;
  }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
  std::istream& operator>>(std::istream& is, 
      auto_ptr<HermBandMatrix<T,U,S,I> >& m)
  { 
    std::string sh;
#ifdef NOSTL
    sh = "  ";
    is >> sh[0];
    if (is) is >> sh[1];
#else
    is >> sh;
#endif
    if (!is)
#ifdef NOTHROW
    { std::cerr<<"HermBandMatrix ReadError: !is \n"; exit(1); }
#else
    throw HermBandMatrixReadError<T>(is);
#endif
    if (IsReal(T())) {
      if (sh != "sB" && sh != "hB") 
#ifdef NOTHROW
      { std::cerr<<"HermBandMatrix ReadError: "<<sh<<" != hB\n"; exit(1); }
#else
      throw HermBandMatrixReadError<T>(is,"hB",sh);
#endif
    } else {
      if (sh != "hB") 
#ifdef NOTHROW
      { std::cerr<<"HermBandMatrix ReadError: "<<sh<<" != hB\n"; exit(1); }
#else
      throw HermBandMatrixReadError<T>(is,"hB",sh);
#endif
    }
    size_t size;
    int nlo;
    is >> size >> nlo;
    if (!is) 
#ifdef NOTHROW
    { std::cerr<<"HermBandMatrix ReadError: !is \n"; exit(1); }
#else
    throw HermBandMatrixReadError<T>(is);
#endif
    m.reset(new HermBandMatrix<T,U,S,I>(size,nlo));
    m->View().Read(is); 
    return is;
  }

  template <class T> std::istream& operator>>(
      std::istream& is, const SymBandMatrixView<T>& m)
  { 
    std::string sh;
#ifdef NOSTL
    sh = "  ";
    is >> sh[0];
    if (is) is >> sh[1];
#else
    is >> sh;
#endif
    if (!is) {
      if (m.issym())
#ifdef NOTHROW
      { std::cerr<<"SymBandMatrix ReadError: !is \n"; exit(1); }
#else
      throw SymBandMatrixReadError<T>(is);
#endif
      else
#ifdef NOTHROW
      { std::cerr<<"HermBandMatrix ReadError: !is \n"; exit(1); }
#else
      throw HermBandMatrixReadError<T>(is);
#endif
    }
    if (IsReal(T())) {
      if (sh != "sB" && sh != "hB") {
        if (m.issym())
#ifdef NOTHROW
        { std::cerr<<"SymBandMatrix ReadError: "<<sh<<" != sB\n"; exit(1); }
#else
        throw SymBandMatrixReadError<T>(is,"sB",sh);
#endif
        else
#ifdef NOTHROW
        { std::cerr<<"HermBandMatrix ReadError: "<<sh<<" != hB\n"; exit(1); }
#else
        throw HermBandMatrixReadError<T>(is,"hB",sh);
#endif
      }
    } else if (m.issym()) {
      if (sh != "sB") 
#ifdef NOTHROW
      { std::cerr<<"SymBandMatrix ReadError: "<<sh<<" != sB\n"; exit(1); }
#else
      throw SymBandMatrixReadError<T>(is,"sB",sh);
#endif
    } else {
      if (sh != "hB") 
#ifdef NOTHROW
      { std::cerr<<"HermBandMatrix ReadError: "<<sh<<" != hB\n"; exit(1); }
#else
      throw HermBandMatrixReadError<T>(is,"hB",sh);
#endif
    }
    size_t s;
    int nlo;
    is >> s >> nlo;
    if (!is) {
      if (m.issym())
#ifdef NOTHROW
      { std::cerr<<"SymBandMatrix ReadError: !is \n"; exit(1); }
#else
      throw SymBandMatrixReadError<T>(is);
#endif
      else
#ifdef NOTHROW
      { std::cerr<<"HermBandMatrix ReadError: !is \n"; exit(1); }
#else
      throw HermBandMatrixReadError<T>(is);
#endif
    }
    if (s != m.size() || nlo != m.nlo()) {
      if (m.issym())
#ifdef NOTHROW
      { std::cerr<<"SymBandMatrix ReadError: Wrong size\n"; exit(1); }
#else
      throw SymBandMatrixReadError<T>(m,is,s,nlo);
#endif
      else
#ifdef NOTHROW
      { std::cerr<<"HermBandMatrix ReadError: Wrong size \n"; exit(1); }
#else
      throw HermBandMatrixReadError<T>(m,is,s,nlo);
#endif
    }
    TMVAssert(m.size() == s);
    TMVAssert(m.nlo() == nlo);
    m.Read(is);
    return is;
  }

#undef RT
#undef CT

#define InstFile "TMV_SymBandMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


