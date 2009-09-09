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



#include "TMV_Blas.h"
#include "TMV_VIt.h"
#include "TMV_BandMatrix.h"
#include "TMV_BandLUDiv.h"
#include "TMV_BandQRDiv.h"
#include "TMV_BandSVDiv.h"
#include "TMV_BandMatrixArith.h"
#include <iostream>

namespace tmv {

  // 
  // Access
  //

  template <class T> T GenBandMatrix<T>::cref(size_t i, size_t j) const
  {
    TMVAssert(i<colsize() && j<rowsize());
    TMVAssert(okij(i,j));
    const T* mi = 
      isrm() ? cptr() + int(i)*stepi() + j :
      iscm() ? cptr() + i + int(j)*stepj() :
	cptr()+int(i)*stepi()+int(j)*stepj();
    return isconj() ? CONJ(*mi) : *mi;
  }

  template <class T, IndexStyle I> RefType(T) BandMatrixView<T,I>::ref(
      size_t i, size_t j) const
  {
    TMVAssert(this->okij(i,j));
    T* mi = 
      isrm() ? ptr() + int(i)*stepi() + j :
      iscm() ? ptr() + i + int(j)*stepj() :
	ptr()+int(i)*stepi()+int(j)*stepj();
#ifdef TMVFLDEBUG
    TMVAssert(mi >= first);
    TMVAssert(mi < last);
#endif
    return REF(mi,ct());
  }

  template <class T> void GenBandMatrix<T>::NewDivider() const
  {
    switch (this->GetDivType()) {
      case LU : 
	this->SetDiv(new BandLUDiv<T>(*this,this->IsDivInPlace())); break;
      case QR : 
	this->SetDiv(new BandQRDiv<T>(*this,this->IsDivInPlace())); break;
      case SV : 
	this->SetDiv(new BandSVDiv<T>(*this,true,true)); break;
      case SVS :
	this->SetDiv(new BandSVDiv<T>(*this,false,false)); break;
      case SVU :
	this->SetDiv(new BandSVDiv<T>(*this,true,false)); break;
      case SVV :
	this->SetDiv(new BandSVDiv<T>(*this,false,true)); break;
      default : TMVAssert(FALSE);
    }
  }

  size_t BandStorageLength(StorageType s, size_t cs, size_t rs, int lo, int hi)
  {
    TMVAssert(s == RowMajor || s == ColMajor || s == DiagMajor);
    if (s == RowMajor) {
      if (cs == rs) return (cs-1)*(lo+hi+1)+1;
      else if (cs > rs) return std::min(rs+lo,cs)*(lo+hi)+rs;
      else return (cs-1)*(lo+hi)+std::min(cs+hi,rs);
    } else if (s == ColMajor) {
      if (cs == rs) return (cs-1)*(lo+hi+1)+1;
      else if (cs > rs) return (rs-1)*(lo+hi)+std::min(rs+lo,cs);
      else return std::min(cs+hi,rs)*(lo+hi)+cs;
    } else {
      if (cs == rs) return (cs-1)*(lo+hi+1)+1;
      else if (cs > rs) return rs*(lo+hi+1);
      else return std::min(cs+hi,rs)*(lo+hi)+cs;
    }
  }


  //
  // OK? (SubMatrix, etc.)
  //

  template <class T> bool GenBandMatrix<T>::OKSubMatrix(
      int i1, int i2, int j1, int j2, int istep, int jstep) const
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
      std::cout<<"row range ("<<j2-j1<<") must be multiple of jstep (";
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
    if (!okij(i1,j2-jstep)) {
      ok = false;
      std::cout<<"Upper right corner ("<<i1<<','<<j2-jstep;
      std::cout<<") must be in band\n";
    }
    if (!okij(i2-istep,j1)) {
      ok = false;
      std::cout<<"Lower left corner ("<<i2-istep<<','<<j1;
      std::cout<<") must be in band\n";
    }
    if (!okij(i2-istep,j2-jstep)) {
      ok = false;
      std::cout<<"Lower right corner ("<<i2-istep<<','<<j2-jstep;
      std::cout<<") must be in band\n";
    }
    return ok;
  }

  template <class T> bool GenBandMatrix<T>::OKSubVector(
      int i, int j, int istep, int jstep, size_t size) const
  {
    if (size==0) return true;
    bool ok = true;
    if (istep == 0 && jstep == 0) {
      ok = false;
      std::cout<<"istep ("<<istep<<") and jstep ("<<jstep;
      std::cout<<") can not both be 0\n";
    }
    if (i<0 || i >= int(colsize())) {
      ok = false;
      std::cout<<"i ("<<i<<") must be in 0 -- "<<colsize()-1<<std::endl;
    }
    if (j<0 || j >= int(rowsize())) {
      ok = false;
      std::cout<<"j ("<<j<<") must be in 0 -- "<<rowsize()-1<<std::endl;
    }
    int i2 = int(i)+istep*int(size-1);
    int j2 = int(j)+jstep*int(size-1);
    if (i2 < 0 || i2 >= int(colsize())) {
      ok = false;
      std::cout<<"last element's i ("<<i2<<") must be in 0 -- ";
      std::cout<<colsize()-1<<std::endl;
    }
    if (j2 < 0 || j2 >= int(rowsize())) {
      ok = false;
      std::cout<<"last element's j ("<<j2<<") must be in 0 -- ";
      std::cout<<rowsize()-1<<std::endl;
    }
    if (!okij(i,j)) {
      ok = false;
      std::cout<<"First element ("<<i<<','<<j<<") must be in band\n";
    }
    if (!okij(i2,j2)) {
      ok = false;
      std::cout<<"Last element ("<<i2<<','<<j2<<") must be in band\n";
    }
    return ok;
  }
  
  template <class T> bool GenBandMatrix<T>::OKSubBandMatrix(
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
    return ok;
  } 

  template <class T> bool ConstBandMatrixView<T,FortranStyle>::OKSubMatrix(
      int i1, int i2, int j1, int j2, int istep, int jstep) const
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
    if (j1 < 0 || j1 >= int(this->rowsize())) {
      ok = false;
      std::cout<<"first row element ("<<j1<<") must be in 1 -- ";
      std::cout<<this->rowsize()<<std::endl;
    }
    if (j2 < 0 || j2 >= int(this->rowsize())) {
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
    if (!this->okij(i1-1,j2-1)) {
      ok = false;
      std::cout<<"Upper right corner ("<<i1<<','<<j2<<") must be in band\n";
    }
    if (!this->okij(i2-1,j1-1)) {
      ok = false;
      std::cout<<"Lower left corner ("<<i2<<','<<j1<<") must be in band\n";
    }
    if (!this->okij(i2-1,j2-1)) {
      ok = false;
      std::cout<<"Lower right corner ("<<i2<<','<<j2<<") must be in band\n";
    }
    return ok;
  }

  template <class T> bool ConstBandMatrixView<T,FortranStyle>::OKSubVector(
      int i, int j, int istep, int jstep, size_t size) const
  {
    if (size==0) return true;
    bool ok = true;
    if (istep == 0 && jstep == 0) {
      ok = false;
      std::cout<<"istep ("<<istep<<") and jstep ("<<jstep;
      std::cout<<") can not both be 0\n";
    }
    if (i < 1 || i > int(this->colsize())) {
      ok = false;
      std::cout<<"i ("<<i<<") must be in 1 -- "<<this->colsize()<<std::endl;
    }
    if (j < 1 || j > int(this->rowsize())) {
      ok = false;
      std::cout<<"j ("<<j<<") must be in 1 -- "<<this->rowsize()<<std::endl;
    }
    int i2 = int(i)+istep*int(size-1);
    int j2 = int(j)+jstep*int(size-1);
    if (i2 < 1 || i2 > int(this->colsize())) {
      ok = false;
      std::cout<<"last element's i ("<<i2<<") must be in 1 -- ";
      std::cout<<this->colsize()<<std::endl;
    }
    if (j2 < 1 || j2 > int(this->rowsize())) {
      ok = false;
      std::cout<<"last element's j ("<<j2<<") must be in 1 -- ";
      std::cout<<this->rowsize()<<std::endl;
    }
    if (!this->okij(i-1,j-1)) {
      ok = false;
      std::cout<<"First element ("<<i<<','<<j<<") must be in band\n";
    }
    if (!this->okij(i2-1,j2-1)) {
      ok = false;
      std::cout<<"Last element ("<<i2<<','<<j2<<") must be in band\n";
    }
    return ok;
  }
  
  template <class T> bool ConstBandMatrixView<T,FortranStyle>::OKSubBandMatrix(
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
    return ok;
  } 

  template <class T, IndexStyle I> 
    bool ConstBandMatrixView<T,I>::CanLinearize() const
    {
      if (linsize == 1 && (rowsize() != 1 || colsize() != 1)) {
	TMVAssert(this->isrm() || this->iscm());
	if (this->isrm()) linsize = rowsize() + (colsize()-1)*stepi();
	else linsize = colsize() + (rowsize()-1)*stepj();
      }
      return linsize > 0;
    }

  template <class T, IndexStyle I> 
    bool BandMatrixView<T,I>::CanLinearize() const
    {
      if (linsize == 1 && (rowsize() != 1 || colsize() != 1)) {
	TMVAssert(isrm() || iscm());
	if (isrm()) linsize = rowsize() + (colsize()-1)*stepi();
	else linsize = colsize() + (rowsize()-1)*stepj();
      }
      return linsize > 0;
    }

  template <class T> auto_ptr<BaseMatrix<T> > GenBandMatrix<T>::NewCopy() const
  {
    if (isrm()) return auto_ptr<BaseMatrix<T> >(
	new BandMatrix<T,RowMajor>(*this));
    else if (iscm()) return auto_ptr<BaseMatrix<T> >(
	new BandMatrix<T,ColMajor>(*this));
    else return auto_ptr<BaseMatrix<T> >(
	new BandMatrix<T,DiagMajor>(*this));
  }

  template <class T> auto_ptr<BaseMatrix<T> > GenBandMatrix<T>::NewView() const
  {
    return auto_ptr<BaseMatrix<T> >(new ConstBandMatrixView<T>(View()));
  }

  template <class T> 
    auto_ptr<BaseMatrix<T> > GenBandMatrix<T>::NewTranspose() const
    {
      return auto_ptr<BaseMatrix<T> >(new ConstBandMatrixView<T>(Transpose()));
    }

  template <class T> 
    auto_ptr<BaseMatrix<T> > GenBandMatrix<T>::NewConjugate() const
    {
      return auto_ptr<BaseMatrix<T> >(new ConstBandMatrixView<T>(Conjugate()));
    }

  template <class T> 
    auto_ptr<BaseMatrix<T> > GenBandMatrix<T>::NewAdjoint() const
    {
      return auto_ptr<BaseMatrix<T> >(new ConstBandMatrixView<T>(Adjoint()));
    }

  template <class T> 
    auto_ptr<BaseMatrix<T> > GenBandMatrix<T>::NewInverse() const
    {
      auto_ptr<Matrix<T,ColMajor> > minv(
	  new Matrix<T,ColMajor>(rowsize(),colsize()));
      Inverse(minv->View());
      BaseMatrix<T>* ret1 = minv.release();
      auto_ptr<BaseMatrix<T> > ret(ret1);
      return ret;
    }

  template <class T> QuotXB<T,T> GenBandMatrix<T>::QInverse() const
  { return QuotXB<T,T>(T(1),*this); }

  //
  // Norms
  //

  template <class T> RealType(T) GenBandMatrix<T>::NormSq() const
  { 
    RealType(T) sum = 0;
    if (colsize() > 0 && rowsize() > 0) {
      if (isrm()) {
	size_t j1=0;
	size_t j2=nhi()+1;
	size_t k=nlo();
	for(size_t i=0;i<colsize();++i) {
	  sum += row(i,j1,j2).NormSq();
	  if (k>0) --k; else ++j1;
	  if (j2<rowsize()) ++j2;
	  else if (j1==rowsize()) break;
	}
      } else if (iscm()) {
	size_t i1=0;
	size_t i2=nlo()+1;
	size_t k=nhi();
	for(size_t j=0;j<rowsize();++j) {
	  sum += col(j,i1,i2).NormSq();
	  if (k>0) --k; else ++i1;
	  if (i2<colsize()) ++i2;
	  else if (i1==colsize()) break;
	}
      } else {
	for(int i=-nlo();i<=nhi();++i) sum += diag(i).NormSq();
      }
    }
    return sum;
  }

  template <class T> inline RealType(T) NonLapMaxAbsElement(
      const GenBandMatrix<T>& m)
  { 
    RealType(T) max = 0;
    if (m.colsize() > 0 && m.rowsize() > 0) {
      if (m.isrm()) {
	size_t j1=0;
	size_t j2=m.nhi()+1;
	size_t k=m.nlo();
	for(size_t i=0;i<m.colsize();++i) {
	  RealType(T) temp = m.row(i,j1,j2).MaxAbsElement();
	  if (temp > max) max = temp;
	  if (k>0) --k; else ++j1;
	  if (j2<m.rowsize()) ++j2;
	  else if (j1==m.rowsize()) break;
	}
      } else if (m.iscm()) {
	size_t i1=0;
	size_t i2=m.nlo()+1;
	size_t k=m.nhi();
	for(size_t j=0;j<m.rowsize();++j) {
	  RealType(T) temp = m.col(j,i1,i2).MaxAbsElement();
	  if (temp > max) max = temp;
	  if (k>0) --k; else ++i1;
	  if (i2<m.colsize()) ++i2;
	  else if (i1==m.colsize()) break;
	}
      } else {
	for(int i=-m.nlo();i<=m.nhi();++i) {
	  RealType(T) temp = m.diag(i).MaxAbsElement();
	  if (temp > max) max = temp;
	}
      }
    }
    return max;
  }

  // 1-Norm = max_j (sum_i |a_ij|)
  template <class T> inline RealType(T) NonLapNorm1(
      const GenBandMatrix<T>& m)
  { 
    RealType(T) max = 0;
    if (m.colsize() > 0 && m.rowsize() > 0) {
      size_t i1=0;
      size_t i2=m.nlo()+1;
      size_t k=m.nhi();
      for(size_t j=0;j<m.rowsize();++j) {
	RealType(T) temp = m.col(j,i1,i2).Norm1();
	if (temp > max) max = temp;
	if (k>0) --k; else ++i1;
	if (i2<m.colsize()) ++i2;
	else if (i1==m.colsize()) break;
      }
    }
    return max;
  }

  template <class T> inline RealType(T) NonLapNormF(
      const GenBandMatrix<T>& m)
  { return tmv::SQRT(m.NormSq()); }

  template <class T> inline RealType(T) NonLapNormInf(
      const GenBandMatrix<T>& m)
  { return NonLapNorm1(m.Transpose()); }

#ifdef XLAP
  template <class T> inline RealType(T) LapNorm(
      const char c, const GenBandMatrix<T>& m)
  { 
    switch(c) {
      case 'M' : return NonLapMaxAbsElement(m);
      case '1' : return NonLapNorm1(m);
      case 'F' : return NonLapNormF(m);
      case 'I' : return NonLapNormInf(m);
      default : TMVAssert(FALSE); 
    }
    return RealType(T)(0);
  }
#ifdef INST_DOUBLE
  template <> inline double LapNorm(
      const char c, const GenBandMatrix<double>& m)
  { 
    TMVAssert(m.iscm() || (m.isdm() && m.nlo()==1 && m.nhi()==1));
    TMVAssert(m.IsSquare());
    char cc = c;
    int n = m.colsize();
    double norm;
    if (m.iscm()) {
      int kl = m.nlo();
      int ku = m.nhi();
      int lda = m.stepj()+1;
#ifndef LAPNOWORK
      auto_array<double> work(c == 'I' ? new double[n] : 0);
#endif
      norm = LAPNAMEX(dlangb) (LAPCM LAPV(cc),LAPV(n),LAPV(kl),LAPV(ku),
	  LAPP(m.cptr()-m.nhi()),LAPV(lda) LAPWK(work.get()));
    } else {
      norm = LAPNAMEX(dlangt) (LAPCM LAPV(cc),LAPV(n),
	  LAPP(m.cptr()+m.stepi()),LAPP(m.cptr()),LAPP(m.cptr()+m.stepj()));
    }
    return norm;
  }
  template <> inline double LapNorm(
      const char c, const GenBandMatrix<std::complex<double> >& m)
  { 
    TMVAssert(m.iscm() || (m.isdm() && m.nlo()==1 && m.nhi()==1));
    TMVAssert(m.IsSquare());
    char cc = c;
    int n = m.colsize();
    double norm;
    if (m.iscm()) {
      int kl = m.nlo();
      int ku = m.nhi();
      int lda = m.stepj()+1;
#ifndef LAPNOWORK
      auto_array<double> work(c == 'I' ? new double[n] : 0);
#endif
      norm = LAPNAMEX(zlangb) (LAPCM LAPV(cc),LAPV(n),LAPV(kl),LAPV(ku),
	  LAPP(m.cptr()-m.nhi()),LAPV(lda) LAPWK(work.get()));
    } else {
      norm = LAPNAMEX(zlangt) (LAPCM LAPV(cc),LAPV(n),
	  LAPP(m.cptr()+m.stepi()),LAPP(m.cptr()),LAPP(m.cptr()+m.stepj()));
    }
    return norm;
  }
#endif
#ifdef INST_FLOAT
  template <> inline float LapNorm(
      const char c, const GenBandMatrix<float>& m)
  { 
    TMVAssert(m.iscm() || (m.isdm() && m.nlo()==1 && m.nhi()==1));
    TMVAssert(m.IsSquare());
    char cc = c;
    int n = m.colsize();
    float norm;
    if (m.iscm()) {
      int kl = m.nlo();
      int ku = m.nhi();
      int lda = m.stepj()+1;
#ifndef LAPNOWORK
      auto_array<float> work(c == 'I' ? new float[n] : 0);
#endif
      norm = LAPNAMEX(slangb) (LAPCM LAPV(cc),LAPV(n),LAPV(kl),LAPV(ku),
	  LAPP(m.cptr()-m.nhi()),LAPV(lda) LAPWK(work.get()));
    } else {
      norm = LAPNAMEX(slangt) (LAPCM LAPV(cc),LAPV(n),
	  LAPP(m.cptr()+m.stepi()),LAPP(m.cptr()),LAPP(m.cptr()+m.stepj()));
    }
    return norm;
  }
  template <> inline float LapNorm(
      const char c, const GenBandMatrix<std::complex<float> >& m)
  { 
    TMVAssert(m.iscm() || (m.isdm() && m.nlo()==1 && m.nhi()==1));
    TMVAssert(m.IsSquare());
    char cc = c;
    int n = m.colsize();
    float norm;
    if (m.iscm()) {
      int kl = m.nlo();
      int ku = m.nhi();
      int lda = m.stepj()+1;
#ifndef LAPNOWORK
      auto_array<float> work(c == 'I' ? new float[n] : 0);
#endif
      norm = LAPNAMEX(clangb) (LAPCM LAPV(cc),LAPV(n),LAPV(kl),LAPV(ku),
	  LAPP(m.cptr()-m.nhi()),LAPV(lda) LAPWK(work.get()));
    } else {
      norm = LAPNAMEX(clangt) (LAPCM LAPV(cc),LAPV(n),
	  LAPP(m.cptr()+m.stepi()),LAPP(m.cptr()),LAPP(m.cptr()+m.stepj()));
    }
    return norm;
  }
#endif
#endif // XLAP

  template <class T> RealType(T) GenBandMatrix<T>::MaxAbsElement() const
  {
#ifdef XLAP
    if (isrm() && this->IsSquare()) return LapNorm('M',Transpose());
    else if (iscm() && this->IsSquare()) return LapNorm('M',*this);
    else if (isdm() && this->IsSquare() && nlo()==1 && nhi()==1)
      return LapNorm('M',*this);
    else
#endif
      return NonLapMaxAbsElement(*this);
  }
  template <class T> RealType(T) GenBandMatrix<T>::Norm1() const
  {
#ifdef XLAP
    if (isrm() && this->IsSquare()) return LapNorm('I',Transpose());
    else if (iscm() && this->IsSquare()) return LapNorm('1',*this);
    else if (isdm() && this->IsSquare() && nlo()==1 && nhi()==1)
      return LapNorm('1',*this);
    else
#endif
      return NonLapNorm1(*this);
  }
  template <class T> RealType(T) GenBandMatrix<T>::NormF() const
  {
#ifdef XLAP
    if (isrm() && this->IsSquare()) return LapNorm('F',Transpose());
    else if (iscm() && this->IsSquare()) return LapNorm('F',*this);
    else if (isdm() && this->IsSquare() && nlo()==1 && nhi()==1)
      return LapNorm('F',*this);
    else
#endif
      return NonLapNormF(*this);
  }

  //
  // Modifying Functions
  //

  template <class T, IndexStyle I>
    const BandMatrixView<T,I>& BandMatrixView<T,I>::Clip(
	RealType(T) thresh) const
    {
      if (colsize() > 0 && rowsize() > 0) {
	if (isrm()) {
	  size_t j1=0;
	  size_t j2=nhi()+1;
	  size_t k=nlo();
	  for(size_t i=0;i<colsize();++i) {
	    row(i,j1,j2).Clip(thresh);
	    if (k>0) --k; else ++j1;
	    if (j2<rowsize()) ++j2;
	    else if (j1==rowsize()) break;
	  }
	} else if (iscm()) {
	  size_t i1=0;
	  size_t i2=nlo()+1;
	  size_t k=nhi();
	  for(size_t j=0;j<rowsize();++j) {
	    col(j,i1,i2).Clip(thresh);
	    if (k>0) --k; else ++i1;
	    if (i2<colsize()) ++i2;
	    else if (i1==colsize()) break;
	  }
	} else {
	  for(int i=-nlo();i<=nhi();++i) diag(i).Clip(thresh);
	}
      }
      return *this;
    }

  template <class T, IndexStyle I> 
    const BandMatrixView<T,I>& BandMatrixView<T,I>::SetAllTo(T x) const
    {
      if (colsize() > 0 && rowsize() > 0) {
	if (isrm()) {
	  size_t j1=0;
	  size_t j2=nhi()+1;
	  size_t k=nlo();
	  for(size_t i=0;i<colsize();++i) {
	    row(i,j1,j2).SetAllTo(x);
	    if (k>0) --k; else ++j1;
	    if (j2<rowsize()) ++j2;
	    else if (j1==rowsize()) break;
	  }
	} else if (iscm()) {
	  size_t i1=0;
	  size_t i2=nlo()+1;
	  size_t k=nhi();
	  for(size_t j=0;j<rowsize();++j) {
	    col(j,i1,i2).SetAllTo(x);
	    if (k>0) --k; else ++i1;
	    if (i2<colsize()) ++i2;
	    else if (i1==colsize()) break;
	  }
	} else {
	  for(int i=-nlo();i<=nhi();++i) diag(i).SetAllTo(x);
	}
      }
      return *this;
    }

  template <class T, IndexStyle I> 
    void BandMatrixView<T,I>::DoTransposeSelf() const
    {
      TMVAssert(colsize() == rowsize());
      TMVAssert(nlo() == nhi());
      for(int i=1;i<=nhi();++i) Swap(diag(-i),diag(i));
    }

  template <class T, IndexStyle I> 
    const BandMatrixView<T,I>& BandMatrixView<T,I>::ConjugateSelf() const
    {
      if (IsComplex(T()) && colsize() > 0 && rowsize() > 0) {
	if (isrm()) {
	  size_t j1=0;
	  size_t j2=nhi()+1;
	  size_t k=nlo();
	  for(size_t i=0;i<colsize();++i) {
	    row(i,j1,j2).ConjugateSelf();
	    if (k>0) --k; else ++j1;
	    if (j2<rowsize()) ++j2;
	    else if (j1==rowsize()) break;
	  }
	} else if (iscm()) {
	  size_t i1=0;
	  size_t i2=nlo()+1;
	  size_t k=nhi();
	  for(size_t j=0;j<rowsize();++j) {
	    col(j,i1,i2).ConjugateSelf();
	    if (k>0) --k; else ++i1;
	    if (i2<colsize()) ++i2;
	    else if (i1==colsize()) break;
	  }
	} else {
	  for(int i=-nlo();i<=nhi();++i) diag(i).ConjugateSelf();
	}
      }
      return *this;
    }


  //
  // Special Constructors
  //

  template <class T> BandMatrix<T,DiagMajor> UpperBiDiagMatrix(
      const GenVector<T>& v1, const GenVector<T>& v2)
  {
    if (v1.size() == v2.size()) {
      BandMatrix<T,DiagMajor> temp(v1.size(),v1.size()+1,0,1);
      temp.diag(0) = v1;
      temp.diag(1) = v2;
      return temp;
    } else {
      TMVAssert2(v2.size() == v1.size()-1);
      BandMatrix<T,DiagMajor> temp(v1.size(),v1.size(),0,1);
      temp.diag(0) = v1;
      temp.diag(1) = v2;
      return temp;
    }
  }

  template <class T> BandMatrix<T,DiagMajor> LowerBiDiagMatrix(
      const GenVector<T>& v1, const GenVector<T>& v2)
  {
    if (v1.size() == v2.size()) {
      BandMatrix<T,DiagMajor> temp(v2.size()+1,v2.size(),1,0);
      temp.diag(-1) = v1;
      temp.diag(0) = v2;
      return temp;
    } else {
      TMVAssert2(v1.size() == v2.size()-1);
      BandMatrix<T,DiagMajor> temp(v2.size(),v2.size(),1,0);
      temp.diag(-1) = v1;
      temp.diag(0) = v2;
      return temp;
    }
  }

  template <class T> BandMatrix<T,DiagMajor> TriDiagMatrix(
      const GenVector<T>& v1, const GenVector<T>& v2,
      const GenVector<T>& v3)
  {
    if (v1.size() == v2.size()) {
      TMVAssert2(v3.size() == v2.size()-1);
      BandMatrix<T,DiagMajor> temp(v2.size()+1,v2.size(),1,1);
      temp.diag(-1) = v1;
      temp.diag(0) = v2;
      temp.diag(1) = v3;
      return temp;
    } else if (v2.size() == v3.size()) {
      TMVAssert2(v1.size() == v2.size()-1);
      BandMatrix<T,DiagMajor> temp(v2.size(),v2.size()+1,1,1);
      temp.diag(-1) = v1;
      temp.diag(0) = v2;
      temp.diag(1) = v3;
      return temp;
    } else {
      TMVAssert2(v1.size() == v2.size()-1);
      TMVAssert2(v3.size() == v2.size()-1);
      BandMatrix<T,DiagMajor> temp(v2.size(),v2.size(),1,1);
      temp.diag(-1) = v1;
      temp.diag(0) = v2;
      temp.diag(1) = v3;
      return temp;
    }
  }

  template <class T> ConstBandMatrixView<T> BandMatrixViewOf(
      const T* m, size_t cs, size_t rs, int nlo, int nhi,
      StorageType stor)
  {
    TMVAssert2(stor==RowMajor || stor==ColMajor || stor==DiagMajor);
    TMVAssert2(cs>0);
    TMVAssert2(rs>0);
    TMVAssert2(nlo<int(cs));
    TMVAssert2(nhi<int(rs));
    if (stor == DiagMajor) {
      int stepi = rs >= cs ? -int(cs)+1 : -int(rs);
      int stepj = rs >= cs ? int(cs) : int(rs)+1;
      return ConstBandMatrixView<T>(m-nlo*stepi,cs,rs,nlo,nhi,stepi,stepj,1,
	  DiagMajor,NonConj);
    } else {
      int lohi = nlo+nhi;
      if (stor == RowMajor)
	return ConstBandMatrixView<T>(m,cs,rs,nlo,nhi,lohi,1,lohi+1,
	    RowMajor,NonConj);
      else 
	return ConstBandMatrixView<T>(m,cs,rs,nlo,nhi,1,lohi,lohi+1,
	    ColMajor,NonConj);
    }
  }

  template <class T> BandMatrixView<T> BandMatrixViewOf(
      T* m, size_t cs, size_t rs, int nlo, int nhi,
      StorageType stor)
  {
    TMVAssert2(stor==RowMajor || stor==ColMajor || stor==DiagMajor);
    TMVAssert2(cs>0);
    TMVAssert2(rs>0);
    TMVAssert2(nlo<int(cs));
    TMVAssert2(nhi<int(rs));
    if (stor == DiagMajor) {
      int stepi = rs >= cs ? -int(cs)+1 : -int(rs);
      int stepj = rs >= cs ? int(cs) : int(rs)+1;
      return BandMatrixView<T>(m-nlo*stepi,cs,rs,nlo,nhi,stepi,stepj,1,
	  DiagMajor,NonConj
	  FIRSTLAST1(m,m+BandStorageLength(DiagMajor,cs,rs,nlo,nhi)));
    } else {
      int lohi = nlo+nhi;
      if (stor == RowMajor)
	return BandMatrixView<T>(m,cs,rs,nlo,nhi,lohi,1,lohi+1,RowMajor,NonConj
	    FIRSTLAST1(m,m+BandStorageLength(RowMajor,cs,rs,nlo,nhi)));
      else
	return BandMatrixView<T>(m,cs,rs,nlo,nhi,1,lohi,lohi+1,ColMajor,NonConj
	    FIRSTLAST1(m,m+BandStorageLength(ColMajor,cs,rs,nlo,nhi)));
    }
  }



  //
  // Swap
  //

  template <class T> void Swap(const BandMatrixView<T>& m1,
      const BandMatrixView<T>& m2)
  {
    TMVAssert2(m1.colsize() == m2.colsize());
    TMVAssert2(m1.rowsize() == m2.rowsize());
    TMVAssert2(m1.nlo() == m2.nlo());
    TMVAssert2(m1.nhi() == m2.nhi());
    if (m1.colsize() > 0 && m1.rowsize() > 0) {
      if (m1.isrm() && m2.isrm()) {
	size_t j1=0;
	size_t j2=m1.nhi()+1;
	size_t k=m1.nlo();
	for(size_t i=0;i<m1.colsize();++i) {
	  Swap(m1.row(i,j1,j2),m2.row(i,j1,j2));
	  if (k>0) --k; else ++j1;
	  if (j2<m1.rowsize()) ++j2;
	  else if (j1==m1.rowsize()) break;
	}
      } else if (m1.iscm() && m2.iscm()) {
	size_t i1=0;
	size_t i2=m1.nlo()+1;
	size_t k=m1.nhi();
	for(size_t j=0;j<m1.rowsize();++j) {
	  Swap(m1.col(j,i1,i2),m2.col(j,i1,i2));
	  if (k>0) --k; else ++i1;
	  if (i2<m1.colsize()) ++i2;
	  else if (i1==m1.colsize()) break;
	}
      } else {
	for(int i=-m1.nlo();i<=m1.nhi();++i) Swap(m1.diag(i),m2.diag(i));
      }
    }
  }

  //
  // m1 == m2
  //

  template <class T1, class T2> bool operator==(const GenBandMatrix<T1>& m1,
      const GenBandMatrix<T2>& m2)
  {
    if (m1.colsize() != m2.colsize()) return false;
    if (m1.rowsize() != m2.rowsize()) return false;
    if (m1.nlo() != m2.nlo()) return false;
    if (m1.nhi() != m2.nhi()) return false;
    if (m1.SameAs(m2)) return true;
    for(int i=-m1.nlo();i<=m1.nhi();++i) 
      if (m1.diag(i) != m2.diag(i)) return false;
    return true;
  }

  //
  // I/O
  //
 
  template <bool conj, bool rm, bool compact, bool th, class T> 
    inline void DoWrite(std::ostream& os, const GenBandMatrix<T>& m,
	RealType(T) thresh)
    {
      size_t j1=0;
      size_t len=m.nhi()+1;
      const T* mrowi = m.cptr();
      const int ds = m.diagstep();
      const int si = m.stepi();
      const int sj = m.stepj();
      size_t k = m.nlo();

      if (compact) {
	os << "B "<<m.colsize()<<' '<<m.rowsize();
	os <<' '<<m.nlo()<<' '<<m.nhi()<<std::endl;
      }
      else 
	os << m.colsize()<<' '<<m.rowsize()<<std::endl;

      for(size_t i=m.colsize();i>0;--i) {
	os << "( ";
	if (!compact)
	  for(size_t j=j1;j>0;--j) os << ' '<<T(0)<<' ';

	const T* mij = mrowi;
	for(size_t j=len;j>0;--j,rm?++mij:mij+=sj) 
	  if (conj) 
	    if (th)
	      os << ' ' <<(ABS(*mij)<thresh ? T(0) : CONJ(*mij)) << ' ';
	    else
	      os << ' ' << CONJ(*mij) << ' ';
	  else 
	    if (th)
	      os << ' ' <<(ABS(*mij)<thresh ? T(0) : *mij) << ' ';
	    else
	      os << ' ' << *mij << ' ';

	if (!compact)
	  for(size_t j=m.rowsize()-len-j1;j>0;--j) os << ' '<<T(0)<<' ';
	os << " )\n";
	if (k>0) { --k; ++len; mrowi+=si; } 
	else { 
	  mrowi+=ds; 
	  if (j1<m.rowsize()) ++j1; 
	}
	if (j1+len>m.rowsize()) --len;
      }
    }

  template <bool rm, bool compact, bool th, class T> inline void DoWrite1(
      std::ostream& os, const GenBandMatrix<T>& m, T thresh)
  { DoWrite<false,rm,compact,th>(os,m,thresh); }

  template <bool rm, bool compact, bool th, class T> inline void DoWrite1(
      std::ostream& os, const GenBandMatrix<std::complex<T> >& m, T thresh)
  {
    if (m.isconj())
      DoWrite<true,rm,compact,th>(os,m,thresh);
    else
      DoWrite<false,rm,compact,th>(os,m,thresh);
  }

  template <class T> void GenBandMatrix<T>::Write(std::ostream& os) const
  {
    if (isrm())
      DoWrite1<true,false,false>(os,*this,RealType(T)(0));
    else
      DoWrite1<false,false,false>(os,*this,RealType(T)(0));
  }

  template <class T> void GenBandMatrix<T>::Write(std::ostream& os, 
      RealType(T) thresh) const
  {
    if (isrm())
      DoWrite1<true,false,true>(os,*this,thresh);
    else
      DoWrite1<false,false,true>(os,*this,thresh);
  }

  template <class T> void GenBandMatrix<T>::WriteCompact(std::ostream& os) const
  {
    if (isrm())
      DoWrite1<true,true,false>(os,*this,RealType(T)(0));
    else
      DoWrite1<false,true,false>(os,*this,RealType(T)(0));
  }

  template <class T> void GenBandMatrix<T>::WriteCompact(std::ostream& os,
      RealType(T) thresh) const
  {
    if (isrm())
      DoWrite1<true,true,true>(os,*this,thresh);
    else
      DoWrite1<false,true,true>(os,*this,thresh);
  }

  template <class T> class BandMatrixReadError :
    public ReadError
  {
    public :
      size_t i,j;
      mutable auto_ptr<BandMatrix<T> > m;
      char exp,got;
      size_t cs,rs;
      int lo,hi;
      bool is, iseof, isbad;

      inline BandMatrixReadError(
	  size_t _i, size_t _j, const GenBandMatrix<T>& _m,
	  std::istream& _is) throw() :
	ReadError("BandMatrix"),
	i(_i), j(_j), m(new BandMatrix<T>(_m)), exp(0), got(0),
	cs(_m.colsize()), rs(_m.rowsize()), lo(_m.nlo()), hi(_m.nhi()),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline BandMatrixReadError(std::istream& _is) throw() :
	ReadError("BandMatrix"),
	i(0), j(0), m(0), exp(0), got(0),
	cs(0), rs(0), lo(0), hi(0),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline BandMatrixReadError(
	  size_t _i, size_t _j, const GenBandMatrix<T>& _m,
	  std::istream& _is, char _e, char _g) throw() :
	ReadError("BandMatrix"),
	i(_i), j(_j), m(new BandMatrix<T>(_m)), exp(_e), got(_g),
	cs(_m.colsize()), rs(_m.rowsize()), lo(_m.nlo()), hi(_m.nhi()),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline BandMatrixReadError(std::istream& _is, char _e, char _g) throw() :
	ReadError("BandMatrix"),
	i(0), j(0), m(0), exp(_e), got(_g),
	cs(0), rs(0), lo(0), hi(0),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline BandMatrixReadError(const GenBandMatrix<T>& _m,
	  std::istream& _is, size_t _cs, size_t _rs, int _lo, int _hi) throw() :
	ReadError("BandMatrix"),
	i(0), j(0), m(new BandMatrix<T>(_m)), exp(0), got(0),
	cs(_cs), rs(_rs), lo(_lo), hi(_hi),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

      inline BandMatrixReadError(const BandMatrixReadError<T>& rhs) :
	i(rhs.i), j(rhs.j), m(rhs.m), exp(rhs.exp), got(rhs.got),
	cs(rhs.cs), rs(rhs.rs), lo(rhs.lo), hi(rhs.hi),
	is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
      virtual inline ~BandMatrixReadError() throw() {}

      virtual void Write(std::ostream& os) const throw();
  };

  template <class T> void BandMatrixReadError<T>::Write(
      std::ostream& os) const throw()
  {
    os<<"TMV Read Error: Reading istream input for BandMatrix\n";
    if (exp != got) {
      os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
    }
    if (m.get() && cs != m->colsize()) {
      os<<"Wrong colsize: expected "<<m->colsize()<<", got "<<cs<<".\n";
    }
    if (m.get() && rs != m->rowsize()) {
      os<<"Wrong rowsize: expected "<<m->rowsize()<<", got "<<rs<<".\n";
    }
    if (m.get() && lo != m->nlo()) {
      os<<"Wrong nlo: expected "<<m->nlo()<<", got "<<lo<<".\n";
    }
    if (m.get() && hi != m->nhi()) {
      os<<"Wrong nhi: expected "<<m->nhi()<<", got "<<hi<<".\n";
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
      os<<"The portion of the Bandatrix which was successfully read is: \n";
      ConstBandMatrixView<T> mm = m->View();
      for(size_t ii=0;ii<i;++ii) {
	os<<"( ";
	for(size_t jj=0;jj<mm.rowsize();++jj)
	  os<<' '<<mm(ii,jj)<<' ';
	os<<" )\n";
      }
      os<<"( ";
      for(size_t jj=0;jj<j;++jj)
	os<<' '<<mm(i,jj)<<' ';      
      os<<" )\n";
    }
  }

  template <class T, IndexStyle I> void BandMatrixView<T,I>::Read(
      std::istream& is) const
  {
    char paren;
    size_t j1=0;
    size_t len=nhi()+1;
    T* mrowi = ptr();
    const int ds = diagstep();
    const int si = stepi();
    const int sj = stepj();
    size_t k = nlo();
    for(size_t i=colsize();i>0;--i) {
      is >> paren;
      if (!is || paren != '(') 
	throw BandMatrixReadError<T>(colsize()-i,0,*this,is,
	    '(',is?paren:'(');
      T* mij = mrowi;
      if (isrm()) 
	for(size_t j=len;j>0;--j,++mij) {
	  is >> *mij;
	  if (!is)
	    throw BandMatrixReadError<T>(colsize()-i,j1+len-j,*this,is);
	}
      else 
	for(size_t j=len;j>0;--j,mij+=sj) {
	  is >> *mij;
	  if (!is)
	    throw BandMatrixReadError<T>(colsize()-i,j1+len-j,*this,is);
	}
      is >> paren;
      if ((!is && i>1) || paren != ')') 
	throw BandMatrixReadError<T>(colsize()-i,rowsize(),*this,is,
	    ')',is?paren:')');
      if (k>0) { --k; ++len; mrowi+=si; } 
      else { 
	mrowi+=ds; 
	if (j1<rowsize()) ++j1; 
      }
      if (j1+len>rowsize()) --len;
    }
    if (isconj()) ConjugateSelf();
  }

  template <class T, StorageType S, IndexStyle I> std::istream& operator>>(
      std::istream& is, auto_ptr<BandMatrix<T,S,I> >& m)
  { 
    char b;
    is >> b;
    if (!is)
      throw BandMatrixReadError<T>(is);
    if (b != 'B') 
      throw BandMatrixReadError<T>(is,'B',b);
    size_t cs,rs;
    int lo,hi;
    is >> cs >> rs >> lo >> hi;
    if (!is) 
      throw BandMatrixReadError<T>(is);
    m.reset(new BandMatrix<T,S,I>(cs,rs,lo,hi));
    m->View().Read(is); 
    return is;
  }

  template <class T> std::istream& operator>>(std::istream& is,
      const BandMatrixView<T>& m)
  { 
    char b;
    is >> b;
    if (!is)
      throw BandMatrixReadError<T>(is);
    if (b != 'B') 
      throw BandMatrixReadError<T>(is,'B',b);
    size_t cs,rs;
    int lo,hi;
    is >> cs >> rs >> lo >> hi;
    if (!is)
      throw BandMatrixReadError<T>(is);
    if (cs != m.colsize() || rs != m.rowsize() ||
	lo != m.nlo() || hi != m.nhi())
      throw BandMatrixReadError<T>(m,is,cs,rs,lo,hi);
    TMVAssert(m.colsize() == cs);
    TMVAssert(m.rowsize() == rs);
    TMVAssert(m.nlo() == lo);
    TMVAssert(m.nhi() == hi);
    m.Read(is);
    return is;
  }

#define InstFile "TMV_BandMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


