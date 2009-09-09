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



#define XXDEBUG1

#include "TMV_Blas.h"
#include "TMV_VIt.h"
#include "TMV_Matrix.h"
#include "TMV_SymBandMatrix.h"
#include "TMV_SymBandCHDiv.h"
#include "TMV_BandLUDiv.h"
#include "TMV_SymBandSVDiv.h"
#include "TMV_SymBandMatrixArith.h"
#include <algorithm>
#include <iostream>
#include <string>

//#define XDEBUG

#ifdef XDEBUG
using std::cerr;
using std::endl;
#endif

namespace tmv {

  //
  // Access
  //

  template <class T> T GenSymBandMatrix<T>::cref(size_t i, size_t j) const
  {
    TMVAssert(i < size());
    TMVAssert(j < size());
    TMVAssert(okij(i,j));
    const T* mi;
    if ((uplo() == Upper && i<=j) || (uplo() == Lower && i>=j)) {
      if (isrm()) mi = cptr() + int(i)*stepi() + j;
      else if (iscm()) mi = cptr() + i + int(j)*stepj();
      else mi = cptr() + int(i)*stepi() + int(j)*stepj();
      return isconj() ? CONJ(*mi) : *mi;
    } else {
      if (isrm()) mi = cptr() + int(j)*stepi() + i;
      else if (iscm()) mi = cptr() + j + int(i)*stepj();
      else mi = cptr() + int(j)*stepi() + int(i)*stepj();
      return issym() != isconj() ? *mi : CONJ(*mi);
    }
  }

  template <class T, IndexStyle I> RefType(T) SymBandMatrixView<T,I>::ref(
      size_t i, size_t j) const
  {
    TMVAssert(i < size());
    TMVAssert(j < size());
    TMVAssert(okij(i,j));
    T* mi;
    if ((uplo() == Upper && i<=j) || (uplo() == Lower && i>=j)) {
      if (this->isrm()) mi = ptr() + int(i)*stepi() + j;
      else if (this->iscm()) mi = ptr() + i + int(j)*stepj();
      else mi = ptr() + int(i)*stepi() + int(j)*stepj();
    } else {
      if (this->isrm()) mi = ptr() + int(j)*stepi() + i;
      else if (this->iscm()) mi = ptr() + j + int(i)*stepj();
      else mi = ptr() + int(j)*stepi() + int(i)*stepj();
    }
#ifdef TMVFLDEBUG
    TMVAssert(mi>=first);
    TMVAssert(mi<last);
#endif
    if ((uplo() == Upper && i<=j) || (uplo() == Lower && i>=j)) 
      return REF(mi,ct());
    else
      return this->issym() != this->isconj() ? REF(mi,NonConj) : REF(mi,Conj);
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
		  if (isherm()) {
		    this->SetDiv(new HermBandSVDiv<T>(*this,true)); 
		  } else {
		    this->SetDiv(new SymBandSVDiv<T>(*this,true,true)); 
		  }
		  break; 
		}
      case SVU : {
		   if (isherm()) {
		     this->SetDiv(new HermBandSVDiv<T>(*this,true)); 
		   } else {
		     this->SetDiv(new SymBandSVDiv<T>(*this,true,false)); 
		   }
		   break; 
		 }
      case SVV : {
		   if (isherm()) {
		     this->SetDiv(new HermBandSVDiv<T>(*this,true)); 
		   } else {
		     this->SetDiv(new SymBandSVDiv<T>(*this,false,true)); 
		   }
		   break; 
		 }
      case SVS : {
		   if (isherm()) {
		     this->SetDiv(new HermBandSVDiv<T>(*this,false)); 
		   } else {
		     this->SetDiv(new SymBandSVDiv<T>(*this,false,false)); 
		   }
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
      int i, int j, int istep, int jstep, size_t n) const 
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
      int i, int j, int istep, int jstep, size_t n) const 
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

  template <class T> RealType(T) GenSymBandMatrix<T>::NormSq() const
  { 
    RealType(T) ans = diag().NormSq();
    if (size() > 0 && nlo() > 0) 
      ans += RealType(T)(2) * UpperBandOff().NormSq();
    return ans;
  }

  template <class T> inline RealType(T) NonLapNorm1(
      const GenSymBandMatrix<T>& m) 
  { 
    if (m.nlo() > 0) {
      RealType(T) max(0);
      if (m.size() > 0) {
	size_t i1=0;
	size_t i2=m.nlo()+1;
	size_t k=m.nlo();
	for(size_t j=0;j<m.size();++j) {
	  RealType(T) temp = m.col(j,i1,j).Norm1();
	  temp += m.col(j,j,i2).Norm1();
	  if (temp > max) max = temp;
	  if (k>0) --k; else ++i1;
	  if (i2 < m.size()) ++i2;
	}
      }
      return max;
    } else {
      return m.diag().NormInf();
    }
  } 

  template <class T> inline RealType(T) NonLapNormF(
      const GenSymBandMatrix<T>& m)
  { return SQRT(m.NormSq()); }

  template <class T> inline RealType(T) NonLapMaxAbsElement(
      const GenSymBandMatrix<T>& m)
  { return m.UpperBand().MaxAbsElement(); }

#ifdef XLAP
  template <class T> inline RealType(T) LapNorm(const char c, 
      const GenSymBandMatrix<T>& m)
  {
    switch(c) {
      case 'M' : return NonLapMaxAbsElement(m);
      case '1' : return NonLapNorm1(m);
      case 'F' : return NonLapNormF(m);
      default : TMVAssert(FALSE);
    }
    return RealType(T)(0);
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
      return LAPNAMEX(dlansb) (LAPCM LAPV(cc),
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
	return LAPNAMEX(zlanhb) (LAPCM LAPV(cc),
	    (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO,
	    LAPV(N),LAPV(noff),LAPP(mp),LAPV(lda) LAPWK(work.get()) 
	    LAP1 LAP1);
      else 
	return LAPNAMEX(zlansb) (LAPCM LAPV(cc),
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
      return LAPNAMEX(slansb) (LAPCM LAPV(cc),
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
	return LAPNAMEX(clanhb) (LAPCM LAPV(cc),
	    (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO,
	    LAPV(N),LAPV(noff),LAPP(mp),LAPV(lda) LAPWK(work.get()) 
	    LAP1 LAP1);
      else 
	return LAPNAMEX(clansb) (LAPCM LAPV(cc),
	    (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO,
	    LAPV(N),LAPV(noff),LAPP(mp),LAPV(lda) LAPWK(work.get()) 
	    LAP1 LAP1);
    }
  }
#endif
#endif // XLAP

  template <class T> RealType(T) GenSymBandMatrix<T>::MaxAbsElement() const
  {
#ifdef XLAP
    if ((isrm() && stepi()>0) || (iscm() && stepj()>0))
      return LapNorm('M',*this);
    else
#endif
      return NonLapMaxAbsElement(*this);
  }
  template <class T> RealType(T) GenSymBandMatrix<T>::Norm1() const
  {
#ifdef XLAP
    if ((isrm() && stepi()>0) || (iscm() && stepj()>0))
      return LapNorm('1',*this);
    else
#endif
      return NonLapNorm1(*this);
  }
  template <class T> RealType(T) GenSymBandMatrix<T>::NormF() const
  {
#ifdef XLAP
    if ((isrm() && stepi()>0) || (iscm() && stepj()>0))
      return LapNorm('F',*this);
    else
#endif
      return NonLapNormF(*this);
  }

  template <UpLoType U, class T> 
    SymBandMatrix<T,U,DiagMajor> SymTriDiagMatrix(
	const GenVector<T>& v1, const GenVector<T>& v2)
    {
      TMVAssert2(v2.size() == v1.size()-1);
      SymBandMatrix<T,U,DiagMajor> temp(v1.size(),1);
      temp.diag() = v1;
      temp.diag(1) = v2;
      return temp;
    }

  template <UpLoType U, class T> 
    HermBandMatrix<T,U,DiagMajor> HermTriDiagMatrix(
	const GenVector<T>& v1, const GenVector<T>& v2)
    {
      TMVAssert2(v2.size() == v1.size()-1);
      HermBandMatrix<T,U,DiagMajor> temp(v1.size(),1);
      temp.diag() = v1;
      if (U == Upper)
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
      return ConstSymBandMatrixView<T>(m,s,nlo,2*nlo,1,2*nlo+1,
	  Sym,uplo,RowMajor,NonConj);
    else
      return ConstSymBandMatrixView<T>(m,s,nlo,1,2*nlo,2*nlo+1,
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
      return SymBandMatrixView<T>(m,s,nlo,2*nlo,1,2*nlo+1,
	  Sym,uplo,RowMajor,NonConj
	  FIRSTLAST1(m,m+BandStorageLength(RowMajor,s,s,nlo,nlo)));
    else
      return SymBandMatrixView<T>(m,s,nlo,1,2*nlo,2*nlo+1,
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
      return ConstSymBandMatrixView<T>(m,s,nlo,2*nlo,1,2*nlo+1,
	  Herm,uplo,RowMajor,NonConj);
    else
      return ConstSymBandMatrixView<T>(m,s,nlo,1,2*nlo,2*nlo+1,
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
      return SymBandMatrixView<T>(m,s,nlo,2*nlo,1,2*nlo+1,
	  Herm,uplo,RowMajor,NonConj
	  FIRSTLAST1(m,m+BandStorageLength(RowMajor,s,s,nlo,nlo)));
    else
      return SymBandMatrixView<T>(m,s,nlo,1,2*nlo,2*nlo+1,
	  Herm,uplo,ColMajor,NonConj
	  FIRSTLAST1(m,m+BandStorageLength(ColMajor,s,s,nlo,nlo)));
  }


  //
  // I/O
  //

  template <bool sym, bool conj, bool rm, bool cm, bool compact, bool th, class T>
    inline void DoWrite(std::ostream& os, const GenSymBandMatrix<T>& m, 
	RealType(T) thresh)
    {
      TMVAssert(m.uplo() == Lower);
      size_t gaplen1 = 0;
      size_t rowlen1 = compact ? 1 : 0;
      size_t rowlen2 = m.nlo()+1;
      size_t gaplen2 = m.size()-rowlen2;
      const int si = m.stepi();
      const int sj = rm?1:m.stepj();
      const int sd = m.diagstep();
      const T* mrowi = m.cptr();
      size_t k=m.nlo();

      if (compact) 
	os << (sym?"sB ":"hB ")<< m.size()<<' '<<m.nlo() << std::endl;
      else
	os << m.size() <<' '<< m.size() << std::endl;

      while (rowlen2>0) {
	os << "( ";
	if (!compact) {
	  for(size_t j=gaplen1;j>0;--j) os << ' '<<T(0)<<' ';
	}

	const T* mij = mrowi;
	for(size_t j=rowlen1;j>0;--j,rm?++mij:mij+=sj) {
	  if (conj)
	    if (th)
	      os << ' '<<(ABS(*mij)<thresh ? T(0) : CONJ(*mij))<<' ';
	    else
	      os << ' '<<CONJ(*mij)<<' ';
	  else
	    if (th)
	      os << ' '<<(ABS(*mij)<thresh ? T(0) : *mij)<<' ';
	    else
	      os << ' '<<*mij<<' ';
	}

	if (!compact) {
	  for(size_t j=rowlen2;j>0;--j,cm?++mij:mij+=si) {
	    if (sym == conj)
	      if (th)
		os << ' '<<(ABS(*mij)<thresh ? T(0) : CONJ(*mij))<<' ';
	      else
		os << ' '<<CONJ(*mij)<<' ';
	    else
	      if (th)
		os << ' '<<(ABS(*mij)<thresh ? T(0) : *mij)<<' ';
	      else
		os << ' '<<*mij<<' ';
	  }

	  for(size_t j=gaplen2;j>0;--j) os << ' '<<T(0)<<' ';
	}
	os << " )\n";

	if (k>0) { --k; ++rowlen1; mrowi += si; }
	else { if(!compact) ++gaplen1; mrowi += sd; }
	if (gaplen2 > 0) --gaplen2; 
	else --rowlen2;
      }
    }

  template <bool rm, bool cm, bool compact, bool th, class T> inline void DoWrite1(
      std::ostream& os, const GenSymBandMatrix<T>& m, T thresh)
  { DoWrite<true,false,rm,cm,compact,th>(os,m,thresh); }

  template <bool rm, bool cm, bool compact, bool th, class T> inline void DoWrite1(
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
	DoWrite1<true,false,false,false>(os,*this,RealType(T)(0));
      else if (iscm()) 
	DoWrite1<false,true,false,false>(os,*this,RealType(T)(0));
      else
	DoWrite1<false,false,false,false>(os,*this,RealType(T)(0));
    }
  }

  template <class T> void GenSymBandMatrix<T>::Write(std::ostream& os,
      RealType(T) thresh) const
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
	DoWrite1<true,false,true,false>(os,*this,RealType(T)(0));
      else if (iscm())
	DoWrite1<false,true,true,false>(os,*this,RealType(T)(0));
      else
	DoWrite1<false,false,true,false>(os,*this,RealType(T)(0));
  }

  template <class T> void GenSymBandMatrix<T>::WriteCompact(std::ostream& os,
      RealType(T) thresh) const
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

  template <class T> class SymBandMatrixReadError :
    public ReadError
  {
    public :
      size_t i,j;
      mutable auto_ptr<SymBandMatrix<T> > m;
      std::string exp,got;
      size_t s;
      int lo;
      bool is, iseof, isbad;

      inline SymBandMatrixReadError(
	  size_t _i, size_t _j, const GenSymBandMatrix<T>& _m,
	  std::istream& _is) throw() :
	ReadError("SymBandMatrix"),
	i(_i), j(_j), m(new SymBandMatrix<T>(_m)), exp(""), got(""), 
	s(_m.size()), lo(_m.nlo()),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline SymBandMatrixReadError(std::istream& _is) throw() :
	ReadError("SymBandMatrix"),
	i(0), j(0), m(0), exp(""), got(""), s(0), lo(0),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline SymBandMatrixReadError(
	  size_t _i, size_t _j, const GenSymBandMatrix<T>& _m,
	  std::istream& _is, std::string _e, std::string _g) throw() :
	ReadError("SymBandMatrix"),
	i(_i), j(_j), m(new SymBandMatrix<T>(_m)), 
	exp(_e), got(_g), s(_m.size()), lo(_m.nlo()),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline SymBandMatrixReadError(std::istream& _is, 
	  std::string _e, std::string _g) throw() :
	ReadError("SymBandMatrix"),
	i(0), j(0), m(0), exp(_e), got(_g), s(0), lo(0),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline SymBandMatrixReadError(
	  const GenSymBandMatrix<T>& _m, std::istream& _is, 
	  size_t _s, int _lo) throw() :
	ReadError("SymBandMatrix"),
	i(0), j(0), m(new SymBandMatrix<T>(_m)), exp(""), got(""), 
	s(_s), lo(_lo),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

      inline SymBandMatrixReadError(const SymBandMatrixReadError<T>& rhs) :
	ReadError("SymBandMatrix"),
	i(rhs.i), j(rhs.j), m(rhs.m), exp(rhs.exp), got(rhs.got), 
	s(rhs.s), lo(rhs.lo),
	is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}

      virtual inline ~SymBandMatrixReadError() throw() {}

      virtual void Write(std::ostream& os) const throw();
  };

  template <class T> void SymBandMatrixReadError<T>::Write(
      std::ostream& os) const throw()
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
      size_t rowlen = 1;
      size_t rowstart = 0;
      for(size_t ii=0;ii<i;++ii) {
	os<<"( ";
	for(size_t jj=0;jj<rowlen;++jj)
	  os<<' '<<mm(ii,jj+rowstart)<<' ';
	os<<" )\n";
	if (k > 0) { ++rowlen; --k; }
	else ++rowstart;
      }
      os<<"( ";
      for(size_t jj=0;jj<j;++jj)
	os<<' '<<mm(i,jj+rowstart)<<' ';      
      os<<" )\n";
    }
  }

  template <class T> class HermBandMatrixReadError :
    public ReadError
  {
    public :
      size_t i,j;
      mutable auto_ptr<HermBandMatrix<T> > m;
      std::string exp,got;
      size_t s;
      int lo;
      T dv;
      bool is, iseof, isbad;

      inline HermBandMatrixReadError(
	  size_t _i, size_t _j, const GenSymBandMatrix<T>& _m,
	  std::istream& _is) throw() :
	ReadError("HermBandMatrix"),
	i(_i), j(_j), m(new HermBandMatrix<T>(_m)), exp(""), got(""), 
	s(_m.size()), lo(_m.nlo()),
	dv(0), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline HermBandMatrixReadError(std::istream& _is) throw() :
	ReadError("HermBandMatrix"),
	i(0), j(0), m(0), exp(""), got(""), s(0), lo(0),
	dv(0), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline HermBandMatrixReadError(
	  size_t _i, size_t _j, const GenSymBandMatrix<T>& _m,
	  std::istream& _is, std::string _e, std::string _g) throw() :
	ReadError("HermBandMatrix"),
	i(_i), j(_j), m(new HermBandMatrix<T>(_m)), exp(_e), got(_g), 
	s(_m.size()), lo(_m.nlo()),
	dv(0), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline HermBandMatrixReadError(
	  size_t _i, size_t _j, const GenSymBandMatrix<T>& _m,
	  std::istream& _is, T _dv) throw() :
	ReadError("HermBandMatrix"),
	i(_i), j(_j), m(new HermBandMatrix<T>(_m)), exp(""), got(""), 
	s(_m.size()), lo(_m.nlo()),
	dv(_dv), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline HermBandMatrixReadError(std::istream& _is, 
	  std::string _e, std::string _g) throw() :
	ReadError("HermBandMatrix"),
	i(0), j(0), m(0), exp(_e), got(_g), s(0), lo(0),
	dv(0), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline HermBandMatrixReadError(
	  const GenSymBandMatrix<T>& _m, std::istream& _is, 
	  size_t _s, int _lo) throw() :
	ReadError("HermBandMatrix"),
	i(0), j(0), m(new HermBandMatrix<T>(_m)), exp(""), got(""), 
	s(_s), lo(_lo),
	dv(0), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

      inline HermBandMatrixReadError(const HermBandMatrixReadError<T>& rhs) :
	ReadError("HermBandMatrix"),
	i(rhs.i), j(rhs.j), m(rhs.m), exp(rhs.exp), got(rhs.got),
	s(rhs.s), lo(rhs.lo),
	dv(rhs.dv), is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}

      virtual inline ~HermBandMatrixReadError() throw() {}

      virtual void Write(std::ostream& os) const throw();
  };

  template <class T> void HermBandMatrixReadError<T>::Write(
      std::ostream& os) const throw()
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
      for(size_t ii=0;ii<i;++ii) {
	os<<"( ";
	for(size_t jj=0;jj<(ii<j?i+1:i);++jj)
	  os<<' '<<mm(ii,jj)<<' ';
	os<<" )\n";
      }
      os<<"( ";
      for(size_t jj=0;jj<j;++jj)
	os<<' '<<mm(i,jj)<<' ';      
      os<<" )\n";
    }
  }

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
      size_t rowlen = 1;
      int k = nlo();
      const int si = stepi();
      const int sj = stepj();
      const int sd = si+sj;
      for(size_t i=size();i>0;--i) {
#ifdef NOSTL
	char p;
	is >> p;
	paren = p;
#else
	is >> paren;
#endif
	if (!is || paren != "(") {
	  if (this->issym())
	    throw SymBandMatrixReadError<T>(size()-i,0,*this,is,
		"(",is?paren:"(");
	  else
	    throw HermBandMatrixReadError<T>(size()-i,0,*this,is,
		"(",is?paren:"(");
	}
	T* mij = mrowi;
	if (this->isrm()) {
	  for(size_t j=rowlen;j>0;--j,++mij) {
	    is >> *mij;
	    if (!is)
	      if (this->issym())
		throw SymBandMatrixReadError<T>(size()-i,rowlen-j,*this,is);
	      else
		throw HermBandMatrixReadError<T>(size()-i,rowlen-j,*this,is);
	    if (this->isherm() && j==1 && IMAG(*mij) != RealType(T)(0))
	      throw HermBandMatrixReadError<T>(size()-i,rowlen-j,*this,is,*mij);
	  }
	}
	else {
	  for(size_t j=rowlen;j>0;--j,mij+=sj) {
	    is >> *mij;
	    if (!is)
	      if (this->issym())
		throw SymBandMatrixReadError<T>(size()-i,rowlen-j,*this,is);
	      else
		throw HermBandMatrixReadError<T>(size()-i,rowlen-j,*this,is);
	    if (this->isherm() && j==1 && IMAG(*mij) != RealType(T)(0))
	      throw HermBandMatrixReadError<T>(size()-i,rowlen-j,*this,is,*mij);
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
	    throw SymBandMatrixReadError<T>(size()-i,size()-i+1,*this,is,
		")",is?paren:")");
	  else
	    throw HermBandMatrixReadError<T>(size()-i,size()-i+1,*this,is,
		")",is?paren:")");
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
	throw SymBandMatrixReadError<T>(is);
      if (IsReal(T())) {
	if (sh != "sB" && sh != "hB") 
	  throw SymBandMatrixReadError<T>(is,"sB",sh);
      } else {
	if (sh != "sB") 
	  throw SymBandMatrixReadError<T>(is,"sB",sh);
      }
      size_t size;
      int nlo;
      is >> size >> nlo;
      if (!is) 
	throw SymBandMatrixReadError<T>(is);
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
	throw HermBandMatrixReadError<T>(is);
      if (IsReal(T())) {
	if (sh != "sB" && sh != "hB") 
	  throw HermBandMatrixReadError<T>(is,"hB",sh);
      } else {
	if (sh != "hB") 
	  throw HermBandMatrixReadError<T>(is,"hB",sh);
      }
      size_t size;
      int nlo;
      is >> size >> nlo;
      if (!is) 
	throw HermBandMatrixReadError<T>(is);
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
    if (!is)
      if (m.issym())
	throw SymBandMatrixReadError<T>(is);
      else
	throw HermBandMatrixReadError<T>(is);
    if (IsReal(T())) {
      if (sh != "sB" && sh != "hB") 
	if (m.issym())
	  throw SymBandMatrixReadError<T>(is,"sB",sh);
	else
	  throw HermBandMatrixReadError<T>(is,"hB",sh);
    } else if (m.issym()) {
      if (sh != "sB") 
	throw SymBandMatrixReadError<T>(is,"sB",sh);
    } else {
      if (sh != "hB") 
	throw HermBandMatrixReadError<T>(is,"hB",sh);
    }
    size_t s;
    int nlo;
    is >> s >> nlo;
    if (!is) 
      if (m.issym())
	throw SymBandMatrixReadError<T>(is);
      else
	throw HermBandMatrixReadError<T>(is);
    if (s != m.size() || nlo != m.nlo())
      if (m.issym())
	throw SymBandMatrixReadError<T>(m,is,s,nlo);
      else
	throw HermBandMatrixReadError<T>(m,is,s,nlo);
    TMVAssert(m.size() == s);
    TMVAssert(m.nlo() == nlo);
    m.Read(is);
    return is;
  }


#define InstFile "TMV_SymBandMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


