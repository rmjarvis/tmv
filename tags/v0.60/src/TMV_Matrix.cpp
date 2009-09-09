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
#include "TMV_Matrix.h"
#include "TMV_Divider.h"
#include "TMV_LUDiv.h"
#include "TMV_QRDiv.h"
#include "TMV_QRPDiv.h"
#include "TMV_SVDiv.h"
#include "TMV_MatrixArith.h"
#include "TMV_VIt.h"
#include <iostream>

namespace tmv {

#ifdef TMV_BLOCKSIZE
  const size_t PERM_BLOCKSIZE = TMV_BLOCKSIZE/2;
#else
  const size_t PERM_BLOCKSIZE = 32;
#endif

  //
  // Constructors
  //

#define NEW_SIZE(cs,rs) \
  linsize((cs)*(rs)), \
  itsm(new T[linsize]), itscs(cs), itsrs(rs) \
  DEFFIRSTLAST(itsm.get(),itsm.get()+ls())

  template <class T, StorageType S, IndexStyle I> Matrix<T,S,I>::Matrix(
      const std::vector<std::vector<T> >& vv) :
    NEW_SIZE(vv.size(),(vv.size()>0?vv[0].size():0))
    {
      TMVAssert(S==RowMajor || S==ColMajor);
      T* vi=itsm.get();
      if (S == RowMajor) {
	for(size_t i=0;i<colsize();++i) {
	  TMVAssert(vv[i].size() == rowsize());
	  typename std::vector<T>::const_iterator vvi = vv[i].begin();
	  for(size_t j=0;j<rowsize();++j,++vi,++vvi) 
	    *vi = *vvi;
	}
      } else { 
	for(size_t j=0;j<rowsize();++j) {
	  for(size_t i=0;i<colsize();++i,++vi) {
	    TMVAssert(vv[i].size() == rowsize());
	    *vi = vv[i][j];
	  }
	}
      }
    }

#undef NEW_SIZE


  //
  // Access
  //

  template <class T> T GenMatrix<T>::cref(size_t i, size_t j) const
  {
    TMVAssert(i<colsize() && j<rowsize());
    const T* mi;
    if (isrm()) mi = cptr() + int(i)*stepi() + j;
    else if (iscm()) mi = cptr() + i + int(j)*stepj();
    else mi = cptr() + int(i)*stepi() + int(j)*stepj();
    return isconj() ? CONJ(*mi) : *mi;
  }

  template <class T, IndexStyle I> RefType(T) MatrixView<T,I>::ref(
      size_t i, size_t j) const
  {
    TMVAssert(i<colsize());
    TMVAssert(j<rowsize());
    T* mi;
    if (this->isrm()) mi = ptr() + int(i)*stepi() + int(j);
    else if (this->iscm()) mi = ptr() + i + int(j)*stepj();
    else mi = ptr() + int(i)*itssi + int(j)*stepj();
#ifdef TMVFLDEBUG
    TMVAssert(mi >= first);
    TMVAssert(mi < last);
#endif
    return REF(mi,ct());
  }

  template <class T> void GenMatrix<T>::NewDivider() const
  {
    switch (this->GetDivType()) {
      case tmv::LU : this->SetDiv(new LUDiv<T>(*this,
			   this->IsDivInPlace())); 
		     break;
      case tmv::QR : this->SetDiv(new QRDiv<T>(*this,
			   this->IsDivInPlace())); 
		     break;
      case tmv::QRP : this->SetDiv(new QRPDiv<T>(*this,
			    this->IsDivInPlace())); 
		      break;
      case tmv::SV : this->SetDiv(new SVDiv<T>(*this,
			   this->IsDivInPlace(),true,true)); 
		     break;
      case tmv::SVS : this->SetDiv(new SVDiv<T>(*this,
			    this->IsDivInPlace(),false,false)); 
		      break;
      case tmv::SVU : this->SetDiv(new SVDiv<T>(*this,
			    this->IsDivInPlace(),true,false)); 
		      break;
      case tmv::SVV : this->SetDiv(new SVDiv<T>(*this,
			    this->IsDivInPlace(),false,true)); 
		      break;
      default : std::cerr<<"dt = "<<Text(this->GetDivType())<<std::endl; 
		TMVAssert(FALSE);
    }
  }

  //
  // OK? (SubMatrix, SubVector)
  //

  template <class T> bool GenMatrix<T>::OKSubMatrix(
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
      std::cout<<"row range ("<<j2-j1<<") must be multiple of istep (";
      std::cout<<jstep<<")\n";
    }
    if ((j2-j1)/jstep < 0) {
      ok = false;
      std::cout<<"n row elements ("<<(j2-j1)/jstep<<") must be nonnegative\n";
    }
    return ok;
  }
  
  template <class T> bool GenMatrix<T>::OKSubVector(
      int i, int j, int istep, int jstep, size_t size) const 
  {
    //cerr<<"SubVector "<<tmv::Type(*this)<<" "<<i<<','<<j<<','<<istep<<','<<jstep<<','<<size<<endl;
    if (size==0) return true;
    bool ok = true;
    if (istep == 0 && jstep == 0) {
      ok = false;
      std::cout<<"istep ("<<istep<<") and jstep ("<<jstep;
      std::cout<<") can not both be 0\n";
    }
    if (i < 0 || i >= int(colsize())) {
      ok = false;
      std::cout<<"i ("<<i<<") must be in 0 -- "<<colsize()-1<<std::endl;
    }
    if (j < 0 || j >= int(rowsize())) {
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
    return ok;
  }

  template <class T> bool ConstMatrixView<T,FortranStyle>::OKSubMatrix(
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
    return ok;
  }

  template <class T> bool ConstMatrixView<T,FortranStyle>::OKSubVector(
      int i, int j, int istep, int jstep, size_t size) const 
  {
    //cerr<<"SubVectorF "<<tmv::Type(*this)<<" "<<i<<','<<j<<','<<istep<<','<<jstep<<','<<size<<endl;
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
    return ok;
  }

  //
  // Norms
  //

  template <class T> RealType(T) GenMatrix<T>::NormSq() const
  {
    if (CanLinearize()) return ConstLinearView().NormSq();
    else {
      RealType(T) sum(0);
      if (isrm()) 
	for(size_t i=0;i<colsize();++i) 
	  sum += row(i).NormSq();
      else 
	for(size_t j=0;j<rowsize();++j) 
	  sum += col(j).NormSq();
      return sum;
    }
  }

  template <class T> inline RealType(T) NonLapMaxAbsElement(
      const GenMatrix<T>& m)
  {
    if (m.CanLinearize()) return m.ConstLinearView().MaxAbsElement();
    else {
      RealType(T) max(0);
      if (m.iscm())
	for(size_t j=0;j<m.rowsize();++j) {
	  RealType(T) temp = m.col(j).NormInf();
	  if (temp > max) max = temp;
	}
      else 
	for(size_t i=0;i<m.colsize();++i) {
	  RealType(T) temp = m.row(i).NormInf();
	  if (temp > max) max = temp;
	}
      return max;
    }
  }

  template <class T> inline RealType(T) NonLapNorm1(const GenMatrix<T>& m)
  { 
    RealType(T) max(0);
    for(size_t j=0;j<m.rowsize();++j) {
      RealType(T) temp = m.col(j).Norm1();
      if (temp > max) max = temp;
    }
    return max;
  }

  template <class T> inline RealType(T) NonLapNormF(const GenMatrix<T>& m)
  { return tmv::SQRT(m.NormSq()); }

  template <class T> inline RealType(T) NonLapNormInf(const GenMatrix<T>& m)
  { return NonLapNorm1(m.Transpose()); }

#ifdef XLAP
  template <class T> inline RealType(T) LapNorm(const char c,
      const GenMatrix<T>& m)
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
  template <> inline double LapNorm(const char c, const GenMatrix<double>& m)
  { 
    TMVAssert(m.iscm());
    char cc = c;
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
#ifndef LAPNOWORK
    auto_array<double> work(c == 'I' ? new double[M] : 0);
#endif
    double norm = LAPNAMEX(dlange) (LAPCM LAPV(cc),LAPV(M),LAPV(N),
	LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()));
    return norm;
  }
  template <> inline double LapNorm(const char c, 
      const GenMatrix<std::complex<double> >& m)
  { 
    TMVAssert(m.iscm());
    char cc = c;
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
#ifndef LAPNOWORK
    auto_array<double> work(c == 'I' ? new double[M] : 0);
#endif
    double norm = LAPNAMEX(zlange) (LAPCM LAPV(cc),LAPV(M),LAPV(N),
	LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()));
    return norm;
  }
#endif
#ifdef INST_FLOAT
  template <> inline float LapNorm(const char c, const GenMatrix<float>& m)
  { 
    TMVAssert(m.iscm());
    char cc = c;
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
#ifndef LAPNOWORK
    auto_array<float> work(c == 'I' ? new float[M] : 0);
#endif
    double norm = LAPNAMEX(slange) (LAPCM LAPV(cc),LAPV(M),LAPV(N),
	LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()));
    return norm;
  }
  template <> inline float LapNorm(const char c, 
      const GenMatrix<std::complex<float> >& m)
  { 
    TMVAssert(m.iscm());
    char cc = c;
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
#ifndef LAPNOWORK
    auto_array<float> work(c == 'I' ? new float[M] : 0);
#endif
    double norm = LAPNAMEX(clange) (LAPCM LAPV(cc),LAPV(M),LAPV(N),
	LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()));
    return norm;
  }
#endif
#endif // XLAP

  template <class T> RealType(T) GenMatrix<T>::MaxAbsElement() const
  {
#ifdef XLAP
    if (isrm() && stepi() > 0) return LapNorm('M',Transpose());
    else if (iscm() && stepj() > 0) return LapNorm('M',*this);
    else
#endif
      return NonLapMaxAbsElement(*this);
  }
  template <class T> RealType(T) GenMatrix<T>::Norm1() const
  {
#ifdef XLAP
    if (isrm() && stepi() > 0) return LapNorm('I',Transpose());
    else if (iscm() && stepj() > 0) return LapNorm('1',*this);
    else
#endif
      return NonLapNorm1(*this);
  }
  template <class T> RealType(T) GenMatrix<T>::NormF() const
  {
#ifdef XLAP
    if (isrm() && stepi() > 0) return LapNorm('F',Transpose());
    else if (iscm() && stepj() > 0) return LapNorm('F',*this);
    else
#endif
      return NonLapNormF(*this);
  }

  template <class T> QuotXM<T,T> GenMatrix<T>::QInverse() const
  { return QuotXM<T,T>(T(1),*this); }

  template <class T> auto_ptr<BaseMatrix<T> > GenMatrix<T>::NewCopy() const
  {
    auto_ptr<BaseMatrix<T> > a;
    if (isrm()) a.reset(new Matrix<T,RowMajor>(*this));
    else a.reset(new Matrix<T,ColMajor>(*this));
    return a;
  }

  template <class T> auto_ptr<BaseMatrix<T> > GenMatrix<T>::NewView() const
  {
    auto_ptr<BaseMatrix<T> > a(new ConstMatrixView<T>(View()));
    return a;
  }

  template <class T> auto_ptr<BaseMatrix<T> > GenMatrix<T>::NewTranspose() const
  {
    auto_ptr<BaseMatrix<T> > a(new ConstMatrixView<T>(Transpose()));
    return a;
  }

  template <class T> auto_ptr<BaseMatrix<T> > GenMatrix<T>::NewConjugate() const
  {
    auto_ptr<BaseMatrix<T> > a(new ConstMatrixView<T>(Conjugate()));
    return a;
  }

  template <class T> auto_ptr<BaseMatrix<T> > GenMatrix<T>::NewAdjoint() const
  {
    auto_ptr<BaseMatrix<T> > a(new ConstMatrixView<T>(Adjoint()));
    return a;
  }

  template <class T> auto_ptr<BaseMatrix<T> > GenMatrix<T>::NewInverse() const
  {
    auto_ptr<Matrix<T,ColMajor> > minv(
	new Matrix<T,ColMajor>(rowsize(),colsize()));
    Inverse(minv->View());
    BaseMatrix<T>* ret1 = minv.release();
    auto_ptr<BaseMatrix<T> > ret(ret1);
    return ret;
  }
  
  //
  // Modifying Functions
  //

  template <class T, IndexStyle I> 
    const MatrixView<T,I>& MatrixView<T,I>::Clip(RealType(T) thresh) const
    { 
      TMVAssert(I==CStyle);
      if (this->CanLinearize()) LinearView().Clip(thresh);
      else {
	if (this->isrm())
	  for(size_t i=0;i<colsize();++i) row(i).Clip(thresh);
	else 
	  for(size_t j=0;j<rowsize();++j) col(j).Clip(thresh);
      }
      return *this; 
    }

  template <class T, IndexStyle I> 
    const MatrixView<T,I>& MatrixView<T,I>::SetAllTo(T x) const
    { 
      TMVAssert(I==CStyle);
      if (this->CanLinearize()) LinearView().SetAllTo(x);
      else {
	if (this->isrm())
	  for(size_t i=0;i<colsize();++i) row(i).SetAllTo(x); 
	else 
	  for(size_t j=0;j<rowsize();++j) col(j).SetAllTo(x); 
      }
      return *this; 
    }

  template <class T, IndexStyle I> 
    const MatrixView<T,I>& MatrixView<T,I>::TransposeSelf() const
    {
      TMVAssert(I==CStyle);
      TMVAssert(colsize() == rowsize());
      for(size_t i=1;i<colsize();++i) Swap(row(i,0,i),col(i,0,i));
      return *this; 
    }

  template <class T, IndexStyle I> 
    const MatrixView<T,I>& MatrixView<T,I>::ConjugateSelf() const
    { 
      TMVAssert(I==CStyle);
      if (IsComplex(T())) {
	if (this->CanLinearize()) LinearView().ConjugateSelf();
	else {
	  if (this->isrm()) 
	    for(size_t i=0;i<colsize();++i) row(i).ConjugateSelf();
	  else 
	    for(size_t j=0;j<rowsize();++j) col(j).ConjugateSelf();
	}
      }
      return *this; 
    }

  template <class T, IndexStyle I> 
    const MatrixView<T,I>& MatrixView<T,I>::SetToIdentity(T x) const 
    {
      TMVAssert(I==CStyle);
      TMVAssert(colsize() == rowsize());
      Zero();
      diag().SetAllTo(x);
      return *this;
    }

  template <class T, IndexStyle I> 
    const MatrixView<T,I>& MatrixView<T,I>::PermuteRows(
	const size_t* p, size_t i1, size_t i2) const
    { 
      TMVAssert(I==CStyle);
      TMVAssert(i2<=colsize());
      TMVAssert(i1<=i2);
      // This idea of doing the permutations a block at a time
      // is cribbed from the LAPack code.  It does speed things up 
      // quite a bit for large matrices.  On my machine where BLOCKSIZE=64
      // is optimal for most routines, blocks of 32 were optimal here,
      // so I use BLOCKSIZE/2 in general.
      const size_t Nx = rowsize()/PERM_BLOCKSIZE*PERM_BLOCKSIZE;
      if (Nx != 0) {
	for(size_t j=0;j<Nx;) {
	  size_t j2 = j+PERM_BLOCKSIZE;
	  const size_t* pi = p+i1;
	  for(size_t i=i1;i<i2;++i,++pi) {
            TMVAssert(*pi < colsize());
	    Cols(j,j2).SwapRows(i,*pi);
	  }
	  j = j2;
	}
      }
      if (Nx != rowsize()) {
	const size_t* pi = p+i1;
	for(size_t i=i1;i<i2;++i,++pi) {
	  TMVAssert(*pi < colsize());
	  Cols(Nx,rowsize()).SwapRows(i,*pi);
	}
      }
      return *this;
    }

  template <class T, IndexStyle I> 
    const MatrixView<T,I>& MatrixView<T,I>::ReversePermuteRows(
	const size_t* p, size_t i1, size_t i2) const
    { 
      TMVAssert(I==CStyle);
      TMVAssert(i2<=colsize());
      TMVAssert(i1<=i2);
      const size_t Nx = rowsize()/PERM_BLOCKSIZE*PERM_BLOCKSIZE;
      if (Nx != 0) {
	for(size_t j=0;j<Nx;) {
	  size_t j2 = j+PERM_BLOCKSIZE;
	  const size_t* pi = p+i2;
	  for(size_t i=i2;i>i1;) {
	    --i; --pi;
	    TMVAssert(*pi < colsize());
	    Cols(j,j2).SwapRows(i,*pi);
	  }
	  j = j2;
	}
      }
      if (Nx != rowsize()) {
	const size_t* pi = p+i2;
	for(size_t i=i2;i>i1;) {
	  --i; --pi;
	  TMVAssert(*pi < colsize());
	  Cols(Nx,rowsize()).SwapRows(i,*pi);
	}
      }
      return *this;
    }

  //
  // Copy Matrices
  //

  template <bool c1, class T> inline void NonLapCopy(
      const GenMatrix<T>& m1, const MatrixView<T>& m2)
  {
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert(m2.ct()==NonConj);
    TMVAssert(!(m1.isrm() && m2.isrm()));
    TMVAssert(c1 == m1.isconj());
    if (m1.iscm() && m2.iscm()) 
      if (c1)
	for(size_t j=0;j<m2.rowsize();++j) 
	  DoCopySameType<c1>(m1.col(j),m2.col(j));
      else
	for(size_t j=0;j<m2.rowsize();++j) 
	  memmove(m2.col(j).ptr(),m1.col(j).cptr(),m2.colsize()*sizeof(T));
    else if (m2.colsize() > m2.rowsize())
      if (ShouldReverse(m1.stepi(),m2.stepi()))
	for(size_t j=0;j<m2.rowsize();++j) 
	  DoCopySameType<c1>(m1.col(j).Reverse(),m2.col(j).Reverse());
      else
	for(size_t j=0;j<m2.rowsize();++j) 
	  DoCopySameType<c1>(m1.col(j),m2.col(j));
    else 
      if (ShouldReverse(m1.stepj(),m2.stepj()))
	for(size_t i=0;i<m2.colsize();++i) 
	  DoCopySameType<c1>(m1.row(i).Reverse(),m2.row(i).Reverse());
      else
	for(size_t i=0;i<m2.colsize();++i) 
	  DoCopySameType<c1>(m1.row(i),m2.row(i));
  }
#ifdef ELAP
  template <class T> inline void LapCopy(
      const GenMatrix<T>& m1, const MatrixView<T>& m2)
  { NonLapCopy<false>(m1,m2); }
#ifdef INST_DOUBLE
  template <> inline void LapCopy(
      const GenMatrix<double>& m1, const MatrixView<double>& m2)
  { 
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert(m1.iscm());
    TMVAssert(m2.iscm());
    char c = 'A';
    int m = m1.colsize();
    int n = m1.rowsize();
    int ld1 = m1.stepj();
    int ld2 = m2.stepj();
    LAPNAMEX(dlacpy) (LAPCM LAPV(c),LAPV(m),LAPV(n),LAPP(m1.cptr()),LAPV(ld1),
	LAPP(m2.ptr()),LAPV(ld2));
  }
  template <> inline void LapCopy(
      const GenMatrix<std::complex<double> >& m1,
      const MatrixView<std::complex<double> >& m2)
  { 
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert(m1.iscm());
    TMVAssert(m2.iscm());
    TMVAssert(m1.ct() == NonConj);
    TMVAssert(m2.ct() == NonConj);
    char c = 'A';
    int m = m1.colsize();
    int n = m1.rowsize();
    int ld1 = m1.stepj();
    int ld2 = m2.stepj();
    LAPNAMEX(zlacpy) (LAPCM LAPV(c),LAPV(m),LAPV(n),LAPP(m1.cptr()),LAPV(ld1),
	LAPP(m2.ptr()),LAPV(ld2));
  }
#endif
#ifdef INST_FLOAT
  template <> inline void LapCopy(
      const GenMatrix<float>& m1, const MatrixView<float>& m2)
  { 
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert(m1.iscm());
    TMVAssert(m2.iscm());
    char c = 'A';
    int m = m1.colsize();
    int n = m1.rowsize();
    int ld1 = m1.stepj();
    int ld2 = m2.stepj();
    LAPNAMEX(slacpy) (LAPCM LAPV(c),LAPV(m),LAPV(n),LAPP(m1.cptr()),LAPV(ld1),
	LAPP(m2.ptr()),LAPV(ld2));
  }
  template <> inline void LapCopy(
      const GenMatrix<std::complex<float> >& m1,
      const MatrixView<std::complex<float> >& m2)
  { 
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert(m1.iscm());
    TMVAssert(m2.iscm());
    TMVAssert(m1.ct() == NonConj);
    TMVAssert(m2.ct() == NonConj);
    char c = 'A';
    int m = m1.colsize();
    int n = m1.rowsize();
    int ld1 = m1.stepj();
    int ld2 = m2.stepj();
    LAPNAMEX(clacpy) (LAPCM LAPV(c),LAPV(m),LAPV(n),LAPP(m1.cptr()),LAPV(ld1),
	LAPP(m2.ptr()),LAPV(ld2));
  }
#endif
#endif
  template <bool c1, class T> void DoCopySameType(const GenMatrix<T>& m1, 
      const MatrixView<T>& m2)
  {
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert(m2.rowsize() > 0);
    TMVAssert(m2.colsize() > 0);
    TMVAssert(m2.ct() == NonConj);
    TMVAssert(!(m1.isrm() && m2.isrm()));
    TMVAssert(!m2.SameAs(m1));
    TMVAssert(c1 == m1.isconj());

#ifdef ELAP
    if (!c1 && m1.iscm() && m2.iscm() && m1.stepj()>0 && m2.stepj()>0) 
      LapCopy(m1,m2);
    else
#endif
      NonLapCopy<c1>(m1,m2);
  }

  // 
  // Swap
  //

  template <class T> void Swap(
      const MatrixView<T>& m1, const MatrixView<T>& m2)
  {
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    if (m1.stor() == m2.stor() && m1.CanLinearize() && m2.CanLinearize()) {
      TMVAssert(m1.LinearView().size() == m2.LinearView().size());
      TMVAssert(m1.stepi()==m2.stepi() && m1.stepj()==m2.stepj());
      Swap(m1.LinearView(),m2.LinearView());
    }
    else {
      if (m1.isrm() && m2.isrm())
	for(size_t i=0;i<m1.colsize();++i) Swap(m1.row(i),m2.row(i)); 
      else 
	for(size_t j=0;j<m1.rowsize();++j) Swap(m1.col(j),m2.col(j)); 
    }
  }

  //
  // m1 == m2
  //

  template <class T1, class T2> bool operator==(
      const GenMatrix<T1>& m1, const GenMatrix<T2>& m2)
  {
    if (m1.colsize() != m2.colsize()) return false;
    if (m1.rowsize() != m2.rowsize()) return false;
    if (m1.SameAs(m2)) return true;
    if (m1.stepi()==m2.stepi() && m1.stepj()==m2.stepj() &&
	m1.CanLinearize() && m2.CanLinearize())
      return m1.ConstLinearView() == m2.ConstLinearView();
    else {
      for(size_t i=0;i<m1.colsize();++i) 
	if (m1.row(i) != m2.row(i)) return false;
      return true;  
    }
  }

  //
  // I/O
  //

  template <bool conj, bool rm, bool th, class T> 
    inline void DoWrite(std::ostream& os, const GenMatrix<T>& m,
	RealType(T) thresh)
  {
    const T* mrowi = m.cptr();
    const int sj = rm ? 1 : m.stepj();
    os << m.colsize() <<"  "<<m.rowsize()<<std::endl;
    for(size_t i=m.colsize();i>0;--i,mrowi+=m.stepi()) {
      os << "( ";
      const T* mij = mrowi;
      for(size_t k=m.rowsize();k>0;--k,rm?++mij:mij+=sj) 
	if (conj) 
	  if (th)
	    os << ' '<<(ABS(*mij) < thresh ? T(0) : CONJ(*mij))<<' ';
	  else
	    os << ' '<<CONJ(*mij)<<' ';
	else 
	  if (th)
	    os << ' '<<(ABS(*mij) < thresh ? T(0) : *mij)<<' ';
	  else
	    os << ' '<<*mij<<' ';
      os << " )\n";
    }
  }

  template <bool rm, bool th, class T> inline void DoWrite1(
      std::ostream& os, const GenMatrix<T>& m, T thresh)
  { DoWrite<false,rm,th>(os,m,thresh); }

  template <bool rm, bool th, class T> inline void DoWrite1(
      std::ostream& os, const GenMatrix<std::complex<T> >& m, T thresh)
  {
    if (m.isconj())
      DoWrite<true,rm,th>(os,m,thresh);
    else
      DoWrite<false,rm,th>(os,m,thresh);
  }

  template <class T> void GenMatrix<T>::Write(std::ostream& os) const
  {
    if (isrm())
      DoWrite1<true,false>(os,*this,RealType(T)(0));
    else
      DoWrite1<false,false>(os,*this,RealType(T)(0));
  }

  template <class T> void GenMatrix<T>::Write(std::ostream& os,
      RealType(T) thresh) const
  {
    if (isrm())
      DoWrite1<true,true>(os,*this,thresh);
    else
      DoWrite1<false,true>(os,*this,thresh);
  }

  template <class T> class MatrixReadError :
    public ReadError
  {
    public :
      size_t i,j;
      mutable auto_ptr<Matrix<T> > m;
      char exp,got;
      size_t cs,rs;
      bool is,iseof,isbad;

      inline MatrixReadError(size_t _i, size_t _j, const GenMatrix<T>& _m,
	  std::istream& _is) throw() :
	ReadError("Matrix"),
	i(_i), j(_j), m(new Matrix<T>(_m)), exp(0), got(0),
	cs(_m.colsize()), rs(_m.rowsize()),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline MatrixReadError(std::istream& _is) throw() :
	ReadError("Matrix"),
	i(0), j(0), m(0), exp(0), got(0), cs(0), rs(0),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline MatrixReadError(size_t _i, size_t _j, const GenMatrix<T>& _m,
	  std::istream& _is, char _e, char _g) throw() :
	ReadError("Matrix"),
	i(_i), j(_j), m(new Matrix<T>(_m)), exp(_e), got(_g),
	cs(_m.colsize()), rs(_m.rowsize()),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline MatrixReadError(const GenMatrix<T>& _m,
	  std::istream& _is, size_t _cs, size_t _rs) throw() :
	ReadError("Matrix"),
	i(0), j(0), m(new Matrix<T>(_m)), exp(0), got(0), cs(_cs), rs(_rs),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

      inline MatrixReadError(const MatrixReadError<T>& rhs) :
	i(rhs.i), j(rhs.j), m(rhs.m), exp(rhs.exp), got(rhs.got),
	cs(rhs.cs), rs(rhs.rs),
	is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
      virtual inline ~MatrixReadError() throw() {}

      virtual void Write(std::ostream& os) const throw();
  };

  template <class T> void MatrixReadError<T>::Write(
      std::ostream& os) const throw()
  {
    os<<"TMV Read Error: Reading istream input for Matrix\n";
    if (exp != got) {
      os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
    }
    if (m.get() && cs != m->colsize()) {
      os<<"Wrong column size: expected "<<m->colsize()<<", got "<<cs<<".\n";
    }
    if (m.get() && rs != m->rowsize()) {
      os<<"Wrong row size: expected "<<m->rowsize()<<", got "<<rs<<".\n";
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
      os<<"The portion of the Matrix which was successfully read is: \n";
      for(size_t ii=0;ii<i;++ii) {
	os<<"( ";
	for(size_t jj=0;jj<m->rowsize();++jj)
	  os<<' '<<(*m)(ii,jj)<<' ';
	os<<" )\n";
      }
      os<<"( ";
      for(size_t jj=0;jj<j;++jj)
	os<<' '<<(*m)(i,jj)<<' ';
      os<<" )\n";
    }
  }

  template <class T, IndexStyle I> void MatrixView<T,I>::Read(
      std::istream& is) const
  {
    TMVAssert(I==CStyle);
    T* mrowi = ptr();
    const int sj = stepj();
    char paren;
    for(size_t i=0;i<colsize();++i,mrowi+=stepi()) {
      is >> paren;
      if (!is || paren != '(') 
	throw MatrixReadError<T>(i,0,*this,is,'(',is?paren:'(');
      T* mij = mrowi;
      if (this->isrm()) 
	for(size_t k=rowsize();k>0;--k,++mij) {
	  is >> *mij;
	  if (!is) 
	    throw MatrixReadError<T>(i,rowsize()-k,*this,is);
	}
      else
	for(size_t k=rowsize();k>0;--k,mij+=sj) {
	  is >> *mij;
	  if (!is) 
	    throw MatrixReadError<T>(i,rowsize()-k,*this,is);
	}
      is >> paren;
      if ((!is && i+1<colsize())  || paren != ')') 
	throw MatrixReadError<T>(i,rowsize(),*this,is,')',is?paren:')');
    }
    if (this->isconj()) ConjugateSelf();
  }

  template <class T, StorageType S, IndexStyle I> std::istream& operator>>(
      std::istream& is, auto_ptr<Matrix<T,S,I> >& m)
  { 
    size_t cs,rs;
    is >> cs >> rs;
    if (!is) 
      throw MatrixReadError<T>(is);
    m.reset(new Matrix<T,S,I>(cs,rs));
    m->View().Read(is); 
    return is;
  }

  template <class T> std::istream& operator>>(std::istream& is, 
      const MatrixView<T>& m)
  { 
    size_t cs,rs;
    is >> cs >> rs;
    if (!is) 
      throw MatrixReadError<T>(is);
    if (m.colsize() != cs || m.rowsize() != rs) 
      throw MatrixReadError<T>(m,is,cs,rs);
    TMVAssert(m.colsize() == cs);
    TMVAssert(m.rowsize() == rs);
    m.Read(is);
    return is;
  }


#define InstFile "TMV_Matrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


