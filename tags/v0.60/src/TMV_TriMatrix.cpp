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
#include "TMV_TriMatrix.h"
#include <iostream>

namespace tmv {

  //
  // Access
  //

  template <class T> T GenUpperTriMatrix<T>::cref(size_t i, size_t j) const
  {
    TMVAssert(i < size());
    TMVAssert(j < size());
    if (i==j && isunit()) return T(1);
    else if (i>j) return T(0);
    else {
      const T* mi;
      if (isrm()) mi = cptr() + int(i)*stepi() + j;
      else if (iscm()) mi = cptr() + i + int(j)*stepj();
      else mi = cptr() + int(i)*stepi() + int(j)*stepj();
      return (isconj() ? CONJ(*mi) : *mi);
    }
  }

  template <class T, IndexStyle I> RefType(T) UpperTriMatrixView<T,I>::ref(
      size_t i, size_t j) const
  {
    TMVAssert(i < size());
    TMVAssert(j < size());
    TMVAssert(this->okij(i,j));
    T* mi;
    if (isrm()) mi = ptr() + int(i)*stepi() + j;
    else if (iscm()) mi = ptr() + i + int(j)*stepj();
    else mi = ptr() + int(i)*stepi() + int(j)*stepj();
#ifdef TMVFLDEBUG
    TMVAssert(mi>=first);
    TMVAssert(mi<last);
#endif
    return REF(mi,ct());
  }

  template <class T> T GenLowerTriMatrix<T>::cref(size_t i, size_t j) const
  {
    TMVAssert(i < size());
    TMVAssert(j < size());
    if (i==j && isunit()) return T(1);
    else if (i<j) return T(0);
    else {
      const T* mi;
      if (isrm()) mi = cptr() + int(i)*stepi() + j;
      else if (iscm()) mi = cptr() + i + int(j)*stepj();
      else mi = cptr() + int(i)*stepi() + int(j)*stepj();
      return (isconj() ? CONJ(*mi) : *mi);
    }
  }

  template <class T, IndexStyle I> RefType(T) LowerTriMatrixView<T,I>::ref(
      size_t i, size_t j) const
  {
    TMVAssert(i < size());
    TMVAssert(j < size());
    TMVAssert(this->okij(i,j));
    T* mi;
    if (isrm()) mi = ptr() + int(i)*stepi() + j;
    else if (iscm()) mi = ptr() + i + int(j)*stepj();
    else mi = ptr() + int(i)*stepi() + int(j)*stepj();
#ifdef TMVFLDEBUG
    TMVAssert(mi>=first);
    TMVAssert(mi<last);
#endif
    return REF(mi,ct());
  }

  //
  // OK? (SubMatrix, etc.)
  //

  template <class T> bool GenUpperTriMatrix<T>::OKSubMatrix(
      int i1, int i2, int j1, int j2, int istep, int jstep) const
  { 
    if (i1==i2 || j1==j2) return true; // no elements, so whatever...
    bool ok = true;
    if (istep == 0) {
      ok = false;
      std::cout<<"istep ("<<istep<<") can not be 0\n";
    }
    if (i1 < 0 || i1 >= int(size())) {
      ok = false;
      std::cout<<"first col element ("<<i1<<") must be in 0 -- ";
      std::cout<<size()-1<<std::endl;
    }
    if (i2-istep < 0 || i2-istep >= int(size())) {
      ok = false;
      std::cout<<"last col element ("<<i2-istep<<") must be in 0 -- ";
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
      std::cout<<"row range ("<<j2-j1<<") must be multiple of istep (";
      std::cout<<jstep<<")\n";
    }
    if ((j2-j1)/jstep < 0) {
      ok = false;
      std::cout<<"n row elements ("<<(j2-j1)/jstep<<") must be nonnegative\n";
    }
    if (!this->okij(i1,j1)) {
      ok = false;
      std::cout<<"Upper left corner ("<<i1<<','<<j1;
      std::cout<<") must be in Upper Triangle\n";
    }
    if (!this->okij(i1,j2-jstep)) {
      ok = false;
      std::cout<<"Upper right corner ("<<i1<<','<<j2-jstep;
      std::cout<<") must be in Upper Triangle\n";
    }
    if (!this->okij(i2-istep,j1)) {
      ok = false;
      std::cout<<"Lower left corner ("<<i2-istep<<','<<j1;
      std::cout<<") must be in Upper Triangle\n";
    }
    if (!this->okij(i2-istep,j2-jstep)) {
      ok = false;
      std::cout<<"Lower right corner ("<<i2-istep<<','<<j2-jstep;
      std::cout<<") must be in Upper Triangle\n";
    }
    return ok;
  }

  template <class T> bool GenUpperTriMatrix<T>::OKSubVector(
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
    if (!this->okij(i,j)) {
      ok = false;
      std::cout<<"First element ("<<i<<','<<j<<") must be in Triangle\n";
    }
    if (!this->okij(i2,j2)) {
      ok = false;
      std::cout<<"Last element ("<<i2<<','<<j2<<") must be in Triangle\n";
    }
    return ok;
  }

  template <class T> bool GenUpperTriMatrix<T>::OKSubTriMatrix(
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
    return ok;
  }

  template <class T> bool ConstUpperTriMatrixView<T,FortranStyle>::OKSubMatrix(
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
      std::cout<<"last col element ("<<i2<<") must be in 1 -- ";
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
      std::cout<<"last row element ("<<j2<<") must be in 1 -- ";
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
    if (!this->okij(i1-1,j1-1)) {
      ok = false;
      std::cout<<"Upper left corner ("<<i1<<','<<j1;
      std::cout<<") must be in Upper Triangle\n";
    }
    if (!this->okij(i1-1,j2-1)) {
      ok = false;
      std::cout<<"Upper right corner ("<<i1<<','<<j2;
      std::cout<<") must be in Upper Triangle\n";
    }
    if (!this->okij(i2-1,j1-1)) {
      ok = false;
      std::cout<<"Lower left corner ("<<i2<<','<<j1;
      std::cout<<") must be in Upper Triangle\n";
    }
    if (!this->okij(i2-1,j2-1)) {
      ok = false;
      std::cout<<"Lower right corner ("<<i2<<','<<j2;
      std::cout<<") must be in Upper Triangle\n";
    }
    return ok;
  }

  template <class T> bool ConstUpperTriMatrixView<T,FortranStyle>::OKSubVector(
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
    if (!this->okij(i-1,j-1)) {
      ok = false;
      std::cout<<"First element ("<<i<<','<<j<<") must be in Triangle\n";
    }
    if (!this->okij(i2-1,j2-1)) {
      ok = false;
      std::cout<<"Last element ("<<i2<<','<<j2<<") must be in Triangle\n";
    }
    return ok;
  }

  template <class T> 
    bool ConstUpperTriMatrixView<T,FortranStyle>::OKSubTriMatrix(
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
      if (i2<1 || i2 > int(this->size())) {
	ok = false;
	std::cout<<"last diag element ("<<i2<<") must be in 1 -- ";
	std::cout<<this->size()<<std::endl;
      }
      if ((i2-i1)%istep != 0) {
	ok = false;
	std::cout<<"range ("<<i2-i1<<") must be multiple of istep (";
	std::cout<<istep<<")\n";
      }
      if ((i2-i1)/istep < 0) {
	ok = false;
	std::cout<<"n diag elements ("<<(i2-i1)/istep+1;
	std::cout<<") must be positive\n";
      }
      return ok;
    }

  //
  // Norms
  //

  template <class T> RealType(T) GenUpperTriMatrix<T>::NormSq() const
  {
    RealType(T) sum(0);
    if (isrm()) 
      if (isunit())
	for(size_t i=0;i<size();++i) 
	  sum += row(i,i+1,size()).NormSq();
      else
	for(size_t i=0;i<size();++i) 
	  sum += row(i,i,size()).NormSq();
    else
      if (isunit())
	for(size_t j=0;j<size();++j) 
	  sum += col(j,0,j).NormSq();
      else
	for(size_t j=0;j<size();++j) 
	  sum += col(j,0,j+1).NormSq();
    if (isunit()) sum += size();
    return sum;
  }

  template <class T> inline RealType(T) NonLapNormF(
      const GenUpperTriMatrix<T>& m)
  { return SQRT(m.NormSq()); }

  template <class T> inline RealType(T) NonLapMaxAbsElement(
      const GenUpperTriMatrix<T>& m)
  {
    RealType(T) max(0);
    if (m.isrm())
      for(size_t i=0;i<m.size();++i) {
	RealType(T) temp;
	if (m.isunit())
	  temp = m.row(i,i+1,m.size()).NormInf();
	else 
	  temp = m.row(i,i,m.size()).NormInf();
	if (temp > max) max = temp;
      }
    else
      for(size_t j=0;j<m.size();++j) {
	RealType(T) temp;
	if (m.isunit())
	  temp = m.col(j,0,j).NormInf();
	else 
	  temp = m.col(j,0,j+1).NormInf();
	if (temp > max) max = temp;
      }
    if (m.isunit() && max < RealType(T)(1)) max = RealType(T)(1);
    return max;
  }

  template <class T> inline RealType(T) NonLapNorm1(
      const GenUpperTriMatrix<T>& m)
  { 
    RealType(T) max(0);
    for(size_t j=0;j<m.size();++j) {
      RealType(T) temp;
      if (m.isunit()) {
	temp = m.col(j,0,j).Norm1();
	temp += RealType(T)(1);
      } else temp = m.col(j,0,j+1).Norm1();
      if (temp > max) max = temp;
    }
    return max;
  } 

  template <class T> inline RealType(T) NonLapNormInf(
      const GenUpperTriMatrix<T>& m)
  { 
    RealType(T) max(0);
    for(size_t j=0;j<m.size();++j) {
      RealType(T) temp;
      if (m.isunit()) {
	temp = m.row(j,j+1,m.size()).Norm1();
	temp += RealType(T)(1);
      } else temp = m.row(j,j,m.size()).Norm1();
      if (temp > max) max = temp;
    }
    return max;
  }
#ifdef XLAP
  template <class T> inline RealType(T) LapNorm(
      const char c, const GenUpperTriMatrix<T>& m)
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
      const char c, const GenUpperTriMatrix<double>& m)
  {
    TMVAssert(m.iscm() || m.isrm());
    char cc = c;
    if (m.isrm()) {
      if (c == '1') cc = 'I';
      else if (c == 'I') cc = '1';
    }
    int M = m.size();
    int N = M;
    int lda = m.iscm() ? m.stepj() : m.stepi();
#ifndef LAPNOWORK
    auto_array<double> work(cc == 'I' ? new double[M] : 0);
#endif
    double norm = LAPNAMEX(dlantr) (LAPCM LAPV(cc),
	m.iscm() ? LAPCH_UP : LAPCH_LO, m.isunit() ? LAPCH_U : LAPCH_NU,
	LAPV(M),LAPV(N),LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
	LAP1 LAP1 LAP1);
    return norm;
  }
  template <> inline double LapNorm(
      const char c, const GenUpperTriMatrix<std::complex<double> >& m)
  {
    TMVAssert(m.iscm() || m.isrm());
    char cc = c;
    if (m.isrm()) {
      if (c == '1') cc = 'I';
      else if (c == 'I') cc = '1';
    }
    int M = m.size();
    int N = M;
    int lda = m.iscm() ? m.stepj() : m.stepi();
#ifndef LAPNOWORK
    auto_array<double> work(cc == 'I' ? new double[M] : 0);
#endif
    double norm = LAPNAMEX(zlantr) (LAPCM LAPV(cc),
	m.iscm() ? LAPCH_UP : LAPCH_LO, m.isunit() ? LAPCH_U : LAPCH_NU,
	LAPV(M),LAPV(N),LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
	LAP1 LAP1 LAP1);
    return norm;
  }
#endif
#ifdef INST_FLOAT
  template <> inline float LapNorm(
      const char c, const GenUpperTriMatrix<float>& m)
  {
    TMVAssert(m.iscm() || m.isrm());
    char cc = c;
    if (m.isrm()) {
      if (c == '1') cc = 'I';
      else if (c == 'I') cc = '1';
    }
    int M = m.size();
    int N = M;
    int lda = m.iscm() ? m.stepj() : m.stepi();
#ifndef LAPNOWORK
    auto_array<float> work(cc == 'I' ? new float[M] : 0);
#endif
    float norm = LAPNAMEX(slantr) (LAPCM LAPV(cc),
	m.iscm() ? LAPCH_UP : LAPCH_LO, m.isunit() ? LAPCH_U : LAPCH_NU,
	LAPV(M),LAPV(N),LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
	LAP1 LAP1 LAP1);
    return norm;
  }
  template <> inline float LapNorm(
      const char c, const GenUpperTriMatrix<std::complex<float> >& m)
  {
    TMVAssert(m.iscm() || m.isrm());
    char cc = c;
    if (m.isrm()) {
      if (c == '1') cc = 'I';
      else if (c == 'I') cc = '1';
    }
    int M = m.size();
    int N = M;
    int lda = m.iscm() ? m.stepj() : m.stepi();
#ifndef LAPNOWORK
    auto_array<float> work(cc == 'I' ? new float[M] : 0);
#endif
    float norm = LAPNAMEX(clantr) (LAPCM LAPV(cc),
	m.iscm() ? LAPCH_UP : LAPCH_LO, m.isunit() ? LAPCH_U : LAPCH_NU,
	LAPV(M),LAPV(N),LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
	LAP1 LAP1 LAP1);
    return norm;
  }
#endif
#endif // XLAP

  template <class T> RealType(T) GenUpperTriMatrix<T>::MaxAbsElement() const
  {
#ifdef XLAP
    return LapNorm('M',*this);
#else
    return NonLapMaxAbsElement(*this);
#endif
  }
  template <class T> RealType(T) GenUpperTriMatrix<T>::Norm1() const
  {
#ifdef XLAP
    return LapNorm('1',*this);
#else
    return NonLapNorm1(*this);
#endif
  }
  template <class T> RealType(T) GenUpperTriMatrix<T>::NormInf() const
  {
#ifdef XLAP
    return LapNorm('I',*this);
#else
    return NonLapNormInf(*this);
#endif
  }
  template <class T> RealType(T) GenUpperTriMatrix<T>::NormF() const
  {
#ifdef XLAP
    return LapNorm('F',*this);
#else
    return NonLapNormF(*this);
#endif
  }

  template <class T> 
    auto_ptr<BaseMatrix<T> > GenUpperTriMatrix<T>::NewCopy() const
  {
    auto_ptr<BaseMatrix<T> > a;
    if (isunit()) {
      if (isrm()) a.reset(new UpperTriMatrix<T,UnitDiag,RowMajor>(*this));
      else a.reset(new UpperTriMatrix<T,UnitDiag,ColMajor>(*this));
    } else {
      if (isrm()) a.reset(
	  new UpperTriMatrix<T,NonUnitDiag,RowMajor>(*this));
      else a.reset(new UpperTriMatrix<T,NonUnitDiag,ColMajor>(*this));
    }
    return a;
  }

  template <class T> 
    auto_ptr<BaseMatrix<T> > GenUpperTriMatrix<T>::NewView() const
  {
    auto_ptr<BaseMatrix<T> > a(new ConstUpperTriMatrixView<T>(View()));
    return a;
  }

  template <class T> 
    auto_ptr<BaseMatrix<T> > GenUpperTriMatrix<T>::NewTranspose() const
  {
    auto_ptr<BaseMatrix<T> > a(
	new ConstLowerTriMatrixView<T>(Transpose()));
    return a;
  }

  template <class T> 
    auto_ptr<BaseMatrix<T> > GenUpperTriMatrix<T>::NewConjugate() const
  {
    auto_ptr<BaseMatrix<T> > a(
	new ConstUpperTriMatrixView<T>(Conjugate()));
    return a;
  }

  template <class T> 
    auto_ptr<BaseMatrix<T> > GenUpperTriMatrix<T>::NewAdjoint() const
  {
    auto_ptr<BaseMatrix<T> > a(new ConstLowerTriMatrixView<T>(Adjoint()));
    return a;
  }

  template <class T> 
    auto_ptr<BaseMatrix<T> > GenUpperTriMatrix<T>::NewInverse() const
    {
      if (isunit()) {
	if (isrm()) {
	  auto_ptr<UpperTriMatrix<T,UnitDiag,ColMajor> > minv(
	      new UpperTriMatrix<T,UnitDiag,ColMajor>(*this));
	  minv->InvertSelf();
	  BaseMatrix<T>* ret1 = minv.release();
	  auto_ptr<BaseMatrix<T> > ret(ret1);
	  return ret;
	} else {
	  auto_ptr<UpperTriMatrix<T,UnitDiag,ColMajor> > minv(
	      new UpperTriMatrix<T,UnitDiag,ColMajor>(*this));
	  minv->InvertSelf();
	  BaseMatrix<T>* ret1 = minv.release();
	  auto_ptr<BaseMatrix<T> > ret(ret1);
	  return ret;
	}
      } else {
	if (isrm()) {
	  auto_ptr<UpperTriMatrix<T,NonUnitDiag,ColMajor> > minv(
	      new UpperTriMatrix<T,NonUnitDiag,ColMajor>(*this));
	  minv->InvertSelf();
	  BaseMatrix<T>* ret1 = minv.release();
	  auto_ptr<BaseMatrix<T> > ret(ret1);
	  return ret;
	} else {
	  auto_ptr<UpperTriMatrix<T,NonUnitDiag,ColMajor> > minv(
	      new UpperTriMatrix<T,NonUnitDiag,ColMajor>(*this));
	  minv->InvertSelf();
	  BaseMatrix<T>* ret1 = minv.release();
	  auto_ptr<BaseMatrix<T> > ret(ret1);
	  return ret;
	}
      }
    }

  template <class T> 
    auto_ptr<BaseMatrix<T> > GenLowerTriMatrix<T>::NewCopy() const
  {
    auto_ptr<BaseMatrix<T> > a;
    if (isunit()) {
      if (isrm()) a.reset(new LowerTriMatrix<T,UnitDiag,RowMajor>(*this));
      else a.reset(new LowerTriMatrix<T,UnitDiag,ColMajor>(*this));
    } else {
      if (isrm()) a.reset(
	  new LowerTriMatrix<T,NonUnitDiag,RowMajor>(*this));
      else a.reset(new LowerTriMatrix<T,NonUnitDiag,ColMajor>(*this));
    }
    return a;
  }

  template <class T> 
    auto_ptr<BaseMatrix<T> > GenLowerTriMatrix<T>::NewView() const
  {
    auto_ptr<BaseMatrix<T> > a(new ConstLowerTriMatrixView<T>(View()));
    return a;
  }

  template <class T> 
    auto_ptr<BaseMatrix<T> > GenLowerTriMatrix<T>::NewTranspose() const
  {
    auto_ptr<BaseMatrix<T> > a(
	new ConstUpperTriMatrixView<T>(Transpose()));
    return a;
  }

  template <class T> 
    auto_ptr<BaseMatrix<T> > GenLowerTriMatrix<T>::NewConjugate() const
  {
    auto_ptr<BaseMatrix<T> > a(
	new ConstLowerTriMatrixView<T>(Conjugate()));
    return a;
  }

  template <class T> 
    auto_ptr<BaseMatrix<T> > GenLowerTriMatrix<T>::NewAdjoint() const
  {
    auto_ptr<BaseMatrix<T> > a(new ConstUpperTriMatrixView<T>(Adjoint()));
    return a;
  }

  template <class T> 
    auto_ptr<BaseMatrix<T> > GenLowerTriMatrix<T>::NewInverse() const
    {
      if (isunit()) {
	if (isrm()) {
	  auto_ptr<LowerTriMatrix<T,UnitDiag,ColMajor> > minv(
	      new LowerTriMatrix<T,UnitDiag,ColMajor>(*this));
	  minv->InvertSelf();
	  BaseMatrix<T>* ret1 = minv.release();
	  auto_ptr<BaseMatrix<T> > ret(ret1);
	  return ret;
	} else {
	  auto_ptr<LowerTriMatrix<T,UnitDiag,ColMajor> > minv(
	      new LowerTriMatrix<T,UnitDiag,ColMajor>(*this));
	  minv->InvertSelf();
	  BaseMatrix<T>* ret1 = minv.release();
	  auto_ptr<BaseMatrix<T> > ret(ret1);
	  return ret;
	}
      } else {
	if (isrm()) {
	  auto_ptr<LowerTriMatrix<T,NonUnitDiag,ColMajor> > minv(
	      new LowerTriMatrix<T,NonUnitDiag,ColMajor>(*this));
	  minv->InvertSelf();
	  BaseMatrix<T>* ret1 = minv.release();
	  auto_ptr<BaseMatrix<T> > ret(ret1);
	  return ret;
	} else {
	  auto_ptr<LowerTriMatrix<T,NonUnitDiag,ColMajor> > minv(
	      new LowerTriMatrix<T,NonUnitDiag,ColMajor>(*this));
	  minv->InvertSelf();
	  BaseMatrix<T>* ret1 = minv.release();
	  auto_ptr<BaseMatrix<T> > ret(ret1);
	  return ret;
	}
      }
    }


  //
  // Modifying Functions
  //

  template <class T, IndexStyle I> 
    const UpperTriMatrixView<T,I>& UpperTriMatrixView<T,I>::SetAllTo(T x) const
    { 
      const size_t N = size();

      if (isrm())
	if (isunit())
	  for(size_t i=0;i<N;++i) row(i,i+1,N).SetAllTo(x); 
	else
	  for(size_t i=0;i<N;++i) row(i,i,N).SetAllTo(x); 
      else 
	if (isunit())
	  for(size_t j=0;j<N;++j) col(j,0,j).SetAllTo(x); 
	else
	  for(size_t j=0;j<N;++j) col(j,0,j+1).SetAllTo(x); 
      return *this; 
    }

  template <class T, IndexStyle I> 
    const UpperTriMatrixView<T,I>& UpperTriMatrixView<T,I>::Clip(
	RealType(T) thresh) const
    { 
      const size_t N = size();

      if (isrm())
	if (isunit())
	  for(size_t i=0;i<N;++i) row(i,i+1,N).Clip(thresh);
	else
	  for(size_t i=0;i<N;++i) row(i,i,N).Clip(thresh);
      else 
	if (isunit())
	  for(size_t j=0;j<N;++j) col(j,0,j).Clip(thresh);
	else
	  for(size_t j=0;j<N;++j) col(j,0,j+1).Clip(thresh);
      return *this; 
    }

  template <class T, IndexStyle I> 
    const UpperTriMatrixView<T,I>& UpperTriMatrixView<T,I>::ConjugateSelf() const
    { 
      const size_t N = size();

      if (IsComplex(T())) {
	if (isrm())
	  if (isunit())
	    for(size_t i=0;i<N;++i) row(i,i+1,N).ConjugateSelf();
	  else
	    for(size_t i=0;i<N;++i) row(i,i,N).ConjugateSelf();
	else
	  if (isunit())
	    for(size_t j=0;j<N;++j) col(j,0,j).ConjugateSelf();
	  else
	    for(size_t j=0;j<N;++j) col(j,0,j+1).ConjugateSelf();
      }
      return *this; 
    }

  template <class T, IndexStyle I> 
    const UpperTriMatrixView<T,I>& UpperTriMatrixView<T,I>::SetToIdentity(
	T x) const 
    {
      TMVAssert(!isunit() || x==T(1));
      Zero();
      if (!isunit()) diag().SetAllTo(x);
      return *this;
    }

  //
  // Swap
  //

  template <class T> void Swap(
      const UpperTriMatrixView<T>& m1, const UpperTriMatrixView<T>& m2)
  {
    TMVAssert(m1.size() == m2.size());
    TMVAssert(m1.dt() == m2.dt());
    const size_t N = m1.size();

    if (m1.isrm() && m2.isrm())
      if (m1.isunit()) 
	for(size_t i=0;i<N;++i) 
	  Swap(m1.row(i,i+1,N), m2.row(i,i+1,N));
      else
	for(size_t i=0;i<N;++i) 
	  Swap(m1.row(i,i,N), m2.row(i,i,N));
    else
      if (m1.isunit()) 
	for(size_t j=0;j<N;++j) 
	  Swap(m1.col(j,0,j), m2.col(j,0,j));
      else
	for(size_t j=0;j<N;++j) 
	  Swap(m1.col(j,0,j+1), m2.col(j,0,j+1));
  }

  //
  // m1 == m2
  //

  template <class T1, class T2> bool operator==(
      const GenUpperTriMatrix<T1>& m1, const GenUpperTriMatrix<T2>& m2)
  {
    if (m1.size() != m2.size()) return false;
    if (m1.SameAs(m2)) return true;
    for(size_t j=0;j<m1.size();++j) {
      if (m1.col(j,0,j) != m2.col(j,0,j)) return false;
    }

    if (m1.isunit() && !m2.isunit()) {
      for(size_t i=0;i<m1.size();++i) if (m2(i,i) != T2(1)) return false;
    }
    else if (m2.isunit() && !m1.isunit()) {
      for(size_t i=0;i<m1.size();++i) if (m1(i,i) != T1(1)) return false;
    }
    else if (!m1.isunit() && !m2.isunit()) {
      if (m1.diag() != m2.diag()) return false;
    }

    return true;  
  }

  //
  // I/O
  //

  template <bool conj, bool rm, bool compact, bool th, class T> 
    inline void DoWrite(std::ostream& os, const GenUpperTriMatrix<T>& m,
	RealType(T) thresh)
  {
    const T* mrowi = m.cptr();
    const int sj = rm?1:m.stepj();
    const int ds = m.stepi()+sj;
    size_t len = m.size();

    if (m.isunit()) {
      mrowi += sj;
      --len;
    }

    if (compact)
      os << "U " << m.size() << std::endl;
    else
      os << m.size() <<' '<< m.size() << std::endl;

    for(size_t i=0;i<m.size();++i,--len,mrowi+=ds) {
      os << "( ";
      if (!compact) {
	for(size_t j=0;j<i;j++) os << ' '<<T(0)<<' ';
      }
      if (m.isunit()) os << ' '<<T(1)<<' ';
      const T* mij = mrowi;
      for(size_t k=len;k>0;--k,rm?++mij:mij+=sj) 
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

  template <bool rm, bool compact, bool th, class T> inline void DoWrite1(
      std::ostream& os, const GenUpperTriMatrix<T>& m, T thresh)
  { DoWrite<false,rm,compact,th>(os,m,thresh); }

  template <bool rm, bool compact, bool th, class T> inline void DoWrite1(
      std::ostream& os, const GenUpperTriMatrix<std::complex<T> >& m,
      T thresh)
  { 
    if (m.isconj())
      DoWrite<true,rm,compact,th>(os,m,thresh); 
    else
      DoWrite<false,rm,compact,th>(os,m,thresh); 
  }

  template <class T> void GenUpperTriMatrix<T>::Write(std::ostream& os) const
  {
    if (isrm())
      DoWrite1<true,false,false>(os,*this,RealType(T)(0));
    else
      DoWrite1<false,false,false>(os,*this,RealType(T)(0));
  }

  template <class T> void GenUpperTriMatrix<T>::Write(std::ostream& os,
      RealType(T) thresh) const
  {
    if (isrm())
      DoWrite1<true,false,true>(os,*this,thresh);
    else
      DoWrite1<false,false,true>(os,*this,thresh);
  }

  template <class T> void GenUpperTriMatrix<T>::WriteCompact(
      std::ostream& os) const
  {
    if (isrm())
      DoWrite1<true,true,false>(os,*this,RealType(T)(0));
    else
      DoWrite1<false,true,false>(os,*this,RealType(T)(0));
  }

  template <class T> void GenUpperTriMatrix<T>::WriteCompact(std::ostream& os,
      RealType(T) thresh) const
  {
    if (isrm())
      DoWrite1<true,true,true>(os,*this,thresh);
    else
      DoWrite1<false,true,true>(os,*this,thresh);
  }

  template <class T> class UpperTriMatrixReadError :
    public ReadError
  {
    public :
      size_t i,j;
      mutable auto_ptr<UpperTriMatrix<T> > m;
      char exp,got;
      T unitgot;
      size_t s;
      bool is, iseof, isbad;

      inline UpperTriMatrixReadError(size_t _i, size_t _j,
	  const GenUpperTriMatrix<T>& _m, std::istream& _is) throw() :
	ReadError("UpperTriMatrix"),
	i(_i), j(_j), m(new UpperTriMatrix<T>(_m)),
	exp(0), got(0), unitgot(T(1)), s(_m.size()),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline UpperTriMatrixReadError(std::istream& _is) throw() :
	ReadError("UpperTriMatrix"),
	i(0), j(0), m(0), exp(0), got(0), unitgot(T(1)), s(0),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline UpperTriMatrixReadError(size_t _i, size_t _j,
	  const GenUpperTriMatrix<T>& _m, std::istream& _is,
	  char _e, char _g) throw() :
	ReadError("UpperTriMatrix"),
	i(_i), j(_j), m(new UpperTriMatrix<T>(_m)),
	exp(_e), got(_g), unitgot(T(1)), s(_m.size()),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline UpperTriMatrixReadError(std::istream& _is, 
	  char _e, char _g) throw() :
	ReadError("UpperTriMatrix"),
	i(0), j(0), m(0), exp(_e), got(_g), unitgot(T(1)), s(0),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline UpperTriMatrixReadError(size_t _i, size_t _j,
	  const GenUpperTriMatrix<T>& _m, std::istream& _is, T _u) throw() :
	ReadError("UpperTriMatrix"),
	i(_i), j(_j), m(new UpperTriMatrix<T>(_m)),
	exp(0), got(0), unitgot(_u), s(_m.size()),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline UpperTriMatrixReadError(const GenUpperTriMatrix<T>& _m,
	  std::istream& _is, size_t _s) throw() :
	ReadError("UpperTriMatrix"),
	i(0), j(0), m(new UpperTriMatrix<T>(_m)),
	exp(0), got(0), unitgot(T(1)), s(_s),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

      inline UpperTriMatrixReadError(const UpperTriMatrixReadError<T>& rhs) :
	i(rhs.i), j(rhs.j), m(rhs.m), exp(rhs.exp), got(rhs.got),
	unitgot(rhs.unitgot), s(rhs.s),
	is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
      virtual inline ~UpperTriMatrixReadError() throw() {}

      virtual void Write(std::ostream& os) const throw();
  };

  template <class T> void UpperTriMatrixReadError<T>::Write(
      std::ostream& os) const throw()
  {
    os<<"TMV Read Error: Reading istream input for UpperTriMatrix\n";
    if (exp != got) {
      os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
    }
    if (unitgot != T(1)) {
      os<<"Wrong format: expected 1 on the diagonal, got "<<unitgot<<".\n";
    }
    if (m.get() && s != m->size()) {
      os<<"Wrong size: expected "<<m->size()<<", got "<<s<<".\n";
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
      os<<"The portion of the UpperTriMatrix which was successfully read is:\n";
      ConstUpperTriMatrixView<T> mm = m->View();
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

  template <class T, IndexStyle I> void UpperTriMatrixView<T,I>::Read(
      std::istream& is) const
  {
    T* mrowi = ptr();
    const int sj = stepj();
    const int ds = stepi()+sj;
    size_t len = size();
    if (isunit()) {
      mrowi += sj;
      --len;
    }
    char paren;
    for(size_t i=0;i<size();++i,--len,mrowi+=ds) {
      is >> paren;
      if (!is || paren != '(') 
	throw UpperTriMatrixReadError<T>(i,0,*this,is,'(',is?paren:'(');
      if (isunit()) {
	T unit;
	is >> unit;
	if (!is || unit != T(1))
	  throw UpperTriMatrixReadError<T>(i,i,*this,is,is?unit:T(1));
      }
      T* mij = mrowi;
      if (isrm()) 
	for(size_t k=len;k>0;--k,++mij) {
	  is >> *mij;
	  if (!is) 
	    throw UpperTriMatrixReadError<T>(i,size()-k,*this,is);
	}
      else 
	for(size_t k=len;k>0;--k,mij+=sj) {
	  is >> *mij;
	  if (!is) 
	    throw UpperTriMatrixReadError<T>(i,size()-k,*this,is);
	}
      is >> paren;
      if ((!is && i+1<size())  || paren != ')') 
	throw UpperTriMatrixReadError<T>(i,size(),*this,is,')',is?paren:')');
    }
    if (isconj()) ConjugateSelf();
  }

  template <class T, DiagType D, StorageType S, IndexStyle I> 
    std::istream& operator>>(std::istream& is, 
	auto_ptr<UpperTriMatrix<T,D,S,I> >& m)
    { 
      char ul;
      is >> ul;
      if (!is)
	throw UpperTriMatrixReadError<T>(is);
      if (ul != 'U')
	throw UpperTriMatrixReadError<T>(is,'U',ul);
      size_t size;
      is >> size;
      if (!is) 
	throw UpperTriMatrixReadError<T>(is);
      m.reset(new UpperTriMatrix<T,D,S,I>(size));
      m->View().Read(is); 
      return is;
    }

  template <class T> std::istream& operator>>(
      std::istream& is, const UpperTriMatrixView<T>& m)
  { 
    char ul;
    is >> ul;
    if (!is)
      throw UpperTriMatrixReadError<T>(is);
    if (ul != 'U')
      throw UpperTriMatrixReadError<T>(is,'U',ul);
    size_t s;
    is >> s;
    if (!is)
      throw UpperTriMatrixReadError<T>(is);
    if (m.size() != s)
      throw UpperTriMatrixReadError<T>(m,is,s);
    TMVAssert(m.size() == s);
    m.Read(is);
    return is;
  }

  template <bool conj, bool rm, bool compact, bool th, class T> 
    inline void DoWrite(std::ostream& os, const GenLowerTriMatrix<T>& m,
	RealType(T) thresh)
  {
    const T* mrowi = m.cptr();
    const int sj = rm?1:m.stepj();
    const int si = m.stepi();
    size_t len = m.isunit() ? 0 : 1;

    if (compact)
      os << "L " << m.size() << std::endl;
    else
      os << m.size()<<' '<<m.size() << std::endl;

    for(size_t i=m.size();i>0;--i,++len,mrowi+=si) {
      os << "( ";
      const T* mij = mrowi;
      for(size_t k=len;k>0;--k,rm?++mij:mij+=sj) 
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

      if (m.isunit()) os << ' '<<T(1)<<' ';
      if (!compact)
	for(size_t j=m.isunit()?len+1:len;j<m.size();j++) 
	  os << ' '<<T(0)<<' ';
      os << " )\n";
    }
  }

  template <bool rm, bool compact, bool th, class T> inline void DoWrite1(
      std::ostream& os, const GenLowerTriMatrix<T>& m, T thresh)
  { DoWrite<false,rm,compact,th>(os,m,thresh); }

  template <bool rm, bool compact, bool th, class T> inline void DoWrite1(
      std::ostream& os, const GenLowerTriMatrix<std::complex<T> >& m,
      T thresh)
  { 
    if (m.isconj())
      DoWrite<true,rm,compact,th>(os,m,thresh); 
    else
      DoWrite<false,rm,compact,th>(os,m,thresh); 
  }

  template <class T> void GenLowerTriMatrix<T>::Write(std::ostream& os) const
  {
    if (isrm())
      DoWrite1<true,false,false>(os,*this,RealType(T)(0));
    else
      DoWrite1<false,false,false>(os,*this,RealType(T)(0));
  }

  template <class T> void GenLowerTriMatrix<T>::Write(std::ostream& os,
      RealType(T) thresh) const
  {
    if (isrm())
      DoWrite1<true,false,true>(os,*this,thresh);
    else
      DoWrite1<false,false,true>(os,*this,thresh);
  }

  template <class T> void GenLowerTriMatrix<T>::WriteCompact(
      std::ostream& os) const
  {
    if (isrm())
      DoWrite1<true,true,false>(os,*this,RealType(T)(0));
    else
      DoWrite1<false,true,false>(os,*this,RealType(T)(0));
  }

  template <class T> void GenLowerTriMatrix<T>::WriteCompact(std::ostream& os,
      RealType(T) thresh) const
  {
    if (isrm())
      DoWrite1<true,true,true>(os,*this,thresh);
    else
      DoWrite1<false,true,true>(os,*this,thresh);
  }

  template <class T> class LowerTriMatrixReadError :
    public ReadError
  {
    public :
      size_t i,j;
      mutable auto_ptr<LowerTriMatrix<T> > m;
      char exp,got;
      T unitgot;
      size_t s;
      bool is, iseof, isbad;

      inline LowerTriMatrixReadError(size_t _i, size_t _j,
	  const GenLowerTriMatrix<T>& _m, std::istream& _is) :
	ReadError("LowerTriMatrix"),
	i(_i), j(_j), m(new LowerTriMatrix<T>(_m)),
	exp(0), got(0), unitgot(T(1)), s(_m.size()),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline LowerTriMatrixReadError(std::istream& _is) :
	ReadError("LowerTriMatrix"),
	i(0), j(0), m(0), exp(0), got(0), unitgot(T(1)), s(0),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline LowerTriMatrixReadError(size_t _i, size_t _j,
	  const GenLowerTriMatrix<T>& _m, std::istream& _is,
	  char _e, char _g) :
	ReadError("LowerTriMatrix"),
	i(_i), j(_j), m(new LowerTriMatrix<T>(_m)),
	exp(_e), got(_g), unitgot(T(1)), s(_m.size()),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline LowerTriMatrixReadError(std::istream& _is, char _e, char _g) :
	ReadError("LowerTriMatrix"),
	i(0), j(0), m(0), exp(_e), got(_g), unitgot(T(1)), s(0),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline LowerTriMatrixReadError(size_t _i, size_t _j,
	  const GenLowerTriMatrix<T>& _m, std::istream& _is, T _u) :
	ReadError("LowerTriMatrix"),
	i(_i), j(_j), m(new LowerTriMatrix<T>(_m)),
	exp(0), got(0), unitgot(_u), s(_m.size()),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
      inline LowerTriMatrixReadError(const GenLowerTriMatrix<T>& _m,
	  std::istream& _is, size_t _s) :
	ReadError("LowerTriMatrix"),
	i(0), j(0), m(new LowerTriMatrix<T>(_m)),
	exp(0), got(0), unitgot(T(1)), s(_s),
	is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

      inline LowerTriMatrixReadError(const LowerTriMatrixReadError<T>& rhs) :
	i(rhs.i), j(rhs.j), m(rhs.m), exp(rhs.exp), got(rhs.got),
	unitgot(rhs.unitgot), s(rhs.s),
	is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
      virtual inline ~LowerTriMatrixReadError() throw() {}

      void Write(std::ostream& os) const throw();
  };
  
  template <class T> void LowerTriMatrixReadError<T>::Write(
      std::ostream& os) const throw()
  {
    os<<"TMV Read Error: Reading istream input for LowerTriMatrix\n";
    if (exp != got) {
      os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
    }
    if (unitgot != T(1)) {
      os<<"Wrong format: expected 1 on the diagonal, got "<<unitgot<<".\n";
    }
    if (m.get() && s != m->size()) {
      os<<"Wrong size: expected "<<m->size()<<", got "<<s<<".\n";
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
      os<<"The portion of the LowerTriMatrix which was successfully read is:\n";
      ConstLowerTriMatrixView<T> mm = m->View();
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

  template <class T, IndexStyle I> void LowerTriMatrixView<T,I>::Read(
      std::istream& is) const
  {
    T* mrowi = ptr();
    const int sj = stepj();
    size_t len = isunit() ? 0 : 1;
    char paren;
    for(size_t i=0;i<size();++i,++len,mrowi+=stepi()) {
      is >> paren;
      if (!is || paren != '(') 
	throw LowerTriMatrixReadError<T>(i,0,*this,is,'(',is?paren:'(');
      T* mij = mrowi;
      if (isrm()) 
	for(size_t k=len;k>0;--k,++mij) {
	  is >> *mij;
	  if (!is) 
	    throw LowerTriMatrixReadError<T>(i,len-k,*this,is);
	}
      else 
	for(size_t k=len;k>0;--k,mij+=sj) {
	  is >> *mij;
	  if (!is) 
	    throw LowerTriMatrixReadError<T>(i,len-k,*this,is);
	}
      if (isunit()) {
	T unit;
	is >> unit;
	if (!is || unit != T(1))
	  throw LowerTriMatrixReadError<T>(i,i,*this,is,is?unit:T(1));
      }
      is >> paren;
      if ((!is && i+1<size())  || paren != ')') 
	throw LowerTriMatrixReadError<T>(i,size(),*this,is,')',is?paren:')');
    }
    if (isconj()) ConjugateSelf();
  }

  template <class T, DiagType D, StorageType S, IndexStyle I>
    std::istream& operator>>(std::istream& is, 
	auto_ptr<LowerTriMatrix<T,D,S,I> >& m)
    { 
      char ul;
      is >> ul;
      if (!is)
	throw LowerTriMatrixReadError<T>(is);
      if (ul != 'L')
	throw LowerTriMatrixReadError<T>(is,'L',ul);
      size_t size;
      is >> size;
      if (!is) 
	throw LowerTriMatrixReadError<T>(is);
      m.reset(new LowerTriMatrix<T,D,S,I>(size));
      m->View().Read(is); 
      return is;
    }

  template <class T> std::istream& operator>>(
      std::istream& is, const LowerTriMatrixView<T>& m)
  { 
    char ul;
    is >> ul;
    if (!is)
      throw LowerTriMatrixReadError<T>(is);
    if (ul != 'L')
      throw LowerTriMatrixReadError<T>(is,'L',ul);
    size_t s;
    is >> s;
    if (!is)
      throw LowerTriMatrixReadError<T>(is);
    if (m.size() != s)
      throw LowerTriMatrixReadError<T>(m,is,s);
    TMVAssert(m.size() == s);
    m.Read(is);
    return is;
  }


#define InstFile "TMV_TriMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


