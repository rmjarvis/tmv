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
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_Divider.h"
#include "tmv/TMV_LUD.h"
#include "tmv/TMV_QRD.h"
#include "tmv/TMV_QRPD.h"
#include "tmv/TMV_SVD.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_VIt.h"
#include <iostream>
#include "portable_platform.h"

namespace tmv {

#define RT RealType(T)

#ifdef TMV_BLOCKSIZE
#define PERM_BLOCKSIZE TMV_BLOCKSIZE/2
#else
#define PERM_BLOCKSIZE 32
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
      const int M = colsize();
      const int N = rowsize();
      for(int i=0;i<M;++i) {
        TMVAssert(vv[i].size() == rowsize());
        typename std::vector<T>::const_iterator vvi = vv[i].begin();
        for(int j=0;j<N;++j,++vi,++vvi) {
#ifdef TMVFLDEBUG
          TMVAssert(vi >= first);
          TMVAssert(vi < last);
#endif
          *vi = *vvi;
        }
      }
    } else { 
      const int M = colsize();
      const int N = rowsize();
      for(int j=0;j<N;++j) {
        for(int i=0;i<M;++i,++vi) {
          TMVAssert(vv[i].size() == rowsize());
#ifdef TMVFLDEBUG
          TMVAssert(vi >= first);
          TMVAssert(vi < last);
#endif
          *vi = vv[i][j];
        }
      }
    }
  }

#undef NEW_SIZE


  //
  // Access
  //

  template <class T> T GenMatrix<T>::cref(int i, int j) const
  {
    const T* mi = cptr() + int(i)*stepi() + int(j)*stepj();
    return isconj() ? CONJ(*mi) : *mi;
  }

  template <class T, IndexStyle I> RefType(T) MatrixView<T,I>::ref(
      int i, int j) const
  {
    T* mi = ptr() + int(i)*itssi + int(j)*stepj();
    return REF(mi,ct());
  }

  template <class T> void GenMatrix<T>::NewDivider() const
  {
    switch (this->GetDivType()) {
      case LU : this->SetDiv(new LUDiv<T>(*this,
                      this->IsDivInPlace())); 
                break;
      case QR : this->SetDiv(new QRDiv<T>(*this,
                      this->IsDivInPlace())); 
                break;
      case QRP : this->SetDiv(new QRPDiv<T>(*this,
                       this->IsDivInPlace())); 
                 break;
      case SV : this->SetDiv(new SVDiv<T>(*this,
                      this->IsDivInPlace())); 
                break;
      default : TMVAssert(FALSE);
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
      int i, int j, int istep, int jstep, int size) const 
  {
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
      int i, int j, int istep, int jstep, int size) const 
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
    return ok;
  }

  //
  // Norms
  //

  template <class T> RT GenMatrix<T>::NormSq(const RT scale) const
  {
    if (CanLinearize()) return ConstLinearView().NormSq(scale);
    else {
      RT sum(0);
      if (isrm()) {
        const int M = colsize();
        for(int i=0;i<M;++i) 
          sum += row(i).NormSq(scale);
      } else {
        const int N = rowsize();
        for(int j=0;j<N;++j) 
          sum += col(j).NormSq(scale);
      }
      return sum;
    }
  }

  template <class T> static RT NonLapMaxAbsElement(
      const GenMatrix<T>& m)
  {
    if (m.CanLinearize()) return m.ConstLinearView().MaxAbsElement();
    else {
      RT max(0);
      if (m.iscm()) {
        const int N = m.rowsize();
        for(int j=0;j<N;++j) {
          RT temp = m.col(j).NormInf();
          if (temp > max) max = temp;
        }
      } else {
        const int M = m.colsize();
        for(int i=0;i<M;++i) {
          RT temp = m.row(i).NormInf();
          if (temp > max) max = temp;
        }
      }
      return max;
    }
  }

  template <class T> static RT NonLapNorm1(const GenMatrix<T>& m)
  { 
    RT max(0);
    const int N = m.rowsize();
    for(int j=0;j<N;++j) {
      RT temp = m.col(j).Norm1();
      if (temp > max) max = temp;
    }
    return max;
  }

  template <class T> static RT NonLapNormF(
      const GenMatrix<T>& m)
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

  template <class T> static inline RT NonLapNormInf(
      const GenMatrix<T>& m)
  { return NonLapNorm1(m.Transpose()); }

#ifdef XLAP
  template <class T> static RT LapNorm(const char c,
      const GenMatrix<T>& m)
  { 
    switch(c) {
      case 'M' : return NonLapMaxAbsElement(m);
      case '1' : return NonLapNorm1(m);
      case 'F' : return NonLapNormF(m);
      case 'I' : return NonLapNormInf(m);
      default : TMVAssert(FALSE); 
    }
    return RT(0);
  }
#ifdef INST_DOUBLE
  template <> double LapNorm(const char c, const GenMatrix<double>& m)
  { 
    TMVAssert(m.iscm());
    char cc = c;
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
#ifndef LAPNOWORK
    auto_array<double> work(c == 'I' ? new double[M] : 0);
#endif
    double norm = LAPNAME(dlange) (LAPCM LAPV(cc),LAPV(M),LAPV(N),
        LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()));
    return norm;
  }
  template <> double LapNorm(const char c, 
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
    double norm = LAPNAME(zlange) (LAPCM LAPV(cc),LAPV(M),LAPV(N),
        LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()));
    return norm;
  }
#endif
#ifdef INST_FLOAT
  template <> float LapNorm(const char c, const GenMatrix<float>& m)
  { 
    TMVAssert(m.iscm());
    char cc = c;
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
#ifndef LAPNOWORK
    auto_array<float> work(c == 'I' ? new float[M] : 0);
#endif
    double norm = LAPNAME(slange) (LAPCM LAPV(cc),LAPV(M),LAPV(N),
        LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()));
    return norm;
  }
  template <> float LapNorm(const char c, 
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
    double norm = LAPNAME(clange) (LAPCM LAPV(cc),LAPV(M),LAPV(N),
        LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()));
    return norm;
  }
#endif
#endif // XLAP

  template <class T> RT GenMatrix<T>::MaxAbsElement() const
  {
#ifdef XLAP
    if (isrm() && stepi() > 0) return LapNorm('M',Transpose());
    else if (iscm() && stepj() > 0) return LapNorm('M',*this);
    else
#endif
      return NonLapMaxAbsElement(*this);
  }
  template <class T> RT GenMatrix<T>::Norm1() const
  {
#ifdef XLAP
    if (isrm() && stepi() > 0) return LapNorm('I',Transpose());
    else if (iscm() && stepj() > 0) return LapNorm('1',*this);
    else
#endif
      return NonLapNorm1(*this);
  }
  template <class T> RT GenMatrix<T>::NormF() const
  {
#ifdef XLAP
    if (isrm() && stepi() > 0) return LapNorm('F',Transpose());
    else if (iscm() && stepj() > 0) return LapNorm('F',*this);
    else
#endif
      return NonLapNormF(*this);
  }

  template <class T> RT GenMatrix<T>::DoNorm2() const
  {
    if (this->colsize() < this->rowsize()) return Transpose().DoNorm2();
    if (this->rowsize() == 0) return RT(0);
    Matrix<T> m = *this;
    DiagMatrix<RT> S(this->rowsize());
    SV_Decompose(m.View(),S.View(),false);
    return S(0);
  }

  template <class T> RT GenMatrix<T>::DoCondition() const
  {
    if (this->colsize() < this->rowsize()) return Transpose().DoNorm2();
    if (this->rowsize() == 0) return RT(1);
    Matrix<T> m = *this;
    DiagMatrix<RT> S(this->rowsize());
    SV_Decompose(m.View(),S.View(),false);
    return S(0)/S(S.size()-1);
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
  const MatrixView<T,I>& MatrixView<T,I>::Clip(RT thresh) const
  { 
    TMVAssert(I==CStyle);
    if (this->CanLinearize()) LinearView().Clip(thresh);
    else {
      if (this->isrm()) {
        const int M = colsize();
        for(int i=0;i<M;++i) row(i).Clip(thresh);
      }
      else {
        const int N = rowsize();
        for(int j=0;j<N;++j) col(j).Clip(thresh);
      }
    }
    return *this; 
  }

  template <class T, IndexStyle I> 
  const MatrixView<T,I>& MatrixView<T,I>::Zero() const
  { 
    TMVAssert(I==CStyle);
    if (this->CanLinearize()) LinearView().Zero();
    else {
      if (this->isrm()) {
        const int M = colsize();
        for(int i=0;i<M;++i) row(i).Zero();
      }
      else {
        const int N = rowsize();
        for(int j=0;j<N;++j) col(j).Zero();
      }
    }
    return *this; 
  }

  template <class T, IndexStyle I> 
  const MatrixView<T,I>& MatrixView<T,I>::SetAllTo(T x) const
  { 
    TMVAssert(I==CStyle);
    if (this->CanLinearize()) LinearView().SetAllTo(x);
    else {
      if (this->isrm()) {
        const int M = colsize();
        for(int i=0;i<M;++i) row(i).SetAllTo(x); 
      } else  {
        const int N = rowsize();
        for(int j=0;j<N;++j) col(j).SetAllTo(x); 
      }
    }
    return *this; 
  }

  template <class T, IndexStyle I> 
  const MatrixView<T,I>& MatrixView<T,I>::TransposeSelf() const
  {
    TMVAssert(I==CStyle);
    TMVAssert(colsize() == rowsize());
    const int M = colsize();
    for(int i=1;i<M;++i) Swap(row(i,0,i),col(i,0,i));
    return *this; 
  }

  template <class T, IndexStyle I> 
  const MatrixView<T,I>& MatrixView<T,I>::ConjugateSelf() const
  { 
    TMVAssert(I==CStyle);
    if (IsComplex(T())) {
      if (this->CanLinearize()) LinearView().ConjugateSelf();
      else {
        if (this->isrm()) {
          const int M = colsize();
          for(int i=0;i<M;++i) row(i).ConjugateSelf();
        } else  {
          const int N = rowsize();
          for(int j=0;j<N;++j) col(j).ConjugateSelf();
        }
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
      const int* p, int i1, int i2) const
  { 
    TMVAssert(I==CStyle);
    TMVAssert(i2<=int(colsize()));
    TMVAssert(i1<=i2);
    // This idea of doing the permutations a block at a time
    // is cribbed from the LAPack code.  It does speed things up 
    // quite a bit for large matrices.  On my machine where BLOCKSIZE=64
    // is optimal for most routines, blocks of 32 were optimal here,
    // so I use BLOCKSIZE/2 in general.
    const int N = rowsize();
    const int Nx = N/PERM_BLOCKSIZE*PERM_BLOCKSIZE;
    if (Nx != 0) {
      for(int j=0;j<Nx;) {
        int j2 = j+PERM_BLOCKSIZE;
        const int* pi = p+i1;
        for(int i=i1;i<i2;++i,++pi) {
          TMVAssert(*pi < int(colsize()));
          Cols(j,j2).SwapRows(i,*pi);
        }
        j = j2;
      }
    }
    if (Nx != N) {
      const int* pi = p+i1;
      for(int i=i1;i<i2;++i,++pi) {
        TMVAssert(*pi < int(colsize()));
        Cols(Nx,N).SwapRows(i,*pi);
      }
    }
    return *this;
  }

  template <class T, IndexStyle I> 
  const MatrixView<T,I>& MatrixView<T,I>::ReversePermuteRows(
      const int* p, int i1, int i2) const
  { 
    TMVAssert(I==CStyle);
    TMVAssert(i2<=int(colsize()));
    TMVAssert(i1<=i2);
    const int N = rowsize();
    const int Nx = N/PERM_BLOCKSIZE*PERM_BLOCKSIZE;
    if (Nx != 0) {
      for(int j=0;j<Nx;) {
        int j2 = j+PERM_BLOCKSIZE;
        const int* pi = p+i2;
        for(int i=i2;i>i1;) {
          --i; --pi;
          TMVAssert(*pi < int(colsize()));
          Cols(j,j2).SwapRows(i,*pi);
        }
        j = j2;
      }
    }
    if (Nx != N) {
      const int* pi = p+i2;
      for(int i=i2;i>i1;) {
        --i; --pi;
        TMVAssert(*pi < int(colsize()));
        Cols(Nx,N).SwapRows(i,*pi);
      }
    }
    return *this;
  }

  //
  // Copy Matrices
  //

  template <class T> static void NonLapCopy(
      const GenMatrix<T>& m1, const MatrixView<T>& m2)
  {
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert(m1.ct()==NonConj);
    TMVAssert(m2.ct()==NonConj);
    TMVAssert(!(m1.isrm() && m2.isrm()));
    const int M = m2.colsize();
    const int N = m2.rowsize();

    if (m1.iscm() && m2.iscm()) {
      const T* p1 = m1.cptr();
      T* p2 = m2.ptr();
      const int s1 = m1.stepj();
      const int s2 = m2.stepj();
      for(int j=0;j<N;++j,p1+=s1,p2+=s2) {
        std::copy(p1,p1+M,p2);
      }
    } else if (M > N) {
      if (ShouldReverse(m1.stepi(),m2.stepi()))
        for(int j=0;j<N;++j) 
          DoCopySameType(m1.col(j).Reverse(),m2.col(j).Reverse());
      else
        for(int j=0;j<N;++j) 
          DoCopySameType(m1.col(j),m2.col(j));
    } else {
      if (ShouldReverse(m1.stepj(),m2.stepj()))
        for(int i=0;i<M;++i) 
          DoCopySameType(m1.row(i).Reverse(),m2.row(i).Reverse());
      else
        for(int i=0;i<M;++i) 
          DoCopySameType(m1.row(i),m2.row(i));
    }
  }
#ifdef ELAP
  template <class T> static inline void LapCopy(
      const GenMatrix<T>& m1, const MatrixView<T>& m2)
  { NonLapCopy(m1,m2); }
#ifdef INST_DOUBLE
  template <> void LapCopy(
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
    LAPNAME(dlacpy) (LAPCM LAPV(c),LAPV(m),LAPV(n),LAPP(m1.cptr()),LAPV(ld1),
        LAPP(m2.ptr()),LAPV(ld2));
  }
  template <> void LapCopy(
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
    LAPNAME(zlacpy) (LAPCM LAPV(c),LAPV(m),LAPV(n),LAPP(m1.cptr()),LAPV(ld1),
        LAPP(m2.ptr()),LAPV(ld2));
  }
#endif
#ifdef INST_FLOAT
  template <> void LapCopy(
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
    LAPNAME(slacpy) (LAPCM LAPV(c),LAPV(m),LAPV(n),LAPP(m1.cptr()),LAPV(ld1),
        LAPP(m2.ptr()),LAPV(ld2));
  }
  template <> void LapCopy(
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
    LAPNAME(clacpy) (LAPCM LAPV(c),LAPV(m),LAPV(n),LAPP(m1.cptr()),LAPV(ld1),
        LAPP(m2.ptr()),LAPV(ld2));
  }
#endif
#endif
  template <class T> void DoCopySameType(const GenMatrix<T>& m1, 
      const MatrixView<T>& m2)
  {
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert(m2.rowsize() > 0);
    TMVAssert(m2.colsize() > 0);
    TMVAssert(m1.ct() == NonConj);
    TMVAssert(m2.ct() == NonConj);
    TMVAssert(!m2.isrm());
    TMVAssert(!m2.SameAs(m1));

#ifdef ELAP
    if (m1.iscm() && m2.iscm() && m1.stepj()>0 && m2.stepj()>0) 
      LapCopy(m1,m2);
    else
#endif
      NonLapCopy(m1,m2);
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
      if (m1.isrm() && m2.isrm()) {
        const int M = m1.colsize();
        for(int i=0;i<M;++i) Swap(m1.row(i),m2.row(i)); 
      } else {
        const int N = m1.rowsize();
        for(int j=0;j<N;++j) Swap(m1.col(j),m2.col(j)); 
      }
    }
  }

  //
  // m1 == m2
  //

  template <class T1, class T2> bool operator==(
      const GenMatrix<T1>& m1, const GenMatrix<T2>& m2)
  {
    if (m1.colsize() != m2.colsize()) return false;
    else if (m1.rowsize() != m2.rowsize()) return false;
    else if (m1.SameAs(m2)) return true;
    else if (m1.stepi()==m2.stepi() && m1.stepj()==m2.stepj() &&
        m1.CanLinearize() && m2.CanLinearize())
      return m1.ConstLinearView() == m2.ConstLinearView();
    else {
      const int M = m1.colsize();
      for(int i=0;i<M;++i) 
        if (m1.row(i) != m2.row(i)) return false;
      return true;  
    }
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

  template <bool conj, bool rm, bool th, class T> 
  static void DoWrite(std::ostream& os, const GenMatrix<T>& m,
      RT thresh)
  {
    const T* mrowi = m.cptr();
    const int sj = rm ? 1 : m.stepj();
    os << m.colsize() <<"  "<<m.rowsize()<<std::endl;
    for(int i=m.colsize();i>0;--i,mrowi+=m.stepi()) {
      os << "( ";
      const T* mij = mrowi;
      for(int k=m.rowsize();k>0;--k,rm?++mij:mij+=sj) 
        if (conj) 
          if (th)
            os << ' '<<Value(ABS(*mij) < thresh ? T(0) : CONJ(*mij))<<' ';
          else
            os << ' '<<Value(CONJ(*mij))<<' ';
        else 
          if (th)
            os << ' '<<Value(ABS(*mij) < thresh ? T(0) : *mij)<<' ';
          else
            os << ' '<<Value(*mij)<<' ';
      os << " )\n";
    }
  }

  template <bool rm, bool th, class T> static inline void DoWrite1(
      std::ostream& os, const GenMatrix<T>& m, T thresh)
  { DoWrite<false,rm,th>(os,m,thresh); }

  template <bool rm, bool th, class T> static inline void DoWrite1(
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
      DoWrite1<true,false>(os,*this,RT(0));
    else
      DoWrite1<false,false>(os,*this,RT(0));
  }

  template <class T> void GenMatrix<T>::Write(std::ostream& os,
      RT thresh) const
  {
    if (isrm())
      DoWrite1<true,true>(os,*this,thresh);
    else
      DoWrite1<false,true>(os,*this,thresh);
  }

#ifndef NOTHROW
  template <class T> class MatrixReadError :
    public ReadError
  {
  public :
    int i,j;
    mutable auto_ptr<Matrix<T> > m;
    char exp,got;
    size_t cs,rs;
    bool is,iseof,isbad;

    MatrixReadError(int _i, int _j, const GenMatrix<T>& _m,
        std::istream& _is) throw() :
      ReadError("Matrix."),
      i(_i), j(_j), m(new Matrix<T>(_m)), exp(0), got(0),
      cs(_m.colsize()), rs(_m.rowsize()),
      is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
    MatrixReadError(std::istream& _is) throw() :
      ReadError("Matrix."),
      i(0), j(0), m(0), exp(0), got(0), cs(0), rs(0),
      is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
    MatrixReadError(int _i, int _j, const GenMatrix<T>& _m,
        std::istream& _is, char _e, char _g) throw() :
      ReadError("Matrix."),
      i(_i), j(_j), m(new Matrix<T>(_m)), exp(_e), got(_g),
      cs(_m.colsize()), rs(_m.rowsize()),
      is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
    MatrixReadError(const GenMatrix<T>& _m,
        std::istream& _is, size_t _cs, size_t _rs) throw() :
      ReadError("Matrix."),
      i(0), j(0), m(new Matrix<T>(_m)), exp(0), got(0), cs(_cs), rs(_rs),
      is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

    MatrixReadError(const MatrixReadError<T>& rhs) :
      i(rhs.i), j(rhs.j), m(rhs.m), exp(rhs.exp), got(rhs.got),
      cs(rhs.cs), rs(rhs.rs),
      is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
    virtual ~MatrixReadError() throw() {}

    virtual void Write(std::ostream& os) const throw()
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
        const int N = m->rowsize();
        os<<"The portion of the Matrix which was successfully read is: \n";
        for(int ii=0;ii<i;++ii) {
          os<<"( ";
          for(int jj=0;jj<N;++jj)
            os<<' '<<(*m)(ii,jj)<<' ';
          os<<" )\n";
        }
        os<<"( ";
        for(int jj=0;jj<j;++jj)
          os<<' '<<(*m)(i,jj)<<' ';
        os<<" )\n";
      }
    }
  };
#endif

  template <class T, IndexStyle I> void MatrixView<T,I>::Read(
      std::istream& is) const
  {
    TMVAssert(I==CStyle);
    T* mrowi = ptr();
    const int sj = stepj();
    char paren;
    const int M = colsize();
    for(int i=0;i<M;++i,mrowi+=stepi()) {
      is >> paren;
      if (!is || paren != '(') 
#ifdef NOTHROW
      { std::cerr<<"Matrix ReadError: "<<paren<<" != (\n"; exit(1); }
#else
      throw MatrixReadError<T>(i,0,*this,is,'(',is?paren:'(');
#endif
      T* mij = mrowi;
      if (this->isrm()) 
        for(int k=rowsize();k>0;--k,++mij) {
          is >> *mij;
          if (!is) 
#ifdef NOTHROW
          { std::cerr<<"Matrix ReadError: !is\n"; exit(1); }
#else
          throw MatrixReadError<T>(i,rowsize()-k,*this,is);
#endif
        }
      else
        for(int k=rowsize();k>0;--k,mij+=sj) {
          is >> *mij;
          if (!is) 
#ifdef NOTHROW
          { std::cerr<<"Matrix ReadError: !is\n"; exit(1); }
#else
          throw MatrixReadError<T>(i,rowsize()-k,*this,is);
#endif
        }
      is >> paren;
      if ((!is && i+1<M)  || paren != ')') 
#ifdef NOTHROW
      { std::cerr<<"Matrix ReadError: "<<paren<<" != )\n"; exit(1); }
#else
      throw MatrixReadError<T>(i,rowsize(),*this,is,')',is?paren:')');
#endif
    }
    if (this->isconj()) ConjugateSelf();
  }

  template <class T, StorageType S, IndexStyle I> std::istream& operator>>(
      std::istream& is, auto_ptr<Matrix<T,S,I> >& m)
  { 
    size_t cs,rs;
    is >> cs >> rs;
    if (!is) 
#ifdef NOTHROW
    { std::cerr<<"Matrix ReadError: !is\n"; exit(1); }
#else
    throw MatrixReadError<T>(is);
#endif
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
#ifdef NOTHROW
    { std::cerr<<"Matrix ReadError: !is\n"; exit(1); }
#else
    throw MatrixReadError<T>(is);
#endif
    if (m.colsize() != cs || m.rowsize() != rs) 
#ifdef NOTHROW
    { std::cerr<<"Matrix ReadError: Wrong size\n"; exit(1); }
#else
    throw MatrixReadError<T>(m,is,cs,rs);
#endif
    TMVAssert(m.colsize() == cs);
    TMVAssert(m.rowsize() == rs);
    m.Read(is);
    return is;
  }

#undef RT

#define InstFile "TMV_Matrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv

