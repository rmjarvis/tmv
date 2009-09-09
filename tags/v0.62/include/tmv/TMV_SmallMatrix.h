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


//---------------------------------------------------------------------------
//
// This file defines the TMV SmallMatrix class.
//
// Constructors:
//
//    SmallMatrix<T,M,N,stor,I>()
//        Makes a SmallMatrix with column size = M and row size = N
//        with _uninitialized_ values
//
//    SmallMatrix<T,M,N,stor,I>(T x)
//        Makes a SmallMatrix of size n with all values = x
//
//    SmallMatrix<T,M,N,stor,I>(const vector<vector<T> >& m)
//        Makes a SmallMatrix with a_ij = m[i][j]
//
//    SmallMatrix<T,M,N,stor,I>(const T* m)
//    SmallMatrix<T,M,N,stor,I>(const vector<T>& m)
//    SmallMatrix<T,M,N,stor,I>(const GenMatrix<T>& m)
//        Make a SmallMatrix which copies the elements of m.
//
// SmallMatrix doesn't have views like a regular Matrix.
// All the normal viewing kinds of routines just return a regular MatrixView.
// It is mostly useful for fast element access and simple tasks
// like multiplication and addition.  For most division routines,
// SmallMatrix just sends the task to a regular matrix.  The exception
// is 2x2, which has specialization for Det, Inverse, etc.


#ifndef TMV_SmallMatrix_H
#define TMV_SmallMatrix_H

#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_SmallVector.h"

namespace tmv {

#define MIN(M,N) (M<N ? M : N)

  template <class T1, class T2, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
  inline void Copy(const SmallMatrix<T1,M,N,S1,I1>& m1,
      SmallMatrix<T2,M,N,S2,I2>& m2);

  template <class T, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
  inline void Copy(const SmallMatrix<std::complex<T>,M,N,S1,I1>& ,
      SmallMatrix<T,M,N,S2,I2>& )
  { TMVAssert(FALSE); }

  template <class T, class Tm, int M, int N, StorageType S, IndexStyle I> 
  class QuotXm_1;

  template <class T, int M, int N> 
  class SmallMatrixComposite;

  template <class T, int M, int N, StorageType S, IndexStyle I> 
  inline T DoDet(const SmallMatrix<T,M,N,S,I>& m);

  template <class T, class T2, int M, int N, StorageType S, StorageType S2, IndexStyle I, IndexStyle I2> 
  inline void DoInverse(const SmallMatrix<T,M,N,S,I>& m,
      SmallMatrix<T2,N,M,S2,I2>& minv);

  template <class T, int M, int N, StorageType S, StorageType S2, IndexStyle I, IndexStyle I2> 
  inline void DoInverseATA(const SmallMatrix<T,M,N,S,I>& m, 
      SmallMatrix<T,N,N,S2,I2>& ata);

#ifdef XTEST
#ifdef TMVDEBUG
#define XTEST_DEBUG
#endif
#endif

#define Si (S==RowMajor?N:1)
#define Sj (S==RowMajor?1:M)
  // defaults S=ColMajor and I=CStyle are set in TMV_BaseMatrix.h
  template <class T, int M, int N, StorageType S, IndexStyle I> 
  class SmallMatrix 
  {

  public:

    //
    // Constructors
    //

    inline SmallMatrix() 
    {
      TMVAssert(M>0);
      TMVAssert(N>0);
      TMVAssert(S==RowMajor || S==ColMajor);
#ifdef TMVDEBUG
      SetAllTo(T(888));
#endif
    }

    explicit inline SmallMatrix(T x) 
    { 
      TMVAssert(M>0);
      TMVAssert(N>0);
      TMVAssert(S==RowMajor || S==ColMajor);
      if (x == T(0)) Zero();
      else SetAllTo(x);
    }

    inline SmallMatrix(const T* vv) 
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(M>0);
      TMVAssert(N>0);
      TMVAssert(S==RowMajor || S==ColMajor);
      for(int i=0;i<M*N;++i) itsm[i] = vv[i];
    }

    inline SmallMatrix(const std::vector<T>& vv) 
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(M>0);
      TMVAssert(N>0);
      TMVAssert(S==RowMajor || S==ColMajor);
      TMVAssert(vv.size() == M*N);
      for(int i=0;i<M*N;++i) itsm[i] = vv[i];
    }

    explicit inline SmallMatrix(
        const std::vector<std::vector<T> >& vv) 
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(M>0);
      TMVAssert(N>0);
      TMVAssert(S==RowMajor || S==ColMajor);
      for(int i=0;i<M;++i) {
        TMVAssert(vv[i].size() == N);
        for(int j=0;j<N;++j) ref(i,j) = vv[i][j];
      }
    }

    inline SmallMatrix(const SmallMatrix<T,M,N,S,I>& m2) 
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(M>0);
      TMVAssert(N>0);
      TMVAssert(S==RowMajor || S==ColMajor);
      for(int i=0;i<M*N;++i) itsm[i] = m2.itsm[i];
    }

    template <IndexStyle I2> 
    inline SmallMatrix(const SmallMatrix<T,M,N,S,I2>& m2) 
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(M>0);
      TMVAssert(N>0);
      TMVAssert(S==RowMajor || S==ColMajor);
      for(int i=0;i<M*N;++i) itsm[i] = m2.cptr()[i];
    }

    template <class T2, StorageType S2, IndexStyle I2> 
    inline SmallMatrix(const SmallMatrix<T2,M,N,S2,I2>& m2) 
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(M>0);
      TMVAssert(N>0);
      TMVAssert(S==RowMajor || S==ColMajor);
      TMVAssert(IsComplex(T()) || IsReal(T2()));
      Copy(m2,*this);
    }

    inline SmallMatrix(const GenMatrix<T>& m2) 
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(M>0);
      TMVAssert(N>0);
      TMVAssert(S==RowMajor || S==ColMajor);
      TMVAssert(m2.colsize() == M);
      TMVAssert(m2.rowsize() == N);
      View() = m2;
    }

    template <class T2> 
    inline SmallMatrix(const GenMatrix<T2>& m2) 
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(M>0);
      TMVAssert(N>0);
      TMVAssert(IsComplex(T()) || IsReal(T2()));
      TMVAssert(S==RowMajor || S==ColMajor);
      TMVAssert(m2.colsize() == M);
      TMVAssert(m2.rowsize() == N);
      View() = m2;
    }

    inline SmallMatrix(const AssignableToMatrix<RealType(T)>& m2) 
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(M>0);
      TMVAssert(N>0);
      TMVAssert(S==RowMajor || S==ColMajor);
      TMVAssert(m2.colsize() == M);
      TMVAssert(m2.rowsize() == N);
      View() = m2;
    }

    inline SmallMatrix(const AssignableToMatrix<ComplexType(T)>& m2) 
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(M>0);
      TMVAssert(N>0);
      TMVAssert(IsComplex(T()));
      TMVAssert(S==RowMajor || S==ColMajor);
      TMVAssert(m2.colsize() == M);
      TMVAssert(m2.rowsize() == N);
      View() = m2;
    }

    inline SmallMatrix(
        const SmallMatrixComposite<RealType(T),M,N>& m2) 
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(M>0);
      TMVAssert(N>0);
      TMVAssert(S==RowMajor || S==ColMajor);
      m2.AssignTom(*this);
    }

    inline SmallMatrix(
        const SmallMatrixComposite<ComplexType(T),M,N>& m2) 
    { 
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(M>0);
      TMVAssert(N>0);
      TMVAssert(S==RowMajor || S==ColMajor);
      m2.AssignTom(*this);
    }

    inline ~SmallMatrix()
    {
#ifdef TMVDEBUG
      SetAllTo(T(999));
#endif
    }


    //
    // Op=
    //

    inline SmallMatrix<T,M,N,S,I>& operator=(
        const SmallMatrix<T,M,N,S,I>& m2)
    { 
      if (&m2 != this) Copy(m2,*this);
      return *this; 
    }

    template <IndexStyle I2> 
    inline SmallMatrix<T,M,N,S,I>& operator=(
        const SmallMatrix<T,M,N,S,I2>& m2)
    { 
      Copy(m2,*this);
      return *this; 
    }

    template <class T2, StorageType S2, IndexStyle I2> 
    inline SmallMatrix<T,M,N,S,I>& operator=(
        const SmallMatrix<T2,M,N,S2,I2>& m2)
    { 
      TMVAssert(IsComplex(T()) || IsReal(T2()));
      Copy(m2,*this);
      return *this; 
    }

    inline SmallMatrix<T,M,N,S,I>& operator=(T x) 
    { return SetToIdentity(x); }

    template <class T2> 
    inline SmallMatrix<T,M,N,S,I>& operator=(const GenMatrix<T2>& m2)
    {
      TMVAssert(m2.colsize() == M);
      TMVAssert(m2.rowsize() == N);
      TMVAssert(IsComplex(T()) || IsReal(T2()));
      View() = m2;
      return *this;
    }

    inline SmallMatrix<T,M,N,S,I>& operator=(
        const AssignableToMatrix<RealType(T)>& m2)
    {
      TMVAssert(m2.colsize() == M);
      TMVAssert(m2.rowsize() == N);
      View() = m2;
      return *this;
    }

    inline SmallMatrix<T,M,N,S,I>& operator=(
        const AssignableToMatrix<ComplexType(T)>& m2)
    {
      TMVAssert(m2.colsize() == M);
      TMVAssert(m2.rowsize() == N);
      TMVAssert(IsComplex(T()));
      View() = m2;
      return *this;
    }

    inline SmallMatrix<T,M,N,S,I>& operator=(
        const SmallMatrixComposite<RealType(T),M,N>& m2)
    {
      m2.AssignTom(*this);
      return *this;
    }

    inline SmallMatrix<T,M,N,S,I>& operator=(
        const SmallMatrixComposite<ComplexType(T),M,N>& m2)
    {
      m2.AssignTom(*this);
      return *this;
    }

    ListAssigner<T,VIt<T,Unit,NonConj> > inline operator=(ListInitClass)
    { 
      return ListAssigner<T,VIt<T,Unit,NonConj> >(
          VIt<T,Unit,NonConj>(ptr(),1),M*N); 
    }

    //
    // Access
    //

    typedef T& reference;

    inline T operator()(int i,int j) const
    { 
      if (I == CStyle) {
        TMVAssert(i>=0 && i<M);
        TMVAssert(j>=0 && j<N);
        return cref(i,j); 
      } else {
        TMVAssert(i>=1 && i<=M);
        TMVAssert(j>=1 && j<=N);
        return cref(i-1,j-1); 
      }
    }

    inline T& operator()(int i,int j) 
    { 
      if (I == CStyle) {
        TMVAssert(i>=0 && i<M);
        TMVAssert(j>=0 && j<N);
        return ref(i,j);
      } else {
        TMVAssert(i>=1 && i<=M);
        TMVAssert(j>=1 && j<=N);
        return ref(i-1,j-1);
      }
    }

    inline ConstVectorView<T,I> row(int i) const 
    { return View().row(i); }

    inline ConstVectorView<T,I> row(int i, int j1, int j2) const 
    { return View().row(i,j1,j2); }

    inline ConstVectorView<T,I> operator[](int i) const
    { return View().row(i); }

    inline ConstVectorView<T,I> col(int j) const
    { return View().col(j); }

    inline ConstVectorView<T,I> col(int j, int i1, int i2) const
    { return View().col(j,i1,i2); }

    inline ConstVectorView<T,I> diag() const
    { return View().diag(); }

    inline ConstVectorView<T,I> diag(int i) const
    { return View().diag(i); }

    inline ConstVectorView<T,I> diag(int i, int j1, int j2) const
    { return View().diag(i,j1,j2); }

    inline VectorView<T,I> row(int i)
    { return View().row(i); }

    inline VectorView<T,I> row(int i, int j1, int j2)
    { return View().row(i,j1,j2); }

    inline VectorView<T,I> operator[](int i)
    { return View().row(i); }

    inline VectorView<T,I> col(int j)
    { return View().col(j); }

    inline VectorView<T,I> col(int j, int i1, int i2)
    { return View().col(j,i1,i2); }

    inline VectorView<T,I> diag()
    { return View().diag(); }

    inline VectorView<T,I> diag(int i)
    { return View().diag(i); }

    inline VectorView<T,I> diag(int i, int j1, int j2) 
    { return View().diag(i,j1,j2); }

    //
    // Functions of Matrix
    //

    inline T Trace() const
    { 
      TMVAssert(M == N);
      T sum(0);
      for(int i=0; i<M*N; i+=M+1) sum += itsm[i];
      return sum;
    }

    inline RealType(T) Norm() const 
    { return NormF(); }

    inline RealType(T) NormF() const
    { return SQRT(NormSq()); }

    // NormF()^2
    inline RealType(T) NormSq(RealType(T) scale=RealType(T)(1)) const
    { 
      RealType(T) sum(0);
      if (scale == RealType(T)(1))
        for(int i=0;i<M*N; ++i) sum += NORM(itsm[i]);
      else
        for(int i=0;i<M*N; ++i) sum += NORM(itsm[i]*scale);
      return sum;
    }

    // 1-Norm = max_j (sum_i |a_ij|)
    inline RealType(T) Norm1() const
    {
      RealType(T) max(0);
      for(int j=0;j<N;++j) {
        RealType(T) temp(0);
        for(int i=0;i<M;++i) temp += ABS(cref(i,j));
        if (temp > max) max = temp;
      }
      return max;
    }

    // inf-Norm = max_i (sum_j |a_ij|)
    inline RealType(T) NormInf() const
    {
      RealType(T) max(0);
      for(int i=0;i<M;++i) {
        RealType(T) temp(0);
        for(int j=0;j<N;++j) temp += ABS(cref(i,j));
        if (temp > max) max = temp;
      }
      return max;
    }

    // = max_i,j (|a_ij|)
    inline RealType(T) MaxAbsElement() const
    { 
      RealType(T) max(0);
      for(int i=0;i<M*N; ++i) {
        RealType(T) temp = ABS(itsm[i]);
        if (temp > max) max = temp;
      }
      return max;
    }

    inline T Det() const
    {
      TMVAssert(M == N);
      return DoDet(*this);
    }
    inline RealType(T) LogDet(T* sign=0) const
    {
      TMVAssert(M==N);
      if (M <= 3) {
        T det = Det();
        RealType(T) absdet = ABS(det);
        if (sign) {
          if (IsReal(T())) *sign = REAL(det) > 0 ?
            RealType(T)(1) : RealType(T)(-1);
          else *sign = det / absdet;
        }
        return LOG(absdet);
      }
      else return View().LogDet(sign);
    }
    inline RealType(T) Norm2() const
    { return View().DoNorm2(); }
    inline RealType(T) DoNorm2() const
    { return View().DoNorm2(); }
    inline bool Singular() const
    { return View().Singular(); }
    inline RealType(T) Condition() const
    { return View().DoCondition(); }
    inline RealType(T) DoCondition() const
    { return View().DoCondition(); }

    // 
    // Division Control
    //

    inline QuotXm_1<T,T,M,N,S,I> Inverse() const
    { return QuotXm_1<T,T,M,N,S,I>(T(1),*this); }

    inline void Inverse(const MatrixView<T>& minv) const
    { View().Inverse(minv); }

    template <class T1> 
    inline void Inverse(const MatrixView<T1>& minv) const
    { View().Inverse(minv); }

    template <class T2, StorageType S2, IndexStyle I2> 
    inline void Inverse(Matrix<T2,S2,I2>& minv) const
    { View().Inverse(minv); }

    inline void InverseATA(const MatrixView<T>& ata) const
    { View().InverseATA(ata); }

    template <StorageType S2, IndexStyle I2> 
    inline void InverseATA(Matrix<T,S2,I2>& ata) const
    { View().InverseATA(ata); }

    template <class T2, StorageType S2, IndexStyle I2> 
    inline void Inverse(SmallMatrix<T2,N,M,S2,I2>& minv) const
    { DoInverse(*this,minv); }

    template <class T2, StorageType S2, IndexStyle I2> 
    inline void InverseATA(SmallMatrix<T2,N,N,S2,I2>& ata) const
    { DoInverseATA(*this,ata); }

    //
    // Modifying Functions
    //

    inline SmallMatrix<T,M,N,S,I>& Zero() 
    { for(int i=0;i<M*N;++i) itsm[i] = T(0); }

    inline SmallMatrix<T,M,N,S,I>& Clip(RealType(T) thresh)
    { 
      for(int i=0; i<M*N; ++i)
        if (ABS(itsm[i]) < thresh) itsm[i] = T(0);
      return *this;
    }

    inline SmallMatrix<T,M,N,S,I>& SetAllTo(T x) 
    {
      for(int i=0; i<M*N; ++i) itsm[i] = x;
      return *this;
    }

    inline SmallMatrix<T,M,N,S,I>& TransposeSelf() 
    {
      TMVAssert(M == N);
      for(int i=1; i<M; ++i) 
        for(int j=0; j<i; ++j) __TMV_SWAP(ref(i,j),ref(j,i));
      return *this;
    }

    inline SmallMatrix<T,M,N,S,I>& ConjugateSelf() 
    {
      if (IsComplex(T())) {
        RealType(T)* itsmi = reinterpret_cast<RealType(T)*>(itsm)+1;
        for(int i=0;i<2*M*N;i+=2) itsmi[i] = -itsmi[i];
      }
      return *this;
    }

    inline SmallMatrix<T,M,N,S,I>& SetToIdentity(T x=T(1)) 
    { 
      TMVAssert(M == N);
      Zero();
      for(int i=0; i<M*N; i+=M+1) itsm[i] = x;
      return *this;
    }

    inline SmallMatrix<T,M,N,S,I>& SwapRows(int i1, int i2)
    {
      if (I == CStyle) {
        TMVAssert(i1 >= 0 && i1 < M);
        TMVAssert(i2 >= 0 && i2 < M);
        if (i1 != i2)
          for(int j=0; j<N; ++j) __TMV_SWAP(ref(i1,j),ref(i2,j));
      } else {
        TMVAssert(i1 >= 1 && i1 <= M);
        TMVAssert(i2 >= 1 && i2 <= M);
        if (i1 != i2)
          for(int j=0; j<N; ++j) __TMV_SWAP(ref(i1-1,j),ref(i2-1,j));
      }
      return *this;
    }

    inline SmallMatrix<T,M,N,S,I>& SwapCols(int j1, int j2)
    {
      if (I == CStyle) {
        TMVAssert(j1 >= 0 && j1 < N);
        TMVAssert(j2 >= 0 && j2 < N);
        if (j1 != j2)
          for(int i=0; i<M; ++i) __TMV_SWAP(ref(i,j1),ref(i,j2));
      } else {
        TMVAssert(j1 >= 1 && j1 <= N);
        TMVAssert(j2 >= 1 && j2 <= N);
        if (j1 != j2)
          for(int i=0; i<M; ++i) __TMV_SWAP(ref(i,j1-1),ref(i,j2-1));
      }
      return *this;
    }

    inline SmallMatrix<T,M,N,S,I>& PermuteRows(
        const int* p, int i1, int i2)
    {
      if (I == CStyle) {
        TMVAssert(i1 >= 0 && i1 <= i2 && i2 <= M);
        for(int i=i1;i<i2;++i) SwapRows(i,p[i]);
      } else {
        TMVAssert(i1 >= 1 && i1 <= i2 && i2 <= M);
        for(int i=i1-1;i<i2;++i) SwapRows(i,p[i]);
      }
      return *this;
    }

    inline SmallMatrix<T,M,N,S,I>& PermuteRows(const int* p)
    { PermuteRows(p,I==CStyle?0:1,M); return *this; }

    inline SmallMatrix<T,M,N,S,I>& PermuteCols(
        const int* p, int j1, int j2)
    {
      if (I == CStyle) {
        TMVAssert(j1 >= 0 && j1 <= j2 && j2 <= N);
        for(int j=j1;j<j2;++j) SwapCols(j,p[j]);
      } else {
        TMVAssert(j1 >= 1 && j1 <= j2 && j2 <= N);
        for(int j=j1-1;j<j2;++j) SwapCols(j,p[j]);
      }
      return *this;
    }

    inline SmallMatrix<T,M,N,S,I>& PermuteCols(const int* p)
    { PermuteCols(p,I==CStyle?0:1,N); return *this; }

    inline SmallMatrix<T,M,N,S,I>& ReversePermuteRows(
        const int* p, int i1, int i2)
    {
      if (I == CStyle) {
        TMVAssert(i1 >= 0 && i1 <= i2 && i2 <= M);
        for(int i=i2-1;i>=i1;--i) SwapRows(i,p[i]);
      } else {
        TMVAssert(i1 >= 1 && i1 <= i2 && i2 <= M);
        for(int i=i2-1;i>=i1-1;--i) SwapRows(i,p[i]);
      }
      return *this;
    }

    inline SmallMatrix<T,M,N,S,I>& ReversePermuteRows(const int* p)
    { ReversePermuteRows(p,I==CStyle?0:1,M); return *this; }

    inline SmallMatrix<T,M,N,S,I>& ReversePermuteCols(
        const int* p, int j1, int j2)
    {
      if (I == CStyle) {
        TMVAssert(j1 >= 0 && j1 <= j2 && j2 <= N);
        for(int j=j2-1;j>=j1;--j) SwapCols(j,p[j]);
      } else {
        TMVAssert(j1 >= 1 && j1 <= j2 && j2 <= N);
        for(int j=j2-1;j>=j1-1;--j) SwapCols(j,p[j]);
      }
      return *this;
    }

    inline SmallMatrix<T,M,N,S,I>& ReversePermuteCols(const int* p)
    { ReversePermuteCols(p,I==CStyle?0:1,N); return *this; }

    //
    // SubMatrix
    //

    inline ConstMatrixView<T,I> SubMatrix(
        int i1, int i2, int j1, int j2) const
    { return View().SubMatrix(i1,i2,j1,j2); }

    inline ConstMatrixView<T,I> SubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) const
    { return View().SubMatrix(i1,i2,j1,j2,istep,jstep); }

    inline ConstVectorView<T,I> SubVector(
        int i, int j, int istep, int jstep, int s) const
    { return View().SubVector(i,j,istep,jstep,s); }

    inline ConstMatrixView<T,I> ColPair(int j1, int j2) const
    { return View().ColPair(j1,j2); }

    inline ConstMatrixView<T,I> RowPair(int i1, int i2) const
    { return View().RowPair(i1,i2); }

    inline ConstMatrixView<T,I> Cols(int j1, int j2) const
    { return View().Cols(j1,j2); }

    inline ConstMatrixView<T,I> Rows(int i1, int i2) const
    { return View().Rows(i1,i2); }

    inline ConstMatrixView<RealType(T)> Real() const
    { return View().Real(); }

    inline ConstMatrixView<RealType(T)> Imag() const
    { return View().Imag(); }

    inline MatrixView<T,I> SubMatrix(int i1, int i2, int j1, int j2)
    { return View().SubMatrix(i1,i2,j1,j2); }

    inline MatrixView<T,I> SubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) 
    { return View().SubMatrix(i1,i2,j1,j2,istep,jstep); }

    inline VectorView<T,I> SubVector(
        int i, int j, int istep, int jstep, int s) 
    { return View().SubVector(i,j,istep,jstep,s); }

    inline MatrixView<T,I> ColPair(int j1, int j2) 
    { return View().ColPair(j1,j2); }

    inline MatrixView<T,I> RowPair(int i1, int i2) 
    { return View().RowPair(i1,i2); }

    inline MatrixView<T,I> Cols(int j1, int j2) 
    { return View().Cols(j1,j2); }

    inline MatrixView<T,I> Rows(int i1, int i2) 
    { return View().Rows(i1,i2); }

    inline MatrixView<RealType(T)> Real() 
    { return View().Real(); }

    inline MatrixView<RealType(T)> Imag() 
    { return View().Imag(); }


    //
    // Views
    //

    inline ConstMatrixView<T,I> View() const
    { return ConstMatrixView<T,I>(cptr(),M,N,Si,Sj,S,NonConj,M*N); }

    inline ConstMatrixView<T,I> Transpose() const
    {
      return ConstMatrixView<T,I>(cptr(),N,M,Sj,Si,TransOf(S),NonConj,M*N); 
    }

    inline ConstMatrixView<T,I> Conjugate() const
    {
      return ConstMatrixView<T,I>(cptr(),M,N,Si,Sj,S,
          IsReal(T())?NonConj:Conj,M*N); 
    }

    inline ConstMatrixView<T,I> Adjoint() const
    {
      return ConstMatrixView<T,I>(cptr(),N,M,Sj,Si,TransOf(S),
          IsReal(T())?NonConj:Conj,M*N); 
    }

    inline ConstUpperTriMatrixView<T,I> UpperTri(
        DiagType dt=NonUnitDiag) const
    { return ConstUpperTriMatrixView<T,I>(cptr(),N,Si,Sj,dt,S,NonConj); }

    inline ConstLowerTriMatrixView<T,I> LowerTri(
        DiagType dt=NonUnitDiag) const
    { return ConstLowerTriMatrixView<T,I>(cptr(),M,Si,Sj,dt,S,NonConj); }

    inline ConstVectorView<T,I> ConstLinearView() const
    { return ConstVectorView<T,I>(cptr(),M*N,1,NonConj); }

    inline MatrixView<T,I> View()
    { return MatrixView<T,I>(ptr(),M,N,Si,Sj,S,NonConj,M*N); }

    inline MatrixView<T,I> Transpose() 
    {
      return MatrixView<T,I>(ptr(),N,M,Sj,Si,TransOf(S),NonConj,M*N); 
    }

    inline MatrixView<T,I> Conjugate()
    { 
      return MatrixView<T,I>(ptr(),M,N,Si,Sj,S,
          IsReal(T())?NonConj:Conj,M*N); 
    }

    inline MatrixView<T,I> Adjoint()
    { 
      return MatrixView<T,I>(ptr(),N,M,Sj,Si,TransOf(S),
          IsReal(T())?NonConj:Conj,M*N); 
    }

    inline UpperTriMatrixView<T,I> UpperTri(
        DiagType dt=NonUnitDiag) 
    { return UpperTriMatrixView<T,I>(ptr(),N,Si,Sj,dt,S,NonConj); }

    inline LowerTriMatrixView<T,I> LowerTri(
        DiagType dt=NonUnitDiag)
    { return LowerTriMatrixView<T,I>(ptr(),M,Si,Sj,dt,S,NonConj); }

    inline VectorView<T,I> LinearView()
    { return VectorView<T,I>(ptr(),M*N,1,NonConj); }


    //
    // I/O
    //

    inline void Write(std::ostream& os) const
    { View().Write(os); }

    inline void Write(std::ostream& os, RealType(T) thresh) const
    { View().Write(os,thresh); }

    inline size_t colsize() const { return M; }
    inline size_t rowsize() const { return N; }
    inline bool IsSquare() const { return M == N; }
    inline int stepi() const { return Si; }
    inline int stepj() const { return Sj; }
    inline bool isrm() const { return S == RowMajor; }
    inline bool iscm() const { return S == ColMajor; }
    inline bool isconj() const { return false; }
    inline StorageType stor() const { return isrm() ? RowMajor : ColMajor; }
    inline ConjType ct() const { return NonConj; }
    inline size_t ls() const { return M*N; }
    inline const T* cptr() const { return itsm; }
    inline T* ptr() { return itsm; }

    inline T cref(int i, int j) const
    { return S == RowMajor ? itsm[i*N+j] : itsm[j*M+i]; }

    inline T& ref(int i, int j)
    { return S == RowMajor ? itsm[i*N+j] : itsm[j*M+i]; }

  protected :

    T itsm[M*N];

  }; // SmallMatrix

#undef Si
#undef Sj

  //---------------------------------------------------------------------------

  //
  // Copy Matrices
  //

  template <int M, int N, StorageType S1, StorageType S2, class T1, class T2> 
  struct DoCopy2 {};
  template <int M, int N, StorageType S, class T1, class T2> 
  struct DoCopy2<M,N,S,S,T1,T2>
  {
    inline DoCopy2(const T1* m1, T2* m2)
    { for(int i=0;i<M*N;++i) m2[i] = m1[i]; }
  };
  template <int M, int N, class T1, class T2> 
  struct DoCopy2<M,N,ColMajor,RowMajor,T1,T2>
  {
    inline DoCopy2(const T1* m1, T2* m2)
    { for(int i=0;i<M;++i) for(int j=0;j<N;++j) m2[i*N+j] = m1[j*M+i]; }
  };
  template <int M, int N, class T1, class T2> 
  struct DoCopy2<M,N,RowMajor,ColMajor,T1,T2>
  {
    inline DoCopy2(const T1* m1, T2* m2)
    { for(int i=0;i<M;++i) for(int j=0;j<N;++j) m2[j*M+i] = m1[i*N+j]; }
  };

  template <int M, int N, StorageType S, class T>
  struct DoCopy2<M,N,S,S,std::complex<T>,T>
  { inline DoCopy2(const std::complex<T>*, T*) { TMVAssert(FALSE); } };
  template <int M, int N, class T>
  struct DoCopy2<M,N,ColMajor,RowMajor,std::complex<T>,T>
  { inline DoCopy2(const std::complex<T>*, T*) { TMVAssert(FALSE); } };
  template <int M, int N, class T>
  struct DoCopy2<M,N,RowMajor,ColMajor,std::complex<T>,T>
  { inline DoCopy2(const std::complex<T>*, T*) { TMVAssert(FALSE); } };

  template <int M, int N, StorageType S1, StorageType S2, class T1, class T2>
  inline void DoCopy(const T1* m1, T2* m2)
  { DoCopy2<M,N,S1,S2,T1,T2>(m1,m2); }

  template <class T1, class T2, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
  inline void Copy(const SmallMatrix<T1,M,N,S1,I1>& m1, 
      SmallMatrix<T2,M,N,S2,I2>& m2)
  { 
    TMVAssert(IsComplex(T2()) || IsReal(T1()));
    DoCopy<M,N,S1,S2>(m1.cptr(),m2.ptr()); 
  }

  //
  // Swap Matrices
  //

  template <class T, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
  inline void Swap(SmallMatrix<T,M,N,S1,I1>& m1, 
      SmallMatrix<T,M,N,S2,I2>& m2) 
  {
    if (S1==S2)
      for(int i=0;i<M*N;++i) __TMV_SWAP(m1.ptr()[i],m2.ptr()[i]);
    else
      for(int i=0;i<M;++i) for(int j=0;j<N;++j) 
        __TMV_SWAP(m1.ref(i,j),m2.ref(i,j));
  }

  //
  // Functions of Matrices:
  //

  template <class T, int M, int N, StorageType S, IndexStyle I>
  inline T Det(const SmallMatrix<T,M,N,S,I>& m)
  { return m.Det(); }

  template <class T, int M, int N, StorageType S, IndexStyle I>
  inline RealType(T) LogDet(const SmallMatrix<T,M,N,S,I>& m)
  { return m.LogDet(); }

  template <class T, int M, int N, StorageType S, IndexStyle I>
  inline T Trace(const SmallMatrix<T,M,N,S,I>& m)
  { return m.Trace(); }

  template <class T, int M, int N, StorageType S, IndexStyle I>
  inline RealType(T) Norm(const SmallMatrix<T,M,N,S,I>& m)
  { return m.Norm(); }

  template <class T, int M, int N, StorageType S, IndexStyle I>
  inline RealType(T) NormSq(const SmallMatrix<T,M,N,S,I>& m)
  { return m.NormSq(); }

  template <class T, int M, int N, StorageType S, IndexStyle I>
  inline RealType(T) NormF(const SmallMatrix<T,M,N,S,I>& m)
  { return m.NormF(); }

  template <class T, int M, int N, StorageType S, IndexStyle I>
  inline RealType(T) Norm1(const SmallMatrix<T,M,N,S,I>& m)
  { return m.Norm1(); }

  template <class T, int M, int N, StorageType S, IndexStyle I>
  inline RealType(T) Norm2(const SmallMatrix<T,M,N,S,I>& m)
  { return m.Norm2(); }

  template <class T, int M, int N, StorageType S, IndexStyle I>
  inline RealType(T) NormInf(const SmallMatrix<T,M,N,S,I>& m)
  { return m.NormInf(); }

  template <class T, int M, int N, StorageType S, IndexStyle I>
  inline RealType(T) MaxAbsElement(const SmallMatrix<T,M,N,S,I>& m)
  { return m.MaxAbsElement(); }

  template <class T, int M, int N, StorageType S, IndexStyle I> 
  inline ConstMatrixView<T,I> Transpose(const SmallMatrix<T,M,N,S,I>& m)
  { return m.Transpose(); }

  template <class T, int M, int N, StorageType S, IndexStyle I> 
  inline MatrixView<T,I> Transpose(SmallMatrix<T,M,N,S,I>& m)
  { return m.Transpose(); }

  template <class T, int M, int N, StorageType S, IndexStyle I> 
  inline ConstMatrixView<T,I> Conjugate(const SmallMatrix<T,M,N,S,I>& m)
  { return m.Conjugate(); }

  template <class T, int M, int N, StorageType S, IndexStyle I> 
  inline MatrixView<T,I> Conjugate(SmallMatrix<T,M,N,S,I>& m)
  { return m.Conjugate(); }

  template <class T, int M, int N, StorageType S, IndexStyle I> 
  inline ConstMatrixView<T,I> Adjoint(const SmallMatrix<T,M,N,S,I>& m)
  { return m.Adjoint(); }

  template <class T, int M, int N, StorageType S, IndexStyle I> 
  inline MatrixView<T,I> Adjoint(SmallMatrix<T,M,N,S,I>& m)
  { return m.Adjoint(); }

  template <class T, int M, int N, StorageType S, IndexStyle I>
  inline QuotXm_1<T,T,M,N,S,I> Inverse(const SmallMatrix<T,M,N,S,I>& m)
  { return m.Inverse(); }

  //
  // Matrix ==, != Matrix
  //

  template <class T1, class T2, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
  inline bool operator==(const SmallMatrix<T1,M,N,S1,I1>& m1, 
      const SmallMatrix<T2,M,N,S2,I2>& m2)
  { 
    if (S1==S2)
      for(int i=0;i<M*N;++i) {
        if (m1.cptr()[i] != m2.cptr()[i]) return false;
      }
    else
      for(int i=0;i<M;++i) for(int j=0;j<N;++j) {
        if (m1(i,j) != m2(i,j)) return false;
      }
    return true;
  }

  template <class T1, class T2, int M, int N, StorageType S1, StorageType S2, IndexStyle I1, IndexStyle I2> 
  inline bool operator!=(const SmallMatrix<T1,M,N,S1,I1>& m1, 
      const SmallMatrix<T2,M,N,S2,I2>& m2)
  { return !(m1 == m2); }


  //
  // I/O
  //

  template <class T, int M, int N, StorageType S, IndexStyle I> 
  inline std::ostream& operator<<(std::ostream& os,
      const SmallMatrix<T,M,N,S,I>& m)
  { m.Write(os); return os; }

  template <class T, int M, int N, StorageType S, IndexStyle I> 
  inline std::istream& operator>>(std::istream& is, SmallMatrix<T,M,N,S,I>& m)
  { return is>>m.View(); }

  template <class T, int M, int N, StorageType S, IndexStyle I> 
  inline std::string TypeText(const SmallMatrix<T,M,N,S,I>& )
  { 
    std::ostringstream s;
    s << std::string("SmallMatrix<")<<TypeText(T())<<','<<M<<','<<N;
    s <<','<<Text(S)<<','<<Text(I)<<">"; 
    return s.str();
  }

} // namespace tmv

#endif
