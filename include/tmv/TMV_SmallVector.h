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


//-----------------------------------------------------------------------------
//
// This file defines the TMV SmallVector class.
//
// Constructors:
//
//    explicit SmallVector<T,N,I>()  
//        Makes a Vector of size N with _uninitialized_ values
//
//    SmallVector<T,N,I>(T x)
//        Makes a Vector of size N with all values = x
//
//    SmallVector<T,N,I>(const T* vv)
//    SmallVector<T,N,I>(const vector<T>& vv)
//    SmallVector<T,N,I>(const GenVector<T>& vv)
//        Makes a SmallVector which copies the elements of vv.
//
// 
// SmallVector doesn't have views like a regular Vector.
// All the normal viewing kinds of routines just return a regular VectorView.
// It is mostly useful for fast element access and simple tasks
// like multiplication and addition.  All the calculations are done
// inline, so the compiler can optimize the calculation for the particular
// value of N.  It is ually faster, but may take longer to compile
// depending on what calculations you are doing with them.
// 

#ifndef TMV_SmallVector_H
#define TMV_SmallVector_H

#include "tmv/TMV_Vector.h"
#include "tmv/TMV_VIt.h"
#include <sstream>
#include <algorithm>

namespace tmv {

  template <class T, int N> 
  class SmallVectorComposite;

  template <class T1, class T2, int N, IndexStyle I1, IndexStyle I2> 
  inline void Copy(const SmallVector<T1,N,I1>& v1, SmallVector<T2,N,I2>& v2);

#ifdef XTEST
#ifdef TMVDEBUG
#define XTEST_DEBUG
#endif
#endif

  template <class T, int N, IndexStyle I> 
  class SmallVector 
  {
  public:

    //
    // Constructors
    //

    inline SmallVector() 
    {
      TMVAssert(N > 0);
#ifdef TMVDEBUG
      SetAllTo(T(888));
#endif
    }

    explicit inline SmallVector(T x) 
    {
      TMVAssert(N > 0);
      if (x == T(0)) Zero();
      else SetAllTo(x); 
    }

    explicit inline SmallVector(const T* vv) 
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(N > 0);
      for(int i=0;i<N;++i) itsv[i] = vv[i];
    }

    explicit inline SmallVector(const std::vector<T>& vv) 
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(N > 0);
      for(int i=0;i<N;++i) itsv[i] = vv[i];
    }

    inline SmallVector(const SmallVector<T,N,I>& v2) 
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(N > 0);
      for(int i=0;i<N;++i) itsv[i] = v2.cref(i);
    }

    template <IndexStyle I2> 
    inline SmallVector(const SmallVector<T,N,I2>& v2) 
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(N > 0);
      TMVAssert(v2.size() == N);
      for(int i=0;i<N;++i) itsv[i] = v2.cref(i);
    }

    template <class T2, IndexStyle I2> 
    inline SmallVector(const SmallVector<T2,N,I2>& v2) 
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(N > 0);
      TMVAssert(v2.size() == N);
      Copy(v2,*this); 
    }

    template <IndexStyle I2> 
    inline SmallVector(const Vector<T,I2>& v2) 
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(N > 0);
      TMVAssert(v2.size() == N);
      for(int i=0;i<N;++i) itsv[i] = v2.cref(i);
    }

    template <class T2> 
    inline SmallVector(const GenVector<T2>& v2) 
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(N > 0);
      TMVAssert(IsReal(T2()) || IsComplex(T()));
      TMVAssert(v2.size() == N);
      View() = v2;
    }

    inline SmallVector(const AssignableToVector<RealType(T)>& v2) 
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(N > 0);
      TMVAssert(v2.size() == N);
      View() = v2;
    }

    inline SmallVector(const AssignableToVector<ComplexType(T)>& v2) 
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(N > 0);
      TMVAssert(IsComplex(T()));
      TMVAssert(v2.size() == N);
      View() = v2;
    }

    inline SmallVector(const SmallVectorComposite<RealType(T),N>& v2) 
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(N > 0);
      v2.AssignTov(*this);
    }

    inline SmallVector(const SmallVectorComposite<ComplexType(T),N>& v2) 
    {
#ifdef XTEST_DEBUG
      SetAllTo(T(888));
#endif
      TMVAssert(IsComplex(T()));
      TMVAssert(N > 0);
      v2.AssignTov(*this);
    }

    virtual inline ~SmallVector()
    {
#ifdef TMVDEBUG
      SetAllTo(T(999));
#endif
    }

    //
    // Op =
    //

    inline SmallVector<T,N,I>& operator=(SmallVector<T,N,I>& v2)
    { 
      if (&v2 != this) 
        for(int i=0;i<N;++i) itsv[i] = v2.cref(i);
      return *this; 
    }

    template <IndexStyle I2> 
    inline SmallVector<T,N,I>& operator=(SmallVector<T,N,I2>& v2)
    { 
      if (&v2 != this) 
        for(int i=0;i<N;++i) itsv[i] = v2.cref(i);
      return *this; 
    }

    template <class T2, IndexStyle I2> 
    inline SmallVector<T,N,I>& operator=(const SmallVector<T2,N,I2>& v2) 
    { 
      TMVAssert(IsReal(T2()) || IsComplex(T()));
      TMVAssert(v2.size() == N);
      Copy(v2,*this);
      return *this; 
    }

    inline SmallVector<T,N,I>& operator=(const Vector<T>& v2) 
    {
      TMVAssert(v2.size() == N);
      for(int i=0;i<N;++i) itsv[i] = v2.cref(i);
      return *this; 
    }

    template <class T2> 
    inline SmallVector<T,N,I>& operator=(const GenVector<T2>& v2) 
    {
      TMVAssert(IsReal(T2()) || IsComplex(T()));
      TMVAssert(v2.size() == N);
      View() = v2;
      return *this; 
    }

    inline SmallVector<T,N,I>& operator=(
        const AssignableToVector<RealType(T)>& v2) 
    {
      TMVAssert(v2.size() == N);
      View() = v2;
      return *this; 
    }

    inline SmallVector<T,N,I>& operator=(
        const AssignableToVector<ComplexType(T)>& v2) 
    {
      TMVAssert(IsComplex(T()));
      TMVAssert(v2.size() == N);
      View() = v2;
      return *this; 
    }

    inline SmallVector<T,N,I>& operator=(
        const SmallVectorComposite<RealType(T),N>& v2) 
    {
      v2.AssignTov(*this);
      return *this; 
    }

    inline SmallVector<T,N,I>& operator=(
        const SmallVectorComposite<ComplexType(T),N>& v2) 
    {
      TMVAssert(IsComplex(T()));
      v2.AssignTov(*this);
      return *this; 
    }


    //
    // Access Functions
    //

    typedef T value_type;
    typedef VIt<T,Unit,NonConj> iterator;
    typedef CVIt<T,Unit,NonConj> const_iterator;
    typedef VIt<T,Step,NonConj> reverse_iterator;
    typedef CVIt<T,Step,NonConj> const_reverse_iterator;
    typedef T& reference;

    inline const_iterator begin() const 
    { return const_iterator(itsv,1); }
    inline const_iterator end() const 
    { return begin()+N; }
    inline const_reverse_iterator rbegin() const 
    { return const_reverse_iterator(itsv+(N-1),-1); }
    inline const_reverse_iterator rend() const 
    { return rbegin()+N; }

    ListAssigner<T,iterator> inline operator=(ListInitClass)
    { return ListAssigner<T,iterator>(begin(),size()); }

    inline T operator[](int i) const 
    { 
      if (I == CStyle) {
        TMVAssert(i>=0 && i<N);
        return cref(i); 
      } else {
        TMVAssert(i>=1 && i<=N);
        return cref(i-1); 
      }
    }
    inline T operator()(int i) const 
    { 
      if (I == CStyle) {
        TMVAssert(i>=0 && i<N);
        return cref(i); 
      } else {
        TMVAssert(i>=1 && i<=N);
        return cref(i-1); 
      }
    }

    inline iterator begin() 
    { return iterator(itsv,1 FIRSTLAST ); }
    inline iterator end() 
    { return begin() + N; }

    inline reverse_iterator rbegin() 
    { return reverse_iterator(itsv+N-1,-1 FIRSTLAST ); }
    inline reverse_iterator rend() 
    { return rbegin()+N; }

    inline T& operator[](int i) 
    { 
      if (I == CStyle) {
        TMVAssert(i>=0 && i<N);
        return ref(i); 
      } else {
        TMVAssert(i>=1 && i<=N);
        return ref(i-1); 
      }
    }
    inline T& operator()(int i) 
    { 
      if (I == CStyle) {
        TMVAssert(i>=0 && i<N);
        return ref(i); 
      } else {
        TMVAssert(i>=1 && i<=N);
        return ref(i-1); 
      }
    }

    //
    // Functions of Vector
    //

    inline RealType(T) Norm() const // = Norm2
    { return Norm2(); }

    inline RealType(T) NormSq(RealType(T) scale = RealType(T)(1)) const
    {
      RealType(T) sum(0);
      if (scale == RealType(T)(1))
        for(int i=0;i<N;++i) sum += NORM(itsv[i]);
      else
        for(int i=0;i<N;++i) sum += NORM(scale*itsv[i]);
      return sum;
    }

    inline RealType(T) Norm1() const // sum_i |v_i|
    { return SumAbsElements(); }

    inline RealType(T) Norm2() const // sqrt( sum_i |v_i|^2 )
    { return SQRT(NormSq()); }

    inline RealType(T) NormInf() const // max_i |v_i|
    { return size() > 0 ? MaxAbsElement() : RealType(T)(0); }

    inline T SumElements() const
    {
      T sum(0);
      for(int i=0;i<N;++i) sum += itsv[i];
      return sum;
    }

    inline RealType(T) SumAbsElements() const
    {
      RealType(T) sum(0);
      for(int i=0;i<N;++i) sum += ABS(itsv[i]);
      return sum;
    }

    inline T MinElement(int* iminout=0) const
    {
      T min = N>0 ? itsv[0] : T(0);
      if (iminout) *iminout = 0;
      for(int i=1;i<N;++i) {
        if (REAL(itsv[i]) < REAL(min)) {
          min = itsv[i];
          if (iminout) *iminout = i;
        }
      }
      if (I == FortranStyle && iminout) ++(*iminout);
      return min;
    }

    inline T MaxElement(int* imaxout=0) const
    {
      T max = N>0 ? itsv[0] : T(0);
      if (imaxout) *imaxout = 0;
      for(int i=1;i<N;++i) {
        if (REAL(itsv[i]) > REAL(max)) {
          max = itsv[i];
          if (imaxout) *imaxout = i;
        }
      }
      if (I == FortranStyle && imaxout) ++(*imaxout);
      return  max;
    }

    inline RealType(T) MinAbsElement(int* iminout=0) const
    {
      RealType(T) min = N>0 ? 
      (IsReal(T()) ? ABS(itsv[0]) : NORM(itsv[0])) : RealType(T)(0);
      if (iminout) *iminout = 0;
      for(int i=1;i<N;++i) {
        RealType(T) absvi = IsReal(T()) ? ABS(itsv[i]) : NORM(itsv[i]);
        if (absvi < min) {
          min = absvi;
          if (iminout) *iminout = i;
        }
      }
      if (I == FortranStyle && iminout) ++(*iminout);
      return IsReal(T()) ? min : SQRT(min);
    }

    inline RealType(T) MaxAbsElement(int* imaxout=0) const
    {
      RealType(T) max = N>0 ? 
      (IsReal(T()) ? ABS(itsv[0]) : NORM(itsv[0])) : RealType(T)(0);
      if (imaxout) *imaxout = 0;
      for(int i=1;i<N;++i) {
        RealType(T) absvi = IsReal(T()) ? ABS(itsv[i]) : NORM(itsv[i]);
        if (absvi > max) {
          max = absvi;
          if (imaxout) *imaxout = i;
        }
      }
      if (I == FortranStyle && imaxout) ++(*imaxout);
      return IsReal(T()) ? max : SQRT(max);
    }


    //
    // Modifying Functions
    //

    inline SmallVector<T,N,I>& Zero() 
    { 
      for(int i=0;i<N;++i) itsv[i] = T(0);
      return *this;
    }

    inline SmallVector<T,N,I>& Clip(RealType(T) thresh)
    {
      for(int i=0; i<N; ++i)
        if (std::abs(itsv[i]) < thresh) itsv[i] = T(0);
      return *this;
    }

    inline SmallVector<T,N,I>& SetAllTo(T x)
    { 
      for(int i=0;i<N;++i) itsv[i] = x;
      return *this;
    }

    inline SmallVector<T,N,I>& AddToAll(T x)
    {
      for(int i=0;i<N;++i) itsv[i] += x;
      return *this;
    }

    inline SmallVector<T,N,I>& ConjugateSelf()
    { 
      if (IsComplex(T())) {
        RealType(T)* itsvi = reinterpret_cast<RealType(T)*>(itsv)+1;
        for(int i=0;i<2*N;i+=2) itsvi[i] = -itsvi[i];
      }
      return *this; 
    }

    inline SmallVector<T,N,I>& MakeBasis(int i)
    { 
      if (I == CStyle) {
        TMVAssert(i>=0 && i<N);
        Zero(); itsv[i] = T(1);
      } else {
        TMVAssert(i>=1 && i<=N);
        Zero(); itsv[i-1] = T(1);
      }
      return *this; 
    }

    inline SmallVector<T,N,I>& Swap(int i1, int i2)
    {
      if (I == CStyle) {
        TMVAssert(i1>=0 && i1<N);
        TMVAssert(i2>=0 && i2<N);
        if (i1 != i2) __TMV_SWAP(itsv[i1],itsv[i2]);
      } else {
        TMVAssert(i1>=1 && i1<=N);
        TMVAssert(i2>=1 && i2<=N);
        if (i1 != i2) __TMV_SWAP(itsv[i1-1],itsv[i2-1]);
      }
      return *this;
    }

    inline SmallVector<T,N,I>& Permute(const int* p, int i1, int i2) 
    {
      if (I == CStyle) {
        TMVAssert(i1>=0 && i1<=i2 && i2<=N);
        for(int i=i1;i<i2;++i) Swap(i,p[i]);
      } else {
        TMVAssert(i1>=1 && i1<=i2 && i2<=N);
        for(int i=i1-1;i<i2;++i) Swap(i,p[i]);
      }
      return *this;
    }

    inline SmallVector<T,N,I>& Permute(const int* p) 
    { return Permute(p,I==CStyle?0:1,N); }

    inline SmallVector<T,N,I>& ReversePermute(const int* p, 
        int i1, int i2)
    {
      if (I == CStyle) {
        TMVAssert(i1>=0 && i1<=i2 && i2<=N);
        for(int i=i2-1;i>=i1;--i) Swap(i,p[i]);
      } else {
        TMVAssert(i1>=1 && i1<=i2 && i2<=N);
        for(int i=i2-1;i>=i1-1;--i) Swap(i,p[i]);
      }
      return *this;
    }

    inline SmallVector<T,N,I>& ReversePermute(const int* p)
    { return ReversePermute(p,I==CStyle?0:1,N); }

    inline SmallVector<T,N,I>& ReverseSelf()
    {
      for(int i1=0,i2=N-1;i1<i2;++i1,--i2) __TMV_SWAP(itsv[i1],itsv[i2]);
      return *this;
    }

    inline SmallVector<T,N,I>& Sort(
        int* P, ADType ad=ASCEND, COMPType comp=REAL_COMP)
    { View().Sort(P,ad,comp); return *this; }

    //
    // SubVector
    //

    inline ConstVectorView<T,I> SubVector(int i1, int i2) const
    { return View().SubVector(i1,i2); }

    inline VectorView<T,I> SubVector(int i1, int i2)
    { return View().SubVector(i1,i2); }

    inline ConstVectorView<T,I> SubVector(int i1, int i2, int istep) const
    { return View().SubVector(i1,i2,istep); }

    inline VectorView<T,I> SubVector(int i1, int i2, int istep)
    { return View().SubVector(i1,i2,istep); }

    inline ConstVectorView<T,I> Reverse() const
    { return View().Reverse(); }

    inline VectorView<T,I> Reverse()
    { return View().Reverse(); }

    inline ConstVectorView<T,I> View() const
    { return ConstVectorView<T,I>(itsv,N,1,NonConj); }

    inline VectorView<T,I> View()
    { return VectorView<T,I>(itsv,N,1,NonConj); }

    inline ConstVectorView<T,I> Conjugate() const
    { return ConstVectorView<T,I>(itsv,N,1,IsReal(T())?NonConj:Conj); }

    inline VectorView<T,I> Conjugate()
    { return VectorView<T,I>(itsv,N,1,IsReal(T())?NonConj:Conj); }

    inline ConstVectorView<RealType(T),I> Real() const
    { return View().Real(); }

    inline ConstVectorView<RealType(T),I> Imag() const
    { return View().Imag(); }

    inline ConstVectorView<RealType(T),I> Flatten() const
    { return View().Flatten(); }

    inline VectorView<RealType(T),I> Real()
    { return View().Real(); }

    inline VectorView<RealType(T),I> Imag()
    { return View().Imag(); }

    inline VectorView<RealType(T),I> Flatten()
    { return View().Flatten(); }

    // 
    // I/O
    //

    inline void Write(std::ostream& os) const
    { View().Write(os); }

    inline void Write(std::ostream& os, RealType(T) thresh) const
    { View().Write(os,thresh); }

    inline size_t size() const { return N; }
    inline const T* cptr() const { return itsv; }
    inline T* ptr() { return itsv; }
    inline int step() const { return 1; }
    inline ConjType ct() const { return NonConj; }
    inline bool isconj() const { return false; }

    inline T cref(int i) const
    { return itsv[i]; }

    inline T& ref(int i)
    { return itsv[i]; }

  protected :

    T itsv[N];

  }; // SmallVector


  //
  // Copy SmallVectors
  //

  template <int N, class T1, class T2> 
  struct DoCopy1
  {
    DoCopy1(const T1* v1, T2* v2)
    { for(int i=0;i<N;++i) v2[i] = v1[i]; }
  };

  template <int N, class T> 
  struct DoCopy1<N,std::complex<T>,T>
  {
    DoCopy1(const std::complex<T>* , T* )
    { TMVAssert(FALSE); }
  };

  template <int N, class T1, class T2> 
  inline void DoCopy(const T1* v1, T2* v2)
  { DoCopy1<N,T1,T2>(v1,v2); }

  template <class T1, class T2, int N, IndexStyle I1, IndexStyle I2> 
  inline void Copy(const SmallVector<T1,N,I1>& v1, SmallVector<T2,N,I2>& v2)
  { DoCopy<N>(v1.cptr(),v2.ptr()); }

  //
  // Swap SmallVectors
  //

  template <class T, int N, IndexStyle I1, IndexStyle I2> 
  inline void Swap(SmallVector<T,N,I1>& v1, SmallVector<T,N,I2>& v2)
  {
    for(int i=0;i<N;++i) __TMV_SWAP(v1.ref(i),v2.ref(i));
  }

  //
  // Functions of Vectors
  //

  template <class T, int N, IndexStyle I> 
  inline RealType(T) Norm(const SmallVector<T,N,I>& v)
  { return v.Norm(); }

  template <class T, int N, IndexStyle I> 
  inline RealType(T) Norm1( const SmallVector<T,N,I>& v)
  { return v.Norm1(); }

  template <class T, int N, IndexStyle I> 
  inline RealType(T) NormSq( const SmallVector<T,N,I>& v)
  { return v.NormSq(); }

  template <class T, int N, IndexStyle I> 
  inline RealType(T) Norm2( const SmallVector<T,N,I>& v)
  { return v.Norm2(); }

  template <class T, int N, IndexStyle I> 
  inline RealType(T) NormInf( const SmallVector<T,N,I>& v)
  { return v.NormInf(); }

  template <class T, int N, IndexStyle I> 
  inline T SumElements( const SmallVector<T,N,I>& v)
  { return v.SumElements(); }

  template <class T, int N, IndexStyle I> 
  inline RealType(T) SumAbsElements(const SmallVector<T,N,I>& v)
  { return v.SumAbsElements(); }

  template <class T, int N, IndexStyle I> 
  inline T MinElement(const SmallVector<T,N,I>& v)
  { return v.MinElement(); }

  template <class T, int N, IndexStyle I> 
  inline T MaxElement(const SmallVector<T,N,I>& v)
  { return v.MaxElement(); }

  template <class T, int N, IndexStyle I> 
  inline RealType(T) MinAbsElement(const SmallVector<T,N,I>& v)
  { return v.MinAbsElement(); }

  template <class T, int N, IndexStyle I> 
  inline RealType(T) MaxAbsElement(const SmallVector<T,N,I>& v)
  { return v.MaxAbsElement(); }

  template <class T, int N, IndexStyle I> 
  inline VectorView<T,I> Conjugate(const SmallVector<T,N,I>& v)
  { return v.Conjugate(); }

  template <class T, int N, IndexStyle I> 
  inline VectorView<T,I> Conjugate(SmallVector<T,N,I>& v)
  { return v.Conjugate(); }


  //
  // Vector ==, != Vector
  //

  template <class T1, class T2, int N, IndexStyle I1, IndexStyle I2> 
  inline bool operator==(
      const SmallVector<T1,N,I1>& v1, const SmallVector<T2,N,I2>& v2)
  {
    for(int i=0;i<N;++i) if (v1.cref(i) != v2.cref(i)) return false;
    return true;
  }

  template <class T1, class T2, int N, IndexStyle I1, IndexStyle I2> 
  inline bool operator!=(
      const SmallVector<T1,N,I1>& v1, const SmallVector<T2,N,I2>& v2)
  { return !(v1 == v2); }

  template <class T1, class T2, int N, IndexStyle I> 
  inline bool operator==(
      const GenVector<T1>& v1, const SmallVector<T2,N,I>& v2)
  { return v1 == v2.View(); }

  template <class T1, class T2, int N, IndexStyle I> 
  inline bool operator==(
      const SmallVector<T1,N,I>& v1, const GenVector<T2>& v2)
  { return v1.View() == v2; }

  template <class T1, class T2, int N, IndexStyle I> 
  inline bool operator!=(
      const GenVector<T1>& v1, const SmallVector<T2,N,I>& v2)
  { return v1 != v2.View(); }

  template <class T1, class T2, int N, IndexStyle I> 
  inline bool operator!=(
      const SmallVector<T1,N,I>& v1, const GenVector<T2>& v2)
  { return v1.View() != v2; }


  //
  // I/O
  //

  template <class T, int N, IndexStyle I> 
  inline std::ostream& operator<<(
      std::ostream& os, const SmallVector<T,N,I>& v)
  { v.Write(os); return os; }

  template <class T, int N, IndexStyle I> 
  inline std::istream& operator>>(
      std::istream& is, SmallVector<T,N,I>& v)
  { return is >> v.View(); }

  template <class T, int N, IndexStyle I> 
  inline std::string TypeText(const SmallVector<T,N,I>& )
  { 
    std::ostringstream s;
    s << "SmallVector<"<<TypeText(T())<<","<<N<<","<<Text(I)<<">";
    return s.str();
  }

} // namespace tmv

#endif
