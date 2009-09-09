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


//-----------------------------------------------------------------------------
//
// This file defines the TMV SmallVector class.
//
// Constructors:
//
//    explicit SmallVector<T,N>()  
//        Makes a Vector of size N with _uninitialized_ values
//
//    SmallVector<T,N>(T x)
//        Makes a Vector of size N with all values = x
//
//    SmallVector<T,N>(const T* vv)
//    SmallVector<T,N>(const vector<T>& vv)
//        Makes a SmallVector which copies the elements of vv.
//        vv must have length N.
//
//    SmallVectorViewOf<N>(T* vv)
//    SmallVectorViewOf<N>(const T* vv)
//        Makes a SmallVectorView which refers to the exact
//        elements of vv, not copying them to new storage.
//        The first one returns a SmallVectorView, 
//        the second a ConstSmallVectorView.
//
// 
// A SmallVector _is a_ GenVector.  ie. it inherits from GenVector.  
// So everything you can do with a Vector, you can also do with a SmallVector.
// The advantage is that many access functions are doable at compile 
// time, since N is a template parameter.
//
// An advantage that may also be a disadvantage is that every calculation 
// is inline, since it doesn't make sense to instantiate every possible N 
// that someone might want.  This can speed up some calculations.  But
// it may take significantly longer to compile, depending on what 
// calculations you are doing with them.
// 

#ifndef TMV_SmallVector_H
#define TMV_SmallVector_H

#include "TMV_Vector.h"
#include "TMV_VIt.h"
#include <sstream>
#include <algorithm>

namespace tmv {

  template <class T, int N, int S, bool C> class GenSmallVector;
  template <class T, int N, int S=1, bool C=false, IndexStyle I=CStyle> 
    class ConstSmallVectorView;
  template <class T, int N, int S=1, bool C=false, IndexStyle I=CStyle>
    class SmallVectorView;
  template <class T, int N, IndexStyle I=CStyle> class SmallVector;
  template <class T, int N> class SmallVectorComposite;

  template <class T1, class T2, int N, int S1, int S2, bool C1>
    inline void Copy(
	const GenSmallVector<T1,N,S1,C1>& v1, 
	const SmallVectorView<T2,N,S2,false>& v2);
  template <class T1, class T2, int N, int S1, int S2, bool C1>
    inline void Copy(
	const GenSmallVector<T1,N,S1,C1>& v1, 
	const SmallVectorView<T2,N,S2,true>& v2)
    { Copy(v1.Conjugate(),v2.Conjugate()); }
  template <class T1, class T2, int N, int S1, bool C1>
    inline void Copy(
	const GenSmallVector<T1,N,S1,C1>& v1, T2* v2p, int S2);
#define CT std::complex<T>
  template <class T, int N, int S1, int S2, bool C1>
    inline void Copy(
	const GenSmallVector<CT,N,S1,C1>& , 
	const SmallVectorView<T,N,S2,false>& )
    { TMVAssert(FALSE); }
  template <class T, int N, int S1, bool C1>
    inline void Copy(
	const GenSmallVector<CT,N,S1,C1>& , T* , int )
    { TMVAssert(FALSE); }
#undef CT

  template <class T, int N> struct AssignableToSmallVector :
    virtual public AssignableToVector<T>
  {
    virtual void AssignTov(
	const SmallVectorView<RealType(T),N,1,false>& v2) const = 0;
    inline void AssignTov(
	const SmallVectorView<RealType(T),N,1,true>& v2) const
    { AssignTov(v2.Conjugate()); }
    virtual void AssignTov(
	const SmallVectorView<ComplexType(T),N,1,false>& v2) const = 0;
    virtual void AssignTov(
	const SmallVectorView<ComplexType(T),N,1,true>& v2) const = 0;

    virtual void DoAssignTov(RealType(T)* v2p, const int S2) const = 0;
    virtual void DoAssignTov(ComplexType(T)* v2p, 
	const int S2, const bool C2) const = 0;
    template <int S2, bool C2> inline void AssignTov(
	const SmallVectorView<RealType(T),N,S2,C2>& v2) const
    { DoAssignTov(v2.ptr(),S2); }
    template <int S2, bool C2> inline void AssignTov(
	const SmallVectorView<ComplexType(T),N,S2,C2>& v2) const
    { DoAssignTov(v2.ptr(),S2,C2); }
  };

  template <int N, IndexStyle I, class T> 
    inline SmallVectorView<T,N,1,false,I> SmallVectorViewOf(T* v)
    { return SmallVectorView<T,N,1,false,I>(v FIRSTLAST1(v,v+N)); }

  template <int N, IndexStyle I, class T> 
    inline ConstSmallVectorView<T,N,1,false,I> SmallVectorViewOf(const T* v)
    { return ConstSmallVectorView<T,N,1,false,I>(v); }

  template <int N, class T> 
    inline SmallVectorView<T,N,1,false> SmallVectorViewOf(T* v)
    { return SmallVectorView<T,N,1,false>(v FIRSTLAST1(v,v+N)); }

  template <int N, class T> 
    inline ConstSmallVectorView<T,N,1,false> SmallVectorViewOf(const T* v)
    { return ConstSmallVectorView<T,N,1,false>(v); }

#define GENBASE

  template <class T, int N, int S, bool C> class GenSmallVector : 
#ifdef GENBASE
    public GenVector<T>,
#endif
    virtual public AssignableToSmallVector<T,N>
  {
    public:

      //
      // Constructor
      //

      inline GenSmallVector() {}
      inline GenSmallVector(const GenSmallVector<T,N,S,C>& ) {}
      virtual inline ~GenSmallVector() {}

      //
      // Access Functions
      //

      inline size_t size() const { return N; }

      inline void AssignTov(const SmallVectorView<RealType(T),
	  N,1,false>& v2) const
      { TMVAssert(IsReal(T())); Copy(*this,v2); }
      inline void AssignTov(const SmallVectorView<ComplexType(T),
	  N,1,false>& v2) const
      { Copy(*this,v2); }
      inline void AssignTov(const SmallVectorView<ComplexType(T),
	  N,1,true>& v2) const
      { Copy(*this,v2); }
      inline void DoAssignTov(RealType(T)* v2p, const int S2) const
      { 
	TMVAssert(IsReal(T()));
	if (S2 == 1) Copy(*this,SmallVectorViewOf<N>(v2p));
	else Copy(*this,v2p,S2);
      }
      inline void DoAssignTov(ComplexType(T)* v2p, 
	  const int S2, const bool C2) const
      { 
	if (C2) Conjugate().DoAssignTov(v2p,S,false);
	else if (S2 == 1) Copy(*this,SmallVectorViewOf<N>(v2p));
	else Copy(*this,v2p,S2);
      }
      inline void AssignToV(const VectorView<RealType(T)>& v2) const
      {
	TMVAssert(IsReal(T()));
	TMVAssert(v2.size() == N);
	DoAssignTov(v2.ptr(),v2.step()); 
      }
      inline void AssignToV(const VectorView<ComplexType(T)>& v2) const
      { 
	TMVAssert(v2.size() == N);
	DoAssignTov(v2.ptr(),v2.step(),v2.isconj()); 
      }

      inline T operator[](int i) const 
      { 
	TMVAssert(i>=0 && i<N);
	return cref(i); 
      }
      inline T operator()(int i) const 
      {
	TMVAssert(i>=0 && i<N);
	return cref(i); 
      }

      typedef CVIt<T,S==1?Unit:Step,C?Conj:NonConj> const_iterator;
      typedef CVIt<T,-S==1?Unit:Step,C?Conj:NonConj> const_reverse_iterator;

      inline const_iterator begin() const 
      { return const_iterator(cptr(),S); }
      inline const_iterator end() const 
      { return begin()+N; }
      inline const_reverse_iterator rbegin() const 
      { return const_reverse_iterator(cptr()+S*(N-1),-S); }
      inline const_reverse_iterator rend() const 
      { return rbegin()+N; }

      //
      // SubVector
      //

      inline ConstVectorView<T> SubVector(int i1, int i2) const
      {
	TMVAssert(RegView().OKSubVector(i1,i2,1));
	return RegView().SubVector(i1,i2);
      }
      inline ConstVectorView<T> SubVector(int i1, int i2, int istep) const
      {
	TMVAssert(RegView().OKSubVector(i1,i2,istep));
	return RegView().SubVector(i1,i2,istep);
      }

      inline ConstSmallVectorView<T,N,-S,C> Reverse() const
      { return ConstSmallVectorView<T,N,-S,C>(cptr()+(N-1)*S); }

      inline ConstSmallVectorView<T,N,S,C> View() const
      { return ConstSmallVectorView<T,N,S,C>(cptr()); }

      inline ConstSmallVectorView<T,N,S,!C> Conjugate() const
      { return ConstSmallVectorView<T,N,S,!C>(cptr()); }

      inline ConstVectorView<T> RegView() const
      { return ConstVectorView<T>(cptr(),N,S,C?Conj:NonConj); }


      //
      // Functions of Vector
      //

      inline RealType(T) Norm() const // = Norm2
      { return Norm2(); }

      inline RealType(T) NormSq(RealType(T) scale = RealType(T)(1)) const
      {
	RealType(T) sum(0);
	const T* itsv = cptr();
	for(int i=0,ii=0;i<N;++i,ii+=S) sum += NORM(scale*itsv[ii]);
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
	const T* itsv = cptr();
	for(int i=0,ii=0;i<N;++i,ii+=S) sum += itsv[ii];
	return C ? CONJ(sum) : sum;
      }

      inline RealType(T) SumAbsElements() const
      {
	RealType(T) sum(0);
	const T* itsv = cptr();
	for(int i=0,ii=0;i<N;++i,ii+=S) sum += ABS(itsv[ii]);
	return sum;
      }

      inline T MinElement(int* iminout=0) const
      {
	const T* itsv = cptr();
	T min = N>0 ? itsv[0] : T(0);
	if (iminout) *iminout = 0;
	for(int i=1,ii=S;i<N;++i,ii+=S) {
	  if (REAL(itsv[ii]) < REAL(min)) {
	    min = itsv[ii];
	    if (iminout) *iminout = i;
	  }
	}
	return C ? CONJ(min) : min;
      }

      inline T MaxElement(int* imaxout=0) const
      {
	const T* itsv = cptr();
	T max = N>0 ? itsv[0] : T(0);
	if (imaxout) *imaxout = 0;
	for(int i=1,ii=S;i<N;++i,ii+=S) {
	  if (REAL(itsv[ii]) > REAL(max)) {
	    max = itsv[ii];
	    if (imaxout) *imaxout = i;
	  }
	}
	return  C ? CONJ(max) :max;
      }

      inline RealType(T) MinAbsElement(int* iminout=0) const
      {
	const T* itsv = cptr();
	RealType(T) min = N>0 ? 
	  (IsReal(T()) ? ABS(itsv[0]) : NORM(itsv[0])) : RealType(T)(0);
	if (iminout) *iminout = 0;
	for(int i=1,ii=S;i<N;++i,ii+=S) {
	  RealType(T) absvi = IsReal(T()) ? ABS(itsv[ii]) : NORM(itsv[ii]);
	  if (absvi < min) {
	    min = absvi;
	    if (iminout) *iminout = i;
	  }
	}
	return IsReal(T()) ? min : SQRT(min);
      }

      inline RealType(T) MaxAbsElement(int* imaxout=0) const
      {
	const T* itsv = cptr();
	RealType(T) max = N>0 ? 
	  (IsReal(T()) ? ABS(itsv[0]) : NORM(itsv[0])) : RealType(T)(0);
	if (imaxout) *imaxout = 0;
	for(int i=1,ii=S;i<N;++i,ii+=S) {
	  RealType(T) absvi = IsReal(T()) ? ABS(itsv[ii]) : NORM(itsv[ii]);
	  if (absvi > max) {
	    max = absvi;
	    if (imaxout) *imaxout = i;
	  }
	}
	return IsReal(T()) ? max : SQRT(max);
      }


      // 
      // I/O
      //

      inline void Write(std::ostream& os) const
      {
	os << N << " ( ";
	const T* itsv = cptr();
	for(int i=0,ii=0;i<N;++i,ii+=S) os << " " << itsv[ii] << " ";
	os << " )";
      }

      virtual const T* cptr() const =0;

      inline int step() const { return S; }
      inline ConjType ct() const { return isconj() ? Conj : NonConj; }
      inline bool isconj() const { return C && IsComplex(T()); }

    protected:

      inline T cref(int i) const
      {
	TMVAssert(i>=0 && i < N);
	const T* vi = S > 0 ? cptr() + i*S : S < 0 ? cptr() - i*(-S) : cptr();
	return C ? CONJ(*vi) : *vi;
      }

    private:

      inline GenSmallVector<T,N,S,C>& operator=(
	  const GenSmallVector<T,N,S,C>&) 
      { TMVAssert(FALSE); return *this; }

  }; // GenSmallVector

  template <class T, int N, int S, bool C, IndexStyle I> 
    class ConstSmallVectorView : 
    public GenSmallVector<T,N,S,C>
  {
    public:

      inline ConstSmallVectorView(const ConstSmallVectorView<T,N,S,C,I>& rhs) : 
	itsv(rhs.itsv) {}
      inline ConstSmallVectorView(const GenSmallVector<T,N,S,C>& rhs) : 
	itsv(rhs.cptr()) {}
      inline ConstSmallVectorView(const GenVector<T>& rhs) : 
	itsv(rhs.cptr()) 
      {
	TMVAssert(rhs.size() == N);
	TMVAssert(rhs.step() == S);
	TMVAssert(rhs.isconj() == C);
      }
      inline ConstSmallVectorView(const T* inv) : itsv(inv) {}
      virtual inline ~ConstSmallVectorView() {}

      inline const T* cptr() const { return itsv; }

    private:

      const T*const itsv;

      inline ConstSmallVectorView<T,N,S,C,I>& operator=(
	  const ConstSmallVectorView<T,N,S,C,I>&) 
      { TMVAssert(FALSE); return *this; }

  }; // ConstSmallVectorView

  template <class T, int N, int S, bool C> 
    class ConstSmallVectorView<T,N,S,C,FortranStyle> : 
    public ConstSmallVectorView<T,N,S,C,CStyle>
    {
      public:

	inline ConstSmallVectorView(
	    const ConstSmallVectorView<T,N,S,C,FortranStyle>& rhs) : 
	  ConstSmallVectorView<T,N,S,C,CStyle>(rhs) {}
	inline ConstSmallVectorView(
	    const ConstSmallVectorView<T,N,S,C,CStyle>& rhs) : 
	  ConstSmallVectorView<T,N,S,C,CStyle>(rhs) {}
	inline ConstSmallVectorView(const GenSmallVector<T,N,S,C>& rhs) : 
	  ConstSmallVectorView<T,N,S,C,CStyle>(rhs) {}
	inline ConstSmallVectorView(const GenVector<T>& rhs) : 
	  ConstSmallVectorView<T,N,S,C,CStyle>(rhs) {}
	inline ConstSmallVectorView(const T* inv) :
	  ConstSmallVectorView<T,N,S,C,CStyle>(inv) {}
	virtual inline ~ConstSmallVectorView() {}

	//
	// Access Functions
	//

	inline T operator[](int i) const 
	{
	  TMVAssert(i>0 && i<=N);
	  return GenSmallVector<T,N,S,C>::cref(i-1); 
	}
	inline T operator()(int i) const 
	{
	  TMVAssert(i>0 && i<=N);
	  return GenSmallVector<T,N,S,C>::cref(i-1); 
	}

	//
	// SubVector
	//

	inline ConstVectorView<T,FortranStyle> SubVector(int i1, int i2) const
	{
	  TMVAssert(RegView().OKSubVector(i1,i2,1));
	  return GenVector<T>::SubVector(i1-1,i2);
	}

	inline ConstVectorView<T,FortranStyle> SubVector(
	    int i1, int i2, int istep) const
	{
	  TMVAssert(RegView().OKSubVector(i1,i2,istep));
	  return GenVector<T>::SubVector(i1-1,i2-1+istep,istep);
	}

	inline ConstSmallVectorView<T,N,-S,C,FortranStyle> Reverse() const
	{ return GenSmallVector<T,N,S,C>::Reverse(); }

	inline ConstSmallVectorView<T,N,S,C,FortranStyle> View() const
	{ return GenSmallVector<T,N,S,C>::View(); }

	inline ConstSmallVectorView<T,N,S,!C,FortranStyle> Conjugate() const
	{ return GenSmallVector<T,N,S,C>::Conjugate(); }

	inline ConstVectorView<T,FortranStyle> RegView() const
	{ return GenSmallVector<T,N,S,C>::RegView(); }

	inline ConstVectorView<RealType(T),FortranStyle> Real() const
	{ return GenSmallVector<T,N,S,C>::Real(); }

	inline ConstVectorView<RealType(T),FortranStyle> Imag() const
	{ return GenSmallVector<T,N,S,C>::Imag(); }

	inline ConstVectorView<RealType(T),FortranStyle> Flatten() const
	{ return GenSmallVector<T,N,S,C>::Flatten(); }

	inline T MinElement(int* iminout=0) const
	{
	  T ret = View().MinElement(iminout);
	  if (iminout) (*iminout)--;  
	  return ret;
	}
	inline T MaxElement(int* imaxout=0) const
	{
	  T ret = View().MaxElement(imaxout);
	  if (imaxout) (*imaxout)--;  
	  return ret;
	}
	inline RealType(T) MinAbsElement(int* iminout=0) const
	{
	  RealType(T) ret = View().MinElement(iminout);
	  if (iminout) (*iminout)--;  
	  return ret;
	}
	inline RealType(T) MaxAbsElement(int* imaxout=0) const
	{
	  RealType(T) ret = View().MaxElement(imaxout);
	  if (imaxout) --(*imaxout);  
	  return ret;
	}

      private:

	inline ConstSmallVectorView<T,N,S,C,FortranStyle>& operator=(
	    const ConstSmallVectorView<T,N,S,C,FortranStyle>&) 
	{ TMVAssert(FALSE); return *this; }

    }; // FortranStyle ConstSmallVectorView

  template <class T, int N, int S, bool C, IndexStyle I> 
    class SmallVectorView : 
    public GenSmallVector<T,N,S,C>
  {
    public:

      //
      // Constructors 
      //

      inline SmallVectorView(const SmallVectorView<T,N,S,C,I>& rhs) : 
	itsv(rhs.itsv) DEFFIRSTLAST(rhs.first,rhs.last) {}

      inline SmallVectorView(T* inv PARAMFIRSTLAST(T) ) :
	itsv(inv) DEFFIRSTLAST(_first,_last) {}

      inline SmallVectorView(const VectorView<T>& rhs) : 
	itsv(rhs.ptr()) DEFFIRSTLAST(rhs.first,rhs.last)
      {
	TMVAssert(rhs.size() == N);
	TMVAssert(rhs.step() == S);
	TMVAssert(rhs.isconj() == C);
      }

      virtual inline ~SmallVectorView() {}

      inline operator VectorView<T,I>() const
      { return VectorView<T,I>(ptr(),N,S,C?Conj:NonConj FIRSTLAST); }

      //
      // Op =
      //

      inline const SmallVectorView<T,N,S,C,I>& operator=(
	  const SmallVectorView<T,N,S,C,I>& v2) const
      { 
	Copy(v2,*this); 
	return *this; 
      }

      inline const SmallVectorView<T,N,S,C,I>& operator=(
	  const GenSmallVector<RealType(T),N,S,C>& v2) const
      { 
	v2.AssignTov(*this);
	return *this; 
      }

      inline const SmallVectorView<T,N,S,C,I>& operator=(
	  const GenSmallVector<ComplexType(T),N,S,C>& v2) const
      { 
	TMVAssert(IsComplex(T()));
	v2.AssignTov(*this);
	return *this; 
      }

      template <class T2, int S2, bool C2> 
	inline const SmallVectorView<T,N,S,C,I>& operator=(
	  const GenSmallVector<T2,N,S2,C2>& v2) const
	{ 
	  TMVAssert(IsReal(T2()) || IsComplex(T()));
	  Copy(v2,*this);
	  return *this; 
	}

      inline const SmallVectorView<T,N,S,C,I>& operator=(
	  const GenVector<RealType(T)>& v2) const
      { 
	TMVAssert(v2.size() == N);
	v2.AssignToV(RegView());
	return *this; 
      }

      inline const SmallVectorView<T,N,S,C,I>& operator=(
	  const GenVector<ComplexType(T)>& v2) const
      { 
	TMVAssert(v2.size() == N);
	TMVAssert(IsComplex(T()));
	v2.AssignToV(RegView());
	return *this; 
      }

      template <class T2> inline const SmallVectorView<T,N,S,C,I>& operator=(
	  const GenVector<T2>& v2) const
      { 
	TMVAssert(v2.size() == N);
	v2.AssignToV(RegView());
	return *this; 
      }

      inline const SmallVectorView<T,N,S,C,I>& operator=(
	  const AssignableToSmallVector<RealType(T),N>& v2) const
      { v2.AssignTov(*this); return *this; }

      inline const SmallVectorView<T,N,S,C,I>& operator=(
	  const AssignableToSmallVector<ComplexType(T),N>& v2) const
      { v2.AssignTov(*this); return *this; }

      //
      // Access Functions
      //
      
      typedef VIt<T,S==1?Unit:Step,C?Conj:NonConj> iterator;
      typedef CVIt<T,S==1?Unit:Step,C?Conj:NonConj> const_iterator;
      typedef VIt<T,-S==1?Unit:Step,C?Conj:NonConj> reverse_iterator;
      typedef CVIt<T,-S==1?Unit:Step,C?Conj:NonConj> const_reverse_iterator;
      typedef RefType(T) reference;

      inline reference operator[](int i) const 
      { 
	TMVAssert(i>=0 && i<N);
	return ref(i); 
      }
      inline reference operator()(int i) const 
      {
	TMVAssert(i>=0 && i<N);
	return ref(i); 
      }

      inline iterator begin() const 
      { return iterator(ptr(),S FIRSTLAST ); }
      inline iterator end() const 
      { return begin() + N; }
      inline reverse_iterator rbegin() const 
      { return iterator(ptr()+S*(N-1),-S FIRSTLAST ); }
      inline reverse_iterator rend() const 
      { return rbegin()+N; }

      //
      // Modifying Functions
      //

      inline const SmallVectorView<T,N,S,C,I>& Zero() const 
      { 
	if (S == 1) memset(ptr(),0,N*sizeof(T));
	else for(int i=0,ii=0; i<N; ++i,ii+=S) itsv[ii] = T(0);
	return *this;
      }

      inline const SmallVectorView<T,N,S,C,I>& Clip(RealType(T) thresh) const
      {
	for(int i=0,ii=0; i<N; ++i,ii+=S)
	  if (std::abs(itsv[ii]) < thresh) itsv[ii] = T(0);
	return *this;
      }

      inline const SmallVectorView<T,N,S,C,I>& SetAllTo(T x) const
      {
	if (C) x = CONJ(x);
	for(int i=0,ii=0; i<N; ++i,ii+=S) itsv[ii] = x;
	return *this;
      }

      inline const SmallVectorView<T,N,S,C,I>& AddToAll(T x) const
      {
	if (C) x = CONJ(x);
	for(int i=0,ii=0; i<N; ++i,ii+=S) itsv[ii] += x;
	return *this;
      }

      inline const SmallVectorView<T,N,S,C,I>& ConjugateSelf() const
      {
	if (IsComplex(T())) {
	  RealType(T)* itsvi = reinterpret_cast<RealType(T)*>(itsv)+1;
	  for(int i=0,ii=0; i<N; ++i,ii+=2*S) itsvi[ii] = -itsvi[ii];
	}
	return *this;
      }

      inline const SmallVectorView<T,N,S,C,I>& MakeBasis(int i) const
      { 
	TMVAssert(i>=0 && i < N); 
	Zero(); 
	itsv[i*S] = T(1); 
	return *this; 
      }

      inline const SmallVectorView<T,N,S,C,I>& Swap(int i1, int i2) const
      {
	TMVAssert(i1>=0 && i1<N); 
	TMVAssert(i2>=0 && i2<N); 
	if (i1 != i2) SWAP(itsv[i1*S],itsv[i2*S]); 
	return *this; 
      }

      inline const SmallVectorView<T,N,S,C,I>& Permute(const int* p, 
	  int i1, int i2) const
      {
	TMVAssert(i1>=0 && i1<=i2 && i2<=N);
	for(int i=i1;i<i2;++i) Swap(i,p[i]);
	return *this;
      }

      inline const SmallVectorView<T,N,S,C,I>& Permute(const int* p) const
      { return Permute(p,0,N); }

      inline const SmallVectorView<T,N,S,C,I>& ReversePermute(const int* p, 
	  int i1, int i2) const
      {
	TMVAssert(i1>=0 && i1<=i2 && i2<=N);
	for(int i=i2;i>i1;) { --i; Swap(i,p[i]); }
	return *this;
      }

      inline const SmallVectorView<T,N,S,C,I>& ReversePermute(
	  const int* p) const
      { return ReversePermute(p,0,N); }

      inline const SmallVectorView<T,N,S,C,I>& ReverseSelf() const
      {
	for(int i1=0,i2=N-1;i1<i2;++i1,--i2) Swap(i1,i2);
	return *this;
      }

      inline const SmallVectorView<T,N,S,C,I>& Sort(
	  int* P, ADType ad=ASCEND, COMPType comp=REAL_COMP) const
      { RegView().Sort(P,ad,comp); return *this; }

      //
      // SubVector
      //

      inline VectorView<T,I> SubVector(int i1, int i2) const
      {
	TMVAssert(RegView().OKSubVector(i1,i2,1));
	return RegView().SubVector(i1,i2);
      }

      inline VectorView<T,I> SubVector(int i1, int i2, int istep) const
      {
	TMVAssert(RegView().OKSubVector(i1,i2,istep));
	return RegView().SubVector(i1,i2,istep);
      }

      inline SmallVectorView<T,N,-S,C,I> Reverse() const
      { return SmallVectorView<T,N,-S,C,I>(ptr()+(N-1)*S FIRSTLAST ); }

      inline SmallVectorView<T,N,S,C,I> View() const
      { return *this; }

      inline SmallVectorView<T,N,S,!C,I> Conjugate() const
      { return SmallVectorView<T,N,S,!C,I>(ptr() FIRSTLAST ); }

      inline VectorView<T,I> RegView() const
      { return VectorView<T,I>(ptr(),N,S,C?Conj:NonConj FIRSTLAST ); }

      inline VectorView<RealType(T),I> Real() const
      { return RegView().Real(); }

      inline VectorView<RealType(T),I> Imag() const
      { return RegView().Imag(); }

      inline VectorView<RealType(T),I> Flatten() const
      { return RegView().Flatten(); }


      //
      // I/O
      //

      inline void Read(std::istream& fin) const
      { RegView().Read(fin); }

      // 
      // Iterator Typedefs
      //

      inline const T* cptr() const { return itsv; }
      inline T* ptr() const { return itsv; }

    protected:

      inline reference ref(int i) const
      {
	TMVAssert(i>=0 && i < N);
	// This S>0? thing shouldn't be necessary, but g++-4 -O3 on my
	// MacBook seems to optimize this wrong when S is negative.
	// This should be as easily optimized as the simpler statement of
	// T* vi = itsv + i*S;
	// so there shouldn't be any performance hit from this.
	T* vi = S > 0 ? itsv + i*S : S < 0 ? itsv - i*(-S) : itsv;
	return REF(vi,this->ct());
      }

    private:

      T*const itsv;

#ifdef TMVFLDEBUG
    public:
      const T*const first;
      const T*const last;
#endif

  }; // SmallVectorView

  template <class T, int N, int S, bool C> 
    class SmallVectorView<T,N,S,C,FortranStyle> : 
    public SmallVectorView<T,N,S,C,CStyle>
    {
      public:

	//
	// Constructors 
	//

	inline SmallVectorView(
	    const SmallVectorView<T,N,S,C,FortranStyle>& rhs) : 
	  SmallVectorView<T,N,S,C,CStyle>(rhs) {}

	inline SmallVectorView(const SmallVectorView<T,N,S,C,CStyle>& rhs) : 
	  SmallVectorView<T,N,S,C,CStyle>(rhs) {}

	inline SmallVectorView(const VectorView<T>& rhs) : 
	  SmallVectorView<T,N,S,C,CStyle>(rhs) {}

	inline SmallVectorView(T* inv PARAMFIRSTLAST(T) ) :
	  SmallVectorView<T,N,S,C,CStyle>(inv FIRSTLAST1(_first,_last) ) {}

	virtual inline ~SmallVectorView() {}

	//
	// Op =
	//

	inline const SmallVectorView<T,N,S,C,FortranStyle>& operator=(
	    const SmallVectorView<T,N,S,C,FortranStyle>& v2) const
	{ SmallVectorView<T,N,S,C,CStyle>::operator=(v2); return *this; }

	inline const SmallVectorView<T,N,S,C,FortranStyle>& operator=(
	    const GenSmallVector<T,N,S,C>& v2) const
	{ SmallVectorView<T,N,S,C,CStyle>::operator=(v2); return *this; }

	template <class T2, int S2, bool C2> 
	  inline const SmallVectorView<T,N,S,C,FortranStyle>& operator=(
	      const GenSmallVector<T2,N,S2,C2>& v2) const
	  { SmallVectorView<T,N,S,C,CStyle>::operator=(v2); return *this; }

	inline const SmallVectorView<T,N,S,C,FortranStyle>& operator=(
	    const GenVector<T>& v2) const
	{ SmallVectorView<T,N,S,C,CStyle>::operator=(v2); return *this; }

	template <class T2> 
	  inline const SmallVectorView<T,N,S,C,FortranStyle>& operator=(
	      const GenVector<T2>& v2) const
	  { SmallVectorView<T,N,S,C,CStyle>::operator=(v2); return *this; }

	inline const SmallVectorView<T,N,S,C,FortranStyle>& operator=(
	    const SmallVectorComposite<T,N>& vcomp) const
	{ SmallVectorView<T,N,S,C,CStyle>::operator=(vcomp); return *this; }

	//
	// Access Functions
	//

	inline RefType(T) operator[](int i) const 
	{ 
	  TMVAssert(i>0 && i<=N);
	  return SmallVectorView<T,N,S,C,CStyle>::ref(i-1); 
	}
	inline RefType(T) operator()(int i) const 
	{
	  TMVAssert(i>0 && i<=N);
	  return SmallVectorView<T,N,S,C,CStyle>::ref(i-1); 
	}

	//
	// Modifying Functions
	//

	inline const SmallVectorView<T,N,S,C,FortranStyle>& Zero() const 
	{ SmallVectorView<T,N,S,C,CStyle>::Zero(); return *this; }

	inline const SmallVectorView<T,N,S,C,FortranStyle>& Clip(
	    RealType(T) thresh) const
	{ SmallVectorView<T,N,S,C,CStyle>::Clip(thresh); return *this; }

	inline const SmallVectorView<T,N,S,C,FortranStyle>& SetAllTo(T x) const
	{ SmallVectorView<T,N,S,C,CStyle>::SetAllTo(x); return *this; }

	inline const SmallVectorView<T,N,S,C,FortranStyle>& AddToAll(T x) const
	{ SmallVectorView<T,N,S,C,CStyle>::AddToAll(x); return *this; }

	inline const SmallVectorView<T,N,S,C,FortranStyle>& ConjugateSelf() const
	{ SmallVectorView<T,N,S,C,CStyle>::ConjugateSelf(); return *this; }

	inline const SmallVectorView<T,N,S,C,FortranStyle>& MakeBasis(
	    int i) const
	{ 
	  TMVAssert(i>0 && i<=N); 
	  SmallVectorView<T,N,S,C,CStyle>::MakeBasis(i-1);
	  return *this; 
	}

	inline const SmallVectorView<T,N,S,C,FortranStyle>& Swap(
	    int i1, int i2) const
	{ 
	  TMVAssert(i1>0 && i1<=N);
	  TMVAssert(i2>0 && i2<=N);
	  if (i1!=i2) SmallVectorView<T,N,S,C,CStyle>::Swap(i1-1,i2-1);
	  return *this; 
	}

	inline const SmallVectorView<T,N,S,C,FortranStyle>& Permute(
	    const int* p, int i1, int i2) const
	{
	  TMVAssert(i1>0);
	  SmallVectorView<T,N,S,C,CStyle>::Permute(p,i1-1,i2);
	  return *this;
	}

	inline const SmallVectorView<T,N,S,C,FortranStyle>& Permute(
	    const int* p) const
	{ SmallVectorView<T,N,S,C,CStyle>::Permute(p); return *this; }

	inline const SmallVectorView<T,N,S,C,FortranStyle>& ReversePermute(
	    const int* p, int i1, int i2) const
	{
	  TMVAssert(i1>0);
	  SmallVectorView<T,N,S,C,CStyle>::ReversePermute(p,i1-1,i2);
	  return *this;
	}

	inline const SmallVectorView<T,N,S,C,FortranStyle>& ReversePermute(
	    const int* p) const
	{ SmallVectorView<T,N,S,C,CStyle>::ReversePermute(p); return *this; }

	inline const SmallVectorView<T,N,S,C,FortranStyle>& ReverseSelf() const
	{ SmallVectorView<T,N,S,C,CStyle>::ReverseSelf(); return *this; }

	inline const SmallVectorView<T,N,S,C,FortranStyle>& Sort(
	    int* P, ADType ad=ASCEND, COMPType comp=REAL_COMP) const
	{ SmallVectorView<T,N,S,C,CStyle>::Sort(P,ad,comp); return *this; }

	//
	// SubVector
	//

	inline VectorView<T,FortranStyle> SubVector(int i1, int i2) const
	{
	  TMVAssert(RegView().OKSubVector(i1,i2,1));
	  return RegView().SubVector(i1,i2);
	}

	inline VectorView<T,FortranStyle> SubVector(
	    int i1, int i2, int istep) const
	{
	  TMVAssert(RegView().OKSubVector(i1,i2,istep));
	  return RegView().SubVector(i1,i2,istep);
	}

	inline SmallVectorView<T,N,-S,C,FortranStyle> Reverse() const
	{ return SmallVectorView<T,N,S,C,CStyle>::Reverse(); }

	inline SmallVectorView<T,N,S,C,FortranStyle> View() const
	{ return SmallVectorView<T,N,S,C,CStyle>::View(); }

	inline SmallVectorView<T,N,S,!C,FortranStyle> Conjugate() const
	{ return SmallVectorView<T,N,S,C,CStyle>::Conjugate(); }

	inline VectorView<T,FortranStyle> RegView() const
	{ return SmallVectorView<T,N,S,C,CStyle>::RegView(); }

	inline VectorView<RealType(T),FortranStyle> Real() const
	{ return SmallVectorView<T,N,S,C,CStyle>::Real(); }

	inline VectorView<RealType(T),FortranStyle> Imag() const
	{ return SmallVectorView<T,N,S,C,CStyle>::Imag(); }

	inline VectorView<RealType(T),FortranStyle> Flatten() const
	{ return SmallVectorView<T,N,S,C,CStyle>::Flatten(); }

    }; // FortranStyle SmallVectorView

#ifdef XTEST
#ifdef TMVDEBUG
#define XTEST_DEBUG
#endif
#endif

  template <class T, int N, IndexStyle I> class SmallVector : 
    public GenSmallVector<T,N,1,false>
  {
    public:

      //
      // Constructors
      //

#ifdef TMVFLDEBUG
#define NEW_SIZE \
      : DEFFIRSTLAST(itsv,itsv+N)
#else
#define NEW_SIZE 
#endif

      inline SmallVector() NEW_SIZE
      {
	TMVAssert(N > 0);
#ifdef TMVDEBUG
	SetAllTo(T(888));
#endif
      }

      explicit inline SmallVector(T val) NEW_SIZE
      {
	TMVAssert(N > 0);
	SetAllTo(val); 
      }

      explicit inline SmallVector(const T* vv) NEW_SIZE
      {
#ifdef XTEST_DEBUG
	SetAllTo(T(888));
#endif
	TMVAssert(N > 0);
	memmove(itsv,vv,N*sizeof(T)); 
      }

      explicit inline SmallVector(const std::vector<T>& vv) NEW_SIZE
      {
#ifdef XTEST_DEBUG
	SetAllTo(T(888));
#endif
	TMVAssert(N > 0);
	for(int i=0;i<N;i++) itsv[i] = vv[i];
      }

      inline SmallVector(const SmallVector<T,N,I>& rhs) NEW_SIZE
      {
#ifdef XTEST_DEBUG
	SetAllTo(T(888));
#endif
	TMVAssert(N > 0);
	memmove(itsv,rhs.cptr(),N*sizeof(T)); 
      }

      inline SmallVector(const GenSmallVector<T,N,1,false>& rhs) NEW_SIZE
      {
#ifdef XTEST_DEBUG
	SetAllTo(T(888));
#endif
	TMVAssert(N > 0);
	memmove(itsv,rhs.cptr(),N*sizeof(T)); 
      }

      template <IndexStyle I2> inline SmallVector(const Vector<T,I2>& rhs) 
	NEW_SIZE
      {
#ifdef XTEST_DEBUG
	SetAllTo(T(888));
#endif
	TMVAssert(N > 0);
	TMVAssert(rhs.size() == N);
	memmove(itsv,rhs.cptr(),N*sizeof(T)); 
      }

      template <class T2, int S2, bool C2> inline SmallVector(
	  const GenSmallVector<T2,N,S2,C2>& rhs) NEW_SIZE
      {
#ifdef XTEST_DEBUG
	SetAllTo(T(888));
#endif
	TMVAssert(N > 0);
	TMVAssert(IsReal(T2()) || IsComplex(T()));
	Copy(rhs,View()); 
      }

      inline SmallVector(const GenVector<RealType(T)>& v2) NEW_SIZE
      {
#ifdef XTEST_DEBUG
	SetAllTo(T(888));
#endif
	TMVAssert(N > 0);
	TMVAssert(v2.size() == N); v2.AssignToV(RegView()); 
      }

      inline SmallVector(const GenVector<ComplexType(T)>& v2) NEW_SIZE
      {
#ifdef XTEST_DEBUG
	SetAllTo(T(888));
#endif
	TMVAssert(N > 0);
	TMVAssert(IsComplex(T()));
	TMVAssert(v2.size() == N); 
	v2.AssignToV(RegView()); 
      }

      template <class T2> inline SmallVector(const GenVector<T2>& rhs) 
	NEW_SIZE
      {
#ifdef XTEST_DEBUG
	SetAllTo(T(888));
#endif
	TMVAssert(N > 0);
	TMVAssert(IsReal(T2()) || IsComplex(T()));
	TMVAssert(rhs.size() == N);
	Copy(rhs,RegView()); 
      }

      inline SmallVector(const AssignableToSmallVector<RealType(T),N>& v2) 
	NEW_SIZE
      {
#ifdef XTEST_DEBUG
	SetAllTo(T(888));
#endif
	TMVAssert(N > 0);
	v2.AssignTov(View()); 
      }

      inline SmallVector(const AssignableToSmallVector<ComplexType(T),N>& v2) 
	NEW_SIZE
      {
#ifdef XTEST_DEBUG
	SetAllTo(T(888));
#endif
	TMVAssert(N > 0);
	TMVAssert(IsComplex(T())); v2.AssignTov(View()); 
      }

#undef NEW_SIZE

      virtual inline ~SmallVector() {}


      //
      // Op =
      //

      inline SmallVector<T,N>& operator=(SmallVector<T,N>& v2)
      { 
	if (&v2 != this) memmove(itsv,v2.cptr(),N*sizeof(T)); 
	return *this; 
      }

      template <class T2, int S, bool C> inline SmallVector<T,N,I>& operator=(
	  const GenSmallVector<T2,N,S,C>& v2) 
      { 
	TMVAssert(IsReal(T2()) || IsComplex(T()));
	Copy(v2,View()); 
	return *this; 
      }

      inline SmallVector<T,N,I>& operator=(const GenVector<RealType(T)>& v2) 
      {
	TMVAssert(v2.size() == N);
	v2.AssignToV(RegView());
	return *this; 
      }

      inline SmallVector<T,N,I>& operator=(const GenVector<ComplexType(T)>& v2) 
      {
	TMVAssert(v2.size() == N);
	TMVAssert(IsComplex(T()));
	v2.AssignToV(RegView());
	return *this; 
      }

      template <class T2> inline SmallVector<T,N,I>& operator=(
	  const GenVector<T2>& v2) 
      {
	TMVAssert(IsReal(T2()) || IsComplex(T()));
	TMVAssert(v2.size() == N);
	Copy(v2,RegView());
	return *this; 
      }

      inline SmallVector<T,N,I>& operator=(
	  const AssignableToSmallVector<RealType(T),N>& v2)
      { v2.AssignTov(View()); return *this; }

      inline SmallVector<T,N,I>& operator=(
	  const AssignableToSmallVector<ComplexType(T),N>& v2)
      {
	TMVAssert(IsComplex(T()));
	v2.AssignTov(View()); 
	return *this; 
      }

      //
      // Access Functions
      //

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

      inline T operator[](int i) const 
      { 
	if (I == CStyle) { TMVAssert(i>=0 && i<N); return cref(i); }
	else { TMVAssert(i>0 && i<=N); return cref(i-1); }
      }
      inline T operator()(int i) const 
      { 
	if (I == CStyle) { TMVAssert(i>=0 && i<N); return cref(i); }
	else { TMVAssert(i>0 && i<=N); return cref(i-1); }
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
	if (I == CStyle) { TMVAssert(i>=0 && i<N); return ref(i); }
	else { TMVAssert(i>0 && i<=N); return ref(i-1); }
      }
      inline T& operator()(int i) 
      { 
	if (I == CStyle) { TMVAssert(i>=0 && i<N); return ref(i); }
	else { TMVAssert(i>0 && i<=N); return ref(i-1); }
      }

      //
      // Modifying Functions
      //

      inline SmallVector<T,N,I>& Zero() 
      { 
	memset(itsv,0,N*sizeof(T)); 
	return *this;
      }

      inline SmallVector<T,N,I>& Clip(RealType(T) thresh)
      { View().Clip(thresh); return *this; }

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
	}
	else {
	  TMVAssert(i>0 && i<=N);
	  Zero(); itsv[i-1] = T(1);
	}
	return *this; 
      }

      inline SmallVector<T,N,I>& Swap(int i1, int i2)
      {
	if (I == CStyle) {
	  TMVAssert(i1>=0 && i1<N);
	  TMVAssert(i2>=0 && i2<N);
	  if (i1 != i2) SWAP(itsv[i1],itsv[i2]);
	} else {
	  TMVAssert(i1>0 && i1<=N);
	  TMVAssert(i2>0 && i2<=N);
	  if (i1 != i2) SWAP(itsv[i1-1],itsv[i2-1]);
	}
	return *this;
      }

      inline SmallVector<T,N,I>& Permute(const int* p, int i1, int i2)
      {
	if (I == CStyle) TMVAssert(i1>=0 && i1<=i2 && i2<=N);
	else TMVAssert(i1>0 && i1<=i2 && i2<=N);
	View().Permute(p,i1,i2); return *this; 
      }
      inline SmallVector<T,N,I>& Permute(const int* p) 
      { return Permute(p,0,N); }
      inline SmallVector<T,N,I>& ReversePermute(
	  const int* p, int i1, int i2)
      {
	if (I == CStyle) TMVAssert(i1>=0 && i1<=i2 && i2<=N);
	else TMVAssert(i1>0 && i1<=i2 && i2<=N);
	View().ReversePermute(p,i1,i2); return *this; 
      }
      inline SmallVector<T,N,I>& ReversePermute(const int* p) 
      { return ReversePermute(p,0,N); }

      inline SmallVector<T,N,I>& ReverseSelf()
      { View().ReverseSelf(); return *this; }
      inline SmallVector<T,N,I>& Sort(int* P, ADType ad=ASCEND, 
	  COMPType comp=REAL_COMP) 
      { View().Sort(P,ad,comp); return *this; }

      //
      // SubVector
      //

      inline ConstVectorView<T> SubVector(int i1, int i2) const
      {
	TMVAssert(RegView().OKSubVector(i1,i2,1));
	return View().SubVector(i1,i2);
      }

      inline VectorView<T> SubVector(int i1, int i2)
      {
	TMVAssert(RegView().OKSubVector(i1,i2,1));
	return View().SubVector(i1,i2);
      }

      inline ConstVectorView<T> SubVector(int i1, int i2, int istep) const
      {
	TMVAssert(RegView().OKSubVector(i1,i2,istep));
	return View().SubVector(i1,i2,istep);
      }

      inline VectorView<T> SubVector(int i1, int i2, int istep)
      {
	TMVAssert(RegView().OKSubVector(i1,i2,istep));
	return View().SubVector(i1,i2,istep);
      }

      inline ConstSmallVectorView<T,N,-1,false,I> Reverse() const
      { return ConstSmallVectorView<T,N,-1,false,I>(itsv+N-1); }

      inline SmallVectorView<T,N,-1,false,I> Reverse()
      { return SmallVectorView<T,N,-1,false,I>(itsv+N-1 FIRSTLAST ); }

      inline ConstSmallVectorView<T,N,1,false,I> View() const
      { return ConstSmallVectorView<T,N,1,false,I>(itsv); }

      inline SmallVectorView<T,N,1,false,I> View()
      { return SmallVectorView<T,N,1,false,I>(itsv FIRSTLAST ); }

      inline ConstVectorView<T,I> RegView() const
      { return ConstVectorView<T,I>(itsv,N,1,NonConj); }

      inline VectorView<T,I> RegView()
      { return VectorView<T,I>(itsv,N,1,NonConj FIRSTLAST ); }

      inline ConstSmallVectorView<T,N,1,true,I> Conjugate() const
      { return ConstSmallVectorView<T,N,1,true,I>(itsv); }

      inline SmallVectorView<T,N,1,true,I> Conjugate()
      { return SmallVectorView<T,N,1,true,I>(itsv FIRSTLAST ); }

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

      inline size_t size() const { return N; }
      inline const T* cptr() const { return itsv; }
      inline T* ptr() { return itsv; }
      inline int step() const { return 1; }
      inline ConjType ct() const { return NonConj; }
      inline bool isconj() const { return false; }

    private:

      T itsv[N];

#ifdef TMVFLDEBUG
    public:
      const T*const first;
      const T*const last;
    private:
#endif

      inline T cref(int i) const
      { 
	TMVAssert(i>=0 && i < N); 
	return itsv[i];
      }

      inline T& ref(int i)
      { 
	TMVAssert(i>=0 && i < N); 
#ifdef TMVFLDEBUG
	TMVAssert(itsv+i >= first);
	TMVAssert(itsv+i < last);
#endif
	return itsv[i];
      }

  }; // SmallVector


  //
  // Copy SmallVectors
  //

  template <class T1, class T2, int N, int S1, int S2, bool C1>
    inline void Copy(
	const GenSmallVector<T1,N,S1,C1>& v1, 
	const SmallVectorView<T2,N,S2,false>& v2)
  {
    TMVAssert(IsComplex(T2()) || IsReal(T1()));
    for(int i=0;i<N;++i) v2[i] = v1[i];
  }
  template <class T1, class T2, int N, int S1, bool C1>
    inline void Copy(
	const GenSmallVector<T1,N,S1,C1>& v1, T2* v2p, const int S2)
    {
      TMVAssert(IsComplex(T2()) || IsReal(T1()));
      for(int i=0,i2=0;i<N;++i,i2+=S2) v2p[i2] = v1[i];
    }


  //
  // Swap SmallVectors
  //

  template <class T, int N, int S1, int S2, bool C1, bool C2> 
    inline void Swap(
      const SmallVectorView<T,N,S1,C1>& v1,
      const SmallVectorView<T,N,S2,C2>& v2)
  {
    TMVAssert(S1!=0);
    TMVAssert(S2!=0);
    for(int i=0;i<N;++i) SWAP(v1[i],v2[i]);
  }
  template <class T, int N, int S1, bool C1, IndexStyle I2> 
    inline void Swap(
      const SmallVectorView<T,N,S1,C1>& v1, SmallVector<T,N,I2>& v2)
  { Swap(v1,v2.View()); }
  template <class T, int N, IndexStyle I1, int S2, bool C2> 
    inline void Swap(
      SmallVector<T,N,I1>& v1, const SmallVectorView<T,N,S2,C2>& v2) 
  { Swap(v1.View(),v2); }
  template <class T, int N, IndexStyle I1, IndexStyle I2> 
    inline void Swap(
      SmallVector<T,N,I1>& v1, SmallVector<T,N,I2>& v2)
  { Swap(v1.View(),v2.View()); }

  //
  // Functions of Vectors
  //

  template <class T, int N, int S, bool C> inline RealType(T) Norm(
      const GenSmallVector<T,N,S,C>& v)
  { return v.Norm(); }

  template <class T, int N, int S, bool C> inline RealType(T) Norm1(
      const GenSmallVector<T,N,S,C>& v)
  { return v.Norm1(); }

  template <class T, int N, int S, bool C> inline RealType(T) NormSq(
      const GenSmallVector<T,N,S,C>& v)
  { return v.NormSq(); }

  template <class T, int N, int S, bool C> inline RealType(T) Norm2(
      const GenSmallVector<T,N,S,C>& v)
  { return v.Norm2(); }

  template <class T, int N, int S, bool C> inline RealType(T) NormInf(
      const GenSmallVector<T,N,S,C>& v)
  { return v.NormInf(); }

  template <class T, int N, int S, bool C> inline T SumElements(
      const GenSmallVector<T,N,S,C>& v)
  { return v.SumElements(); }

  template <class T, int N, int S, bool C> inline RealType(T) SumAbsElements(
      const GenSmallVector<T,N,S,C>& v)
  { return v.SumAbsElements(); }

  template <class T, int N, int S, bool C> inline T MinElement(
      const GenSmallVector<T,N,S,C>& v)
  { return v.MinElement(); }

  template <class T, int N, int S, bool C> inline T MaxElement(
      const GenSmallVector<T,N,S,C>& v)
  { return v.MaxElement(); }

  template <class T, int N, int S, bool C> inline RealType(T) MinAbsElement(
      const GenSmallVector<T,N,S,C>& v)
  { return v.MinAbsElement(); }

  template <class T, int N, int S, bool C> inline RealType(T) MaxAbsElement(
      const GenSmallVector<T,N,S,C>& v)
  { return v.MaxAbsElement(); }

  template <class T, int N, int S, bool C> 
    inline ConstSmallVectorView<T,N,S,!C> Conjugate(
	const GenSmallVector<T,N,S,C>& v)
    { return v.Conjugate(); }

  template <class T, int N, int S, bool C> 
    inline SmallVectorView<T,N,S,!C> Conjugate(
	const SmallVectorView<T,N,S,C>& v)
    { return v.Conjugate(); }

  template <class T, int N, IndexStyle I> 
    inline SmallVectorView<T,N,1,true> Conjugate(SmallVector<T,N,I>& v)
    { return v.Conjugate(); }


  //
  // Vector ==, != Vector
  //

  template <class T1, class T2, int N, int S1, int S2, bool C1, bool C2> 
    inline bool operator==(
      const GenSmallVector<T1,N,S1,C1>& v1,
      const GenSmallVector<T2,N,S2,C2>& v2)
    {
      for(int i=0;i<N;++i) if (v1[i] != v2[i]) return false;
      return true;
    }

  template <class T1, class T2, int N, int S1, int S2, bool C1, bool C2> 
    inline bool operator!=(
      const GenSmallVector<T1,N,S1,C1>& v1, 
      const GenSmallVector<T2,N,S2,C2>& v2)
  { return !(v1 == v2); }

  //
  // I/O
  //

  template <class T, int N, int S, bool C> inline std::ostream& operator<<(
      std::ostream& os, const GenSmallVector<T,N,S,C>& v)
  { v.Write(os); return os;}

  template <class T, int N, int S, bool C> inline std::istream& operator>>(
      std::istream& is, const SmallVectorView<T,N,S,C>& v)
  { is >> v.RegView(); return is; }

  template <class T, int N, IndexStyle I> inline std::istream& operator>>(
      std::istream& is, SmallVector<T,N,I>& v)
  { return is >> v.RegView(); }

#ifdef XDEBUG
  template <class T, int N, IndexStyle I> inline std::string Type(
      const SmallVector<T,N,I>& )
  { 
    std::ostringstream s;
    s << "SmallVector<"<<Type(T())<<","<<N<<","<<Text(I)<<">";
    return s.str();
  }
  template <class T, int N, int S, bool C> inline std::string Type(
      const GenSmallVector<T,N,S,C>& )
  { 
    std::ostringstream s;
    s << "GenSmallVector<"<<Type(T())<<","<<N<<","<<S<<",";
    s << (C ? "true" : "false") << ">";
    return s.str();
  }
  template <class T, int N, int S, bool C, IndexStyle I> 
    inline std::string Type(
	const ConstSmallVectorView<T,N,S,C,I>& )
    { 
      std::ostringstream s;
      s << "ConstSmallVectorView<"<<Type(T())<<","<<N<<","<<S<<",";
      s << (C ? "true" : "false") <<","<<Text(I)<< ">";
      return s.str();
    }
  template <class T, int N, int S, bool C, IndexStyle I> 
    inline std::string Type(
	const SmallVectorView<T,N,S,C,I>& )
    { 
      std::ostringstream s;
      s << "SmallVectorView<"<<Type(T())<<","<<N<<","<<S<<",";
      s << (C ? "true" : "false") <<","<<Text(I)<< ">";
      return s.str();
    }
#endif

} // namespace tmv

#endif
