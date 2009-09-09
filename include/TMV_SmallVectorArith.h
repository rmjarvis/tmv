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


#ifndef TMV_SmallVectorArith_H
#define TMV_SmallVectorArith_H

#include "TMV_SmallVectorArithFunc.h"
#include "TMV_VectorArithFunc.h"

#define CT std::complex<T>
#define CCT ConjRef<std::complex<T> >
#define VCT VarConjRef<std::complex<T> >

namespace tmv {

  template <class T, size_t N> class SmallVectorComposite : 
    public GenSmallVector<T,N,1,false>
  {
    public:

      inline SmallVectorComposite() : itsv(0) {}
      inline SmallVectorComposite(const SmallVectorComposite<T,N>&) :
	itsv(0) {}
      virtual inline ~SmallVectorComposite() {}

      inline const T* cptr() const
      {
	if (!itsv.get()) {
	  itsv.reset(new T[N]);
	  AssignTov(SmallVectorView<T,N,1,false>(itsv.get()
	      FIRSTLAST1(itsv.get(),itsv.get()+N) ) );
	}
	return itsv.get();
      }

    private:
      mutable auto_array<T> itsv;
  };

  // These are what we want to do no matter what type Tx is:
  template <class T, size_t N, IndexStyle I, class Tx> 
    inline SmallVector<T,N,I>& operator+=(
	SmallVector<T,N,I>& v, const Tx& x)
    { v.View() += x; return v; }

  template <class T, size_t N, IndexStyle I, class Tx> 
    inline SmallVector<T,N,I>& operator-=(
	SmallVector<T,N,I>& v, const Tx& x)
    { v.View() -= x; return v; }

  template <class T, size_t N, IndexStyle I, class Tx> 
    inline SmallVector<T,N,I>& operator*=(
	SmallVector<T,N,I>& v, const Tx& x)
    { v.View() *= x; return v; }

  template <class T, size_t N, IndexStyle I, class Tx> 
    inline SmallVector<T,N,I>& operator/=(
	SmallVector<T,N,I>& v, const Tx& x)
    { v.View() /= x; return v; }

  template <class T, size_t N, IndexStyle I, class Tx> 
    inline SmallVector<T,N,I>& operator%=(
	SmallVector<T,N,I>& v, const Tx& x)
    { v.View() %= x; return v; }


  //
  // Vector * / Scalar
  //

  template <class T, class T1, size_t N, int S1, bool C1> class ProdXv : 
    public SmallVectorComposite<T,N> 
  {
    public:
      inline ProdXv(T _x, const GenSmallVector<T1,N,S1,C1>& _v) : 
	x(_x), v(_v) {}
      inline size_t size() const { return v.size(); }
      inline T GetX() const { return x; }
      inline const GenSmallVector<T1,N,S1,C1>& GetV() const { return v; }
      inline void AssignTov(
	  const SmallVectorView<RealType(T),N,1,false>& v0) const
      { TMVAssert(IsReal(T())); MultXV(x,v,v0); }
      inline void AssignTov(
	  const SmallVectorView<ComplexType(T),N,1,false>& v0) const
      { MultXV(x,v,v0); }
      inline void AssignTov(
	  const SmallVectorView<ComplexType(T),N,1,true>& v0) const
      { MultXV(x,v,v0); }
      inline void DoAssignTov(RealType(T)* vp, const int S) const
      { TMVAssert(IsReal(T())); MultXV(x,v,vp,S); }
      inline void DoAssignTov(ComplexType(T)* vp,
	  const int S, const bool C) const
      {
	if (C) MultXV(CONJ(x),v.Conjugate(),vp,S);
	else MultXV(x,v,vp,S);
      }
    private:
      const T x;
      const GenSmallVector<T1,N,S1,C1>& v;
  };

  template <class T, size_t N, int S, bool C> 
    inline const SmallVectorView<T,N,S,C>& operator*=(
	const SmallVectorView<T,N,S,C>& v1, T x2)
    { MultXV(x2,v1); return v1; }

  template <class T, size_t N, int S, bool C> 
    inline const SmallVectorView<CT,N,S,C>& operator*=(
	const SmallVectorView<CT,N,S,C>& v1, T x2)
    { MultXV(x2,v1); return v1; }

  template <class T, size_t N, int S, bool C> 
    inline const SmallVectorView<CT,N,S,C>& operator*=(
	const SmallVectorView<CT,N,S,C>& v1, CCT x2)
    { MultXV(CT(x2),v1); return v1; }

  template <class T, size_t N, int S, bool C> 
    inline const SmallVectorView<CT,N,S,C>& operator*=(
	const SmallVectorView<CT,N,S,C>& v1, VCT x2)
    { MultXV(CT(x2),v1); return v1; }

  template <class T, size_t N, int S, bool C> 
    inline const SmallVectorView<T,N,S,C>& operator/=(
	const SmallVectorView<T,N,S,C>& v1, T x2)
    { MultXV(T(1)/x2,v1); return v1; }

  template <class T, size_t N, int S, bool C> 
    inline const SmallVectorView<CT,N,S,C>& operator/=(
	const SmallVectorView<CT,N,S,C>& v1, T x2)
    { MultXV(T(1)/x2,v1); return v1; }

  template <class T, size_t N, int S, bool C> 
    inline const SmallVectorView<CT,N,S,C>& operator/=(
	const SmallVectorView<CT,N,S,C>& v1, CCT x2)
    { MultXV(T(1)/CT(x2),v1); return v1; }

  template <class T, size_t N, int S, bool C> 
    inline const SmallVectorView<CT,N,S,C>& operator/=(
	const SmallVectorView<CT,N,S,C>& v1, VCT x2)
    { MultXV(T(1)/CT(x2),v1); return v1; }

#define GENVECTOR GenSmallVector
#define PRODXV ProdXv
#define X ,N,S1,C1
#define Y ,size_t N,int S1,bool C1
#include "TMV_AuxProdXV.h"
#undef GENVECTOR
#undef PRODXV

  //
  // Vector + Vector
  //

  template <class T, class T1, class T2, size_t N, int S1, bool C1, int S2, bool C2> class Sumvv_1_1 : 
    public SmallVectorComposite<T,N> 
  {
    public:
      inline Sumvv_1_1(
	  T DEBUGPARAM(_x1), const GenSmallVector<T1,N,S1,C1>& _v1, 
	  T DEBUGPARAM(_x2), const GenSmallVector<T2,N,S2,C2>& _v2) :
	v1(_v1),v2(_v2)
      { TMVAssert(_x1 == T(1) && _x2 == T(1)); }
      inline size_t size() const { return v1.size(); }
      inline const GenSmallVector<T1,N,S1,C1>& GetV1() const { return v1; }
      inline const GenSmallVector<T2,N,S2,C2>& GetV2() const { return v2; }
      inline void AssignTov(
	  const SmallVectorView<RealType(T),N,1,false>& v0) const
      { TMVAssert(IsReal(T())); AddVV_1_1(v1,v2,v0); }
      inline void AssignTov(
	  const SmallVectorView<ComplexType(T),N,1,false>& v0) const
      { AddVV_1_1(v1,v2,v0); }
      inline void AssignTov(
	  const SmallVectorView<ComplexType(T),N,1,true>& v0) const
      { AddVV_1_1(v1,v2,v0); }
      inline void DoAssignTov(RealType(T)* vp, const int S) const
      { TMVAssert(IsReal(T())); AddVV_1_1(v1,v2,vp,S); }
      inline void DoAssignTov(ComplexType(T)* vp, 
	  const int S, const bool C) const
      {
	if (C) AddVV_1_1(v1.Conjugate(),v2.Conjugate(),vp,S); 
	else AddVV_1_1(v1,v2,vp,S); 
      }
    private:
      const GenSmallVector<T1,N,S1,C1>& v1;
      const GenSmallVector<T2,N,S2,C2>& v2;
  };

  template <class T, class T1, class T2, size_t N, int S1, bool C1, int S2, bool C2> class Sumvv_1_m1 : 
    public SmallVectorComposite<T,N> 
  {
    public:
      inline Sumvv_1_m1(
	  T DEBUGPARAM(_x1), const GenSmallVector<T1,N,S1,C1>& _v1, 
	  T DEBUGPARAM(_x2), const GenSmallVector<T2,N,S2,C2>& _v2) :
	v1(_v1),v2(_v2)
      { TMVAssert(_x1 == T(1) && _x2 == T(-1)); }
      inline size_t size() const { return v1.size(); }
      inline const GenSmallVector<T1,N,S1,C1>& GetV1() const { return v1; }
      inline const GenSmallVector<T2,N,S2,C2>& GetV2() const { return v2; }
      inline void AssignTov(
	  const SmallVectorView<RealType(T),N,1,false>& v0) const
      { TMVAssert(IsReal(T())); AddVV_1_m1(v1,v2,v0); }
      inline void AssignTov(
	  const SmallVectorView<ComplexType(T),N,1,false>& v0) const
      { AddVV_1_m1(v1,v2,v0); }
      inline void AssignTov(
	  const SmallVectorView<ComplexType(T),N,1,true>& v0) const
      { AddVV_1_m1(v1,v2,v0); }
      inline void DoAssignTov(RealType(T)* vp, const int S) const
      { TMVAssert(IsReal(T())); AddVV_1_m1(v1,v2,vp,S); }
      inline void DoAssignTov(ComplexType(T)* vp, 
	  const int S, const bool C) const
      {
	if (C) AddVV_1_m1(v1.Conjugate(),v2.Conjugate(),vp,S); 
	else AddVV_1_m1(v1,v2,vp,S); 
      }
    private:
      const GenSmallVector<T1,N,S1,C1>& v1;
      const GenSmallVector<T2,N,S2,C2>& v2;
  };

  template <class T, class T1, class T2, size_t N, int S1, bool C1, int S2, bool C2> class Sumvv_1_x : 
    public SmallVectorComposite<T,N> 
  {
    public:
      inline Sumvv_1_x(
	  T DEBUGPARAM(_x1), const GenSmallVector<T1,N,S1,C1>& _v1, 
	  T _x2, const GenSmallVector<T2,N,S2,C2>& _v2) :
	v1(_v1),x2(_x2), v2(_v2)
      { TMVAssert(_x1==T(1)); }
      inline size_t size() const { return v1.size(); }
      inline const GenSmallVector<T1,N,S1,C1>& GetV1() const { return v1; }
      inline T GetX2() const { return x2; }
      inline const GenSmallVector<T2,N,S2,C2>& GetV2() const { return v2; }
      inline void AssignTov(
	  const SmallVectorView<RealType(T),N,1,false>& v0) const
      { TMVAssert(IsReal(T())); AddVV_1_x(v1,x2,v2,v0); }
      inline void AssignTov(
	  const SmallVectorView<ComplexType(T),N,1,false>& v0) const
      { AddVV_1_x(v1,x2,v2,v0); }
      inline void AssignTov(
	  const SmallVectorView<ComplexType(T),N,1,true>& v0) const
      { AddVV_1_x(v1,x2,v2,v0); }
      inline void DoAssignTov(RealType(T)* vp, const int S) const
      { TMVAssert(IsReal(T())); AddVV_1_x(v1,x2,v2,vp,S); }
      inline void DoAssignTov(ComplexType(T)* vp, 
	  const int S, const bool C) const
      {
	if (C) AddVV_1_x(v1.Conjugate(),CONJ(x2),v2.Conjugate(),vp,S); 
	else AddVV_1_x(v1,x2,v2,vp,S); 
      }
    private:
      const GenSmallVector<T1,N,S1,C1>& v1;
      const T x2;
      const GenSmallVector<T2,N,S2,C2>& v2;
  };

  template <class T, class T1, class T2, size_t N, int S1, bool C1, int S2, bool C2> class Sumvv_x_1 : 
    public SmallVectorComposite<T,N> 
  {
    public:
      inline Sumvv_x_1(
	  T _x1, const GenSmallVector<T1,N,S1,C1>& _v1, 
	  T DEBUGPARAM(_x2), const GenSmallVector<T2,N,S2,C2>& _v2) :
	x1(_x1),v1(_v1),v2(_v2)
      { TMVAssert(_x2 == T(1)); }
      inline size_t size() const { return v1.size(); }
      inline T GetX1() const { return x1; }
      inline const GenSmallVector<T1,N,S1,C1>& GetV1() const { return v1; }
      inline const GenSmallVector<T2,N,S2,C2>& GetV2() const { return v2; }
      inline void AssignTov(
	  const SmallVectorView<RealType(T),N,1,false>& v0) const
      { TMVAssert(IsReal(T())); AddVV_x_1(x1,v1,v2,v0); }
      inline void AssignTov(
	  const SmallVectorView<ComplexType(T),N,1,false>& v0) const
      { AddVV_x_1(x1,v1,v2,v0); }
      inline void AssignTov(
	  const SmallVectorView<ComplexType(T),N,1,true>& v0) const
      { AddVV_x_1(x1,v1,v2,v0); }
      inline void DoAssignTov(RealType(T)* vp, const int S) const
      { TMVAssert(IsReal(T())); AddVV_x_1(x1,v1,v2,vp,S); }
      inline void DoAssignTov(ComplexType(T)* vp, 
	  const int S, const bool C) const
      {
	if (C) AddVV_x_1(CONJ(x1),v1.Conjugate(),v2.Conjugate(),vp,S); 
	else AddVV_x_1(x1,v1,v2,vp,S); 
      }
    private:
      const T x1;
      const GenSmallVector<T1,N,S1,C1>& v1;
      const GenSmallVector<T2,N,S2,C2>& v2;
  };

  template <class T, class T1, class T2, size_t N, int S1, bool C1, int S2, bool C2> class Sumvv_x_m1 : 
    public SmallVectorComposite<T,N> 
  {
    public:
      inline Sumvv_x_m1(
	  T _x1, const GenSmallVector<T1,N,S1,C1>& _v1, 
	  T DEBUGPARAM(_x2), const GenSmallVector<T2,N,S2,C2>& _v2) :
	x1(_x1),v1(_v1),v2(_v2)
      { TMVAssert(_x2 == T(-1)); }
      inline size_t size() const { return v1.size(); }
      inline T GetX1() const { return x1; }
      inline const GenSmallVector<T1,N,S1,C1>& GetV1() const { return v1; }
      inline const GenSmallVector<T2,N,S2,C2>& GetV2() const { return v2; }
      inline void AssignTov(
	  const SmallVectorView<RealType(T),N,1,false>& v0) const
      { TMVAssert(IsReal(T())); AddVV_x_m1(x1,v1,v2,v0); }
      inline void AssignTov(
	  const SmallVectorView<ComplexType(T),N,1,false>& v0) const
      { AddVV_x_m1(x1,v1,v2,v0); }
      inline void AssignTov(
	  const SmallVectorView<ComplexType(T),N,1,true>& v0) const
      { AddVV_x_m1(x1,v1,v2,v0); }
      inline void DoAssignTov(RealType(T)* vp, const int S) const
      { TMVAssert(IsReal(T())); AddVV_x_m1(x1,v1,v2,vp,S); }
      inline void DoAssignTov(ComplexType(T)* vp, 
	  const int S, const bool C) const
      {
	if (C) AddVV_x_m1(CONJ(x1),v1.Conjugate(),v2.Conjugate(),vp,S); 
	else AddVV_x_m1(x1,v1,v2,vp,S); 
      }
    private:
      const T x1;
      const GenSmallVector<T1,N,S1,C1>& v1;
      const GenSmallVector<T2,N,S2,C2>& v2;
  };

  template <class T, class T1, class T2, size_t N, int S1, bool C1, int S2, bool C2> class Sumvv : 
    public SmallVectorComposite<T,N> 
  {
    public:
      inline Sumvv(T _x1, const GenSmallVector<T1,N,S1,C1>& _v1, 
	  T _x2, const GenSmallVector<T2,N,S2,C2>& _v2) :
	x1(_x1),v1(_v1),x2(_x2), v2(_v2)
      { TMVAssert(v1.size() == v2.size()); }
      inline size_t size() const { return v1.size(); }
      inline T GetX1() const { return x1; }
      inline const GenSmallVector<T1,N,S1,C1>& GetV1() const { return v1; }
      inline T GetX2() const { return x2; }
      inline const GenSmallVector<T2,N,S2,C2>& GetV2() const { return v2; }
      inline void AssignTov(
	  const SmallVectorView<RealType(T),N,1,false>& v0) const
      { TMVAssert(IsReal(T())); AddVV(x1,v1,x2,v2,v0); }
      inline void AssignTov(
	  const SmallVectorView<ComplexType(T),N,1,false>& v0) const
      { AddVV(x1,v1,x2,v2,v0); }
      inline void AssignTov(
	  const SmallVectorView<ComplexType(T),N,1,true>& v0) const
      { AddVV(x1,v1,x2,v2,v0); }
      inline void DoAssignTov(RealType(T)* vp, const int S) const
      { TMVAssert(IsReal(T())); AddVV(x1,v1,x2,v2,vp,S); }
      inline void DoAssignTov(ComplexType(T)* vp, 
	  const int S, const bool C) const
      {
	if (C) AddVV(CONJ(x1),v1.Conjugate(),CONJ(x2),v2.Conjugate(),vp,S); 
	else AddVV(x1,v1,x2,v2,vp,S); 
      }
    private:
      const T x1;
      const GenSmallVector<T1,N,S1,C1>& v1;
      const T x2;
      const GenSmallVector<T2,N,S2,C2>& v2;
  };

  // v+=v
  template <class T, size_t N, int S1, bool C1, int S2, bool C2> 
    inline const SmallVectorView<T,N,S1,C1>& operator+=(
	const SmallVectorView<T,N,S1,C1>& v1, 
	const GenSmallVector<T,N,S2,C2>& v2) 
    { AddVV_1(v2,v1); return v1; }

  template <class T, size_t N, int S1, bool C1, int S2, bool C2> 
    inline const SmallVectorView<CT,N,S1,C1>& operator+=(
	const SmallVectorView<CT,N,S1,C1>& v1, 
	const GenSmallVector<T,N,S2,C2>& v2) 
    { AddVV_1(v2,v1); return v1; }

  // v-=v
  template <class T, size_t N, int S1, bool C1, int S2, bool C2> 
    inline const SmallVectorView<T,N,S1,C1>& operator-=(
	const SmallVectorView<T,N,S1,C1>& v1,
	const GenSmallVector<T,N,S2,C2>& v2)
    { AddVV_m1(v2,v1); return v1; }

  template <class T, size_t N, int S1, bool C1, int S2, bool C2> 
    inline const SmallVectorView<CT,N,S1,C1>& operator-=(
	const SmallVectorView<CT,N,S1,C1>& v1, 
	const GenSmallVector<T,N,S2,C2>& v2) 
    { AddVV_m1(v2,v1); return v1; }

  // v+=(x*v)
  template <class T, class T2, size_t N, int S1, bool C1, int S2, bool C2> 
    inline const SmallVectorView<T,N,S1,C1>& operator+=(
	const SmallVectorView<T,N,S1,C1>& v, const ProdXv<T,T2,N,S2,C2>& pxv)
    { AddVV(pxv.GetX(),pxv.GetV(),v); return v; }

  template <class T, size_t N, int S1, bool C1, int S2, bool C2> 
    inline const SmallVectorView<CT,N,S1,C1>& operator+=(
	const SmallVectorView<CT,N,S1,C1>& v, const ProdXv<T,T,N,S2,C2>& pxv)
    { AddVV(pxv.GetX(),pxv.GetV(),v); return v; }

  // v-=(x*v)
  template <class T, class T2, size_t N, int S1, bool C1, int S2, bool C2> 
    inline const SmallVectorView<T,N,S1,C1>& operator-=(
	const SmallVectorView<T,N,S1,C1>& v, const ProdXv<T,T2,N,S2,C2>& pxv)
    { AddVV(-pxv.GetX(),pxv.GetV(),v); return v; }

  template <class T, size_t N, int S1, bool C1, int S2, bool C2> 
    inline const SmallVectorView<CT,N,S1,C1>& operator-=(
	const SmallVectorView<CT,N,S1,C1>& v, const ProdXv<T,T,N,S2,C2>& pxv)
    { AddVV(-pxv.GetX(),pxv.GetV(),v); return v; }

  // Mix with Vector
  template <class T, size_t N, int S1, bool C1> 
    inline const SmallVectorView<T,N,S1,C1>& operator+=(
	const SmallVectorView<T,N,S1,C1>& v1, const GenVector<T>& v2) 
    {
      TMVAssert(v2.size() == N);
      AddVV(T(1),v2,v1.RegView()); 
      return v1; 
    }

  template <class T, size_t N, int S1, bool C1> 
    inline const SmallVectorView<CT,N,S1,C1>& operator+=(
	const SmallVectorView<CT,N,S1,C1>& v1, const GenVector<T>& v2) 
    {
      TMVAssert(v2.size() == N);
      AddVV(T(1),v2,v1.RegView()); 
      return v1; 
    }

  template <class T, size_t N, int S1, bool C1> 
    inline const SmallVectorView<T,N,S1,C1>& operator-=(
	const SmallVectorView<T,N,S1,C1>& v1, const GenVector<T>& v2) 
    {
      TMVAssert(v2.size() == N);
      AddVV(T(-1),v2,v1.RegView()); 
      return v1; 
    }

  template <class T, size_t N, int S1, bool C1> 
    inline const SmallVectorView<CT,N,S1,C1>& operator-=(
	const SmallVectorView<CT,N,S1,C1>& v1, const GenVector<T>& v2) 
    {
      TMVAssert(v2.size() == N);
      AddVV(T(-1),v2,v1.RegView()); 
      return v1; 
    }

  template <class T, class Tv> class ProdXV;

  template <class T, class T2, size_t N, int S1, bool C1> 
    inline const SmallVectorView<T,N,S1,C1>& operator+=(
	const SmallVectorView<T,N,S1,C1>& v, const ProdXV<T,T2>& pxv) 
    {
      TMVAssert(pxv.size() == N);
      AddVV(pxv.GetX(),pxv.GetV(),v.RegView()); 
      return v; 
    }

  template <class T, size_t N, int S1, bool C1> 
    inline const SmallVectorView<CT,N,S1,C1>& operator+=(
	const SmallVectorView<CT,N,S1,C1>& v, const ProdXV<T,T>& pxv) 
    {
      TMVAssert(pxv.size() == N);
      AddVV(pxv.GetX(),pxv.GetV(),v.RegView()); 
      return v; 
    }

  template <class T, class T2, size_t N, int S1, bool C1> 
    inline const SmallVectorView<T,N,S1,C1>& operator-=(
	const SmallVectorView<T,N,S1,C1>& v, const ProdXV<T,T2>& pxv) 
    {
      TMVAssert(pxv.size() == N);
      AddVV(-pxv.GetX(),pxv.GetV(),v.RegView()); 
      return v; 
    }

  template <class T, size_t N, int S1, bool C1> 
    inline const SmallVectorView<CT,N,S1,C1>& operator-=(
	const SmallVectorView<CT,N,S1,C1>& v, const ProdXV<T,T>& pxv) 
    {
      TMVAssert(pxv.size() == N);
      AddVV(-pxv.GetX(),pxv.GetV(),v.RegView()); 
      return v; 
    }

#define GENMATRIX1 GenSmallVector
#define GENMATRIX2 GenSmallVector
#define SUMMM Sumvv
#define SUMMM_1_1 Sumvv_1_1
#define SUMMM_1_m1 Sumvv_1_m1
#define SUMMM_1_x Sumvv_1_x
#define SUMMM_x_1 Sumvv_x_1
#define SUMMM_x_m1 Sumvv_x_m1
#define PRODXM1 ProdXv
#define PRODXM2 ProdXv
#define X1 ,N,S1,C1
#define X2 ,N,S2,C2
#define X3 ,N,S1,C1,S2,C2
#define Y ,size_t N,int S1,bool C1,int S2,bool C2
#include "TMV_AuxSumMM.h"

#define SUMMM_1_1 Sumvv_1_1
#define SUMMM_1_m1 Sumvv_1_m1
#define SUMMM_1_x Sumvv_1_x
#define SUMMM_x_1 Sumvv_x_1
#define SUMMM_x_m1 Sumvv_x_m1
#define X3 ,N,S1,C1,S2,C2
#define Y ,size_t N,int S1,bool C1,int S2,bool C2
#include "TMV_AuxSumMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef SUMMM
#undef PRODXM1
#undef PRODXM2

  //
  // Vector * Vector
  //

  template <class T, size_t N, int S1, bool C1, int S2, bool C2> 
    inline T operator*(
      const GenSmallVector<T,N,S1,C1>& v1, const GenSmallVector<T,N,S2,C2>& v2) 
  { return MultVV(v1,v2); }

  template <class T, size_t N, int S1, bool C1, int S2, bool C2> 
    inline CT operator*(
      const GenSmallVector<CT,N,S1,C1>& v1, const GenSmallVector<T,N,S2,C2>& v2)
  { return MultVV(v1,v2); }

  template <class T, size_t N, int S1, bool C1, int S2, bool C2> 
    inline CT operator*(
      const GenSmallVector<T,N,S1,C1>& v1, const GenSmallVector<CT,N,S2,C2>& v2)
  { return MultVV(v2,v1); }

  // v * (x*v)

  template <class T, class T2, size_t N, int S1, bool C1, int S2, bool C2> 
    inline T operator*(
      const GenSmallVector<T,N,S1,C1>& v1, const ProdXv<T,T2,N,S2,C2>& v2) 
  { return v2.GetX()*MultVV(v1,v2.GetV()); }

  template <class T, size_t N, int S1, bool C1, int S2, bool C2> 
    inline CT operator*(
      const GenSmallVector<CT,N,S1,C1>& v1, const ProdXv<T,T,N,S2,C2>& v2)
  { return v2.GetX()*MultVV(v1,v2.GetV()); }

  template <class T, class T2, size_t N, int S1, bool C1, int S2, bool C2> 
    inline CT operator*(
      const GenSmallVector<T,N,S1,C1>& v1, const ProdXv<CT,T2,N,S2,C2>& v2)
  { return v2.GetX()*MultVV(v2.GetV(),v1); }

  // (x*v) * v

  template <class T, class T1, size_t N, int S1, bool C1, int S2, bool C2> 
    inline T operator*(
      const ProdXv<T,T1,N,S1,C1>& v1, const GenSmallVector<T,N,S2,C2>& v2)
  { return v1.GetX()*MultVV(v1.GetV(),v2); }

  template <class T, size_t N, int S1, bool C1, int S2, bool C2> 
    inline CT operator*(
      const ProdXv<T,T,N,S1,C1>& v1, const GenSmallVector<CT,N,S2,C2>& v2)
  { return v1.GetX()*MultVV(v2,v1.GetV()); }

  template <class T, class T1, size_t N, int S1, bool C1, int S2, bool C2> 
    inline CT operator*(
      const ProdXv<CT,T1,N,S1,C1>& v1, const GenSmallVector<T,N,S2,C2>& v2)
  { return v1.GetX()*MultVV(v1.GetV(),v2); }

  // (x*v) * (x*v)

  template <class T, class T1, class T2, size_t N, int S1, bool C1, int S2, bool C2> 
    inline T operator*(
      const ProdXv<T,T1,N,S1,C1>& v1, const ProdXv<T,T2,N,S2,C2>& v2)
  { return v1.GetX()*v2.GetX()*MultVV(v1.GetV(),v2.GetV()); }

  template <class T, class T1, size_t N, int S1, bool C1, int S2, bool C2> 
    inline CT operator*(
      const ProdXv<CT,T1,N,S1,C1>& v1, const ProdXv<T,T,N,S2,C2>& v2)
  { return v1.GetX()*v2.GetX()*MultVV(v1.GetV(),v2.GetV()); }

  template <class T, class T2, size_t N, int S1, bool C1, int S2, bool C2> 
    inline CT operator*(
      const ProdXv<T,T,N,S1,C1>& v1, const ProdXv<CT,T2,N,S2,C2>& v2)
  { return v1.GetX()*v2.GetX()*MultVV(v1.GetV(),v2.GetV()); }

} // namespace tmv

#undef CT
#undef CCT
#undef VCT

#endif 
